"""Multi-tone BCC: N parallel BFSK (or OOK) tones over a measured frequency-selective body channel.

Purpose: each tone carries its own stream of HBC packets (HBCProtocol.PacketGenerator, spreading off, 13125kbps,
pilot_info '010', burst '0', 255-byte PSDUs of different text per tone, 1000 idle chips before and after each packet as
in an earlier script (v2_packet_run); the 64-chip SFD trailer that the reference chip file appends to each
packet is omitted, so CPB = 8800 / 2040 as in) on BFSK tones at f_i -/+ R_c
(PHY_layer_BCC_simulation_v2.PHY2 m = 1: tone separation 2 R_c, modulation index h = 2). Each tone is received by the
PHY2 chain (band-pass 4 R_c centred on f_i, non-coherent energy detection per chip at f_i -/+ R_c) and scored with
the repository PacketParser exactly as that script's score() does (detected, decoded, error-free; PSR = error-free / transmitted). OOK variant: tones
2 R_c (1 + guard) apart, band-pass 2 R_c, energy at f_i against the PHY2 adaptive threshold.
Channel: per-tone gains at f_i -/+ R_c (OOK: at f_i) interpolated in log10(f) from a curve of channel_curves.csv
(PHY2 gain_db_tones semantics). Channel phase and group delay are NOT modelled (non-coherent detection), except in
'isi', which filters the wave with the minimum-phase response of the measured magnitude.
Noise: one-sided N0 = 4 k T R 10^(NF/10), T = 290 K, R = 50 Ohm (8.0e-19 V^2/Hz before NF)
for every load, the 1 MOhm curves included (not a like-for-like noise reference for a high-impedance input),
or N0 set so that the best tone reaches an Eb/N0 target. Eb = mean chip energy a^2 / (2 R_c) (PHY2 convention;
BFSK: mean over the two tones, OOK: a 'one' chip). Drive: per-tone amplitude A, split none | power (A / sqrt(N)) |
peak (A / N). Tones: spacing 4 R_c (1 + guard), centred in [f_lo, f_hi] (default: from the lowest to the highest
edge of the 3 dB usable bins of the curve, so bins outside the tolerance inside that span also carry tones).
Async (default): random fractional-chip timing offset and random packet start offset per tone; --sync aligns both.
Transmitter phase: each tone is continuous-phase FSK (carrier f_i runs continuously, the -/+ R_c deviation restarts
every chip; with h = 2 the phase is continuous at every chip boundary); OOK gates a continuous carrier.
This equals PHY2, whose tone template restarts every chip, whenever f_i / R_c is an integer (every PHY2 validation);
PHY2's restart at an arbitrary tone centre would add a phase jump at every chip boundary, an artefact.
Simulation: complex-baseband equivalent per victim tone with its M nearest neighbours on each side (2; OOK 4),
an even integer number of samples per chip (>= 32, raised so neighbours do not alias into the receive band), the
PHY2 FIR band-pass as its low-pass equivalent, complex noise of variance N0 fs per component (the envelope of
real noise of one-sided density N0). 'check' proves the equivalence with PHY2 and a direct passband simulation.
A single tone whose receive band f -/+ 2 R_c reaches 0 Hz is simulated in real passband instead: there the negative-frequency image of the real signal, which the baseband model lacks, enters the band.
Aggregates: payload per tone R_c / CPB, i.e. back-to-back packets (idle chips excluded, as in the
estimate; the simulated streams carry 2000 idle chips per 8800-chip packet, duty 0.815). In-band columns count the
tones whose receive band f_i -/+ 2 R_c (OOK -/+ R_c) lies inside the 3 dB bins of the curve.
Inputs:  channel_curves.csv, ../HBCProtocol.py, ../PHY_layer_BCC_simulation_v2.py
Outputs: run   -> per-tone table and summary on stdout
         sweep -> results/bcc_fdma_tones.csv (one row per configuration and tone), results/bcc_fdma_summary.csv
         check -> results/bcc_fdma_equivalence.csv (and a markdown table on stdout)
         rates -> results/bcc_rate_stability.csv (one BFSK tone at channel_curves.tone_centre, 50 Ohm, 1 V
                  at NF 10 / 20 / 30 and the Eb/N0 target KNEE_DB, every measured path and the cohort curve, per chip
                  rate of RATES; psr_model is the engine's frame model, phy.frame_outcome of the as-built 255-octet
                  frame at the one chip Ec/N0 of channel.bfsk_ec_n0_db) and its summary table on stdout
         isi   -> results/bcc_rate_isi.csv and its summary table on stdout
argv:    run [--curve C] [--load L] [--rc R] [--guard G] [--nf NF | --ebn0 E] [--drive V] [--split S] [--mod M]
             [--sync] [--packets P] [--neighbours M] [--band LO HI] [--n-tones N] [--seed S] [--workers W]
         sweep [--workers W] [--grid-packets P] [--headline-packets P] [--aci-packets P]
         check [--chips N] [--workers W]
         rates [--packets P] [--chunks C] [--workers W]
         isi [--packets P] [--chunks C] [--workers 1|2]
"""
from __future__ import annotations
import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")                    # one BLAS thread per worker process
import argparse
import contextlib
import csv
import functools
import io
import sys
import time
from dataclasses import dataclass
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.signal import fftconvolve, firwin

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import HBCProtocol as protocol                        # noqa: E402
import PHY_layer_BCC_simulation_v2 as v2              # noqa: E402
import channel_curves as cc                           # noqa: E402

RESULTS = HERE / "results"
N0_THERMAL = 4 * 1.380649e-23 * 290.0 * 50.0          # 8.0e-19 V^2/Hz, thermal noise of 50 Ohm at 290 K
MAX_WORKERS = 8
M_DEFAULT = {"fsk": 2, "ook": 4}                      # neighbours per side; OOK at 2 R_c needs 4 (check (c))
IDLE, PSDU_BYTES, F_FLAT = 1000, 255, 10e6
ALPHABET = np.array(list("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ ,."))
PARSER = protocol.PacketParser(protocol.PacketGenerator(spreading_enabled=False))
PLCP = np.array(protocol.PacketGenerator.PLCP_CODE) * 2 - 1
SFD = np.array(protocol.PacketGenerator.SFD_CODE) * 2 - 1


def text_bits(s):
    return "".join(format(ord(c), "08b") for c in s)


def make_packet(gen, text):
    b = text_bits(text)
    return np.array(gen.generate_packet(data=b, data_rate="13125kbps", pilot_info="010", burst_mode="0",
                                        psdu_length=len(b) // 8), dtype=np.uint8)


CHIPS_PER_PACKET = len(make_packet(protocol.PacketGenerator(spreading_enabled=False), "x" * PSDU_BYTES))
CPB = CHIPS_PER_PACKET / (PSDU_BYTES * 8)            # chips per payload bit (8800 / 2040 = 4.31)
PERIOD = CHIPS_PER_PACKET + 2 * IDLE                  # idle, packet, idle


@functools.lru_cache(maxsize=64)
def tone_stream(seed, tone, n_packets, sync):
    """(bits uint8, texts) of one tone: [random idle start if async] + n x (idle, packet, idle) + idle pad.
    Bit 0 = lower BFSK tone (chip -1, also sent while idle) or OOK off; bit 1 = upper tone or OOK on."""
    rng = np.random.default_rng([seed, tone, 1])
    gen = protocol.PacketGenerator(spreading_enabled=False)
    texts = tuple("".join(rng.choice(ALPHABET, PSDU_BYTES)) for _ in range(n_packets))
    off = 0 if sync else int(rng.integers(PERIOD))
    bits = np.zeros(n_packets * PERIOD + (0 if sync else PERIOD), np.uint8)
    for i, t in enumerate(texts):
        s = off + i * PERIOD + IDLE
        bits[s:s + CHIPS_PER_PACKET] = make_packet(gen, t)
    bits.flags.writeable = False
    return bits, texts


def score(rx_bits, texts):
    """an earlier script (v2_packet_run) score(): (detected, decoded, error-free)."""
    stream = rx_bits.astype(int) * 2 - 1
    with contextlib.redirect_stdout(io.StringIO()):
        pk = PARSER.extract_packets(stream, PLCP, SFD, 50, 50)
    bers = []
    for p in pk:
        if len(p) < 1000:
            continue
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                data, *_, crc_ok, rawc = PARSER.unpack_packet(((p + 1) // 2).astype(int))
            b = text_bits(protocol.find_best_match(PARSER.bin_to_str(data), texts)); n = min(len(b), len(data))
            bers.append(min(1.0, (sum(x != y for x, y in zip(b[:n], data[:n])) + abs(len(b) - len(data))) / len(b)))
        except Exception:
            pass
    return len(pk), len(bers), int(sum(b == 0 for b in bers))


# ---------------------------------------------------------------- geometry
@dataclass(frozen=True)
class Geometry:
    series: str
    curve: str                    # 'cohort', a path 'TX-RX', or 'flat' (0 dB at every frequency)
    load: str                     # '50' or '1M' ('' for flat)
    Rc: float
    guard: float                  # spacing = 4 R_c (1 + guard) for BFSK, 2 R_c (1 + guard) for OOK
    mod: str = "fsk"
    sync: bool = False
    n_packets: int = 6
    M: int = 2                    # neighbours simulated on each side of a victim
    seed: int = 1
    band: tuple | None = None     # (f_lo, f_hi) Hz; None = 3 dB usable band of the curve
    n_tones: int | None = None    # None = as many as fit the band
    spc: int | None = None        # samples per chip; None = automatic


@dataclass
class Plan:
    f: np.ndarray                 # tone centres (Hz)
    gdb: np.ndarray               # (n, 2) gain dB at f - Rc, f + Rc (OOK: at f twice)
    amps: np.ndarray              # (n, 2) received amplitude per 1 V drive for bit 0 / bit 1
    e_chip: np.ndarray            # mean chip energy per 1 V drive, a^2 / (2 Rc)
    lo: float
    hi: float
    spc: int
    in_band: int | None           # tones whose receive band lies inside the 3 dB bins (None for 'flat')


@functools.lru_cache(maxsize=None)
def plan(g: Geometry) -> Plan:
    spacing = (4.0 if g.mod == "fsk" else 2.0) * g.Rc * (1.0 + g.guard)
    if g.curve == "flat":
        n = g.n_tones or 1
        lo, hi = F_FLAT - n * spacing / 2, F_FLAT + n * spacing / 2
        gain, bins = np.zeros_like, ()
    else:
        fr, gd = cc.read_curve(g.curve, g.load)
        ub = cc.usable_band(fr, gd)
        lo, hi = g.band if g.band else ub[2:4]
        bins = ub[4]
        n = g.n_tones or int(np.floor((hi - lo) / spacing + 1e-9))
        gain = lambda x: np.interp(np.log10(x), np.log10(fr), gd)          # noqa: E731
    f = (lo + hi) / 2 + (np.arange(n) - (n - 1) / 2) * spacing
    if g.mod == "fsk":
        gdb = np.stack([gain(f - g.Rc), gain(f + g.Rc)], axis=1); amps = 10 ** (gdb / 20)
        e = (amps ** 2).mean(axis=1) / (2 * g.Rc)
    else:
        gdb = np.stack([gain(f)] * 2, axis=1); amps = 10 ** (gdb / 20); amps[:, 0] = 0.0
        e = amps[:, 1] ** 2 / (2 * g.Rc)
    fmax = min(g.M, n - 1) * spacing + g.Rc            # farthest simulated component from a victim centre
    spc = g.spc or max(32, 2 * int(np.ceil(fmax / g.Rc)) + 8)
    half = (2.0 if g.mod == "fsk" else 1.0) * g.Rc     # half the receive band
    in_band = int(sum(any(a - 1e-3 <= x - half and x + half <= b + 1e-3 for a, b in bins) for x in f)) if bins else None
    return Plan(f, gdb, amps, e, float(lo), float(hi), spc + spc % 2, in_band)


def offsets(g: Geometry, n, spc):
    """Timing offset of each tone in samples (fraction of a chip x spc)."""
    return [0 if g.sync else int(np.random.default_rng([g.seed, j, 3]).integers(spc)) for j in range(n)]


# ---------------------------------------------------------------- complex-baseband core
def rx_kernel(Rc, mod, spc):
    """Receive filter (PHY2 band-pass as its low-pass equivalent) and chip correlators folded into one kernel,
    returned as (m, spc, n_templates) plus the filter half length h (samples)."""
    fs = spc * Rc
    bw = 4.0 * Rc if mod == "fsk" else 2.0 * Rc
    ntaps = int(4 * fs / bw) | 1
    taps = firwin(ntaps, bw / 2, fs=fs)                 # symmetric, unit gain at DC (= at f_i in passband)
    r = np.arange(spc)
    tpl = [np.exp(-2j * np.pi * s * Rc * r / fs) for s in ((-1, 1) if mod == "fsk" else (0,))]
    Q = np.stack([np.convolve(t, taps) for t in tpl], axis=1)       # length spc + ntaps - 1
    m = -(-Q.shape[0] // spc)
    Q = np.vstack([Q, np.zeros((m * spc - Q.shape[0], Q.shape[1]))])
    return Q.reshape(m, spc, -1), (ntaps - 1) // 2


def victim_corr(v, f, amps, bits, d, Rc, mod, spc, M, noise_seed, block=4096):
    """Chip correlator outputs of victim tone v: signal for 1 V drive and noise for N0 = 1 V^2/Hz, each (L, n_tpl).
    f, amps: tone centres and per-bit amplitudes; bits[j]: chip bits of tone j (only |j - v| <= M are read);
    d[j]: timing offset of tone j in samples. Chip k of tone j starts at sample k spc + d[j] = t0 fs; its passband
    wave is a cos(2 pi f_j t + 2 pi s R_c (t - t0)) (CPFSK, s = -/+1; OOK: s = 0 and a = 0 for an 'off' chip), whose
    envelope about f_v is a exp(j 2 pi (f_j - f_v) t0) exp(j 2 pi (f_j + s R_c - f_v)(t - t0)). When f_j / R_c is an
    integer this equals PHY2's template restarted every chip, a cos(2 pi (f_j + s R_c)(t - t0)), up to one constant
    phase per tone."""
    n, L, dv = len(f), len(bits[v]), d[v]
    fs = spc * Rc
    Q, h = rx_kernel(Rc, mod, spc)
    m = Q.shape[0]
    dev = np.array((-Rc, Rc) if mod == "fsk" else (0.0, 0.0))
    r = np.arange(spc)
    J = range(max(0, v - M), min(n, v + M + 1))
    rows = {j: np.exp(2j * np.pi * np.outer(f[j] - f[v] + dev, r) / fs) for j in J}
    Cs = np.empty((L, Q.shape[2]), complex); Cn = np.empty_like(Cs)
    rng = np.random.default_rng(noise_seed)
    sd = np.sqrt(fs)                                     # per-component sd for N0 = 1: sqrt(N0 fs)
    ext = (m - 1) * spc
    tail = (rng.standard_normal(ext) + 1j * rng.standard_normal(ext)) * sd
    for k0 in range(0, L, block):
        B = min(block, L - k0)
        n0, nlen = k0 * spc + dv - h, (B + m - 1) * spc
        z = np.zeros(nlen, complex)
        for j in J:
            klo, khi = (n0 - d[j]) // spc, (n0 + nlen - 1 - d[j]) // spc + 1
            kk = np.arange(klo, khi)
            b = bits[j][np.clip(kk, 0, L - 1)]
            cyc = np.mod((f[j] - f[v]) / Rc * (kk + d[j] / spc), 1.0)              # (f_j - f_v) t0 in cycles
            c = np.where((kk >= 0) & (kk < L), amps[j][b], 0.0) * np.exp(2j * np.pi * cyc)
            s = n0 - (klo * spc + d[j])
            z += (rows[j][b] * c[:, None]).ravel()[s:s + nlen]
        Z = z.reshape(-1, spc)
        Cs[k0:k0 + B] = sum(Z[i:i + B] @ Q[i] for i in range(m))
        zn = np.concatenate([tail, (rng.standard_normal(B * spc) + 1j * rng.standard_normal(B * spc)) * sd])
        tail = zn[len(zn) - ext:]
        Zn = zn.reshape(-1, spc)
        Cn[k0:k0 + B] = sum(Zn[i:i + B] @ Q[i] for i in range(m))
    return Cs, Cn


def passband_corr(v, f, amps, bits, d, Rc, mod, spc, M, noise_seed, block=4096, chans=()):
    """victim_corr in real passband, one BFSK tone. The complex-baseband model has no negative-frequency
    image, so it fails where the receive band f -/+ 2 R_c reaches 0 Hz (run_victim switches there): a lower tone at
    f - R_c << R_c makes a real chip a slowly turning DC level whose energy fades with the carrier phase. Transmitter:
    the continuous-phase chips of passband_errors (phy2_tx False), chip k from sample k spc + d; receiver: the PHY2
    band-pass [max(f - 2 R_c, 1 Hz), f + 2 R_c] (its taps) and correlators at f -/+ R_c over the victim's chips; real
    noise of one-sided density 1 V^2/Hz. Outputs as victim_corr: signal for 1 V drive, noise for N0 = 1.
    chans: (channel impulse response at fs, receiver lag in samples[, amps]) entries; then the signal output
    is a list, one Cs per entry, the transmitted wave (tone amplitudes amps unless the entry has its own) through the
    channel and the band-pass, chip k read from sample k spc + lag (the receiver's chip timing follows the channel's
    delay; noise enters after the channel)."""
    assert len(f) == 1 and mod == "fsk", "passband_corr simulates one BFSK tone"
    f, fs, bits, d = f[0], spc * Rc, bits[v], d[v]
    taps = firwin(spc | 1, [max(f - 2 * Rc, 1.0), min(f + 2 * Rc, fs / 2 - 1.0)], pass_zero=False, fs=fs)
    h = (len(taps) - 1) // 2
    filt = [(np.convolve(c[0], taps), c[1], c[2] if len(c) > 2 else amps) for c in chans] or [(taps, 0, amps)]
    tpl = np.exp(-2j * np.pi * np.outer(np.arange(spc), (f - Rc, f + Rc)) / fs)
    L = len(bits)
    rng = np.random.default_rng(noise_seed)
    sd = np.sqrt(fs / 2)
    tail = rng.standard_normal(2 * h) * sd
    Cs, Cn = [np.empty((L, 2), complex) for _ in filt], np.empty((L, 2), complex)
    for k0 in range(0, L, block):
        B = min(block, L - k0)
        s0 = [k0 * spc + lag + h - (len(g) - 1) for g, lag, _ in filt]     # first sample each channel reads
        n = np.arange(min(s0), (k0 + B) * spc + max(lag for _, lag, _ in filt) + h)     # samples from chip 0
        k, r = n // spc, n % spc
        b = bits[np.clip(k, 0, L - 1)]
        wave = np.cos(2 * np.pi * (f * (n + d) + (2.0 * b - 1) * Rc * r) / fs) * ((k >= 0) & (k < L))
        for (g, lag, am), C, i in zip(filt, Cs, s0):
            j = slice(i - n[0], i - n[0] + B * spc + len(g) - 1)
            C[k0:k0 + B] = fftconvolve(am[0][b[j]] * wave[j], g, mode="valid").reshape(B, spc) @ tpl
        zn = np.concatenate([tail, rng.standard_normal(B * spc) * sd])
        tail = zn[len(zn) - 2 * h:]
        Cn[k0:k0 + B] = fftconvolve(zn, taps, mode="valid").reshape(B, spc) @ tpl
    return (Cs if chans else Cs[0]), Cn


# ---------------------------------------------------------------- minimum-phase channel
MP_DF, MP_HALF = 1e3, 2 ** 18         # dense grid: 0 to 262.144 MHz in 1 kHz steps (cepstrum FFT of 2^19 points)
ARTEFACT_BINS = (8e6, 9e6, 12.5e6, 15e6, 20e6)   # measurement-arrangement effects , not body
EXTENSIONS = {"A": ("slope", "hold"), "B": ("hold", "slope")}     # (below 0.05 MHz, above 20 MHz)


def minphase(db):
    """Minimum-phase response exp(ln|H| + j phi), phi = -Hilbert{ln|H|}, of a magnitude db (20 log10 |H|) given on the
    uniform grid k MP_DF, k = 0 .. len(db) - 1 (the last point is the Nyquist frequency of the grid): the real cepstrum
    of ln|H| on the full circle, folded onto n >= 0 (numpy's FFT sign: H = sum h[n] exp(-j w n), so a causal cepstrum
    gives a causal, minimum-phase h and a phase lag, phi < 0 for a low-pass)."""
    n = 2 * (len(db) - 1)
    c = np.fft.irfft(np.asarray(db) * np.log(10) / 20, n)
    c[1:n // 2] *= 2
    c[n // 2 + 1:] = 0
    return np.exp(np.fft.rfft(c))


def curve_db(curve, variant="measured", ext="A", f=None):
    """20 log10 |H| of a channel_curves.csv 50 Ohm curve at f (default: the minphase grid), a monotone cubic (PCHIP)
    in log10 f through the measured points: C1, so the group delay stays finite at the points, where a line in log f
    would put a log singularity. variant 'smoothed' removes the arrangement artefacts: 8 and 9 MHz on the log-f line
    between 7 and 10 MHz (the notch), 12.5 to 20 MHz held at the 10 MHz value (the peak and the 20 MHz drop).
    ext (EXTENSIONS), each end: 'slope' continues the end segment's dB per decade for one decade (to 5 kHz, to
    200 MHz), 'hold' does not; beyond, the curve is held (a repeated knot a factor 2 further out makes the PCHIP slope
    0 there, so the hold is C1 too); 0 Hz takes the lowest-frequency value."""
    fr, gd = cc.read_curve(curve, "50")
    lf, gd = np.log10(fr), gd.copy()
    if variant == "smoothed":
        i7, i10 = (int(np.flatnonzero(fr == x)[0]) for x in (7e6, 10e6))
        gd[i7 + 1:i10] = np.interp(lf[i7 + 1:i10], lf[[i7, i10]], gd[[i7, i10]])
        gd[i10 + 1:] = gd[i10]
    lo, hi = EXTENSIONS[ext]
    if lo == "slope":
        lf, gd = np.r_[lf[0] - 1, lf], np.r_[gd[0] - (gd[1] - gd[0]) / (lf[1] - lf[0]), gd]
    if hi == "slope":
        lf, gd = np.r_[lf, lf[-1] + 1], np.r_[gd, gd[-1] + (gd[-1] - gd[-2]) / (lf[-1] - lf[-2])]
    lf, gd = np.r_[lf[0] - np.log10(2), lf, lf[-1] + np.log10(2)], np.r_[gd[0], gd, gd[-1]]
    f = np.arange(MP_HALF + 1) * MP_DF if f is None else np.asarray(f, float)
    return PchipInterpolator(lf, gd)(np.clip(np.log10(np.maximum(f, 1.0)), lf[0], lf[-1]))


@functools.lru_cache(maxsize=None)
def channel_response(curve, variant="measured", ext="A"):
    """(f, H) on the minphase grid: the minimum-phase response of a 50 Ohm curve (curve_db)."""
    return np.arange(MP_HALF + 1) * MP_DF, minphase(curve_db(curve, variant, ext))


def group_delay(H):
    """-d phi / d omega (s) of H on the minphase grid."""
    return -np.gradient(np.unwrap(np.angle(H)), 2 * np.pi * MP_DF)


def channel_fir(curve, variant, ext, fs, n=2 ** 16, keep=1 - 1e-8):
    """Minimum-phase impulse response at sample rate fs (causal taps): the cepstrum of curve_db on the n-point FFT
    grid of fs itself, so the discrete filter is minimum phase with |H| exact on 0 to fs / 2 (interpolating the
    continuous-time phase of channel_response onto a grid that stops at fs / 2 would leave a precursor), truncated
    after the sample from which the tail holds 1 - keep of the energy (-80 dB)."""
    h = np.fft.irfft(minphase(curve_db(curve, variant, ext, np.arange(n // 2 + 1) * fs / n)), n)
    e = np.cumsum(h ** 2)
    return h[:int(np.searchsorted(e, keep * e[-1])) + 1]


def smear(curve, variant, ext="A", fs=40e6):
    """(rms delay spread, 90 % energy duration) in ns of the minimum-phase impulse response over the measured band
    (channel_fir at fs = 40 MHz: 0 to 20 MHz): the least smear any channel with this magnitude can have, since
    minimum phase has the least energy delay (its partial energy sum over the first n taps is the largest for every
    n). The 90 % duration runs from the 5 % to the 95 % point of the energy of the taps interpolated 8 times
    (zero-padded spectrum, read circularly: the interpolation puts a small precursor before t = 0); the rms is taken
    on the taps themselves (the interpolation rings at the truncated end, which the t^2 weight would amplify). Both
    sit near the time resolution of a 20 MHz band, 1 / 40 MHz = 25 ns."""
    h = channel_fir(curve, variant, ext, fs)
    e, t = h ** 2 / np.sum(h ** 2), np.arange(len(h)) / fs * 1e9
    rms = float(np.sqrt(np.sum(t ** 2 * e) - np.sum(t * e) ** 2))
    u = np.roll(np.fft.irfft(np.fft.rfft(h, 4 * len(h)), 32 * len(h)) ** 2, 16 * len(h))
    c = np.cumsum(u) / np.sum(u)
    return rms, float((np.searchsorted(c, 0.95) - np.searchsorted(c, 0.05)) / (8 * fs) * 1e9)


def fir_delay(hc, f, fs, n=2 ** 16):
    """Group delay (samples) of the taps hc at frequency f."""
    ph = np.unwrap(np.angle(np.fft.rfft(hc, n)))
    k = int(round(f / fs * n))
    return -(ph[k + 1] - ph[k - 1]) / (2 * 2 * np.pi / n)


def isi_spc(fc, rc):
    """Samples per chip of the passband run: even, >= 32, fs >= 3 (f + 2 R_c), so what the sampled wave folds
    over fs / 2 lands in the receive band only from 2 (f + 2 R_c) up, >= 28 R_c from the tone (sidelobes near
    -60 dB)."""
    return max(32, 2 * int(np.ceil(3 * (fc + 2 * rc) / (2 * rc))))


def rx_lag(hc, amps, fc, rc, spc, n=3000, seed=7):
    """(lag, margin): the receiver's chip timing over the channel hc (samples), the lag within -/+ half a chip of the
    channel's group delay at the tone centre that maximises the smallest normalised decision margin
    (E_right - E_wrong) / (E_right + E_wrong) over n noise-free random chips (tone amplitudes amps), and that margin
    (1: no energy on the wrong tone; 0: the eye is shut). An ideal timing recovery: a pure delay costs nothing."""
    lag0 = int(round(fir_delay(hc, fc, spc * rc))) if len(hc) > 1 else 0
    bits = np.random.default_rng(seed).integers(0, 2, n).astype(np.uint8)
    lags = lag0 + np.arange(-(spc // 2), spc // 2 + 1)
    Cs, _ = passband_corr(0, np.array([fc]), amps, [bits], [0], rc, "fsk", spc, 0, seed, chans=[(hc, x) for x in lags])
    k, e = np.arange(50, n - 50), bits[50:n - 50]
    marg = [((E[k, e] - E[k, 1 - e]) / (E[k, e] + E[k, 1 - e])).min() for E in (np.abs(C) ** 2 for C in Cs)]
    return int(lags[int(np.argmax(marg))]), float(max(marg))


@functools.lru_cache(maxsize=None)
def isi_channel(curve, variant, ext, rc):
    """(channel taps, tone amplitudes (1, 2), receiver lag, margin, mean chip energy per 1 V) of one curve variant at
    the sample rate of the run at chip rate rc. variant 'tones': the model without dispersion, H = 1
    (one tap) and the curve's gains at f -/+ R_c on the tones (plan(), a line in log f); otherwise the minimum-phase
    channel of curve_db and a 1 V drive on both tones."""
    fc = cc.tone_centre(rc)
    spc = isi_spc(fc, rc)
    if variant == "tones":
        p = plan(Geometry("rates", curve, "50", rc, 0.0, band=(fc, fc), n_tones=1))
        hc, amps, e_chip = np.ones(1), p.amps, float(p.e_chip[0])
    else:
        hc, amps = channel_fir(curve, variant, ext, spc * rc), np.ones((1, 2))
        e_chip = float(np.mean(10 ** (curve_db(curve, variant, ext, [fc - rc, fc + rc]) / 10)) / (2 * rc))
    return (hc, amps, *rx_lag(hc, amps, fc, rc, spc), e_chip)


def run_isi(task):
    """One chunk of packets of one curve at every chip rate through each channel variant -> rows per rate,
    variant and setting: NF 10 / 20 / 30 at 1 V (chip errors and packets), and the Eb/N0 target KNEE_DB on the mean
    chip energy of the two tones (chip errors only, for the dB cost of the dispersion). Every variant sees the same
    chips, carrier offset and noise. The chip stream does not depend on the rate, so the error-free score is computed
    once."""
    curve, seed, n_packets, rcs, chans, nfs = task
    bits, texts = tone_stream(seed, 0, n_packets, False)
    out, clean = [], None
    for rc in rcs:
        fc = cc.tone_centre(rc)
        spc = isi_spc(fc, rc)
        d = int(np.random.default_rng([seed, 0, 3]).integers(spc))          # as offsets()
        ch = [isi_channel(curve, v, e, rc) for v, e in chans]
        Cs, Cn = passband_corr(0, np.array([fc]), None, [bits], [d], rc, "fsk", spc, 0, [seed, 0, 2],
                               chans=[(hc, lag, amps) for hc, amps, lag, *_ in ch])
        for (v, e), (_, _, lag, marg, e_chip), C in zip(chans, ch, Cs):
            for nf in nfs + (None,):
                N0 = N0_THERMAL * 10 ** (nf / 10) if nf else e_chip / 10 ** (KNEE_DB / 10)
                rx = detect(C, Cn, 1.0, N0, "fsk")
                errs = int(np.count_nonzero(rx != bits))
                ef = None
                if nf and errs == 0:
                    clean = clean or score(bits, texts)
                    ef = clean[2]
                elif nf:
                    ef = score(rx, texts)[2]
                out.append(dict(rc_cps=rc, curve=curve, variant=v, ext=e, nf_db=nf,
                                ebn0_target_db=None if nf else KNEE_DB, lag_ns=lag / (spc * rc) * 1e9, margin=marg,
                                packets=n_packets, chips=len(rx), chip_errors=errs, error_free=ef))
    return out


def gd_spread(curve, variant, ext, rc):
    """Group delay of the minimum-phase curve (ns) at chip rate rc, tone at channel_curves.tone_centre: the spread
    (max - min) over each BFSK tone's main lobe f -/+ R_c -/+ R_c within the measured 0.05 to 20 MHz, lower and upper
    tone, counting only where |H| is within 10 dB of its value at the tone (the lobe's received energy; below 0.1 MHz
    the extension's corner puts a large delay where the channel is 25 dB down), and the delay of the upper tone minus
    that of the lower one."""
    f, H = channel_response(curve, variant, ext)
    tau, db, fc = group_delay(H) * 1e9, 20 * np.log10(np.abs(H)), cc.tone_centre(rc)

    def spread(ft):
        m = (f >= max(ft - rc, 50e3)) & (f <= min(ft + rc, 20e6)) & (db >= db[round(ft / MP_DF)] - 10)
        return float(np.ptp(tau[m]))

    return spread(fc - rc), spread(fc + rc), float(tau[round((fc + rc) / MP_DF)] - tau[round((fc - rc) / MP_DF)])


def detect(Cs, Cn, alpha, N0, mod):
    """PHY2 detectors on the correlator outputs: BFSK energy comparison, OOK adaptive threshold."""
    X = alpha * Cs + np.sqrt(N0) * Cn
    E = X.real ** 2 + X.imag ** 2
    if mod == "fsk":
        return (E[:, 1] > E[:, 0]).astype(np.uint8)
    lo, hi = np.percentile(E[:, 0], [25, 75])
    return (E[:, 0] > 0.5 * (lo + hi)).astype(np.uint8)


def resolve(s, p: Plan):
    """(per-tone drive alpha V, N0 V^2/Hz) of a setting dict."""
    n = len(p.f)
    alpha = s["drive_v"] / {"none": 1.0, "power": np.sqrt(n), "peak": float(n)}[s["split"]]
    if s["nf_db"] is not None:
        return alpha, N0_THERMAL * 10 ** (s["nf_db"] / 10)
    return alpha, alpha ** 2 * p.e_chip.max() / 10 ** (s["ebn0_target_db"] / 10)


def run_victim(task):
    """One victim tone of one geometry under several settings -> list of per-tone result rows."""
    g, v, settings = task
    p = plan(g)
    n = len(p.f)
    bits = [None] * n
    for j in range(max(0, v - g.M), min(n, v + g.M + 1)):
        bits[j] = tone_stream(g.seed, j, g.n_packets, g.sync)[0]
    texts = tone_stream(g.seed, v, g.n_packets, g.sync)[1]
    corr = victim_corr if p.f[v] > 2 * g.Rc else passband_corr     # the receive band reaches 0 Hz: real passband
    Cs, Cn = corr(v, p.f, p.amps, bits, offsets(g, n, p.spc), g.Rc, g.mod, p.spc, g.M, [g.seed, v, 2])
    clean, out = None, []
    for s in settings:
        alpha, N0 = resolve(s, p)
        rx = detect(Cs, Cn, alpha, N0, g.mod)
        errs = int(np.count_nonzero(rx != bits[v]))
        if errs == 0:                                    # same input as the transmitted stream: score once
            clean = clean or score(bits[v], texts)
        det, dec, ef = clean if errs == 0 else score(rx, texts)
        out.append(dict(series=g.series, curve=g.curve, load=g.load, Rc=g.Rc, guard=g.guard, nf_db=s["nf_db"],
                        drive_v=s["drive_v"], split=s["split"], mod=g.mod, sync=int(g.sync),
                        ebn0_target_db=s["ebn0_target_db"], seed=g.seed, fc_hz=int(round(p.f[v])), gain_m_db=p.gdb[v, 0],
                        gain_p_db=p.gdb[v, 1], ebn0_db=10 * np.log10(alpha ** 2 * p.e_chip[v] / N0),
                        chips=len(rx), chip_errors=errs, cer=errs / len(rx), packets_tx=g.n_packets,
                        detected=det, decoded=dec, error_free=ef, psr=ef / g.n_packets, _g=g))
    return out


# ---------------------------------------------------------------- tables
TONE_COLS = ["series", "curve", "load", "Rc", "guard", "nf_db", "drive_v", "split", "mod", "sync", "ebn0_target_db",
             "seed", "fc_hz", "gain_m_db", "gain_p_db", "ebn0_db", "chips", "chip_errors", "cer", "packets_tx", "detected",
             "decoded", "error_free", "psr"]
KEY = TONE_COLS[:12]
SUM_COLS = KEY + ["n_tones", "band_lo_hz", "band_hi_hz", "payload_per_tone_kbps", "nominal_aggregate_kbps",
                  "delivered_aggregate_kbps", "in_band_tones", "in_band_nominal_kbps", "min_psr", "median_psr", "tones_with_psr_1", "packets_per_tone"]


def sort_key(r):
    return tuple(-1 if r[k] is None else r[k] for k in KEY) + (r.get("fc_hz", 0),)


def summarise(rows):
    groups = {}
    for r in rows:
        groups.setdefault(tuple(r[k] for k in KEY), []).append(r)
    out = []
    for key, rs in groups.items():
        p, psr = plan(rs[0]["_g"]), np.array([r["psr"] for r in rs])
        pay = key[3] / CPB / 1e3
        out.append(dict(zip(KEY, key), n_tones=len(rs), band_lo_hz=int(round(p.lo)), band_hi_hz=int(round(p.hi)),
                        payload_per_tone_kbps=pay, nominal_aggregate_kbps=len(rs) * pay,
                        delivered_aggregate_kbps=float(psr.sum()) * pay, in_band_tones=p.in_band,
                        in_band_nominal_kbps=None if p.in_band is None else p.in_band * pay, min_psr=psr.min(), median_psr=float(np.median(psr)),
                        tones_with_psr_1=int((psr == 1).sum()), packets_per_tone=rs[0]["packets_tx"]))
    return sorted(out, key=sort_key)


def fmt(x):
    if x is None:
        return ""
    if isinstance(x, (float, np.floating)):
        return str(int(x)) if float(x).is_integer() else f"{x:.6g}"
    return x


def write_csv(path, rows, cols):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(cols)
        w.writerows([[fmt(r[c]) for c in cols] for r in rows])


def run_tasks(tasks, workers):
    """Run victim tasks, largest first, on at most MAX_WORKERS processes."""
    cost = lambda t: len(tone_stream(t[0].seed, 0, t[0].n_packets, t[0].sync)[0]) * plan(t[0]).spc * (2 * t[0].M + 1)
    tasks = sorted(tasks, key=cost, reverse=True)
    rows, t0 = [], time.time()
    workers = max(1, min(MAX_WORKERS, workers))
    with Pool(workers, maxtasksperchild=50) as pool:
        for i, r in enumerate(pool.imap_unordered(run_victim, tasks), 1):
            rows += r
            if i % 200 == 0 or i == len(tasks):
                print(f"  {i}/{len(tasks)} victim tasks, {time.time() - t0:.0f} s", flush=True)
    return rows


def settings(nf=(None,), drives=((1.0, "none"),), targets=(None,)):
    return [dict(nf_db=a, drive_v=dv, split=sp, ebn0_target_db=e) for a in nf for dv, sp in drives for e in targets]


# ---------------------------------------------------------------- modes
CURVES = (("cohort", "50"), ("RHW1-LHW1", "50"), ("cohort", "1M"))
RCS = (50e3, 100e3, 200e3, 500e3, 1e6, 2e6)
ACI_FACTORS = (0.5, 0.75, 1.0, 1.25, 1.5, 2.0)
ACI_SEEDS = (1, 2, 3, 4, 5, 6)                        # async ACI: six draws of timing and packet offsets
# the as-built rates, the HBC Walsh-chip rates, then upward until the lower tone f - R_c nears 0 Hz (1 MHz steps,
# 0.5 MHz at the end; 9.95 M is the last 50 k step whose two tones fit the measured 0.05 to 20 MHz)
RATES = (50e3, 200e3, 500e3, 656.25e3, 1.3125e6, 2e6, 2.625e6, 5.25e6, 7e6, 8e6, 9e6, 9.5e6, 9.95e6)
KNEE_DB = 7.0                                         # Eb/N0 target on the slope of the frame model's PSR (0.56)
STAB_COLS = ["rc_cps", "tone_hz", "band_lo_hz", "band_hi_hz", "curve", "nf_db", "ebn0_target_db", "gain_m_db",
             "gain_p_db", "ec_n0_m_db", "ec_n0_p_db", "ec_n0_db", "packets", "chips", "chip_errors", "error_free", "psr",
             "psr_model"]


# the rates of the dispersion check and the
# channel variants (curve_db variant, extension); 'tones' is the model without dispersion, run alongside
ISI_RATES = (0.5e6, 1.3125e6, 2e6, 2.625e6, 5.25e6, 7e6, 9e6, 9.5e6)
ISI_CHANS = (("tones", ""), ("measured", "A"), ("smoothed", "A"), ("measured", "B"))
ISI_COLS = ["rc_cps", "tone_hz", "curve", "variant", "ext", "nf_db", "ebn0_target_db", "gain_m_db", "gain_p_db",
            "gd_m_ns", "gd_p_ns", "gd_dt_ns", "rms_ns", "e90_ns", "lag_ns", "margin", "packets", "chips",
            "chip_errors", "cer", "error_free", "psr"]


def sweep(a):
    tasks = []

    def add(g, st):
        tasks.extend((g, v, st) for v in range(len(plan(g).f)))

    grid = settings(nf=(10.0, 20.0, 30.0), drives=((1.0, "power"), (1.0, "peak"), (0.1, "power"), (0.016, "power")))
    for curve, load in CURVES:
        for Rc in RCS:
            for guard in (0.0, 0.25, 0.5):
                add(Geometry("grid", curve, load, Rc, guard, n_packets=a.grid_packets), grid)
    head = settings(nf=(20.0,), drives=((1.0, "power"), (1.0, "peak")))
    for mod in ("fsk", "ook"):
        for curve, load in CURVES[:2]:
            for Rc in RCS:
                for guard in (0.0, 0.25):
                    add(Geometry("headline", curve, load, Rc, guard, mod=mod, n_packets=a.headline_packets,
                                 M=M_DEFAULT[mod]), head)
    for fct in ACI_FACTORS:                              # 5 flat tones, every tone simulated against all others
        for seed in ACI_SEEDS:
            add(Geometry("aci", "flat", "", 500e3, fct - 1.0, n_packets=a.aci_packets, M=4, n_tones=5, seed=seed),
                settings(targets=(8.0, 10.0, 12.0)))
        add(Geometry("aci", "flat", "", 500e3, fct - 1.0, sync=True, n_packets=a.aci_packets, M=4, n_tones=5),
            settings(targets=(10.0,)))
    print(f"sweep: {len(tasks)} victim tasks on {min(MAX_WORKERS, a.workers)} workers", flush=True)
    t0 = time.time()
    rows = sorted(run_tasks(tasks, a.workers), key=sort_key)
    write_csv(RESULTS / "bcc_fdma_tones.csv", rows, TONE_COLS)
    write_csv(RESULTS / "bcc_fdma_summary.csv", summarise(rows), SUM_COLS)
    print(f"wrote {RESULTS / 'bcc_fdma_tones.csv'} ({len(rows)} rows) and bcc_fdma_summary.csv in {time.time() - t0:.0f} s")


def run(a):
    g = Geometry("run", a.curve, a.load if a.curve != "flat" else "", a.rc, a.guard, mod=a.mod, sync=a.sync,
                 n_packets=a.packets, M=a.neighbours if a.neighbours is not None else M_DEFAULT[a.mod], seed=a.seed, band=tuple(a.band) if a.band else None,
                 n_tones=a.n_tones)
    st = [dict(nf_db=None if a.ebn0 is not None else a.nf, drive_v=a.drive, split=a.split, ebn0_target_db=a.ebn0)]
    p = plan(g)
    print(f"{len(p.f)} tones, band {p.lo / 1e6:.3f}-{p.hi / 1e6:.3f} MHz, {p.spc} samples per chip")
    rows = sorted(run_tasks([(g, v, st) for v in range(len(p.f))], a.workers), key=sort_key)
    print("fc (MHz) | gain -/+ (dB) | Eb/N0 (dB) | chip errors | detected/decoded/error-free | PSR")
    for r in rows:
        print(f"{r['fc_hz'] / 1e6:8.3f} | {r['gain_m_db']:6.1f} {r['gain_p_db']:6.1f} | {r['ebn0_db']:6.1f} | "
              f"{r['chip_errors']:6d} ({r['cer']:.2e}) | {r['detected']}/{r['decoded']}/{r['error_free']} of {r['packets_tx']} | {r['psr']:.2f}")
    s = summarise(rows)[0]
    print(f"payload per tone {s['payload_per_tone_kbps']:.1f} kbit/s, nominal {s['nominal_aggregate_kbps']:.0f}, "
          f"delivered {s['delivered_aggregate_kbps']:.0f} kbit/s, min PSR {s['min_psr']:.2f}, tones with PSR 1: {s['tones_with_psr_1']}")


# ---------------------------------------------------------------- equivalence proofs
def passband_errors(f, amps, bits, frac, Rc, N0, seed, mod="fsk", phy2_tx=True, chunk=20000):
    """Direct real-passband reference: all tones summed at one sample rate, real noise of one-sided density N0, the
    PHY2 band-pass and detector per tone. amps: (n, 2) received amplitude for bit 0 / bit 1; frac: timing offsets in
    chips (multiples of 1/32). Transmitter: phy2_tx -> PHY2 modulate and apply_channel (BFSK; its template restarts
    every chip, equal to continuous phase up to one phase per tone when f_j / R_c is an integer); otherwise an
    independent continuous-phase modulator valid at any f_j: chip k of tone j starts at t0 = (k + frac_j) / R_c and
    is a cos(2 pi f_j t + 2 pi s R_c (t - t0)), s = -/+1 (OOK: s = 0, a = 0 when off), t absolute across chunks.
    Returns chip errors per tone. Chunks restart the receive filter (edge effect: one chip per chunk)."""
    spc = 32 * int(np.ceil(4 * (max(f) + Rc) / Rc / 32))
    fs = spc * Rc
    d = [int(round(x * spc)) for x in frac]
    assert all(abs(dj - x * spc) < 1e-9 for dj, x in zip(d, frac)), "offsets must be multiples of 1/32 chip"
    phys = [v2.PHY2(fc=fj, Rc=Rc, mod=mod, A=1.0, N0=N0, fs=fs,
                    gain_db_tones=tuple(20 * np.log10(a)) if phy2_tx else None) for fj, a in zip(f, amps)]
    r = np.arange(spc)
    rng = np.random.default_rng(seed)
    err = np.zeros(len(f), int)
    for c0 in range(0, len(bits[0]), chunk):
        ch = [b[c0:c0 + chunk] for b in bits]
        nc = len(ch[0])
        x = np.zeros((nc + 1) * spc)
        for j, (phy, c, dj) in enumerate(zip(phys, ch, d)):
            if phy2_tx:
                pm = c.astype(int) * 2 - 1
                w = phy.apply_channel(phy.modulate(pm), pm)[0]
            else:
                sgn = (c * 2.0 - 1.0) if mod == "fsk" else np.zeros(nc)
                t = ((c0 + np.arange(nc))[:, None] * spc + dj + r) / fs
                w = amps[j][c][:, None] * np.cos(2 * np.pi * f[j] * t + 2 * np.pi * sgn[:, None] * Rc * r / fs)
            x[dj:dj + nc * spc] += w.ravel()
        x += rng.normal(0.0, np.sqrt(N0 * fs / 2), x.size)
        for j, (phy, c, dj) in enumerate(zip(phys, ch, d)):
            y = phy.rx_filter(x)[dj:dj + nc * spc].reshape(nc, spc)
            err[j] += int(np.count_nonzero((phy.detect(y) > 0) != (c > 0)))
    return err


def baseband_errors(f, amps, bits, d, Rc, N0, alpha, M, spc, seed, mod="fsk", victims=None):
    """Chip errors of the complex-baseband model per victim tone (random-chip streams)."""
    out = []
    for v in (range(len(f)) if victims is None else victims):
        Cs, Cn = victim_corr(v, f, amps, bits, d, Rc, mod, spc, M, [seed, v, 2])
        out.append(int(np.count_nonzero(detect(Cs, Cn, alpha, N0, mod) != bits[v])))
    return np.array(out)


def theory_ber(ebn0_db):
    return 0.5 * np.exp(-10 ** (np.asarray(ebn0_db) / 10) / 2)


def offset_db(ber, ebn0_db):
    """Simulated Eb/N0 minus the Eb/N0 at which 0.5 exp(-Eb/(2 N0)) equals ber (positive = worse than theory)."""
    return ebn0_db - 10 * np.log10(2 * np.log(1 / (2 * ber)))


def z_diff(e1, e2, n):
    p1, p2 = e1 / n, e2 / n
    se = np.sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / n)
    return (p2 - p1) / se if se > 0 else 0.0


def check_a(n, ebs, seed=1):
    """(a) single tone, random chips: baseband vs theory and vs PHY2 (fc 6 MHz, Rc 200 kchip/s, 0 dB, A 1 V)."""
    Rc, rows = 200e3, []
    bits = np.random.default_rng(seed).integers(0, 2, n).astype(np.uint8)
    for e in ebs:
        N0 = 0.5 / Rc / 10 ** (e / 10)
        bb = baseband_errors(np.array([6e6]), np.ones((1, 2)), [bits], [0], Rc, N0, 1.0, 2, 32, seed)[0]
        pb = 0
        for c0 in range(0, n, 40000):                    # PHY2 in chunks (5 M samples each)
            c = bits[c0:c0 + 40000].astype(int) * 2 - 1
            pb += int(np.count_nonzero(v2.PHY2(fc=6e6, Rc=Rc, mod="fsk", A=1.0, N0=N0, seed=seed + c0).run(c) != c))
        rows.append(dict(check="a", setting="single tone vs PHY2", tone=0, ebn0_db=e, chips=n, errors_ref=pb,
                         cer_ref=pb / n, errors_bb=int(bb), cer_bb=bb / n, z=z_diff(pb, bb, n),
                         theory=float(theory_ber(e)), offset_db=offset_db(bb / n, e) if bb else float("nan")))
    return rows


def check_b(n, seed=2):
    """(b) 3 tones at Rc 500 kchip/s, 4 Rc spacing, cohort 50 Ohm gains at 8, 10, 12 MHz, async offsets 0, 1/4, 3/4
    chip, best tone Eb/N0 10 dB: direct passband (PHY2 signals) vs baseband, per tone."""
    Rc = 500e3
    p = plan(Geometry("check", "cohort", "50", Rc, 0.0, band=(7e6, 13e6)))
    bits = [np.random.default_rng([seed, j]).integers(0, 2, n).astype(np.uint8) for j in range(3)]
    frac = (0.0, 0.25, 0.75)
    N0 = p.e_chip.max() / 10 ** 1.0
    pb = passband_errors(p.f, p.amps, bits, frac, Rc, N0, seed)
    bb = baseband_errors(p.f, p.amps, bits, [int(x * p.spc) for x in frac], Rc, N0, 1.0, 2, p.spc, seed)
    return [dict(check="b", setting=f"3 tones passband vs baseband, gains {p.gdb[i, 0]:.1f}/{p.gdb[i, 1]:.1f} dB",
                 tone=int(round(p.f[i])), ebn0_db=10 * np.log10(p.e_chip[i] / N0), chips=n, errors_ref=int(pb[i]),
                 cer_ref=pb[i] / n, errors_bb=int(bb[i]), cer_bb=bb[i] / n, z=z_diff(pb[i], bb[i], n),
                 theory=float(theory_ber(10 * np.log10(p.e_chip[i] / N0))), offset_db=float("nan")) for i in range(3)]


def check_b2(n, seed=2, mod="fsk"):
    """(b2) independent continuous-phase passband transmitter vs baseband, timing offsets that are not quarter chips
    (7/32, 21/32, 29/32, 3/32, 13/32), so the carrier-phase term (f_j - f_v) t0 of victim_corr is not a whole number
    of cycles. BFSK: 3 tones as (b), cohort 50 Ohm gains, centred 0.37 R_c off the integer grid (8.185, 10.185 and
    12.185 MHz, f / R_c = 16.37 to 24.37), best tone Eb/N0 10 dB, every tone. OOK: 5 flat tones at 2 R_c on the
    integer grid (8 to 12 MHz), 'on' chip Eb/N0 12 dB, centre tone."""
    Rc = 500e3
    if mod == "fsk":
        g, ebn0, victims = Geometry("check", "cohort", "50", Rc, 0.0, band=(7e6 + 0.37 * Rc, 13e6 + 0.37 * Rc)), 10.0, [0, 1, 2]
    else:
        g, ebn0, victims = Geometry("check", "flat", "", Rc, 0.0, mod="ook", M=4, n_tones=5), 12.0, [2]
    p = plan(g)
    assert p.spc % 32 == 0, p.spc
    nt = len(p.f)
    frac = (7 / 32, 21 / 32, 29 / 32, 3 / 32, 13 / 32)[:nt]
    bits = [np.random.default_rng([seed, j]).integers(0, 2, n).astype(np.uint8) for j in range(nt)]
    N0 = p.e_chip.max() / 10 ** (ebn0 / 10)
    pb = passband_errors(p.f, p.amps, bits, frac, Rc, N0, seed, mod, phy2_tx=False)
    bb = baseband_errors(p.f, p.amps, bits, [int(round(x * p.spc)) for x in frac], Rc, N0, 1.0, g.M, p.spc, seed, mod,
                         victims)
    label = ("3 tones continuous-phase passband off the grid vs baseband" if mod == "fsk"
             else "OOK centre of 5 tones at 2 R_c, gated-carrier passband vs baseband")
    return [dict(check="b2", setting=f"{label}, f/R_c {p.f[v] / Rc:.2f}, offset {frac[v] * 32:.0f}/32 chip",
                 tone=int(round(p.f[v])), ebn0_db=10 * np.log10(p.e_chip[v] / N0), chips=n, errors_ref=int(pb[v]),
                 cer_ref=pb[v] / n, errors_bb=int(e), cer_bb=e / n, z=z_diff(pb[v], e, n),
                 theory=float(theory_ber(10 * np.log10(p.e_chip[v] / N0))) if mod == "fsk" else float("nan"),
                 offset_db=float("nan")) for v, e in zip(victims, bb)]


def check_c(n, mod="fsk", ebn0=8.0, m_ref=2, m_new=4, seed=3):
    """(c) 2 m_new + 1 flat tones (4 R_c spacing BFSK, 2 R_c OOK), async: centre tone with m_ref vs m_new neighbours
    on each side, same samples per chip and the same noise, so the difference is the far neighbours only."""
    Rc, n_t = 500e3, 2 * m_new + 1
    p = plan(Geometry("check", "flat", "", Rc, 0.0, mod=mod, M=m_new, n_tones=n_t))
    bits = [np.random.default_rng([seed, j]).integers(0, 2, n).astype(np.uint8) for j in range(n_t)]
    d = offsets(Geometry("check", "flat", "", Rc, 0.0, seed=seed), n_t, p.spc)
    N0 = p.e_chip.max() / 10 ** (ebn0 / 10)
    e1, e2 = (baseband_errors(p.f, p.amps, bits, d, Rc, N0, 1.0, M, p.spc, seed, mod, victims=[m_new])[0]
              for M in (m_ref, m_new))
    return [dict(check="c", setting=f"{mod} M={m_ref} (ref) vs M={m_new} (bb), centre of {n_t} tones", tone=m_new,
                 ebn0_db=ebn0, chips=n, errors_ref=int(e1), cer_ref=e1 / n, errors_bb=int(e2), cer_bb=e2 / n,
                 z=z_diff(e1, e2, n), theory=float(theory_ber(ebn0)) if mod == "fsk" else float("nan"),
                 offset_db=float("nan"))]


def _check_job(job):
    name, kw = job
    return {"a": check_a, "b": check_b, "b2": check_b2, "c": check_c}[name](**kw)


def check(a):
    n = a.chips
    jobs = [("a", dict(n=n, ebs=(4.0, 6.0, 8.0, 10.0, 12.0)))] + [("b", dict(n=n, seed=s)) for s in (2, 5, 7, 11)] + [
        ("b2", dict(n=n, seed=s, mod=m)) for m in ("fsk", "ook") for s in (2, 5, 7, 11)] + [
        ("c", dict(n=n, mod="fsk", ebn0=8.0)), ("c", dict(n=n, mod="fsk", ebn0=30.0)),
        ("c", dict(n=n, mod="ook", ebn0=12.0)), ("c", dict(n=n, mod="ook", ebn0=12.0, m_ref=4, m_new=6)),
        ("c", dict(n=n, mod="ook", ebn0=30.0, m_ref=4, m_new=6))]
    with Pool(max(1, min(MAX_WORKERS, a.workers, len(jobs)))) as pool:
        rows = [r for rs in pool.map(_check_job, jobs) for r in rs]
    pooled = {}                                          # (b), (b2): pool the four seeds per tone
    for r in (r for r in rows if r["check"] in ("b", "b2")):
        q = pooled.setdefault((r["check"], r["tone"]), dict(r, chips=0, errors_ref=0, errors_bb=0))
        q["chips"] += r["chips"]; q["errors_ref"] += r["errors_ref"]; q["errors_bb"] += r["errors_bb"]
    for q in pooled.values():
        q.update(cer_ref=q["errors_ref"] / q["chips"], cer_bb=q["errors_bb"] / q["chips"],
                 z=z_diff(q["errors_ref"], q["errors_bb"], q["chips"]), setting=q["setting"] + ", 4 seeds pooled")
    rows = [r for r in rows if r["check"] == "a"] + list(pooled.values()) + [r for r in rows if r["check"] == "c"]
    cols = ["check", "setting", "tone", "ebn0_db", "chips", "errors_ref", "cer_ref", "errors_bb", "cer_bb", "z", "theory", "offset_db"]
    write_csv(RESULTS / "bcc_fdma_equivalence.csv", rows, cols)
    print("| " + " | ".join(cols) + " |\n|" + "---|" * len(cols))
    for r in rows:
        print("| " + " | ".join(str(fmt(r[c])) for c in cols) + " |")


def rates(a):
    """ stability study: waveform PSR of one BFSK tone per measured path and chip rate (chunks of packets with
    their own seeds, pooled), next to the engine's frame model at the same Ec/N0; at 1 V and NF 10 / 20 / 30, and at
    the Eb/N0 target KNEE_DB, where the frame model is checked on the slope of its PSR curve."""
    import channel as ch                                 # the engine's link and frame model, only for this mode
    import phy
    import protocols as pr
    with open(cc.CURVES_CSV, encoding="utf-8") as fh:
        curves = sorted({r["curve"] for r in csv.DictReader(fh) if r["load"] == "50"})     # 36 paths + cohort
    nfs = (10.0, 20.0, 30.0)
    st = settings(nf=nfs) + settings(targets=(KNEE_DB,))
    tasks = [(Geometry("rates", c, "50", rc, 0.0, n_packets=a.packets, seed=1 + 100 * i + j,
                       band=(cc.tone_centre(rc),) * 2, n_tones=1), 0, st)
             for rc in RATES for i, c in enumerate(curves) for j in range(a.chunks)]
    print(f"rates: {len(tasks)} victim tasks on {min(MAX_WORKERS, a.workers)} workers", flush=True)
    pooled = {}
    for r in run_tasks(tasks, a.workers):
        q = pooled.setdefault((r["Rc"], r["curve"], r["nf_db"] or 0.0, r["ebn0_target_db"] or 0.0),
                              dict(r, packets=0, chips=0, chip_errors=0, error_free=0))
        for x, y in (("packets", "packets_tx"), ("chips", "chips"), ("chip_errors", "chip_errors"),
                     ("error_free", "error_free")):
            q[x] += r[y]
    spec = pr.get_spec("bcc", "thesis")                  # one segment over the frame: the same outcome at any R_c
    lay = pr.frame_layout(spec, PSDU_BYTES)
    rows = []
    for (rc, curve, nf, target), q in sorted(pooled.items()):
        fc, g = cc.tone_centre(rc), np.array([q["gain_m_db"], q["gain_p_db"]])
        e = (ch.ec_n0_db(1.0, g, rc, nf) if nf else target + g - 10 * np.log10(np.mean(10 ** (g / 10)))).tolist()
        eq = float(ch.bfsk_ec_n0_db(e))
        model = phy.frame_outcome(None, [(0.0, lay[-1]["t1_us"], 10 ** (eq / 10), np.inf)], lay, spec.phy,
                                  sampled=False)["success"]
        rows.append(dict(rc_cps=rc, tone_hz=fc, band_lo_hz=fc - 2 * rc, band_hi_hz=fc + 2 * rc, curve=curve,
                         nf_db=q["nf_db"], ebn0_target_db=q["ebn0_target_db"], gain_m_db=g[0], gain_p_db=g[1],
                         ec_n0_m_db=e[0], ec_n0_p_db=e[1], ec_n0_db=eq, packets=q["packets"], chips=q["chips"],
                         chip_errors=q["chip_errors"], error_free=q["error_free"], psr=q["error_free"] / q["packets"],
                         psr_model=float(model)))
    write_csv(RESULTS / "bcc_rate_stability.csv", rows, STAB_COLS)
    print("| R_c (Mchip/s) | tone (MHz) | f -/+ 2 R_c (MHz) | weakest path | Ec/N0 at f - R_c, f + R_c (dB, NF 20) | "
          f"PSR min / median NF 10 | NF 20 | NF 30 | median at {KNEE_DB:g} dB: waveform / model |\n|" + "---|" * 9)
    for rc in RATES:
        by = {x: [r for r in rows if r["rc_cps"] == rc and x in (r["nf_db"], r["ebn0_target_db"])]
              for x in nfs + (KNEE_DB,)}
        w = min(by[20.0], key=lambda r: r["ec_n0_db"])
        pm = [f"{min(r['psr'] for r in by[nf]):.4g} / {np.median([r['psr'] for r in by[nf]]):.4g}" for nf in nfs]
        print(f"| {rc / 1e6:g} | {w['tone_hz'] / 1e6:.3f} | {w['band_lo_hz'] / 1e6:.2f} to "
              f"{w['band_hi_hz'] / 1e6:.2f} | {w['curve']} | {w['ec_n0_m_db']:.1f}, {w['ec_n0_p_db']:.1f} | "
              f"{' | '.join(pm)} | {np.median([r['psr'] for r in by[KNEE_DB]]):.3f} / "
              f"{np.median([r['psr_model'] for r in by[KNEE_DB]]):.3f} |")


def isi(a):
    """: the one-tone stability study of rates() with the minimum-phase response of each 50 Ohm curve applied to
    the transmitted wave, per variant of ISI_CHANS, in real passband at every rate, 1 V at NF 10 / 20 / 30 and the
    Eb/N0 target KNEE_DB; the same chunk seeds and chip streams as rates()."""
    with open(cc.CURVES_CSV, encoding="utf-8") as fh:
        curves = sorted({r["curve"] for r in csv.DictReader(fh) if r["load"] == "50"})     # 36 paths + cohort
    nfs = (10.0, 20.0, 30.0)
    tasks = [(c, 1 + 100 * i + j, a.packets, ISI_RATES, ISI_CHANS, nfs) for i, c in enumerate(curves)
             for j in range(a.chunks)]
    print(f"isi: {len(tasks)} tasks on {a.workers} workers", flush=True)
    pooled, t0 = {}, time.time()
    with Pool(a.workers) as pool:                       # a curve's chunks on one worker: its channels built once
        for i, rs in enumerate(pool.imap_unordered(run_isi, tasks, chunksize=a.chunks), 1):
            for r in rs:
                q = pooled.setdefault(tuple(r[k] or 0 for k in ("rc_cps", "curve", "variant", "ext", "nf_db")),
                                      dict(r, packets=0, chips=0, chip_errors=0, error_free=0))
                for x in ("packets", "chips", "chip_errors", "error_free"):
                    q[x] += r[x] or 0
            print(f"  {i}/{len(tasks)} tasks, {time.time() - t0:.0f} s", flush=True)
    rows = []
    for key, q in sorted(pooled.items()):
        rc, curve, v, e = key[:4]
        fc = cc.tone_centre(rc)
        g = (isi_channel(curve, v, e, rc)[1][0] if v == "tones" else
             10 ** (curve_db(curve, v, e, [fc - rc, fc + rc]) / 20))
        gd = gd_spread(curve, v, e, rc) + smear(curve, v, e) if v != "tones" else (None,) * 5
        knee = q["nf_db"] is None
        rows.append(dict(q, tone_hz=fc, gain_m_db=20 * np.log10(g[0]), gain_p_db=20 * np.log10(g[1]), gd_m_ns=gd[0],
                         gd_p_ns=gd[1], gd_dt_ns=gd[2], rms_ns=gd[3], e90_ns=gd[4], cer=q["chip_errors"] / q["chips"],
                         error_free=None if knee else q["error_free"],
                         psr=None if knee else q["error_free"] / q["packets"]))
    write_csv(RESULTS / "bcc_rate_isi.csv", rows, ISI_COLS)
    print(f"| R_c (Mchip/s) | channel | PSR min / median NF 10 | NF 20 | NF 30 | most chip errors of a path at NF 30 | "
          f"smallest margin | cost at {KNEE_DB:g} dB, median / largest (dB) |\n|" + "---|" * 8)
    for rc in ISI_RATES:
        ref = {r["curve"]: r["cer"] for r in rows
               if r["rc_cps"] == rc and r["variant"] == "tones" and r["ebn0_target_db"]}
        for v, e in ISI_CHANS:
            sel = [r for r in rows if (r["rc_cps"], r["variant"], r["ext"]) == (rc, v, e)]
            by = [[r for r in sel if r["nf_db"] == nf] for nf in nfs]
            pm = " | ".join(f"{min(r['psr'] for r in x):.4g} / {np.median([r['psr'] for r in x]):.4g}" for x in by)
            cost = [10 * np.log10(np.log(0.5 / ref[r["curve"]]) / np.log(0.5 / r["cer"])) for r in sel
                    if r["ebn0_target_db"]]
            print(f"| {rc / 1e6:g} | {v} {e} | {pm} | {max(r['chip_errors'] for r in by[2])} | "
                  f"{min(r['margin'] for r in sel):.3f} | {np.median(cost):.2f} / {max(cost):.2f} |")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="mode", required=True)
    r = sub.add_parser("run", help="one configuration")
    r.add_argument("--curve", default="cohort", help="'cohort', a path 'TX-RX' or 'flat'")
    r.add_argument("--load", default="50", choices=("50", "1M"))
    r.add_argument("--rc", type=float, default=500e3, help="chip rate (chip/s)")
    r.add_argument("--guard", type=float, default=0.0)
    r.add_argument("--nf", type=float, default=20.0, help="receiver noise figure (dB)")
    r.add_argument("--ebn0", type=float, default=None, help="Eb/N0 target of the best tone (dB), replaces --nf")
    r.add_argument("--drive", type=float, default=1.0, help="per-tone amplitude A (V) before the split")
    r.add_argument("--split", default="power", choices=("none", "power", "peak"))
    r.add_argument("--mod", default="fsk", choices=("fsk", "ook"))
    r.add_argument("--sync", action="store_true", help="align chip timing and packet starts of all tones")
    r.add_argument("--packets", type=int, default=6, help="packets per tone")
    r.add_argument("--neighbours", type=int, default=None, help="neighbour tones simulated on each side (default 2 BFSK, 4 OOK)")
    r.add_argument("--band", type=float, nargs=2, default=None, metavar=("LO", "HI"), help="tone band (Hz)")
    r.add_argument("--n-tones", type=int, default=None, help="number of tones (default: as many as fit)")
    r.add_argument("--seed", type=int, default=1)
    r.add_argument("--workers", type=int, default=MAX_WORKERS)
    s = sub.add_parser("sweep", help="the grid, headline rows and ACI series")
    s.add_argument("--workers", type=int, default=MAX_WORKERS)
    s.add_argument("--grid-packets", type=int, default=6)
    s.add_argument("--headline-packets", type=int, default=30)
    s.add_argument("--aci-packets", type=int, default=30)
    c = sub.add_parser("check", help="equivalence proofs (a), (b), (b2), (c)")
    c.add_argument("--chips", type=int, default=200000)
    c.add_argument("--workers", type=int, default=MAX_WORKERS)
    t = sub.add_parser("rates", help=" chip-rate stability study, one tone per measured path")
    t.add_argument("--packets", type=int, default=100, help="packets per chunk")
    t.add_argument("--chunks", type=int, default=4, help="chunks (own seeds) per path and rate, pooled")
    t.add_argument("--workers", type=int, default=MAX_WORKERS)
    i = sub.add_parser("isi", help=" dispersion check: rates() through the minimum-phase channel")
    i.add_argument("--packets", type=int, default=100, help="packets per chunk")
    i.add_argument("--chunks", type=int, default=4, help="chunks (own seeds) per path, pooled")
    i.add_argument("--workers", type=int, default=2, choices=(1, 2), help="worker processes (at most 2)")
    a = ap.parse_args()
    {"run": run, "sweep": sweep, "check": check, "rates": rates, "isi": isi}[a.mode](a)


if __name__ == "__main__":
    main()

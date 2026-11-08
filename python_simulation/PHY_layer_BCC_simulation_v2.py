"""PHY_layer_BCC_simulation_v2.py

Vectorised BCC physical layer with a physically defined noise model, a receive filter, near-optimal detectors,
optional HackRF-like front-end impairments and an electrode-contact dropout model.

Differences from PHY_layer_BCC_simulation.py
--------------------------------------------
* Noise is specified as a one-sided spectral density N0 (V^2/Hz); the per-sample standard deviation is
  sqrt(N0 * fs / 2), so results do not depend on the sample-rate convention and Eb/N0 is defined.
* A receive band-pass filter (FIR, default bandwidth 4 * Rc around fc for BFSK, 2 * Rc for OOK) precedes detection.
* BFSK detector 'energy': non-coherent tone-energy comparison over each chip (quadrature correlators at f_m and f_p).
  BFSK detector 'discriminator': the v1 Hilbert/instantaneous-frequency method, kept for comparison.
* OOK detector: chip energy at fc against a threshold; 'adaptive' threshold = midpoint between the two energy clusters,
  or a fixed absolute threshold (as the GNU Radio receiver uses).
* Channel gain may differ per tone (two-point frequency selectivity) and may be time varying through a two-state
  contact model (connected / dropout) with geometric dwell times, to represent partial electrode contact.
* Front end: 8-bit quantisation with a full-scale setting, DC offset, hard clipping.

Usage (library):
    phy = PHY2(fc=8e6, Rc=50e3, mod='fsk', A=1.0, N0=1e-14, gain_db=-52.4)
    rx_chips = phy.run(chips)             # chips in {-1,+1} (fsk) or {0,1} (ook)
    phy.ebn0_db()                         # Eb/N0 implied by A, gain and N0
"""
from __future__ import annotations
import numpy as np
from scipy.signal import firwin, fftconvolve


class Contact:
    """Two-state contact model on the channel gain. States: 0 = connected (nominal gain), 1 = dropout (gain - depth_db).
    Dwell times are geometric with means given in chips."""
    def __init__(self, mean_good_chips: float, mean_bad_chips: float, depth_db: float, start_bad: bool = False, seed: int | None = None):
        self.p_gb = 1.0 / max(mean_good_chips, 1.0); self.p_bg = 1.0 / max(mean_bad_chips, 1.0)
        self.depth_db = depth_db; self.start_bad = start_bad; self.rng = np.random.default_rng(seed)

    def states(self, n: int) -> np.ndarray:
        s = np.zeros(n, dtype=np.int8); cur = 1 if self.start_bad else 0
        u = self.rng.random(n)
        for i in range(n):
            s[i] = cur
            if cur == 0 and u[i] < self.p_gb: cur = 1
            elif cur == 1 and u[i] < self.p_bg: cur = 0
        return s


class PHY2:
    def __init__(self, fc: float, Rc: float, mod: str = "fsk", A: float = 1.0, N0: float | None = None, sigma: float | None = None,
                 gain_db: float = 0.0, gain_db_tones: tuple | None = None, m: float = 1.0, fs: float | None = None,
                 rx_bw: float | None = None, detector: str = "energy", ook_threshold: str | float = "adaptive",
                 adc_bits: int | None = None, adc_fullscale: float = 1.0, dc_offset: float = 0.0, clip: float | None = None,
                 contact: Contact | None = None, seed: int | None = None):
        self.fc, self.Rc, self.mod, self.A, self.m = fc, Rc, mod, A, m
        self.fm, self.fp = fc - m * Rc, fc + m * Rc
        # sample rate: integer samples per chip, at least 4 x highest tone
        if fs is None:
            spc = int(np.ceil(4.0 * (fc + m * Rc) / Rc)); fs = spc * Rc
        self.fs = fs; self.spc = int(round(fs / Rc))
        assert abs(self.spc * Rc - fs) < 1e-6, "fs must be an integer multiple of Rc"
        if N0 is None and sigma is None: raise ValueError("give N0 (V^2/Hz) or sigma (V per sample)")
        self.N0 = N0 if N0 is not None else 2.0 * sigma ** 2 / fs
        self.sigma = float(np.sqrt(self.N0 * fs / 2.0))
        self.gain_db = gain_db
        self.gain_tones = gain_db_tones if gain_db_tones is not None else (gain_db, gain_db)   # (at f_m, at f_p)
        self.rx_bw = rx_bw if rx_bw is not None else (4.0 * m * Rc if mod == "fsk" else 2.0 * Rc)
        self.detector, self.ook_threshold = detector, ook_threshold
        self.adc_bits, self.adc_fullscale, self.dc_offset, self.clip = adc_bits, adc_fullscale, dc_offset, clip
        self.contact = contact; self.rng = np.random.default_rng(seed)
        t = np.arange(self.spc) / fs
        self.tpl_m = np.cos(2 * np.pi * self.fm * t); self.tpl_p = np.cos(2 * np.pi * self.fp * t); self.tpl_c = np.cos(2 * np.pi * fc * t)
        self.q_m = np.sin(2 * np.pi * self.fm * t); self.q_p = np.sin(2 * np.pi * self.fp * t); self.q_c = np.sin(2 * np.pi * fc * t)
        ntaps = int(4 * fs / self.rx_bw) | 1
        self.taps = firwin(ntaps, [max(fc - self.rx_bw / 2, 1.0), min(fc + self.rx_bw / 2, fs / 2 - 1.0)], pass_zero=False, fs=fs)

    # ---------------------------------------------------------------- energy bookkeeping
    def rx_amplitude(self, gain_db: float | None = None) -> float:
        g = self.gain_db if gain_db is None else gain_db
        return self.A * 10 ** (g / 20)

    def ebn0_db(self, gain_db: float | None = None) -> float:
        """Energy per chip (one BFSK tone or an OOK 'one' chip) over N0."""
        a = self.rx_amplitude(gain_db); Eb = a ** 2 / 2.0 / self.Rc
        return 10 * np.log10(Eb / self.N0)

    # ---------------------------------------------------------------- transmitter and channel
    def modulate(self, chips: np.ndarray) -> np.ndarray:
        chips = np.asarray(chips)
        if self.mod == "fsk":
            sel = (chips > 0)[:, None]
            wave = np.where(sel, self.tpl_p[None, :], self.tpl_m[None, :]) * self.A
        else:
            wave = (chips > 0)[:, None] * self.tpl_c[None, :] * self.A
        return wave.astype(np.float64)                     # shape (n_chips, spc)

    def apply_channel(self, wave: np.ndarray, chips: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        n = wave.shape[0]
        if self.mod == "fsk":
            g = np.where(chips > 0, self.gain_tones[1], self.gain_tones[0]).astype(float)
        else:
            g = np.full(n, self.gain_db, dtype=float)
        states = np.zeros(n, dtype=np.int8)
        if self.contact is not None:
            states = self.contact.states(n); g = g - states * self.contact.depth_db
        return wave * (10 ** (g / 20))[:, None], states

    def add_noise(self, wave: np.ndarray) -> np.ndarray:
        return wave + self.rng.normal(0.0, self.sigma, size=wave.shape)

    def frontend(self, x: np.ndarray) -> np.ndarray:
        x = x + self.dc_offset
        if self.clip is not None: x = np.clip(x, -self.clip, self.clip)
        if self.adc_bits is not None:
            q = self.adc_fullscale / (2 ** (self.adc_bits - 1))
            x = np.clip(np.round(x / q) * q, -self.adc_fullscale, self.adc_fullscale - q)
        return x

    def rx_filter(self, x_flat: np.ndarray) -> np.ndarray:
        return fftconvolve(x_flat, self.taps, mode="same")

    # ---------------------------------------------------------------- detection
    def detect(self, x: np.ndarray) -> np.ndarray:
        """x: (n_chips, spc) filtered samples -> chips."""
        if self.mod == "fsk":
            if self.detector == "energy":
                Em = (x @ self.tpl_m) ** 2 + (x @ self.q_m) ** 2; Ep = (x @ self.tpl_p) ** 2 + (x @ self.q_p) ** 2
                return np.where(Ep > Em, 1, -1)
            elif self.detector == "discriminator":
                from scipy.signal import hilbert
                out = np.empty(x.shape[0], dtype=int); t = np.arange(self.spc) / self.fs
                for i in range(x.shape[0]):
                    z = hilbert(x[i]) * np.exp(-1j * 2 * np.pi * self.fc * t)
                    f_inst = np.diff(np.unwrap(np.angle(z))) * self.fs / (2 * np.pi)
                    out[i] = 1 if np.mean(f_inst) > 0 else -1
                return out
            raise ValueError(self.detector)
        E = (x @ self.tpl_c) ** 2 + (x @ self.q_c) ** 2
        if self.ook_threshold == "adaptive":
            lo, hi = np.percentile(E, 25), np.percentile(E, 75); thr = 0.5 * (lo + hi)
        else:
            thr = float(self.ook_threshold) ** 2 * (self.spc / 2) ** 2       # threshold given as amplitude in V
        return (E > thr).astype(int)

    # ---------------------------------------------------------------- end to end
    def run(self, chips: np.ndarray, return_states: bool = False):
        chips = np.asarray(chips)
        w = self.modulate(chips); w, states = self.apply_channel(w, chips); w = self.add_noise(w)
        x = self.frontend(w.ravel()); x = self.rx_filter(x).reshape(w.shape)
        rx = self.detect(x)
        return (rx, states) if return_states else rx


def ber_vs_ebn0(fc=6e6, Rc=200e3, ebn0_list=(4, 6, 8, 10, 12, 14), n_chips=20000, detector="energy", rx_bw=None, seed=1, **kw):
    """Textbook check: returns list of (Eb/N0 dB, simulated BER) for random BFSK chips at 0 dB gain, A = 1 V."""
    rng = np.random.default_rng(seed); chips = rng.choice([-1, 1], size=n_chips); out = []
    for e in ebn0_list:
        Eb = 0.5 / Rc; N0 = Eb / 10 ** (e / 10)
        phy = PHY2(fc=fc, Rc=Rc, mod="fsk", A=1.0, N0=N0, gain_db=0.0, detector=detector, rx_bw=rx_bw, seed=seed, **kw)
        rx = phy.run(chips); out.append((e, float(np.mean(rx != chips))))
    return out


if __name__ == "__main__":
    for det, bw in (("energy", None), ("discriminator", None), ("discriminator", 1e9)):
        print(det, "rx_bw", bw, ber_vs_ebn0(detector=det, rx_bw=bw, n_chips=4000))

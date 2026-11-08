"""Physical-layer abstraction of the coexistence simulator: bit, chip, block and frame error models.

Purpose: the bit or chip error probability of each modulation (PARAMETERS.md 3.1 to 3.8 "BER model"), the per-bit
SNR of a segment from its SINR in the receiver noise bandwidth (5.4), bounded-distance block decoding (BCH of 3.4
and 3.5), the Walsh-16 symbol of the BCC and HBC stacks (3.8), and the outcome of one frame from the SINR segments
it saw, with the sync threshold and the co-channel C/I floor (5.5, DESIGN.md 5.2), lowered per block by its
processing and BCH coding gain (ci_gain).
Inputs:  linear SNR or SINR values; a frame layout (protocols.frame_layout) and the phy dict of a ProtocolSpec.
         phy fields read here: mod, noise_bw_hz, bit_rate_bps (on-air binary symbol or chip rate, after FEC and
         repetition: 600e3 for 802.15.6 NB at every S), rep (repetition or spreading factor S, combined ideally;
         per block override 'rep'), payload_crc, ci_co_lin (None or 0: no floor), optional ci_pg (processing gain of
         the payload against the floor beyond rep: HBC FSC SF / 8), optional ci_credit (False: no block lowers the
         floor), sync_sinr_db (None: no threshold),
         optional h (modulation index of 'gfsk_marcum', default 0.5). BLE Coded: bit_rate_bps 1e6 (symbol rate),
         rep 1, the S=2 / S=8 gain of PARAMETERS 3.1 (model gain minus extra loss) rides on ('conv', gain_db) blocks.
Outputs: probabilities (ber, block_fail_prob, walsh16_symbol_ok_prob), SNRs (gamma_from_sinr, sinr_at_ber) and
         frame_outcome dicts.
argv:    none (library); tests in test_phy.py.
"""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import comb, factorial, i0e
from scipy.stats import binom, ncx2

# PHY_layer_BCC_simulation_v2 runs 0.29..0.40 dB better than 0.5 exp(-g/2) (PARAMETERS 3.8)
V2_OFFSET_DB = json.loads(Path(__file__).with_name("parameters.json").read_text(encoding="utf-8"))["bcc"]["ber"]["sim_offset_dB"]
# Credit of a block against the co-channel C/I floor, calibrated on the BLE specification's own coded C/I
# values (PARAMETERS 3.1: 21 / 17 / 12 dB for 1M, S=2, S=8; the co-channel gains are 4 and 9 dB where the sensitivity
# gains are 5 and 12 dB). Repetition or spreading: S=2 -> S=8 is the P = 4 pattern mapper alone and buys 17 - 12 = 5 dB,
# so a factor S buys 10 log10(S) x 5 / 6.02 dB (4 -> 5.0 dB, 8 -> 7.5 dB). FEC: the S=2 code buys 4 dB against
# interference for 5 dB of sensitivity, 4 / 5 of its AWGN gain; the one BCH gain PARAMETERS gives is the ETSI
# informative sensitivity of SmartBAN at the 1 Msym/s channel rate (3.5, [44] Table 13 with 3.5's rate correction:
# Eb/N0 8.8 against 7.4 dB, 1.4 dB), so a BCH block buys 0.8 x 1.4 = 1.12 dB. Applied to every t = 2 BCH block (the
# 802.15.6 NB BCH(63,51) PSDU and BCH(31,19) header, the SmartBAN BCH(36,22) header) as the analogue (A).
REP_CI_SLOPE = 5 / (10 * np.log10(4))
BCH_CI_GAIN_DB = 0.8 * 1.4
# P(Walsh-16 symbol decoded | e chip errors), e = 0..6+: dmin 8, all 3-error patterns corrected, 4 and 5 errors on
# average 73.1 % and 32.7 % (exhaustive run of the HBCProtocol.py decoder, PARAMETERS 3.8)
WALSH_OK = np.array([1.0, 1.0, 1.0, 1.0, 0.731, 0.327, 0.0])
# bit errors of a wrongly decoded 16-ary symbol, the wrong symbol uniform over the other 15 labels:
# Hamming weights 1, 2, 3, 4 occur 4, 6, 4, 1 times (mean 32/15)
WALSH_BIT_W = np.array([4, 6, 4, 1]) / 15
_K16 = np.arange(2, 17)
_C16 = (-1.0) ** _K16 * comb(16, _K16)


def _nc_marcum(a, b):
    """Q1(a, b) - 0.5 exp(-(a^2 + b^2) / 2) I0(ab), the non-coherent correlated-FSK and DQPSK form (Proakis)."""
    return ncx2.sf(b * b, 2, a * a) - 0.5 * np.exp(-(a - b) ** 2 / 2) * i0e(a * b)


def ber(mod, gamma, h=0.5, n_comb=1):
    """Bit (chip) error probability at linear per-bit (per-chip) SNR gamma; 'oqpsk_dsss' takes the 2 MHz SINR.

    'gfsk', 'bfsk_chip': 0.5 exp(-g/2) (PARAMETERS 3.1 curve a, 3.8); 'bfsk_chip_v2': same with the +0.35 dB of the
    v2 simulator; 'ook_chip': 0.5 exp(-g/4); 'gfsk_marcum': correlated non-coherent FSK, |rho| = sinc(h) (3.1 curve b,
    3.6); 'oqpsk_dsss': 802.15.4 Annex E formula (3.2); 'dqpsk': pi/4-DQPSK Gray, g = Eb/N0 (3.4); 'dbpsk':
    n_comb-branch post-detection combining, g = total SNR over the branches, n_comb = 1 gives 0.5 exp(-g) (3.4).
    """
    g = np.asarray(gamma, dtype=float)
    if mod in ('gfsk', 'bfsk_chip'):
        p = 0.5 * np.exp(-g / 2)
    elif mod == 'bfsk_chip_v2':
        p = 0.5 * np.exp(-g * 10 ** (V2_OFFSET_DB / 10) / 2)
    elif mod == 'ook_chip':
        p = 0.5 * np.exp(-g / 4)   # textbook, not checked against the adaptive threshold (PARAMETERS 3.8, A)
    elif mod == 'gfsk_marcum':
        r = np.sqrt(1 - np.sinc(h) ** 2)   # np.sinc(h) = sin(pi h) / (pi h)
        p = _nc_marcum(np.sqrt(g / 2 * (1 - r)), np.sqrt(g / 2 * (1 + r)))
    elif mod == 'dqpsk':
        p = _nc_marcum(np.sqrt(g * (2 - np.sqrt(2))), np.sqrt(g * (2 + np.sqrt(2))))
    elif mod == 'dbpsk':
        L = int(n_comb)
        c = [sum(comb(2 * L - 1, n) for n in range(L - k)) / factorial(k) for k in range(L)]
        p = 2.0 ** (1 - 2 * L) * np.exp(-g) * np.polyval(c[::-1], g)
    elif mod == 'oqpsk_dsss':
        p = (8 / 15) / 16 * (_C16 * np.exp(20 * g[..., None] * (1 / _K16 - 1))).sum(-1)
    else:
        raise ValueError(f'unknown modulation {mod!r}')
    p = np.clip(p, 0.0, 0.5)
    p = np.where(p < 1e-300, 0.0, p)    # near the double floor scipy's binom.pmf overflows (seen at 1.7e-308)
    return p if p.ndim else float(p)


def gamma_from_sinr(sinr_lin, phy, rep=None, gain_db=0.0):
    """Per-bit (per-chip) SNR: SINR x B_noise / R x S x coded gain; 'oqpsk_dsss' keeps the SINR (its formula wants it).

    rep None uses phy['rep']; gain_db is the coded-PHY gain of a ('conv', gain_db) block (the layout carries it,
    because BLE Coded FEC block 1 always runs at S=8 while block 2 may run at S=2).
    """
    s = np.asarray(sinr_lin, dtype=float)
    g = s if phy['mod'] == 'oqpsk_dsss' else s * phy['noise_bw_hz'] / phy['bit_rate_bps']
    return g * (phy['rep'] if rep is None else rep) * 10 ** (gain_db / 10)


def sinr_at_ber(phy, p=1e-3):
    """Segment SINR (dB) at which the bit (chip) error probability of phy's modulation is p; p = 1e-3 is the BLE
    sensitivity criterion (0.1 % BER; common GFSK curve: 10.94 dB)."""
    def f(x):
        return ber(phy['mod'], gamma_from_sinr(10 ** (x / 10), phy), phy.get('h', 0.5), phy['rep']) - p
    return brentq(f, -40.0, 60.0)


def block_fail_prob(n, t, p):
    """Bounded-distance decoder failure of an (n, k, t) block: P(more than t errors among n) at independent p."""
    return binom.sf(t, n, p)


def walsh16_symbol_ok_prob(p):
    """P(Walsh-16 symbol decoded) at chip error p: sum_{e<=3} B(e) + 0.731 B(4) + 0.327 B(5), B = Binomial(16, p)."""
    p = np.asarray(p, dtype=float)
    r = (WALSH_OK[:6] * binom.pmf(np.arange(6), 16, p[..., None])).sum(-1)
    return r if r.ndim else float(r)


def ci_gain(block, phy):
    """Linear factor by which a block lowers the co-channel C/I floor ci_co_lin: its processing gain S, the
    block's repetition or spreading factor (rep, 802.15.6 NB S; HBC: the FSC SF over SF 8 via phy['ci_pg'], since SF 8
    is spectral shaping only, PARAMETERS 3.7), credited as S ** REP_CI_SLOPE, times BCH_CI_GAIN_DB for a BCH block
    (not the sync correlators ('bch', n, 0, t)). BLE Coded carries its specification C/I in ci_co_lin (12 dB S=8,
    17 dB S=2, PARAMETERS 5.5), which already holds its coding: ('conv', gain) blocks add nothing. phy['ci_credit']
    False (BCC full): no credit, the 0 dB chip rule holds for every block."""
    if not phy.get('ci_credit', True):
        return 1.0
    code = block['code'] or ('none',)
    fec = 10 ** (BCH_CI_GAIN_DB / 10) if code[0] == 'bch' and code[2] > 0 else 1.0
    return (block.get('rep', phy['rep']) * phy.get('ci_pg', 1.0)) ** REP_CI_SLOPE * fec


def decisions(block):
    """Bits a block carries before FEC, the basis of raw_bit_errors: coded bits of a BCH or uncoded block, bits after
    despreading (802.15.4, whose formula is per bit), 4 bits per 16-chip Walsh symbol (BCC, HBC: the chips are not
    bits and there is no FEC, so raw = post there), none for a 64-chip sync or pilot correlator (('bch', 64, 0, 7))."""
    code = block['code'] or ('none',)
    if code[0] == 'walsh16':
        return int(block['n_bits']) // 4
    return 0 if code[0] == 'bch' and code[2] == 0 else int(block['n_bits'])


def frame_outcome(rng, segments, layout, phy, sampled=True):
    """Outcome of one frame from the (t0_us, t1_us, sinr_lin, sir_lin) segments that cover it in time order.

    segments[0] starts at the frame start (absolute times are shifted to the layout's frame-relative times).
    Blocks: dict(t0_us, t1_us, kind 'sync'|'header'|'payload', n_bits, info_bits, code, optional rep). n_bits counts
    the hard decisions (coded bits after repetition combining, or chips), spread uniformly over [t0_us, t1_us).
    code None: all n_bits must be right; ('bch', n, k, t): codewords of n (the last one shortened) fail above t
    errors (('bch', 64, 0, 7) is the BCC 64-chip sync correlator); ('walsh16',): 16-chip symbols, WALSH_OK;
    ('conv', gain_db): SNR x gain, then an uncoded decision per bit (n_bits = info bits; stated approximation).
    Errors per codeword = sum over the segments it overlaps of Binomial(its decisions there, segment p). A segment
    with sir_lin below a block's floor ci_co_lin / ci_gain(block) fails the frame and the block's decisions in it are
    guesses (p = 0.5); a sync block needs
    SINR >= sync_sinr_db in every segment it overlaps. A failed codeword leaves Hypergeometric(k, n - k, e) info-bit
    errors (mean e k / n), a failed Walsh symbol 1..4 bit errors (mean 32/15). A missed pilot (a 'payload' block
    coded ('bch', 64, 0, 7)) shifts the rest of the frame (PARAMETERS 3.8): every later payload decision and info
    bit is a guess (Binomial(n, 0.5); expected mode weighted by the probability that a pilot before was missed).
    Returns success (frame delivered correct), sync_ok, header_ok, raw_bit_errors (payload bits in error before FEC,
    on the bit basis of decisions()), post_fec_bit_errors (payload info bits in error after decoding), payload_bits
    (info), payload_coded_bits (decisions() summed over the payload blocks), delivered_corrupt (sync and header
    passed, frame not correct, no payload CRC). sampled=False returns the probabilities and expected counts instead
    of one draw.
    """
    seg = np.asarray(segments, dtype=float).reshape(-1, 4)
    ts, te = seg[:, 0] - seg[0, 0], seg[:, 1] - seg[0, 0]
    ci_co = phy['ci_co_lin'] or 0.0
    floor_ok = np.ones(len(seg), dtype=bool)   # a segment clears the floor of every block it overlaps
    thr = -np.inf if phy['sync_sinr_db'] is None else 10 ** (phy['sync_sinr_db'] / 10)
    ok = {'sync': [1.0], 'header': [1.0], 'payload': [1.0]}   # per block: 0/1 (sampled) or probability
    raw = post = pay_bits = pay_coded = 0
    found = 1.0     # every pilot so far found: 0/1 (sampled) or probability
    guess = (lambda x: rng.binomial(int(x), 0.5)) if sampled else (lambda x: 0.5 * x)
    for b in layout:
        n, kind = int(b['n_bits']), b['kind']
        code = b['code'] or ('none',)
        over = (ts < b['t1_us']) & (te > b['t0_us'])
        ci_ok = seg[:, 3] >= ci_co / ci_gain(b, phy)
        floor_ok &= ci_ok | ~over
        if kind == 'sync':
            ok['sync'].append(float((seg[over, 2] >= thr).all()))
        if kind == 'payload':
            pay_bits += b['info_bits']
            pay_coded += decisions(b)
        if n == 0:
            continue
        rep = b.get('rep', phy['rep'])
        g = gamma_from_sinr(seg[:, 2], phy, rep, code[1] if code[0] == 'conv' else 0.0)
        p_seg = np.where(ci_ok, ber(phy['mod'], g, h=phy.get('h', 0.5), n_comb=rep), 0.5)
        # decisions per segment, then pieces that lie in one segment and one codeword.
        # ponytail: decisions evenly spread in time, BCH shortening on the last codeword only (NB spreads it),
        # no BCH miscorrection; exact per-codeword Binomial sums otherwise
        dur = b['t1_us'] - b['t0_us']
        e0, e1 = (np.clip(np.rint((x - b['t0_us']) / dur * n), 0, n).astype(int) for x in (ts, te))
        assert e0[0] == 0 and e1[-1] == n, 'segments must cover the frame'
        size = code[1] if code[0] == 'bch' else 16 if code[0] == 'walsh16' else n
        assert code[0] != 'walsh16' or n % 16 == 0, 'Walsh block of whole symbols'
        gedge = np.append(np.arange(0, n, size), n)
        gsize = np.diff(gedge)
        m = len(gsize)
        info = (np.maximum(gsize - (code[1] - code[2]), 0) if code[0] == 'bch'
                else np.full(m, 4) if code[0] == 'walsh16' else np.array([b['info_bits']]))
        tol = 5 if code[0] == 'walsh16' else code[3] if code[0] == 'bch' else 0
        bp = np.unique(np.concatenate((e0, e1, gedge)))
        start, plen = bp[:-1], np.diff(bp)
        pp = p_seg[np.searchsorted(e0, start, 'right') - 1]   # p of the segment holding each piece
        gi = np.searchsorted(gedge, start, 'right') - 1       # codeword holding each piece
        if sampled:
            e = np.bincount(gi, weights=rng.binomial(plen, pp), minlength=m).astype(int)
            if code[0] == 'walsh16':
                gok = rng.random(m) < WALSH_OK[np.minimum(e, 6)]
                res = (rng.choice(4, size=int((~gok).sum()), p=WALSH_BIT_W) + 1).sum()
            else:
                gok = e <= tol
                res = rng.hypergeometric(info[~gok], gsize[~gok] - info[~gok], e[~gok]).sum()
            blk_ok, blk_raw = float(gok.all()), int(e.sum())
        else:
            pmf = binom.pmf(np.arange(tol + 1), plen[:, None], pp[:, None])   # error pmf per piece, 0..tol
            acc = np.zeros((m, tol + 1))
            acc[:, 0] = 1.0
            for j in range(len(plen)):   # pmf of a codeword's errors truncated at tol (exact up to tol)
                acc[gi[j]] = np.convolve(acc[gi[j]], pmf[j])[:tol + 1]
            mean_e = np.bincount(gi, weights=plen * pp, minlength=m)
            if code[0] == 'walsh16':
                gok = acc @ WALSH_OK[:6]
                res = (1 - gok).sum() * 32 / 15
            else:
                gok = acc.sum(1)
                res = (info / gsize * (mean_e - acc @ np.arange(tol + 1))).sum()   # E[e 1{e > t}] k / n
            blk_ok, blk_raw = float(np.prod(gok)), float(mean_e.sum())
        ok[kind].append(blk_ok)
        if kind == 'payload':
            # raw on the bit basis of decisions(): demapped Walsh bits, nothing for a pilot correlator
            r = res if code[0] == 'walsh16' else blk_raw if decisions(b) else 0
            if code[0] in ('bch', 'walsh16') and 0 < b['info_bits'] < info.sum():
                # the codewords also carry MAC header and FCS bits: keep the residual errors that fall on the
                # application payload (info_bits of the block), so post_fec_bit_errors / payload_bits is a BER
                share = b['info_bits'] / info.sum()
                res = rng.binomial(int(res), share) if sampled else res * share
            if found < 1:   # after a missed pilot the demapper runs on shifted chips
                r = found * r + (1 - found) * guess(decisions(b))
                res = found * res + (1 - found) * guess(b['info_bits'])
            raw += r
            post += res
            if code[0] == 'bch' and code[2] == 0:   # a pilot correlator
                found *= blk_ok
    sync, hdr = float(np.prod(ok['sync'])), float(np.prod(ok['header']))
    succ = float(floor_ok.all()) * sync * hdr * float(np.prod(ok['payload']))
    corrupt = 0.0 if phy['payload_crc'] else max(sync * hdr - succ, 0.0)
    cast, count = (bool, int) if sampled else (float, float)
    return dict(success=cast(succ), sync_ok=cast(sync), header_ok=cast(hdr), raw_bit_errors=count(raw),
                post_fec_bit_errors=count(post), payload_bits=pay_bits, payload_coded_bits=pay_coded,
                delivered_corrupt=cast(corrupt))

"""Plain-assert tests of phy.py, no pytest needed.

Purpose: the BER formulas of PARAMETERS.md 3.x at hand-derived points, BCH failure against the binomial tail, the
Walsh-16 symbol and the BCC packet check values of PARAMETERS 3.8, and frame_outcome: sampled draws against the
expected values, the C/I floor, the sync threshold, repetition S=4 as 6.02 dB, the coded-PHY gain shift and the
missed pilot (PARAMETERS 3.8: the rest of the frame is guessed).
Every expected value is derived by hand (or read in PARAMETERS.md) in the comment above its assert.
Inputs:  phy.py.
Outputs: one line per passed test on stdout; an AssertionError on failure.
argv:    none (python python_simulation/coexistence/test_phy.py), runtime about ten seconds.
"""
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parent))
import phy   # noqa: E402

DB = lambda x: 10 ** (x / 10)   # noqa: E731


def blk(t0, t1, kind, n, info, code=None, **kw):
    return dict(t0_us=t0, t1_us=t1, kind=kind, n_bits=n, info_bits=info, code=code, **kw)


def mk_phy(mod, bw, rate, rep=1, crc=True, ci_db=None, sync_db=None):
    return dict(mod=mod, noise_bw_hz=bw, bit_rate_bps=rate, rep=rep, payload_crc=crc,
                ci_co_lin=None if ci_db is None else DB(ci_db), sync_sinr_db=sync_db)


# BCC frame of PARAMETERS 3.8 at L = 255, 1 chip = 1 us: 3 preamble repetitions not in the PSR formula
# [0, 192), 4th preamble and SFD as 64-chip correlators (<= 7 errors), header 8 Walsh symbols, PSDU 8160 chips
# with a 64-chip pilot after each full 2048 PSDU chips (3 pilots): 448 + 8160 + 192 = 8800 chips.
BCC_255 = [blk(192, 256, 'sync', 64, 0, ('bch', 64, 0, 7)), blk(256, 320, 'sync', 64, 0, ('bch', 64, 0, 7)),
           blk(320, 448, 'header', 128, 32, ('walsh16',))]
for _t in (448, 2560, 4672):
    BCC_255 += [blk(_t, _t + 2048, 'payload', 2048, 512, ('walsh16',)),
                blk(_t + 2048, _t + 2112, 'sync', 64, 0, ('bch', 64, 0, 7))]
BCC_255.append(blk(6784, 8800, 'payload', 2016, 504, ('walsh16',)))
BCC_PHY = mk_phy('bfsk_chip', 1e6, 1e6, crc=False)   # B = Rc: the segment SINR is Ec/N0


def test_ber_points():
    # non-coherent BFSK / common GFSK curve at 10 dB: 0.5 e^-5 = 0.5 x 6.7379e-3 = 3.3690e-3
    for mod in ('gfsk', 'bfsk_chip'):
        assert abs(phy.ber(mod, DB(10)) - 3.3690e-3) < 1e-7
    # v2 offset: the curve shifted by 0.35 dB, so 9.65 dB on 'bfsk_chip_v2' is 10 dB on 'bfsk_chip'
    assert abs(phy.ber('bfsk_chip_v2', DB(9.65)) - 3.3690e-3) < 1e-7
    # OOK chips: 0.5 exp(-g/4) at 10 dB = 0.5 e^-2.5 = 0.5 x 0.082085 = 0.041042
    assert abs(phy.ber('ook_chip', DB(10)) - 0.041042) < 1e-6
    # common GFSK curve at BER 1e-3: 0.5 exp(-g/2) = 1e-3 -> g = 2 ln 500 = 12.4292 = 10.944 dB (PARAMETERS 10.94)
    assert abs(phy.ber('gfsk', 2 * np.log(500)) - 1e-3) < 1e-12
    # DBPSK: 0.5 exp(-g) at g = 3: 0.5 e^-3 = 0.024894. Two-branch combining (L = 2): c0 = C(3,0) + C(3,1) = 4,
    # c1 = C(3,0) / 1! = 1, Pb = 2^-3 e^-g (4 + g); at g = 4: (1/8) e^-4 x 8 = e^-4 = 0.0183156
    assert abs(phy.ber('dbpsk', 3.0) - 0.024894) < 1e-6
    assert abs(phy.ber('dbpsk', 4.0, n_comb=2) - 0.0183156) < 1e-7
    # correlated non-coherent FSK at h = 1: rho = sin(pi)/pi = 0, so a = 0, b = sqrt(g), Q1(0, b) = exp(-b^2/2) and
    # I0(0) = 1: Pb = e^-g/2 - 0.5 e^-g/2 = 0.5 e^-g/2, the orthogonal curve (3.3690e-3 at 10 dB)
    assert abs(phy.ber('gfsk_marcum', DB(10), h=1.0) - 3.3690e-3) < 1e-7
    # PARAMETERS 3.1 curve b: BER 1e-3 at 14.31 dB for h = 0.5; 3.6: 17.80 dB for h = 0.32 (both rounded to 0.01 dB)
    for h, db in ((0.5, 14.31), (0.32, 17.80)):
        x = brentq(lambda d: np.log10(phy.ber('gfsk_marcum', DB(d), h=h)) + 3, 5, 25)
        assert abs(x - db) < 0.011, (h, x)
    # DQPSK at g = 0: a = b = 0, Q1(0, 0) = 1, minus 0.5 I0(0) = 0.5: a coin flip
    assert abs(phy.ber('dqpsk', 0.0) - 0.5) < 1e-12
    # arrays in, arrays out
    assert phy.ber('gfsk', np.array([DB(10), DB(10)])).shape == (2,)
    print('ber points ok')


def test_oqpsk_dsss():
    # Pb = (8/15)(1/16) sum_{k=2..16} (-1)^k C(16,k) exp(20 SINR (1/k - 1)) = (1/30) sum.
    # SINR = 0: every exp is 1; sum_{k=0..16} (-1)^k C(16,k) = (1-1)^16 = 0, so sum_{k=2..16} = 0 - 1 + 16 = 15,
    # Pb = 15/30 = 0.5.
    assert abs(phy.ber('oqpsk_dsss', 0.0) - 0.5) < 1e-12
    # SINR = 1 (0 dB), terms C(16,k) exp(20/k - 20), k = 2..16:
    # +120 e^-10 = 5.44799e-3, -560 e^-13.333 = -9.06974e-4, +1820 e^-15 = 5.56742e-4, -4368 e^-16 = -4.91554e-4,
    # +8008 e^-16.667 = 4.62682e-4, -11440 e^-17.143 = -4.10561e-4, +12870 e^-17.5 = 3.23166e-4,
    # -11440 e^-17.778 = -2.17588e-4, +8008 e^-18 = 1.21962e-4, -4368 e^-18.182 = -5.54650e-5,
    # +1820 e^-18.333 = 1.98612e-5, -560 e^-18.462 = -5.37581e-6, +120 e^-18.571 = 1.03208e-6,
    # -16 e^-18.667 = -1.25109e-7, +1 e^-18.75 = 7.19413e-9; sum = 4.84580e-3; Pb = 4.84580e-3 / 30 = 1.61527e-4
    # (PARAMETERS 3.2: 1.6e-4 at 0 dB)
    assert abs(phy.ber('oqpsk_dsss', 1.0) - 1.61527e-4) < 2e-9
    print('802.15.4 formula ok')


def test_gamma_from_sinr():
    # BLE 2M: SINR in 2 MHz, 2 Mb/s -> Eb/N0 = SINR x 2e6 / 2e6 = SINR
    assert phy.gamma_from_sinr(7.0, mk_phy('gfsk', 2e6, 2e6)) == 7.0
    # 802.15.6 NB pi/4-DQPSK: SINR in the 600 kHz symbol band, 1.2 Mb/s raw -> gb = gs / 2 (PARAMETERS 3.4)
    assert abs(phy.gamma_from_sinr(7.0, mk_phy('dqpsk', 600e3, 1.2e6)) - 3.5) < 1e-12
    # 802.15.4: the formula takes the 2 MHz SINR itself, not 8 SINR (PARAMETERS 3.2, corrected)
    assert phy.gamma_from_sinr(7.0, mk_phy('oqpsk_dsss', 2e6, 250e3)) == 7.0
    # repetition S = 4 is a 6.02 dB gain: 10 log10 4 = 6.0206 dB
    p1, p4 = mk_phy('gfsk', 1e6, 1e6, rep=1), mk_phy('gfsk', 1e6, 1e6, rep=4)
    assert abs(phy.gamma_from_sinr(2.0, p4) - phy.gamma_from_sinr(2.0 * DB(6.0206), p1)) < 1e-4
    print('gamma_from_sinr ok')


def test_block_codes():
    # BCH(63,51,t=2) at p = 1e-2: 0.99^63 = exp(63 ln 0.99) = exp(-0.633172) = 0.530906; 63 x 0.01 x 0.99^62 =
    # 0.337849; C(63,2) = 1953, 1953 x 1e-4 x 0.99^61 = 0.105791; P(<= 2) = 0.974546, failure 0.025454
    assert abs(phy.block_fail_prob(63, 2, 1e-2) - 0.025454) < 2e-6
    # BCH(127,113,t=2) at p = 1e-2: 0.99^127 = exp(-1.276393) = 0.279042; 127 x 0.01 x 0.99^126 = 0.357963;
    # C(127,2) = 8001, 8001 x 1e-4 x 0.99^125 = 0.227795; P(<= 2) = 0.864800, failure 0.135200
    assert abs(phy.block_fail_prob(127, 2, 1e-2) - 0.135200) < 2e-6
    # Walsh-16 at p = 5.5 %, q = 0.945: B(0) = q^16 = 0.404491, B(1) = 16 p q^15 = 0.376669, B(2) = 120 p^2 q^14 =
    # 0.164419, B(3) = 560 p^3 q^13 = 0.044657, B(4) = 1820 p^4 q^12 = 0.008447, B(5) = 4368 p^5 q^11 = 0.001180;
    # ok = 0.990236 + 0.731 x 0.008447 + 0.327 x 0.001180 = 0.996797
    assert abs(phy.walsh16_symbol_ok_prob(0.055) - 0.996797) < 2e-6
    # p = 1.5 %: B(4) = 1820 x 5.0625e-8 x 0.985^12 = 7.6855e-5, B(5) = 4368 x 7.59375e-10 x 0.985^11 = 2.8089e-6,
    # B(>= 6) = 7.9e-8; 1 - ok = 0.269 B(4) + 0.673 B(5) + B(>= 6) = 2.0674e-5 + 1.890e-6 + 7.9e-8 = 2.264e-5
    assert abs(phy.walsh16_symbol_ok_prob(0.015) - (1 - 2.264e-5)) < 2e-8
    # packet PSR = Q64^(2 + P) x ok^(8 + 2L), L = 255, P = floor(32 x 255 / 2048) = 3: PARAMETERS 3.8 check values
    # 0.168 at p = 5.5 % (measured 2/12) and 0.99 at p = 1.5 % (12/12); Q64 = P(<= 7 errors of 64)
    for p, psr, tol in ((0.055, 0.168, 1e-3), (0.015, 0.99, 5e-3)):
        got = (1 - phy.block_fail_prob(64, 7, p)) ** 5 * phy.walsh16_symbol_ok_prob(p) ** 518
        assert abs(got - psr) < tol, (p, got)
        # the same through frame_outcome on the chip layout; p = 0.5 exp(-g/2) -> g = -2 ln(2p)
        seg = [(0, 8800, -2 * np.log(2 * p), np.inf)]
        out = phy.frame_outcome(None, seg, BCC_255, BCC_PHY, sampled=False)
        assert abs(out['success'] - got) < 1e-9, (out['success'], got)
        # raw BER on the bit basis: 8160 PSDU chips = 510 Walsh symbols = 2040 bits, pilots carry none; the stack has
        # no FEC, so raw = post = the bit errors of the wrong symbols, 510 x (1 - ok) x 32/15 (WALSH_BIT_W mean)
        assert out['payload_bits'] == 2040 and out['payload_coded_bits'] == 2040
        want = 510 * (1 - phy.walsh16_symbol_ok_prob(p)) * 32 / 15
        assert abs(out['raw_bit_errors'] - want) < 1e-9 and abs(out['post_fec_bit_errors'] - want) < 1e-9, out
    print('BCH, Walsh-16 and BCC packet ok')


def agree(rng, segs, layout, ph, n):
    """Mean of n sampled outcomes within 4.5 standard errors of the expected values."""
    exp = phy.frame_outcome(None, segs, layout, ph, sampled=False)
    draws = [phy.frame_outcome(rng, segs, layout, ph) for _ in range(n)]
    for key in ('success', 'sync_ok', 'header_ok', 'delivered_corrupt', 'raw_bit_errors', 'post_fec_bit_errors'):
        x = np.array([d[key] for d in draws], dtype=float)
        # Bernoulli sd for probabilities; for counts the sample sd, floored at sqrt(mean) (bursty counts have
        # variance >= mean), so a rare count that never occurred in n draws still passes
        sd = np.sqrt(exp[key] * (1 - exp[key])) if isinstance(draws[0][key], bool) else max(x.std(), np.sqrt(exp[key]))
        assert abs(x.mean() - exp[key]) <= 4.5 * sd / np.sqrt(n) + 1e-12, (key, x.mean(), exp[key], sd)
    return exp


def test_sampled_vs_expected():
    rng = np.random.default_rng(48)
    # 802.15.6 NB-like frame, pi/2-DBPSK 600 kb/s: 90 sync bits, BCH(31,19) header at S = 4 (4-branch combining),
    # PSDU of 166 coded bits = BCH(63,51) x 2 + a codeword shortened to 40 (28 info bits); the segment edge at
    # 450 us cuts the first PSDU codeword (bits 0..62 over 356.7..461.7 us)
    ph = mk_phy('dbpsk', 600e3, 600e3)
    lay = [blk(0, 150, 'sync', 90, 0), blk(150, 356.7, 'header', 31, 19, ('bch', 31, 19, 2), rep=4),
           blk(356.7, 633.3, 'payload', 166, 130, ('bch', 63, 51, 2))]
    segs = [(0, 200, 6.5, np.inf), (200, 450, 4.0, 30.0), (450, 633.3, 3.2, np.inf)]
    e = agree(rng, segs, lay, ph, 20000)
    assert 0.2 < e['success'] < 0.9, e   # a mid-range case, so the comparison has power
    # BCC-like frame without payload CRC, L = 20: delivered_corrupt must follow too
    lay = [blk(0, 64, 'sync', 64, 0, ('bch', 64, 0, 7)), blk(64, 192, 'header', 128, 32, ('walsh16',)),
           blk(192, 832, 'payload', 640, 160, ('walsh16',))]
    segs = [(1000, 1300, DB(6.0), np.inf), (1300, 1832, DB(7.0), np.inf)]   # absolute times: shifted to 0
    e = agree(rng, segs, lay, BCC_PHY, 20000)
    assert 0.2 < e['success'] < 0.9 and e['delivered_corrupt'] > 0.05, e
    print('sampled vs expected ok')


def test_ci_floor_and_sync():
    rng = np.random.default_rng(7)
    # BLE 1M-like frame: preamble + AA 40 bits, then header + 20-octet payload + CRC = 200 bits, 1 us per bit
    ph = mk_phy('gfsk', 1e6, 1e6, ci_db=21.0)
    lay = [blk(0, 40, 'sync', 40, 0), blk(40, 240, 'payload', 200, 160)]
    # SINR 30 to 40 dB: 0.5 exp(-500) = 0 per bit, so only the floor can fail the frame
    bad = [(0, 100, 1e4, np.inf), (100, 150, 1e3, DB(20.0)), (150, 240, 1e4, np.inf)]
    good = [(0, 100, 1e4, np.inf), (100, 150, 1e3, DB(25.0)), (150, 240, 1e4, np.inf)]
    assert phy.frame_outcome(None, bad, lay, ph, sampled=False)['success'] == 0.0
    assert not any(phy.frame_outcome(rng, bad, lay, ph)['success'] for _ in range(200))
    assert phy.frame_outcome(None, good, lay, ph, sampled=False)['success'] == 1.0
    # the 50 payload bits under SIR 20 dB < 21 dB are guesses: expected raw errors 50 x 0.5 = 25
    assert phy.frame_outcome(None, bad, lay, ph, sampled=False)['raw_bit_errors'] == 25.0
    # sync threshold: 802.15.4 SHR (40 bits) partly at SINR 2 dB, where Pb is about 5e-7 (PARAMETERS 3.2: 1e-6 at
    # 1.81 dB); a 3 dB threshold loses the frame, a 1 dB threshold keeps it (40 x 5e-7 = 2e-5 < 1e-3 loss)
    lay = [blk(0, 160, 'sync', 40, 0), blk(160, 192, 'header', 8, 8), blk(192, 832, 'payload', 160, 160)]
    segs = [(0, 100, DB(2.0), np.inf), (100, 832, DB(10.0), np.inf)]
    for thr, want in ((3.0, 0.0), (1.0, 1.0)):
        out = phy.frame_outcome(None, segs, lay, mk_phy('oqpsk_dsss', 2e6, 250e3, sync_db=thr), sampled=False)
        assert abs(out['sync_ok'] - want) < 1e-3 and abs(out['success'] - want) < 1e-3, (thr, out)
    print('C/I floor and sync threshold ok')


def test_repetition_and_coded_gain():
    lay = [blk(0, 200, 'payload', 200, 200)]
    # repetition S = 4 at SINR x equals S = 1 at x + 6.0206 dB (ideal combining on the GFSK curve)
    s4 = phy.frame_outcome(None, [(0, 200, 2.48, np.inf)], lay, mk_phy('gfsk', 1e6, 1e6, rep=4), sampled=False)
    s1 = phy.frame_outcome(None, [(0, 200, 2.48 * DB(6.0206), np.inf)], lay, mk_phy('gfsk', 1e6, 1e6), sampled=False)
    assert abs(s4['success'] - s1['success']) < 1e-6 and 0.2 < s1['success'] < 0.9, (s4, s1)
    # coded-PHY gain: a ('conv', 7.99) block at SINR 10 / 10^0.799 sees g = 10 per bit: Pb = 0.5 e^-5 = 3.36897e-3,
    # 200 bits all right: (1 - 3.36897e-3)^200 = exp(200 x -3.374661e-3) = exp(-0.674932) = 0.509191
    lay = [blk(0, 200, 'payload', 200, 200, ('conv', 7.99))]
    out = phy.frame_outcome(None, [(0, 200, 10 / DB(7.99), np.inf)], lay, mk_phy('gfsk', 1e6, 1e6), sampled=False)
    assert abs(out['success'] - 0.509191) < 1e-6, out
    # expected raw errors 200 x 3.36897e-3 = 0.673795, uncoded decisions: post-FEC errors the same
    assert abs(out['raw_bit_errors'] - 0.673795) < 1e-6 and abs(out['post_fec_bit_errors'] - 0.673795) < 1e-6
    print('repetition and coded gain ok')


def test_post_fec_counts_payload_only():
    # BCH(63,51) payload block of 120 decisions: 96 info positions in 2 codewords (51, and 45 after shortening), of
    # which 24 are application bits (3 octets; the rest MAC header and FCS). DBPSK at SINR 0: p = 0.5 e^0 = 0.5, a
    # coin flip per decision, and a t = 2 decoder fails but for P(e <= 2 of 63) ~ 2e-16. Expected raw errors 60;
    # residual info errors 96 x 0.5 = 48, of which the application share 24 / 96 gives 12: post / payload_bits = 0.5
    lay = [blk(0, 200, 'payload', 120, 24, ('bch', 63, 51, 2))]
    out = phy.frame_outcome(None, [(0, 200, 0.0, np.inf)], lay, mk_phy('dbpsk', 6e5, 6e5), sampled=False)
    assert abs(out['raw_bit_errors'] - 60) < 1e-9 and abs(out['post_fec_bit_errors'] - 12) < 1e-9, out
    # 0.5 e^-707 = 4.5e-308 sits below the 1e-300 floor and is flushed to 0 (scipy's binom.pmf overflowed there),
    # so the expected mode at that SNR succeeds with certainty instead of raising
    assert phy.ber('dbpsk', 707.0) == 0.0
    out = phy.frame_outcome(None, [(0, 200, 707.0, np.inf)], lay, mk_phy('dbpsk', 6e5, 6e5), sampled=False)
    assert out['success'] == 1.0
    print('post-FEC errors on the payload ok')


def test_missed_pilot():
    # PARAMETERS 3.8: 'a missed pilot shifts the rest'. BCC-like frame, 1 chip = 1 us (Ec/N0 = SINR): sync correlator,
    # 8-symbol header, 128 Walsh symbols (512 bits), a pilot inside the payload, 40 symbols (160 bits). SINR 1e4 gives
    # p = 0.5 e^-5000 = 0 everywhere but the pilot, which sits at SINR 0 (linear): p = 0.5 e^0 = 0.5, so it is found
    # with P(<= 7 errors of 64) = (1 + 64 + 2016 + 41664 + 635376 + 7624512 + 74974368 + 621216192) / 2^64 =
    # 704494193 / 2^64 = 3.8e-11.
    # The 160 bits after it are guesses: raw = post = 0.5 x 160 = 80 expected, the 512 before it are clean; the
    # sync and header passed and there is no payload CRC, so the frame is delivered corrupt.
    lay = [blk(0, 64, 'sync', 64, 0, ('bch', 64, 0, 7)), blk(64, 192, 'header', 128, 32, ('walsh16',)),
           blk(192, 2240, 'payload', 2048, 512, ('walsh16',)), blk(2240, 2304, 'payload', 64, 0, ('bch', 64, 0, 7)),
           blk(2304, 2944, 'payload', 640, 160, ('walsh16',))]
    segs = [(0, 2240, 1e4, np.inf), (2240, 2304, 0.0, np.inf), (2304, 2944, 1e4, np.inf)]
    out = phy.frame_outcome(None, segs, lay, BCC_PHY, sampled=False)
    assert abs(out['raw_bit_errors'] - 80) < 1e-6 and abs(out['post_fec_bit_errors'] - 80) < 1e-6, out
    assert out['success'] < 1e-10 and out['delivered_corrupt'] > 1 - 1e-10, out
    # sampled draws follow the expected values when the pilot is found about half the time: chip error 0.1 on the
    # pilot (g = -2 ln 0.2), 0.05 elsewhere (g = -2 ln 0.1), so the Walsh errors and the guesses both count
    segs = [(0, 2240, -2 * np.log(0.1), np.inf), (2240, 2304, -2 * np.log(0.2), np.inf),
            (2304, 2944, -2 * np.log(0.1), np.inf)]
    e = agree(np.random.default_rng(3), segs, lay, BCC_PHY, 10000)
    assert 0.2 < e['delivered_corrupt'] < 0.9, e
    print('missed pilot ok')


if __name__ == '__main__':
    test_ber_points()
    test_oqpsk_dsss()
    test_gamma_from_sinr()
    test_block_codes()
    test_sampled_vs_expected()
    test_ci_floor_and_sync()
    test_repetition_and_coded_gain()
    test_post_fec_counts_payload_only()
    test_missed_pilot()
    print('all phy tests passed')

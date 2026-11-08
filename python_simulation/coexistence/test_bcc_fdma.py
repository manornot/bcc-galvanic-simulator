"""Plain-assert tests of channel_curves.py and bcc_fdma.py, no pytest needed.

Purpose: the equivalence checks (a), (b2) and (c) of bcc_fdma at reduced size, packets at Eb/N0 = 8 dB on one tone,
far-spaced tones against a single tone, the chips per payload bit, the tone placement of plan(), the one-tone placement
of channel_curves.tone_centre and the real-passband tone near 0 Hz, the usable-band widths of
channel_curves.csv, and the minimum-phase channel of (RC and RC-ladder reconstruction, the curve variants, a pure
delay that the receiver's timing absorbs).
Every expected value is derived by hand in the comment above its assert.
Inputs:  channel_curves.csv, the repository HBCProtocol and PHY2 modules.
Outputs: one line per passed test on stdout; an AssertionError on failure.
argv:    none (python python_simulation/coexistence/test_bcc_fdma.py), runtime well under 2 minutes.
"""
import sys
import time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bcc_fdma as b          # noqa: E402
import channel_curves as cc   # noqa: E402


def test_usable_band_widths():
    # 3 dB rule on the grid 50k ... 20 MHz, bin edges at geometric midpoints clipped to 0.6 and 20 MHz, as in
    # an earlier script (fdma_ceiling). Edges used below (MHz):
    # sqrt(1*2)=1.414, sqrt(2*3)=2.449, sqrt(6*7)=6.481, sqrt(7*8)=7.483, sqrt(8*9)=8.485, sqrt(9*10)=9.487,
    # sqrt(10*12.5)=11.180, sqrt(12.5*15)=13.693, sqrt(15*20)=17.321, top 20 (clipped).
    # cohort 50 Ohm, bins 7, 10, 12.5, 15, 20: (7.483-6.481)+(11.180-9.487)+(20-11.180) = 1.002+1.693+8.820 = 11.5
    # cohort 1 MOhm, bins 2..7, 9, 10, 12.5, 15: (7.483-1.414)+(11.180-8.485)+(17.321-11.180) = 6.069+2.695+6.141 = 14.9
    # RHW1-LHW1 50 Ohm, bins 3..7, 10..20: (7.483-2.449)+(20-9.487) = 5.034+10.513 = 15.5
    # RHW1-LHW1 1 MOhm, bins 2..7, 10, 12.5: 6.069+(13.693-9.487) = 6.069+4.206 = 10.3
    # (the published widths of the 2026-09-17 script: 11.5, 14.9, 15.5, 10.3 MHz)
    for curve, load, width in (("cohort", "50", 11.5), ("cohort", "1M", 14.9), ("RHW1-LHW1", "50", 15.5),
                               ("RHW1-LHW1", "1M", 10.3)):
        got = cc.usable_band(*cc.read_curve(curve, load))[1] / 1e6
        assert abs(got - width) < 0.05, (curve, load, got, width)


def test_chips_per_payload_bit():
    # preamble 4 x 64 = 256, SFD 64, header 32 bits x 4 chips = 128, PSDU 2040 bits x 4 = 8160 chips with a
    # 64-chip pilot after each full 2048-chip chunk (8160 // 2048 = 3 pilots, 192 chips): 256+64+128+8160+192 = 8800
    assert b.CHIPS_PER_PACKET == 8800, b.CHIPS_PER_PACKET
    assert abs(b.CPB - 8800 / 2040) < 1e-12 and abs(b.CPB - 4.314) < 1e-3


def test_single_tone_ber():
    # (a) reduced: non-coherent BFSK theory 0.5 exp(-Eb/(2 N0)): 6 dB -> Eb/N0 = 3.981 -> 0.5 exp(-1.991) = 0.0683;
    # 8 dB -> 6.310 -> 0.5 exp(-3.155) = 0.0213. The simulated Eb/N0 offset from theory must be within 0.5 dB
    # (PHY2 itself is 0.3 to 0.4 dB better than theory, v2_ber_vs_theory.md), and the
    # baseband error count must match PHY2 on the same chips within 3 binomial standard errors
    # (40 000 chips at p = 0.016: one error count has sd sqrt(40000 x 0.016) = 25, the difference sd 36).
    rows = b.check_a(40000, (6.0, 8.0))
    for r in rows:
        assert abs(r["offset_db"]) < 0.5, r
        assert abs(r["z"]) < 3.0, r
    assert abs(rows[0]["theory"] - 0.0683) < 1e-3 and abs(rows[1]["theory"] - 0.0213) < 1e-3


def test_neighbour_count():
    # (c) reduced: centre of 9 tones at 4 R_c spacing, 8 dB. Neighbours 3 and 4 sit 12 and 16 R_c away, their
    # nearest tone component 11 R_c from the victim centre, 9 R_c beyond the 2 R_c filter edge, so adding them
    # (same noise, same samples per chip) must change the error count by less than one binomial standard error
    # sqrt(n p (1 - p)) ~= sqrt(30000 x 0.016) = 22 errors.
    r = b.check_c(30000, "fsk", 8.0)[0]
    sd = np.sqrt(r["chips"] * r["cer_ref"] * (1 - r["cer_ref"]))
    assert r["errors_ref"] > 100, r
    assert abs(r["errors_bb"] - r["errors_ref"]) < sd, r


def test_single_tone_packets_8db():
    # One tone, 3 packets of the 3-packet chip file layout (3 x (1000 + 8800 + 1000) = 32 400 chips), Eb/N0 8 dB.
    # Chip error ~1.6 % (theory 2.1 % less the 0.4 dB filter gain of PHY2). A 16-chip Walsh symbol (minimum distance
    # 8) needs >= 4 chip errors to fail: C(16,4) 0.016^4 0.984^12 = 1820 x 6.6e-8 x 0.82 = 9.8e-5, and only about a
    # third of those land on a wrong codeword, so ~3e-5 per symbol, 518 symbols per packet -> ~1.6 % packet error
    # (the ACI series measures 2 to 3 %). So 8 dB is NOT zero packet error: with p = 0.03 per packet, 2 or 3 of 3
    # lost has probability 3 p^2 (1 - p) + p^3 = 0.0026, so at least 2 of 3 must be error-free (any seed). The 256-chip
    # preamble and 64-chip SFD are found at this chip error rate, so all 3 packets are detected.
    g = b.Geometry("test", "flat", "", 50e3, 0.0, sync=True, n_packets=3, n_tones=1)
    r = b.run_victim((g, 0, b.settings(targets=(8.0,))))[0]
    assert r["chips"] == 32400, r["chips"]
    assert abs(r["ebn0_db"] - 8.0) < 1e-9
    assert 0.005 < r["cer"] < 0.03, r["cer"]
    assert r["detected"] == 3 and r["decoded"] >= r["error_free"] >= 2, r


def test_plan_tone_placement():
    # cohort 50 Ohm, 50 kchip/s, guard 0: band lo = sqrt(6 x 7) = 6.4807 MHz, hi = 20 MHz, span 13.5193 MHz, spacing
    # 4 x 50 k = 0.2 MHz -> floor(67.60) = 67 tones (68 x 0.2 = 13.6 > 13.519 does not fit). Centred:
    # f[0] = 13.2404 - 33 x 0.2 = 6.6404 MHz, f[0] - 2 R_c = 6.5404 >= 6.4807; f[-1] + 2 R_c = 19.9404 <= 20.
    # In-band (3 dB bins [6.4807, 7.4833] and [9.4868, 20] MHz, tone i spans 6.5404 + 0.2 i to 6.7404 + 0.2 i):
    # i = 0..3 in the first (6.7404 + 0.6 = 7.3404 <= 7.4833, i = 4 gives 7.5404), i >= 15 in the second
    # (6.5404 + 3.0 = 9.5404 >= 9.4868, i = 14 gives 9.3404), i.e. 4 + 52 = 56 tones.
    p = b.plan(b.Geometry("test", "cohort", "50", 50e3, 0.0))
    assert len(p.f) == 67, len(p.f)
    assert abs(p.lo - 6.4807e6) < 100 and p.hi == 20e6, (p.lo, p.hi)
    assert p.f[0] - 100e3 >= p.lo and p.f[-1] + 100e3 <= p.hi, (p.f[0], p.f[-1])
    assert abs(p.f[0] - 6.6404e6) < 100 and np.allclose(np.diff(p.f), 200e3)
    assert p.in_band == 56, p.in_band


def test_tone_centre():
    # one-tone placement on the cohort 50 Ohm curve (csv: 5 MHz -45.8622, 6 -46.0428, 7 -45.2582, 9 -46.2123,
    # 10 -43.6213, 12.5 -42.3660, 15 -42.8687, 20 -44.7671 dB), f a multiple of 1 kHz with f -/+ R_c in [0.05, 20] MHz:
    # - 7 Mchip/s: f in [7.05, 13] MHz, lower tone f - 7 in [0.05, 6], upper f + 7 in [14.05, 20]. Below 6 MHz the
    #   curve peaks at its 5 MHz point (-45.86; 4 MHz -48.06, 6 MHz -46.04, log-f segments monotonic), while the upper
    #   tone never falls below -44.77 (20 MHz): the lower tone is the weaker everywhere, best at 5 MHz: f = 12.000 MHz
    # - 5.25 Mchip/s: f in [5.30, 14.75], lower tone in [0.05, 9.5] MHz. Its best gain is at 9.5 MHz, -46.2123 +
    #   2.591 log(9.5 / 9) / log(10 / 9) = -44.88 (7 MHz gives -45.26), still below the upper tone's -44.77 at 20 MHz:
    #   f = 20 - 5.25 = 14.750 MHz
    # - 9.95 Mchip/s: f in [10.00, 10.05], lower tone 50 to 100 kHz (-77.6 to -73.6 dB) against about -44.8 at 20 MHz:
    #   f = 10.050 MHz; 10.5 Mchip/s: 2 R_c = 21 MHz exceeds the 19.95 MHz span, no placement
    assert [cc.tone_centre(r) for r in (7e6, 5.25e6, 9.95e6, 10.5e6)] == [12.0e6, 14.75e6, 10.05e6, None]


def test_off_grid_passband():
    # (b2) reduced, BFSK: 3 tones 0.37 R_c off the integer grid with offsets 7/32, 21/32, 29/32 chip, an independent
    # continuous-phase passband transmitter against the baseband model on the same chips. The two must agree within
    # 3 standard errors of the difference on every tone. A carrier-phase term that changes from chip to chip (f_j in
    # place of f_j - f_v, the PHY2 per-chip restart artefact) fails this with z = 9 to 38 at this size.
    for r in b.check_b2(30000, 2, "fsk"):
        assert r["errors_ref"] > 50 and abs(r["z"]) < 3.0, r


def test_passband_near_dc():
    # One BFSK tone in real passband against the complex-baseband model on the same chips,
    # timing and noise seed, tone at channel_curves.tone_centre, 0 dB gains:
    # - 5.25 Mchip/s, 8 dB, 2 x 100 000 chips, the highest rate run_victim keeps on the baseband (f - 2 R_c = 4.25 MHz
    #   > 0). The real band-pass also passes the lower tone's image at -9.5 MHz, which the f - R_c correlator sees
    #   2 x 9.5 / 5.25 = 3.62 cycles per chip off: amplitude |sin(3.62 pi)| / (3.62 pi) = 0.08, energy (1 -/+ 0.08)^2,
    #   -/+ 0.7 dB with the carrier phase, a mean loss well below 0.2 dB: the Eb/N0 offsets from theory agree to 0.2 dB
    # - 9.95 Mchip/s, 30 dB on both tones, 40 000 chips: the lower tone sits at 100 kHz, 0.01 cycle per chip, so a real
    #   lower chip is a DC level a cos(phi) whose correlator energy is 4 cos^2(phi) x nominal (the image adds to it); it
    #   errs with 0.5 exp(-2 cos^2(phi) x 1000), over a slowly turning phi 0.5 exp(-1000) I0(1000) = 0.5 / sqrt(2 pi
    #   1000) = 0.0063, on the half of the chips that are lower: 0.0032, 126 of 40 000 (asserted within a factor of 2).
    #   The baseband model has no image: 0.5 exp(-500), no error.
    for rc, e, n, seeds in ((5.25e6, 8.0, 100000, (5, 6)), (9.95e6, 30.0, 40000, (4,))):
        fc, err = np.array([cc.tone_centre(rc)]), [0, 0]
        for s in seeds:
            bits = [np.random.default_rng(s).integers(0, 2, n).astype(np.uint8)]
            for i, corr in enumerate((b.victim_corr, b.passband_corr)):
                Cs, Cn = corr(0, fc, np.ones((1, 2)), bits, [s], rc, "fsk", 32, 2, [s, 0, 2])
                err[i] += int(np.count_nonzero(b.detect(Cs, Cn, 1.0, 0.5 / rc / 10 ** (e / 10), "fsk") != bits[0]))
        if rc == 5.25e6:
            d = b.offset_db(err[1] / n / 2, e) - b.offset_db(err[0] / n / 2, e)
            assert 0 <= d < 0.2, (err, d)
        else:
            assert err[0] == 0 and 63 <= err[1] <= 252, err


def test_far_spaced_tones():
    # 3 tones 32 R_c apart (guard 7): a neighbour's nearest component is 31 R_c from the victim centre, deep in the
    # filter stop band, so each tone must behave like a single tone: the same tone simulated alone (M = 0, same chips,
    # noise and samples per chip) differs by less than one binomial standard error (20 000 chips at p ~= 0.056 at
    # 6 dB: sqrt(20000 x 0.056 x 0.944) = 33 errors), and each tone is within 0.5 dB of theory (6 dB -> 0.0683).
    n, Rc = 20000, 500e3
    g = b.Geometry("test", "flat", "", Rc, 7.0, n_tones=3)
    p = b.plan(g)
    N0 = p.e_chip[0] / 10 ** 0.6
    bits = [np.random.default_rng([9, j]).integers(0, 2, n).astype(np.uint8) for j in range(3)]
    d = b.offsets(g, 3, p.spc)
    far, alone = (b.baseband_errors(p.f, p.amps, bits, d, Rc, N0, 1.0, M, p.spc, 9) for M in (2, 0))
    for e, e0 in zip(far, alone):
        assert abs(int(e) - int(e0)) < np.sqrt(n * 0.056 * 0.944), (far, alone)
        assert abs(b.offset_db(e / n, 6.0)) < 0.5, (e, n)


def test_minphase_rc():
    # b.minphase from the magnitude alone against filters whose phase is known exactly, x = f / 1 MHz:
    # - first-order RC, H = 1 / (1 + j x): phase -atan(x), -pi/4 at 1 MHz; group delay RC / (1 + x^2), RC = 1 / (2 pi
    #   1 MHz) = 159.15 ns at 0 Hz;
    # - two equal RC sections (unloaded ladder), H = 1 / (1 - x^2 + 3 j x): phase -pi/2 at x = 1 (H = 1 / 3j), group
    #   delay 3 RC = 477.46 ns at 0 Hz (phase -3x for small x).
    # The magnitude is kept exactly (the folded cepstrum keeps its even part). The error is the grid's end at
    # F = 262.144 MHz, where both fall at k x 20 dB per decade (k = 1, 2): by Bode's gain-phase integral the missing
    # tail beyond F adds (2k / pi)(f / F) of phase for f << F, and the grid's mirror image beyond F (rising again)
    # about as much, so about (4k / pi)(f / F): 0.0049 k rad at 1 MHz, a pure delay 2k / (pi^2 F) = 0.77 k ns, constant
    # to order (f / F)^2 (< 1 % of it up to 20 MHz). Asserted: 1.3 times these leading-order values, and the delay error
    # constant within 0.05 ns over 0.05 to 20 MHz (a constant delay adds no ISI).
    f = np.arange(b.MP_HALF + 1) * b.MP_DF
    x, F = f / 1e6, b.MP_HALF * b.MP_DF
    band = (f >= 50e3) & (f <= 20e6)
    rc = 1 / (2 * np.pi * 1e6)
    for k, H, ph in ((1, 1 / (1 + 1j * x), -np.pi / 4), (2, 1 / (1 - x ** 2 + 3j * x), -np.pi / 2)):
        assert abs(b.group_delay(H)[0] - (1 if k == 1 else 3) * rc) < 0.1e-9, k
        Hm = b.minphase(20 * np.log10(np.abs(H)))
        assert np.allclose(np.abs(Hm), np.abs(H), rtol=1e-9), k
        assert abs(np.angle(Hm[1000]) - ph) < 1.3 * 4 * k / np.pi * 1e6 / F, (k, np.angle(Hm[1000]))
        err = (b.group_delay(Hm) - b.group_delay(H))[band]
        assert np.abs(err).max() < 1.3 * 2 * k / (np.pi ** 2 * F) and np.ptp(err) < 0.05e-9, (k, err.min(), err.max())


def test_curve_variants():
    # curve_db, cohort 50 Ohm (csv: 50 kHz -77.5781, 100 kHz -73.5576, 7 MHz -45.2582, 10 MHz -43.6213 dB):
    # - extension A below 50 kHz: the first segment's slope 4.0205 dB per log10(2) = 0.30103 decade = 13.356 dB per
    #   decade, one decade down: 5 kHz at -77.5781 - 13.356 = -90.934 dB, held below (0 Hz reads 1 Hz: the same)
    # - smoothed: 8 MHz on the log-f line from 7 to 10 MHz, t = log(8/7) / log(10/7) = 0.13353 / 0.35667 = 0.37438:
    #   -45.2582 + 0.37438 x 1.6369 = -44.6454; 9 MHz, t = log(9/7) / log(10/7) = 0.25131 / 0.35667 = 0.70461:
    #   -45.2582 + 0.70461 x 1.6369 = -44.1048; 12.5 MHz and above held at the 10 MHz -43.6213 dB (1 GHz too)
    # - as measured: PCHIP passes through every measured point (10 MHz -43.6213)
    got = b.curve_db("cohort", "smoothed", "A", [5e3, 0.0, 8e6, 9e6, 12.5e6, 20e6, 1e9])
    want = [-90.934, -90.934, -44.6454, -44.1048, -43.6213, -43.6213, -43.6213]
    assert np.allclose(got, want, atol=2e-3), got
    assert abs(b.curve_db("cohort", "measured", "A", [10e6])[0] + 43.6213) < 1e-9


def test_pure_delay_channel():
    # a channel that only delays (7 samples, 24.3 ns at 9 Mchip/s x 32) must cost nothing. Delay commutes with
    # the band-pass, so with the receiver's timing moved by the same 7 samples every correlator output equals the one
    # without a channel (up to FFT rounding): rx_lag must find the undelayed optimum plus 7, with the same margin.
    rc = 9e6
    fc = cc.tone_centre(rc)
    spc = b.isi_spc(fc, rc)
    ones, hd = np.ones((1, 2)), np.r_[np.zeros(7), 1.0]
    (l0, m0), (l7, m7) = b.rx_lag(np.ones(1), ones, fc, rc, spc), b.rx_lag(hd, ones, fc, rc, spc)
    assert l7 == l0 + 7 and abs(m7 - m0) < 1e-9, (l0, l7, m0, m7)
    bits = [np.random.default_rng(3).integers(0, 2, 5000).astype(np.uint8)]
    (C0, C7), _ = b.passband_corr(0, np.array([fc]), ones, bits, [5], rc, "fsk", spc, 0, 3,
                                  chans=[(np.ones(1), l0), (hd, l0 + 7)])
    assert np.allclose(C0, C7, rtol=0, atol=1e-9 * np.abs(C0).max())


if __name__ == "__main__":
    t0 = time.time()
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            t = time.time(); fn(); print(f"ok   {name} ({time.time() - t:.1f} s)")
    print(f"all tests passed in {time.time() - t0:.1f} s")

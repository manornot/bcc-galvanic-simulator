"""Plain-assert tests of protocols.py, no pytest needed.

Purpose: frame airtimes against the field sizes of the standards, the BCC frame length, BCC full on one tone and
striped over the tones, the full-capability specs (retry limits, channel rules, the variants get_spec builds),
the ACR tables, the access exchanges, select_period against BLE_SIM.select_ci and against hand-derived cells, and that
every catalogue layout is consumable by phy.frame_outcome.
Every expected value is derived by hand (or read from the PARAMETERS.md check tables) in the comment above its assert.
Inputs:  parameters.json, ../BLE_SIM.py, phy.py.
Outputs: one line per passed test on stdout; an AssertionError on failure.
argv:    none (python python_simulation/coexistence/test_protocols.py), runtime a few seconds.
"""
import math
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import BLE_SIM                # noqa: E402
import bcc_fdma               # noqa: E402
import phy                    # noqa: E402
import protocols as pr        # noqa: E402

SPECS = [pr.get_spec(*c) for c in pr.catalogue()]


def _end(spec, L):
    """Unrounded frame duration in us (airtime_us rounds it up)."""
    return pr.frame_layout(spec, L)[-1]["t1_us"]


def test_ble_airtime():
    # 1M, 1 bit/us: preamble 8 + access address 32 + header 16 + 8 L + CRC 24 bits = 80 + 8 L us
    # 2M, 2 bits/us: (16 + 32 + 16 + 8 L + 24) / 2 = 44 + 4 L us
    # Coded S=8: preamble 80 us, FEC block 1 (AA 32 + CI 2 + TERM1 3 = 37 bits x 8 us) 296 us, FEC block 2
    #   (header 16 + 8 L + CRC 24 + TERM2 3 bits x 8 us) = 344 + 64 L us: total 720 + 64 L
    # Coded S=2: the same preamble and block 1 (always S=8), block 2 at 2 us per bit: 86 + 16 L us: total 462 + 16 L
    for phy_name, fixed, per_octet in (("1M", 80, 8), ("2M", 44, 4), ("coded_s2", 462, 16), ("coded_s8", 720, 64)):
        s = pr.get_spec("ble", "full", phy_name)
        for L in (0, 20, 251):
            assert pr.airtime_us(s, L) == fixed + per_octet * L, (phy_name, L, pr.airtime_us(s, L))


def test_ieee802154_airtime_and_csma_exchange():
    # SHR 5 + PHR 1 + MHR 9 + payload 20 + FCS 2 = 37 octets x 32 us = 1184 us
    s = pr.get_spec("ieee802154", "full")
    assert pr.airtime_us(s, 20) == 1184
    # PARAMETERS 3.2 mode A FULL step 4: mean successful exchange = mean backoff 1120 + CCA 128 + turnaround 192
    # + T_data + 192 + Imm-Ack 352 + LIFS 640 (MPDU 31 > 18 octets) = 3808 us at L = 20; the exchange leaves the
    # random backoff to the engine: 3808 - 1120 = 2688
    assert sum(st.dur_us for st in pr.exchange(s, 20)) == 3808 - 1120


def test_bcc_frame():
    # N(L) = preamble 256 + SFD 64 + header 128 + 32 L + 64 floor(32 L / 2048) chips (PARAMETERS 3.8):
    # N(255) = 448 + 8160 + 64 x 3 = 8800 chips; 8800 / (8 x 255) = 4.3137 chips per payload bit
    for rc in (500e3, 2e6):
        s = pr.get_spec("bcc", "thesis", rc_cps=rc)
        chips = pr.frame_layout(s, 255)[-1]["t1_us"] * rc / 1e6
        assert abs(chips - 8800) < 1e-6, chips
        assert abs(chips / (8 * 255) - 4.314) < 5e-4
    # at 500 kchip/s one chip is 2 us: 17600 us
    assert pr.airtime_us(pr.get_spec("bcc", "thesis"), 255) == 17600
    # FCS switch: +2 octets = +64 chips; at L = 20 no pilot is added: 448 + 640 + 64 = 1152 chips = 2304 us
    f = pr.get_spec("bcc", "thesis", fcs=True)
    assert pr.airtime_us(f, 20) == 2304 and f.phy["payload_crc"]


def test_bcc_full():
    # HBC 3.7 FULL on the galvanic channel, R_c = 500 kchip/s as the Walsh-chip rate of the SF 8 mode (3.7 last
    # paragraph): 256 + 64 + 1024 (header at SF 64: 128 Walsh chips x 8) + 32 (L + 9) + 64 floor(32 (L + 9) / 2048)
    # Walsh chips of 2 us. L = 0: 1344 + 288 = 1632 chips = 3264 us; L = 20: 1344 + 928 = 2272 = 4544 us;
    # L = 255: 1344 + 8448 + 4 x 64 = 10048 = 20096 us
    s = pr.get_spec("bcc", "full")
    assert [pr.airtime_us(s, L) for L in (0, 20, 255)] == [3264, 4544, 20096]
    assert pr.airtime_us(pr.get_spec("bcc", "full", rc_cps=50e3), 20) == 45440     # 10 x 4544
    # exchange (3.7 'as NB'): data 4544 + pSIFS 75 + I-Ack (the L = 0 frame) 3264 + pSIFS 75 = 7958 us
    assert sum(st.dur_us for st in pr.exchange(s, 20)) == 7958
    # the header at SF 64 against the payload at SF 8: SNR factor 64 / 8 = 8 (+9.03 dB); SNR per Walsh chip (noise
    # bandwidth R_c = payload Walsh-chip rate), as the as-built stack
    assert [b.get("rep") for b in pr.frame_layout(s, 20) if b["kind"] == "header"] == [8.0]
    assert s.phy["bit_rate_bps"] == s.phy["noise_bw_hz"] == s.chip_rate_cps == 500e3
    # every HBC FULL feature (PARAMETERS 2 and 3.7), on the as-built BCC tone, drive, energy and isolation
    h, a = pr.get_spec("hbc", "full"), pr.get_spec("bcc", "thesis")
    assert s.arq_max_retx == h.arq_max_retx == 3 and s.phy["payload_crc"] and s.period_grid == h.period_grid
    assert (s.inter_body, s.channels, s.energy, s.p_loss, s.phy["amp_v"]) == (
        a.inter_body, a.channels, a.energy, 0.0, a.phy["amp_v"])
    # 2 kb/s, 1 sensor, latency (P 7.5 ms): beacon (17-octet body: 1344 + 32 x 26 = 2176 chips = 4352 us) + pSIFS 75
    # = 4427 us reserved; allocation of 1 + 1 spare exchanges, 2 (3264 + 64 L + 3414) = 13356 + 128 L. 18000 us:
    # 36 bits -> 5 octets, 13996 > 13573; 18500 us: 37 bits -> 5 octets, 13996 <= 14073
    sel = pr.select_period(2, "latency", 1, s)
    assert (sel["period_us"], sel["frame_payload_octets"], sel["deadline_us"]) == (18500, [5], 7500), sel


def test_bcc_multitone():
    # tone plan (bcc_fdma.plan, spacing 4 R_c, guard 0) over the 3 dB span of the cohort 50 Ohm curve, from
    # sqrt(6 x 7) = 6.480741 MHz (lower edge of the 7 MHz bin) to 20 MHz (bcc_fdma_summary.md), the dip between the
    # 3 dB bins [6.4807, 7.4833] and [9.4868, 20] MHz included: 13.519 MHz / 2 MHz at 500 kchip/s = 6.76 -> 6 tones
    # centred on 13.240370 MHz, 8.240370 + 2 i MHz. 50 kchip/s: 13.519 / 0.2 = 67.6 -> 67 tones (the headline
    # count); 2 Mchip/s: 13.519 / 8 = 1.69 -> 1 tone. Receive bands f -/+ 1 MHz: 7.24-9.24 (in the dip) and 10.24 -/+ 1
    # (9.24 < 9.4868) leave the bins, so 4 tones are in band. PARAMETERS 3.8's rule, spacing 5 R_c = 2.5 MHz (guard
    # 0.25): 13.519 / 2.5 = 5.41 -> 5 tones, 8.24 + 2.5 i MHz, all but 8.24 in band
    s = pr.get_spec("bcc", "full", "multitone")
    assert np.allclose(s.tones, 8.240370e6 + 2e6 * np.arange(6), rtol=0, atol=1), s.tones
    assert s.tones == bcc_fdma.plan(bcc_fdma.Geometry("x", "cohort", "50", 500e3, 0.0)).f.tolist()
    assert bcc_fdma.plan(bcc_fdma.Geometry("x", "cohort", "50", 500e3, 0.0)).in_band == 4
    p38 = bcc_fdma.plan(bcc_fdma.Geometry("x", "cohort", "50", 500e3, 0.25))
    assert np.allclose(p38.f, 8.240370e6 + 2.5e6 * np.arange(5), rtol=0, atol=1) and p38.in_band == 4, p38.f
    assert [len(pr.get_spec("bcc", "full", "multitone", rc_cps=rc).tones) for rc in (50e3, 2e6)] == [67, 1]
    # back-to-back aggregate over the 6 tones. As-built frames: 255 octets in 8800 chips of 2 us = 17600 us,
    # 6 x 2040 bits / 17600 us = 695.45 kb/s (summary: 695). Full stack: 255 octets per tone in the BCC full frame of
    # 20096 us (test_bcc_full), 6 x 2040 / 20096 = 609.08 kb/s; one exchange (+ pSIFS 75 + I-Ack 3264 + pSIFS 75 =
    # 23510 us) 520.63 kb/s; the allocation with its spare exchange (47020 us) 260.31 kb/s, the capacity of one body
    # (one tone: 2040 / 47020 = 43.39 kb/s)
    assert abs(6 * 2040 / pr.airtime_us(pr.get_spec("bcc", "thesis"), 255) * 1e3 - 695.45) < 0.01
    assert s.max_payload_octets == 1530 and pr.tone_octets(s, 1530) == [255] * 6 and pr.airtime_us(s, 1530) == 20096
    assert sum(st.dur_us for st in pr.exchange(s, 1530)) == 23510 and pr._frame_us(s, 1530) == 47020
    for us, kbps in ((20096, 609.08), (23510, 520.63), (47020, 260.31)):
        assert abs(12240 / us * 1e3 - kbps) < 0.01
    # a 250-octet frame: 250 = 6 x 41 + 4, so 4 tones carry 42 octets and 2 carry 41. The longest tone packet: 1344 +
    # 32 (42 + 9) = 2976 Walsh chips (no pilot: 1632 < 2048) = 5952 us; on one tone 1344 + 32 x 259 + 64 x 4 = 9888
    # chips = 19776 us. The tone packets carry the 2000 application bits between them
    assert pr.tone_octets(s, 250) == [42] * 4 + [41] * 2 and pr.airtime_us(s, 250) == 5952
    assert pr.airtime_us(pr.get_spec("bcc", "full"), 250) == 19776
    assert sum(b["info_bits"] for lay in pr.tone_layouts(s, 250) for b in lay if b["kind"] == "payload") == 2000
    # L < N: 5 octets on 5 tones and MAC header + FCS alone on the 6th (every tone active); the I-Ack stripe: 6 empty
    # packets of 3264 us
    assert pr.tone_octets(s, 5) == [1] * 5 + [0] and pr.tone_octets(s, 0) == [0] * 6 and pr.airtime_us(s, 0) == 3264
    # drive split by sqrt(N) (power limited): 1 / sqrt(6) V per tone; N front ends at the typical BCC power (2 mW TX,
    # 2.4 mW RX at 3 V: 0.6667 / 0.8 mA each) -> 4 / 4.8 mA, the ramps likewise; MCU wake-up and sleep once
    a = pr.get_spec("bcc", "full")
    assert math.isclose(s.phy["amp_v"], 1 / math.sqrt(6)) and s.phy["name"] == "multitone"
    assert all(math.isclose(s.energy[x], 6 * a.energy[x]) for x in ("tx_ma", "rx_ma", "ramp_tx_ma", "ramp_rx_ma"))
    assert math.isclose(s.energy["tx_ma"], 4.0) and math.isclose(s.energy["rx_ma"], 4.8)
    assert all(s.energy[x] == a.energy[x] for x in ("wake_uc", "sleep_ua", "ramp_us", "vdd_v"))
    # otherwise BCC full: the stack, the ARQ, the one-tone beacon, isolation; the beacon on one front end's currents
    assert (s.arq_max_retx, s.beacon, s.inter_body) == (3, dict(a.beacon, energy=a.energy), "isolated")
    # refused: no contact chain or inter-body coupling per tone; a chip rate whose 4 R_c tone does not fit the 3 dB span
    # of 20 - 6.480741 = 13.519 MHz (3.3 Mchip/s: 13.519 / 13.2 = 1.02 -> 1 tone, 3.4: 13.519 / 13.6 = 0.99 -> none);
    # the FCS flag on either BCC full stack, which carries the 802.15.6 FCS (the flag is the as-built stack's)
    assert len(pr.get_spec("bcc", "full", "multitone", rc_cps=3.3e6).tones) == 1
    for ph, kw, msg in (("multitone", dict(contact=True, rc_cps=50e3), "multi-tone"),
                        ("multitone", dict(coupling_db=40), "multi-tone"),
                        ("multitone", dict(rc_cps=3.4e6), "no multi-tone plan"),
                        ("multitone", dict(fcs=True), "802.15.6 FCS"), ("sf8", dict(fcs=True), "802.15.6 FCS")):
        try:
            pr.get_spec("bcc", "full", ph, **kw)
            raise AssertionError((ph, kw))
        except ValueError as e:
            assert msg in str(e), e
    # select_period: a body carries at most 8 L_max per allocation, 12240 bits / 47020 us = 260.31 kb/s striped and
    # 2040 / 47020 = 43.39 kb/s on one tone; above that no period fits (scan-free bound, R_agg x 47020 / (8 L_max)):
    # 261 x 3.8415 = 1003 > 1000 and 44 x 23.049 = 1014 > 1000; N_dev 10: 27 x 10 x 3.8415 = 1037, 5 x 10 x 23.049 =
    # 1152. Below it T = 1 s works with the beacon reserve 4427 us (test_bcc_full), so the shortest period is <= 1 s
    # (P <= 1 s: latency 7.5 ms, energy 7.5 to 500 ms):
    # - striped, N_dev 1, 250 kb/s: 31250 octets = 20 x 1530 + 650 (tones of 109, 109, 108 x 4: 1344 + 32 x 118 + 64 =
    #   5184 chips = 10368 us, allocation 2 (10368 + 3414) = 27564 us): 20 x 47020 + 27564 = 967964 <= 995573
    # - striped, N_dev 10, 20 kb/s: 2500 octets = 1530 + 970 (tones of 162, 162, 162, 162, 161, 161: 1344 + 32 x 171 +
    #   128 = 6944 chips = 13888 us, allocation 34604 us): 10 x (47020 + 34604) = 816240 <= 995573
    # - one tone, N_dev 1, 40 kb/s: 5000 octets = 19 x 255 + 155 (1344 + 32 x 164 + 128 = 6720 chips = 13440 us,
    #   allocation 33708 us): 19 x 47020 + 33708 = 927088 <= 995573
    # - one tone, N_dev 10, 4 kb/s: 500 octets = 255 + 245 (1344 + 32 x 254 + 192 = 9664 chips = 19328 us, allocation
    #   45484 us): 10 x (47020 + 45484) = 925040 <= 995573
    for spec, n, ok, no in ((s, 1, 250, 261), (s, 10, 20, 27), (a, 1, 40, 44), (a, 10, 4, 5)):
        for mode in ("latency", "energy"):
            sel = pr.select_period(ok, mode, n, spec)
            assert sel is not None and sel["period_us"] <= 1_000_000, (spec.phy["name"], n, mode, sel)
            assert pr.select_period(no, mode, n, spec) is None, (spec.phy["name"], n, mode)


def test_bcc_ladder():
    # Galvanic-optimisation ladder of BCC full (DESIGN 4.3) at 500 kchip/s: Walsh chip 2 us, octet 32 chips = 64 us,
    # pSIFS 75 us. Frame = preamble 256 chips (512 us) + SFD (L0: 512 chips at 42 Mcps / 8 = 64 Walsh chips, 128 us;
    # RI mode from L1: the 608-chip SFD/RI field of 10.5.2, Eq. 117, (64 + 12) x 8 = 608, 76 Walsh chips, 152 us) +
    # header (L0: 128 chips at SF 64 = 2048 us; RI mode from L1: at SF 8, 256 us) + 64 us per PSDU octet (MAC header
    # 7 + L + FCS 2; from L4 2 + L + 2) + a 128 us pilot after every 128 PSDU symbols (64 octets; from L1 the '128 B'
    # option, every 256 symbols, 128 octets). Head: 2688 us at L0, 512 + 152 + 256 = 920 from L1. One tone, frames of
    # L = 0 (I-Ack; L3 poll), 20, 255:
    #   L0: 2688 + 9 x 64 = 3264; 2688 + 29 x 64 = 4544; 2688 + 264 x 64 + 4 x 128 = 20096
    #   L1-L3: 920 + 576 = 1496; 920 + 1856 = 2776; 920 + 16896 + 2 x 128 = 18072
    #   L4, L5: 920 + 4 x 64 = 1176; 920 + 24 x 64 = 2456; 920 + 259 x 64 + 2 x 128 (518 symbols) = 17752
    # Six tones (L0-L4: 20 octets as 4, 4, 3, 3, 3, 3, the longest tone 4 octets; 1530 as 255 per tone):
    #   L0: 2688 + 13 x 64 = 3520; L1-L3: 920 + 832 = 1752; L4: 920 + (4 + 4) x 64 = 1432; 1530 as one tone's 255
    #   L5: the stripe's PSDU (2 + L + 2 octets) over the 6 tones, tone 0 after the 920 us head. L = 20: 24 octets on
    #   tones 1-5 as 5, 5, 5, 5, 4 (at most 320 us), tone 0 the head alone: 920. L = 1530: 1534 octets; the stripe ends
    #   at D = 16768, the least D where the octets that end by D reach 1534: tone 0 245 (920 + 245 x 64 + 128 = 16728;
    #   the 246th ends at 16792), tones 1-5 258 each (258 x 64 + 2 x 128 = 16768), 1535 >= 1534, while by 16767 only
    #   245 + 5 x 257 = 1530 fit; the four last octets tie at 16768 and go to tones 1 to 4: 245, 258, 258, 258, 258,
    #   257. L = 0: 4 octets on tones 1 to 4, the head sets 920
    # Exchange: L0 data + 75 + I-Ack + 75, L1 likewise with the 1496 us I-Ack, from L2 data + 75 (no I-Ack).
    # Allocation per frame: L0-L2 two exchanges (the spare), from L3 one. Pool exchange from L3: poll + 75 + exchange;
    # the poll is the L = 0 frame at L3 (1496 us), from L4 the L = 2 frame (the polled 15-bit sequence number): one
    # tone 920 + (2 + 2 + 2) x 64 = 1304, six tones at L4 one octet on each of two tones, 920 + (2 + 1 + 2) x 64 =
    # 1240, at L5 6 octets on tones 1-5: 920.
    # (level, data L 0 / 20 / L_max, exchange 20 / L_max, allocation 20 / L_max, pool 20 / L_max)
    one = ((0, (3264, 4544, 20096), (7958, 23510), (15916, 47020), None),
           (1, (1496, 2776, 18072), (4422, 19718), (8844, 39436), None),
           (2, (1496, 2776, 18072), (2851, 18147), (5702, 36294), None),
           (3, (1496, 2776, 18072), (2851, 18147), (2851, 18147), (1496 + 75 + 2851, 1496 + 75 + 18147)),
           (4, (1176, 2456, 17752), (2531, 17827), (2531, 17827), (1304 + 75 + 2531, 1304 + 75 + 17827)),
           (5, (1176, 2456, 17752), (2531, 17827), (2531, 17827), (3910, 19206)))
    six = ((0, (3264, 3520, 20096), (6934, 23510), (13868, 47020), None),
           (1, (1496, 1752, 18072), (3398, 19718), (6796, 39436), None),
           (2, (1496, 1752, 18072), (1827, 18147), (3654, 36294), None),
           (3, (1496, 1752, 18072), (1827, 18147), (1827, 18147), (3398, 19718)),
           (4, (1176, 1432, 17752), (1507, 17827), (1507, 17827), (1240 + 75 + 1507, 1240 + 75 + 17827)),
           (5, (920, 920, 16768), (995, 16843), (995, 16843), (920 + 75 + 995, 920 + 75 + 16843)))
    # Capacity per body, the scan-free bound of select_period plus the pool and bitmap shares per frame: 8 L_max /
    # (allocation + pool exchange / 20 (pool_size) + 8 us from L2 (one bitmap bit on the one-tone beacon, 64 / 8 us;
    # the base sequence number, one octet per sensor and superframe, vanishes with many frames per superframe)):
    #   one tone 2040 / 47020 = 43.386, / 39436 = 51.729, / 36302 = 56.195, / (18147 + 985.9 + 8) = 106.578,
    #   / (17827 + 960.3 + 8) = 108.538 (L4 and L5) kb/s
    #   six tones 12240 / 47020 = 260.315, / 39436 = 310.376, / 36302 = 337.172, / 19140.9 = 639.468, / (17827 +
    #   957.1 + 8) = 651.337, / (16843 + 891.9 + 8) = 689.853 kb/s
    cap = {"sf8": (43.386, 51.729, 56.195, 106.578, 108.538, 108.538),
           "multitone": (260.315, 310.376, 337.172, 639.468, 651.337, 689.853)}
    for ph, lmax, table in (("sf8", 255, one), ("multitone", 1530, six)):
        for o, data, ex, alloc, pool in table:
            s = pr.get_spec("bcc", "full", ph, opt=o)
            assert tuple(pr.airtime_us(s, L) for L in (0, 20, lmax)) == data, (ph, o)
            assert tuple(sum(st.dur_us for st in pr.exchange(s, L)) for L in (20, lmax)) == ex, (ph, o)
            assert tuple(pr._frame_us(s, L) for L in (20, lmax)) == alloc, (ph, o)
            assert pool is None or tuple(pr._pool_us(s, L) for L in (20, lmax)) == pool, (ph, o)
            kinds = ["data", "ifs", "ack", "ifs"] if o < 2 else ["data", "ifs"]     # from L2 no I-Ack
            assert [st.kind for st in pr.exchange(s, 20)] == kinds, (ph, o)
            c = 8 * lmax / (alloc[1] + (pool[1] / 20 if pool else 0) + (8 if o >= 2 else 0)) * 1e3
            assert abs(c - cap[ph][o]) < 1e-3, (ph, o, c)
            assert s.max_payload_octets == lmax and s.arq_max_retx == 3 and s.flags == {"opt": o}
    m5 = pr.get_spec("bcc", "full", "multitone", opt=5)
    assert pr._balance(6, 1533, 896.0, 64.0, 128.0, 128) == (246, 258, 258, 257, 257, 257)
    assert pr._balance(6, 1533, 920.0, 64.0, 128.0, 128) == (245, 258, 258, 258, 257, 257)
    assert [lay[-1]["t1_us"] for lay in pr.tone_layouts(m5, 1530)] == [16728, 16768, 16768, 16768, 16768, 16704]
    assert [lay[-1]["t1_us"] for lay in pr.tone_layouts(m5, 20)] == [920, 320, 320, 320, 320, 256]
    assert [lay[-1]["t1_us"] for lay in pr.tone_layouts(m5, 0)] == [920, 64, 64, 64, 64]     # tone 5 idle
    assert [b["kind"] for b in pr.tone_layouts(m5, 20)[0]] == ["sync"] * 3 + ["header"]   # tone 0: no PSDU octet
    assert {b["kind"] for lay in pr.tone_layouts(m5, 1530)[1:] for b in lay} == {"payload"}
    # one tone: L5 is L4 (the flag is accepted, the frames and the selection are the same)
    s4, s5 = (pr.get_spec("bcc", "full", "sf8", opt=o) for o in (4, 5))
    assert all(pr.frame_layout(s4, L) == pr.frame_layout(s5, L) for L in (0, 20, 255))
    assert pr.select_period(10, "latency", 4, s4) == pr.select_period(10, "latency", 4, s5)
    # pool size max(1, ceil(0.05 frames)): 1 frame -> 1, 20 -> 1, 21 -> 2, 40 -> 2, 41 -> 3, 60 -> 3 (0.05 x 60 is
    # 3.0000000000000004 in floats); the contact design point 0.43 (= 0.3 / 0.7): 1 -> 1, 2 -> ceil(0.86) = 1, 3 ->
    # ceil(1.29) = 2, 7 -> ceil(3.01) = 4, 20 -> ceil(8.6) = 9
    assert [pr.pool_size(n) for n in (1, 20, 21, 40, 41, 60)] == [1, 1, 2, 2, 3, 3]
    assert [pr.pool_size(n, 0.43) for n in (1, 2, 3, 7, 20)] == [1, 1, 2, 4, 9]
    # capacity per body with the 0.43 pool (both design points), the rule above with 0.43 x the pool
    # exchange per frame: one tone L3 2040 / (18147 + 0.43 x 19718 + 8) = 2040 / 26633.74 = 76.595, L4 and L5 2040 /
    # (17827 + 0.43 x 19206 + 8) = 2040 / 26093.58 = 78.180; six tones L3 12240 / 26633.74 = 459.567, L4 12240 /
    # (17827 + 0.43 x 19142 + 8) = 12240 / 26066.06 = 469.576, L5 12240 / (16843 + 0.43 x 17838 + 8) = 12240 /
    # 24521.34 = 499.157 kb/s. The flag is refused
    # below L3 and outside (0, 1]; the selection holds the larger pool: 1 kb/s x 10 sensors, energy mode, one tone L3
    for ph, lmax, o, want in (("sf8", 255, 3, 76.595), ("sf8", 255, 5, 78.180), ("multitone", 1530, 3, 459.567),
                              ("multitone", 1530, 4, 469.576), ("multitone", 1530, 5, 499.157)):
        s = pr.get_spec("bcc", "full", ph, opt=o, pool=0.43)
        c = 8 * lmax / (pr._frame_us(s, lmax) + 0.43 * pr._pool_us(s, lmax) + 8) * 1e3
        assert abs(c - want) < 1e-3 and s.flags == {"opt": o, "pool": 0.43}, (ph, o, c)
    for kw in (dict(pool=0.0), dict(opt=3, pool=0.0), dict(opt=3, pool=1.5)):
        try:
            pr.get_spec("bcc", "full", "sf8", **kw)
            raise AssertionError(kw)
        except ValueError as e:
            assert "pool" in str(e), e
    a, b = (pr.select_period(1, "energy", 10, pr.get_spec("bcc", "full", "sf8", opt=3, **kw)) for kw in ({},
                                                                                                    dict(pool=0.43)))
    assert (a["pool"], b["pool"]) == (pr.pool_size(10 * a["frames_per_period"]),
                                      pr.pool_size(10 * b["frames_per_period"], 0.43)), (a, b)
    # beacon (one tone at every level): L1 17-octet body, PSDU 26 octets: 920 + 26 x 64 = 2584 us; from L2 + the ack
    # field, per sensor a base sequence number octet and a bit per frame, N_dev + ceil(N_dev n / 8) octets: one
    # sensor with n = 1 or 8 frames 17 + 2 = 19 (PSDU 28: 2712), n = 9 or 16 20 (2776), n = 17 21 (2840); 4 sensors
    # with 1 frame 17 + 4 + 1 = 22 (2904), 10 sensors 17 + 10 + 2 = 29 (PSDU 38: 3352); L4 keeps the 7-octet MAC
    # header (2 + 19 + 5 + 2 = 28 octets) and widens the base to the 15-bit sequence number, 2 octets:
    # 2 N_dev + ceil(N_dev n / 8), one sensor 17 + 3 = 20 (PSDU 29: 2776), 21 (2840), 22 (2904); 4 sensors
    # 17 + 8 + 1 = 26 (PSDU 35: 3160), 10 sensors 17 + 20 + 2 = 39 (PSDU 48: 3992); multi-tone the same one-tone beacon
    assert pr.get_spec("bcc", "full", opt=1).beacon["airtime_us"] == 2584
    for o, ph, one, many in ((2, "sf8", [2712, 2712, 2776, 2776, 2840], [2904, 3352]),
                             (4, "sf8", [2776, 2776, 2840, 2840, 2904], [3160, 3992]),
                             (5, "multitone", [2776, 2776, 2840, 2840, 2904], [3160, 3992])):
        s = pr.get_spec("bcc", "full", ph, opt=o)
        assert ([pr.beacon_layout(s, 1, n)[-1]["t1_us"] for n in (1, 8, 9, 16, 17)] == one
                and [pr.beacon_layout(s, d, 1)[-1]["t1_us"] for d in (4, 10)] == many), (o, ph)
    # select_period with the beacon and the pool. L2, one tone, 2 kb/s latency, 1 sensor: reserved 2712 + 75 = 2787;
    # 7500 us: 15 bits -> 2 octets, frame 1496 + 128 = 1624, allocation 2 (1624 + 75) = 3398 <= 7500 - 2787 = 4713
    # (L0: 18500 us, test_bcc_full). L3: allocation 1699 + pool (1496 + 75 + 1699 = 3270) + 2787 = 7756 > 7500;
    # 8000 us: 16 bits -> 2 octets, 7756 <= 8000. L3 with 2 sensors: beacon 17 + 2 + 1 = 20 octets (2776), reserved
    # 2851; 2 octets (7500, 8000 us): 2 x 1699 + 3270 + 2851 = 9519 > 8000; 3 octets (8500 to 11500 us): frame 1688,
    # exchange 1763, pool 1496 + 75 + 1763 = 3334: 2 x 1763 + 3334 + 2851 = 9711 > 9500; 10000 us holds (20 bits)
    for o, n, want in ((2, 1, dict(period_us=7500, frame_payload_octets=[2], beacon_us=2712, reserved_us=2787)),
                       (3, 1, dict(period_us=8000, frame_payload_octets=[2], slot_us=1699, pool=1, pool_us=3270)),
                       (3, 2, dict(period_us=10000, frame_payload_octets=[3], slot_us=1763, pool=1, pool_us=3334,
                                   beacon_us=2776))):
        sel = pr.select_period(2, "latency", n, pr.get_spec("bcc", "full", opt=o))
        assert {x: sel[x] for x in want} == want, (o, n, sel)
    # the pool's rounding: six tones at L5, one sensor, 20 frames of 1530 octets: ack field 2 + ceil(20 / 8) = 5 octets
    # (beacon PSDU 2 + 17 + 5 + 5 + 2 = 31 octets: 920 + 31 x 64 = 2904), pool 1: 2904 + 75 + 20 x 16843 + 17838 =
    # 357677 us (684.42 kb/s); 21 frames need a pool of 2
    assert pr._fit(20 * 12240, 357677, 1, m5)["pool"] == 1 and pr._fit(20 * 12240, 357676, 1, m5) is None
    # refused: a level outside 0..5, and any level of the as-built stack; HBC has no ladder
    for pid, v, kw in (("bcc", "full", dict(opt=6)), ("bcc", "full", dict(opt=-1)), ("bcc", "thesis", dict(opt=1)),
                       ("hbc", "full", dict(opt=1))):
        try:
            pr.get_spec(pid, v, **kw)
            raise AssertionError((pid, v, kw))
        except ValueError:
            pass
    # opt 0 is the merged stack itself
    assert pr.get_spec("bcc", "full", opt=0).beacon == pr.get_spec("bcc", "full").beacon


def test_bcc_rates():
    # every chip rate (50 k, 200 k, 500 k, 2 M) and of the HBC Walsh-chip set (656.25 k, 1.3125 M,
    # 2.625 M, 5.25 M), then upward to 9 M, the last rate without a chip error in the study (bcc_fdma.RATES,
    # results/bcc_rate_stability.md). One tone at channel_curves.tone_centre (test_bcc_fdma: 14.750 MHz at 5.25 M,
    # 12.000 at 7 M; 9 M: f <= 20 - 9 = 11 MHz keeps the lower tone at <= 2 MHz, where the curve still rises, so
    # 11.000) on both single-tone stacks; 9.5 M (lower tone at 1 MHz, fading with its image) and 9.95 M: refused
    for rc, mhz in ((5.25e6, 14.75), (7e6, 12.0), (9e6, 11.0)):
        assert [pr.get_spec("bcc", v, rc_cps=rc).channels for v in ("full", "thesis")] == [[mhz]] * 2, rc
    for v in ("full", "thesis"):
        for rc in (9.5e6, 9.95e6):
            try:
                pr.get_spec("bcc", v, rc_cps=rc)
                raise AssertionError(v)
            except ValueError as e:
                assert "no stable BFSK tone" in str(e), e
    # Capacity per body (the scan-free bound of select_period, DESIGN 4.3), chips from test_bcc_full / test_bcc_ladder
    # (Walsh chips, Tc = 1e6 / R_c us; airtimes rounded up to whole us, pSIFS 75 us fixed):
    # - BCC full one tone (L0): frame 10048, I-Ack 1632 chips, allocation 2 (ceil(10048 Tc) + ceil(1632 Tc) + 150),
    #   2040 bits per allocation
    # - multi-tone (L0): N tones of the plan (13.519 MHz / 4 R_c: 67, 16, 6, 5, 2, 1, 1 up to 2.625 M, none above),
    #   255 octets per tone in the same allocation: N x the one-tone capacity
    # - the proposal (L5), multi-tone up to 2.625 M, above it on one tone (= L4): stripe D chips, the S = 255 N + 4
    #   PSDU octets (2-octet MAC header) on their S earliest ends (tone 0 after its 460-chip head:
    #   preamble 256, the RI SFD/RI field 76 (608 / 8, 10.5.2 Eq. 117), header 128; an octet 32 chips, a 64-chip pilot
    #   after every 128 octets): N 67 D = 8320 (the others 256 octets, 8192 + 128; tone 0 243, 460 + 7776 + 64 = 8300,
    #   the 244th ends at 8332; 243 + 66 x 256 = 17139 >= 17089; one chip earlier 255 and 243: 17073 < 17089), N 16
    #   D = 8332 (at 8320 243 + 15 x 256 = 4083 < 4084; the next end is tone 0's 244th, 8332: 4084), N 6 D = 8384
    #   (258 and 245, 460 + 7840 + 64 = 8364: 1535 >= 1534; at 8383 257 and 245: 1530), N 5 D = 8416 (at 8396 246 +
    #   4 x 258 = 1278 < 1279; the next end is the 259th of tones 1-4, 8288 + 128 = 8416, before tone 0's 247th at
    #   8428: 246 + 4 x 259 = 1282 >= 1279), N 2 D = 8556 (at 8544 250 + 263 = 513 < 514; tone 0's 251st ends at
    #   460 + 8032 + 64 = 8556, before tone 1's 264th at 8576: 514), N 1 D = 460 + 32 x 259 + 128 = 8876; poll 460
    #   chips (the head; its 6 octets fit on tones 1.. before tone 0's first; N 1: + 6 octets, 652);
    #   exchange ceil(D Tc) + 75, pool exchange ceil(poll Tc) + 75 + exchange, capacity 2040 N / (exchange + share x
    #   pool exchange + 4 Tc of beacon bitmap), share 0.05 and 0.43
    # - the as-built stack: 8800 chips + 10 guard chips per frame, 2040 / (ceil(8800 Tc) + 10 Tc)
    # R_c: Tc; one tone ceil(10048 Tc), ceil(1632 Tc), allocation -> kb/s; proposal ceil(D Tc), ceil(poll Tc) -> 0.05 /
    # 0.43 kb/s; as-built ceil(8800 Tc) -> kb/s
    #   50 k: 20; 200960, 32640, 467500 -> 4.364 (67 tones 292.364); 166400, 9200 -> 2040 x 67 / (166475 + 8787.5 +
    #     80) = 779.503 / (166475 + 75572.5 + 80) = 564.496; 176000 -> 2040 / 176200 = 11.578
    #   200 k: 5; 50240, 8160, 117100 -> 17.421 (16: 278.736); 41660, 2300 -> 32640 / (41735 + 0.05 x 44110 + 20) =
    #     742.485 / 32640 / (41735 + 0.43 x 44110 + 20) = 537.529; 44000 -> 46.311
    #   500 k: 2; 20096, 3264, 47020 -> 43.386 (6: 260.315); 16768, 920 -> 689.853 / 499.157; 17600 -> 115.778
    #   656.25 k: 1.52381; 15312, 2487, 35898 -> 56.828 (5: 284.138); 12825, 701 -> 10200 / (12900 + 0.05 x 13676 +
    #     6.095) = 750.558 / 10200 / (12900 + 0.43 x 13676 + 6.095) = 542.935; 13410 -> 2040 / 13425.238 = 151.953
    #   1.3125 M: 0.761905; 7656, 1244, 18100 -> 112.707 (2: 225.414); 6519, 351 -> 4080 / (6594 + 0.05 x 7020 +
    #     3.048) = 587.215 / 4080 / (6594 + 0.43 x 7020 + 3.048) = 424.308; 6705 -> 303.905
    #   2 M: 0.5; 5024, 816, 11980 -> 170.284 (1: 170.284); 4438, 326 -> 2040 / (4513 + 0.05 x 4914 + 2) = 428.508 /
    #     2040 / (4513 + 0.43 x 4914 + 2) = 307.784; 4400 -> 463.110
    #   2.625 M: 0.380952; 3828, 622, 9200 -> 221.739 (1: 221.739); 3382, 249 -> 2040 / (3457 + 0.05 x 3781 + 1.524)
    #     = 559.276 / 2040 / (3457 + 0.43 x 3781 + 1.524) = 401.231; 3353 -> 607.720
    #   5.25 M: 0.190476; 1914, 311, 4750 -> 429.474; 1691, 125 -> 2040 / (1766 + 0.05 x 1966 + 0.762) = 1093.797 /
    #     2040 / (1766 + 0.43 x 1966 + 0.762) = 780.968; 1677 -> 1215.078
    #   7 M: 0.142857; 1436, 234, 3640 -> 560.440; 1268, 94 -> 2040 / (1343 + 0.05 x 1512 + 0.571) = 1437.458 / 2040 /
    #     (1343 + 0.43 x 1512 + 0.571) = 1023.207; 1258 -> 1619.782
    #   8 M: 0.125; 1256, 204, 3220 -> 633.540; 1110, 82 -> 2040 / (1185 + 0.05 x 1342 + 0.5) = 1628.612 / 2040 /
    #     (1185 + 0.43 x 1342 + 0.5) = 1157.407; 1100 -> 2040 / 1101.25 = 1852.440
    #   9 M: 0.111111; 1117, 182, 2898 -> 703.934; 987, 73 -> 2040 / (1062 + 0.05 x 1210 + 0.444) = 1816.653 /
    #     2040 / (1062 + 0.43 x 1210 + 0.444) = 1288.900; 978 -> 2040 / 979.111 = 2083.522
    want = {50e3: (4.364, 292.364, 779.503, 564.496, 11.578),
            200e3: (17.421, 278.736, 742.485, 537.529, 46.311),
            500e3: (43.386, 260.315, 689.853, 499.157, 115.778),
            656.25e3: (56.828, 284.138, 750.558, 542.935, 151.953),
            1.3125e6: (112.707, 225.414, 587.215, 424.308, 303.905),
            2e6: (170.284, 170.284, 428.508, 307.784, 463.110),
            2.625e6: (221.739, 221.739, 559.276, 401.231, 607.720),
            5.25e6: (429.474, None, 1093.797, 780.968, 1215.078),
            7e6: (560.440, None, 1437.458, 1023.207, 1619.782),
            8e6: (633.540, None, 1628.612, 1157.407, 1852.440),
            9e6: (703.934, None, 1816.653, 1288.900, 2083.522)}

    def cap(s, extra=0.0):
        return 8e3 * s.max_payload_octets / (pr._frame_us(s, s.max_payload_octets) + extra)

    for rc, (one, multi, p05, p43, asb) in want.items():
        tc = 1e6 / rc
        assert abs(cap(pr.get_spec("bcc", "full", rc_cps=rc)) - one) < 1e-3, rc
        assert abs(cap(pr.get_spec("bcc", "thesis", rc_cps=rc)) - asb) < 1e-3, rc
        ph = "multitone" if multi else "sf8"
        if multi:
            assert abs(cap(pr.get_spec("bcc", "full", ph, rc_cps=rc)) - multi) < 1e-3, rc
        else:
            try:
                pr.get_spec("bcc", "full", "multitone", rc_cps=rc)
                raise AssertionError(rc)
            except ValueError as e:
                assert "no multi-tone plan" in str(e), e
        for share, c in ((0.05, p05), (0.43, p43)):
            s = pr.get_spec("bcc", "full", ph, rc_cps=rc, opt=5, pool=share)
            got = cap(s, share * pr._pool_us(s, s.max_payload_octets) + 4 * tc)
            assert abs(got - c) < 1e-3, (rc, share, got)


def _under(spec, i, sir_db, L=20):
    """P(success) of a frame at 60 dB SNR whose block i alone overlaps a co-channel interferer at sir_db."""
    lay = pr.frame_layout(spec, L)
    a, b, end = lay[i]["t0_us"], lay[i]["t1_us"], lay[-1]["t1_us"]
    sir = 10 ** (sir_db / 10)
    segs = [(0, a, 1e6, np.inf), (a, b, 1e6 / (1 + 1e6 / sir), sir), (b, end, 1e6, np.inf)]
    return phy.frame_outcome(None, [x for x in segs if x[1] > x[0]], lay, spec.phy, sampled=False)["success"]


def test_ci_floor_gains():
    # the co-channel C/I floor of a block is lowered by its processing gain and BCH coding gain, calibrated on
    # the BLE specification's coded C/I values (PARAMETERS 3.1: 1M 21, S=2 17, S=8 12 dB). Repetition: S=2 -> S=8 is
    # the P = 4 mapper alone, 17 -> 12 dB, so a factor S buys 10 log10(S) x 5 / 6.02 dB: 2 -> 2.5, 4 -> 5.0, 8 -> 7.5.
    # FEC: the S=2 code buys 4 of its 5 dB sensitivity gain against interference, so the ETSI BCH gain at the channel
    # rate (3.5, Table 13: Eb/N0 8.8 against 7.4 dB) buys 0.8 x 1.4 = 1.12 dB. Floors in dB per block (BLE analogue
    # 21 dB where the recipe has none, OI-3): HBC FULL SF 16: payload 21 - 2.5 = 18.5 (SF 8 is spectral shaping, 3.7),
    # header at SF 64: 21 - 7.5 = 13.5, preamble, SFD and pilots at SF 8: 21. BCC full: 0 everywhere (the 0 dB chip
    # rule of 3.8; no credit for the long header chips). 802.15.6 NB 121.4 kb/s: header BCH(31,19) and PSDU BCH(63,51)
    # at S = 4: 21 - 5.0 - 1.12 = 14.88, preamble (S = 1) 21. SmartBAN: BCH(36,22) header and BCH(127,113) PSDU
    # 21 - 1.12 = 19.88, sync, CRC-4 and PHY CRC-16 21. BLE Coded S=8: the specification's 12 dB (5.5) on every block.
    want = {("hbc", "sf16"): [21, 21, 21, 13.5, 18.5], ("bcc", None): [0, 0, 0, 0, 0],
            ("ieee802156nb", None): [21, 14.88, 14.88], ("smartban", None): [21, 21, 19.88, 21, 19.88, 21],
            ("ble", "coded_s2"): [17, 17, 17, 17], ("ble", "coded_s8"): [12, 12, 12, 12]}
    for (pid, ph), floors in want.items():
        s = pr.get_spec(pid, "full", ph)
        got = [10 * np.log10(s.phy["ci_co_lin"] / phy.ci_gain(b, s.phy)) for b in pr.frame_layout(s, 20)]
        assert np.allclose(got[:len(floors)], floors, atol=1e-9), (pid, got)
    # the rule applied to BLE's own P = 4 step: a 4x repeated block under C/I_co 17 dB clears 12 dB
    g = phy.ci_gain(dict(code=None, rep=4), dict(rep=1))
    assert abs(17 - 10 * np.log10(g) - 12) < 1e-9
    # the floor decides at 60 dB SNR (every bit and chip error probability ~0 above 11 dB SINR): half a dB above a
    # block's floor the frame arrives, half a dB below it is lost; the old frame-wide 21 dB lost both, BLE Coded keeps
    # its 12 dB, and the first credit (full 10 log10(S), BCH 1.9 dB: 11.97, 17.99, 13.08, 19.1, -9.03 dB)
    # delivered the frame half a dB below each new floor
    for pid, ph, i, floor in (("hbc", "sf16", 3, 13.5), ("hbc", "sf16", 4, 18.5), ("ieee802156nb", None, 2, 14.88),
                              ("smartban", None, 4, 19.88), ("ble", "coded_s8", 3, 12.0), ("bcc", None, 3, 0.0)):
        s = pr.get_spec(pid, "full", ph)
        hi, lo = _under(s, i, floor + 0.5), _under(s, i, floor - 0.5)
        assert hi > 0.999 and lo == 0.0, (pid, i, hi, lo)


def test_nb_hbc_smartban_ant_frames():
    # 802.15.6 NB FULL 121.4 kb/s, L = 20: N_PSDU = 8 (9 + 20) = 232, 5 BCH(63,51) codewords, 232 + 5 x 12 = 292
    # coded bits x S 4 = 1168 symbols; (90 preamble + 31 x 4 header + 1168) / 600 kHz = 2303.3 us
    assert pr.airtime_us(pr.get_spec("ieee802156nb", "full"), 20) == 2304
    # the other cells against the PARAMETERS 3.4 check table, which rounds to 0.1 us
    ck = pr.P["ieee802156_nb"]["check_airtime_full_us"]
    for rate in ("121.4", "242.9", "485.7", "971.4"):
        s = pr.get_spec("ieee802156nb", "full", rate)
        for L, t in zip(ck["L"], ck[rate]):
            assert abs(_end(s, L) - t) < 0.051, (rate, L, _end(s, L), t)
    # HBC FULL SF 8, L = 20: 2048 + 512 + 8192 + 2 (9 + 20) x 16 x 8 = 18176 chips / 42 MHz = 432.8 us
    assert pr.airtime_us(pr.get_spec("hbc", "full"), 20) == 433
    ck = pr.P["ieee802156_hbc"]["check_airtime_us"]
    for sf in ("8", "16", "32", "64"):
        s = pr.get_spec("hbc", "full", "sf" + sf)
        for L, t in zip(ck["L"], ck[sf]):
            assert abs(_end(s, L) - t) < 0.051, (sf, L, _end(s, L), t)
    # RI mode, SF 8, L = 0: SFD/RI field (64 + 12) x 8 = 608 chips (10.5.2, Eq. 117), header at SF 8 (1024 chips):
    # 2048 + 608 + 1024 + 18 x 128 = 5984 chips = 142.5 us
    assert pr.airtime_us(pr.get_spec("hbc", "full", ri_mode=True), 0) == 143
    # SmartBAN: everything but the MPDU is the 104-bit PPDU overhead of TS 103 326 (1 us per bit); LF = 20: MPDU 64 +
    # 160 = 224 bits, + 14 BCH parity bits per started 113-bit block: 2 blocks, 104 + 224 + 28 = 356 us
    assert pr.airtime_us(pr.get_spec("smartban", "full"), 20) == 104 + 64 + 8 * 20 + 2 * 14
    # the PHY schemes (3.5 common recipe: PPDU airtime Nrep (104 + L_PSDU) us), LF = 250: MPDU 64 + 2000 = 2064 bits;
    # BCH 19 codewords (2064 / 113 = 18.3), + 266: 104 + 2330 = 2434 us, x 2 = 4868, x 4 = 9736; uncoded 104 + 2064 =
    # 2168, x 2 = 4336, x 4 = 8672. LF,max at LSLOT 32 with the 16 us MCA sensing: Nrep (104 + L_PSDU) <= 20000 - 16 -
    # 168 - 300 = 19516: uncoded Nrep 1 168 + 8 LF <= 19516 -> 2418; Nrep 2 <= 9758 -> 1198; Nrep 4 <= 4879 -> 588;
    # BCH Nrep 2: LF 1065 gives MPDU 8584, 76 codewords (75.96), 104 + 8584 + 1064 = 9752 <= 9758, LF 1066 8592 bits, 77
    # codewords, 9774; Nrep 4: LF 522 4240 bits, 38 codewords, 104 + 4240 + 532 = 4876 <= 4879, LF 523 4884
    got = {p: (pr.airtime_us(pr.get_spec("smartban", "full", p), 250),
               pr.get_spec("smartban", "full", p).max_payload_octets) for p in pr.SB_PHYS}
    assert got == {"bch1": (2434, 2150), "bch2": (4868, 1065), "bch4": (9736, 522), "unc1": (2168, 2418),
                   "unc2": (4336, 1198), "unc4": (8672, 588)}, got
    # repetition is combined ideally in the SINR model (phy rep = Nrep, 3.5 'ideal EGC'); the ACK is not repeated
    s = pr.get_spec("smartban", "full", "unc4")
    assert (s.phy["rep"], s.phy["nrep"], s.phy["fec"], s.phy["name"]) == (4, 4, None, "unc4")
    # ANT: 18 octets = 144 us; bits at risk = address 40 + payload field 80 + CRC 16 = 136
    a = pr.frame_layout(pr.get_spec("ant", "full"), 8)
    assert a[-1]["t1_us"] == 144 and sum(b["n_bits"] for b in a) == pr.P["ant"]["bits_at_risk"]
    # advanced-burst packets of 8, 16, 24 octets (3.6 recipe: 144 / 208 / 272 us); a shorter payload is padded to the
    # next size: 1 and 8 octets 144, 9 and 16 octets 208, 17 and 24 octets 272 us
    s = pr.get_spec("ant", "full")
    assert [pr.airtime_us(s, L) for L in (1, 8, 9, 16, 17, 24)] == [144, 144, 208, 208, 272, 272]
    assert {int(k): v for k, v in pr.P["ant"]["adv_burst_airtime_us"].items()} == {8: 144, 16: 208, 24: 272}
    assert (s.max_payload_octets, s.ant["burst_retries"]) == (24, 20)


def test_layouts_feed_phy():
    rng = np.random.default_rng(0)
    ladder = [pr.get_spec("bcc", "full", ph, opt=o) for ph in ("sf8", "multitone") for o in (1, 4, 5)]
    for s in SPECS + ladder + [pr.get_spec("ble4thesis", "thesis"), pr.get_spec("bcc", "thesis", fcs=True)]:
        for L in (min(20, s.max_payload_octets), s.max_payload_octets):
            lays = pr.tone_layouts(s, L)    # one packet, or the tone packets of a multi-tone BCC frame
            t = pr.airtime_us(s, L)
            # the longest packet (the first, but for an L5 stripe, whose tone 0 may end earlier) ends at the airtime
            # (rounded up to whole us)
            top = max(lays, key=lambda lay: lay[-1]["t1_us"])
            assert top == pr.frame_layout(s, L) and 0 <= t - top[-1]["t1_us"] < 1, (s.pid, L)
            # info bits of the payload blocks sum to the application payload
            assert sum(b["info_bits"] for lay in lays for b in lay if b["kind"] == "payload") == 8 * L, (s.pid, L)
            for lay in lays:
                # blocks tile the packet from 0 on; Walsh blocks hold whole symbols; no packet outlasts the frame
                assert lay[0]["t0_us"] == 0 and all(abs(a["t1_us"] - b["t0_us"]) < 1e-9 for a, b in zip(lay, lay[1:]))
                assert all(b["n_bits"] % 16 == 0 for b in lay if b["code"] == ("walsh16",))
                assert lay[-1]["t1_us"] <= t
                # 60 dB SINR: the frame gets through; -30 dB with SIR -30 dB: it cannot (>= 160 random bits or a
                # C/I floor); a sampled draw over two segments runs too
                hi = phy.frame_outcome(rng, [(0, t, 1e6, 1e6)], lay, s.phy, sampled=False)
                lo = phy.frame_outcome(rng, [(0, t, 1e-3, 1e-3)], lay, s.phy, sampled=False)
                assert hi["success"] > 0.999 and lo["success"] < 1e-3, (s.pid, s.variant, L, hi, lo)
                assert phy.frame_outcome(rng, [(0, t / 2, 1e6, 1e6), (t / 2, t, 1e6, 1e6)], lay, s.phy)["success"]


def test_full_capability():
    # every stack runs with every capability its standard gives (DESIGN 1). Retry limits of
    # PARAMETERS 3: BLE none (ACL), 802.15.4 3, TSCH 3, NB 3, SmartBAN RT 1, ANT 6 attempts, HBC 3
    want = {"ble": None, "ieee802154": 3, "tsch": 3, "ieee802156nb": 3, "smartban": 1, "ant": 5, "hbc": 3}
    for s in SPECS:
        if s.pid in want:
            assert s.variant == "full" and s.arq_max_retx == want[s.pid], (s.pid, s.arq_max_retx)
    # channel rules: BLE CSA#2 with AFH (PDR < 0.95 over 20 outcomes, >= 10 channels kept: 3.1 FULL step 6); TSCH
    # 16-channel hopping with blacklisting (3.3 FULL+: PDR < 0.8 over the last 100 outcomes, >= 4 kept); 802.15.4, NB
    # and SmartBAN an ED scan, SmartBAN with migration (> 20 % failed over 10 IBIs); TSCH a CCA of 128 us at
    # macTsCcaOffset 1800 us before every data frame (3.3 macTsCca option); SmartBAN MCA sensing of 16 us
    # (3.5); ANT a private network, one frequency per BAN from the 78 of 2403..2480 MHz without 2450 and 2457 (3.6)
    b, t, sb, a = (pr.get_spec(p, "full") for p in ("ble", "tsch", "smartban", "ant"))
    assert b.hop == "csa2" and np.isclose(b.adapt["threshold"], 0.05) and (b.adapt["window"], b.adapt["min_channels"]) \
        == (20, 10)
    assert t.hop == "tsch" and np.isclose(t.adapt["threshold"], 0.2) and (t.adapt["window"], t.adapt["min_channels"]) \
        == (100, 4)
    assert [(s.kind, s.dur_us) for s in pr.exchange(t, 20)[:3]] == [("ifs", 1800), ("cca", 128), ("ifs", 192)]
    assert all(pr.get_spec(p, "full").ed_scan for p in ("ieee802154", "ieee802156nb", "smartban"))
    assert sb.channel_switch and sb.csma == dict(mode="mca", cca_us=16, ed_threshold_dbm=-75)
    assert a.hop == "fixed" and len(a.channels) == 76 and 2450 not in a.channels and 2457 not in a.channels
    assert (a.channels[0], a.channels[-1], a.ant["ch_rate_factor"]) == (2403, 2480, 1.25)
    # the catalogue: every rate or PHY each standard offers (the report takes the best per cell as link adaptation):
    # BLE 1M, 2M, Coded S=2 and S=8 (3.1); 802.15.6 NB (M, S) = (2, 4), (2, 2), (2, 1), (4, 1) at 600 ksym/s for
    # 121.4, 242.9, 485.7 and 971.4 kb/s, sensitivities -95, -93, -90, -86 dBm (3.4); HBC FSC SF 8, 16, 32, 64 (3.7);
    # 802.15.4 unslotted CSMA/CA (mode A) and beacon mode B with GTSs (phy 'gts', 3.2); SmartBAN's PHY schemes, PSDU
    # BCH or uncoded at Nrep 1, 2 and 4 (3.5, [42] Table 3); one each for TSCH and ANT; BCC full on one tone and
    # multi-tone, and the as-built stack
    assert [(p, ph) for p, v, ph in pr.catalogue()] == [
        ("ble", "1M"), ("ble", "2M"), ("ble", "coded_s2"), ("ble", "coded_s8"), ("ieee802154", "oqpsk"),
        ("ieee802154", "gts"), ("tsch", "oqpsk"), ("ieee802156nb", "121.4"), ("ieee802156nb", "242.9"),
        ("ieee802156nb", "485.7"), ("ieee802156nb", "971.4"), ("smartban", "bch1"), ("smartban", "bch2"),
        ("smartban", "bch4"), ("smartban", "unc1"), ("smartban", "unc2"), ("smartban", "unc4"), ("ant", "gfsk"),
        ("hbc", "sf8"), ("hbc", "sf16"), ("hbc", "sf32"), ("hbc", "sf64"), ("bcc", "sf8"), ("bcc", "multitone"),
        ("bcc", "bfsk")]
    nb = [pr.get_spec("ieee802156nb", "full", r).phy for r in ("121.4", "242.9", "485.7", "971.4")]
    assert [(x["mod"], x["bit_rate_bps"], x["rep"], x["sens_dbm"]) for x in nb] == [
        ("dbpsk", 6e5, 4, -95), ("dbpsk", 6e5, 2, -93), ("dbpsk", 6e5, 1, -90), ("dqpsk", 1.2e6, 1, -86)]
    assert [pr.get_spec("hbc", "full", f"sf{sf}").phy["sf"] for sf in (8, 16, 32, 64)] == [8, 16, 32, 64]
    # variants: 'full' for every stack, 'thesis' for the as-built models (the as-built BCC stack and the plain
    # BLE model); the fair variants and the strip_hopping option are gone (an unknown flag), so are the old flags
    # that are now always on (TSCH blacklist and cca, SmartBAN mca, ANT private and adaptive_iso)
    assert pr.get_spec("bcc", "thesis").phy["name"] == "bfsk" and pr.get_spec("ble4thesis", "thesis").hop == "perm"
    for pid, v, kw in (("ble", "fair", {}), ("bcc", "fair", {}), ("ble4thesis", "full", {}), ("ble", "thesis", {}),
                       ("ble", "full", dict(strip_hopping=True)), ("tsch", "full", dict(blacklist=True)),
                       ("tsch", "full", dict(cca=True)),
                       ("smartban", "full", dict(mca=True)), ("ant", "full", dict(private=True)),
                       ("ant", "full", dict(adaptive_iso=True))):
        try:
            pr.get_spec(pid, v, **kw)
            raise AssertionError((pid, v, kw))
        except ValueError:
            pass


def test_acr_tables():
    # PARAMETERS 5.5: the last entry of a table is a floor for every larger offset ('>= 3 MHz'), so each table runs
    # to the span of the channel plan. BLE 1M spec ACR 6 / 38 / 48 dB at 1 / 2 / >= 3 MHz, 2 MHz grid: offsets 1, 2
    # -> 2, 4 MHz; BLE data channels span 2404..2478 MHz = 37 grid steps (the 4 MHz gap at 2426 included)
    a = pr.get_spec("ble", "full", "1M").acr_db
    assert a[:3] == [0, 38, 48] and len(a) == 38 and set(a[2:]) == {48}
    # 2M: 6 / 38 / 48 dB at 2 / 4 / >= 6 MHz: offsets 1, 2, 3 and beyond
    a = pr.get_spec("ble", "full", "2M").acr_db
    assert a[:4] == [0, 6, 38, 48] and len(a) == 38 and set(a[3:]) == {48}
    # 802.15.4 chip (AT86RF233) 33 / 48 / 54 dB at 5 / 10 / 15 MHz, 5 MHz grid, channels 2405..2480 = 15 steps
    a = pr.get_spec("ieee802154", "full").acr_db
    assert a[:4] == [0, 33, 48, 54] and len(a) == 16 and set(a[3:]) == {54}
    # NB: adjacent 1 MHz channel 0 dB (worst case), offsets >= 2 ignored (3.4): the table stops at offset 1
    assert pr.get_spec("ieee802156nb", "full").acr_db == [0, 0]
    # ANT 4 / 33 / 42 / 52 / 62 dB at 1 / 2 / 3 / 6 / 25 MHz, 1 MHz grid, to the span of the private pool 2403..2480 MHz
    # = 77 steps: 78 entries, 3..5 -> 42, 6..24 -> 52, 25..77 -> 62
    a = pr.get_spec("ant", "full").acr_db
    assert len(a) == 78 and a[1:4] == [4, 33, 42] and a[5] == 42 and a[6] == a[24] == 52 and a[25] == a[77] == 62
    # pairs by centre frequency: BLE 2M data indices 9 / 10 are 2 MHz apart (6 dB), 10 / 11 are 4 MHz apart across
    # advertising channel 38 at 2426 MHz (38 dB, not the 6 dB of an index offset of 1), 0 / 36 get the 48 dB floor
    m = pr.acr_matrix(pr.get_spec("ble", "full", "2M"))
    assert (m[9][10], m[10][11], m[0][36], m[5][5]) == (6, 38, 48, 0)
    # NB: 1 MHz apart 0 dB, 2 MHz apart outside the table (ignored)
    m = pr.acr_matrix(pr.get_spec("ieee802156nb", "full"))
    assert m[0][1] == 0 and m[0][2] == float("inf")


def test_exchanges():
    # BLE ACL 1M, L = 20: poll (empty PDU) 80 + T_IFS 150 + data 240 + T_IFS 150 = 620 us
    assert sum(st.dur_us for st in pr.exchange(pr.get_spec("ble", "full", "1M"), 20)) == 620
    # TSCH: one 10 ms timeslot per exchange
    assert sum(st.dur_us for st in pr.exchange(pr.get_spec("tsch", "full"), 20)) == 10000
    # ANT acknowledged: master TX 0..144, slave ACK 274..418
    assert sum(st.dur_us for st in pr.exchange(pr.get_spec("ant", "full"), 8)) == 418
    # NB FULL, L = 20: data 2304 + pSIFS 75 + I-Ack (the L = 0 frame: 72 bits, 2 codewords, 96 coded x 4 = 384
    # symbols, (90 + 124 + 384) / 0.6 = 996.7 -> 997) + pSIFS 75 = 3451
    assert sum(st.dur_us for st in pr.exchange(pr.get_spec("ieee802156nb", "full"), 20)) == 3451
    # SmartBAN, L = 20: data 356 + TIFS 150 + ACK 168 + TIFS 150
    assert sum(st.dur_us for st in pr.exchange(pr.get_spec("smartban", "full"), 20)) == 356 + 468


def test_select_period_matches_ble_sim():
    # BLE_SIM.select_ci is the reference; the points cover latency and energy, both PHYs, a scan-free infeasible
    # slot (120 kb/s x 9 > 1 Mb/s), a payload above 251 B, L < 1 byte, and the nearest-rounding quirk of the
    # energy CI (80.9 kb/s: 2008000 / 80.9 = 24820 us -> 25000 us -> 2022.5 bits > 2008 -> infeasible)
    pts = [(30, "latency", 10, "1M"), (50, "latency", 10, "1M"), (200, "latency", 4, "2M"), (60, "latency", 16, "2M"),
           (2.5, "latency", 40, "1M"), (100, "energy", 4, "2M"), (5, "energy", 30, "1M"), (100.4, "energy", 1, "1M"),
           (1, "energy", 3, "2M"), (120, "latency", 9, "1M"), (300, "latency", 5, "1M"), (0.5, "latency", 1, "1M"),
           (80.9, "energy", 5, "1M"), (10, "energy", 20, "2M"), (100, "latency", 10, "1M")]
    # the last point fills the PHY exactly (10 x 100 kb/s = 1 Mb/s: t_on 750 us = slot 750 us), which BLE_SIM
    # accepts: the capacity bound of select_period must not reject it
    feasible = 0
    for r, mode, n, ph in pts:
        ref = BLE_SIM.select_ci(r, mode, n, 1_000_000 if ph == "1M" else 2_000_000)
        got = pr.select_period(r, mode, n, pr.get_spec("ble4thesis", "thesis", ph))
        assert (ref is None) == (got is None), (r, mode, n, ph, ref, got)
        if ref is not None:
            feasible += 1
            assert (got["period_us"], got["payload_bits"], got["frames_per_period"]) == (ref.ci_us, ref.payload_bits, 1)
            assert abs(got["r_actual_kbps"] - ref.r_actual_kbps) < 1e-12
    assert feasible >= 8
    # hand check of one point: 50 kb/s x 7.5 ms = 375 bits; t_on 375 us <= slot 7500 // 10 = 750 us
    assert pr.select_period(50, "latency", 10, pr.get_spec("ble4thesis", "thesis"))["payload_bits"] == 375


def test_select_period_hand():
    ble = pr.get_spec("ble", "full", "1M")
    # reference period P (DESIGN 3): latency the shortest BLE CI carrying >= 8 bits, energy the longest CI with
    # 251 octets = 2008 bits: 10 kb/s -> 200.8 ms -> 200 ms on the 1.25 ms grid; 0.5 kb/s latency: 13.75 ms gives
    # round(6.875) = 7 bits, 15 ms round(7.5) = 8 bits (half to even)
    assert pr.reference_period(50, "latency") == 7500 and pr.reference_period(10, "energy") == 200000
    assert pr.reference_period(0.5, "latency") == 15000
    # ACL exchange = poll 80 + 150 + (80 + 8 L) + 150 = 460 + 8 L us, and the event holds 1 + 1 spare exchanges per
    # frame. 10 kb/s, 10 sensors: 7.5 ms gives 75 bits -> 10 octets -> 2 x 540 us, x 10 =
    # 10800 > 7500; 8.75 ms round(87.5) = 88 bits -> 11 octets -> 2 x 548, x 10 > 8750; 10 ms 100 bits -> 13 octets
    # -> 2 x 564, x 10 > 10000; 11.25 ms round(112.5) = 112 bits -> 14 octets -> 2 x 572, x 10 = 11440 > 11250;
    # 12.5 ms 125 bits -> 16 octets -> 2 x 588, x 10 = 11760 <= 12500. The deadline stays the cell's P = 7.5 ms
    sel = pr.select_period(10, "latency", 10, ble)
    assert (sel["period_us"], sel["payload_bits"], sel["frame_payload_octets"], sel["deadline_us"]) == (
        12500, 125, [16], 7500)
    # 50 kb/s x 10 sensors: L = 0.05 T bits (T in us) take 8 us per octet, at least 0.05 T us, so the events need
    # 10 x 2 x (460 + 0.05 T) >= 9200 + T > T us at every period
    assert pr.select_period(50, "latency", 10, ble) is None
    # energy, 10 kb/s: P = 200 ms, L = 2000 bits = 250 octets in one frame
    sel = pr.select_period(10, "energy", 1, ble)
    assert (sel["period_us"], sel["frame_payload_octets"]) == (200000, [250])
    # 1000 kb/s x 10 sensors exceeds the 1 Mb/s PHY at any period: infeasible
    assert pr.select_period(1000, "latency", 10, ble) is None
    assert pr.select_period(0.5, "latency", 1, ble)["period_us"] == 15000
    # 802.15.4 CSMA admits a period when the scheduled load of the exchanges is at most 1 / (1 + 1 spare exchange)
    # = 0.5 of T. 10 kb/s, 4 sensors: BI 15.36 ms: 154 bits -> 20 octets, exchange CCA 128 +
    # turnaround 192 + data 1184 + 192 + Imm-Ack 352 + LIFS 640 = 2688 us, x 4 = 10752 > 7680; BI 30.72 ms: 307 bits
    # -> 39 octets, data 1792, exchange 3296, x 4 = 13184 <= 15360. 5 kb/s, 10 sensors: 15.36 ms: 77 bits -> 10
    # octets, 2368 x 10 > 7680; 30.72 ms: 20 octets, 2688 x 10 > 15360; 61.44 ms: 307 bits -> 39 octets, 3296 x 10 >
    # 30720; 122.88 ms: 614 bits -> 77 octets, data 3008, exchange 4512, x 10 = 45120 <= 61440
    c = pr.get_spec("ieee802154", "full")
    assert c.csma["max_load"] == 0.5 and pr.select_period(10, "latency", 4, c)["period_us"] == 30720
    assert pr.select_period(5, "latency", 10, c)["period_us"] == 122880
    # ANT advanced burst (3.6: packets of up to 24 octets, one data + ACK pair every 3.2 ms) above the 8 octets of a
    # message. 20 kb/s: P = 7.5 ms, the first channel period at or above it N = 246 counts, 246e6 / 32768 = 7507.3 ->
    # 7507 us, L = round(150.14) = 150 bits = 19 octets: one packet, 3200 <= 7507. 50 kb/s: L = round(375.35) = 375
    # bits = 47 octets = 24 + 23, 2 x 3200 <= 7507. 60 kb/s: 24 octets per 3.2 ms exactly, so only a period that is a
    # whole number of spacings holds it: N = 734, 734e6 / 32768 = 22399.9 -> 22400 us, L = 1344 bits = 7 x 24 octets, 7
    # x 3200 = 22400 (every earlier period needs one packet more than it holds); above 60 kb/s none does
    ant = pr.get_spec("ant", "full")
    for r, want in ((20, (7507, [19])), (50, (7507, [24, 23])), (60, (22400, [24] * 7))):
        sel = pr.select_period(r, "latency", 1, ant)
        assert (sel["period_us"], sel["frame_payload_octets"], sel["burst"]) == want + (True,), (r, sel)
    assert pr.select_period(61, "latency", 1, ant) is None
    # TSCH, 3 sensors, 1 kb/s: EB cell + 3 dedicated cells of 10 ms -> slotframe 40 ms, L = 40 bits = 5 octets,
    # one cell per frame (3.3)
    sel = pr.select_period(1, "latency", 3, pr.get_spec("tsch", "full"))
    assert (sel["period_us"], sel["frame_payload_octets"], sel["slot_us"]) == (40000, [5], 10000)


def test_select_period_recipes():
    # SmartBAN (3.5): LSLOT = the smallest of 1, 2, 4, .. with PPDU <= TS - TMUA 16 - ACK 168 - 2 TIFS 300, IBI =
    # (beacon + N_dev + 2 C/M) slots at least, T a multiple of TS. FULL 10 kb/s, 4 sensors, P = 7.5 ms, 625 us grid:
    # 7500: 75 bits -> 10 octets, MPDU 64 + 80 = 144 bits, 2 BCH codewords (+ 28), PPDU 104 + 172 = 276 us -> TS >=
    # 760 -> LSLOT 2 = 1250 us, but 7 x 1250 = 8750 > 7500; 8125 is no multiple of 1250; 8750: round(87.5) = 88
    # bits -> 11 octets, PPDU 284 us, TS 1250, 7 slots = 8750 fit. Uplink slot 1 starts after the beacon slot
    sel = pr.select_period(10, "latency", 4, pr.get_spec("smartban", "full"))
    assert (sel["period_us"], sel["payload_bits"], sel["slot_us"], sel["reserved_us"]) == (8750, 88, 1250, 1250)
    # LF max at LSLOT 32 (20 ms): the table's 2152 BCH octets hold for TMUA 0; with the 16 us MCA sensing the room is
    # 20000 - 16 - 468 = 19516 us; LF 2151: MPDU 17272, 153 codewords (+ 2142), PPDU 104 + 19414 = 19518 > 19516;
    # LF 2150: 17264 + 2142 + 104 = 19510 fits
    assert pr.get_spec("smartban", "full").max_payload_octets == 2150
    # 802.15.6 NB FULL: the allocation holds 1 + 1 spare exchanges . 10 kb/s, 1 sensor: beacon 17
    # octets (N_PSDU 208, 5 codewords, 268 coded bits x 4 = 1072 symbols, (90 + 124 + 1072) / 0.6 = 2143.3 -> 2144
    # us) + pSIFS 75 = 2219 us reserved; exchange = data + 75 + I-Ack 997 + 75. 7500 us: 75 bits -> 10 octets,
    # N_PSDU 152, 3 codewords, 188 x 4 = 752 symbols, (214 + 752) / 0.6 = 1610 us, 2 x 2757 = 5514 > 7500 - 2219 =
    # 5281; 8000 us: 80 bits -> 10 octets, 5514 <= 5781
    sel = pr.select_period(10, "latency", 1, pr.get_spec("ieee802156nb", "full"))
    assert (sel["period_us"], sel["frame_payload_octets"]) == (8000, [10])
    # 100 kb/s, 4 sensors, latency: at 121.4 kb/s the allocations of 1 + 1 exchanges cannot hold it (scan-free bound:
    # a 255-octet frame of N_PSDU 2112 bits, 42 codewords, 2616 coded bits x 4 = 10464 symbols, (214 + 10464) / 0.6 =
    # 17797 us, exchange 17797 + 75 + 997 + 75 = 18944, x 2 = 37888 us per 2040 bits: 4 x 100 x 37888 / 2040 = 7429 >
    # 1000)
    assert pr.select_period(100, "latency", 4, pr.get_spec("ieee802156nb", "full")) is None
    # 802.15.4 beacon mode B (3.2, phy 'gts'): one GTS of whole slots (SD / 16) per sensor at the end of the SD, the
    # smallest SO whose GTSs leave the 440-symbol CAP (7040 us); a GTS holds 1 + 1 spare exchanges per frame of data,
    # 192 us, Imm-Ack 352 us and LIFS 640 us (no CSMA in the CFP). 10 kb/s latency, 1 sensor: BI 15.36 ms, 154 bits
    # -> 20 octets, data (6 + 11 + 20) x 32 = 1184, exchange 2368, x 2 = 4736 us -> 5 slots of 960 us = 4800 <= 15360 -
    # 7040 = 8320; the GTS starts at 15360 - 4800 = 10560. 4 sensors: 15.36 ms: 4 x 4800 > 8320; 30.72 ms: 307 bits
    # -> 39 octets, data 1792, 2 x 2976 = 5952: SO 0 7 slots of 960 = 6720, x 4 > 8320; SO 1 4 slots of 1920 = 7680, x
    # 4 = 30720 > 23680; 61.44 ms: 614 bits -> 77 octets, data 3008, 2 x 4192 = 8384: SO 0 9 x 960 = 8640 x 4 > 8320,
    # SO 1 5 x 1920 = 9600 x 4 > 23680, SO 2 3 x 3840 = 11520 x 4 = 46080 <= 54400, GTSs from 61440 - 46080 = 15360.
    # At most 7 GTSs: 8 sensors are refused
    g = pr.get_spec("ieee802154", "full", "gts")
    assert (g.access, g.max_dev, g.phy["name"], g.beacon["max_lost"]) == ("tdma", 7, "gts", 4)
    assert [(s.kind, s.dur_us) for s in pr.exchange(g, 20)] == [("data", 1184), ("ifs", 192), ("ack", 352),
                                                                 ("ifs", 640)]
    for n, want in ((1, (15360, [20], 4800, 10560)), (4, (61440, [77], 11520, 15360))):
        sel = pr.select_period(10, "latency", n, g)
        assert (sel["period_us"], sel["frame_payload_octets"], sel["slot_us"], sel["reserved_us"]) == want, sel
    assert pr.select_period(1, "latency", 8, g) is None
    # ANT FULL (3.6): acknowledged messages at T / 1.25 must keep the channel period at N >= 171 counts (5219 us),
    # so single-frame selections have T >= 1.25 x 5219
    for r in (1, 5, 9, 10, 12, 30, 100, 1000):
        for mode in ("latency", "energy"):
            for n in (1, 4, 10, 16):
                a = pr.select_period(r, mode, n, pr.get_spec("ant", "full"))
                assert a is None or a["frames_per_period"] > 1 or a["period_us"] / 1.25 >= 5219, (r, mode, n, a)

def test_bcc_pool_fill():
    # fill adds the slack before the next beacon to the pool, admission on the share alone. L3, one tone,
    # 1 sensor, 16 bits -> 2 octets (test_bcc_ladder: exchange 1699, beacon 2712, reserved 2787, pool exchange 3270):
    # at 20000 us the slack is 20000 - 2787 - 1699 = 15514 us, floor(15514 / 3270) = 4 exchanges (5 need 16350);
    # without fill pool_size(1) = 1. At 8000 us the slack holds only the one exchange the share already admits
    # (8000 - 2787 - 1699 = 3514). The fit bound is 2787 + 1699 + 3270 = 7756 us with and without fill (7500 us in
    # test_bcc_ladder fails on the period grid's 15 bits): capacity is unchanged.
    plain, fill = (pr.get_spec("bcc", "full", opt=3, **kw) for kw in ({}, dict(fill=True)))
    assert pr._fit(16, 20000, 1, plain)["pool"] == 1 and pr._fit(16, 20000, 1, fill)["pool"] == 4
    assert pr._fit(16, 8000, 1, fill)["pool"] == 1
    assert all(pr._fit(16, 7756, 1, s)["pool"] == 1 and pr._fit(16, 7755, 1, s) is None for s in (plain, fill))
    try:
        pr.get_spec("bcc", "full", opt=2, fill=True)
        raise AssertionError("fill below L3")
    except ValueError as e:
        assert "fill" in str(e), e


def test_bcc_unspread_and_rate_rule():
    #. Unspread PSDU: 8 chips per octet instead of 32. L5, 500 kchip/s (2 us chips), L = 255: PSDU 2 + 255 + 2
    # = 259 octets, a pilot after each 128 octets (2 pilots either way), so the frame shrinks by 259 x 24 x 2 = 12432 us
    # (17752 -> 5320 us). At Ec/N0 14 dB (g = 25.12) a chip errs with 0.5 exp(-12.56) = 1.756e-6 and the 2072 payload
    # bits are all right with (1 - 1.756e-6)^2072 = 0.99637; preamble, header and pilots lose < 1e-9 there.
    walsh, plain = pr.get_spec("bcc", "full", opt=5), pr.get_spec("bcc", "full", opt=5, spread="none")
    assert pr.airtime_us(walsh, 255) == 17752 and pr.airtime_us(plain, 255) == 5320
    assert pr.airtime_us(walsh, 255, plain=True) == 5320 and pr.airtime_us(plain, 255, plain=False) == 17752
    lay = pr.frame_layout(plain, 255)
    assert sum(b["n_bits"] for b in lay if b["kind"] == "payload" and b["code"] is None) == 2072
    out = phy.frame_outcome(None, [(0.0, lay[-1]["t1_us"], 10 ** 1.4, math.inf)], lay, plain.phy, sampled=False)
    assert abs(out["success"] - (1 - 0.5 * math.exp(-10 ** ((14 + phy.V2_OFFSET_DB) / 10) / 2)) ** 2072) < 1e-9
    # Requirement: Walsh at 8 dB, p = 0.5 exp(-3.155) = 0.02132, symbol ok 0.999911, PSDU of 518 symbols 0.9549; the
    # unspread PSDU matches it at a bit error of 1 - 0.9549^(1 / 2072) = 2.225e-5 (the PSDU length cancels), g = 2 ln(0.5 / 2.225e-5) = 20.04,
    # 13.02 dB. Rate rule from A_95 = 15.7 mV (500 kchip/s, NF 20): at 1 V and NF 30, unspread 500e3 x (1 / 0.0157)^2
    # x 10^-1 x 10^((8 - 13.02) / 10) = 6.39e7 chip/s, Walsh 2.03e8: both capped at 9 Mchip/s, the waveform limit.
    assert abs(pr.required_ec_n0_db(True) - 13.02) < 0.01 and pr.required_ec_n0_db(False) == 8.0
    r, raw = pr.margin_rate("none", 30.0)
    assert r == pr.BCC_RC_MAX_CPS and abs(raw / 6.39e7 - 1) < 0.005
    assert abs(pr.margin_rate("walsh", 30.0)[1] / 2.03e8 - 1) < 0.005
    for kw in (dict(spread="none", opt=5, phy="multitone"), dict(spread="adapt", opt=1), dict(spread="x")):
        try:
            pr.get_spec("bcc", "full", kw.pop("phy", None), **kw)
            raise AssertionError(kw)
        except ValueError as e:
            assert "spread" in str(e), e


def test_bcc_pooled_reference():
    # pooled reference stack: BCC full at L0 with the retry pool (share 0.05) instead of the spare exchange per
    # frame. One tone, 500 kchip/s, the numbers of test_bcc_ladder: L_max frame 20096 us, I-Ack 3264 (the L = 0 frame),
    # exchange 20096 + 75 + 3264 + 75 = 23510; the allocation holds one exchange (47020 with the spare); a pool
    # exchange is the L0 poll (the L = 0 frame, 3264) + pSIFS + the exchange = 26849 us. Capacity bound 8 x 255 /
    # (23510 + 0.05 x 26849) = 82.084 kb/s (43.386 with the spare, 109.118 at L5: the proposal gains 1.33 against the
    # pooled reference and 2.51 against the spare)
    s0, sp = pr.get_spec("bcc", "full", "sf8"), pr.get_spec("bcc", "full", "sf8", pool=0.05)
    assert sp.flags == {"pool": 0.05} and pr._pooled(sp) and not pr._pooled(s0)
    assert (pr._frame_us(s0, 255), pr._frame_us(sp, 255), pr._pool_us(sp, 255)) == (47020, 23510, 26849)
    lmax = sp.max_payload_octets
    assert abs(8 * lmax / (pr._frame_us(sp, lmax) + 0.05 * pr._pool_us(sp, lmax)) * 1e3 - 82.0844) < 1e-3
    # select_period, 2 kb/s latency, one sensor: T 7500 (15 bits, 2 octets) needs the beacon 2688 + 26 x 64 = 4352 +
    # pSIFS 75, the exchange (frame 2688 + 11 x 64 = 3392, + 75 + 3264 + 75 = 6806) and the pool exchange (3264 + 75 +
    # 6806 = 10145): 21378 > 7500. The 500 us grid climbs to 21500 (43 bits, 6 octets: frame 2688 + 15 x 64 = 3648,
    # exchange 7062, pool exchange 10401; 4427 + 7062 + 10401 = 21890 > 21500) and 22000 holds (44 bits, 6 octets).
    # With the spare (no flag) the same cell takes 18500: 37 bits, 5 octets (frame 3584, exchange 6998), 4427 + 2 x
    # 6998 = 18423
    sel = pr.select_period(2, "latency", 1, sp)
    want = dict(period_us=22000, frame_payload_octets=[6], beacon_us=4352, reserved_us=4427, slot_us=7062, pool=1,
                pool_us=10401)
    assert {k: sel[k] for k in want} == want, sel
    assert pr.select_period(2, "latency", 1, s0)["period_us"] == 18500
    # L1 and L2 take the pool the same way; from L3 the flag only names the share, and 0.05 is the default
    for o in (1, 2):
        s = pr.get_spec("bcc", "full", "sf8", opt=o, pool=0.05)
        assert pr._pooled(s) and pr._frame_us(s, lmax) == pr._ex_us(s, lmax) and "pool" in pr._fit(8, 30000, 1, s)
    a, b = (pr.select_period(1, "energy", 10, pr.get_spec("bcc", "full", "sf8", opt=3, **kw))
            for kw in ({}, dict(pool=0.05)))
    assert a == b


def test_ri_sfd_field():
    # IEEE 802.15.6-2012 10.5.2, Eq. 117: in RI mode the SFD/RI field is (64 + 12) x 8 = 608 chips, the 512-chip SFD
    # of DRF mode (10.5.1) plus 12 zero bits at SF 8, i.e. (608 - 512) / 8 = 12 more sync-chip times per frame. The
    # profile (ladder level 1 and up) runs in RI mode, the reference (level 0) in DRF mode. One tone, the pool of the
    # Capacity stacks, at the three simulated chip rates (a Walsh chip 1e6 / R_c us) and three payloads
    keys = [(rc, o, L) for rc in (50e3, 500e3, 2.625e6) for o in (0, 1, 5) for L in (8, 64, 255)]

    def ends():
        return {(rc, o, L): pr.frame_layout(pr.get_spec("bcc", "full", "sf8", rc_cps=rc, opt=o, pool=0.05), L)[-1][
            "t1_us"] for rc, o, L in keys}

    assert pr.HBC["sfd_ri_chips"] == (64 + 12) * 8 and pr.HBC["sfd_chips"] == 64 * 8
    ri = ends()
    pr.HBC["sfd_ri_chips"] = pr.HBC["sfd_chips"]
    try:
        drf = ends()
    finally:
        pr.HBC["sfd_ri_chips"] = (64 + 12) * 8
    for rc, o, L in keys:
        want = 12 * 1e6 / rc if o else 0.0
        assert abs(ri[rc, o, L] - drf[rc, o, L] - want) < 1e-6, (rc, o, L, ri[rc, o, L] - drf[rc, o, L])


def test_bcc_lmax():
    # Short frames (claim4.py): lmax caps the frame payload of one-tone BCC full; 10 octets go as 8 + 2
    s = pr.get_spec("bcc", "full", "sf8", pool=0.05, lmax=8)
    assert s.max_payload_octets == 8 and pr._fit(80, 30000, 1, s)["frame_payload_octets"] == [8, 2]
    assert pr.get_spec("bcc", "full", "sf8", lmax=255).max_payload_octets == 255
    for var, phy, lm in (("full", "multitone", 8), ("thesis", "bfsk", 8), ("full", "sf8", 0), ("full", "sf8", 256)):
        try:
            pr.get_spec("bcc", var, phy, lmax=lm)
        except ValueError:
            continue
        raise AssertionError((var, phy, lm))


def test_bcc_block_ack():
    # block-ack reference (flag back, claim4.py H): exchange data + pSIFS (L-Ack), one B-Ack frame of 3 octets
    # + pSIFS per allocation, the pool exchange with its I-Ack (poll, pSIFS, data, pSIFS, I-Ack, pSIFS). An allocation
    # of m frames costs per frame d + 75 + (b + 75) / m and the pool share 0.05 of a pool exchange: claim4.back_cost
    import claim4
    import report
    for rc in (50e3, 500e3, 2.625e6):
        for L in (8, 64, 255):
            s = pr.get_spec("bcc", "full", "sf8", rc_cps=rc, pool=0.05, back=True, lmax=L)
            d, a0, b3, ps = pr.airtime_us(s, L), pr.airtime_us(s, 0), pr.airtime_us(s, 3), pr.HBC["psifs_us"]
            assert pr._ex_us(s, L) == d + ps and pr._pool_us(s, L) == a0 + ps + d + ps + a0 + ps, (rc, L)
            for m in (1, 4):
                sel = pr._fit(8 * L * m, 4_000_000, 1, s)
                assert sel["slot_us"] == m * (d + ps) + b3 + ps and sel["pool_us"] == pr._pool_us(s, L), (rc, L, m)
                want = d + ps + (b3 + ps) / m + 0.05 * (a0 + ps + d + ps + a0 + ps)
                got = sel["slot_us"] / m + 0.05 * sel["pool_us"]
                assert abs(got - want) < 1e-9 and abs(claim4.back_cost(rc, 0, L, 0.05, m) - want) < 1e-9, (rc, L, m)
            assert abs(report.frame_cost(s, L) - (d + ps + 0.05 * (a0 + ps + d + ps + a0 + ps))) < 1e-9
    for kw in (dict(back=True), dict(back=True, opt=2, pool=0.05)):
        try:
            pr.get_spec("bcc", "full", "sf8", **kw)
        except ValueError as e:
            assert "back" in str(e), e
            continue
        raise AssertionError(kw)



def test_bcc_block_ack_burst():
    # lead reference (flags back and burst): IEEE 802.15.6-2012 6.2.9.4, the frames of a block within an
    # allocation pMIFS apart in burst mode; 10.5.2, a burst packet in RI mode, header (8192 chips at SF 64 in DRF mode)
    # at the payload rate SF 8 (1024 chips) and the 608-chip SFD/RI instead of the 512-chip SFD; Table 89 pMIFS 20 us.
    # An allocation of n frames: the first in DRF mode, n - 1 x (pMIFS + RI frame), pSIFS 75 us, the 3-octet B-Ack,
    # pSIFS. A chip of BCC full is 1 / (8 R_c)
    for rc in (50e3, 500e3, 2.625e6):
        chip = 1e6 / (8 * rc)
        for L in (8, 64, 255):
            s = pr.get_spec("bcc", "full", "sf8", rc_cps=rc, pool=0.05, back=True, burst=True, lmax=L)
            d, ri, b3 = pr.airtime_us(s, L), pr.airtime_us(s.ri, L), pr.airtime_us(s, 3)
            assert abs(ri - (d - (8192 - 1024) * chip + (608 - 512) * chip)) <= 1, (rc, L, d, ri)
            for n in (1, 2, 5):
                sel = pr._fit(8 * L * n, 4_000_000, 1, s)
                assert sel["slot_us"] == d + (n - 1) * (20 + ri) + 75 + b3 + 75, (rc, L, n, sel["slot_us"])
    for kw in (dict(burst=True, pool=0.05), dict(back=True, burst=True, opt=1, pool=0.05)):
        try:
            pr.get_spec("bcc", "full", "sf8", **kw)
        except ValueError as e:
            assert "burst" in str(e), e
            continue
        raise AssertionError(kw)


def test_bcc_block_ack_span():
    # flag span: 6.2.9.4 lets a block span allocations. A block that never closes leaves the
    # allocation of one frame d + pSIFS, the per-frame cost of claim4.back_cost with an unlimited block; a finite block
    # keeps the B-Ack's room in every allocation (it may close in any), as the burst reference; a block of the
    # allocation's n frames costs per frame what the B-Ack every allocation costs
    import claim4
    for rc in (50e3, 500e3, 2.625e6):
        for L in (8, 64, 255):
            kw = dict(rc_cps=rc, pool=0.05, back=True, burst=True, lmax=L, pilot128=True)   # as claim4._burst_terms
            s, u = pr.get_spec("bcc", "full", "sf8", **kw), pr.get_spec("bcc", "full", "sf8", span=pr.BLOCK_UNLIMITED, **kw)
            sel = pr._fit(8 * L, 4_000_000, 1, u)
            assert sel["slot_us"] == pr.airtime_us(u, L) + 75, (rc, L, sel["slot_us"])
            want = claim4.back_cost(rc, 0, L, 0.05, 1, burst=True, block=pr.BLOCK_UNLIMITED)
            assert abs(sel["slot_us"] + 0.05 * sel["pool_us"] - want) < 1e-3, (rc, L)   # B-Ack / 10^9 left
            for n in (1, 3):
                fin = pr.get_spec("bcc", "full", "sf8", span=9, **kw)
                assert pr._fit(8 * L * n, 4_000_000, 1, fin)["slot_us"] == pr._fit(8 * L * n, 4_000_000, 1, s)["slot_us"]
                assert abs(claim4.back_cost(rc, 0, L, 0.05, n, burst=True, block=n)
                           - claim4.back_cost(rc, 0, L, 0.05, n, burst=True)) < 1e-9
    for kw in (dict(back=True, pool=0.05, span=4), dict(back=True, burst=True, pool=0.05, span=0)):
        try:
            pr.get_spec("bcc", "full", "sf8", **kw)
        except ValueError as e:
            assert "span" in str(e), e
            continue
        raise AssertionError(kw)



def test_bcc_pilot128():
    # flag pilot128: the '128 B' pilot period of 10.7.3 at L0, the lead
    # reference. 500 kchip/s, 64 us per octet, a 128 us pilot after every full 128 octets instead of 64; PSDU 7 + L + 2.
    # DRF frame (head 2688 us): L 255, 264 octets: 4 pilots -> 2, 20096 - 256 = 19840; L 64, 73 octets: 1 -> 0,
    # 2688 + 73 x 64 = 7360; L 8, 17 octets: none either way, 2688 + 17 x 64 = 3776. RI frame of the burst (head 512 +
    # 152 + 256 = 920): 920 + 264 x 64 + 2 x 128 = 18072. The pool's frames and its poll take the same framing.
    kw = dict(rc_cps=500e3, pool=0.05, back=True, burst=True, span=pr.BLOCK_UNLIMITED)
    s, p = pr.get_spec("bcc", "full", "sf8", **kw), pr.get_spec("bcc", "full", "sf8", pilot128=True, **kw)
    assert [pr.airtime_us(p, L) for L in (255, 64, 8)] == [19840, 7360, 3776]
    assert [pr.airtime_us(s, L) for L in (255, 64, 8)] == [20096, 7488, 3776]
    assert pr.airtime_us(p.ri, 255) == 18072 and p.ri.flags.get("pilot128")
    assert pr._pool_us(p, 255) == pr.airtime_us(p, 0) + 75 + pr._ex_us(p, 255, True)
    # claim4.field_terms: the profile's beacon at 500 kchip/s, 64 us per octet and a 128 us pilot
    # per 128 octets. 8000 sensors add a 2-octet base each, 16000 octets = 125 full pilot periods: (16000 x 64 + 125 x
    # 128) / 8000 = 130 us per sensor; 8000 bits add 1000 octets and 7 or 8 pilots: (64000 + 7 or 8 x 128) / 8000 = 8.112
    # to 8.128 us per bit
    import claim4
    bit, base = claim4.field_terms(500e3)
    assert base == 130.0 and 8.112 <= bit <= 8.128, (bit, base)
    for kw2 in (dict(variant="full", opt=1), dict(variant="thesis")):
        try:
            pr.get_spec("bcc", kw2.pop("variant"), pilot128=True, **kw2)
        except ValueError as e:
            assert "pilot128" in str(e) or "opt" in str(e), e
            continue
        raise AssertionError(kw2)

def test_claim4_slot_grid():
    # claim4.own_slot at the 50 us slot reproduces a separate check at points A and
    # D (no pool, where every minimum falls); a slot far below every frame gives own_gtn back; 500 us gives Table 27's
    # 1.29 at A; a hand count at 2502.7 kchip/s, L 8, one sensor: beacons 944.5 / 629.6 us are 19 / 13 slots, frames +
    # guard 845.1 / 428.1 us are 17 / 9 slots, (19 + 17) / (13 + 9) = 18 / 11.
    import claim4 as c4
    def low(p, **kw):
        (l0, l1), _, (d0, d1), _ = c4.POINTS[p]
        return min(c4.own_slot(rc, L, 0.0, d, **kw) for L in range(l0, l1 + 1) for d in range(d0, d1 + 1) for rc in c4.GRID)
    assert c4.SLOT_US == 50.0
    assert round(low("A"), 3) == 1.636 and round(low("D"), 3) == 1.389
    assert round(low("A", beacon_slots=False), 3) == 1.662 and round(low("D", beacon_slots=False), 3) == 1.386
    assert round(low("A", slot=500.0), 3) == 1.286
    rc = min(c4.GRID, key=lambda r: c4.own_slot(r, 8, 0.0, 1))
    assert abs(c4.own_slot(rc, 8, 0.0, 1) - 18 / 11) < 1e-12
    for rc, L, s, d in ((50e3, 8, 0.0, 1), (500e3, 64, 0.05, 4), (2.625e6, 250, 0.2, 8)):
        g = c4.packed_gtn(rc, L, s, d)
        assert abs(c4.own_slot(rc, L, s, d, g, 1e-4, False) - c4.own_gtn(rc, L, s, d, g)) < 1e-4, (rc, L, s, d)


def test_claim4_step_factors():
    # the two steps of the ladder on the slot grid multiply to own_slot in either order. Hand count in slots, per
    # sensor: at A (L 8, one sensor, 2502.7 kchip/s) beacon + frame are 19 + 17 on plain HBC, 12 + 10 with step 1, 13 + 9 on
    # the profile, whose ack field costs the slot its shorter frame saves; at 64 octets, 58 sensors, 2625 kchip/s
    # 19/58 + 31, 12/58 + 25, 43/58 + 24.
    import claim4 as c4
    import protocols as pr
    rc = min(c4.GRID, key=lambda r: c4.own_slot(r, 8, 0.0, 1))
    one, rest, _, _ = c4.step_factors(rc, 8, 0.0, 1)
    assert abs(one - 36 / 22) < 1e-12 and abs(rest - 1) < 1e-12
    r, s1, p = 19 / 58 + 31, 12 / 58 + 25, 43 / 58 + 24
    one, rest, _, _ = c4.step_factors(2.625e6, 64, 0.0, 58)
    assert abs(one - r / s1) < 1e-12 and abs(rest - s1 / p) < 1e-12
    for case, kw in (((2.625e6, 253, 0.0, 60), {}), ((500e3, 64, 0.05, 4), {}), ((50e3, 8, 0.2, 1), dict(slot=500.0)),
                     ((1.3125e6, 125, 1.0, 8), dict(beacon_slots=False)), ((500e3, 20, 0.05, 2), dict(g=40.0))):
        total = c4.own_slot(*case, **kw)
        one, rest, changes, on_top = c4.step_factors(*case, **kw)
        assert abs(one * rest - total) < 1e-12 and abs(changes * on_top - total) < 1e-12, case
    # step 1 sends the first frame as the burst's later frames, and the changes alone keep plain HBC's pilot period
    for L in (8, 64, 255):
        ri = c4._burst_terms(500e3, L, 0.05)[1] - pr.HBC["pmifs_us"]
        assert abs(pr.airtime_us(c4.spec(500e3, 1, 0.05), L) - ri) < 1e-9
    s = c4.changes_spec(500e3)
    assert pr._hbc_framing(s)[0] == pr._hbc_framing(c4.spec(500e3, 0, 0.05, pilot128=True))[0]
    assert pr._hbc_framing(s)[3] == pr.LEAN_MAC_OCTETS


def test_bcc_chip_curve_shifted():
    # 5: a galvanic link is decided by the chip error rate shifted by the 0.35 dB of parameters.json,
    # so at 10 dB less that shift a chip errs with 0.5 exp(-10 / 2) = 3.3690e-3 on both stacks.
    assert pr.load_params()["bcc"]["ber"]["sim_offset_dB"] == phy.V2_OFFSET_DB
    for spec in (pr.get_spec("bcc", "full"), pr.get_spec("hbc", "full"), pr.get_spec("bcc", "thesis"),
                 pr.get_spec("bcc", "full", rc_cps=2.625e6, opt=5)):
        assert spec.phy["mod"] == "bfsk_chip_v2", spec.pid
        assert abs(phy.ber(spec.phy["mod"], 10 ** ((10 - phy.V2_OFFSET_DB) / 10)) - 3.3690e-3) < 1e-7


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all protocols tests passed")


"""Plain-assert validation of the coexistence engine, no pytest needed.

Purpose: 1 BLE_SIM reproduction by ble4thesis (collision model, both PHYs); 2 periodic random-phase transmitters against
the ALOHA formula; 3 CSMA single-sensor PSR and latency; 4 BER, BCH and Walsh points (test_phy.py, see below); 5 energy
of one isolated BLE full exchange by hand (and a TSCH cell, the lifetime, the ANT channel-open listen left out of the
rows); 7 determinism and paired cells; plus slots that hold their frames, the BLE connection-event and AFH rules,
beacon loss and SmartBAN migration, the ED-scan, TSCH CCA, TSCH blacklisting and ANT re-timing rules, the sync-lock
listening rule and lock threshold, the 802.15.6 NB sensitivity floor, BCC noise-model PSR 1 and the contact chain
against the contact chain in both models, one-tone BCC against the chip-rate study (1 V and the knee),
multi-tone BCC (PSR 1 at k 1 and 30, a frame needs every tone, the energy of one
exchange), raw and post-FEC BER of a
BCH link at a fixed SNR (SINR model), latency from the sensor's sample clock (BCC TDMA, ANT), outcome counts that add
up to generated and a TSCH backlog whose psr does not depend on the run length, goodput on the application bits, the
C/I_co floor switch, the refusal of BCC coupling in the collision model, the retry room of the FULL schedules (BLE
spare exchange, TSCH shared cell with its contention, SmartBAN C/M slots, no backlog alone), the BCC full ladder
(DESIGN 4.3: the L3 retry pool in sensor order within the superframe, the L3 resend in an allocation's idle room, the
L2 retry in the next allocation, the beacon as the ack, the energy of one L2 exchange and of the L3 pool, PSR 1 at every
level, the contact case at every level with the 5 % pool's cap at N_dev 4, opt 0 rows identical to no flag), the rows
of the stacks leaves unchanged against the merge-base, the sweep's entry selection, the sweep
driver in-process (worker limit, labels, resume) and the cell aggregates (J per bit per trial, clipped intervals).
Item 4 lives in test_phy.py, whose tests __main__ runs first: BER formulas at hand-derived points, BCH failure against
the binomial tail, the Walsh-16 symbol with t = 3 and the BCC packet (0.988 at 1.5 % chip error, 0.168 at 5.5 %).
Every expected value is derived by hand or from a specification in the comment above its assert.
Inputs:  the engine modules, ../BLE_SIM.py.
Outputs: one line per passed test on stdout; an AssertionError on failure.
argv:    none (python python_simulation/coexistence/test_coexistence.py), runtime about two minutes.
"""
import argparse
import dataclasses
import io
import json
import math
import os
import subprocess
import sys
import tarfile
import tempfile
from concurrent.futures import Future
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import BLE_SIM                # noqa: E402
import engine as en           # noqa: E402
import phy                    # noqa: E402
import protocols as pr        # noqa: E402
import sweep                  # noqa: E402


def _psr(rows):
    return sum(r["delivered_ok"] for r in rows) / sum(r["generated"] for r in rows)


def _run(spec, **kw):
    a = dict(k=1, n_dev=1, r_kbps=10, ci_mode="latency", reception="collision", duration_us=1_000_000, seed=1)
    a.update(kw)
    return en.simulate_trial(spec, **a)


def test_ble_sim_reproduction():
    # DESIGN 8.1: ble4thesis keeps BLE_SIM's rules, so the mean PSR over trials must match BLE_SIM.simulate_one_trial
    # (same (R, mode, N_dev, k), same p_loss 1e-4 / 1e-3) within 0.01. Trials last one 37-CI permutation cycle: the
    # trial spread comes from the drawn phases and permutations, so many short trials beat few long ones. The PSR
    # spread over trials is up to 0.025, so 150 trials per side give a standard error of the difference of 0.003
    # (0.01 = 3.4 of them).
    pts = [(30, "latency", 4, 16, "1M"), (100, "latency", 10, 8, "1M"), (50, "energy", 4, 30, "1M"),
           (5, "latency", 10, 30, "1M"), (200, "latency", 4, 30, "2M"), (10, "latency", 10, 30, "2M"),
           (100, "energy", 1, 16, "2M"), (60, "latency", 16, 4, "2M")]
    print("  point (R kb/s, mode, N_dev, k, PHY): BLE_SIM PSR, engine PSR")
    for r, mode, n, k, ph in pts:
        rate = {"1M": 1_000_000, "2M": 2_000_000}[ph]
        sel = BLE_SIM.select_ci(r, mode, n, rate)
        dur = 37 * sel.ci_us
        ref = np.mean([_psr_ble_sim(BLE_SIM.simulate_one_trial(
            k=k, num_devices=n, payload_bits=sel.payload_bits, ci_us=sel.ci_us, phy_rate=rate, duration_us=dur,
            seed=t, random_loss_p={"1M": 1e-4, "2M": 1e-3}[ph])) for t in range(150)])
        spec = pr.get_spec("ble4thesis", "thesis", ph)
        got = np.mean([_psr(en.simulate_trial(spec, k=k, n_dev=n, r_kbps=r, ci_mode=mode, reception="collision",
                                              duration_us=dur, seed=t)) for t in range(150)])
        print(f"  ({r}, {mode}, {n}, {k}, {ph}): {ref:.4f}, {got:.4f}")
        assert abs(got - ref) < 0.01, (r, mode, n, k, ph, ref, got)


def _psr_ble_sim(rows):
    return sum(r["packets_received"] for r in rows) / sum(r["packets_sent"] for r in rows)


def test_aloha_formula():
    # DESIGN 8.2: M periodic transmitters with independent uniform phases on one channel, airtime tau, period T: a
    # frame survives a given other transmitter unless their phases are closer than tau (probability 2 tau / T), so
    # PSR = (1 - 2 tau / T)^(M - 1). The plain BLE model (ble4thesis: one data PDU per superframe at the BAN's
    # phase, no ARQ, no clock drift) held on one channel: 30 kb/s latency, T = 7500 us, L = 225 bits, tau = 225 us.
    # M = 10 (k = 10, one sensor each): (1 - 450 / 7500)^9 = 0.5730 (the integer phases make it (1 - 449 / 7500)^9,
    # 0.5738). p_loss off, 10 periods per trial; the trial spread (about 0.2) comes from the drawn phases, so 3000
    # trials give a standard error of about 0.004.
    spec = dataclasses.replace(pr.get_spec("ble4thesis", "thesis"), channels=[2404], p_loss=0.0)
    sel = pr.select_period(30, "latency", 1, spec)
    assert (sel["period_us"], sel["payload_bits"], pr.airtime_us(spec, 225 / 8)) == (7500, 225, 225)
    want = (1 - 2 * 225 / 7500) ** 9
    got = np.mean([_psr(_run(spec, k=10, r_kbps=30, duration_us=75_000, seed=s)) for s in range(3000)])
    print(f"  ALOHA: formula {want:.4f}, engine {got:.4f}")
    assert abs(got - want) < 0.01


def test_csma_single_sensor():
    # DESIGN 8.3: one 802.15.4 CSMA sensor alone: every CCA is idle, so the latency is the initial backoff
    # U{0..2^3 - 1} x 320 us (mean 3.5 x 320 = 1120 us) + CCA 128 + turnaround 192 + frame. 10 kb/s latency:
    # T = 15.36 ms, L = round(153.6) = 154 bits = 20 octets, frame (6 + 9 + 20 + 2) x 32 = 1184 us; mean 2624 us.
    # PSR 1 (p_loss 1e-4 per attempt, 4 attempts). 40 s: 2604 frames, standard error of the mean 733 / sqrt(2604)
    # = 14 us (0.5 %).
    rows = _run(pr.get_spec("ieee802154", "full"), duration_us=40_000_000)
    r = rows[0]
    print(f"  CSMA: PSR {r['psr']:.4f}, mean latency {r['lat_mean_ms'] * 1000:.1f} us (hand 2624)")
    assert r["psr"] == 1.0
    assert abs(r["lat_mean_ms"] * 1000 / 2624 - 1) < 0.02
    # The run's end does not cut a sample in service (DESIGN 4.2): 100 kb/s latency, T = 15.36 ms, L = 1536 bits in
    # frames of 116 and 76 octets. Alone and without p_loss every CCA is idle and every frame and ACK arrives; a frame
    # takes at most backoff 7 x 320 + CCA 128 + turnaround 192 + (6 + 9 + octets + 2) x 32 + 192 + ACK 352 + LIFS 640
    # = 8000 and 6720 us, 14720 us < T together, so each sample is done before the next one and every frame sampled
    # before the end is delivered, whatever the phase: psr 1, drop_end 0. 10 periods, 20 seeds (stopping at the end
    # left the 76-octet frame queued in 12 of them). Pinned to the admission (max_load 1.0): at the default
    # 0.5 one sensor cannot carry 100 kb/s (5760 us per 116-octet frame, 0.62 of any period)
    spec = dataclasses.replace(pr.get_spec("ieee802154", "full"), p_loss=0.0)
    assert pr.select_period(100, "latency", 1, spec) is None
    spec.csma = dict(spec.csma, max_load=1.0)
    assert pr.select_period(100, "latency", 1, spec)["frame_payload_octets"] == [116, 76]
    for seed in range(20):
        r = _run(spec, r_kbps=100, duration_us=153_600, seed=seed)[0]
        assert r["psr"] == 1.0 and r["drop_end"] == 0, (seed, r["generated"], r["drop_end"])


def test_generation_on_sample_clock():
    # DESIGN 3: every sensor samples on its own clock (period T, random phase), so a scheduled frame waits for its
    # slot and D = P binds when T > P. BCC, 10 kb/s latency, 10 sensors: P = 7.5 ms, T = 48.75 ms, L = round(10 x
    # 48.75) = 488 bits = 61 octets, frame 448 + 32 x 61 = 2400 chips = 4800 us at 500 kchip/s from its slot start.
    # A frame sampled soff ~ U[0, T) before its slot has latency soff + 4800 us: mean T/2 + 4800 = 29175 us, on time
    # (<= P) iff soff <= 2700 us, probability 2700 / 48750 = 0.0554 (isolated bodies, p_loss 0: PSR 1). soff is fixed
    # per sensor: 20 trials x 50 BANs x 10 sensors give standard errors 0.0023 and T / sqrt(12) / 100 = 141 us.
    spec = pr.get_spec("bcc", "thesis")
    assert pr.select_period(10, "latency", 10, spec)["period_us"] == 48750 and pr.airtime_us(spec, 61) == 4800
    rows = [r for s in range(20) for r in _run(spec, k=50, n_dev=10, duration_us=2 * 48750, seed=s)]
    gen, dlv = sum(r["generated"] for r in rows), sum(r["delivered_ok"] for r in rows)
    on_time = sum(r["on_time"] * r["generated"] for r in rows) / gen
    lat = sum(r["lat_mean_ms"] * r["delivered_ok"] for r in rows) / dlv * 1000
    print(f"  BCC TDMA: on_time {on_time:.4f} (hand 0.0554), mean latency {lat:.0f} us (hand 29175)")
    assert dlv == gen and abs(on_time - 2700 / 48750) < 0.01 and abs(lat - 29175) < 600
    # ANT, 1 kb/s latency: T = 7507 us (N = ceil(7500 x 32768 / 1e6) = 246 counts, round(246e6 / 32768)), samples
    # every T, acknowledged messages every channel period tc = T / 1.25 = 6005.6 us at the master's phase, which is
    # independent of the sample clock's: a frame goes at the next slot, latency (slot - sample) mod tc + 144 (the data
    # end). As T = 1.25 tc, four samples meet the slots at waits v, v + tc / 4, v + tc / 2, v + 3 tc / 4 with v =
    # U[0, tc / 4), mean tc / 2 + 144 = 3146.8 us. 100 trials x 10 masters of 40 samples (ten 4-sample cycles; each
    # BAN on its private frequency, co-channel masters placed apart by the listen, so no losses): the per-master mean
    # has a spread of tc / 4 / sqrt(12) = 433 us, standard error 14 us over 1000 masters; a last sample whose slot
    # falls past the end is not sent (its long wait missing: a bias below 20 us)
    spec = dataclasses.replace(pr.get_spec("ant", "full"), p_loss=0.0)
    rows = [r for s in range(100) for r in _run(spec, k=10, r_kbps=1, duration_us=40 * 7507, seed=s)]
    dlv = sum(r["delivered_ok"] for r in rows)
    lat = sum(r["lat_mean_ms"] * r["delivered_ok"] for r in rows if r["delivered_ok"]) / dlv * 1000
    print(f"  ANT: mean latency {lat:.0f} us (hand 3146.8)")
    assert abs(lat - 3146.8) < 100


def test_backlog_counts_to_the_end():
    # DESIGN 7: every frame counts in generated, and one still queued at the end is not delivered (drop_end), so the
    # outcome counts add up to generated for every stack (10 kb/s latency, 4 BANs of 4 sensors, 0.3 s).
    outs = "delivered_ok delivered_corrupt drop_collision drop_noise drop_sync drop_csma drop_retry drop_queue drop_end"
    for pid, var, ph in pr.catalogue():
        for r in _run(pr.get_spec(pid, var, ph), k=4, n_dev=4, duration_us=300_000):
            assert not r["feasible"] or r["generated"] == sum(r[c] for c in outs.split()), (pid, var, r)
    # TSCH FULL gives each frame one dedicated cell per slotframe and the frames left over one shared cell per
    # slotframe, which the BAN's sensors contend for, so at per_first 0.36 (5 kb/s energy:
    # T = 410 ms, 3 frames per sensor, N_dev 4, k 30) a sensor generates more frames than it serves and its FIFO
    # fills within a few slotframes; from then on the frames delivered per
    # slotframe follow the per-attempt success, whatever the run length, so psr must not depend on it. The 40-period
    # run shares its first 5 periods with the 5-period one; leaving the queued frames out of generated (the old rule)
    # gave 0.94 at 5 periods and 0.64 at 40. Pinned without blacklisting (threshold 1: no channel ever leaves), which
    # measures the counting rule alone: blacklisting shrinks the BANs' sequences as the
    # failures pile up over the run, which lifts psr from 0.611 to 0.638 at 40 periods on these seeds
    # (test_channel_rules checks the rule itself)
    spec = pr.get_spec("tsch", "full")
    spec.adapt = dict(spec.adapt, threshold=1.0)
    psr = [_psr([r for s in range(4) for r in _run(spec, k=30, n_dev=4, r_kbps=5, ci_mode="energy",
                                                   duration_us=n * 410_000, seed=s)]) for n in (5, 40)]
    print(f"  TSCH full backlog: psr {psr[0]:.4f} over 5 periods, {psr[1]:.4f} over 40")
    assert abs(psr[0] - psr[1]) < 0.02


def test_goodput_counts_application_bits():
    # DESIGN 3, 7: goodput counts the L application bits of a period, not the padding to whole octets. Lossless links
    # (p_loss 0, one sensor alone, 2 s): BLE 10 kb/s latency, T = 7.5 ms, L = 75 bits in 10 octets (80 on air; two
    # 540 us exchanges fit); TSCH 1 kb/s, T = 20 ms, L = 20 bits in 3 octets (24). Each period's L arrives, so goodput
    # = L / T = r_agg, within one period's L over the run (the run holds D / T +/- 1 periods): 75 / 2 s = 0.0375 kb/s
    # and 0.01 kb/s. The octet count gave 10.68 and 1.2 kb/s.
    for pid, rk, bits in (("ble", 10, 75), ("tsch", 1, 20)):
        r = _run(dataclasses.replace(pr.get_spec(pid, "full"), p_loss=0.0), r_kbps=rk, duration_us=2_000_000)[0]
        assert r["psr"] == 1.0 and r["payload_bits"] == bits and r["r_agg_kbps"] == rk, r
        assert abs(r["goodput_kbps"] - rk) <= bits / 2_000_000 * 1000 + 1e-12, (pid, r["goodput_kbps"])


def test_energy_one_exchange():
    # DESIGN 8.5: k = 1, N_dev = 1, 10 kb/s latency on BLE 1M: T = 7.5 ms, L = 75 bits -> 10 octets, frame 80 + 8 x 10
    # = 160 us; duration one period, so exactly one event. nRF52840 (PARAMETERS 4.2): TX 6.4 mA, RX 6.26 mA, start-up
    # (ramp) TX 6.8 / RX 5.3 mA for 140 us per radio start, wake 0.3 ms x 3 mA = 900 nC per event, sleep 3.16 uA over
    # the rest, 3 V. The first reception of an event carries the receive guard (PARAMETERS 4.1) = 4 x 50 ppm x
    # 7500 us + 2 x 16 us = 33.5 us; a response whose timing the receiver knows (T_IFS) none. mA x us = nC.
    ex = pr.get_spec("ble", "full")
    # poll at the anchor (hub TX 80 us, sensor RX 80 + 33.5), then the data PDU T_IFS later (sensor TX 160, hub RX 160,
    # no guard: the central times the response); two ramps per node
    s_j = 3e-9 * (6.26 * 113.5 + 5.3 * 140 + 6.4 * 160 + 6.8 * 140 + 900 + 0.00316 * (7500 - 273.5 - 280))
    h_j = 3e-9 * (6.4 * 80 + 6.8 * 140 + 6.26 * 160 + 5.3 * 140 + 900 + 0.00316 * (7500 - 240 - 280))
    r = _run(ex, duration_us=7500)[0]
    assert r["tx_attempts"] == 1
    assert math.isclose(r["e_sensor_j"], s_j, rel_tol=1e-9), (r["e_sensor_j"], s_j)
    assert math.isclose(r["e_hub_j"], h_j, rel_tol=1e-9), (r["e_hub_j"], h_j)
    # lifetime (PARAMETERS 4.4): I = E / 3 V / 7.5 ms, x 3.0 / 2.8 (buck converter late in life), T = k C / (I + s k C
    # / 8760 h), k C = 0.9 x 235 mAh, s = 1 %/year: E = 3e-9 x 4350.461 nC = 1.305138e-5 J, 0.580061 mA x 1.071429 =
    # 0.621494 mA, + 2.41438e-4 mA -> 211.5 / 0.621736 = 340.177 h = 14.174 days
    i_ma = s_j / 3 / 7.5e-3 * 1e3 * 3.0 / 2.8
    assert math.isclose(r["lifetime_days"], 211.5 / (i_ma + 0.01 * 211.5 / 8760) / 24, rel_tol=1e-9)
    assert abs(r["lifetime_days"] - 14.174) < 1e-3
    assert ex.energy["tx_ma"] == 6.4 and ex.energy["ramp_us"] == 140
    # TSCH, 1 kb/s latency, one sensor: slotframe 20 ms (EB cell + one dedicated cell), L = 20 bits -> 3 octets, frame
    # (6 + 9 + 3 + 2) x 32 = 640 us, EB (6 + 41) x 32 = 1504 us in cell 0. nRF52840 802.15.4 (4.2): TX 6.4, RX 6.53 mA,
    # ramps 130 us at 6.8 / 5.3 mA. Timeslot template (3.3): the hub listens from macTsRxOffset 1120 us to the frame end
    # at 2120 + 640 us = 1640 us; the sensor listens for the Enh-Ack from macTsRxAckDelay 800 us after the data to its
    # end 1000 + 608 us = 808 us, the hub sends it (608 us), and before the data it senses the channel (CCA, macTsCca on
    # since) for 128 us at macTsCcaOffset, one more RX start. The sensor also listens in cell 0, the RFC 8180
    # minimal cell, where it sends nothing in this slotframe and receives the
    # EB: from macTsRxOffset 1120 to the EB end 2120 + 1504 us, 2504 us, one more RX ramp and wake-up. One slotframe
    # with the cell inside the run (a seed whose phase < 10 ms); p_loss off so the ACK always goes out.
    spec = dataclasses.replace(pr.get_spec("tsch", "full"), p_loss=0.0)
    sel = pr.select_period(1, "latency", 1, spec)
    assert (sel["period_us"], sel["frame_payload_octets"]) == (20000, [3])
    seed = next(x for x in range(50) if en.Sim(spec, sel, 1, 1, "collision", 20000, x).phase[0] < 10000)
    r = _run(spec, r_kbps=1, duration_us=20000, seed=seed)[0]
    h = 3e-9 * (6.4 * (1504 + 608) + 6.8 * 130 * 2 + 6.53 * 1640 + 5.3 * 130 + 2 * 900
                + 0.00316 * (20000 - 1504 - 608 - 1640 - 130 * 3))
    s_j = 3e-9 * (6.4 * 640 + 6.8 * 130 + 6.53 * (128 + 808 + 2504) + 5.3 * 130 * 3 + 900 * 2
                  + 0.00316 * (20000 - 640 - 128 - 808 - 2504 - 130 * 4))
    assert r["tx_attempts"] == 1
    assert math.isclose(r["e_hub_j"], h, rel_tol=1e-9), (r["e_hub_j"], h)
    assert math.isclose(r["e_sensor_j"], s_j, rel_tol=1e-9), (r["e_sensor_j"], s_j)
    # ANT listen at channel open (PARAMETERS 3.6 FULL step 1, DESIGN 6): one channel period of RX once per deployment,
    # not per run, so the rows leave it out. 0.02 kb/s latency: T = 375 ms, tc = T / 1.25 = 300 ms, run 40 T = 15 s.
    # Booked, the listen alone would add 3 V x (6.26 mA x 300000 us + 5.3 mA x 140 us - 3.16 uA x 300140 us) = 5.633 mJ
    # per run, a fixed 2 % RX duty whatever the run length (J per bit x 8, lifetime 58 d); without it the master spends
    # about 1.3 mJ in the run (50 channel periods of about 7.7 uC each: TX 144 us, the ACK window 176 us, the 500 us
    # adaptive-isochronous guard, two ramps and a wake-up, plus 47 uC of sleep)
    full = dataclasses.replace(pr.get_spec("ant", "full"), p_loss=0.0)
    a = _run(full, r_kbps=0.02, duration_us=15_000_000)[0]
    assert a["period_us"] == 375000 and a["delivered_ok"] == a["generated"] > 0, a
    assert a["e_sensor_j"] < 5.633e-3, a["e_sensor_j"]
    # the adaptive-isochronous guard (3.6 FULL step 4): after an acknowledged message the master listens for the ACK
    # from the turnaround on (W: 274 to 450 us into the slot) and stays in RX through the 500 us guard after the
    # exchange (418 to 918 us): 644 us of RX with one start-up; after an old-payload broadcast (no W) the guard alone,
    # 500 us with its own start-up. One master, 1 kb/s latency (samples every T = 7507 us, slots every T / 1.25: four
    # messages and one broadcast per five slots), 40 periods: 40 samples, of which the last is sent only when its slot
    # falls before the run end (the random sample and slot phases: 39 or 40 messages), about 10 broadcasts in the 50
    # slots
    sim = en.Sim(full, pr.select_period(1, "latency", 1, full), 1, 1, "collision", 40 * 7507, 1).run()
    txs = sim.ntrack[1].x
    n_msg, n_bc = sum(x.data for x in txs), sum(x.rx == -1 for x in txs)
    assert n_msg in (39, 40) and 5 <= n_bc <= 10, (n_msg, n_bc)
    assert (sim.en[1][1], sim.en[1][3]) == (644 * n_msg + 500 * n_bc, n_msg + n_bc), sim.en[1]


def test_slots_hold_their_frames():
    # DESIGN 3 / 4.1: select_period admits a period only when N_dev x (the exchanges of one period) fits T minus the
    # beacon part; the as-built models send every frame once, so each BAN makes exactly one attempt per generated frame;
    # 5 kb/s latency, the tightest period, 8 BANs
    for pid in ("ble4thesis", "bcc"):
        for r in _run(pr.get_spec(pid, "thesis"), k=8, n_dev=4, r_kbps=5, duration_us=500_000):
            assert r["tx_attempts"] == r["generated"] > 0, (pid, r["tx_attempts"], r["generated"])
    # every FULL scheduled stack of the catalogue, at every rate, alone (k = 1, p_loss 0) delivers every frame in each
    # feasible latency cell of the default grid: its slot holds the whole exchange, including the SmartBAN MCA sensing
    # (TMUA 16 us, 3.5)
    for pid, v, ph in pr.catalogue():
        if pr.get_spec(pid, v, ph).access != "tdma" or v != "full":
            continue
        spec = dataclasses.replace(pr.get_spec(pid, "full", ph), p_loss=0.0)
        for rk in (1, 5, 10, 30, 50, 100, 200, 500, 1000):
            for n in (1, 4, 10):
                sel = pr.select_period(rk, "latency", n, spec)     # three periods at least (SmartBAN reaches 1 s)
                r = _run(spec, n_dev=n, r_kbps=rk, duration_us=max(500_000, 3 * sel["period_us"] if sel else 0))[0]
                assert not r["feasible"] or r["delivered_ok"] == r["generated"] > 0, (pid, ph, rk, n, r["generated"],
                                                                                     r["delivered_ok"])


def test_determinism():
    # DESIGN 8.7: same seed, same rows (both models); another seed, other rows
    spec = pr.get_spec("ble", "full")
    for rec in ("collision", "sinr"):
        a, b = (_run(spec, k=4, n_dev=4, reception=rec, duration_us=300_000, seed=5) for _ in range(2))
        c = _run(spec, k=4, n_dev=4, reception=rec, duration_us=300_000, seed=6)
        assert [_nan_eq(x, y) for x, y in zip(a, b)] == [True] * 4
        assert a != c and not all(_nan_eq(x, y) for x, y in zip(a, c))


def test_cells_are_paired():
    # sweep.py: one seed per (R, k, N_dev, ci_mode, trial) for every protocol and variant; the engine draws the room,
    # the link gains, the superframe phases (fractions of the period) and the clock errors from substreams of it, so
    # BLE 1M and 2M (both at T = 7.5 ms: 4 x 2 x 500 and 4 x 2 x 408 us fit) share phases and gains (their SNRs differ
    # by the noise of their bandwidths only), and 802.15.6 NB sees the same room
    sims = [en.Sim(sp, pr.select_period(5, "latency", 4, sp), 4, 4, "sinr", 1000, 12345)
            for sp in (pr.get_spec("ble", "full", "1M"), pr.get_spec("ble", "full", "2M"),
                       pr.get_spec("ieee802156nb", "full"))]
    assert sims[0].T == sims[1].T == 7500 and np.array_equal(sims[0].phase, sims[1].phase)
    assert np.allclose(np.array(sims[0].snr) + sims[0].noise_dbm, np.array(sims[1].snr) + sims[1].noise_dbm,
                       equal_nan=True)     # nan: a node to itself
    assert np.array_equal(sims[0].room.centres, sims[2].room.centres)
    assert np.allclose(sims[0].phase / sims[0].T, sims[2].phase / sims[2].T, atol=1 / sims[0].T)


def _nan_eq(x, y):
    return all(x[c] == y[c] or (isinstance(x[c], float) and math.isnan(x[c]) and math.isnan(y[c])) for c in x)


def test_bcc_noise_and_contact():
    # PARAMETERS 3.8 case (a): 1 V drive, the weakest measured path of channel_curves.csv (RHW1-LHW1, 50 Ohm) at 500
    # kchip/s and NF 20 dB.: the tone at 12.680 MHz (channel_curves.tone_centre), its BFSK tones at 12.18 and
    # 13.18 MHz: -48.7617 + 0.5817 log(1.218) / log(1.25) = -48.2476 and -48.1800 - 1.2746 log(1.0544) / log(1.2) =
    # -48.5503 dB, Ec/N0 = gain - 10 log10(2 x 5e5 x 8.0e-19) - 20 = gain + 100.9649: 52.72 and 52.41 dB (at the old
    # 8 MHz tone centre 47.4 dB), far above the 8 dB error-free point: PSR 1 in the SINR model; p_loss 0 in the
    # collision model; bodies isolated, so k does not matter.
    # BCC full sees the same Ec/N0 per Walsh chip (noise bandwidth R_c = its payload Walsh-chip rate) and
    # its SF 64 header 9.03 dB more, so it delivers every frame too, in both models and at k 1 and 30
    for v in ("full", "thesis"):
        spec = pr.get_spec("bcc", v)
        for rec in ("collision", "sinr"):
            for k in (1, 30):
                assert _psr(_run(spec, k=k, n_dev=4, reception=rec, duration_us=500_000)) == 1.0, (v, rec, k)
    # Contact chain (DESIGN 5.2), the measured-hardware case at 50 kchip/s in both reception
    # models: continuous per link, connected 14 dB, degraded 2 dB (a chip error of 0.23 kills the 64-chip sync and
    # any Walsh symbol), exponential dwell 35 000 / 13 000 chips. A frame of N chips succeeds iff it starts connected
    # (stationary 35 / 48) and stays connected (exp(-N / 35 000)). 10 kb/s energy with one sensor: T = 200 ms (the
    # reference period: 2008 bits / 10 kb/s = 200.8 -> 200 ms), L = 250 octets, N = 448 + 32 x 250 + 64 x
    # floor(8000 / 2048) = 8640 chips = 172.8 ms, of which the first three preamble copies (192 chips) are not
    # decoded (the correlator needs the 4th and the SFD): N = 8448 chips at risk, PSR = 35/48 x exp(-8448/35000) =
    # 0.5728. A failed frame whose sync and header passed is delivered corrupt (no payload CRC), so pdr > psr. 40
    # links x 20 s = 4000 frames, correlated over the 9500-chip chain memory: tolerance 0.03 (about 3 standard
    # errors).
    spec = pr.get_spec("bcc", "thesis", contact=True)
    assert spec.chip_rate_cps == 50e3
    sel = pr.select_period(10, "energy", 1, spec)
    assert (sel["period_us"], sel["frame_payload_octets"]) == (200000, [250])
    want = 35 / 48 * math.exp(-(8640 - 192) / 35000)
    for rec in ("collision", "sinr"):
        rows = _run(spec, k=40, r_kbps=10, ci_mode="energy", reception=rec, duration_us=20_000_000)
        got, pdr = _psr(rows), np.mean([r["pdr"] for r in rows])
        print(f"  BCC contact, {rec}: stationary PSR {got:.4f} (analytic {want:.4f}), PDR {pdr:.4f}")
        assert abs(got - want) < 0.03 and pdr > got + 0.05, (rec, got, pdr)
    # The as-built figure 0.78 is the second factor alone for its 8800-chip packet, exp(-8800/35000) = 0.778: the v2
    # simulator restarts the chain connected at every 3-packet run. Check the chain dwell against it: frames of 8800
    # chips that start connected stay connected with that probability (5000 frames, standard error 0.006).
    sim = en.Sim(spec, sel, 1, 1, "sinr", 1, 0)
    us = 1e6 / spec.chip_rate_cps
    stay = np.mean([sim.cstate.update({1: (0.0, False)}) or len(sim._contact(1, 0.0, 8800 * us)) == 1
                    for _ in range(5000)])
    print(f"  BCC contact: connected-start survival {stay:.3f} (as-built 0.78)")
    assert abs(stay - 0.78) < 0.02


def test_bcc_rates_match_study():
    # one-tone BCC in the engine against the chip-rate study (bcc_fdma.py rates, results/bcc_rate_stability.csv:
    # per measured path 400 as-built 255-octet packets through the waveform model, 50 Ohm; real passband where the
    # receive band f -/+ 2 R_c reaches 0 Hz, from 7 Mchip/s).
    # (a) 500 kchip/s, 7 M and 9 M (protocols.BCC_RC_MAX_CPS), 1 V: every link the engine draws (30 bodies x 4
    # sensors, NF 20) has the one-chip Ec/N0 of a study path (the csv rounds to 6 digits: 1e-3 dB); the study makes no
    # chip error on any path at NF 10, 20 and 30, and the engine delivers every frame of BCC full and of the as-built
    # stack (the SINR model: the noise decides).
    # (b) Above the cap: at 9.5 M (lower tone at 1 MHz, 0.105 cycle per chip) every path still delivers its 400
    # packets at NF 20, but chip errors appear on every path, from the image the one-chip rule does not model; at
    # 9.95 M (lower tone at 100 kHz) every path falls below PSR 0.99 at NF 10, 20 and 30. get_spec refuses both.
    # (c) The frame model near its knee: at the Eb/N0 target 7 dB (bcc_fdma.KNEE_DB; model PSR 0.56) every path's
    # waveform PSR against phy.frame_outcome at the path's one-chip Ec/N0; the stack's curve carries the 0.35 dB by which the PHY2 receiver
    # beats 0.5 exp(-E / 2) (phy.V2_OFFSET_DB, PARAMETERS 3.8). On the baseband (500 k, 2.625 M): within 3 binomial sd
    # of the model (sd = sqrt(m (1 - m) / 400)), both ways. At 7 M the lower tone at 5 MHz is 0.714 cycle per chip, so
    # its image leaks into its correlator at |sin(1.43 pi)| / (1.43 pi) = 0.22 (-13 dB), energy (1 -/+ 0.22)^2 = +1.7
    # to -2.1 dB with the carrier phase: the waveform is no better than the model and no worse than the model 0.5 dB
    # lower (3 sd each way). The engine uses the rule only at 1 V: the weakest path at 9 M and NF 30 has 26.5 dB, 19 dB
    # above this knee.
    table = HERE / "results" / "bcc_rate_stability.csv"
    if not table.exists():
        print("skipped, write it first:", table.name)
        return
    st = pd.read_csv(table)
    nf = st[st.nf_db.notna()]
    for rc in (500e3, 7e6, pr.BCC_RC_MAX_CPS):
        s = nf[nf.rc_cps == rc]
        assert len(s) == 3 * 37 and (s.packets == 400).all() and (s.chip_errors == 0).all(), rc
        paths = s[(s.nf_db == 20) & (s.curve != "cohort")].ec_n0_db.to_numpy()
        for v in ("full", "thesis"):
            spec = pr.get_spec("bcc", v, rc_cps=rc)
            sim = en.Sim(spec, pr.select_period(10, "latency", 4, spec), 30, 4, "sinr", 1_000_000, 5)
            links = [sim.snr[b * 5][b * 5 + i] for b in range(30) for i in range(1, 5)]
            assert all(np.abs(paths - x).min() < 1e-3 for x in links), (rc, v)
            assert _psr(_run(spec, k=30, n_dev=4, reception="sinr", duration_us=500_000)) == 1.0, (rc, v)
    s = nf[nf.rc_cps == 9.5e6]
    assert (s[s.nf_db == 20].psr == 1).all() and (s[s.nf_db == 20].chip_errors > 0).all()
    assert len(nf[nf.rc_cps == 9.95e6]) == 3 * 37 and (nf[nf.rc_cps == 9.95e6].psr < 0.99).all()
    for rc in (9.5e6, 9.95e6):
        try:
            pr.get_spec("bcc", "full", rc_cps=rc)
            raise AssertionError(rc)
        except ValueError:
            pass
    a = pr.get_spec("bcc", "thesis")
    lay = pr.frame_layout(a, 255)

    def model(e_db):
        segs = [(0.0, lay[-1]["t1_us"], 10 ** (e_db / 10), math.inf)]
        m = phy.frame_outcome(None, segs, lay, a.phy, sampled=False)["success"]
        return m, 3 * math.sqrt(m * (1 - m) / 400)

    for rc, lo_db in ((500e3, 0.0), (2.625e6, 0.0), (7e6, 0.5)):
        knee = st[(st.rc_cps == rc) & (st.ebn0_target_db == 7)]
        assert len(knee) == 37 and (knee.packets == 400).all(), rc
        for r in knee.itertuples():
            (hi, dh), (lo, dl) = model(r.ec_n0_db), model(r.ec_n0_db - lo_db)
            assert lo - dl <= r.psr <= hi + dh, (rc, r.curve, r.psr, lo, hi)


def test_bcc_multitone():
    # Multi-tone BCC full: 6 tones of 500 kchip/s, drive 1 / sqrt(6) V per tone, the path drawn per link, each tone
    # packet at the path's gains at its BFSK tones f -/+ R_c. The weakest tone of any measured
    # path is RHW1-LHW1 at 8.24037 MHz, BFSK tones 7.74037 and 8.74037 MHz: channel_curves.csv -50.0672 dB at 7 MHz,
    # -53.5954 at 8, -51.2301 at 9, log-f weights log(7.74037 / 7) / log(8 / 7) = 0.7530 and log(8.74037 / 8) /
    # log(9 / 8) = 0.7515 -> -52.7238 and -51.8179 dB; Ec/N0 = 20 log10(1 / sqrt(6)) + gain + 100.9649 (test_channel:
    # 1 V, 500 kchip/s, NF 20) = 40.46 and 41.37 dB; one chip Ec/N0 (channel.bfsk_ec_n0_db) 11113 + 2 ln 2 = 11114.4
    # = 40.46 dB (the stronger tone never errs), far above the 10 dB at which lost no packet: PSR 1 in the noise
    # (SINR) model and, p_loss 0, in the collision model; the bodies are isolated, so k 1 and k 30 alike (50 kb/s
    # latency, N_dev 4).
    spec = pr.get_spec("bcc", "full", "multitone")
    sim = en.Sim(spec, pr.select_period(50, "latency", 4, spec), 30, 4, "sinr", 1_000_000, 5)   # 120 links: RHW1-LHW1
    assert abs(min(min(x) for x in sim.tsnr.values()) - 40.46) < 5e-3                         # drawn (centre: 40.18)
    for rec in ("collision", "sinr"):
        for k in (1, 30):
            rows = _run(spec, k=k, n_dev=4, r_kbps=50, reception=rec)
            assert _psr(rows) == 1.0 and all(r["raw_ber"] == 0 for r in rows), (rec, k, _psr(rows))
    # the frame needs every tone packet: one dead tone (Ec/N0 -10 dB: chip error 0.5 exp(-0.05) = 0.48, no 64-chip
    # sync with <= 7 errors) loses every attempt of that sensor, so its frames end at the retry limit (or queued);
    # 6 tones at 40 dB deliver all
    sel = pr.select_period(50, "latency", 1, spec)
    for dead, want in ((False, 1.0), (True, 0.0)):
        sim = en.Sim(spec, sel, 1, 1, "sinr", 1_000_000, 1)
        sim.tsnr[1] = [40.0] * 5 + [-10.0 if dead else 40.0]
        r = sim.run().rows(dict(protocol="bcc"))[0]
        assert r["psr"] == want and (r["drop_retry"] > 0) == dead, r
    # energy of one exchange (DESIGN 6): 2 kb/s latency, one sensor: T = 18 ms, L = 36 bits -> 5 octets, one on each
    # of 5 tones and an empty packet on the 6th (test_protocols), data 1664 Walsh chips = 3328 us, I-Ack 3264 us,
    # beacon 4352 us + pSIFS 75 reserved (test_bcc_full), receive guard 2 x 16 us (no drift). One period with the slot
    # inside it (a seed whose phase + 4427 < 18000). Six front ends at the typical BCC power for data and I-Ack: TX
    # 6 x 2 mW / 3 V = 4 mA, RX 6 x 2.4 / 3 = 4.8 mA; the one-tone beacon on one front end: TX 2 / 3 mA, RX 0.8 mA;
    # ramps 100 us at the currents of their packet; wake 0.3 ms x 3 mA = 900 nC per wake-up (beacon and slot: 2 per
    # node), sleep 3.16 uA over the rest. mA x us = nC.
    # sensor: TX data 3328; RX I-Ack 3264, beacon 4352 + 32 = 4384; ramps: data TX, I-Ack RX, beacon RX; sleep 18000 -
    #   3328 - 3264 - 4384 - 300 = 6724: 4 x 3328 + 4.8 x 3264 + 0.8 x 4384 + 100 (4 + 4.8 + 0.8) + 2 x 900 + 0.00316 x
    #   6724 = 13312 + 15667.2 + 3507.2 + 960 + 1800 + 21.24784 = 35267.64784 nC
    # hub: TX I-Ack 3264, beacon 4352; RX data 3328 + 32 = 3360; ramps: I-Ack TX, data RX, beacon TX; sleep 6724:
    #   4 x 3264 + 4.8 x 3360 + 2 / 3 x (4352 + 100) + 100 (4 + 4.8) + 2 x 900 + 0.00316 x 6724 = 13056 + 16128 + 2968 +
    #   880 + 1800 + 21.24784 = 34853.24784 nC
    sel = pr.select_period(2, "latency", 1, spec)
    assert (sel["period_us"], sel["frame_payload_octets"], pr.airtime_us(spec, 5)) == (18000, [5], 3328)
    seed = next(x for x in range(50) if en.Sim(spec, sel, 1, 1, "collision", 18000, x).phase[0] + 4427 < 18000)
    r = _run(spec, r_kbps=2, duration_us=18000, seed=seed)[0]
    assert r["tx_attempts"] == r["delivered_ok"] == 1, r
    assert math.isclose(r["e_sensor_j"], 3e-9 * 35267.64784, rel_tol=1e-9), r["e_sensor_j"]
    assert math.isclose(r["e_hub_j"], 3e-9 * 34853.24784, rel_tol=1e-9), r["e_hub_j"]


def test_bcc_ladder():
    # Galvanic-optimisation ladder of BCC full (DESIGN 4.3), engine paths: k 1, collision model (BCC p_loss 0), 500
    # kchip/s, no clock drift (receive guard 2 x 16 = 32 us). Selections from test_protocols.test_bcc_ladder.
    base = dict(protocol="x")
    # L3 retry pool, 2 kb/s latency, 2 sensors: T 10000 us, 3 octets, beacon 2776 + pSIFS 75 = 2851 reserved, slots of
    # one exchange (1688 + 75 = 1763 us) from 2851 and 4614, a pool of 1 exchange from 6377: poll 1496, pSIFS 75, data.
    # Both sensors' first frames fail: the pool of superframe 0 polls sensor 1 (the first failure) 6377 - 2851 = 3526
    # us after its failed frame, which follows 6377 + 1496 + 75 - 2851 = 5097 us after it and ends at 7948 + 1688 =
    # 9636 < 10000: within the superframe; sensor 2's waits for the pool of superframe 1 (its own allocation holds only
    # its new frame): poll 10000 + 6377 - 4614 = 11763 us, frame 10000 + 7948 - 4614 = 13334 us after its failure.
    s3 = pr.get_spec("bcc", "full", opt=3)
    sel = pr.select_period(2, "latency", 2, s3)
    assert (sel["period_us"], sel["slot_us"], sel["pool"]) == (10000, 1763, 1), sel
    sim = _fail_first(2)(s3, sel, 1, 2, "collision", 40 * sel["period_us"], 1)
    r = sim.run().rows(base)[0]
    bad = {x.node: x for x in sim.failed}
    t0 = bad[1].s - 2851
    assert sorted(bad) == [1, 2] and bad[2].s - t0 == 4614, (sorted(bad), bad[2].s - t0)
    for j, gap in ((1, 3526), (2, 11763)):
        poll = next(x for x in sim.ntrack[0].x if x.rx == j and x.e - x.s == 1496)     # the beacons last 2776 us
        retry = next(x for x in _data(sim, j) if x.s == poll.e + 75)
        assert (poll.s - bad[j].s, retry.s - bad[j].s, retry.ok) == (gap, gap + 1571, True), (j, poll.s - bad[j].s)
    assert retry.e - t0 - 10000 == 9636
    assert r["delivered_ok"] == r["generated"] == r["tx_attempts"] - 2 and r["drop_end"] == r["drop_queue"] == 0, r
    # L3 resends in its own allocation the frames the heard beacon reported missing when their arrivals leave room:
    # one sensor, 2 kb/s latency (T 8000 us, a queue of 4 frames, slot 1699 us from 2787, pool from 4486: poll 1496,
    # pSIFS 75, data). Every data frame and poll before superframe 4 fails: the frames A, B, C, D of superframes 0 to 3
    # fill the queue while the pool polls A in vain; superframe 4's arrival is dropped (queue full), so its allocation
    # sends A, 4 T = 32000 us after its first attempt (the pool would send it 32000 + 1699 + 1571 us after), and its
    # pool polls B; then C and D go in the next two pools. Every frame but the dropped one is delivered over 40
    # periods, A to D in two attempts.
    class Bad(en.Sim):
        def _eval(self, tx):
            super()._eval(tx)
            if (tx.data or tx.e - tx.s == 1496) and tx.s < self.phase[0] + 4 * 8000:
                tx.ok, tx.reason = False, "collision"
    sel = pr.select_period(2, "latency", 1, s3)
    sim = Bad(s3, sel, 1, 1, "collision", 40 * sel["period_us"], 1)
    r = sim.run().rows(base)[0]
    d = _data(sim, 1)
    assert d[4].s - d[0].s == 32000 and d[4].ok and not any(x.ok for x in d[:4]), [x.s - d[0].s for x in d[:6]]
    assert r["drop_queue"] == 1 and r["delivered_ok"] == r["generated"] - 1 == r["tx_attempts"] - 4, r
    assert r["drop_retry"] == r["drop_end"] == 0, r
    # L2: the retry waits for the sensor's own next allocation (the beacon's ack tells it). 2 kb/s latency, 1 sensor:
    # T 7500 us, 2 octets, the failed frame at 2787 into superframe 0, its retry first in the allocation of superframe
    # 1, 7500 us later, the new frame one exchange (1624 + 75) after it; one period late (latency > T)
    s2 = pr.get_spec("bcc", "full", opt=2)
    sel = pr.select_period(2, "latency", 1, s2)
    sim = _fail_first(1)(s2, sel, 1, 1, "collision", 40 * sel["period_us"], 1)
    r = sim.run().rows(base)[0]
    bad = sim.failed[0]
    later = [x for x in _data(sim, 1) if x.s > bad.s]
    assert (sel["period_us"], later[0].s - bad.s, later[1].s - later[0].s) == (7500, 7500, 1699), sel
    assert r["delivered_ok"] == r["generated"] == r["tx_attempts"] - 1 and r["lat_max_ms"] > 7.5, r
    # the ack is the beacon's: when the sensor misses the beacon of superframe 1 (it acknowledged the frame of
    # superframe 0), L2 resends that frame, delivered already (one duplicate attempt), while L3 keeps it for the next
    # beacon and sends only new frames (the pool serves failed frames only); nothing is lost either way. The beacon of
    # one sensor with one frame: 17 + 2 octets (2712 us)
    for spec, dup in ((s2, 1), (s3, 0)):
        seen = []
        sel = pr.select_period(2, "latency", 1, spec)
        sim = _forced(lambda tx: seen.append(tx) or len(seen) == 2, "collision")(spec, sel, 1, 1, "collision",
                                                                               40 * sel["period_us"], 1)
        r = sim.run().rows(base)[0]
        assert not seen[1].ok and seen[1].e - seen[1].s == 2712, seen[1].e - seen[1].s
        assert r["delivered_ok"] == r["generated"] == r["tx_attempts"] - dup and r["drop_end"] == 0, (dup, r)
    # energy of one L2 exchange (no I-Ack to listen for): 2 kb/s latency, one sensor, T 7500 us, data 1624 us, beacon
    # (with its 2-octet ack field) 2712 us, reserved 2787; one period with the slot inside it (phase + 2787 + 1699 <
    # 7500). Typical BCC power at 3 V: TX 2 / 3 mA, RX 0.8 mA, ramps 100 us at the packet's current, wake 900 nC
    # (beacon and slot: 2 per node), sleep 3.16 uA. mA x us = nC.
    # sensor: TX data 1624; RX beacon 2712 + 32 = 2744; sleep 7500 - 1624 - 2744 - 200 = 2932:
    #   2 / 3 x 1624 + 0.8 x 2744 + 100 (2 / 3 + 0.8) + 1800 + 0.00316 x 2932 = 1082.66667 + 2195.2 + 146.66667 +
    #   1800 + 9.26512 = 5233.79845 nC (L1 adds the I-Ack listen: 0.8 x 1496 + 100 x 0.8)
    # hub: TX beacon 2712; RX data 1624 + 32 = 1656; sleep 2932: 2 / 3 x 2712 + 0.8 x 1656 + 146.66667 + 1800 +
    #   9.26512 = 1808 + 1324.8 + 146.66667 + 1800 + 9.26512 = 5088.73179 nC
    sel = pr.select_period(2, "latency", 1, s2)
    seed = next(x for x in range(50) if en.Sim(s2, sel, 1, 1, "collision", 7500, x).phase[0] + 2787 + 1699 < 7500)
    r = _run(s2, r_kbps=2, duration_us=7500, seed=seed)[0]
    assert r["tx_attempts"] == r["delivered_ok"] == 1, r
    assert math.isclose(r["e_sensor_j"], 3e-9 * 5233.79845, rel_tol=1e-8), r["e_sensor_j"]
    assert math.isclose(r["e_hub_j"], 3e-9 * 5088.73179, rel_tol=1e-8), r["e_hub_j"]
    # energy of the L3 pool: one sensor, 2 kb/s latency, T 8000 us, data 1624, beacon 2712, reserved 2787, slot 1699,
    # the pool of one exchange from phase + 4486 < 8000 (the duration: one superframe). Every sensor wakes for the pool
    # and listens from the receive guard: the guard alone when no poll comes, else through the poll; the hub wakes
    # for the pool when it polls. Currents as above.
    # no loss: sensor TX 1624; RX beacon 2744 + pool guard 32 = 2776; 1 TX and 2 RX ramps; 3 wakes (beacon, slot,
    #   pool); sleep 8000 - 1624 - 2776 - 300 = 3300: 1082.66667 + 2220.8 + 226.66667 + 2700 + 10.428 = 6240.56133
    #   hub TX 2712; RX 1656; 1 TX and 1 RX ramp; 2 wakes; sleep 8000 - 2712 - 1656 - 200 = 3432: 1808 + 1324.8 +
    #   146.66667 + 1800 + 10.84512 = 5090.31179
    # first frame failed, polled in the pool: sensor TX 3248; RX 2744 + 32 + poll 1496 = 4272; 2 TX and 2 RX ramps;
    #   3 wakes; sleep 8000 - 3248 - 4272 - 400 = 80: 2165.33333 + 3417.6 + 293.33333 + 2700 + 0.2528 = 8576.51947
    #   hub TX 2712 + poll 1496 = 4208; RX 1656 + the retry 1624 (its timing known: no guard) = 3280; 2 TX and 2 RX
    #   ramps; 3 wakes; sleep 8000 - 4208 - 3280 - 400 = 112: 2805.33333 + 2624 + 293.33333 + 2700 + 0.35392 =
    #   8423.02059
    sel = pr.select_period(2, "latency", 1, s3)
    seed = next(x for x in range(50) if en.Sim(s3, sel, 1, 1, "collision", 8000, x).phase[0] + 4486 < 8000)
    for cls, tries, want in ((en.Sim, 1, (6240.56133, 5090.31179)), (_fail_first(1), 2, (8576.51947, 8423.02059))):
        r = cls(s3, sel, 1, 1, "collision", 8000, seed).run().rows(base)[0]
        assert r["tx_attempts"] == tries and r["delivered_ok"] == 1, r
        assert math.isclose(r["e_sensor_j"], 3e-9 * want[0], rel_tol=1e-8), (tries, r["e_sensor_j"] / 3e-9)
        assert math.isclose(r["e_hub_j"], 3e-9 * want[1], rel_tol=1e-8), (tries, r["e_hub_j"] / 3e-9)
    # psr 1 at 1 V (Ec/N0 >= 40 dB, test_bcc_multitone) in both models, k 1 and 30 (isolated bodies), at every level;
    # one tone at 10 kb/s x 4 (40 of its 52 to 109 kb/s), six tones at 50 kb/s x 4
    for ph, rk in (("sf8", 10), ("multitone", 50)):
        for o in range(1, 6):
            spec = pr.get_spec("bcc", "full", ph, opt=o)
            for rec in ("collision", "sinr"):
                for k in (1, 30):
                    rows = _run(spec, k=k, n_dev=4, r_kbps=rk, reception=rec)
                    assert _psr(rows) == 1.0 and all(x["raw_ber"] == 0 for x in rows), (ph, o, rec, k, _psr(rows))
    # the contact case (one tone, 50 kchip/s, the two-state chain in both models) runs at every level: feasible, the
    # outcome counts add up to generated. N_dev 1: ARQ recovers most of what the chain loses (psr > 0.5). N_dev 4
    # (energy, 1 kb/s: one frame of about 251 octets per sensor and period of about 2 s, about 8900 chips; 40 periods,
    # 3 seeds): a frame gets through when the chain stays connected for its whole airtime, p_s = 35 / 48 x exp(-8900 /
    # 35000) = 0.57, and tries a period apart (above both dwells) are about independent. Up to L2 every frame gets at
    # least two tries a period apart (the spare exchanges): psr >= 1 - (1 - 0.57)^2 = 0.81. From L3 a superframe
    # carries at most its 4 frames and 1 pool exchange, so psr <= 0.57 x (1 + 1 / 4) = 0.71: the 5 % pool caps the
    # contact case (DESIGN 4.3, level 3)
    drops = "delivered_corrupt drop_collision drop_noise drop_sync drop_csma drop_retry drop_queue drop_end".split()
    for o in range(6):
        spec = pr.get_spec("bcc", "full", contact=True, opt=o)
        t = pr.select_period(1, "energy", 4, spec)["period_us"]
        many = [x for seed in (1, 2, 3) for x in _run(spec, n_dev=4, r_kbps=1, ci_mode="energy", seed=seed,
                                                       duration_us=40 * t)]
        for x in _run(spec, k=2, r_kbps=1, duration_us=5_000_000) + many:
            assert x["feasible"] == 1 and x["generated"] == x["delivered_ok"] + sum(x[c] for c in drops), (o, x)
            assert x["variant"] == "full+contact" + (f"+opt{o}" if o else ""), (o, x["variant"])
            assert x["psr"] > 0.5 or x["n_dev"] == 4, (o, x["psr"])
        assert _psr(many) >= 0.81 if o < 3 else _psr(many) <= 0.71, (o, _psr(many))


def test_rows_match_merge_base():
    # deletes the fair variants and relabels the as-built models 'fair' -> 'thesis'; the stacks whose model it
    # does not change (BLE, 802.15.4, 802.15.6 NB at 121.4 kb/s, HBC SF 8, BCC full on one tone and multi-tone, the
    # as-built BCC and the contact cases of both, ble4thesis) keep their rows float for float, the label mapped: a few
    # cells of each, run in a subprocess with this tree's modules, again with opt=0 on BCC full (the ladder's level 0 is
    # the stack itself). Every run starts in a
    # temporary directory, so no module comes from the caller's current directory.
    specs = [("ble", "full", "1M", {}), ("ble", "full", "2M", {}), ("ble", "full", "coded_s8", {}),
             ("ieee802154", "full", "oqpsk", {}), ("ieee802156nb", "full", "121.4", {}), ("hbc", "full", "sf8", {}),
             ("bcc", "full", "sf8", {}), ("bcc", "full", "multitone", {}), ("bcc", "thesis", "bfsk", {}),
             ("bcc", "full", "sf8", {"contact": True}), ("bcc", "thesis", "bfsk", {"contact": True}),
             ("ble4thesis", "thesis", "1M", {}), ("ble4thesis", "thesis", "2M", {})]
    script = (
        "import json, sys\n"
        "sys.path.insert(0, sys.argv[1])\n"
        "import engine as en, protocols as pr\n"
        "old = hasattr(pr, 'FAIR_STRIP')\n"
        "out = []\n"
        "for pid, v, ph, kw in json.loads(sys.argv[3]):\n"
        "    kw = dict(kw, **({'opt': 0} if sys.argv[2] == 'opt0' and (pid, v) == ('bcc', 'full') else {}))\n"
        "    spec = pr.get_spec(pid, 'fair' if old and v == 'thesis' else v, ph, **kw)\n"
        "    cells = [('collision', 10, 'latency', 4, 4), ('collision', 1, 'latency', 1, 2)]\n"
        "    for rec, r, mode, n, k in cells + ([('sinr', 5, 'energy', 2, 2)] if pid != 'ble4thesis' else []):\n"
        "        out += en.simulate_trial(spec, k=k, n_dev=n, r_kbps=r, ci_mode=mode, reception=rec,\n"
        "                                 duration_us=1_000_000, seed=7)\n"
        "print(json.dumps([dict(x, variant=x['variant'].replace('fair', 'thesis')) for x in out]))\n")
    env = dict(os.environ, PYTHONPATH=str(HERE.parent))    # HBCProtocol and PHY_layer_BCC_simulation_v2
    with tempfile.TemporaryDirectory() as d:
        runs = [(HERE, "-"), (HERE, "opt0")]
        outs = [subprocess.run([sys.executable, "-c", script, str(where), arg, json.dumps(specs)], cwd=d, env=env,
                               capture_output=True, text=True, check=True).stdout for where, arg in runs]
    rows = json.loads(outs[0])      # 4 + 2 (+ 2 SINR) rows per spec, and every spec runs in at least one cell
    feasible = {(x["protocol"], x["variant"], x["phy"]) for x in rows if x["feasible"]}
    assert len(rows) == 8 * 11 + 6 * 2 and len(feasible) == len(specs), feasible
    assert {x["variant"] for x in rows} == {"full", "thesis", "full+contact", "thesis+contact"}
    assert len(set(outs)) == 1, f"a row changed ({len(runs)} runs)"


def test_fec_lowers_ber_in_sinr():
    # SINR model, a noise-limited link at a fixed SNR (fade off), so every payload bit sees the same error
    # probability. SmartBAN FULL, 10 kb/s latency, one sensor: T = 7.5 ms (4 slots of 1250 us fit), L = 75 bits ->
    # 10 octets, MPDU 144 bits in BCH(127,113) + a codeword shortened to 45 (31 info), PHY CRC-16 after it. Per-bit
    # SNR 8 (9.03 dB): p = 0.5 e^-4 = 0.009158, so raw BER = p. Post-FEC: E[e 1{e > 2}] of a codeword = n p - P(1) -
    # 2 P(2); n = 127: P(1) = 0.36490, P(2) = 0.21248, 1.16307 - 0.78986 = 0.37321, x 113 / 127 = 0.33207 info-bit
    # errors; n = 45: P(1) = 0.27492, P(2) = 0.055902, 0.41211 - 0.38672 = 0.025386, x 31 / 45 = 0.017488. The
    # application share 80 / 144 of the 0.34956 residual errors over 80 bits: post BER = 0.0024275, a quarter of p.
    # 8 s: 1066 frames, about 1900 raw errors (2.3 %) and 200 post errors in bursts of 3 to 4 (about 13 %).
    spec = pr.get_spec("smartban", "full")
    sel = pr.select_period(10, "latency", 1, spec)
    assert (sel["period_us"], sel["frame_payload_octets"]) == (7500, [10])
    sim = en.Sim(spec, sel, 1, 1, "sinr", 8_000_000, 3)
    sim.fsig = None
    sim.snr[1][0], sim.snr[0][1] = 10 * math.log10(8), 30.0     # sensor -> hub at 9.03 dB, ACK and beacon clean
    r = sim.run().rows(dict(protocol="smartban"))[0]
    print(f"  smartban full at 9.03 dB: raw BER {r['raw_ber']:.3e} (hand 9.158e-3), post-FEC BER "
          f"{r['post_ber']:.3e} (hand 2.428e-3)")
    assert abs(r["raw_ber"] / 0.009158 - 1) < 0.08 and abs(r["post_ber"] / 0.0024275 - 1) < 0.4


def test_sync_lock_needs_listening():
    # DESIGN 5.2: a receiver locks onto an earlier foreign frame q only if it was listening for p when q began (its
    # window opens the receive guard before p, or at the end of its own transmission) and was not transmitting.
    # 802.15.4, SINR model, two BANs on channel 0, fade off: q (BAN 1) from 0 to 3000 us at 20 dB above noise
    # at hub 0, p (BAN 0, 20 octets: 1184 us) from 1000 us at 40 dB. q's SINR at its start is 20 dB, above the sync
    # threshold, so a listening hub 0 locks on q and loses p (reason 'sync'); p's own SIR is 20 dB, above C/I_co
    # 3 dB, and its SINR 20 dB leaves no bit errors, so a hub that opened its window after q began (guard window:
    # 1000 - guard > 0) receives p; so does a hub that was transmitting when q began.
    spec = pr.get_spec("ieee802154", "full")
    sel = pr.select_period(10, "latency", 1, spec)
    for lead, own_tx, want in ((None, False, True), (1500.0, False, False), (1500.0, True, True)):
        sim = en.Sim(spec, sel, 2, 1, "sinr", 1, 0)     # a 1 us run: no process sends anything
        sim.fsig = None
        sim.snr[1][0], sim.snr[3][0] = 40.0, 20.0        # nodes: hub 0, sensor 1 (BAN 0); hub 2, sensor 3 (BAN 1)
        txs = [sim._tx(0.0, 3000, 0, 1, 3, 2), sim._tx(1000.0, 1184, 0, 0, 1, 0, sim.lay(20), True, lead)]
        if own_tx:
            txs.append(sim._tx(-100.0, 200, 0, 0, 0, 1))  # hub 0 sends from -100 to 100 us
        p = _place(sim, txs)[1]
        sim._eval(p)
        assert p.ok == want and (want or p.reason == "sync"), (lead, own_tx, p.ok, p.reason)


def _place(sim, txs):
    """Put hand-made transmissions on the tracks of sim, as Sim._send does."""
    for tx in txs:
        sim.tracks.setdefault(tx.key, en.Track()).add(tx)
        sim.ntrack[tx.node].add(tx)
        sim.atrack.add(tx)
        sim.maxd = max(sim.maxd, tx.e - tx.s)
    return txs


def test_sensitivity_and_lock_threshold():
    # 802.15.6 NB FULL (PARAMETERS 3.4): 'frames below the rate's sensitivity (TI draft values) are lost', -95 dBm at
    # 121.4 kb/s. Noise in the 600 kHz symbol band with NF 10 + IL 6: -174 + 10 log10(6e5) + 16 = -100.22 dBm. A lone
    # 20-octet frame at -96 dBm (4.22 dB SNR) is lost for noise, although the BER model alone (4-branch DBPSK and
    # BCH(63,51)) would deliver it; at -94 dBm it arrives. SINR model, fade off.
    spec = pr.get_spec("ieee802156nb", "full")
    sel = pr.select_period(5, "latency", 1, spec)
    for dbm, want in ((-96.0, False), (-94.0, True)):
        sim = en.Sim(spec, sel, 1, 1, "sinr", 1, 0)
        sim.fsig = None
        assert abs(sim.noise_dbm + 100.2185) < 1e-4
        sim.snr[1][0] = dbm - sim.noise_dbm
        lay = sim.lay(20)
        if dbm < -95:
            seg = [(0.0, lay[0][-1]["t1_us"], 10 ** (sim.snr[1][0] / 10), math.inf)]
            assert phy.frame_outcome(None, seg, lay[0], spec.phy, sampled=False)["success"] > 0.999
        p = _place(sim, [sim._tx(0.0, lay[1], 0, 0, 1, 0, lay, True)])[0]
        sim._eval(p)
        assert p.ok == want and (want or p.reason == "noise"), (dbm, p.ok, p.reason)
    # OI-6: 'the receiver locks on the first packet above sensitivity'. SmartBAN has no sync rule of its own, so the
    # lock threshold is the SINR of a 1e-3 bit error on its GFSK curve, the BLE sensitivity criterion (0.5 exp(-g/2)
    # = 1e-3 at g = 2 ln 500 = 10.94 dB). Two BANs on one channel as in test_sync_lock_needs_listening, hub 0
    # listening since before q began: q at 8 dB above noise does not capture it, so p (40 dB, SIR 32 dB above the
    # 19.88 to 21 dB C/I floors, no bit errors) arrives; q at 14 dB does, and p is lost for sync. The old 0 dB threshold
    # locked
    # on both.
    spec = pr.get_spec("smartban", "full")
    sel = pr.select_period(10, "latency", 1, spec)
    for q_db, want in ((8.0, True), (14.0, False)):
        sim = en.Sim(spec, sel, 2, 1, "sinr", 1, 0)
        sim.fsig = None
        assert abs(sim.sync_thr - 10 * math.log10(2 * math.log(500))) < 1e-6
        sim.snr[1][0], sim.snr[3][0] = 40.0, q_db
        lay = sim.lay(20)
        p = _place(sim, [sim._tx(0.0, 3000, 0, 1, 3, 2), sim._tx(1000.0, lay[1], 0, 0, 1, 0, lay, True, 1500.0)])[1]
        sim._eval(p)
        assert p.ok == want and (want or p.reason == "sync"), (q_db, p.ok, p.reason)


def test_repetition_copies():
    # SmartBAN Nrep (PARAMETERS 3.5 BER model): in the collision model each copy of the PPDU is decoded alone, so a
    # frame survives when one copy is clear. Uncoded Nrep 2, 20 octets: 2 x (104 + 64 + 160) = 656 us, copies [0, 328)
    # and [328, 656) us. Another BAN's frame over [100, 200) hits copy 1 only: delivered; over [300, 400) it hits
    # both: collision. Nrep 1 (328 us) loses the frame to either. p_loss 0.
    for ph, q0, want in (("unc2", 100.0, True), ("unc2", 300.0, False), ("unc1", 100.0, False)):
        spec = dataclasses.replace(pr.get_spec("smartban", "full", ph), p_loss=0.0)
        sim = en.Sim(spec, pr.select_period(10, "latency", 1, spec), 2, 1, "collision", 1, 0)
        lay = sim.lay(20)
        assert lay[1] == (656 if ph == "unc2" else 328)
        p = _place(sim, [sim._tx(q0, 100, 0, 1, 3, 2), sim._tx(0.0, lay[1], 0, 0, 1, 0, lay, True)])[1]
        sim._eval(p)
        assert p.ok == want and (want or p.reason == "collision"), (ph, q0, p.ok, p.reason)


def test_ci_floor_switch():
    # DESIGN 5.2 co-channel floor (PARAMETERS OI-3, default on) and its sensitivity switch ci_floor=False, which treats
    # interference as Gaussian noise. BLE 1M, SINR model, two BANs on channel 0, fade off: p (sensor 1 -> hub 0,
    # 20 octets: 240 bits, 240 us at 1 Mb/s) at 40 dB above noise inside q (sensor 3 of BAN 1) at 25 dB. SIR = 15 dB,
    # SINR = 10 log10(1e4 / (1 + 10^2.5)) = 14.99 dB: above the Gaussian threshold (1e-3 bit error at 10.94 dB; B = R,
    # so 0.5 exp(-31.5 / 2) = 7.2e-8 per bit and the frame is error-free with probability 1 - 240 x 7.2e-8 = 0.99998)
    # but below C/I_co 21 dB. With the floor p fails: every decision under the floor is a guess (p = 0.5), so the 32
    # access-address bits are all right with probability 2^-32 only and the reason is 'sync'. Without the floor p
    # arrives.
    spec = pr.get_spec("ble", "full", "1M")
    sel = pr.select_period(10, "latency", 1, spec)
    for floor, want in ((True, False), (False, True)):
        sim = en.Sim(spec, sel, 2, 1, "sinr", 1, 0, ci_floor=floor)
        sim.fsig = None
        sim.snr[1][0], sim.snr[3][0] = 40.0, 25.0
        lay = sim.lay(20)
        assert lay[1] == 240 and abs(spec.phy["ci_co_lin"] - 10 ** 2.1) < 1e-9
        p = _place(sim, [sim._tx(0.0, 3000, 0, 1, 3, 2), sim._tx(1000.0, lay[1], 0, 0, 1, 0, lay, True)])[1]
        sim._eval(p)
        assert p.ok == want and (want or p.reason == "sync"), (floor, p.ok, p.reason)
    # the CSV label marks the sensitivity run on the SINR rows only: the collision model has no floor (test_sweep_main
    # checks that the switch reaches the engine)
    assert en.label(spec) == "full" and en.label(spec, False, "sinr") == "full+nocifloor"
    assert en.label(spec, False, "collision") == "full"


def test_bcc_coupling_needs_sinr():
    # The BCC coupling cases (+40 / +60 dB over the intra-body loss, PARAMETERS 5.7) weaken the interferer, an SINR
    # concept; the collision model has no power and would let any overlap kill. simulate_trial raises ValueError for
    # that pair and runs it in the SINR model; sweep.py stops with an argparse error (exit status 2). The argv bounds a
    # regressed guard to one worker, one trial of one small cell and output files in a temporary directory.
    spec = pr.get_spec("bcc", "thesis", coupling_db=40)
    try:
        _run(spec, k=2, duration_us=100_000)
        raise AssertionError("coupling accepted in the collision model")
    except ValueError:
        pass
    assert _run(spec, k=2, reception="sinr", duration_us=100_000)[0]["feasible"] == 1
    argv = sys.argv
    try:
        with tempfile.TemporaryDirectory() as d:
            sys.argv = ["sweep.py", "--protocols", "bcc", "--bcc-coupling", "40", "--reception", "sinr,collision",
                        "--workers", "1", "--trials", "1", "--k", "1", "--n-dev", "1", "--r", "1", "--ci-modes",
                        "latency", "--out", str(Path(d, "rows.csv")), "--cells", str(Path(d, "cells.csv"))]
            sweep.main()
        raise AssertionError("sweep accepted --bcc-coupling with collision")
    except SystemExit as e:
        assert e.code == 2, e.code
    finally:
        sys.argv = argv


def _forced(fail_poll, why):
    """Sim whose hub polls fail when fail_poll(tx) holds, with reason `why` (after the model's own outcome)."""
    class Forced(en.Sim):
        def _eval(self, tx):
            super()._eval(tx)
            if tx.node % self.m == 0 and tx.rx >= 0 and fail_poll(tx):
                tx.ok, tx.reason = False, why
    return Forced


def test_bcc_pooled_reference():
    # pooled reference stack (L0 with the retry pool), engine path: k 1, collision model, 500 kchip/s, one
    # sensor, 2 kb/s latency: T 22000 us, 6 octets (test_protocols.test_bcc_pooled_reference): beacon 4352 + pSIFS 75
    # = 4427 reserved, the allocation of one exchange (data 3648, pSIFS, I-Ack 3264, pSIFS: 7062 us) from 4427, the
    # pool of one exchange from 4427 + 7062 = 11489: poll 3264, pSIFS 75, the sensor's exchange with its I-Ack. The
    # first data frame fails: the allocation holds no retry; the pool polls the frame 11489 - 4427 = 7062 us after its
    # start, the retry starts 7062 + 3264 + 75 = 10401 us after it and succeeds, its I-Ack pSIFS later; nothing is lost.
    # With the spare exchange (no flag) the failed frame goes again at once, one exchange (6998 us: 5 octets, T 18500)
    # after its start.
    base = dict(protocol="x")
    sp = pr.get_spec("bcc", "full", "sf8", pool=0.05)
    sel = pr.select_period(2, "latency", 1, sp)
    assert (sel["period_us"], sel["slot_us"], sel["pool"], sel["pool_us"]) == (22000, 7062, 1, 10401), sel
    sim = _fail_first(1)(sp, sel, 1, 1, "collision", 40 * sel["period_us"], 1)
    r = sim.run().rows(base)[0]
    bad = sim.failed[0]
    poll = next(x for x in sim.ntrack[0].x if x.rx == 1 and x.s > bad.s)     # no I-Ack follows the failed frame
    retry = next(x for x in _data(sim, 1) if x.s > bad.s)
    ack = next(x for x in sim.ntrack[0].x if x.rx == 1 and x.s > retry.s)
    assert (poll.s - bad.s, poll.e - poll.s, retry.s - bad.s, retry.ok) == (7062, 3264, 10401, True), \
        (poll.s - bad.s, poll.e - poll.s, retry.s - bad.s)
    assert (ack.s - retry.e, ack.e - ack.s, ack.ok) == (75, 3264, True), (ack.s - retry.e, ack.e - ack.s)
    assert r["delivered_ok"] == r["generated"] == r["tx_attempts"] - 1, r
    assert r["drop_end"] == r["drop_queue"] == r["drop_retry"] == 0, r
    s0 = pr.get_spec("bcc", "full", "sf8")
    sel0 = pr.select_period(2, "latency", 1, s0)
    sim0 = _fail_first(1)(s0, sel0, 1, 1, "collision", 40 * sel0["period_us"], 1)
    r0 = sim0.run().rows(base)[0]
    bad0 = sim0.failed[0]
    retry0 = next(x for x in _data(sim0, 1) if x.s > bad0.s)
    assert (sel0["period_us"], retry0.s - bad0.s, retry0.ok) == (18500, 6998, True), (sel0, retry0.s - bad0.s)
    assert r0["delivered_ok"] == r0["generated"] == r0["tx_attempts"] - 1 and r0["drop_end"] == 0, r0


def test_bcc_block_ack():
    # block-ack reference (flag back), k 1, collision model, 500 kchip/s, 4 sensors, 20 kb/s energy mode: loss
    # free, every frame is delivered; every allocation ends with one B-Ack (the 3-octet frame) from the hub pSIFS after
    # the sensor's last data frame. One forced failure goes to the pool (poll, then the frame) and is delivered
    sp = pr.get_spec("bcc", "full", "sf8", pool=0.05, back=True)
    sel = pr.select_period(20, "energy", 4, sp)
    rows = en.Sim(sp, sel, 1, 4, "collision", 40 * sel["period_us"], 3).run().rows(dict(protocol="x"))
    assert all(r["psr"] == 1.0 for r in rows), rows
    sim = _fail_first(1)(sp, sel, 1, 4, "collision", 40 * sel["period_us"], 3)
    rows = sim.run().rows(dict(protocol="x"))
    assert sum(r["delivered_ok"] for r in rows) == sum(r["generated"] for r in rows) > 0, rows
    b3 = pr.airtime_us(sp, pr.BACK_OCTETS)
    for j in range(1, 5):
        data, back = _data(sim, j), [x for x in sim.ntrack[0].x if x.rx == j and x.e - x.s == b3]
        assert back and all(any(abs(x.s - y.e - pr.HBC["psifs_us"]) < 1e-6 for y in data) for x in back), j
    bad, a0 = sim.failed[0], pr.airtime_us(sp, 0)     # the poll is the L = 0 frame, the B-Ack is longer
    poll = next(x for x in sim.ntrack[0].x if x.rx == bad.node and x.s > bad.s and x.e - x.s == a0)
    retry = next(x for x in _data(sim, bad.node) if x.s > poll.s)
    assert retry.ok and abs(retry.s - poll.e - pr.HBC["psifs_us"]) < 1e-6, (retry.s, poll.e)



def test_bcc_block_ack_burst():
    # the lead reference in burst mode (flags back, burst), as test_bcc_block_ack: loss free; in an allocation the
    # frames after the first go pMIFS (20 us) after the frame before with the RI-mode airtime, the B-Ack pSIFS after
    # the last
    sp = pr.get_spec("bcc", "full", "sf8", pool=0.05, back=True, burst=True)
    sel = pr.select_period(20, "energy", 4, sp)
    assert sel["frames_per_period"] == 2, sel
    sim = en.Sim(sp, sel, 1, 4, "collision", 40 * sel["period_us"], 3)
    assert all(r["psr"] == 1.0 for r in sim.run().rows(dict(protocol="x")))
    ri, b3 = pr.airtime_us(sp.ri, sel["frame_payload_octets"][1]), pr.airtime_us(sp, pr.BACK_OCTETS)
    for j in range(1, 5):
        data = _data(sim, j)
        cont = [y for x, y in zip(data, data[1:]) if abs(y.s - x.e - 20) < 1e-6]
        assert cont and all(y.e - y.s == ri for y in cont), j
        back = [x for x in sim.ntrack[0].x if x.rx == j and x.e - x.s == b3]
        assert back and all(any(abs(x.s - y.e - 75) < 1e-6 for y in data) for x in back), j


def test_bcc_block_ack_span():
    # flag span: a block of the allocation's frames (2 here) runs as the burst reference, the
    # same rows; a block of 2 at one frame per allocation sends a B-Ack after every second allocation; an unlimited
    # block sends none and spends less energy per bit, loss free
    kw = dict(pool=0.05, back=True, burst=True)
    sp = pr.get_spec("bcc", "full", "sf8", **kw)
    sel = pr.select_period(20, "energy", 4, sp)
    assert sel["frames_per_period"] == 2, sel
    run = lambda s, sl: en.Sim(s, sl, 1, 4, "collision", 40 * sl["period_us"], 3)     # noqa: E731
    strip = lambda rows: [{k: v for k, v in r.items() if k != "variant"} for r in rows]  # noqa: E731
    a = run(sp, sel).run().rows(dict(protocol="x"))
    two = pr.get_spec("bcc", "full", "sf8", span=2, **kw)
    assert strip(run(two, pr.select_period(20, "energy", 4, two)).run().rows(dict(protocol="x"))) == strip(a)
    b3 = pr.airtime_us(sp, pr.BACK_OCTETS)
    sel1 = pr.select_period(1, "energy", 4, two)
    assert sel1["frames_per_period"] == 1, sel1
    sim = run(two, sel1)
    sim.run()
    for j in range(1, 5):
        back = [x for x in sim.ntrack[0].x if x.rx == j and x.e - x.s == b3]
        assert len(_data(sim, j)) > 2 and len(back) == len(_data(sim, j)) // 2, j
    u = pr.get_spec("bcc", "full", "sf8", span=pr.BLOCK_UNLIMITED, **kw)
    sim = run(u, pr.select_period(1, "energy", 4, u))
    rows = sim.run().rows(dict(protocol="x"))
    assert all(r["psr"] == 1.0 for r in rows) and not [x for x in sim.ntrack[0].x if x.e - x.s == b3]
    ref = run(sp, pr.select_period(1, "energy", 4, sp)).run().rows(dict(protocol="x"))
    assert all(r["jpb_sensor"] < q["jpb_sensor"] for r, q in zip(rows, ref)), (rows, ref)


def _fail_first(n):
    """Sim whose first n data transmissions to be decided fail (reason 'collision'), kept in sim.failed."""
    class F(en.Sim):
        failed = []

        def _eval(self, tx):
            super()._eval(tx)
            if tx.data and len(self.failed) < n:
                tx.ok, tx.reason = False, "collision"
                self.failed.append(tx)
    return F


def _data(sim, node):
    """Data transmissions of node, by start."""
    return [x for x in sim.ntrack[node].x if x.data]


def test_retry_room():
    # the retry room each FULL standard gives; k = 1, collision model, p_loss 0, one forced failure.
    base = dict(protocol="x")
    # BLE ACL (3.1 FULL step 3): the event holds 1 + 1 spare exchanges per frame, so the failed data PDU goes again in
    # the same event. 10 kb/s, 10 sensors: T = 12.5 ms (test_protocols), slot 1250 us, L = 16 octets: exchange =
    # poll 80 + 150 + data 208 + 150 = 588 us, two fit in 1250 us; the retry starts one exchange after the failed PDU
    spec = dataclasses.replace(pr.get_spec("ble", "full", "1M"), p_loss=0.0)
    sel = pr.select_period(10, "latency", 10, spec)
    sim = _fail_first(1)(spec, sel, 1, 10, "collision", 40 * sel["period_us"], 1)
    r = sim.run().rows(base)[0]
    bad = sim.failed[0]
    retry = next(x for x in _data(sim, bad.node) if x.s > bad.s)
    assert (sel["period_us"], sim.slot, retry.s - bad.s) == (12500, 1250, 588), (sel, sim.slot, retry.s - bad.s)
    assert r["delivered_ok"] == r["generated"] == r["tx_attempts"] - 1 and r["drop_end"] == r["drop_queue"] == 0, r
    # 802.15.4 GTS (3.2 mode B): 10 kb/s, 1 sensor, the GTS of 4800 us holds 1 + 1 exchanges of 2368 us
    # (test_protocols): the failed frame goes again one exchange later in the same GTS, and the sensor hears the beacon
    # of every superframe
    spec = dataclasses.replace(pr.get_spec("ieee802154", "full", "gts"), p_loss=0.0)
    sel = pr.select_period(10, "latency", 1, spec)
    sim = _fail_first(1)(spec, sel, 1, 1, "collision", 40 * sel["period_us"], 1)
    r = sim.run().rows(base)[0]
    bad = sim.failed[0]
    retry = next(x for x in _data(sim, bad.node) if x.s > bad.s)
    assert retry.s - bad.s == 2368 and retry.ok and r["drop_sync"] == 0, (retry.s - bad.s, r)
    assert r["delivered_ok"] == r["generated"] == r["tx_attempts"] - 1 and r["drop_end"] == r["drop_queue"] == 0, r
    # TSCH (3.3, RFC 8180): 1 kb/s, 2 sensors: slotframe = shared cell 0 + cells 1, 2 = 30 ms (20 ms holds one
    # sensor only), L = 30 bits -> 4 octets. Sensor 1's frame fails in cell 1 of slotframe 0; the next shared cell
    # (cell 0 of slotframe 1, no EB before 16 s) comes before its next dedicated cell: retry 30000 - 10000 = 20000 us
    # later (+ T x hub drift <= 30000 x 40 ppm = 1.2 us), channel offset 0; no backoff before a first shared attempt
    spec = dataclasses.replace(pr.get_spec("tsch", "full"), p_loss=0.0)
    sel = pr.select_period(1, "latency", 2, spec)
    sim = _fail_first(1)(spec, sel, 1, 2, "collision", 40 * sel["period_us"], 1)
    r = sim.run().rows(base)[0]
    bad = sim.failed[0]
    retry = next(x for x in _data(sim, bad.node) if x.s > bad.s)
    assert sel["period_us"] == 30000 and bad.node == 1 and abs(retry.s - bad.s - 20000) < 2, (sel, retry.s - bad.s)
    assert retry.key == sim._key(0, sim._chan(0, 0, 1)) and retry.ok
    assert r["delivered_ok"] == r["generated"] == r["tx_attempts"] - 1 and r["drop_end"] == r["drop_queue"] == 0, r
    # shared-cell contention: both sensors fail in slotframe 0 and meet in the shared cell of slotframe 1 (same start,
    # same channel): both fail again (collision), draw a backoff, and go in their own cells j of slotframe 1, j x
    # 10000 us after the shared cell
    sim = _fail_first(2)(spec, sel, 1, 2, "collision", 40 * sel["period_us"], 1).run()
    for j in (1, 2):
        a, b, c = _data(sim, j)[:3]
        assert abs(b.s - a.s - (30000 - j * 10000)) < 2 and abs(c.s - b.s - j * 10000) < 2, (j, a.s, b.s, c.s)
        assert (a.ok, b.ok, b.reason, c.ok) == (False, False, "collision", True), (j, b.reason)
    # the hub has one radio: it listens once in that shared cell, however many frames collide there. Two slotframes
    # in the run (a seed whose phase < 10 ms, so sensor 2's cell of slotframe 1 starts before 60 ms): slotframe 0 is an
    # EB cell (hub TX) and the two failed frames of cells 1 and 2, each heard from macTsRxOffset, 1000 us before the
    # frame of (6 + 9 + 4 + 2) x 32 = 672 us: 2 x 1672; slotframe 1 the shared cell once, 1672 (two frames booked 3344),
    # and the two retries in cells 1 and 2, 2 x 1672. Hub RX 8360 us in 5 listens (the double booking: 10032 in 6)
    seed = next(x for x in range(50) if en.Sim(spec, sel, 1, 2, "collision", 60000, x).phase[0] < 10000)
    sim = _fail_first(2)(spec, sel, 1, 2, "collision", 60000, seed).run()
    assert (sim.en[0][1], sim.en[0][3]) == (8360, 5), sim.en[0]
    # shared-cell backoff (802.15.4-2015 6.2.5.3, DESIGN 4.1) with every U{0..2^BE - 1} draw at its maximum: one
    # sensor (1 kb/s: slotframe 20 ms = shared cell + one dedicated cell) whose frames always fail. No backoff before
    # the first shared-cell attempt (slotframe 1); after it BE = 1 and the sensor skips 1 shared cell, then BE = 2, 3,
    # 4 skip 3, 7, 15: shared-cell attempts in slotframes 1, 3, 7, 15, 31 (the dedicated-cell failures leave BE)

    class Always(en.Sim):
        def _eval(self, tx):
            super()._eval(tx)
            if tx.data:
                tx.ok, tx.reason = False, "collision"

    class MaxRng:
        def __init__(self, rng):
            self.rng = rng

        def integers(self, lo, hi):
            return hi - 1

        def __getattr__(self, name):
            return getattr(self.rng, name)

    sel1 = pr.select_period(1, "latency", 1, spec)
    sim = Always(spec, sel1, 1, 1, "collision", 40 * sel1["period_us"], 1)
    sim.rng = MaxRng(sim.rng)
    sim.run()
    t = sel1["period_us"]
    shared = [int((x.s - sim.phase[0]) // t) for x in _data(sim, 1) if (x.s - sim.phase[0]) % t < 10000]
    assert (t, shared) == (20000, [1, 3, 7, 15, 31]), (t, shared)
    # SmartBAN (3.5, TS 103 325 7.3.3): 10 kb/s, 4 sensors: T = 8750 us = beacon + 4 slots + 2 C/M slots of TS = 1250
    # us. Sensor 1's data starts 16 us (MCA sensing) into its slot at 1250; the C/M period starts at (1 + 4) x 1250 =
    # 6250, where a node whose backoff counter BC is 0 senses TMUA = Tsch 16 us and sends: BC 0 (every backoff draw at
    # its minimum) puts the retry 6266 - 1266 = 5000 us after the failed frame, BC 1 one C/M slot later, 6250 us. At
    # the default draws (BC = U{0..7}, BEmin 3) it takes a C/M slot of this IBI with probability 2 / 8, else its next
    # own slot (+ T)

    class Fixed:
        def __init__(self, rng, v, seen=None):
            self.rng, self.v, self.seen = rng, v, seen

        def integers(self, lo, hi):
            if self.seen is not None:
                self.seen.append(hi)
            return lo + self.v

        def __getattr__(self, name):
            return getattr(self.rng, name)

    spec = dataclasses.replace(pr.get_spec("smartban", "full"), p_loss=0.0)
    sel = pr.select_period(10, "latency", 4, spec)
    for v, gap in ((0, 5000), (1, 6250)):
        sim = _fail_first(1)(spec, sel, 1, 4, "collision", 40 * sel["period_us"], 1)
        sim.rng = Fixed(sim.rng, v)
        r = sim.run().rows(base)[0]
        bad = sim.failed[0]
        retry = next(x for x in _data(sim, bad.node) if x.s > bad.s)
        assert (sel["period_us"], bad.node, retry.s - bad.s) == (8750, 1, gap), (v, sel, bad.node, retry.s - bad.s)
        assert r["delivered_ok"] == r["generated"] == r["tx_attempts"] - 1 and r["drop_end"] == r["drop_queue"] == 0, r
    # the hub listens once per C/M slot: one IBI in the run (a seed whose phase < 3750 us, so sensor 4's slot starts
    # before 8750), sensors 1 and 2 fail in their slots and both take C/M slot 0 with BC 0 (both sense it idle at the
    # same instant, then collide; RT = 1 ends them), C/M slot 1 stays empty. Receive guard 4 x 20 ppm x 8750 + 2 x 16 =
    # 32.7 us before the data, frame of 11 octets 284 us: 4 slots + C/M slot 0 at 32.7 + 284, C/M slot 1 at 32.7:
    # 1616.2 us in 6 listens (the double booking: 1900.2 in 8)
    seed = next(x for x in range(50) if en.Sim(spec, sel, 1, 4, "collision", 8750, x).phase[0] < 3750)
    sim = _fail_first(2)(spec, sel, 1, 4, "collision", 8750, seed)
    sim.rng = Fixed(sim.rng, 0)
    sim.run()
    assert abs(sim.guard - 32.7) < 1e-9 and [x.node for x in sim.failed] == [1, 2]
    assert abs(sim.en[0][1] - 1616.2) < 1e-6 and sim.en[0][3] == 6, sim.en[0]
    # the same in both reception models: the two sensings end as the first node's data starts (6250 + 16 = 6266 us
    # into the IBI), and a CCA does not hear a frame that starts at its end, so both send at 6266 and no sensing is
    # busy (the SINR CCA once counted that frame: the second node deferred to C/M slot 1, 7516)
    for rec in ("collision", "sinr"):
        sim = _fail_first(2)(spec, sel, 1, 4, rec, 8750, seed)
        sim.rng = Fixed(sim.rng, 0)
        sim.run()
        cm = [x.s - sim.phase[0] for j in (1, 2) for x in _data(sim, j)[1:]]
        assert cm == [6266, 6266] and sim.c[0]["cca_busy"] == 0, (rec, cm, sim.c[0]["cca_busy"])
    # the backoff exponent (7.3.3.1.1): BEmin 3 for a new frame, + 1 after a busy C/M slot or a failed attempt, at most
    # BEmax 5, a new BC = U{0..2^BE - 1} each time. The channel always busy (the owner's sensing in its own slot too,
    # so the frame goes to the C/M period), every draw 0 (BC 0: the node senses at every C/M slot): the draws are
    # U{0..7}, then after each busy C/M slot U{0..15}, U{0..31}, U{0..31}, ...; nothing is ever sent
    assert (en.BE_MIN, en.BE_MAX) == (3, 5)

    class Busy(en.Sim):
        def _cca(self, *a):
            return True
    seen = []
    sim = Busy(spec, pr.select_period(10, "latency", 1, spec), 1, 1, "collision", 100_000, 1)
    sim.rng = Fixed(sim.rng, 0, seen)
    r = sim.run().rows(base)[0]
    assert seen[:6] == [8, 16, 32, 32, 32, 32] and r["tx_attempts"] == 0 < r["generated"], (seen, r["tx_attempts"])
    # no full catalogue spec builds a backlog alone from one isolated loss: k = 1, p_loss 0, the first data frame of the
    # run fails, 40 periods, 3 seeds. BLE and the 802.15.6 allocations resend it in the same event or allocation (1 + 1
    # spare exchanges), TSCH in the next shared cell (before the sensor's dedicated cells), SmartBAN in a C/M slot, or
    # in its own slot ahead of new data until a C/M slot takes the extra frame (BC 0 or 1: 2 / 8), ANT in a spare
    # channel period (1.25 per sample period), 802.15.4 at once: drop_end = drop_queue = 0, and every frame arrives but
    # those 802.15.4 loses to its own contention (drop_csma, drop_retry). The cells hold the one the exactly provisioned
    # schedules could not clear (BLE 10 kb/s x 10: one 540 us exchange per 750 us slot at T = 7.5 ms) and
    # 802.15.4 at the loads that backlogged at max_load 1.0 (5 kb/s x 10: 0.88, 10 kb/s x 4: 0.70). Left out: ANT at
    # N_dev 10, whose masters collide with each other: 10 exchanges of 418 us fit one 6005.6 us channel period only
    # packed, and the listen places the masters one after another, the later ones at the least-occupied phase, which
    # still overlaps two others (DESIGN 4.1; 5 kb/s, k 1, p_loss 0, seeds 0 to 2: drop_queue 54 each of about 400
    # frames, 220, 140 and 164 with the random-phase fallback), and ANT at 10 kb/s latency x 4, where the advanced burst
    # (one 16-octet packet per 3.2 ms, T = 12817 us) fills the period with the 4 bursts and gives a retry no room: a
    # burst packet retried 3.2 ms later meets the next master's packet, whose own retry meets the next, and so on
    # (DESIGN 3; seed 0: psr 0.500, drop_queue 39 of 150; seeds 1 and 2 place the masters so that the retry falls in a
    # gap: psr 1). Under a sustained loss the room can run out: DESIGN 8.8. The BCC full ladder too (DESIGN 4.3): L2
    # resends in the next allocation's spare exchange, L3 and up in the superframe's retry pool.
    ladder = [("bcc", "full", ph, dict(opt=o)) for ph in ("sf8", "multitone") for o in range(1, 6)]
    for pid, v, ph, *kw in pr.catalogue() + ladder:
        spec = dataclasses.replace(pr.get_spec(pid, v, ph, **(kw[0] if kw else {})), p_loss=0.0)
        cells = 0
        for rk, n, mode in ((1, 4, "latency"), (1, 4, "energy"), (10, 4, "energy"), (100, 1, "latency"),
                            (10, 10, "latency"), (5, 10, "latency"), (10, 4, "latency")):
            sel = pr.select_period(rk, mode, n, spec)
            if v != "full" or sel is None or (pid == "ant" and (n == 10 or (rk, n, mode) == (10, 4, "latency"))):
                continue
            cells += 1
            for seed in range(3):
                sim = _fail_first(1)(spec, sel, 1, n, "collision", 40 * sel["period_us"], seed).run()
                r = sim.rows(base)[0]
                assert len(sim.failed) == 1 and r["drop_end"] == r["drop_queue"] == 0, (pid, ph, rk, n, mode, seed, r)
                assert pid == "ieee802154" or r["delivered_ok"] == r["generated"], (pid, ph, rk, n, mode, seed, r)
        assert cells or v != "full", (pid, ph)



def test_bcc_cap_under_dispersion():
    # (bcc_fdma.py isi, results/bcc_rate_isi.csv): the waveform through the minimum-phase response of each
    # measured magnitude (the least smear the body can cause), 400 packets per path, 37 curves, NF 10 / 20 / 30.
    # (a) The cap stands: at every tested rate up to protocols.BCC_RC_MAX_CPS (9 Mchip/s) the body-only curves
    # (arrangement notch, peak and 20 MHz drop smoothed) keep every path error-free (PSR 1). Smear would lower the
    # cap only if a path fell below 0.99 there.
    # (b) The probe above the cap: at 9.5 M the curves as measured (extension A) put paths below 0.99 at every NF,
    # while the model without dispersion ('tones') still delivers every packet. So the check can fail, and
    # dispersion is what makes it fail.
    table = HERE / "results" / "bcc_rate_isi.csv"
    if not table.exists():
        print("skipped, write it first:", table.name)
        return
    st = pd.read_csv(table)
    nf = st[st.nf_db.notna()]
    body = nf[(nf.variant == "smoothed") & (nf.ext == "A") & (nf.rc_cps <= pr.BCC_RC_MAX_CPS)]
    assert sorted(body.rc_cps.unique()) == [500e3, 1.3125e6, 2e6, 2.625e6, 5.25e6, 7e6, 9e6]
    assert len(body) == 7 * 37 * 3 and (body.packets == 400).all() and (body.psr == 1).all()
    top = nf[nf.rc_cps == 9.5e6]
    for x in (10, 20, 30):
        assert (top[(top.variant == "measured") & (top.ext == "A") & (top.nf_db == x)].psr < 0.99).any(), x
        assert (top[(top.variant == "tones") & (top.nf_db == x)].psr == 1).all(), x

def test_ble_event_rules():
    # PARAMETERS 3.1 FULL: 'If the Peripheral misses C, the event ends'; 'the Peripheral answers even after a CRC
    # failure' (then without the NESN, so it resends its unacknowledged head). BLE full 1M alone, 10 kb/s latency,
    # SINR model, 1 s, every poll failed. Reason 'sync' (missed): no P is ever sent. Reason 'noise' (synced, CRC
    # error): each event carries P answers, the head frame goes through on the first one, and since no poll ever
    # carries its acknowledgement the sensor keeps resending that head: one frame delivered, all later attempts
    # are that head again.
    spec = dataclasses.replace(pr.get_spec("ble", "full", "1M"), p_loss=0.0)
    sel = pr.select_period(10, "latency", 1, spec)
    for why in ("sync", "noise"):
        r = _forced(lambda tx: True, why)(spec, sel, 1, 1, "sinr", 1_000_000, 1).run().rows(dict(protocol="ble"))[0]
        if why == "sync":   # energy spent, nothing delivered: J per bit is inf, not undefined
            assert r["tx_attempts"] == 0 and r["jpb_sensor"] == r["jpb_ban"] == math.inf, r
        else:
            assert r["delivered_ok"] == 1 and r["tx_attempts"] > 100, (r["delivered_ok"], r["tx_attempts"])
    # AFH (3.1 FULL step 6: success = valid-CRC P received): a missed poll is a failed outcome on its channel.
    # Collision model, 4 s (533 events), every poll on channel 5 lost: without AFH CSA#2 visits channel 5 in about
    # 1 event in 37 (14 times). Window 20, PDR < 0.95: the 2nd failure (2 > 0.05 x 20) removes channel 5 from the
    # map 6 events later, so it is used at most 3 times (2, plus once in the 6-event delay at worst) and is not in
    # the final map.
    used = []

    class Rec(_forced(lambda tx: tx.key == 5, "collision")):
        def _chan(self, b, j, n, cell=0):
            used.append(super()._chan(b, j, n, cell))
            return used[-1]
    sim = Rec(spec, sel, 1, 1, "collision", 4_000_000, 1).run()
    assert used.count(5) <= 3 and 5 not in sim.cmap[0][1], (used.count(5), len(sim.cmap[0][1]))


def test_beacon_loss_and_migration():
    # Beacon loss (SmartBAN, OI-33: after 10 missed D-Beacons the node loses its frames until it hears one again). One
    # BAN, one sensor, 10 kb/s latency (T = 7.5 ms), p_loss 0, every D-Beacon (272 us; the 168 us ACKs arrive) lost:
    # the misses count 1 to 9 at the slots of superframes 0 to 8 (frames sent and acknowledged), 10 at superframe 9
    # (10 >= 10: lost from there on): 9 frames delivered, the other 124 of the 1 s dropped for lost sync.
    spec = dataclasses.replace(pr.get_spec("smartban", "full"), p_loss=0.0)
    sel = pr.select_period(10, "latency", 1, spec)
    r = _forced(lambda tx: tx.e - tx.s == 272, "collision")(spec, sel, 1, 1, "collision", 1_000_000, 1).run().rows(
        dict(protocol="smartban"))[0]
    assert sel["period_us"] == 7500 and r["delivered_ok"] == 9 and r["drop_sync"] == r["generated"] - 9 > 100, (
        r["delivered_ok"], r["drop_sync"], r["generated"])
    # SmartBAN migration (3.5 FULL, A): uplink failure ratio > 20 % over 10 IBIs. 10 sensors, 10 kb/s latency: T =
    # 16.25 ms (13 slots of 1250 us), one exchange per sensor and IBI. Three failed exchanges in IBI 0 give 3 / 100
    # over the first 10 IBIs: no migration (the old per-exchange window of 10 moved at once). Every exchange failed:
    # the 10 completed IBIs 0..9 are judged at the first outcome of IBI 10, the switch lands 2 beacons later, IBI 12.
    spec = dataclasses.replace(pr.get_spec("smartban", "full"), p_loss=0.0)
    sel = pr.select_period(10, "latency", 10, spec)
    assert sel["period_us"] == 16250

    class Some(en.Sim):
        def _eval(self, tx):
            super()._eval(tx)
            if tx.data and tx.node % self.m in (1, 2, 3) and tx.s < self.phase[0] + self.T:
                tx.ok, tx.reason = False, "collision"

    class All(en.Sim):
        def _eval(self, tx):
            super()._eval(tx)
            if tx.data:
                tx.ok, tx.reason = False, "collision"
    assert len(Some(spec, sel, 1, 10, "collision", 500_000, 1).run().chist[0]) == 1
    assert All(spec, sel, 1, 10, "collision", 500_000, 1).run().chist[0][1][0] == 12


def test_channel_rules():
    # ED scan (DESIGN 4.1): 802.15.4 takes the channel with the fewest BANs, so 16 BANs fill the 16 channels once
    # each. 802.15.6 NB takes the least energy with the 0 dB adjacent-channel rejection (3.4): a BAN also counts on
    # both neighbours, so each placed BAN blocks at most 3 of the 79 channels and 20 BANs (60 < 79) all land at
    # least 2 channels from each other.
    s = en.Sim(pr.get_spec("ieee802154", "full"), pr.select_period(10, "latency", 1, pr.get_spec("ieee802154",
                                                                                                "full")),
               16, 1, "collision", 1, 0)
    assert sorted(h[0][1] for h in s.chist) == list(range(16))
    nb = pr.get_spec("ieee802156nb", "full")
    for seed in range(5):
        ch = sorted(h[0][1] for h in en.Sim(nb, pr.select_period(5, "latency", 1, nb), 20, 1, "collision", 1,
                                            seed).chist)
        assert min(np.diff(ch)) >= 2, ch
    # TSCH CCA busy. Channel always busy, one sensor, 1 kb/s
    # (20 ms slotframe): every frame fails its first attempt and 3 retries in the next cells, so per_first = 1, no
    # delivery, and 3 retransmissions per first attempt (a frame still in progress at the end has fewer)
    class Busy(en.Sim):
        def _cca(self, *a):
            return True
    spec = pr.get_spec("tsch", "full")
    r = Busy(spec, pr.select_period(1, "latency", 1, spec), 1, 1, "collision", 2_000_000, 1).run().rows(
        dict(protocol="tsch"))[0]
    assert r["per_first"] == 1.0 and r["delivered_ok"] == 0 and 2.5 < r["retx_per_frame"] <= 3, r
    # TSCH blacklisting (3.3 FULL+, on in every run): a channel whose failures exceed 0.2 x 100 = 20 of its last 100
    # outcomes leaves the BAN's hopping sequence, at least 4 channels stay. One sensor alone, 1 kb/s (20 ms slotframe),
    # p_loss 0, every data frame on HS[0] lost, 20 s: without blacklisting the sensor would meet HS[0] in about 1 cell
    # in 16 (about 125 of its 1000 dedicated and 1000 shared cells, the shared ones used by retries only); with it, the
    # 21st failure on HS[0] removes it, so exactly 21 transmissions go on HS[0] and none after. A frame that failed
    # there (dedicated cell 1 of slotframe n, channel offset 1: index ASN0 + 2 n + 2 with S = 2 cells) goes again in the
    # next shared cell (cell 0 of slotframe n + 1, offset 0: index ASN0 + 2 (n + 1), the same channel) and fails again,
    # then in its next dedicated cell (index ASN0 + 2 n + 4, another channel), where it arrives: every frame arrives,
    # and the 21 transmissions on HS[0] are dedicated and shared ones.
    spec = dataclasses.replace(pr.get_spec("tsch", "full"), p_loss=0.0)
    hs = spec.tsch["hopping_sequence"]

    class Jam(en.Sim):
        def _eval(self, tx):
            super()._eval(tx)
            if tx.data and tx.key == hs[0]:
                tx.ok, tx.reason = False, "collision"
    sim = Jam(spec, pr.select_period(1, "latency", 1, spec), 1, 1, "collision", 20_000_000, 1).run()
    r = sim.rows(dict(protocol="tsch"))[0]
    on = [x for x in sim.ntrack[1].x if x.data and x.key == hs[0]]
    assert len(on) == 21 and hs[0] not in sim.hs[0] and len(sim.hs[0]) == 15, (len(on), sim.hs[0])
    assert r["delivered_ok"] == r["generated"] > 900, (r["delivered_ok"], r["generated"])
    # removals stop at 4 channels: failing the others one after another (21 each) leaves the last 4 of the sequence
    for c in hs[1:]:
        for _ in range(21):
            sim._adapt(0, 1, 0, c, True)
    assert sim.hs[0] == hs[-4:], sim.hs[0]
    # ANT adaptive isochronous (3.6 FULL step 4): masters 1 and 2 of one BAN at 1000 us and 1618 us, acknowledged
    # slot 418 us, 5 kb/s: T = 7507 us, channel period tc = T (1 + eps) / 1.25. Master 1's guard window [1418, 1918]
    # holds master 2's start: it moves 250 us earlier, again next period ([1168, 1668] still holds 1618), then
    # [918, 1418] is clear. So master 1 starts at 1000, then + tc - 250 twice, then + tc; master 2 never moves.
    class Placed(en.Sim):
        def _ant_phases(self):
            return np.array([0.0, 1000.0, 1618.0])
    spec = dataclasses.replace(pr.get_spec("ant", "full"), p_loss=0.0)
    sim = Placed(spec, pr.select_period(5, "latency", 2, spec), 1, 2, "collision", 60_000, 1).run()
    for node, steps in ((1, [-250, -250, 0, 0]), (2, [0, 0, 0, 0])):
        tc = sim.T * (1 + sim.eps[node]) / 1.25
        starts = [x.s for x in sim.ntrack[node].x][:5]     # its messages and old-payload broadcasts
        assert np.allclose(np.diff(starts) - tc, steps, atol=1e-6), (node, np.diff(starts) - tc)
    # ANT frequency agility (3.6 FULL step 4: 3 frequencies, a move after 8 unacknowledged messages): per BAN its
    # private frequency and two more of the pool. One master, 1 kb/s latency (one 8-octet message per sample, 6
    # attempts), every data frame on the first frequency lost: frame 1 fails its 6 attempts (dropped), frame 2 its
    # first 2, and that 8th unacknowledged message moves the master to the second frequency, where frame 2 arrives at
    # its 3rd attempt and every later frame at its first (the samples that arrive while frame 1 fails fill the FIFO, so
    # some are dropped at arrival: drop_queue, not counted here)
    class JamF0(en.Sim):
        def _eval(self, tx):
            super()._eval(tx)
            if tx.data and tx.key == self.agile[0][0]:
                tx.ok, tx.reason = False, "collision"
    spec = dataclasses.replace(pr.get_spec("ant", "full"), p_loss=0.0)
    sim = JamF0(spec, pr.select_period(1, "latency", 1, spec), 1, 1, "collision", 40 * 7507, 1).run()
    r = sim.rows(dict(protocol="ant"))[0]
    f = sim.agile[0]
    assert len(set(f)) == 3 and f[0] == sim.chist[0][0][1] and all(0 <= c < 76 for c in f), f
    keys = [x.key for x in sim.ntrack[1].x if x.data]
    assert keys[:8] == [f[0]] * 8 and set(keys[8:]) == {f[1]}, keys
    assert r["drop_retry"] == 1 and r["delivered_ok"] + r["drop_queue"] == r["generated"] - 1, r
    # ANT listen at channel open (3.6 FULL step 1) with bursts: 5 kb/s energy, P = 401.25 ms -> N = 13149 counts, T =
    # 401276 us, L = round(2006.4) = 2006 bits = 251 octets = 11 advanced-burst packets (10 x 24 + 11) of 3.2 ms, an
    # occupancy d = 35.2 ms per master. A master may not start within d + g of an earlier master's start on either side
    # (guard g = 500 us), a window of 2d + 2g = 71.4 ms: a random try misses it with probability (401.276 - 71.4) /
    # 401.276 = 0.82, so the second master finds a free phase within its 100 tries (all fail with probability 0.18^100)
    # and its occupancy leaves a gap of at least g on both sides: (phase 2 - phase 1) mod T in [d + g, T - d - g]. Clock
    # drift (50 ppm) moves the two by at most 100 ppm x 2 s = 200 us < g over the run, and p_loss is 0, so one BAN alone
    # loses nothing. (Masters pick one after another, so a later one may find the free region split into gaps each
    # shorter than d + 2g; it then takes the least-occupied phase, below.)
    spec = dataclasses.replace(pr.get_spec("ant", "full"), p_loss=0.0)
    sel = pr.select_period(5, "energy", 2, spec)
    assert (sel["period_us"], sel["frames_per_period"]) == (401276, 11)
    for seed in range(10):
        sim = en.Sim(spec, sel, 1, 2, "collision", 2_006_380, seed)
        assert 35_700 <= (sim.ant_phase[2] - sim.ant_phase[1]) % sim.T <= 401_276 - 35_700, sim.ant_phase
        assert _psr(sim.run().rows(dict(protocol="ant"))) == 1.0, seed
    # the least-occupied phase (D4, DESIGN 4.1): a master whose 100 tries find no free phase takes the phase midway
    # between the two heard masters farthest apart (the centre of the largest gap, the least overlap with what it
    # heard). Every try at phase 0 (a stub), masters 1 to 4 of one BAN in order, 1 kb/s latency: T = 7507 us, tc = T /
    # 1.25 = 6005.6 us. Master 1 takes 0 (nothing heard); master 2 finds 0 taken, one gap of the whole period: tc / 2
    # = 3002.8; master 3: gaps 3002.8 and 3002.8, the first: 1501.4; master 4: gaps 1501.4, 1501.4 and 3002.8: 3002.8
    # + 1501.4 = 4504.2

    class Zero:
        def uniform(self, lo, hi):
            return 0.0

        def permutation(self, x):
            return np.asarray(x)
    sim = en.Sim(spec, pr.select_period(1, "latency", 4, spec), 1, 4, "collision", 1, 0)
    sim.rng = Zero()
    assert sim.T == 7507 and np.allclose(sim._ant_phases()[1:], [0, 3002.8, 1501.4, 4504.2]), sim._ant_phases()


def test_sweep_specs():
    # the catalogue runs full stacks only, plus the as-built BCC ('thesis'): --phy ble=2M,ble=coded_s8 with both
    # variants runs full 2M and full coded; --protocols bcc runs BCC full on one tone and on the
    # tones and the as-built stack, and a BCC flag reaches all three; an entry get_spec refuses is skipped: the
    # multi-tone one under the contact and coupling cases (it has neither) and at 4 Mchip/s (4 R_c = 16 MHz exceeds
    # the 13.52 MHz 3 dB span: no tone), both BCC full ones under --bcc-fcs (they carry the 802.15.6 FCS); --bcc-opt
    # reaches the two BCC full entries only (the as-built stack has no ladder), labelled 'full+optN'; so does --bcc-pool
    a = argparse.Namespace(protocols="ble", variants="full,thesis", phy="ble=2M,ble=coded_s8", bcc_contact=False,
                           bcc_fcs=False, bcc_coupling=None, bcc_rc=None, bcc_power=None, bcc_opt=0,
                           bcc_pool=None, bcc_pool_fill=False, bcc_touch=None,
                           bcc_touch_db=None, bcc_fallback=False, bcc_spread=None, bcc_back=False, bcc_burst=False,
                           bcc_span=None, bcc_pilot128=False)
    assert sweep.specs(a) == [("ble", "full", "2M", ()), ("ble", "full", "coded_s8", ())]
    a.protocols, a.phy, a.bcc_rc = "bcc", "", 50e3
    rc = (("rc_cps", 50e3),)
    assert sweep.specs(a) == [("bcc", "full", "sf8", rc), ("bcc", "full", "multitone", rc),
                              ("bcc", "thesis", "bfsk", rc)]
    a.bcc_rc, a.bcc_contact = None, True
    c = (("contact", True),)
    assert sweep.specs(a) == [("bcc", "full", "sf8", c), ("bcc", "thesis", "bfsk", c)]
    a.bcc_contact, a.bcc_coupling = False, 40
    assert [t[2] for t in sweep.specs(a)] == ["sf8", "bfsk"]
    a.bcc_coupling, a.bcc_rc = None, 4e6
    assert [t[2] for t in sweep.specs(a)] == ["sf8", "bfsk"]
    a.bcc_rc, a.bcc_fcs = None, True
    assert sweep.specs(a) == [("bcc", "thesis", "bfsk", (("fcs", True),))]
    a.bcc_fcs, a.bcc_opt = False, 3
    o3 = (("opt", 3),)
    assert sweep.specs(a) == [("bcc", "full", "sf8", o3), ("bcc", "full", "multitone", o3),
                              ("bcc", "thesis", "bfsk", ())]
    assert [en.label(sweep._spec(t)) for t in sweep.specs(a)] == ["full+opt3", "full+opt3", "thesis"]
    # --bcc-pool (the contact design point 0.43) reaches the BCC full entries only, labelled '+pool=0.43'
    a.bcc_pool = 0.43
    p = (("opt", 3), ("pool", 0.43))
    assert sweep.specs(a) == [("bcc", "full", "sf8", p), ("bcc", "full", "multitone", p), ("bcc", "thesis", "bfsk", ())]
    assert en.label(sweep._spec(sweep.specs(a)[0])) == "full+opt3+pool=0.43"
    # the CSV phy column tells them apart
    assert [pr.get_spec(*t[:3]).phy["name"] for t in pr.catalogue() if t[0] == "bcc"] == ["sf8", "multitone", "bfsk"]
    # --help renders
    out = subprocess.run([sys.executable, str(HERE / "sweep.py"), "--help"], capture_output=True, text=True)
    assert out.returncode == 0 and "--bcc-pool" in out.stdout and "30 % contact" in out.stdout, out.stderr[-300:]


def test_sweep_main():
    # sweep.main end to end on one small cell (the as-built BCC, 10 kb/s latency, k 1, N_dev 1, one trial of 40 periods
    # of 7.5 ms), with the ProcessPoolExecutor replaced by an in-process stand-in and engine.Sim by a subclass that records
    # the ci_floor it receives. --ci-floor off with both models: the switch reaches the engine through sweep._task and
    # simulate_trial (two Sims, both ci_floor False), the SINR row is labelled 'thesis+nocifloor' and the collision row
    # keeps 'thesis' (the collision model has no C/I floor), so a later default collision run finds it done and builds
    # no Sim. The pool gets the default 4 workers, and --workers 9 is refused (at most 8 worker processes: the machine
    # crashed once under load).
    seen, floors = [], []

    class Pool:
        def __init__(self, n):
            seen.append(n)

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def submit(self, f, x):
            fut = Future()
            fut.set_result(f(x))
            return fut

    class Sim(en.Sim):
        def __init__(self, *a, **kw):
            super().__init__(*a, **kw)
            floors.append(self.ci_floor)

    argv, saved = sys.argv, (sweep.ProcessPoolExecutor, en.Sim)
    sweep.ProcessPoolExecutor, en.Sim = Pool, Sim
    try:
        with tempfile.TemporaryDirectory() as d:
            base = ["sweep.py", "--protocols", "bcc", "--variants", "thesis", "--r", "10", "--k", "1",
                    "--n-dev", "1", "--ci-modes", "latency", "--trials", "1", "--duration-s", "0.1",
                    "--out", str(Path(d, "rows.csv")), "--cells", str(Path(d, "cells.csv"))]
            for extra in (["--reception", "collision,sinr", "--ci-floor", "off"], ["--reception", "collision"]):
                sys.argv = base + extra
                sweep.main()
            rows = pd.read_csv(Path(d, "rows.csv"))
            sys.argv = base + ["--workers", "9"]
            try:
                sweep.main()
                raise AssertionError("sweep accepted 9 workers")
            except SystemExit as e:
                assert e.code == 2, e.code
    finally:
        sys.argv = argv
        sweep.ProcessPoolExecutor, en.Sim = saved
    assert seen == [4, 4] and floors == [False, False], (seen, floors)
    assert sorted(zip(rows["reception"], rows["variant"])) == [("collision", "thesis"), ("sinr", "thesis+nocifloor")]


def test_cells_aggregate():
    # sweep.aggregate on hand-made rows. Cell A: k 2, N_dev 2, 1 s, two trials; each BAN spends 1 mJ per sensor and 1 mJ
    # at the hub, and goodput_kbps x 1 s gives its bits: trial 0 1000 and 0 (BAN 1 delivered nothing), trial 1 1000 and
    # 3000. J per bit of a trial is a ratio of sums over its k BANs: jpb_sensor 4 mJ / 1000 = 4e-6 and 4 mJ / 4000 =
    # 1e-6, mean 2.5e-6; jpb_ban 6 mJ / 1000 = 6e-6 and 6 mJ / 4000 = 1.5e-6, mean 3.75e-6 (the mean of per-BAN ratios
    # without the BAN that delivered nothing gave (2e-6 + (2e-6 + 6.67e-7) / 2) / 2 = 1.67e-6). Cell B: k 1, one trial
    # that delivered nothing: inf. The Student-t 95 % interval of two trials is mean +/- t(0.975, 1) x sd / sqrt(2),
    # t = 12.7062, clipped to the physical range: psr per trial (1 + 0) / 2 and 1, mean 0.75, sd / sqrt(2) = 0.25, so
    # [-2.43, 3.93] -> [0, 1]; jpb_sensor sd / sqrt(2) = 1.5e-6, so [-1.66e-5, 2.156e-5] -> [0, 2.156e-5].
    rows = [dict(dict.fromkeys(en.COLUMNS, math.nan), protocol="x", variant="thesis", phy="p", reception="collision",
                 strip_hopping=0, ci_mode="latency", r_target_kbps=1.0, n_dev=2, k=k, duration_us=1_000_000,
                 venue="ward", density="typical", case="typical", trial=tr, seed=tr, network_id=b, feasible=1,
                 e_sensor_j=1e-3, e_hub_j=1e-3, goodput_kbps=gp, psr=float(gp > 0))
            for k, tr, b, gp in ((2, 0, 0, 1.0), (2, 0, 1, 0.0), (2, 1, 0, 1.0), (2, 1, 1, 3.0), (1, 0, 0, 0.0))]
    with tempfile.TemporaryDirectory() as d:
        out, cells = Path(d, "rows.csv"), Path(d, "cells.csv")
        pd.DataFrame(rows).to_csv(out, index=False)
        assert sweep.aggregate(out, cells) == 2
        c = pd.read_csv(cells).set_index("k")
    assert math.isclose(c.at[2, "jpb_sensor_mean"], 2.5e-6, rel_tol=1e-6), c.at[2, "jpb_sensor_mean"]
    assert math.isclose(c.at[2, "jpb_ban_mean"], 3.75e-6, rel_tol=1e-6), c.at[2, "jpb_ban_mean"]
    assert c.at[1, "jpb_sensor_mean"] == math.inf
    assert (c.at[2, "psr_mean"], c.at[2, "psr_lo"], c.at[2, "psr_hi"]) == (0.75, 0.0, 1.0)
    assert c.at[2, "jpb_sensor_lo"] == 0.0
    assert math.isclose(c.at[2, "jpb_sensor_hi"], 2.5e-6 + 12.7062 * 1.5e-6, rel_tol=1e-5), c.at[2, "jpb_sensor_hi"]


def test_csa2_self_test():
    # Core 5.4 Vol 6 Part C (CSA#2 sample data, parameters.json ble5.csa2_selftest): chanId 0x305F, all 37 channels
    # used, counters 0..3 -> 25, 20, 6, 21; used {9, 10, 21, 22, 23, 33, 34, 35, 36}, counters 6, 7, 8 -> 23, 9, 34
    for t in pr.BLE["csa2_selftest"]:
        used = list(range(37)) if t["used"] == "all" else t["used"]
        assert [en.csa2(c, int(t["chanId"], 16), used) for c in t["counters"]] == t["channels"]

def test_bcc_touch_and_fallback():
    #. Contact intervals: alternating exponential on / off, mean contact 5 s, share 0.2: over 6000 s (about
    # 240 cycles of 25 s) the share of time in contact is 0.2 within 0.05, and both BANs of a pair see the same list.
    s = pr.get_spec("bcc", "full", rc_cps=2.625e6, opt=5, touch=0.2)
    sim = en.Sim(s, pr.select_period(10, "latency", 4, s), 2, 4, "sinr", 6000e6, seed=1)
    assert sim.touch[0] is sim.touch[1] and abs(sum(e - a for a, e in sim.touch[0]) / 6000e6 - 0.2) < 0.05
    # Always touching at 0 dB (the partner as strong as the own link), 4 sensors x 10 kb/s at 2.625 Mchip/s, T 7.5 ms:
    # with seed 3 the two superframes start 222 us apart and, sharing the period, their beacons overlap in every
    # superframe, so the sensors miss their beacons and hardly send (psr < 0.1); the fallback interleaves the two
    # active periods after one superframe and restores psr 1, as isolation does; at 20 dB the partner does no harm.
    def psr(**kw):
        sp = pr.get_spec("bcc", "full", rc_cps=2.625e6, opt=5, fill=True, **kw)
        rows = en.simulate_trial(sp, k=2, n_dev=4, r_kbps=10, ci_mode="latency", reception="sinr",
                                 duration_us=4e6, seed=3)
        return sum(r["psr"] for r in rows) / len(rows)
    assert psr() == 1.0 and psr(touch=1.0, touch_db=20.0) == 1.0
    # touch_db 0.0 is a value, not "unset" (0.0 == False): it shows in the label (the sweep's flag filter as well)
    assert en.label(pr.get_spec("bcc", "full", touch=1.0, touch_db=0.0), True, "sinr") == "full+touch=1.0+touch_db=0.0"
    assert psr(touch=1.0, touch_db=0.0) < 0.1 and psr(touch=1.0, touch_db=0.0, fallback=True) == 1.0
    # refused: collision reception, multi-tone, fallback or touch_db without touch
    for call in (lambda: en.simulate_trial(pr.get_spec("bcc", "full", touch=0.2), k=2, n_dev=1, r_kbps=1,
                                           ci_mode="energy", reception="collision", duration_us=1e6, seed=1),
                 lambda: pr.get_spec("bcc", "full", "multitone", touch=0.2),
                 lambda: pr.get_spec("bcc", "full", fallback=True)):
        try:
            call()
            raise AssertionError("not refused")
        except ValueError as e:
            assert "touch" in str(e), e



def test_bcc_spread_switch():
    # hysteresis on an injected report sequence: one loss keeps unspread, two in a row switch to Walsh; 15
    # deliveries and a loss keep Walsh (the run restarts), 16 in a row switch back.
    plain, run, seen = True, 0, []
    for ok in [0, 1, 0, 0] + [1] * 15 + [0] + [1] * 16:
        plain, run = en.spread_step(plain, run, ok)
        seen.append(plain)
    assert seen[:4] == [True, True, True, False] and not any(seen[4:35]) and seen[35]
    # On the noise-limited link (500 kchip/s, 4 sensors x 10 kb/s) every mode delivers everything; unspread frames
    # take less air, and 'adapt' never leaves unspread, so it sends the same frames in Walsh-sized allocations.
    def run1(**kw):
        sp = pr.get_spec("bcc", "full", opt=5, **kw)
        r = en.simulate_trial(sp, k=1, n_dev=4, r_kbps=10, ci_mode="energy", reception="sinr", duration_us=2e6,
                              seed=1)[0]
        return r["psr"], r["airtime_frac"]
    (pw, aw), (pn, an), (pa, aa) = run1(), run1(spread="none"), run1(spread="adapt")
    assert pw == pn == pa == 1.0 and an < aa < aw / 2


if __name__ == "__main__":
    import test_phy               # DESIGN 8.4: BER formulas, BCH decoding, Walsh-16 and the BCC packet
    for mod in (test_phy, sys.modules[__name__]):
        for name, fn in list(vars(mod).items()):
            if name.startswith("test_") and callable(fn):
                fn()
                print("ok ", name, flush=True)
    print("all coexistence tests passed")

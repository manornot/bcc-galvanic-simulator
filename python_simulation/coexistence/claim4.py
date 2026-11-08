"""Capacity numbers: the galvanic profile (BCC
full at ladder level 5, pool 0.05) against the pooled reference (the HBC stack with hub-polled retries, level 0, pool
0.05), one tone unless named.

Usage, from python_simulation/coexistence/:   python claim4.py [--workers N] [--raw DIR]

  A  air time per frame beyond the payload and the capacity ratio by frame length (8, 64, 255 octets), closed form
     8 L / report.frame_cost(L) (DESIGN 4.3);
  B  the capacity ratio at maximum frames by pool share (the same share on both hubs, 0.01 to 1.0, and the limit share
     -> 0) and chip rate (50 kchip/s to 2.625 Mchip/s), the spare-exchange policy on both hubs (L2 against L0), six
     tones at 500 kchip/s; unequal pools (reference share -> 0 and 0.01, profile 0.05);
  C  attribution by frame length: the two options (L1 against L0, both with the pool) and the three changes (L5
     against L1 with the pool), and the other order: the three changes on L0 and the two options on top;
  D  simulated at the admitted capacity (report.capacity): both stacks, N_dev 1 and 4, k 1, collision model, 2 seeds,
     both ci modes, at the capacity and 0.5 kb/s above it; maximum frames at 50 kchip/s, 500 kchip/s (one and six
     tones) and 2.625 Mchip/s, and the short frames (payload limit 64 and 8 octets, protocols flag lmax) at 500
     kchip/s and 2.625 Mchip/s;
  E  energy per delivered bit, paired per cell (reference / proposal) over the cells feasible for both: results/
     cells_lossfree.csv,
     cells_ladder.csv (500 kchip/s), cells_rates.csv (2.625 Mchip/s); the sensor and the BAN (hub included); the
     transceiver power presets at 50 kchip/s (cells_lossfree.csv, run_sweeps.sh pooled);
  F  engine reruns, k 1, collision, 2 seeds: every one-tone (rate, N_dev) point of E with the proposal at the
     reference's own period (both ci modes), and the power presets at 500 kchip/s and 2.625 Mchip/s;
  G  the reference with the standard's block acknowledgement (L-Ack, one B-Ack per block) instead of the I-Ack,
     closed form: ratios at 50, 500, 2625 kchip/s, step 1 (the two options) with block ack, equal pool
     shares and B-Ack block sizes over 461 chip rates. --only-g recomputes G alone (seconds, no engine run), --closed
     A, B, C, E and G (the parts without an engine run);
  H  the block-ack reference in the engine (protocols flag back): H1 simulated at its admitted
     capacity, the cases of D (k 1, collision, 2 seeds, both ci modes) beside the profile and the closed form of G;
     H2 energy per delivered bit, reference with block ack / profile, paired per cell as E: results/cells_back.csv
     (run_sweeps.sh back) against the profile cells of E. --only-h recomputes H alone (minutes);
  G2, I  the lead reference: the block-ack reference in the HBC burst mode (protocols flags back
     and burst), G2 closed form by frames per allocation n (also with the profile sending its own frames pMIFS
     apart), I1 at the admitted capacity as H1, I2 energy against results/cells_burst.csv (run_sweeps.sh burst), also
     against the best of the three references per cell. --only-i recomputes I alone;
  J  the lead reference: the burst reference of I with its block spanning allocations
     (protocols flag span, 6.2.9.4: a block may span more than one allocation interval), an unlimited block (no B-Ack;
     the limit of an unlimited block, below every block the standard allows, 6.2.9.4); G3 gives its closed form
     (block M), G4 its typical value over the working points A-H and their J2 energy cells. J1 its admitted
     capacity and packet success ratio as I1, J2 energy against results/cells_span.csv (run_sweeps.sh span), also against the best
     of the four references per cell. --only-j recomputes J alone.
  G5 the standard guard time GTn (IEEE 802.15.6-2012 6.11.1.1, Equations (6) to (8), Table 24; 6.11.1.3) charged to
     both stacks, also in the pool's polled allocations: G4 with
     the guard in the capacity reading and with the beacon period equal to the longest
     connection interval, a zero-guard self-check against G4, the long-burst floor with a guard per allocation, the
     'any frame' bound, and the J2 energy with the 6.11.1.2 beacon receive window on both stacks. Part of G (--only-g).
  G6 G5 on the profile's allocation slot: each scheduled allocation of both stacks in whole slots, the
     beacon in slots or continuous, the standard's 500 us beside it, the 'any frame' bound. Part of G (--only-g).
  G7 the ratios of G6 split into step 1 and steps 2 to 4 at the cell where each lowest value falls. Part of G.

When DIR/t_l091/lossfree_cells.csv exists (default DIR results/raw, git-ignored) it rewrites results/cells_lossfree.csv
with DIR/t_l091/power_cells.csv (6 significant digits), DIR/back/back_cells.csv results/cells_back.csv and
DIR/burst/burst_cells.csv results/cells_burst.csv, DIR/span/span_cells.csv results/cells_span.csv. Output:
results/tables/t15_claim4.md, printed as well.
"""
import argparse
import dataclasses
import math
import re
import time
from functools import lru_cache
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

import engine as en
import protocols as pr
import report as rp

RES = Path(__file__).resolve().parent / "results"
REF, PROP = (0, 0.05), (5, 0.05)         # (ladder level, pool share)
RATES = (50e3, 200e3, 500e3, 656.25e3, 1.3125e6, 2e6, 2.625e6)
HEAD = (50e3, 500e3, 2.625e6)
SHARES = (0.01, 0.05, 0.1, 0.2, 0.43, 1.0)
FRAMES = (8, 64, 255)
CAP_CASES = [("bcc-sf8", 50e3, None), ("bcc-sf8", 500e3, None), ("bcc-multitone", 500e3, None),
             ("bcc-sf8", 2.625e6, None), ("bcc-sf8", 500e3, 64), ("bcc-sf8", 500e3, 8), ("bcc-sf8", 2.625e6, 64),
             ("bcc-sf8", 2.625e6, 8)]
CELL = ["r_target_kbps", "n_dev", "k"]
KEEP = ["protocol", "variant", "phy", "reception", "ci_mode"] + CELL + [
    "period_us", "feasible", "n_trials", "psr_mean", "on_time_mean", "jpb_sensor_mean", "jpb_sensor_lo",
    "jpb_sensor_hi", "jpb_ban_mean"]


def spec(rc, lv, pool, lmax=None, sid="bcc-sf8", back=False, burst=False, span=None, pilot128=False):
    return rp.bcc_spec(sid, rc, lv, pool, None, lmax, back, burst, span, pilot128)


def cost(rc, lv, pool, octets, sid="bcc-sf8"):
    return rp.frame_cost(spec(rc, lv, pool, sid=sid), octets)


def payload_us(rc, octets):
    return 32 * octets / rc * 1e6       # Walsh-16: 2 symbols of 16 chips per octet


def md(header, rows):
    return "\n".join(["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
                     + ["| " + " | ".join(map(str, r)) + " |" for r in rows]) + "\n"


def part_a():
    rows, ov = [], []
    for rc in HEAD:
        for L in FRAMES:
            p, cr, cp = payload_us(rc, L), cost(rc, *REF, L), cost(rc, *PROP, L)
            ov.append((cp - p) / (cr - p))
            rows.append([f"{rc / 1e3:g}", L, f"{cr:.0f} / {cp:.0f}", f"{(cr - p) / 1e3:.2f} / {(cp - p) / 1e3:.2f}",
                         f"{(cp - p) / (cr - p):.3f}", f"{p / cr:.2f} / {p / cp:.2f}", f"{cr / cp:.3f}"])
    return ("A. Air time per frame, pooled reference / profile, pool 0.05, one tone (closed form, DESIGN 4.3)\n\n"
            + md(["chip rate (kchip/s)", "payload (octets)", "frame cost (us)", "beyond the payload (ms)",
                  "beyond the payload, profile / reference", "payload share of the air time", "capacity ratio"], rows)
            + f"\nAir time beyond the payload, profile / reference: {min(ov):.3f} to {max(ov):.3f}.\n")


def part_b():
    out, lo = [], {}
    for L in FRAMES:
        rows = []
        for rc in RATES:
            r = [cost(rc, 0, s, L) / prof_cost(rc, L, s, n=1) for s in SHARES]
            lim = pr._frame_us(spec(rc, *REF), L) / prof_cost(rc, L, 0, n=1)
            spare = cost(rc, 0, None, L) / cost(rc, 2, None, L)
            lo.setdefault(L, []).extend([lim, *r])
            rows.append([f"{rc / 1e3:g}", f"{lim:.3f}"] + [f"{x:.3f}" for x in r] + [f"{spare:.3f}"])
        out.append(f"\nPayload {L} octets:\n\n" + md(["chip rate (kchip/s)", "share -> 0"] + [f"share {s}" for s in SHARES]
                                                    + ["spare exchange on both (L2 / L0)"], rows))
    # unequal pools: frame_cost is linear and rising in the share, so the ratio rises with the reference's share and
    # falls with the profile's; a reference pool not smaller than the profile's keeps the ratio within the range above
    asym = "; ".join(
        f"{L} octets, reference share -> 0 {min(a):.3f} to {max(a):.3f}, 0.01 {min(b):.3f} to {max(b):.3f}"
        for L in FRAMES
        for a, b in [([pr._frame_us(spec(rc, *REF), L) / prof_cost(rc, L, 0.05, n=1) for rc in RATES],
                      [cost(rc, 0, 0.01, L) / prof_cost(rc, L, 0.05, n=1) for rc in RATES])])
    six = [cost(500e3, 0, s, 1530, "bcc-multitone") / cost(500e3, 5, s, 1530, "bcc-multitone") for s in SHARES]
    short = min(min([cost(rc, 0, s, L) / prof_cost(rc, L, s, n=1) for s in SHARES]
                    + [pr._frame_us(spec(rc, *REF), L) / prof_cost(rc, L, 0, n=1)])
                for rc in RATES for L in range(1, 65))
    return ("B. Capacity ratio profile / pooled reference at the same pool share on both hubs (closed form; the"
            " profile with its share of the beacon's ack field at one frame per superframe, the least favourable, a bit"
            " and the 2-octet base per frame, claim4.prof_cost n 1; six tones without it)\n"
            + "".join(out) + "\nRange over every share and chip rate above: "
            + "; ".join(f"{L} octets {min(v):.3f} to {max(v):.3f}" for L, v in lo.items())
            + f"\nSmallest ratio for every payload of 1 to 64 octets, every share and chip rate above: {short:.3f}"
            + f"\nUnequal pools, profile at share 0.05, every chip rate above: {asym}"
            + "\nSix tones, 500 kchip/s, maximum frame (1530 octets), by share: "
            + ", ".join(f"{s} {x:.3f}" for s, x in zip(SHARES, six)) + "\n")


class _NoOptions(int):
    """Ladder level 5 for every '>= n' test but falsy, so neither the RI SFD nor the doubled pilot period applies
    (protocols._hbc_framing)."""
    def __bool__(self):
        return False


def changes_only_cost(rc, octets):
    """Frame cost of the three changes (L2, L4, L5) without the two options of L1: the L5 spec with the L0 phy (no RI
    mode) and the pilot period of L0."""
    s = spec(rc, *PROP)
    s = dataclasses.replace(s, phy=spec(rc, *REF).phy, flags=dict(s.flags, opt=_NoOptions(5)))
    return rp.frame_cost(s, octets)


def part_c():
    rows = []
    for rc in HEAD:
        for L in FRAMES:
            c0, c1, c5 = cost(rc, *REF, L), cost(rc, 1, 0.05, L), cost(rc, *PROP, L)
            ch = changes_only_cost(rc, L)
            rows.append([f"{rc / 1e3:g}", L, f"{c0 / c1:.3f}", f"{c1 / c5:.3f}", f"{c0 / ch:.3f}", f"{ch / c5:.3f}",
                         f"{c0 / c5:.3f}"])
    return ("C. Attribution, pool 0.05 on every level (closed form); the ladder is cumulative, so each pair of columns"
            " multiplies to the total: the two options first and the three changes on top, or the three changes first"
            " (on L0) and the two options on top\n\n"
            + md(["chip rate (kchip/s)", "payload (octets)", "two options (L1 / L0)", "three changes (L5 / L1)",
                  "three changes on L0", "two options on top", "profile / reference"], rows))


GRID = [50e3 + 250 * i for i in range(401)] + [float(x) for x in np.geomspace(150e3, 2625e3, 61)[1:]]
UNLIMITED = pr.BLOCK_UNLIMITED


def back_cost(rc, lv, L, pool, m=UNLIMITED, grow=False, burst=False, block=None):
    """Frame cost with the standard's block acknowledgement instead of the I-Ack (L-Ack frames, one B-Ack frame of the
    level's framing per block of m frames, pSIFS after each; PARAMETERS 3.7b, 6.2.9.4): 3 B-Ack payload octets, or
    2 + ceil(m / 8) with a bitmap that grows with the block (the 2012 format is not in the repository). The pool
    exchange keeps its I-Ack, which favours the profile only at large pools; pool 0 is the limit of no pool.
    burst (L0): the block of m frames of one allocation in the HBC burst mode (6.2.9.4, 10.5.2,
    10.11.3): the first frame in DRF mode, m - 1 x (pMIFS + RI-mode frame), pSIFS, B-Ack, pSIFS. block: one B-Ack per block of block >= m frames spanning allocations (6.2.9.4); None: m, a B-Ack every
    allocation."""
    if burst:
        a0, nxt, b, pu = _burst_terms(rc, L, pool or 0.05)
        return (a0 + (m - 1) * nxt + pr.HBC["psifs_us"]) / m + b / (block or m) + pool * pu
    s = spec(rc, lv, pool or 0.05)
    ps = pr.HBC["psifs_us"]
    b = pr.airtime_us(s, 2 + math.ceil(m / 8) if grow else 3) + ps
    return pr.airtime_us(s, L) + ps + b / m + pool * pr._pool_us(s, L)


@lru_cache(maxsize=None)
def _burst_terms(rc, L, pool):
    """(DRF frame, pMIFS + RI frame, B-Ack + pSIFS, pool exchange) of the burst reference, us; with the '128 B' pilot period (10.7.3, Table 86; protocols flag pilot128), the strongest the standard offers."""
    s = spec(rc, 0, pool, back=True, burst=True, pilot128=True)
    return (pr.airtime_us(s, L), pr.airtime_us(s.ri, L) + pr.HBC["pmifs_us"],
            pr.airtime_us(s, 3) + pr.HBC["psifs_us"], pr._pool_us(s, L))


@lru_cache(maxsize=None)
def field_terms(rc):
    """(one bit, one base sequence number) of the ack field in the profile's beacon, us, measured from
    protocols.beacon_layout as the slope over 8000 bits and 8000 sensors (pilots included)."""
    sp = spec(rc, *PROP)
    t = lambda d, n: pr.beacon_layout(sp, d, n)[-1]["t1_us"]     # noqa: E731
    return (t(1, 8000) - t(1, 0)) / 8000, (t(8000, 0) - t(0, 0)) / 8000


def prof_frame(rc, L, pool):
    """Profile frame cost without the beacon (report.frame_cost; pool 0: the limit of no pool)."""
    return cost(rc, 5, pool, L) if pool else pr._frame_us(spec(rc, *PROP), L)


@lru_cache(maxsize=None)
def prof_cost(rc, L, pool, m=None, n=None):
    """Profile frame cost with its share of the beacon's ack field: one bit per frame and the
    sensor's base sequence number over its n frames per superframe (default m, else the limit of many frames). The
    rest of both beacons is left out: the reference's is about 1.6 to 1.7 times the profile's, so leaving it out lowers
    every ratio below that and raises those above it . m: the profile
    also sends the m frames of one allocation pMIFS apart instead of pSIFS (the standard's burst spacing, an option of
    the author's)."""
    bit, base = field_terms(rc)
    c = prof_frame(rc, L, pool) + bit + base / (n or m or UNLIMITED)
    return c - (1 - 1 / m) * (pr.HBC["psifs_us"] - pr.HBC["pmifs_us"]) if m else c


def part_g():
    rows = []
    for rc in HEAD:
        for L in FRAMES:
            ref, prof, opt = cost(rc, *REF, L), prof_cost(rc, L, 0.05, n=1), back_cost(rc, 1, L, 0.05)
            rows.append([f"{rc / 1e3:g}", L, f"{ref / prof:.4f}", f"{back_cost(rc, 0, L, 0.05) / prof:.4f}",
                         f"{ref / opt:.3f}", f"{opt / prof:.3f}"])
    shares = (0.0,) + SHARES
    share = []
    for L in FRAMES:
        r = [(s, [back_cost(rc, 0, L, s) / prof_cost(rc, L, s, n=1) for rc in GRID]) for s in shares]
        share.append(f"{L} octets: " + "; ".join(f"{s if s else '-> 0'} {min(v):.4f} to {max(v):.4f}" for s, v in r))
    short = min(back_cost(rc, 0, L, s) / prof_cost(rc, L, s, n=1) for s in shares for rc in GRID for L in range(1, 65))
    blocks = []
    for L in FRAMES:
        f = [f"m {m if m < UNLIMITED else 'unlimited'}{', bitmap grows' if g else ''} "
             + f"{min(back_cost(rc, 0, L, 0.05, m, g) / prof_cost(rc, L, 0.05, n=m) for rc in GRID):.4f}"
             for m in (4, 15, 39, 280, UNLIMITED) for g in (False, True) if not (g and m == UNLIMITED)]
        blocks.append(f"{L} octets: " + "; ".join(f))
    burst = []
    ns = (1, 2, 3, 4, 8, 15, 39, UNLIMITED)
    for L in FRAMES:
        cells = []
        for m in ns:
            a = min(back_cost(rc, 0, L, 0.05, m, burst=True) / prof_cost(rc, L, 0.05, n=m) for rc in GRID)
            b = min(back_cost(rc, 0, L, 0.05, m, burst=True) / prof_cost(rc, L, 0.05, m) for rc in GRID)
            z = min(back_cost(rc, 0, L, s, m, burst=True) / prof_cost(rc, L, s, n=m) for s in (0.0,) + SHARES
                    for rc in GRID)
            cells.append(f"{a:.3f} / {b:.3f} / {z:.3f}")
        burst.append([L] + cells)
    bshort = [min(back_cost(rc, 0, L, 0.05, m, burst=True) / prof_cost(rc, L, 0.05, n=m) for rc in GRID
                 for L in range(1, 65)) for m in ns]
    sh = (0.0,) + SHARES
    ball = [(min(back_cost(rc, 0, L, s, m, burst=True) / prof_cost(rc, L, s, n=m) for s in sh for rc in GRID
                 for L in range(1, 65)),
             min(back_cost(rc, 0, L, s, m, burst=True) / prof_cost(rc, L, s, m) for s in sh for rc in GRID
                 for L in range(1, 65)),
             min(back_cost(rc, 0, 255, s, m, burst=True) / prof_cost(rc, 255, s, m) for s in sh for rc in GRID))
            for m in ns]
    span = []
    for m in ns:
        cells = []
        for M in sorted({m, max(m, 9), max(m, 255), UNLIMITED}):
            r = [min(back_cost(rc, 0, L, s, m, burst=True, block=M) / prof_cost(rc, L, s, n=m) for s in sh for rc in GRID
                     for L in Ls) for Ls in ((255,), range(1, 65))]
            cells.append(f"M {M if M < UNLIMITED else 'unlimited'} {r[0]:.4f} / {r[1]:.4f}")
        span.append(f"n {m if m < UNLIMITED else 'unlimited'}: " + "; ".join(cells))
    head = [f"{rc / 1e3:g} kchip/s: " + ", ".join(
        f"{L} octets {back_cost(rc, 0, L, 0.05, 1, burst=True) / prof_cost(rc, L, 0.05, n=1):.3f} / "
        f"{back_cost(rc, 0, L, 0.05, burst=True) / prof_cost(rc, L, 0.05):.3f}" for L in FRAMES) for rc in HEAD]
    return ("G. Against the standard's block acknowledgement (closed form, no engine run): the"
            " reference with L-Ack and one B-Ack frame per block of m frames (back_cost), the pool exchange with its"
            " I-Ack. The profile carries its share of the beacon's ack field (prof_cost: one"
            " bit per frame and the 2-octet base sequence number over the sensor's n frames per superframe), and the"
            " burst references of G2 to G4 the '128 B' pilot period (_burst_terms). The table and the share lines below"
            " take the profile at one frame per superframe (n 1), the least favourable for it\n\n"
            + md(["chip rate (kchip/s)", "payload (octets)", "profile / I-Ack reference",
                  "profile / block-ack reference (unlimited block)", "step 1 (two options) with block ack / I-Ack reference",
                  "profile on top of step 1"], rows)
            + f"\nSame pool share on both hubs, unlimited block, {len(GRID)} chip rates 50 kchip/s to 2.625 Mchip/s"
            + " (0.25 kchip/s apart up to 150 kchip/s, 60 geometric above), by share:\n\n"
            + "\n".join(f"- {x}" for x in share)
            + "\n\nSmallest ratio for every payload of 1 to 64 octets, every share above, every chip rate of the grid: "
            + f"{short:.4f}\n\nLowest ratio over the grid by B-Ack block size m, pool 0.05 on both:\n\n"
            + "\n".join(f"- {x}" for x in blocks)
            + "\n\nG2. The lead reference: block ack in the HBC burst mode (back_cost burst; 6.2.9.4,"
            " 10.5.2, Table 89): per allocation of n frames the first in DRF mode, the other n - 1 pMIFS (20 us) apart in"
            " RI mode (header at the payload rate, 608-chip SFD/RI), pSIFS, the 3-octet B-Ack, pSIFS; pool exchanges"
            " with their I-Ack. Lowest ratio over the grid by frames per allocation n: profile as it is at pool 0.05 /"
            " profile also with pMIFS between its own n frames at pool 0.05 / profile as it is, lowest over the shares"
            " (-> 0 and " + ", ".join(map(str, SHARES)) + ", the same on both hubs)\n\n"
            + md(["payload (octets)"] + [f"n {m if m < UNLIMITED else 'unlimited'}" for m in ns], burst)
            + "\nSmallest ratio for every payload of 1 to 64 octets, pool 0.05, by n: "
            + ", ".join(f"n {m if m < UNLIMITED else 'unlimited'} {x:.3f}" for m, x in zip(ns, bshort))
            + "\n\nEvery share (-> 0 and the shares above, the same on both hubs), by n: every payload of 1 to 64 octets,"
            " profile as it is / profile also with pMIFS; 255 octets, profile also with pMIFS: "
            + ", ".join(f"n {m if m < UNLIMITED else 'unlimited'} {a:.4f} / {b:.4f}; {c:.4f}"
                        for m, (a, b, c) in zip(ns, ball))
            + "\n\nAt the three simulated chip rates, pool 0.05, n 1 / unlimited: " + "; ".join(head)
            + "\n\nG3. The lead reference: the burst reference of G2 with one B-Ack per block of"
            " M >= n frames spanning allocations (6.2.9.4; back_cost block), per frame (a0 + (n - 1) (pMIFS + RI frame) +"
            " pSIFS) / n + (B-Ack + pSIFS) / M + share x pool exchange. M = n is G2 (a B-Ack every allocation), M = 9"
            " the largest block whose losses one B-Ack bitmap always reports (mBAckLimit 8), M = 255 the largest Next"
            " Block Size of a one-octet field, M unlimited the lead reference. Lowest ratio, profile as it is, over every share (-> 0 and the shares above, the same on both"
            " hubs) and the grid, 255 octets / every payload of 1 to 64 octets, by frames per allocation n:\n\n"
            + "\n".join(f"- {x}" for x in span) + "\n" + part_g4() + part_g5() + part_g6())


# working points A-H : payload L (octets), R_agg (kb/s), N_dev, k
POINTS = {"A": ((8, 8), (0.006, 0.048), (1, 3), (1, 10)), "B": ((8, 20), (0.06, 1.6), (1, 5), (1, 20)),
          "C": ((12, 32), (1.9, 38), (1, 3), (1, 15)), "D": ((20, 50), (16, 107), (1, 2), (1, 10)),
          "E": ((32, 64), (341, 683), (10, 10), (1, 5)), "F": ((32, 128), (26, 137), (1, 1), (1, 5)),
          "G": ((125, 250), (100, 267), (1, 1), (1, 3)), "H": ((125, 250), (400, 2133), (4, 8), (3, 10))}


@lru_cache(maxsize=None)
def beacons(rc, n_dev):
    """Beacon + pSIFS per superframe of the lead reference and of the profile (with its ack field for n_dev sensors
    of one frame each), us (protocols.beacon_layout)."""
    sr = spec(rc, 0, 0.05, back=True, burst=True, span=UNLIMITED, pilot128=True)
    return tuple(pr.beacon_layout(s, n_dev, 1)[-1]["t1_us"] + pr.HBC["psifs_us"] for s in (sr, spec(rc, *PROP)))


def own(rc, L, s, n_dev):
    """Ratio of the superframe capacities, one frame per sensor and superframe, n_dev sensors, both beacons (the
    profile's with its ack field) and the share s of the pool on both."""
    br, bp = beacons(rc, n_dev)
    return ((br / n_dev + back_cost(rc, 0, L, s, 1, burst=True, block=UNLIMITED))
            / (bp / n_dev + prof_frame(rc, L, s)))


@lru_cache(maxsize=None)
def g4_lowest():
    """G4's lowest ratio per point: {point: (ratio, L, N_dev)} over its L and N_dev, the grid and every share."""
    return {p: min((own(rc, L, s, d), L, d) for L in range(l0, l1 + 1) for d in range(d0, d1 + 1) for rc in GRID
                   for s in (0.0,) + SHARES) for p, ((l0, l1), _, (d0, d1), _) in POINTS.items()}


def j2_cells():
    """The J2 energy cells (50 kchip/s, energy mode, lead reference _r against the profile _p) and the mask of those
    inside the R_agg, N_dev and k ranges of a working point."""
    span = pd.read_csv(RES / "cells_span.csv")
    m = paired(span[span.variant == f"full+back+burst+pilot128+pool=0.05+rc_cps=50000.0+span={UNLIMITED}"], cases()[0][3],
               "energy")
    m["ra"] = m.n_dev * m.r_target_kbps      # POINTS holds R_agg = N_dev x r_target (report.py Rules)
    inside = np.logical_or.reduce([m.ra.between(*r) & m.n_dev.between(*n) & m.k.between(*k)
                                   for _, r, n, k in POINTS.values()])
    wps = rp.read_working_points()
    assert [(r, n, k) for _, r, n, k in POINTS.values()] == [((w[2], w[3]), (w[6], w[7]), (w[4], w[5])) for w in wps]
    assert set(m.index[inside]) == set().union(*(rp.wp_cells(m, w)[0].index for w in wps)), "not report.wp_cells"
    return m, inside


def part_g4():
    """G4: per point
    the lowest ratio against the lead reference (one frame per allocation, as one frame of L per connection interval)
    over its payloads, its N_dev, the grid and every share, at the point's own N_dev with both beacons; and the J2 energy cells inside a point's R_agg, N_dev and k ranges."""
    sh = (0.0,) + SHARES
    lo = g4_lowest()
    ns = (1, 2, 3, 4, 8, 15, 39, 255, UNLIMITED)
    every = [min((back_cost(rc, 0, L, s, n, burst=True, block=UNLIMITED) / prof_cost(rc, L, s, n=n), L)
                 for L in range(1, 256) for rc in GRID for s in sh) for n in ns]
    m, inside = j2_cells()
    m["q"] = m.jpb_sensor_mean_r / m.jpb_sensor_mean_p
    w = m.loc[inside]
    inside = [m[m.ra.between(*r) & m.n_dev.between(*n) & m.k.between(*k)] for _, r, n, k in POINTS.values()]
    return ("\nG4. Typical value over the working points A-H : one"
            " frame of L octets per connection interval, i.e. one frame per allocation (n 1), the lead reference of G3"
            " (block unlimited), the ratio of the superframe capacities at the point's own N_dev with both beacons (own)."
            " Per point the lowest ratio over its payload and N_dev ranges, the grid and every share (at L, N_dev): "
            + ", ".join(f"{p} {v:.4f} ({L}, {d})" for p, (v, L, d) in lo.items())
            + f". Median over the eight points {np.median([v[0] for v in lo.values()]):.4f}, range"
            f" {min(lo.values())[0]:.4f} to {max(lo.values())[0]:.4f}.\n\nFloor over every payload of 1 to 255 octets,"
            " the grid and every share, by frames per allocation n, beacons left out but the profile's ack field"
            " (prof_cost) (lowest ratio, at payload): "
            + ", ".join(f"n {n if n < UNLIMITED else 'unlimited'} {v:.4f} ({L})" for n, (v, L) in zip(ns, every))
            + ".\n\nEnergy, 50 kchip/s energy mode (the J2 cells, lead reference / profile):"
            f" {len(m)} cells {m.q.min():.4f} to {m.q.max():.4f}, median {m.q.median():.4f}; the cells inside a"
            " point's R_agg, N_dev and k ranges: "
            + ", ".join(f"{p} {len(c)}" for p, c in zip(POINTS, inside) if len(c))
            + f", together {len(w)} cells {w.q.min():.4f} to {w.q.max():.4f}, median {w.q.median():.4f}\n")


# IEEE 802.15.6-2012 6.11.1.1, Equations (6) to (8), Table 24: GTn = pSIFS + pExtraIFS + mClockResolution + 2 SIn ppm,
# SIn = mNominalSynchInterval = 8 beacon periods, mHubClockPPMLimit 40 ppm; HBC pExtraIFS 10 us. The model
# books pSIFS only (PARAMETERS OI-24), so the guard charged below is GTn - pSIFS.
PEXTRA_US, MCLK_US, SYNCH_PERIODS, PPM = 10.0, 4.0, 8, 40e-6
# connection interval of each working point, s 
CI_S = {"A": (4, 10), "B": (0.5, 1.0), "C": (0.02, 0.05), "D": (0.0075, 0.01), "E": (0.0075, 0.0075),
        "F": (0.0075, 0.01), "G": (0.0075, 0.01), "H": (0.0075, 0.01)}


def gtn(period_us):
    """GTn - pSIFS (us) for a beacon period."""
    return PEXTRA_US + MCLK_US + 2 * SYNCH_PERIODS * PPM * period_us


def packed_gtn(rc, L, s, d):
    """GTn - pSIFS (us) in the capacity reading: the beacon period is the lead reference's superframe packed with its d
    allocations of one frame, T = beacon + d (allocation + (1 + s) gtn(T)), solved for T."""
    k = 2 * SYNCH_PERIODS * PPM
    a = back_cost(rc, 0, L, s, 1, burst=True, block=UNLIMITED)
    return gtn((beacons(rc, d)[0] + d * (a + (1 + s) * gtn(0))) / (1 - d * (1 + s) * k))


def own_gtn(rc, L, s, d, g=None):
    """own() with the guard g (GTn - pSIFS, us) in every scheduled allocation of both stacks and in the polled ones
    (6.11.1.3: s polled allocations per frame at the share s); g None: packed_gtn, the reference's (larger) guard
    charged to the profile as well."""
    br, bp = beacons(rc, d)
    g = (packed_gtn(rc, L, s, d) if g is None else g) * (1 + s)
    return (br / d + back_cost(rc, 0, L, s, 1, burst=True, block=UNLIMITED) + g) / (bp / d + prof_frame(rc, L, s) + g)


# The profile's allocation slot for the galvanic PHY, us: pAllocationSlotMin = pAllocationSlotResolution. The
# 500 us of Table 27 belongs to the clause 10 HBC PHY, which neither stack runs here.
SLOT_US = 50.0
SLOT_SHARES = (0.0,) + SHARES


def on_slots(t, slot):
    return math.ceil(t / slot - 1e-9) * slot


def own_slot(rc, L, s, d, g=None, slot=SLOT_US, beacon_slots=True):
    """own_gtn on the slot grid, one frame per allocation: frame + guard of each scheduled allocation of both
    stacks rounded up to whole slots, and with beacon_slots the beacon + pSIFS too, since allocations start on a slot
    boundary (5.7.7.2) and slot 0 starts at the beacon (5.3.1.2.1). Pool exchanges and polled guards stay continuous."""
    br, bp = beacons(rc, d)
    g = packed_gtn(rc, L, s, d) if g is None else g
    fr = back_cost(rc, 0, L, 0.0, 1, burst=True, block=UNLIMITED)
    fp = prof_frame(rc, L, 0.0)
    pool_r = back_cost(rc, 0, L, s, 1, burst=True, block=UNLIMITED) - fr
    pool_p = prof_frame(rc, L, s) - fp
    return (slot_cost(br, fr, pool_r, s, d, g, slot, beacon_slots)
            / slot_cost(bp, fp, pool_p, s, d, g, slot, beacon_slots))


def slot_cost(beacon, frame, pool, s, d, g, slot, beacon_slots):
    """Time per sensor and superframe of one stack on the slot grid, us (own_slot)."""
    return ((on_slots(beacon, slot) if beacon_slots else beacon) / d + on_slots(frame + g, slot) + pool + s * g)


def changes_spec(rc):
    """Steps 2 to 4 without step 1: the profile's framing on plain HBC's phy and pilot period (DRF mode, no RI field)."""
    s = spec(rc, *PROP)
    return dataclasses.replace(s, phy=spec(rc, *REF).phy, flags=dict(s.flags, opt=_NoOptions(5), pilot128=True))


def step_factors(rc, L, s, d, g=None, slot=SLOT_US, beacon_slots=True):
    """own_slot split by the ladder, (step 1, steps 2 to 4 on top, steps 2 to 4 first, step 1 on top); each pair
    multiplies to own_slot. Step 1 is plain HBC with the rate indicator on every frame, beacon and poll (ladder level 1,
    which has the 128-byte pilot period) and its acknowledgement otherwise unchanged; the guard is own_slot's on all."""
    g = packed_gtn(rc, L, s, d) if g is None else g
    ps, one, ch = pr.HBC["psifs_us"], spec(rc, 1, 0.05), changes_spec(rc)
    br, bp = beacons(rc, d)
    fr, fp = back_cost(rc, 0, L, 0.0, 1, burst=True, block=UNLIMITED), prof_frame(rc, L, 0.0)
    costs = [slot_cost(*x, s, d, g, slot, beacon_slots) for x in (
        (br, fr, back_cost(rc, 0, L, s, 1, burst=True, block=UNLIMITED) - fr),
        (pr.beacon_layout(one, d, 1)[-1]["t1_us"] + ps, pr.airtime_us(one, L) + ps + (pr.airtime_us(one, 3) + ps) / UNLIMITED,
         s * pr._pool_us(one, L)),
        (pr.beacon_layout(ch, d, 1)[-1]["t1_us"] + ps, pr._frame_us(ch, L), s * pr._pool_us(ch, L)),
        (bp, fp, prof_frame(rc, L, s) - fp))]
    r, a, c, p = costs
    return r / a, a / p, r / c, c / p


def part_g5():
    """G5: the G4
    working points with GTn charged to both stacks in the capacity reading and with the beacon period equal to the
    point's longest connection interval; the long-burst floor with a guard per allocation; the 'any frame' bound; the J2
    energy with the 6.11.1.2 beacon receive window (GTn - GT0 = 2 SIn ppm early per period) on both stacks."""
    sh = (0.0,) + SHARES
    zero, cap, ci = {}, {}, {}
    for p, ((l0, l1), _, (d0, d1), _) in POINTS.items():
        gci = gtn(CI_S[p][1] * 1e6)
        cells = [(rc, L, s, d) for L in range(l0, l1 + 1) for d in range(d0, d1 + 1) for rc in GRID for s in sh]
        zero[p] = min((own_gtn(rc, L, s, d, 0.0), L, d) for rc, L, s, d in cells)
        cap[p] = min((own_gtn(rc, L, s, d), L, d, rc, s) for rc, L, s, d in cells)
        ci[p] = min((own_gtn(rc, L, s, d, gci), L, d, rc, s) for rc, L, s, d in cells)
    g4 = g4_lowest()
    assert all(abs(zero[p][0] - g4[p][0]) < 1e-12 and zero[p][1:] == g4[p][1:] for p in POINTS), "zero guard is not G4"
    g10 = gtn(10e6)
    ns = (1, 2, 3, 4, 8, 15, 39, 255, UNLIMITED)
    floor = min(((back_cost(rc, 0, L, s, n, burst=True, block=UNLIMITED) + g10 / n)
                 / (prof_cost(rc, L, s, n=n) + g10 / n), L, n) for n in ns for L in range(1, 256) for rc in GRID for s in sh)
    anyf = []
    for lo_, hi_ in ((1, 64), (65, 255)):
        v = min((own_gtn(rc, L, s, d), L, d, rc) for L in range(lo_, hi_ + 1) for d in (1, 64) for rc in GRID for s in sh)
        # check at that payload and rate that no sensor count between 1 and 64 is lower
        assert all(own_gtn(v[3], v[1], s, d) >= v[0] - 1e-12 for d in range(1, 65) for s in sh), v
        anyf.append(v)
    m, inside = j2_cells()
    p_rx = pr.load_params()["energy"]["bcc_mW"]["typical"]["rx"] * 1e-3
    win = p_rx * 2 * SYNCH_PERIODS * PPM / (m.r_target_kbps * 1e3)     # J per bit, the window of every period
    q = (m.jpb_sensor_mean_r + win) / (m.jpb_sensor_mean_p + win)
    vals = lambda dct: [v[0] for v in dct.values()]                    # noqa: E731
    row = lambda dct: (", ".join(f"{p} {v[0]:.4f}" for p, v in dct.items()) + f"; median {np.median(vals(dct)):.4f},"  # noqa: E731
                       f" range {min(vals(dct)):.4f} to {max(vals(dct)):.4f}")
    return ("\nG5. The standard guard time on both stacks (IEEE 802.15.6-2012 6.11.1.1,"
            " Equations (6) to (8), Table 24, 6.11.1.3): GTn - pSIFS = pExtraIFS 10 us + mClockResolution 4 us + 2 x 8"
            " beacon periods x 40 ppm in every scheduled allocation of both stacks and in every polled allocation of the"
            " pool (own_gtn). Working points as G4 (lowest per point over its L and N_dev, the grid and every share):\n\n"
            f"- zero guard (self-check, equal to G4): {row(zero)}\n"
            f"- capacity reading, the beacon period the reference's packed superframe: {row(cap)}; minima at "
            + ", ".join(f"{p} L {v[1]} N {v[2]} {v[3] / 1e3:g} kchip/s share {v[4]}" for p, v in cap.items()) + "\n"
            f"- beacon period the point's longest connection interval (a sparse network): {row(ci)}\n\n"
            f"Long-burst floor with a guard of {g10:.0f} us per allocation (the 10 s beacon period) on both, every payload of"
            " 1 to 255 octets, the grid, every share, frames per allocation n " + ", ".join(
                str(n) if n < UNLIMITED else "unlimited" for n in ns)
            + f": {floor[0]:.6f} at {floor[1]} octets, n {floor[2] if floor[2] < UNLIMITED else 'unlimited'}"
            " (the guard per frame vanishes in an unlimited burst).\n\n'Any frame' bound, capacity reading, N_dev 1 and"
            " 64 (checked for 1 to 64 at the minimum), the grid, every share: "
            + "; ".join(f"payload {a}-{b} octets {v[0]:.4f} at {v[1]} octets, N_dev {v[2]}, {v[3] / 1e3:g} kchip/s"
                        for (a, b), v in zip(((1, 64), (65, 255)), anyf))
            + f".\n\nEnergy, the J2 cells with the beacon receive window on both stacks ({win.min() * 1e9:.2f} to"
            f" {win.max() * 1e9:.2f} nJ per bit at the typical receive power): {len(q)} cells {q.min():.4f} to {q.max():.4f},"
            f" median {q.median():.4f}; the {int(inside.sum())} cells inside a working point {q[inside].min():.4f} to"
            f" {q[inside].max():.4f}, median {q[inside].median():.4f}; by N_dev: "
            + ", ".join(f"{n} {int((m.n_dev == n).sum())} cells {q[m.n_dev == n].min():.4f} to {q[m.n_dev == n].max():.4f}"
                        f" (median {q[m.n_dev == n].median():.4f})" for n in sorted(m.n_dev.unique())) + "\n")


def _g6_point(p):
    """The slot-grid lowest ratios of one working point, (name, lowest + where) for each reading of part_g6."""
    (l0, l1), _, (d0, d1), _ = POINTS[p]
    cells = [(rc, L, s, d) for L in range(l0, l1 + 1) for d in range(d0, d1 + 1) for rc in GRID for s in SLOT_SHARES]
    gci = gtn(CI_S[p][1] * 1e6)
    readings = {"slots": {}, "air": dict(beacon_slots=False), "ci": dict(g=gci), "ci air": dict(g=gci, beacon_slots=False),
                "500": dict(slot=500.0), "500 air": dict(slot=500.0, beacon_slots=False)}
    return {k: min((own_slot(rc, L, s, d, **kw), L, d, rc, s) for rc, L, s, d in cells) for k, kw in readings.items()}


def _g6_any(L):
    """Lowest ratio of the payload L over 1 to 64 sensors, the grid and every share, {beacon in slots: (ratio, L, N_dev,
    rate, share)}. Exhaustive: on the slot grid the ratio is neither monotonic in N_dev nor lowest at no pool."""
    return {bs: min((own_slot(rc, L, s, d, beacon_slots=bs), L, d, rc, s) for d in range(1, 65) for rc in GRID
                    for s in SLOT_SHARES) for bs in (True, False)}


def _g7_scan(L):
    """Over the sensor counts, the grid and every share at the payload L: the smallest margin of step 1 over steps 2 to 4
    and the smallest steps 2 to 4, in the ladder's order, each with its cell."""
    fs = [(step_factors(rc, L, s, d), rc, s, d) for d in range(1, 65) for rc in GRID for s in SLOT_SHARES]
    return (min((f[0] / f[1], L, d, rc, s) for f, rc, s, d in fs), min((f[1], L, d, rc, s) for f, rc, s, d in fs))


def part_g7(cells):
    """G7: the G6 values split at the cell where each falls, on the slot and with no rounding at all."""
    rows = []
    for name, (v, L, d, rc, s) in cells.items():
        one, rest, ch, top = step_factors(rc, L, s, d)
        g = packed_gtn(rc, L, s, d)
        co, cr, _, _ = step_factors(rc, L, s, d, g, 1e-4, False)
        lo = [step_factors(r, L, s, d) for r in GRID]
        rows.append([name, f"{rc / 1e3:.1f}", L, d, s, f"{one:.4f}", f"{rest:.4f}", f"{one * rest:.4f}", f"{v:.4f}",
                     f"{ch:.4f} x {top:.4f}", f"{co:.4f} x {cr:.4f}", f"{min(x[0] for x in lo):.3f} to {max(x[0] for x in lo):.3f}",
                     f"{min(x[1] for x in lo):.3f} to {max(x[1] for x in lo):.3f}"])
    with Pool(8) as pool:
        scan = pool.map(_g7_scan, range(1, 65))
    bm, br = min(x[0] for x in scan), min(x[1] for x in scan)
    return ("\nG7. The G6 values split by the ladder, one tone: step 1 (plain HBC with the rate indicator on every frame,"
            " beacon and poll, step_factors) and steps 2 to 4 on top (on one tone step 4 changes nothing), the guard of the"
            " cell on all stacks, at the cell (chip rate, L, N_dev, share) where the G6 value falls. Rounded to the slot, the"
            " factors of a cell belong to that cell only: the minima of the total, step 1 and steps 2 to 4 fall at different"
            " rates\n\n"
            + md(["value", "chip rate (kchip/s)", "L (octets)", "N_dev", "share", "step 1", "steps 2 to 4", "product",
                  "G6 value", "other order: steps 2 to 4 x step 1", "no rounding: step 1 x steps 2 to 4",
                  "step 1 over the grid, same L, N_dev, share", "steps 2 to 4 over the grid"], rows)
            + "\nEvery payload of 1 to 64 octets, N_dev 1 to 64, the grid, every share: lowest step 1 / steps 2 to 4 "
            f"{bm[0]:.4f} (L {bm[1]}, N_dev {bm[2]}, {bm[3] / 1e3:.1f} kchip/s, share {bm[4]}); lowest steps 2 to 4 {br[0]:.4f}"
            f" (L {br[1]}, N_dev {br[2]}, {br[3] / 1e3:.1f} kchip/s, share {br[4]}).\n")


def part_g6():
    """G6: G5 on the allocation-slot grid, both stacks, beacon in slots and beacon as continuous air time; the
    packed superframe and the connection interval; the standard's 500 us for comparison; the 'any frame' bound."""
    with Pool(4) as pool:
        res = dict(zip(POINTS, pool.map(_g6_point, list(POINTS))))
    row = lambda k: (", ".join(f"{p} {res[p][k][0]:.4f}" for p in POINTS)  # noqa: E731
                     + f"; median {np.median([res[p][k][0] for p in POINTS]):.4f}")
    with Pool(8) as pool:
        by_l = pool.map(_g6_any, range(1, 256))
    anyf = {(bs, a): min(x[bs] for x in by_l[a - 1:b]) for bs in (True, False) for a, b in ((1, 64), (65, 255))}
    return (f"\nG6. G5 on an allocation slot of {SLOT_US:g} us, both stacks: frame + guard of each scheduled allocation"
            " rounded up to whole slots (own_slot); pool exchanges and polled guards continuous. Lowest per point over its L"
            " and N_dev, the grid and every share:\n\n"
            f"- capacity reading, beacon in slots: {row('slots')}; minima at "
            + ", ".join(f"{p} L {res[p]['slots'][1]} N {res[p]['slots'][2]} {res[p]['slots'][3] / 1e3:.1f} kchip/s share"
                        f" {res[p]['slots'][4]}" for p in POINTS) + "\n"
            f"- capacity reading, beacon continuous: {row('air')}\n"
            f"- beacon period the longest connection interval, beacon in slots: {row('ci')}\n"
            f"- beacon period the longest connection interval, beacon continuous: {row('ci air')}\n"
            f"- the standard's 500 us, beacon in slots: {row('500')}\n"
            f"- the standard's 500 us, beacon continuous: {row('500 air')}\n\n"
            "'Any frame' bound, capacity reading, N_dev 1 to 64, the grid, every share, payload 1-64 / 65-255 octets: "
            + "; ".join(f"beacon {'in slots' if bs else 'continuous'} {anyf[bs, 1][0]:.4f} at {anyf[bs, 1][1]} octets, N_dev"
                        f" {anyf[bs, 1][2]}, {anyf[bs, 1][3] / 1e3:.1f} kchip/s share {anyf[bs, 1][4]} / {anyf[bs, 65][0]:.4f} at {anyf[bs, 65][1]}"
                        f" octets, N_dev {anyf[bs, 65][2]}, {anyf[bs, 65][3] / 1e3:.1f} kchip/s share {anyf[bs, 65][4]}" for bs in (True, False))
            + ".\n\nThese hold SLOT_US fixed over the 461 rates and the whole of a packed superframe in one beacon period;"
            " The fields and a denser rate scan were checked separately.\n"
            + part_g7({**{f"point {p}": res[p]["slots"] for p in POINTS}, "8 octets": by_l[7][True],
                       "1-64 octets": anyf[True, 1], "65-255 octets": anyf[True, 65], "255 octets": by_l[254][True]}))


def at_capacity(task):
    sid, rc, lmax, st, n_dev = task
    # (level, share), + block ack, burst, block, '128 B' pilot period
    lv, pool, back, burst, span, p128 = (*st, False, False, None, False)[:6]
    s = spec(rc, lv, pool, lmax, sid, back, burst, span, p128)
    cap = rp.capacity(sid, rc, n_dev, lv, pool, None, lmax, back, burst, span, p128)
    bound = 8000 * s.max_payload_octets / rp.frame_cost(s, s.max_payload_octets)
    psr, periods = [], []
    for mode in ("energy", "latency"):
        sel = pr.select_period(cap / n_dev, mode, n_dev, s)
        assert sel, (task, cap)
        periods.append(sel["period_us"] / 1e3)
        for seed in (0, 1):
            psr += [r["psr"] for r in en.Sim(s, sel, 1, n_dev, "collision", max(2e6, 3 * sel["period_us"]), seed)
                    .run().rows(dict(network_id=0))]
    above = any(pr.select_period((cap + 0.5) / n_dev, m, n_dev, s) for m in ("energy", "latency"))
    return dict(sid=sid, rc=rc, lmax=s.max_payload_octets, stack="profile" if lv else "reference", n_dev=n_dev,
                cap=cap, bound=bound, psr_min=min(psr), T_ms=periods, above=above)


def check_pool_bound():
    """capacity() scans down from the pool bound 8 L / frame_cost(L); for the short frames as for 255 octets
    (report.self_check) the one short frame of a period saves less than the beacon's share per sensor."""
    for rc in (500e3, 2.625e6):
        for lmax in (8, 64):
            for lv, pool in (REF, PROP):
                s, n = spec(rc, lv, pool, lmax), 10
                save = max(o * rp.frame_cost(s, lmax) / lmax - rp.frame_cost(s, o) for o in range(1, lmax + 1))
                bcn = math.ceil(pr.beacon_layout(s, n, 1)[-1]["t1_us"] - 1e-6) + pr.HBC["psifs_us"]
                assert save < bcn / n, (rc, lmax, lv, save, bcn / n)


def part_d(workers):
    check_pool_bound()
    tasks = [(sid, rc, lmax, st, n) for sid, rc, lmax in CAP_CASES for st in (REF, PROP) for n in (1, 4)]
    with Pool(workers) as p:
        res = p.map(at_capacity, tasks)
    rows = []
    for i in range(0, len(res), 4):
        r, q = res[i:i + 2], res[i + 2:i + 4]    # reference N_dev 1, 4; profile N_dev 1, 4
        a = r + q
        rows.append([a[0]["sid"].replace("bcc-", ""), f"{a[0]['rc'] / 1e3:g}", a[0]["lmax"],
                     f"{r[0]['cap']:g} / {q[0]['cap']:g}", f"{r[1]['cap']:g} / {q[1]['cap']:g}",
                     f"{q[0]['cap'] / r[0]['cap']:.3f}", f"{q[0]['bound'] / r[0]['bound']:.3f}",
                     f"{min(x['psr_min'] for x in a):.4f}", "none" if not any(x["above"] for x in a) else "FEASIBLE"])
    return ("D. Simulated at the admitted capacity (k 1, collision, 2 seeds, both ci modes)\n\n"
            + md(["tones", "chip rate (kchip/s)", "payload limit (octets)", "capacity N_dev 1, ref / profile (kb/s)",
                  "capacity N_dev 4, ref / profile (kb/s)", "ratio N_dev 1", "closed-form ratio", "lowest psr",
                  "schedule 0.5 kb/s above"], rows))


def compact(raw, back_only=False):
    src = Path(raw) / "t_l091"
    for name in ("back", "burst", "span"):
        b = Path(raw) / name / f"{name}_cells.csv"
        if b.exists():
            pd.read_csv(b)[KEEP].to_csv(RES / f"cells_{name}.csv", index=False, float_format="%.6g",
                                        lineterminator="\n")
            print(f"claim4: {RES / f'cells_{name}.csv'} rewritten from {b}")
    if back_only or not (src / "lossfree_cells.csv").exists():
        return
    d = pd.concat([pd.read_csv(f) for f in (src / "lossfree_cells.csv", src / "power_cells.csv") if f.exists()])
    d[KEEP].to_csv(RES / "cells_lossfree.csv", index=False, float_format="%.6g", lineterminator="\n")
    print(f"claim4: {RES / 'cells_lossfree.csv'} rewritten from {src}")


def cases():
    """(label, chip rate, cells of the reference, cells of the profile), collision model, one tone unless named."""
    lf = pd.read_csv(RES / "cells_lossfree.csv")
    lad = pd.read_csv(RES / "cells_ladder.csv", dtype={"option": str}, keep_default_na=False, na_values=[""])
    rat = pd.read_csv(RES / "cells_rates.csv", dtype={"option": str}, keep_default_na=False, na_values=[""])
    lad = lad[(lad.sweep == "ladder") & (lad.reception == "collision")]
    rat = rat[(rat.phy == "sf8") & (rat.reception == "collision")]
    one, six = lad[lad.phy == "sf8"], lad[lad.phy == "multitone"]
    v = lf.variant
    return [("50 kchip/s", 50e3, lf[v == "full+pool=0.05+rc_cps=50000.0"], lf[v == "full+fill+opt5+rc_cps=50000.0"]),
            ("500 kchip/s", 500e3, one[one.option == "pool=0.05"], one[one.option == "opt5"]),
            ("2.625 Mchip/s", 2.625e6, rat[rat.option == "pool=0.05+rc_cps=2625000.0"],
             rat[rat.option == "opt5+rc_cps=2625000.0"]),
            ("500 kchip/s, six tones", None, six[six.option == "pool=0.05"], six[six.option == "opt5"]),
            ("50 kchip/s, low power", None, lf[v == "full+pool=0.05+power=low+rc_cps=50000.0"],
             lf[v == "full+fill+opt5+power=low+rc_cps=50000.0"]),
            ("50 kchip/s, high power", None, lf[v == "full+pool=0.05+power=high+rc_cps=50000.0"],
             lf[v == "full+fill+opt5+power=high+rc_cps=50000.0"])]


def paired(a, b, mode):
    a, b = a[(a.ci_mode == mode) & (a.feasible == 1)], b[(b.ci_mode == mode) & (b.feasible == 1)]
    return a.merge(b, on=CELL, suffixes=("_r", "_p"))


EHEAD = ["case", "ci mode", "cells", "min", "median of ratios", "max", "ratio of medians", "medians ref / profile (nJ)",
         "equal-period cells: n, range", "hub included: min / median"]


def erows(pairs):
    """E's rows for (label, reference cells, profile cells), and the trial-limit floor of the 50 kchip/s energy mode."""
    rows, extra = [], []
    for label, a, b in pairs:
        for mode in ("energy", "latency"):
            m = paired(a, b, mode)
            if m.empty:
                continue
            q = m.jpb_sensor_mean_r / m.jpb_sensor_mean_p
            h = m.jpb_ban_mean_r / m.jpb_ban_mean_p
            e = q[m.period_us_r == m.period_us_p]
            rows.append([label, mode, len(m), f"{q.min():.3f}", f"{q.median():.3f}", f"{q.max():.3f}",
                         f"{m.jpb_sensor_mean_r.median() / m.jpb_sensor_mean_p.median():.3f}",
                         f"{m.jpb_sensor_mean_r.median() * 1e9:.4g} / {m.jpb_sensor_mean_p.median() * 1e9:.4g}",
                         f"{len(e)}: {e.min():.3f} to {e.max():.3f}" if len(e) else "0",
                         f"{h.min():.3f} / {h.median():.3f}"])
            if label == "50 kchip/s" and mode == "energy":
                w = m.jpb_sensor_lo_r / m.jpb_sensor_hi_p
                extra.append(f"50 kchip/s energy mode: lowest per-cell ratio of the trial limits (reference lower 95 % "
                             f"limit / profile upper limit) {w.min():.3f}.")
    return rows, extra


def part_e():
    rows, extra = erows([(label, a, b) for label, _, a, b in cases()])
    # the filling pool against the fixed pool, 50 kchip/s (cells_rates.csv carries L5 with the 5 % pool at k 1 and 8)
    rat = pd.read_csv(RES / "cells_rates.csv", dtype={"option": str}, keep_default_na=False, na_values=[""])
    fixed = rat[(rat.option == "opt5+rc_cps=50000.0") & (rat.phy == "sf8") & (rat.reception == "collision")]
    fill = cases()[0][3]
    for mode in ("energy", "latency"):
        m = paired(fixed, fill, mode)
        d = (m.jpb_sensor_mean_r / m.jpb_sensor_mean_p - 1).abs().max()
        extra.append(f"50 kchip/s {mode} mode: profile with the 5 % pool against the filling pool, {len(m)} cells, "
                     f"largest relative difference of J per bit {d:.4f}.")
    return ("E. Energy per delivered bit, reference / profile, paired per cell (k, N_dev, R), collision model\n\n"
            + md(EHEAD, rows)
            + "\n" + "\n".join(extra) + "\n")


def jpb(s, sel, n_dev, dur):
    return float(np.mean([np.mean([r["jpb_sensor"] for r in en.Sim(s, sel, 1, n_dev, "collision", dur, seed).run()
                                   .rows(dict(network_id=0))]) for seed in (0, 1)]))


def rerun(task):
    """Reference at its own period, the profile at its own and at the reference's period (protocols._fit)."""
    rc, r, n_dev, mode, power = task
    f = {"power": power} if power else {}
    ref = pr.get_spec("bcc", "full", "sf8", rc_cps=rc, pool=0.05, **f)
    pro = pr.get_spec("bcc", "full", "sf8", rc_cps=rc, opt=5, pool=0.05, **f)
    a, b = pr.select_period(r, mode, n_dev, ref), pr.select_period(r, mode, n_dev, pro)
    c = dict(pr._fit(a["payload_bits"], a["period_us"], n_dev, pro), deadline_us=a["deadline_us"])
    dur = max(2e6 if rc >= 5e5 else 15e6, 3 * a["period_us"])
    ja, jb, jc = jpb(ref, a, n_dev, dur), jpb(pro, b, n_dev, dur), jpb(pro, c, n_dev, dur)
    return dict(rc=rc, r=r, n_dev=n_dev, mode=mode, power=power or "typical", Tr=a["period_us"] / 1e3,
                Tp=b["period_us"] / 1e3, own=ja / jb, equal=ja / jc)


def part_f(workers):
    tasks = []
    for _, rc, a, b in cases()[:3]:
        for mode in ("energy", "latency"):
            m = paired(a, b, mode)
            tasks += [(rc, r, n, mode, None) for r, n in sorted(set(zip(m.r_target_kbps, m.n_dev)))]
    pw = [(rc, r, n, "energy", p) for rc in (500e3, 2.625e6) for p in ("low", "high")
          for r, n in ((1.0, 1), (1.0, 4), (10.0, 4))]
    with Pool(workers) as p:
        res = p.map(rerun, tasks + pw)
    d = pd.DataFrame(res)
    rows = []
    for (rc, mode, power), g in d.groupby(["rc", "mode", "power"], sort=False):
        u = g[g.Tr != g.Tp]
        rows.append([f"{rc / 1e3:g}", mode, power, len(g), f"{g.own.min():.3f} to {g.own.max():.3f}",
                     f"{g.equal.min():.3f} to {g.equal.max():.3f}", f"{g.equal.median():.3f}",
                     ", ".join(f"{x.n_dev} x {x.r:g} kb/s: T {x.Tr:g} / {x.Tp:g} ms, {x.own:.3f} -> {x.equal:.3f}"
                               for x in u.itertuples()) if mode == "energy" else f"{len(u)} of {len(g)}"])
    return ("F. Engine reruns (k 1, collision, 2 seeds): J per bit reference / profile at each stack's own period and with"
            " the profile at the reference's period\n\n"
            + md(["chip rate (kchip/s)", "ci mode", "power", "(R, N_dev) points", "own periods", "reference's period",
                  "median at the reference's period", "points with unequal periods"], rows))


def part_h(workers, burst=False):
    """H: the block-ack reference (flag back); burst: I, the lead reference, the same in the HBC burst mode."""
    check_pool_bound()
    back = (0, 0.05, True, burst)
    tasks = [(sid, rc, lmax, st, n) for sid, rc, lmax in CAP_CASES for st in (back, PROP) for n in (1, 4)]
    with Pool(workers) as p:
        res = p.map(at_capacity, tasks)
    rows = []
    for i in range(0, len(res), 4):
        r, q = res[i:i + 2], res[i + 2:i + 4]    # block-ack reference N_dev 1, 4; profile N_dev 1, 4
        a = r + q
        rows.append([a[0]["sid"].replace("bcc-", ""), f"{a[0]['rc'] / 1e3:g}", a[0]["lmax"],
                     f"{r[0]['cap']:g} / {q[0]['cap']:g}", f"{r[1]['cap']:g} / {q[1]['cap']:g}",
                     f"{q[0]['cap'] / r[0]['cap']:.3f}", f"{q[1]['cap'] / r[1]['cap']:.3f}",
                     f"{q[0]['bound'] / r[0]['bound']:.3f}", f"{min(x['psr_min'] for x in a):.4f}",
                     "none" if not any(x["above"] for x in a) else "FEASIBLE"])
    name = "burst" if burst else "back"
    cb = pd.read_csv(RES / f"cells_{name}.csv")
    c = {x[0]: x for x in cases()}
    v = lambda rc, d=cb, f=name: d[d.variant == f"full+back{'+burst' if f == 'burst' else ''}+pool=0.05+rc_cps={rc}"]  # noqa: E731
    pairs = [(lab, v(rc), c[lab][2], c[lab][3]) for lab, rc in (("50 kchip/s", 50000.0), ("500 kchip/s", 500000.0),
                                                               ("2.625 Mchip/s", 2625000.0))]
    e, extra = erows([(lab, a, b) for lab, a, _, b in pairs])
    if burst:   # per cell the reference that spends least energy per bit: burst, block ack without burst, I-Ack
        nb = pd.read_csv(RES / "cells_back.csv")
        for (lab, a, i, b), rc in zip(pairs, (50000.0, 500000.0, 2625000.0)):
            m = paired(a, b, "energy")
            for d, sfx in ((v(rc, nb, "back"), "_b"), (i, "_i")):
                m = m.merge(paired(d, b, "energy")[CELL + ["jpb_sensor_mean_r"]], on=CELL, suffixes=("", sfx))
            best = m[["jpb_sensor_mean_r", "jpb_sensor_mean_r_b", "jpb_sensor_mean_r_i"]].min(axis=1)
            q = best / m.jpb_sensor_mean_p
            extra.append(f"{lab} energy mode, against the best of the three references in each cell ({len(m)} cells; "
                         f"burst best in {int((m.jpb_sensor_mean_r == best).sum())}): {q.min():.3f} to {q.max():.3f}, "
                         f"median {q.median():.3f}.")
    else:
        for lab, a, i, b in pairs:     # per cell the reference, block ack or I-Ack, that spends less energy per bit
            m = paired(a, b, "energy").merge(paired(i, b, "energy")[CELL + ["jpb_sensor_mean_r"]], on=CELL,
                                             suffixes=("", "_i"))
            q = np.minimum(m.jpb_sensor_mean_r, m.jpb_sensor_mean_r_i) / m.jpb_sensor_mean_p
            extra.append(f"{lab} energy mode, against the better of the two references in each cell ({len(m)} cells, "
                         f"the block-ack one in {int((m.jpb_sensor_mean_r < m.jpb_sensor_mean_r_i).sum())}): "
                         f"{q.min():.3f} to {q.max():.3f}, median {q.median():.3f}.")
    if burst:
        title = ("I. The lead reference, the block-ack reference of H in the HBC burst mode, in the engine (protocols"
                 " flags back and burst: the sensor's first frame of its allocation in DRF mode, the others pMIFS 20 us"
                 " apart in RI mode, the B-Ack pSIFS after the last)\n\nI1. Simulated at the"
                 " admitted capacity (k 1, collision, 2 seeds, both ci modes); the closed form is 8 L /"
                 " report.frame_cost(L), for the burst reference G2's unlimited burst\n\n")
        e_title = ("\nI2. Energy per delivered bit, burst reference / profile, paired per cell (k, N_dev, R) as E,"
                   " collision model (results/cells_burst.csv, run_sweeps.sh burst, against the profile cells of E)\n\n")
    else:
        title = ("H. The reference with the standard's block acknowledgement in the engine (protocols flag back: L-Ack"
                 " data frames, one B-Ack frame of 3 octets pSIFS after the sensor's last frame of its allocation, the"
                 " pool exchanges with their I-Ack)\n\nH1. Simulated at the admitted capacity (k 1,"
                 " collision, 2 seeds, both ci modes); the closed form is 8 L / report.frame_cost(L), for the block-ack"
                 " reference G's unlimited block\n\n")
        e_title = ("\nH2. Energy per delivered bit, block-ack reference / profile, paired per cell (k, N_dev, R) as E,"
                   " collision model (results/cells_back.csv, run_sweeps.sh back, against the profile cells of E)\n\n")
    return (title
            + md(["tones", "chip rate (kchip/s)", "payload limit (octets)",
                  "capacity N_dev 1, block-ack ref / profile (kb/s)", "capacity N_dev 4, block-ack ref / profile (kb/s)",
                  "ratio N_dev 1", "ratio N_dev 4", "closed-form ratio", "lowest psr", "schedule 0.5 kb/s above"], rows)
            + e_title + md(EHEAD, e) + "\n" + "\n".join(extra) + "\n")


def part_j(workers):
    """J: the lead reference, the burst reference of I with an unlimited block spanning
    allocations (flags back, burst, span): J1 at its admitted capacity beside the profile and I1's burst reference
    (capacity scans; the profile's packet success is I1's), J2 energy as I2, also against the best of the four
    references per cell."""
    check_pool_bound()
    lead = (0, 0.05, True, True, UNLIMITED, True)     # '128 B' pilot period
    with Pool(workers) as p:
        res = p.map(at_capacity, [(sid, rc, lmax, lead, n) for sid, rc, lmax in CAP_CASES for n in (1, 4)])
    rows = []
    for (sid, rc, lmax), i in zip(CAP_CASES, range(0, len(res), 2)):
        r, lm = res[i:i + 2], res[i]["lmax"]
        q = [rp.capacity(sid, rc, n, *PROP, None, lmax) for n in (1, 4)]
        b = [rp.capacity(sid, rc, n, *REF, None, lmax, True, True) for n in (1, 4)]
        rows.append([sid.replace("bcc-", ""), f"{rc / 1e3:g}", lm, f"{r[0]['cap']:g} / {q[0]:g}",
                     f"{r[1]['cap']:g} / {q[1]:g}", f"{b[0]:g} / {b[1]:g}", f"{q[0] / r[0]['cap']:.3f}",
                     f"{q[1] / r[1]['cap']:.3f}", f"{min(x['psr_min'] for x in r):.4f}",
                     "none" if not any(x["above"] for x in r) else "FEASIBLE"])
    cs, cb, nb = (pd.read_csv(RES / f"cells_{x}.csv") for x in ("span", "burst", "back"))
    c = {x[0]: x for x in cases()}
    v = lambda d, rc, f, s="": d[d.variant == f"full+back{f}+pool=0.05+rc_cps={rc}{s}"]     # noqa: E731
    lead = lambda rc: v(cs, rc, "+burst+pilot128", f"+span={UNLIMITED}")    # noqa: E731
    rcs = (("50 kchip/s", 50000.0), ("500 kchip/s", 500000.0), ("2.625 Mchip/s", 2625000.0))
    e, extra = erows([(lab, lead(rc), c[lab][3]) for lab, rc in rcs])
    for lab, rc in rcs:     # per cell the reference that spends least energy per bit: span, burst, block ack, I-Ack
        b = c[lab][3]
        m = paired(lead(rc), b, "energy")
        for d, sfx in ((v(cb, rc, "+burst"), "_u"), (v(nb, rc, ""), "_b"), (c[lab][2], "_i")):
            m = m.merge(paired(d, b, "energy")[CELL + ["jpb_sensor_mean_r"]], on=CELL, suffixes=("", sfx))
        best = m[["jpb_sensor_mean_r", "jpb_sensor_mean_r_u", "jpb_sensor_mean_r_b", "jpb_sensor_mean_r_i"]].min(axis=1)
        q = best / m.jpb_sensor_mean_p
        lo = m.loc[q.idxmin()]
        extra.append(f"{lab} energy mode, against the best of the four references in each cell ({len(m)} cells; span "
                     f"best in {int((m.jpb_sensor_mean_r == best).sum())}): {q.min():.4f} to {q.max():.4f}, median "
                     f"{q.median():.4f}; lowest at N_dev {lo.n_dev:g} x {lo.r_target_kbps:g} kb/s; span alone / "
                     f"burst alone {(m.jpb_sensor_mean_r / m.jpb_sensor_mean_r_u).min():.4f} to "
                     f"{(m.jpb_sensor_mean_r / m.jpb_sensor_mean_r_u).max():.4f}.")
    return ("J. The lead reference: the burst reference of I with one block spanning all its"
            " allocations (protocols flags back, burst, span unlimited: no B-Ack, the pSIFS after the burst stays;"
            " 6.2.9.4, 5.4.2.4), the pool exchanges with their I-Ack\n\nJ1. Simulated at the admitted capacity (k 1,"
            " collision, 2 seeds, both ci modes; the profile's capacities by report.capacity, its packet success in"
            " I1)\n\n"
            + md(["tones", "chip rate (kchip/s)", "payload limit (octets)", "capacity N_dev 1, lead ref / profile (kb/s)",
                  "capacity N_dev 4, lead ref / profile (kb/s)", "I1 burst ref N_dev 1 / 4 (kb/s)", "ratio N_dev 1",
                  "ratio N_dev 4", "lowest psr, lead ref", "schedule 0.5 kb/s above"], rows)
            + "\nJ2. Energy per delivered bit, lead reference / profile, paired per cell (k, N_dev, R) as E, collision"
            " model (results/cells_span.csv, run_sweeps.sh span, against the profile cells of E)\n\n"
            + md(EHEAD, e) + "\n" + "\n".join(extra) + "\n")


def replace_part(text, key, new):
    """text with its part `key` ('G', 'H') replaced by new, or new appended."""
    parts = re.split(r"\n(?=[A-Z]\. )", text)
    i = next((i for i, x in enumerate(parts) if x.startswith(f"{key}. ")), None)
    parts[i if i is not None else len(parts):] = [new] + (parts[i + 1:] if i is not None else [])
    return "\n".join(parts)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--workers", type=int, default=4, help="worker processes (at most 4)")
    ap.add_argument("--raw", default=RES / "raw", type=Path)
    ap.add_argument("--only-g", action="store_true", help="recompute part G only and replace it in t15_claim4.md")
    ap.add_argument("--only-h", action="store_true", help="recompute part H only and replace it in t15_claim4.md")
    ap.add_argument("--only-i", action="store_true", help="recompute part I (the burst reference) only")
    ap.add_argument("--only-j", action="store_true", help="recompute part J (the lead reference, block spanning "
                    "allocations) only")
    ap.add_argument("--closed", action="store_true", help="recompute the parts without an engine run (A, B, C, E, G) "
                    "and replace them in t15_claim4.md")
    a = ap.parse_args()
    t0, w = time.time(), max(1, min(a.workers, 4))
    out = RES / "tables" / "t15_claim4.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    if not out.exists():
        out.write_text("# T15: the galvanic profile against the pooled reference (claim4.py)\n\n")
    if a.only_g:
        text = replace_part(out.read_text(), "G", part_g())
    elif a.closed:
        text = out.read_text()
        for key, f in (("A", part_a), ("B", part_b), ("C", part_c), ("E", part_e), ("G", part_g)):
            text = replace_part(text, key, f())
    elif a.only_j:
        compact(a.raw, back_only=True)
        text = replace_part(out.read_text(), "J", part_j(w))
    elif a.only_h or a.only_i:
        compact(a.raw, back_only=True)
        text = replace_part(out.read_text(), "I" if a.only_i else "H", part_h(w, a.only_i))
    else:
        compact(a.raw)
        parts = [part_a(), part_b(), part_c(), part_e(), part_d(w), part_f(w), part_g(), part_h(w), part_h(w, True), part_j(w)]
        text = ("# T15: the galvanic profile against the pooled reference (claim4.py)\n\n"
                + "\n".join(parts))
    out.write_text(text)
    print(text)
    print(f"claim4: {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()

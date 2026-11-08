"""Sweep driver of the coexistence simulator: engine.simulate_trial over a (protocol, R, k, N_dev) grid.

Purpose: run every selected (protocol, variant, PHY) x reception x ci_mode x R x N_dev x k cell for a number of trials
in parallel (at most 8 worker processes, default 4), append one row per BAN and trial to --out, and write per-cell aggregates to
--cells. A trial runs for --duration-s or 40 periods of the cell (MIN_PERIODS), whichever is longer (energy-mode periods
reach 4 s): a frame still queued at the end counts as not delivered, so the run must be long against the retry and
queue horizon.
Resumable: a (cell, trial, seed) with all its k rows in --out is skipped; a partial last line left by a crash is cut
off at start, and each trial's rows go out in one write. The seed of a trial depends on (--seed, R, k, N_dev,
ci_mode, trial) only, and the engine draws the room, the link gains, the superframe and sample-clock phases (as
fractions of the period) and the clock errors from substreams of it, so every protocol and variant of a cell sees the
same ones; the other draws (channels, backoffs, per-packet fades) are the protocol's own.
The catalogue runs every rate or PHY of each stack at its full capability (protocols.catalogue). BCC runs three
entries: full on one tone (phy 'sf8': the IEEE 802.15.6 HBC stack on the galvanic channel), full on the N tones of the
 plan (phy 'multitone': each frame striped over the tones, one BCC full packet per tone, I-Ack on every tone, N
front ends' power) and the as-built stack as coded (variant 'thesis', phy 'bfsk'); --bcc-rc and
--bcc-power reach all three. An entry protocols.get_spec refuses is skipped with a printed line: the multi-tone entry
under --bcc-contact and --bcc-coupling (no per-tone contact chain or coupling) and at a --bcc-rc without a tone (above
3.38e6), every BCC entry above 9e6, both BCC full
entries under --bcc-fcs (they carry the 802.15.6 FCS; the flag is the as-built stack's).
--bcc-opt N runs both BCC full entries at level N of the galvanic-optimisation ladder (DESIGN 4.3, variant label
'full+optN'; levels 1 and 3 802.15.6 options, 2, 4 and 5 proposals); the as-built stack keeps no flag.
Inputs:  the engine modules; an existing --out CSV when resuming.
Outputs: --out CSV (appended, header written once), one row per BAN and trial, columns engine.COLUMNS:
           protocol, variant (full | thesis, plus '+flag' for non-default flags, e.g. 'thesis+contact', 'full+opt3',
           and '+nocifloor' on the sinr rows with --ci-floor off), phy (bcc: sf8, multitone or bfsk),
           reception (collision | sinr), strip_hopping (always 0: the v1 schema), ci_mode (latency | energy),
           period_us, r_target_kbps,
           r_actual_kbps, r_agg_kbps (= n_dev x r_actual), n_dev, k, payload_bits (L per period), frames_per_period,
           duration_us, trial, seed, network_id (BAN 0..k-1), feasible (0: no period fits, metrics empty);
           frame counts: generated (frames whose first access opportunity falls within the run), tx_attempts (data
           transmissions, and TSCH cells whose CCA found the channel busy: a failed attempt, PARAMETERS 3.3),
           delivered_ok, delivered_corrupt (residual errors, no payload CRC: BCC), drop_collision, drop_noise,
           drop_sync (final attempt lost to overlap / noise / sync; frames sent once), drop_csma (channel access
           failure), drop_retry (retry limit), drop_queue (arrival to a full queue), drop_end (still queued at the
           end: waiting for a retry or behind a backlog); the outcome counts add up to generated;
           psr = delivered_ok / generated, pdr = (delivered_ok + delivered_corrupt) / generated, per_first (failed
           first attempts / first attempts), retx_per_frame (retransmissions / first attempts), lat_mean_ms,
           lat_p50_ms, lat_p95_ms, lat_p99_ms, lat_max_ms (end of the delivering frame minus its sample instant),
           on_time (delivered_ok within the cell's reference period, the same deadline for every protocol /
           generated), raw_ber (payload bit errors before FEC / payload bits on air; Walsh stacks: demapped bits,
           so raw = post), post_ber (payload info-bit errors after FEC / info bits), goodput_kbps
           (delivered_ok application bits / duration, padding excluded), airtime_frac (the BAN's air time /
           duration), e_sensor_j (mean sensor radio energy), e_hub_j, p_sensor_uw, jpb_sensor and jpb_ban (J per
           delivered application bit, sensors alone and sensors + hub; inf when the BAN delivered nothing),
           lifetime_days (CR2032, mean sensor current), collision_frac (data attempts overlapping another
           same-channel transmission), cca_busy_frac; per trial, repeated on each BAN row: jain_fairness (of psr
           over the k BANs), channel_occupancy (mean busy fraction over channels);
           venue, density, case (room and channel case of the SINR model).
         --cells CSV (rewritten from the whole --out at the end), one row per cell: protocol, variant, phy,
           reception, strip_hopping, ci_mode, r_target_kbps, n_dev, k, duration_us, venue, density, case,
           period_us, r_actual_kbps, payload_bits, frames_per_period, feasible, n_trials, and for each of psr, pdr,
           on_time, goodput_kbps, lat_p50_ms, lat_p95_ms, lat_p99_ms, jpb_sensor, jpb_ban, lifetime_days, raw_ber,
           post_ber, collision_frac, jain_fairness: <m>_mean, <m>_lo, <m>_hi (mean over trials of the per-trial
           value, Student-t 95 % interval over trials, clipped to [0, 1] for the fractions (FRAC) and to >= 0
           for the others; empty with one trial or an inf trial). The per-trial value
           is the mean over the k BANs, except J per bit: the ratio of sums over the k BANs, sum(e_sensor_j x n_dev)
           (+ e_hub_j for jpb_ban) / sum(goodput_kbps x duration_us / 1000) bits, so a BAN that delivered nothing
           counts; inf when the trial delivered nothing.
argv:    see --help; example: python sweep.py --protocols ble,ieee802156nb --k 1,8 --n-dev 4 --r 10,100 --trials 2
"""
from __future__ import annotations
import argparse
import csv
import io
import os
import sys
import time
import zlib
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import engine as en            # noqa: E402
import protocols as pr         # noqa: E402

R_COARSE = "1,5,10,30,50,100,200,500,1000"    # BLE_SIM.py coarse grid (kb/s per link)
# periods per trial at least: frames awaiting a retry at the end (drop_end) bias psr down by about per_first / 40
MIN_PERIODS = 40
CELL = ["protocol", "variant", "phy", "reception", "strip_hopping", "ci_mode", "r_target_kbps", "n_dev", "k",
        "duration_us", "venue", "density", "case"]
AGG = ["psr", "pdr", "on_time", "goodput_kbps", "lat_p50_ms", "lat_p95_ms", "lat_p99_ms", "jpb_sensor", "jpb_ban",
       "lifetime_days", "raw_ber", "post_ber", "collision_frac", "jain_fairness"]
FRAC = {"psr", "pdr", "on_time", "raw_ber", "post_ber", "collision_frac", "jain_fairness"}   # intervals within [0, 1]


def _floats(s):
    return [float(x) for x in s.split(",")]


def _ints(s):
    return [int(x) for x in s.split(",")]


def specs(args):
    """(pid, variant, phy, flags) of the run: catalogue entries of the chosen protocols and variants, the PHYs replaced
    by the --phy entries of that protocol when there are any; ble4thesis runs '1M' and '2M'; --bcc-fcs runs the
    as-built BCC only (BCC full carries the HBC FCS)."""
    cat = pr.catalogue() + [("ble4thesis", "thesis", p) for p in ("1M", "2M")]
    pids = args.protocols.split(",") if args.protocols else sorted({c[0] for c in pr.catalogue()})
    out = []
    for pid, var, ph in cat:
        if pid not in pids or var not in args.variants.split(","):
            continue
        chosen = [x.split("=", 1)[1] for x in args.phy.split(",") if x.split("=", 1)[0] == pid] if args.phy else []
        for p in chosen or [ph]:
            flags = {}
            if pid == "bcc" and var == "full" and args.bcc_fcs:
                continue    # the FCS case is the as-built stack's: BCC full carries the HBC FCS already
            if pid == "bcc":    # only non-default flags, so the variant label stays 'thesis' for the reference
                flags = {k: v for k, v in (("contact", args.bcc_contact), ("fcs", args.bcc_fcs),
                                           ("coupling_db", args.bcc_coupling), ("rc_cps", args.bcc_rc),
                                           ("power", args.bcc_power), ("opt", args.bcc_opt if var == "full" else 0),
                                           ("pool", args.bcc_pool if var == "full" else None),
                                           ("fill", args.bcc_pool_fill and var == "full"),
                                           ("touch", args.bcc_touch if var == "full" else None),
                                           ("touch_db", args.bcc_touch_db if var == "full" else None),
                                           ("fallback", args.bcc_fallback and var == "full"),
                                           ("spread", args.bcc_spread if var == "full" else None),
                                           ("back", args.bcc_back and var == "full"),
                                           ("burst", args.bcc_burst and var == "full"),
                                           ("span", args.bcc_span if var == "full" else None),
                                           ("pilot128", args.bcc_pilot128 and var == "full"))
                         if v not in (None, False) or k == "touch_db" and v is not None}   # opt 0 (== False) is the
                # default: no flag; touch_db 0.0 (== False) is a real value
                try:    # get_spec refuses the cases a stack has no model for; the other entries still run
                    pr.get_spec(pid, var, p, **flags)
                except ValueError as e:
                    print(f"bcc {var} {p} skipped: {e}")
                    continue
            t = (pid, var, p, tuple(sorted(flags.items())))
            if t not in out:
                out.append(t)
    return out


def _spec(t):
    return pr.get_spec(t[0], t[1], t[2], **dict(t[3]))


def _task(a):
    t, r, k, n, mode, rec, dur, trial, seed, venue, dens, case, floor = a
    return en.simulate_trial(_spec(t), k=k, n_dev=n, r_kbps=r, ci_mode=mode, reception=rec, duration_us=dur,
                             seed=seed, venue=venue, density=dens, case=case, trial=trial, ci_floor=floor)


def _key(row):
    return tuple(str(row[c]) for c in CELL) + (str(row["trial"]), str(row["seed"]))


def _cut_partial_line(path):
    """Truncate a file after its last newline (a crash inside a write can leave half a row)."""
    with open(path, "rb+") as f:
        end = f.seek(0, 2)
        f.seek(max(0, end - 65536))
        tail = f.read()
        f.truncate(end - len(tail) + tail.rfind(b"\n") + 1)


def _norm(v):
    """The CSV spelling of a cell value (so resumed keys compare equal)."""
    return f"{v:g}" if isinstance(v, float) else str(v)


def aggregate(out_csv, cells_csv):
    """Per-cell means and 95 % intervals over trials; reads only the columns it needs (strings as categories,
    metrics as float32: every block of same-type columns stays below 70 MB for the 0.95 M rows of the default grid)
    and keeps the last copy of a row that a resumed crash wrote twice."""
    once = ["period_us", "r_actual_kbps", "payload_bits", "frames_per_period", "feasible"]
    text = ["protocol", "variant", "phy", "reception", "ci_mode", "venue", "density", "case"]
    num = AGG + ["e_sensor_j", "e_hub_j"]
    df = pd.read_csv(out_csv, usecols=CELL + ["trial", "seed", "network_id"] + once + num,
                     dtype={**{c: "category" for c in text}, **{c: "float32" for c in num}})
    df = df.drop_duplicates(CELL + ["trial", "seed", "network_id"], keep="last")
    # J per bit of a trial: a ratio of sums over its k BANs, so a BAN that delivered nothing still counts (inf when
    # the whole trial delivered nothing); a mean of per-BAN ratios would drop it and follow the BANs with few bits
    df["e_s"] = df["e_sensor_j"].astype(float) * df["n_dev"]
    df["e_b"] = df["e_s"] + df["e_hub_j"]
    df["bits"] = df["goodput_kbps"].astype(float) * df["duration_us"] / 1000
    rows = []
    for key, g in df.groupby(CELL, sort=False, dropna=False, observed=True):
        per = g.groupby(["trial", "seed"])
        tot = per[["e_s", "e_b", "bits"]].sum(min_count=1)
        jpb = {"jpb_sensor": tot["e_s"] / tot["bits"], "jpb_ban": tot["e_b"] / tot["bits"]}
        first = g.iloc[0]
        row = dict(zip(CELL, key), period_us=first["period_us"], r_actual_kbps=first["r_actual_kbps"],
                   payload_bits=first["payload_bits"], frames_per_period=first["frames_per_period"],
                   feasible=int(first["feasible"]), n_trials=per.ngroups)
        for m in AGG:
            x = (jpb[m] if m in jpb else per[m].mean()).dropna().to_numpy(float)
            mean = float(x.mean()) if x.size else np.nan
            half = (float(stats.t.ppf(0.975, x.size - 1) * x.std(ddof=1) / np.sqrt(x.size))
                    if x.size > 1 and np.isfinite(x).all() else np.nan)
            top = 1.0 if m in FRAC else np.inf     # the symmetric interval clipped to the physical range
            row.update({f"{m}_mean": mean, f"{m}_lo": float(np.clip(mean - half, 0.0, top)),
                        f"{m}_hi": float(np.clip(mean + half, 0.0, top))})
        rows.append(row)
    pd.DataFrame(rows).to_csv(cells_csv, index=False, lineterminator="\n")
    return len(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0], formatter_class=argparse.RawTextHelpFormatter,
                                 epilog=__doc__[__doc__.index("Outputs:"):])
    g = ap.add_argument_group("protocol selection")
    g.add_argument("--protocols", default="", help="comma list of ids: ble, ieee802154, tsch, ieee802156nb, "
                   "smartban, ant, hbc, bcc, ble4thesis (default: every catalogue protocol, ble4thesis excluded)")
    g.add_argument("--variants", default="full,thesis", help="full, thesis or both (thesis: the as-built BCC stack as\n"
                   "coded and ble4thesis; bcc full = the HBC stack on the galvanic channel, on one tone (phy sf8)\n"
                   "and striped over the N tones of the plan (phy multitone))")
    g.add_argument("--phy", default="", help="pid=name pairs replacing the catalogue PHYs of that protocol, e.g.\n"
                   "ble=2M,ble=coded_s8,ieee802156nb=242.9,hbc=sf16 (names: protocols.get_spec); default: the\n"
                   "catalogue's (ble 1M, 2M, coded_s2, coded_s8; ieee802154 oqpsk (CSMA) and gts; 802.15.6 NB\n"
                   "121.4, 242.9, 485.7, 971.4; smartban bch1, bch2, bch4, unc1, unc2, unc4; HBC sf8, sf16, sf32,\n"
                   "sf64; BCC full sf8 and multitone, BCC thesis bfsk)")
    g.add_argument("--bcc-contact", action="store_true", help="BCC measured-hardware case: two-state contact chain "
                   "in both reception models, at 50 kchip/s (other --bcc-rc values are refused); not for the\n"
                   "multitone entry (skipped)")
    g.add_argument("--bcc-fcs", action="store_true", help="as-built BCC with a 2-octet FCS (F7 parity); not for\n"
                   "the BCC full entries, which carry the 802.15.6 FCS (skipped)")
    g.add_argument("--bcc-coupling", type=int, choices=(40, 60), default=None, help="BCC inter-body coupling case: "
                   "intra-body loss + 40 or + 60 dB (default: isolated bodies); --reception sinr only; not for\n"
                   "the multitone entry (skipped)")
    g.add_argument("--bcc-rc", type=float, default=None, help="BCC chip rate in chip/s, every BCC entry (default "
                   "500e3; 50e3 with --bcc-contact);\nmultitone: also the tone plan (spacing 4 R_c: 6 tones at "
                   "500e3, 67 at 50e3, 1 at 2e6,\nnone above 3.38e6: skipped); one tone: none above 9e6 "
                   "(skipped)")
    g.add_argument("--bcc-power", choices=("low", "typical", "high"), default=None, help="BCC transceiver power "
                   "per state (multitone: times its N tones)")
    g.add_argument("--bcc-opt", type=int, choices=range(6), default=0, help="galvanic-optimisation ladder level of "
                   "the BCC full entries (DESIGN 4.3; variant 'full+optN'),\ncumulative, default 0 = the HBC stack:\n"
                   "1 RI mode (PHY header at the payload rate) and a pilot every 4096 chips (802.15.6 HBC options)\n"
                   "2 no I-Ack: the next beacon acknowledges the superframe's frames (proposal)\n"
                   "3 retry pool: one exchange per frame, hub-polled retries at the superframe end (802.15.6 polls)\n"
                   "4 lean framing: 1-octet MAC header (proposal)\n"
                   "5 one preamble, SFD and header per multi-tone stripe (proposal; one tone: as 4)")
    g.add_argument("--bcc-pool", type=float, default=None, help="retry-pool share of the superframe's frames, "
                   "BCC full (default 0.05 from --bcc-opt 3; 0.43 = 0.3 / 0.7\nsizes it for the measured ~30 %% contact "
                   "loss; below level 3 it replaces the spare exchange per frame:\nthe pooled reference stack of; "
                   "variant '+pool=S')")
    g.add_argument("--bcc-touch", type=float, default=None, help="share of time each pair of bodies touches and "
                   "shares the galvanic medium, one-tone BCC full, SINR only\n(variant '+touch=S')")
    g.add_argument("--bcc-touch-db", type=float, default=None, help="loss of the partner body's signal against its "
                   "own link while touching, dB (default 10; with --bcc-touch)")
    g.add_argument("--bcc-fallback", action="store_true", help="interleave the active periods of a touching pair from "
                   "one superframe after contact starts (802.15.6 6.13.3; with --bcc-touch)")
    g.add_argument("--bcc-spread", choices=("none", "adapt"), default=None, help="payload spreading of one-tone "
                   "BCC full (default Walsh-16): none = unspread BFSK,\nadapt = per link, Walsh after 2 "
                   "reported losses, unspread after 16 deliveries (from --bcc-opt 2; variant '+spread=S')")
    g.add_argument("--bcc-pool-fill", action="store_true", help="the retry pool also takes the superframe's slack "
                   "before the next beacon, BCC full from --bcc-opt 3\n(admitted on the pool share; variant "
                   "'+fill')")
    g.add_argument("--bcc-back", action="store_true", help="block acknowledgement instead of the I-Ack: L-Ack frames "
                   "and one B-Ack per allocation, the pool exchanges with\ntheir I-Ack; BCC full at --bcc-opt 0 or 1 with "
                   "--bcc-pool (reference, claim4.py H; variant '+back')")
    g.add_argument("--bcc-burst", action="store_true", help="with --bcc-back at --bcc-opt 0: the block in the HBC burst "
                   "mode, frames 2.. of an allocation pMIFS apart in RI mode\n(lead reference; variant '+burst')")
    g.add_argument("--bcc-span", type=int, default=None, help="with --bcc-burst: one B-Ack per block of this many "
                   "frames, spanning allocations; 1000000000 (protocols.BLOCK_UNLIMITED) never\ncloses the block "
                   "(lead reference; variant '+span=N')")
    g.add_argument("--bcc-pilot128", action="store_true", help="the '128 B' pilot period (802.15.6 10.7.3, Table 86) "
                   "at --bcc-opt 0\n(lead reference; variant '+pilot128')")
    g = ap.add_argument_group("model and grid")
    g.add_argument("--reception", default="collision", help="collision, sinr or both (comma list)")
    g.add_argument("--ci-floor", default="on", choices=("on", "off"), help="C/I_co floor of the sinr model "
                   "(PARAMETERS OI-3); off: interference as Gaussian noise\nfor every stack, variant label "
                   "'+nocifloor' on the sinr rows (the collision model has no floor) (default on)")
    g.add_argument("--r", default=R_COARSE, help=f"per-link rate grid, kb/s (default {R_COARSE})")
    g.add_argument("--k", default="1,2,4,8,16,30", help="coexisting BANs (default 1,2,4,8,16,30)")
    g.add_argument("--n-dev", default="1,4,10", help="sensors per BAN (default 1,4,10)")
    g.add_argument("--ci-modes", default="latency,energy", help="period choice: latency, energy or both")
    g.add_argument("--trials", type=int, default=8)
    g.add_argument("--duration-s", type=float, default=2.0, help="simulated time per trial, s (at least "
                   f"{MIN_PERIODS} periods of the cell)")
    g.add_argument("--seed", type=int, default=0, help="base seed")
    g.add_argument("--venue", default="ward", choices=("ward", "gym", "office"), help="room of the sinr model")
    g.add_argument("--density", default="typical", choices=("low", "typical", "high"))
    g.add_argument("--case", default="typical", choices=("low", "typical", "high"),
                   help="channel case (on-body and body-to-body loss sets; HBC coupling)")
    g = ap.add_argument_group("run and output")
    g.add_argument("--workers", type=int, default=4, help="worker processes (1..8, default 4)")
    g.add_argument("--out", default="coexistence_rows.csv", help="per-BAN rows, appended (resumable)")
    g.add_argument("--cells", default="coexistence_cells.csv", help="per-cell aggregates, rewritten")
    a = ap.parse_args()
    if not 1 <= a.workers <= 8:
        ap.error("--workers must be 1..8 (16 CPUs, 23 GB without swap; the machine crashed once under load)")
    if a.bcc_coupling and "collision" in a.reception.split(","):
        ap.error("--bcc-coupling needs --reception sinr: the coupling cases weaken the interferer by 40 / 60 dB, "
                 "the collision model would let any overlap kill")
    floor = a.ci_floor == "on"
    dur = int(round(a.duration_s * 1e6))
    done = set()
    if os.path.exists(a.out) and os.path.getsize(a.out):
        _cut_partial_line(a.out)
        with open(a.out, newline="", encoding="utf-8") as f:
            n_rows = Counter(_key(r) for r in csv.DictReader(f))
        done = {x for x, n in n_rows.items() if n >= int(x[CELL.index("k")])}    # all k rows (a rerun appends k)
    tasks = []
    for t in specs(a):
        sp = _spec(t)
        for rec in a.reception.split(","):
            for mode in a.ci_modes.split(","):
                for r in _floats(a.r):
                    for n in _ints(a.n_dev):
                        sel = pr.select_period(r, mode, n, sp)
                        d = max(dur, MIN_PERIODS * sel["period_us"]) if sel else dur
                        for k in _ints(a.k):
                            for tr in range(a.trials):
                                row = dict(protocol=sp.pid, variant=en.label(sp, floor, rec),
                                           phy=sp.phy["name"], reception=rec, strip_hopping=0, ci_mode=mode,
                                           r_target_kbps=r, n_dev=n, k=k, duration_us=d, venue=a.venue,
                                           density=a.density, case=a.case, trial=tr)
                                seed = zlib.crc32(f"{a.seed}|{r:g}|{k}|{n}|{mode}|{tr}".encode())
                                if tuple(_norm(row[c]) for c in CELL + ["trial"]) + (str(seed),) in done:
                                    continue
                                tasks.append((t, r, k, n, mode, rec, d, tr, seed, a.venue, a.density, a.case,
                                              floor))
    print(f"{len(tasks)} trials to run ({len(done)} already in {a.out}), {a.workers} workers", flush=True)
    new = not (os.path.exists(a.out) and os.path.getsize(a.out))
    t0 = time.time()
    # a worker killed from outside (earlyoom) raises BrokenProcessPool, where multiprocessing.Pool would wait forever
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "a", newline="", encoding="utf-8") as f, ProcessPoolExecutor(a.workers) as pool:
        w = csv.DictWriter(f, fieldnames=en.COLUMNS, lineterminator="\n")
        if new:
            w.writeheader()
        for i, fut in enumerate(as_completed([pool.submit(_task, t) for t in tasks]), 1):
            rows = fut.result()
            buf = io.StringIO()     # the trial's rows in one write
            csv.DictWriter(buf, fieldnames=en.COLUMNS, lineterminator="\n").writerows(
                {c: _norm(v) if isinstance(v, float) else v for c, v in r.items()} for r in rows)
            f.write(buf.getvalue())
            f.flush()
            if i % 50 == 0 or i == len(tasks):
                print(f"  {i}/{len(tasks)} trials, {time.time() - t0:.0f} s", flush=True)
    print(f"{aggregate(a.out, a.cells)} cells -> {a.cells}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Fill the factor-2 grid datarate 10..1280 kbit/s x sensors 1..128 x networks 1..128 for the honeycomb maps
(psr_map.py --lockmap): the lock steps above 1u were only run up to 320 kbit/s and 32 sensors / networks.

A cell is run only when none of its lower neighbours (one factor-2 step down in one parameter) is known to fail:
PSR below --floor or no schedule fits. More load never raises the PSR, so such a cell cannot pass either. Cells go
cheapest first (sensors x networks, then datarate), so a cell's lower neighbours are measured before it.

Known PSR comes from --map (psr_map.py --build output) and from the cells files this script writes:
results/triangle/cells_L_<protocol>.csv, and cells_L_prof<chip rate>.csv for galvanic BCC (the profile of
run_zoomout.sh: --bcc-opt 5 --bcc-pool-fill at both chip rates). psr_map.py --build reads them with the others.

--pilot N runs N needed cells spread over the cost range, reports the seconds per cell against the cost proxy
sensors x networks x datarate (offered frames) of the cells that had a fitting schedule, and extrapolates an
upper bound for all needed cells.

Survives earlyoom (which prefers python processes here): a sweep that loses a worker exits with BrokenProcessPool;
the cell is rerun with half the workers once --reserve-gb of memory is available (sweep.py resumes from its rows
file, so finished trials are kept).
After each clean cell one worker is added back while the memory allows, up to --workers. run_lockgrid.sh restarts
this script if it is killed itself.
"""
import argparse
import itertools as it
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE / "results" / "triangle"
R = [10 * 2 ** i for i in range(8)]
N = [2 ** i for i in range(8)]
STACKS = {"BLE": "ble", "802.15.4 CSMA/CA": "ieee802154", "802.15.4 TSCH": "tsch", "802.15.6 NB": "ieee802156nb",
          "SmartBAN": "smartban", "ANT": "ant", "802.15.6 HBC capacitive": "hbc", "Galvanic BCC (profile)": "bcc"}
BCC_RC = ("500e3", "2.625e6")


def runs(proto):
    """(rows file, cells file, extra sweep.py arguments) per sweep.py call of one cell."""
    if proto == "bcc":
        return [(OUT / f"rows_L_prof{rc}.csv", OUT / f"cells_L_prof{rc}.csv",
                 ["--bcc-opt", "5", "--bcc-pool-fill", "--bcc-rc", rc]) for rc in BCC_RC]
    return [(OUT / f"rows_L_{proto}.csv", OUT / f"cells_L_{proto}.csv", [])]


def known_psr(stack, map_csv):
    """{(datarate, sensors, networks): PSR} from the map and the cells this script has written; infeasible = 0."""
    d = pd.read_csv(map_csv)
    d = d[d["stack"] == stack]
    known = {(r, n, k): (p if f else 0.0) for r, n, k, p, f in
             zip(d.datarate, d.n_sensors, d.n_networks, d.psr.fillna(0.0), d.feasible)}
    proto = STACKS[stack]
    files = [cells for _, cells, _ in runs(proto) if cells.exists()]
    if files:
        # best configuration per cell (PHY, period mode, BCC chip rate), as psr_map.py --build takes it
        c = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
        c = c[(c.protocol == proto) & (c.reception == "sinr")]
        c = c.assign(p=np.where(c.feasible == 1, c.psr_mean.fillna(0.0), 0.0))
        known.update(c.groupby(["r_target_kbps", "n_dev", "k"]).p.max().to_dict())
    return known


def lower(cell):
    r, n, k = cell
    return [(r // 2, n, k)] * (r > R[0]) + [(r, n // 2, k)] * (n > 1) + [(r, n, k // 2)] * (k > 1)


def needed(known, floor):
    """Missing cells none of whose lower neighbours is known to fail, cheapest first."""
    miss = [c for c in it.product(R, N, N) if c not in known]
    keep = [c for c in miss if all(known.get(l, 1.0) >= floor for l in lower(c))]
    return sorted(keep, key=lambda c: (c[1] * c[2], c[0]))


def cost(cell):
    r, n, k = cell
    return r * n * k


def mem_available_gb():
    with open("/proc/meminfo") as f:
        for line in f:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / 2 ** 20
    return 0.0


def wait_for_memory(reserve_gb):
    while mem_available_gb() < reserve_gb:
        print(f"  waiting for memory: {mem_available_gb():.1f} GB available, {reserve_gb} GB wanted", flush=True)
        time.sleep(30)


def run_sweep(cmd):
    """True when sweep.py ends cleanly; a killed worker makes it exit with BrokenProcessPool, a killed sweep with
    a signal: both False."""
    return subprocess.run(cmd, cwd=HERE, stdout=subprocess.DEVNULL).returncode == 0


def run_cell(proto, cell, trials, workers):
    r, n, k = cell
    for rows, cells, extra in runs(proto):
        cmd = [sys.executable, str(HERE / "sweep.py"), "--protocols", proto, "--variants", "full", *extra,
               "--reception", "sinr", "--ci-modes", "latency,energy", "--trials", str(trials), "--workers",
               str(workers), "--out", str(rows), "--cells", str(cells), "--r", str(r), "--n-dev", str(n), "--k", str(k)]
        if not run_sweep(cmd):
            return False
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stacks", default="Galvanic BCC (profile),BLE,SmartBAN,802.15.6 NB",
                    help="comma list of psr_map stack names, run in this order")
    ap.add_argument("--map", default=str(OUT / "psr_map_fixed.csv"))
    ap.add_argument("--floor", type=float, default=0.5, help="lowest PSR target the maps show")
    ap.add_argument("--trials", type=int, default=8)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--pilot", type=int, default=0, help="time this many cells per stack and extrapolate, then stop")
    ap.add_argument("--dry", action="store_true", help="list the needed cells only")
    ap.add_argument("--reserve-gb", type=float, default=4.0, help="available memory to wait for before a rerun "
                    "and to keep free when adding a worker back")
    a = ap.parse_args()
    for stack in a.stacks.split(","):
        proto = STACKS[stack]
        todo = needed(known_psr(stack, a.map), a.floor)
        print(f"{stack}: {len(todo)} cells needed now (upper bound; failing cells prune the rest)", flush=True)
        if a.dry:
            continue
        if a.pilot:
            pick = sorted(todo, key=cost)
            pick = [pick[int(round(q * (len(pick) - 1)))] for q in np.linspace(0, 1, a.pilot)]
            secs = []
            for cell in dict.fromkeys(pick):
                t0 = time.time(); run_cell(proto, cell, a.trials, a.workers); dt = time.time() - t0
                secs.append((cost(cell), dt))
                p = known_psr(stack, a.map).get(cell)
                print(f"  pilot {cell} (r, sensors, networks): {dt:7.1f} s, PSR {p}", flush=True)
            # a cell without a fitting schedule ends in seconds: the rate comes from the cells that ran trials
            ran = [(x, y) for (x, y), c in zip(secs, dict.fromkeys(pick)) if known_psr(stack, a.map).get(c, 0) > 0]
            if not ran:
                print("  no pilot cell had a fitting schedule: no rate", flush=True)
                continue
            x, y = np.array(ran).T
            per = np.sum(y) / np.sum(x)
            est = per * sum(cost(c) for c in todo)
            print(f"  {per * 1e3:.3f} ms per unit of r x sensors x networks (cells with a schedule); all {len(todo)} "
                  f"needed cells at that rate ~ {est / 3600:.1f} h, an upper bound: cells without a schedule and "
                  "pruned cells cost seconds", flush=True)
            continue
        workers = a.workers
        while todo:
            cell = todo[0]
            wait_for_memory(a.reserve_gb)
            t0 = time.time()
            if not run_cell(proto, cell, a.trials, workers):
                workers = max(1, workers // 2)
                print(f"  {cell}: a worker died or the sweep failed after {time.time() - t0:.0f} s; "
                      f"rerun with {workers} workers", flush=True)
                continue
            p = known_psr(stack, a.map).get(cell)
            todo = needed(known_psr(stack, a.map), a.floor)
            print(f"  {cell}: {time.time() - t0:6.1f} s, PSR {p}, {workers} workers; {len(todo)} left", flush=True)
            if workers < a.workers and mem_available_gb() > a.reserve_gb * 2:
                workers += 1
    print("LOCKGRID_DONE", flush=True)


if __name__ == "__main__":
    main()

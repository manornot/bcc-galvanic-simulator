"""
BLE coexistence simulator (event-driven, throughput-axis edition).

Sweep axes
----------
- Primary x-axis: per-link throughput R = L / CI (kbps), application-natural.
- For each R, two CI choices are simulated:
      latency-optimal: CI = 7.5 ms (BLE Core minimum)
      energy-optimal:  CI = max(7.5 ms, min(L_max / R, 4 s))
                       -- longest CI consistent with the BLE max payload.
  These two cases bracket the real-world deployment envelope.
- y-axis: number of coexisting networks k.
- Sliced over: devices-per-network N_dev (from --n-dev).
- Run for both PHYs (1M with random_loss_p = 1e-4 by default,
  2M with random_loss_p = 1e-3 by default).

Adaptive refinement
-------------------
For each (PHY, ci_mode, N_dev, k) cell, the throughput axis is sampled
adaptively. A coarse R-grid is run first; intervals whose mean PSRs
straddle any of the contour levels (0.99, 0.95, 0.90) are bisected
until either |dPSR| < tolerance or a recursion-depth limit is hit.
Coarse-only outside contour-straddling intervals; fine inside them.

Collision rule (verified against Ancans et al. 2019): any time-overlap
on the same channel kills all overlapping packets. No capture, no partial
recovery, no link-layer ACK, no retransmission. Surviving packets are
additionally lost with independent probability `random_loss_p` to model
the physical-layer PER from thermal noise; this defaults to 1e-4 for
1M and 1e-3 for 2M (the ~10x BER penalty of doubled symbol rate).

Output schema
-------------
One CSV row per network per trial, with full parameters in every row.
Columns: phy_rate, ci_mode, ci_us, r_target_kbps, r_actual_kbps,
num_devices, k, payload_bits, duration_us, trial_id, seed,
network_id, packets_sent, packets_received, success_rate, random_loss_p
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import random
import sys
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

NUM_BLE_CHANNELS = 37
DEFAULT_CONTOUR_LEVELS = (0.99, 0.95, 0.90)

# BLE Core constants
CI_MIN_US = 7_500           # 7.5 ms minimum
CI_MAX_US = 4_000_000       # 4 s maximum
CI_GRANULARITY_US = 1_250   # 1.25 ms quantum
PAYLOAD_MAX_BYTES = 251     # BLE 4.2+ Data PDU max
PAYLOAD_MIN_BYTES = 1


def quantise_ci_us(ci_us: int) -> int:
    """Round CI to the nearest 1.25 ms multiple, clamped to [7.5 ms, 4 s]."""
    q = round(ci_us / CI_GRANULARITY_US) * CI_GRANULARITY_US
    return max(CI_MIN_US, min(CI_MAX_US, q))


# ---------------------------------------------------------------------------
# Single-trial event-driven engine (unchanged from v3)
# ---------------------------------------------------------------------------
def simulate_one_trial(
    *,
    k: int,
    num_devices: int,
    payload_bits: int,
    ci_us: int,
    phy_rate: int,
    duration_us: int,
    seed: int,
    random_loss_p: float = 0.0,
) -> List[Dict[str, int]]:
    rng = random.Random(seed)

    t_on_us = max(1, math.ceil(payload_bits * 1_000_000 / phy_rate))
    slot_us = ci_us // num_devices
    if slot_us < t_on_us:
        raise ValueError(
            f"Slot time ({slot_us} us) shorter than packet air-time ({t_on_us} us)."
        )

    hop_sequences: List[List[int]] = []
    ci_phase_us: List[int] = []
    for _ in range(k):
        seq = list(range(NUM_BLE_CHANNELS))
        rng.shuffle(seq)
        hop_sequences.append(seq)
        ci_phase_us.append(rng.randint(0, ci_us - 1))

    events: List[Tuple[int, int, int, int]] = []
    for net_id in range(k):
        seq = hop_sequences[net_id]
        phase = ci_phase_us[net_id]
        n_ci = (duration_us - phase + ci_us - 1) // ci_us
        if n_ci <= 0:
            continue
        for cycle in range(n_ci):
            ch = seq[cycle % NUM_BLE_CHANNELS]
            ci_start = phase + cycle * ci_us
            for node_idx in range(num_devices):
                start = ci_start + node_idx * slot_us
                if start >= duration_us:
                    break
                events.append((start, start + t_on_us, net_id, ch))

    packets_sent = [0] * k
    for _, _, net_id, _ in events:
        packets_sent[net_id] += 1

    failed = [False] * len(events)

    by_channel: Dict[int, List[int]] = {}
    for idx, (_, _, _, ch) in enumerate(events):
        by_channel.setdefault(ch, []).append(idx)

    for idxs in by_channel.values():
        idxs.sort(key=lambda i: events[i][0])
        n = len(idxs)
        i = 0
        while i < n:
            j = i + 1
            cluster_end = events[idxs[i]][1]
            cluster = [idxs[i]]
            while j < n and events[idxs[j]][0] < cluster_end:
                cluster.append(idxs[j])
                if events[idxs[j]][1] > cluster_end:
                    cluster_end = events[idxs[j]][1]
                j += 1
            if len(cluster) > 1:
                for c in cluster:
                    failed[c] = True
            i = j

    if random_loss_p > 0.0:
        for idx in range(len(events)):
            if not failed[idx] and rng.random() < random_loss_p:
                failed[idx] = True

    packets_received = [0] * k
    for idx, (_, _, net_id, _) in enumerate(events):
        if not failed[idx]:
            packets_received[net_id] += 1

    return [
        {"network_id": nid,
         "packets_sent": packets_sent[nid],
         "packets_received": packets_received[nid]}
        for nid in range(k)
    ]


# ---------------------------------------------------------------------------
# CI selection from throughput
# ---------------------------------------------------------------------------
@dataclass
class CISelection:
    """A (CI, payload) pair that delivers the requested throughput."""
    ci_us: int
    payload_bits: int
    r_actual_kbps: float


def select_ci(r_kbps: float, mode: str, n_dev: int, phy_rate: int) -> Optional[CISelection]:
    """Pick (CI, L) for the given throughput and mode, or None if infeasible.

    Latency-optimal: CI = 7.5 ms (BLE minimum). L = R * CI.
    Energy-optimal:  CI = longest CI such that L <= L_max. Bounded by 4 s.

    Returns None when the resulting payload would be smaller than 1 byte
    or larger than 251 bytes, or when the per-network slot is too small
    for the resulting packet on-air time.
    """
    if mode == "latency":
        ci_us = CI_MIN_US
        # L (in bits) = R * CI = r_kbps * 1000 * CI_us / 1e6 = r_kbps * CI_us / 1000
        payload_bits = round(r_kbps * ci_us / 1000.0)
    elif mode == "energy":
        # Choose the longest CI such that payload <= 251 B.
        # L_bits = r_kbps * ci_us / 1000  ==>  ci_us = L_bits * 1000 / r_kbps
        max_payload_bits = PAYLOAD_MAX_BYTES * 8
        if r_kbps <= 0:
            return None
        ci_us_raw = int(max_payload_bits * 1000.0 / r_kbps)
        ci_us = quantise_ci_us(ci_us_raw)
        payload_bits = round(r_kbps * ci_us / 1000.0)
    else:
        raise ValueError(f"Unknown ci_mode: {mode}")

    # Feasibility checks.
    payload_bytes = payload_bits / 8
    if payload_bytes < PAYLOAD_MIN_BYTES or payload_bytes > PAYLOAD_MAX_BYTES:
        return None

    # Slot must accommodate the packet's on-air time.
    t_on_us = math.ceil(payload_bits * 1_000_000 / phy_rate)
    slot_us = ci_us // n_dev
    if slot_us < t_on_us:
        return None

    r_actual_kbps = payload_bits * 1000.0 / ci_us
    return CISelection(ci_us=ci_us, payload_bits=int(payload_bits),
                       r_actual_kbps=r_actual_kbps)


# ---------------------------------------------------------------------------
# Cached evaluator and adaptive refinement
# ---------------------------------------------------------------------------
CSV_HEADER = [
    "phy_rate", "ci_mode", "ci_us", "r_target_kbps", "r_actual_kbps",
    "num_devices", "k", "payload_bits", "duration_us",
    "trial_id", "seed",
    "network_id", "packets_sent", "packets_received", "success_rate",
    "random_loss_p",
]


@dataclass
class PointResult:
    mean_psr: float
    n_trials: int


@dataclass
class SweepConfig:
    phy_rate: int
    ci_modes: List[str]
    n_dev_values: List[int]
    k_values: List[int]
    r_coarse_kbps: List[float]
    n_trials: int
    duration_us: int
    random_loss_p: float
    contour_levels: Tuple[float, ...] = DEFAULT_CONTOUR_LEVELS
    refine_tol: float = 0.01
    refine_max_depth: int = 5
    base_seed: int = 0


class CachedEvaluator:
    def __init__(self, cfg: SweepConfig, writer, flush_fn):
        self.cfg = cfg
        self.writer = writer
        self.flush_fn = flush_fn
        # cache key: (phy, ci_mode, r_target_round, n_dev, k)
        self._cache: Dict[Tuple[int, str, float, int, int], PointResult] = {}
        self.point_count = 0
        self.cache_hits = 0
        self.skipped_infeasible = 0

    @staticmethod
    def _key(phy, ci_mode, r, n_dev, k):
        return (phy, ci_mode, round(r, 4), n_dev, k)

    def evaluate(self, r_kbps: float, ci_mode: str, n_dev: int, k: int) -> Optional[PointResult]:
        cfg = self.cfg
        key = self._key(cfg.phy_rate, ci_mode, r_kbps, n_dev, k)
        if key in self._cache:
            self.cache_hits += 1
            return self._cache[key]

        sel = select_ci(r_kbps, ci_mode, n_dev, cfg.phy_rate)
        if sel is None:
            self.skipped_infeasible += 1
            return None

        per_trial_psr: List[float] = []
        for trial_id in range(cfg.n_trials):
            seed = (cfg.base_seed
                    + 1000 * trial_id
                    + 17 * k
                    + 31 * n_dev
                    + (1 if ci_mode == "latency" else 2) * 7
                    + int(r_kbps * 100) % 9973)
            results = simulate_one_trial(
                k=k,
                num_devices=n_dev,
                payload_bits=sel.payload_bits,
                ci_us=sel.ci_us,
                phy_rate=cfg.phy_rate,
                duration_us=cfg.duration_us,
                seed=seed,
                random_loss_p=cfg.random_loss_p,
            )
            trial_sent = 0
            trial_recv = 0
            for r in results:
                sent = r["packets_sent"]
                recv = r["packets_received"]
                rate = (recv / sent) if sent else 0.0
                self.writer.writerow([
                    cfg.phy_rate, ci_mode, sel.ci_us,
                    f"{r_kbps:.4f}", f"{sel.r_actual_kbps:.4f}",
                    n_dev, k, sel.payload_bits, cfg.duration_us,
                    trial_id, seed,
                    r["network_id"], sent, recv, f"{rate:.6f}",
                    f"{cfg.random_loss_p:.6g}",
                ])
                trial_sent += sent
                trial_recv += recv
            per_trial_psr.append(trial_recv / trial_sent if trial_sent else 0.0)
        self.flush_fn()

        mean_psr = sum(per_trial_psr) / len(per_trial_psr)
        result = PointResult(mean_psr=mean_psr, n_trials=cfg.n_trials)
        self._cache[key] = result
        self.point_count += 1
        return result


def _interval_straddles(p_lo, p_hi, contours):
    a, b = sorted((p_lo, p_hi))
    return any(a < c < b for c in contours)


def adaptive_r_sweep(evaluator: CachedEvaluator, cfg: SweepConfig,
                     ci_mode: str, n_dev: int, k: int, log) -> None:
    coarse: List[Tuple[float, Optional[PointResult]]] = []
    for r in cfg.r_coarse_kbps:
        coarse.append((r, evaluator.evaluate(r, ci_mode, n_dev, k)))

    queue: List[Tuple[float, float, float, float, int]] = []
    for (r_a, res_a), (r_b, res_b) in zip(coarse, coarse[1:]):
        if res_a is None or res_b is None:
            continue
        if _interval_straddles(res_a.mean_psr, res_b.mean_psr, cfg.contour_levels):
            queue.append((r_a, r_b, res_a.mean_psr, res_b.mean_psr, 0))

    while queue:
        r_lo, r_hi, p_lo, p_hi, depth = queue.pop()
        if depth >= cfg.refine_max_depth:
            continue
        if abs(p_hi - p_lo) < cfg.refine_tol:
            continue
        r_mid = 0.5 * (r_lo + r_hi)
        res_mid = evaluator.evaluate(r_mid, ci_mode, n_dev, k)
        if res_mid is None:
            continue
        log(f"      refine {ci_mode} n_dev={n_dev} k={k} "
            f"R=[{r_lo:.2f},{r_hi:.2f}] depth={depth} "
            f"mid={r_mid:.2f} -> PSR={res_mid.mean_psr:.4f}")
        if _interval_straddles(p_lo, res_mid.mean_psr, cfg.contour_levels):
            queue.append((r_lo, r_mid, p_lo, res_mid.mean_psr, depth + 1))
        if _interval_straddles(res_mid.mean_psr, p_hi, cfg.contour_levels):
            queue.append((r_mid, r_hi, res_mid.mean_psr, p_hi, depth + 1))


# ---------------------------------------------------------------------------
# Top-level driver
# ---------------------------------------------------------------------------
def run_landscape(cfg: SweepConfig, out_path: str, progress: bool = True) -> None:
    file_exists = os.path.exists(out_path)
    f = open(out_path, "a", newline="", buffering=1)
    writer = csv.writer(f)
    if not file_exists:
        writer.writerow(CSV_HEADER)

    def _flush():
        f.flush()

    def _log(msg):
        if progress:
            sys.stdout.write(msg + "\n")
            sys.stdout.flush()

    evaluator = CachedEvaluator(cfg, writer, _flush)
    t0 = time.perf_counter()

    total_cells = len(cfg.ci_modes) * len(cfg.n_dev_values) * len(cfg.k_values)
    cell_idx = 0

    try:
        for ci_mode in cfg.ci_modes:
            for n_dev in cfg.n_dev_values:
                for k in cfg.k_values:
                    cell_idx += 1
                    _log(f"[{cell_idx}/{total_cells}] phy={cfg.phy_rate//1_000_000}M "
                         f"ci_mode={ci_mode} n_dev={n_dev} k={k}")
                    adaptive_r_sweep(evaluator, cfg, ci_mode, n_dev, k, _log)
                    elapsed = time.perf_counter() - t0
                    rate = evaluator.point_count / elapsed if elapsed > 0 else 0.0
                    _log(f"    cumulative: {evaluator.point_count} unique points, "
                         f"{evaluator.cache_hits} cache hits, "
                         f"{evaluator.skipped_infeasible} infeasible, "
                         f"{rate:.1f} pts/s, {elapsed:.0f}s")
    finally:
        f.close()
        if progress:
            sys.stdout.write("\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _parse_float_list(s):
    return [float(x) for x in s.split(",")]


def _parse_int_list(s):
    if "-" in s and "," not in s:
        a, b = s.split("-", 1)
        return list(range(int(a), int(b) + 1))
    return [int(x) for x in s.split(",")]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="BLE coexistence landscape sweep over (R, ci_mode, N_dev, k)."
    )
    parser.add_argument("--out", default="ble_landscape.csv")
    parser.add_argument("--phy", choices=["1M", "2M", "both"], default="1M")
    parser.add_argument("--ci-modes", default="latency,energy",
                        help="Subset of {latency, energy}.")
    parser.add_argument("--n-dev", default="1,2,4,6,10")
    parser.add_argument("--k", default="1,2,4,8,16,30")
    parser.add_argument("--r-coarse", default="1,5,10,30,50,100,200,500,1000",
                        help="Initial coarse throughput grid in kbps.")
    parser.add_argument("--trials", type=int, default=8)
    parser.add_argument("--duration-s", type=float, default=10.0)
    parser.add_argument("--refine-tol", type=float, default=0.01)
    parser.add_argument("--refine-max-depth", type=int, default=5)
    parser.add_argument("--random-loss-p", type=float, default=None,
                        help="Default 1e-4 for 1M, 1e-3 for 2M.")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    phys: List[Tuple[int, float]] = []
    if args.phy in ("1M", "both"):
        phys.append((1_000_000, 1e-4))
    if args.phy in ("2M", "both"):
        phys.append((2_000_000, 1e-3))

    ci_modes = [m.strip() for m in args.ci_modes.split(",") if m.strip()]
    for m in ci_modes:
        if m not in ("latency", "energy"):
            parser.error(f"Unknown ci_mode '{m}'.")

    n_dev_values = _parse_int_list(args.n_dev)
    k_values = _parse_int_list(args.k)
    r_coarse = _parse_float_list(args.r_coarse)
    duration_us = int(args.duration_s * 1_000_000)

    print(f"PHYs:        {[f'{p[0]//1_000_000}M' for p in phys]}")
    print(f"ci_modes:    {ci_modes}")
    print(f"N_dev:       {n_dev_values}")
    print(f"k:           {k_values[:8]}{'...' if len(k_values) > 8 else ''}")
    print(f"R coarse:    {r_coarse} kbps")
    print(f"trials/cell: {args.trials}, duration: {args.duration_s} s\n")

    for phy_rate, default_loss in phys:
        loss = args.random_loss_p if args.random_loss_p is not None else default_loss
        cfg = SweepConfig(
            phy_rate=phy_rate,
            ci_modes=ci_modes,
            n_dev_values=n_dev_values,
            k_values=k_values,
            r_coarse_kbps=r_coarse,
            n_trials=args.trials,
            duration_us=duration_us,
            random_loss_p=loss,
            refine_tol=args.refine_tol,
            refine_max_depth=args.refine_max_depth,
            base_seed=args.seed,
        )
        print(f"=== sweeping PHY {phy_rate//1_000_000}M  random_loss_p={loss} ===")
        run_landscape(cfg, args.out, progress=not args.quiet)


if __name__ == "__main__":
    main()
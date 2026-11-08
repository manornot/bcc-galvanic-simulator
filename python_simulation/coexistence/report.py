"""Tables, figures and compact CSVs of the full-stack protocol comparison.

Usage, from python_simulation/coexistence/:   python report.py [--raw DIR]

Inputs
  DIR/*_cells.csv (--raw, default results/raw)  per-cell means and 95 % intervals of every sweep of run_sweeps.sh
                                 (git-ignored). When present they rewrite the compact CSVs below (4 significant
                                 digits; a compact CSV none of whose sweeps is in DIR is kept, and the sweeps read
                                 are printed); every table and figure is then computed from the compact CSVs, so a clone
                                 without the raw sweep regenerates the same output. A last digit can differ from the
                                 raw values.
  DIR/ladder_contact_rows.csv    retx_per_frame per BAN of the ladder's contact case (T10): its cell mean goes into
                                 cells_ladder.csv (the cells files carry no retransmission count).
  protocols.py                   BCC frame air times, exchanges, allocations and capacity per body (select_period
                                 scanned in per-body steps of 0.5 kb/s) at every chip rate (T6a, T12, F8).
  report.WORKING_POINTS          working points A-H.
  ble_landscape_cells.csv         BLE_SIM cells (T9), not part of this repository; T9 is skipped without it.

Series (sid): each stack at its full capability, one line or row each: BLE (1M, 2M, Coded S=2, S=8), 802.15.4 (CSMA,
  GTS), TSCH, 802.15.6 NB (4 rates), SmartBAN (6 PHY schemes), ANT, IEEE HBC capacitive (SF 8 to 64), BCC full (one
  tone or multi-tone, any chip rate), BCC galvanic (the proposal: BCC full at ladder level 5, retry pool 0.05, one tone
  or multi-tone, any chip rate; DESIGN 4.3) and BCC as built (variant 'thesis', any chip rate). Link adaptation: per
  cell a stack takes its best configuration (rate, PHY, tone plan, chip rate) and ci mode, the highest PSR, PSRs within
  TIE tie and the lower sensor J/bit wins (pick). The configuration of each row stays in cid (protocol-phy) and rc
  (BCC chip rate); T11 counts which configuration wins where. Chip rates: 500 kchip/s from the main sweep (the proposal
  from the ladder sweep, level 5, pool 0.05), the others from the rates sweep. The plain BLE model (ble4thesis,
  S0) is used by T9 only. Rows of any other variant are dropped at load.
  BCC pooled reference: the HBC stack with one exchange per frame in the allocation and its
  retries in a hub-polled pool of 5 % of the superframe's frames, a hub policy the standard permits (BCC full with the
  pool flag at level 0; the lead reference). It is a reference line of T6a, T10, T12, F8 and F9 only (500 kchip/s
  from the ladder sweeps, 2.625 Mchip/s from the rates sweep), not a series of the landscape, dominance or scenario
  tables; BCC full stays the spare-exchange hub (one exchange + one spare per frame).
  Gateway limit (GATEWAY): a BLE cell with more sensors than the nRF52840 S140 SoftDevice central keeps connections
  (20) is infeasible; 802.15.4 GTS (7), SmartBAN (16) and ANT (15 channels) carry theirs in protocols max_dev; the
  802.15.6 BAN size (mMaxBANSize 64) lies above every grid N_dev. A BLE mesh or a second hub is not modelled.
  Ladder (T10, F9): BCC full one tone and multi-tone at levels L0 (the main sweep) to L5, pool 0.05 and, from L3,
  0.43 (sweeps ladder and ladder_contact), at 500 kchip/s; the pooled reference as its own row at L0, pool 0.05.

Outputs (all under results/)
  cells_main.csv, cells_sensitivity.csv, cells_ladder.csv, cells_rates.csv, cells_dense.csv  the cell metrics the
                                         tables and figures use, 4 significant digits (a missing file is skipped).
  landscape.csv                          PSR landscape value per series, reception, R_agg band and k.
  tables/t1..t14_*.md, t3b_*.md          markdown tables, included verbatim by REPORT.md (T3b: per R_agg split).
  figs/f1a..f10b_*.pdf and .png          figures; F1 and F2 on two pages each (a, b), 6.3 in wide (A4 text).

Rules
  R_agg of a cell = N_dev x r_target_kbps (the same grid cells for every protocol); bands [lo, hi).
  PSR landscape (T1, T2, F3, landscape.csv), rule: per (N_dev, k, R) the highest PSR over the stack's
    configurations and the two ci modes, then the median over N_dev at each R_agg grid point, then the median over
    the points of a band. Intervals: the same medians of the chosen cells' 95 % limits. F1, F2: the class per R_agg
    grid point (the first two steps, no band median), so a tile can differ from its band in T1. F3, F6: a line at half
    its width marks a band series with fewer feasible grid points than the band has (the * of T1).
  Chip rates: T1-T3b, T5, T11, T13, T14 and F1-F6, F10 take BCC at its best chip rate up to RC_HEAD; the faster rates only in T12, F8 and the upper-bound view of T13 (UPPER).
  Other metrics (T3-T6, T11, T13, T14, F4, F5, F10): per cell the link-adaptation row (pick 'better'), BCC at its best
    chip rate (RATE_NOTE); T4, T6 and the 500 kchip/s view of T13 and F4 hold BCC at 500 kchip/s. T7: BCC
    configurations paired at the same ci mode; T8: each stack's best configuration per ci mode, 500 kchip/s for BCC.
    F6: the latency mode (P = 7.5 ms). F7 and the T6 contact case: every feasible row (frame length matters).
  Working points (T3, T6b, F5, T14): the cells whose R_agg and k lie in the point's ranges (inclusive) and whose N_dev
    is the grid value inside its N_dev range (all such values), else the nearest grid value(s); the median over the
    feasible matching cells, with the count.
  Dominance (T13, F10): per cell (model, N_dev, R, k) of the main and S6 dense sweeps the stacks with psr >= 0.99
    (F10b, T13: also on_time >= 0.9), and among them the lowest sensor J/bit; T13 also with BCC at 500 kchip/s only
    and with every chip rate to 9 Mchip/s (the upper bound).
  Scenarios (T14): the smallest grid cell at or above the need; no winner when its load is below the need or above
    1.5 times it; J/bit at the sweep's ci-mode periods, not at a period chosen per deadline.
  Ladder (T10, F9) and chip rates (T12, F8): collision model, the rows of each ci mode as run (no better-mode pick); at
    a point (N_dev x R) the median over k of the feasible cells (the isolated BCC bodies do not interact).
"""
import argparse
import math
import re
import time
from functools import lru_cache
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter  # noqa: E402

import protocols  # noqa: E402

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
TAB, FIG = RES / "tables", RES / "figs"
WORKING_POINTS = r"""A & Overnight skin patch                          & $8$        & $4$--$10$\,s     & $1$--$3$  & $0.006$--$0.048$  & $1$--$10$ & $0.005$--$0.013$ \\
B & Ward vital-signs patch                        & $8$--$20$  & $0.5$--$1.0$\,s  & $1$--$5$  & $0.06$--$1.6$     & $1$--$20$ & $0.05$--$0.12$   \\
C & Elderly fall monitor                          & $12$--$32$ & $20$--$50$\,ms   & $1$--$3$  & $1.9$--$38$       & $1$--$15$ & $1.1$--$3.6$     \\
D & Fitness band, $100$\,Hz inertial              & $20$--$50$ & $7.5$--$10$\,ms  & $1$--$2$  & $16$--$107$       & $1$--$10$ & $6.2$--$11.5$    \\
E & Motion-capture suit                           & $32$--$64$ & $7.5$\,ms        & $10$      & $341$--$683$      & $1$--$5$  & $9.5$--$13.0$    \\
F & 8-channel electromyography                    & $32$--$128$& $7.5$--$10$\,ms  & $1$       & $26$--$137$       & $1$--$5$  & $7.2$--$19.8$    \\
G & 16-channel electroencephalography             & $125$--$250$& $7.5$--$10$\,ms & $1$       & $100$--$267$      & $1$--$3$  & $14.6$--$32.8$   \\
H & Multi-bed intensive-care suite                & $125$--$250$& $7.5$--$10$\,ms & $4$--$8$  & $400$--$2133$     & $3$--$10$ & as in G          \\
"""
THESIS_CELLS = HERE / "ble_landscape_cells.csv"

KS = [1, 2, 4, 8, 16, 30]
TIE = 0.005    # PSRs this close count as equal: 95 % of the k 1 gaps between the ci modes (p_loss noise)
EDGES = [0, 1, 10, 30, 100, 300, 500, 1000, 3000, 1e9]
BANDS = ["0-1", "1-10", "10-30", "30-100", "100-300", "300-500", "500-1000", "1000-3000", "3000-10000"]
MAIN_SWEEPS = ["main", "s0_ble4thesis"]
DOM_SWEEPS = ["main", "s6_dense"]
SK = ["sweep", "sid", "option", "reception"]                       # one series of one sweep
CK = SK + ["r_target_kbps", "n_dev", "k"]                          # one cell (every configuration and ci mode)
KEEP = ["period_us", "r_actual_kbps", "payload_bits", "feasible", "psr_mean", "psr_lo", "psr_hi", "on_time_mean",
        "on_time_lo", "on_time_hi", "goodput_kbps_mean", "lat_p95_ms_mean", "lat_p95_ms_lo", "lat_p95_ms_hi",
        "jpb_sensor_mean", "jpb_ban_mean", "lifetime_days_mean", "raw_ber_mean", "post_ber_mean",
        "collision_frac_mean", "jain_fairness_mean", "n_trials"]
COMPACT = {"main": "cells_main.csv", "s0_ble4thesis": "cells_main.csv", "rates": "cells_rates.csv",
           "s6_dense": "cells_dense.csv"}       # else cells_ladder.csv (ladder*) or cells_sensitivity.csv

# series id -> label, in table order: each stack at full capability (link adaptation over what the label lists)
SERIES = [("ble", "BLE"), ("ieee802154", "802.15.4"), ("tsch", "TSCH"), ("ieee802156nb", "802.15.6 NB"),
          ("smartban", "SmartBAN"), ("ant", "ANT"), ("hbc", "IEEE HBC (capacitive)"), ("bcc-full", "BCC full"),
          ("bcc-opt", "BCC galvanic (proposal)"), ("bcc-thesis", "BCC as built")]
# the pooled reference is a reference line of T6a, T10, T12, F8 and F9, not a series of SERIES: in SERIES it
# would enter every landscape, dominance and scenario table as one more BCC row
LABEL = {**dict(SERIES), "bcc-pooled": "BCC pooled reference"}
SID = [s for s, _ in SERIES]
CONTEND = [s for s in SID if s != "bcc-thesis"]     # full stacks: the as-built BCC is a reference line only
BCC = ["bcc-full", "bcc-opt", "bcc-thesis"]
ADAPT = {"ble": "1M, 2M, Coded S=2, S=8", "ieee802154": "CSMA, GTS", "tsch": "one PHY",
         "ieee802156nb": "121.4, 242.9, 485.7, 971.4 kb/s", "smartban": "BCH and uncoded, Nrep 1, 2, 4",
         "ant": "one PHY", "hbc": "SF 8, 16, 32, 64", "bcc-full": "one tone or multi-tone, every chip rate",
         "bcc-opt": "one tone or multi-tone at ladder level 5, pool 0.05, every chip rate",
         "bcc-pooled": "one tone or multi-tone at ladder level 0, pool 0.05, 500 kchip/s and 2.625 Mchip/s",
         "bcc-thesis": "every chip rate"}
# the headline and every verdict take BCC at most at the fastest IEEE 802.15.6 HBC Walsh-chip rate;
# faster rates (to protocols.BCC_RC_MAX_CPS) only in T12, F8 and the upper-bound view of T13
RC_HEAD = 2.625e6
UPPER = "model upper bound: no phase data above what bounds; channel above 20 MHz unmeasured"
# what the chip-rate model holds (protocols.get_spec, results/bcc_rate_stability.md)
RATE_MODEL = ("the SNR per chip falls as 1/R_c at a fixed 1 V drive (no power increase), the noise band is 4 R_c, the "
              "tone gains are the measured path gains at f -/+ R_c, and the minimum-phase smear check shows "
              "no loss up to 2.625 Mchip/s (<= 0.03 dB; 0.1-0.27 dB at 5.25 M, about 1 dB at 9 M)")
# the chip-rate basis of every table and figure that takes a BCC series at its best chip rate per cell
RATE_NOTE = ("BCC series at their best chip rate per cell up to 2.625 Mchip/s (the fastest IEEE 802.15.6 HBC "
             "Walsh-chip rate, the headline limit; 50 k to 2.625 Mchip/s in the main grid, 0.5 to 2.625 Mchip/s in "
             "S6 dense, 500 kchip/s in the other sweeps; the faster rates to 9 Mchip/s only in T12, F8 and T13's "
             "upper bound); the BCC transceiver-IC power is held constant over the chip rate, so BCC J/bit falls as "
             "the rate rises; " + RATE_MODEL)
# profile version of the T3 / T3b working-point tables : the sweep CSVs predate the
# later profile changes; claim4.py uses with the current profile
T3_PROFILE_NOTE = ("The proposal series is the profile of the report v3 run (6b5b6b5, 2026-09-25: BCC full at "
                   "ladder level 5, pool share 0.05). It predates the slack-filling pool, the 608-chip SFD/RI "
                   "field and the 2-octet MAC header with the beacon ack field, and has no GTn. The "
                   "Capacity figures come from claim4.py with the current profile. Not regenerated.")
RATE_FIG = "BCC at its best chip rate per cell up to 2.625 Mchip/s (headline limit; IC power constant over the rate)"
# gateway limit not in protocols max_dev: the modelled BLE controller (nRF52840, S140 SoftDevice) keeps at most 20
# concurrent connections (S140 SoftDevice Specification v7, multilink); above it BLE needs a second hub or a mesh
GATEWAY = {"ble": 20}
# BCC configurations of T6a, T12, F8: (series, cid, ladder level, pool share, label), in series order: the
# spare-exchange hub (BCC full), the pooled reference, the proposal, the as-built stack
BCFG = [("bcc-full", "bcc-sf8", 0, None, "BCC full, one tone"),
        ("bcc-full", "bcc-multitone", 0, None, "BCC full, multi-tone"),
        ("bcc-pooled", "bcc-sf8", 0, 0.05, "pooled reference, one tone"),
        ("bcc-pooled", "bcc-multitone", 0, 0.05, "pooled reference, multi-tone"),
        ("bcc-opt", "bcc-sf8", 5, None, "proposal, one tone"),
        ("bcc-opt", "bcc-multitone", 5, None, "proposal, multi-tone"),
        ("bcc-thesis", "bcc-bfsk", 0, None, "BCC as built")]
LCID = ["bcc-sf8", "bcc-multitone"]        # the ladder's stacks (T10, F9)
# colour per series (dataviz reference palette, light mode): one hue per radio stack; BCC a neutral family (full ink,
# the as-built stack light grey dashed) and the proposal a dark wine with star marks; the BCC configurations of F8
# and F9 (cid) ink (one tone) and dark grey (multi-tone)
COLOUR = {"ble": "#3987e5", "ieee802154": "#eb6834", "tsch": "#1baf7a", "ieee802156nb": "#eda100",
          "smartban": "#e87ba4", "ant": "#008300", "hbc": "#4a3aa7", "bcc-full": "#0b0b0b", "bcc-opt": "#6e1846",
          "bcc-pooled": "#6e6d68", "bcc-thesis": "#9a9993", "bcc-sf8": "#0b0b0b", "bcc-multitone": "#5c5b57",
          "bcc-bfsk": "#9a9993"}
MARKER = {"bcc-full": "D", "bcc-opt": "*", "bcc-pooled": "s", "bcc-thesis": "D", "bcc-sf8": "D", "bcc-multitone": "v",
          "bcc-bfsk": "D"}
DEADLINE = "#d03b3b"
PANELS = [("BLE", ["ble"]), ("802.15.4", ["ieee802154"]), ("TSCH", ["tsch"]), ("802.15.6 NB", ["ieee802156nb"]),
          ("SmartBAN", ["smartban"]), ("ANT", ["ant"]), ("IEEE HBC", ["hbc"])]
CLASS_COL = ["#83bf76", "#f0c46d", "#cf6262", "#e4e3de"]                           # class colours
CLASS_LAB = [r"PSR $\geq$ 0.99", "0.90 - 0.99", "< 0.90", "infeasible"]
INK, MUTED, GRID, MEAS = "#0b0b0b", "#898781", "#e1e0d9", "#e6e5df"
# BCC cases: (sweep, reception model, label); the S5 sweeps run the BCC entries get_spec builds for the case
CASES = [("main", "collision", "typical, 500 kchip/s"), ("main", "sinr", "typical, 500 kchip/s"),
         ("s5_bcc_contact", "collision", "contact chain, 50 kchip/s"),
         ("s5_bcc_contact", "sinr", "contact chain, 50 kchip/s"), ("s5_bcc_fcs", "collision", "payload FCS"),
         ("s5_bcc_fcs", "sinr", "payload FCS"), ("s5_bcc_coupling40", "sinr", "inter-body coupling +40 dB"),
         ("s5_bcc_coupling60", "sinr", "inter-body coupling +60 dB"),
         ("s5_bcc_power_low", "collision", "transceiver power low"),
         ("s5_bcc_power_high", "collision", "transceiver power high")]


def sty(sid, lw=1.4, ms=3.2):
    """Line and marker style of a series: the as-built BCC dashed with hollow markers."""
    c = COLOUR[sid]
    ref = sid in ("bcc-thesis", "bcc-bfsk")
    ms *= 1.6 if sid == "bcc-opt" else 1        # a star reads smaller than a disc of the same size
    return dict(color=c, ls="--" if ref else "-", marker=MARKER.get(sid, "o"), ms=ms, mew=0.9,
                mfc="white" if ref else c, lw=lw)


def handle(sid, label=None, **kw):
    return Line2D([], [], label=label or LABEL[sid], **{**sty(sid), **kw})


def band_of(ra):
    return pd.Categorical(np.array(BANDS)[np.searchsorted(EDGES, ra, side="right") - 1], BANDS, ordered=True)


# ---------------------------------------------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------------------------------------------
def build_compact(raw):
    """Raw cells files -> results/cells_main.csv, cells_sensitivity.csv and cells_ladder.csv (4 significant digits)."""
    files = sorted(Path(raw).glob("*_cells.csv"))
    if not files:
        return False
    parts = []
    for f in files:
        d = pd.read_csv(f)
        rows = f.with_name("ladder_contact_rows.csv")
        if f.name == "ladder_contact_cells.csv" and rows.exists():     # T10: retransmissions are a row column only
            keys = ["protocol", "variant", "phy", "reception", "ci_mode", "r_target_kbps", "n_dev", "k"]
            r = pd.read_csv(rows, usecols=keys + ["retx_per_frame"]).groupby(keys).retx_per_frame.mean()
            d = d.merge(r.rename("retx_per_frame_mean").reset_index(), on=keys, how="left")
        d.insert(0, "sweep", f.name[:-len("_cells.csv")])
        vo = d.variant.str.split("+", n=1, expand=True)
        d["variant"], d["option"] = vo[0], (vo[1] if vo.shape[1] > 1 else None)
        parts.append(d)
    d = pd.concat(parts, ignore_index=True)
    d["option"] = d.option.fillna("")
    retx = ["retx_per_frame_mean"] if "retx_per_frame_mean" in d else []
    cols = ["sweep", "protocol", "variant", "option", "phy", "reception", "ci_mode", "r_target_kbps",
            "n_dev", "k"] + KEEP
    d = d[cols + retx].sort_values(cols[:10], kind="mergesort")
    print("report: raw sweeps read:", ", ".join(d.sweep.unique()))
    name = d.sweep.map(lambda w: COMPACT.get(w, "cells_ladder.csv" if w.startswith("ladder") else
                                             "cells_sensitivity.csv"))
    for n in sorted(set(COMPACT.values()) | {"cells_ladder.csv", "cells_sensitivity.csv"}):
        m = name == n
        if m.any():         # a raw directory filled in stages must not blank the other file
            d[m].drop(columns=[] if n == "cells_ladder.csv" else retx).to_csv(
                RES / n, index=False, float_format="%.4g", lineterminator="\n")
        else:
            print(f"report: no sweep of {n} in {raw}, the existing file is kept")
    return True


def load():
    """Every configuration row; sid the stack series, cid the configuration, rc the BCC chip rate, src the sweep
    file. The rates sweep and the proposal rows of the ladder (level 5, pool 0.05) join the main sweep."""
    dt = {"phy": str, "option": str}
    names = [n for n in ["cells_main.csv", "cells_sensitivity.csv", "cells_ladder.csv", "cells_rates.csv",
                         "cells_dense.csv"] if (RES / n).exists()]
    d = pd.concat([pd.read_csv(RES / n, keep_default_na=False, na_values=[""], dtype=dt) for n in names],
                  ignore_index=True)
    d = d[d.variant.isin(["full", "thesis"])].reset_index(drop=True)     # the variants since
    d["option"] = d.option.fillna("")
    bcc = d.protocol == "bcc"
    rc = d.option.str.extract(r"rc_cps=([^+]+)", expand=False).astype(float)
    d["rc"] = rc.fillna(pd.Series(np.where(d.option.str.contains("contact"), 50e3, 500e3), d.index)).where(bcc)
    d["option"] = d.option.str.replace(r"\+?rc_cps=[^+]+", "", regex=True).str.strip("+")
    # the dense sweep runs BCC at 500 kchip/s twice (default and --bcc-rc 500e3), identical rows; keep one
    d = d.drop_duplicates(["sweep", "protocol", "variant", "option", "phy", "reception", "ci_mode", "r_target_kbps",
                           "n_dev", "k", "rc"]).reset_index(drop=True)
    bcc = d.protocol == "bcc"
    d["cid"] = d.protocol + "-" + d.phy
    opt5 = bcc & (d.variant == "full") & (d.option == "opt5")
    # the pooled reference stack: the HBC stack with the retry pool instead of the spare exchange per frame
    pooled = bcc & (d.variant == "full") & d.option.str.contains("pool=") & ~d.option.str.contains("opt")
    d["sid"] = np.select([d.protocol == "ble4thesis", opt5, pooled, bcc & (d.variant == "full"), bcc],
                         [d.cid, "bcc-opt", "bcc-pooled", "bcc-full", "bcc-thesis"], d.protocol)
    d["src"] = d.sweep
    d.loc[d.sweep == "rates", "sweep"] = "main"
    d.loc[opt5 & (d.sweep != "ladder"), "option"] = ""
    d.loc[pooled & (d.sweep == "main"), "option"] = ""
    # the proposal and the pooled reference at 500 kchip/s ran in the ladder sweep: their rows join the main sweep
    lad = d[(opt5 | pooled) & (d.sweep == "ladder")].assign(sweep="main", option="")
    d = pd.concat([d, lad], ignore_index=True)
    for pid, n in GATEWAY.items():
        m = (d.protocol == pid) & (d.n_dev > n)
        d.loc[m, [c for c in KEEP if c != "feasible"]] = np.nan
        d.loc[m, "feasible"] = 0
    d["ra"] = d.n_dev * d.r_target_kbps
    return d


def pick(d, rule, keys=CK):
    """One row per cell (keys): link adaptation over the stack's configurations and ci modes. rule 'psr': the highest
    PSR; 'better': PSRs within TIE tie, the lower sensor J/bit wins; 'latency' (F6): the lower lat_p95."""
    d = d.assign(_f=d.feasible)
    if rule == "psr":
        d = d.sort_values(["_f", "psr_mean", "ci_mode"], ascending=[False, False, True], kind="mergesort")
    else:
        mx = d.groupby(keys).psr_mean.transform("max")
        d = d.assign(_c=(d.psr_mean >= mx - TIE).astype(int))
        m = "lat_p95_ms_mean" if rule == "latency" else "jpb_sensor_mean"
        d = d.sort_values(["_f", "_c", m, "ci_mode"], ascending=[False, False, True, True], kind="mergesort")
    return d.groupby(keys, sort=False).head(1).drop(columns=[c for c in ["_f", "_c"] if c in d])


def points(sel, cols):
    """Median over N_dev at each (series, k, R_agg) point; feasible = any N_dev feasible."""
    g = sel.groupby(SK + ["k", "ra"], observed=True)
    out = g[cols].median()
    out["feasible"] = g.feasible.max()
    return out.reset_index()


def bands(pts, cols):
    """Median over the R_agg points of each band, with the counts of feasible and grid points."""
    p = pts.assign(band=band_of(pts.ra))
    g = p.groupby(SK + ["band", "k"], observed=True)
    out = g[cols].median()
    out["n_feas"], out["n_pts"] = g.feasible.sum(), g.size()
    return out.reset_index()


def rows_of(df):
    """Series present in df, in SERIES order."""
    have = set(df.sid)
    return [s for s, _ in SERIES if s in have]


# ---------------------------------------------------------------------------------------------------------------
# working points
# ---------------------------------------------------------------------------------------------------------------
def read_working_points():
    """Rows of the working-point table: (id, application, ra_lo, ra_hi, k_lo, k_hi, n_lo, n_hi)."""
    wps = []
    for line in WORKING_POINTS.splitlines():
        m = re.match(r"\s*([A-H]) & (.+?)\s*&", line)
        if not m:
            continue
        c = [re.findall(r"\d+(?:\.\d+)?", x) for x in line.split("&")]
        n, ra, k = c[4], c[5], c[6]
        app = re.sub(r"\\,", " ", m.group(2)).replace("$", "")
        wps.append((m.group(1), app, float(ra[0]), float(ra[-1]), int(k[0]), int(k[-1]), int(n[0]),
                    int(n[-1])))
    assert [w[0] for w in wps] == list("ABCDEFGH"), wps
    return wps


def wp_ndev(grid, n_lo, n_hi):
    inside = [n for n in grid if n_lo <= n <= n_hi]
    if inside:
        return inside
    dist = {n: max(n_lo - n, n - n_hi) for n in grid}
    return [n for n in grid if dist[n] == min(dist.values())]


def wp_cells(sel, wp):
    """The rows of sel that match working point wp, and the grid N_dev values used."""
    _, _, ra_lo, ra_hi, k_lo, k_hi, n_lo, n_hi = wp
    ns = wp_ndev(sorted(sel.n_dev.unique()), n_lo, n_hi)
    return sel[sel.n_dev.isin(ns) & sel.ra.between(ra_lo, ra_hi) & sel.k.between(k_lo, k_hi)], ns


def no_inf(m):
    """J/bit columns with inf (no delivered bit) as NaN, so a median runs over the cells that deliver."""
    jc = [c for c in m if c.startswith("jpb")]
    return m.assign(**{c: m[c].replace([np.inf, -np.inf], np.nan) for c in jc})


def wp_table(sel, wp, cols):
    """Median over the feasible matching cells per series (sel: one row per cell); J/bit over the cells with delivered
    bits (n_zero: feasible cells that deliver none, J/bit inf)."""
    m, ns = wp_cells(sel, wp)
    m = no_inf(m.assign(lat=((m.ci_mode == "latency") & (m.feasible == 1)).astype(int),
                        n_zero=((m.feasible == 1) & np.isinf(m.jpb_sensor_mean)).astype(int)))
    g = m.groupby(SK, observed=True)
    out = g[cols].median()
    out["n_feas"], out["n_cells"], out["n_lat"], out["n_zero"] = g.feasible.sum(), g.size(), g.lat.sum(), g.n_zero.sum()
    return out.reset_index(), ns


def self_check():
    """Working-point rule on a constructed frame (values derived by hand)."""
    rows = []
    for n in [1, 4, 10]:
        for r in [5, 10, 100]:
            for k in [1, 2, 4]:
                # n 1: PSR 1.0 (latency) vs 0.9995 (energy) tie within TIE -> energy (J/bit r < 100)
                # n 4: PSR 1.0 vs 0.98 -> latency (J/bit 100); energy J/bit is r
                e_psr = 0.9995 if n == 1 else 0.98
                rows += [dict(ci_mode="latency", psr_mean=1.0, jpb_sensor_mean=100.0),
                         dict(ci_mode="energy", psr_mean=e_psr, jpb_sensor_mean=float(r))]
                for x in rows[-2:]:
                    x.update(sweep="t", sid="x", option="", reception="collision",
                             r_target_kbps=r, n_dev=n, k=k, feasible=1, ra=n * r)
    t = pd.DataFrame(rows)
    t.loc[(t.n_dev == 1) & (t.r_target_kbps == 5) & (t.k == 2), "feasible"] = 0      # one infeasible cell
    t.loc[t.feasible == 0, ["psr_mean", "jpb_sensor_mean"]] = np.nan
    sel = pick(t, "better")
    assert len(sel) == 27
    # N_dev range 2-3: grid 1 and 4 are both 1 away -> both; R_agg 5-50: n1 r 5, 10; n4 r 5, 10; k 1-2 -> 8 cells
    out, ns = wp_table(sel, ("T", "test", 5, 50, 1, 2, 2, 3), ["jpb_sensor_mean"])
    assert ns == [1, 4]
    o = out.iloc[0]
    assert (o.n_cells, o.n_feas, o.n_lat) == (8, 7, 4)
    # n1 (energy): r5 k1 -> 5, r5 k2 infeasible, r10 k1 and k2 -> 10, 10; n4 (latency): 100 x 4
    # sorted [5, 10, 10, 100, 100, 100, 100] -> median 100, 7 feasible of 8, 4 in latency mode
    assert o.jpb_sensor_mean == 100.0
    out, ns = wp_table(sel, ("T", "test", 0, 20, 1, 1, 1, 1), ["jpb_sensor_mean"])    # n1: r 5, 10 -> 5, 10
    assert ns == [1] and out.iloc[0].jpb_sensor_mean == 7.5 and out.iloc[0].n_lat == 0
    assert wp_ndev([1, 4, 10], 5, 8) == [4] and wp_ndev([1, 4, 10], 10, 10) == [10]
    # class limits: 0.98996 would print as 0.990 and hide that it is below 0.99; 0.9901 is at or above
    assert (fpsr(0.98996), fpsr(0.9899), fpsr(0.9901), fpsr(0.8999)) == ("0.98996", "0.9899", "0.990", "0.8999")
    # capacity: the scan down from the bound finds the maximum of the definition (every feasible per-body rate of
    # the 0.5 kb/s grid, scanned upward) at every N_dev; as-built BCC at 500 kchip/s, one tone at 50 kchip/s
    for sid, rc in [("bcc-bfsk", 500e3), ("bcc-sf8", 50e3)]:
        s = bcc_spec(sid, rc)
        for n in [1, 4, 10]:
            ok = [i * 0.5 for i in range(1, 400) if any(protocols.select_period(i * 0.5 / n, m, n, s)
                                                        for m in ("latency", "energy"))]
            assert capacity(sid, rc, n) == max(ok) < 199, (sid, n, max(ok))
    # the ladder scan (pool bound, memoised layouts) against select_period itself: its capacity is admitted, the next
    # 0.5 kb/s is not (one tone, 50 kchip/s, L3, pool 0.43; the full upward scan takes 15 s)
    s, c = bcc_spec("bcc-sf8", 50e3, 3, 0.43), capacity("bcc-sf8", 50e3, 1, 3, 0.43)
    ok = [any(protocols.select_period(x, m, 1, s) for m in ("latency", "energy")) for x in (c, c + 0.5)]
    assert ok == [True, False], (c, ok)
    # hub limits of T14's 'why': GATEWAY (BLE), else protocols max_dev (802.15.4 CSMA: none)
    assert [hub_limit(*c) for c in [("ble", "1M"), ("smartban", "bch1"), ("ant", "gfsk"), ("ieee802154", "gts"),
                                    ("ieee802154", "oqpsk")]] == [20, 16, 15, 7, math.inf]
    # the pool bound of capacity() at 500 kchip/s: the one short frame of a period saves (against L_max's cost per
    # octet) less than the beacon's share per sensor, whose ack field is smallest with one frame per sensor; every
    # level with a pool, the pooled reference (L0, pool 0.05) included
    for sid in LCID:
        for lv, pool in [(lv, p) for lv, p in LADDER if lv >= 3 or p]:
            s, n = bcc_spec(sid, 500e3, lv, float(pool) if pool else None), 10
            lmax = s.max_payload_octets
            save = max(o * frame_cost(s, lmax) / lmax - frame_cost(s, o) for o in range(1, lmax + 1))
            bcn = math.ceil(protocols.beacon_layout(s, n, 1)[-1]["t1_us"] - 1e-6) + protocols.HBC["psifs_us"]
            assert save < bcn / n, (sid, lv, pool, save, bcn / n)


# ---------------------------------------------------------------------------------------------------------------
# formatting
# ---------------------------------------------------------------------------------------------------------------
def md(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(x) for x in r) + " |" for r in rows]
    return "\n".join(out) + "\n"


def f3(x):
    return "" if pd.isna(x) else f"{x:.3f}"


def fe(x, d=2):
    if pd.isna(x):
        return ""
    return "inf" if np.isinf(x) else ("0" if x == 0 else f"{x:.{d}e}")


def fg(x):
    """3 significant digits, trailing zeros kept (39.0, not 39); whole numbers from 1000 up."""
    if pd.isna(x):
        return ""
    return "inf" if np.isinf(x) else (f"{x:.0f}" if abs(x) >= 1000 else f"{x:#.3g}".rstrip("."))


def fpsr(x):
    """PSR with 3 decimals, more where rounding would hide that it lies below a class limit (0.99, 0.90)."""
    s = f3(x)
    if s in ("0.990", "0.900") and x < float(s):
        s = next(t for t in (f"{x:.{n}f}" for n in (4, 5, 6, 12)) if float(t) < float(s) or len(t) == 14)
    return s


def fd(x):
    return "" if pd.isna(x) else f"{x:+.3f}".replace("-0.000", "+0.000")


def ratio(p, col):
    """Median over paired cells of the ratio a / b (finite pairs)."""
    r = (p[col + "_a"] / p[col + "_b"]).replace([np.inf, -np.inf], np.nan).dropna()
    return f"{r.median():.2f}" if len(r) else ""


def write(name, text):
    (TAB / name).write_text(text, encoding="ascii", newline="\n")


# ---------------------------------------------------------------------------------------------------------------
# BCC stacks from protocols.py
# ---------------------------------------------------------------------------------------------------------------
@lru_cache(maxsize=None)
def bcc_spec(sid, rc, opt=0, pool=None, spread=None, lmax=None, back=False, burst=False, span=None, pilot128=False):
    """BCC stack of a configuration (cid) at chip rate rc; opt, pool: the ladder level and retry-pool share (None:
    0.05) of BCC full (DESIGN 4.3); bcc-opt is multi-tone at level 5; spread: the payload spreading; lmax:
    the frame payload limit of one-tone BCC full (claim4.py); back: the block-ack reference (claim4.py H), burst: in
    the HBC burst mode, span: its block size across allocations (claim4.py J); pilot128: the '128 B' pilot period at
    level 0 (the lead reference)."""
    if sid == "bcc-opt":
        sid, opt = "bcc-multitone", 5
    var, phy = {"bcc-sf8": ("full", "sf8"), "bcc-multitone": ("full", "multitone"), "bcc-bfsk": ("thesis", None)}[sid]
    flags = {k: v for k, v in (("opt", opt), ("pool", pool), ("spread", spread), ("lmax", lmax), ("back", back),
                                       ("burst", burst), ("span", span), ("pilot128", pilot128)) if v}
    return protocols.get_spec("bcc", var, phy, rc_cps=rc, **flags)


def frame_cost(s, octets):
    """Schedule time (us) of a frame in capacity()'s bound: its allocation, with the pool (from ladder level 3, or the
    pooled reference stack) plus the pool's share of a pool exchange (the pool holds at least that share of the
    superframe's frames). Under the flag back the allocation's B-Ack is left out (claim4.back_cost with an unlimited
    block), so the bound stays above what select_period admits with one B-Ack per allocation; under burst the frame
    is a later frame of the burst (pMIFS, RI mode), the limit of an unlimited burst."""
    x = protocols._frame_us(s, octets, cont=True)
    if protocols._pooled(s):
        x += s.flags.get("pool", protocols.POOL_SHARE) * protocols._pool_us(s, octets)
    return x


@lru_cache(maxsize=None)
def capacity(sid, rc, n_dev, opt=0, pool=None, spread=None, lmax=None, back=False, burst=False, span=None,
             pilot128=False):
    """Largest per-body rate n_dev x R (kb/s) in 0.5 kb/s steps that select_period admits in either ci mode (0: none).
    The scan runs down from a bound, so the first feasible rate is the largest: 8 L_max / frame_cost(L_max), up to
    ladder level 2 select_period's own bound (N_dev max-size frames' share of the medium <= 1). From level 3 it adds
    the pool share, which spares the scan the rates between the two bounds (each a full scan of the period grid).
    The cost per octet of a frame is a sawtooth in its length (some frames shorter than L_max cost less), but a period
    carries only one short frame per sensor, the last, and what it saves stays below the beacon and pSIFS share per
    sensor (self_check asserts this at every level with a pool and N_dev 10, the largest grid value)."""
    s = bcc_spec(sid, rc, opt, pool, spread, lmax, back, burst, span, pilot128)
    lmax = s.max_payload_octets
    i = math.floor(8000 * lmax / frame_cost(s, lmax) / 0.5)
    # an infeasible rate lays out the frames of every period of the grid (multi-tone: a stripe balance each): the
    # scan memoises the layouts of s by octets while it runs
    real, memo = {f: getattr(protocols, f) for f in ("_ex_us", "airtime_us")}, {}

    def memoised(f):
        def g(spec, octets, *a):
            if spec is not s:
                return real[f](spec, octets, *a)
            if (f, octets, a) not in memo:
                memo[f, octets, a] = real[f](s, octets, *a)
            return memo[f, octets, a]
        return g
    try:
        for f in real:
            setattr(protocols, f, memoised(f))
        while i > 0 and not any(protocols.select_period(i * 0.5 / n_dev, m, n_dev, s) for m in ("energy", "latency")):
            i -= 1
    finally:
        for f, fn in real.items():
            setattr(protocols, f, fn)
    return i * 0.5


# ---------------------------------------------------------------------------------------------------------------
# tables
# ---------------------------------------------------------------------------------------------------------------
def t1_t2(land):
    for rec in ["collision", "sinr"]:
        L = land[land.reception == rec]
        rows = rows_of(L)
        parts = [f"PSR landscape, {rec} model: median over the R_agg grid points of the band of (median over N_dev "
                 "of (the highest PSR over the stack's configurations and the two ci modes: link adaptation)). "
                 "pts = grid points with a feasible cell / grid points in "
                 "the band; * = not every point feasible; n/f = no feasible point. A value is printed with 4 "
                 "decimals where 3 would round it up onto a class limit (0.99, 0.90) it lies below. " + RATE_NOTE +
                 ".\n"]
        for b in BANDS:
            B = L[L.band == b]
            tab = []
            for s in rows:
                x = B[B.sid == s].set_index("k")
                if x.empty:
                    continue
                nf, npt = int(x.n_feas.iloc[0]), int(x.n_pts.iloc[0])
                star = "*" if 0 < nf < npt else ""
                tab.append([LABEL[s], f"{nf}/{npt}"] +
                           [("n/f" if nf == 0 else fpsr(x.psr_mean.get(k)) + star) for k in KS])
            if tab:
                parts.append(f"\nR_agg {b} kb/s\n\n" + md(["series", "pts"] + [f"k = {k}" for k in KS], tab))
        write(f"t1_psr_landscape_{rec}.md", "".join(parts))
        tab = []
        for s in rows:
            cells = []
            for b in BANDS:
                x = L[(L.sid == s) & (L.band == b)].set_index("k")
                if x.empty or x.n_feas.iloc[0] == 0:
                    cells.append("n/f")
                    continue
                star = "*" if x.n_feas.iloc[0] < x.n_pts.iloc[0] else ""
                ks = [next((str(k) for k in KS if x.psr_mean.get(k, np.nan) < t), "-") for t in (0.99, 0.90)]
                cells.append(f"{ks[0]} / {ks[1]}{star}")
            tab.append([LABEL[s]] + cells)
        write(f"t2_k_thresholds_{rec}.md",
              f"First grid k (1, 2, 4, 8, 16, 30) at which the band PSR of T1 ({rec} model) falls below 0.99 / "
              "below 0.90; - = not below up to k 30; n/f = no feasible point; * = not every grid point of the band "
              "feasible. " + RATE_NOTE + ".\n\n" + md(["series"] + [f"{b}" for b in BANDS], tab))


def t3(good, wps):
    cols = ["psr_mean", "on_time_mean", "lat_p95_ms_mean", "jpb_sensor_mean", "jpb_ban_mean", "lifetime_days_mean"]
    res = {}
    for rec in ["collision", "sinr"]:
        sel = good[(good.sweep == "main") & (good.reception == rec)]
        summ, parts, parts_b = {}, [], []
        for wp in wps:
            out, ns = wp_table(sel, wp, cols)
            res[(rec, wp[0])] = out
            wid, app, ra_lo, ra_hi, k_lo, k_hi, n_lo, n_hi = wp
            m = wp_cells(sel, wp)[0].query("feasible == 1")
            k_top = max(k for k in KS if k_lo <= k <= k_hi)
            top = m[m.k == k_top].groupby("sid").psr_mean.median()
            tab = []
            for s in rows_of(out):
                o = out[out.sid == s].iloc[0]
                if o.n_feas == 0:
                    tab.append([LABEL[s], "n/f", "", "", "", "", "", "", f"0/{o.n_cells}", ""])
                    continue
                z = f" [{o.n_zero} none]" if o.n_zero else ""
                tab.append([LABEL[s], f3(o.psr_mean), f3(top.get(s)) or "n/f", f3(o.on_time_mean),
                            fg(o.lat_p95_ms_mean), fe(o.jpb_sensor_mean) + z, fe(o.jpb_ban_mean) + z,
                            fg(o.lifetime_days_mean), f"{o.n_feas}/{o.n_cells}", f"{o.n_lat}/{o.n_feas - o.n_lat}"])
                summ.setdefault(s, {})[wp[0]] = fe(o.jpb_sensor_mean) + z
            nr = f"{n_lo}" if n_lo == n_hi else f"{n_lo}-{n_hi}"
            title = (f"\n{wid} {app}: R_agg {ra_lo:g}-{ra_hi:g} kb/s, k {k_lo}-{k_hi}, N_dev {nr} "
                     f"(grid N_dev {', '.join(map(str, ns))})\n\n")
            parts.append(title + md(["series", "psr", f"psr k {k_top}", "on_time", "lat_p95 (ms)", "J/bit sensor",
                                     "J/bit BAN", "lifetime (d)", "cells", "ci lat/en"], tab))
            # T3b: the same cells split by R_agg grid value (the medians above mix per-link rates)
            mm = no_inf(wp_cells(sel, wp)[0])
            ras = sorted(mm.ra.unique())
            heads = [f"R_agg {ra:g} (" + ", ".join(f"{n} x {r:g}" for n, r in
                                                    sorted(set(zip(*[mm[mm.ra == ra][c] for c in
                                                                     ["n_dev", "r_target_kbps"]])))) + ")"
                     for ra in ras]
            tab = []
            for s in rows_of(mm):
                row = [LABEL[s]]
                for ra in ras:
                    x = mm[(mm.sid == s) & (mm.ra == ra)]
                    f = x[x.feasible == 1]
                    row.append(f"n/f (0/{len(x)})" if f.empty else
                               f"{f3(f.psr_mean.median())} / {f3(f.on_time_mean.median())} / "
                               f"{fe(f.jpb_sensor_mean.median())} / {fg(f.lifetime_days_mean.median())} "
                               f"({len(f)}/{len(x)})")
                tab.append(row)
            parts_b.append(title + md(["series"] + heads, tab))
        head = (f"Working points, {rec} model. Rule: cells whose R_agg = N_dev x R and k lie in the "
                "point's ranges and whose N_dev is the grid value (1, 4, 10) inside the point's N_dev range, else the "
                "nearest; per cell the stack's best configuration and ci mode (link adaptation: higher PSR; within "
                "0.005 the lower sensor J/bit); median over "
                "the feasible matching cells. Each value is the median of its own column, not the values of one "
                "cell, and each series runs over its own feasible cells (cells = feasible / matching): a series "
                "infeasible in the harder cells of a point is judged on the easier ones (T3b splits the cells by "
                "R_agg). psr k N: the median over the feasible cells at the top grid k of the point's k range. ci "
                "lat/en = cells in latency / energy mode; n/f = no feasible cell. J/bit per delivered application "
                "bit, the median over the feasible cells that deliver bits; [n none]: n feasible cells deliver none "
                f"(J/bit infinite, left out of the median; T3b likewise). {RATE_NOTE}; T6b: the BCC series at 500 "
                f"kchip/s. {T3_PROFILE_NOTE}\n\nSummary: sensor J per delivered bit\n\n")
        tab = [[LABEL[s]] + [summ[s].get(w[0], "n/f") for w in wps]
               for s in rows_of(sel) if s in summ]
        write(f"t3_working_points_{rec}.md", head + md(["series"] + [w[0] for w in wps], tab) + "".join(parts))
        write(f"t3b_working_points_by_ragg_{rec}.md",
              f"T3b, {rec} model: the cells of T3 (same matching rule, better ci mode per cell) split by the R_agg "
              "grid value (column head: R_agg (N_dev x R per link, kb/s)). Per series and column: psr / on_time / "
              "sensor J per delivered bit (over the cells that deliver bits) / lifetime (d), each the median of its "
              "own column over the feasible cells "
              "(the k values of the point's range), then (feasible / matching cells); n/f = no feasible cell. " +
              RATE_NOTE + ". " + T3_PROFILE_NOTE + "\n" +
              "".join(parts_b))
    return res


def t4(good):
    s1r = sorted(good[good.sweep == "s1_nocifloor"].r_target_kbps.unique())
    run = good[(good.reception == "sinr") & good.sweep.isin(["main", "s1_nocifloor"]) & good.r_target_kbps.isin(s1r)]
    sel = run[run.feasible == 1]
    tab = []
    for s in rows_of(sel[sel.sweep == "main"]):
        r = [LABEL[s]]
        for k in [8, 30]:
            for sw in ["main", "s1_nocifloor"]:
                x = sel[(sel.sid == s) & (sel.sweep == sw) & (sel.k == k)]
                r += [f"{fe(x.raw_ber_mean.mean(), 1)} / {fe(x.post_ber_mean.mean(), 1)}" if len(x) else
                      "n/f" if ((run.sid == s) & (run.sweep == sw) & (run.k == k)).any() else "n/a (not run)"]
        tab.append(r)
    write("t4_ber_sinr.md",
          f"Payload BER raw / post-FEC in the SINR model, mean over the feasible cells of the S1 grid (R "
          f"{', '.join(f'{r:g}' for r in s1r)} kb/s, every N_dev, both ci modes reduced to the better one), with the "
          "C/I floor of each receiver (main) and without it (S1); each stack at its best configuration per cell "
          "(BCC at 500 "
          "kchip/s in both, the rates sweep ran the main grid only). FEC: 802.15.6 NB and SmartBAN BCH (the "
          "uncoded SmartBAN schemes have none), BLE Coded (convolutional, modelled as a coding gain on the SNR, so "
          "raw = post by construction); raw = post for the others (802.15.4 after despreading, HBC and the BCC "
          "stacks after Walsh demapping). n/f: no feasible cell; n/a (not run): the sweep has no row of the series "
          "(the proposal: the ladder sweep ran the main grid only).\n\n" +
          md(["series", "k 8 floor", "k 8 no floor", "k 30 floor", "k 30 no floor"], tab))


def t5(good):
    cols = ["collision_frac_mean", "jain_fairness_mean", "goodput_kbps_mean"]
    sel = good[good.sweep.isin(MAIN_SWEEPS)]
    B = bands(points(sel, cols), cols)
    parts = ["Medium metrics, band medians (better ci mode per cell, median over N_dev, median over the band's "
             "R_agg points): collision_frac (attempts overlapping another same-channel transmission), Jain fairness "
             "of PSR over the k BANs, goodput per body (kb/s). *: fewer feasible grid points than the band has (as in "
             "T1), so the goodput median sits below the band's offered load. " + RATE_NOTE + ".\n"]
    for rec in ["collision", "sinr"]:
        for b in ["30-100", "300-500"]:
            x = B[(B.reception == rec) & (B.band == b) & (B.n_feas > 0)]
            tab = []
            for s in rows_of(x):
                y = x[x.sid == s].set_index("k")
                star = "*" if (y.n_feas < y.n_pts).any() else ""
                tab.append([LABEL[s] + star] + [f(y[c].get(k)) for k in [8, 30] for c, f in zip(cols, [f3, f3, fg])])
            parts.append(f"\n{rec} model, R_agg {b} kb/s\n\n" + md(
                ["series", "k 8 coll_frac", "k 8 Jain", "k 8 goodput", "k 30 coll_frac", "k 30 Jain", "k 30 goodput"],
                tab))
    write("t5_medium.md", "".join(parts))


def bcc_rows(d, sid, cid, rc=None):
    """Collision-model rows of one BCC configuration in the main sweep (rc None: every chip rate)."""
    x = d[(d.sweep == "main") & (d.sid == sid) & (d.cid == cid) & (d.reception == "collision")]
    return x if rc is None else x[x.rc == rc]


def t6(d, good5, wps):
    tab = []
    for sid, cid, lv, pool, lab in BCFG:
        s = bcc_spec(cid, 500e3, lv, pool)
        a = bcc_rows(d, sid, cid, 500e3)
        f = a[a.feasible == 1]
        bcn = s.beacon["reserved_us"]
        if lv >= 2 or protocols._pooled(s):      # select_period's reserve: the ladder beacon (ack field from L2)
            bcn = math.ceil(protocols.beacon_layout(s, 1, 1)[-1]["t1_us"] - 1e-6) + protocols.HBC["psifs_us"]
        cap = [capacity(cid, 500e3, n, lv, pool) for n in (1, 4, 10)]
        tab.append([lab, len(s.tones) if s.tones else 1,
                    " / ".join(fg(protocols.airtime_us(s, n) / 1e3) for n in (20, 255)),
                    " / ".join(fg(sum(st.dur_us for st in protocols.exchange(s, n)) / 1e3) for n in (20, 255)),
                    fg(protocols._frame_us(s, 255) / 1e3), fg(bcn / 1e3),
                    " / ".join(f"{c:g}" if c else "< 0.5" for c in cap),
                    f"{f.ra.max():g}" if len(f) else "n/f" if len(a) else "not run"])
    out = ["BCC: what the IEEE 802.15.6 HBC features, the parallel tones and the ladder cost and buy the as-built stack, "
           "at 500 kchip/s (T12: every chip rate). BCC full carries the HBC stack (SF 64 header, MAC header and FCS, "
           "I-Ack, 3 retransmissions, allocations of one exchange + one spare per frame; the spare: this model's "
           "provisioning rule of 2026-09-24, not an HBC feature); the pooled reference (, the reference of Claim "
           "4) is the same HBC stack with one exchange per frame in the allocation and its retries in a hub-polled "
           "pool of 5 % of the superframe's frames after the allocations, a hub policy the standard permits; "
           "multi-tone stripes each "
           "frame over the tones of the plan (one BCC full packet per tone, I-Ack on every tone; bcc_fdma.plan: "
           "spacing 4 R_c over the 3 dB span of the cohort 50 Ohm curve, DESIGN 4.1); the proposal is BCC full at "
           "ladder level 5 with the 0.05 retry pool (DESIGN 4.3, T10), one tone or multi-tone; the as-built stack is "
           "the as-built stack as coded (no FCS, no ACK, hub sync, TDMA slot with guard). The series BCC full and the "
           "proposal take, per cell, the better of their one-tone and multi-tone configurations at any chip rate.\n\n"
           "(a) From protocols.py. Air time of a frame with 20 and 255 application octets (multi-tone: striped over "
           "its tones, the longest tone packet), one exchange (frame, pSIFS, I-Ack, pSIFS; as built: frame and slot "
           "guard; the proposal: data and pSIFS, the beacon acknowledges), the allocation per 255-octet frame (the "
           "pooled reference and the proposal: one exchange, a shared retry pool follows the allocations), the beacon "
           "(as built: sync; the "
           "proposal: with its acknowledgement field for N_dev 1 and one frame) reserved per period; capacity per body "
           "= the largest N_dev x R that select_period admits in either ci mode (0.5 kb/s steps; < 0.5: none); largest "
           "feasible grid R_agg of the sweep, collision model (main; the pooled reference and the proposal: the ladder "
           "sweep; not run: no row).\n\n",
           md(["configuration", "tones", "frame air time L 20 / 255 (ms)", "exchange L 20 / 255 (ms)",
               "allocation per 255-octet frame (ms)", "beacon (ms)", "capacity per body N_dev 1 / 4 / 10 (kb/s)",
               "largest feasible grid R_agg (kb/s)"], tab)]
    tab = []
    sel = good5[(good5.sweep == "main") & (good5.reception == "collision")]
    for wp in wps:
        o = wp_table(sel, wp, ["psr_mean", "on_time_mean", "lat_p95_ms_mean", "jpb_sensor_mean"])[0]
        row = [f"{wp[0]} ({wp[2]:g}-{wp[3]:g} kb/s, k {wp[4]}-{wp[5]})"]
        for sid in BCC:
            x = o[(o.sid == sid) & (o.n_feas > 0)]
            row.append("n/f" if x.empty else
                       f"{f3(x.psr_mean.iloc[0])} / {f3(x.on_time_mean.iloc[0])} / {fg(x.lat_p95_ms_mean.iloc[0])} / "
                       f"{fe(x.jpb_sensor_mean.iloc[0])} ({x.n_feas.iloc[0]}/{x.n_cells.iloc[0]})")
        tab.append(row)
    out += ["\n(b) Working points at 500 kchip/s, collision model (T3 rule: the series' best configuration and ci "
            "mode per cell at 500 kchip/s, median over the feasible matching cells; T3 takes each series at its best "
            "chip rate): psr / on_time / lat_p95 (ms) / sensor J per delivered bit "
            "(feasible / matching cells); n/f = no feasible cell.\n\n", md(["point"] + [LABEL[s] for s in BCC], tab)]
    cols = ["jpb_sensor_mean", "on_time_mean", "lat_p95_ms_mean"]
    B = bands(points(good5[good5.sid.isin(BCC)], cols), cols)
    B = B[(B.k == 1) & (B.n_feas > 0)]
    allb = d[d.sid.isin(BCC) & (d.src != "rates")]         # every configuration and ci mode, 500 kchip/s in main
    tab = []
    for sw, rec, lab in CASES:
        for sid in BCC:
            a = allb[(allb.sweep == sw) & (allb.reception == rec) & (allb.sid == sid)]
            if a.empty:
                continue
            f = a[a.feasible == 1]
            bv = B[(B.sweep == sw) & (B.reception == rec) & (B.sid == sid)].set_index("band")
            row = [lab, rec, LABEL[sid], f"{len(f)}/{len(a)}"]
            if f.empty:
                tab.append(row + ["n/f"] + [""] * 8)
                continue
            w = f.loc[f.psr_mean.idxmin()]
            at_ = (f"{w.ci_mode}, {w.phy}, R {w.r_target_kbps:g}, N_dev {w.n_dev}, k {w.k} [{f3(w.psr_lo)}, "
                   f"{f3(w.psr_hi)}]" if w.psr_mean < 1 else "every row")
            tab.append(row + [f3(w.psr_mean), at_, str(int((f.psr_mean < 0.995).sum())),
                              f3(f[f.k == 30].psr_mean.median()) or "no k 30"] +
                       [fe(bv.jpb_sensor_mean.get(b)) for b in ["1-10", "10-30", "100-300"]] +
                       [f3(bv.on_time_mean.get("10-30")), fg(bv.lat_p95_ms_mean.get("10-30"))])
    out += ["\n(c) Cases per BCC series (S5; typical = the main sweep at 500 kchip/s). Feasible rows: "
            "configuration and "
            "ci-mode rows with a feasible period / all rows of the series in that sweep. min psr: every feasible row, "
            "configuration, k and ci mode (at: where, with its 95 % interval); rows < 0.995: feasible rows below "
            "0.995; psr k 30: median over the feasible rows at k 30; J/bit, on_time, lat_p95: band medians at k 1 "
            "(best configuration and ci mode per cell, median over N_dev at each R_agg point, then over the band's "
            "points; they mix sensors whose per-link rates differ tenfold). A missing row: get_spec refuses the case "
            "for that configuration (multi-tone: contact chain and coupling; BCC full: the FCS switch, it carries the "
            "802.15.6 FCS), or, for the proposal, the case was not run (the ladder sweep ran the main grid). The "
            "coupling cases give the BCC receivers a co-channel C/I floor of 0 dB (a chip under SIR < 0 dB is a guess, "
            "PARAMETERS 3.8).\n\n",
            md(["case", "model", "series", "feasible rows", "min psr", "at", "rows < 0.995", "psr k 30", "J/bit 1-10",
                "J/bit 10-30", "J/bit 100-300", "on_time 10-30", "lat_p95 10-30 (ms)"], tab)]
    c = d[(d.sweep == "s5_bcc_contact") & (d.feasible == 1)]          # both ci modes: frame length matters
    tab = []
    for sid in BCC:
        for rec in ["collision", "sinr"]:
            x = c[(c.sid == sid) & (c.reception == rec)]
            if x.empty:
                continue
            lng, sht = x[x.payload_bits >= 2000], x[x.payload_bits < 2000]
            tab.append([LABEL[sid], rec, f3(sht.psr_mean.median()),
                        f"{sht.psr_mean.min():.3f}-{sht.psr_mean.max():.3f}" if len(sht) else "",
                        f3(lng.psr_mean.median()) or "n/f",
                        f"{35 / 48 * np.exp(-8448 / 35000):.3f}" if sid == "bcc-thesis" else "", "0.66-0.72"])
    out += ["\n(d) Contact case (measured-hardware two-state contact chain, 50 kchip/s), all feasible cells, k and ci "
            "modes, the single-tone configurations (no per-tone chain for multi-tone). The chain and the retry pool "
            "0.43 are the pessimistic case, calibrated on the HackRF proof of concept: the "
            "author's "
            "real-time DBFSK audio stream over BCC (>= 44.1 kb/s, no audible errors) shows the HackRF packet rates are "
            "not the channel limit. Short frames: payload < 2000 bits; 250-octet frames: payload >= 2000 bits (energy "
            "mode). Chain formula 35/48 exp(-N / 35000) with N = 8448 decoded chips of an as-built 250-octet frame; "
            "measured: HackRF on the arm.\n\n",
            md(["series", "model", "psr short frames (median)", "range", "psr 250-octet frames (median)",
                "chain formula 250 octets", "measured"], tab)]
    write("t6_bcc.md", "".join(out))


def pair(a, b, keys):
    m = a.merge(b, on=keys, suffixes=("_a", "_b"))
    return m[(m.feasible_a == 1) & (m.feasible_b == 1)]


def t7(d):
    keys = ["reception", "ci_mode", "r_target_kbps", "n_dev", "k"]
    main = d[(d.sweep == "main") & (d.rc == 500e3)]
    lab = {(c[0], c[1]): c[-1] for c in BCFG}
    tab = []
    for a, b in [(("bcc-full", "bcc-sf8"), ("bcc-thesis", "bcc-bfsk")),
                 (("bcc-full", "bcc-multitone"), ("bcc-full", "bcc-sf8")),
                 (("bcc-opt", "bcc-sf8"), ("bcc-full", "bcc-sf8")),
                 (("bcc-opt", "bcc-multitone"), ("bcc-full", "bcc-multitone"))]:
        for rec in ["collision", "sinr"]:
            for mode in ["latency", "energy"]:
                x = main[(main.reception == rec) & (main.ci_mode == mode)]
                A, B = [x[(x.sid == s) & (x.cid == c)] for s, c in (a, b)]
                p = pair(A, B, keys)
                tab.append([f"{lab[a]} - {lab[b]}", rec, mode,
                            f"{int(A.feasible.sum())} / {int(B.feasible.sum())} of {max(len(A), len(B))}", str(len(p)),
                            fd((p.psr_mean_a - p.psr_mean_b).mean()), fd((p.on_time_mean_a - p.on_time_mean_b).mean()),
                            ratio(p, "lat_p95_ms_mean"), ratio(p, "jpb_sensor_mean"), ratio(p, "jpb_ban_mean")])
    write("t7_bcc_pairs.md",
          "Paired cells of BCC configurations, main sweep at 500 kchip/s: BCC full on one tone minus the as-built "
          "stack (what the IEEE 802.15.6 HBC features cost and buy), multi-tone minus one tone (the same stack on the "
          "tones of the plan) and the proposal (ladder level 5, pool 0.05) minus BCC full on the same tones "
          "(T10). A pair: the same model, ci mode, R, N_dev and k, feasible in both configurations, each at its own "
          "period. feasible cells: first / second out of the grid (the capacity cost); dPSR and don_time: mean over "
          "the pairs of first - second; lat_p95 and J/bit: median over the pairs of the ratio first / second (finite "
          "pairs). Every k pooled: the isolated BCC bodies do not interact.\n\n" +
          md(["pair", "model", "ci mode", "feasible cells", "pairs", "dPSR", "don_time", "lat_p95 ratio",
              "J/bit sensor ratio", "J/bit BAN ratio"], tab))


def t8(d):
    keys = ["sid", "reception", "ci_mode", "r_target_kbps", "n_dev", "k"]
    d = pick(d[d.src != "rates"], "better", CK + ["ci_mode"])     # best configuration per ci mode, BCC 500 kchip/s
    main = d[d.sweep == "main"]
    runs = [("s1_nocifloor", "sinr", "S1 no C/I floor"), ("s3_gym", "sinr", "S3 gym, high density"),
            ("s4_case_low", "sinr", "S4 low inter-body loss"), ("s4_case_high", "sinr", "S4 high inter-body loss")]
    res = {(sw, rec): pair(d[(d.sweep == sw) & (d.reception == rec)], main[main.reception == rec], keys)
           for sw, rec, _ in runs + CASES[2:]}
    tab = []
    for s in rows_of(main):
        r = [LABEL[s]]
        for sw, rec, _ in runs:
            p = res[(sw, rec)]
            p = p[p.sid == s]
            r.append("n/a" if p.empty else " / ".join(fd((q.psr_mean_a - q.psr_mean_b).mean())
                                                   for q in [p[p.k == 8], p[p.k == 30]]))
        tab.append(r)
    out = ["Sensitivity runs against the main sweep: pairs of cells (same series, model, ci mode, R, N_dev, k; "
           "feasible in both), each series at its best configuration per cell and ci mode (BCC at 500 kchip/s, "
           "the chip rate of the sensitivity sweeps).\n\n(a) S1, S3, S4: mean over the pairs of PSR(run) - "
           "PSR(main), at k 8 / k 30. S4 "
           "'low' is the optimistic channel set of PARAMETERS 5 (low inter-body loss), which for coexistence means "
           "stronger interference; 'high' the reverse. n/a = the series is not in that run.\n\n",
           md(["series"] + [x[2] for x in runs], tab)]
    tab = []
    for sw, rec, lab in CASES[2:]:
        row = [lab, rec]
        for s in BCC:
            p = res[(sw, rec)]
            p = p[p.sid == s]
            row.append("n/a" if p.empty else
                       f"{fd((p.psr_mean_a - p.psr_mean_b).mean())} / {ratio(p, 'jpb_sensor_mean')} ({len(p)})")
        tab.append(row)
    out += ["\n(b) S5 BCC cases: mean over the pairs of PSR(case) - PSR(main) / median over the pairs of the sensor "
            "J/bit ratio case / main (pairs), every k. The contact case changes the chip rate as well (50 kchip/s), "
            "so its pairs hold both effects. n/a = the case does not run that stack or has "
            "no pair.\n\n", md(["case", "model"] + [LABEL[s] for s in BCC], tab)]
    write("t8_sensitivity.md", "".join(out))


def read_ble_degradation():
    out = {}
    for line in (HERE / "ble_degradation_table.tex").read_text().splitlines():
        m = re.match(r"\s*(\d+)--(\d+) & (.+?)\\\\", line)
        if m and "tab:ble" not in line:
            out[f"{m.group(1)}-{m.group(2)}"] = [float(x) for x in m.group(3).split("&")]
    assert list(out) == ["1-10", "10-30", "30-100", "100-300", "300-500", "500-1000", "1000-3000"], out
    return out


def t9(d):
    printed = read_ble_degradation()
    th = pd.read_csv(THESIS_CELLS)
    th = th[th.phy_rate == 1000000]
    d = d[d.cid.isin(["ble-1M", "ble4thesis-1M"]) & d.sweep.isin(MAIN_SWEEPS)].assign(sid=lambda x: x.cid)
    common_r = sorted(set(th.r_target_kbps) & set(d[d.sid == "ble-1M"].r_target_kbps))
    th_c = th[th.r_target_kbps.isin(common_r) & th.num_devices.isin([1, 4, 10])]
    b = th_c.groupby(["num_devices", "k", "r_target_kbps"]).mean_psr.max().reset_index()
    b["ra"] = b.num_devices * b.r_target_kbps
    p = b.groupby(["k", "ra"]).mean_psr.median().reset_index()
    p["band"] = band_of(p.ra)
    thesis_common = p.groupby(["band", "k"], observed=True).mean_psr.median()

    x = d[(d.reception == "collision") & d.r_target_kbps.isin(common_r)]
    L = bands(points(pick(x, "psr"), ["psr_mean"]), ["psr_mean"])
    s0 = d[(d.sid == "ble4thesis-1M")].merge(th, left_on=["ci_mode", "n_dev", "k", "r_target_kbps"],
                                             right_on=["ci_mode", "num_devices", "k", "r_target_kbps"])
    s0 = s0[s0.feasible == 1]
    diff = (s0.psr_mean - s0.mean_psr).abs()
    srcs = [("printed table", None), ("BLE_SIM cells, common grid", None),
            ("this engine, ble4thesis 1M", "ble4thesis-1M"), ("this engine, BLE 1M full", "ble-1M")]
    tab = []
    for band in printed:
        for name, sid in srcs:
            if sid is None and name.startswith("printed"):
                vals = [f3(v) for v in printed[band]]
            elif sid is None:
                vals = [f3(thesis_common.get((band, k))) for k in KS]
            else:
                y = L[(L.sid == sid) & (L.band == band) & (L.n_feas > 0)].set_index("k")
                vals = [f3(y.psr_mean.get(k)) for k in KS] if len(y) else ["n/f"] * len(KS)
            tab.append([band if name.startswith("printed") else "", name] + vals)
    trials = " or ".join(str(t) for t in sorted(th.n_trials.unique()))
    write("t9_ble_chapter6.md",
          "BLE 1M PSR, collision model, rule (higher PSR of the two ci modes, median over N_dev, median "
          f"over the band). Printed: the BLE table (BLE_SIM, {th.r_target_kbps.nunique()} R values, N_dev "
          f"{', '.join(str(n) for n in sorted(th.num_devices.unique()))}, {trials} trials per cell, "
          f"{', '.join(str(t) for t in sorted(th_c.n_trials.unique()))} on the common grid). The other rows use the "
          f"common grid of both sweeps: R {', '.join(f'{r:g}' for r in common_r)} "
          "kb/s and N_dev 1, 4, 10. n/f = no feasible cell. Cell by cell, ble4thesis reproduces BLE_SIM: "
          f"{len(s0)} matched cells, mean |dPSR| {diff.mean():.4f}, max {diff.max():.3f}.\n\n" +
          md(["R_agg (kb/s)", "source"] + [f"k = {k}" for k in KS], tab) +
          "\nReading the rows: printed vs BLE_SIM common grid is the cell mix of the band (the as-built sweep has "
          "many more R values and N_dev 2 and 6). BLE_SIM vs ble4thesis is the engine (the same rules). ble4thesis "
          "vs BLE 1M full is the protocol model: BLE full carries 10 octets of framing per PDU (80 us at 1M, so a "
          "10-octet payload takes 160 us instead of 80 us, which doubles the vulnerable window at small payloads; "
          "4 % at 250 octets), a hub poll before each data PDU, ARQ without a limit, AFH and CSA#2 instead of one "
          "permutation per BAN.\n")


# ladder levels and pools of T10 / F9 (pool '': 0.05, from L3; (0, '0.05'): the pooled reference), what each level
# changes (DESIGN 4.3), and the points
LADDER = [(0, ""), (0, "0.05"), (1, ""), (2, ""), (3, ""), (3, "0.43"), (4, ""), (4, "0.43"), (5, ""), (5, "0.43")]
STEP = ["BCC full: the IEEE 802.15.6 HBC stack", "RI mode (header at the payload rate), a pilot per 4096 chips",
        "the next beacon acknowledges, no I-Ack", "the retry pool as at the pooled reference (L2 with the pool)",
        "1-octet MAC header (the BAN ID in the beacon)", "one preamble, SFD and header per stripe (one tone: as L4)"]
STEP_POOLED = ("the HBC stack with its retries in a hub-polled pool (5 %) instead of a spare exchange per frame; the "
               "the lead reference since")
ROLE = ["802.15.6 HBC stack (multi-tone: on the as-built PHY)", "802.15.6 option", "proposal", "hub policy within 802.15.6",
        "proposal", "proposal"]
ROLE_POOLED = "hub policy within 802.15.6"
PT_E, PT_L, PT_C = (4, 10), [(4, 10), (4, 50)], (4, 1)      # (N_dev, R kb/s): J/bit, latency mode, contact case


def ladder(d):
    """Collision-model rows of BCC full one tone and multi-tone with their ladder level and pool ('' = 0.05 or none):
    L0 from the main sweep (contact case: the 'contact' rows of ladder_contact), L1 to L5 from the ladder sweeps."""
    x = d[d.cid.isin(LCID) & (d.reception == "collision") &
          (((d.src == "main") & (d.sid == "bcc-full")) | d.sweep.isin(["ladder", "ladder_contact"]))]
    x = x.assign(sid=x.cid)
    return x.assign(level=x.option.str.extract(r"opt(\d)", expand=False).fillna("0").astype(int),
                    pool=x.option.str.extract(r"pool=([\d.]+)", expand=False).fillna(""))


def at(x, pt, mode):
    """Median over k of the feasible cells at pt = (N_dev, R) in ci mode (NaN: none feasible); None: none ran."""
    c = x[(x.n_dev == pt[0]) & (x.r_target_kbps == pt[1]) & (x.ci_mode == mode)]
    return None if c.empty else c[c.feasible == 1].median(numeric_only=True)


def cell(v, text):
    """A point of at(): not run, n/f or text(v)."""
    return "not run" if v is None else "n/f" if pd.isna(v.period_us) else text(v)


def latency(v):
    return f"{v.period_us / 1e3:g} / {fg(v.lat_p95_ms_mean)} / {f3(v.on_time_mean)}"


def energy(v):
    return f"{fe(v.jpb_sensor_mean)} ({v.period_us / 1e3:g})"


def t10(d):
    """T10, the galvanic-optimisation ladder; returns one record per stack, level and pool for F9."""
    L = ladder(d)
    grid, con = L[L.sweep != "ladder_contact"], L[L.sweep == "ladder_contact"]

    def contact(v):
        return f"{f3(v.psr_mean)} [{f3(v.psr_lo)}, {f3(v.psr_hi)}] / {fg(v.get('retx_per_frame_mean'))}"

    recs, tab = [], []
    for sid in LCID:
        n_t = len(bcc_spec(sid, 500e3).tones or [0])
        cap0, capp = capacity(sid, 500e3, 1), capacity(sid, 500e3, 1, 0, 0.05)   # the spare hub, the pooled reference
        for lv, pool in LADDER:
            x = grid[(grid.sid == sid) & (grid.level == lv) & (grid.pool == pool)]
            c = con[(con.sid == sid) & (con.level == lv) & (con.pool == pool)]
            f = x[x.feasible == 1]
            cap = [capacity(sid, 500e3, n, lv, float(pool) if pool else None) for n in (1, 4, 10)]
            ref = (lv, pool) == (0, "0.05")
            r = dict(sid=sid, level=lv, pool=pool, cap=cap[0], e=at(x, PT_E, "energy"),
                     lat=[at(x, pt, "latency") for pt in PT_L], ce=at(c, PT_C, "energy"), cl=at(c, PT_C, "latency"))
            recs.append(r)
            tab.append([("one tone" if n_t == 1 else f"multi-tone ({n_t} tones)") if (lv, pool) == (0, "") else "",
                        f"L{lv}", pool or ("-" if lv < 3 else "0.05"),
                        STEP_POOLED if ref else STEP[lv] if not pool else "as above, pool 0.43",
                        ROLE_POOLED if ref else ROLE[lv],
                        " / ".join(f"{v:g}" for v in cap), f"{cap[0] / cap0:.2f}", f"{cap[0] / capp:.2f}",
                        f"{f.ra.max():g}" if len(f) else "n/f" if len(x) else "not run",
                        f"{f3(f.psr_mean.min())} ({f3(f.psr_mean.median())})" if len(f) else "",
                        cell(r["e"], energy)] +
                       [cell(v, latency) for v in r["lat"]] +
                       [cell(v, contact) if sid == "bcc-sf8" else "no case" for v in (r["ce"], r["cl"])])
    te = {r["e"].period_us for r in recs if r["e"] is not None and pd.notna(r["e"].period_us)}
    same = (f"the energy-mode period is {te.pop() / 1e3:g} ms at every level and stack there, so the J/bit steps are "
            "the stack's own" if len(te) == 1 else "the energy-mode period differs across the rows (in parentheses)")
    (n_e, r_e), (n_c, r_c) = PT_E, PT_C
    write("t10_ladder.md",
          "The galvanic-optimisation ladder of BCC full: each level switches off one more layer "
          "of the IEEE 802.15.6 HBC stack that the isolated galvanic link makes redundant, L1 with 802.15.6 options, L3 as "
          "the hub policy of the pooled reference (held from L0 in its own row), the others with 802.15.6 "
          "options, L2, L4 and L5 as proposals for a galvanic amendment (L3 builds on the L2 beacon acknowledgement). "
          "The retry pool (from L3) at two design points: share 0.05 of the superframe's frames (the headline: the "
          "isolated link as modelled) and 0.43 = 0.3 / 0.7 (the pessimistic case, sized for the ~30 % contact loss of "
          "the HackRF measurement; the author's real-time DBFSK audio stream over BCC, >= 44.1 kb/s with no audible "
          "errors, shows the HackRF packet rates are not the channel limit). The series BCC galvanic (proposal) is "
          "L5, pool 0.05, on one tone or multi-tone at any chip rate (T12). The row 'L0, pool 0.05' is the pooled "
          "reference stack: the HBC stack with one exchange per frame in the allocation and its retries in "
          "the hub-polled pool of 5 % of the superframe's frames, a hub policy the standard permits, and the reference "
          "of the capacity comparison; L0 without the pool is the spare-exchange hub (one exchange + one spare per frame, this model's "
          "provisioning rule of 2026-09-24, not an HBC feature). L2 with the pool equals L3, so against the pooled "
          "reference the ladder has four steps (L1, L2, L4, L5). L0 is the main sweep, the pooled reference and L1 to "
          "L5 the ladder sweep (same grid, 500 kchip/s), collision model, every ci mode as run.\n\n"
          "Columns: capacity per body at 500 kchip/s from protocols.select_period (as T6a, N_dev 1 / 4 / 10, 0.5 kb/s "
          "steps, the beacon charged per period; DESIGN 4.3 and quote the scan-free bound, up to about 1 kb/s "
          "higher: 689.9 kb/s at multi-tone L5, pool 0.05), and its "
          "ratios at N_dev 1 to L0 (the spare-exchange hub) and to the pooled reference; largest feasible grid R_agg "
          "(any ci mode); psr: the minimum over the feasible rows "
          f"(their median); J/bit per sensor in energy mode at N_dev {n_e} x R {r_e:g} kb/s (T in ms): {same}; latency "
          "mode (reference period 7.5 ms): T (ms) / lat_p95 (ms) / on_time at "
          f"{' and '.join(f'N_dev {n} x R {r:g} kb/s' for n, r in PT_L)}; contact case (one tone, 50 kchip/s, the "
          f"measured two-state contact chain, sweep ladder_contact): psr [95 % interval over trials] / retransmissions "
          f"per frame at N_dev {n_c} x R {r_c:g} kb/s, energy mode (one frame of about 250 octets per period) and "
          "latency mode (short frames); measured on hardware: psr 0.66-0.72. Each point is the median over k of its "
          "feasible cells (the isolated BCC bodies do not interact; the interval: the medians of the cells' limits). "
          "n/f: no feasible cell; not run: the sweep has no row there; no case: multi-tone has no contact chain.\n\n" +
          md(["stack", "level", "pool", "what the level changes", "802.15.6 option or proposal",
              "capacity per body N_dev 1 / 4 / 10 (kb/s)", "x L0", "x pooled reference",
              "largest feasible grid R_agg (kb/s)",
              "psr min (median)", f"J/bit sensor, energy, {n_e} x {r_e:g} (T ms)"] +
             [f"latency {n} x {r:g}: T / lat_p95 / on_time" for n, r in PT_L] +
             [f"contact energy {n_c} x {r_c:g}: psr [95 %] / retx", f"contact latency {n_c} x {r_c:g}: psr [95 %] / "
              "retx"], tab))
    return recs


def cfg_label(pid, phy, rc):
    if pid == "bcc":
        return {"sf8": "one tone", "multitone": "multi-tone", "bfsk": "as built"}[phy] + f", {rc / 1e6:g} Mchip/s"
    if pid == "ieee802156nb":
        return f"{phy} kb/s"
    return {"oqpsk": "CSMA" if pid == "ieee802154" else "O-QPSK", "gts": "GTS"}.get(
        phy, phy.replace("coded_s", "Coded S=").replace("sf", "SF "))


@lru_cache(maxsize=None)
def tones(rc):
    """Tones of the multi-tone plan at chip rate rc (1 from 2 Mchip/s; none above 3.38 Mchip/s counts as 1)."""
    try:
        return len(bcc_spec("bcc-multitone", rc).tones or [0])
    except ValueError:
        return 1


def t11(good):
    """Which configuration each stack's link adaptation takes, per R_agg band (main and S6 dense sweeps)."""
    x = good[good.sweep.isin(DOM_SWEEPS) & good.sid.isin(SID) & (good.feasible == 1)]
    one = (x.phy == "multitone") & (x.rc.map(lambda r: tones(r) if pd.notna(r) else 2) < 2)
    x = x.assign(band=band_of(x.ra), phy=x.phy.where(~one, "sf8"))     # a one-tone plan runs the one-tone stack
    parts = ["Link adaptation: the configuration (rate, PHY, access, tone plan, chip rate) each stack takes per cell "
             "(the highest PSR over its configurations and ci modes; PSRs within 0.005 tie and the lower sensor J/bit "
             "wins), main and S6 dense sweeps, every N_dev, k and ci mode. Per configuration: the cells it wins in "
             "each R_agg band (empty: none), then its cells / the stack's feasible cells. The BCC series and the "
             "proposal choose among one tone and multi-tone at every chip rate run (multi-tone: one tone from 2 "
             "Mchip/s, the same rows as one tone, and no tone plan above 3.38 Mchip/s): from 2 Mchip/s the multi-tone "
             "rows are counted as one tone (the proposal there is ladder level 4 on one tone). " + RATE_NOTE + ".\n"]
    for rec in ["collision", "sinr"]:
        y = x[x.reception == rec]
        bs = [b for b in BANDS if (y.band == b).any()]
        tab = []
        for sid in rows_of(y):
            z = y[y.sid == sid]
            for i, ((ph, rc), g) in enumerate(z.groupby(["phy", "rc"], dropna=False, sort=True)):
                tab.append([LABEL[sid] if i == 0 else "", cfg_label(g.protocol.iloc[0], ph, rc)] +
                           [str(int((g.band == b).sum())) if (g.band == b).any() else "" for b in bs] +
                           [f"{len(g)} / {len(z)}"])
        parts.append(f"\n{rec} model\n\n" + md(["series", "configuration"] + [f"{b}" for b in bs] +
                                                ["cells / feasible"], tab))
    write("t11_configurations.md", "".join(parts))


def t12(d):
    """T12, BCC against chip rate; returns one record per configuration and rate for F8."""
    rcs = sorted(d.loc[(d.sweep == "main") & d.rc.notna(), "rc"].unique())
    recs, tab = [], []
    for sid, cid, lv, pool, lab in BCFG:
        first = True
        for rc in rcs:
            try:
                s = bcc_spec(cid, rc, lv, pool)
            except ValueError:          # multi-tone: no tone plan above 3.38 Mchip/s
                continue
            if cid == "bcc-multitone" and len(s.tones or []) < 2:     # one tone: the one-tone configuration's rows
                continue
            x = bcc_rows(d, sid, cid, rc)
            f = x[x.feasible == 1]
            r = dict(sid=sid, cid=cid, lab=lab, rc=rc, cap=capacity(cid, rc, 1, lv, pool),
                     ra=f.ra.max() if len(f) else np.nan, psr=f.psr_mean.min() if len(f) else np.nan,
                     e=at(x, PT_E, "energy"), lat=at(x, PT_L[0], "latency"))
            recs.append(r)
            tab.append([lab if first else "", f"{rc / 1e6:g}", "headline" if rc <= RC_HEAD else "upper bound",
                        len(s.tones) if s.tones else 1, f"{r['cap']:g}",
                        f"{r['ra']:g}" if len(f) else "n/f" if len(x) else "not run",
                        f"{f3(r['psr'])} ({f3(f.psr_mean.median())})" if len(f) else "",
                        cell(r["e"], energy), cell(r["lat"], latency)])
            first = False
    (n_e, r_e), (n_l, r_l) = PT_E, PT_L[0]
    mt = max(r["rc"] for r in recs if r["cid"] == "bcc-multitone") / 1e6
    write("t12_bcc_chip_rates.md",
          "BCC against chip rate: BCC full (the spare-exchange hub) on one tone and multi-tone, the pooled reference "
          "(the HBC stack with its retries in a hub-polled pool of 5 % of the superframe's frames, a hub "
          "policy the standard permits; the lead reference; run at 500 kchip/s and 2.625 Mchip/s, its capacity "
          "at every rate) on one tone and multi-tone, the proposal (ladder level 5, pool 0.05) on "
          "one tone and multi-tone, and the as-built stack, at every chip rate run (500 kchip/s: the main sweep, the "
          "pooled reference and the proposal the ladder sweep; the others: the rates sweep: the as-built sweep, "
          "the IEEE 802.15.6 HBC "
          "Walsh-chip rates and upward). Headline limit: 2.625 Mchip/s, the Walsh-chip rate of the 656 kb/s HBC "
          "mode and the highest HBC rate whose occupied band f -/+ 2 R_c stays in the measured 0.05-20 MHz ;"
          " every verdict of the report (T1-T3b, T13, T14, F1-F6, F10) uses BCC up to it. Model limit: one "
          f"BFSK tone stays error-free on every measured path up to {protocols.BCC_RC_MAX_CPS / 1e6:g} Mchip/s "
          "(results/bcc_rate_stability.md; at 9.5 Mchip/s the lower tone near 0 Hz errs on every path, at 9.95 every "
          f"path falls below PSR 0.99), so faster rates are refused. The rows above 2.625 Mchip/s ('upper bound') are "
          f"a {UPPER}: from 5.25 Mchip/s up the occupied band f -/+ 2 R_c leaves the measured 0.05-20 MHz channel. The "
          f"model: {RATE_MODEL}; it leaves out dispersion beyond the bound, the channel above 20 MHz and "
          "receiver impairments. Multi-tone rows only where the plan has 2 or more tones "
          f"(among the rates run, to {mt:g} Mchip/s): above it the plan holds one tone (the one-tone rows) up to "
          "3.38 Mchip/s and none beyond (spacing 4 R_c). Energy basis: the BCC transceiver-IC "
          "power per radio state is held constant over the chip rate (PARAMETERS 4.3, static bias dominates; typical "
          "case from a 10 Mb/s IC), so J/bit falls as the chip rate rises. "
          "Collision model. Columns: capacity per body at N_dev 1 from protocols.select_period (either ci mode, 0.5 "
          "kb/s steps); largest feasible grid R_agg of the sweep; psr: minimum over the feasible rows (their median); "
          f"J/bit per sensor in energy mode at N_dev {n_e} x R {r_e:g} kb/s (T ms); latency mode at N_dev {n_l} x R "
          f"{r_l:g} kb/s: T (ms) / lat_p95 (ms) / on_time (deadline 7.5 ms); each point the median over k of its "
          "feasible cells. n/f: no feasible cell; not run: no row.\n\n" +
          md(["configuration", "chip rate (Mchip/s)", "basis", "tones", "capacity per body N_dev 1 (kb/s)",
              "largest feasible grid R_agg (kb/s)", "psr min (median)", f"J/bit sensor, energy, {n_e} x {r_e:g} (T ms)",
              f"latency {n_l} x {r_l:g}: T / lat_p95 / on_time"], tab))
    return recs


def qualified(x, ok):
    """The rows of x that meet ok, cheapest first: lower sensor J/bit, at equal J/bit a radio stack before BCC (a BCC
    win must be strict), then by series and configuration (deterministic)."""
    return x[ok].assign(_b=x.sid.isin(BCC)).sort_values(["jpb_sensor_mean", "_b", "sid", "cid", "rc", "ci_mode"],
                                                        kind="mergesort")


def dominance(d, sids=CONTEND):
    """Per cell (model, N_dev, R, k) of the main and S6 dense sweeps and criterion ('psr': psr >= 0.99; 'on_time':
    also on_time >= 0.9): a stack qualifies when one of its configurations and ci modes meets the criterion, and takes
    the cheapest of those (link adaptation to the requirement). The winner (lowest sensor J/bit, 'none'), the best
    radio and the best BCC series with their J/bit (the winner's chip rate win_rc), and margin = radio J/bit / BCC J/bit
    (NaN: no radio qualifies)."""
    x = d[d.sweep.isin(DOM_SWEEPS) & d.sid.isin(sids)]
    key = ["reception", "n_dev", "r_target_kbps", "k"]
    cells = x.loc[~x.sid.isin(BCC), key].drop_duplicates().sort_values(key)    # cells the radio stacks ran
    x = x.merge(cells)
    out = []
    for crit in ("psr", "on_time"):
        ok = (x.feasible == 1) & (x.psr_mean >= 0.99) & ((x.on_time_mean >= 0.9) if crit == "on_time" else True)
        q = qualified(x, ok)

        def best(y, name):
            return y.groupby(key).head(1)[key + ["sid", "rc", "jpb_sensor_mean"]].rename(
                columns={"sid": name, "rc": name + "_rc", "jpb_sensor_mean": name + "_jpb"})
        m = cells.merge(best(q, "win"), how="left").merge(best(q[~q._b], "radio"), how="left")
        out.append(m.merge(best(q[q._b], "bcc"), how="left").assign(crit=crit))
    dm = pd.concat(out, ignore_index=True)
    return dm.assign(win=dm.win.fillna("none"), ra=dm.n_dev * dm.r_target_kbps, margin=dm.radio_jpb / dm.bcc_jpb)


def rates(rc):
    """Counts of the BCC chip rates in rc, most frequent first ('9M 10, 0.5M 2')."""
    return ", ".join(f"{r / 1e6:g}M {n}" for r, n in rc.value_counts().items())


def t13(dm, dm5, dmu, ref, dense_rc, trials):
    """dm: BCC at its best chip rate up to RC_HEAD per cell (the headline); dm5: BCC at 500 kchip/s only; dmu: every
    chip rate to 9 Mchip/s (the model upper bound); the radio stacks as in dm; the contenders CONTEND. ref: the three
    views with the as-built BCC a contender too (its wins, for reference). dense_rc: the BCC chip rates of S6."""
    crits = [(rec, c) for rec in ["collision", "sinr"] for c in ["psr", "on_time"]]
    head = {"psr": "psr >= 0.99", "on_time": "psr >= 0.99, on_time >= 0.9"}
    views = (dm, dm5, dmu)

    def won(x, w):
        return int((x.win == w).sum())
    sel = [[x[(x.reception == rec) & (x.crit == c)] for rec, c in crits] for x in views]
    tab = [[LABEL.get(w, "no stack serves")] + [" / ".join(str(won(v, w)) for v in vs) for vs in zip(*sel)]
           for w in SID + ["none"] if any((v.win == w).any() for v in views)]
    tab.append(["cells"] + [str(len(a)) for a in sel[0]])
    sref = [[x[(x.reception == rec) & (x.crit == c)] for rec, c in crits] for x in ref]
    tab.append([f"reference: {LABEL['bcc-thesis']} as a contender"] +
               [" / ".join(str(won(v, "bcc-thesis")) for v in vs) for vs in zip(*sref)])
    nb = [int(v.win.isin(BCC).sum()) for v in views]
    hd = ", ".join(f"{r / 1e6:g}" for r in dense_rc if r <= RC_HEAD)
    lo = dm[(dm.reception == "sinr") & (dm.margin < 1)]
    hbc = (f" Of the SINR-model cells where both qualify and the margin is below 1 (the other stack spends less per "
           "bit), IEEE HBC (capacitive, not a 2.4 GHz radio) is the best other stack in "
           f"{(lo.radio == 'hbc').mean():.0%}"
           f" ({len(lo)} cells, both criteria)." if len(lo) else "")
    out = ["Dominance: per cell (model, N_dev, R per sensor, k) of the main sweep (N_dev 1, 4, 10) and the S6 dense "
           "sweep (N_dev 20, 30, 50) the stacks that qualify: one of their configurations and ci modes has psr >= 0.99 "
           "(second criterion also on_time >= 0.9), and each takes the cheapest of those (link adaptation to the "
           "need); among them the lowest sensor J per delivered bit wins (F10), a BCC series only when strictly "
           "lower than every other stack. The contenders are the full stacks: the radio stacks, IEEE HBC, BCC full "
           "(one tone or multi-tone) and the galvanic proposal; the as-built BCC (no FCS, no ACK) is a reference only, "
           "and the last row of (a) counts the cells it would win as a contender. Gateway limits: BLE above 20 "
           "sensors per hub is infeasible (the nRF52840 S140 SoftDevice central "
           "keeps at most 20 connections; a BLE mesh or a second hub is not modelled), 802.15.4 GTS above 7, "
           "SmartBAN above 16 (nodes per hub) and ANT above 15 (channels of the S212 hub) are infeasible in "
           "protocols; the 802.15.6 BAN size (64) lies above the grid. Three views of BCC, the radio stacks "
           f"unchanged: the headline, {RATE_NOTE}; BCC at 500 kchip/s only (the as-built chip rate); and the {UPPER} "
           "(every chip rate run, to 9 Mchip/s). 'Other stack': the radio stacks and IEEE HBC (capacitive), every "
           "contender that is not BCC. " + trials + " The S6 dense sweep ran BCC at "
           f"{', '.join(f'{r / 1e6:g}' for r in dense_rc)} Mchip/s (BCC full and the proposal; the headline takes "
           f"{hd} Mchip/s at N_dev 20 to 50). Over every model, criterion and cell BCC wins {nb[0]} of {len(dm)} "
           "in the headline, "
           f"{nb[1]} at 500 kchip/s and {nb[2]} at the upper bound: the wins above 500 kchip/s rest on the "
           "constant-IC-power basis." + hbc + "\n\n(a) Cells won per series: headline (to 2.625 Mchip/s) / 500 "
           "kchip/s only / "
           "upper bound (to 9 Mchip/s).\n\n",
           md(["winner"] + [f"{rec}, {head[c]}" for rec, c in crits], tab)]

    def wins(x):
        """BCC wins of x, their chip rates, R_agg, margin median [min] and the other stack most often second."""
        w = x[x.win.isin(BCC)]
        return [f"{len(w)} ({int(w.margin.isna().sum())})", rates(w.win_rc),
                f"{w.ra.min():g}-{w.ra.max():g}" if len(w) else "",
                f"{w.margin.median():.2f} [{w.margin.min():.2f}]" if w.margin.notna().any() else "",
                LABEL[w.radio.mode().iloc[0]] if w.radio.notna().any() else ""]

    def both(x):
        return f"{x.margin.median():.2f}" if x.margin.notna().any() else "-"
    tab, tab_u = [], []
    for rec, c in crits:
        for n in sorted(dm.n_dev.unique()):
            x, x5, xu = [y[(y.reception == rec) & (y.crit == c) & (y.n_dev == n)] for y in views]
            w5 = x5[x5.win.isin(BCC)]
            tab.append([rec, head[c], str(n), str(len(x)), str(int((x.win != "none").sum()))] + wins(x) +
                       [f"{len(w5)} ({int(w5.margin.isna().sum())})", f"{both(x)} / {both(x5)}"])
            tab_u.append([rec, head[c], str(n), str(int((xu.win != "none").sum()))] + wins(xu) + [both(xu)])
    out += ["\n(b) Per model, criterion and N_dev: cells, cells some stack serves, BCC wins in the headline view "
            "(BCC full or the proposal; in parentheses: cells only a BCC series serves), the chip "
            "rates of those wins (Mchip/s and cells), their R_agg, their J/bit margin = best other stack J/bit / BCC "
            "J/bit, "
            "median [minimum] over the wins with another stack qualifying, the other stack most often second there; "
            "then the BCC wins with BCC at 500 kchip/s only, and the median margin over every cell where both another "
            "stack and a BCC series qualify (winner or not; > 1: BCC spends less per bit), headline / 500 kchip/s.\n\n",
            md(["model", "criterion", "N_dev", "cells", "served", "BCC wins (BCC only)", "chip rates of the BCC wins",
                "R_agg of BCC wins (kb/s)", "margin median [min]", "other stack second", "BCC wins at 500 k (BCC only)",
                "margin, all cells with both: headline / 500 k"], tab),
            f"\n(c) Upper bound, BCC at every chip rate to 9 Mchip/s ({UPPER}); columns as (b): cells served, BCC "
            "wins (BCC only), their chip rates, R_agg, margin median [min], other stack second, and the median margin "
            "over "
            "every cell with both.\n\n",
            md(["model", "criterion", "N_dev", "served", "BCC wins (BCC only)", "chip rates of the BCC wins",
                "R_agg of BCC wins (kb/s)", "margin median [min]", "other stack second", "margin, all cells with both"],
               tab_u)]
    write("t13_dominance.md", "".join(out))


# realistic body networks: (name, N_dev, R per sensor kb/s, k range, deadline ms, sweep (None: by N_dev), source)
SCENARIOS = [
    ("motion capture suit", 17, 25.0, (1, 4), 16.7, None,
     "17 inertial trackers, 60 Hz full-body update (Xsens MVN Awinda, Schepers et al. 2018, Xsens white paper); "
     "52 octets per sample assumed (quaternion, acceleration, rate of turn, field; 32-bit values); D one sample"),
    ("EEG, 21 wireless electrodes", 21, 6.1, (1, 8), 100.0, None,
     "the 10-20 system, 21 electrodes (Jasper 1958; ACNS Guideline 2, Acharya et al. 2016; verify); 256 Hz "
     "(ACNS Guideline 1, Sinha et al. 2016), 24 bit; "
     "D assumed"),
    ("EEG cap, 32 channels on one node", 1, 384.0, (1, 4), 100.0, None,
     "32 channels x 500 Hz x 24 bit on one amplifier node; D assumed"),
    ("12-lead ECG, 10 electrode nodes", 10, 8.0, (1, 16), 100.0, None,
     "10 electrodes, 500 Hz (AHA/ACCF/HRS, Kligfield et al. 2007), 16 bit; D assumed"),
    ("audio stream, author's demonstration", 1, 44.1, (1, 8), 150.0, None,
     "the author's real-time DBFSK .wav stream over BCC (>= 44.1 kb/s, no audible errors); D: ITU-T G.114 one-way "
     "delay 150 ms"),
    ("audio, CD-quality mono", 1, 705.6, (1, 8), 150.0, None,
     "44.1 kHz x 16 bit (IEC 60908); D: ITU-T G.114 150 ms"),
    ("smart textile, 50 sensors", 50, 0.8, (1, 8), 100.0, None,
     "50 strain, pressure or temperature sensors at 50 Hz x 16 bit (assumed)"),
    ("hospital ward", 4, 1.0, (16, 30), 1000.0, None,
     "point B (ward vital-signs patch, CI 0.5-1 s), k 16 to 30 patients in radio range (assumed)"),
    ("crowded gym", 4, 10.0, (30, 30), 100.0, "s3_gym",
     "S3 venue gym, high density; 4 fitness sensors at 10 kb/s (points C, D; assumed); D assumed"),
]


@lru_cache(maxsize=None)
def hub_limit(pid, phy):
    """Sensors per hub of a configuration: GATEWAY, else protocols max_dev (inf: no cap)."""
    if pid in GATEWAY:
        return GATEWAY[pid]
    m = None if pid == "bcc" else protocols.get_spec(pid, "full", phy).max_dev
    return math.inf if m is None else m


def rank(m, ok, n, dl=None):
    """Per series over the matched cells (m: every configuration and ci-mode row, ok: the row meets the need): sid,
    the share of cells where one of its rows meets the need, the median over those cells of the cheapest such row's
    J/bit and lat_p95, why a series serves no cell, the most frequent chip rate of those rows (BCC) and whether most
    of them are latency-mode periods below dl / 2 (their J/bit is not that of a period near the deadline dl); sorted
    by share, then J/bit."""
    key = ["r_target_kbps", "n_dev", "k"]
    q = qualified(m, ok).groupby(["sid"] + key).head(1)
    rows = []
    for sid in rows_of(m):
        y, c = m[m.sid == sid], q[q.sid == sid]
        lim = max(hub_limit(p, ph) for p, ph in set(zip(y.protocol, y.phy)))
        why = ("" if len(c) else
               f"above its gateway limit {lim}" if n > lim else
               "no period fits" if (y.feasible == 0).all() else
               "psr < 0.99" if not (y.psr_mean >= 0.99).any() else "deadline")
        short = dl is not None and ((c.ci_mode == "latency") & (c.period_us < dl * 500)).mean() > 0.5
        rows.append((sid, len(c) / y.groupby(key).ngroups, c.jpb_sensor_mean.median(), c.lat_p95_ms_mean.median(),
                     why, c.rc.mode().iloc[0] if c.rc.notna().any() else None, short))
    return sorted(rows, key=lambda r: (-r[1], r[2] if r[1] > 0 else np.inf))


def name(x):
    """Series label of a rank() row, with its chip rate (BCC) and ^ where its J/bit is that of a short period."""
    return LABEL[x[0]] + (f" ({x[5] / 1e6:g} Mchip/s)" if x[5] else "") + (" ^" if x[6] else "")


def verdict(r, over=1.0):
    """Winner, its served share, runner-up, best other stack and why, from rank() (sorted by served share, then
    J/bit: the winner is the cheapest of the stacks with the highest share); over: grid load / need, outside
    [1, 1.5] the winner is indicative."""
    srv = [x for x in r if x[1] > 0]
    if not r:
        return "", "", "", "", "no grid cell in the range"
    if not srv:
        return "none", "", "", "", "no stack serves: " + "; ".join(f"{LABEL[x[0]]} {x[4]}" for x in r)
    w = srv[0]
    ru = srv[1] if len(srv) > 1 else None
    assert ru is None or ru[1] < w[1] or ru[2] >= w[2], (w, ru)       # never dearer than an equal server
    rad = next((x for x in srv if x[0] not in BCC), None)
    why = f"serves {w[1]:.0%} of the cells, J/bit {fe(w[2])}, lat_p95 {fg(w[3])} ms"
    if ru:
        why += f"; runner-up {ru[1]:.0%}, J/bit x {ru[2] / w[2]:.2f} of the winner's"
    lost = [f"{LABEL[x[0]]} {x[4]}" for x in r if x[1] == 0]
    why += "; none: " + ", ".join(lost) if lost else ""
    win = name(w)
    if not 1 <= over <= 1.5:
        win, why = f"{win} (indicative: grid x {over:.2f} of the need)", "indicative, at the grid cell: " + why
    return (win, f"{w[1]:.0%}", f"{name(ru)} ({fe(ru[2])}, {ru[1]:.0%})" if ru else "",
            f"{name(rad)} ({fe(rad[2])}, x {rad[2] / w[2]:.2f}, {rad[1]:.0%})" if rad else "no other stack serves", why)


def trials_note(d):
    """on_time granularity: the sample phase is fixed per trial, so at N_dev 1 on_time moves in steps of 1/trials."""
    n = {w: ", ".join(f"{t:g}" for t in sorted(d[d.sweep == w].n_trials.dropna().unique())) for w in DOM_SWEEPS}
    return ("on_time: the sensor sampling phase is fixed per trial, so at N_dev 1 on_time moves in steps of about "
            f"1 / trials ({n['main']} trials per cell in the main sweep, {n['s6_dense']} in S6); an on_time >= 0.9 "
            "verdict at N_dev 1 (T14 points D, F, G) lies within that granularity.")


def t14(d, wps):
    tab = []
    for name, n, r, (k_lo, k_hi), dl, sw, src in SCENARIOS:
        sw = sw or ("main" if n <= 10 else "s6_dense")
        shown = False
        for rec in ["collision", "sinr"]:
            x = d[(d.sweep == sw) & (d.reception == rec) & d.sid.isin(CONTEND)]
            if x.empty:
                continue
            ns = sorted(x.n_dev.unique())
            n_g = next((v for v in ns if v >= n), ns[-1])
            x = x[x.n_dev == n_g]
            rs = sorted(x.r_target_kbps.unique())
            r_g = next((v for v in rs if v >= r), rs[-1])
            m = x[(x.r_target_kbps == r_g) & x.k.between(k_lo, k_hi)]
            ok = (m.feasible == 1) & (m.psr_mean >= 0.99) & (m.lat_p95_ms_mean <= dl)
            ks = sorted(m.k.unique())
            over = n_g * r_g / (n * r)
            grid = f"{n_g} x {r_g:g} (x {over:.2f} of the need), k {', '.join(map(str, ks))}"
            first, shown = not shown, True
            tab.append([name if first else "", src if first else "",
                        f"{n} x {r:g}, k {k_lo}-{k_hi}, D {dl:g} ms" if first else "", grid, rec,
                        *verdict(rank(m, ok, n_g, dl), over)])
    for wp in wps:
        for rec in ["collision", "sinr"]:
            m, ns = wp_cells(d[(d.sweep == "main") & (d.reception == rec) & d.sid.isin(CONTEND)], wp)
            ok = (m.feasible == 1) & (m.psr_mean >= 0.99) & (m.on_time_mean >= 0.9)
            tab.append([f"{wp[0]} {wp[1]}" if rec == "collision" else "",
                        "working-point table" if rec == "collision" else "",
                        f"R_agg {wp[2]:g}-{wp[3]:g}, N_dev {wp[6]}-{wp[7]}, k {wp[4]}-{wp[5]}" if rec == "collision"
                        else "", f"N_dev {', '.join(map(str, ns))}, "
                        f"{len(m[['r_target_kbps', 'n_dev', 'k']].drop_duplicates())} cells", rec,
                        *verdict(rank(m, ok, max(ns)))])
    write("t14_scenarios.md",
          "Realistic body networks: per scenario the matched grid cells (the smallest grid N_dev and R per sensor at "
          "or above the need, from the main sweep up to N_dev 10 and the S6 dense sweep above; the gym from S3, SINR "
          "model only; the k grid values in the range), with the grid cell's load against the need (x: grid N_dev x R "
          "/ need N_dev x R). A cell below the need or more than 1.5 times above it names its winner as indicative: "
          "a larger network counts against the stacks with slot or connection limits, and "
          "grid points at the need are a sweep change. A stack serves a cell when one of its configurations and ci "
          "modes (link adaptation to the need) has psr >= 0.99 and lat_p95 <= D (the working points A-H: "
          "on_time >= 0.9, the deadline D = P of each cell), and takes the cheapest of those. J/bit is taken at the "
          "sweep's two ci-mode periods (energy: up to 200 ms; latency: the shortest period that fits), not at the "
          "cheapest period that meets D: where the energy period misses D a stack is charged its latency-mode J/bit; "
          "^ marks a series whose J/bit comes mostly from latency-mode periods below D / 2, so its J/bit ratio "
          "overstates its cost at a period near D. The contenders are the full stacks: the radio stacks, IEEE HBC, "
          "BCC full and the galvanic proposal; the as-built BCC (no FCS, no ACK) is a reference only (T3, T12). "
          "Winner: among the stacks that serve every matched cell (served 100 %; when none does, those with the "
          "highest served share) the lowest median sensor J per delivered bit over its served cells, so a winner is "
          "never dearer than a runner-up that serves as many cells; runner-up: the next by served share, then J/bit; "
          "served: the share of the matched cells a stack serves; best other stack: the "
          "first non-BCC stack in that order (a radio stack or IEEE HBC, capacitive), its J/bit and its ratio to "
          "the winner's; why: the winner's numbers, the "
          "runner-up's J/bit against the winner's and the stacks that serve no cell (gateway limit: BLE above 20 "
          "sensors per hub, a second hub or a mesh not modelled, 802.15.4 GTS above 7, SmartBAN above 16, ANT above "
          "15; no period fits: the stack cannot carry the load; psr or deadline missed). BCC names its chip rate: "
          f"{RATE_NOTE}. {trials_note(d)} Sources name the figures; 'assumed' and 'verify' mark a figure to "
          "confirm.\n\n" +
          md(["scenario", "source", "need (N_dev x R kb/s, k, D)", "grid cells", "model", "winner", "served",
              "runner-up (J/bit, served)", "best other stack (J/bit, x winner, served)", "why"], tab))


# ---------------------------------------------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------------------------------------------
def style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 7.5, "axes.titlesize": 7.5, "axes.labelsize": 7.5,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7, "axes.spines.top": False,
        "axes.spines.right": False, "axes.edgecolor": "#c3c2b7", "axes.linewidth": 0.6, "xtick.color": "#52514e",
        "ytick.color": "#52514e", "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4, "axes.labelcolor": INK, "text.color": INK, "grid.color": GRID,
        "grid.linewidth": 0.5, "lines.linewidth": 1.4, "pdf.fonttype": 42, "axes.titlelocation": "left",
        "legend.frameon": False, "figure.dpi": 100})


def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf", metadata={"CreationDate": None})
    fig.savefig(FIG / f"{name}.png", dpi=200)
    plt.close(fig)


def kaxis(ax):
    ax.set_xscale("log", base=2)
    ax.set_xticks(KS, [str(k) for k in KS])
    ax.minorticks_off()


def note(ax, text, y=0.4):
    ax.text(0.5, y, text, transform=ax.transAxes, ha="center", va="center", color=MUTED, fontsize=6.5)


def f1_f2(lp):
    xs = np.array(sorted(lp.ra.unique()))
    mids = np.sqrt(xs[1:] * xs[:-1])
    edges = np.concatenate([[xs[0] ** 2 / mids[0]], mids, [xs[-1] ** 2 / mids[-1]]])
    cmap = ListedColormap(CLASS_COL)
    pages = {"a": BCC + ["ble", "ieee802154"], "b": ["tsch", "ieee802156nb", "smartban", "ant", "hbc"]}
    for rec, name in [("collision", "f1{}_landscape_collision"), ("sinr", "f2{}_landscape_sinr")]:
        for page, sids in pages.items():
            fig, axs = plt.subplots(len(sids), 1, figsize=(6.3, 0.95 * len(sids) + 1.45), sharex=True)
            for ax, sid in zip(axs, sids):
                ax.set_xscale("log")
                ax.set_yticks(np.arange(6) + 0.5, [str(k) for k in KS])
                ax.set_ylim(0, 6)
                ax.tick_params(length=2)
                ax.set_ylabel("k")
                ax.set_title(LABEL[sid])
                x = lp[(lp.sid == sid) & (lp.reception == rec)]
                if x.empty:
                    note(ax, "no cell of this series in the sweep", 0.5)
                    continue
                z = np.full((6, len(xs)), np.nan)          # blank: no cell of the series at that point
                for _, r in x.iterrows():
                    z[KS.index(r.k), np.searchsorted(xs, r.ra)] = 3 if r.feasible == 0 else (
                        0 if r.psr_mean >= 0.99 else (1 if r.psr_mean >= 0.90 else 2))
                ax.pcolormesh(edges, np.arange(7), np.ma.masked_invalid(z), cmap=cmap, vmin=-0.5, vmax=3.5,
                              edgecolors="#fcfcfb", linewidth=0.3)
            axs[-1].set_xlabel("R_agg per body (kb/s)")
            axs[-1].set_xlim(edges[0], edges[-1])
            fig.legend(handles=[Patch(color=c, label=lab) for c, lab in zip(CLASS_COL, CLASS_LAB)],
                       loc="lower center", ncol=4, bbox_to_anchor=(0.5, 0.0))
            rest = ("BCC full, BCC galvanic (proposal), BCC as built, BLE, 802.15.4" if page == "a"
                    else "TSCH, 802.15.6 NB, SmartBAN, ANT, IEEE HBC")
            fig.suptitle(f"PSR class per R_agg grid point, {rec} model, page {'1' if page == 'a' else '2'} of 2: the "
                         "highest PSR over the stack's\nconfigurations and ci modes (link adaptation), median over "
                         f"N_dev (no band median, unlike T1); blank: no cell.\n{rest}.\n{RATE_FIG}.", x=0.02, ha="left",
                         fontsize=8)
            h = fig.get_figheight()
            fig.tight_layout(rect=(0, 0.25 / h, 1, 1), h_pad=0.4)
            save(fig, name.format(page))


def partial(x, lw=1.4):
    """Line width of a band series: half lw where not every grid point of the band is feasible (the * of T1)."""
    return lw / 2 if (x.n_feas < x.n_pts).any() else lw


BCC_LW, BCC_MS = 1.1, 2.4           # the BCC overlay of F3 and F6, a little thinner than the panel's own series
BCC_DX = {"bcc-thesis": 0.35, "bcc-full": 0.15, "bcc-opt": -0.35}    # F3, F6: marks at k x 2 ** dx, in draw order


def bcc_lines(ax, B, xcol, ycol, marks=True):
    """The BCC series as reference lines under the panel's own series (the proposal uppermost: it often shares BCC
    full's values; drawn in BCC_DX order).
    marks (F3, F6, x = k): marks at k x 2 ** BCC_DX, off the radio marks at k, width by partial(); else (F4) thin lines
    without marks. '' when at least one has a point, else a title suffix saying none has."""
    drawn = False
    for sid in BCC_DX:
        x = B[B.sid == sid].sort_values(xcol)
        if len(x):
            if marks:
                ax.plot(x[xcol] * 2 ** BCC_DX[sid], x[ycol], zorder=2.5, **sty(sid, partial(x, BCC_LW), BCC_MS))
            else:
                ax.plot(x[xcol], x[ycol], zorder=2.5, **{**sty(sid, 0.8), "marker": None})
            drawn = True
    return "" if drawn else "; BCC n/f"


def f3_psr_vs_k(land, contact):
    L = land[(land.reception == "collision") & (land.n_feas > 0)]
    bands_ = ["1-10", "30-100", "300-500"]
    fig, axs = plt.subplots(len(PANELS), 3, figsize=(6.3, 9.3), sharex=True, sharey=True)
    for i, (title, sids) in enumerate(PANELS):
        for j, b in enumerate(bands_):
            ax = axs[i, j]
            kaxis(ax)
            ax.grid(axis="y")
            B = L[L.band == b]
            if b == "1-10":      # 50 kchip/s: the rate of the measurement and the contact chain
                ax.axhspan(0.66, 0.72, fc=MEAS, ec=MUTED, lw=0.4, zorder=0)      # no hatch: it hid the as-built value
                for sid, v in contact.items():
                    ax.axhline(v, color=COLOUR[sid], lw=1.2, ls=":", zorder=4)
            no_bcc = bcc_lines(ax, B, "k", "psr_mean")
            drawn = False
            for sid in sids:
                x = B[B.sid == sid].sort_values("k")
                if x.empty:
                    continue
                ax.fill_between(x.k, x.psr_lo, x.psr_hi, color=COLOUR[sid], alpha=0.13, lw=0, zorder=2)
                drawn = True
                ax.plot(x.k, x.psr_mean, zorder=3, **sty(sid, partial(x)))
            if not drawn:
                note(ax, "no feasible cell")
            ax.set_title(f"{title}, {b} kb/s{no_bcc}", fontsize=7)
            ax.set_ylim(0, 1.04)
            if j == 0:
                ax.set_ylabel("PSR")
    for ax in axs[-1]:
        ax.set_xlabel("k (coexisting BANs)")
    hs = [handle(s, lw=BCC_LW, ms=BCC_MS) for s in BCC]
    hs += [Line2D([], [], color="#52514e", lw=0.7, label="half width: not every grid point of the band feasible"),
           Patch(color="#52514e", alpha=0.13, label="median of the cells' 95 % limits")]
    hs += [Line2D([], [], color=COLOUR[s], lw=1.2, ls=":", label=f"{LABEL[s]}, contact chain, 50 kchip/s "
                  f"({contact[s]:.3f}),\npessimistic case (HackRF-calibrated)") for s in contact]
    hs += [Patch(fc=MEAS, ec=MUTED, lw=0.4, label="BCC measured on hardware, 50 kchip/s OOK (0.66-0.72)")]
    fig.legend(handles=hs, loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("PSR against k in three R_agg bands, collision model (rule, link adaptation).\nEvery "
                 "panel: the three BCC series where feasible (BCC n/f: none is; their marks set apart along k);\n"
                 f"{RATE_FIG}.\nDotted, 1-10 kb/s panels: the S5 contact-chain case (50 kchip/s, the pessimistic case) "
                 "by the\nsame rule and band, then the median over k",
                 x=0.02, ha="left", fontsize=8)
    fig.tight_layout(rect=(0, 0.1, 1, 1), h_pad=0.4)
    save(fig, "f3_psr_vs_k")


def f4_jpb(good, good5, metric, name, what):
    """good: BCC at its best chip rate per cell; good5: at 500 kchip/s (the power band and BCC full's line in it)."""
    def pts(g, sweeps):
        sel = g[(g.reception == "collision") & (g.n_dev == 4) & g.sweep.isin(sweeps)]
        p = points(sel, [metric])
        return p[np.isfinite(p[metric])]
    P = pts(good, ["main"])
    P5 = pts(good5, ["main", "s5_bcc_power_low", "s5_bcc_power_high"])
    fig, axs = plt.subplots(2, 4, figsize=(6.3, 5.0), sharex=True, sharey=True)

    def line(sweep, sid, k, p=P):
        return p[(p.sweep == sweep) & (p.sid == sid) & (p.k == k)].sort_values("ra")

    pw = line("s5_bcc_power_low", "bcc-full", 1, P5).merge(line("s5_bcc_power_high", "bcc-full", 1, P5), on="ra")
    f5 = line("main", "bcc-full", 1, P5)
    B = P[(P.sweep == "main") & (P.k == 1)]
    axs.flat[-1].axis("off")         # 7 panels on a 2 x 4 grid
    for ax, (title, sids) in zip(axs.flat, PANELS):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.grid(axis="y", which="major")
        ax.fill_between(pw.ra, pw[metric + "_x"], pw[metric + "_y"], color="#c3c2b7", alpha=0.6, lw=0, zorder=0)
        ax.plot(f5.ra, f5[metric], color=INK, lw=0.8, ls=":", zorder=2.4)
        bcc_lines(ax, B, "ra", metric, marks=False)     # markers: the radio stacks at k 16
        for sid in sids:
            x1, x16 = line("main", sid, 1), line("main", sid, 16)
            ax.plot(x1.ra, x1[metric], zorder=3, **{**sty(sid, 1.1), "marker": None})
            ax.plot(x16.ra, x16[metric], zorder=4, **{**sty(sid), "ls": "none"})
            if x1.empty and x16.empty:
                note(ax, "no finite value")
        ax.set_title(title, fontsize=7)
    for ax in axs[-1]:
        ax.set_xlabel("R_agg per body (kb/s)")
    for ax in axs[:, 0]:
        ax.set_ylabel(f"J per delivered bit, {what}")
    hs = [Line2D([], [], color="#52514e", label="radio stack, k 1"),
          Line2D([], [], color="#52514e", ls="none", marker="o", ms=3.2, label="radio stack, k 16")]
    hs += [handle(s, f"{LABEL[s]}, best chip rate to 2.625 M, k 1", marker=None, lw=0.8) for s in BCC]
    hs += [Line2D([], [], color=INK, lw=0.8, ls=":", label=f"{LABEL['bcc-full']} at 500 kchip/s, k 1"),
           Patch(color="#c3c2b7", alpha=0.6, label=f"{LABEL['bcc-full']} at 500 kchip/s, transceiver power low - "
                 "high, k 1")]
    fig.legend(handles=hs, loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle(f"J per delivered bit ({what}) against R_agg, collision model, N_dev 4, each stack at its best "
                 "configuration\nand ci mode per cell; every panel: the three BCC series at k 1 (BCC does not depend "
                 f"on k),\n{RATE_FIG},\nand BCC full at 500 kchip/s", x=0.02, ha="left", fontsize=8)
    fig.tight_layout(rect=(0, 0.17, 1, 1))
    save(fig, name)


def f5_lifetime(wpres, wps):
    order = SID
    fig, axs = plt.subplots(4, 2, figsize=(6.3, 9.0), sharex=True, sharey=True)
    for ax, wp in zip(axs.flat, wps):
        out = wpres[("collision", wp[0])]
        ax.set_xscale("log")
        ax.grid(axis="x")
        if not (out.n_feas > 0).any():
            note(ax, "no series has a feasible cell", 0.5)
        for i, sid in enumerate(order):
            x = out[(out.sid == sid) & (out.n_feas > 0)]
            if x.empty:
                continue
            o = x.iloc[0]
            ax.plot(o.lifetime_days_mean, [i], **{**sty(sid, ms=4.5), "ls": "none"})
            if o.n_feas < o.n_cells:        # median over a subset of the point's cells
                ax.annotate(f"{o.n_feas}/{o.n_cells}", (o.lifetime_days_mean, i), xytext=(5, 0),
                            textcoords="offset points", fontsize=5.5, color="#52514e", va="center", ha="left")
        ax.set_yticks(range(len(order)), [LABEL[s] for s in order])
        ax.set_ylim(len(order) - 0.5, -0.5)
        ax.set_title(f"{wp[0]} {wp[1]}\nR_agg {wp[2]:g}-{wp[3]:g} kb/s, k {wp[4]}-{wp[5]}", fontsize=7)
        ax.tick_params(axis="y", length=0)
    for ax in axs[-1]:
        ax.set_xlabel("sensor lifetime (days)")
        ax.set_xlim(0.8, 4000)
    hs = [Line2D([], [], ls="none", label="n/m: median over n of the m matching cells (the rest infeasible)")]
    fig.legend(handles=hs, loc="lower center", bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("Sensor lifetime at the working points (T3 rule, collision model); missing dot = no feasible cell;"
                 f"\n{RATE_FIG}",
                 x=0.02, ha="left", fontsize=8)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    save(fig, "f5_lifetime_working_points")


def f6_latency(d, band):
    sel = pick(d[(d.sweep == "main") & (d.reception == "collision") & (d.ci_mode == "latency")], "latency")
    cols = ["lat_p95_ms_mean", "lat_p95_ms_lo", "lat_p95_ms_hi", "on_time_mean", "on_time_lo", "on_time_hi"]
    B = bands(points(sel, cols), cols)
    B = B[(B.band == band) & (B.n_feas > 0)]
    lat_min = 0.5
    B["lat_p95_ms_lo"] = B.lat_p95_ms_lo.where(B.lat_p95_ms_lo >= lat_min)   # a 0 limit has no place on a log axis
    fig, axs = plt.subplots(len(PANELS), 2, figsize=(6.3, 9.3), sharex=True)
    for i, (title, sids) in enumerate(PANELS):
        for j, m in enumerate(["lat_p95_ms", "on_time"]):
            ax = axs[i, j]
            kaxis(ax)
            ax.grid(axis="y")
            no_bcc = bcc_lines(ax, B, "k", m + "_mean")
            head = f"{title} (no feasible cell)"      # in the title: a note in the plot would cross the BCC lines
            for sid in sids:
                x = B[B.sid == sid].sort_values("k")
                if x.empty:
                    continue
                ax.fill_between(x.k, x[m + "_lo"], x[m + "_hi"], color=COLOUR[sid], alpha=0.13, lw=0, zorder=2)
                ax.plot(x.k, x[m + "_mean"], zorder=3, **sty(sid, partial(x)))
                head = title
            if m == "on_time":
                ax.set_ylim(0, 1.04)
            else:
                ax.set_yscale("log")
                ax.set_ylim(lat_min, 300)
                ax.set_yticks([1, 3, 10, 30, 100, 300], ["1", "3", "10", "30", "100", "300"])
                ax.yaxis.set_minor_formatter(NullFormatter())
                ax.axhline(7.5, color=DEADLINE, lw=0.9, ls=(0, (2, 1.5)), zorder=1.8)     # under every series
            ax.set_title(f"{head}: {'latency p95 (ms)' if j == 0 else 'on_time (deadline 7.5 ms)'}{no_bcc}",
                         fontsize=7)
    for ax in axs[-1]:
        ax.set_xlabel("k (coexisting BANs)")
    hs = [handle(s, lw=BCC_LW, ms=BCC_MS) for s in BCC]
    hs += [Line2D([], [], color="#52514e", lw=0.7, label="half width: not every grid point of the band feasible"),
           Patch(color="#52514e", alpha=0.13, label="median of the cells' 95 % limits (not below 0.5 ms; BCC: none)"),
           Line2D([], [], color=DEADLINE, lw=0.9, ls=(0, (2, 1.5)), label="deadline D = P = 7.5 ms")]
    fig.legend(handles=hs, loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle(f"Latency p95 and on_time against k, R_agg {band} kb/s, latency ci mode (P = 7.5 ms), collision "
                 "model,\nper cell the configuration with the highest PSR, then the lowest lat_p95; every panel: the "
                 f"three\nBCC series where feasible (BCC n/f: none is; marks set apart along k);\n{RATE_FIG}", x=0.02,
                 ha="left",
                 fontsize=8)
    fig.tight_layout(rect=(0, 0.085, 1, 1), h_pad=0.4)
    save(fig, "f6_latency_on_time")


def f7_ber(d):
    floor = 1e-7
    s1r = sorted(d[d.sweep == "s1_nocifloor"].r_target_kbps.unique())
    x = d[(d.reception == "sinr") & (d.feasible == 1) & d.r_target_kbps.isin(s1r) & (d.src != "rates")]
    fig, axs = plt.subplots(1, 2, figsize=(6.3, 3.9), sharex=True, sharey=True)
    for ax, sw, title in zip(axs, ["main", "s1_nocifloor"], ["with the C/I floor (main)", "without the floor (S1)"]):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.plot([floor, 1], [floor, 1], color=MUTED, lw=0.8, zorder=0)
        for sid, lab in [("ieee802156nb", "NB"), ("smartban", "SmartBAN"), ("ble", "BLE"), ("bcc-full", "BCC full"),
                         ("bcc-opt", "proposal"), ("bcc-thesis", "BCC as built")]:
            y = x[(x.sweep == sw) & (x.sid == sid)]
            ax.plot(y.raw_ber_mean.clip(lower=floor), y.post_ber_mean.clip(lower=floor), alpha=0.55,
                    **{**sty(sid), "ls": "none", "mew": 0.9 if sid == "bcc-thesis" else 0},
                    label=f"{lab}: {int((y.raw_ber_mean == 0).sum())} / "
                          f"{int((y.post_ber_mean < 0.8 * y.raw_ber_mean).sum())} / {len(y)}")
        ax.legend(loc="upper left", fontsize=6, title="cells: raw BER 0 /\npost < 0.8 raw / all",
                  title_fontsize=6, alignment="left", borderaxespad=0.2)
        ax.set_title(title)
        ax.set_xlabel("raw BER (before FEC)")
        ax.grid(which="major")
    axs[0].set_ylabel("post-FEC BER")
    fig.suptitle(f"Payload BER per cell, SINR model, R {', '.join(f'{r:g}' for r in s1r)} kb/s, every "
                 "configuration, k and ci mode:\nthe FEC stacks and BCC at 500 kchip/s (no FEC; the proposal from the "
                 "ladder sweep,\n"
                 f"main only); 0 plotted at {floor:g}; grey line: post = raw (FEC corrects "
                 "nothing)", x=0.02, ha="left", fontsize=8)
    fig.tight_layout()          # tight_layout makes room for the suptitle itself
    save(fig, "f7_ber_fec")


def cfg_sty(sid, cid):
    """F8 style of a BCC configuration: BCC full by its tone plan (ink, grey), the pooled reference mid-grey squares
    and the proposal wine stars (both: multi-tone dotted), the as-built stack dashed, hollow and thin under the others
    (it shares values with the proposal)."""
    st = {**sty(cid), "color": COLOUR[cid if sid == "bcc-full" else sid], "zorder": 3}
    if sid == "bcc-thesis":
        st.update(lw=0.9, zorder=2.5)
    if sid in ("bcc-opt", "bcc-pooled"):
        st.update(marker=MARKER[sid], ms=5 if sid == "bcc-opt" else 3.2, mfc=COLOUR[sid],
                  ls="-" if cid == "bcc-sf8" else ":")
    return st


F8_TICKS = {"0.05", "0.2", "0.5", "1.3125", "2.625", "5.25", "9"}      # labelled chip-rate ticks of F8


def f8_chip_rate(recs):
    R = pd.DataFrame(recs)
    for m, col in [("e", "jpb_sensor_mean"), ("lat", "lat_p95_ms_mean"), ("lat", "on_time_mean")]:
        R[col] = R[m].map(lambda v: np.nan if v is None else v[col])
    R["jpb"] = R.jpb_sensor_mean * 1e9
    mt = R[R.cid == "bcc-multitone"].rc.max() / 1e6
    (n_e, r_e), (n_l, r_l) = PT_E, PT_L[0]
    panels = [("cap", "(a) capacity per body, N_dev 1 (kb/s)"), ("psr", "(b) psr, minimum over the feasible rows"),
              ("jpb", f"(c) sensor nJ per delivered bit, energy, {n_e} x {r_e:g} kb/s"),
              ("lat_p95_ms_mean", f"(d) latency p95 (ms), latency mode, {n_l} x {r_l:g} kb/s"),
              ("on_time_mean", f"(e) on_time (deadline 7.5 ms), latency, {n_l} x {r_l:g} kb/s")]
    rcs = sorted(R.rc.unique())
    fig, axs = plt.subplots(3, 2, figsize=(6.3, 7.6))
    groups = list(R.groupby(["sid", "cid"], sort=False))
    for ax, (key, title) in zip(axs.flat, panels):
        ax.set_xscale("log")
        ax.set_xticks(rcs, [f"{r / 1e6:g}" if f"{r / 1e6:g}" in F8_TICKS else "" for r in rcs], fontsize=6)
        ax.minorticks_off()
        ax.grid(axis="y", which="major")
        ax.set_title(title, fontsize=7)
        ax.axvspan(RC_HEAD * 1.15, rcs[-1] * 1.3, color="#f1f0ec", lw=0, zorder=0)
        ax.axvline(RC_HEAD, color=INK, lw=0.9, ls="--", zorder=1)
        ax.axvline(protocols.BCC_RC_MAX_CPS, color=DEADLINE, lw=0.9, ls=(0, (2, 1.5)), zorder=1)
        for (sid, cid), x in groups:
            ax.plot(x.rc, x[key], **cfg_sty(sid, cid))
            if key == "cap":
                ax.plot(x.rc, x.ra, ls="none", marker="o", ms=6, mfc="none", mec=cfg_sty(sid, cid)["color"],
                        mew=0.8, zorder=4)
    for a in (axs[0, 0], axs[1, 0], axs[1, 1]):
        a.set_yscale("log")
        a.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    axs[0, 1].set_ylim(min(0.9, R.psr.min() - 0.01), 1.005)
    axs[1, 1].axhline(7.5, color=MUTED, lw=0.7, ls=":")
    axs[2, 0].set_ylim(0, 1.04)
    for a in (axs[2, 0], axs[1, 1]):
        a.set_xlabel("chip rate per tone (Mchip/s; ticks at every rate run)")
    axs[2, 1].axis("off")
    labs = {(c[0], c[1]): c[-1] for c in BCFG}
    hs = [Line2D([], [], label=labs[k], **cfg_sty(*k)) for k, _ in groups]
    hs += [Line2D([], [], color=INK, lw=0.9, ls="--", label=f"headline limit {RC_HEAD / 1e6:g} Mchip/s (highest HBC "
                  "Walsh-chip rate whose\nband stays in the measured channel; every verdict uses BCC up to it)"),
           Patch(color="#f1f0ec", label="model upper bound: no phase data above what\n bounds; channel above "
                 "20 MHz unmeasured"),
           Line2D([], [], color=DEADLINE, lw=0.9, ls=(0, (2, 1.5)),
                  label=f"model limit {protocols.BCC_RC_MAX_CPS / 1e6:g} Mchip/s (one tone error-free\non every "
                        "measured path, bcc_rate_stability.md)"),
           Line2D([], [], ls="none", marker="o", ms=6, mfc="none", mec="#52514e", mew=0.8,
                  label="(a) simulated: largest feasible grid R_agg"),
           Line2D([], [], color=MUTED, lw=0.7, ls=":", label="(d) 7.5 ms")]
    axs[2, 1].legend(handles=hs, loc="center", fontsize=6.3)
    fig.suptitle("BCC against chip rate, collision model: capacity (protocols.select_period) and simulated metrics\n"
                 f"(median over k of the feasible cells at the point); multi-tone drawn where the plan has 2 or more\n"
                 f"tones (to {mt:g} Mchip/s, one tone above). IC power held constant over the chip rate (PARAMETERS\n"
                 "4.3), so J/bit falls as the rate rises. SNR per chip falls as 1/R_c at a fixed 1 V drive (no power\n"
                 "increase), noise band 4 R_c, tone gains at f -/+ R_c; the minimum-phase smear check shows\n"
                 "no loss up to 2.625 Mchip/s (<= 0.03 dB); above it the band leaves the measured 0.05-20 MHz",
                 x=0.02, ha="left", fontsize=8)
    fig.tight_layout(h_pad=0.6)
    save(fig, "f8_bcc_chip_rate")


def f10_dominance(dm):
    order = SID + ["none"]
    cmap = ListedColormap([COLOUR[s] for s in SID] + ["white"])      # no stack serves: white, hatched
    dark = {"ble", "hbc", "ant", "bcc-full", "bcc-opt", "tsch"}          # white text on these tiles
    for crit, name, what in [("psr", "f10a_dominance_psr", "psr >= 0.99"),
                             ("on_time", "f10b_dominance_on_time", "psr >= 0.99 and on_time >= 0.9")]:
        D = dm[dm.crit == crit]
        ns = sorted(D.n_dev.unique())
        hr = [D[D.n_dev == n].k.nunique() + 2 for n in ns]         # k rows plus the title
        fig, axs = plt.subplots(len(ns), 2, figsize=(6.3, 0.2 * sum(hr) + 1.9), squeeze=False,
                                gridspec_kw={"height_ratios": hr})
        for i, n in enumerate(ns):
            for j, rec in enumerate(["collision", "sinr"]):
                ax = axs[i, j]
                x = D[(D.n_dev == n) & (D.reception == rec)]
                rc0 = x.win_rc[x.win.isin(BCC)].mode()      # the panel's usual BCC chip rate, named in the title
                rc0 = rc0.iloc[0] if len(rc0) else None
                ax.set_title(f"{rec} model, N_dev {n}" + (f"; BCC at {rc0 / 1e6:g} Mchip/s" if rc0 else ""),
                             fontsize=7)
                if x.empty:
                    ax.axis("off")
                    continue
                rs, ks = sorted(x.r_target_kbps.unique()), sorted(x.k.unique())
                z = np.full((len(ks), len(rs)), np.nan)
                for r in x.itertuples():
                    iy, ix = ks.index(r.k), rs.index(r.r_target_kbps)
                    z[iy, ix] = order.index(r.win)
                    if r.win == "none":
                        ax.add_patch(Rectangle((ix - 0.5, iy - 0.5), 1, 1, fill=False, hatch="////", ec=MUTED, lw=0))
                    if r.win in BCC:
                        rc = f"{r.win_rc / 1e6:g}" if len(f"{r.win_rc / 1e6:g}") < 6 else f"{r.win_rc / 1e6:.2g}"
                        ax.text(ix, iy, ("only" if pd.isna(r.margin) else f"{r.margin:.2g}") +
                                ("" if r.win_rc == rc0 else f"\n{rc}M"), ha="center", va="center", fontsize=5.5,
                                linespacing=0.9, color="white" if r.win in dark else INK)
                ax.imshow(z, cmap=cmap, vmin=-0.5, vmax=len(order) - 0.5, aspect="auto", origin="lower",
                          interpolation="nearest")
                ax.set_xticks(range(len(rs)), [f"{v:g}" for v in rs], fontsize=5.5)
                ax.set_yticks(range(len(ks)), [str(v) for v in ks], fontsize=5.5)
                ax.tick_params(length=0)
                for sp in ax.spines.values():
                    sp.set_visible(False)
                if j == 0:
                    ax.set_ylabel("k")
        for ax in axs[-1]:
            ax.set_xlabel("R per sensor (kb/s)")
        have = [w for w in order if (D.win == w).any()]
        fig.legend(handles=[Patch(fc=cmap(order.index(w)), ec=MUTED if w == "none" else "none",
                                  hatch="////" if w == "none" else None, label=LABEL.get(w, "no stack serves"))
                            for w in have],
                   loc="lower center", ncol=4, bbox_to_anchor=(0.5, 0.0), fontsize=6.3)
        fig.suptitle(f"Dominance: per cell the stack with the lowest sensor J/bit among those with {what},\neach at "
                     "its cheapest qualifying configuration (link adaptation); BLE above 20 sensors, SmartBAN above "
                     "16,\nANT above 15 infeasible (gateway limits; no BLE mesh or second hub). On a BCC tile: best "
                     "other stack J/bit / BCC\nJ/bit (radio or IEEE HBC; only: none qualifies), and the chip rate "
                     "(Mchip/s) where it is "
                     "not the one in the panel title.\nBCC at its best chip rate per cell up to 2.625 Mchip/s "
                     "(headline limit),\nIC power constant over the rate, SNR per chip falls as 1/R_c at 1 V;\n"
                     "T13: also at 500 kchip/s only and the model "
                     "upper bound to 9 Mchip/s; BCC as built: a reference, not a contender", x=0.02, ha="left",
                     fontsize=7.5)
        fig.tight_layout(rect=(0, 0.3 / fig.get_figheight(), 1, 1), h_pad=0.5)
        save(fig, name)


def f9_ladder(recs):
    R = pd.DataFrame(recs)
    R["lat0"] = R.lat.str[0]

    def val(key, col):
        return R[key].map(lambda v: np.nan if v is None else v[col])
    R["jpb"], R["t_ms"] = val("e", "jpb_sensor_mean") * 1e9, val("lat0", "period_us") / 1e3
    R["lat95"], R["on"] = val("lat0", "lat_p95_ms_mean"), val("lat0", "on_time_mean")
    for m in ("ce", "cl"):
        for col in ("psr_mean", "psr_lo", "psr_hi"):
            R[f"{m}_{col}"] = val(m, col)
    (n_e, r_e), (n_l, r_l), (n_c, r_c) = PT_E, PT_L[0], PT_C
    panels = [("cap", "(a) capacity per body at 500 kchip/s, N_dev 1 (kb/s)"),
              ("jpb", f"(b) sensor nJ per delivered bit, energy mode, {n_e} x {r_e:g} kb/s"),
              ("t_ms", f"(c) latency-mode period T (ms), {n_l} x {r_l:g} kb/s"),
              ("lat95", f"(d) latency p95 (ms), latency mode, {n_l} x {r_l:g} kb/s"),
              ("on", f"(e) on_time (deadline 7.5 ms), latency mode, {n_l} x {r_l:g} kb/s"),
              ("contact", f"(f) contact case psr (one tone, 50 kchip/s), {n_c} x {r_c:g} kb/s")]
    fig, axs = plt.subplots(3, 2, figsize=(6.3, 7.7), sharex=True)
    for ax, (key, title) in zip(axs.flat, panels):
        for lv in (2, 4, 5):
            ax.axvspan(lv - 0.5, lv + 0.5, color="#f1f0ec", lw=0, zorder=0)
        ax.set_title(title, fontsize=7)
        ax.grid(axis="y", which="major")
        for sid in LCID if key != "contact" else ["bcc-sf8"]:
            for pool in ("", "0.43"):      # pool 0.43 branches off at L2, the last level without a pool
                x = R[(R.sid == sid) & ((R.pool == pool) | ((R.level == 2) & (pool != "")))].sort_values("level")
                # the branch under pool 0.05, its hollow marks from L3: where both pools agree the filled marks show
                st = {**sty(sid), **({"ls": "--", "mfc": "white", "markevery": slice(1, None)} if pool else {}),
                      "zorder": 2.9 if pool else 3}
                if key != "contact":
                    ax.plot(x.level, x[key], **st)
                    continue
                for m, dx, mk in (("ce", -0.1, "D"), ("cl", 0.1, "o")):
                    y = x[f"{m}_psr_mean"]
                    ax.errorbar(x.level + dx, y, yerr=[y - x[f"{m}_psr_lo"], x[f"{m}_psr_hi"] - y],
                                elinewidth=0.7, capsize=1.5,
                                **{**st, "marker": mk, "ls": ":" if m == "cl" else st["ls"]})
            # the pooled reference: a square at L0, in the stack's colour; (f) hollow for the latency mode
            ref = R[(R.sid == sid) & (R.level == 0) & (R.pool == "0.05")]
            st = {**sty(sid, ms=3.6), "marker": "s", "ls": "none", "zorder": 3.5}
            if key != "contact":
                ax.plot(ref.level, ref[key], **st)
                continue
            r0 = ref.iloc[0]        # one record per stack; scalars: errorbar warns on single-element Series
            for m, dx in (("ce", -0.1), ("cl", 0.1)):
                y, lo, hi = (r0[f"{m}_psr_{c}"] for c in ("mean", "lo", "hi"))
                ax.errorbar([r0.level + dx], [y], yerr=[[y - lo], [hi - y]], elinewidth=0.7, capsize=1.5,
                            **{**st, "mfc": "white" if m == "cl" else st["mfc"]})
        if key != "contact":
            opt = R[(R.sid == "bcc-multitone") & (R.level == 5) & (R.pool == "")]
            ax.plot(opt.level, opt[key], zorder=4, **{**sty("bcc-opt", ms=5), "ls": "none"})
    for a in (axs[0, 0], axs[1, 0], axs[1, 1]):
        a.set_yscale("log")
        a.yaxis.set_major_locator(LogLocator(subs=(1, 2, 5)))
        a.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        a.yaxis.set_minor_formatter(NullFormatter())
    axs[0, 1].set_ylim(bottom=0)       # (b) spans less than a decade: a log axis left the one-tone values unlabelled
    te = val("e", "period_us").dropna().unique()        # T10's check: the J/bit steps are the stack's own
    p = R[R.pool == "0.43"].merge(R[R.pool == ""], on=["sid", "level"])
    axs[0, 1].text(0.97, 0.95, (f"T {te[0] / 1e3:g} ms at every level" if len(te) == 1 else "T differs by level (T10)")
                   + ("\npool 0.43: as pool 0.05" if np.allclose(p.jpb_x, p.jpb_y) else ""),
                   transform=axs[0, 1].transAxes, ha="right", va="top", fontsize=6.5, color="#52514e")
    a = axs[0, 0]
    for _, r in R[R.level.isin([0, 5]) & (R.pool == "")].iterrows():
        a.annotate(f"{r.cap:g}", (r.level, r.cap), xytext=(0, 5), textcoords="offset points", ha="center",
                   fontsize=6, color="#52514e")
    for _, r in R[((R.level == 5) & (R.pool != "")) | ((R.level == 0) & (R.pool == "0.05"))].iterrows():
        a.annotate(f"{r.cap:g}", (r.level, r.cap), xytext=(0, -9), textcoords="offset points", ha="center",
                   fontsize=6, color="#52514e")
    axs[1, 1].axhline(7.5, color=DEADLINE, lw=0.9, ls=(0, (2, 1.5)), zorder=1.8)
    axs[2, 0].set_ylim(0, 1.04)
    axs[2, 1].axhspan(0.66, 0.72, fc=MEAS, ec=MUTED, lw=0.4, zorder=0.5)
    axs[2, 1].set_ylim(0, 1.04)
    for a in axs[-1]:
        a.set_xticks(range(6), [f"L{i}\n" + ("HBC" if i == 0 else ROLE[i].split()[-1]) for i in range(6)], fontsize=6.5)
        a.set_xlim(-0.5, 5.5)
    hs = [handle("bcc-sf8", "BCC full one tone and its ladder, pool 0.05"),
          handle("bcc-sf8", "one tone, pool 0.43 (from L3), pessimistic case", ls="--", mfc="white"),
          handle("bcc-multitone", "BCC full multi-tone (6 tones) and its ladder, pool 0.05"),
          handle("bcc-multitone", "multi-tone, pool 0.43 (from L3), pessimistic case", ls="--", mfc="white"),
          handle("bcc-opt", "BCC galvanic (proposal), multi-tone: L5, pool 0.05", ls="none", ms=8),
          Line2D([], [], color="#52514e", marker="s", ms=3.6, ls="none", label="pooled reference at L0 (the HBC stack "
                 "with the 5 % retry pool),\nthe lead reference; colour by stack"),
          Patch(color="#f1f0ec", label="proposal level (galvanic amendment); L1, L3: 802.15.6 options"),
          Line2D([], [], color=INK, marker="D", ms=3.2, ls="-",
                 label="(f) energy mode, one ~250-octet frame per period"),
          Line2D([], [], color=INK, marker="o", ms=3.2, ls=":", label="(f) latency mode, short frames; bars: 95 % over "
                 "trials"),
          Patch(fc=MEAS, ec=MUTED, lw=0.4, label="(f) BCC measured on hardware, 50 kchip/s (0.66-0.72)"),
          Line2D([], [], color=DEADLINE, lw=0.9, ls=(0, (2, 1.5)), label="(d) deadline D = P = 7.5 ms")]
    fig.legend(handles=hs, loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.0), fontsize=6.3)
    fig.suptitle("The galvanic-optimisation ladder of BCC full (DESIGN 4.3): one layer of the 802.15.6 HBC stack "
                 "off per level,\nL1 RI mode and pilot period, L2 beacon acknowledgement, L3 retry pool (an 802.15.6 "
                 "option, but a sensor\nlearns its failures from the L2 acknowledgement; poll and pool details "
                 "assumed), L4 1-octet MAC header,\nL5 one head per stripe. Collision model; (b) to (f): median over k "
                 "of the feasible cells at the point (T10).\nPessimistic case, calibrated on the HackRF measurement: "
                 "pool 0.43 and the contact case (f)", x=0.02, ha="left", fontsize=8)
    fig.tight_layout(rect=(0, 0.11, 1, 1), h_pad=0.6)
    save(fig, "f9_ladder")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--raw", default=RES / "raw", type=Path, help="directory of the sweeps' *_cells.csv "
                    "(default results/raw); without any, the committed compact CSVs are used")
    a = ap.parse_args()
    t0 = time.time()
    self_check()
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    built = build_compact(a.raw)
    d = load()
    dh = d[~d.sid.isin(BCC) | (d.rc <= RC_HEAD)]        # the headline: BCC up to 2.625 Mchip/s
    best, good = pick(dh, "psr"), pick(dh, "better")
    good5 = pick(d[d.src != "rates"], "better")        # BCC at 500 kchip/s: the chip rate of the S1 to S5 sweeps
    lp = points(best[best.sweep.isin(MAIN_SWEEPS) & best.sid.isin(SID)], ["psr_mean", "psr_lo", "psr_hi"])
    land = bands(lp, ["psr_mean", "psr_lo", "psr_hi"])
    land = land.assign(cls=np.select([land.n_feas == 0, land.psr_mean >= 0.99, land.psr_mean >= 0.90],
                                     ["infeasible", ">=0.99", "0.90-0.99"], "<0.90"))
    land.drop(columns=["sweep", "option"]).to_csv(RES / "landscape.csv", index=False, float_format="%.4g",
                                                  lineterminator="\n")
    wps = read_working_points()
    t1_t2(land)
    wpres = t3(good, wps)
    t4(good5)
    t5(good)
    t6(d, good5, wps)
    t7(d)
    t8(d)
    if THESIS_CELLS.exists():    # the BLE landscape of the BLE model, not in this repository
        t9(d)
    recs = t10(d) if (d.sweep == "ladder").any() else None
    t11(good)
    rrecs = t12(d)
    dm = dominance(dh)
    d5 = d[~d.sid.isin(BCC) | (d.rc == 500e3)]
    t13(dm, dominance(d5), dominance(d), [dominance(x, SID) for x in (dh, d5, d)],
        sorted(d.loc[(d.src == "s6_dense") & d.rc.notna(), "rc"].unique()), trials_note(d))
    t14(dh, wps)
    style()
    f1_f2(lp)
    c = best[(best.sweep == "s5_bcc_contact") & (best.reception == "collision")]
    c = bands(points(c, ["psr_mean"]), ["psr_mean"])
    c = c[c.band == "1-10"].groupby("sid").psr_mean.median()
    f3_psr_vs_k(land, {s: c[s] for s in BCC if s in c})
    f4_jpb(good, good5, "jpb_sensor_mean", "f4a_jpb_sensor", "sensor")
    f4_jpb(good, good5, "jpb_ban_mean", "f4b_jpb_ban", "BAN")
    f5_lifetime(wpres, wps)
    f6_latency(dh, "30-100")
    f7_ber(d)
    f8_chip_rate(rrecs)
    if recs:
        f9_ladder(recs)
    else:
        print("report: no ladder sweep in the compact CSVs, T10 and F9 skipped")
    f10_dominance(dm)
    print(f"report: raw {a.raw} {'read' if built else 'absent, committed CSVs used'}; {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()

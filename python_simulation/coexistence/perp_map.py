"""PSR region of one stack on a triangle with one axis per side, read by perpendiculars.

Each side carries one parameter on a log2 scale: 1 at the side's midpoint, 128 at one corner. A cell (sensors,
datarate, networks) sits where the perpendiculars from its three axis values meet. The three feet, measured from the
midpoints in one rotational sense, sum to zero, so the perpendiculars meet only for cells with
log(sensors) - log(datarate) +/- log(networks) = 0; the sign depends on which corner the networks axis grows to, and
both choices are drawn. Only simulated cells on that surface are plotted (results/triangle/psr_map_fixed.csv); the
region is the convex hull of the cells with PSR >= --psr.

usage: python perp_map.py [--stack BLE] [--psr 0.99] [--out results/triangle/perp_map_ble.png]
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull, QhullError

HALF_LOG2 = 7                                   # midpoint 1 .. corner 128
R = 1 / np.sqrt(3)                              # side 1, centroid at the origin
T, L, Rt = (np.array([R * np.cos(a), R * np.sin(a)]) for a in np.radians([90, 210, 330]))


def unit(a, b):
    return (b - a) / np.linalg.norm(b - a)


def axes(net_to):
    """(name, side start, side end, direction of growth) per parameter."""
    return [("sensors per network", T, L, unit(T, L)),
            ("datarate per sensor, x 10 kbit/s", Rt, L, unit(Rt, L)),
            ("number of networks", T, Rt, unit(T, Rt) if net_to == "R" else unit(Rt, T))]


def place(cells, spec):
    """Intersection point of the perpendiculars per cell, NaN where they do not meet in one point."""
    lv = np.log2(cells) / (2 * HALF_LOG2)
    P = lv[:, :2] @ np.linalg.inv(np.array([spec[0][3], spec[1][3]])).T
    P[np.abs(P @ spec[2][3] - lv[:, 2]) > 1e-9] = np.nan
    return P


def panel(ax, b, net_to, psr, title):
    spec = axes(net_to)
    cells = np.c_[b["n_sensors"], b["datarate"] / 10, b["n_networks"]].astype(float)
    P = place(cells, spec)
    on = ~np.isnan(P[:, 0]) & (cells >= 1).all(1) & (cells <= 2 ** HALF_LOG2).all(1)
    good = on & b["feasible"].to_numpy() & (b["psr"].fillna(0).to_numpy() >= psr)
    ax.plot(*np.array([T, L, Rt, T]).T, color="0.3", lw=1.2)
    for name, a, e, d in spec:                  # ticks from the midpoint to the corner the axis grows to
        m = (a + e) / 2
        nrm = np.array([d[1], -d[0]])
        if nrm @ m < 0:
            nrm = -nrm
        for k in range(HALF_LOG2 + 1):
            p = m + d * k / (2 * HALF_LOG2)
            ax.plot(*np.array([p, p + 0.015 * nrm]).T, color="0.3", lw=1)
            ax.text(*(p + 0.04 * nrm), f"{2 ** k}", ha="center", va="center", fontsize=8)
            g, = ax.plot(*np.array([p, p - 0.9 * nrm]).T, color="0.75", lw=0.4, ls=":", zorder=0)
            g.set_clip_path(plt.Polygon([T, L, Rt], transform=ax.transData))
        ang = np.degrees(np.arctan2(d[1], d[0]))
        ang = ang - 180 if abs(ang) > 90 else ang
        ax.text(*(m + 0.11 * nrm + d * 0.25), name, rotation=ang, ha="center", va="center", fontsize=9)
    pts = P[good]
    try:
        hull = pts[ConvexHull(pts).vertices]
        ax.fill(*hull.T, color="#2b7bd6", alpha=0.25, zorder=1)
        ax.plot(*np.vstack([hull, hull[:1]]).T, color="#2b7bd6", lw=2, zorder=2)
    except (QhullError, ValueError):            # fewer than three points or all on one line
        ax.plot(*pts[np.lexsort(pts.T[::-1])].T, color="#2b7bd6", lw=2, zorder=2)
    feas = b["feasible"].to_numpy()
    ax.plot(*P[good].T, "o", color="#08306b", ms=6, zorder=3)
    ax.plot(*P[on & feas & ~good].T, "o", mfc="w", mec="#08306b", ms=6, zorder=3)
    ax.plot(*P[on & ~feas].T, "x", color="0.4", ms=6, zorder=3)
    for c, p in zip(cells[on], P[on]):
        ax.annotate(f"{c[0]:g},{c[1]:g},{c[2]:g}", p, xytext=(4, 3), textcoords="offset points", fontsize=6,
                    color="0.35")
    ax.set_title(title, fontsize=10)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(-0.65, 0.65); ax.set_ylim(-0.45, 0.72)


def main():
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser()
    p.add_argument("--stack", default="BLE")
    p.add_argument("--psr", type=float, default=0.99)
    p.add_argument("--csv", default=str(here / "results/triangle/psr_map_fixed.csv"))
    p.add_argument("--out", default=str(here / "results/triangle/perp_map_ble.png"))
    a = p.parse_args()
    df = pd.read_csv(a.csv)
    b = df[df["stack"] == a.stack]
    fig, axs = plt.subplots(1, 2, figsize=(14, 7.5))
    panel(axs[0], b, "R", a.psr, f"{a.stack}: networks axis grows to the lower right corner\n"
                                 "cells that can be drawn: sensors = datarate x networks")
    panel(axs[1], b, "T", a.psr, f"{a.stack}: networks axis grows to the top corner\n"
                                 "cells that can be drawn: datarate = sensors x networks")
    h = [plt.Line2D([], [], ls="", marker="o", color="#08306b", label=f"PSR >= {a.psr:g}"),
         plt.Line2D([], [], ls="", marker="o", mfc="w", mec="#08306b", label=f"PSR < {a.psr:g}"),
         plt.Line2D([], [], ls="", marker="x", color="0.4", label="no schedule fits"),
         plt.Rectangle((0, 0), 1, 1, fc="#2b7bd6", alpha=0.25, ec="#2b7bd6", label="region (outer passing cells)")]
    fig.legend(handles=h, loc="lower center", ncol=4, frameon=False,
               title="labels: sensors, datarate (x 10 kbit/s), networks")
    fig.savefig(a.out, dpi=150, bbox_inches="tight")
    print("saved", a.out)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Maps of the PSR = 0.99 boundary of the BAN stacks over three parameters.

Three parameters: datarate per sensor, number of networks (co-located BANs), sensors per network. The map consists
of three 4-sided regions; each is a 2D slice of the parameter space with the third parameter locked (LOCK).

Layouts:
  center        - minimum of all scales in the triangle center, the three lines from the center to the side
                  midpoints are the axes, scale labels on the outer sides (two half scales per side).
  corner        - minimum of the scales in the corners, each side one parameter, maximum at the side midpoint.
  circle        - the corner layout with the sides bent to arcs.
  circle_center - the center layout on a circle: minimum in the center, the three radii are the axes.
  polar         - minimum in the center, the axes on the three radii; inside a slice the angle is the mix of its two
                  parameters and the radius their sum in log units (the combined load), so equal-load lines are
                  arcs; the corner zone where both are large (u + v > 1) is pressed into an outer ring.
  disk          - polar with radius 1 - (1 - u)(1 - v): the axes keep their scale, every outer edge lies on one
                  circle, and the load arcs of polar bend only where both parameters are large.

Data: --csv with columns stack, datarate, n_networks, n_sensors, psr, psr_lo, feasible (one row per grid point),
or --build to write that CSV from the slice sweeps in results/triangle/ (run_triangle.sh). Per stack and grid point
the best configuration (PHY, period mode, BCC chip rate) by mean PSR is taken, with its lower 95 % limit.

Drawing: the coloured line closes the region where the stack works: a schedule fits and the PSR (or its lower
limit, --metric psr_lo) is at least 0.99. Earlier stacks are drawn wider, so coinciding lines stay visible.
--fill bcc shades the galvanic region, --fill all shades every region with its own hatch so overlaps stay readable.
"""
import argparse
import re
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from contourpy import contour_generator

PSR_TARGET = 0.99

# All scales logarithmic. Every axis shows the same numbers (lo/scale .. hi/scale); 'scale' divides the value for
# the labels, so datarate 10..1280 kbit/s reads 1..128 in units of 10 kbit/s.
AXES = {
    "datarate":   dict(label="datarate per sensor", unit="x 10 kbit/s", lo=10, hi=1280, scale=10, real_unit="kbit/s"),
    "n_networks": dict(label="number of networks", unit="", lo=1, hi=128, scale=1, real_unit=""),
    "n_sensors":  dict(label="sensors per network", unit="", lo=1, hi=128, scale=1, real_unit=""),
}
LOCK = {"datarate": 10, "n_networks": 1, "n_sensors": 1}
TICK_FACTOR = 2
GRID_FACTOR = 2

STYLES = [
    ("#2a78d6", "-"), ("#eb6834", (0, (6, 3))), ("#1baf7a", (0, (2, 2))),
    ("#eda100", (0, (8, 3, 2, 3))), ("#e87ba4", (0, (1, 3))),
    ("#008300", (0, (10, 4))), ("#6250d6", (0, (4, 2, 1, 2, 1, 2))),
]
GALVANIC_STYLE = ("#c0392b", "-")      # own colour so its fill is never read as a grey "no data" area
HATCHES = ["////", "\\\\\\\\", "||||", "----", "++", "xx", "..", "oo"]
GALVANIC = "Galvanic BCC"

NAMES = {"ble": "BLE", "ieee802154": "802.15.4 CSMA/CA", "tsch": "802.15.4 TSCH", "ieee802156nb": "802.15.6 NB",
         "smartban": "SmartBAN", "ant": "ANT", "hbc": "802.15.6 HBC capacitive", "bcc": GALVANIC + " (profile)"}

# ----------------------------------------------------------------- geometry
T = np.array([0.5, np.sqrt(3) / 2]); BL = np.array([0.0, 0.0]); BR = np.array([1.0, 0.0])
G = (T + BL + BR) / 3
ML, MB, MR = (BL + T) / 2, (BL + BR) / 2, (BR + T) / 2
RADIUS = np.hypot(*(T - G))


def unit(d):
    return d / np.hypot(*d)
LAYOUTS = ["center", "corner", "circle", "circle_center", "hex", "hex_center", "polar", "disk"]
POLAR_RING = 0.25


def hexed(p):
    """Point pushed out from G to the circle radius: the hexagon vertices of the hex layouts."""
    return G + RADIUS * unit(p - G)


def outer_corner(region, layout):
    return hexed(region[5]) if layout.startswith("hex") else region[5]
REGIONS = [
    ("datarate", "n_sensors", "n_networks", ML, MB, BL),
    ("n_sensors", "n_networks", "datarate", MB, MR, BR),
    ("datarate", "n_networks", "n_sensors", ML, MR, T),
]


def to_xy(u, v, region, layout):
    a, b, corner = region[3:]
    if layout.startswith("hex"):
        # three rhombi G, a, corner, b with all four outer points on the hexagon: the "cube" picture
        a, b, corner = hexed(a), hexed(b), hexed(corner)
        layout = "center" if layout == "hex_center" else "corner"
    u = np.asarray(u, float)[..., None]; v = np.asarray(v, float)[..., None]
    if layout in ("polar", "disk"):
        th_a, th_b = (np.arctan2(p[1] - G[1], p[0] - G[0]) for p in (a, b))
        sweep = (th_b - th_a + np.pi) % (2 * np.pi) - np.pi
        s = u + v
        th = th_a + sweep * np.divide(v, s, out=np.full_like(s, 0.5), where=s > 0)
        r = RADIUS * (np.where(s <= 1, s, 1 + POLAR_RING * (s - 1)) if layout == "polar" else 1 - (1 - u) * (1 - v))
        return G + r * np.concatenate([np.cos(th), np.sin(th)], axis=-1)
    if layout.startswith("circle"):
        th0, tha, thb = (np.arctan2(p[1] - G[1], p[0] - G[0]) for p in (corner, a, b))

        def arc(th1, t):
            th = th0 + ((th1 - th0 + np.pi) % (2 * np.pi) - np.pi) * t
            return G + RADIUS * np.concatenate([np.cos(th), np.sin(th)], axis=-1)
        if layout == "circle_center":
            return v * arc(thb, 1 - u) + u * arc(tha, 1 - v) + (1 - u) * (1 - v) * G - u * v * corner
        return (1 - v) * arc(tha, u) + (1 - u) * arc(thb, v) - (1 - u) * (1 - v) * corner + u * v * G
    origin, far = (G, corner) if layout == "center" else (corner, G)
    return origin * (1 - u) * (1 - v) + a * u * (1 - v) + b * (1 - u) * v + far * u * v


def rotate(p, pivot, M):
    """p mapped by the 2x2 matrix M about pivot (M None: unchanged)."""
    if M is None:
        return p
    q = np.asarray(p, float) - pivot
    return pivot + np.stack([M[0, 0] * q[..., 0] + M[0, 1] * q[..., 1],
                             M[1, 0] * q[..., 0] + M[1, 1] * q[..., 1]], axis=-1)


def mirror(phi_deg):
    """Reflection across the line through the origin at angle phi."""
    c, s_ = np.cos(np.radians(2 * phi_deg)), np.sin(np.radians(2 * phi_deg))
    return np.array([[c, s_], [s_, -c]])


def panels(layout, fan):
    """The drawn 2D slices: (region, value of its locked parameter, map about its outer corner).

    fan: datarates for the networks x sensors slice, smallest first. The first replaces the locked value of
    that slice; each further one is the previous slice mirrored across its clockwise edge at the outer vertex,
    so neighbouring slices share an edge that carries the same parameter on both sides, and the copies close
    the circle: triangle layouts (60-degree kites) take up to six, hex layouts (120-degree rhombi) up to three."""
    out = [(r, LOCK[r[2]], None) for r in REGIONS]
    if fan:
        if layout.startswith("circle"):
            raise SystemExit("--fan works with the triangle and hex layouts only")
        region = REGIONS[1]
        pivot = outer_corner(region, layout)
        ends = [panel_xy0(1.0, 0.0, region, layout), panel_xy0(0.0, 1.0, region, layout)]  # edges at the pivot
        ang = sorted(np.degrees(np.arctan2(e[1] - pivot[1], e[0] - pivot[0])) % 360 for e in ends)
        wedge = ang[1] - ang[0] if ang[1] - ang[0] <= 180 else 360 - (ang[1] - ang[0])
        n = int(round(360 / wedge))
        if len(fan) > n:
            raise SystemExit(f"--fan: at most {n} datarates in layout {layout}")
        out[1] = (region, fan[0], None)
        M = np.eye(2)
        edge = ang[0] if ang[1] - ang[0] <= 180 else ang[1]     # clockwise edge of the first slice
        for val in fan[1:]:
            M = mirror(edge) @ M
            out.append((region, val, M))
            edge -= wedge
    return out


def panel_xy0(u, v, region, layout):
    return to_xy(u, v, region, layout)


def panel_xy(u, v, panel, layout):
    region, _, M = panel
    return rotate(to_xy(u, v, region, layout), outer_corner(region, layout), M)


def log_pos(param, value):
    lo, hi = AXES[param]["lo"], AXES[param]["hi"]
    return np.log(np.asarray(value, float) / lo) / np.log(hi / lo)


def scale_values(param, factor):
    lo, hi = AXES[param]["lo"], AXES[param]["hi"]
    k = int(np.floor(np.log(hi / lo) / np.log(factor) + 1e-9))
    return [lo * factor ** i for i in range(k + 1)]


# -------------------------------------------------------------------- data
def build_csv(src, out, agg="best", bcc="profile"):
    """One row per stack and grid point from the slice sweeps (cells_*.csv in src).

    Every cell holds the mean PSR over its Monte-Carlo trials and the lower 95 % limit of that mean.
    agg = best:  per grid point the configuration (PHY, period mode, BCC chip rate) with the highest mean PSR;
                 neighbouring points may use different configurations.
    agg = fixed: one configuration per stack for the whole map, the one with the highest mean PSR over all
                 grid points of the three slices (infeasible points count as PSR 0).
    """
    import pandas as pd
    parts = []
    for f in sorted(Path(src).glob("cells_*.csv")):
        c = pd.read_csv(f)
        m = re.search(r"(500e3|2\.625e6)$", f.stem)
        c["rc"] = m.group(1) if m else ""
        c["src"] = str(f)
        parts.append(c)
    if not parts:
        raise SystemExit(f"no cells_*.csv in {Path(src).resolve()}: run the sweeps of results/triangle/ first, "
                         "or point --csv into that folder")
    d = pd.concat(parts, ignore_index=True)
    d = d[d.reception == "sinr"]
    # galvanic BCC is drawn as the profile (ladder level 5, filling pool); the plain HBC
    # stack on the galvanic channel (the simulator's default 'full' entry) is left out
    if bcc == "profile":
        d = d[(d.protocol != "bcc") | d.variant.astype(str).str.contains("opt5")]
    else:
        d = d[(d.protocol != "bcc") | ~d.variant.astype(str).str.contains("opt")]
    # phase-2 cells (cells_P2_*.csv, 30 trials on the 0.99 boundary) replace the phase-1 cells of the same configuration
    d["p2"] = [Path(f).name.startswith("cells_P2_") for f in d["src"]]
    same = ["protocol", "variant", "phy", "ci_mode", "r_target_kbps", "n_dev", "k"]
    d = d.sort_values("p2").drop_duplicates(same + ["rc"], keep="last")
    d["psr_eff"] = np.where(d.feasible == 1, d.psr_mean, -1.0)
    # the BCC chip rate is not part of the fixed configuration: the hub sets it for the load (as BLE picks 1M or
    # 2M), so with agg = fixed BCC keeps one PHY and period mode and takes the better chip rate at each point
    d["cfg"] = (d.variant.astype(str).str.replace(r"\+rc_cps=[0-9.e+]+", "", regex=True) + "/" + d.phy.astype(str)
                + "/" + d.ci_mode.astype(str))
    key = ["protocol", "r_target_kbps", "n_dev", "k"]
    if agg == "fixed":
        score = d.assign(s=d.psr_eff.clip(lower=0)).groupby(["protocol", "cfg"]).s.mean()
        keep = score.groupby(level=0).idxmax().tolist()
        d = d[[tuple(x) in set(keep) for x in d[["protocol", "cfg"]].values]]
        print("fixed configuration per stack:", ", ".join(f"{p} {c}" for p, c in keep))
    best = d.loc[d.groupby(key).psr_eff.idxmax()].copy()
    best["feasible"] = best.psr_eff >= 0
    best["stack"] = best["protocol"].map(NAMES)
    best = best.rename(columns={"r_target_kbps": "datarate", "k": "n_networks", "n_dev": "n_sensors",
                                "psr_mean": "psr"})
    best.loc[~best.feasible, ["psr", "psr_lo"]] = np.nan
    order = [NAMES[p] for p in NAMES]
    best["o"] = best["stack"].map(order.index)
    best = best.sort_values(["o", "datarate", "n_networks", "n_sensors"])
    best[["stack", "datarate", "n_networks", "n_sensors", "psr", "psr_lo", "feasible", "n_trials"]].to_csv(out, index=False)
    print("wrote", out, len(best), "rows")


def _pava_dec(y, w):
    """Least-squares non-increasing fit of y with weights w (pool adjacent violators)."""
    vals, wts, lens = [], [], []
    for yi, wi in zip(y, w):
        vals.append(yi); wts.append(wi); lens.append(1)
        while len(vals) > 1 and vals[-2] < vals[-1]:
            v2, w2, l2 = vals.pop(), wts.pop(), lens.pop()
            v1, w1, l1 = vals.pop(), wts.pop(), lens.pop()
            ww = w1 + w2
            vals.append((v1 * w1 + v2 * w2) / ww); wts.append(ww); lens.append(l1 + l2)
    return np.repeat(vals, lens)


def monotone(z, iters=50):
    """PSR made non-increasing along both axes of a slice (more load never helps) by isotonic regression:
    alternating least-squares fits along rows and columns. A single low cell among high ones is averaged with
    its neighbours instead of being copied to every cell of higher load (a running minimum would do that)."""
    m = np.ma.filled(z.astype(float), np.nan)
    ok = ~np.isnan(m)
    out = m.copy()
    for _ in range(iters):
        prev = out.copy()
        for axis in (0, 1):
            o = out if axis == 1 else out.T
            k = ok if axis == 1 else ok.T
            for i in range(o.shape[0]):
                idx = np.where(k[i])[0]
                if len(idx) > 1:
                    o[i, idx] = _pava_dec(o[i, idx], np.ones(len(idx)))
        if np.nanmax(np.abs(out - prev)) < 1e-6:
            break
    return np.ma.masked_invalid(out)


def slice_from_csv(df, stack, region, metric, lock=None):
    u_par, v_par, locked = region[:3]
    sl = df[(df["stack"] == stack) & np.isclose(df[locked], LOCK[locked] if lock is None else lock)]
    for par in (u_par, v_par):
        sl = sl[(sl[par] >= AXES[par]["lo"] * 0.999) & (sl[par] <= AXES[par]["hi"] * 1.001)]
    z = sl.pivot_table(index=v_par, columns=u_par, values=metric, dropna=False)
    f = sl.pivot_table(index=v_par, columns=u_par, values="feasible", aggfunc="max", dropna=False)
    # grid lines that only the neighbouring slices sampled (their shared edge has a finer grid): drop them
    f = f.reindex(index=z.index, columns=z.columns)
    while True:                    # drop the sparsest grid line, one at a time, until no line has a gap
        have = f.notna()
        r, c = have.mean(axis=1), have.mean(axis=0)
        if min(r.min(), c.min()) >= 1.0 or min(len(r), len(c)) < 2:
            break
        if r.min() <= c.min():
            keep = r.index != r.idxmin(); z, f = z.loc[keep], f.loc[keep]
        else:
            keep = c.index != c.idxmin(); z, f = z.loc[:, keep], f.loc[:, keep]
    # a slice needs its interior: the edge rows that other slices share (one network, one sensor) are not enough
    if z.shape[0] < 2 or z.shape[1] < 2 or f.reindex(index=z.index, columns=z.columns).notna().values.mean() < 0.5:
        print(f"skipped: {stack}, slice {u_par} x {v_par}")
        return None
    f = f.reindex(index=z.index, columns=z.columns)
    return (log_pos(u_par, z.columns.values), log_pos(v_par, z.index.values),
            np.ma.masked_invalid(z.values.astype(float)), f.values.astype(float))


# ------------------------------------------------------------------ drawing
def text_angle(d):
    a = np.degrees(np.arctan2(d[1], d[0]))
    return a - 180 if a > 90 else a + 180 if a <= -90 else a


def axis_title(par):
    a = AXES[par]
    return a["label"] + (", " + a["unit"] if a["unit"] else "")


def label(par, val):
    return f"{val / AXES[par]['scale']:g}"


def draw_frame(ax, layout, fan=None, lock_note=""):
    grid = dict(color="0.85", lw=0.5, zorder=1); frame = dict(color="0.45", lw=1.1, zorder=5)
    centered = layout in ("center", "circle_center", "hex_center"); round_ = layout.startswith("circle")
    edge = 1 if centered else 0
    s = np.linspace(0, 1, 60); seen = set()

    def curve(p, **kw):
        ax.plot(p[:, 0], p[:, 1], **kw)

    pts = []; edges_done = set()
    for panel in panels(layout, fan):
        region, lockval, deg = panel
        u_par, v_par, locked, a, b, corner = region
        pivot = outer_corner(region, layout)
        Gp = rotate(G, pivot, deg)          # centre of the figure as seen from this panel
        pts += [panel_xy(x, y, panel, layout) for x in (0, 1) for y in (0, 1)]
        for par, is_u in ((u_par, True), (v_par, False)):
            xy = (lambda u, v, P=panel: panel_xy(u, v, P, layout)) if is_u else \
                 (lambda u, v, P=panel: panel_xy(v, u, P, layout))
            for val in scale_values(par, GRID_FACTOR)[1:]:
                curve(xy(log_pos(par, val), s), **grid)
            curve(xy(s, edge), **frame); curve(xy(s, 1 - edge), **frame)
            e0, e1 = xy(0, edge), xy(1, edge)
            ekey = tuple(sorted((tuple(np.round(e0, 3)), tuple(np.round(e1, 3)))))
            if ekey in edges_done:
                continue                  # an edge shared by two fan slices carries one scale
            edges_done.add(ekey)

            def normal(p, e0=e0, e1=e1, Gp=Gp):
                if round_:
                    return unit(p - Gp)
                n = unit(np.array([e1[1] - e0[1], e0[0] - e1[0]]))
                return n if np.dot(n, (e0 + e1) / 2 - Gp) > 0 else -n

            for val in scale_values(par, TICK_FACTOR):
                t = log_pos(par, val); p = xy(t, edge); n = normal(p)
                curve(np.array([p, p + 0.012 * n]), **frame)
                q = p + 0.035 * n
                if layout != "corner" and t == 0:
                    q = q + 0.018 * unit(xy(0.02, edge) - e0)
                if layout in ("circle_center", "hex_center", "hex", "polar", "disk") and np.isclose(t, 1):
                    q = q - 0.03 * unit(e1 - xy(0.98, edge))
                key = (round(q[0], 3), round(q[1], 3), label(par, val))
                if fan and region is REGIONS[1] and np.isclose(t, 1):
                    key = ("fan pivot", label(par, val))     # every fan slice ends at the pivot: one label
                    q = pivot + 0.045 * unit(pivot - G)
                if key not in seen:
                    seen.add(key)
                    ax.text(q[0], q[1], label(par, val), ha="center", va="center", fontsize=8, color="0.4")
            anchor = xy(0.5, edge) if centered else e1
            n = normal(anchor); q = anchor + 0.085 * n
            direction = np.array([-n[1], n[0]]) if round_ else e1 - e0
            ha = "center"
            if layout in ("polar", "disk"):       # beside the outer part of the spoke, on the side away from its tick labels
                direction = e1 - e0; q = xy(0.72, edge) - 0.05 * normal(xy(0.72, edge))
            if centered or par not in seen:
                seen.add(par)
                ax.text(q[0], q[1], axis_title(par), ha=ha, va="center", rotation=text_angle(direction),
                        rotation_mode="anchor", fontsize=9, color="0.25")
        corner = outer_corner(region, layout)
        d = unit(corner - G)
        text = f"{AXES[locked]['label']} = {lockval:g} {AXES[locked]['real_unit']}".strip() + lock_note
        if fan and region is REGIONS[1]:
            # the fan's rhombi meet at their outer vertex: name each inside, between that vertex and its centre
            q = np.mean([panel_xy(x, y, panel, layout) for x in (0, 1) for y in (0, 1)], axis=0)
            ax.text(q[0], q[1], text.replace(" = ", "\n"), ha="center", va="center", fontsize=10,
                    fontweight="bold", color="0.35", alpha=0.75, zorder=6)
            continue
        if round_ or layout.startswith(("hex", "polar", "disk")):
            q = corner + (0.17 + (0.1 if layout == "polar" else 0.05 if layout == "disk" else 0)) * d; ha = "center" if abs(d[0]) < 0.1 else "right" if d[0] < 0 else "left"
        else:
            q = corner + 0.075 * d + np.array([0, 0.035 if corner is T else -0.085]); ha = "center"
        ax.text(q[0], q[1], text, ha=ha, va="center", fontsize=9, fontweight="bold")
    return np.array(pts)


def draw(layout, stacks, df, metric, fill, out, mono=False, fan=None):
    fig, ax = plt.subplots(figsize=(10, 10) if not fan else (10.5, 9))
    frame_pts = draw_frame(ax, layout, fan)
    handles = []
    for i, stack in enumerate(stacks):
        galv = stack.startswith(GALVANIC)
        color, ls = GALVANIC_STYLE if galv else STYLES[i % len(STYLES)]
        # earlier stacks drawn wider, so a line that coincides with a later one shows as an outline around it
        lw = 1.5 + 0.3 * (len(stacks) - 1 - i)
        for panel in panels(layout, fan):
            region, lockval, _ = panel
            data = slice_from_csv(df, stack, region, metric, lockval)
            if data is None:
                continue
            u, v, z, feas = data
            if mono:
                z = monotone(np.ma.masked_where(np.nan_to_num(feas, nan=0.0) <= 0, z))
            uu, vv = np.meshgrid(u, v)
            xy = panel_xy(uu, vv, panel, layout)
            feas = np.nan_to_num(feas, nan=0.0)
            zv = np.ma.filled(z.astype(float), np.nan)
            # region where the stack works: feasible and PSR >= target (infeasible counts as PSR 0)
            zwork = np.where(feas > 0, np.nan_to_num(zv, nan=0.0), 0.0)
            if (fill == "all" or (fill == "bcc" and galv)) and (zwork >= PSR_TARGET).any():
                ax.contourf(xy[..., 0], xy[..., 1], zwork, levels=[PSR_TARGET, 1.01], colors=[color],
                            alpha=0.18 if galv else 0.07, zorder=2)
                if not galv:
                    cs = ax.contourf(xy[..., 0], xy[..., 1], zwork, levels=[PSR_TARGET, 1.01], colors="none",
                                     hatches=[HATCHES[i % len(HATCHES)]], zorder=2)
                    cs.set_edgecolor(color); cs.set_linewidth(0)
            # one closed line per region: PSR drop and capacity wall alike (infeasible counts as PSR 0)
            for line in contour_generator(x=u, y=v, z=zwork).lines(PSR_TARGET):
                p = panel_xy(line[:, 0], line[:, 1], panel, layout)
                ax.plot(p[:, 0], p[:, 1], color=color, ls=ls, lw=lw, zorder=4)
        handles.append(plt.Line2D([], [], color=color, ls=ls, lw=lw, label=stack))
    round_ = layout.startswith(("circle", "hex"))
    place = dict(loc="upper right", bbox_to_anchor=(0.0, 1.0)) if round_ else dict(loc="upper left")
    title = f"PSR = {PSR_TARGET:g}" + (" (lower 95 % limit)" if metric == "psr_lo" else "")
    ax.legend(handles=handles, frameon=False, fontsize=9, handlelength=4.5, title=title, title_fontsize=9,
              alignment="left", **place)
    ax.set_aspect("equal"); ax.axis("off")
    if fan:
        lo, hi = frame_pts.min(axis=0) - 0.2, frame_pts.max(axis=0) + 0.2
        ax.set_xlim(min(lo[0], -0.32), hi[0]); ax.set_ylim(lo[1], max(hi[1], 1.0))
    else:
        ax.set_xlim(-0.32, 1.32); ax.set_ylim(-0.46 if round_ else -0.22, 1.0)
    fig.savefig(out, dpi=200, bbox_inches="tight"); plt.close(fig); print("saved", out)


LEVELS = [(1e-9, "any packet delivered"), (0.5, "PSR >= 0.5"), (0.9, "PSR >= 0.9"), (0.99, "PSR >= 0.99")]
LEVEL_COLORS = ["#dbe9f6", "#9ecae1", "#4292c6", "#08306b"]


def draw_levels(layout, stacks, df, metric, out, mono=True, fan=None, ncol=4):
    """One map per stack, the regions where its PSR reaches each level of LEVELS filled from light to dark:
    the widest region (anything delivered) first, the 0.99 region (the strictest use) last and smallest."""
    nrow = int(np.ceil(len(stacks) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 4.9 * nrow), squeeze=False)
    for ax, stack in zip(axes.flat, stacks):
        draw_frame(ax, layout, fan)
        for panel in panels(layout, fan):
            region, lockval, _ = panel
            data = slice_from_csv(df, stack, region, metric, lockval)
            if data is None:
                continue
            u, v, z, feas = data
            if mono:
                z = monotone(np.ma.masked_where(np.nan_to_num(feas, nan=0.0) <= 0, z))
            feas = np.nan_to_num(feas, nan=0.0)
            zwork = np.where(feas > 0, np.nan_to_num(np.ma.filled(z.astype(float), np.nan), nan=0.0), -1.0)
            uu, vv = np.meshgrid(u, v)
            xy = panel_xy(uu, vv, panel, layout)
            bounds = [lv for lv, _ in LEVELS] + [1.01]
            ax.contourf(xy[..., 0], xy[..., 1], zwork, levels=bounds, colors=LEVEL_COLORS, zorder=2)
            for line in contour_generator(x=u, y=v, z=feas).lines(0.5):     # capacity wall
                q = panel_xy(line[:, 0], line[:, 1], panel, layout)
                ax.plot(q[:, 0], q[:, 1], color="0.2", lw=1.0, zorder=3)
        ax.set_title(stack, fontsize=12, fontweight="bold")
        ax.set_aspect("equal"); ax.axis("off")
        ax.set_xlim(-0.32, 1.32); ax.set_ylim(-0.22, 1.0)
    for ax in axes.flat[len(stacks):]:
        ax.axis("off")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in LEVEL_COLORS] + [plt.Line2D([], [], color="0.2", lw=1.0)]
    labels = [t for _, t in LEVELS] + ["no schedule fits beyond"]
    fig.legend(handles, labels, loc="lower center", ncol=5, frameon=False, fontsize=11)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close(fig); print("saved", out)


# Locked parameter in units u of its LOCK value (1 sensor, 1 network, 10 kbit/s), up to the end of its axis (128u).
# Steps above 1u come from the factor-2 3D grid (run_fan*.sh, lockgrid.py).
LOCK_STEPS = [1, 2, 4, 8, 16, 32, 64, 128]


def with_pruned(df, floor=0.5):
    """df plus the factor-2 grid cells lockgrid.py skipped because a lower neighbour (one step down in one
    parameter) fails below floor. More load never raises the PSR, so such a cell's PSR is at most the lowest PSR
    of its lower neighbours; that bound is entered as its PSR (no schedule if the bound is 0). Every PSR target of
    the lock maps is >= floor, so the cell fails all of them, as it would with its true value."""
    import itertools as it
    import pandas as pd
    from lockgrid import N, R, lower
    rows = []
    for stack, s in df.groupby("stack", sort=False):
        psr = {(r, n, k): (p if f else 0.0) for r, n, k, p, f in
               zip(s.datarate, s.n_sensors, s.n_networks, s.psr.fillna(0.0), s.feasible)}
        idx = {v: i for i, v in enumerate(N)} | {v: i for i, v in enumerate(R)}
        for c in sorted(it.product(R, N, N), key=lambda c: idx[c[0]] + idx[c[1]] + idx[c[2]]):
            lows = [psr[l] for l in lower(c) if l in psr]
            if c in psr or not lows or min(lows) >= floor:
                continue
            psr[c] = b = min(lows)
            rows.append(dict(stack=stack, datarate=c[0], n_sensors=c[1], n_networks=c[2],
                             psr=b if b > 0 else np.nan, psr_lo=b if b > 0 else np.nan, feasible=b > 0))
    return pd.concat([df, pd.DataFrame(rows, columns=df.columns)], ignore_index=True)


def lock_count(df, stack, region, metric, thr, mono, n=161):
    """Per point of an n x n grid over the slice: how many lock steps in a row, from 1u up, reach PSR >= thr
    (a step that fails ends the count, so the regions of higher steps nest inside the lower ones).
    Outside the grid a step was simulated on, that step counts as failed; extent[i] is its (u_max, v_max)."""
    from scipy.interpolate import RegularGridInterpolator
    s = np.linspace(0, 1, n)
    uu, vv = np.meshgrid(s, s)
    count = np.zeros_like(uu); alive = np.ones_like(uu, bool); extent = []
    for step in LOCK_STEPS:
        data = slice_from_csv(df, stack, region, metric, LOCK[region[2]] * step)
        if data is None:
            break
        u, v, z, feas = data
        feas = np.nan_to_num(feas, nan=0.0)
        if mono:
            z = monotone(np.ma.masked_where(feas <= 0, z))
        zwork = np.where(feas > 0, np.nan_to_num(np.ma.filled(z.astype(float), np.nan), nan=0.0), 0.0)
        zi = RegularGridInterpolator((v, u), zwork, bounds_error=False, fill_value=np.nan)(np.stack([vv, uu], -1))
        alive &= zi >= thr
        count += alive
        extent.append((u.max(), v.max()))
    return uu, vv, count, extent


def draw_coverage(ax, panel, layout, extent):
    """Dotted outline of the area where the steps above 1u were simulated."""
    if len(extent) < 2 or min(extent[1]) >= 1 - 1e-9:
        return
    um, vm = extent[1]; s = np.linspace(0, 1, 40)
    for p in (panel_xy(um, s * vm, panel, layout), panel_xy(s * um, vm, panel, layout)):
        ax.plot(p[:, 0], p[:, 1], color="0.35", ls=":", lw=1.0, zorder=5)


def step_label(i):
    return f"{LOCK_STEPS[i]}u"


def draw_lockmap(stacks, df, metric, out, mono, thrs=(PSR_TARGET,), layout="hex_center"):
    """Two figures on the isometric layout (hex_center, or polar), one hexagon per stack.

    out_lock:  the PSR >= thrs[0] region at 1u filled, the same boundary at 2u, 4u ... dashed.
    out_steps: one hexagon per stack and threshold, coloured by how many lock steps from 1u up reach it."""
    K = len(LOCK_STEPS)
    cache = {}

    def counts(stack, panel, thr):
        key = (stack, panel[0][2], thr)
        if key not in cache:
            cache[key] = lock_count(df, stack, panel[0], metric, thr, mono)
        return cache[key]

    def finish(ax, title):
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.set_aspect("equal"); ax.axis("off")
        ax.set_xlim(-0.62, 1.62); ax.set_ylim(-0.42 - 0.12 * (layout == "polar"), 1.1 + 0.08 * (layout == "polar"))

    fig, axes = plt.subplots(1, len(stacks), figsize=(10 * len(stacks), 7.5), squeeze=False)
    for i, (ax, stack) in enumerate(zip(axes.flat, stacks)):
        color = GALVANIC_STYLE[0] if stack.startswith(GALVANIC) else STYLES[i % len(STYLES)][0]
        draw_frame(ax, layout, lock_note=" (1u)")
        for panel in panels(layout, None):
            uu, vv, count, extent = counts(stack, panel, thrs[0])
            xy = panel_xy(uu, vv, panel, layout)
            X, Y = xy[..., 0], xy[..., 1]
            if (count >= 1).any():
                ax.contourf(X, Y, count, levels=[0.5, K + 1], colors=[color], alpha=0.3, zorder=2)
                ax.contour(X, Y, count, levels=[0.5], colors=[color], linewidths=2.0, zorder=4)
            lv = [k + 0.5 for k in range(1, K) if (count >= k + 1).any()]
            if lv:
                cs = ax.contour(X, Y, count, levels=lv, colors=[color], linewidths=1.3, linestyles=[(0, (5, 3))],
                                alpha=0.6, zorder=4)
                ax.clabel(cs, fmt={l: step_label(int(l)) for l in lv}, fontsize=8, inline=True)
            draw_coverage(ax, panel, layout, extent)
        finish(ax, stack)
    handles = [plt.Rectangle((0, 0), 1, 1, color="0.5", alpha=0.3),
               plt.Line2D([], [], color="0.3", ls=(0, (5, 3)), alpha=0.6),
               plt.Line2D([], [], color="0.35", ls=":")]
    labels = [f"PSR >= {thrs[0]:g} with the corner parameter at 1u",
              "same boundary at 2u, 4u ... (u = 1 sensor, 1 network, 10 kbit/s)",
              "steps above 1u simulated inside this line only"]
    fig.legend(handles, labels, loc="lower center", ncol=1, frameon=False, fontsize=10)
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(f"{out}_lock.png", dpi=200, bbox_inches="tight"); plt.close(fig); print("saved", f"{out}_lock.png")

    colors = plt.get_cmap("Blues")(np.linspace(0.25, 1.0, K))
    fig, axes = plt.subplots(len(stacks), len(thrs), figsize=(10 * len(thrs), 7 * len(stacks)), squeeze=False)
    for row, stack in zip(axes, stacks):
        for ax, thr in zip(row, thrs):
            draw_frame(ax, layout, lock_note=" (1u)")
            for panel in panels(layout, None):
                uu, vv, count, extent = counts(stack, panel, thr)
                xy = panel_xy(uu, vv, panel, layout)
                if (count >= 1).any():
                    ax.contourf(xy[..., 0], xy[..., 1], count, levels=np.arange(K + 1) + 0.5, colors=colors,
                                zorder=2)
                draw_coverage(ax, panel, layout, extent)
            finish(ax, f"{stack}, PSR >= {thr:g}")
    from matplotlib.colors import BoundaryNorm, ListedColormap
    sm = plt.cm.ScalarMappable(cmap=ListedColormap(colors), norm=BoundaryNorm(np.arange(K + 1) + 0.5, K))
    cb = fig.colorbar(sm, ax=axes.ravel().tolist(), ticks=np.arange(1, K + 1), orientation="horizontal",
                      shrink=0.45, pad=0.0, aspect=30)
    cb.ax.set_xticklabels([step_label(k) for k in range(K)])
    cb.set_label("PSR target met at every step from 1u up to\n"
                 "(dotted line: steps above 1u simulated inside it only)", fontsize=10)
    fig.savefig(f"{out}_steps.png", dpi=200, bbox_inches="tight"); plt.close(fig); print("saved", f"{out}_steps.png")


def draw_psr_variants(stack, df, metric, out, mono, layout, hi=PSR_TARGET, lo=0.9):
    """One IPS of one stack, four ways to show two PSR targets on it (the lock-step counts of lock_count)."""
    from matplotlib.colors import BoundaryNorm, ListedColormap
    K = len(LOCK_STEPS)
    lv = np.arange(K + 1) + 0.5
    blues = plt.get_cmap("Blues")(np.linspace(0.25, 1.0, K))
    oranges = plt.get_cmap("Oranges")(np.linspace(0.2, 0.85, K))
    ORANGE = "#c2410c"
    fmt = {l: step_label(int(l)) for l in lv[:-1]}
    data = [(panel, lock_count(df, stack, panel[0], metric, hi, mono), lock_count(df, stack, panel[0], metric, lo, mono)[2])
            for panel in panels(layout, None)]
    titles = [f"A  fill: steps at {hi:g};  dashed lines: steps at {lo:g}",
              f"B  orange fill: steps at {lo:g};  blue fill on top: steps at {hi:g}",
              f"C  fill: steps at {hi:g};  hatch: extra steps when {lo:g} is enough (/ +1, // +2, xx +3 or more)",
              f"D  lines only: {hi:g} blue, {lo:g} orange, 1u regions lightly filled"]
    fig, axes = plt.subplots(2, 2, figsize=(20, 15))
    for i, ax in enumerate(axes.flat):
        draw_frame(ax, layout, lock_note=" (1u)")
        for panel, (uu, vv, c_hi, extent), c_lo in data:
            xy = panel_xy(uu, vv, panel, layout)
            X, Y = xy[..., 0], xy[..., 1]
            if i == 0:
                ax.contourf(X, Y, c_hi, levels=lv, colors=blues, zorder=2)
                cs = ax.contour(X, Y, c_lo, levels=lv[:-1], colors=ORANGE, linewidths=1.1,
                                linestyles=[(0, (5, 3))], zorder=4)
                ax.clabel(cs, fmt=fmt, fontsize=7, inline=True)
            elif i == 1:
                ax.contourf(X, Y, c_lo, levels=lv, colors=oranges, zorder=2)
                ax.contourf(X, Y, np.ma.masked_less(c_hi, 0.5), levels=lv, colors=blues, zorder=3)
            elif i == 2:
                ax.contourf(X, Y, c_hi, levels=lv, colors=blues, zorder=2)
                cs = ax.contourf(X, Y, c_lo - c_hi, levels=[0.5, 1.5, 2.5, K + 1], colors="none",
                                 hatches=["/", "//", "xx"], zorder=3)
                cs.set_edgecolor(ORANGE); cs.set_linewidth(0)
            else:
                ax.contourf(X, Y, c_lo, levels=[0.5, K + 1], colors=[ORANGE], alpha=0.12, zorder=2)
                ax.contourf(X, Y, c_hi, levels=[0.5, K + 1], colors=[STYLES[0][0]], alpha=0.25, zorder=2)
                for c, col, w in ((c_lo, ORANGE, 1.1), (c_hi, STYLES[0][0], 1.4)):
                    cs = ax.contour(X, Y, c, levels=lv[:-1], colors=col, linewidths=w, zorder=4,
                                    linestyles=["-"] + [(0, (5, 3))] * (K - 1))
                    ax.clabel(cs, fmt=fmt, fontsize=7, inline=True)
            draw_coverage(ax, panel, layout, extent)
        ax.set_title(titles[i], fontsize=11, loc="left")
        ax.set_aspect("equal"); ax.axis("off")
        ax.set_xlim(-0.62, 1.62); ax.set_ylim(-0.42 - 0.12 * (layout == "polar"), 1.1 + 0.08 * (layout == "polar"))
    for cols, label, x in ((blues, f"PSR >= {hi:g} met at every step from 1u up to", 0.12),
                           (oranges, f"PSR >= {lo:g} met at every step from 1u up to", 0.56)):
        cax = fig.add_axes([x, 0.035, 0.32, 0.012])
        sm = plt.cm.ScalarMappable(cmap=ListedColormap(cols), norm=BoundaryNorm(lv, K))
        cb = fig.colorbar(sm, cax=cax, orientation="horizontal", ticks=np.arange(1, K + 1))
        cb.ax.set_xticklabels([step_label(k) for k in range(K)]); cb.set_label(label, fontsize=10)
    fig.suptitle(f"{stack}: two PSR targets on one isometric Pareto surface ({layout})", fontsize=14, fontweight="bold")
    fig.subplots_adjust(left=0.02, right=0.98, top=0.95, bottom=0.08, wspace=0.18, hspace=0.08)
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig); print("saved", out)


PSR_BINS = [0.5, 0.9, 0.99, 0.999]
PSR_COLORS = ["#b2182b", "#ef8a62", "#fddbc7", "#67a9cf", "#2166ac"]     # < 0.5 ... >= 0.999
LEVEL_CMAPS = ["Oranges", "Blues", "Greens", "Purples"]


def node_psr(df, stack, region, metric, mono):
    """PSR at the factor-2 grid nodes of a slice for every lock step: array [step, v node, u node], NaN where
    that step was not simulated; a configuration without a fitting schedule counts as PSR 0."""
    u_par, v_par = region[:2]
    un, vn = log_pos(u_par, scale_values(u_par, 2)), log_pos(v_par, scale_values(v_par, 2))
    P = np.full((len(LOCK_STEPS), len(vn), len(un)), np.nan)
    for k, step in enumerate(LOCK_STEPS):
        data = slice_from_csv(df, stack, region, metric, LOCK[region[2]] * step)
        if data is None:
            continue
        u, v, z, feas = data
        feas = np.nan_to_num(feas, nan=0.0)
        if mono:
            z = monotone(np.ma.masked_where(feas <= 0, z))
        zw = np.where(feas > 0, np.nan_to_num(np.ma.filled(z.astype(float), np.nan), nan=0.0), 0.0)
        iu = [np.flatnonzero(np.isclose(u, x)) for x in un]
        iv = [np.flatnonzero(np.isclose(v, y)) for y in vn]
        for j, a in enumerate(iv):
            for i, b in enumerate(iu):
                if len(a) and len(b):
                    P[k, j, i] = zw[a[0], b[0]]
    return un, vn, P


def cell_poly(u0, u1, v0, v1, panel, layout, n=6):
    """Outline of the (u, v) box [u0, u1] x [v0, v1] on the map, edges sampled so curved layouts stay curved."""
    t = np.linspace(0, 1, n)
    uu = np.concatenate([u0 + (u1 - u0) * t, np.full(n, u1), u1 - (u1 - u0) * t, np.full(n, u0)])
    vv = np.concatenate([np.full(n, v0), v0 + (v1 - v0) * t, np.full(n, v1), v1 - (v1 - v0) * t])
    return panel_xy(uu, vv, panel, layout)


def draw_cells(stack, df, metric, out, mono, layout, mode, thrs):
    """Every grid node of a slice drawn as a cell, the cell split into sub-cells.

    mode levels: sub-cells interleave the PSR targets thrs (2 targets: 3 x 3 checker, n targets: n x n Latin
    square); a sub-cell shows in its target's colour ramp how many lock steps from 1u up reach that target.
    mode steps:  2 rows of sub-cells, one per lock step (1u top left ... 128u bottom right), coloured by PSR band."""
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import BoundaryNorm, ListedColormap
    K = len(LOCK_STEPS)
    if mode == "levels":
        n = len(thrs)
        pattern = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]]) if n == 2 else (np.add.outer(np.arange(n), np.arange(n)) % n)
        ramps = [plt.get_cmap(LEVEL_CMAPS[i % len(LEVEL_CMAPS)])(np.linspace(0.2, 1.0, K)) for i in range(n)]
    else:
        pattern = np.arange(K).reshape(2, -1)
    rows, cols = pattern.shape
    fig, ax = plt.subplots(figsize=(10, 9))
    draw_frame(ax, layout, lock_note=" (1u)")
    polys, faces = [], []
    for panel in panels(layout, None):
        un, vn, P = node_psr(df, stack, panel[0], metric, mono)
        du, dv = un[1] - un[0], vn[1] - vn[0]
        for j, v in enumerate(vn):
            for i, u in enumerate(un):
                u0, u1 = max(u - du / 2, 0), min(u + du / 2, 1)
                v0, v1 = max(v - dv / 2, 0), min(v + dv / 2, 1)
                col = P[:, j, i]
                for r in range(rows):
                    for c in range(cols):
                        idx = pattern[r, c]
                        if mode == "levels":
                            cnt = np.cumprod(np.nan_to_num(col, nan=-1) >= thrs[idx]).sum()
                            face = ramps[idx][cnt - 1] if cnt else (1, 1, 1, 1)
                        else:
                            p = col[idx]
                            face = (0.93, 0.93, 0.93, 1) if np.isnan(p) else PSR_COLORS[np.searchsorted(PSR_BINS, p, side="right")]
                        a0, a1 = u0 + (u1 - u0) * c / cols, u0 + (u1 - u0) * (c + 1) / cols
                        b0, b1 = v1 - (v1 - v0) * (r + 1) / rows, v1 - (v1 - v0) * r / rows
                        polys.append(cell_poly(a0, a1, b0, b1, panel, layout)); faces.append(face)
    ax.add_collection(PolyCollection(polys, facecolors=faces, edgecolors="none", zorder=2))
    edges = []
    for panel in panels(layout, None):
        un, vn = log_pos(panel[0][0], scale_values(panel[0][0], 2)), log_pos(panel[0][1], scale_values(panel[0][1], 2))
        du, dv = un[1] - un[0], vn[1] - vn[0]
        for v in vn:
            for u in un:
                edges.append(cell_poly(max(u - du / 2, 0), min(u + du / 2, 1), max(v - dv / 2, 0), min(v + dv / 2, 1),
                                       panel, layout))
    ax.add_collection(PolyCollection(edges, facecolors="none", edgecolors="0.55", linewidths=0.5, zorder=3))
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(-0.62, 1.62); ax.set_ylim(-0.42 - 0.12 * (layout == "polar"), 1.1 + 0.08 * (layout == "polar"))
    lv = np.arange(K + 1) + 0.5
    if mode == "levels":
        for i, thr in enumerate(thrs):
            cax = fig.add_axes([0.2, 0.06 - 0.045 * i, 0.6, 0.012])
            sm = plt.cm.ScalarMappable(cmap=ListedColormap(ramps[i]), norm=BoundaryNorm(lv, K))
            cb = fig.colorbar(sm, cax=cax, orientation="horizontal", ticks=np.arange(1, K + 1))
            cb.ax.set_xticklabels([step_label(k) for k in range(K)], fontsize=8)
            cb.ax.set_ylabel(f">= {thr:g}", rotation=0, ha="right", va="center", fontsize=9)
        pat = "checker 3 x 3" if len(thrs) == 2 else f"Latin square {len(thrs)} x {len(thrs)}"
        title = f"{stack}: cells split by PSR target ({pat}); colour = steps from 1u up meeting the target, white = none"
    else:
        cax = fig.add_axes([0.2, 0.06, 0.6, 0.014])
        bounds = [0] + PSR_BINS + [1.0001]
        sm = plt.cm.ScalarMappable(cmap=ListedColormap(PSR_COLORS), norm=BoundaryNorm(bounds, len(PSR_COLORS)))
        cb = fig.colorbar(sm, cax=cax, orientation="horizontal", ticks=[0.25, 0.7, 0.945, 0.9945, 0.9995])
        cb.ax.set_xticklabels(["< 0.5", "0.5-0.9", "0.9-0.99", "0.99-0.999", ">= 0.999"])
        cb.set_label("PSR (grey: not simulated)", fontsize=9)
        half = K // 2
        title = (f"{stack}: cells split by lock step ({' '.join(map(step_label, range(half)))} / "
                 f"{' '.join(map(step_label, range(half, K)))}), colour = PSR")
    fig.suptitle(title, fontsize=11, fontweight="bold", y=0.97)
    fig.savefig(out, dpi=200, bbox_inches="tight"); plt.close(fig); print("saved", out)


def spiral_centres(n, d):
    """Hex cell centres along a spiral: the middle, then ring after ring anticlockwise from the right; the first
    cell of a new ring touches the last and the first cell of the ring before it."""
    dirs = [d * np.array([np.cos(np.radians(a)), np.sin(np.radians(a))]) for a in range(0, 360, 60)]
    out, ring = [np.zeros(2)], 1
    while len(out) < n:
        p = ring * dirs[0]
        walk = []
        for k in range(6):
            for _ in range(ring):
                walk.append(p.copy()); p = p + dirs[(k + 2) % 6]
        if ring > 1:            # start at the cell just clockwise of the right one: it touches the previous ring's end
            walk = walk[-1:] + walk[:-1]
        out += walk
        ring += 1
    return out[:n]


def stack_ramp(color, K):
    from matplotlib.colors import LinearSegmentedColormap, to_rgb
    c = np.array(to_rgb(color))
    return LinearSegmentedColormap.from_list("", [0.82 + 0.18 * c, c, 0.45 * c])(np.linspace(0, 1, K))


def draw_honeycomb(stacks, df, metric, out, mono, thrs=(PSR_TARGET,), hues=False, width_in=None):
    """IPS on a hex spiral: all stacks at the first PSR target (the first stack in the middle), then all stacks
    again in the same order at the next target in the next colour, and so on. Opposite edges of an IPS carry the
    same parameter with the same scale direction, so touching edges of neighbours match without rotation.
    hues (one target only): every stack in its own colour, light to dark = more lock steps.
    width_in: figure width in inches (default about 14), so the fonts keep their size in print at that width."""
    from matplotlib.colors import BoundaryNorm, ListedColormap
    K = len(LOCK_STEPS)
    lv = np.arange(K + 1) + 0.5
    ramps = [plt.get_cmap(["Blues", "Oranges", "Greens", "Purples"][i % 4])(np.linspace(0.25, 1.0, K))
             for i in range(len(thrs))]
    items = [(stack, i) for i in range(len(thrs)) for stack in stacks]
    layout = "hex_center"
    d = np.sqrt(3) * RADIUS
    centres = spiral_centres(len(items), d)

    def occupied(p):
        return any(np.hypot(*(p - c)) < 1e-6 for c in centres)

    cs = np.array(centres)
    span = cs.max(axis=0) - cs.min(axis=0) + 2 * RADIUS + 0.5
    w = width_in or 15 * span[0] / 3.0
    fig, ax = plt.subplots(figsize=(w, w * span[1] / span[0]))
    pt = span[0] / w / 72                   # data units per point of the saved figure
    big = width_in is None
    fs = dict(tick=7, title=8, name=10, bar=9) if big else dict(tick=5.5, title=5.5, name=7, bar=6)
    short = {"datarate": "data rate (10 kbit/s)", "n_networks": "networks", "n_sensors": "sensors per network"}
    s = np.linspace(0, 1, 40)
    for (stack, li), c in zip(items, centres):
        off = c - G
        ramp, thr = ramps[li], thrs[li]
        if hues:
            ramp = stack_ramp(GALVANIC_STYLE[0] if stack.startswith(GALVANIC)
                              else STYLES[list(NAMES.values()).index(stack) % len(STYLES)][0], K)
        for panel in panels(layout, None):
            region = panel[0]
            uu, vv, count, extent = lock_count(df, stack, region, metric, thr, mono)
            xy = panel_xy(uu, vv, panel, layout) + off
            if (count >= 1).any():
                ax.contourf(xy[..., 0], xy[..., 1], count, levels=lv, colors=ramp, zorder=2)
            for par, is_u in ((region[0], True), (region[1], False)):
                for val in scale_values(par, GRID_FACTOR)[1:]:
                    t = log_pos(par, val)
                    p = (panel_xy(t, s, panel, layout) if is_u else panel_xy(s, t, panel, layout)) + off
                    ax.plot(p[:, 0], p[:, 1], color="0.85", lw=0.4, zorder=1)
            for p in (panel_xy(s, 0, panel, layout), panel_xy(0, s, panel, layout)):
                p = p + off
                ax.plot(p[:, 0], p[:, 1], color="0.5", lw=0.8, zorder=3)
            # outer edges: u = 1 carries the v parameter, v = 1 the u parameter; scales only where no neighbour
            for par, edge in ((region[1], lambda t: panel_xy(1.0, t, panel, layout)),
                              (region[0], lambda t: panel_xy(t, 1.0, panel, layout))):
                e = edge(s) + off
                ax.plot(e[:, 0], e[:, 1], color="0.3", lw=1.2, zorder=4)
                n = unit((e[0] + e[-1]) / 2 - c)
                if occupied(c + d * n):
                    continue
                for val in scale_values(par, 4):
                    q = edge(log_pos(par, val)) + off
                    ax.plot(*np.array([q, q + 0.6 * fs["tick"] * pt * n]).T, color="0.3", lw=1.0, zorder=4)
                    ax.text(*(q + 1.55 * fs["tick"] * pt * n), label(par, val), ha="center", va="center",
                            fontsize=fs["tick"], color="0.4")
                title = axis_title(par) if big else short[par]
                m = (e[0] + e[-1]) / 2 + (3.4 if big else 2.9) * fs["title"] * pt * n
                if not big:     # move the title away from an inner corner, where the neighbour's title ends
                    inner = [int(sum(np.hypot(*(cc - v)) < RADIUS + 1e-6 for cc in centres) > 1) for v in (e[0], e[-1])]
                    m = m + 0.2 * (inner[0] - inner[1]) * (e[-1] - e[0])
                ax.text(m[0], m[1], title, ha="center", va="center",
                        fontsize=fs["title"], color="0.25", rotation=text_angle(e[-1] - e[0]), rotation_mode="anchor")
        name = stack.replace(" (profile)", "")
        name = name + f"\nPSR >= {thr:g}" if big else textwrap.fill(name, 14)
        ax.text(c[0], c[1] - 0.07, name, ha="center", va="center", fontsize=fs["name"],
                fontweight="bold", zorder=6, bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.6", alpha=0.85))
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(cs[:, 0].min() - RADIUS - 0.25, cs[:, 0].max() + RADIUS + 0.25)
    ax.set_ylim(cs[:, 1].min() - RADIUS - 0.25, cs[:, 1].max() + RADIUS + 0.25)
    for i, thr in enumerate(thrs):
        cax = fig.add_axes([0.3, (0.075 if big else 0.12) - 0.03 * i, 0.4, 0.008 if big else 0.012])
        bar = plt.get_cmap("Greys")(np.linspace(0.2, 0.9, K)) if hues else ramps[i]
        sm = plt.cm.ScalarMappable(cmap=ListedColormap(bar), norm=BoundaryNorm(lv, K))
        cb = fig.colorbar(sm, cax=cax, orientation="horizontal", ticks=np.arange(1, K + 1))
        cb.ax.set_xticklabels([step_label(k) for k in range(K)], fontsize=fs["bar"])
        cb.ax.set_ylabel(f"PSR >= {thr:g}", rotation=0, ha="right", va="center", fontsize=fs["bar"] + 1)
    if width_in is None:    # the note is about 10 in long at 10 pt; a page-width figure leaves it to the caption
        fig.text(0.5, 0.075 - 0.03 * len(thrs) - 0.005,
                 "colour depth: target met at every step from 1u up to; rhombus at the top: sensors per network locked, "
                 "lower left: networks, lower right: datarate (1u = 1 sensor, 1 network, 10 kbit/s)",
                 ha="center", va="top", fontsize=10)
    fig.savefig(out, dpi=200, bbox_inches="tight"); plt.close(fig); print("saved", out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--layout", choices=LAYOUTS + ["all"], default="all")
    ap.add_argument("--csv", default="results/triangle/psr_map.csv")
    ap.add_argument("--build", action="store_true", help="write --csv from results/triangle/cells_*.csv first")
    ap.add_argument("--metric", choices=("psr", "psr_lo"), default="psr")
    ap.add_argument("--fill", choices=("none", "bcc", "all"), default="bcc")
    ap.add_argument("--agg", choices=("best", "fixed"), default="best", help="with --build: configuration per grid point or one per stack")
    ap.add_argument("--bcc", choices=("profile", "stack"), default="profile", help="with --build: the galvanic profile or the plain HBC stack")
    ap.add_argument("--monotone", action="store_true", help="PSR made non-increasing in every parameter per slice")
    ap.add_argument("--fan", default="", help="datarates (kbit/s) of the networks x sensors slices, e.g. 10,40,160; "
                    "in the hex layouts the 2nd and 3rd are drawn as rotated copies around the slice's outer vertex")
    ap.add_argument("--levels", action="store_true", help="one map per stack with nested PSR levels instead of 0.99 lines")
    ap.add_argument("--lockmap", action="store_true", help="isometric maps of the locked parameter stepped 1u..32u")
    ap.add_argument("--lock-psr", default=str(PSR_TARGET), help="with --lockmap: PSR targets, one hexagon each, e.g. 0.99,0.9")
    ap.add_argument("--cells", choices=("levels", "steps"), help="with --lockmap: grid nodes drawn as split cells")
    ap.add_argument("--honeycomb", action="store_true", help="with --lockmap: first stack in the middle, the others on a spiral")
    ap.add_argument("--hues", action="store_true", help="with --honeycomb and one PSR target: one colour per stack")
    ap.add_argument("--stacks", default="", help="comma-separated stack names to draw (default: all)")
    ap.add_argument("--out", default="results/triangle/psr_map")
    ap.add_argument("--fmt", default="png")
    ap.add_argument("--width-in", type=float, help="with --honeycomb: figure width in inches")
    args = ap.parse_args()
    # relative paths are taken from this script's folder, so the script runs from any working directory
    here = Path(__file__).resolve().parent
    args.csv = str(here / args.csv) if not Path(args.csv).is_absolute() else args.csv
    args.out = str(here / args.out) if not Path(args.out).is_absolute() else args.out
    import pandas as pd
    if args.build:
        build_csv(Path(args.csv).parent, args.csv, args.agg, args.bcc)
    df = pd.read_csv(args.csv)
    stacks = list(dict.fromkeys(df["stack"]))
    if args.stacks:
        stacks = [s for s in args.stacks.split(",") if s in stacks]
    if args.lockmap:
        df = with_pruned(df)
        tag = "_".join(s.split()[0].lower() for s in stacks) if len(stacks) < 3 else "all"
        layout = "hex_center" if args.layout == "all" else args.layout
        out = f"{args.out}_{tag}" + ("" if layout == "hex_center" else f"_{layout}")
        thrs = [float(x) for x in args.lock_psr.split(",")]
        if args.honeycomb:
            draw_honeycomb(stacks, df, args.metric,
                           f"{args.out}_honeycomb{len(stacks)}_" + "_".join(f"{t:g}" for t in thrs)
                           + ("_hues" if args.hues else "") + "." + args.fmt, args.monotone, thrs, args.hues,
                           args.width_in)
            return
        if args.cells:
            draw_cells(stacks[0], df, args.metric, f"{out}_cells_{args.cells}{len(thrs) if args.cells == 'levels' else ''}.png",
                       args.monotone, layout, args.cells, thrs)
            return
        if len(thrs) == 2 and len(stacks) == 1:
            draw_psr_variants(stacks[0], df, args.metric, f"{out}_variants.png", args.monotone, layout, *thrs)
        else:
            draw_lockmap(stacks, df, args.metric, out, args.monotone, thrs, layout)
        return
    if args.levels:
        fan = [float(x) for x in args.fan.split(",")] if args.fan else None
        for layout in (LAYOUTS if args.layout == "all" else [args.layout]):
            draw_levels(layout, stacks, df, args.metric, f"{args.out}_levels_{layout}_{args.metric}.{args.fmt}", args.monotone, fan)
        return
    for layout in (LAYOUTS if args.layout == "all" else [args.layout]):
        fan = [float(x) for x in args.fan.split(",")] if args.fan else None
        tag = f"{layout}_{args.fill}_{args.metric}" + ("_mono" if args.monotone else "") + ("_fan" if fan else "")
        draw(layout, stacks, df, args.metric, args.fill, f"{args.out}_{tag}.{args.fmt}", args.monotone, fan)


if __name__ == "__main__":
    main()

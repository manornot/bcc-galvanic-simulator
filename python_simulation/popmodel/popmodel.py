"""Draw bodies from the population channel model of the galvanic channel at 50 Ohm.

Gain (dB) of link l of body p between sites a and b at the 15 tabulated frequencies (0.6-20 MHz):
    y = M + (b0 + b1 log10 f_MHz)(lc - lc_ref) + u_p + e_l,   lc = log10(c_a c_b / 272)
M: population-mean table (the population mean, valid at the site size lc_ref); c_a, c_b: typical circumferences of the
two sites (cm). The first two terms are the mean layer, the site term carries the dependence on the sites; u_p = W_u z_u + s_u n_u (once per
body) and e_l = W_e z_e + s_e n_e (once per link) are the deviation layer, z ~ N(0, I_K), n ~ N(0, I_15). A reference
mean measured elsewhere at the same site size replaces M through mean_db. Table: params.json, from
the fit of the population model on the cleaned 50 Ohm curves of the open dataset. A draw is static for a body and
session and reciprocal; between frequencies interpolate linearly in log f.
"""
import json
from pathlib import Path
import numpy as np

PARAMS = json.loads((Path(__file__).resolve().parent / "params.json").read_text(encoding="ascii"))
FREQS = np.array(PARAMS["freqs_hz"])
PATHS = PARAMS["selected_paths"]
SITES = sorted(PARAMS["site_circumference_cm"])


def path_mean(path, mean_db=None):
    """Mean layer of a TX-RX path between any two sites of SITES (dB at FREQS); mean_db: reference mean at lc_ref."""
    a, b = path.split("-")
    c, t = PARAMS["site_circumference_cm"], PARAMS["site_term"]
    lc = np.log10(c[a] * c[b] / t["ref_c2_cm2"])
    mean = np.asarray(PARAMS["mean_db"] if mean_db is None else mean_db, float)
    return mean + (t["b0_db"] + t["b1_db"] * np.log10(FREQS / 1e6)) * (lc - t["lc_ref"])


def _term(t, rng):
    W = np.asarray(t["W_db"])
    return W @ rng.standard_normal(W.shape[1]) + np.sqrt(t["nugget_db2"]) * rng.standard_normal(len(W))


def draw_body(rng, paths=PATHS, mean_db=None):
    """One body: {path: gain curve in dB at FREQS} for the given paths (each path one link)."""
    u = _term(PARAMS["body"], rng)
    return {k: path_mean(k, mean_db) + u + _term(PARAMS["link"], rng) for k in paths}


def gain_at(curve, f_hz):
    """Gain of a drawn curve at f_hz (0.6-20 MHz), linear in log f."""
    return np.interp(np.log10(f_hz), np.log10(FREQS), curve)


def covariances():
    """Body and link covariances of the table (15 x 15, dB^2)."""
    return tuple(np.asarray(PARAMS[t]["W_db"]) @ np.asarray(PARAMS[t]["W_db"]).T + PARAMS[t]["nugget_db2"] * np.eye(len(FREQS))
                 for t in ("body", "link"))

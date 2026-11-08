"""Measured body-channel gain curves for the coexistence simulations.

Purpose: turn the cleaned measurement set into aggregate gain-vs-frequency curves that a simulation can
interpolate. The cleaning rules are those of an earlier script (fdma_ceiling):
curves kept by the cleaning-v2 manifest, repaired bins re-interpolated in log f, subject-path means,
subject-balanced grand mean (mean over subjects of each subject's mean over paths), path means over subjects.
Inputs:  PhD-DataProcessing/data/all_measurements.csv from the dataset, and the cleaning manifest (--flags),
         which is not part of this repository; channel_curves.csv is the stored result
Outputs: channel_curves.csv next to this file, columns curve,load,freq_hz,gain_db. curve is 'cohort' (grand
         mean) or a path 'TX-RX' (mean over subjects of the subject-path mean). Means only: no subject-level row
         and no per-subject spread is written.
         Prints the usable band (3 dB rule) of the cohort and RHW1-LHW1 curves.
         Library: read_curve, usable_band, tone_centre.
argv:    [--csv PATH] [--flags PATH] [--out PATH]
"""
from __future__ import annotations
import argparse
from functools import lru_cache
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
NAMES = ["LE1", "RE1", "LE2", "RE2", "F", "NOSE", "NECK", "LSB", "RSB", "LHD", "RHD", "LHF", "LHP", "RHP",
         "LHW1", "RHW1", "LHW2", "RHW2", "LHE", "LHT", "CHEST", "LOIN", "LLQ", "LLB", "LLK", "LLC", "LLAT",
         "RLAT", "LLAJ", "RHT", "RHE"]
FREQS = np.array([50e3, 100e3, 200e3, 400e3, 600e3, 800e3, 1e6, 2e6, 3e6, 4e6, 5e6, 6e6, 7e6, 8e6, 9e6, 10e6,
                  12.5e6, 15e6, 20e6])
LOGF = np.log10(FREQS)
F_LO, F_HI = 0.6e6, 20e6
CURVES_CSV = HERE / "channel_curves.csv"


def load_curves(csv_path, flags_path):
    """{load: (grand_mean[19], subject-path means DataFrame indexed (subject, path))} for loads '50' and '1M'."""
    df = pd.read_csv(csv_path)
    df["path"] = df.tx_point.map(dict(enumerate(NAMES))) + "-" + df.rx_point.map(dict(enumerate(NAMES)))
    df["subject"] = df.subject_id.astype(str)
    flags = pd.read_csv(flags_path)
    out = {}
    for load in ("50", "1M"):
        G = df[[f"rx_gain_{load}_f_{int(f)}" for f in FREQS]].to_numpy(float)
        fl = flags[flags.load == load].reset_index(drop=True)
        assert len(fl) == len(df)
        for i, r in fl.iterrows():
            if isinstance(r.repaired_bins_hz, str) and r.repaired_bins_hz:
                bad = [int(np.where(FREQS == float(h))[0][0]) for h in r.repaired_bins_hz.split(";")]
                good = np.setdiff1d(np.arange(len(FREQS)), bad)
                G[i, bad] = np.interp(LOGF[bad], LOGF[good], G[i, good])
        keep = fl.keep.to_numpy(bool)
        d = df[keep].reset_index(drop=True)
        sp = pd.DataFrame(G[keep]).groupby([d.subject, d.path]).mean()       # subject-path means
        grand = sp.groupby(level=0).mean().mean(axis=0).to_numpy()           # subject-balanced grand mean
        out[load] = (grand, sp)
    return out


def usable_band(freqs, gain_db, tol=3.0):
    """3 dB rule of: bins within tol dB of the maximum over 0.6-20 MHz, edges at geometric
    midpoints of the grid clipped to 0.6 and 20 MHz. Returns (bin freqs, total width Hz, lowest edge, highest edge,
    the kept bins merged into contiguous (lo, hi) intervals)."""
    freqs, gain_db = np.asarray(freqs, float), np.asarray(gain_db, float)
    lo = np.clip(np.r_[freqs[0] ** 2 / freqs[1], np.sqrt(freqs[:-1] * freqs[1:])], F_LO, F_HI)
    hi = np.clip(np.r_[np.sqrt(freqs[:-1] * freqs[1:]), freqs[-1] ** 2 / freqs[-2]], F_LO, F_HI)
    band = (freqs >= F_LO) & (freqs <= F_HI)
    q = np.where(band & (gain_db >= gain_db[band].max() - tol))[0]
    iv = []
    for a, b in zip(lo[q], hi[q]):
        if iv and a <= iv[-1][1] + 1e-6:
            iv[-1][1] = b
        else:
            iv.append([a, b])
    return freqs[q], float(np.sum(hi[q] - lo[q])), float(lo[q[0]]), float(hi[q[-1]]), [(float(a), float(b)) for a, b in iv]


def read_curve(curve, load, path=CURVES_CSV):
    """(freq_hz, gain_db) of one curve of channel_curves.csv."""
    c = pd.read_csv(path, dtype={"load": str})
    c = c[(c.curve == curve) & (c.load == load)].sort_values("freq_hz")
    if c.empty:
        raise KeyError(f"no curve {curve} at load {load} in {path}")
    return c.freq_hz.to_numpy(float), c.gain_db.to_numpy(float)


@lru_cache(maxsize=None)
def tone_centre(rc, curve="cohort", load="50"):
    """Centre (Hz) of one BFSK tone at chip rate rc: the multiple of 1 kHz that maximises the weaker of the
    two tone gains g(f - rc), g(f + rc) of the curve (log-f interpolation; ties: the lowest f), both tones inside the
    measured 0.05 to 20 MHz; None when 2 rc exceeds that span.
    ponytail: one placement per rate, from the cohort curve, serves every link (one BAN tone); a per-path placement
    would lift the weak paths, which every rate the engine offers keeps far above the error knee."""
    fr, gd = read_curve(curve, load)
    f = np.arange(np.ceil((fr[0] + rc) / 1e3), np.floor((fr[-1] - rc) / 1e3) + 1) * 1e3
    if not len(f):
        return None
    g = lambda x: np.interp(np.log10(x), np.log10(fr), gd)          # noqa: E731
    return float(f[np.argmax(np.minimum(g(f - rc), g(f + rc)))])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", default=ROOT / "PhD-DataProcessing/data/all_measurements.csv")
    ap.add_argument("--flags", default=ROOT / "cleaning_manifest.csv")
    ap.add_argument("--out", default=CURVES_CSV)
    a = ap.parse_args()
    rows = []
    for load, (grand, sp) in load_curves(a.csv, a.flags).items():
        rows += [dict(curve="cohort", load=load, freq_hz=int(f), gain_db=g) for f, g in zip(FREQS, grand)]
        for path, grp in sp.groupby(level=1):
            rows += [dict(curve=path, load=load, freq_hz=int(f), gain_db=g) for f, g in zip(FREQS, grp.mean(axis=0))]
    pd.DataFrame(rows).to_csv(a.out, index=False, float_format="%.4f", lineterminator="\n")
    print(f"wrote {a.out} ({len(rows)} rows)")
    for curve, load in (("cohort", "50"), ("cohort", "1M"), ("RHW1-LHW1", "50"), ("RHW1-LHW1", "1M")):
        f, g = read_curve(curve, load, a.out)
        bins, width, lo, hi, _ = usable_band(f, g)
        print(f"{curve:10s} {load:>2s} Ohm: usable band {width / 1e6:.1f} MHz, edges {lo / 1e6:.2f}-{hi / 1e6:.2f} MHz, "
              f"bins {', '.join(f'{b / 1e6:g}' for b in bins)}")


if __name__ == "__main__":
    main()

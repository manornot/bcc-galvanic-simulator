"""Channel of the coexistence simulator: room geometry, 2.4 GHz path loss, noise, BCC and HBC links.

Purpose: PARAMETERS.md sections 5.1 to 5.7 as functions for the SINR reception model (DESIGN.md 5.2).
- make_room / Room: square room of area k a_p (a_p by venue and density), body centres uniform with >= 0.6 m
  spacing (random sequential placement), facing uniform; hub on the chest (node 0) and sensors 1..n_dev on the
  body sites of PARAMETERS 5.3 in the parameters.json order head, belt, left wrist, right wrist, left ankle,
  right ankle, back (repeated beyond 7 sensors). Body frame x forward, y left, z up; room frame z up.
- link_gain_db / gain_matrix_db: static gain (dB, negative) per unordered node pair, drawn once per room (trial)
  and cached in it, so it is reciprocal. Same body: -(median of the site pair (5.1) + N(0, static sigma of the
  case)). Different bodies: -(PL0 + 10 n log10(max(d, 0.3 m) / 1 m) + X + L_sh min(m_blk, 2)) of 5.2 with d the
  3-D node distance, X ~ N(0, sigma of the case), m_blk the torso cylinders the straight node segment crosses
  (own torsos included). fade_db: per-packet fade in dB, zero mean, added to the static gain.
- noise_dbm: -174 dBm/Hz + 10 log10(B) + NF + L_impl (5.4); dbm_to_mw, mw_to_dbm.
- ec_n0_db, bcc_link: galvanic link Ec/N0 = (A |H|)^2 / (2 Rc N0), N0 = 4 k T 50 Ohm x NF at 290 K (link
  budget, an earlier script (link_budget)), |H| of a path mean of channel_curves.csv
  interpolated in log10 f; N0 stays the 50 Ohm value for the 1 MOhm curves too (as in bcc_fdma.py). bfsk_ec_n0_db:
  the one Ec/N0 of a BFSK chip whose two tones (f -/+ R_c) arrive at different Ec/N0.
- hbc_coupling_db: body-to-body loss in excess of the intra-body loss, HBC standard (5.7) or BCC sensitivity case.
Inputs:  parameters.json (sections common, channel, ieee802156_hbc, bcc), channel_curves.csv
Outputs: none (library)
"""
from __future__ import annotations
import json
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
_P = json.loads((HERE / "parameters.json").read_text(encoding="utf-8"))
CH = _P["channel"]
_HBC = _P["ieee802156_hbc"]
CURVES_CSV = HERE / "channel_curves.csv"
SITES = tuple(CH["room"]["nodes_body_frame_m"])                 # 'hub_chest' first, then the sensor sites
_XYZ = np.array(list(CH["room"]["nodes_body_frame_m"].values()), float)
DENSITIES = ("low", "typical", "high")                          # order of area_per_person_m2
TORSO_R = CH["room"]["torso_radius_m"]
TORSO_Z = (0.8, 1.5)    # torso cylinder height range (m), assumption: crotch to shoulder of a 1.75 m adult (A)
D_MIN = CH["inter_body"]["d_min_m"]
N0_THERMAL = 4 * 1.380649e-23 * 290.0 * 50.0                    # 8.0e-19 V^2/Hz, 50 Ohm at 290 K
HBC_INTRA_DB = _HBC["intra_body_loss_dB"]                       # 60 dB, HBC standard intra-body loss

# Site pair -> PARAMETERS 5.1 per-link key. The hub is taken on the right of the chest (the NICTA set of [84] has
# the right wrist and ankle 10-30 dB stronger than the left ones), so right = same side (A). The belt link to the
# chest uses the hip-chest value (no chest-belt value exists) (A). Unlisted pairs use the typical median 68 dB.
_ON_BODY_KEY = {frozenset(k): v for k, v in {
    ("hub_chest", "right_wrist"): "chest-wrist_same", ("hub_chest", "left_wrist"): "chest-wrist_opp",
    ("hub_chest", "right_ankle"): "chest-ankle_same", ("hub_chest", "left_ankle"): "chest-ankle_opp",
    ("hub_chest", "back"): "chest-back", ("hub_chest", "belt"): "hip-chest", ("head", "belt"): "head-belt",
    ("belt", "left_wrist"): "belt-wrist", ("belt", "right_wrist"): "belt-wrist",
    ("head", "left_wrist"): "head-wrist", ("head", "right_wrist"): "head-wrist"}.items()}


@dataclass
class Room:
    """k bodies: centres (k, 2) m, facing (k,) rad (body x axis in the room), n_dev sensors per body.
    Derived: sites (node index -> site name, 0 = 'hub_chest'), pos (k, n_dev + 1, 3) node positions in m.
    gains caches the static gain matrix per case (gain_matrix_db)."""
    centres: np.ndarray
    facing: np.ndarray
    n_dev: int
    side_m: float = float("nan")
    gains: dict = field(default_factory=dict, repr=False)

    def __post_init__(self):
        self.centres = np.asarray(self.centres, float).reshape(-1, 2)
        self.facing = np.asarray(self.facing, float).reshape(-1)
        # ponytail: sites repeat beyond 7 sensors (co-located sensors), add sites if N_dev > 7 matters
        self.sites = ("hub_chest",) + tuple(SITES[1 + i % (len(SITES) - 1)] for i in range(self.n_dev))
        xyz = _XYZ[[SITES.index(s) for s in self.sites]]
        c, s = np.cos(self.facing)[:, None], np.sin(self.facing)[:, None]
        self.pos = np.stack([self.centres[:, :1] + xyz[:, 0] * c - xyz[:, 1] * s,
                             self.centres[:, 1:] + xyz[:, 0] * s + xyz[:, 1] * c,
                             np.broadcast_to(xyz[:, 2], c.shape[:1] + xyz[:, 2].shape)], axis=-1)

    @property
    def k(self):
        return len(self.centres)


def make_room(rng, k, n_dev, venue="ward", density="typical"):
    """Room of one trial (PARAMETERS 5.3). ponytail: static bodies, the 1 s mobility epochs are not modelled."""
    side = float(np.sqrt(k * CH["room"]["area_per_person_m2"][venue][DENSITIES.index(density)]))
    pts = np.empty((0, 2))
    for _ in range(1000 * k):
        if len(pts) == k:
            break
        p = rng.uniform(0.0, side, 2)
        if np.all(np.hypot(*(pts - p).T) >= CH["room"]["min_spacing_m"]):
            pts = np.vstack([pts, p])
    if len(pts) < k:
        raise RuntimeError(f"could not place {k} bodies with {CH['room']['min_spacing_m']} m spacing in {side:.2f} m")
    return Room(pts, rng.uniform(0.0, 2 * np.pi, k), n_dev, side)


def on_body_median_db(site_a, site_b):
    """Median on-body path loss (dB, positive) of a site pair, PARAMETERS 5.1 per-link values."""
    key = _ON_BODY_KEY.get(frozenset((site_a, site_b)))
    return float(CH["on_body"]["per_link_dB"][key] if key else CH["on_body"]["median_dB"]["typical"])


def inter_body_pl_db(d_m, m_blk, case="typical"):
    """PARAMETERS 5.2 path loss without shadowing and fade: PL0 + 10 n log10(max(d, 0.3) / 1 m) + L_sh min(m_blk, 2)."""
    c = CH["inter_body"][case]
    return c["PL0"] + 10 * c["n"] * np.log10(np.maximum(d_m, D_MIN)) + c["L_sh"] * np.minimum(m_blk, 2)


def torso_blocks(room, a, b):
    """Number of torso cylinders (radius TORSO_R, height TORSO_Z, one per body) that each straight segment a[i] ->
    b[i] (points (N, 3) m) passes through; a segment that only touches a cylinder (node on the torso surface
    pointing away) does not count."""
    a, b = np.atleast_2d(a), np.atleast_2d(b)
    d = b - a
    A = d[:, 0] ** 2 + d[:, 1] ** 2
    flat = A < 1e-12                                            # vertical or zero-length segment
    n = np.zeros(len(a), int)
    for cx, cy in room.centres:
        ex, ey = a[:, 0] - cx, a[:, 1] - cy
        B = 2 * (ex * d[:, 0] + ey * d[:, 1])
        C = ex ** 2 + ey ** 2 - TORSO_R ** 2
        with np.errstate(divide="ignore", invalid="ignore"):
            sq = np.sqrt(np.maximum(B * B - 4 * A * C, 0.0))
            t1 = np.where(flat, 0.0, np.clip((-B - sq) / (2 * A), 0, 1))
            t2 = np.where(flat, (C < 0).astype(float), np.clip((-B + sq) / (2 * A), 0, 1))
        z1, z2 = a[:, 2] + t1 * d[:, 2], a[:, 2] + t2 * d[:, 2]
        n += (t2 - t1 > 1e-9) & (np.maximum(z1, z2) > TORSO_Z[0]) & (np.minimum(z1, z2) < TORSO_Z[1])
    return n


def gain_matrix_db(room, rng, case="typical"):
    """Static gain (dB, negative) between every pair of nodes, flat index body * (n_dev + 1) + node, symmetric,
    NaN on the diagonal. Drawn with rng on the first call per case, then returned from room.gains."""
    if case not in room.gains:
        m = room.n_dev + 1
        P = room.pos.reshape(-1, 3)
        body, node = np.repeat(np.arange(room.k), m), np.tile(np.arange(m), room.k)
        i, j = np.triu_indices(len(P), 1)
        same = body[i] == body[j]
        med = np.array([[on_body_median_db(x, y) for y in room.sites] for x in room.sites])
        loss = np.empty(len(i))
        loss[same] = med[node[i[same]], node[j[same]]] + rng.normal(
            0.0, CH["on_body"]["static_sigma_dB"][case], int(same.sum()))
        ii, jj = i[~same], j[~same]
        loss[~same] = inter_body_pl_db(np.linalg.norm(P[ii] - P[jj], axis=1), torso_blocks(room, P[ii], P[jj]),
                                       case) + rng.normal(0.0, CH["inter_body"][case]["sigma"], len(ii))
        G = np.full((len(P), len(P)), np.nan)
        G[i, j] = G[j, i] = -loss
        room.gains[case] = G
    return room.gains[case]


def link_gain_db(room, rng, tx, rx, case="typical"):
    """Static gain (dB, negative) from node tx to node rx, nodes (body_index, node_index), node_index 0 = hub."""
    m = room.n_dev + 1
    return float(gain_matrix_db(room, rng, case)[tx[0] * m + tx[1], rx[0] * m + rx[1]])


def fade_db(rng, size, kind, case="typical"):
    """Per-packet fade (dB, zero mean, add to the gain): kind 'on_body' (sigma 2.8 / 3.5 / 4.3 dB by case) or
    'inter_body' (sigma 5 dB)."""
    sigma = {"on_body": CH["on_body"]["fade_sigma_dB"][case], "inter_body": CH["inter_body"]["fade_sigma_dB"]}[kind]
    return rng.normal(0.0, sigma, size)


def noise_dbm(bw_hz, nf_db, l_impl_db=0.0):
    """Receiver noise power (dBm) in the bandwidth bw_hz, PARAMETERS 5.4."""
    return _P["common"]["thermal_noise_dBm_per_Hz"] + 10 * np.log10(bw_hz) + nf_db + l_impl_db


def dbm_to_mw(x):
    return 10.0 ** (np.asarray(x, float) / 10)


def mw_to_dbm(x):
    return 10 * np.log10(x)


def ec_n0_db(amp_v, gain_db, rc, nf_db):
    """Chip energy over noise density (dB) of a link with drive amplitude amp_v (V), amplitude gain gain_db (20 log10
    |H|), chip rate rc (chip/s): Ec = (A |H|)^2 / (2 Rc), N0 = 4 k T 50 Ohm x 10^(NF / 10)."""
    return 20 * np.log10(amp_v) + gain_db - 10 * np.log10(2 * rc * N0_THERMAL) - nf_db


def bfsk_ec_n0_db(e_db):
    """One Ec/N0 (dB) for a BFSK chip whose two tones arrive at the Ec/N0 of the last axis of e_db: a chip
    on tone i errs with 0.5 exp(-E_i / 2) (non-coherent BFSK), so equiprobable chips err with the mean of the two,
    which 0.5 exp(-E / 2) gives at E = E_min - 2 ln((1 + exp(-(E_max - E_min) / 2)) / 2): E_min to E_min + 2 ln 2.
    ponytail: the two tones as independent channels; the waveform also leaks the strong tone into the weak one's
    chips through the 4 R_c filter (tones tens of dB apart) and, from 7 Mchip/s, a lower tone near 0 Hz fades with
    its own image (up to 0.5 dB at the knee at 7 M); the engine offers only rates where the waveform makes no chip
    error at 1 V (results/bcc_rate_stability.md)."""
    x = 10 ** (np.asarray(e_db, float) / 10)
    lo, hi = x.min(axis=-1), x.max(axis=-1)
    return 10 * np.log10(lo - 2 * np.log((1 + np.exp(-(hi - lo) / 2)) / 2))


@lru_cache(maxsize=None)
def _curve_gains_db(load, f_hz):
    """{curve: gain dB at f_hz} for every curve of channel_curves.csv at one load, interpolated in log10 f."""
    c = pd.read_csv(CURVES_CSV, dtype={"load": str})
    return {name: float(np.interp(np.log10(f_hz), np.log10(g.freq_hz), g.gain_db))
            for name, g in c[c.load == load].sort_values("freq_hz").groupby("curve")}


def bcc_link(rng, n_links, amp_v, rc, nf_db, f_hz, load="50", path=None):
    """Ec/N0 (dB) of n_links galvanic links at the frequency f_hz: path = a curve name of channel_curves.csv
    ('RHW1-LHW1', 'cohort', ...) for every link, or None to draw one of the 36 measured paths per link (uniform).
    f_hz a sequence (the BFSK tones f -/+ R_c of every BCC tone): an (n_links, len(f_hz)) array, all frequencies of a
    link on its path."""
    g = [_curve_gains_db(str(load), float(f)) for f in np.atleast_1d(f_hz)]
    names = sorted(p for p in g[0] if p != "cohort") if path is None else [path]
    pick = rng.integers(len(names), size=n_links) if path is None else np.zeros(n_links, int)
    e = ec_n0_db(amp_v, np.array([[gf[p] for gf in g] for p in names])[pick], rc, nf_db)
    return e if np.ndim(f_hz) else e[:, 0]


def hbc_coupling_db(room, case="typical"):
    """(k, k) body-to-body loss in excess of the intra-body loss (dB, diagonal 0), body centre distance d.
    HBC standard (PARAMETERS 5.7): case 'low' / 'typical' / 'high' = 10 / 15 / 20 dB at 1 m + 20 log10(max(d, 0.3)
    / 1 m); absolute HBC loss = HBC_INTRA_DB + this. BCC sensitivity: case 'bcc40' / 'bcc60' = 40 / 60 dB flat over
    the intra-body loss of the measured link (the 0.3 m clip is an assumption; the generator keeps d >= 0.6 m)."""
    d = np.hypot(*(room.centres[:, None, :] - room.centres[None, :, :]).transpose(2, 0, 1))
    if case.startswith("bcc"):
        if int(case[3:]) not in _P["bcc"]["inter_body_sensitivity_dB"]:
            raise KeyError(case)
        x = np.full(d.shape, float(case[3:]))
    else:
        x = _HBC["inter_body_extra_dB"][case] + _HBC["inter_body_slope_dB_per_decade"] * np.log10(np.maximum(d, D_MIN))
    np.fill_diagonal(x, 0.0)
    return x

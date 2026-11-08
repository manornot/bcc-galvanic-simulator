"""Plain-assert tests of channel.py, no pytest needed.

Purpose: room size and spacing, node placement, inter-body path loss (1 m, 3 m, 0.3 m clip, torso blocking and
its cap), on-body medians per site pair, static shadowing and fade sigmas, reciprocity, noise, the BCC link Ec/N0
of RHW1-LHW1 at 8 MHz, the one Ec/N0 of a BFSK chip from its two tones, and the HBC / BCC inter-body coupling.
Every expected value is derived by hand in the comment next to its assert; statistical checks use tolerances of
about 4 standard errors.
Inputs:  parameters.json, channel_curves.csv (through channel.py).
Outputs: one line per passed test on stdout; an AssertionError on failure.
argv:    none (python python_simulation/coexistence/test_channel.py), runtime a few seconds.
"""
import inspect
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import channel as ch   # noqa: E402

# node indices with n_dev = 7 (parameters.json order): 0 hub_chest, 1 head, 2 belt, 3 left_wrist, 4 right_wrist,
# 5 left_ankle, 6 right_ankle, 7 back
HUB, HEAD, BELT, LW, RW, LA, RA, BACK = range(8)


def test_room_area_and_spacing():
    # A = k a_p (PARAMETERS 5.3): ward typical 16.0 m2 x 5 = 80 m2; gym high 2.3 x 30 = 69 m2; office low 18.6 x 1
    rng = np.random.default_rng(1)
    for k, venue, dens, area in ((5, "ward", "typical", 80.0), (30, "gym", "high", 69.0), (1, "office", "low", 18.6)):
        for _ in range(50):
            r = ch.make_room(rng, k, 2, venue, dens)
            assert abs(r.side_m ** 2 - area) < 1e-9 and r.k == k and r.pos.shape == (k, 3, 3)
            assert np.all((r.centres >= 0) & (r.centres <= r.side_m))
            dist = np.hypot(*(r.centres[:, None] - r.centres[None]).transpose(2, 0, 1))
            assert np.all(dist[np.triu_indices(k, 1)] >= 0.6), dist.min()     # spacing >= 0.6 m


def test_node_positions():
    # body at (2, 3) facing +y (theta = pi/2):
    # x_room = cx + x cos - y sin = cx - y, y_room = cy + x sin + y cos = cy + x
    # hub (0.15, 0, 1.3) -> (2, 3.15, 1.3); left wrist (0.05, 0.30, 0.9) -> (1.70, 3.05, 0.9) (left of +y is -x);
    # back (-0.15, 0, 1.3) -> (2, 2.85, 1.3)
    r = ch.Room([[2.0, 3.0]], [np.pi / 2], 7)
    assert r.sites[HUB] == "hub_chest" and r.sites[RW] == "right_wrist" and r.sites[BACK] == "back"
    assert np.allclose(r.pos[0, HUB], [2.0, 3.15, 1.3])
    assert np.allclose(r.pos[0, LW], [1.70, 3.05, 0.9])
    assert np.allclose(r.pos[0, BACK], [2.0, 2.85, 1.3])
    # sites cycle beyond 7 sensors: sensor 8 is again the head
    assert ch.Room([[0, 0]], [0], 8).sites[8] == "head"


def test_inter_body_pl():
    # typical PL0 48, n 2, L_sh 15: 1 m -> 48; 3 m -> 48 + 20 log10 3 = 48 + 9.5424 = 57.5424
    assert abs(ch.inter_body_pl_db(1.0, 0) - 48.0) < 1e-9
    assert abs(ch.inter_body_pl_db(3.0, 0) - 57.5424) < 1e-4
    # 0.3 m clip: 0.1 m and 0 m give 48 + 20 log10 0.3 = 48 - 10.4576 = 37.5424
    assert abs(ch.inter_body_pl_db(0.1, 0) - 37.5424) < 1e-4 and abs(ch.inter_body_pl_db(0.0, 0) - 37.5424) < 1e-4
    # blocking: 1 torso +15 -> 63; 3 torsos capped at 2 -> 48 + 30 = 78
    assert ch.inter_body_pl_db(1.0, 1) == 63.0 and ch.inter_body_pl_db(1.0, 3) == 78.0
    # low case 1 m one torso: 42 + 6 = 48; high case 3 m: 55 + 30 log10 3 = 55 + 14.3136 = 69.3136
    assert ch.inter_body_pl_db(1.0, 1, "low") == 48.0 and abs(ch.inter_body_pl_db(3.0, 0, "high") - 69.3136) < 1e-4


def _stats(make_room, tx, rx, case="typical", n=4000, seed=7):
    rng = np.random.default_rng(seed)
    g = np.array([ch.link_gain_db(make_room(), rng, tx, rx, case) for _ in range(n)])
    return g.mean(), g.std()


def test_inter_body_gain_geometry():
    # A at (0, 0) facing +x, B at (1.3, 0) facing -x: chests at x 0.15 and 1.15, distance 1.0 m, the segment leaves
    # A's torso surface outwards and meets B's surface from outside -> no blocking; gain = -(48 + X), X ~ N(0, 5).
    # 4000 draws: standard error of the mean 5 / sqrt(4000) = 0.079 dB, of the std about 5 / sqrt(8000) = 0.056 dB
    m, s = _stats(lambda: ch.Room([[0, 0], [1.3, 0]], [0, np.pi], 1), (0, HUB), (1, HUB))
    assert abs(m + 48.0) < 0.3 and abs(s - 5.0) < 0.25, (m, s)
    # B at (3.3, 0): 3.0 m -> -57.5424
    m, s = _stats(lambda: ch.Room([[0, 0], [3.3, 0]], [0, np.pi], 1), (0, HUB), (1, HUB))
    assert abs(m + 57.5424) < 0.3, m
    # side by side 0.6 m apart, both facing +x: A's left wrist (0.05, 0.30, 0.9) and B's right wrist
    # (0.05, 0.6 - 0.30, 0.9) coincide, d = 0 -> clipped to 0.3 m: -37.5424, no torso on a zero-length segment
    m, s = _stats(lambda: ch.Room([[0, 0], [0, 0.6]], [0, 0], 7), (0, LW), (1, RW))
    assert abs(m + 37.5424) < 0.3, m
    # A's back to B's chest (1.3 m apart as above): back (-0.15, 0) -> (1.15, 0), 1.3 m, crosses A's torso:
    # 48 + 20 log10 1.3 + 15 = 48 + 2.2789 + 15 = 65.2789
    m, s = _stats(lambda: ch.Room([[0, 0], [1.3, 0]], [0, np.pi], 7), (0, BACK), (1, HUB))
    assert abs(m + 65.2789) < 0.3, m


def test_torso_blocks():
    # line of bodies on y = 0: A (0, 0) facing +x, C (1.5, 0) facing +x, B (3, 0) facing -x, all 7 sensors
    r = ch.Room([[0, 0], [1.5, 0], [3, 0]], [0, 0, np.pi], 7)
    A, C, B = 0, 1, 2
    P = r.pos
    cases = [
        (P[A, HUB], P[C, BACK], 0),   # x 0.15 -> 1.35: leaves A's surface outwards, stops on C's surface
        (P[A, HUB], P[C, HUB], 1),    # x 0.15 -> 1.65: through C's torso
        (P[A, BACK], P[C, BACK], 1),  # x -0.15 -> 1.35: through A's own torso
        (P[A, BACK], P[B, BACK], 3),  # x -0.15 -> 3.15: through A, C and B (the PL caps this at 2)
        (P[A, HEAD], P[B, HEAD], 0),  # z = 1.65 above the torso top 1.5 m: passes over C, heads inside A/B discs
        (P[A, RA], P[B, LA], 0),      # ankles at y = -0.12, z = 0.1 below the torso bottom 0.8 m
    ]
    got = ch.torso_blocks(r, np.array([c[0] for c in cases]), np.array([c[1] for c in cases]))
    assert list(got) == [c[2] for c in cases], got
    # lateral miss and hit: C moved to y = 0.16 (> 0.15 radius) or y = 0.14 (< 0.15) beside the A-B chest line
    for y, want in ((0.16, 0), (0.14, 1)):
        r = ch.Room([[0, 0], [1.5, y], [3, 0]], [0, 0, np.pi], 1)
        assert ch.torso_blocks(r, r.pos[0, HUB], r.pos[2, HUB])[0] == want, y


def test_on_body_medians():
    # PARAMETERS 5.1 per-link values: right side = hub side (same), belt-chest read as hip-chest, unlisted -> 68
    want = {("hub_chest", "right_wrist"): 64, ("hub_chest", "left_wrist"): 70, ("hub_chest", "right_ankle"): 75,
            ("hub_chest", "left_ankle"): 82, ("hub_chest", "back"): 78, ("hub_chest", "belt"): 62.5,
            ("head", "belt"): 54, ("belt", "right_wrist"): 56, ("left_wrist", "belt"): 56, ("head", "left_wrist"): 60,
            ("hub_chest", "head"): 68, ("left_wrist", "right_wrist"): 68, ("head", "head"): 68}
    for (a, b), v in want.items():
        assert ch.on_body_median_db(a, b) == v and ch.on_body_median_db(b, a) == v, (a, b)
    # drawn gain hub -> right wrist: -(64 + N(0, 4.5)); SE of the mean 4.5 / sqrt(4000) = 0.071 dB
    m, s = _stats(lambda: ch.Room([[0, 0]], [0], 7), (0, HUB), (0, RW))
    assert abs(m + 64.0) < 0.3 and abs(s - 4.5) < 0.25, (m, s)
    # high case: static sigma 6.9 dB, same median
    m, s = _stats(lambda: ch.Room([[0, 0]], [0], 7), (0, BACK), (0, HUB), case="high")
    assert abs(m + 78.0) < 0.45 and abs(s - 6.9) < 0.35, (m, s)


def test_reciprocity_and_cache():
    # one draw per unordered pair and trial: both directions and repeated calls give the same value
    rng = np.random.default_rng(3)
    r = ch.make_room(rng, 6, 4)
    g1 = ch.link_gain_db(r, rng, (2, 0), (4, 3))
    assert g1 == ch.link_gain_db(r, rng, (4, 3), (2, 0)) == ch.link_gain_db(r, rng, (2, 0), (4, 3))
    G = ch.gain_matrix_db(r, rng)
    assert G.shape == (30, 30) and np.all(np.isnan(np.diag(G))) and np.allclose(G, G.T, equal_nan=True)
    assert G[2 * 5 + 0, 4 * 5 + 3] == g1                        # flat index body * (n_dev + 1) + node


def test_fade():
    # sigmas of PARAMETERS 5.1 (on-body typical 3.5, high 4.3) and 5.2 (5 dB); 200 000 draws: SE of std ~ 0.16 %
    rng = np.random.default_rng(5)
    for args, sigma in ((("on_body",), 3.5), (("on_body", "high"), 4.3), (("inter_body",), 5.0)):
        x = ch.fade_db(rng, 200_000, *args)
        assert abs(x.std() - sigma) < 0.02 and abs(x.mean()) < 0.05, (args, x.std())


def test_noise():
    # -174 + 10 log10(1e6) + NF = -174 + 60 + 6 = -108 dBm; with NF 10 and L_impl 5: -174 + 60 + 15 = -99 dBm
    assert abs(ch.noise_dbm(1e6, 6.0) + 108.0) < 1e-9 and abs(ch.noise_dbm(1e6, 10.0, 5.0) + 99.0) < 1e-9
    assert abs(ch.dbm_to_mw(-30.0) - 1e-3) < 1e-15 and abs(ch.mw_to_dbm(1.0)) < 1e-12


def test_bcc_link():
    # RHW1-LHW1, 50 Ohm, 8 MHz: channel_curves.csv row gain_db = -53.5954 (= 20 log10 |H|).
    # N0 thermal = 4 x 1.380649e-23 x 290 x 50 = 8.00776e-19 V^2/Hz -> 10 log10 = -180.9649 dB; NF 20 -> -160.9649;
    # 2 Rc = 1e6 Hz (+60 dB) -> 2 Rc N0 = -100.9649 dB; (1 V x |H|)^2 -> -53.5954 dB;
    # Ec/N0 = -53.5954 + 100.9649 = 47.3695 dB
    rng = np.random.default_rng(0)
    e = ch.bcc_link(rng, 3, 1.0, 500e3, 20.0, 8e6, "50", "RHW1-LHW1")
    assert e.shape == (3,) and np.all(np.abs(e - 47.3695) < 1e-3), e
    # path drawn per link: the csv holds 6 TX sites (CHEST, LOIN, RHD, RHE, RHT, RHW1) x 6 RX sites (LHW1, LLAJ,
    # LLB, LLC, LLK, LLQ) = 36 paths besides 'cohort'; 5000 uniform draws hit each (P(miss one) = 36 (35/36)^5000 ~ 0)
    e = ch.bcc_link(rng, 5000, 1.0, 500e3, 20.0, 8e6)
    assert len(np.unique(np.round(e, 9))) == 36 and np.any(np.abs(e - 47.3695) < 1e-3)
    # cohort at 8 MHz: -48.2794 dB -> 47.3695 + (53.5954 - 48.2794) = 52.6855 dB
    assert abs(ch.bcc_link(rng, 1, 1.0, 500e3, 20.0, 8e6, path="cohort")[0] - 52.6855) < 1e-3
    # several tones (multi-tone BCC): one column per tone. RHW1-LHW1 at 12.5 MHz: csv -48.1800 dB -> 100.9649 - 48.1800
    # = 52.7849 dB next to the 47.3695 dB at 8 MHz. With a drawn path every tone of a link is on the path the
    # single-tone call draws from the same stream (one draw per link), so column 0 repeats it.
    e = ch.bcc_link(rng, 2, 1.0, 500e3, 20.0, [8e6, 12.5e6], "50", "RHW1-LHW1")
    assert e.shape == (2, 2) and np.allclose(e, [[47.3695, 52.7849]] * 2, atol=1e-3), e
    a, b = (ch.bcc_link(np.random.default_rng(1), 50, 1.0, 500e3, 20.0, f_hz=f) for f in ([8e6, 12.5e6], 8e6))
    assert np.array_equal(a[:, 0], b)
    # no default frequency: the old 8 MHz default was the one-tone centre removed (it lies in the dip)
    assert inspect.signature(ch.bcc_link).parameters["f_hz"].default is inspect.Parameter.empty


def test_bfsk_ec_n0():
    # One Ec/N0 for the two BFSK tones: the chip error is the mean of 0.5 exp(-E_i / 2) over the two tones.
    # Equal tones: E itself. 0 and 3 dB (1 and 1.995262 linear): 0.5 exp(-0.5) = 0.303265 and 0.5 exp(-0.997631) =
    # 0.184376, mean 0.243821 = 0.5 exp(-E / 2) -> E = -2 ln(0.487641) = 1.43635 = 1.5726 dB. 10 dB against 100 dB (the
    # strong tone never errs): mean 0.25 exp(-5) = 0.5 exp(-E / 2) -> E = 10 + 2 ln 2 = 11.38629 = 10.5638 dB
    got = ch.bfsk_ec_n0_db([[20.0, 20.0], [0.0, 3.0], [3.0, 0.0], [10.0, 100.0]])
    assert got.shape == (4,) and np.allclose(got, [20.0, 1.5726, 1.5726, 10.5638], atol=1e-4), got


def test_hbc_coupling():
    # centres (0, 0), (1, 0), (10, 0): typical 15 + 20 log10 d: d 1 -> 15; d 10 -> 35; d 9 -> 15 + 19.0849 = 34.0849
    r = ch.Room([[0, 0], [1, 0], [10, 0]], [0, 0, 0], 1)
    x = ch.hbc_coupling_db(r)
    assert np.allclose(np.diag(x), 0) and np.allclose(x, x.T)
    assert abs(x[0, 1] - 15) < 1e-9 and abs(x[0, 2] - 35) < 1e-9 and abs(x[1, 2] - 34.0849) < 1e-4
    # low / high at 1 m: 10 / 20 dB; 0.6 m typical: 15 + 20 log10 0.6 = 15 - 4.4370 = 10.5630
    assert ch.hbc_coupling_db(r, "low")[0, 1] == 10 and ch.hbc_coupling_db(r, "high")[0, 1] == 20
    r6 = ch.Room([[0, 0], [0.6, 0]], [0, 0], 1)
    assert abs(ch.hbc_coupling_db(r6)[0, 1] - 10.5630) < 1e-4
    # BCC sensitivity cases: flat +40 / +60 dB off the diagonal
    assert np.array_equal(ch.hbc_coupling_db(r, "bcc40"), 40 * (1 - np.eye(3)))
    assert np.array_equal(ch.hbc_coupling_db(r, "bcc60"), 60 * (1 - np.eye(3)))
    assert ch.HBC_INTRA_DB == 60


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print(f"ok {name}")
    print("all channel tests passed")

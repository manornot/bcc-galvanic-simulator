"""Checks of the model checker (python -m pytest test_formal.py, or python test_formal.py): textbook protocols,
superframes traced by hand, and mutations whose verdict is known in advance."""
import unittest

from mc import check, init, step
from faults import check_L
from textbook import run


def bad(r):
    return sorted(e for e in r["events"] if e.isupper() and e != "DONE")


def age_cfg(A, M, **kw):
    return dict(dict(M=M, Q=A + 1, R=2, P=1, B=1, A=A, CAP=A + 2, N=1, K=None, var="age"), **kw)


class Textbook(unittest.TestCase):
    def test_selective_repeat_needs_2W(self):
        # alternating bit (W = 1, M = 2) and selective repeat: safe exactly when M >= 2W (Tanenbaum and Wetherall, selective repeat)
        for W in (1, 2, 3):
            ok, short = run(W, 2 * W), run(W, 2 * W - 1)
            self.assertEqual(bad(ok), [])
            self.assertIn("DONE", ok["events"])
            self.assertIn("ALIAS", short["events"])


class HandTrace(unittest.TestCase):
    cfg = age_cfg(2, 5, Q=2)

    def test_delivered_frame_and_base(self):
        ns, ev, used = step(init(self.cfg), [False, False], self.cfg)
        self.assertEqual(ns, (1, ((1, ()),), (((0, 1, True, 0),),), ()))
        self.assertEqual((ev, used), ([], 0))

    def test_lost_frame_expires_at_age_A(self):
        s, ev, used = step(init(self.cfg), [False, True, False, True], self.cfg)   # frame and its retry lost
        self.assertEqual(s, (1, ((0, ()),), (((0, 2, False, 0),),), ((0, 0, 1),)))
        s, ev, used = step(s, [True, True, True], self.cfg)                         # beacon, new frame, poll lost
        self.assertEqual(s, (2, ((0, ()),), (((0, 1, False, 0), (1, 2, False, 0)),), ((0, 0, 2), (0, 1, 0))))
        s, ev, used = step(s, [False, True], self.cfg)                              # poll lost; queue full
        self.assertEqual(ev, ["drop_queue"])
        self.assertEqual(s[1], ((1, ()),))           # the hub's base passes the expired frame 0
        self.assertEqual(s[3], ((0, 1, 0), (0, 2, 0)))
        s, ev, used = step(s, [False] * 4, self.cfg)
        self.assertIn("drop_retry", ev)               # the sensor gives frame 0 up at the same age


class Mutations(unittest.TestCase):
    def test_age_rule_needs_2A_plus_1(self):
        for A in (2, 3):
            self.assertEqual(bad(check(age_cfg(A, 2 * A + 1))), [])
            self.assertIn("ALIAS_LOSS", bad(check(age_cfg(A, 2 * A))))

    def test_without_hub_expiry_the_base_stalls(self):
        self.assertNotEqual(bad(check(age_cfg(2, 5, hub_expiry=False))), [])

    def test_bitmap_is_not_needed_for_soundness(self):
        self.assertEqual(bad(check(age_cfg(2, 5, B=0))), [])

    def test_retry_limit_by_attempts_acknowledges_a_lost_frame(self):
        cfg = dict(M=8, Q=4, R=2, P=1, B=1, A=4, CAP=12, N=1, var="count")
        self.assertIn("FALSE_ACK", check_L(cfg, 3)[1])

    def test_queue_of_A_frames_refuses_a_frame_within_the_deadline(self):
        # the age rule keeps a frame A superframes plus the one it was sent in: a queue of A frames refuses a new one
        for A in (2, 3):
            self.assertIn("drop_queue", check_L(age_cfg(A, 2 * A + 1, Q=A), A)[1])
            self.assertNotIn("drop_queue", check_L(age_cfg(A, 2 * A + 1), A)[1])

    def test_two_frames_per_superframe_need_2n_A_plus_1_minus_1(self):
        n, A = 2, 2
        self.assertEqual(bad(check(age_cfg(A, 2 * n * (A + 1) - 1, n=n, Q=n * (A + 1), B=n))), [])
        self.assertIn("ALIAS_LOSS", bad(check(age_cfg(A, 2 * n * (A + 1) - 2, n=n, Q=n * (A + 1), B=n))))


if __name__ == "__main__":
    unittest.main()

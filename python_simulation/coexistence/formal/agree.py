"""Python state counts of the two TLC runs that sweep_*.py and faults.py do not cover: the bitmap-free mutation
(tlc/mut_B0_A2_M5) and two sensors under any losses (tlc/sound_N2_A2_M5). Same constants as the TLC runs."""
import time
from mc import check

for name, cfg in (("mut_B0_A2_M5", dict(M=5, Q=3, R=2, P=1, B=0, A=2, CAP=4, N=1, K=None, var="age")),
                  ("sound_N2_A2_M5", dict(M=5, Q=3, R=2, P=1, B=1, A=2, CAP=4, N=2, K=None, var="age"))):
    t0 = time.time()
    r = check(cfg, limit=1_500_000)
    bad = sorted(e for e in r["events"] if e.isupper())
    print(name, "states", r["states"], "trunc", r["truncated"], "bad", bad, "%.1fs" % (time.time() - t0), flush=True)

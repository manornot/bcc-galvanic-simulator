"""Design alternative: the hub takes every decoded number as a new frame (no half window behind the base; under the
age rule no frame reaches the hub behind its base). Smallest sound M, N = 1, any loss pattern. Also the retry limit
counted in attempts ('count', M = 8, a queue of four frames as in the engine, OI-18): the smallest loss budget that
reaches FALSE_ACK, for R = 2 and for R = 3 retries (the engine's HBC retry limit, protocols.py, parameters.json)."""
import time
from mc import check
from faults import check_L

for n, A in ((1, 2), (1, 3), (1, 4), (2, 2), (2, 3)):
    for M in range(n * (A + 1), 2 * n * (A + 1)):
        cfg = dict(M=M, Q=n * (A + 1), R=2, P=1, B=n, A=A, CAP=A + 2, N=1, K=None, var="age", n=n, hub_window=False)
        t0 = time.time()
        r = check(cfg, limit=1_500_000)
        bad = sorted(e for e in r["events"] if e.isupper())
        print("hub_window off n", n, "A", A, "M", M, "2nA+1", 2 * n * A + 1, "states", r["states"], "trunc",
              r["truncated"], "bad", bad, "%.1fs" % (time.time() - t0), flush=True)
        if not bad and not r["truncated"]:
            break
for R in (2, 3):
    for L in (0, 1, 2, 3, 4):
        n_, found, tr = check_L(dict(M=8, Q=4, R=R, P=1, B=1, A=4, CAP=12, N=1, var="count"), L)
        print("count M 8 R", R, "budget", L, "states", n_, "trunc", tr, "FALSE_ACK at depth", found.get("FALSE_ACK"),
              "ALIAS_LOSS at depth", found.get("ALIAS_LOSS"), flush=True)

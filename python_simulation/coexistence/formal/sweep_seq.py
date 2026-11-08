"""Smallest sequence space M without FALSE_ACK / ALIAS_LOSS (age variant, any loss pattern, N = 1), and whether a
B-bit bitmap acknowledges every delivered frame before its deadline (no drop_delivered)."""
import time
from mc import check

for A in (2, 3, 4):
    for B in (1, A):
        for M in range(3, 2 * A + 5):
            cfg = dict(M=M, Q=A + 1, R=2, P=1, B=B, A=A, CAP=A + 2, N=1, K=None, var="age")
            t0 = time.time()
            r = check(cfg, limit=300_000)
            bad = sorted(e for e in r["events"] if e.isupper())
            print("A", A, "B", B, "M", M, "states", r["states"], "trunc", r["truncated"], "bad", bad,
                  "drop_delivered", "drop_delivered" in r["events"], "%.1fs" % (time.time() - t0), flush=True)
            if not bad and not r["truncated"]:
                break

"""Smallest sequence space M without FALSE_ACK / ALIAS_LOSS for n frames per sensor and superframe (age variant, any
loss pattern, N = 1, a bitmap of B = n or 1 bits). With n = 1 the bound is 2A + 1 (sweep_seq.py); does 2nA + 1 hold
for n > 1, or 2n(A + 1) - 1 (the frames alive at a delivery span n(A + 1) numbers: A superframes and the current one)?"""
import time
from mc import check

for n, A in ((1, 2), (1, 3), (2, 2), (3, 2), (2, 3)):
    for B in (n, 1) if n > 1 else (1,):
        for M in range(2 * n * A, 2 * n * (A + 1) + 1):
            cfg = dict(M=M, Q=n * (A + 1), R=2, P=1, B=B, A=A, CAP=A + 2, N=1, K=None, var="age", n=n)
            t0 = time.time()
            r = check(cfg, limit=1_500_000)
            bad = sorted(e for e in r["events"] if e.isupper())
            print("n", n, "A", A, "B", B, "M", M, "2nA+1", 2 * n * A + 1, "2n(A+1)-1", 2 * n * (A + 1) - 1, "states",
                  r["states"], "trunc", r["truncated"], "bad", bad, "%.1fs" % (time.time() - t0), flush=True)
            if not bad and not r["truncated"]:
                break

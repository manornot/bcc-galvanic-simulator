"""State counts for N = 1..10 (age variant, P = 1): (A 3, M 7, L 3) and (A 4, M 9, L 4); stop a row past 150 s."""
import time
from faults import check_L

for A, L in ((3, 3), (4, 4)):
    for N in range(1, 11):
        cfg = dict(M=2 * A + 1, Q=A + 1, R=2, P=1, B=1, A=A, CAP=A + 2, N=N, var="age")
        t0 = time.time()
        n, found, tr = check_L(cfg, L, limit=1_500_000)
        print("A", A, "L", L, "N", N, "states", n, "truncated", tr, "miss", "drop_retry" in found or "drop_queue" in found,
              "bad", sorted(e for e in found if e.isupper()), "%.1fs" % (time.time() - t0), flush=True)
        if time.time() - t0 > 150 or tr:
            break

"""k-fault tolerance of the deadline (age variant, M = 2A + 1): the largest number L of arbitrary losses (beacon copies,
data frames, polls, anywhere in the run) after which no frame misses its deadline A (no undelivered 'drop_retry' and no
frame refused by a full queue, 'drop_queue'; the sensor buffers the n(A + 1) frames of the last A + 1 superframes),
for N sensors and a pool of P exchanges. Also the profile as written ('count') for comparison: FALSE_ACK reachable?"""
import time
from collections import deque
from mc import step, Need


def check_L(cfg, L, limit=300_000):
    N = cfg["N"]
    s0 = ((0, tuple((0, ()) for _ in range(N)), tuple(() for _ in range(N)), ()), L)
    seen, dq, found = {s0}, deque([(s0, 0)]), {}
    while dq:
        (s, rem), d = dq.popleft()
        stack = [[]]
        while stack:
            ch = stack.pop()
            try:
                ns, ev, u = step(s, ch, cfg)
            except Need:
                stack.append(ch + [False])
                if sum(ch) < rem:
                    stack.append(ch + [True])
                continue
            for e in ev:
                found.setdefault(e, d + 1)
            x = (ns, rem - u)
            if x not in seen:
                seen.add(x)
                dq.append((x, d + 1))
                if len(seen) > limit:
                    return len(seen), found, True
    return len(seen), found, False


if __name__ == "__main__":
    for N in (1, 2, 3):
        for P in (1, 2):
            for A in (2, 3, 4):
                row = []
                for L in range(0, 6):
                    cfg = dict(M=2 * A + 1, Q=A + 1, R=2, P=P, B=1, A=A, CAP=A + 2, N=N, var="age")
                    t0 = time.time()
                    n, found, tr = check_L(cfg, L)
                    bad = sorted(e for e in found if e.isupper())
                    miss = "drop_retry" in found or "drop_queue" in found
                    row.append((L, n, tr, miss, bad, round(time.time() - t0, 1)))
                    if miss or tr:
                        break
                print("N", N, "P", P, "A", A, "max L without deadline miss:",
                      max([r[0] for r in row if not r[3] and not r[2]], default=None), row, flush=True)

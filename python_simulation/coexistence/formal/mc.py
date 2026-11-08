"""Explicit-state model checker of the galvanic profile's medium access (DESIGN 4.3, L2-L4;).

Superframe-granular model, N sensors, n periodic frames per sensor and superframe (cfg 'n', default 1), frame j of
superframe c carrying sequence number (n c + j) mod M. Distributed: the sensor acts only on the beacon ack field it
heard (base seq + B-bit bitmap); the hub only on frames it decoded and on its own poll count. Ghost variables (true
delivery, frame age) check the properties. Losses (beacon copies, data frames, polls) are nondeterministic: at most K
per superframe (K = None: unbounded) here, a total budget in faults.py.
Variants: 'count' = the profile as first written (retry limit by attempts, each side counts its own);
          'age'   = repair: both sides give a frame up at age A superframes (the beacon numbers superframes).
Events in capitals are property violations: FALSE_ACK (a sensor drops as acknowledged a frame the hub never decoded),
ALIAS_LOSS (the hub takes a new frame for an old one), STUCK (a frame outlives CAP superframes).
See test_formal.py.
"""
from collections import deque


class Need(Exception):
    pass


def step(s, choices, cfg):
    """One superframe from state s; choices: the loss outcomes consumed in order (True = lost). Returns (s', events,
    losses used)."""
    N, M, Q, R, P, B, A, var, CAP = (cfg[k] for k in "N M Q R P B A var CAP".split())
    n, W = cfg.get("n", 1), cfg.get("W", M // 2)
    it = iter(choices)
    used = [0]

    def lost():
        try:
            c = next(it)
        except StopIteration:
            raise Need
        used[0] += c
        return c

    c, hub, sens, pool = s            # c: superframe count mod M; hub[i] = (base, held) held: ((seq, ghost_age),..)
    hub = [[b, dict(h)] for b, h in hub]
    sens = [list(map(list, q)) for q in sens]   # queue entries [age, tries, delivered, j]
    pool = [list(e) for e in pool]              # [sensor, seq, hub_polls]
    ev = []
    seq_of = lambda f: (n * (c - 1 - f[0]) + f[3]) % M   # frame j generated `age` superframes ago
    # 1 beacon: base + bitmap per sensor; each sensor may miss it
    for i in range(N):
        base, held = hub[i]
        heard = not lost()
        for f in list(sens[i]):
            sq = seq_of(f)
            d = (base - sq) % M
            acked = heard and ((0 < d <= W) or (1 <= (sq - base) % M <= B and sq in held))
            if var == "age" and f[0] >= A:          # given up first: past its age the hub has passed it too
                sens[i].remove(f)
                ev.append("drop_retry" if not f[2] else "drop_delivered")
            elif acked:
                if not f[2]:
                    ev.append("FALSE_ACK")
                sens[i].remove(f)
            elif var == "count" and f[1] > R:
                sens[i].remove(f)
                ev.append("drop_retry" if not f[2] else "drop_delivered")
    # 2 generation (seq n c + j), queue cap Q frames
    for i in range(N):
        for j in range(n):
            if len(sens[i]) < Q:
                sens[i].append([-1, 0, False, j])     # age -1 now, 0 after this superframe's aging
            else:
                ev.append("drop_queue")
    c1 = (c + 1) % M
    new_seqs = [(n * c + j) % M for j in range(n)]

    def deliver(i, f):
        sq = seq_of(f)
        base, held = hub[i]
        if sq in held and held[sq] != f[0]:
            ev.append("ALIAS_LOSS")              # hub takes a different frame with this seq as a duplicate
        elif cfg.get("hub_window", True) and (base - sq) % M and (base - sq) % M <= W and not f[2]:
            ev.append("ALIAS_LOSS")              # seq looks older than base: discarded as a duplicate
        else:
            held[sq] = f[0]
            f[2] = True

    # 3 allocation: one exchange per new frame, in order; the hub lists every expected frame it did not decode
    for i in range(N):
        for f in [x for x in sens[i] if x[1] == 0]:
            f[1] += 1
            if not lost():
                deliver(i, f)
        for sq in new_seqs:
            if sq not in hub[i][1] and not any(e[0] == i and e[1] == sq for e in pool):
                pool.append([i, sq, 0])          # hub: expected frame not decoded -> retry list
    # 4 pool: first P entries, at most once each
    for e in pool[:P]:
        i, sq, _ = e
        e[2] += 1
        if lost():
            continue                              # sensor missed the poll: no attempt (hub cannot tell)
        f = next((x for x in sens[i] if seq_of(x) == sq), None)
        if f is None:
            ev.append("poll_no_frame")
            continue
        f[1] += 1
        if not lost():
            deliver(i, f)
    # hub list clean-up: delivered, or given up (own poll count / age at the next beacon)
    keep = []
    for e in pool:
        i, sq, k = e
        age_next = ((n * c1 - 1 - sq) % M) // n
        if sq in hub[i][1]:
            continue
        if (var == "count" and k > R) or (var == "age" and cfg.get("hub_expiry", True) and age_next >= A):
            continue
        keep.append(e)
    pool = keep
    # 5 aging; hub base advances over held seqs (and, 'age', over expired ones, never past the next new seq)
    nxt = (n * c1) % M
    for i in range(N):
        for f in sens[i]:
            f[0] += 1
            if f[0] >= CAP:
                ev.append("STUCK")
                f[0] = CAP
        base, held = hub[i]
        held = {k: v + 1 for k, v in held.items()}
        while True:
            if base in held:
                del held[base]
            elif var == "age" and cfg.get("hub_expiry", True) and base != nxt and ((nxt - 1 - base) % M) // n >= A:
                pass
            else:
                break
            base = (base + 1) % M
        hub[i] = [base, held]
    for i in range(N):   # hub forgets ghost records of very old frames (keeps the state finite)
        hub[i][1] = {k: min(v, CAP) for k, v in hub[i][1].items()}
    ns = (c1, tuple((b, tuple(sorted(h.items()))) for b, h in hub),
          tuple(tuple(sorted(tuple(x) for x in q)) for q in sens), tuple(tuple(e) for e in pool))
    return ns, ev, used[0]


def init(cfg):
    N = cfg["N"]
    return (0, tuple((0, ()) for _ in range(N)), tuple(() for _ in range(N)), ())


def successors(s, cfg):
    K = cfg["K"]
    stack, out = [[]], []
    while stack:
        ch = stack.pop()
        try:
            ns, ev, u = step(s, ch, cfg)
        except Need:
            stack.append(ch + [False])
            if K is None or sum(ch) < K:
                stack.append(ch + [True])
            continue
        out.append((ns, ev))
    return out


def check(cfg, limit=400_000, succ=successors, s0=None):
    """Breadth-first search of every state reachable from s0 (default: the model's initial state) by succ; records
    the depth at which each event is first seen."""
    s0 = init(cfg) if s0 is None else s0
    seen, dq, found = {s0: None}, deque([(s0, 0)]), {}
    trans = 0
    while dq:
        s, d = dq.popleft()
        for ns, ev in succ(s, cfg):
            trans += 1
            for e in ev:
                if e not in found:
                    found[e] = d + 1
            if ns not in seen:
                seen[ns] = s
                dq.append((ns, d + 1))
                if len(seen) > limit:
                    return dict(states=len(seen), trans=trans, events=found, truncated=True)
    return dict(states=len(seen), trans=trans, events=found, truncated=False)

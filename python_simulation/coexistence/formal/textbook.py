"""Validation of the checker of mc.py on textbook protocols: selective repeat with window W and sequence numbers mod M
over two lossy FIFO channels (Tanenbaum and Wetherall, Computer Networks, 5th ed., the selective-repeat protocol: safe exactly when M >= 2W); the alternating-bit
protocol is the case W = 1, M = 2. The same breadth-first search (mc.check) runs a different successor function.

State: (sbase, snext, sacked, chan_d, chan_a, rbase, rbuf), absolute frame indices < T as ghosts. The sender may
(re)send any unacknowledged frame of its window at any time (every timeout pattern); either channel may lose any
message; each channel holds at most C messages. Events: ALIAS (the receiver takes a frame for another one), ACK_ALIAS
(the sender takes an acknowledgement for another frame), DONE (all T frames acknowledged).
"""
import sys
from mc import check


def succ(s, cfg):
    W, M, T, C = cfg["W"], cfg["M"], cfg["T"], cfg["C"]
    sbase, snext, sacked, cd, ca, rbase, rbuf = s
    out = []
    # sender: send or resend a frame of the window
    if len(cd) < C:
        for k in range(sbase, min(sbase + W, T)):
            if k not in sacked:
                out.append(((sbase, max(snext, k + 1), sacked, cd + ((k % M, k),), ca, rbase, rbuf), []))
    # channels: lose any message
    for i in range(len(cd)):
        out.append(((sbase, snext, sacked, cd[:i] + cd[i + 1:], ca, rbase, rbuf), []))
    for i in range(len(ca)):
        out.append(((sbase, snext, sacked, cd, ca[:i] + ca[i + 1:], rbase, rbuf), []))
    # receiver: the head of the data channel
    if cd:
        (sq, g), rest = cd[0], cd[1:]
        ev, nb, buf, ack = [], rbase, set(rbuf), None
        d = (sq - rbase) % M
        if d < W:                                   # in the receive window
            k = rbase + d
            if k != g:
                ev.append("ALIAS")
            else:
                buf.add(k)
                while nb in buf:                    # deliver in order
                    buf.discard(nb)
                    nb += 1
                ack = (sq, g)
        elif 1 <= (rbase - sq) % M <= W:            # previous window: acknowledge again
            ack = (sq, g)
        nca = ca + (ack,) if ack and len(ca) < C else ca
        out.append(((sbase, snext, sacked, rest, nca, nb, frozenset(buf)), ev))
    # sender: the head of the acknowledgement channel
    if ca:
        (sq, g), rest = ca[0], ca[1:]
        ev, acked, nb = [], set(sacked), sbase
        d = (sq - sbase) % M
        if d < W and sbase + d < snext:
            k = sbase + d
            if k != g:
                ev.append("ACK_ALIAS")
            else:
                acked.add(k)
                while nb in acked:
                    acked.discard(nb)
                    nb += 1
        if nb == T:
            ev.append("DONE")
        out.append(((nb, snext, frozenset(acked), cd, rest, rbase, rbuf), ev))
    return out


def run(W, M, T=None, C=None, limit=2_000_000):
    cfg = dict(W=W, M=M, T=T or 2 * W + 1, C=C or W)
    r = check(cfg, limit=limit, succ=succ, s0=(0, 0, frozenset(), (), (), 0, frozenset()))
    return r


if __name__ == "__main__":
    import time
    for W in range(1, 5):
        for M in range(max(1, W), 2 * W + 2):
            t0 = time.time()
            r = run(W, M)
            bad = sorted(e for e in r["events"] if e != "DONE")
            print("W", W, "M", M, "T", 2 * W + 1, "C", W, "states", r["states"], "trunc", r["truncated"], "bad", bad,
                  "done", "DONE" in r["events"], "%.1fs" % (time.time() - t0), flush=True)

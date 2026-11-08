------------------------------ MODULE GalvanicMAC ------------------------------
(* Medium access of the galvanic profile (DESIGN.md 4.3 L2-L4) at superframe granularity.
   N sensors send n frames per superframe each; frame j of superframe c carries the sequence number (n c + j) mod M.
   The beacon acknowledges, per sensor, by a base sequence number (every earlier frame) and a B-bit bitmap (the B
   numbers after the base). In the allocations every new frame gets one exchange; the hub lists each expected frame
   it did not decode and polls the first P frames of that list in a retry pool after the allocations. An adversary
   loses any beacon copy, data frame or poll, within a total budget of Budget losses (Bounded = FALSE: unbounded).
   Variant "age": hub and sensor give a frame up at the same age of A superframes; "count": each side gives up after
   R retries that it counts itself (a missed poll costs the sensor no attempt, the hub cannot tell).
   Ghost state: each frame's age (superframes since it was generated) and whether the hub decoded it. The model
   transcribes mc.py of this folder action for action, so that TLC and the Python checker can be compared state for
   state; a frame is the tuple <<age, tries, decoded, j>>. *)
EXTENDS Integers, Sequences, FiniteSets

CONSTANTS N, n, M, A, P, B, R, CAP, Q, Bounded, Budget, Variant, Track, HubExpiry, HubWindow

Absent == -9
Sensors == 1..N
W == M \div 2
Mod(x) == ((x % M) + M) % M
Min(a, b) == IF a < b THEN a ELSE b

\* the places where the adversary can lose something in one superframe
Slots == {<<"b", i>> : i \in Sensors} \cup {<<"a", i, j>> : i \in Sensors, j \in 0..n-1}
         \cup {<<"p", k>> : k \in 1..P} \cup {<<"d", k>> : k \in 1..P}
\* every set of at most k slots (computed once: a constant definition)
Upto[k \in 0..Budget] ==
    IF k = 0 THEN {{}} ELSE Upto[k - 1] \cup {s \cup {x} : s \in Upto[k - 1], x \in Slots}

VARIABLES c, base, held, sens, pool, rem, falseAck, alias, miss, stuck
vars == <<c, base, held, sens, pool, rem, falseAck, alias, miss, stuck>>

Init == /\ c = 0
        /\ base = [i \in Sensors |-> 0]
        /\ held = [i \in Sensors |-> [s \in 0..M-1 |-> Absent]]   \* hub: decoded numbers above base, ghost age
        /\ sens = [i \in Sensors |-> {}]
        /\ pool = <<>>                                           \* hub's retry list <<sensor, seq, polls>>
        /\ rem = Budget
        /\ falseAck = FALSE /\ alias = FALSE /\ miss = FALSE /\ stuck = FALSE

SeqOf(f) == Mod(n * (c - 1 - f[1]) + f[4])
c1 == Mod(c + 1)
Nxt == Mod(n * c1)

\* 1 beacon
Heard(i, L) == <<"b", i>> \notin L
Acked(i, f, L) == /\ Heard(i, L)
                  /\ LET sq == SeqOf(f) IN
                       \/ 0 < Mod(base[i] - sq) /\ Mod(base[i] - sq) <= W
                       \/ 1 <= Mod(sq - base[i]) /\ Mod(sq - base[i]) <= B /\ held[i][sq] # Absent
Expired(f) == Variant = "age" /\ f[1] >= A
CountOut(f) == Variant = "count" /\ f[2] > R
Kept(i, L) == {f \in sens[i] : ~(Expired(f) \/ Acked(i, f, L) \/ CountOut(f))}
FalseAckIn(L) == \E i \in Sensors : \E f \in sens[i] : ~Expired(f) /\ Acked(i, f, L) /\ ~f[3]
MissIn(L) == \/ \E i \in Sensors : \E f \in sens[i] : ~f[3] /\ (Expired(f) \/ (~Acked(i, f, L) /\ CountOut(f)))
             \/ \E i \in Sensors : Cardinality(Kept(i, L)) + n > Q      \* a full queue refuses a new frame

\* 2 generation: up to n new frames while the queue holds fewer than Q
NewFrames(k) == {<<-1, 0, FALSE, j>> : j \in 0..(Min(n, Q - k) - 1)}

\* the hub decodes frame f (number sq) of sensor i: a different frame held under sq, or a number in the half window
\* behind the base for a frame not decoded before, makes it discard the frame as a duplicate (alias)
Dlv(h, bs, sq, f) ==
    IF h[sq] # Absent /\ h[sq] # f[1] THEN [h |-> h, ok |-> FALSE, al |-> TRUE]
    ELSE IF HubWindow /\ Mod(bs - sq) # 0 /\ Mod(bs - sq) <= W /\ ~f[3] THEN [h |-> h, ok |-> FALSE, al |-> TRUE]
    ELSE [h |-> [h EXCEPT ![sq] = f[1]], ok |-> TRUE, al |-> FALSE]

Send(S, i, f, lost) ==   \* one attempt of frame f; lost: the frame is lost
    LET d == Dlv(S.held[i], base[i], SeqOf(f), f) IN
    IF lost THEN [S EXCEPT !.sens[i] = (@ \ {f}) \cup {<<f[1], f[2] + 1, f[3], f[4]>>}, !.used = @ + 1]
    ELSE [S EXCEPT !.sens[i] = (@ \ {f}) \cup {<<f[1], f[2] + 1, f[3] \/ d.ok, f[4]>>},
                   !.held[i] = d.h, !.alias = @ \/ d.al]

\* 3 allocations, sensor by sensor, frame j by frame j
AllocOne(S, i, j, L) ==
    LET fs == {f \in S.sens[i] : f[2] = 0 /\ f[4] = j} IN
    IF fs = {} THEN S ELSE Send(S, i, CHOOSE f \in fs : TRUE, <<"a", i, j>> \in L)
RECURSIVE Alloc(_, _, _)
Alloc(S, k, L) == IF k > N * n THEN S ELSE Alloc(AllocOne(S, (k - 1) \div n + 1, (k - 1) % n, L), k + 1, L)

ListOne(S, i, j) ==      \* the hub lists an expected frame it did not decode
    LET sq == Mod(n * c + j) IN
    IF S.held[i][sq] = Absent /\ ~\E k \in 1..Len(S.pool) : S.pool[k][1] = i /\ S.pool[k][2] = sq
    THEN [S EXCEPT !.pool = Append(@, <<i, sq, 0>>)] ELSE S
RECURSIVE List(_, _)
List(S, k) == IF k > N * n THEN S ELSE List(ListOne(S, (k - 1) \div n + 1, (k - 1) % n), k + 1)

\* 4 retry pool: the first P entries of the list, once each
PollOne(S, k, L) ==
    LET e == S.pool[k]
        S1 == [S EXCEPT !.pool[k] = <<e[1], e[2], e[3] + 1>>]
        fs == {f \in S1.sens[e[1]] : SeqOf(f) = e[2]}
    IN IF <<"p", k>> \in L THEN [S1 EXCEPT !.used = @ + 1]      \* missed poll: no attempt
       ELSE IF fs = {} THEN S1
       ELSE Send(S1, e[1], CHOOSE f \in fs : TRUE, <<"d", k>> \in L)
RECURSIVE Pool(_, _, _)
Pool(S, k, L) == IF k > Min(P, Len(S.pool)) THEN S ELSE Pool(PollOne(S, k, L), k + 1, L)

\* 5 clean-up of the list, aging, base advance over held numbers (age variant: and over expired ones)
AgeNext(sq) == Mod(n * c1 - 1 - sq) \div n
KeepEntry(S, e) == /\ S.held[e[1]][e[2]] = Absent
                   /\ ~(Variant = "count" /\ e[3] > R)
                   /\ ~(Variant = "age" /\ HubExpiry /\ AgeNext(e[2]) >= A)
RECURSIVE Adv(_, _)
Adv(bs, h) == IF h[bs] # Absent THEN Adv(Mod(bs + 1), [h EXCEPT ![bs] = Absent])
              ELSE IF Variant = "age" /\ HubExpiry /\ bs # Nxt /\ Mod(Nxt - 1 - bs) \div n >= A
                   THEN Adv(Mod(bs + 1), h)
                   ELSE <<bs, h>>
Older(h) == [s \in 0..M-1 |-> IF h[s] = Absent THEN Absent ELSE h[s] + 1]
Capped(h) == [s \in 0..M-1 |-> IF h[s] = Absent THEN Absent ELSE Min(h[s], CAP)]

Choices == IF Bounded THEN Upto[rem] ELSE SUBSET Slots

Next ==
  \E L \in Choices :
    LET S0 == [held |-> held, sens |-> [i \in Sensors |-> Kept(i, L) \cup NewFrames(Cardinality(Kept(i, L)))],
               pool |-> pool, used |-> Cardinality(L \cap {<<"b", i>> : i \in Sensors}), alias |-> FALSE]
        S1 == List(Alloc(S0, 1, L), 1)
        S == Pool(S1, 1, L)
        adv == [i \in Sensors |-> Adv(base[i], Older(S.held[i]))]
    IN /\ c' = c1
       /\ base' = [i \in Sensors |-> adv[i][1]]
       /\ held' = [i \in Sensors |-> Capped(adv[i][2])]
       /\ sens' = [i \in Sensors |-> {<<Min(f[1] + 1, CAP), f[2], f[3], f[4]>> : f \in S.sens[i]}]
       /\ pool' = SelectSeq(S.pool, LAMBDA e : KeepEntry(S, e))
       /\ rem' = IF Bounded THEN rem - S.used ELSE rem
       /\ falseAck' = (falseAck \/ ("falseAck" \in Track /\ FalseAckIn(L)))
       /\ alias' = (alias \/ ("alias" \in Track /\ S.alias))
       /\ miss' = (miss \/ ("miss" \in Track /\ MissIn(L)))
       /\ stuck' = (stuck \/ ("stuck" \in Track /\ \E i \in Sensors : \E f \in S.sens[i] : f[1] + 1 >= CAP))

Spec == Init /\ [][Next]_vars

NoFalseAck == ~falseAck   \* no sensor drops as acknowledged a frame that the hub never decoded
NoAlias == ~alias         \* the hub never discards a new frame as a duplicate of another
NoMiss == ~miss           \* no frame is given up undecoded or refused by a full queue (deadline A in the age variant)
NoStuck == ~stuck
=============================================================================

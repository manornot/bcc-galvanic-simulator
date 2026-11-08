# Model checking of the galvanic profile's medium access

Two checkers explore the same superframe-level model: `mc.py` (explicit-state, breadth first, Python) and
`GalvanicMAC.tla` with TLC. The run details are kept here; the numbers are those of the last run, with the sensor queue Q = n(A + 1) under the age rule and
Q = 4 under the attempt rule, with a frame refused by a full queue counted as a deadline miss.

## Tools

- Python checker: `mc.py`; drivers `sweep_seq.py`, `sweep_n.py`, `sweep_hub.py`, `faults.py`, `scale2.py`,
  `textbook.py`, `agree.py`; unit tests `test_formal.py` (9 tests). Each driver writes the `.out` file beside it.
- TLC: TLC2 version 2.19 of 8 August 2024 (rev 5a47802), OpenJDK 21.0.12.1, from
  `tla2tools.jar` (set `TLA2TOOLS` to its path), run by `tlc.sh` (soundness) and `tlc_deadline.sh` (deadline) with
  `-workers 2 -Xmx3g`. Each run writes `tlc/NAME.cfg` and `tlc/NAME.out`. A TLC run stops at the first violation
  and then gives no complete state count; the Python checker explores every reachable state also after a violation.

## Validation of the checkers

- `textbook.out`: selective repeat over two lossy channels, W = 1 to 3 with 2W + 1 frames and channels that hold W
  messages, W = 4 and 5 with W + 2 frames and channels of two; ALIAS for every M < 2W, none at M = 2W; W = 1, M = 2
  (alternating bit) correct, M = 1 aliased.
- `test_formal.out`: four hand-traced superframes as tests of the step function; mutations `hub_expiry False`
  (alias follows), `M = 2A` (alias), `B = 0` (sound, acknowledgement later), `var 'count'` (false acknowledgement),
  and a queue of A frames (`drop_queue` for A 2 and 3 with L = A, none with A + 1).
- `agree.out`: 16 configurations explored completely by both tools with equal distinct-state counts, 9 failing
  configurations found by both.

## Configurations and state counts (Python / TLC)

Equal counts: `age_A2_M5` 586, `age_A3_M7` 3886, `age_A4_M9` 21547, `n2_A2_M11` 19504, `count_M8_L2` 105,
`count_M8_R3_L3` 305, `dl_A3_L3_N10` 20280, `dl_A3_L4_N2` 2668, `dl_A4_L4_N1` to `N6` 892, 3448, 8695, 17731,
31816, 52354, `mut_B0_A2_M5` 516, `sound_N2_A2_M5` 13906.

Violations found by both: `age_A2_M4`, `age_A3_M6`, `age_A4_M8`, `n2_A2_M9`, `n2_A2_M10` (NoAlias),
`count_M8_L3`, `count_M8_R3_L4` (NoFalseAck), `dl_A3_L5_N2` (NoMiss), `mut_noexpiry_A2_M5` (NoAlias).
Python counts of the failing configurations: 393, 2995, 18001, 15471 (M 10), 425, 929, 8667.

Not run with TLC: `A 4, L 4, 10 sensors` (Python 232156 states, no miss).

Run times: `sound_N2_A2_M5` 2 min 30 s; `dl_A3_L3_N10` 4 min 39 s; `dl_A4_L4_N6` 5 min 24 s;
`count_M8_L3_workers1` gives the shortest trace (7 states = 6 superframes) with one worker.

In an earlier model, with Q = nA and only `drop_retry` counted as a miss, the counts were 536, 4572, 29368, 15104,
20420 and 919 to 52516 for the age and deadline runs, and N 2 and 3 with P 1 missed at L = A + 1 (the refused
frame stayed on the hub's retry list and took the poll).

## Trace of the attempt-rule failure (count_M8_L3)

Frame 0, its pool retry and its second pool retry are lost; the sensor's tries exceed R = 2, then the hub's polls
exceed R; the base stays at 0; frame 4 is discarded as a duplicate since (0 - 4) mod 8 = 4 <= M / 2, and the next
beacon acknowledges it. With R = 3 the same needs four losses. With M = 2^15 the first false acknowledgement comes
after half the space, 16384 frames.

## Sequence-space bound

`sweep_seq.out`: smallest sound M 5, 7, 9 for A 2, 3, 4 with n 1 (B 1 and B A). `sweep_n.out`: n 2, A 2: ALIAS at
M 8, 9, 10, none at 11; n 3, A 2: 17; n 2, A 3: 15 (B n and B 1). All follow M >= 2n(A + 1) - 1. A hub without
the half-window duplicate rule is sound with M 4, 5, 7 for A 2, 3, 4 (`sweep_hub.out`). With 15 bits and A = 4:
2 x 3276 x 5 - 1 = 32759 <= 32768 < 2 x 3277 x 5 - 1; at A = 8 (nine periods held) 1820 frames.

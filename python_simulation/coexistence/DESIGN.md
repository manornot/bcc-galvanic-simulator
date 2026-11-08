# Coexistence simulator: design contract

This file fixes the semantics of `python_simulation/coexistence/`. Parameter values and the per-protocol
simulation recipes live in `PARAMETERS.md` (machine copy `parameters.json`, read by `protocols.py`).
Precedence: for what a protocol does on air (frame, timing, retries, channel rule) the recipe of
`PARAMETERS.md` wins; for how the engine decides outcomes and what a metric means, this file wins. Where the
code disagrees with either, the code is wrong or the document is updated in the same commit.

## 1 Question the simulator answers

k bodies share one room. Each body carries one body-area network (BAN): one hub and N_dev sensors. Each
sensor produces L payload bits every period T, so the per-link rate is R = L / T and the aggregate per-body
rate is R_agg = N_dev x R (the axis of `viz_headline.py`). All k BANs run the same protocol and
variant. The simulator reports, per BAN and trial, how much of the offered traffic arrives, how late, with how
many bit errors, at what energy, and how the medium was used.

The comparison runs stacks at their full capability : every protocol with every capability
its standard gives, at every rate or PHY it offers, and BCC with every feature of the IEEE 802.15.6 HBC standard stack on the galvanic channel (`bcc` full), on one tone (phy `sf8`) and striped over the N parallel tones of the plan, with the galvanic-optimisation ladder as the proposal (4.3). Variant `thesis` marks the as-built models:
the stack as coded (`bcc` thesis), the as-built reference, and the plain BLE model (`ble4thesis`). The fair
variants are gone (section 2). Energy: the radio stacks on their chips' system currents (nRF52840 compounded), BCC and
HBC on transceiver-IC figures; no nRF overhead is added to BCC (section 6).

Protocols (ids used in code and CSV):

| id | protocol | variants |
|---|---|---|
| `ble` | Bluetooth LE 5, PHY 1M, 2M, Coded S=2, Coded S=8 | full (ACL: poll, LL ACK, retransmission in the same event, CSA#2, AFH) |
| `ieee802154` | IEEE 802.15.4 O-QPSK 2.4 GHz | full: phy `oqpsk` mode A (unslotted CSMA/CA, ACK, retries, ED-scan channel choice), phy `gts` beacon mode B (one GTS per sensor, ACK, retries inside the GTS, ED-scan channel; 4.1) |
| `tsch` | IEEE 802.15.4 TSCH | full (ACK, retries in dedicated and shared cells, blacklisting, CCA before every frame: macTsCca, 4.1) |
| `ieee802156nb` | IEEE 802.15.6 NB 2.4 GHz, 121.4, 242.9, 485.7 and 971.4 kb/s | full (BCH(63,51), spreading, I-ACK, retries, scan, the rate's draft sensitivity; not modelled: the optional re-scan, beacon shifting, superframe interleaving and hub hopping, 4.1) |
| `smartban` | ETSI SmartBAN, PHY schemes `bch1`, `bch2`, `bch4`, `unc1`, `unc2`, `unc4` (PSDU BCH or uncoded, Nrep 1, 2, 4) | full (header BCH, ACK, retries in own and C/M slots, the C/M ones by the MCA backoff, MCA sensing, channel switching) |
| `ant` | ANT (private network) | full (acknowledged messages, retries, advanced burst, listen at channel open with the least-occupied phase as fallback, adaptive isochronous, one private frequency per BAN with frequency agility over it and two more, 4.1) |
| `hbc` | IEEE 802.15.6 HBC standard (capacitive, 21 MHz, 42 Mchip/s), FSC SF 8, 16, 32 and 64 | full (FSC, FCS, I-ACK, 3 retransmissions); coupled between bodies (PARAMETERS 5.7) |
| `bcc` | BCC on the galvanic channel | full (the HBC stack of PARAMETERS 3.7 FULL with R_c as the Walsh-chip rate of the SF 8 mode: SF 64 header, MAC header + FCS, I-ACK, 3 retransmissions, NB MAC allocations; phy `sf8`), full multi-tone (phy `multitone`: the same stack on the N tones of `bcc_fdma.plan`, spacing 4 R_c over the 3 dB span of the cohort 50 Ohm curve (6.48 to 20 MHz, dip included), 6 tones at 500 kchip/s, 4 of them inside the 3 dB bins (4.1); a frame striped over the tones, one BCC full packet per tone, I-Ack on every tone, drive split by sqrt(N), N front ends), full at ladder level 1 to 5 (flag `opt`, `--bcc-opt`, label `full+optN`, one tone and multi-tone: the HBC layers galvanic isolation makes redundant switched off one at a time, levels 1 and 3 802.15.6 options, 2, 4 and 5 proposals; 4.3), thesis (the stack as coded, `HBCProtocol.py`: spreading off, no FCS, no ACK; the as-built reference); options for all three: `--bcc-rc`, `--bcc-power`; for the single-tone stacks: `--bcc-contact` (measured-hardware case), inter-body coupling +40 / +60 dB sensitivity cases (`get_spec` refuses them for multi-tone); `--bcc-fcs` (+2-octet FCS) for the as-built stack only (`get_spec` refuses it for BCC full, which carries the 802.15.6 FCS); `get_spec` also refuses multi-tone at a `--bcc-rc` above 3.38 Mchip/s (no 4 R_c tone fits the 3 dB span) and all three above 9 Mchip/s; `sweep.py` skips every entry `get_spec` refuses, with a printed line |
| `ble4thesis` | `BLE_SIM.py` rules inside this engine (no framing overhead, data packets only) | thesis: the plain BLE model, for the reproduction test and the side-by-side with |

## 2 Fair variant (removed)

The fair variants (each stack stripped of the features the BCC stack lacks: FEC, repetition and spreading, ARQ,
carrier sensing, adaptive channel selection, link adaptation) are removed from the simulator: the comparison is between stacks at their full capability (section 1). Their definition is this section in the
git history before (commit 97e0720), and their results stay in `results/raw/v1`.

## 3 Traffic and period selection

- Each sensor generates L bits per period T. If L exceeds the protocol's maximum application payload per
  frame, the sensor enqueues ceil(L / L_max) frames per period (last one shorter). Each frame carries the
  protocol's framing overhead.
- One reference period P per cell, the same for every protocol (PARAMETERS OI-4: payload and period from the
  `BLE_SIM` axis R = L / CI; OI-19 and 3.2 mode B: every protocol runs at P rounded to its own grid). It is
  `BLE_SIM.select_ci`'s CI rule on the BLE CI grid (`reference_period`): `latency` the shortest CI carrying
  L >= 8 bits (7.5 ms from 1.07 kb/s up), `energy` the longest CI with L within one 251-octet PDU (at least
  7.5 ms, at most 4 s).
- `select_period(r_kbps, mode, n_dev, spec)`: T = the shortest period of the protocol's grid at or above P whose
  schedule holds the N_dev sensors (scanned upward), L = round(R T) bits (at least one byte). The schedule holds when
  the time of one period of all N_dev sensors (all frames, polls, IFS, ACKs, beacon) fits T for scheduled access and
  ANT, or `csma_max_load` x T for CSMA (0.5, section 4.2), with the recipe's slot rules: TSCH one cell per frame after
  the shared cell, SmartBAN one LSLOT slot per frame with TMUA counted, the 2 C/M slots and T a multiple of TS, NB, HBC
  and BCC FULL allocations, 802.15.4 GTSs (of whole superframe slots) and BLE FULL connection events of 1 + 1 spare
  exchanges per frame (`protocols.SPARE_EXCHANGES` replacing PARAMETERS 3.4 "sized for (R + 1) full
  exchanges"; a provisioning rule of this model, not a feature of any of the standards: see the pooled reference stack
  of 4.3; the GTS by the same rule, A, where PARAMETERS 3.2 says 'retry inside the GTS if it fits'; BLE from:
  a failed data PDU goes again in the same event; later retries go to the next allocation, GTS or event), ANT
  acknowledged messages at T / 1.25 >= the N >= 171 clamp. Every rate of a stack scans its own grid. Infeasible cells
  are written with `feasible = 0` and no metrics.
- Deadline of `on_time`: D = P, the same absolute deadline for every protocol of the cell (one sample period,
  OI-4), even where the protocol's grid or schedule makes T longer than P: a frame then waits up to T for its
  access opportunity, so about P / T of the frames at most are on time.
- Generation time of a frame: the sensor's own sample clock, period T with a random phase per sensor (drawn from
  the cell's substream, so every protocol and variant of a cell sees the same phases), for every access type. A
  scheduled sensor sends the frame at its next slot or connection event, a CSMA sensor starts its backoff at once,
  an ANT master sends it in its next channel period (its slot phase is independent of the sample clock). ponytail:
  a scheduled sensor samples on its hub's superframe clock, so sampling never slips against the slots. Latency is
  delivery time (end of the successful frame at the hub) minus generation time, so it includes the wait for the
  access opportunity (`BLE_SIM.py` stamps a frame at its slot start; only its PSR is reproduced).
- Retry room of the FULL schedules:
  every FULL stack with ARQ gets retry room where its standard provides a place for a retry. BLE events, the 802.15.6
  allocations (NB, HBC, BCC full) and the 802.15.4 GTSs hold one spare exchange per frame: the standards permit a retry
  in the same allocation or event when time remains, but none sizes the allocation for it (802.15.6 sizes allocations
  in slots, 5.7.2.5 and 6.7.1, defers a retry to the next allocation interval when time is short, Table 22, or grants
  a polled allocation, 6.6.1), so the per-frame spare is this model's provisioning rule and not a standard feature
  (the pooled reference stack of 4.3 provisions the retries of BCC full by the pool instead); TSCH sends the
  frames left over from an earlier
  slotframe in the slotframe's shared cell (the RFC 8180 minimal cell, before the dedicated cells) with the shared-cell
  backoff; SmartBAN sends them in the 2 C/M slots of the IBI by the MCA backoff (4.1); 802.15.4 CSMA retries at once and
  its admission holds the same room, a scheduled load of at most 1 / (1 + 1) of T (4.2). Later retries go to the next
  event, allocation, shared cell or own slot, ahead of new data. Where that room is contended (the BAN's sensors share
  the one TSCH shared cell and the SmartBAN C/M slots, entered after a backoff) or used up, a lost exchange still delays
  the later frames of its sensor by one period until its FIFO (4 periods) is full, after which new frames are dropped at
  arrival, so on_time, latency and even psr can follow a backlog rather than the ARQ. measured the exactly
  provisioned schedules (one exchange per frame; collision model, 2 s, 6 seeds): BLE 1M, 10 kb/s, N_dev 10 (T 7.5 ms): k
  4 psr 0.917, on_time 0.030, drop_queue 0.075; TSCH 1 kb/s, N_dev 4, k 8 (T 50 ms): psr 0.921, drop_queue 0.056; the
  2-exchange BLE event (T 12.5 ms) gave k 4 psr 1.000, on_time 0.54 (D = P punishes the longer T). Not taken: letting a
  backlogged sensor use the idle rest of the period (in no recipe). ANT FULL bursts lag the same way without a slot: a
  burst carries the
  frames sampled before it starts, and one that runs past the channel period resumes at the next whole period. At 30
  kb/s latency (T 7.507 ms, 29 octets = 2 packets = 6.4 ms per period) a single retry stretches a burst to 9.6 ms; from
  then on the master sends 4 packets every 2 periods (12.8 of 15.0 ms) and never catches up, so every later frame is
  about one period late. N_dev 1, k 1, collision, p_loss 0, seed 0: latency 1.4 / 4.5 ms and on_time 1.000; with the
  first packet of the run lost, 7.8 to 12.0 ms after it and on_time 0.002. Where the N_dev bursts fill the period (10
  kb/s latency x 4: one 16-octet packet per master, 4 x 3.2 ms in T = 12.817 ms) a retry 3.2 ms later meets the next
  master's packet, whose retry meets the next, and the losses cascade (seed 0, the first packet lost: psr 0.500,
  drop_queue 39 of 150; seeds 1 and 2 place the retry in a gap). The recipe gives a burst no retry room, and the hub's
  pre-emption of its other channels by a burst (3.6 common recipe) is not modelled: the masters collide instead.

## 4 Time structure per access type

Times in microseconds, integer arithmetic where the standard allows it.

### 4.1 Scheduled access (TDMA family)

Per BAN: superframe of length T with random phase phi_j ~ U[0, T) against the other BANs (as in `BLE_SIM.py`).
Layout: optional hub beacon at the superframe start (airtime `beacon_bits / phy_rate` plus IFS), then
N_dev equal slots, or the recipe's quantised slots where `select_period` returns them (`slot_us`,
`reserved_us`): TSCH the shared cell (slot 0: the EB, and in FULL the frames left over), one 10 ms dedicated cell
per frame after it, the rest of the slotframe idle; SmartBAN the D-Beacon slot, one TS slot per frame and node, then
the 2 C/M slots (FULL retries) and the inactive slots; 802.15.4 GTSs contiguous at the end of the SD.
Inside its slot a sensor runs exchanges back to back while the next exchange fits before the slot end (one per
quantised slot). An exchange is the protocol's sequence, for example:

- BLE full (ACL): hub poll (empty PDU) + T_IFS + sensor data PDU + T_IFS. The acknowledgement of the data
  PDU travels in the next hub packet (NESN): a failed data PDU is retransmitted in the next exchange of the
  same event, which holds one spare exchange per frame (section 3), else in the next event. A missed packet closes
  the event (a missed poll: no answer);
  a poll received with a CRC error is still answered, with the unacknowledged head PDU; the event also closes
  after two consecutive CRC errors at either side (PARAMETERS 3.1 closing rules). The SINR model tells a miss
  (sync) from a CRC error; the collision model treats a failed poll as missed and a failed data PDU as a CRC
  error. AFH scores a missed poll as a failure on its channel. No retry limit (ACL); the queue per sensor holds
  `queue_len` frames and a frame that arrives to a full queue is dropped (`drop_queue`).
- Beacon loss (SmartBAN, 802.15.4 GTS): the beacon is decided at every sensor like a frame; a sensor that has missed
  `max_lost` consecutive beacons (SmartBAN 10, OI-33; 802.15.4 aMaxLostBeacons 4) loses the frames of that superframe
  (`drop_sync`) until it hears a beacon again.
- 802.15.4 GTS: the hub's beacon (PSDU 13 + 1 + 3 N_dev octets) at the superframe start, the CAP of at least
  440 symbols, then one GTS per sensor (N_dev <= 7), contiguous at the end of the SD, of the whole slots that hold 1 + 1
  exchanges per frame; SO the smallest that leaves the CAP. Exchange: data + 192 us + Imm-Ack + SIFS / LIFS, no CCA (no
  CSMA in the CFP); a failed frame goes again in the same GTS, else in the next superframe's (up to 3 retries). The
  CAP stays idle (the sensors send nothing by CSMA). The spare exchange costs capacity: 10 kb/s latency, N_dev 4 runs
  at T = 61.44 ms (CSMA 30.72 ms).
- TSCH: every exchange sits in its own timeslot of `timeslot_us`, one dedicated cell per frame; full: CCA of 128 us at
  macTsCcaOffset (1800 us), data at `ts_tx_offset`, ACK after `ts_tx_ack_delay`; a busy CCA counts as a failed attempt
  (PARAMETERS 3.3, OI-16). The CCA is the standard's macTsCca option, off in Contiki-NG and in the PARAMETERS default,
  on since (full capability, as blacklisting): in the collision model it lifts psr at 10 kb/s latency, N_dev 4, k
  8 from 0.657 to 0.745; in the SINR model it is about neutral (10 kb/s latency N_dev 4 k 30: 0.781 with it, 0.788 without). A
  failed frame retries in the sensor's next dedicated cell of the slotframe if it has one, else in the next shared cell
  (RFC 8180 minimal cell: slot 0, channel offset 0, where the hub also sends its EB and so cannot receive), which comes
  before its next dedicated cell (up to `max_retries`). Every frame still queued at the shared cell (failed, or behind a
  failed one) goes there; ponytail: the period's new sample joins the queue after the shared cell (samples on the
  superframe clock). Shared-cell backoff (802.15.4-2015 6.2.5.3, PARAMETERS 3.3): none before a first shared-cell
  attempt; after a failed one BE = macMinBE (1), then min(BE + 1, macMaxBE (7)), and the sensor skips U{0..2^BE - 1}
  shared cells; a success in a shared cell, or one in a dedicated cell that empties the queue, resets it; a failure in a
  dedicated cell leaves it (Contiki-NG). The sensors of a BAN contend in the one shared cell (same cell and channel:
  they collide). Channels hop with blacklisting (below).
- 802.15.6 NB coexistence options not modelled (PARAMETERS 3.4 FULL and MAC row, OI-25): the re-scan when the BAN's
  PER exceeds 10 % over 10 superframes, beacon shifting, active superframe interleaving and hub channel hopping. All
  are optional and marked A in PARAMETERS, and each is hub-level coordination whose control traffic (beacon fields,
  negotiation with neighbouring hubs) the engine does not carry; NB keeps the channel of its ED scan (4.1 below). The
  re-scan could reuse the SmartBAN migration rule.
- 802.15.6 NB, HBC and BCC full: data + pSIFS + I-ACK; retransmission in the same allocation (it holds one spare
  exchange per frame), else in the next one, up to `max_retries` (multi-tone BCC: the whole stripe, below). BCC full
  from ladder level 2: no I-Ack, the beacon acknowledges; from level 3 a retry pool replaces the spare exchanges (4.3).
- SmartBAN PHY schemes (PARAMETERS 3.5, TS 103 326 Table 3; the hub assigns one statically, so each is a
  catalogue entry and the report takes the best per cell): PSDU BCH(127,113) or uncoded, the PPDU repeated Nrep = 1, 2
  or 4 times, airtime Nrep (104 + L_PSDU) us, LF,max per scheme from the LSLOT 32 slot rule (BCH 2150, 1065, 522;
  uncoded 2418, 1198, 588 octets). SINR model: ideal EGC (phy rep = Nrep: SINR x Nrep per decision, and the C/I floor
  lowered as for any repetition, 5.2); ponytail: each block is stretched over the Nrep copies, so an overlap of one copy
  hits every block in proportion instead of that copy's blocks. Collision model: each copy decided alone (the 3.5 BER
  model: PER = the product over the copies), the frame survives when one copy is clear, then p_loss^(clear copies).
  The ACK (168 us) and the D-Beacon (272 us) are not repeated. Nrep 1 BCH is the stack as before.
- SmartBAN full: data + T_IFS + ACK + T_IFS in the node's TS slot after the MCA sensing (TMUA 16 us). A frame still
  queued after the node's slot (failed, or its slot went to a retry or was sensed busy; TS 103 325 7.3.2.0: C/M
  access 'when it has not been allocated sufficient time slots') is sent in the IBI's 2 C/M slots by the MCA backoff
  (D2; PARAMETERS 3.5 FULL: 'BE/BC applies to RCA, C/M and PCA'; TS 103 325 7.3.3.2.2, 7.3.3.1.1:
  with the MCA bit set in the D-Beacon the C/M period uses the backoff algorithm): per head frame BE = BEmin and a
  backoff counter BC = U{0..2^BE - 1} of C/M slots, counted down one per C/M slot; at BC 0 the node senses TMUA =
  Tsch (16 us, Table 2 'C&M access, any traffic') at the slot start and, idle, sends its data after it, the ACK in the
  same slot; busy, or sent and not acknowledged, BE = min(BE + 1, BEmax) and a new BC. BEmin 3, BEmax 5 are not in
  the TS (the 802.15.4 defaults, A in PARAMETERS 3.5). The BAN's nodes contend in the same slots: they sense at the
  same instant, so two with BC 0 in one slot both find it idle and collide, in both reception models; the sensing detects frames of other BANs on the DCH. A new frame takes a C/M slot of its IBI
  with probability 2 / 8 (BC 0 or 1; slotted ALOHA at user priority 0, the reading, took one with 1 -
  (7/8)^2 = 0.23). A frame not sent there waits for the node's
  next own slot, ahead of new data. RT = 1 (the recipe's retry limit, not the MCA's RTmax 3): one retransmission,
  then the frame is dropped. Not modelled: re-use channel access (RCA) of other nodes' idle scheduled slots (a
  'may').
- ANT: no common TDMA. Each sensor is a master with its own channel period and phase; all sensors of all BANs on the
  same RF channel collide when they overlap. Acknowledged messages at T / 1.25 (the hub answers in the same period), a
  failed message is resent in the next period (up to `max_retries`), advanced burst when L exceeds 8 octets, adaptive isochronous
  re-timing: the master senses a 500 us guard window after its slot and
  moves 250 us per period earlier while a co-channel frame starts in it. After a message the master is already listening
  for the ACK (W, 274 to 450 us into the slot) and stays in RX to the guard's end (918 us: 644 us of RX, one start-up);
  after an old-payload broadcast the guard is 500 us with its own start-up. The guard is the feature's price: at k 1
  (N_dev 4, collision) J per bit 728 against 421 nJ without it at 5 kb/s latency, 3460 against 2000 nJ at 1 kb/s (booked
  with a start-up after every message, as first built, 803 and 3814). Private network: each BAN
  on one frequency drawn from the 76 of 2403..2480 MHz without 2450 and 2457 (3.6 recipe), so masters of other BANs meet
  only when their BANs drew the same frequency. Listen at channel open (3.6 FULL step 1): masters open one after
  another; each tries up to 100 random phases for one whose occupancy, widened by the 500 us guard, overlaps no master
  it heard, and when none is free takes the least-occupied phase (D4): midway between the two heard
  masters farthest apart, the centre of the largest gap, which overlaps the least heard airtime and leaves the most
  clearance for the drift. The greedy placement cannot pack a crowded period: N_dev 10 at T = 7.5 ms (tc = 6.0 ms, 10 x
  418 us) still overlaps its own masters: PSR 0.795 without other BANs (1 kb/s latency, collision, 3 seeds of 2 s),
  0.384 with the recipe's random-phase fallback; N_dev 4 1.000 with either. Frequency agility: each BAN draws two
  more frequencies of the pool beside its private one (the default 2403 / 2439 / 2475 MHz would put every BAN on the
  same three), and a master moves to the next of the three after 8 consecutive unacknowledged messages or burst packets.
  ponytail: the hub's channel follows at once and the master keeps its phase (no new listen); in the recipe a slave
  follows only after 8 s without reception (4 times a 2 s run) and then searches, which the engine does not model, so
  the move is optimistic by that gap.
- BCC and HBC: TDMA of `MAC_layer.TDMA`: hub sync sequence at the frame start, slots with guard, one tone per BAN. BCC
  full multi-tone runs every exchange on the BAN's N tones in parallel: `spec.tones` from
  `bcc_fdma.plan`. A frame of L
  octets is striped (`protocols.tone_octets`): ceil(L / N) octets on the first L mod N tones, floor(L / N) on the
  others, each tone one BCC full packet (preamble, SFD, SF 64 header, MAC header + FCS, pilots; with L < N the last
  tones send MAC header and FCS alone, so every tone is active), so the frame airtime is the longest tone packet (250
  octets: 5952 us against 19776 us on one tone); the hub answers with an I-Ack on every tone in parallel (one empty
  packet's airtime). The frame is delivered when every tone packet is, the ACK when every I-Ack arrives. ARQ resends the
  whole stripe (ponytail: the engine's ARQ is per frame; resending only the unacknowledged tone packets would take at
  most the same airtime and fewer front ends, so the shortcut costs energy after a failure, not capacity; at 1 V no tone
  packet fails, 5.2). Allocation: one exchange + SPARE_EXCHANGES, as BCC full; frames of up to N x 255 octets (one full
  packet per tone). The beacon stays the single-tone 17-octet BCC full beacon, 4352 us on one tone, sent and heard on
  one front end (section 6; striped over the N tones it would take 3456 us on N front ends). One body carries at most 8
  x 1530 bits per 2 x 23510 us = 260.3 kb/s at 500 kchip/s, against 43.4 kb/s on one tone. PARAMETERS 3.8's FDMA option asks
  for a tone spacing of at least 4 R_c x 1.25 = 5 R_c within the usable band; taken instead: 4 R_c over the whole 3 dB
  span, as the ticket asks and found free of adjacent-channel loss. Under the 3.8 rule the plan holds 5 tones (4
  inside the 3 dB bins): 5 x 2040 bits / 47020 us = 217 kb/s per body, 173.5 kb/s with the in-band tones only.
- BCC at every chip rate: the as-built rates (50 k, 200 k, 500 k, 2 M), the HBC Walsh-chip rates (656.25 k,
  1.3125 M, 2.625 M, 5.25 M) and upward to 9 M, the last rate at which one tone makes no chip error on any measured
  path at NF 10 to 30 (`results/bcc_rate_stability.md`, real passband from 7 M: at 9.5 M the lower tone at 1 MHz
  fades with its image and errs, although every packet still arrives; at 9.95 M every path falls below PSR 0.99 at
  every NF; both refused). The Walsh chip counts do not change with R_c, pSIFS stays 75 us,
  airtimes round up to whole us. Capacity per body (kb/s; the scan-free bound of `select_period`, 4.3's rule for the
  proposal, without the beacon; hand-derived in test_protocols.test_bcc_rates). The proposal (L5) is multi-tone up to
  2.625 Mchip/s and one tone (= L4) above:

| R_c (chip/s) | 50 k | 200 k | 500 k | 656.25 k | 1.3125 M | 2 M | 2.625 M | 5.25 M | 7 M | 8 M | 9 M |
|---|---|---|---|---|---|---|---|---|---|---|---|
| tones of the multi-tone plan | 67 | 16 | 6 | 5 | 2 | 1 | 1 | none | none | none | none |
| BCC full, one tone | 4.36 | 17.42 | 43.39 | 56.83 | 112.71 | 170.28 | 221.74 | 429.47 | 560.44 | 633.54 | 703.93 |
| BCC full, multi-tone | 292.36 | 278.74 | 260.31 | 284.14 | 225.41 | 170.28 | 221.74 | refused | refused | refused | refused |
| proposal, pool 0.05 | 779.50 | 743.55 | 689.85 | 752.36 | 588.02 | 430.17 | 561.41 | 1097.89 | 1442.23 | 1634.62 | 1824.04 |
| proposal, pool 0.43 | 564.50 | 538.29 | 499.16 | 544.22 | 424.88 | 309.50 | 403.45 | 785.23 | 1028.38 | 1163.45 | 1296.05 |
| as-built stack | 11.58 | 46.31 | 115.78 | 151.95 | 303.91 | 463.11 | 607.72 | 1215.08 | 1619.78 | 1852.44 | 2083.52 |

  Multi-tone carries 260 to 292 kb/s up to 656.25 kchip/s (N falls as 1 / R_c); with 2 and 1 tones the integer
  packing leaves span unused (225, 170, 222 kb/s). One tone grows less than R_c: each exchange keeps its 150 us of
  pSIFS.

Channel per exchange:

- `perm`: per BAN a random permutation of the data channels, one channel per superframe, cycling
  (`ble4thesis`, identical to `BLE_SIM.py`).
- `csa2`: per BAN and superframe a pseudo-random channel from the used channel map (BLE full).
- `tsch`: channel = hopping_sequence[(ASN + channel_offset_j) mod n_channels] per timeslot.
- `fixed`: one channel per BAN for the whole run: uniform random (ANT, its private frequency; HBC and BCC have one
  channel); with ED scan (802.15.4, NB, SmartBAN) BANs are formed in random order and each picks the channel with the
  fewest BANs already on it (ties random); 802.15.6 NB picks the least energy instead (PARAMETERS 3.4), BANs counted
  on their channel and, with the 0 dB adjacent rejection, on both neighbours.
- AFH (BLE), blacklisting: the hub keeps the failures per channel
  over its last `adapt_window` outcomes; a channel with more than `adapt_threshold` x window failures is removed from
  the map (BLE: PDR < 0.95 over 20, at least 10 channels kept, else the map resets to all; TSCH: PDR < 0.8 over 100, so
  the 21st failure in 100 removes it, at least 4 kept). Against other TSCH BANs, whose interference is spread over all
  channels, it removes the channels a BAN fails on most, and a heavily hit BAN shrinks to 4 (5 kb/s energy, N_dev 4, k
  30, 40 periods, seeds 0 and 1: 3 of 30 BANs at 4 channels, 22 and 24 keep all 16; psr 0.638 over 4 seeds against 0.611
  without it). At 10 kb/s energy (N_dev 4, k 30, 3 frames per sensor per 200 ms slotframe of S = 20 cells) 174 of 300
  BANs shrink to 4 channels (10 trials), and because S is a multiple of the sequence length the channel of a cell then
  repeats every slotframe (with 16 channels it cycles through 16 / gcd(20, 16) = 4): a cell that meets another BAN's
  cell on its channel meets it every slotframe. on_time falls from 0.107 to 0.048 (k 8: 0.579 to 0.507), psr stays
  (0.346 to 0.353; k 8: 0.808 to 0.850), 20 paired trials. The recipe's slotframe is P rounded to 10 ms cells, not a
  length coprime with the sequence. (These numbers predate the TSCH CCA.) Where it acts:
  the rule needs 21 failures on one channel among its last 100 outcomes (ponytail: a young window needs the full-window
  count, so a channel seen fewer than 21 times is never removed, even at PDR 0), which a 2 s run reaches only where a
  BAN sends many frames per slotframe. BANs of 30 with a removed channel, k 30, SINR, the sweep seed of trial 0, with
  CCA: N_dev 1 none at any R, both modes; latency N_dev 4 none but 1 at 0.1 kb/s, N_dev 10 2 to 3 up to 1 kb/s
  and 11 at 5 kb/s; energy N_dev 4 0 or 1 up to 5 kb/s and 10 at 10 kb/s, N_dev 10 0 to 13 up to 1 kb/s and 29 at 5
  kb/s. In the collision model at 10 kb/s latency, N_dev 4, k 30, a channel with 16 failures in 17 outcomes stays. So
  most TSCH rows are TSCH without a removal, and the report credits blacklisting (D15) only in the cells above where
  BANs lose channels. Channel switching (SmartBAN full): the uplink failure ratio of the BAN over its last 10
  IBIs above 20 % moves it to the channel with the lowest load, 2 beacons later (PARAMETERS 3.5). The adaptation cost
  (control messages) is neglected and stated.

### 4.2 CSMA/CA (802.15.4 full)

Unslotted CSMA/CA as in the standard: NB = 0, BE = macMinBE; wait random(0 .. 2^BE - 1) x aUnitBackoffPeriod;
CCA of 8 symbols; busy -> NB += 1, BE = min(BE + 1, macMaxBE), channel access failure when NB > macMaxCSMABackoffs
(frame dropped, `drop_csma`); idle -> aTurnaroundTime -> frame; the hub replies with an ACK after
aTurnaroundTime; the sensor waits macAckWaitDuration; missing ACK -> retry from CSMA up to macMaxFrameRetries
(`drop_retry` after the last). The hub is half-duplex: a frame that arrives while the hub transmits an ACK is
lost. The hub (coordinator) keeps its receiver on all the time (non-beacon PAN). At the end of the run a sensor, like
a scheduled slot, finishes the sample it is serving (its later fragments and retries run past the end); frames of a
later sample still queued then count as `drop_end`.

Admission: a CSMA cell is feasible when the scheduled load, N_dev
x the exchanges of one period (CCA, turnaround, frame, ACK, IFS; no backoff), is at most csma_max_load = 1 / (1 +
SPARE_EXCHANGES) = 0.5 of T (`protocols.CSMA_MAX_LOAD`), the room of the one spare exchange per frame that the 802.15.6
allocations and the BLE events hold (section 3). admitted up to 1.0, where one BAN alone loses to its own
contention (two sensors whose CCAs end within a few microseconds both find the channel idle, a sensor whose CCA ended
just before the hub's ACK to another sensor sends into that ACK, the hub being half duplex, and sensors that find the
channel busy five times drop the frame) and backlogs. k 1, p_loss 0, collision, 40 periods, seeds 0 to 2: at 1.0,
5 kb/s x 10 (T 30.72 ms, load 0.88) psr 0.78 / 0.73 / 0.75, on_time 0.29 / 0.26 / 0.26, frames left queued 2 / 5 / 0;
10 kb/s x 4 (T 15.36 ms, load 0.70) psr 0.956 / 0.969 / 0.938. At 0.5 these cells move to T 122.88 and 30.72 ms: psr
1.000 / 0.990 / 0.998 and 0.994 / 1.000 / 1.000, on_time 0.76 and 0.78 to 1.00, nothing left queued. The price: 802.15.4
carries at most half the scheduled load it did (10 kb/s x 10 is now infeasible at every period: at 5760 us per
116-octet frame the 10 sensors need 0.62 of any period). Not taken: 1.0, the admission (reachable through
`CSMA_MAX_LOAD`).

CCA busy: collision model: any transmission on the same channel active at the end of the CCA window;
SINR model: received power on the channel (adjacent channels weighted by the ACR) above the ED threshold. In both a
frame that starts at the end of the window is not heard.

### 4.3 Galvanic-optimisation ladder

BCC full carries every layer of the IEEE 802.15.6 HBC stack, several of which serve a shared, capacitively coupled
channel: other BANs' frames, collisions, a weak header. On the galvanic channel the bodies are assumed isolated
(PARAMETERS 5.7, A: no measured galvanic body-to-body number exists) and the link margin at 1 V is 26 to 56 dB
(PARAMETERS 3.8), so those layers cost airtime and protect against no other BAN. The proposal holds only under that
assumption (a touch between two bodies breaks it: level 4). Losses remain (the measured contact case, PARAMETERS 3.8:
PSR 0.66 to 0.72), so the FCS and the retries stay. The flag `opt` of `get_spec('bcc', 'full', 'sf8' | 'multitone')`
(`sweep.py --bcc-opt N`, variant label `full+optN`) switches the layers off one at a time, cumulatively; `opt` 0
(default) is BCC full as merged in /, row for row (test_coexistence checks it against the rows without the
flag). Level 1 uses 802.15.6 options, level 3 a hub policy within the standard (the pooled reference below holds it
from L0); 2, 4 and 5 are proposals for a galvanic amendment. Marks
as in PARAMETERS 1 (S secondary source, A assumed).

| level | what changes | 802.15.6 option or proposal | what isolation makes redundant | what still needs the layer |
|---|---|---|---|---|
| 1 | RI mode: the PHY header at the payload rate (SF 8) instead of SF 64: 2048 -> 256 us per frame and I-Ack at 500 kchip/s; a pilot every 256 symbols = 4096 Walsh chips (2048): 2 pilots instead of 4 in a 255-octet frame, 256 us | options: SFD / rate indicator and the pilot period '128 B' (PARAMETERS 3.7, S) | the header's 9 dB spreading gain: the link has 26 to 56 dB of margin and no interferer (A, 5.7). The pilot period is a clock question, not an isolation one: a 40 ppm drift is 0.16 chip over 4096 chips per clock, 0.33 chip between a sensor and the hub at 40 ppm each (A), a fraction of a chip between pilots | the contact case: in the degraded state (2 dB) the SF 64 header would survive, but the payload fails there anyway; a pilot per 4096 chips |
| 2 | no I-Ack or its pSIFS: the exchange is data + pSIFS; the hub acknowledges in its next beacon, per sensor a base sequence number (1 octet: every earlier frame of the sensor is held by the hub or given up; hub and sensor give a frame up at the same age of A = 4 superframes, the beacon numbering the superframes, so the base passes a frame the hub gave up only when the sensor has given it up too, formal/) and one bit for each of the n sequence numbers after it (n frames per sensor and superframe), N_dev + ceil(N_dev n / 8) octets; a failed frame goes again in the sensor's next allocation, which keeps the spare exchange per frame | proposal (A). The 802.15.6 B-Ack policy (PARAMETERS 3.4, S) answers with a B-Ack frame pSIFS after the frame that asks for it (Castalia BaselineBANMac times it as the I-Ack, S): with one frame per allocation it costs at least an I-Ack (1560 us with a 1-octet bitmap, against 1496 us, at 500 kchip/s); but at the admitted capacity the reference sends 15, 39 and 280 frames per allocation, so a reference with block acknowledgement (L-Ack, one B-Ack per block, PARAMETERS 3.7b) saves most of the I-Ack: against it the proposal's capacity ratio is 1.12 to 1.27 at 255 octets and at least 1.35 for 1 to 64 octets (claim4.py, T15 G, closed form: the 2-octet header and the beacon's ack field at one frame per superframe), against 1.306 to 1.363 and 1.987 against the I-Ack reference (T15 B). Precedent for the group acknowledgement in the beacon: the 802.15.4e DSME GACK, sent within the adjacent beacon interval (Kurunathan et al., IEEE 802.15.4e in a Nutshell, IEEE COMST 2018, S) | an immediate answer after every frame, the fast recovery from collisions and hidden nodes: on the isolated link losses are rare, one acknowledgement per superframe suffices | the acknowledgement itself (contact losses); a retry now waits one period |
| 3 | one exchange per frame (no spare); a pool of max(1, ceil(share x frames per superframe)) hub-polled exchanges, share = flag `pool` (default 0.05; 0.43 for the measured contact loss) (poll, pSIFS, the frame, pSIFS) after the allocations; failed frames polled in the order they failed, each at most once per pool, the others wait for the next pool; a sensor also resends the frames the beacon reported missing in the room its allocation has left; retry limit 3 kept in the engine, whose acknowledgement is an oracle (`_acked` reads the hub's record); the protocol gives a frame up at age A on both sides instead: a limit counted in attempts by each side lets the counts disagree after a missed poll and acknowledges a frame the hub never received after three losses (formal/, model checking) | option: improvised access, a poll grants a polled allocation starting pSIFS after it (Ullah 2.3.2, S); the poll as the empty frame, MAC header + FCS (Castalia BaselineBANMac, S), whose sequence-number field names the polled frame (A); a retry in the next allocation (3.4, S). Pool share 0.05 and 0.43 = 0.3 / 0.7 (the measured contact loss); both design points are reported. It builds on the level 2 proposal: the sensor learns its failures from the next beacon | the spare exchange per frame, room for a loss of every frame: on the isolated link losses are rare and independent across sensors, so a shared 5 % suffices | a pool sized for the loss: the measured contact case loses about a third of the tries, which the 5 % pool cannot carry (below) |
| 4 | lean framing: MAC header 2 octets (7): 1 bit frame type (data or poll), then a data frame's 15-bit sequence number or a poll's recipient NID (6 bits, mMaxBANSize 64) and spare bits; a poll carries the polled frame's 15-bit sequence number in a 2-octet body (L3: in its header), and the beacon's base sequence number per sensor (L2) takes 2 octets from L4. Under the age rule a sensor holds at most A + 1 = 5 periods of frames (the current one and A = 4 behind it; a buffer of 4 periods refuses a new frame after A lost beacons: formal/test_formal.py; the engine keeps its FIFO of 4 frames, OI-18), so 2^15 numbers keep up to 3276 frames per sensor and period unambiguous (the age rule with the hub's half-window duplicate test is sound for M >= 2n(A + 1) - 1 sequence numbers, n frames per sensor and superframe, 2A + 1 for n = 1; A = 4: 10 n - 1 <= 32768: formal/sweep_seq.out, sweep_n.out and the TLC runs in formal/tlc/); the admitted capacities of claim4.py reach at most 2020. The beacon body is not capped at the standard's pMaxFrameBodyLength of 255 octets: above 1888 frames of one sensor per superframe (460 each for four sensors) the ack field would need a second beacon frame, which is not modelled. The 1-octet header with a 7-bit number of earlier rounds held only 32 frames per period, below what the capacity cases send, and a resend could not be told from a new frame | proposal (A) | the BAN ID in every frame (1 octet): it tells apart frames of other BANs, which cannot arrive. The rest is a proposal for scheduled access, not a galvanic one: the TDMA slot names the sender and the hub is the recipient (2 octets) in any scheduled 802.15.6 BAN, capacitive HBC and NB too, and frame control 32 -> 16 bits (type and sequence number) is compression | the BAN ID: a touch between two bodies breaks the isolation, so the beacon keeps the 7-octet header with the BAN ID and a sensor that syncs to a foreign hub rejects it; the hub rejects a foreign sensor's frame at no airtime when the FCS starts from a value seeded with the BAN ID, as BLE seeds its CRC per connection (A, not modelled); the FCS (2 octets) stays |
| 5 | multi-tone: one preamble, SFD and header per stripe, on tone 0; the stripe's PSDU (MAC header, L, FCS) split so the stripe ends earliest, the other tones carry PSDU from the stripe start; the header's length field widened to 11 bits (a stripe carries up to 1534 PSDU octets) within the 8 reserved bits of the 32-bit header (PARAMETERS 3.8: rate 3, pilot 3, burst 1, seed 1, length 8, CRC 8), so the one header costs no airtime, and each tone's share follows from the length; one tone: as 4 | proposal (A: the tones share one transmitter clock and one galvanic path, so tone 0's sync gives the others their chip timing) | a preamble, SFD and header per tone: they synchronise tones that come from one transmitter | the one preamble, SFD and header (the stripe fails with tone 0); pilots on every tone |

Pooled reference stack. The spare exchange per frame of levels 0 to 2 is this model's provisioning rule
(section 3), not a feature of IEEE 802.15.6-2012, which sizes an allocation in slots requested by the node and granted
by the hub (5.7.2.5, 6.7.1), defers a retry to the next allocation interval when time is short and lets the
hub grant polled allocations for retries (6.6.1; Castalia's BaselineBANMac reserves no retry room and polls the free
slots). The flag `pool` below level 3 (`sweep.py --bcc-pool 0.05`, label `full+pool=0.05`) therefore builds the
reference stack with that hub policy: the level's own framing and exchange (below L2 with the I-Ack), one exchange per
frame in the allocation, and the pool of `pool_size` hub-polled exchanges after the allocations, each a poll, pSIFS and
the exchange (at L0 the poll is the L = 0 frame of 3264 us and a pool exchange 3264 + 75 + 23510 = 26849 us). The
engine runs it on the L2/L3 allocation path (`_late`, `_pool`) and `_acked` treats an I-Ack stack as informed at every
allocation (test_protocols and test_coexistence `test_bcc_pooled_reference`). Capacity per body at 500 kchip/s with the
pool at every level (share 0.05; the scan-free bound, the admitted 0.5 kb/s step in brackets), one tone: 82.1 (81.5),
98.2 (98), 106.6 (106.5), 106.6, 109.0 (108.5), 109.0 kb/s; six tones: 492.5 (491.5), 589.0 (588.5), 639.7 (639), 639.7,
654.1 (653), 690.2 (689) kb/s. Level 3 adds nothing to a pooled reference: L2 with the pool is L3. Against the pooled
reference the proposal (L5) gains 1.33 on one tone and 1.40 on six tones at 500 kchip/s, 1.34 at 2.625 Mchip/s (419 ->
561 kb/s) and 1.31 at 50 kchip/s (8 -> 10.5), against 2.52, 2.65, 2.53 and 2.6 with the spare-exchange reference: the
spare is exactly a factor 2 of those ratios (47020 -> 23510 us per 255-octet frame).

Airtime at 500 kchip/s (test_protocols.test_bcc_ladder, hand-derived; from L1 the RI mode's 608-chip SFD/RI field,
76 Walsh chips, 152 us, instead of the 64-chip SFD of L0, 802.15.6-2012 10.5.2 Eq. 117): frame of 255 octets on one
tone 20096, 18072 (L1-L3), 17688 us (L4-L5); exchange 23510, 19718, 18147 (L2, L3), 17763 us (L4); allocation per frame
47020, 39436, 36294, 18147, 17763 us; poll 1496 (L3), 1176 us (L4, with the sequence-number octet); pool exchange 19718
(L3), 19014 us (L4). Six tones, frames of 1530 octets: the same per tone, and at L5 a stripe of 16768 us (245 octets on
tone 0 after its 920 us head, 258, 258, 258, 257, 257 on the others), exchange 16843, poll 920, pool exchange 17838 us;
a 20-octet frame takes 920 us at L5 (its 23 PSDU octets fit tones 1 to 5 within tone 0's head). Capacity per body, 8 L_max /
(allocation per frame + pool exchange / 20 + 8 us of beacon bitmap from L2; the scan-free bound of select_period with
the pool and bitmap shares; the base sequence number, one octet per sensor and superframe, vanishes with many frames per
superframe):

| level | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| one tone (kb/s) | 43.39 | 51.73 | 56.20 | 106.58 | 108.96 | 108.96 |
| six tones (kb/s) | 260.3 | 310.4 | 337.2 | 639.5 | 653.8 | 689.9 |

RI mode (the header at SF 8, with the longer SFD/RI field) gives most of L1's gain, the pilot option the rest. The pool
rounds up per superframe: 20 frames of 1530 octets at L5 need 357613 us (684.5 kb/s), 21 frames a pool of 2.

The measured contact case (one tone, 50 kchip/s, the two-state chain: a try succeeds with p_s = 35 / 48 x exp(-frame
chips / 35000), 0.57 for a 251-octet frame of about 8900 chips) is where levels 3 to 5 pay: a superframe carries at most
its frames and the pool, so psr <= p_s (1 + pool / frames) (test_coexistence.test_bcc_ladder: N_dev 4, energy 1 kb/s: <=
0.71, while L0 to L2 reach >= 0.81). psr L0 / L2 / L3 / L5 (collision model, k 1, 3 seeds of max(2 s, 40 T)): N_dev 4
energy 1 kb/s 0.92 / 0.91 / 0.63 / 0.62, N_dev 4 latency 0.5 kb/s 0.94 / 0.95 / 0.84 / 0.78, N_dev 10 latency 0.1 kb/s
0.95 / 0.96 / 0.73 / 0.74; at N_dev 1 the levels are within noise. The resends in an allocation's idle room lift L3 / L5
there from 0.35 / 0.50, 0.83 / 0.68 and 0.54 / 0.59; the rest is the pool size: at 40 % of the frames L3 / L5 reach 0.80
/ 0.82, 0.96 / 0.89 and 0.94 / 0.93. So the 0.05 pool assumes the isolated, nearly lossless link. the ladder is reported at two design points of the pool share (flag `pool`, `sweep.py --bcc-pool`): 0.05, the link as
modelled, and 0.43 = p / (1 - p) for the measured per-try loss p = 0.3 (the pool then carries the expected retries of
a superframe), which trades capacity for contact-case psr: six tones 639.5 / 653.8 / 689.9 kb/s at L3 / L4 / L5 with
0.05, 459.6 / 471.7 / 499.2 with 0.43 (one tone 106.6 / 109.0 against 76.6 / 78.6). A pool the hub resizes from the loss
it measures, announced in the beacon, is the runtime form of the same rule (not modelled).

Slack-filling pool: the schedule is admitted on the pool
share alone, as above, so capacity is that of the share; at run time the pool also takes every pool exchange that fits
between the allocations' end and the next beacon, max(pool, floor((T - beacon - pSIFS - N_dev x allocation) / pool
exchange)) (test_protocols.test_bcc_pool_fill). The superframe's idle time thus carries retries. A loss estimate (the
runtime rule above) was planned and dropped (dev-flow 6.3 premise check): in this engine an unused pool exchange costs a
sensor only its receive guard, since polls stop when the failed frames are served, so a pool as long as the slack
costs next to nothing and the hub has no size to choose; an estimate would matter only for admission (how much share
to reserve), which the fixed period of a cell does not model. The standard leaves the pool policy open (PARAMETERS 3.7b:
polled allocations exist, no pool rule, no HBC retry limit).

Touching bodies and the fallback. The BANs are paired
(0-1, 2-3, ..); each pair touches during alternating exponential intervals, mean contact 5 s (A) and the share S of time
in contact, drawn from their own random substream. While a pair touches, a frame of the partner BAN reaches a receiver at
the partner node's own-link Ec/N0 less X dB (default 10; the contact case: the chain's connected level less X), with no
sync lock across BANs. No loss figure for touching bodies exists, so X is a sweep. Two BANs with the same period keep the offset of their superframes, so a pair
whose active periods overlap once overlaps in every superframe: its beacons fail, the sensors that miss a beacon do not
send, and packet success collapses for as long as the contact lasts (test_coexistence.test_bcc_touch_and_fallback: psr <
0.1 at X = 0 when the superframes start 222 us apart). The fallback is the hub's reaction: on hearing a foreign beacon it
interleaves the pair's active periods (802.15.6 6.13.3, active superframe interleaving; PARAMETERS 3.7b), from one
superframe after the contact starts, when both active periods (beacon, allocations and the admitted pool, the slack that
`fill` takes being given up) fit in one period. ponytail: the interleaving is modelled as "no partner frame reaches the
receiver after detection", not as a moved schedule; a real shift costs one superframe of latency. Foreign frames are
never delivered as the BAN's own: the engine attributes every frame to its BAN, which is what an FCS seeded with the BAN
ID gives (REPORT 5.4); an unseeded FCS could accept a foreign frame, which is not modelled.

Payload spreading and the rate rule.
'walsh' (default) is the HBC PSDU: 2 Walsh-16 symbols per octet, 32 chips. 'none' sends the PSDU unspread, one BFSK
chip per bit, 8 chips per octet, with the chip error 0.5 exp(-g / 2) of `phy.ber('bfsk_chip')` and no decoding: all
bits must be right (the FCS finds a wrong one). Preamble, SFD, header and pilots stay as they are (the header carries
the mode, like the RI field), and a pilot still follows every 128 PSDU octets. The PSDU of a 255-octet frame shrinks
from 8256 to 2064 chips, so the frame at 500 kchip/s shrinks from 17688 to 5304 us
(test_protocols.test_bcc_unspread_and_rate_rule). 'adapt' (from L2) is the per-link switch. A sensor sends unspread,
and turns to Walsh after ADAPT_DOWN = 2 consecutive frames that the beacon reports lost. It turns back to unspread
after ADAPT_UP = 16 consecutive frames reported delivered (A: the hysteresis keeps one bad contact dwell from
toggling it). The sensor learns the outcomes from the next beacon it hears, one superframe late, and the hub keeps the
same count. The allocations are sized for Walsh frames, so a switched link still fits: 'adapt' has the capacity of
'walsh'; an unspread frame only ends earlier (test_coexistence.test_bcc_spread_switch, with an injected report
sequence).
Rate rule (`protocols.margin_rate`): the highest chip rate at which the 95th-percentile participant on the weakest
path of the link budget (RHW1-LHW1, A_95 = 15.7 mV at 500 kchip/s and NF 20 dB for Ec/N0 8 dB) still
meets the mode's requirement at the drive A and noise figure NF, capped at BCC_RC_MAX_CPS. The rate is
R_c = 500e3 (A / 15.7 mV)^2 10^((20 - NF) / 10) 10^((8 - E_req) / 10). E_req is 8 dB for Walsh. For unspread it is the
Ec/N0 at which the 255-octet PSDU is delivered as often as the Walsh PSDU at 8 dB (0.955): 13.02 dB
(`required_ec_n0_db`). At the 1 V drive and NF 30 dB the rule gives 64 Mchip/s unspread and 203 Mchip/s Walsh. The
population margin therefore never binds below the waveform's 9 Mchip/s, and the rule selects unspread at every
offered rate. It binds from a drive of 0.375 V (unspread, NF 30, 9 Mchip/s).

Scheduling (engine): from L2 the beacon is `protocols.beacon_layout` (17-octet body, A, plus the ack field; from L4 with
the 7-octet MAC header, which names the hub and the BAN) and is decided at every sensor like a frame, on its one-tone
layout (multi-tone: on tone 0). At its allocation, after the beacon, a sensor drops the frames the hub holds when it
heard the beacon (by sequence number, so a duplicate is acknowledged too; ponytail: the model acknowledges every frame
the hub holds, also one more than n sequence numbers past the sensor's base, which the field cannot name; that takes
consecutive losses) and the frames past the retry limit (`drop_retry`); a sensor that missed the beacon keeps its frames
and, at L2, resends them (duplicates of delivered frames count as attempts); A: it keeps its slot timing, as at L0. L2
allocation: every unacknowledged frame once, oldest first, while its exchange fits the slot of 2 exchanges per frame
(the previous superframe's failed frames go first, into the room of the spare exchanges). From L3 the slots hold one
exchange per frame, packed from the beacon on (`slot_us`), the pool follows the last slot (`pool`, `pool_us`), the rest
of the period is idle; a sensor sends first the frames it has not sent, then, when it heard the beacon, the frames the
beacon reported missing while room is left (only when arrivals were dropped at a full queue); a frame that fails for the
first time joins the hub's list for the pool (ordered by failure: sensor order within a superframe, earlier superframes
first). A poll the sensor misses costs the frame no attempt: it keeps its place in the queue and the list. Pool energy:
the hub wakes when it polls; every sensor wakes for the pool and listens from the receive guard, poll after poll, the
guard alone when no poll comes (A: before the next beacon a sensor does not know whether its frame failed; the 802.15.6
future poll, a B-Ack or poll that tells a sensor when it will be polled (Castalia handlePoll, S), would wake only the
polled sensors, at the cost of that frame: not modelled); the hub waits the frame's airtime after a poll the sensor
missed (A, as a sensor for a missing I-Ack). Not taken: resending only the failed tone packets of a stripe (4.1).

## 5 Reception models

Every transmission has start, end, channel, transmitter, intended receiver, power. The outcome of a
transmission is decided at its end, when every overlapping transmission is known. A receiver that
transmits during any part of the packet loses it (half duplex).

### 5.1 `collision` (rule)

A packet fails if any other transmission on the same channel overlaps it in time by any amount. A packet that survives
is lost with probability `p_loss` (PHY noise): 1e-4 for 1M-class PHYs, 1e-3 for 2M, and the value in `PARAMETERS.md` for
the others. No capture, adjacent channels ignored. FEC and repetition have no effect in this model (they neither prevent
a collision nor lower `p_loss`, which already represents a link with margin), except the SmartBAN whole-PPDU repetition,
whose copies are decided alone (4.1; PARAMETERS 3.5); ARQ, CSMA, adaptive channel selection and hopping do. BER in this
model: bits under an overlap count with error probability 0.5, other bits with the BER that gives `p_loss` for the
frame.

### 5.2 `sinr` (physical)

Geometry per trial: the room generator of PARAMETERS 5.3 (square room of area k x a_p, body centres uniform
with 0.6 m spacing, random facing, torso cylinders for blocking, node positions on the body frame); hub on
the chest and sensors on the listed body sites. Path loss: on-body (sensor <-> own hub) drawn per link from the
on-body model; body-to-body (any node of body j -> any node of body l) = log-distance model at the body
distance plus log-normal shadowing per (transmitter, receiver) pair, drawn once per trial. Transmit power
per protocol (default 0 dBm). Noise N = -174 dBm/Hz + 10 log10(B_rx) + NF.

For a packet p at receiver r: split [start, end) at every start and end of overlapping transmissions q on
any channel. Per segment SINR = S / (N + sum_q I_q x 10^(-ACR(df)/10)), df the centre-frequency offset of the two
channels (`protocols.acr_matrix`; BLE data indices 10 and 11 are 4 MHz apart); the last entry of an ACR table holds
for every larger offset (PARAMETERS 5.5 floors such as '>= 3 MHz'), except 802.15.6 NB, which ignores |dnc| >= 2.
- Co-channel floor: p also fails when the signal-to-interference ratio of a segment falls below the floor of a block the
  segment overlaps, because treating interference as Gaussian noise overstates capture; the block's decisions in that
  segment are guesses. The floor of a block is the receiver's C/I_co (PARAMETERS 5.5) lowered by the block's processing
  and coding gain, calibrated on the BLE specification's own coded C/I values (1M 21, S=2 17,
  S=8 12 dB; PARAMETERS 3.1: the co-channel gains are 4 and 9 dB where the sensitivity gains are 5 and 12 dB), A: a
  repetition or spreading factor S buys 10 log10(S) x 5 / 6.02 dB (BLE S=2 -> S=8 is the P = 4 mapper alone, 17 -> 12
  dB), a BCH block 4 / 5 of its AWGN gain (the S=2 code buys 4 of its 5 dB), 0.8 x 1.4 = 1.12 dB (the ETSI BCH gain of
  3.5 at the 1 Msym/s channel rate, Table 13: Eb/N0 8.8 against 7.4 dB; applied to the 802.15.6 NB codes as the
  analogue). S is the 802.15.6 NB S = 4 on header and PSDU and, for HBC, the SF over SF 8, which is spectral shaping
  only (3.7): the SF 16 / 32 / 64 modes and the SF 64 header gain 2.5 / 5.0 / 7.5 dB. BLE Coded takes the
  specification's own C/I (17 dB S=2, 12 dB S=8) for the whole packet, which holds its coding. Floors: NB 121.4 kb/s
  header and PSDU 14.88 dB, SmartBAN coded header and PSDU 19.88 dB (uncoded PSDU 21 dB; Nrep 2 and 4: every block a
  further 2.5 and 5.0 dB lower, preamble and sync included, as they are repeated too), HBC SF 64 header 13.5 dB. The
  first credit (the full 10 log10(S) and the BCH gain 1.9 dB at the information rate: 13.08, 19.1, 11.97 dB) was
  1 to 2 dB more lenient; the review of found the same psr and on_time to 4 decimals with both (HBC SF 8 and NB,
  10 kb/s latency, N_dev 4, k 8 and 30, SINR, 5 trials). BCC full gets no credit (`ci_credit` False): its SF 64 header
  chip integrates 8 R_c chips, and a same-tone interferer whose 8 chips all sit on the other tone still wins below 0 dB
  (at -9 dB, 5 of 8 on the other tone suffice, probability 93/256), so the 0 dB chip rule of 3.8 holds for every block
  of both BCC stacks (no BCC case puts an interferer within 40 dB of the own link anyway). Blocks without repetition,
  spreading or BCH keep C/I_co: CRC fields outside the BCH, the HBC SF 8 payload, and every preamble and sync block (the
  802.15.6 NB m-sequence preamble at S = 1, the HBC preamble, SFD and pilots at SF 8, BLE 1M / 2M, SmartBAN, ANT): a
  frame whose preamble is overlapped is judged at C/I_co, so the gain acts only on overlaps that miss the preamble
  (pessimistic for NB and HBC, whose preamble correlators are not modelled; BLE Coded's specification value covers its
  whole packet, preamble included). C/I_co itself is on a mixed basis (PARAMETERS 5.5, OI-3): the BLE specification's 21
  dB for BLE 1M and 2M and as the analogue for SmartBAN, NB and HBC, the chip values for 802.15.4 (3 dB) and ANT
  (nRF24AP2, 12 dB), so uncoded ANT GFSK clears a floor 9 dB below that of BLE 1M and of BCH-coded SmartBAN at the same
  1 Msym/s. ANT full, 1 kb/s x 4 latency, SINR, the sweep seeds of trials 0 to 2, every BAN on 2457 MHz: psr 0.189 at 12 dB against 0.122 at 21 dB (k 8), 0.035 against 0.011 (k 30). To decide
 ; kept as
  PARAMETERS gives it. Default taken 2026-09-24: the floor stays the default (PARAMETERS OI-3);
  `simulate_trial(..., ci_floor=False)` and `sweep.py --ci-floor off` drop it for every stack (interference as Gaussian
  noise, variant label '+nocifloor'), a sensitivity run in which FEC can show its effect against weak interference. The
  collision model has no floor, so its rows keep the default label (a rerun finds them done instead of repeating them).
- Fading: on-body and body-to-body gains carry a static per-link shadowing drawn per trial and a per-packet
  fade (PARAMETERS 5.1, 5.2), redrawn for a retry only after the coherence time.
- Sync: the synchronisation part of p (preamble + access address / SFD) must see SINR >= `sync_sinr_db` in all
  its segments. Protocols with a common preamble (`sync_on_foreign = True`: 802.15.4, 802.15.6, SmartBAN, ANT)
  also lose p when r locked onto an earlier foreign frame on the same channel that is still in the air at
  p's start and whose own SINR at r at its start was above the lock threshold, provided r was listening for p when
  that frame began: its listening for p had started (DESIGN 6: the receive guard before a first reception, p's
  start for a known response, the TSCH template offsets, the CSMA sensor's data end; the CSMA coordinator always
  listens) and r was not transmitting. BLE (access address) and TSCH (addressed, scheduled) do not lock onto
  foreign frames. Lock threshold (OI-6: 'the receiver locks on the first packet above sensitivity'):
  `sync_sinr_db` where the recipe gives one (802.15.4: -5 dB), else the sensitivity (802.15.6 NB FULL, below),
  else the SINR of a 1e-3 bit or chip error, the BLE sensitivity criterion (10.94 dB on the GFSK curve, which is the
  calibrated sensitivity of the SmartBAN and ANT chips; the same rule for the HBC and BCC chips, A).
- Sensitivity: a frame whose received power is below the protocol's sensitivity is lost (reason `noise`,
  PARAMETERS 5.6, 3.4 FULL). Only 802.15.6 NB FULL carries one, the TI draft value of its rate (-95 / -93 / -90 /
  -86 dBm): its BER model (post-detection combining, BCH) alone would put the 121.4 kb/s threshold about 3 dB lower.
  The other RF stacks' BER models are calibrated to their chip's sensitivity (PARAMETERS 5.4), so a floor would
  repeat them; BCC and HBC have none.
- Bits: per segment the bit (or chip) error probability of the modulation at the segment SINR
  (`phy.bit_error_prob(mod, sinr_lin, ...)`, with repetition or spreading factor S as SINR x S). The frame
  is cut into its coded blocks (header block and payload codewords of the protocol's FEC, or single bits when
  uncoded); the number of errors in a block is the sum over segments of Binomial(bits of the block in the
  segment, BER of the segment), drawn with the trial RNG. Block code (n, k, t): the block fails when errors
  > t. Convolutional code (BLE coded): SINR shift by the coding gain of `PARAMETERS.md`, then uncoded
  decision per bit (stated approximation). CRC: a frame with any residual error is discarded when the
  protocol has a payload CRC; the BCC stack has none, so its frame is delivered and counted
  `delivered_corrupt`.
- BCC and HBC use the same machinery with the 16-chip Walsh symbol as block (n = 16 chips, k = 4 bits,
  t = 3) and the chip error probability of non-coherent BFSK, P = 0.5 exp(-Ec / (2 N0)). The
  sync part and every pilot are 64-chip correlators (<= 7 chip errors); a missed pilot shifts the rest of the frame
  (PARAMETERS 3.8), so every payload bit after it counts as a guess (error probability 0.5) in raw and post BER.
  They see no inter-body interference (isolation assumption). Default taken 2026-09-24: the +40 / +60 dB BCC coupling cases are an SINR concept (the interferer 40 to 60 dB weaker), which the
  collision model would turn into 'any overlap kills', so `simulate_trial` raises ValueError and `sweep.py` refuses
  `--bcc-coupling` with `--reception collision`. Their S / N comes from the
  measured path gain, drive amplitude, chip rate and
  N0 = 4 k T R NF. One tone: its centre is `channel_curves.tone_centre(R_c)`, the
  multiple of 1 kHz that maximises the weaker of the two BFSK tone gains g(f - R_c), g(f + R_c) of the cohort 50 Ohm
  curve with both tones inside the measured 0.05 to 20 MHz (12.517 MHz at 50 kchip/s, 12.680 at 500 k, 13.230 at 2 M,
  14.750 at 5.25 M, 12.000 at 7 M and 8 M, 11.000 at 9 M; refused above 9 M, see 4.1); a link takes its drawn path's
  gains at f - R_c and f + R_c, and its chips decide at the one Ec/N0 whose error 0.5 exp(-E / 2) is the mean of the
  two tones' (`channel.bfsk_ec_n0_db`, between E_min and E_min + 2 ln 2), the rule the chip-rate study checks against
  the waveform model (`results/bcc_rate_stability.md`, test_coexistence.test_bcc_rates_match_study). Before the
  tone sat at 8 MHz, in the dip, with the path gain at its centre: at 500 kchip/s the weakest link (RHW1-LHW1) moves
  from 47.4 to 52.4 dB. The BCC rows at 500 kchip/s and those of the 50 kchip/s contact case stay bit-identical (976
  rows of BCC full, the as-built stack, multi-tone, the coupling cases, the contact case and ladder levels 3 and 5,
  collision and SINR, k 4 to 30, two seeds): every link stays tens of dB above the error-free point, the coupled
  interferer 40 / 60 dB below its own link, and the contact chain carries its own Ec/N0. One gain per BFSK tone: the
  gain slope across a tone and the inter-chip interference it causes are not simulated in the engine, and from 5.25
  Mchip/s on the band f -/+ 2 R_c leaves the measured range (the results file lists what this misses). The
  repository holds no phase of the measured paths.
  checked the dispersion with the minimum-phase response of each measured magnitude, which is the least
  smear the body can cause. The waveform was filtered in real passband at 0.5 to 9.5 Mchip/s, NF 10 to 30, every
  path, with ideal chip timing (`bcc_fdma.py isi`, `results/bcc_rate_isi.csv`, section "Dispersion" of
  `results/bcc_rate_stability.md`).
  - On the body-only curves every path keeps PSR 1 up to 9.5 Mchip/s. These curves have the arrangement's 8 to 9
    MHz notch, 12.5 to 15 MHz peak and 20 MHz drop smoothed out. Their least smear is about 20 ns rms, with 80 to
    110 ns holding 90 % of the energy.
  - So smear sets no bound below the image limit, and `BCC_RC_MAX_CPS` stays 9e6.
  - The curves as measured fail at 9.5 M (up to 26 paths below 0.99), and at 9 M a few chips err at NF 30.
  - At 9 M the engine's one-chip rule is optimistic by up to 1 dB near the error knee. At 1 V the weakest path stays
    19 dB above that knee at NF 30.
  - Every result above 2.625 Mchip/s is an upper bound: it holds for this receiver without an equaliser, and
    for the measured magnitude to 20 MHz with the stated extension. Excess phase would only lower it. Option
  `--bcc-contact` (the measured-hardware case, at the measured 50 kchip/s only): the two-state contact model of Section
  5.8 as a time-varying gain per link (connected / degraded dwell exponential with means 35 000 and 13 000 chips, depth
  12 dB, Eb/N0 14 dB when connected), continuous per link and used by both reception models (BCC is isolated, so the
  model choice must not matter); a failure after a good sync and header is delivered corrupt. The stationary chain gives
  35/48 exp(-N / 35 000) = 0.573 for a 250-octet frame (8640 chips, of which N = 8448 are decoded: the first three
  preamble copies are not), below the measured 0.66-0.72 on the arm: the as-built runs start connected (0.78 for a first
  packet, about 0.67 over its 3-packet runs). Default taken 2026-09-24: not recalibrated; the
  report prints the chain result (0.573) next to the measured 0.66-0.72.
- Multi-tone BCC full: the path drawn per link (the same draw as for one tone) gives the gains at the BFSK
  tones f_i -/+ R_c of every tone (`channel.bcc_link` with the frequency list, log-f interpolation), combined per tone
  by `channel.bfsk_ec_n0_db` as for one tone, and the drive is split by sqrt(N), so each tone packet has its own Ec/N0: at 1 V, 6 tones, the weakest tone of the weakest path
  (RHW1-LHW1 at 8.24 MHz: -52.72 and -51.82 dB at 7.74 and 8.74 MHz) has 40.46 dB. `Sim._eval_tones` decides every tone packet with
  `phy.frame_outcome`; the stripe succeeds only when all do, and raw and post BER add up over the tones. No inter-tone
  interference and none inside the isolated BAN (its TDMA slots never overlap; ponytail: an
  overlap would fail the whole stripe). The contact chain and the coupling cases would need a per-tone model (new
  physics): `get_spec` refuses them for multi-tone and `sweep.py` skips that entry.

HBC FULL between bodies: both models put HBC far below the retired fair variant. Its frames are
over ten times longer (SF 8 data after the 164 kb/s header: 372 us for 10 octets and a 415 us beacon, against 26 and 31
us at SF 1; SF 16 to 64 longer still), so on the single 21 MHz channel almost every exchange overlaps a frame of another
body. The collision model gives no credit for spreading; in the SINR model the C/I_co floor (21 dB, the BLE value as
analogue, OI-3) holds for the SF 8 payload, preamble, SFD and pilots, while a body 1 m away couples only 15 dB below the own link (typical case, PARAMETERS 5.7), so such an overlap
fails, sync included., 10 kb/s latency, 3 trials, full / fair: collision, N_dev 1, k 8: psr 0.042 / 0.958,
collision_frac 0.99 / 0.04; SINR, N_dev 4, k 8: 0.678 / 1.000, and 1.000 / 1.000 with `--ci-floor off`.

## 6 Energy

Energy basis : every stack at its own figures, the radio stacks on their chips' system currents
(nRF52840 compounded, PARAMETERS 4.2), BCC and HBC on transceiver-IC figures; no nRF system overhead is added to BCC or
HBC.

Every node accumulates time in TX, RX (receiving, listening for a poll, ACK, beacon, CCA), ramp-up (count of wake-ups x
ramp time) and sleep. Energy = V (I_tx t_tx + I_rx t_rx + I_ramp t_ramp + I_sleep t_sleep) with the radio of
`PARAMETERS.md` (default chip per protocol; BCC low / typical / high from the transceiver literature, as power per
state; multi-tone BCC: N parallel front ends, each at the single-tone power in every radio state, so TX, RX and ramp
currents x N for the data packets, the I-Acks and their receive guards, the one-tone beacon on one front end
(`spec.beacon['energy']`), the MCU wake-up and the sleep current once; ponytail: conservative, a shared wideband front
end may draw less and a frame with fewer octets than tones still powers all N; the `--bcc-power` cases scale alike).
Default taken 2026-09-24: the beacon on one front end, the consistent reading of a
one-tone packet. Booked at N front ends it was a third of the sensor energy and made multi-tone J per bit 26 to 28 %
higher at every R in energy mode and 33 to 35 % in latency mode (collision, N_dev 1, 0.02 to 30 kb/s); multi-tone
against single-tone J per bit is now 1.9 to 2.1 (energy) and 3.2 to 4.3 (latency), was 2.5 to 3.0 and 4.7 to 6.4). Radio
only: MCU and sensor energy are excluded. Receive time per frame: the first reception of a scheduled event (beacon, BLE
poll at the anchor, the first frame of a slot, an ANT message) starts the receive guard early (window widening of
PARAMETERS 4.1: 4 x drift x T + 2 x 16 us); a response whose timing the receiver knows (ACK, the BLE T_IFS turn, later
exchanges of a slot) has no guard; TSCH follows its timeslot template (hub from macTsRxOffset, macTsRxWait when nothing
is sent; sensor from macTsRxAckDelay, macTsAckWait when no Enh-Ack comes); in the TSCH FULL shared cell, the RFC 8180
minimal cell (TX | RX | SHARED for every node), every node of the BAN that sends nothing there listens and wakes for it:
from macTsRxOffset to the end of the longest own-BAN frame that starts there (the EB for the sensors, a sensor's frame
for the hub and the other sensors), macTsRxWait when none does; the SmartBAN FULL hub opens its receive guard at both C/M slots of every IBI, the longest
frame there adds its airtime, and a sensor that sends there wakes for it. In such a contention cell or slot a node
listens once, however many frames collide there (booking each colliding frame counted 1.5 to 3.6 % of the TSCH
hub's RX time twice at k 8 and 30, up to 0.45 % of SmartBAN's); ponytail: a frame of another BAN in the cell is not
received. The 802.15.4 CSMA sensor listens until the ACK ends or macAckWaitDuration (864 us) runs out; the ANT master
listens for its RX window W (176 us) and, with adaptive isochronous re-timing, on through its 500 us guard window after
the exchange: 644 us of RX with one start-up per message, 500 us with its own start-up after an old-payload broadcast
(4.1; the guard size is assumed, A in PARAMETERS 3.6 FULL step 4, and it raises ANT's J per bit by about 73 % at k 1:
728 against 421 nJ at 5 kb/s latency, N_dev 4); the non-beacon coordinator all the time. The ANT master's listen at
channel open (one channel period tc of RX, PARAMETERS 3.6 FULL step 1) happens once per deployment, not per run, and is
left out, like the hub's search as a slave: booked in each run of 40 periods it was a fixed
2 % RX duty (tc = T / 1.25 in 40 T, about 125 uA) that made 88 % of the sensor energy at 0.02 kb/s (J per bit 2.0e-5
against 2.5e-6, lifetime 58 d against 476 d); spread over the battery it is at most 6.26 mA x 4 s = 7e-3 mAh, 3e-5 of
the 211.5 mAh. Outputs: energy per sensor and per hub (J), average power (uW), energy per delivered payload bit for the
sensor alone and for the whole BAN (sensors + hub), and sensor lifetime on the battery of `PARAMETERS.md` (4.4, with the
buck-converter factor 3.0 / V_cell taken at 2.8 V over the whole life).

## 7 Metrics per BAN and trial (CSV columns)

protocol, variant (full or thesis, plus '+flag' for a non-default flag and '+nocifloor'), phy, reception, strip_hopping
, ci_mode, period_us, r_target_kbps, r_actual_kbps,
r_agg_kbps, n_dev, k, payload_bits, frames_per_period, duration_us, trial, seed, network_id, feasible,
generated (every frame whose first access opportunity falls before the end of the run: its slot, connection
event or ANT channel period; CSMA its sample instant), tx_attempts, delivered_ok, delivered_corrupt, drop_collision,
drop_noise, drop_sync, drop_csma, drop_retry, drop_queue, drop_end (still queued at the end, waiting for a retry or
behind a backlog: counted as not delivered, so the outcome counts add up to generated and a backlog cannot hide in
the denominator; the sweep runs at least 40 periods, so frames awaiting a retry at the end bias psr down by about
per_first / 40 at most), psr (= delivered_ok / generated), pdr (= (delivered_ok + delivered_corrupt) / generated),
per_first (failed first attempts / first attempts), retx_per_frame, lat_mean_ms, lat_p50_ms, lat_p95_ms,
lat_p99_ms, lat_max_ms, on_time (delivered_ok within the deadline D = the cell's reference period P of
section 3, PARAMETERS OI-4, over generated), raw_ber (bit errors before FEC over the payload bits on air: coded bits
of BCH codewords, bits after DSSS despreading (802.15.4) or Walsh demapping (BCC, HBC: no FEC there, so raw = post);
`phy.decisions`),
post_ber (after FEC, over transmitted payload bits), goodput_kbps (delivered_ok application bits / duration: each
frame carries its share of L, the last frame of a period the rest; the padding to whole octets is airtime only),
airtime_frac (BAN airtime / duration), e_sensor_j, e_hub_j, p_sensor_uw, jpb_sensor (J per delivered application
bit; inf when the BAN delivered nothing), jpb_ban, lifetime_days,
collision_frac (attempts that overlapped another same-channel transmission), cca_busy_frac.

Per trial (all BANs): jain_fairness of psr over the k BANs, channel_occupancy (mean over channels of the
busy time fraction), written as columns repeated on each BAN row of that trial.

Cells (`sweep.py --cells`): per metric the mean over trials of a per-trial value and a Student-t 95 % interval over
trials, clipped to the physical range ([0, 1] for psr, pdr, on_time, the BERs, collision_frac and jain_fairness, >= 0
for the others). The per-trial value is the mean over the k BANs, except J per bit, a ratio of sums over the k BANs
(sensor energy, plus the hub's for jpb_ban, over delivered application bits): a BAN that delivered nothing counts, and
a trial that delivered nothing gives inf (the cell mean is then inf, without an interval). A mean of per-BAN ratios
would drop the BANs that delivered nothing and follow those with a handful of deliveries.

## 8 Validation (test_coexistence.py)

1. `ble4thesis` in the collision model reproduces `BLE_SIM.simulate_one_trial` mean PSR within 0.01 on six or
   more (R, k, N_dev) points, both PHYs.
2. Periodic random-phase transmitters on one channel (ble4thesis held on one channel, collision model, no p_loss):
   PSR within 0.01 of the analytic (1 - 2 tau / T)^(M - 1) for M transmitters with airtime tau.
3. CSMA single sensor, no contention: PSR 1 and mean latency within 2 % of the hand-derived value (mean
   initial backoff (2^macMinBE - 1) / 2 x aUnitBackoffPeriod + CCA + turnaround + frame airtime).
4. BER functions at hand-derived points; BCH block failure against the binomial tail; Walsh symbol with
   t = 3 at chip error 1.5 % gives a packet error-free probability consistent with.
5. Energy of one isolated BLE full exchange and of a TSCH slotframe against the hand sum of currents x times.
6. Full capability: the retry limits, channel rules and rates of every catalogue stack, and the variants
   `get_spec` builds (`test_full_capability` in test_protocols); TSCH blacklisting removes a jammed channel at its 21st
   failure in 100 and stops at 4 channels; an ANT listen that finds no free phase takes the least-occupied one (0, tc /
   2, tc / 4, 3 tc / 4 for four masters whose tries all hit phase 0); SmartBAN C/M retries by the MCA backoff (the retry
   5000 / 6250 us after the failed frame at BC 0 / 1, BE 3, 4, 5, 5 on a busy channel; two sensors with BC 0 collide in
   both reception models); the 802.15.4 GTS selections by hand and a failed frame retried one exchange later in its GTS;
   the SmartBAN PHY schemes' airtimes and LF,max by hand, and in the collision model a frame with one clear copy arrives
   (`test_repetition_copies`); ANT advanced burst selections by hand (20, 50 and 60 kb/s, none above), and frequency
   agility moves a master jammed on its first frequency to the second after its 8th unacknowledged message.
7. Same seed gives the same rows; different trials give different rows.
8. Retry room: a failed BLE data PDU goes again one exchange later in the same event, a TSCH
   frame in the next shared cell (two sensors collide there, and the hub listens there once), the shared-cell backoff
   skips 1, 3, 7, 15 shared cells when every draw is at its maximum, a SmartBAN frame goes in the first C/M slot at BC 0
   (the hub listens once per C/M slot), an 802.15.4 GTS frame one exchange later in its GTS. Alone (k = 1, p_loss 0)
   with one forced loss no full catalogue spec (nor BCC full at ladder levels 1 to 5, one tone and multi-tone) leaves a
   frame queued at the end or dropped at arrival over 40 periods, in cells that include those the exactly provisioned
   schedules and the CSMA admission could not clear; except ANT at N_dev 10, whose masters still overlap
   each other at the least-occupied phases (4.1; 5 kb/s, k 1, p_loss 0, seeds 0 to 2: drop_queue 54 of about 400 frames
   each, 220, 140 and 164 with the random-phase fallback), and ANT at 10 kb/s latency x 4, whose advanced bursts fill
   the period (3; seed 0: psr 0.500). The room is sized for isolated losses: under a sustained loss it can run out. TSCH
   FULL (before its CCA) at 10 kb/s energy (3 frames per sensor per 200 ms slotframe, N_dev 4, k 1, collision, p_loss
   0.02, 200 periods, seeds 0 to 4) drops 51 frames at arrival and leaves 65 queued of 11994; SmartBAN FULL in the same
   cell leaves 4 of 4000 queued (awaiting their retry at the end).
9. BCC full ladder: frame, exchange, allocation, poll and
   pool airtimes, beacon sizes, selections and capacity per level at 500 kchip/s, hand-derived; the L3 pool polls the
   failures in their order, within the superframe and in the next; a sensor resends a reported failure in its
   allocation's idle room; the L2 retry in the next allocation; the beacon as the ack (after a missed beacon L2 resends
   a duplicate, L3 waits); the energy of one L2 exchange and of the L3 pool with and without a poll; psr 1 at every
   level in both models at k 1 and 30; the contact case at every level (N_dev 1: psr > 0.5; N_dev 4: L0 to L2 >= 0.81,
   from L3 <= 0.71, the cap of the 5 % pool); opt 0 prints the rows of no flag.
10. One-tone BCC against the chip-rate study: at 500 kchip/s, 7 and 9 Mchip/s
   every link the engine draws has the one-chip Ec/N0 of a study path, the study makes no chip error at NF 10 to 30 and
   the engine delivers every frame of both single-tone stacks; 9.5 M (chip errors on every path at NF 20) and 9.95 M
   (every path below PSR 0.99 at every NF) are refused. Near the knee (Eb/N0 7 dB) the frame model plus the 0.35 dB of
   the PHY2 receiver (`phy.V2_OFFSET_DB`) holds every path's waveform PSR within 3 binomial sd at 500 k and 2.625 M;
   at 7 M the image of the 5 MHz lower tone costs the waveform up to 0.5 dB, so the model is optimistic there by at
   most that much (`results/bcc_rate_stability.md`). Under dispersion, at
   every tested rate up to the cap, every path of the body-only minimum-phase curves keeps PSR 1. At 9.5 M the curves
   as measured fail, while the model without dispersion does not.

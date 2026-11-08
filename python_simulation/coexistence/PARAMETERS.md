# Parameter record of the BAN coexistence comparison

Ticket (epic E-L04). Written 2026-09-23. This record feeds (simulator) and.

## 1. Purpose, marks and use

This file lists every constant the coexistence simulator uses: frame formats, timing, PHY rates and modulation, channel plans, FEC and repetition, ARQ, channel access, coexistence features, BER models, radio current per state and the propagation model. It covers BLE 5, IEEE 802.15.4 (unslotted CSMA/CA and beacon/GTS), IEEE 802.15.4 TSCH, IEEE 802.15.6 NB, ETSI SmartBAN, ANT (with Nordic ESB as an optional comparator), the IEEE 802.15.6 HBC standard stack and the as-built BCC stack of `python_simulation/HBCProtocol.py`. Each protocol runs in two variants, FULL and FAIR (section 2).

Every value has a source and a mark:

| mark | meaning |
|---|---|
| V (verified) | Read in the primary source (standard, datasheet, paper) or in a verbatim quote of it. Plain arithmetic from V field sizes (for example an airtime formula whose endpoints reproduce the spec) is also V, with "derived" in the note. |
| S (secondary) | Read in a secondary source: review, textbook, patent of a draft proposal, ns-3, Contiki-NG or Zephyr code, application note, forum, Wikipedia. |
| A (assumed) | No source found, or a model output (fit, union bound, numerical evaluation), or a value that only one draft-era or secondary source gives and the adversarial check could not confirm. The note gives the reasoning or names that single source. |

Every research result was checked by an adversarial verifier. Where the verifier corrected a value with evidence, the corrected value is used and the note says "corrected". Where the verifier could not confirm a value, it keeps its value and drops one mark (V to S when a secondary source shows it, otherwise to A; S to A unless the verifier re-read the same secondary source and found the value there).

Sources are numbered; `[n] locator` in a table points to section 7, which also names the `biblio.bib` key when one exists.

How the simulator uses this file:
- reads the constants from `parameters.json` (the machine copy written with this record; if the two differ, this file wins and the JSON is corrected).
- Symbols: L = payload octets, k = bodies (BANs) in the room, N_dev = sensors per BAN, CI = BLE connection interval, Rc = BCC chip rate, γ = linear SNR or SINR per symbol unless stated.
- Units: times in µs unless stated, powers in dBm, gains and losses in dB.
- Two interference modes, chosen per run and identical for FULL and FAIR:
  - collision mode: any time overlap of two frames on the same channel destroys all of them (the `BLE_SIM.py` rule, checked against Ancans et al. [102]); surviving frames are lost with a residual PER taken from the BER model at the link's nominal SNR;
  - SINR mode: per-chunk SINR from all overlapping transmissions weighted by adjacent-channel rejection (section 5.6), then the protocol's BER and packet-error model, plus the co-channel C/I floor of the receiver.

## 2. Fair-variant rule and feature matrix

Rule (decided 2026-09-23). The FULL variant runs the protocol with the robustness features its standard or reference stack normally uses. The FAIR variant removes every robustness feature the BCC stack does not implement:
- removed in FAIR: F1 FEC, F2 repetition or spreading beyond the base modulation, F3 ARQ, F4 carrier sensing, F5 adaptive channel selection, F11 link adaptation;
- kept in FAIR: F6 blind hopping (with an option to strip it), F7 payload CRC, F8 header protection, F9 base modulation including orthogonal symbol mapping or DSSS, F10 scheduled access, F13 whitening or scrambling;
- off in both variants unless intrinsic: F12 encryption, F14 interleaving.

F7 asymmetry: the BCC stack has no payload CRC (checked in `HBCProtocol.py`, section 3.8), while every other protocol keeps its CRC in FAIR. In FAIR the CRC changes only the airtime, because a corrupted frame counts as lost in both cases. A switch adds a 2-octet FCS (64 chips) to the BCC frame so all stacks carry the same detection.

Matrix cell format: `presence; FULL / FAIR`. Presence is yes, no or opt (optional in the standard). "-" means nothing to do.

| id | BCC stack | HBC standard | BLE 5 | 802.15.4 CSMA | 802.15.4 TSCH | 802.15.6 NB | SmartBAN | ANT |
|---|---|---|---|---|---|---|---|---|
| F1 FEC | no; - / - | no (A); - / - | opt (LE Coded conv. K=4 r 1/2); Coded only as FULL+ case / off | no; - / - | no; - / - | yes BCH(63,51) PSDU, BCH(31,19) hdr; on / off | yes BCH(127,113) PSDU, BCH(36,22) hdr; on / off | no; - / - |
| F2 repetition, spreading | no (FSC path in code, disabled); - / - | yes FSC SF 8..64; SF 8 (sweep 16..64) / SF 1 | opt (S=8 mapper, ISO NSE>BN); with Coded only / off | no (DSSS is F9); - / - | no; - / - | yes S=4/2/1; S=4 (sweep) / S=1 | yes Nrep 2/4; Nrep 1 (2/4 only via F11 flag) / Nrep 1 | no; - / - |
| F3 ARQ | no; - / - | yes I-Ack; R=3 retx / off | yes SN/NESN, no retry limit; on / off | yes, 3 retries; on / AR=0 | yes, 3 retx (RFC 8180); on / off | yes I-Ack; R=3 retx / off | yes ACK/NACK; RT=1 (sweep 3) / off | opt (ack msg without auto retx; burst 5 retries); host retry, 6 attempts / broadcast only |
| F4 carrier sensing | no; - / - | no (slotted Aloha, not CS); - / - | no; - / - | yes CSMA/CA; on / off (ALOHA) | opt CCA; off (flag on) / off | yes CSMA/CA in contention phases; unused (scheduled data) / off | opt MCA; flag on, owner senses Tsch only / off | opt listen at channel open; on / off (fast channel initiation) |
| F5 adaptive channel | no (static FDMA tones); - / - | no (single 21 MHz channel); - / - | yes AFH; PDR blacklisting / off, all 37 | opt ED scan (upper layer); least-used channel (A) / random fixed | opt, not normative; off (FULL+ flag) / off | yes hub channel choice; ED min-energy channel (A) / random fixed | yes DCH choice and migration; scan + migration (A) / random fixed DCH | yes adaptive isochronous (patent), freq. agility; both flags off / off |
| F6 blind hopping | no; - / - | no; - / - | yes CSA#1/#2; CSA#2 over AFH map / CSA#1 over 37 (strip option) | no; - / - | yes; on / on (strip option) | opt hub hopping; off / off (hop option) | no; - / - | opt (adv. burst 6 freq.); only with adv. burst / same (strip option) |
| F7 payload CRC | no; - / - (FCS switch) | yes CRC-16 FCS; on / on | yes CRC-24; on / on | yes FCS-16; on / on | yes FCS-16; on / on | yes FCS-16; on / on | yes CRC-16 MAC + CRC-16 PHY; on / on | yes CRC-16 (S); on / on |
| F8 header check | yes CRC-8 (dmin 2); on / on | yes CRC-8 HCS; on / on | no (CRC covers header); - / - | no; - / - | no; - / - | yes CRC-4 HCS; on / on | yes CRC-4 + FCS-8; on / on | no (CRC covers address); - / - |
| F9 orthogonal / DSSS | yes Walsh-16; on / on | yes Walsh-16; on / on | no (binary GFSK); - / - | yes 16-ary DSSS; on / on | yes 16-ary DSSS; on / on | no (DBPSK/DQPSK); - / - | no (GFSK); - / - | no (GFSK); - / - |
| F10 scheduled access | opt (TDMA class, unused in results); - / hub TDMA | yes superframe; on / on | yes connection events; on / on | opt (beacon GTS, mode B); mode B on / on | yes slotframe; on / on | yes superframe; on / on | yes IBI TDMA; on / on | yes per-link period; on / on |
| F11 link adaptation | no (rate field inert); - / - | opt (4 rates, no algorithm); rate rule (A) / fixed | opt PHY update, power control; off (flag) / off | opt TX power; fixed 0 dBm / fixed | opt TX power; fixed / fixed | opt rate field; fixed per run (swept) / fixed | partial (static PHY scheme); static, adaptive Nrep flag (A) / fixed | no; - / - |
| F12 encryption | no; - / - | opt; off / off | opt MIC 4 B; off / off | opt; off / off | opt; off / off | opt; off / off | opt flag; off / off | opt; off / off |
| F13 whitening, scrambling | yes 32-stage scrambler; on / on | yes 32nd order; on / on | yes whitening; on / on | no; - / - | no; - / - | yes scrambler; on / on | opt scrambler; on / on | no (S); - / - |
| F14 interleaving | no; - / - | no; - / - | no; - / - | no; - / - | no; - / - | yes, only for S=2/4; on with S>1 / identity at S=1 | no; - / - | no; - / - |

What each variant does, one line per protocol:
- BCC stack. FULL: identical to FAIR (the stack has none of F1 to F5, F11); optional ablation "FULL-BCC" turns the FSC spreading on (F2). FAIR: the stack as coded: Walsh-16 chips, header CRC-8, scrambler, no payload CRC, no ACK, hub TDMA inside the body, one fixed tone.
- HBC standard. FULL: FSC spreading SF 8 (sweep 16/32/64), header at 164 kbps, 7-octet MAC header + FCS, I-Ack with 3 retransmissions, scheduled allocations, optional slotted-Aloha contention, rate rule. FAIR: same frame with SF 1 everywhere and the header at the PSDU rate, MAC header + FCS kept, no ACK, scheduled allocations, fixed rate.
- BLE 5. FULL: one ACL connection per sensor, CSA#2, SN/NESN ARQ with no retry limit, PDR-based AFH, 2M (and 1M) PHY, Coded as an extra FULL+ case. FAIR: same connection interval and anchors, one transmission per event without ACK or retry, CSA#1 blind hopping over all 37 channels (strip option: one fixed channel), fixed 1M or 2M, 0 dBm.
- 802.15.4 CSMA. FULL: unslotted CSMA/CA (BE 3..5, 4 backoffs), Imm-Ack with 3 retries, least-used channel from an ED scan at BAN start (mode A); in beacon mode B the same ARQ inside GTSs. FAIR: unslotted ALOHA (mode A) or GTS (mode B), no ACK, one random fixed channel.
- 802.15.4 TSCH. FULL: 10 ms slots, default 16-channel hopping, dedicated cells, Enh-Ack with 3 retransmissions, CCA off (flag), blacklisting off (FULL+ flag). FAIR: same schedule and hopping, no ACK, no CCA (strip-F6 option).
- 802.15.6 NB. FULL: BCH on header and PSDU, 121.4 kb/s (S=4, sweep other rates), I-Ack with 3 retransmissions in scheduled allocations, ED-scan channel choice at start. FAIR: uncoded π/2-DBPSK S=1 (600 kb/s raw), no ACK, same superframe, random fixed channel.
- SmartBAN. FULL: BCH on header and PSDU, Nrep 1, ACK policy 0 with RT=1 retransmission, MCA flag on, DCH chosen by CCH scan and migrated on failure. FAIR: uncoded, Nrep 1, no ACK, no MCA, same IBI layout, random fixed DCH.
- ANT. FULL: acknowledged data with host retry in the next slot (6 attempts), listen-before-open phase choice, burst with 5 retries above the broadcast rate; adaptive isochronous re-timing and frequency agility only as flags. FAIR: broadcast only, random phase drifting with the clock, fixed frequency, no retries (fast channel initiation).

## 3. Protocols

### 3.1 BLE 5 (Core 5.4, with 6.0 and 6.2 notes): LE 1M, LE 2M, LE Coded S=2 and S=8

| group | parameter | value | unit | source | mark |
|---|---|---|---|---|---|
| PHY | modulation | GFSK, BT 0.5, h 0.45..0.55; 1 Msym/s (1M, Coded), 2 Msym/s (2M) | - | [1] A 3.1 | V |
| PHY | TX power | power-class maximum between -20 and +20 (not a control range); nRF52840: +8 max, settings down to -40 | dBm | [1] A 3; [6] | V (corrected) |
| PHY | channels | 40 x 2 MHz, 2402..2480 MHz. Data index 0..10 at 2404 + 2 idx, index 11..36 at 2428 + 2 (idx - 11); advertising 37/38/39 at 2402/2426/2480 | MHz | [1] B 1.4.1 | V |
| frame | preamble | 1 octet (1M, 8 µs), 2 octets (2M, 8 µs), 80 µs uncoded on Coded | - | [1] B 2.1, 2.2 | V |
| frame | access address | 4 octets (32 µs 1M, 16 µs 2M, 256 µs Coded) | - | [1] B 2.1 | V |
| frame | data PDU header | 16 bits (24 with CTEInfo): LLID 2, NESN, SN, MD, CP, Length 8 | bit | [1] B 2.4 | V |
| frame | payload L | 0..27 octets without DLE, 0..251 with DLE; Length counts payload + MIC | octet | [1] B 2.4, 4.5.10 | V |
| frame | MIC | 4 octets, encrypted non-empty PDUs only; off in both variants | octet | [1] B 2.4 | V |
| frame | CRC | 24 bits over header + payload + MIC; no separate header check | bit | [1] B 2.1.4 | V |
| airtime | LE 1M | 80 + 8L | µs | [1] B 2.1 (derived, range 44..2128 reproduced) | V |
| airtime | LE 2M | 44 + 4L | µs | [1] B 2.1 (derived) | V |
| airtime | LE Coded S=2 | 462 + 16L | µs | [1] B 2.2, Table 2.1 (derived) | V |
| airtime | LE Coded S=8 | 720 + 64L; the 2704 µs of B 4.5.10 is 27 octets + 4-octet MIC | µs | [1] B 2.2, 4.5.10 | V |
| airtime | empty PDU | 80 (1M), 44 (2M), 462 (S=2), 720 (S=8) | µs | derived | V |
| coding | convolutional code (Coded) | rate 1/2, K=4, G0 = 1+x+x^2+x^3, G1 = 1+x^2+x^3, zero start, 3-bit termination | - | [1] B 3.3.1 | V |
| coding | pattern mapper | P=1 (S=2); P=4 (S=8: 0 -> 0011, 1 -> 1100) | symbol/bit | [1] B 3.3.2 | V |
| coding | coding indicator | CI 00 = S=8, 01 = S=2 for FEC block 2; FEC block 1 (AA, CI, TERM1) always S=8 | - | [1] B 2.2.3 | V |
| coding | distance spectrum | dfree 6; A_d (d=6..14) 1,3,5,11,25,55,121,267,589; B_d 2,7,18,49,130,333,836,2069,5060 | - | trellis search, reproduced by the verifier | A |
| coding | whitening | LFSR x^7+x^4+1 seeded with the channel index, over PDU + CRC | - | [1] B 3.2 | V |
| timing | T_IFS | 150 (Core 5.4). From 6.0 T_IFS_ACL_CP/PC and T_IFS_CIS are negotiable (default 150, FS_Max <= 10 ms); T_MCES default 150 | µs | [1] B 4.1.1; [3] B 4.1.3 | V |
| timing | connInterval | 7.5 ms..4 s in 1.25 ms steps (5.4). Core 6.2: 375 µs..4 s in 125 µs steps, creation and update use values >= 7.5 ms | - | [1] B 4.5.1; [4] B 4.5.1 | V |
| timing | supervision timeout | 100 ms..32 s in 10 ms steps, > (1 + latency) x subrate x CI x 2; not established after 6 events means lost | - | [1] B 4.5.2 | V |
| timing | clock accuracy | active +/-50 ppm, 2 µs jitter; sleep +/-500 ppm, 16 µs jitter | ppm | [1] B 4.2.1, 4.2.2 | V |
| timing | receive window | windowWidening = (SCA_central + SCA_peripheral) x time since last anchor; listen from windowWidening before to windowWidening after the expected start | - | [3] B 4.2.4 | V (corrected) |
| timing | nRF52840 ramp-up | tTXEN = tRXEN = 140 (fast 40); tTXDIS 6 (1M), 4 (2M); tRXDIS 0 | µs | [6] 6.20.15.8 | V |
| event | connection event | Central sends at the anchor; the two sides alternate with T_IFS on one channel; the Peripheral answers even after a CRC failure; the event ends >= T_IFS (6.0: T_MCES) before the next anchor | - | [1] B 4.5.1 | V |
| event | closing | neither side has MD = 1; or a missed packet; or 2 consecutive CRC failures | - | [1] B 4.5.6 | V |
| ARQ | SN/NESN | 1-bit ack in the peer's next packet; a PDU with failed CRC is not acked and is resent; duplicates dropped by SN | - | [1] B 4.5.9 | V |
| ARQ | retry limit (ACL) | none; LE-U data is non-flushable; only the supervision timeout ends it | - | [1] B 4.5.9; [2] HCI 5.4.2 | V |
| hopping | CSA#1 | unmapped = (last + hopIncrement) mod 37, last = 0 at the first event; unused: remap = unmapped mod N_used into the ascending used table | - | [1] B 4.5.8.2 | V |
| hopping | hopIncrement | uniform 5..16 (same in Core 4.2) | - | [1] B 2.3.3.1; [5] | V |
| hopping | CSA#2 | chanId = AA[31:16] XOR AA[15:0]; prn = counter XOR chanId; 3 x {bit-reverse each byte; prn = (17 prn + chanId) mod 2^16}; prn_e = prn XOR chanId; ch = prn_e mod 37; unused: ch = used[floor(N_used prn_e / 2^16)] | - | [1] B 4.5.8.3 | V |
| hopping | self tests | CSA#2, chanId 0x305F, all 37 used: counters 0..3 -> 25, 20, 6, 21. CSA#2, used {9,10,21,22,23,33,34,35,36}: counters 6, 7, 8 -> 23, 9, 34. CSA#1, hop 10, channel 11 removed at instant 100: events 99, 100, 101 -> 1, 12, 21 | - | [8]; [1] B 5.1.2 | S (Zephyr) / V (spec example) |
| hopping | CSA choice | CSA#1 for legacy connections where a side has ChSel = 0, else CSA#2; CIS and BIS always CSA#2 | - | [1] B 4.5, 4.5.13.6 | V |
| AFH | minimum used channels | 2 (LE); LL_MIN_USED_CHANNELS_IND 2..37; the value 20 is BR/EDR AFH | channel | [1] B 4.5.8.1; [5]; [2] Baseband 2.3.1 | V |
| AFH | classification | vendor specific; only the Central changes the map (LL_CHANNEL_MAP_IND, 8-octet control PDU) | - | [1] B 4.5.8.1 | V |
| AFH | instant | absolute connEventCount written at the first transmission; the Central "should" allow >= 6 events; if the Peripheral has not received it by the instant, the procedure fails | event | [1] B 5.5.1 | V |
| AFH | practical rule | window 20 outcomes per channel, blacklist when PDR < 95 %, new map at N + 6, re-enable all 37 when used < C_min; during re-enabling the CI is divided by OWL = ceil(37/C_min) | - | [7] Sec. 4.1-4.3 | V |
| AFH | C_min | 10 (worked example in [7], not a design constant) | channel | [7] Sec. 4.3 | A |
| ISO | CIS | ISO_Interval 5 ms..4 s (1.25 ms steps, < supervision timeout / 2); BN 0..15; NSE 1..31; FT 1..255; Max_PDU 0..251; needs an associated ACL | - | [1] B 4.5.13 | V |
| ISO | BIS | ISO_Interval 5 ms..4 s; BN 1..7; NSE 1..31 multiple of BN; IRC 1..15; PTO 0..15; Max_PDU 1..251; no ACK; needs a periodic advertising train carrying BIGInfo | - | [1] B 4.4.6 | V |
| receiver | spec sensitivity | uncoded <= -70, S=2 <= -75, S=8 <= -82; BER 0.1 % (payload 1-37 B), 0.064 % (38-63), 0.034 % (64-127), 0.017 % (128-255) | dBm | [1] A 4.1, Tables 4.1, 4.2 | V |
| receiver | spec C/I 1M | co-channel 21; 1 MHz 15; 2 MHz -17; >= 3 MHz -27; image -9; image +/-1 MHz -15 (5 spurious channels may relax to -17) | dB | [1] A Table 4.3 | V |
| receiver | spec C/I 2M | co 21; 2 MHz 15; 4 MHz -17; >= 6 MHz -27; image -9; image +/-2 MHz -15 | dB | [1] A Table 4.4 | V |
| receiver | spec C/I Coded S=2 | co 17; 1 MHz 11; 2 MHz -21; >= 3 MHz -31; image -13; image +/-1 MHz -19 | dB | [1] A Table 4.6 | V |
| receiver | spec C/I Coded S=8 | co 12; 1 MHz 6; 2 MHz -26; >= 3 MHz -36; image -18; image +/-1 MHz -24 | dB | [1] A Table 4.5 | V |
| receiver | max usable input | > -10 at BER <= 0.1 % | dBm | [1] A 4.5 | V |
| receiver | nRF52840 sensitivity | 1M -95 (<= 37 B, BER 1e-3), -94 (>= 128 B, BER 1e-4); 2M -92; S=2 (500k) -99; S=8 (125k) -103 | dBm | [6] 6.20.15.5 | V |
| receiver | nRF52840 C/I 1M | co 6; -1 MHz -2; +1 MHz -9; -2 MHz -22; +2 MHz -46; >= 3 MHz -50; image -22 | dB | [6] 6.20.15.6 | V |
| receiver | nRF52840 C/I 2M | co 6.8; +/-2 MHz -10; +/-4 MHz -45; >= 6 MHz -48; image -24 | dB | [6] 6.20.15.6 | V |
| receiver | nRF52840 C/I 125k | co 4.4; -1 MHz -4.0; +1 MHz -12; -2 MHz -28; +2 MHz -50; >= 3 MHz -55; image -29; no LE500k table exists in v1.5 | dB | [6] 6.20.15.6 | V |
| legacy | BLE_SIM.py | airtime t_on = payload_bits / phy_rate (no 80/44 µs overhead, no poll or T_IFS); one shuffled 37-channel permutation per body, shared by its N_dev links; residual loss 1e-4 (1M), 1e-3 (2M) with no BER basis | - | [100] BLE_SIM.py lines 14-15, 85, 94-96 | V |
| host | L2CAP + ATT notification overhead | 4 + 3 = 7 octets per notification (optional, count inside L) | octet | not read | A |

Simulation recipe, common to both variants:
- One hub (Central) and N_dev sensors (Peripherals) per body. Data channel frequency from the index rule above.
- Anchor of link j on hub h: anchor_j(n) = φ_h + j CI / N_dev + n CI, with φ_h uniform in [0, CI) per body (the BLE_SIM staggering). Optional clock error per body uniform in +/-50 ppm.
- CI a multiple of 1250 µs in 7.5 ms..4 s. Core 6.2 short intervals only as a sensitivity case.
- Payload: L <= 251 (DLE assumed); larger application frames are split into several PDUs.
- Hub scheduling when two of its links' events collide: the earlier anchor wins and the later event is skipped, its counter still advances (assumed, open issue OI-5).

FULL (F3 + F5 + CSA#2; F1, F2, F11 as flags):
1. One ACL connection per sensor, random 32-bit AA, CSA#2 over the used-channel map, connEventCounter from 0.
2. Each event on channel ch: Central sends C (queued downlink, unacked retransmission, else an Empty PDU). If the Peripheral misses C, the event ends. Otherwise, T_IFS after C, the Peripheral sends P (its unacked PDU, else the next queued PDU, else Empty). SN/NESN as in [1] B 4.5.9: a side updates its state only from packets with a valid CRC.
3. The Central sends another exchange only if the last P was received, fewer than 2 consecutive CRC errors occurred at either side, MD was set by a side (or the Central holds an unacked PDU), and the exchange plus 150 µs fits before the next anchor of this link and before the hub's next scheduled event.
4. Unacked PDUs wait for the next event on the next channel. No retry limit. A payload is delivered when the hub receives it with a valid CRC and a new SN; it is late if not delivered within the application deadline D (default one sample period, open issue OI-4). Report PDR and on-time ratio separately.
5. Supervision timeout TO = max(2 s, 3 CI) rounded up to 10 ms (assumed); when no valid packet arrives for TO, the link drops, its queue counts as lost, and it reconnects after 100 ms (assumed).
6. AFH at the Central per link: success = valid-CRC P received. Sliding window of 20 outcomes per channel; PDR < 0.95 marks the channel bad and queues an LL_CHANNEL_MAP_IND (L = 8, acked like data); the new map applies at instant = n_first_tx + 6 (corrected: the instant is fixed at the first transmission, not after the ack). If a removal would leave fewer than C_min = 10 channels, re-enable all 37 and clear the statistics. The OWL interval speed-up of [7] is left out (declared deviation; flag).
7. PHY: default 2M, second run 1M. FULL+ case: Coded S=2 or S=8 for all packets. F11 flag (off): switch 2M to 1M to Coded when link PER > 10 % over 100 events (assumed). TX power 0 dBm.

FAIR (no F1 to F5, F11; keeps F6, F7, F10, F13):
1. Same CI, anchors and staggering.
2. One transmission per event, never repeated. Default option B (BIS-like, uplink only, BN = NSE = IRC = 1): the sensor sends its PDU once at its anchor, no poll, no ACK; this matches `BLE_SIM.py`. Option A (CIS-like, BN = NSE = FT = 1): hub poll C (Empty), T_IFS, sensor P; use it when the poll should count as interference. Both options omit the standard's side traffic (CIS needs an ACL kept alive; BIS needs a periodic advertising train); this omission is declared.
3. Channel, blind: CSA#1 per link with hopIncrement uniform 5..16 and a random start u0 in 0..36, ch(n) = (u0 + (n+1) hopIncrement) mod 37 (reproduces the spec example with u0 = 0). CSA#1 over ISO transport is not standard (ISO uses CSA#2); option CSA#2 with the full map; strip option: one random fixed channel per BAN.
4. PHY 1M or 2M uncoded (run both), 0 dBm, no PHY update.
5. Delivered iff P is received with a valid CRC. Latency <= CI.
6. Fix for `BLE_SIM.py` alignment: add the 80/44 µs overhead to t_on, and use CSA#1 per link instead of one shuffled permutation per body.

BER model:
- Noise bandwidth = symbol rate: 1 MHz (1M, Coded), 2 MHz (2M). γ_s = P_rx / (N + Σ_j I_j 10^(-ACR(Δf_j)/10)), N = -174 + 10 log10(B) + NF + L_impl.
- ACR(Δf) = C/I_co - C/I(Δf). Spec: 1M 6 / 38 / 48 dB at 1 / 2 / >= 3 MHz; 2M 6 / 38 / 48 dB at 2 / 4 / >= 6 MHz (a 2M neighbour channel 2 MHz away is only 6 dB down); Coded 6 / 38 / 48 dB. nRF52840 1M: 28 to 52 dB at 2 MHz (asymmetric).
- Uncoded bit error, default (a): Pb = 0.5 exp(-γ/2). Corrected by the verifier: this is exact only for orthogonal tones (h = 1, [11]); at h = 0.5 it is an optimistic curve used as a calibrated fit. BER 1e-3 needs 10.94 dB, 1e-4 12.31 dB, 1e-5 13.35 dB. With (a), NF + L_impl = 8.1 dB reproduces the nRF52840 at 1M (-95 dBm) and 2M (-92 dBm).
- Conservative (b): correlated non-coherent FSK, Pb = Q1(a,b) - 0.5 exp(-(a^2+b^2)/2) I0(ab), a,b = sqrt(γ/2 (1 -/+ sqrt(1-ρ^2))), |ρ| = sin(πh)/(πh) (form of [9] eq. 14 with ρ changed); 1e-3 at 15.10 / 14.31 / 13.63 dB for h = 0.45 / 0.50 / 0.55; then NF + L_impl = 4.7 dB fits the nRF52840. Q1(a,b) = scipy.stats.ncx2.sf(b^2, 2, a^2). NF + L_impl is a fit parameter in both cases.
- Threshold (c): lost if SIR < C/I_co anywhere in the packet; C/I_co = 21 dB (1M, 2M), 17 (S=2), 12 (S=8) by the spec, or 6 / 6.8 / 4.4 dB on the nRF52840.
- Packet: bits at risk = AA 32 + header 16 + 8L + CRC 24 = 72 + 8L. P(sync) = Π over AA chunks of (1 - Pb)^n; P(CRC ok | sync) over the other 40 + 8L bits. A lost AA means the peer does not answer.
- Coded PHY: per coded bit γ_c = γ_s (S=2) or 4 γ_s (S=8, ideal combining of the 4 chips). Hard decision: P_d = Σ_{k>d/2} C(d,k) p^k (1-p)^(d-k) + 0.5 C(d,d/2) (p(1-p))^(d/2) for even d. Soft non-coherent: P_d = 2^-(2d-1) e^(-d γ_c / 2) Σ_{k=0..d-1} c_k (d γ_c / 2)^k, c_k = (1/k!) Σ_{n=0..d-1-k} C(2d-1, n). First-event error P_E = min(1, Σ A_d P_d). PER = 1 - (1 - P_E(S=8))^34 (1 - P_E(S))^(40 + 8L).
- Coded gains over 1M at BER 1e-3 (model): S=2 4.67 dB soft, 3.17 dB hard; S=8 10.69 dB soft, 9.19 dB hard. Spec: 5 and 12 dB; nRF52840: 4 and 8 dB. Default: soft model, extra loss 0 dB for S=2 and 2.7 dB for S=8 (hits -103 dBm with NF + L_impl = 8.1 dB). Against interference use the Coded C/I tables: the co-channel gain is only 4 dB (S=2) and 9 dB (S=8), not 5 and 12 dB (corrected).
- ARQ is simulated attempt by attempt (retries in one event share the interferer and are correlated).

### 3.2 IEEE 802.15.4 at 2.4 GHz: O-QPSK PHY, unslotted CSMA/CA (mode A) and beacon-enabled GTS (mode B)

All PHY and MAC values were read in IEEE Std 802.15.4-2006 [13]; they are assumed unchanged in the 2015 and 2020 editions [14] (constants in 2020 confirmed secondarily through Zephyr [25]).

| group | parameter | value | unit | source | mark |
|---|---|---|---|---|---|
| PHY | band and rates | 2400-2483.5 MHz; 2000 kchip/s; 250 kb/s; 62.5 ksym/s; 16-ary orthogonal | - | [13] Table 1, 6.5.1 | V |
| PHY | channels | k = 11..26, fc = 2405 + 5 (k - 11) | MHz | [13] 6.1.2.1 | V |
| PHY | symbol, octet | 16 µs (+/-40 ppm); octet = 2 symbols = 32 µs | µs | [13] 6.5.3.2 | V |
| PHY | spreading | 4 bits -> one of 16 near-orthogonal 32-chip sequences; 8 chips per bit; minimum pairwise Hamming distance 12 | - | [13] 6.5.2.3 Table 24; distance computed | V |
| PHY | chip modulation | O-QPSK, half-sine pulses; no scrambler, FEC or interleaver | - | [13] 6.5.2.4, Figure 18 | V |
| frame | SHR | preamble 4 octets (8 symbols, 128 µs) + SFD 1 octet | octet | [13] 6.3.1, 6.3.2, Table 19 | V |
| frame | PHR | 1 octet (7-bit length + 1 reserved), no header check | octet | [13] 6.3.3 | V |
| frame | aMaxPHYPacketSize | 127 (PSDU); phyMaxFrameDuration 266 symbols = 4256 µs | octet | [13] Tables 22, 23 | V |
| frame | airtime | T = (6 + L_psdu) x 32 | µs | derived from [13] 6.3 | V |
| frame | data frame | L_psdu = 11 + L_app (MHR 9 with short addresses and PAN ID compression, FCS 2); L_app <= 116 | octet | [13] 7.2.1 Figures 41, 42 (sum derived) | V |
| frame | Imm-Ack | 5-octet PSDU (FC 2, Seq 1, FCS 2), 352 µs | octet | [13] 7.2.2.3 Figure 53 | V |
| frame | FCS | 16 bits, x^16 + x^12 + x^5 + 1 | bit | [13] 7.2.1.9 | V |
| frame | beacon | 13 octets without GTS, 35 octets with 7 GTS (608 / 1312 µs) | octet | [13] 7.2.2.1 (sum derived) | V |
| timing | aTurnaroundTime | 12 symbols = 192 µs (RX-to-TX and TX-to-RX maximum) | symbol | [13] Table 22 | V |
| timing | aUnitBackoffPeriod | 20 symbols = 320 µs | symbol | [13] Table 85 | V |
| timing | CCA duration | 8 symbols = 128 µs | symbol | [13] 6.9.9 | V |
| timing | macAckWaitDuration | 54 symbols = 864 µs (20 + 12 + 10 + 6 x 2) | symbol | [13] 7.4.2 Eq. 13 | V |
| timing | ACK start | aTurnaroundTime after the last data symbol (non-beacon, CFP); in the CAP optionally at a backoff boundary; no CSMA | - | [13] 7.5.6.4.2 | V |
| timing | SIFS / LIFS | 12 / 40 symbols (192 / 640 µs); SIFS after MPDU <= 18 octets (aMaxSIFSFrameSize) | symbol | [13] Table 3, Table 85 | V |
| CSMA | macMinBE | 3 (range 0..macMaxBE) | - | [13] Table 86 | V |
| CSMA | macMaxBE | 5 (range 3..8) | - | [13] Table 86 | V |
| CSMA | macMaxCSMABackoffs | 4 (range 0..5) | - | [13] Table 86 | V |
| CSMA | macMaxFrameRetries | 3 (range 0..7) | - | [13] Table 86 | V |
| CSMA | unslotted algorithm | NB = 0, BE = macMinBE; wait U{0..2^BE - 1} x 320 µs; CCA; busy: NB++, BE = min(BE+1, macMaxBE), fail if NB > macMaxCSMABackoffs; idle: transmit. NB and BE are reset for every attempt, retransmissions included | - | [13] 7.5.1.4 Figure 69 | V |
| CSMA | retransmission | no ACK within macAckWaitDuration: repeat up to macMaxFrameRetries times with the same DSN; in beacon mode only inside the same portion (CAP or GTS), else deferred to the same portion of the next superframe | - | [13] 7.5.6.4.3 | V |
| CCA | modes and threshold | mode 1 energy, 2 carrier sense, 3 both; ED threshold at most 10 dB above the specified sensitivity (<= -75 dBm for the -85 dBm requirement) | dBm | [13] 6.9.9; [15] 6.9.9 | V |
| CCA | ED measurement | 8 symbols, 8-bit value, range >= 40 dB, +/-6 dB | - | [13] 6.9.7 | V |
| RX/TX | sensitivity requirement | -85 or better (PER < 1 %, 20-octet PSDU) | dBm | [13] 6.5.3.3, Table 4 | V |
| RX/TX | channel rejection requirement | adjacent (5 MHz) 0 dB, alternate (10 MHz) 30 dB, wanted 3 dB above -85 dBm | dB | [13] 6.5.3.4 Table 26 | V |
| RX/TX | TX power, max input | TX capable of >= -3 dBm; RX max input >= -20 dBm | dBm | [13] 6.9.5, 6.9.6 | V |
| beacon | superframe | aBaseSlotDuration 60 symbols, 16 slots, aBaseSuperframeDuration 960 symbols (15.36 ms); BI = 15.36 ms x 2^BO, SD = 15.36 ms x 2^SO, 0 <= SO <= BO <= 14; BO = 15 means non-beacon | - | [13] Table 85, 7.5.1.1 | V |
| beacon | CAP, GTS | aMinCAPLength 440 symbols (7.04 ms); up to 7 GTS in the CFP, one TX and/or one RX GTS per device, no CSMA in the CFP | - | [13] Table 85, 7.5.7 | V |
| beacon | beacon loss | aMaxLostBeacons 4; aGTSDescPersistenceTime 4 | - | [13] Table 85 | V |
| chip | nRF52840 | sensitivity -100 dBm; no 802.15.4 C/I given; tTXEN/tRXEN 130 µs (fast 40), tTXDIS 21 µs, turnaround <= 40 µs | - | [6] 6.20.15 (local v1.5 text) | V |
| chip | CC2652R | -99 dBm typ (PER 1 %), -100 dBm (front page); ACR 36 dB (+/-5 MHz), 57 dB (+/-10), 59 dB (>= 15), wanted -82 dBm; no co-channel figure | - | [22] 8.12 | V |
| chip | CC2420 | -95 dBm typ (-90 min); ACR +5 MHz 45, -5 MHz 30, +10 MHz 54, -10 MHz 53, >= 15 MHz 62 dB; co-channel rejection -3 dB (about +3 dB C/I needed) | - | [23] 6.3 | V |
| chip | AT86RF233 | -101 dBm (250 kb/s); NF 6 dB; ACR -5 MHz 32, +5 MHz 35, +/-10 MHz 48, +/-15 MHz 54 dB (P_RF -82 dBm) | - | [24] 12.7 | V |

Simulation recipe, common:
- Airtime T(L_psdu) = (6 + L_psdu) x 32 µs. Channels 11..26; in collision mode only same-channel overlaps interact. Half-duplex radios; 192 µs turnaround in protocol timing (chip energy uses 40 µs, section 4).
- Each sensor generates one L_app frame per reporting period P. Each BAN has an independent uniform time origin.

Mode A, unslotted CSMA/CA. FULL:
1. Channel: at BAN start the hub ED-scans 16 channels and picks the one with the fewest active BANs (ties random); BANs start one after another in random order. This is a Zigbee-style upper-layer policy, not 802.15.4 behaviour (A).
2. Per frame: retry = 0; NB = 0, BE = 3; wait U{0..2^BE - 1} x 320 µs; CCA 128 µs (busy if any co-channel frame is on air during the window, assuming every such frame is above the ED threshold); busy: NB++, BE = min(BE+1, 5), drop when NB > 4; idle: 192 µs turnaround, then the data frame.
3. The hub answers with an Imm-Ack 192 µs after the data end. The sender succeeds if the ACK arrives intact within 864 µs of the data end; otherwise retry++ and, if retry <= 3, restart step 2 with NB = 0, BE = 3; else drop.
4. IFS after the ACK (or timeout): 640 µs if MPDU > 18 octets, else 192 µs. Mean successful exchange = 1120 + 128 + 192 + T_data + 192 + 352 + IFS µs; for L_app = 20: 3808 µs.
5. Queue while a frame is in backoff or retry: FIFO depth 4, drop newest (assumed, OI-18). Duplicates from lost ACKs are discarded by DSN.

Mode A FAIR: channel fixed per BAN, uniform from 11..26 (worst-case option: all BANs on one channel); no CCA, no backoff (unslotted ALOHA at generation time); AR = 0; IFS and FCS kept; airtime T(11 + L_app).

Mode B, beacon-enabled with GTS:
- The hub sends a beacon at slot 0 without CSMA (PSDU 13 + 1 + 3 n_gts octets). One TX GTS per sensor, contiguous at the end of the superframe, so N_dev <= 7 per BAN.
- The period is quantised: BI = 15.36 ms x 2^BO (15.36, 30.72, 61.44, 122.88 ms, ...). The simulator rounds P to this grid, and runs the other protocols on the same P for the comparison (corrected).
- Pick the smallest SO whose 16 - n_gts g slots still hold the 440-symbol CAP; g slots per GTS hold one FAIR attempt (T_data + IFS) or the FULL attempts (T_data + 192 + 352 + IFS each).
- A sensor that misses 4 consecutive beacons loses its frames until it hears a beacon again.
- FULL: ED-scan channel (as mode A), AR = 1, retry inside the GTS if it fits, else in the same GTS next superframe, up to 3 retries. FAIR: random channel, AR = 0, CAP used only by the beacon.
- Energy note: mode A keeps the hub in continuous RX; report hub energy separately or treat the hub as phone-powered.

BER model (both modes, both variants, and TSCH):
- Pb(SINR) = (8/15)(1/16) Σ_{k=2..16} (-1)^k C(16,k) exp(20 SINR (1/k - 1)), SINR linear in a 2 MHz noise bandwidth ([13] Annex E.4.1.8 item 7, from IEEE 802.15.2-2003 5.3; same code in ns-3 [21]).
- This is the textbook non-coherent 16-ary orthogonal formula with Es/N0 = 20 SINR, so the implied Eb/N0 = 5 SINR (+6.99 dB). Feed SINR directly; do not use Eb/N0 = 8 SINR (corrected).
- Evaluated: BER 1e-3 / 1e-4 / 1e-5 / 1e-6 at SINR -0.92 / 0.21 / 1.09 / 1.81 dB; 1.6e-4 at 0 dB, 1.6e-2 at -3 dB, 0.12 at -6 dB. PER 1 % at 0.40 dB (20-octet PSDU) or 0.51 dB (26 octets with SHR and PHR) (A, numerical).
- The formula is an ideal receiver: [13] Annex E.3.2.1 says a typical low-cost detector meets 1 % PER at 5 to 6 dB SNR. Calibrate NF + L_impl to the chip: nRF52840 (-100 dBm) 10.6 dB, AT86RF233 (-101 dBm) 9.6 dB, CC2420 (-95 dBm) 15.6 dB, spec (-85 dBm) 25.6 dB.
- PER = 1 - Π_i (1 - Pb(SINR_i))^(bits_i) over constant-SINR chunks (packet form as in [101]), bits = 8 (6 + L_psdu). ns-3 also drops frames whose SINR at preamble start is <= -5 dB (S).
- Co-channel floor in SINR mode: C/I_co = 3 dB (CC2420 co-channel rejection -3 dB); the Gaussian formula alone would accept 0.4 dB.
- ACR: typical chip 33 / 48 / 54 dB at 5 / 10 / 15 MHz (AT86RF233); spec minimum 0 / 30 dB at 5 / 10 MHz (parameter; with near-far between bodies the 0 dB case matters).
- Optional chip-level model (A): pc = Q(sqrt(2 Ec/N0)), hard decision corrects t = 5 chip errors per 32-chip symbol; Ps = Σ_{i=6..32} C(32,i) pc^i (1-pc)^(32-i).

### 3.3 IEEE 802.15.4 TSCH (802.15.4e, now in 802.15.4-2015/2020)

PHY, frame, FCS, BER model and chip data as in section 3.2. The IEEE TSCH clauses were not readable; the timeslot template comes from draft-ietf-6tisch-minimal-05 and Contiki-NG.

| group | parameter | value | unit | source | mark |
|---|---|---|---|---|---|
| slot | macTsTimeslotLength | 10000 | µs | [17] Sec. 3 table; [18] tsch-timeslot-timing.c | S |
| slot | macTsCcaOffset / macTsCca | 1800 / 128 | µs | [17]; [18] | S |
| slot | macTsTxOffset | 2120 (slot start to frame start) | µs | [17]; [18] | S |
| slot | macTsRxOffset | 1120 in [17] (read); the Contiki-NG code uses TxOffset - RxWait/2 = 1020; affects idle-listen energy only | µs | [17]; [18] | S |
| slot | macTsRxAckDelay / macTsTxAckDelay | 800 / 1000 (from data end) | µs | [17]; [18] | S |
| slot | macTsRxWait / macTsAckWait | 2200 / 400 | µs | [17]; [18] | S |
| slot | macTsRxTx / macTsMaxAck / macTsMaxTx | 192 / 2400 / 4256 | µs | [17]; [18] | S |
| hopping | channel | ch = HS[(ASN + chOff) mod len(HS)]; ASN counts slots (5 octets) | - | [16] RFC 7554 A.6-A.7 | S |
| hopping | default sequence (16 ch) | [16, 17, 23, 18, 26, 15, 25, 22, 19, 11, 12, 13, 24, 14, 20, 21] (LFSR x^9 + x^5 + 1, seed 255); RFC 8180 requires macHoppingSequenceID 0 | - | [18] tsch-const.h; [16] RFC 8180 Fig. 1 | S |
| ARQ | retransmissions | 3 (4 attempts); frames without ACK request (EBs) are not retransmitted | - | [16] RFC 8180 Fig. 1, 4.3 | V |
| ARQ | Contiki-NG defaults | TSCH_MAC_MAX_FRAME_RETRIES 7; CCA disabled (TSCH_CCA_ENABLED 0); default sequence 4 channels {15, 25, 26, 20}; queue 4 per neighbour; EB period 16 s | - | [18] tsch-conf.h | S |
| backoff | macMinBE / macMaxBE (TSCH) | 1 / 7 (Chen et al.); Contiki-NG 1 / 5; keep as a parameter | - | [19] Table 2; [18] | S |
| backoff | shared-cell rule | backoff only in shared cells; the first attempt has no backoff; BE = macMinBE on the first failure, then BE = min(BE+1, macMaxBE), wait U{0..2^BE - 1} shared cells; reset on success in a shared cell, or on success in a dedicated cell that empties the queue | - | [20] (802.15.4-2015 6.2.5.3 quoted) | S (corrected) |
| beacon | EB period | not defined by 802.15.4; tunable (RFC 8180 4.5.2); 16 s in Contiki-NG | s | [16]; [18] | S |
| frame | EB | about 41-octet PSDU (30 octets of IEs + 9-octet MHR + FCS) | octet | [16] App. A.1 plus sum | A |
| frame | Enh-Ack | about 13-octet PSDU (FC 2, Seq 1, PAN 2, dst 2, time-correction IE 4, FCS 2), 608 µs | octet | built as in Contiki-NG | A |

Simulation recipe, common: 10 ms slots; each BAN has its own ASN origin and slot phase; slotframe length S >= N_dev + 1 with S x 10 ms = P (P quantised to 10 ms steps). Slot 0: shared advertising cell, chOff 0, EB every 16 s (unacked). Slots 1..N_dev: dedicated TX cell sensor i -> hub, chOff = i mod 16 (A). Frame at 2120 µs into the slot; the receiver listens from 1120 µs for up to 2200 µs.

FULL: AR = 1 with Enh-Ack 1000 µs after the data end on the same channel; up to 3 retransmissions in the sensor's next TX cell (dedicated cells: next slotframe, no backoff); FIFO depth 4 (Contiki-NG) with new frames behind the retried one; CCA off by default (Contiki-NG), flag on: CCA at 1800 µs for 128 µs, busy counts as a failed attempt (A, Contiki-NG behaviour). FULL+ flag (off): blacklist channels with PDR < 0.8 over the last 100 attempts, keep >= 4, HS' = HS minus blacklist (A, non-normative).

FAIR: AR = 0, no retry, no CCA, no blacklist; same template and hopping; strip-F6 option: channel = HS[chOff] constant.

Across unsynchronised BANs the channel match rate is about 1/16 per overlapping pair; blacklisting helps little against other TSCH BANs because their interference is spread over all channels.

### 3.4 IEEE 802.15.6 NB PHY (2400-2483.5 MHz) and the common 802.15.6 MAC

The IEEE 802.15.6-2012 text [26] was not readable. The values below come from the TI patent that reproduces the TG6 NB PHY proposal tables [27] and from reviews and simulation papers. A 2026 revision (802.15.6-2026, approved 2026-06-04, up to 50 Mb/s) supersedes 2012 [26].

| group | parameter | value | unit | source | mark |
|---|---|---|---|---|---|
| PHY | channel plan | fc = 2402 + nc, nc = 0..78 (79 channels, 1 MHz); the ns-3 study configured 16 channels | MHz | [27] FIG. 24 | A (single draft-era source) |
| PHY | occupied bandwidth | 1 MHz at -20 dBr; no adjacent-channel rejection figure found | MHz | [27] FIG. 28 | S |
| PHY | symbol rate | 600 ksym/s at every rate; SRRC pulses, roll-off not found (0.3 used by Nasr) | ksym/s | [27]; [28] Table 3; [30] | S |
| PHY | PSDU rates | 121.4 (π/2-DBPSK, S=4), 242.9 (S=2), 485.7 (S=1), 971.4 (π/4-DQPSK, S=1); all BCH(63,51); π/8-D8PSK only below 1 GHz | kb/s | [27] FIG. 2G; [30] Table I; [28] Table 3 | S |
| PHY | rate formula | Rd = Rs log2(M) / S x k/n | - | [28] Eq. 1 | S |
| PHY | PSDU FEC | systematic BCH(63,51), t = 2; N_CW = ceil(N_PSDU/51); shortening bits spread over codewords, not sent | - | [27]; [30] | S |
| PHY | spreading, interleaving | each coded bit repeated S times, then a bit interleaver on blocks of 2S bits (only for S = 2, 4) | - | [27] FIG. 19 text; [28] | S |
| PHY | preamble | 63-bit m-sequence + 27-bit extension = 90 bits, π/2-DBPSK, 150 µs; one of two sequences by rem(nc, 2) | bit | [30]; [27]; [33] Table 1 | S |
| PHY | PLCP header | 15 PHY bits (rate 3, reserved 1, length 8, reserved 1, burst 1, seed 1) + CRC-4 HCS (1+x+x^4) + 12 BCH parity = 31 bits; BCH(31,19), S=4, 91.9 kb/s, 206.7 µs | bit | [28] Fig. 8; [27]; [30] Table I; [33] | S |
| PHY | scrambler | side-stream, 14th-order polynomial, 1-bit seed toggled per frame | - | [28] 3.1.3; [27] | S |
| PHY | sensitivity (TI draft) | -95 / -93 / -90 / -86 at 121.4 / 242.9 / 485.7 / 971.4 kb/s; assumptions NF 10 dB, implementation loss 6 dB, 0 dBi, AWGN | dBm | [27] FIG. 32 | S |
| PHY | sensitivity (final) | 971.4 kb/s: -83 (end of the range in [28], simulation setting in [31]); other rates -92 / -90 / -87 are a guess (draft + 3 dB) contradicted at 607.1 kb/s; use only as a sweep | dBm | [28] 3.1.3; [31] Table 5 | A |
| PHY | EVM, TX | EVM 20 % (π/2-DBPSK), 12.5 % (π/4-DQPSK); TX capable of >= -10 dBm EIRP; ramps <= 5 symbols | - | [27] FIG. 31 | S |
| PHY | pSIFS | 75 (the TI draft had 50) | µs | [31] Table 5; [32] Table 4; [33] Table 1 | S |
| PHY | pMIFS / pExtraIFS | 20 / 10 | µs | [31] Table 5; [27] FIG. 26 | A (draft-era source) |
| PHY | pCCATime / pEDTime | 63 preamble symbols (105 µs) / 8 symbols (13.3 µs) | symbol | [27] FIG. 26, claim 15; [39] 3.4.1 | S |
| PHY | CSMA slot | 125 = pCCATime 105 + 20 | µs | [32] Table 4; [39] 3.4.1 | S |
| PHY | pChannelSwitchTime | 100 | µs | [27] FIG. 26 | A (draft-era source) |
| MAC | MAC header | 7 octets: frame control 32 bits, recipient 8, sender 8, BAN ID 8 | octet | [28] 2.1; [32] Table 4; [35] | S |
| MAC | FCS | 2 octets CRC-16 (x^16+x^12+x^5+1 assumed); Ullah's "18-bit" is a typo | octet | [32] Table 4; [34] | S |
| MAC | frame body | 0..255 octets (PHY Length field counts the body); PSDU <= 264 octets | octet | [28] 2.1, 3.1.2 | S |
| MAC | I-Ack | MAC header + FCS, empty body: 9-octet PSDU (193 bits with preamble and header); also N-Ack, B-Ack | - | [33] Table 1; [34] Sec. III | S |
| MAC | retry limit | no standard constant ("not specified" [34]); 3 measured best on a prototype [34]; 7 in Markov models [32]; Castalia 2 tries [35] | retx | [34] III-B | A |
| MAC | CSMA/CA windows | CWmin/CWmax UP0 16/64, UP1 16/32, UP2 8/32, UP3 8/16, UP4 4/16, UP5 4/8, UP6 2/8, UP7 1/4; backoff U[1, CW]; CW doubles only after an even number of failures; reset on success | slot | [28] Table 1; [35]; [33] | S |
| MAC | slotted Aloha | CPmax/CPmin UP0 0.125/0.0625, UP1 0.125/0.0937, UP2 0.25/0.0937, UP3 0.25/0.125, UP4 0.375/0.125, UP5 0.375/0.1875, UP6 0.5/0.1875, UP7 1/0.25; CP halves after even failures, floor CPmin | - | [28] Table 1 | S |
| MAC | superframe | beacon period of <= 255 allocation slots; EAP1, RAP1, MAP, EAP2, RAP2, MAP, CAP; scheduled 1- or m-periodic and polled allocations | - | [28] 2.2; [29] III-B | S |
| MAC | allocation slot | pAllocationSlotMin + L x pAllocationSlotResolution, 500 µs + L x 500 µs | µs | [31] Table 5 | A (single source) |
| MAC | mClockResolution / mTimeout | 4 / 30 (Castalia's 30 µs pTIFS is another parameter) | µs | [31] Table 5 | A (single source) |
| MAC | mMaxBANSize | 64 | node | [32] Sec. 3 | S |
| MAC | Castalia defaults | 32 slots of 10 ms, RAP1 8 slots; draft-era contention slot 0.36 ms | - | [35] BaselineBANMac.ned | S |
| MAC | security | levels 0/1/2; levels 1-2 add 1 freshness octet + 32-bit MIC; off in both variants | - | [28] 2.1 | S |
| MAC | coexistence tools | beacon shifting, channel hopping, active superframe interleaving (hub level) | - | [31] intro | S |

Airtime (FULL): T = (90 + 124 + N_sym) / 600 kHz, with N_PSDU = 8 (9 + L), N_CW = ceil(N_PSDU/51), N_coded = N_PSDU + 12 N_CW, N_sym = ceil(N_coded S / log2 M). Check values (µs) for L = 0 / 20 / 100 / 255:

| rate | L=0 | L=20 | L=100 | L=255 |
|---|---|---|---|---|
| 121.4 kb/s | 996.7 | 2303.3 | 7610.0 | 17796.7 |
| 242.9 kb/s | 676.7 | 1330.0 | 3983.3 | 9076.7 |
| 485.7 kb/s | 516.7 | 843.3 | 2170.0 | 4716.7 |
| 971.4 kb/s | 436.7 | 600.0 | 1263.3 | 2536.7 |
| FAIR π/2-DBPSK, S=1, uncoded | 301.7 | 568.3 | 1635.0 | 3701.7 |
| FAIR option π/4-DQPSK | 241.7 | 375.0 | 908.3 | 1941.7 |

FAIR airtime: T = (90 + 19 + ceil(8 (9 + L) / log2 M)) / 600 kHz (header 15 + 4 HCS bits at S=1, no BCH). The I-Ack is the L = 0 entry at the data rate (A).

Simulation recipe, common: one hub + N_dev sensors; superframe of N_slot = 32 allocation slots of 10 ms (Castalia default), beacon at slot 0 (17-octet body, A), random superframe phase per BAN, no inter-BAN sync; one 1-periodic scheduled uplink allocation per sensor sized for one exchange plus one spare exchange per frame or one frame (FAIR); the pooled reference stack of holds one exchange per frame and a hub-polled retry pool after the allocations (DESIGN 4.3). IFS: pSIFS 75, pMIFS 20 (burst off), mTimeout 30 µs.

FULL: BCH on header and PSDU; default 121.4 kb/s, sweep the other three rates. Exchange = DATA, pSIFS, I-Ack; success needs both decoded; a missing ACK start by t_data_end + pSIFS + mTimeout fails the attempt (A); each attempt costs T_data + 2 pSIFS + T_ack; up to R = 3 retransmissions (sweep 0/2/3/7), continuing in the next superframe's allocation; drop after R + 1 attempts. Channel: the hub averages energy on nc = 0..78 over one beacon period and takes the minimum (A: the standard gives the ED primitive, not the algorithm); optional re-scan when BAN PER > 10 % over 10 superframes (A). Optional FULL flag: beacon shifting / superframe interleaving (hub level, A). Frames below the rate's sensitivity (TI draft values, consistent with NF 10 + IL 6) are lost.

FAIR: π/2-DBPSK S=1 without BCH (600 kb/s raw), π/4-DQPSK as option; no ACK, one attempt; same superframe; channel static uniform from 0..78 per BAN; F6 option: blind pseudo-random hop per superframe over 79 channels; strip option: common channel. Packet OK iff all 19 header bits and all PSDU bits are correct.

BER model:
- γ_s = SINR in the 600 kHz symbol bandwidth; interference as Gaussian over the overlap. Adjacent channel (|Δnc| = 1): ACR 0 dB worst case, 20 dB sensitivity case (A); |Δnc| >= 2 ignored (A).
- π/2-DBPSK: p = 0.5 exp(-γb), γb = γs.
- π/4-DQPSK (Gray, differential): p = Q1(a,b) - 0.5 I0(ab) exp(-(a^2+b^2)/2), a = sqrt(γb (2 - sqrt 2)), b = sqrt(γb (2 + sqrt 2)), γb = γs/2 [12]; use scipy.special.i0e for stability.
- Repetition S (corrected): S-branch post-detection combining of DBPSK, Pb = 2^(1-2L) e^(-γt) Σ_{k=0..L-1} c_k γt^k, c_k = (1/k!) Σ_{n=0..L-1-k} C(2L-1, n), with L = S and γt = S γb [10]. This replaces "p(S γb) minus 1 dB".
- BCH t = 2: P_cw(n) = 1 - Σ_{i=0..2} C(n,i) p_c^i (1-p_c)^(n-i); PSDU codewords n_i = k_i + 12; header n = 31 at S=4.
- FULL packet success = (1 - P_cw,hdr) Π_i (1 - P_cw,i). FAIR = (1 - p)^(19 + 8 (9 + L)).
- NF + L_impl: NF 10 dB + IL 6 dB reproduce the TI draft sensitivities [27]; use them as the default link calibration.

### 3.5 ETSI SmartBAN (TS 103 326 PHY, TS 103 325 MAC)

| group | parameter | value | unit | source | mark |
|---|---|---|---|---|---|
| PHY | band, channels | 2401-2481 MHz; fc = 2402 + 2n, n = 0..39; control channels n = 0, 12, 39 (2402, 2426, 2480); data channels n = 1..11 -> DCH 0..10 (2404..2424), n = 13..38 -> DCH 11..36 (2428..2478): the BLE data-channel grid | MHz | [42] 5.1, Table 1 | V |
| PHY | channels per BAN | one CCH and one DCH at a time, chosen by the hub; no hopping | - | [43] 5.0, 5.2; [49] II-A; [51] IV | V |
| PHY | modulation | GFSK, BT 0.5, h 0.5, 1 Msym/s on CCH and DCH | - | [42] 7.2 | V |
| PHY | information rates | uncoded 1.0 / 0.5 / 0.25 Mb/s and BCH 0.89 / 0.44 / 0.22 Mb/s for Nrep 1 / 2 / 4 | Mb/s | [42] Table 3 | V |
| PHY | PPDU (V1.2.1) | preamble 16 + sync 32 + PLCP header 40 + PSDU + PHY CRC-16 = 104 bits fixed overhead (V1.1.1: 16 + 40, no sync, no PHY CRC) | bit | [42] 6.1.0 Figure 1 | V |
| PHY | preamble / sync | 1010101010101010 / 10000111101100101000011110110010 | - | [42] 6.1.1, 6.1.2 | V |
| PHY | PLCP header | length 15, PHY scheme 4 (b0b1 FEC, b2b3 repetition), reserved 3, BCH parity 14, CRC-4 1+x+x^4 over 36 bits; unit of the length field not stated | bit | [42] 6.1.3 | V |
| PHY | header FEC | BCH(36,22,t=2), shortened from BCH(127,113); "shall" in V1.1.1, "should" in V1.2.1 | - | [42] 7.3.3 | V |
| PHY | PSDU FEC | optional BCH(127,113,t=2), g(x) = x^14+x^9+x^8+x^6+x^5+x^4+x^2+x+1; last codeword shortened to (L_MPDU mod 113) + 14 bits | - | [42] 7.3.2 | V |
| PHY | repetition | whole PPDU x2 or x4, treated as one PPDU | - | [42] 7.3.1 | V |
| PHY | PHY CRC | 16 bits over the PSDU, x^16+x^12+x^5+1 | bit | [42] 6.1.5 | V |
| PHY | scrambler | optional, 1 + x^14 + x^15, init 000100100001010 ("shall" in V1.1.1) | - | [42] 7.4 | V |
| PHY | CCA | modes 1-3; ED threshold -75 dBm; detection time <= TMUA; only with MCA | dBm | [42] 8.2 | V |
| PHY | slot formula | TACK = (Lpre + LPLCP + Lhdr + Lparity)/Rsym = 120 µs by eq. (6), 168 µs with the full V1.2.1 PPDU; TTX,max = (TS - TMUA - TACK - 2 TIFS)/Nrep | - | [42] 8.1 eq. 6-10 | V |
| PHY | max frame body per slot | full V1.2.1 accounting, TMUA 0, Nrep 1, uncoded/BCH: LSLOT 1: 0; 2: 76/66; 4: 233/205; 8: 545/484; 16: 1170/1039; 32: 2420/2152 octets. LSLOT 4, Nrep 2: 106/92; Nrep 4: 42/35. Standard accounting: LSLOT 1 about 10 octets with the 48-bit TS header (8 octets with TR 103 711's 7-octet header), LSLOT 4 245/217 | octet | computed from [42] eq. 6-10, reproduced by the verifier | A (corrected note) |
| PHY | sensitivity | none normative | - | [42] | V |
| PHY | informative sensitivity | FER 10 %, 255-octet PSDU, NF 13 dB, IL 6 dB: uncoded Nrep 1/2/4 -86.2 / -92.2 / -98.2; BCH -88.1 / -94.3 / -100.2; Eb/N0 8.8 / 5.8 / 2.8 and 7.4 / 4.3 / 1.4 dB. The repetition rows count the gain twice (per-copy Eb/N0 with the information rate): with the 1 Msym/s channel rate uncoded Nrep 2/4 give -89.2 / -92.2 | dBm | [44] 8.2 eq. 39, Table 13; [48] Table 3 | V |
| PHY | fading Eb/N0 (CM3, 45 cm, 50-octet body, BCH) | FER 10 %: Nrep 1 16 dB without / 9.4 with 1 retx; Nrep 2 8.8 / 5.0; Nrep 4 3.8 / 1.5. FER 1 %: 26.3 / 15.3; 13.9 / 8.6; 7.3 / 3.7 | dB | [44] Table 14 | V |
| PHY | ETSI BER model | Pb = Q(sqrt(2 Eb/N0)) (coherent MLSD, EGC of copies [44] eq. 38); Pb = Q(sqrt(2 ε Eb/N0)), ε = 0.78, Eb/N0 = SNR + 3 dB with P_N = -92.5 dBm (a noise floor, not a reception threshold), TX -10..0 dBm [45] | - | [44] 7.6; [45] 7.2, Table 15 | V |
| PHY | alternative link figures | NF 10 dB, IL 5 dB, TX max 0 dBm; TX -10 dBm in STF 511 | dB | [50] Table 2; [44] 7.4.2 | S / V |
| MAC | topology | star, 1 hub + up to 16 nodes; hub NID 00010101, broadcast 11111111 | - | [43] 5.0, Table 5 | V |
| MAC | C-Beacon | on the CCH every TC (TC "variable"); body 120 bits, MPDU 184 bits, 288 µs uncoded; carries DCH number, interference-mitigation bit, destination channel | - | [43] 5.2.1, 6.2.1, 8 | V |
| MAC | D-Beacon | on the DCH in slot 0 of each IBI; body 104 bits (+40 optional incl. channel migration fields), MPDU 168 bits, 272 µs uncoded | - | [43] 5.2.2, 6.2.2 | V |
| MAC | IBI | TD = LD TS: beacon slot, scheduled TDMA slots, C/M slots (slotted ALOHA), inactive period; 10-bit slot index | - | [43] 5.2.2 | V |
| MAC | slot | TS = LSLOT x 625 µs, LSLOT in {1, 2, 4, 8, 16, 32} (0.625..20 ms) | µs | [43] Table 9, 8 | V |
| MAC | TIFS / Tsch | 150 / 16 | µs | [43] 8 | V |
| MAC | slot structure | data, IFS, ACK/NACK, >= IFS before slot end; MCA slot starts with TMUA sensing | - | [43] 5.2.2.1-5.2.2.3 | V |
| MAC | TMUA | UP3 emergency 0; slot owner Tsch; non-owner in scheduled period 2 Tsch; other C/M traffic Tsch | µs | [43] Table 2 | V |
| MAC | frame | MAC header 48 bits (FC 24, recipient 8, sender 8, BAN ID 8) + FCS-8 (x^8+x^7+x^3+x^2+1) + body + CRC-16: MPDU = 64 + LF bits; sequence number kept on retransmission | bit | [43] 6.1 | V |
| MAC | ACK | ACK/NACK with empty body (MPDU 64 bits); policy 0 = ACK on success, policy 1 = NACK on failure only; 168 µs uncoded (derived) | - | [43] 6.1.1.1.2, 6.3.1 | V |
| MAC | slotted ALOHA CP | UP0 1/8..1/16, UP1 1/4..1/16, UP2 1/2..1/8, UP3 1..1/2; keep CP after odd failures, halve after even failures while CP >= 2 CPmin | - | [43] Tables 3-4, 7.3.2.1; [52] | V |
| MAC | MCA backoff | BE from BEmin per packet; BC uniform [0, 2^BE - 1] slots; busy: BE++ up to BEmax; RT limit RTmax; PCA (UP3): BE = 0 first; RCA: TMUA = 2 Tsch. The slot author's scheduled data is exempt: it senses Tsch (for PCA) and, if busy, waits for its next own slot | - | [43] 7.3.3, Annex A | V (corrected) |
| MAC | BEmin / BEmax / RTmax | not in the TS; 3 / 5 / 3 (802.15.4 defaults) | - | [43] 8 (absence); [13] Table 86 | A |
| MAC | retry limit | none normative; STF 511 simulator retransmits once, no combining; the only hardware implementation resends in each following IBI until acked | retx | [44] 7.5; [48]; [51] IV | V (as simulator setting) |
| MAC | channel choice | hub picks a BAN ID and DCH unused by neighbours; migration via D-Beacon fields (switch at a given D-Beacon sequence number); the decision rule is not normative | - | [43] 7.2.1, 6.2.2.10, 7.6; [49] | V / S |
| MAC | neighbour discovery | hub alternates inactive-period lengths and listens on the CCHs to read neighbour C-Beacons (timer not given) | - | [46] 6.2 | V |
| MAC | PHY scheme | static: node capability in C-Req, hub assigns FEC and repetition in C-Ass; no run-time adaptation or TX power control defined | - | [43] Tables 12, 16; [47] 4.8 | V (corrected: F11 partial) |
| MAC | use-case QoS | latency 10 ms, PRR > 90 %; example IBIs 5.6 / 11.2 / 22.4 ms | - | [45] Tables 14, 15 | V |
| impl. | silicon implementation | nRF54L15, LSLOT 16 (10 ms), IBI 220 ms (1 + 16 + 2 C/M + 3 inactive), fixed DCH, one C-Beacon per IBI after the C/M period, IBI drift 23 ppm, slot jitter P99 < 754 µs; BCH applied to all frames (ACK PSDU 78 bits, C-Beacon 212 bits); 10 missed D-Beacons trigger a rescan | - | [51] Table V, Sec. IV-V | S |

Simulation recipe, common: 37 DCHs on the BLE data-channel indices; PPDU airtime Nrep (104 + L_PSDU) µs with L_PSDU = 64 + 8 LF (uncoded) or + 14 ceil(L_MPDU/113) with BCH. LSLOT = smallest value with Nrep (104 + L_PSDU) <= TS - TMUA - 168 - 300. IBI: slot 0 D-Beacon (272 µs, uncoded, Nrep 1), slots 1..N_dev one uplink slot per node, 2 C/M slots, NI inactive slots with TD = reporting period; fragment when the payload exceeds LF,max. Random IBI phase per BAN, clock error uniform +/-20 ppm. C-Beacons (288 µs, TC = TD, A) only in joining or mixed runs. In FULL the ACK starts 150 µs after the data end and lasts 168 µs; D-Beacons and ACKs collide like data.

FULL: header BCH(36,22) and PSDU BCH(127,113), Nrep 1 (0.89 Mb/s). ACK policy 0; a missing ACK leads to a retransmission in the same IBI's C/M slots (slotted ALOHA with the UP's CP), else in the node's next scheduled slot ahead of new data; RT = 1 (sweep 3; option "until ACK" as in [51]). MCA flag (default on, needs every node and the D-Beacon bit): the owner senses 16 µs at slot start, and if a same-DCH signal above -75 dBm is on air it skips to its next own slot. Channel: a new hub listens to the CCHs for one scan window >= TD, learns all active DCHs from C-Beacons and picks a free DCH uniformly (least-used if all 37 are taken). Migration (A): uplink failure ratio > 20 % over 10 IBIs -> free DCH or lowest measured energy, switch at D-Beacon sequence + 2. F11 flag (off, beyond the standard): raise Nrep 1 -> 2 -> 4 when node PER > 10 % over 20 IBIs, lower when < 1 %.

FAIR: uncoded, Nrep 1 (1 Mb/s); the 40-bit header is kept but must be error-free (CRC-4 kept); FCS-8 and both CRC-16s kept; same IBI and slot layout, the ACK part of each slot stays idle; no retransmission, no C/M recovery, no MCA (TMUA = 0); DCH uniform from 37 at start, never changed.

BER model:
- Common GFSK curve (default, same as BLE for fairness): Pb = 0.5 exp(-γ/2), γ per 1 MHz channel bit, NF + L_impl calibrated to the silicon (8.1 dB for an nRF52840-class radio).
- ETSI mode (to reproduce [44]): γ_b = P_rx + 174 - NF - 60 - IL dB with NF 13, IL 6 (γ_b = P_rx + 95 dB); Pb = Q(sqrt(2 ε γ_b)), ε = 1 ([44]) or 0.78 ([45]).
- Repetition: ideal EGC, γ_comb = Nrep γ_b (reproduces 8.8 / 5.8 / 2.8 dB); in collision mode each copy is decoded alone: PER = Π PER_copy.
- Success = P_sync P_hdr P_psdu (1 - p)^16. P_sync = 1, or (1-p)^32 + 32 p (1-p)^31 [50]. FULL header: [Σ_{i=0..2} C(36,i) p^i (1-p)^(36-i)] (1-p)^4. FAIR header: (1-p)^40. FAIR PSDU: (1-p)^(64 + LF). FULL PSDU: product over m = ceil(L_MPDU/113) codewords (n = 127, last n' = (L_MPDU mod 113) + 14) of Σ_{i=0..2} C(n,i) p^i (1-p)^(n-i).
- Calibration: independent-error BCH gives FER 10 % at 5.8 dB, ETSI MLSD 7.4 dB (burst errors). Subtract 1.1 dB (Eb per information bit) to 1.6 dB (per channel bit) from γ for coded frames, or use the pair-error model (events at Pb/2, 2 bit flips each, t = 1 event per codeword, 6.2 dB).
- SINR mode: no C/I or ACR in the standard; use the BLE spec values (co-channel 21 dB, 1 MHz 15, 2 MHz -17) as the analogue (A).

### 3.6 ANT / ANT+ (with Nordic ESB as an optional comparator)

| group | parameter | value | unit | source | mark |
|---|---|---|---|---|---|
| PHY | RF field | 0..124 -> f = 2400 + RF MHz (protocol field range); default RF 66 = 2466 MHz; a channel stays on one frequency | MHz | [53] 5.2.2, 9.5.2.6 | V |
| PHY | usable set | nRF24AP2: 2403-2480 MHz typ (78 frequencies), band edge 2483.5; non-ANT+ networks avoid 2450 and 2457 MHz | MHz | [56] Table 9; [53] p.16 | V (corrected: not 125) |
| PHY | ANT+ frequency | 2457 MHz (RF 57, HRM); many ANT+ devices on 2450 or 2457 | MHz | [53] 5.2.2; [54] Table 1 | V |
| PHY | modulation | GFSK 1 Mb/s, +/-160 kHz deviation (h = 0.32, derived); nRF52832 nRF mode +/-170 kHz (h = 0.34) | - | [56] Table 9; [57] p.2; [59] 23.15 | V |
| PHY | sensitivity | nRF24AP2 -85 dBm (BER 0.1 %); CC2570 -86 dBm (0.1 %); nRF52832 1 Msps nRF mode -93 dBm (BER not stated; 3 dB worse on ADDR1..7) | dBm | [56]; [57]; [59] 23.15.5 fn 16 | V |
| PHY | C/I, ANT interferer | nRF24AP2 equal modulation, wanted -67 dBm: co 12; +/-1 MHz 8; 2 MHz -21; 3 MHz -30; > 6 MHz -40; > 25 MHz -50 dBc. ETSI table: co 9, 1 MHz 8, 2 MHz -20, 3 MHz -30, > 6 MHz -40, > 25 MHz -47 | dB | [56] p.41 | V (corrected) |
| PHY | C/I nRF52832 1 Msps | co 9; -1 MHz -2; +1 MHz -10; +/-2 MHz -19 / -42 | dB | [59] 23.15.5 | V |
| PHY | TX spectrum | nRF24AP2 20 dB bandwidth 950-1100 kHz, adjacent TX power -20 dBc (1 MHz), -40 dBc (2 MHz); nRF52832 -25 / -50 dBc | - | [56]; [59] | V |
| frame | payload | 8 octets per message | octet | [53] 5.6 | V |
| frame | on-air frame | preamble 1 + address 5 + payload field 10 (2 control + 8 data) + CRC 2 = 18 octets = 144 µs; CRC-16 0x1021, init 0xFFFF over address + payload; no whitening | octet | [60] packets.h, shockburst.c; [63] "< 150 µs" | S |
| frame | ACK / reverse frame | 144 µs (full 18-octet frame, static length) | µs | legacy ShockBurst format | A |
| timing | channel period | T = N / 32768 s, 16-bit N, default 8192 (4 Hz); 0.5 Hz to above 200 Hz; nRF24AP2 5.2 ms..2 s (N >= 171) | - | [53] 5.2.4; [56] p.1 | V |
| timing | ANT+ HRM period | 8070 counts (4.06 Hz); receivers may use 16140 or 32280 (minimum 1.02 Hz) | count | [54] 4.1 | V |
| timing | broadcast/ack ceiling | 8 octets x max rate: 12.8 kb/s at 200 Hz (CC2570), 12.3 kb/s on nRF24AP2 (5.2 ms) | kb/s | derived | V |
| timing | burst | 20 kb/s max; packets timed from each other, ~3 ms each (3.2 ms derived); advanced burst 20 / 40 / 60 kb/s with 8 / 16 / 24-octet packets | - | [53] 5.4.3, 5.4.4, 9.5.7.8 | V |
| timing | turnaround | 130 µs on both chip profiles (nRF24 Tstby2a; nRF52 ESB stack default 130, fast 40; ANT timing must match nRF24AP2 peers) | µs | [58] 7.7; [62] esb.h | A (corrected) |
| timing | master RX window W | 144 + 2 x 16 = 176 µs after the turnaround | µs | none | A |
| ARQ | acknowledged data | ACK in the same slot; no automatic retransmission; host may resend next slot; without new data the old payload goes as broadcast | - | [53] 5.4.2, 9.5.5.2 | V |
| ARQ | burst retries | each packet retried up to 5 times (1 + 5 attempts), transfer aborted after that; single-packet burst has no retries; bursts pre-empt the node's other channels | - | [53] 5.4.3 | V |
| ARQ | advanced burst | retry extension in cycles of 5, default 20 retries; stall default 3210 (~10 s); optional 6-frequency hopping per data/ACK pair (must be requested) | - | [53] 9.5.7.x, 9.5.7.7 | V |
| access | open-time search window | the master listens before its first transmission to avoid overlap; "fast channel initiation" (extended assignment bit 4) skips it | - | [53] 5.1, 5.2.1.4.3 | V |
| access | adaptive isochronous | guard window after the own transmission window detects intruders; phase changed "by any amount"; no numbers; the word does not occur in the protocol documents | - | [61] | V (qualitative) |
| access | EVENT_CHANNEL_COLLISION | two channels of one device drift into each other; one is blocked | - | [53] 9.5.6.1 | V |
| access | frequency agility | 3 frequencies (default RF 3/39/75 = 2403/2439/2475 MHz), switch only on degradation; not on shared or Tx/Rx-only channels | - | [53] 5.2.1.4.1, 9.5.2.21 | V |
| access | agility trigger | master after 8 unacked messages, slaves after 8 s without reception | - | forum thread, search summary only | A |
| topology | channels per device | nRF24AP2-8CH 8; nRF24AP1 4; S212 SoftDevice 15 | - | [56]; [53] 5.5; [64] | V |
| topology | shared channel | 1- or 2-octet address in the payload; a slave replies only when addressed | - | [53] 5.6 | V |
| search | timeouts | low priority 5 s (2 x 2.5 s), high priority 25 s (10 x 2.5 s), 255 infinite | s | [53] 9.5.2.5, 9.5.2.15; [55] 6.2 | V |
| search | drop to search | < 2 Hz: after 4 missed messages; > 2 Hz: after 2 s of misses; high-priority search loses up to 50 % of other channels' messages | - | [55] 3, 6.1 | V |
| ESB | frame | preamble 1, address 3..5, PCF 9 bits, payload 0..32 (nRF5 1..252), CRC 1 or 2 octets | - | [58] 7.3; [62] Kconfig | V / S |
| ESB | timing | T_OA = (8 (1 + A + N + C) + 9) / Rb; 137 µs (A 5, N 8, C 2, 1 Mb/s), empty ACK 73 µs; Tstby2a 130 µs | µs | [58] 7.7 Table 19 | V |
| ESB | ARQ | ARC 0..15, ARD 250..4000 µs in 250 µs steps from the end of the packet; PTX gives up after 250 µs without address match | - | [58] 7.4.2 | V |
| ESB | NCS defaults | 2 Mb/s, CRC-16, ARD 600 µs, ARC 3, payload 32 | - | [62] esb.h | S |
| ESB | C/I (nRF24L01+, modulated carrier) | 1 Mb/s co 12, 1 MHz 8, 2 MHz -21; 2 Mb/s co 11, 2 MHz 4, 4 MHz -18 dBc; sensitivity -85 (1M), -82 (2M), -94 (250k) dBm | dB | [58] p.16, Table 7 | V |

Simulation recipe, common:
- Every frame 144 µs; advanced-burst frames 8 (1 + 5 + L + 2 + 2) µs = 144 / 208 / 272 µs for L = 8 / 16 / 24 (A).
- Each sensor opens one bidirectional master channel; the hub opens N_dev slave channels (N_dev <= 8 on nRF24AP2, <= 15 on S212).
- Rate mapping: message rate f = R / 64 bit; N = round(32768 / f), clamped to 171..65535 (nRF24AP2) or 164..65535; true period T_i = N / 32768 (1 + ε_i), ε_i uniform +/-50 ppm per device (the datasheet crystal requirement); initial phase uniform in [0, T_i).
- Above the broadcast ceiling: a burst every slot, B = ceil(R T_i / 64) packets, one data + ACK pair every 3.2 ms.
- Slot: master TX [t0, t0 + 144]; slave reply or ACK [t0 + 274, t0 + 418]. On air: 144 µs broadcast, 418 µs acknowledged.
- Frequencies: default all BANs on 2457 MHz (ANT+); private-network option: one random fixed frequency per BAN from the 78 usable (2450, 2457 excluded).
- Hub radio constraint: two slave slots of the hub that overlap in time lose the later one (EVENT_CHANNEL_COLLISION). Burst transactions pre-empt the hub's other channels.
- Slave loses sync after the drop-to-search rule, then searches 5 s low priority + 25 s high priority (during high priority search the hub's other channels lose each slot with p = 0.5, worst case) and reacquires on the first good master frame.

FULL:
1. Channel open (F4): the master listens for one period and picks a phase whose occupancy, widened by 500 µs each side, overlaps nothing heard; if none, a random phase.
2. Data (F3), R <= ceiling: acknowledged messages at f_ch = 1.25 f_data (A); FIFO; an unacked payload is resent in the next slot; drop after 6 attempts (1 + 5, mirrors the burst rule). Empty FIFO: old payload as broadcast (energy only).
3. R > ceiling: burst, 1 + 5 attempts per packet; advanced burst 20 retries and, if requested, 6-frequency hopping.
4. Flags, default off: adaptive isochronous re-timing (patent [61]: sense a guard window of 500 µs after the own slot, shift 250 µs per period away from an intruder, A); frequency agility (3 frequencies, switch after 8 unacked messages, A).
5. FULL-ANT+ sub-case (how HRM-class sensors ship): broadcast at f_data with step 1 only.

FAIR: broadcast only at f = f_data, no ACK, no host retry; burst packets sent once (spacing 3.2 ms); random initial phase, no listen (a legal stack setting: fast channel initiation), phases drift freely; fixed frequency 2457 MHz (or random fixed); keep CRC, periodic slots, advanced-burst hopping unless stripped. Two same-period masters that meet stay collided for about 2 T_pkt / (|ε_i - ε_j| T) periods (29 periods, about 7 s, for 40 ppm at 4.06 Hz).

ESB (optional config of the same engine, low priority): unslotted PTX per sensor, one fixed channel per BAN, T_OA and T_ACK as above, retransmit at t + T_OA + ARD. FULL: ARC 3, ARD 600 µs, 2 Mb/s (also 1 Mb/s, and ARC 15). FAIR: NO_ACK (ARC 0).

BER model:
- γ = SNR per bit in B = Rb. Default common GFSK curve Pb = 0.5 exp(-γ/2) with NF + L_impl fitted at BER 1e-3: nRF24AP2 (-85 dBm) 18.1 dB, nRF52832 1 Msps (-93 dBm) 10.1 dB, CC2570 (-86 dBm) 17.1 dB.
- Reference (conservative) curve: non-coherent CPFSK with |ρ| = sin(πh)/(πh): 1e-3 at 17.80 dB (h 0.32), 17.30 dB (h 0.34), 14.31 dB (h 0.5).
- PER = 1 - (1 - Pb)^Nb, Nb = 136 (ANT: address 40 + payload field 80 + CRC 16); ESB Nb = 8 (A + N + C) + 9 (129 for 5/8/2). Pb 1e-3 -> 12.7 %, 1e-4 -> 1.35 %, 1e-5 -> 0.14 %.
- No FEC, no repetition: FULL and FAIR share the per-attempt PER. ARQ success needs data and ACK: p_s = (1 - PER_d)(1 - PER_ack).
- Capture (SINR mode): survive if S/I >= 12 dB co-channel (nRF24AP2 equal modulation); 1 MHz apart lost if S/I < 8 dB; channels >= 2 MHz apart are not independent when the interferer is more than 21 dB stronger (apply the full C/I table).

### 3.7 IEEE 802.15.6 HBC (standard stack)

The HBC clause of IEEE 802.15.6-2012 [26] was not readable; all values are secondary or assumed. The MAC is the common 802.15.6 MAC of section 3.4.

| group | parameter | value | unit | source | mark |
|---|---|---|---|---|---|
| PHY | channel | single channel at 21 MHz, 5.25 MHz wide (3 dB band 18.375-23.625 MHz); the draft's 16 / 27 MHz bands in [29] are superseded | MHz | [37] p.440; [38] Table 1; [41] | S |
| PHY | chip rate | 42 Mcps (baseband digital transmission, no RF carrier) | Mcps | [37]; [38]; [41] | S |
| PHY | symbol mapping | 4 bits -> 16-chip Walsh code, then frequency-selective spreading with an FSC of length SF | - | [37]; [38] Fig. 2; [36] 4.2.4 | S |
| PHY | rates | 164 / 328 / 656 / 1312.5 kb/s; symbol rate 41 / 82 / 164 / 328 ksym/s; FSC SF 64 / 32 / 16 / 8; Walsh-chip rate 0.656 / 1.3125 / 2.625 / 5.25 Mcps | kb/s | [38] Fig. 2 | S (corrected) |
| PHY | preamble | 64-bit Gold code x4, FSC SF 8 = 2048 chips (48.76 µs) | chip | [28] 3.2; [41] | S |
| PHY | SFD / RI | 64-bit Gold code, SF 8 = 512 chips (12.19 µs); optional rate indicator by SFD time offset (RI mode sends the header at the PSDU rate) | chip | [28] 3.2; [36] 4.2.2 | S |
| PHY | PLCP header | 32 bits with CRC-8 (rate, pilot info, burst, scrambler seed, length 8, CRC 8) | bit | [28] 3.2 | S |
| PHY | header rate | 164 kb/s (SF 64, 8192 chips, 195 µs) outside RI mode, from the 2010 proposal; [28] only says "spread at 42 Mcps" | - | [36] 4.2.2 | A |
| PHY | header CRC-8 polynomial | not confirmed: 0x8B (x^8+x^7+x^3+x+1) in [40], 1+x^2+x^3+x^7+x^8 in another summary; the earlier code's 0xD5 is its own choice | - | [40] p.603 | A |
| PHY | pilot | SFD sequence inserted periodically; proposal field 3 bits (000..110: 0.125..4 ms or none; 64 B / 128 B / none); [28] says 2 bits | - | [36] Table 10; [28] | S |
| PHY | scrambler | 32nd-order polynomial (the code's taps z[n-11], z[n-31], z[n-32] are not confirmed) | - | [28] 3.2 | S |
| PHY | FEC | none found in any HBC transmitter chain | - | [38] Fig. 2; [37]; [36]; [28] | A |
| PHY | min symbol distance | 64 at 1312.5 kb/s (128-chip symbol) | chip | [38] Table 1 | S |
| MAC | frame | MAC header 7 + body 0..255 + FCS 2 octets; I-Ack 9-octet PSDU | octet | [28] Fig. 2; [32] | S |
| MAC | random access | prioritized slotted Aloha for HBC (CSMA/CA is for NB and UWB); CP tables of section 3.4 | - | [39] 3.4.1; [28] Table 1 | S (corrected) |
| MAC | pSIFS, slot timing | no HBC values found; NB values reused (pSIFS 75 µs, allocation slot 500 µs + L x 500 µs) | µs | - | A |
| MAC | retry limit | 3 retransmissions (as NB); 7 in [32] is a model assumption; sweep 0/3/7 | retx | [34]; [32] | A |
| energy | HBC transceiver ICs | Cho 2015: 5.2 mW (4.3 mA at 1.2 V), 0.13 µm, 802.15.6 HBC compatible; Lee 2013: 5.5 mW; no TX/RX split | mW | [75]; [74] | S |

Frame length in chips at 42 Mcps: 2048 (preamble) + 512 (SFD) + 8192 (header at SF 64; 1024 in RI mode) + N_s 16 SF + floor(N_s / 128) x 512 (pilots, 128-symbol period assumed from the code's "010"), N_s = ceil(8 (9 + L) / 4). Check values (µs):

| SF (rate) | L=0 | L=20 | L=255 |
|---|---|---|---|
| 8 (1312.5 kb/s) | 310.9 (RI mode 142.5: SFD/RI field 608 chips, 10.5.2) | 432.8 | 1913.9 |
| 16 (656.25 kb/s) | 365.7 | 609.5 | 3523.0 |
| 32 (328.125 kb/s) | 475.4 | 963.0 | 6741.3 |
| 64 (164.0625 kb/s) | 694.9 | 1670.1 | 13177.9 |

The same frame expressed in Walsh-chip times of the BCC comparison (Rc plays the 5.25 Mcps Walsh-chip rate of the 1312.5 kb/s mode): 256 + 64 + 1024 (header 8x longer; RI mode 128) + 32 (L + 9) S + 64 floor(32 (L + 9) S / 2048), S = 1, 2, 4, 8 for 1312.5 / 656 / 328 / 164 kb/s.

Simulation recipe:
- Single 21 MHz channel, no F5, no F6. Scheduled allocation per sensor per superframe (as NB), optional slotted-Aloha contention phase.
- Inter-body coupling applies: the HBC standard is capacitive and all BANs share one channel, so FULL and FAIR HBC PSR depend on k (section 5.7); the BCC-galvanic isolation assumption does not carry over.
- FULL: SF 8 (sweep 16/32/64); header at 164 kb/s (flag: RI mode, header at the PSDU rate); FCS detects every payload error; I-Ack pSIFS after the data at the data rate (A), timeout pSIFS + T_ack + 1 slot (A); up to 3 retransmissions; rate rule (A): the highest rate whose predicted PSR >= 0.99 at the link γ.
- FAIR: SF 1 everywhere, header at the PSDU rate, MAC header + FCS kept, no ACK, scheduled allocations, fixed rate.
- Error model: the BCC Walsh and sync formulas of section 3.8, with γ_hdr = γ + 9.03 dB when the header runs at SF 64 and γ_psdu = γ + 10 log10(S) at the lower rates. FSC spreading at equal Walsh-chip energy gives no further AWGN gain (spectral shaping only).

### 3.7a Checked against IEEE 802.15.6-2012 Clause 10

Source: Clause 10 "Human body communications PHY specification", pages 235-249. Clause 10 is the PHY only; the MAC (acknowledgement policies, polls, improvised
access, MAC header, retry limits) is in Clauses 6 and 7, which were not available. (The MAC was checked later against the full standard: section 3.7b.) The rows of 3.7 marked S or A that
this clause settles are listed here; 3.7 itself is left as it was so the audit trail stays readable.

| item | Clause 10 | simulator / earlier code | result |
|---|---|---|---|
| band | centre 21 MHz, fBW 5.25 MHz (10.1, 10.8.1) | HBC: same; galvanic BCC: BFSK at 12.5-13.3 MHz | HBC confirmed. 21 -/+ 2.6 MHz lies above the measured 0.05-20 MHz galvanic band, so the HBC PHY itself cannot be evaluated on earlier work channel; the galvanic series run the HBC framing and MAC over earlier work BFSK PHY |
| clock, spreading | fCK 42 MHz; 4 bits -> 16-chip orthogonal code, then FSC; header 164 kb/s in DRF mode, PSDU 164 / 328 / 656 / 1312.5 kb/s | as 3.7 | confirmed; the header rate (3.7 mark A) is now V |
| preamble | 64-bit Gold code x4, FSC SF 8 = 2048 chips (10.4) | 2048 | confirmed |
| SFD | 64-bit Gold code, SF 8 = 512 chips (10.5.1) | 512 | confirmed |
| RI | SFD with a time offset of up to 12 zero bits: (64 + 12) x 8 = 608 chips (Eq. 117); RI lets header and PSDU run at one rate "especially only for burst packet"; the MAC sets the burst-mode bit (10.5.2) | ladder L1 and up uses RI on every frame with the 608-chip SFD/RI field (`sfd_ri_chips`, 76 Walsh chips: the model had the 512-chip SFD there, 96 chips short per frame); L0 (DRF mode) keeps 512 | confirmed as modelled; RI on single (non-burst) frames, beacons and polls is an extension of the option, not the option as written: 10.5 assigns RI to burst packets and 10.5.2 has the MAC send a non-burst packet in DRF mode |
| PHY header | rate 0-2; pilot info 3-4 (2 bits in the table, 3-bit values listed, bit 5 not listed); reserved 6-7; burst 8; reserved 9-10; scrambler seed 11; reserved 12-15 marked "PSDU length extension"; length 16-23; CRC-8 24-31 | 32-bit header as coded | confirmed; the standard's own table is inconsistent at bits 3-5. Ladder L5 widens the length to 11 bits inside bits 12-15, which the standard reserves for exactly that |
| header CRC-8 | G(x) = x^8 + x^7 + x^3 + x^2 + 1 (0x8D), register initialised to all ones, output complemented, LSB first (10.6.1, Figure 166) | `HBCProtocol.py`: 0xD5, init 0, one bit per byte (3.8) | the earlier code does not implement the standard's CRC-8; of the two candidates in 3.7 the second one is the standard's |
| scrambler | z[n] = z[n-11] xor z[n-31] xor z[n-32]; two seeds, toggled per frame (10.7.1) | same taps, toggled per frame | taps confirmed; seed values not compared |
| pilot | the SFD sequence, every 64 or 128 octets of PSDU, or none (10.7.3, Table 86) | 2048 / 4096 Walsh chips | confirmed (64 octets = 2048 Walsh chips) |
| FEC | none in the transmitter chain | none | confirmed |
| clock tolerance | +/-20 ppm (10.8.3) | L1 assumes 40 ppm relative drift | consistent |
| sensitivity | -97.35 / -94.34 / -91.33 / -88.32 dBm at 164.1 / 328.1 / 656.3 / 1312.5 kb/s for PER < 1 % at 128 octets, NF 10 dB, implementation loss 6 dB | - | implies Eb/N0 = -97.35 + 174 - 10 log10(164.1e3) - 10 - 6 = 8.5 dB per information bit at every rate (3 dB per rate step); reference for |
| timing | pSIFS 75 us, pMIFS 20 us (burst), pExtraIFS 10 us | pSIFS 75 us (reused from NB, mark A); pMIFS 20 us (parameters.json ieee802156_hbc.pmifs_us) between the frames of a burst of the lead reference (flag burst) | confirmed; now V |
| transmit limit | radiated field <= 30 uV/m at 30 m (10.8.2) | - | not modelled |

### 3.7b Checked against IEEE 802.15.6-2012 MAC, Clauses 5 and 6

Source: the whole IEEE Std 802.15.6-2012. 802.15.6-2026 does not support HBC: its introduction says the revision "does not support deprecated FM-UWB PHY and human body communications (HBC)", so 802.15.6-2012 is the last IEEE standard with a body-channel PHY.

| item | 802.15.6-2012 | galvanic profile / simulator | result |
|---|---|---|---|
| MAC header | Frame Control 4 octets, Recipient ID 1, Sender ID 1, BAN ID 1 = 7 octets, one format (5.2.1, Figure 9) | L4: 1 octet | L4 is not compliant: new |
| Ack policies | N-Ack, G-Ack (both 00), I-Ack (01), B-Ack (10), L-Ack (11) (5.2.1.1.2, 6.2.9) | L0: I-Ack; L2: acknowledgement field in the beacon | see next row |
| G-Ack | a node sends data of subtype mG-AckDataSubtype; the hub acknowledges every such frame received since its last B2 by listing the senders' NIDs in the next B2 frame; the node retries on a missing B2 or a missing NID (6.2.9.2, 5.4.8) | L2: per sensor a base sequence number and one bit per frame, in the beacon | L2 is the G-Ack option with two changes: it rides in the beacon instead of a B2 frame, and it acknowledges frames by bitmap, not senders by NID. Existing option, extended |
| B-Ack / L-Ack | block transmission, B-Ack with a frame status bitmap, retransmission of the frames it marks lost (6.2.9.4) | the lead reference (protocols flag back, L0 with the pool: data frames with L-Ack, one B-Ack of 3 octets pSIFS after the sensor's last frame of its allocation, the frames it marks lost in the hub-polled pool with their I-Ack; claim4.py T15 G closed form, H engine); in the HBC burst mode the lead reference (flags back and burst: the allocation's first frame in DRF mode, the others pMIFS apart in RI mode with the 608-chip SFD/RI, pSIFS, B-Ack; T15 G2 closed form, I engine); the lead reference's block spans its allocations, unlimited, so no B-Ack is sent (flag span; T15 G3 closed form, J engine) | an existing alternative to L2 for one sensor's frames |
| retry | "may retry ... taking into consideration delay requirements, fairness policies, channel conditions"; no retry limit for HBC (6.2.6); a limit of 4 exists only for the UWB HARQ | retry limit 3 (PARAMETERS 3.7, A) | the limit is a model assumption, the standard leaves it open |
| polled allocations | improvised access: the hub sends Poll, T-Poll, I-Ack+Poll or B-Ack+Poll and grants a polled allocation; frames in it may be new or retried (6.6.1, 6.6.1.1) | L3: retry pool of hub-polled exchanges after the allocations | L3's mechanism exists; its sizing is a hub policy the standard leaves open: compliant, new policy |
| rate selection | the PHY header carries the rate per frame; the MAC negotiates acknowledgement rates (5.3.6.6, 5.3.7.11); no rate-selection rule |: rate from the path margin, spreading per link | compliant, new policy (within the HBC PHY's rates; the unspread option is a PHY change) |
| coexistence | beacon shifting, channel hopping, active superframe interleaving, all optional (6.13) |: switch to shared-medium operation when bodies touch | the fallback can reuse beacon shifting and superframe interleaving; the trigger (a foreign beacon on a galvanic, normally isolated medium) and the switch back are new |
| BAN identification | BAN ID in every MAC header | L4 keeps it in the beacon only; FCS seeded with the BAN ID proposed (REPORT 5.4) | new, not compliant |

Summary for the claim (compliance is judged against 2012, the last edition with HBC): of the galvanic profile, L1 (RI, pilot period), L2 (as G-Ack) and L3 (polled allocations) are
existing 802.15.6 mechanisms, used with extensions; L4 and the beacon-only BAN ID are new and not compliant; L5 uses the
header bits that Clause 10 reserves as "PSDU length extension". The contributions that are new and still compliant are
policies the standard leaves open: pool sizing, rate selection, and when to switch coexistence mechanisms on. The PHY is
not the HBC PHY (21 MHz lies above the measured galvanic band, 3.7a); the profile is the 802.15.6 MAC over a galvanic PHY.

### 3.8 As-built BCC stack (`HBCProtocol.py`, `MAC_layer.py`, `PHY_layer_BCC_simulation_v2.py`)

All rows were read in the repository code [100] and cross-checked by running the generator; the feature column of section 2 comes from the code, not earlier work text.

| group | parameter | value | unit | source | mark |
|---|---|---|---|---|---|
| frame | preamble | 64-chip PLCP code x4 = 256 chips (earlier work calls it a Gold code; not checked) | chip | HBCProtocol.py:50-55, 149-150 | V |
| frame | SFD | 64 chips; the same sequence is the pilot | chip | HBCProtocol.py:43-48, 120, 220 | V |
| frame | header | 32 bits: rate 3, pilot info 3, 00, burst 1, 00, seed 1, 0000, length 8, CRC-8 8; 8 Walsh symbols = 128 chips at the PSDU chip rate (the 164 kbps multiplier has no effect with spreading off) | bit | HBCProtocol.py:152-158, 212-213 | V |
| frame | PSDU mapping | 4 bits -> 16-chip Walsh (rows of Sylvester H16, +1 -> 1); 32 chips per octet | chip | HBCProtocol.py:57-75, 171-181 | V |
| frame | pilot | 64 chips after every full 2048 PSDU chips (pilot info 010; 011 = 4096; 110 = none); all runs use 010 | chip | HBCProtocol.py:84-93, 183-202 | V |
| frame | length | N(L) = 256 + 64 + 128 + 32 L + 64 floor(32 L / 2048); N(255) = 8800; 4.314 chips per payload bit (6.35 dB) at L = 255, 5.12 at 50, 6.8 at 20 | chip | generator run; ch05 sec04:53 | V |
| frame | on-air file | `protocol_SDR_test/data_to_send` carries one extra 64-chip pilot after each packet: 8864 chips per 255-octet packet (4.345 chips/bit); 2000 idle chips between packets | chip | bit-for-bit comparison | V |
| frame | max PSDU | 255 octets (8-bit length) | octet | HBCProtocol.py:155 | V |
| frame | payload CRC | none; `unpack_packet` returns descrambled bits without a check | - | HBCProtocol.py:469-495 | V |
| frame | header CRC-8 | poly 0xD5, table driven, init 0, fed one header bit per byte: [32,24] code with dmin 2. On failure the parser re-decodes the inverted header and returns its garbage fields, so a detected header error still loses the packet; undetected errors are the real weakness | - | HBCProtocol.py:13-38, 156-157, 356-365, 477-478; exhaustive run | V (corrected) |
| frame | scrambler | additive 32-stage LFSR z[n] = z[n-11] xor z[n-31] xor z[n-32]; 2 seeds toggled per packet, signalled in header bit 11; no error propagation | - | HBCProtocol.py:109-144 | V |
| frame | FSC spreading | disabled (spreading_enabled = False, no caller enables it); if enabled SF = 2 x multiplier (16/8/4/2, 2-chip code), not the standard's 8-chip FSC (64..8); FSC_MAPPING unused | - | HBCProtocol.py:8-10, 79-82, 102-114, 160-169 | V |
| decode | Walsh decoder | hard ML correlation, ties to the lowest index; dmin 8; all 3-chip errors corrected; 4-chip patterns 73.1 % on average (1.000 for 0000 down to 0.481 for 1111); 5-chip patterns 32.7 % | - | HBCProtocol.py:448-461; exhaustive run | V |
| decode | sync / pilot detection | correlation >= 50 over 64 chips (<= 7 chip errors); needs the 4th preamble and the SFD 64 chips later; pilots removed by correlation search anywhere, so a missed pilot shifts the rest and a false match deletes 64 chips | - | HBCProtocol.py:241-268, 297-322, 392-408 | V |
| PHY v2 | BFSK | tones fc -/+ Rc (m = 1), receive band-pass 4 Rc (BFSK) or 2 Rc (OOK); non-coherent quadrature energy detector; noise σ = sqrt(N0 fs / 2) | - | v2.py:52-72, 120-140 | V |
| PHY v2 | energy reference | "Eb" in code and earlier work = energy per chip, a^2 / (2 Rc); per payload bit add 6.35 dB (L = 255) | - | v2.py:75-82 | V |
| PHY v2 | theory match | Pe = 0.5 exp(-Ec / 2N0); the simulator runs 0.29..0.40 dB better | dB | re-run; ch05 sec03:30 | V |
| PHY v2 | error-free packets | Ec/N0 >= 8 dB (1.5 % chip errors, 12/12); 6 dB: 5.5 %, 2/12; step about 2 dB wide (BFSK 8 MHz, 50 kchip/s) | dB | ch05 sec03:32 | V |
| MAC | TDMA class | slot 128 chips, guard 10, sync 8 (unseeded random word, 0 errors allowed); frame N (128 + 10) + 8 chips; not used in any earlier result (SINGLE_PACKET = True, MA = None) | chip | MAC_layer.py:47-49, 96-135 | V |
| MAC | FDMA option | tones linspace(2 MHz + 2 Rc, 8 MHz - 2 Rc, N), Rc 200 kchip/s, static; uses the v1 PHY (chip 1/(Rc x mod_index) = 2.5 µs, tones +/- 400 kHz) | - | MAC_layer.py:37, 154-161; PHY_layer_BCC_simulation.py:24-29 | V |
| MAC | CDMA option | Walsh codes of length 2^ceil(sqrt N) (too short for N = 9; should be 2^ceil(log2 N)); chip-synchronous sum | - | MAC_layer.py:667-689 | V |
| rates | payload rate | Rc / 4.31: 11.6 kb/s at 50 kchip/s, 46.4 at 200k, 116 at 500k (highest simulated), 464 at 2 Mchip/s | kb/s | ch05 sec04:53 | V |
| rates | airtime (L = 255) | 176 ms (50k), 44 ms (200k), 17.6 ms (500k), 4.4 ms (2M) | - | derived | V |
| rates | tone bandwidth, ceiling | 4 Rc per BFSK tone; aggregate B / 17.3 = 0.6..0.9 Mb/s per body (usable B 11.5 MHz at 50 Ω, 14.9 MHz at 1 MΩ; 0.6 Mb/s used in Ch. 6) | - | ch05 sec03:37-39; fdma_ceiling.py | V |
| link | noise | N0 = 4 k T 50 Ω = 8.0e-19 V^2/Hz x NF (10 / 20 / 30 dB) | V^2/Hz | link_budget.py:17-20 | V |
| link | drive | campaign 1 V amplitude (2 Vpp open circuit); margin over the weakest path 26-56 dB up to 500 kchip/s, >= 20 dB at 2 Mchip/s; weakest path gain -52.4 dB (RHW1-LHW1, 8 MHz, 50 Ω) | - | ch03 sec01:9; ch08b sec:13 | V |
| scenario | contact model | dwell 35000 / 13000 chips (good / degraded), depth 12 dB, Ec/N0 14 dB good; gives PSR 0.78 over 36 packets | - | ch05 sec03 (contact model; set to reproduce a capture share of 0.66-0.78 that 5.8 no longer treats as measured) | V |
| withdrawn | HackRF PSR | OOK 50 kchip/s on the arm 0.66-0.72; OOK 100 / 200 kchip/s 0.36 / 0.13; BFSK 100k 0.36; GFSK 8 MHz 0.23 (50k), 0.43 (200k) | - | ch05 sec08:43-66 | V |

Simulation recipe:
- FAIR = FULL = the stack as coded. Packet N(L) chips, airtime N / Rc; one BFSK tone per link; Ec/N0 = (A |H|)^2 / (2 Rc N0).
- Error model below; one Bernoulli draw per packet; no FEC, no FCS, so a failed packet is lost silently. Switch: +2-octet FCS (64 chips) for F7 parity.
- Access inside a body: hub TDMA, collision-free; slot N(L) + 10 chips per device per round, 8-chip sync per frame. FDMA option: tone spacing >= 4 Rc x 1.25 within the usable band.
- Across bodies: no coupling by default (galvanic isolation, section 5.7), so BCC PSR does not depend on k. Sensitivity cases: coupling loss = intra-body loss + 40 or + 60 dB, applied when same-tone slots overlap.
- Reporting cases: (a) noise model, PSR 1.0 at 1 V drive; (b) scenario (not measured): Bernoulli PSR 0.66-0.72 at 50 kchip/s OOK (11.6 kb/s), or the two-state contact chain.
- FULL-BCC ablation (optional): FSC spreading on, airtime x SF, Ec x SF.

BER model:
- Chip: p = 0.5 exp(-γc / 2), γc = Ec/N0 (add +0.35 dB to match the v2 simulator). OOK: p ~= 0.5 exp(-γc / 4) (textbook, not checked against the adaptive threshold, A).
- Walsh symbol: Ps_ok = Σ_{k=0..3} C(16,k) p^k (1-p)^(16-k) + 0.731 C(16,4) p^4 (1-p)^12 + 0.327 C(16,5) p^5 (1-p)^11.
- Sync and pilots: Q64 = Σ_{k=0..7} C(64,k) p^k (1-p)^(64-k).
- PSR = Q64^(2 + P) x Ps_ok^(8 + 2L), P = floor(32 L / 2048). Check: 0.011 / 0.54 / 0.952 / 0.998 / 1.000 at γc = 6 / 7 / 8 / 9 / 10 dB (L = 255, without the 5-error term); p = 5.5 % gives 0.168 with the 5-error term (measured 2/12 = 0.167), p = 1.5 % gives 0.99 (12/12).
- Contact impairment: γc drops 12 dB in the bad state of a chip-level two-state chain with P(g->b) = 1/35000, P(b->g) = 1/13000.
- Interference, only if coupling is enabled: SIR < 0 dB on a chip -> p = 0.5; otherwise p = 0.5 exp(-SINR/2) with interference as Gaussian noise (A).

## 4. Radio energy

### 4.1 Method

- Radio states per node: SLEEP, WAKE (HFXO and CPU), RAMP_TX, RAMP_RX, TX, RX, IDLE (CPU on, radio off). Charge per event Q = Σ I_state t_state; energy E = V Q with V = 3.0 V (datasheet currents are at 3.0 V with DC/DC). Average current over a period T: I_avg = Q_events / T + I_sleep (1 - t_active / T). Log E_TX, E_RX, E_ramp, E_idle, E_sleep per node (hub and each sensor).
- Use system-level ("compounded") currents: they include HFXO, regulator and DMA. For the nRF52840 the system overhead is additive, about 1.6 mA over the radio-only value (derived from pairs of table rows).
- Ramp: current = the start-up current of the target state (nRF52840 compounded RX 5.3 mA, TX 6.8 mA; corrected, the tables do give it); duration = chip ramp time. Protocol turnarounds (BLE T_IFS 150 µs, 802.15.4 192 µs) set the timing; the chip ramp inside them sets the energy.
- WAKE per radio event: 0.3 ms at 3 mA (A: covers HFXO start and CPU; same for every protocol and for BCC).
- Sleep current: 3.16 µA (nRF52840 System ON, full RAM, RTC); option 1.50 µA. The BCC node uses the same host MCU sleep current.
- Receive guard per event: BLE window widening of section 3.1 (two 50 ppm clocks: 2 x 100 ppm x time since anchor, plus 2 x 16 µs jitter, A); TDMA protocols the same rule with their clocks.
- The chip profiles below use different definitions ("compounded" Nordic, "system current with the MCU in EM1" SiLabs, "system bus off" TI, product-page figures for the nRF54L15): treat the cross-chip spread as about +/-15 %.

### 4.2 Current per chip and state (3.0 V, DC/DC unless stated)

| chip | state | value | unit | source | mark |
|---|---|---|---|---|---|
| nRF52840 | TX 1M BLE, 0 dBm (IRADIO_TX1) | 6.40 | mA | [6] 5.2.1.5 | V |
| nRF52840 | TX +8 / -40 dBm (TX0 / TX2) | 16.40 / 3.83 | mA | [6] 5.2.1.5 | V |
| nRF52840 | TX 0 dBm, LDO (TX3) | 10.80 | mA | [6] 5.2.1.5 | V |
| nRF52840 | TX 802.15.4 0 dBm (TX5) | 6.40 | mA | [6] 5.2.1.5 | V |
| nRF52840 | RX 1M (RX0) / LDO (RX1) | 6.26 / 10.10 | mA | [6] 5.2.1.5 | V |
| nRF52840 | RX 802.15.4 (RX2) | 6.53 | mA | [6] 5.2.1.5 | V |
| nRF52840 | RX 2M, compounded | 6.86 = 5.2 + 1.66 overhead | mA | derived | A (corrected from 7.1) |
| nRF52840 | radio-only TX by power | -40 dBm 2.3; -20 2.7; -8 3.3; -4 3.1; 0 4.8; +4 9.6; +8 14.8 (non-linear: use the table + 1.6 mA, not linear interpolation) | mA | [6] 6.20.15.2 | V |
| nRF52840 | radio-only RX 1M / 2M | 4.6 / 5.2 (2M LDO 11.1) | mA | [6] 6.20.15.3 | V |
| nRF52840 | start-up, radio only | ISTART,RX,1M 3.7 (LDO 6.7); ISTART,TX 5.2 at +4 dBm (LDO 11.0) | mA | [6] 6.20.15.2-3 | V |
| nRF52840 | ramp times | BLE tTXEN = tRXEN 140 (fast 40), tTXDIS 6 / 4 (1M / 2M); 802.15.4 130 (fast 40), tTXDIS 21, turnaround <= 40 | µs | [6] 6.20.15.8 | V |
| nRF52840 | sleep | System ON + RTC: 1.50 (no RAM), 3.16 (256 kB RAM); System OFF 0.40 | µA | [6] 5.2.1 | V |
| nRF52840 | Coded PHY and ANT currents | as 1M (TX 6.40, RX 6.26); BG22 shows +5 % at 125 kb/s | mA | none | A |
| nRF52832 | TX 0 / +4 dBm; RX 1M / 2M | 5.3 / 7.5; 5.4 / 5.8 | mA | [59] 23.15.2-3 | V |
| nRF52832 | start-up TX (+4 dBm) / RX; ramp | 4.0 / 3.5 mA; 140 µs (fast 40) | - | [59] | V |
| nRF54L15 | RX 1M; TX 0 dBm; sleep | 3.4; 4.8; 0.7..2.9 µA (product page, radio-only or system not stated) | mA | [66] | V |
| nRF54L15 | SmartBAN implementation, 1.8 V, 0 dBm | TX 5.83 measured (4.8 nominal), RX 12.43 measured (3.4 nominal), radio off 0.54 | mA | [51] Table VI | S |
| CC2652R | RX; TX 0 / +5 dBm; standby | 6.9; 7.0 / 9.2; 0.94 µA (RCOSC_LF), 1.09 µA (XOSC_LF) | mA | [22] 8.5, 8.6 (rev. J) | V |
| CC2652R | idle; RF core delta; wake standby->active | 675 µA; 210.9 µA; 160 µs | - | [22] 8.5, 8.14.2 | V |
| CC2340R5 | RX 1M; TX 0 / +8 dBm; standby | 5.3; 5.1 / 10.7; 0.71 µA | mA | [67] 7.9, 7.10 | V |
| EFR32BG22 | RX 1M / 2M / 125k (VSCALE1); TX 0 dBm; EM2 | 3.8 / 4.2 / 4.0; 4.3; 1.40 µA | mA | [68] Table 4.8 | V |
| EFR32MG24 | RX 1M / 2M / 802.15.4 / 125k (VSCALE1); TX 0 dBm; EM2 | 4.7 / 5.4 / 5.4 / 4.9; 5.0; 1.3 µA (2.9 µA full RAM) | mA | [69] Table 4.6 | V |
| EFR32MG24 | sensitivity 1M / 2M / 125k / O-QPSK | -97.6 / -94.8 / -105.7 / -105.4 | dBm | [69] p.1 | V |
| nRF24AP2 (ANT) | charge per message: master bidirectional / TX only / ack TX / slave RX / ack RX (byte-sync interface) | 27 / 17 / 35 / 21 / 26 | µC | [56] Table 12 | V |
| nRF24AP2 (ANT) | peak TX (0 dBm) / RX; search; base; burst at 20 kb/s | 15 / 17 mA; 2.8 mA; 3.0 µA; 5.9 mA | - | [56] Tables 11-12 | V |
| nRF52 S212 (ANT) | charge per master slot (ramps + TX 144 µs + RX window 250 µs + 0.5 µC overhead) | 3.3..3.7 | µC | estimate from [59] | A |
| nRF24L01+ (ESB) | TX 0 dBm; RX 1M / 2M; TX / RX settling (130 µs); standby-I | 11.3; 13.1 / 13.5; 8.0 / 8.9 mA; 26 µA | mA | [58] Table 5 | V |
| CC2570 (ANT) | RX; TX 0 / +4 dBm | 23.7; 28.8 / 34.3 | mA | [57] p.2 | V |
| 802.15.6 NB IC (Wong 2013, 1 V, 0.13 µm) | RX (802.15.6 and BLE); 802.15.6 TX at -10 dBm; BLE TX at -10 dBm | 4.8 mA; 5.9 mW; 4.6 mA | - | [70] | S |

Default chip per protocol: BLE, SmartBAN, 802.15.6 NB and ANT use the nRF52840 1M currents (the NB and SmartBAN PHYs run on BLE-class silicon; ANT alternatively uses the nRF24AP2 per-message charges as a legacy case); 802.15.4 CSMA and TSCH use the nRF52840 802.15.4 currents (CC2652R as the alternative profile); HBC uses 5.2 mW total (TX = RX = 2.6 mW, A split); BCC uses section 4.3.

Energy per on-air bit (nRF52840, compounded, 3 V, derived A): BLE 1M TX 19.2 / RX 18.8 nJ/b; 802.15.4 TX 76.8 / RX 78.4 nJ/b; Coded S=8 TX 154 nJ/b.

### 4.3 BCC and HBC transceiver literature, and the simulator picks

| work | coupling | figure | derived | source | mark |
|---|---|---|---|---|---|
| Zimmerman 1995, PAN prototype | capacitive | 1.5 mW TX electrode drive (330 kHz, 30 V, 10 pF; 0.5 C f V^2 = 1.49 mW); 2400 b/s | 625 nJ/b; RX front end 0.4 mW | [71] Sec. 1.5 | V |
| Song et al. 2007, 0.25 µm | capacitive (direct-coupled, carrier-less) | 0.2 mW at 2 Mb/s, BER 1.1e-7 | 0.1 nJ/b | [72] abstract via search | S |
| Cho et al. 2009, 0.18 µm | capacitive, 30-120 MHz AFH | 0.37 nJ/b at 10 Mb/s; body driver 0.7 mW; FSK modulator 1.9 mW; -65 dBm | 3.7 mW | [73] | V |
| Bae et al. Jan 2012 | capacitive, 40-120 MHz double-FSK | TX 2.0 mW, RX 2.4 mW at 10 Mb/s | 0.2 / 0.24 nJ/b | [78] abstract | V (corrected: not bae2012low) |
| Lee et al. 2013 | 802.15.6 HBC | 5.5 mW | 4.2 nJ/b at 1.3125 Mb/s | [74] title | S |
| Cho et al. 2015, 0.13 µm | 802.15.6 HBC | 5.2 mW (4.3 mA at 1.2 V) | 3.96 nJ/b | [75] | S |
| Maity et al. 2019, BodyWire, 65 nm | capacitive EQS | 6.3 pJ/b at 30 Mb/s; TX 3.1, RX 3.2 pJ/b | TX 93, RX 96 µW | [76] | V |
| Maity et al. 2021, Sub-µWRComm | capacitive EQS | 415 nW at 1-10 kb/s | 41.5-415 pJ/b | [80] | S |
| Chatterjee et al. ISSCC 2022 / RFIC 2020 | capacitive | 63.3 µW at 15 Mb/s; 2.24 pJ/b TX (NB-HBC mode) | 4.2 pJ/b | [80] | S |
| Jeon et al. 2019, 0.18 µm | galvanic (implant TX) | TX 4.75 pJ/b, RX 26.8 pJ/b at 100 Mb/s, BER 1e-9 | 0.475 / 2.68 mW | [77] | S |
| Wi-R (Ixana) | capacitive EQS | ~0.2 nJ/b at 5 Mb/s (vendor claim) | ~1 mW | [81] | V (vendor) |
| Wegmueller 2007 | galvanic, 256 kHz BPSK, 64 kb/s | 726 mW | - | (wegmuller2007intra) | S |
| Pitou 2023 | galvanic, discrete | <= 59.4 mW total (earlier table); a search summary gives 111 mW at 268.77 MHz; PDF bot-protected | - | [79]  | A |
| Shi 2022 | galvanic, 200 MHz, 250 Mb/s | 2 pJ/b | - | (Shi2022Galvanic-Coupled) | S |

BCC simulator picks (active power, held constant at the BCC rates because static bias dominates):

| case | P_TX | P_RX | basis | mark |
|---|---|---|---|---|
| low | 0.5 mW | 2.7 mW | Jeon 2019, the only galvanic IC (its TX is an implant) | A |
| typical | 2.0 mW | 2.4 mW | Bae Jan 2012 (capacitive); total consistent with the 5.2-5.5 mW 802.15.6 HBC ICs | A |
| high | 59.4 mW | 59.4 mW | Pitou 2023 discrete galvanic total; split assumed | A |
| wake / settling (low / typical / high) | 20 / 100 / 1000 µs | at the target state's power | no PLL; AFE settling (high = MCU + DDS start) | A |
| sleep | 3.16 µA | 3.16 µA | host MCU with RTC, same as the radio nodes (option 1.50 µA) | A |

At 50 kb/s the typical P_TX gives 40 nJ/b, twice the nRF52840 BLE 1M figure. No verified on-body galvanic IC power figure exists; the labels say which designs are capacitive.

### 4.4 Battery and lifetime

| parameter | value | unit | source | mark |
|---|---|---|---|---|
| CR2032 nominal voltage | 3.0 | V | [82] | V |
| capacity | 235 mAh to 2.0 V at 15 kΩ (0.19 mA), 21 °C | mAh | [82] | V |
| self-discharge | ~1 %/year | - | [82] | V |
| pulse test, temperature | 400 Ω (~6.8 mA at 2.7 V) 2 s x 12/day; -30..60 °C | - | [82] | V |
| peak-current derating | 15 -> 30 mA peak costs 9 % capacity on average; poor cells about 50 % of rated; bulk capacitor ~100 µF helps poor cells by > 40 %; IR limit ~30 Ω (30 mA) and ~60 Ω (15 mA) | - | [83] | V |
| usable fraction k | 0.9 typical, 0.8 conservative, 0.5 poor cells | - | derived from [83] | A |
| mean current reference | BLE 10.1 µA, ZigBee 15.7 µA, ANT 28.2 µA at a 120 s interval, 3.3 V; sleep 0.78 / 4.18 / 3.1 µA | µA | [65] Table II | V |

Lifetime: T_life = k C / (I_avg + s k C / 8760 h), C = 235 mAh, s = 0.01 per year. With a buck converter the input current rises by about 3.0 / V_cell late in life (5-7 % at 2.8 V); add this as a correction. Flag peaks above 15 mA (for example +8 dBm TX, 16.4 mA) and assume a 100 µF bulk capacitor there.

## 5. Channel

### 5.1 On-body links at 2.4 GHz

| model | value | source | mark |
|---|---|---|---|
| CM3A hospital room, 2.4-2.5 GHz | PL = 6.6 log10(d_mm) + 36.1 + N(0, 3.80 dB): 50.5 / 53.6 / 55.9 dB at 150 / 450 / 1000 mm | [84] 8.2.5.A | V |
| CM3A anechoic | PL = 29.3 log10(d_mm) - 16.8 + N(0, 6.89 dB): 47.0 / 60.9 / 71.1 dB | [84] 8.2.5.A | V |
| CM3B (IMEC), 2.45 GHz | PL = -10 log10(P0 exp(-m0 d) + P1) + 3.6 nP, P0 -25.8 dB, m0 2.0 dB/cm, P1 -71.3 dB (linear: P0 0.0026, m 46.09 1/m, P1 7.4683e-8): 45.8 / 55.7 / 64.7 / 71.1 dB at 10 / 15 / 20 / 30 cm, floor 71.3 dB | [84] 8.2.5.B; [85] Table 6 | V |
| CM3B fading | Ricean K_dB = 30.6 - 0.43 PL_dB + N(0, 3.4 dB); trms 6 ns (15 cm), 16 ns (45 cm) | [84] 8.2.5.B | V |
| CM3C NICTA, 2.36 GHz, standing / walking / running | Rx right hip from chest 65.3/59.1/55.9, right wrist 44.5/47.3/36.3, left wrist 74.7/59.8/52.5, right ankle 60.9/53.9/55, left ankle 70.7/58.5/59, back 75.3/67.4/68.5; Rx chest from back 73/72/71.7, right wrist 70.5/64.9/57.4, right ankle 66.3/62.4/63.3 dB | [84] 8.2.5.C | V |
| NICTA everyday gain, mean / median | hip-chest -50.1 / -62.5; chest-right wrist -66.3 / -70; chest-left ankle -78.8 / -82; hip-hip -53.2 / -83.5; overall -56.4 / -69.8 dB. Means are linear-power averages dominated by strong samples; the medians are 3-30 dB lower. Overall set lognormal (µ -7.74, σ 1.1, read as amplitude: ~9.6 dB spread); per-link temporal fading gamma, shape 1.44-2.9 (2.8-4.3 dB); fade duration ~22 ms median; level-crossing rate 0.27-1.46 Hz | [84] 8.2.5.C | V |
| CEA-LETI average gain | head-belt -48.2, head-wrist -53.9, belt-wrist -50.0 dB (Rice fading) | [85] Table 1 | S |
| Reusens, Latre survey | n about 3.3 along arm and torso; LOS 3-4, around the torso 5-6; P0 -32..-41 dB at 10 cm; antenna 5 mm vs 5 cm from the body > 20 dB; limb blocking up to 20-30 dB | [93]; [92] | S |

Recommended on-body set (hub on the chest), corrected to medians (the verifier found the means 12-30 dB optimistic for most links):

| case | median PL | static link shadowing σ | per-packet fade σ | use |
|---|---|---|---|---|
| low | 54 dB | 3.6 dB | 2.8 dB | short same-side links (head-belt) |
| typical | 68 dB | 4.5 dB | 3.5 dB | chest to wrist; overall median |
| high | 82 dB | 6.9 dB | 4.3 dB | chest to opposite ankle or back |

Per link (A): sim value = published median, else mean + 6 dB (mean-to-median gap 3-13 dB on the NICTA links): head-belt 54, belt-wrist 56, head-wrist 60, chest-wrist same side 64 / opposite 70, chest-ankle same side 75 / opposite 82, chest-back 78, hip-chest 62.5, hip-head 70 dB. A single generic link uses median 68 dB with σ_total 9.6 dB. Block fading per packet (coherence tens of ms); redraw the fade on a retry only if more than 20 ms have passed.

### 5.2 Body-to-body and off-body links at 2.4 GHz

| model | value | source | mark |
|---|---|---|---|
| CM4 off-body, chest TX, 2.36 GHz | LOS standing 53.81 / 53.12 / 56.04 / 64.72, LOS walking 44.46 / 51.59 / 52.14 / 60.81, NLOS standing 61.81 / 68.64 / 60.12 / 63.10, NLOS walking 59.48 / 70.30 / 63.98 / 62.56 dB at 1 / 2 / 3 / 4 m | [84] 8.2.8 | V |
| fit to CM4 (d0 1 m) | LOS standing PL0 51.6 dB, n 1.53, sd 4.4; LOS walking 44.0 dB, n 2.40, sd 2.9; NLOS flat ~63-64 dB; NLOS excess 6.5 (standing), 11.8 dB (walking) | least squares on 4 points | A |
| Rosini et al. body-to-body, 2.45 GHz | G = G0(α) - 10 n(α) log10(d / 1 m); head-head one room n 1.28, G0 -41.5 dB; two rooms n 0.94, G0 -46.9; belt-belt one room: n = -0.007α + 1.84 (0-90°), -0.027α + 3.64 (90-130°), 0.011α - 1.3 (130-180°); G0 = -0.19α - 29.5 (0-100°), -0.55α + 6.5 (100-120°), 0.06α - 66.7 (120-180°) dB; belt-belt two rooms n 1.7, G0 -0.068α - 41. Belt-belt n drops to 0.13 at 130°; G0 swings ~26-30 dB with orientation; G0 -29.5 dB at 1 m beats free space, so antenna gain or d0 is unclear | [85] Table 4 quoting [86] | S |
| ITU-R P.1238-8 | L = 20 log10 f_MHz - 28 + N log10 d (39.7 dB at 1 m, 2440 MHz); N 28 residential, 30 office; σ 8-10 dB (1.8-2 GHz); for base-to-portable links through walls, so an upper bound for open rooms | [87] eq. 1, Tables 2, 4 | V |
| body shadowing, 5.8 GHz | NLOS - LOS > 10 dB for most positions; LOS n < 2 indoors; deep fades up to 67 dB | [88] | V |
| people crossing the path | up to 20 dB near a terminal; furniture 5-13 dB | [89] | S |
| inter-BAN SIR | dominated by movement, varies up to 60 dB | [90] (Hanlen 2010) | S |
| body-to-body fading | κ-µ; 8.69 dB diversity gain with 4 antennas | [91] | S |
| body attenuation | body 19.2, head 13.0, clothing 1.7 dB | [90] (Alomainy) | S |

Recommended inter-BAN model (transmitter on body i, receiver on body j): PL = PL0 + 10 n log10(max(d, 0.3) / 1 m) + X + L_sh m_blk + per-packet fade.

| case | PL0 | n | X σ | L_sh per blocking torso | meaning |
|---|---|---|---|---|---|
| low loss (strong interference) | 42 dB | 1.3 | 3 dB | 6 dB | pessimistic for coexistence |
| typical | 48 dB | 2.0 | 5 dB | 15 dB | default |
| high loss | 55 dB | 3.0 | 8 dB | 30 dB | optimistic |

m_blk = number of torso cylinders the straight TX-RX segment crosses (including the own torsos), capped at 2; X drawn per node pair per mobility epoch; per-packet fade lognormal σ 5 dB (or Rayleigh). All sources start at 1 m; the 0.3 m clip is an assumption (A).

### 5.3 Room and bodies

| parameter | value | source | mark |
|---|---|---|---|
| US hospital, multi-bed rooms | >= 100 ft2 = 9.29 m2 per bed (single room 120 ft2 = 11.15 m2; the source prints 10.8 m2) | [94] 7.2.A2 | V |
| UK adult in-patient | bed space 3.6 x 3.7 m = 13.3 m2; 4-bed room 64.0 m2 (16.0 m2/bed), 72.5 m2 with en-suite (18.1 m2/bed); single room 23.5 m2 | [95] 3.1-4.24 | V |
| occupant density | office 18.6 m2/person; aerobics room 2.3; weight room 9.3 | [96] Table 403.3.1.1 | V |
| sports hall | ~43 m2/person (4-court badminton hall 690 m2 with 16 players) | court sizes | A |
| staff and visitors | 1.2-1.5 persons per bed in wards | - | A |

Room generator: area A = k a_p, square room. a_p (m2/person) low / typical / high density: ward 18.1 / 16.0 / 9.3; gym 43 / 9.3 / 2.3; office 18.6 / 12 / 8 (office typical and high A). Body centres uniform with >= 0.6 m spacing; facing angle uniform; torso = vertical cylinder of radius 0.15 m. Node positions (body frame x forward, y left, z up, m): hub chest (0.15, 0, 1.3); head (0.10, 0, 1.65); belt (0.15, 0, 1.0); wrists (0.05, +/-0.30, 0.9); ankles (0.05, +/-0.12, 0.1); back (-0.15, 0, 1.3). Mobility epoch 1 s, static by default. All A.

### 5.4 Noise and receiver calibration

- Thermal noise -174 dBm/Hz (kT0 = -173.98 at 290 K). N = -174 + 10 log10(B) + NF + L_impl.
- Published noise figures: AT86RF233 6 dB [24]; the 802.15.6 NB draft assumes NF 10 dB + 6 dB implementation loss [27]; SmartBAN reports assume NF 13 dB + IL 6 dB [44] or NF 10 + IL 5 [50]. Other chips give sensitivity only, so NF + L_impl is fitted.
- Fitted NF + L_impl (A): BLE 1M/2M nRF52840 8.1 dB with curve (a) or 4.7 dB with curve (b); BLE 1M nRF52832 (-96 dBm) 7.1 dB; 802.15.4 nRF52840 10.6, AT86RF233 9.6, CC2420 15.6, spec 25.6 dB; ANT nRF24AP2 18.1, nRF52832 10.1 dB; NB NF 10 + IL 6; SmartBAN as BLE (common curve) or NF 13 + IL 6 (ETSI mode). Recommended generic NF 4 / 6 / 10 dB with L_impl 0 / 2 / 5 dB (low / typical / high).
- nRF52832 sensitivity: -96 dBm BLE 1M, -89 dBm BLE 2M, -93 dBm ANT 1M [59] (product brief); the full PS may say -93 dBm for 2M (not read).

### 5.5 Co-channel and adjacent-channel rejection used in SINR mode

| protocol | C/I_co | ACR by offset | source | mark |
|---|---|---|---|---|
| BLE 1M | 21 (spec), 6 (nRF52840) | 6 / 38 / 48 dB at 1 / 2 / >= 3 MHz (spec) | [1] Table 4.3; [6] | V |
| BLE 2M | 21 (spec), 6.8 | 6 / 38 / 48 dB at 2 / 4 / >= 6 MHz | [1] Table 4.4 | V |
| BLE Coded S=2 / S=8 | 17 / 12 (spec), 4.4 (nRF52840 125k) | 6 / 38 / 48 dB | [1] Tables 4.5-4.6 | V |
| 802.15.4 | 3 (CC2420) | chip 33 / 48 / 54 dB at 5 / 10 / 15 MHz; spec 0 / 30 dB | [23]; [24]; [13] | V |
| ANT | 12 (nRF24AP2) | 4 / 33 / 42 / 52 / 62 dB at 1 / 2 / 3 / > 6 / > 25 MHz (from C/I 8, -21, -30, -40, -50) | [56] | V |
| SmartBAN | none given; BLE spec 21 as analogue | BLE 1M spec | [47] 4.8 says a mask should have been given | A |
| 802.15.6 NB | none found | channel offset 1: 0 dB worst, 20 dB case; offset >= 2 ignored | - | A |
| BLE vs 802.15.4 (cross-technology) | - | 0 dB if carriers <= 1 MHz apart, 10 dB at 2 MHz, else the >= 3 MHz value | - | A |

Treating interference as Gaussian noise gives too permissive a capture threshold (BLE: 10.9 dB against the spec's 21 dB; 802.15.4: 0.4 dB against about 3 dB); SINR mode therefore also requires SIR >= C/I_co over the packet.

### 5.6 Interference bookkeeping

Keep every transmission as (t_start, t_end, carrier, P_tx, node). For each reception, split the packet at every interferer edge into chunks; SINR_j = S / (N + Σ_i 10^((P_rx,i - ACR(|Δf_i|)) / 10)). Packet success = Π_j (1 - BER_j)^(n_j); a packet below sensitivity is lost; the receiver locks on the first packet above sensitivity (a later stronger packet is only interference; capture-restart flag). The flag destroy_on_overlap reproduces `BLE_SIM.py`. Sanity check (typical sets): on-body 68 dB, 0 dBm, S = -68 dBm, 27 dB above BLE 1M sensitivity; a co-channel interferer 1 m away in LOS (48 dB) gives SIR -20 dB (lost); 2 MHz away on 1M (ACR 38 dB) +18 dB (survives, above the 10.9 dB curve threshold, below the 21 dB spec C/I); on 2M at 2 MHz (ACR 6 dB) -14 dB (lost).

### 5.7 Inter-body coupling of BCC and HBC

| item | value | source | mark |
|---|---|---|---|
| capacitive EQS HBC | at 1 m, intra-body loss 60 dB vs inter-body 70 dB (anechoic) to ~80 dB (open air): inter-body only 10-20 dB weaker; Coupling(d) = C_C(d) / C_body, C_body ~150 pF; TX kept 6-9 dB above sensitivity (link margin, not an SNR requirement) | [97] | V (corrected wording) |
| galvanic vs capacitive leakage (phantom arm, 100-400 kHz) | galvanic: RSS at an off-body adversary at the noise floor at 13 cm, never decodable; capacitive: decodable to 15 cm (~1e-4 BER at 1 cm); > 7x less radiation; BER 1e-6 at -2 dBm | [98] | V |
| capacitive HBC interference range | another user's HBC device disturbed at ~150 cm | [99] via search | S |
| HBC co-channel interference at 21 MHz | modelled between HBC users on the single channel | [37] | S |

Simulator defaults (A):
- BCC galvanic: isolated between bodies (no inter-BAN term), an assumption, supported qualitatively by [98] (no measured galvanic body-to-body number exists). Sensitivity cases: inter-body loss = intra-body loss + 40 dB or + 60 dB for bodies not in contact; bodies in contact (handshake) share one channel.
- HBC standard (capacitive, 21 MHz): intra-body loss 60 dB (typical); inter-body loss = intra + 10 / 15 / 20 dB at 1 m (low / typical / high), scaled by 20 log10(d / 1 m) (A, 1/d coupling capacitance). This makes HBC PSR depend on k.

## 6. Open issues

Union of the research agents' open issues and the verifiers' "missing" items. Each line ends with the default the simulator takes.

General
- OI-1 Primary text not read: IEEE 802.15.4-2015/2020 (TSCH clauses, PIB), IEEE 802.15.6-2012 (NB, HBC, MAC) and the 2026 revision. Default: 802.15.4-2006 values and the secondary sources as marked; the texts are available through IEEE GET or a library.
- OI-2 No primary source gives a BER curve for a GFSK discriminator receiver (BLE, SmartBAN, ANT). Default: common curve 0.5 exp(-γ/2) with NF + L_impl fitted per chip; the correlated Marcum-Q curve as the conservative option; the ETSI coherent curve only to reproduce ETSI tables.
- OI-3 Treating interference as Gaussian noise is too permissive for capture; co-channel C/I is missing for SmartBAN, NB, HBC and for 802.15.4 on the nRF52840. Default: collision mode is the main result; SINR mode adds the C/I_co floor of section 5.5 (BLE spec values as analogue where missing).
- OI-4 Application deadline D and the traffic model (payload and sample period per sensor) are not fixed by any source; FULL loss metrics depend on D. Default: D = one sample period; payload and period from the BLE_SIM throughput axis (R = L / CI); PDR and on-time ratio reported separately.
- OI-5 Hub scheduling when events of its own links collide, and the number of exchanges per event. Default: the earlier anchor wins, the later event is skipped (its counter advances); MD-driven continuation in BLE FULL; one exchange per event in FAIR.
- OI-6 Preamble and sync detection under interference (access address or SFD miss, lock, capture restart). Default: sync bits are part of the BER product; the receiver locks on the first packet above sensitivity; capture-restart flag off.
- OI-7 Cross-technology rejection between BLE and 802.15.4 has no source. Default: 0 dB at <= 1 MHz, 10 dB at 2 MHz, else the >= 3 MHz value; used only in mixed runs.
- OI-8 Clock drift and timing jitter. Default: per-device drift uniform (BLE and ANT +/-50 ppm, SmartBAN +/-20 ppm, 802.15.4 +/-40 ppm), no jitter; SmartBAN measured jitter (P99 < 754 µs) as a sensitivity case.

BLE 5
- OI-9 The FAIR BLE is not a standard mode: ACL cannot switch ARQ off; CIS and BIS single-shot modes hop with CSA#2 and need an ACL or a periodic advertising train. Default: option B (BIS-like) with CSA#1 per link, side traffic omitted and declared in earlier work.
- OI-10 Supervision timeout, reconnection delay, peripheral latency and subrating. Default: TO = max(2 s, 3 CI) rounded up to 10 ms; reconnect after 100 ms; latency 0; subrate 1.
- OI-11 AFH classification is vendor specific; the instant rule is a "should". Default: Spoerk rule with C_min 10, OWL speed-up off, instant = first transmission + 6.
- OI-12 Coded S=8 gain: model 10.7 dB, nRF52840 8 dB, spec 12 dB; the distance spectrum is computed, not read. Default: soft model with 2.7 dB extra loss at S=8.
- OI-13 Maximum simultaneous connections per hub is controller specific. Default: no cap; results with N_dev above 20 carry a note.
- OI-14 L2CAP + ATT overhead (7 octets) unverified; Core 4.0 minimum used channels assumed equal to 4.2; newer nRF52840 PS versions not checked. Default: overhead off (L is the LL payload).
- OI-15 `BLE_SIM.py` has no PHY/LL overhead and one permutation per body. Default: aligns it (section 3.1 FAIR step 6) before FULL/FAIR comparisons.

IEEE 802.15.4 and TSCH
- OI-16 TSCH template, macMaxBE (7 vs Contiki-NG 5), RxOffset (1120 vs 1020 µs), sequence origin, EB and Enh-Ack sizes, CCA-busy handling. Default: values as recorded, macMaxBE 7 as a parameter, RxOffset 1120, CCA-busy counts as a failed attempt.
- OI-17 FULL channel rules (Zigbee-style ED scan, TSCH blacklist) are not 802.15.4 behaviour. Default: ED scan on in modes A and B FULL, labelled upper layer; TSCH blacklist off (FULL+ flag).
- OI-18 Mode A buffer policy while a frame is in backoff or retry; duplicate handling. Default: FIFO depth 4, drop newest; duplicates discarded by DSN.
- OI-19 Beacon mode and TSCH quantise the period (15.36 ms x 2^BO; 10 ms). Default: round P to the grid and run every protocol at the rounded P.
- OI-20 Beacon tracking window, TSCH keep-alive (Contiki-NG 12 s) and desync, guard growth, join time. Default: steady state, synchronised from t = 0; add only for energy sensitivity.
- OI-21 Downlink traffic (indirect transmission, RX cells). Default: uplink only, for every protocol.
- OI-22 CCA threshold used by real chips, and whether another BAN's frame is sensed. Default: -75 dBm applied to the modelled received power.

IEEE 802.15.6 NB and MAC
- OI-23 Final NB sensitivities, PER criterion, SRRC roll-off, 1 MHz adjacent-channel rejection, channel count (79 vs 16 in ns-3). Default: BER model with NF 10 + IL 6; roll-off 0.3; ACR 0 dB (20 dB case); 79 channels.
- OI-24 No retry limit in the standard; I-Ack rate and timeout, guard times, clock accuracy (Castalia 100 ppm). Default: R = 3 (sweep 0/2/3/7); I-Ack at the data rate; timeout pSIFS + mTimeout; no extra guard; drift off. The capacity comparison charges the guard GTn of 6.11.1.1 to both stacks in closed form (claim4.py part G5, REPORT 5.13).
- OI-25 Beacon body length, hopping sequence and period, beacon shifting and interleaving. Default: 17-octet body; hop option = per-BAN random permutation of 79 channels, one hop per superframe; shifting and interleaving flag off.
- OI-26 Draft-era or single-source timing (pMIFS, pExtraIFS, pChannelSwitchTime, mClockResolution, mTimeout, allocation slot 500 µs). Default: values as recorded (A).

IEEE 802.15.6 HBC
- OI-27 Header rate (164 kb/s or RI mode), CRC-8 polynomial, pilot interval and field width, scrambler taps, absence of FEC, sensitivity. Default: header at 164 kb/s with an RI-mode flag; the polynomial does not enter the error model; pilot every 128 symbols; no FEC.
- OI-28 HBC MAC timing and contention (pSIFS, allocation slot, Aloha slot). Default: NB values; scheduled access; slotted Aloha option with 500 µs slots.
- OI-29 HBC link budget and inter-body coupling at 21 MHz. Default: intra-body 60 dB, inter-body +10 / +15 / +20 dB at 1 m, 20 dB per decade; γ on the same Ec/N0 axis as BCC.
- OI-30 The code's FSC (2-chip, SF 2..16) differs from the standard (8-chip, SF 8..64). Default: FULL HBC uses the standard SF in the chip-time accounting; the code path is not reused.

SmartBAN
- OI-31 BEmin, BEmax, RTmax, RTmax,PCA, retry limit and path, TC, DCH choice and migration rule, scan window: none normative. Default: 3 / 5 / 3 / 3; RT = 1 (sweep 3 and "until ACK"); TC = TD; rules of section 3.5.
- OI-32 Informative sensitivity counts repetition twice; TR 103 711's -92.5 dBm is a noise floor; LSLOT 1 accounting; Packet Length unit; coding of beacons and ACK; BCH calibration 1.1 or 1.6 dB; header BCH "shall" (V1.1.1) or "should" (V1.2.1). Default: full V1.2.1 accounting (LSLOT 1 carries no data); beacons and ACK uncoded, Nrep 1 (flag: coded as in [51]); 1.6 dB in ETSI mode; header BCH on in FULL.
- OI-33 Beacon-loss and connection timeouts, C-Beacon position, CCH/DCH switch time. Default: after 10 missed D-Beacons the node rescans [51]; C-Beacon at the start of the inactive period; switch time 0.

ANT and ESB
- OI-34 On-air frame (reverse engineered), ACK length, turnaround, master RX window, slave window widening. Default: 144 µs frames, 130 µs turnaround, W = 176 µs, widening 2 |ε| x elapsed time.
- OI-35 Adaptive isochronous algorithm has no public numbers and is absent from the protocol documents; frequency-agility trigger unverified. Default: both flags off; when on, G = 500 µs after the slot, S = 250 µs per period, agility after 8 unacked messages.
- OI-36 nRF52 (S212) charge per message unpublished; burst spacing on nRF52; advanced-burst frame format; nRF24AP2 burst rate depends on the host interface; asynchronous mode; which channel EVENT_CHANNEL_COLLISION blocks. Default: estimates as recorded, 3.2 ms spacing, the later channel is blocked, asynchronous mode not modelled.
- OI-37 ESB on nRF5: payload limit and NCS defaults secondary; ARD reference point. Default: as recorded; ARD counted from the end of the packet.

Energy
- OI-38 WAKE charge per radio event (HFXO start, CPU) assumed. Default: 0.3 ms at 3 mA; replace with Online Power Profiler numbers if obtained.
- OI-39 Not obtained: Coded PHY RX current, nRF52840 2M TX current, nRF54L15 2M / 802.15.4 / ramp, CC2652R7, nRF5340. Default: 1M values; CC2652R7 = CC2652R.
- OI-40 BCC energy basis: no verified on-body galvanic IC power; Pitou 2023 unverifiable (59.4 vs 111 mW); Jeon 2019 TX is an implant; Bae Jan 2012 is not in `biblio.bib` and gives its numbers to `bae2012low`; the Cho2007The row of Table 2.1 may be mis-attributed; the 1.5 mW of Zimmerman should cite `bcc_first_thesis` (the IBM article was not readable); Song 2007 JSSC is not in the bib. Default: the picks of section 4.3; thesis citation fixes go to a separate ticket.
- OI-41 BCC receiver sensitivity and noise figure, TX drive vs link budget, BCC hub listening model. Default: scheduled hub RX (TDMA slot only); link from earlier work Ec/N0 model.
- OI-42 Battery capacity at the real pulsed load and body temperature; nRF52840 1.7 V minimum vs the 2.0 V cut-off. Default: k = 0.9, cut-off 2.0 V, buck correction 3.0 / V_cell.
- OI-43 BCC rate ambiguity: `MAC_layer.py` runs 200 kchip/s on the v1 PHY (2.5 µs chips); the rate labels in `HBCProtocol.py` are inert. Default: Rc is the explicit sweep variable (50k, 200k, 500k, 2M chip/s).

Channel
- OI-44 On-body medians are published for a few links only. Default: the section 5.1 per-link values (median, else mean + 6 dB).
- OI-45 Primary body-to-body sources not read (Rosini 2014, Reusens 2009 Table II, Hanlen 2010, Wang and Cai 2011, a Cotton-Scanlon path-loss exponent); the CM4 fit uses 4 points; no source below 1 m; G0 antenna-gain convention. Default: the typical set of section 5.2 with a 0.3 m clip.
- OI-46 Sports-hall and office densities assumed; ASHRAE read through the NYC code. Default: as recorded.
- OI-47 Body-shadowing dynamics under mobility. Default: static bodies; 1 s mobility epochs as an option.
- OI-48 Galvanic body-to-body coupling never measured; handshake coupling not quantified. Default: isolated; +40 / +60 dB sensitivity cases.

BCC code
- OI-49 Code issues found in the audit: header CRC-8 fed bit per byte (dmin 2) with the inverted-header fallback; earlier work says the CRC rejects a corrupted header; pilots searched anywhere; TDMA sync word unseeded with zero tolerance; CDMA code length 2^ceil(sqrt N); the FDMA main references undefined names; earlier text says 8800 chips while the chip file carries 8864; thesis calls CDMA central although no MAC result uses it. Default: the simulator uses the PSR formula of section 3.8 (any header or payload symbol error loses the packet); code fixes go to a separate ticket.
- OI-50 Fair-rule asymmetry: BCC has no FCS and no MAC header (9 octets) while NB and HBC keep them. Default: keep as coded; the FCS switch for parity runs.
- OI-51 Missing `biblio.bib` entries: Core 5.4/6.0/6.2, IEEE 802.15.4-2006/2020, RFC 8180, draft-ietf-6tisch-minimal, Contiki-NG, the SmartBAN documents (TS 103 326, TS 103 325, TR 103 395, TR 103 711, TS 103 806), Viittala 2017, Hamalainen 2020, Takabayashi 2019, Yoshimi 2026, Politis 2026, Ullah 2013, Hwang 2015, Kang 2016, Khan 2017, Siddik 2023, the ANT documents (D00000652, D00000693, AN11), nRF24AP2, nRF24L01+, nRF52832, nRF52840, CC2652R, CC2420, AT86RF233, US8774072B2, Dementyev 2013, Bae Jan 2012, Song 2007 JSSC, Cho 2015, BodyWire, IEEE P802.15-08-0780, Alam 2016, Rosini 2014, ITU-R P.1238, Energizer CR2032, TI SWRA349. Default: add them before earlier work cites the values (separate ticket).

## 7. References

"bib:" gives the `biblio.bib` key when earlier work already has an entry for the source (or for the closest edition, as stated). All others are not yet in `biblio.bib` (OI-51).

1. Bluetooth SIG, Core Specification v5.4, Vol 6 Part A (Physical Layer) and Part B (Link Layer), https://www.bluetooth.com/wp-content/uploads/Files/Specification/HTML/Core-54/out/en/low-energy-controller/link-layer-specification.html and .../physical-layer-specification.html. bib: bluetooth_5_0 (v5.0 edition; v5.4 not in bib).
2. Bluetooth SIG, Core Specification v5.4, Vol 4 Part E (HCI) and Vol 2 Part B (BR/EDR Baseband), https://www.bluetooth.com/wp-content/uploads/Files/Specification/HTML/Core-54/out/en/host-controller-interface/host-controller-interface-functional-specification.html.
3. Bluetooth SIG, Core Specification v6.0, Vol 6 Part B, https://www.bluetooth.com/wp-content/uploads/Files/Specification/HTML/Core-60/out/en/low-energy-controller/link-layer-specification.html.
4. Bluetooth SIG, Core Specification v6.2, Vol 6 Part B, https://www.bluetooth.com/wp-content/uploads/Files/Specification/HTML/Core-62/out/en/low-energy-controller/link-layer-specification.html.
5. Bluetooth SIG, Core Specification v4.2, Vol 6 Part B (PDF mirror https://community.infineon.com/gfawx74859/attachments/gfawx74859/psoc4/34073/1/Core_v4.2.pdf). bib: bluetooth_4_0 (v4.0 edition).
6. Nordic Semiconductor, nRF52840 Product Specification v1.5 (4413_417 v1.5, 2021-09-06), https://files.seeedstudio.com/wiki/XIAO-BLE/nRF52840_PS_v1.5.pdf; current tables also at https://docs.nordicsemi.com/r/bundle/ps_nrf52840/.
7. M. Spoerk, J. Classen, C. A. Boano, M. Hollick, K. Roemer, "Improving the Reliability of Bluetooth Low Energy Connections", EWSN 2020, pp. 144-155, http://www.carloalbertoboano.com/documents/spoerk20reliability.pdf.
8. Zephyr Project, subsys/bluetooth/controller/ll_sw/lll_chan.c (CSA#2 unit tests), https://raw.githubusercontent.com/zephyrproject-rtos/zephyr/main/subsys/bluetooth/controller/ll_sw/lll_chan.c.
9. J. Song, S. Shetty, D. Gopalpet, "Coexistence of IEEE 802.11b and Bluetooth: An Integrated Performance Analysis", Mobile Netw. Appl. 12:450-459, 2008, doi:10.1007/s11036-008-0047-3.
10. J. G. Proakis, M. Salehi, Digital Communications, 5th ed., McGraw-Hill, 2008 (non-coherent FSK, M-ary orthogonal and multichannel DPSK formulas; not re-read).
11. University of Victoria ECE450 Lab 3 theory (BFSK tone spacing), https://mistic-lab.github.io/ece-communications-labs/ece450/lab3/theory.
12. A. Andras, A. Baricz, Marcum Q-function bounds and DQPSK error probability, arXiv:1203.2760.
13. IEEE Std 802.15.4-2006, full text (mirror https://people.ece.ubc.ca/edc/7860/data/802.15.4-2006.pdf).
14. IEEE Std 802.15.4-2020, doi:10.1109/IEEESTD.2020.9144691 (not read).
15. IEEE 802.15-05-0087-00-004a, 802.15.4b PHY clause draft, https://grouper.ieee.org/groups/802/15/pub/05/15-05-0087-00-004a-4b-phy-clause-draft.pdf.
16. IETF RFC 8180 (6TiSCH minimal configuration), https://www.rfc-editor.org/rfc/rfc8180.txt; RFC 7554 (TSCH in the IoT, informational).
17. draft-ietf-6tisch-minimal-05 (Jan 2015), https://www.ietf.org/archive/id/draft-ietf-6tisch-minimal-05.txt.
18. Contiki-NG source, os/net/mac/tsch: tsch-conf.h, tsch-timeslot-timing.c, tsch-const.h, tsch-queue.c (github.com/contiki-ng/contiki-ng).
19. Chen et al., "Performance Analysis of IEEE 802.15.4e TSCH for Low-Rate Wireless Networks", KSII TIIS 7(1), 2013, https://koreascience.kr/article/JAKO201308138130612.pdf.
20. IETF 6tisch mailing-list archive, verbatim quote of IEEE 802.15.4-2015 6.2.5.3 (message bzRq1rvKsHVlyK8353tntDeA1rI).
21. ns-3, src/lr-wpan (lr-wpan-error-model.cc, lr-wpan-phy.cc), https://github.com/nsnam/ns-3-dev-git/.
22. Texas Instruments, CC2652R data sheet SWRS207J (Nov 2023), https://www.ti.com/lit/ds/symlink/cc2652r.pdf.
23. Texas Instruments, CC2420 data sheet SWRS041c, https://www.ti.com/lit/ds/symlink/cc2420.pdf.
24. Microchip (Atmel), AT86RF233 data sheet 8351E, https://ww1.microchip.com/downloads/en/DeviceDoc/Atmel-8351-MCU_Wireless-AT86RF233_Datasheet.pdf.
25. Zephyr Project, include/zephyr/net/ieee802154_radio.h (802.15.4-2020 constants).
26. IEEE Std 802.15.6-2012 (not read); IEEE 802.15.6-2026 status page https://standards.ieee.org/ieee/802.15.6/10941/. bib: ieee_802_15_6.
27. Texas Instruments, US Patent 9,510,138 B2 (NB PHY of the TG6 proposal), https://patents.google.com/patent/US9510138B2/en.
28. S. Ullah, M. Mohaisen, M. A. Alnuem, "A Review of IEEE 802.15.6 MAC, PHY, and Security Specifications", IJDSN 2013, 950704, doi:10.1155/2013/950704.
29. K. S. Kwak, S. Ullah, N. Ullah, "An Overview of IEEE 802.15.6 Standard", arXiv:1102.4106 (Draft 1). bib: kwak2010overview.
30. Nasr et al., IEEE 802.15.6 NB PHY study, arXiv:1510.00304.
31. ns-3 IEEE 802.15.6 implementation study, IJERPH 17(11):4007, 2020, PMC7312563.
32. Siddik et al., "Effects Investigation of MAC and PHY Layer Parameters on the Performance of IEEE 802.15.6 CSMA/CA", IJWMN 15(6), 2023, arXiv:2401.07725.
33. Khan et al., "Performance Evaluation of RTS/CTS Scheme in Beacon-Enabled IEEE 802.15.6 MAC", PMC7219348.
34. IEEE 802.15.6 prototype with measured retry counts, arXiv:1701.02421.
35. Castalia simulator, BaselineBANMac.h/.cc/.ned, https://github.com/boulis/Castalia.
36. Samsung and ETRI, "EFC Proposal for HBC PHY", IEEE 802.15-10-0049-01-0006, Jan 2010.
37. Hwang, Kang, Kim, Park, "Analysis on Co-channel Interference of Human Body Communication Supporting IEEE 802.15.6 BAN Standard", ETRI J. 37(3):439-449, 2015, doi:10.4218/etrij.15.0114.1029.
38. Kang et al., "Highly Simplified and Bandwidth-Efficient Human Body Communications Based on IEEE 802.15.6 WBAN Standard", ETRI J. 38(6), 2016, doi:10.4218/etrij.16.2716.0003.
39. Khan et al., "Performance Analysis of Different Backoff Algorithms for WBAN-Based Emerging Sensor Networks", Sensors 17(3):492, 2017, doi:10.3390/s17030492.
40. Nataraju et al., IJEECS 35(1):601-609, 2024, doi:10.11591/ijeecs.v35.i1.pp601-609.
41. IBC IEEE 802.15.6 HBC prototype, Sensors 2020, PMC7764060.
42. ETSI TS 103 326 V1.2.1 (2021-07), SmartBAN enhanced ultra-low-power PHY (and V1.1.1), https://www.etsi.org/deliver/etsi_ts/103300_103399/103326/01.02.01_60/ts_103326v010201p.pdf.
43. ETSI TS 103 325 V1.2.1 (2022-07), SmartBAN low-complexity MAC (and V1.1.1), https://www.etsi.org/deliver/etsi_ts/103300_103399/103325/01.02.01_60/ts_103325v010201p.pdf.
44. ETSI TR 103 395 V1.1.2 (2021-06), SmartBAN system performance, https://www.etsi.org/deliver/etsi_tr/103300_103399/103395/01.01.02_60/tr_103395v010102p.pdf.
45. ETSI TR 103 711 V1.1.1 (2020-10), https://www.etsi.org/deliver/etsi_tr/103700_103799/103711/01.01.01_60/tr_103711v010101p.pdf.
46. ETSI TS 103 806 V1.1.1 (2023-11), SmartBAN hub-to-hub communication.
47. ETSI TR 103 944, SmartBAN.
48. Viittala et al., IEEE Access 5:8175-8182, 2017, doi:10.1109/ACCESS.2017.2697502.
49. Hamalainen et al., "ETSI SmartBAN Architecture: The Global Vision for Smart Body Area Networks", IEEE Access 8:150611-150625, 2020, doi:10.1109/ACCESS.2020.3016705.
50. Takabayashi, Tanaka, Sakakibara, Sensors 19(1):30, 2019, doi:10.3390/s19010030.
51. Yoshimi, Ito, Tanaka, Tanaka, "SmartBAN on Silicon by Structured Behavioral Modeling", arXiv:2608.29036, 2026.
52. Politis, Hilas, "Throughput and Delay Performance of Slotted Aloha in SmartBANs under Saturation Conditions", arXiv:2605.22317, 2026.
53. Dynastream, ANT Message Protocol and Usage, D00000652 Rev 5.1 (2014), mirror https://devzone.nordicsemi.com/cfs-file/__key/communityserver-discussions-components-files/4/D00000652_5F00_ANT_5F00_Message_5F00_Protocol_5F00_and_5F00_Usage_5F00_Rev_5F00_5.1.pdf.
54. ANT+ Device Profile Heart Rate Monitor, D00000693 Rev 1.13, https://err.no/tmp/ANT_Device_Profile_Heart_Rate_Monitor.pdf.
55. ANT AN11, ANT Channel Search and Background Scanning Channel, Rev 2.1.
56. Nordic Semiconductor, nRF24AP2 Product Specification v1.2.
57. Texas Instruments, CC2570/CC2571 data sheet SWRS095A, https://www.ti.com/lit/ds/symlink/cc2570.pdf.
58. Nordic Semiconductor, nRF24L01+ Product Specification v1.0, https://cdn.sparkfun.com/assets/3/d/8/5/1/nRF24L01P_Product_Specification_1_0.pdf.
59. Nordic Semiconductor, nRF52832 Product Specification v1.4 (mirror https://e2e.ti.com/cfs-file/__key/communityserver-discussions-components-files/6/nRF52832_5F00_PS_5F00_v1.4.pdf) and nRF52832 Product Brief v2.2.
60. sghctoma, antfs-poc-defcon24 sniffer, sniff/packets.h and shockburst.c, https://github.com/sghctoma/antfs-poc-defcon24.
61. Garmin (Morris, Rooney, Fung), US Patent 8,774,072 B2, https://patents.google.com/patent/US8774072B2/en.
62. Nordic nRF Connect SDK, include/esb.h and subsys/esb/Kconfig, https://github.com/nrfconnect/sdk-nrf.
63. Wikipedia, "ANT (network)", https://en.wikipedia.org/wiki/ANT_(network).
64. Nordic Semiconductor, S212 ANT SoftDevice product page, https://www.nordicsemi.com/Products/Development-software/S212-ANT.
65. A. Dementyev, S. Hodges, S. Taylor, J. Smith, "Power consumption analysis of Bluetooth Low Energy, ZigBee and ANT sensor nodes in a cyclic sleep scenario", IEEE IWS 2013, doi:10.1109/IEEE-IWS.2013.6616827.
66. Nordic Semiconductor, nRF54L15 product page, https://www.nordicsemi.com/Products/nRF54L15.
67. Texas Instruments, CC2340R2/R5 data sheet SWRS272F, https://www.ti.com/lit/ds/symlink/cc2340r5.pdf.
68. Silicon Labs, EFR32BG22 data sheet Rev 1.3, https://www.silabs.com/documents/public/data-sheets/efr32bg22-datasheet.pdf.
69. Silicon Labs, EFR32MG24 data sheet Rev 1.2, https://www.silabs.com/documents/public/data-sheets/efr32mg24-datasheet.pdf.
70. Wong et al., "A 1 V 5 mA Multimode IEEE 802.15.6/Bluetooth Low-Energy WBAN Transceiver", IEEE JSSC 48(1):186, 2013, via Microwaves & RF summary https://mwrf.com/technologies/semiconductors/article/21844956/1v-transceiver-exceeds-needs-of-ieee-802156-narrowbandble-phy.
71. T. G. Zimmerman, Personal Area Networks (PAN): Near-Field Intra-Body Communication, MS thesis, MIT, 1995, https://cba.mit.edu/docs/theses/95.09.zimmerman.pdf. bib: bcc_first_thesis (the IBM Syst. J. 1996 article bcc_first_article was not read).
72. S.-J. Song, N. Cho, H.-J. Yoo, "A 0.2-mW 2-Mb/s Digital Transceiver Based on Wideband Signaling for Human Body Communications", IEEE JSSC 42(9):2021-2033, 2007, doi:10.1109/JSSC.2007.903080. Not in bib (related: song20062mb, song20070).
73. N. Cho, L. Yan, J. Bae, H.-J. Yoo, "A 60 kb/s-10 Mb/s Adaptive Frequency Hopping Transceiver for Interference-Resilient Body Channel Communication", IEEE JSSC 44(3):708-717, 2009, http://ssl.kaist.ac.kr/2007/data/journal/%5BJSSC2009%5DNJCHO.pdf. bib: cho200860kb (ISSCC 2008 version).
74. H. Lee et al., "A 5.5 mW IEEE-802.15.6 WBAN Standard Transceiver for Multichannel Electro-Acupuncture", ISSCC 2013, https://ieeexplore.ieee.org/document/6487811/.
75. H. Cho, H. Lee, J. Bae, H.-J. Yoo, "A 5.2 mW IEEE 802.15.6 HBC Standard Compatible Transceiver With Power Efficient Delay-Locked-Loop Based BPSK Demodulator", IEEE JSSC 50(11):2549-2559, 2015, doi:10.1109/JSSC.2015.2475179.
76. S. Maity, B. Chatterjee, G. Chang, S. Sen, "BodyWire: A 6.3-pJ/b 30-Mb/s -30-dB SIR-Tolerant Broadband Interference-Robust Human Body Communication Transceiver", IEEE JSSC 54(10):2892, 2019, arXiv:1804.06743.
77. Y. Jeon et al., "A 100 Mb/s Galvanically-Coupled Body-Channel-Communication Transceiver with 4.75 pJ/b TX and 26.8 pJ/b RX for Bionic Arms", Symp. VLSI Circuits 2019, https://ieeexplore.ieee.org/document/8778040. bib: Jeon2019A.
78. J. Bae et al., "A 0.24-nJ/b Wireless Body-Area-Network Transceiver With Scalable Double-FSK Modulation", IEEE JSSC 47(1):310-322, Jan 2012, http://ssl.kaist.ac.kr/2007/data/journal/[2012_JSSC]JoonsungBae.pdf. Not in bib (bae2012low is a different paper, JSSC 47(11), Nov 2012).
79. Pitou et al., "First Step Towards a Fully-Miniaturized Intra-Body Communication Transceiver Based on Galvanic Coupling", 2023 (HAL lirmm-03866901, not read). bib: Pitou2023First.
80. Purdue SparcLab publication list (Maity et al., Sub-µWRComm, JSSC 2021; Chatterjee et al., ISSCC 2022; Chatterjee et al., RFIC 2020), https://engineering.purdue.edu/~shreyas/SparcLab/publications.html.
81. Ixana, Wi-R BAN technology page, https://www.ixana.ai/technology/wi-r-ban.
82. Energizer, CR2032 product data sheet, Form 2032NA0618, https://data.energizer.com/pdfs/cr2032.pdf.
83. M. Jensen, "Coin cells and peak current draw", TI white paper SWRA349, https://www.ti.com/lit/wp/swra349/swra349.pdf.
84. K. Y. Yazdandoost, K. Sayrafian-Pour, "Channel Model for Body Area Network (BAN)", IEEE P802.15-08-0780-12-0006, Nov 2010, https://mentor.ieee.org/802.15/dcn/08/15-08-0780-12-0006-tg6-channel-model.pdf.
85. M. M. Alam, E. Ben Hamida, "Realistic Simulation for Body Area and Body-To-Body Networks", Sensors 16(4):561, 2016, PMC4851075.
86. R. Rosini, R. Verdone, R. D'Errico, "Body-to-Body Indoor Channel Modeling at 2.45 GHz", IEEE TAP 62(11):5807-5819, 2014, doi:10.1109/TAP.2014.2352631 (abstract only).
87. Recommendation ITU-R P.1238-8 (07/2015).
88. S. K. Yoo, S. L. Cotton, P. C. Sofotasios, S. Freear, IEEE TWC 15(8):5231-5244, 2016, doi:10.1109/TWC.2016.2555795.
89. A. Kara, H. L. Bertoni, "Blockage/shadowing and polarization measurements at 2.45 GHz for interference evaluation between Bluetooth and IEEE 802.11 WLAN", IEEE APS 2001, doi:10.1109/APS.2001.960112.
90. S. L. Cotton et al., "A review of radio channel models for body centric communications", Radio Science 49, 2014, PMC4579845 (summarises Hanlen et al. 2010, doi:10.1007/s10776-010-0123-z, and Alomainy et al.).
91. S. L. Cotton, W. G. Scanlon, IEEE TAP 57(4):980-990, 2009, doi:10.1109/TAP.2009.2014576.
92. B. Latre et al., "A survey on wireless body area networks", Wireless Networks 17, 2011, https://biblio.ugent.be/publication/3234782/file/3234793.pdf.
93. E. Reusens et al., IEEE TITB 13(6):933-945, 2009, doi:10.1109/TITB.2009.2033054 (via reviews, including Sensors 18(12):4165, 2018, PMC6308834).
94. AIA/FGI, Guidelines for Construction and Equipment of Hospital and Medical Facilities, 1992-93, https://www.fgiguidelines.org/wp-content/uploads/2016/11/1992-93guidelines.pdf.
95. Department of Health, Health Building Note 04-01, Adult in-patient facilities, 2009, https://www.england.nhs.uk/wp-content/uploads/2009/12/HBN_04-01_Final.pdf.
96. NYC Mechanical Code 2022, Table 403.3.1.1 (from IMC / ASHRAE 62.1), https://up.codes/viewer/new_york_city/nyc-mechanical-code-2022/chapter/4/ventilation.
97. M. Nath, S. Maity, S. Avlani, S. Weigand, S. Sen, "Inter-body coupling in electro-quasistatic human body communication", Scientific Reports 11, 2021, doi:10.1038/s41598-020-79788-9. bib: Nath2020Inter-body.
98. Tomlinson, Banou, Yu, Nogueira, Chowdhury, "Secure on-skin biometric signal transmission using galvanic coupling", IEEE INFOCOM 2019, https://par.nsf.gov/servlets/purl/10130328. Journal companion in bib: Tomlinson2019Body-Guided.
99. Hwang et al., "Analysis of signal interference in human body communication using human body as transmission medium", IEEE APS 2006, pp. 495-498. bib: hwang2006analysis.
100. This repository: `python_simulation/HBCProtocol.py`, `MAC_layer.py`, `PHY_layer_BCC_simulation.py`, `PHY_layer_BCC_simulation_v2.py`, `BLE_SIM.py`, `protocol_SDR_test/data_to_send`, `an earlier script (v2_packet_run)`, `an earlier script (fdma_ceiling)`, `an earlier script (link_budget)`, and earlier work chapters ch03 (sec01_hardware), ch05 (sec03_physical_layer, sec04_data_link, sec07_sdr_gnuradio, sec08_evaluation, sec08b_link_budget).
101. M. Zuniga, B. Krishnamachari, "Analyzing the transitional region in low power wireless links", IEEE SECON 2004, pp. 517-526 (packet form PRR = (1 - BER)^(8f); its radio model is non-coherent FSK).
102. Ancans et al., "Bluetooth low energy throughput in densely deployed radio environment", 2019. bib: ancans2019bluetooth.

Related entries already in `biblio.bib`, not used for a value here: siekkinen2012low, karvonen2019experimental, feng2015analysis (BLE and 802.15.4 energy and performance); song20070, song20062mb, bae2012low, choi200624, choi2006low, lin20100, Cho2007The, Shi2022Galvanic-Coupled, Chen2020Design (BCC transceivers); cao2015hospital, hwang2006simulation, Nath2020Toward, tomlinson2018comprehensive, li2017modeling, Maity2020On (channel and coupling).

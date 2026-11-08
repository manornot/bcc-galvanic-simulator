#!/bin/bash
# Regenerates every sweep behind REPORT.md. Run from python_simulation/coexistence/.
# Raw per-BAN rows go to results/raw/ (git-ignored); every sweep is resumable, so an interrupted run is restarted with
# the same command. usage: bash run_sweeps.sh [main|sens|ladder|rates|dense|pool|touch|spread|lossfree|pooled|back|burst|span|contactn|touch_v2|touch_v2_diag|claim4|all] [workers]
# (default: all, 4 workers; at most 8 in total)
set -e
what=${1:-all}
w=${2:-4}
R_ALL=0.02,0.1,0.5,1,5,10,30,50,100,200,500,1000
R_SENS=1,10,30,100,200,500
K=1,2,4,8,16,30
S="python sweep.py --k $K --n-dev 1,4,10 --ci-modes latency,energy --workers $w"
# v3. The raw results of (fair comparison) are in results/raw/v1/, those of in
# results/raw/v2/. Before each sweep the v2 rows of the series left unchanged are copied (checked row
# for row against the merge-base: BLE 1M, 2M, Coded S=8, 802.15.4 CSMA, 802.15.6 NB 121.4, HBC SF 8, every BCC row at
# 500 kchip/s and the 50 kchip/s contact case, ble4thesis), the as-built BCC and ble4thesis relabelled 'thesis'
#, so the resumable sweep reruns only TSCH, ANT, SmartBAN, the 802.15.4 GTS and the new rates.
seed() {
  [ -f results/raw/v2/$1_rows.csv ] && [ ! -f results/raw/$1_rows.csv ] && awk -F, -v OFS=, 'NR == 1 { print; next }
    ($1 == "ble" && ($3 == "1M" || $3 == "2M" || $3 == "coded_s8")) || ($1 == "ieee802154" && $3 == "oqpsk") ||
    ($1 == "ieee802156nb" && $3 == "121.4") || ($1 == "hbc" && $3 == "sf8") || $1 == "bcc" || $1 == "ble4thesis" {
      sub(/^fair/, "thesis", $2); print }' results/raw/v2/$1_rows.csv > results/raw/$1_rows.csv
  return 0
}

if [ "$what" = main ] || [ "$what" = all ]; then
  # every catalogue protocol and variant, both reception models, ward venue, typical density and channel case
  seed main; seed s0_ble4thesis
  $S --reception collision,sinr --r $R_ALL --trials 5 --out results/raw/main_rows.csv --cells results/raw/main_cells.csv
  # S0: the plain BLE model (BLE_SIM rules inside the engine), for the side-by-side with the as-built table
  $S --protocols ble4thesis --reception collision --r $R_ALL --trials 5 --out results/raw/s0_ble4thesis_rows.csv \
     --cells results/raw/s0_ble4thesis_cells.csv
fi

if [ "$what" = sens ] || [ "$what" = all ]; then
  for n in s1_nocifloor s3_gym s4_case_low s4_case_high s5_bcc_contact s5_bcc_fcs s5_bcc_coupling40 \
           s5_bcc_coupling60 s5_bcc_power_low s5_bcc_power_high; do seed $n; done
  # S1: no co-channel C/I floor (interference as Gaussian noise), so FEC can act on weak interference
  $S --reception sinr --ci-floor off --r $R_SENS --trials 3 --out results/raw/s1_nocifloor_rows.csv \
     --cells results/raw/s1_nocifloor_cells.csv
  # S3: crowded gym instead of a ward
  $S --reception sinr --venue gym --density high --r $R_SENS --trials 3 --out results/raw/s3_gym_rows.csv \
     --cells results/raw/s3_gym_cells.csv
  # S4: optimistic and pessimistic channel sets (on-body and body-to-body loss, HBC coupling)
  for c in low high; do
    python sweep.py --reception sinr --case $c --r $R_SENS --k 1,8,30 --n-dev 1,4,10 --ci-modes latency,energy \
      --trials 3 --workers $w --out results/raw/s4_case_${c}_rows.csv --cells results/raw/s4_case_${c}_cells.csv
  done
  # S5: BCC cases (BCC full on one tone and multi-tone, the as-built stack; sweep.py skips the entries get_spec
  # refuses: multi-tone has no contact or coupling case, BCC full no FCS switch): measured-hardware contact chain
  # (50 kchip/s), payload FCS, inter-body coupling, transceiver power (the chip rates are the 'rates' section)
  B="$S --protocols bcc --r $R_ALL --trials 5"
  $B --bcc-contact --reception collision,sinr --out results/raw/s5_bcc_contact_rows.csv \
     --cells results/raw/s5_bcc_contact_cells.csv
  $B --bcc-fcs --reception collision,sinr --out results/raw/s5_bcc_fcs_rows.csv --cells results/raw/s5_bcc_fcs_cells.csv
  for c in 40 60; do
    $B --bcc-coupling $c --reception sinr --out results/raw/s5_bcc_coupling${c}_rows.csv \
       --cells results/raw/s5_bcc_coupling${c}_cells.csv
  done
  for p in low high; do
    $B --bcc-power $p --reception collision --out results/raw/s5_bcc_power_${p}_rows.csv \
       --cells results/raw/s5_bcc_power_${p}_cells.csv
  done
fi
if [ "$what" = ladder ] || [ "$what" = all ]; then
  for n in ladder ladder_contact; do seed $n; done
  # galvanic-optimisation ladder of BCC full (one tone and multi-tone): levels 1 to 5 (level 0 is the main run),
  # levels 3 to 5 also with the retry pool sized for the measured ~30 % contact loss (share 0.43); the contact case
  # (one tone, 50 kchip/s; multi-tone has none and is skipped) at every level and both pools, 10 trials
  L="$S --protocols bcc --variants full --r $R_ALL"
  for o in 0 1 2 3 4 5; do
    [ $o -gt 0 ] && $L --bcc-opt $o --reception collision,sinr --trials 5 --out results/raw/ladder_rows.csv \
      --cells results/raw/ladder_cells.csv
    $L --bcc-opt $o --bcc-contact --reception collision,sinr --trials 10 --out results/raw/ladder_contact_rows.csv \
      --cells results/raw/ladder_contact_cells.csv
  done
  for o in 3 4 5; do
    $L --bcc-opt $o --bcc-pool 0.43 --reception collision,sinr --trials 5 --out results/raw/ladder_rows.csv \
      --cells results/raw/ladder_cells.csv
    $L --bcc-opt $o --bcc-pool 0.43 --bcc-contact --reception collision,sinr --trials 10 \
      --out results/raw/ladder_contact_rows.csv --cells results/raw/ladder_contact_cells.csv
  done
fi
if [ "$what" = pool ]; then
  # slack-filling retry pool (REPORT 5.8): one tone, collision, k 1 and 8 (isolated bodies); not in 'all'
  P="python sweep.py --k 1,8 --n-dev 1,4,10 --ci-modes latency,energy --workers $w --protocols bcc --variants full
     --phy bcc=sf8 --reception collision --r $R_ALL"
  O=results/raw/t_l074
  for kw in "" "--bcc-opt 5" "--bcc-opt 5 --bcc-pool 0.43" "--bcc-opt 5 --bcc-pool-fill" \
            "--bcc-opt 5 --bcc-pool 0.43 --bcc-pool-fill"; do
    $P $kw --bcc-contact --trials 10 --out $O/contact_rows.csv --cells $O/contact_cells.csv
  done
  for rc in 500e3 2.625e6; do
    for kw in "--bcc-opt 5" "--bcc-opt 5 --bcc-pool-fill"; do
      $P $kw --bcc-rc $rc --trials 5 --out $O/clean_rows.csv --cells $O/clean_cells.csv
    done
  done
fi
if [ "$what" = touch ]; then
  # touching bodies and the fallback (REPORT 5.9): the galvanic profile (L5, fill), SINR, k 2 and 8; not in 'all'
  P="python sweep.py --k 2,8 --n-dev 1,4,10 --ci-modes latency,energy --workers $w --protocols bcc --variants full
     --phy bcc=sf8 --reception sinr --r 0.1,1,10,30,100 --trials 10 --bcc-opt 5 --bcc-pool-fill"
  O=results/raw/t_l076
  for base in "--bcc-rc 2.625e6" "--bcc-contact"; do
    $P $base --out $O/rows.csv --cells $O/cells.csv
    for sh in 0.2 1.0; do for db in 0 10; do for fb in "" "--bcc-fallback"; do
      $P $base --bcc-touch $sh --bcc-touch-db $db $fb --out $O/rows.csv --cells $O/cells.csv
    done; done; done
  done
fi
if [ "$what" = spread ]; then
  # payload spreading (REPORT 5.10): the galvanic profile (L5, fill) with Walsh-16, unspread and the per-link
  # switch; contact chain and noise-limited link as in 'pool', touching bodies as in 'touch'; not in 'all'
  P="python sweep.py --n-dev 1,4,10 --ci-modes latency,energy --workers $w --protocols bcc --variants full
     --phy bcc=sf8 --bcc-opt 5 --bcc-pool-fill"
  O=results/raw/t_l075
  for sp in "" "--bcc-spread none" "--bcc-spread adapt"; do
    $P $sp --k 1,8 --reception collision --r $R_ALL --bcc-contact --trials 10 --out $O/contact_rows.csv \
       --cells $O/contact_cells.csv
    for rc in 500e3 2.625e6 9e6; do
      $P $sp --k 1,8 --reception collision --r $R_ALL --bcc-rc $rc --trials 5 --out $O/clean_rows.csv \
         --cells $O/clean_cells.csv
    done
    for db in 0 10; do
      $P $sp --k 2,8 --reception sinr --r 0.1,1,10,30,100 --bcc-rc 2.625e6 --bcc-touch 1.0 --bcc-touch-db $db \
         --trials 10 --out $O/touch_rows.csv --cells $O/touch_cells.csv
    done
  done
fi
if [ "$what" = lossfree ]; then
  # Energy efficiency without contact losses: 50 kchip/s, the HBC stack against the galvanic profile
  # (L5, filling pool), otherwise as the contact sweep (REPORT 5.8); not in 'all'
  P="python sweep.py --k 1,8 --n-dev 1,4,10 --ci-modes latency,energy --workers $w --protocols bcc --variants full
     --phy bcc=sf8 --reception collision --r $R_ALL --bcc-rc 50e3 --trials 10"
  for kw in "" "--bcc-opt 5 --bcc-pool-fill"; do
    $P $kw --out results/raw/t_l091/lossfree_rows.csv --cells results/raw/t_l091/lossfree_cells.csv
  done
fi
if [ "$what" = pooled ]; then
  # pooled reference stack: the HBC stack with retries by a hub-polled pool (share 0.05, the standard's
  # improvised access) instead of a spare exchange per frame, at ladder level 0; the same cells as the series it
  # joins: the ladder (500 kchip/s, one tone and six tones) and its contact case, the 2.625 Mchip/s rate, the
  # loss-free 50 kchip/s energy sweep of; not in 'all'
  L="$S --protocols bcc --variants full --r $R_ALL --bcc-pool 0.05"
  python sweep.py --k 1,8 --n-dev 1,4,10 --ci-modes latency,energy --workers $w --protocols bcc --variants full \
     --phy bcc=sf8 --reception collision --r $R_ALL --bcc-rc 50e3 --trials 10 --bcc-pool 0.05 \
     --out results/raw/t_l091/lossfree_rows.csv --cells results/raw/t_l091/lossfree_cells.csv
  $L --reception collision,sinr --trials 5 --out results/raw/ladder_rows.csv --cells results/raw/ladder_cells.csv
  $L --bcc-rc 2.625e6 --reception collision,sinr --trials 5 --out results/raw/rates_rows.csv \
     --cells results/raw/rates_cells.csv
  $L --bcc-contact --reception collision,sinr --trials 10 --out results/raw/ladder_contact_rows.csv \
     --cells results/raw/ladder_contact_cells.csv
  # the loss-free 50 kchip/s pair at the low and high transceiver power presets (claim4.py E)
  for p in low high; do
    for kw in "--bcc-pool 0.05" "--bcc-opt 5 --bcc-pool-fill"; do
      python sweep.py --k 1,8 --n-dev 1,4,10 --ci-modes latency,energy --workers $w --protocols bcc --variants full \
         --phy bcc=sf8 --reception collision --r $R_ALL --bcc-rc 50e3 --trials 10 --bcc-power $p $kw \
         --out results/raw/t_l091/power_rows.csv --cells results/raw/t_l091/power_cells.csv
    done
  done
fi
if [ "$what" = back ]; then
  # the lead reference, the pooled reference stack with the standard's block
  # acknowledgement (--bcc-back: L-Ack frames, one B-Ack per allocation, the pool exchanges with their I-Ack), one tone,
  # collision model, in the cells of the pooled series it is paired with (claim4.py H2): the loss-free 50 kchip/s sweep,
  # 500 kchip/s and 2.625 Mchip/s; not in 'all'
  python sweep.py --k 1,8 --n-dev 1,4,10 --ci-modes latency,energy --workers $w --protocols bcc --variants full \
     --phy bcc=sf8 --reception collision --r $R_ALL --bcc-rc 50e3 --trials 10 --bcc-pool 0.05 --bcc-back \
     --out results/raw/back/back_rows.csv --cells results/raw/back/back_cells.csv
  for rc in 500e3 2.625e6; do
    $S --protocols bcc --variants full --phy bcc=sf8 --r $R_ALL --bcc-rc $rc --bcc-pool 0.05 --bcc-back \
       --reception collision --trials 5 --out results/raw/back/back_rows.csv --cells results/raw/back/back_cells.csv
  done
fi
if [ "$what" = burst ]; then
  # the lead reference, the block-ack reference above in the HBC burst mode
  # (--bcc-burst: frames 2.. of an allocation pMIFS apart in RI mode, 802.15.6 6.2.9.4, 10.5.2), the same cells
  # (claim4.py H2); not in 'all'
  python sweep.py --k 1,8 --n-dev 1,4,10 --ci-modes latency,energy --workers $w --protocols bcc --variants full \
     --phy bcc=sf8 --reception collision --r $R_ALL --bcc-rc 50e3 --trials 10 --bcc-pool 0.05 --bcc-back --bcc-burst \
     --out results/raw/burst/burst_rows.csv --cells results/raw/burst/burst_cells.csv
  for rc in 500e3 2.625e6; do
    $S --protocols bcc --variants full --phy bcc=sf8 --r $R_ALL --bcc-rc $rc --bcc-pool 0.05 --bcc-back --bcc-burst \
       --reception collision --trials 5 --out results/raw/burst/burst_rows.csv --cells results/raw/burst/burst_cells.csv
  done
fi
if [ "$what" = span ]; then
  # the lead reference, the burst reference above with one block spanning its
  # allocations (--bcc-span 1000000000: no B-Ack, 802.15.6 6.2.9.4, 5.4.2.4), the same cells (claim4.py J2),
  # with the '128 B' pilot period (--bcc-pilot128, 10.7.3 Table 86); not in 'all'
  python sweep.py --k 1,8 --n-dev 1,4,10 --ci-modes latency,energy --workers $w --protocols bcc --variants full \
     --phy bcc=sf8 --reception collision --r $R_ALL --bcc-rc 50e3 --trials 10 --bcc-pool 0.05 --bcc-back --bcc-burst \
     --bcc-span 1000000000 --bcc-pilot128 --out results/raw/span/span_rows.csv --cells results/raw/span/span_cells.csv
  for rc in 500e3 2.625e6; do
    $S --protocols bcc --variants full --phy bcc=sf8 --r $R_ALL --bcc-rc $rc --bcc-pool 0.05 --bcc-back --bcc-burst \
       --bcc-span 1000000000 --bcc-pilot128 --reception collision --trials 5 --out results/raw/span/span_rows.csv \
       --cells results/raw/span/span_cells.csv
  done
fi
if [ "$what" = contactn ]; then
  # the contact chain (50 kchip/s, REPORT 5.8 settings) by number of sensors,
  # both stacks in their current framing: the galvanic profile (L5, filling pool) and the lead reference
  # (block ack in burst mode, one block across the allocations, '128 B' pilots, pool 0.05; the filling pool needs
  # ladder level 3), base seeds 0 and 1 in separate files; at most 2 workers; REPORT 5.14; not in 'all'
  P="python sweep.py --k 1,8 --n-dev 1,4,10 --ci-modes latency,energy --workers $w --protocols bcc --variants full
     --phy bcc=sf8 --reception collision --r $R_ALL --bcc-contact --trials 10"
  for sd in 0 1; do
    for kw in "--bcc-opt 5 --bcc-pool-fill" \
              "--bcc-pool 0.05 --bcc-back --bcc-burst --bcc-span 1000000000 --bcc-pilot128"; do
      $P $kw --seed $sd --out results/raw/contactn/s${sd}_rows.csv --cells results/contactn/s${sd}_cells.csv
    done
  done
fi
if [ "$what" = touch_v2 ]; then
  # touching bodies at 2.625 Mchip/s on one simulator version: shifted chip curve, current framing, 5 % pool on
  # both stacks, latency mode (28 feasible of 30 cells). Then the interleaving hubs at 0 dB and the unspread payload at
  # 10 dB, the other figures of the same sentence. Not in 'all'.
  P="python sweep.py --k 2,8 --n-dev 1,4,10 --ci-modes latency --workers $w --protocols bcc --variants full --phy bcc=sf8
     --reception sinr --r 0.1,1,10,30,100 --bcc-rc 2.625e6 --bcc-touch 1.0 --trials 10 --bcc-pool 0.05"
  PROF="--bcc-opt 5"
  REF="--bcc-back --bcc-burst --bcc-span 1000000000 --bcc-pilot128"
  mkdir -p results/raw/touch_v2 results/touch_v2
  OUT="--out results/raw/touch_v2/rows.csv --cells results/touch_v2/cells.csv"
  for db in 0 10; do
    for kw in "$REF" "$PROF"; do $P $kw --bcc-touch-db $db $OUT; done
  done
  for kw in "$REF" "$PROF"; do $P $kw --bcc-touch-db 0 --bcc-fallback $OUT; done
  $P $PROF --bcc-spread none --bcc-touch-db 10 $OUT
fi
if [ "$what" = touch_v2_diag ]; then
  # stop: what moved the 10 dB comparison. The same cells with the profile's filling pool, and with base seed 1.
  P="python sweep.py --k 2,8 --n-dev 1,4,10 --ci-modes latency --workers $w --protocols bcc --variants full --phy bcc=sf8
     --reception sinr --r 0.1,1,10,30,100 --bcc-rc 2.625e6 --bcc-touch 1.0 --trials 10"
  REF="--bcc-pool 0.05 --bcc-back --bcc-burst --bcc-span 1000000000 --bcc-pilot128"
  for db in 0 10; do
    $P --bcc-opt 5 --bcc-pool-fill --bcc-touch-db $db --out results/raw/touch_v2/fill_rows.csv --cells results/touch_v2/fill_cells.csv
  done
  for kw in "$REF" "--bcc-opt 5 --bcc-pool 0.05" "--bcc-opt 5 --bcc-pool-fill"; do
    $P $kw --bcc-touch-db 10 --seed 1 --out results/raw/touch_v2/s1_rows.csv --cells results/touch_v2/s1_cells.csv
  done
  # the printed 10 dB figures were the interleaving row: that row on the new version
  for kw in "$REF" "--bcc-opt 5 --bcc-pool 0.05"; do
    $P $kw --bcc-touch-db 10 --bcc-fallback --out results/raw/touch_v2/fb10_rows.csv --cells results/touch_v2/fb10_cells.csv
  done
fi
if [ "$what" = claim4 ]; then
  # Capacity in closed form, the short frames and both stacks simulated at the admitted capacity,
  # the energy pairs and their reruns at equal period (results/tables/t15_claim4.md; E needs 'pooled'); not in 'all'
  python claim4.py --workers $w
fi
if [ "$what" = rates ] || [ "$what" = all ]; then
  # chip rates of BCC (the as-built sweep and the IEEE 802.15.6 HBC Walsh-chip rates, then upward to the 9 Mchip/s
  # stability limit; 500 kchip/s is the main run): BCC full one tone and multi-tone (where the tone plan has tones), the
  # as-built stack, and the galvanic proposal (ladder level 5, pool 0.05)
  for rc in 50e3 200e3 656.25e3 1.3125e6 2e6 2.625e6 5.25e6 7e6 8e6 9e6; do
    $S --protocols bcc --bcc-rc $rc --reception collision,sinr --r $R_ALL --trials 5 \
       --out results/raw/rates_rows.csv --cells results/raw/rates_cells.csv
    $S --protocols bcc --variants full --bcc-opt 5 --bcc-rc $rc --reception collision,sinr --r $R_ALL --trials 5 \
       --out results/raw/rates_rows.csv --cells results/raw/rates_cells.csv
  done
fi
if [ "$what" = dense ] || [ "$what" = all ]; then
  # S6 dense sensor networks: 20, 30 and 50 sensors per body on a reduced grid, every stack, and BCC at 500 kchip/s,
  # 5.25 and 9 Mchip/s with the proposal
  D="python sweep.py --k 1,8,30 --n-dev 20,30,50 --ci-modes latency,energy --reception collision,sinr"
  D="$D --r 0.1,1,5,10,30,50,100 --trials 3 --workers $w"
  D="$D --out results/raw/s6_dense_rows.csv --cells results/raw/s6_dense_cells.csv"
  $D
  for rc in 500e3 1.3125e6 2e6 2.625e6 5.25e6 9e6; do
    $D --protocols bcc --bcc-rc $rc
    $D --protocols bcc --variants full --bcc-opt 5 --bcc-rc $rc
  done
fi

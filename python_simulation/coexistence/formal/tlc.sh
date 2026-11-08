#!/bin/sh
# One TLC run of GalvanicMAC.tla: tlc.sh NAME N n M A P B R CAP Q Budget Variant Track Invariants [HubExpiry HubWindow]
# Budget -1: unbounded losses. e.g. tlc.sh age_A2_M5 1 1 5 2 1 1 2 4 3 -1 age '{"falseAck","alias"}' 'NoFalseAck NoAlias'
# writes tlc/NAME.cfg and tlc/NAME.out (TLC 2.19 from ${TLA2TOOLS:-tla2tools.jar}, Java 21, 2 workers)
set -e
cd "$(dirname "$0")"
mkdir -p tlc
NAME=$1; shift
if [ "${10}" -lt 0 ]; then BOUNDED=FALSE; BUDGET=0; else BOUNDED=TRUE; BUDGET=${10}; fi
cat > tlc/$NAME.cfg <<CFG
SPECIFICATION Spec
CONSTANTS N = $1  n = $2  M = $3  A = $4  P = $5  B = $6  R = $7  CAP = $8  Q = $9  Bounded = $BOUNDED  Budget = $BUDGET
  Variant = "${11}"  Track = ${12}  HubExpiry = ${14:-TRUE}  HubWindow = ${15:-TRUE}
INVARIANTS ${13}
CHECK_DEADLOCK FALSE
CFG
META=${TLC_META:-/tmp/tlc-meta}
mkdir -p "$META"
java -Xmx3g -XX:+UseParallelGC -cp "${TLA2TOOLS:-tla2tools.jar}" tlc2.TLC -workers 2 -metadir "$META/$NAME" \
  -config tlc/$NAME.cfg GalvanicMAC.tla > tlc/$NAME.out 2>&1 || true
grep -E "Version|Model checking completed|violated|distinct states|Finished in|Error|depth of the complete" tlc/$NAME.out | head -8

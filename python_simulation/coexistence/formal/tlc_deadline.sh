#!/bin/sh
# Deadline A superframes under any A losses (age rule, M = 2A + 1, queue Q = A + 1, pool P = 1), N = 1, 2, ... until a run passes
# about 10 minutes: tlc_deadline.sh A
A=$1
for N in 1 2 3 4 5 6 7 8 9 10; do
  S=$(date +%s)
  timeout 660 ./tlc.sh dl_A${A}_L${A}_N$N $N 1 $((2*A+1)) $A 1 1 2 $((A+2)) $((A+1)) $A age '{"falseAck","alias","miss"}' 'NoFalseAck NoAlias NoMiss' | grep -E "completed|violated|distinct states found, 0|Finished in"
  E=$(( $(date +%s) - S )); echo "A $A N $N wall ${E}s"
  [ $E -ge 600 ] && break
done

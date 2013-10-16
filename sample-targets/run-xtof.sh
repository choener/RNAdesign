#!/bin/bash

echo "xtof"

#  --optfun "-Ged" \
#  --optfun "eos(1)+eos(2)+eos(3) - 3*gibbs + 1 * ((eos(1)-eos(2))^2 + (eos(1)-eos(3))^2 + (eos(2)-eos(3))^2)" \
#  --optfun "-sum(ed,all)" \
cat sample-targets/xtof | ./.cabal-sandbox/bin/RNAdesign -n 50 --thin 200 -b 100 --scale 1 \
  --optfun "eos(1)+eos(2)+eos(3) - 3*gibbs + 1 * ((eos(1)-eos(2))^2 + (eos(1)-eos(3))^2 + (eos(2)-eos(3))^2)" \
  > xtof.out
cat xtof.out | head -n 40

for i in `cat xtof.out | head -n 7 | tail -n 5 | awk '{print $1}'`
do
  echo $i | RNAfold -p

  cat dot.ps | ./.cabal-sandbox/bin/RNAdotplot \
    -s "((((....))))....((((....))))........" \
    -s "........((((....((((....))))....))))" \
    -s "((((((((....))))((((....))))....))))" \
  --colorize 1,red --colorize 2,green --colorize 3,blue --colorize 1-2,yellow --colorize 1-3,magenta --colorize 2-3,cyan > xtof.$i.ps

done


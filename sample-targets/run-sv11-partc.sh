#!/bin/bash

cat sample-targets/sv11

cat sample-targets/sv11 | tail -n2 | ./.cabal-sandbox/bin/RNAdesign -n 50 --thin 200 -b 100 --scale 2 \
  --optfun "partc(1) + partc(2) + eos(1)+eos(2) - 2*gibbs + 1 * ((eos(1)-eos(2))^2)" \
  > sv11.out

head -n 40 sv11.out

for i in `cat sv11.out | head -n 7 | tail -n5 | awk '{print $1}'`
do
  echo $i | RNAfold -p

  cat dot.ps | ./.cabal-sandbox/bin/RNAdotplot \
    -s "(((((((((((...)))))))..((((((((((....))))))))))........)))).....((((((((........)))))))).((((((((.....))))))))....." \
    -s "(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..)))." \
    --colorize 1,red --colorize 2,blue --colorize 1-2,magenta > sv11.$i.ps
done

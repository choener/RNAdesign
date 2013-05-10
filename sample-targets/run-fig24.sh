#!/bin/bash

cat tests/fig24 | ./cabal-dev/bin/RNASequenceDesign -n 50 --thin 500 -b 100 --scale 1 --optfun "Pdef" > fig24.out # -logMN(0.25,0.25,0.25,0.25)" > fig24.out
cat fig24.out | head -n 40
cat fig24.out | head -n 3 | tail -n1 | awk '{print $1}' | RNAfold -p

for i in `cat fig24.out | head -n 7 | tail -n 5 | awk '{print $1}'`
do
  echo $i | RNAfold -p

  cat dot.ps | ./cabal-dev/bin/RNAdotplot \
    -s "..((((((((....)))..((((....))))..))))).." \
    -s ".....((((.((((...(((....)))..))))))))..." \
    -s ".(((.(((..(((...)))..((((....))))))))))." \
  --colorize 1,red --colorize 2,green --colorize 3,blue --colorize 1-2,yellow --colorize 1-3,magenta --colorize 2-3,cyan > fig24.$i.ps

done


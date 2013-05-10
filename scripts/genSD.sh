#!/bin/bash

# $1 is the sequence design program
# $2 is where the turner params live

for i in *.constraints
do
  IC=`basename $i .constraints`.candidate
  if [[ ! -f $IC ]]
  then
    cat $i | $1 --turner $2 --number 50 --thin 500 --burnin 1000 --scale 2 \
      --optfun "eos(1)+eos(2)+eos(3) - 3*gibbs + 1 * ((eos(1)-eos(2))^2 + (eos(1)-eos(3))^2 + (eos(2)-eos(3))^2)" \
      | head -n 3 | tail -n 1 > $IC
  fi
done

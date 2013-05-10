#!/bin/bash

DONE=0
until [[ $DONE = 1 ]]
do
  randseq -l $1 > tmp
  NAME=`cat tmp`
  cat tmp | $3 > $NAME.rnashapes
  tail -n+2 $NAME.rnashapes | awk '{print $2}' > $NAME.constraints
  rm tmp
  if [[ `wc -l $NAME.constraints | awk '{print $1}'` = $2 ]]
  then
#    echo "$NAME looks good"
    DONE=1
  else
#    cat $NAME.constraints
    rm $NAME.rnashapes
    rm $NAME.constraints
  fi
done

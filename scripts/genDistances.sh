#!/bin/bash

# $1 is the subopt distance tool

rm score-1
rm score-2
rm score-3

for i in *.candidate
do
  B=`basename $i .candidate`
  cat $i | awk '{print $1}' | RNAsubopt -e 5 > tmp
  rm -f tmp.scores.pre
  for j in `cat $B.constraints`
  do
    cat tmp | $1 $2 $3 "$j" | awk '{printf "%s %s\n", $1, $4}' >> tmp.scores.pre
  done
  cat tmp.scores.pre | sort -n -k 2 > tmp.scores
  rm tmp.scores.pre
  sed -n '1p' tmp.scores >> score-1
  sed -n '2p' tmp.scores >> score-2
  sed -n '3p' tmp.scores >> score-3
  tail -n 1 score-1
  tail -n 1 score-2
  tail -n 1 score-3
  rm tmp
  rm tmp.scores
done

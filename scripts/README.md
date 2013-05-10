HELPER SCRIPTS
==============

This directory contains four helper scripts:

- anon.sh : ``anonymizes'' targets. Just copies targets from the constraints
  directory and numbers them. Afterwards, no trace of the original /sequence/
  remains
- genConstraints.sh : creates artificial constraints, starts with a random
  sequence, sends it into RNAshapes, finally we have a K-structure target
- genDistances.sh : calculate distance of each constraint from mfe
- genSD.sh : runs the sequence designer on each target

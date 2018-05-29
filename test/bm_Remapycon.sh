#!/bin/sh
#
export DYLD_LIBRARY_PATH=/opt/intel/lib:$DYLD_LIBRARY_PATH
CDO=src/cdo
#
#
$CDO -V
#
export CDO_REMAP_RADIUS=.1
for N in 3 4 5 6 7 8 9; do
#for N in 3 4 5; do
  IFILE=topor2b$N
  $CDO -f nc topo,icor2b$N $IFILE
  OFILE=${IFILE}ycon
  time $CDO -P 8 remapycon,global_0.2 $IFILE ${OFILE}
done
#
# CDO version 1.9.5rc2
#
#1 bailung: gcc 7.2 sse4_2
#2   8 cores
#3   spherepart
#4   spherepart 8 cores
#       3   4   5    6    7     8      9     10
#1     24  34  71  194  874  x6680
#2      5   6  11   28  126   2005   6660
#x2*8    3   4   8   22   76   533   1808
#x3     27  29  33   40   59   110    280

#1 mistral: intel 18.0.1 avx2 (interactiv on login node)
#       3   4   5    6    7     8      9     10
#1     15  22  54  161  568  4454
#2      3   4   8   22   76   462   2692
#3     17  19  22   29   48    92    235    774
#4      4   4   4    6   10    25     86    351

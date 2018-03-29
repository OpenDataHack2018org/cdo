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
  $CDO remapycon,global_0.2 $IFILE ${OFILE}
done
#
# CDO version 1.9.4rc1
#
#1 bailung: gcc 7.2 sse4_2
#2 mistral: intel 18.0.1 avx2 (interactiv)
#       3   4   5    6    7     8      9
#1     25  35  72  200  770  6680
#2     15  23  55  162  588  4492
#2*8    3   4   8   22   76   533   1808


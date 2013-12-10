#!/bin/sh
#
CDO=src/cdo
#
RMODS="bil bic nn con"
#
cdo -setrtomiss,0,10000 -topo topo05
cdo -setrtomiss,0,10000 -remapbil,global_0.05 -topo topo005
#
IFILE=topo005
#
for RMOD in $RMODS; do
  echo "remap $RMOD:"
  OTYPE=0
  for GRIDTYPE in " " "-setgridtype,curvilinear"; do
    OTYPE=`expr $OTYPE + 1`
    OFILE=${IFILE}_${RMOD}_${OTYPE}
    $CDO remap${RMOD},global_0.5 ${GRIDTYPE} $IFILE ${OFILE}
  done
  cdo diff ${IFILE}_${RMOD}_?
done
#
# result on bailung:
# =================
#         bil    bic    nn   con
# reg2d   4.2    5.1     7   223
# curv   38.8   39.0   842   223
#
# result on blizzard:
# ==================
#         bil    bic    nn   con
# reg2d    15     19    18   
# curv    101    104         
#
#
#

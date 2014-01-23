#!/bin/sh
#
CDO=src/cdo
#
RMODS="bil bic nn con cons"
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
  cdo diff2 ${IFILE}_${RMOD}_?
done
#
# result on hama:
# =================
#         bil    bic    nn   con   cons
# reg2d   1.0    2.1     2   256    265
# curv   43.8   45.0   982   300    356
#
# result on bailung:
# =================
#         bil    bic    nn   con   cons
# reg2d   0.7    1.6     1   223    227
# curv   36.9   37.9   843   227    305
#
# result on blizzard:
# ==================
#         bil    bic    nn   con   cons
# reg2d    15     19    18   528
# curv    101    104  1936   528
#

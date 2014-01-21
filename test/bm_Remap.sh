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
  cdo diff ${IFILE}_${RMOD}_?
done
#
# result on hama:
# =================
#         bil    bic    nn   con   cons
# reg2d   5.5    6.8     9   256    263
# curv   43.8   45.0   982   300
#
# result on bailung:
# =================
#         bil    bic    nn   con   cons
# reg2d   1.0    1.8     2   223    227
# curv   36.9   37.9   843   227    316
#
# result on blizzard:
# ==================
#         bil    bic    nn   con   cons
# reg2d    15     19    18   528
# curv    101    104  1936   528
#

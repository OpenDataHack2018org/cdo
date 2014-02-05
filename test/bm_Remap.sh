#!/bin/sh
#
CDO=src/cdo
#CDO=cdo-1.6.2
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
####################################################
#
# CDO 1.6.3
#
# result on hama:
# =================
#         bil    bic    nn   con   cons
# reg2d   1.0    2.1     2   238    187
# curv   44.8   47.0  1024   240    280
#
# result on bailung:
# =================
#         bil    bic    nn   con   cons
# reg2d   0.7    1.6     1   204    160
# curv   36.9   37.9   843   210    237
#
# result on blizzard:
# ==================
#         bil    bic    nn   con   cons
# reg2d    15     19    18   528
# curv    101    104  1936   528
#
####################################################
#
# CDO 1.6.2
#
# result on hama:
# =================
#         bil    bic    nn   con
# curv     60     61   873   274
#
# result on bailung:
# =================
#         bil    bic    nn   con
# curv     51     52   754   235
#
# result on blizzard:
# ==================
#         bil    bic    nn   con
# curv    101    104  1936   528
#

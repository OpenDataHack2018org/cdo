#!/bin/sh
#
CDO=src/cdo
#CDO=cdo-1.6.2
#
RMODS="bil bic nn con cons"
#
#cdo -setrtomiss,0,10000 -topo topo05
#cdo -setrtomiss,0,10000 -remapbil,global_0.05 -topo topo005
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
#       con/curv  cons/curv  cons/reg2d
#   1   195s       243s       161s
#   2   129s       130s        90s 
#   4   185s        93s        74s
#      4024MB     5836MB     1024MB
####################################################
#
# CDO 1.6.3
#
# result on hama:
# =================
#         bil    bic    nn   con   cons
# reg2d   0.9    2.1     2   224    184
# curv   43.3   44.3   982   255    275
#
# result on bailung:
# =================
#         bil    bic    nn   con   cons
# reg2d   0.7    1.6     1   194    160
# curv   36.9   37.9   843   200    237
#
# result on blizzard: (interactiv)
# ==================
#         bil    bic    nn   con   cons
# reg2d   2.2    5.5     3   499    342
# curv     88     93  1936   500    553
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

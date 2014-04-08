#!/bin/sh
#
CDO=src/cdo
#CDO=cdo-1.6.2
#
RMODS="bil bic nn con ycon"
#
#cdo -setrtomiss,0,10000 -topo topo05
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
####################################################
#
# CDO 1.6.3
#
#       con/curv  ycon/curv  ycon/reg2d  ycon/reg2d/ICC
#   1   195s       243s       140s        134s
#   2   129s       130s        79s         70s
#   4   185s        93s        73s         39s
#   8                                      24s
#      4024MB     5836MB     1240MB
####################################################
#
# CDO 1.6.4
#
# result on hama: gcc 4.8.2
# =================
#         bil    bic    nn   con   ycon
# reg2d   0.9    2.1   1.3   247    160
# curv     44     46   889   255    348
#
# result on hama: icc 12.1.5
# =================
#         bil    bic    nn   con   ycon
# reg2d   1.0    2.2   1.4   295    131
# curv     51     52   513   305    389
#
# result on bailung: gcc
# =================
#         bil    bic    nn   con   ycon
# reg2d   0.7    1.6     1   205    137
# curv     39     40   703   210    293
#
# result on bailung: icc
# =================
#         bil    bic    nn   con   ycon
# reg2d   0.8    1.7     1   217    131
# curv     45     45   361   204    277
#
# result on blizzard: (interactiv)
# ==================
#         bil    bic    nn   con   ycon
# reg2d   2.2    5.5     3   470    326
# curv     88     93  1682   480    692
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

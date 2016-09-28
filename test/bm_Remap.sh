#!/bin/sh
#
export DYLD_LIBRARY_PATH=/opt/intel/lib:$DYLD_LIBRARY_PATH
CDO=src/cdo
#CDO=cdo-1.6.2
#
RMODS="bil bic nn con ycon"
#
$CDO -V
#$CDO -setrtomiss,0,10000 -topo topo05
$CDO -setrtomiss,0,10000 -remapbil,global_0.05 -topo topo005
#
IFILE=topo005
#
export CDO_REMAP_RADIUS=.1
for RMOD in $RMODS; do
  echo "remap $RMOD:"
  OTYPE=0
  for GRIDTYPE in " " "-setgridtype,curvilinear"; do
    OTYPE=`expr $OTYPE + 1`
    OFILE=${IFILE}_${RMOD}_${OTYPE}
    $CDO remap${RMOD},global_0.5 ${GRIDTYPE} $IFILE ${OFILE}
  done
  $CDO diff ${IFILE}_${RMOD}_?
done
####################################################
#
# CDO 1.6.8 ----------------------------------------
# bailung: gcc 4.9.1 / icc 12.1.5
#
#       con/curv  ycon/curv  ycon/reg2d  ycon/reg2d/ICC
#   1   208s       276s       144s        137s
#   2   142s       149s        72s         70s
#   4   204s        88s        39s         39s
#   8               60s        22s         21s
#  12               55s        16s         15s
#  24                          12s         12s
#
# CDO 1.6.3 ----------------------------------------
#
#       con/curv  ycon/curv  ycon/reg2d  ycon/reg2d/ICC
#   1   195s       243s       140s        134s
#   2   129s       130s        79s         70s
#   4   185s        93s        73s         39s
#   8                                      24s
#      4024MB     5836MB     1240MB
####################################################
#
# CDO 1.8.0 ----------------------------------------
#
# result on hama2: gcc 6.1.0 avx2
# =================
#         bil    bic    nn   con   ycon
# reg2d   0.5    1.1   0.6   138     92
# curv     22     22    58   138    172
#
# result on hama2: icc 16.0.2 avx2
# =================
#         bil    bic    nn   con   ycon
# reg2d   0.6    1.3   0.7   120     90
# curv     21     21    63   120    163
#
# CDO 1.6.8 ----------------------------------------
#
# result on hama2: gcc 5.1.0 avx2
# =================
#         bil    bic    nn   con   ycon
# reg2d   0.6    1.3   0.7   145     90
# curv     24     25    66   146    172
#
# result on hama2: icc 15.0.1 avx2
# =================
#         bil    bic    nn   con   ycon
# reg2d   0.7    1.4   0.8   128     92
# curv     20     21    63   126    161
#
# result on bailung: gcc 4.9.1
# =================
#         bil    bic    nn   con   ycon
# reg2d   0.7    1.6   1.0   202    143
# curv     33     34   103   208    273
#
# result on pitpull: gcc 4.8.2 avx2
# =================
#         bil    bic    nn   con   ycon
# reg2d   0.7    1.4   0.8   175    106
# curv     29     30    69   175    210
#
# result on pitbull: icc 15.0.1 avx2
# =================
#         bil    bic    nn   con   ycon
# reg2d   0.6    1.4   0.8   144     78
# curv     23     24    67   144    154
#
# CDO 1.6.4 ----------------------------------------
#
# result on hama: gcc 4.9.1
# =================
#         bil    bic    nn   con   ycon
# reg2d   1.0    2.1   1.3   236    168
# curv     40     41   164   280    325
#
# result on hama: gcc 4.8.2
# =================
#         bil    bic    nn   con   ycon
# reg2d   0.9    2.1   1.3   247    172
# curv     44     46   167   255    356
#
# result on hama: icc 12.1.5
# =================
#         bil    bic    nn   con   ycon
# reg2d   1.0    2.2   1.3   242    131
# curv     51     52   175   263    295
#
# result on bailung: gcc 4.9.1
# =================
#         bil    bic    nn   con   ycon
# reg2d   0.8    1.6   1.0   204    135
# curv     33     34   103   210    284
#
# result on bailung: gcc 4.8.2
# =================
#         bil    bic    nn   con   ycon
# reg2d   0.7    1.5   1.0   205    140
# curv     38     39   105   210    296
#
# result on bailung: icc
# =================
#         bil    bic    nn   con   ycon
# reg2d   0.8    1.7   1.1   217    139
# curv     44     45   117   204    280
#
# result on bull: gcc 4.8.2
# =================
#         bil    bic    nn   con   ycon
# reg2d   0.8    1.6   1.0   227    134
# curv     40     40    94   227    275
#
#
# result on bull: icc 14.0.2
# =================
#         bil    bic    nn   con   ycon
# reg2d   0.7    1.6   0.9   297    101
# curv     45     46    96   298    214
#
# result on blizzard: (interactiv)
# ==================
#         bil    bic    nn   con   ycon
# reg2d   2.2    5.4   2.8   470    287
# curv     88     93   481   480    626
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

#!/bin/sh
#
CDO=cdo
#
FORMAT="-f srv -b 32"
#
IFILE=EH5_AMIP_1_TSURF_6H_1991-1995.grb
OFILE=ts_6h_5years
$CDO $FORMAT remapnn,lon=55_lat=10 $IFILE $OFILE
#
IFILE=$OFILE
OFILE=ts_1d_5years
$CDO $FORMAT daymean $IFILE $OFILE
#
IFILE=$OFILE
OFILE=ts_mm_5years
$CDO $FORMAT monmean $IFILE $OFILE
#
STATS="min max sum avg mean std std1 var var1"
IFILE=$OFILE
for STAT in $STATS; do
  $CDO $FORMAT tim$STAT $IFILE tim${STAT}_ref
done

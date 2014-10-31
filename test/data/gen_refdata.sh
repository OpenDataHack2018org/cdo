#!/bin/sh
#
CDO=cdo-1.6.4
#
FORMAT="-f srv -b 32"
#
########################################################################
#
# Test netCDF files
#
IFILES="testfile01"
for FILE in $IFILES; do
  IFILE=netcdf_${FILE}.nc
  OFILE=netcdf_${FILE}_sinfo_ref
  $CDO sinfo $IFILE > $OFILE
  OFILE=netcdf_${FILE}_info_ref
  $CDO info $IFILE > $OFILE
done
########################################################################
#
# Timstat
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
########################################################################
#
# Vertint
#
IFILE=ECHAM5_T21L19monavg.grb
OFILE=hl_l19.grb
$CDO fldmean -sp2gp -selcode,129,130,152 $IFILE $OFILE
IFILE=$OFILE
OFILE=ml2pl_ref
$CDO $FORMAT ml2pl,92500,85000,50000,20000 $IFILE $OFILE
########################################################################
#
# Spectral
#
IFILE=ECHAM5_T21L19monavg.grb
OFILE=t21_geosp_tsurf.grb
$CDO selname,geosp,tsurf $IFILE $OFILE
#
IFILE=$OFILE
OFILE=gp2sp_ref
$CDO gp2sp $IFILE $OFILE
OFILE=gp2spl_ref
$CDO gp2spl $IFILE $OFILE
IFILE=gp2sp_ref
OFILE=sp2gp_ref
$CDO sp2gp $IFILE $OFILE
IFILE=gp2spl_ref
OFILE=sp2gpl_ref
$CDO sp2gpl $IFILE $OFILE
########################################################################
#
# Remap
#
cdo -f grb setrtomiss,0,10000  -gridboxmean,8,8 -topo bathy4.grb
#
GRIDS="n16 n32"
RMODS="bil bic nn con laf"
IFILE=bathy4.grb
for GRID in $GRIDS; do
  for RMOD in $RMODS; do
    OFILE=${GRID}_${RMOD}
    $CDO $FORMAT remap${RMOD},$GRID $IFILE ${OFILE}_ref
  done
done
########################################################################
#
# Select
#
IFILE=ECHAM5_T21L19.grb
OFILE=pl_data.grb
$CDO -ml2pl,90000,9000,900,90 -remapbil,global_30 -sp2gp -seltimestep,1/5 -dayavg -selcode,129,130,152 $IFILE $OFILE
IFILE=$OFILE
OFILE=select1_ref
$CDO select,code=130,152 $IFILE $OFILE
OFILE=select2_ref
$CDO select,code=130,152,level=9000,90 $IFILE $OFILE
OFILE=select3_ref
$CDO select,code=130,152,level=9000,90,timestep=2,3,5 $IFILE $OFILE
OFILE=select4_ref
$CDO select,code=130 $IFILE $OFILE
OFILE=select5_ref
$CDO select,level=90000 $IFILE $OFILE
########################################################################
#
# Detrend
#
OFILE=detrend_data
$CDO $FORMAT -sqr -for,1,100 $OFILE
IFILE=$OFILE
OFILE=detrend_ref
$CDO detrend $IFILE $OFILE
########################################################################

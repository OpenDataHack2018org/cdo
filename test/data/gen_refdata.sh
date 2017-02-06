#!/bin/sh
#
CDO=cdo
#
FORMAT="-f srv -b F32"
#
########################################################################
#
# EOF
#
export CDO_FILE_SUFFIX=NULL
#export CDO_SVD_MODE=danielson_lanczos
export CDO_WEIGHT_MODE=off
#
IFILE=psl_DJF_anom.grb
cdo $FORMAT eof,1 $IFILE eval_ref eof_ref
cdo $FORMAT eofcoeff eof_ref $IFILE pcoeff
exit
########################################################################
#
# Vertstat
#
STATS="min max sum avg mean std std1 var var1 int"
IFILE=pl_data.grb
for STAT in $STATS; do
  $CDO $FORMAT vert$STAT $IFILE vert${STAT}_ref
done
exit
########################################################################
#
# Comparision
#
STATS="eqc nec lec ltc gec gtc"
#
IFILE=comptest.srv
cdo $FORMAT copy t21_geosp_tsurf.grb $IFILE
for STAT in $STATS; do
  $CDO $FORMAT $STAT,300 $IFILE comp_${STAT}_ref
done
exit
########################################################################
#
# Ymonstat
#
STATS="min max sum avg mean std std1 var var1"
#
IFILE=ts_mm_5years
for STAT in $STATS; do
  $CDO $FORMAT ymon$STAT $IFILE ymon${STAT}_ref
done
#
IFILE=ts_mm_5years
for STAT in $STATS; do
  $CDO $FORMAT yseas$STAT $IFILE yseas${STAT}_ref
done
########################################################################
#
# Timstat Yearstat Monstat Daystat Runstat
#
IFILE=$HOME/data/cdt/cera/EH5_AMIP_1_TSURF_6H_1991-1995.grb
OFILE=ts_6h_5years
$CDO $FORMAT remapnn,lon=55_lat=10 $IFILE $OFILE
#
IFILE=$OFILE
OFILE=ts_1d_5years
$CDO $FORMAT daymean $IFILE $OFILE
$CDO selmon,1 -selyear,1991 $IFILE ts_6h_1mon
#
IFILE=$OFILE
OFILE=ts_mm_5years
$CDO $FORMAT -settime,12:00:00 -setday,15 -monmean $IFILE $OFILE
$CDO selyear,1991 $IFILE ts_1d_1year
#
STATS="min max sum avg mean std std1 var var1"
#
IFILE=ts_mm_5years
for STAT in $STATS; do
  $CDO $FORMAT seas${STAT} $IFILE seas${STAT}_ref
done
#
IFILE=ts_mm_5years
for STAT in $STATS; do
  $CDO $FORMAT run${STAT},12 $IFILE run${STAT}_ref
done
#
IFILE=ts_mm_5years
for STAT in $STATS; do
  $CDO $FORMAT tim$STAT $IFILE tim${STAT}_ref
done
#
IFILE=ts_mm_5years
for STAT in $STATS; do
  $CDO $FORMAT year$STAT $IFILE year${STAT}_ref
done
#
IFILE=ts_1d_1year
for STAT in $STATS; do
  $CDO $FORMAT mon$STAT $IFILE mon${STAT}_ref
done
#
IFILE=ts_6h_1mon
for STAT in $STATS; do
  $CDO $FORMAT day$STAT $IFILE day${STAT}_ref
done
exit
########################################################################
#
# Zonstat
#
STATS="min max sum avg mean std std1 var var1"
IFILE=t21_geosp_tsurf.grb
for STAT in $STATS; do
  $CDO $FORMAT zon$STAT $IFILE zon${STAT}_ref
done
exit
########################################################################
#
# Merstat
#
STATS="min max sum avg mean std std1 var var1"
IFILE=t21_geosp_tsurf.grb
for STAT in $STATS; do
  $CDO $FORMAT mer$STAT $IFILE mer${STAT}_ref
done
exit
########################################################################
#
# MapReduce
#
for grid  in r18x9 icon_cell; do
  REFGRID="${grid}_grid"
  $CDO -f nc -temp,$REFGRID data_${grid}.nc
  $CDO -f nc -gtc,273.15 -temp,$REFGRID mask_${grid}.nc
  $CDO -f nc reducegrid,mask_${grid}.nc data_${grid}.nc reduced_${grid}.nc
  $CDO griddes reduced_${grid}.nc > griddes.${grid}
  rm -f data_${grid}.nc mask_${grid}.nc reduced_${grid}.nc
done
exit
########################################################################
#
# Arith
#
IFILE=arith1.srv
MASK=arithmask.srv
OPS="add sub mul div"
for OP in $OPS; do
  $CDO $FORMAT $OP $IFILE $MASK arith${OP}_ref
done
exit
########################################################################
#
# Fldpctl
#
IFILE=t21_geosp_tsurf.grb
PCTLS="1 20 25 33 50 66 75 80 99 100"
for PCTL in $PCTLS; do
  $CDO $FORMAT fldpctl,$PCTL $IFILE fldpctl${PCTL}_ref
done
rm -f $IFILE
exit
########################################################################
#
# Fldstat
#
STATS="min max sum avg mean std std1 var var1"
IFILE=t21_geosp_tsurf.grb
for STAT in $STATS; do
  $CDO $FORMAT fld$STAT $IFILE fld${STAT}_ref
done
exit
########################################################################
#
# Enspctl
#
IFILE=ts_mm_5year
export CDO_FILE_SUFFIX=NULL
$CDO splityear $IFILE ts_year
IFILE="ts_year????"
PCTLS="1 20 25 33 50 66 75 80 99 100"
for PCTL in $PCTLS; do
  $CDO $FORMAT enspctl,$PCTL $IFILE enspctl${PCTL}_ref
done
rm -f $IFILE
exit
########################################################################
#
# Test File
#
$CDO $FORMAT cdiwrite,1,global_10,3,3,3 file_F32_srv_ref
#
########################################################################
#
# Test GRIB files
#
IFILES="testfile01 testfile02"
for FILE in $IFILES; do
  IFILE=grib_${FILE}.grb
  OFILE=grib_${FILE}_sinfon_ref
  $CDO sinfo $IFILE > $OFILE
  OFILE=grib_${FILE}_infon_ref
  $CDO info $IFILE > $OFILE
done
#
########################################################################
#
# Test netCDF files
#
IFILES="testfile01 testfile02"
for FILE in $IFILES; do
  IFILE=netcdf_${FILE}.nc
  OFILE=netcdf_${FILE}_sinfon_ref
  $CDO sinfon $IFILE > $OFILE
  OFILE=netcdf_${FILE}_infon_ref
  $CDO infon $IFILE > $OFILE
done
exit
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
RMODS="bil bic dis nn con laf"
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

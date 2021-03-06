#!/bin/sh
#
CDO=cdo
#
FORMAT="-f srv -b F32"
########################################################################
#
# Ymonarith
#
IFILE1=ts_mm_5years
IFILE2=ts_mm_1year
OPS="add sub mul div"
for OP in $OPS; do
  $CDO $FORMAT ymon$OP $IFILE1 $IFILE2 ymon${OP}_ref
done
exit
########################################################################
#
# Remap regional grid
#
GRID=spain.grid
cat > $GRID <<EOF
gridtype=lonlat
xsize=20
ysize=18
xfirst=-13
yfirst=33
xinc=.8
yinc=.8
EOF
RMODS="bil bic dis nn con con2 ycon laf"
IFILE=tsurf_spain.grb
for RMOD in $RMODS; do
  OFILE=tsurf_spain_${RMOD}
  for extra in def off on; do
      EXTRA="$extra"
      if [ "$EXTRA" = "def" ]; then EXTRA=""; fi
      REMAP_EXTRAPOLATE=$EXTRA $CDO $FORMAT remap${RMOD},${GRID} $IFILE ${OFILE}_${extra}_ref
  done
done
exit
########################################################################
#
# Filter
#
IFILE=ts_mm_5years
#
OFILE=lowpass_ref
$CDO $FORMAT lowpass,1 $IFILE $OFILE
OFILE=highpass_ref
$CDO $FORMAT highpass,1 $IFILE $OFILE
OFILE=bandpass_ref
$CDO $FORMAT bandpass,.5,5 $IFILE $OFILE
#
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
# Timstat2
#
IFILE=ts_mm_5years
IFILE1=ts_mm_2years1
$CDO selyear,1991,1992 $IFILE $IFILE1
IFILE2=ts_mm_2years2
$CDO selyear,1993,1994 $IFILE $IFILE2
#
OFILE=timcor_ref
$CDO $FORMAT timcor $IFILE1 $IFILE2 $OFILE
OFILE=timcovar_ref
$CDO $FORMAT timcovar $IFILE1 $IFILE2 $OFILE
#
rm $IFILE1 $IFILE2
#
exit
########################################################################
#
# Remap regional grid
#
RMODS="bil bic dis nn con con2 ycon laf"
#cdo -f grb topo,europe_5 topo_eu5.grb
IFILE=topo_eu5.grb
$CDO -f grb -topo,europe_5 $IFILE
for RMOD in $RMODS; do
  OFILE=topo_eu5_${RMOD}
  for extra in def off on; do
      EXTRA="$extra"
      if [ "$EXTRA" = "def" ]; then EXTRA=""; fi
      REMAP_EXTRAPOLATE=$EXTRA $CDO $FORMAT remap${RMOD},global_5 $IFILE ${OFILE}_${extra}_ref
  done
done
exit
########################################################################
#
# smooth
#
IFILE=t21_geosp_tsurf_sea.grb
#
OFILE=smooth
$CDO $FORMAT smooth,radius=5deg $IFILE ${OFILE}1_ref
$CDO $FORMAT smooth,radius=5deg,maxpoints=3 $IFILE ${OFILE}2_ref
$CDO $FORMAT smooth,radius=5deg,nsmooth=9 $IFILE ${OFILE}3_ref
exit
########################################################################
#
# Timpctl Yearpctl Monstat Daypctl
#
PCTLS="50"
#
IFILE=ts_mm_5years
for PCTL in $PCTLS; do
  $CDO $FORMAT seaspctl,${PCTL} $IFILE seasmin_ref seasmax_ref seaspctl${PCTL}_ref
done
#
IFILE=ts_mm_5years
for PCTL in $PCTLS; do
  $CDO $FORMAT timpctl,$PCTL $IFILE timmin_ref timmax_ref timpctl${PCTL}_ref
done
#
IFILE=ts_mm_5years
for PCTL in $PCTLS; do
  $CDO $FORMAT yearpctl,$PCTL $IFILE yearmin_ref yearmax_ref yearpctl${PCTL}_ref
done
#
IFILE=ts_1d_1year
for PCTL in $PCTLS; do
  $CDO $FORMAT monpctl,$PCTL $IFILE monmin_ref monmax_ref monpctl${PCTL}_ref
done
#
IFILE=ts_6h_1mon
for PCTL in $PCTLS; do
  $CDO $FORMAT daypctl,$PCTL $IFILE daymin_ref daymax_ref daypctl${PCTL}_ref
done
exit
#
########################################################################
#
# runpctl
#
IFILE=ts_mm_5years
PCTLS="1 20 25 33 50 66 75 80 99 100"
for PCTL in $PCTLS; do
  $CDO $FORMAT runpctl,$PCTL,12 $IFILE runpctl${PCTL}_ref
done
exit
#
########################################################################
#
# setmiss
#
IFILE=t21_geosp_tsurf_sea.grb
#
$CDO $FORMAT setmisstoc,0 $IFILE setmisstoc_ref
$CDO $FORMAT setmisstonn $IFILE setmisstonn_ref
$CDO $FORMAT setmisstodis $IFILE setmisstodis_ref
exit
########################################################################
#
# Remap
#
cdo -f grb setrtomiss,0,10000  -gridboxmean,8,8 -topo bathy4.grb
#
GRIDS="n16 n32"
RMODS="bil bic dis nn con con2 laf"
RMODS="con2"
IFILE=bathy4.grb
for GRID in $GRIDS; do
  for RMOD in $RMODS; do
    OFILE=${GRID}_${RMOD}
    $CDO $FORMAT remap${RMOD},$GRID $IFILE ${OFILE}_ref
  done
done
exit
########################################################################
#
# Test GRIB files
#
IFILES="testfile01 testfile02 testfile03"
for FILE in $IFILES; do
  IFILE=grib_${FILE}.grb
  OFILE=grib_${FILE}_sinfo_ref
  $CDO sinfo $IFILE > $OFILE
  OFILE=grib_${FILE}_info_ref
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
# Mergetime
#
$CDO $FORMAT settaxis,2001-01-01,,1year -const,1,r1x1 mergetime_y1
$CDO $FORMAT settaxis,2002-01-01,,1year -const,2,r1x1 mergetime_y2
$CDO $FORMAT mergetime mergetime_y1 mergetime_y2 mergetime_y12
env SKIP_SAME_TIME=1 cdo mergetime mergetime_y2 mergetime_y12  mergetime_ref
env SKIP_SAME_TIME=1 cdo mergetime mergetime_y12 mergetime_y2  mergetime_ref2
exit
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
$CDO selmonth,1 -selyear,1991 $IFILE ts_6h_1mon
#
IFILE=$OFILE
OFILE=ts_mm_5years
$CDO $FORMAT -settime,12:00:00 -setday,15 -monmean $IFILE $OFILE
$CDO selyear,1991 $IFILE ts_1d_1year
#
STATS="min max sum avg mean std std1 var var1 range"
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
# Ydrunstat
#
STATS="min max sum avg mean std std1 var var1"
#
IFILE=ts_mm_5years
for STAT in $STATS; do
  $CDO $FORMAT ydrun$STAT,8 $IFILE ydrun${STAT}_ref
done
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
# Zonstat
#
STATS="min max sum avg mean std std1 var var1 range"
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
STATS="min max sum avg mean std std1 var var1 range"
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

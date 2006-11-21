/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

/*
      MODULE      OPERATOR     INDEX    DESCRIPTION
      
      EcaCfd      eca_cfd      CFD      maximum number of consecutive frost days
      EcaCsu      eca_csu      CSU      maximum number of consecutive summer days
      EcaCwdi     eca_cwdi     CWDI     cold wave duration index
      EcaCwfi     eca_cwfi     CWFI     number of cold-spell days
      EcaEtr      eca_etr      ETR      intra-period extreme temperature range
      EcaFd       eca_fd       FD       number of frost days
      EcaGsl      eca_gsl      GSL      growing season length
      EcaHd       eca_hd       HD       heating degree days
      EcaHwdi     eca_hwdi     HWDI     heat wave duration index
      EcaHwfi     eca_hwfi     HWFI     number of warm-spell days
      EcaId       eca_id       ID       number of ice days
      EcaSu       eca_su       SU       number of summer days
      EcaTg10p    eca_tg10p    TG10p    percent of time TX < 10th percentile of daily mean temperature
      EcaTg90p    eca_tg90p    TG90p    percent of time TX > 90th percentile of daily mean temperature
      EcaTn10p    eca_tn10p    TN10p    percent of time TX < 10th percentile of daily minimum temperature   
      EcaTn90p    eca_tn90p    TN90p    percent of time TX > 90th percentile of daily minimum temperature
      EcaTr       eca_tr       TR       number of tropical nights
      EcaTx10p    eca_tx10p    TX10p    percent of time TX < 10th percentile of daily maximum temperature
      EcaTx90p    eca_tx90p    TX90p    percent of time TX > 90th percentile of daily maximum temperature

      EcaCdd      eca_cdd      CDD      maximum number of consecutive dry days
      EcaCwd      eca_cwd      CWD      maximum number of consecutive wet days
      EcaR10mm    eca_r10mm    R10mm    number of days with precipitation >= 10 mm
      EcaR20mm    eca_r20mm    R20mm    number of days with precipitation >= 20 mm
      EcaR75p     eca_r75p     R75p     Percent of time RR > 75th percentile of daily precipitation amount 
      EcaR75ptot  eca_r75ptot  R75pTOT  Percentage of annual total precipitation due to events wit RR > 75th percentile of daily precipitation amount
      EcaR90p     eca_r90p     R90p     Percent of time RR > 90th percentile of daily precipitation amount
      EcaR90ptot  eca_r90ptot  R90pTOT  Percentage of annual total precipitation due to events wit RR > 90th percentile of daily precipitation amount
      EcaR95p     eca_r95p     R95p     Percent of time RR > 95th percentile of daily precipitation amount
      EcaR95ptot  eca_r95ptot  R95pTOT  Percentage of annual total precipitation due to events wit RR > 95th percentile of daily precipitation amount
      EcaR99p     eca_r99p     R99p     Percent of time RR > 75th percentile of daily precipitation amount
      EcaR99ptot  eca_r99ptot  R99pTOT  Percentage of annual total precipitation due to events wit RR > 99th percentile of daily precipitation amount
      EcaRr1      eca_rr1      RR1      number of wet days
      EcaSdii     eca_sdii     SDII     simple daily intensity index
      
      EcaFdns     eca_fdns              frost days without surface snow 
      EcaStrwind  eca_strwind           number of strong-wind days 
*/

#include <stdio.h>
#include <string.h>

#include "cdo.h"
#include "cdo_int.h"
#include "ecacore.h"
#include "ecautil.h"


#define TO_DEG_CELSIUS(x) ((x) - 273.15)
#define TO_KELVIN(x) ((x) + 273.15)


static const char CFD_NAME[]         = "consecutive_frost_days";
static const char CFD_LONGNAME[]     = "greatest number of consecutive frost days";

static const char CSU_NAME[]         = "consecutive_summer_days";
static const char CSU_LONGNAME[]     = "greatest number of consecutive summer days";

static const char CWDI_NAME[]        = "cold_wave_duration_index";
static const char CWDI_LONGNAME[]    = "number of days with Tmin more than %1.0f degree Celsius below mean of reference period";
static const char CWDI_NAME2[]       = "to_be_defined";
static const char CWDI_LONGNAME2[]   = "to be defined";

static const char CWFI_NAME[]        = "cold_spell_days";
static const char CWFI_LONGNAME[]    = "number of days with Tmean below 10th percentile of reference period";
static const char CWFI_NAME2[]       = "to_be_defined";
static const char CWFI_LONGNAME2[]   = "to be defined";

static const char ETR_NAME[]         = "intra_period_extreme_temperature_range";
static const char ETR_LONGNAME[]     = "extreme temperature range in observation period";
static const char ETR_UNITS[]        = "K";

static const char FD_NAME[]          = "frost_days";
static const char FD_LONGNAME[]      = "number of days with Tmin below 0 degree Celsius";

static const char GSL_NAME[]         = "growing_season_length";
static const char GSL_LONGNAME[]     = "growing season length";
static const char GSL_NAME2[]        = "to_be_defined";
static const char GSL_LONGNAME2[]    = "to be defined";
static const char GSL_NAME3[]        = "to_be_defined";
static const char GSL_LONGNAME3[]    = "to be defined";

static const char HD_NAME[]          = "heating_degree_days";
static const char HD_LONGNAME[]      = "Heating degree days";

static const char HWDI_NAME[]        = "heat_wave_duration_index";
static const char HWDI_LONGNAME[]    = "number of days with Tmax more than %1.0f degree Celsius above mean of reference period";
static const char HWDI_NAME2[]       = "to_be_defined";
static const char HWDI_LONGNAME2[]   = "to be defined";

static const char HWFI_NAME[]        = "warm_spell_days";
static const char HWFI_LONGNAME[]    = "number of days with Tmean above 90th percentile of reference period";
static const char HWFI_NAME2[]       = "to_be_defined";
static const char HWFI_LONGNAME2[]   = "to be defined";

static const char ID_NAME[]          = "ice_days";
static const char ID_LONGNAME[]      = "number of days with Tmax below 0 degree Celsius";

static const char SU_NAME[]          = "summer_days";
static const char SU_LONGNAME[]      = "number of days with Tmax above %1.0f degree Celsius";

static const char TG10P_NAME[]       = "cold_days_wrt_10th_percentile_of_reference_period";
static const char TG10P_LONGNAME[]   = "percentage of time with Tmean below 10th percentile of reference period";

static const char TG90P_NAME[]       = "warm_days_wrt_90th_percentile_of_reference_period";
static const char TG90P_LONGNAME[]   = "percentage of time with Tmean above 90th percentile of reference period";

static const char TN10P_NAME[]       = "cold_nights_wrt_10th_percentile_of_reference_period";
static const char TN10P_LONGNAME[]   = "percentage of time with Tmin below 10th percentile of reference period";

static const char TN90P_NAME[]       = "warm_nights_wrt_90th_percentile_of_reference_period";
static const char TN90P_LONGNAME[]   = "percentage of time with Tmin above 90th percentile of reference period";

static const char TR_NAME[]          = "tropical_nights";
static const char TR_LONGNAME[]      = "number of days with Tmin above %1.0f degree Celsius";

static const char TX10P_NAME[]       = "cold_days_wrt_10th_percentile_of_reference_period";
static const char TX10P_LONGNAME[]   = "percentage of time with Tmax below 10th percentile of reference period";

static const char TX90P_NAME[]       = "warm_days_wrt_90th_percentile_of_reference_period";
static const char TX90P_LONGNAME[]   = "percentage of time with Tmax above 90th percentile of reference period";

static const char CDD_NAME[]         = "consecutive_dry_days";
static const char CDD_LONGNAME[]     = "greatest number of consecutive days with daily precipitation below 1 mm";
static const char CDD_NAME2[]        = "to_be_defined";
static const char CDD_LONGNAME2[]    = "to be defined";

static const char CWD_NAME[]         = "consecutive_wet_days";
static const char CWD_LONGNAME[]     = "greatest number of consecutive days with daily precipitation above 1 mm";
static const char CWD_NAME2[]        = "to_be_defined";
static const char CWD_LONGNAME2[]    = "to be defined";

static const char R10MM_NAME[]       = "heavy_precipitation_days";
static const char R10MM_LONGNAME[]   = "number of days with daily precipitation exceeding 10 mm";

static const char R20MM_NAME[]       = "very_heavy_precipitation_days";
static const char R20MM_LONGNAME[]   = "number of days with daily precipitation exceeding 20 mm";

static const char R75P_NAME[]        = "moderate_wet_days_wrt_75th_percentile_of_reference_period";
static const char R75P_LONGNAME[]    = "percentage of time with precipitation sum above 75th percentile of reference period";

static const char R75PTOT_NAME[]     = "precipitation_fraction_due_to_R75p_days";
static const char R75PTOT_LONGNAME[] = "percentage of total precipitation due to moderate wet days";

static const char R90P_NAME[]        = "very_wet_days_wrt_90th_percentile_of_reference_period";
static const char R90P_LONGNAME[]    = "percentage of time with precipitation sum above 90th percentile of reference period";

static const char R90PTOT_NAME[]     = "precipitation_fraction_due_to_R90p_days";
static const char R90PTOT_LONGNAME[] = "percentage of total precipitation due to very wet days";

static const char R95P_NAME[]        = "very_wet_days_wrt_95th_percentile_of_reference_period";
static const char R95P_LONGNAME[]    = "percentage of time with precipitation sum above 95th percentile of reference period";

static const char R95PTOT_NAME[]     =  "precipitation_fraction_due_to_R95p_days";
static const char R95PTOT_LONGNAME[] =  "percentage of total precipitation due to very wet days";

static const char R99P_NAME[]        = "extremely_wet_days_wrt_99th_percentile_of_reference_period";
static const char R99P_LONGNAME[]    = "percentage of time with precipitation sum above 99th percentile of reference period";

static const char R99PTOT_NAME[]     = "precipitation_fraction_due_to_R99p_days";
static const char R99PTOT_LONGNAME[] = "percentage of total precipitation due to extremely wet days";

static const char RR1_NAME[]         = "wet_days";
static const char RR1_LONGNAME[]     = "number of days with daily precipitation exceeding 1 mm";

static const char RX1DAY_NAME[]      = "highest_one_day_precipitation_amount";
static const char RX1DAY_LONGNAME[]  = "highest precipitation sum for one day interval";
static const char RX1DAY_UNITS[]     = "kg/m2";

static const char RX5DAY_NAME[]      = "highest_five_day_precipitation_amount";
static const char RX5DAY_LONGNAME[]  = "highest precipitation sum for five day interval";
static const char RX5DAY_UNITS[]     = "kg/m2";
static const char RX5DAY_NAME2[]     = "to_be_defined";
static const char RX5DAY_LONGNAME2[] = "to be defined";

static const char SDII_NAME[]        = "simple_daily_intensity_index";
static const char SDII_LONGNAME[]    = "mean precipitation amount on wet days";
static const char SDII_UNITS[]       = "kg/m2";

static const char FDNS_NAME[]        = "frost_days_where_no_snow";       
static const char FDNS_LONGNAME[]    = "number of days with Tmin below 0 degree Celsius and no snowcover";       

static const char STRWIND_NAME[]       = "strong_breeze_days";
static const char STRWIND_LONGNAME[]   = "number of days with maximum wind speed above 10.5 m/s";
static const char STRWIND_NAME2[]      = "to_be_defined";
static const char STRWIND_LONGNAME2[]  = "to be defined";
static const char STRWIND2_NAME[]      = "strong_gale_days";
static const char STRWIND2_LONGNAME[]  = "number of days with maximum wind speed above 20.5 m/s";
static const char STRWIND2_NAME2[]     = "to_be_defined";
static const char STRWIND2_LONGNAME2[] = "to be defined";
static const char STRWIND3_NAME[]      = "hurricane_days";
static const char STRWIND3_LONGNAME[]  = "number of days with maximum wind speed above 32.5 m/s";
static const char STRWIND3_NAME2[]     = "to_be_defined";
static const char STRWIND3_LONGNAME2[] = "to be defined";


/* ECA temperature indices */


void *EcaCfd(void *argument)
{
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_cfd", 0, 17, NULL);
  
  request.var1.name     = CFD_NAME;
  request.var1.longname = CFD_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = farselltc;
  request.var1.f1arg    = TO_KELVIN(0.0);
  request.var1.f2       = farnum2;
  request.var1.f3       = farmax;
  request.var1.mulc     = 0.0;
  request.var1.addc     = 0.0;
  request.var1.epilog   = NONE;
  request.var2.h2       = NULL;
  request.var2.h3       = NULL;
   
  eca1(&request);
  cdoFinish();
  
  return (0);
}


void *EcaCsu(void *argument)
{
  double argT = 25.0;
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_csu", 0, 17, NULL);

  if ( operatorArgc() > 0 ) argT = atof(operatorArgv()[0]);

  request.var1.name     = CSU_NAME;
  request.var1.longname = CSU_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = farselgtc;
  request.var1.f1arg    = TO_KELVIN(argT);
  request.var1.f2       = farnum2;
  request.var1.f3       = farmax;
  request.var1.mulc     = 0.0;
  request.var1.addc     = 0.0;
  request.var1.epilog   = NONE;
  request.var2.h2       = NULL;
  request.var2.h3       = NULL;
  
  eca1(&request);
  cdoFinish();
  
  return (0);
}


void *EcaCwdi(void *argument)
{
  static const char func[] = "EcaCwdi";
  char *longname;
  int argN = 6;
  double argT = 5.0;
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_cwdi", 0, 17, NULL);

  if ( operatorArgc() > 0 ) argN = atoi(operatorArgv()[0]);
  if ( operatorArgc() > 1 ) argT = atof(operatorArgv()[1]);
  
  longname = (char *) malloc(strlen(CWDI_LONGNAME) + 40);
  sprintf(longname, CWDI_LONGNAME, argT);

  request.var1.name     = CWDI_NAME;
  request.var1.longname = longname;
  request.var1.units    = NULL;
  request.var1.f1       = NULL;
  request.var1.f2       = farcsub;
  request.var1.f2arg    = argT;
  request.var1.f3       = farsellt;
  request.var1.f4       = farnum2;
  request.var1.f5       = farnum3;
  request.var1.f5arg    = argN;
  request.var1.epilog   = NONE;
  request.var2.name     = CWDI_NAME2;
  request.var2.longname = CWDI_LONGNAME2;
  request.var2.h1       = farseleqc;
  request.var2.h1arg    = argN;
  request.var2.h2       = farnum;
   
  eca2(&request);
  
  free(longname);
  cdoFinish();
  
  return (0);
}


void *EcaCwfi(void *argument)
{
  int argN = 6;
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_cwfi", 0, 17, NULL);

  if ( operatorArgc() > 0 ) argN = atoi(operatorArgv()[0]);

  request.var1.name     = CWFI_NAME;
  request.var1.longname = CWFI_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = NULL;
  request.var1.f2       = NULL;
  request.var1.f3       = farsellt;
  request.var1.f4       = farnum2;
  request.var1.f5       = farnum3;
  request.var1.f5arg    = argN;
  request.var1.epilog   = NONE;
  request.var2.name     = CWFI_NAME2;
  request.var2.longname = CWFI_LONGNAME2;
  request.var2.h1       = farseleqc;
  request.var2.h1arg    = argN;
  request.var2.h2       = farnum;
   
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaEtr(void *argument)
{
  ECA_REQUEST_3 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_etr", 0, 17, NULL);
  
  request.name     = ETR_NAME;
  request.longname = ETR_LONGNAME;
  request.units    = NULL;
  request.f1       = farmax; 
  request.f2       = farmin;
  request.f3       = farsub;
   
  eca3(&request);
  cdoFinish();
  
  return (0);
}


void *EcaFd(void *argument)
{
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_fd", 0, 17, NULL);
  
  request.var1.name     = FD_NAME;
  request.var1.longname = FD_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = farselltc; 
  request.var1.f1arg    = TO_KELVIN(0.0);
  request.var1.f2       = farnum;
  request.var1.f3       = NULL;
  request.var1.mulc     = 0.0;    
  request.var1.addc     = 0.0;
  request.var1.epilog   = NONE;    
  request.var2.h2       = NULL; 
  request.var2.h3       = NULL; 
   
  eca1(&request);
  cdoFinish();
  
  return (0);
}


void *EcaGsl(void *argument)
{
  int argN = 6; 
  double argT = 5.0;
  ECA_REQUEST_4 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_gsl", 0, 8, NULL);
  
  if ( operatorArgc() > 0 ) argN = atoi(operatorArgv()[0]);
  if ( operatorArgc() > 1 ) argT = atof(operatorArgv()[1]);
  
  request.name      = GSL_NAME;
  request.longname  = GSL_LONGNAME;
  request.units     = NULL;
  request.name2     = GSL_NAME2;
  request.longname2 = GSL_LONGNAME2;
  request.units2    = NULL;
  request.name3     = GSL_NAME3;
  request.longname3 = GSL_LONGNAME3;
  request.units3    = NULL;
  request.s1        = farselgtc; 
  request.s1arg     = TO_KELVIN(argT);
  request.s2        = farselltc;
  request.s2arg     = TO_KELVIN(argT);
  request.consecutiveDays = argN;    
   
  eca4(&request);
  cdoFinish();
  
  return (0);
}


void *EcaHd(void *argument)
{
  double argT1 = 17.0;
  double argT2 = 17.0;
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_hd", 0, 17, NULL);

  if ( operatorArgc() > 0 ) 
    {
      argT1 = atof(operatorArgv()[0]);
      argT2 = argT1;
    }
  if ( operatorArgc() > 1 ) 
    argT2 = atof(operatorArgv()[1]);
  
  request.var1.name     = HD_NAME;
  request.var1.longname = HD_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = farselltc; 
  request.var1.f1arg    = TO_KELVIN(argT2);
  request.var1.f2       = farsum;
  request.var1.f3       = NULL;
  request.var1.mulc     = -1.0;    
  request.var1.addc     = TO_KELVIN(argT1);
  request.var1.epilog   = NONE;    
  request.var2.h2       = NULL; 
  request.var2.h3       = NULL; 
   
  eca1(&request);
  cdoFinish();
  
  return (0);
}


void *EcaHwdi(void *argument)
{
  static const char func[] = "EcaHwdi";
  char *longname;
  int argN = 6;
  double argT = 5.0;
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_hwdi", 0, 17, NULL);

  if ( operatorArgc() > 0 ) argN = atoi(operatorArgv()[0]);
  if ( operatorArgc() > 1 ) argT = atof(operatorArgv()[1]);
  
  longname = (char *) malloc(strlen(HWDI_LONGNAME) + 40);
  sprintf(longname, HWDI_LONGNAME, argT);
  
  request.var1.name     = HWDI_NAME;
  request.var1.longname = longname;
  request.var1.units    = NULL;
  request.var1.f1       = NULL;
  request.var1.f2       = farcadd;
  request.var1.f2arg    = argT;
  request.var1.f3       = farselgt;
  request.var1.f4       = farnum2;
  request.var1.f5       = farnum3;
  request.var1.f5arg    = argN;
  request.var1.epilog   = NONE;
  request.var2.name     = HWDI_NAME2;
  request.var2.longname = HWDI_LONGNAME2;
  request.var2.h1       = farseleqc;
  request.var2.h1arg    = argN;
  request.var2.h2       = farnum;
   
  eca2(&request);
  
  free(longname);
  cdoFinish();
  
  return (0);
}


void *EcaHwfi(void *argument)
{
  int argN = 6;
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_hwfi", 0, 17, NULL);

  if ( operatorArgc() > 0 ) argN = atoi(operatorArgv()[0]);

  request.var1.name     = HWFI_NAME;
  request.var1.longname = HWFI_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = NULL;
  request.var1.f2       = NULL;
  request.var1.f3       = farselgt;
  request.var1.f4       = farnum2;
  request.var1.f5       = farnum3;
  request.var1.f5arg    = argN;
  request.var1.epilog   = NONE;
  request.var2.name     = HWFI_NAME2;
  request.var2.longname = HWFI_LONGNAME2;
  request.var2.h1       = farseleqc;
  request.var2.h1arg    = argN;
  request.var2.h2       = farnum;
   
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaId(void *argument)
{
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_id", 0, 17, NULL);

  request.var1.name     = ID_NAME;
  request.var1.longname = ID_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = farselltc;
  request.var1.f1arg    = TO_KELVIN(0.0);
  request.var1.f2       = farnum;
  request.var1.f3       = NULL;
  request.var1.mulc     = 0.0;    
  request.var1.addc     = 0.0;   
  request.var1.epilog   = NONE; 
  request.var2.h2       = NULL; 
  request.var2.h3       = NULL; 
    
  eca1(&request);
  cdoFinish();
  
  return (0);
}


void *EcaSu(void *argument)
{
  static const char func[] = "EcaSu";
  char *longname;
  double argT = 25.0;
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_su", 0, 17, NULL);

  if ( operatorArgc() > 0 ) argT = atof(operatorArgv()[0]);
  longname = (char *) malloc(strlen(SU_LONGNAME) + 40);
  sprintf(longname, SU_LONGNAME, argT);

  request.var1.name     = SU_NAME;
  request.var1.longname = longname;
  request.var1.units    = NULL;
  request.var1.f1       = farselgtc;
  request.var1.f1arg    = TO_KELVIN(argT);
  request.var1.f2       = farnum;
  request.var1.f3       = NULL;
  request.var1.mulc     = 0.0;    
  request.var1.addc     = 0.0;    
  request.var1.epilog   = NONE; 
  request.var2.h2       = NULL; 
  request.var2.h3       = NULL; 
 
  eca1(&request);
  
  free(longname);
  cdoFinish();
  
  return (0);
}


void *EcaTg10p(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_tg10p", 0, 17, NULL);

  request.var1.name     = TG10P_NAME;
  request.var1.longname = TG10P_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = NULL;
  request.var1.f2       = NULL;
  request.var1.f3       = farsellt;
  request.var1.f4       = farnum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TIME;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaTg90p(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_tg90p", 0, 17, NULL);

  request.var1.name     = TG90P_NAME;
  request.var1.longname = TG90P_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = NULL;
  request.var1.f2       = NULL;
  request.var1.f3       = farselgt;
  request.var1.f4       = farnum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TIME;
  request.var2.h2       = NULL;
  
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaTn10p(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_tn10p", 0, 17, NULL);

  request.var1.name     = TN10P_NAME;
  request.var1.longname = TN10P_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = NULL;
  request.var1.f2       = NULL;
  request.var1.f3       = farsellt;
  request.var1.f4       = farnum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TIME;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaTn90p(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_tn90p", 0, 17, NULL);

  request.var1.name     = TN90P_NAME;
  request.var1.longname = TN90P_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = NULL;
  request.var1.f2       = NULL;
  request.var1.f3       = farselgt;
  request.var1.f4       = farnum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TIME;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaTr(void *argument)
{
  static const char func[] = "EcaTr";
  char *longname;
  double argT = TO_KELVIN(20.0);
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_tr", 0, 17, NULL);

  if ( operatorArgc() > 0 ) argT = atof(operatorArgv()[0]);
  longname = (char *) malloc(strlen(TR_LONGNAME) + 40);
  sprintf(longname, TR_LONGNAME, argT);
 
  request.var1.name     = TR_NAME;
  request.var1.longname = longname;
  request.var1.units    = NULL;
  request.var1.f1       = farselgtc;
  request.var1.f1arg    = TO_KELVIN(argT);
  request.var1.f2       = farnum;
  request.var1.f3       = NULL;
  request.var1.mulc     = 0.0;    
  request.var1.addc     = 0.0;    
  request.var1.epilog   = NONE; 
  request.var2.h2       = NULL; 
  request.var2.h3       = NULL; 
   
  eca1(&request);
  
  free(longname);
  cdoFinish();
  
  return (0);
}


void *EcaTx10p(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_tx10p", 0, 17, NULL);

  request.var1.name     = TX10P_NAME;
  request.var1.longname = TX10P_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = NULL;
  request.var1.f2       = NULL;
  request.var1.f3       = farsellt;
  request.var1.f4       = farnum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TIME;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaTx90p(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_tx90p", 0, 17, NULL);
 
  request.var1.name     = TX90P_NAME;
  request.var1.longname = TX90P_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = NULL;
  request.var1.f2       = NULL;
  request.var1.f3       = farselgt;
  request.var1.f4       = farnum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TIME;
  request.var2.h2       = NULL;
   
  eca2(&request);
  cdoFinish();
  
  return (0);
}


/* ECA precipitation indices */


void *EcaCdd(void *argument)
{
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_cdd", 0, 17, NULL);

  request.var1.name     = CDD_NAME;
  request.var1.longname = CDD_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = farselltc;
  request.var1.f1arg    = 1.0;
  request.var1.f2       = farnum2;
  request.var1.f3       = farmax;
  request.var1.mulc     = 0.0;
  request.var1.addc     = 0.0;
  request.var1.epilog   = NONE;
  request.var2.name     = CDD_NAME2;
  request.var2.longname = CDD_LONGNAME2;
  request.var2.h1       = farseleqc;
  request.var2.h1arg    = 6;
  request.var2.h2       = NULL;
  request.var2.h3       = farnum;
   
  eca1(&request);
  cdoFinish();
  
  return (0);
}


void *EcaCwd(void *argument)
{
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_cwd", 0, 17, NULL);

  request.var1.name     = CWD_NAME;
  request.var1.longname = CWD_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = 1.0;
  request.var1.f2       = farnum2;
  request.var1.f3       = farmax;
  request.var1.mulc     = 0.0;
  request.var1.addc     = 0.0;
  request.var1.epilog   = NONE;
  request.var2.name     = CWD_NAME2;
  request.var2.longname = CWD_LONGNAME2;
  request.var2.h1       = farseleqc;
  request.var2.h1arg    = 6;
  request.var2.h2       = NULL;
  request.var2.h3       = farnum;
   
  eca1(&request);
  cdoFinish();
  
  return (0);
}


void *EcaR10mm(void *argument)
{
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_r10mm", 0, 17, NULL);
  
  request.var1.name     = R10MM_NAME;
  request.var1.longname = R10MM_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = 10.0;
  request.var1.f2       = farnum;
  request.var1.f3       = NULL;
  request.var1.mulc     = 0.0;   
  request.var1.addc     = 0.0;    
  request.var1.epilog   = NONE;
  request.var2.h2       = NULL;
  request.var2.h3       = NULL;
  
  eca1(&request);
  cdoFinish();
  
  return (0);
}


void *EcaR20mm(void *argument)
{
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_r20mm", 0, 17, NULL);

  request.var1.name     = R20MM_NAME;
  request.var1.longname = R20MM_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = 20.0;
  request.var1.f2       = farnum;
  request.var1.f3       = NULL;
  request.var1.mulc     = 0.0;   
  request.var1.addc     = 0.0;    
  request.var1.epilog   = NONE;
  request.var2.h2       = NULL;
  request.var2.h3       = NULL;
    
  eca1(&request);
  cdoFinish();
  
  return (0);
}


void *EcaR75p(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_r75p", 0, 17, NULL);

  request.var1.name     = R75P_NAME;
  request.var1.longname = R75P_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = 1.0;
  request.var1.f2       = NULL;
  request.var1.f3       = farselgt;
  request.var1.f4       = farnum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TIME;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaR75ptot(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_r75ptot", 0, 17, NULL);

  request.var1.name     = R75PTOT_NAME;
  request.var1.longname = R75PTOT_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = 1.0;
  request.var1.f2       = NULL;
  request.var1.f3       = farselgt;
  request.var1.f4       = farsum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TOTAL_AMOUNT;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaR90p(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_r90p", 0, 17, NULL);

  request.var1.name     = R90P_NAME;
  request.var1.longname = R90P_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = 1.0;
  request.var1.f2       = NULL;
  request.var1.f3       = farselgt;
  request.var1.f4       = farnum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TIME;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaR90ptot(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_r90ptot", 0, 17, NULL);

  request.var1.name     = R90PTOT_NAME;
  request.var1.longname = R90PTOT_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = 1.0;
  request.var1.f2       = NULL;
  request.var1.f3       = farselgt;
  request.var1.f4       = farsum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TOTAL_AMOUNT;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaR95p(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_r95p", 0, 17, NULL);

  request.var1.name     = R95P_NAME;
  request.var1.longname = R95P_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = 1.0;
  request.var1.f2       = NULL;
  request.var1.f3       = farselgt;
  request.var1.f4       = farnum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TIME;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaR95ptot(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_r95ptot", 0, 17, NULL);

  request.var1.name     = R95PTOT_NAME;
  request.var1.longname = R95PTOT_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = 1.0;
  request.var1.f2       = NULL;
  request.var1.f3       = farselgt;
  request.var1.f4       = farsum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TOTAL_AMOUNT;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaR99p(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_r99p", 0, 17, NULL);

  request.var1.name     = R99P_NAME;
  request.var1.longname = R99P_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = 1.0;
  request.var1.f2       = NULL;
  request.var1.f3       = farselgt;
  request.var1.f4       = farnum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TIME;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaR99ptot(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_r99ptot", 0, 17, NULL);

  request.var1.name     = R99PTOT_NAME;
  request.var1.longname = R99PTOT_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = 1.0;
  request.var1.f2       = NULL;
  request.var1.f3       = farselgt;
  request.var1.f4       = farsum;
  request.var1.f5       = NULL;
  request.var1.epilog   = PERCENT_OF_TOTAL_AMOUNT;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaRr1(void *argument)
{
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_rr1", 0, 17, NULL);

  request.var1.name     = RR1_NAME;
  request.var1.longname = RR1_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = 1.0;
  request.var1.f2       = farnum;
  request.var1.f3       = NULL;
  request.var1.mulc     = 0.0;    
  request.var1.addc     = 0.0; 
  request.var1.epilog   = NONE;   
  request.var2.h2       = NULL;    
  request.var2.h3       = NULL;    
   
  eca1(&request);
  cdoFinish();
  
  return (0);
}


void *EcaRx1day(void *argument)
{
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  if ( operatorArgc() > 0 && 'm' == operatorArgv()[0][0] )
    cdoOperatorAdd("eca_rx1day", 0, 6,  NULL); /* monthly mode */
  else 
    cdoOperatorAdd("eca_rx1day", 0, 17, NULL);

  request.var1.name     = RX1DAY_NAME;
  request.var1.longname = RX1DAY_LONGNAME;
  request.var1.units    = RX1DAY_UNITS;
  request.var1.f1       = NULL;
  request.var1.f2       = farmax;
  request.var1.f3       = NULL;
  request.var1.mulc     = 0.0;    
  request.var1.addc     = 0.0;    
  request.var1.epilog   = NONE;
  request.var2.h2       = NULL;
  request.var2.h3       = NULL;
   
  eca1(&request);
  cdoFinish();
  
  return (0);
}


void *EcaRx5day(void *argument)
{
  double argX = 50.0;
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  if ( operatorArgc() > 0 ) argX = atof(operatorArgv()[0]);
  cdoOperatorAdd("eca_rx5day", 0, 17, NULL);

  request.var1.name     = RX5DAY_NAME;
  request.var1.longname = RX5DAY_LONGNAME;
  request.var1.units    = RX5DAY_UNITS;
  request.var1.f1       = NULL;
  request.var1.f2       = farmax;
  request.var1.f3       = NULL;
  request.var1.mulc     = 0.0; 
  request.var1.addc     = 0.0; 
  request.var1.epilog   = NONE;
  request.var2.name     = RX5DAY_NAME2;
  request.var2.longname = RX5DAY_LONGNAME2;
  request.var2.h1       = farselgec;
  request.var2.h1arg    = argX;
  request.var2.h2       = farnum;
  request.var2.h3       = NULL;
   
  eca1(&request);
  cdoFinish();
  
  return (0);
}


void *EcaSdii(void *argument)
{
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_sdii", 0, 17, NULL);

  request.var1.name     = SDII_NAME;
  request.var1.longname = SDII_LONGNAME;
  request.var1.units    = SDII_UNITS;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = 1.0;
  request.var1.f2       = farsum;
  request.var1.f3       = NULL;
  request.var1.mulc     = 0.0;    
  request.var1.addc     = 0.0;    
  request.var1.epilog   = MEAN;
  request.var2.h2       = NULL;
  request.var2.h3       = NULL;
   
  eca1(&request);
  cdoFinish();
  
  return (0);
}


void *EcaFdns(void *argument)
{
  ECA_REQUEST_2 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_fdns", 0, 17, NULL);

  request.var1.name     = FDNS_NAME;
  request.var1.longname = FDNS_LONGNAME;
  request.var1.units    = NULL;
  request.var1.f1       = farsellec;
  request.var1.f1arg    = TO_KELVIN(0.0);
  request.var1.f2       = farsellec;
  request.var1.f2arg    = 0.01;
  request.var1.f3       = faradd; /* any f with f(a, b) = miss, if a = miss or b = miss will do here */
  request.var1.f4       = farnum;
  request.var1.f5       = NULL;
  request.var1.epilog   = NONE;
  request.var2.h2       = NULL;
    
  eca2(&request);
  cdoFinish();
  
  return (0);
}


void *EcaStrwind(void *argument)
{
  const char *name, *longname, *name2, *longname2;
  double maxWind;
  int beaufort = 0;
  ECA_REQUEST_1 request;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eca_strwind", 0, 17, NULL);

  name      = STRWIND_NAME;
  longname  = STRWIND_LONGNAME;
  name2     = STRWIND_NAME2;
  longname2 = STRWIND_LONGNAME2;
  maxWind   = 10.5; /* strong breeze */
  
  if ( operatorArgc() > 0 )
    beaufort = atoi(operatorArgv()[0]);
    
  switch (beaufort)
    {
      case 9:  /* strong gale */
        name      = STRWIND2_NAME;
        longname  = STRWIND2_LONGNAME;
        name2     = STRWIND2_NAME2;
        longname2 = STRWIND3_LONGNAME2;
        maxWind   = 20.5;
	break;
	
      case 12: /* hurricane */
        name      = STRWIND3_NAME;
        longname  = STRWIND3_LONGNAME;
        name2     = STRWIND3_NAME2;
        longname2 = STRWIND3_LONGNAME2;
        maxWind   = 32.5;
        break;
        	
      default: /* strong breeze */
        name      = STRWIND_NAME;
        longname  = STRWIND_LONGNAME;
        name2     = STRWIND_NAME2;
        longname2 = STRWIND_LONGNAME2;
        maxWind   = 10.5;
	break;
    }
      
  request.var1.name     = name;
  request.var1.longname = longname;
  request.var1.units    = NULL;
  request.var1.f1       = farselgec;
  request.var1.f1arg    = maxWind;
  request.var1.f2       = farnum;
  request.var1.f3       = NULL;
  request.var1.mulc     = 0.0;    
  request.var1.addc     = 0.0;    
  request.var1.epilog   = NONE;
  request.var2.name     = name2;
  request.var2.longname = longname2;
  request.var2.h1       = farselgec;
  request.var2.h1arg    = maxWind;
  request.var2.h2       = farnum2;
  request.var2.h3       = farmax;
   
  eca1(&request);
  cdoFinish();
  
  return (0);
}

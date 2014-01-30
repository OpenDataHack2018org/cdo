/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2007-2012 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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
   This module contains the following operators:

      Derivepar     geopotheight          geopotential height
*/

#include <ctype.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "vinterp.h"
#include "stdnametable.h"

#define  C_RKBOL         (1.380658e-23)     /* Boltzmann constant in J/K   */
#define  C_RNAVO         (6.0221367e+23)    /* Avogadro constant in 1/mol  */
#define  C_RMD           (28.9644)          /* molecular weight of dry air */
#define  C_RMV           (18.0153)          /* molecular weight of water vapor */
#define  C_R             (C_RKBOL * C_RNAVO)
#define  C_RV            (1000. * C_R / C_RMV)

#define  C_EARTH_GRAV    (9.80665)
#define  C_RKBOL         (1.380658e-23)     /* Boltzmann constant in J/K   */
#define  C_RNAVO         (6.0221367e+23)    /* Avogadro constant in 1/mol  */
#define  C_RMD           (28.9644)          /* molecular weight of dry air */
#define  C_R             (C_RKBOL * C_RNAVO)
#define  C_EARTH_RD      (1000. * C_R / C_RMD)

static double Grav          = C_EARTH_GRAV;
static double RD            = C_EARTH_RD;

static
void MakeGeopotHeight(double *geop, double* gt, double *gq, double *ph, int nhor, int nlev)
{
  int i, j;
  double vtmp;
  double zrg;
  double z2log2;
  double *geopl, *gtl, *gql, *phl;

  z2log2 = 2.0 * log(2.0);
  vtmp   = (C_RV / RD) - 1.0;
  zrg    = 1.0 / Grav;

  if ( gq ) /* Humidity is present */
    {
      for ( j = nlev ; j > 1 ; j-- )
        {
          geopl = geop + nhor*(j-1);
          gtl   = gt   + nhor*(j-1);
          gql   = gq   + nhor*(j-1);
          phl   = ph   + nhor*(j-1);
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(_OPENMP)
#pragma omp parallel for
#endif
          for ( i = 0; i < nhor; i++ )
            geopl[i] = geopl[i+nhor] + RD * gtl[i] * (1.0 + vtmp * gql[i])
                     * log(phl[i+nhor] / phl[i]);
        }

#if defined(SX)
#pragma vdir nodep
#endif
#if defined(_OPENMP)
#pragma omp parallel for
#endif
      for ( i = 0; i < nhor; i++ )
        geop[i] = geop[i+nhor] + RD * gt[i] * (1.0 + vtmp * gq[i]) * z2log2;
    }
  else    /* No humidity */
    {
      for ( j = nlev ; j > 1 ; j-- )
#if defined(SX)
#pragma vdir nodep
#endif
        for ( i = nhor * (j-1) ; i < nhor * j ; i++ )
          geop[i] = geop[i+nhor] + RD * gt[i] * log(ph[i+nhor] / ph[i]);

#if defined(SX)
#pragma vdir nodep
#endif
      for ( i = 0; i < nhor; i++ )
        geop[i] = geop[i+nhor] + RD * gt[i] * z2log2;
    }

#if defined(SX)
#pragma vdir nodep
#endif
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for ( i = 0; i < nhor * (nlev+1); i++ ) geop[i] *= zrg;
}


void minmaxval(long nvals, double *array, int *imiss, double *minval, double *maxval)
{
  long i;
  double xmin =  DBL_MAX;
  double xmax = -DBL_MAX;

  if ( imiss )
    {
      for ( i = 0; i < nvals; i++ )
	{
	  if ( ! imiss[i] )
	    {
	      if      ( array[i] > xmax ) xmax = array[i];
	      else if ( array[i] < xmin ) xmin = array[i];
	    }
	}
    }
  else
    {
      for ( i = 0; i < nvals; i++ )
	{
	  if      ( array[i] > xmax ) xmax = array[i];
	  else if ( array[i] < xmin ) xmin = array[i];
	}
    }

  *minval = xmin;
  *maxval = xmax;
}


void *Derivepar(void *argument)
{
  int GEOPOTHEIGHT, SEALEVELPRESSURE;
  int operatorID;
  int mode;
  enum {ECHAM_MODE, WMO_MODE};
  int geop_code = 0, temp_code = 0, ps_code = 0, lsp_code = 0, hum_code = 0;
  int streamID1, streamID2;
  int vlistID1, vlistID2;
  int gridsize, ngp = 0;
  int recID, nrecs;
  int i, offset;
  int tsID, varID, levelID;
  int nvars;
  int zaxisIDh = -1, nzaxis;
  int ngrids, gridID = -1, zaxisID;
  int nlevel;
  int nvct;
  int geopID = -1, tempID = -1, humID = -1, psID = -1, lnpsID = -1, presID = -1;
  // int clwcID = -1, ciwcID = -1;
  int code, param;
  char paramstr[32];
  char varname[CDI_MAX_NAME], stdname[CDI_MAX_NAME];
  double *single2;
  int taxisID1, taxisID2;
  int lhavevct;
  int nhlevf = 0;
  double *vct = NULL;
  double *geop = NULL, *ps = NULL, *temp = NULL, *hum = NULL;
  // double *lwater = NULL, *iwater = NULL;
  double *geopotheight = NULL;
  int nmiss, nmissout = 0;
  int ltq = FALSE;
  double *array = NULL;
  double *half_press = NULL;
  double minval, maxval;
  int instNum, tableNum;
  int useTable;

  cdoInitialize(argument);

  GEOPOTHEIGHT     = cdoOperatorAdd("geopotheight",   0, 0, NULL);
  SEALEVELPRESSURE = cdoOperatorAdd("sealevelpressure",   0, 0, NULL);

  operatorID = cdoOperatorID();

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);

  ngrids  = vlistNgrids(vlistID1);
  for ( i = 0; i < ngrids; i++ )
    {
      gridID = vlistGrid(vlistID1, i);
      if ( gridInqType(gridID) == GRID_SPECTRAL )
	{
	  cdoAbort("Spectral data unsupported!");
	}
      else
	{
	  ngp = gridInqSize(gridID);
	  break;
	}
    }

  /* check gridsize */
  for ( i = 0; i < ngrids; i++ )
    {
      gridID = vlistGrid(vlistID1, i);
      if ( gridInqType(gridID) != GRID_SPECTRAL )
	{
	  if ( ngp != gridInqSize(gridID) )
	    cdoAbort("Grids have different size!");
	}
    }


  nzaxis  = vlistNzaxis(vlistID1);
  lhavevct = FALSE;

  if ( cdoVerbose )
    cdoPrint("nzaxis: %d", nzaxis);

  for ( i = 0; i < nzaxis; i++ )
    {
      zaxisID = vlistZaxis(vlistID1, i);
      nlevel  = zaxisInqSize(zaxisID);
      if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID )
	{
	  if ( nlevel > 1 )
	    {
	      nvct = zaxisInqVctSize(zaxisID);

              if ( cdoVerbose )
                cdoPrint("i: %d, vct size of zaxisID %d = %d", i, zaxisID, nvct);

	      if ( nlevel == (nvct/2 - 1) )
		{
		  if ( lhavevct == FALSE )
		    {
		      lhavevct = TRUE;
		      zaxisIDh = zaxisID;
		      nhlevf   = nlevel;
	      
                      if ( cdoVerbose )
                        cdoPrint("lhavevct=TRUE  zaxisIDh = %d, nhlevf   = %d", zaxisIDh, nlevel);
 
		      vct = malloc(nvct*sizeof(double));
		      zaxisInqVct(zaxisID, vct);

		      if ( cdoVerbose )
			for ( i = 0; i < nvct/2; ++i )
			  cdoPrint("vct: %5d %25.17f %25.17f", i, vct[i], vct[nvct/2+i]);
		    }
		}
              else 
                {
		  if ( cdoVerbose )
		    cdoPrint("nlevel /= (nvct/2 - 1): nlevel = %d", nlevel);
                }
	    }
	}
    }

  if ( zaxisIDh == -1 )
    cdoAbort("No data on hybrid model level found!");

  nvars = vlistNvars(vlistID1);

  useTable = FALSE;
  for ( varID = 0; varID < nvars; varID++ )
    {
      tableNum = tableInqNum(vlistInqVarTable(vlistID1, varID));

      if ( tableNum > 0  && tableNum != 255 )
	{
	  useTable = TRUE;
	  break;
	}
    }

  if ( cdoVerbose && useTable ) cdoPrint("Using code tables!");

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      zaxisID  = vlistInqVarZaxis(vlistID1, varID);
      nlevel   = zaxisInqSize(zaxisID);
      instNum  = institutInqCenter(vlistInqVarInstitut(vlistID1, varID));
      tableNum = tableInqNum(vlistInqVarTable(vlistID1, varID));

      code     = vlistInqVarCode(vlistID1, varID);
      param    = vlistInqVarParam(vlistID1, varID);

      cdiParamToString(param, paramstr, sizeof(paramstr));

      if ( useTable )
	{
	  if ( tableNum == 2 )
	    {
	      mode = WMO_MODE;
	      geop_code  =   6;
	      temp_code  =  11;
	      hum_code   =  51;
	      ps_code    =   1;
	    }
	  else if ( tableNum == 128 || tableNum == 0 )
	    {
	      mode = ECHAM_MODE;
	      geop_code  = 129;
	      temp_code  = 130;
	      hum_code   = 133;
	      ps_code    = 134;
	      lsp_code   = 152;
	    }
	  else
	    mode = -1;
	}
      else
	{
	  mode = ECHAM_MODE;
	  geop_code  = 129;
	  temp_code  = 130;
	  hum_code   = 133;
	  ps_code    = 134;
	  lsp_code   = 152;
	}

      if ( cdoVerbose )
	cdoPrint("Mode = %d  Center = %d  Param = %s", mode, instNum, paramstr);

      if ( code <= 0 || code == 255 )
	{
	  vlistInqVarName(vlistID1, varID, varname);
	  strtolower(varname);

	  vlistInqVarStdname(vlistID1, varID, stdname);
	  strtolower(stdname);

	  code = echamcode_from_stdname(stdname);

	  if ( code < 0 )
	    {
	      if      ( geopID == -1  && strcmp(varname, "geosp")   == 0 ) code = 129;
	      else if ( psID   == -1  && strcmp(varname, "aps")     == 0 ) code = 134;
	      else if ( psID   == -1  && strcmp(varname, "ps")      == 0 ) code = 134;
	      else if ( lnpsID == -1  && strcmp(varname, "lsp")     == 0 ) code = 152;
	      else if ( tempID == -1  && strcmp(varname, "t")       == 0 ) code = 130;
	      else if ( humID  == -1  && strcmp(varname, "q")       == 0 ) code = 133;
	      // else if ( strcmp(varname, "clwc")    == 0 ) code = 246;
	      // else if ( strcmp(varname, "ciwc")    == 0 ) code = 247;
	    }
	}

      if      ( code == geop_code && nlevel == 1      ) geopID    = varID;
      else if ( code == temp_code && nlevel == nhlevf ) tempID    = varID;
      else if ( code == hum_code  && nlevel == nhlevf ) humID     = varID;
      else if ( code == ps_code   && nlevel == 1      ) psID      = varID;
      else if ( code == lsp_code  && nlevel == 1      ) lnpsID    = varID;
      // else if ( code == 246 ) clwcID    = varID;
      // else if ( code == 247 ) ciwcID    = varID;

      if ( gridInqType(gridID) == GRID_SPECTRAL && zaxisInqType(zaxisID) == ZAXIS_HYBRID )
	cdoAbort("Spectral data on model level unsupported!");

      if ( gridInqType(gridID) == GRID_SPECTRAL )
	cdoAbort("Spectral data unsupported!");
    }

  if ( tempID == -1 ) cdoAbort("Temperature not found!");

  array  = malloc(ngp*sizeof(double));

  geop   = malloc(ngp*sizeof(double));
  ps     = malloc(ngp*sizeof(double));

  temp   = malloc(ngp*nhlevf*sizeof(double));

  if ( humID == -1 )
    cdoWarning("Humidity not found - using algorithm without humidity!");
  else
    hum    = malloc(ngp*nhlevf*sizeof(double));

  // lwater = malloc(ngp*nhlevf*sizeof(double));
  // iwater = malloc(ngp*nhlevf*sizeof(double));

  half_press   = malloc(ngp*(nhlevf+1)*sizeof(double));
  geopotheight = malloc(ngp*(nhlevf+1)*sizeof(double));

  if ( zaxisIDh != -1 && geopID == -1 )
    {
      if ( ltq )
	cdoWarning("Orography (surf. geopotential) not found - using zero orography!");

      memset(geop, 0, ngp*sizeof(double));
    }

  presID = lnpsID;
  if ( zaxisIDh != -1 && lnpsID == -1 )
    {
      presID = psID;
      if ( psID != -1 )
	cdoWarning("LOG surface pressure (lsp) not found - using surface pressure (asp)!");
      else
	cdoAbort("Surface pressure not found!");
    }


  vlistID2 = vlistCreate();

  int var_id = -1;

  if ( operatorID == GEOPOTHEIGHT )
    {
      var_id = geopotential_height;
      varID  = vlistDefVar(vlistID2, gridID, zaxisIDh, TSTEP_INSTANT);
    }
  else if ( operatorID == SEALEVELPRESSURE )
    {
      var_id = air_pressure_at_sea_level;
      varID  = vlistDefVar(vlistID2, gridID, zaxisIDh, TSTEP_INSTANT);
    }
  else
    cdoAbort("Internal problem, invalid operatorID: %d!", operatorID);
  
  vlistDefVarParam(vlistID2, varID, cdiEncodeParam(var_echamcode(var_id), 128, 255));
  vlistDefVarName(vlistID2, varID, var_name(var_id));
  vlistDefVarStdname(vlistID2, varID, var_stdname(var_id));
  vlistDefVarUnits(vlistID2, varID, var_units(var_id));

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  zaxisID  = vlistInqVarZaxis(vlistID1, varID);
	  nlevel   = zaxisInqSize(zaxisID);
	  offset   = gridsize*levelID;
	  streamReadRecord(streamID1, array, &nmiss);

	  if ( zaxisIDh != -1 )
	    {
	      if ( varID == geopID )
		{
		  memcpy(geop, array, ngp*sizeof(double));
		}
	      else if ( varID == presID )
		{
		  if ( lnpsID != -1 )
		    for ( i = 0; i < ngp; ++i ) ps[i] = exp(array[i]);
		  else if ( psID != -1 )
		    memcpy(ps, array, ngp*sizeof(double));
		}
	      else if ( varID == tempID )
		memcpy(temp+offset, array, ngp*sizeof(double));
	      else if ( varID == humID )
		memcpy(hum+offset, array, ngp*sizeof(double));
	      /*
	      else if ( varID == clwcID )
		memcpy(lwater+offset, array, ngp*sizeof(double));
	      else if ( varID == ciwcID )
		memcpy(iwater+offset, array, ngp*sizeof(double));
	      */
	    }
	}

      if ( zaxisIDh != -1 )
	{
	  /* check range of ps_prog */
	  minmaxval(ngp, ps, NULL, &minval, &maxval);
	  if ( minval < MIN_PS || maxval > MAX_PS )
	    cdoWarning("Surface pressure out of range (min=%g max=%g)!", minval, maxval);

	  /* check range of geop */
	  minmaxval(ngp, geop, NULL, &minval, &maxval);
	  if ( minval < MIN_FIS || maxval > MAX_FIS )
	    cdoWarning("Orography out of range (min=%g max=%g)!", minval, maxval);
	}

      varID = tempID;
      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      for ( levelID = 0; levelID < nlevel; levelID++ )
	{
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  offset   = gridsize*levelID;
	  single2  = temp + offset;

	  minmaxval(ngp, single2, NULL, &minval, &maxval);
	  if ( minval < MIN_T || maxval > MAX_T )
	    cdoWarning("Input temperature at level %d out of range (min=%g max=%g)!",
		       levelID+1, minval, maxval);
	}

      if ( humID != -1 )
	{
	  varID = humID;
	  nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	      offset   = gridsize*levelID;
	      single2  = hum + offset;

	      // corr_hum(gridsize, single2, MIN_Q);

	      minmaxval(ngp, single2, NULL, &minval, &maxval);
	      if ( minval < -0.1 || maxval > MAX_Q )
		cdoWarning("Input humidity at level %d out of range (min=%g max=%g)!",
			   levelID+1, minval, maxval);
	    }
	}

      presh(NULL, half_press, vct, ps, nhlevf, ngp);

      memcpy(geopotheight+ngp*nhlevf, geop, ngp*sizeof(double));
      MakeGeopotHeight(geopotheight, temp, hum, half_press, ngp, nhlevf);

      nmissout = 0;
      varID = 0;
      nlevel = nhlevf;
      for ( levelID = 0; levelID < nlevel; levelID++ )
	{
	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, geopotheight+levelID*ngp, nmissout);
	}

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  vlistDestroy(vlistID2);

  free(ps);
  free(geop);
  free(temp);
  free(geopotheight);
  if ( hum ) free(hum);

  if ( half_press ) free(half_press);

  free(array);
  if ( vct ) free(vct);

  cdoFinish();

  return (0);
}

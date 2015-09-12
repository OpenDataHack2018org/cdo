/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2015 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
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

      Vertint    ap2pl           Model air pressure level to pressure level interpolation
*/


#include <ctype.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "after_vertint.h"
#include "list.h"
#include "stdnametable.h"

static
int is_height_axis(int zaxisID, int nlevel)
{
  int isheight = FALSE;
  if ( nlevel > 1 )
    {
      if ( zaxisInqType(zaxisID) == ZAXIS_REFERENCE )
        {
          char units[CDI_MAX_NAME];
          char stdname[CDI_MAX_NAME];
          zaxisInqUnits(zaxisID, units);
          zaxisInqStdname(zaxisID, stdname);
          if ( strcmp(stdname, "height") == 0 && *units == 0 )
            isheight = TRUE;
        }
    }
  return isheight;
}

static
int is_hybrid_axis(int zaxisID, int nlevel)
{
  int ishybrid = FALSE;
  if ( (zaxisInqType(zaxisID) == ZAXIS_HYBRID || zaxisInqType(zaxisID) == ZAXIS_HYBRID_HALF) &&
       nlevel > 1 ) ishybrid = TRUE;
  return ishybrid;
}


void *Vertintap(void *argument)
{
  enum {func_pl, func_hl};
  enum {type_lin, type_log};
  int recID, nrecs;
  int i, k, offset;
  int varID, levelID;
  int zaxisIDp, zaxisIDh = -1, nzaxis;
  int gridID, zaxisID;
  int nhlev = 0, nhlevf = 0, nhlevh = 0, nlevel;
  int *vert_index = NULL;
  int nvct;
  int apressID = -1, dpressID = -1;
  int psID = -1, tempID = -1;
  int param;
  //int sortlevels = TRUE;
  int *pnmiss = NULL;
  char paramstr[32];
  char stdname[CDI_MAX_NAME];
  double minval, maxval;
  double missval;
  double *vct = NULL;
  double *single1, *single2;
  double *ps_prog = NULL, *full_press = NULL, *dpress = NULL;
  double *hyb_press = NULL;
  int Extrapolate = 0;
  int mono_level;
  LIST *flist = listNew(FLT_LIST);

  cdoInitialize(argument);

  int AP2PL     = cdoOperatorAdd("ap2pl",     func_pl, type_lin, "pressure levels in pascal");
  int AP2PLX    = cdoOperatorAdd("ap2plx",    func_pl, type_lin, "pressure levels in pascal");
  int AP2PL_LP  = cdoOperatorAdd("ap2pl_lp",  func_pl, type_log, "pressure levels in pascal");
  int AP2PLX_LP = cdoOperatorAdd("ap2plx_lp", func_pl, type_log, "pressure levels in pascal");

  int operatorID = cdoOperatorID();
  int operfunc   = cdoOperatorF1(operatorID);
  int opertype   = cdoOperatorF2(operatorID);

  if ( operatorID == AP2PL || operatorID == AP2PL_LP )
    {
      char *envstr = getenv("EXTRAPOLATE");

      if ( envstr && isdigit((int) envstr[0]) )
	{
          Extrapolate = atoi(envstr);
          if ( Extrapolate == 1 )
            cdoPrint("Extrapolation of missing values enabled!");
	}
    }
  else if ( operatorID == AP2PLX || operatorID == AP2PLX_LP )
    {
      Extrapolate = 1;
    }

  operatorInputArg(cdoOperatorEnter(operatorID));

  int nplev = 0;
  double *plev = NULL;
  if ( operatorArgc() == 1 && strcmp(operatorArgv()[0], "default") == 0 )
    {
      /*
      double stdlev[] = {100000, 92500, 85000, 77500, 70000, 60000, 50000, 40000, 30000, 25000, 20000,
                          15000, 10000, 7000, 5000, 3000, 2000, 1000, 700, 500, 300, 200, 100, 50, 20, 10};
      */
      double stdlev[] = {100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000, 20000, 15000,
                          10000,  7000,  5000,  3000,  2000, 1000 };
      nplev = sizeof(stdlev)/sizeof(*stdlev);
      plev  = (double *) malloc(nplev*sizeof(double));
      for ( i = 0; i < nplev; ++i ) plev[i] = stdlev[i];
    }
  else
    {
      nplev = args2fltlist(operatorArgc(), operatorArgv(), flist);
      plev  = (double *) listArrayPtr(flist);
    }

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int gridsize = vlist_check_gridsize(vlistID1);

  if ( operfunc == func_hl )
    zaxisIDp = zaxisCreate(ZAXIS_HEIGHT, nplev);
  else
    zaxisIDp = zaxisCreate(ZAXIS_PRESSURE, nplev);

  zaxisDefLevels(zaxisIDp, plev);
  nzaxis  = vlistNzaxis(vlistID1);

  int lheight = FALSE;
  for ( i = 0; i < nzaxis; i++ )
    {
      mono_level = TRUE;
      zaxisID = vlistZaxis(vlistID1, i);
      nlevel  = zaxisInqSize(zaxisID);

      if ( is_height_axis(zaxisID, nlevel) )
        {
	  double *level = (double *) malloc(nlevel*sizeof(double));
	  zaxisInqLevels(zaxisID, level);
	  int l;
	  for ( l = 0; l < nlevel; l++ )
	    {
	      if ( (l+1) != (int) (level[l]+0.5) ) break;
	    }
	  if ( l == nlevel ) mono_level = TRUE; 
	  free(level);          
        }
      
      if ( is_height_axis(zaxisID, nlevel) && mono_level )
        {
          lheight  = TRUE;
          zaxisIDh = zaxisID;
          nhlev    = nlevel;
          nhlevf   = nhlev;
          nhlevh   = nhlevf + 1;

          vlistChangeZaxisIndex(vlistID2, i, zaxisIDp);
          break;
        }
    }

  if ( ! lheight ) { // old version, not needed anymore!!
  int lhavevct = FALSE;
  for ( i = 0; i < nzaxis; i++ )
    {
      /* mono_level = FALSE; */
      mono_level = TRUE;
      zaxisID = vlistZaxis(vlistID1, i);
      nlevel  = zaxisInqSize(zaxisID);

      if ( is_hybrid_axis(zaxisID, nlevel) )
	{
	  double *level = (double *) malloc(nlevel*sizeof(double));
	  zaxisInqLevels(zaxisID, level);
	  int l;
	  for ( l = 0; l < nlevel; l++ )
	    {
	      if ( (l+1) != (int) (level[l]+0.5) ) break;
	    }
	  if ( l == nlevel ) mono_level = TRUE; 
	  free(level);
	}

      if ( is_hybrid_axis(zaxisID, nlevel) && mono_level )
	{
	  nvct = zaxisInqVctSize(zaxisID);
	  if ( nlevel == (nvct/2 - 1) )
	    {
	      if ( lhavevct == FALSE )
		{
		  lhavevct = TRUE;
		  zaxisIDh = zaxisID;
		  nhlev    = nlevel;
		  nhlevf   = nhlev;
		  nhlevh   = nhlevf + 1;
	      
		  vct = (double*) malloc(nvct*sizeof(double));
		  zaxisInqVct(zaxisID, vct);

		  vlistChangeZaxisIndex(vlistID2, i, zaxisIDp);
		}
	      else
		{
		  if ( memcmp(vct, zaxisInqVctPtr(zaxisID), nvct*sizeof(double)) == 0 )
		    vlistChangeZaxisIndex(vlistID2, i, zaxisIDp);
		}
	    }
	  else if ( nlevel == (nvct/2) )
	    {
	      if ( lhavevct == FALSE )
		{
		  lhavevct = TRUE;
		  zaxisIDh = zaxisID;
		  nhlev    = nlevel;
		  nhlevf   = nhlev - 1;
		  nhlevh   = nhlev;
	      
		  vct = (double*) malloc(nvct*sizeof(double));
		  zaxisInqVct(zaxisID, vct);

		  vlistChangeZaxisIndex(vlistID2, i, zaxisIDp);
		}
	      else
		{
		  if ( memcmp(vct, zaxisInqVctPtr(zaxisID), nvct*sizeof(double)) == 0 )
		    vlistChangeZaxisIndex(vlistID2, i, zaxisIDp);
		}
	    }
	}
    }
  }
  
  int nvars = vlistNvars(vlistID1);

  int vars[nvars];
  double *vardata1[nvars];
  double *vardata2[nvars];
  int *varnmiss[nvars];
  int varinterp[nvars];

  int maxlev = nhlevh > nplev ? nhlevh : nplev;

  if ( Extrapolate == 0 )
    pnmiss = (int *) malloc(nplev*sizeof(int));

  // check levels
  if ( zaxisIDh != -1 )
    {
      int nlev = zaxisInqSize(zaxisIDh);
      if ( nlev != nhlev ) cdoAbort("Internal error, wrong number of height level!");
      double levels[nlev];
      zaxisInqLevels(zaxisIDh, levels);

      for ( int ilev = 0; ilev < nlev; ++ilev )
	{
	  if ( (ilev+1) != (int)levels[ilev] )
	    {
	      //sortlevels = FALSE;
	      break;
	    }
	}
    }

  if ( zaxisIDh != -1 && gridsize > 0 )
    {
      vert_index = (int *) malloc(gridsize*nplev*sizeof(int));
      ps_prog    = (double *) malloc(gridsize*sizeof(double));
      full_press = (double *) malloc(gridsize*nhlevf*sizeof(double));
      dpress     = (double *) malloc(gridsize*nhlevf*sizeof(double));
    }
  else
    cdoWarning("No 3D variable with generalized height level found!");

  if ( operfunc == func_hl )
    {
      double phlev[nplev];
      height2pressure(phlev, plev, nplev);

      if ( cdoVerbose )
	for ( i = 0; i < nplev; ++i )
	  cdoPrint("level = %d   height = %g   pressure = %g", i+1, plev[i], phlev[i]);

      memcpy(plev, phlev, nplev*sizeof(double));
    }

  if ( opertype == type_log )
    for ( k = 0; k < nplev; k++ ) plev[k] = log(plev[k]);

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      zaxisID  = vlistInqVarZaxis(vlistID1, varID);
      nlevel   = zaxisInqSize(zaxisID);

      vlistInqVarStdname(vlistID1, varID, stdname);
      strtolower(stdname);

      if      ( strcmp(stdname, var_stdname(surface_air_pressure)) == 0 ) psID = varID; 
      else if ( strcmp(stdname, var_stdname(air_pressure))         == 0 ) apressID = varID; 
      else if ( strcmp(stdname, var_stdname(pressure_thickness))   == 0 ) dpressID = varID; 
      else if ( strcmp(stdname, var_stdname(air_temperature))      == 0 ) tempID = varID; 

      if ( gridInqType(gridID) == GRID_SPECTRAL )
	cdoAbort("Spectral data unsupported!");

      vardata1[varID] = (double *) malloc(gridsize*nlevel*sizeof(double));

      if ( zaxisID == zaxisIDh )
	{
	  varinterp[varID] = TRUE;
	  vardata2[varID]  = (double *) malloc(gridsize*nplev*sizeof(double));
	  varnmiss[varID]  = (int *) malloc(maxlev*sizeof(int));
	  memset(varnmiss[varID], 0, maxlev*sizeof(int));
	}
      else
	{
	  varinterp[varID] = FALSE;
	  vardata2[varID]  = vardata1[varID];
	  varnmiss[varID]  = (int *) malloc(nlevel*sizeof(int));
	}
    }

  if ( cdoVerbose )
    {
      cdoPrint("Found:");
      if ( psID      != -1 ) cdoPrint("  %s", var_stdname(surface_air_pressure));
      if ( apressID  != -1 ) cdoPrint("  %s", var_stdname(air_pressure));
      if ( dpressID  != -1 ) cdoPrint("  %s", var_stdname(pressure_thickness));
      if ( tempID    != -1 ) cdoPrint("  %s", var_stdname(air_temperature));
    }

  if ( apressID == -1 ) cdoAbort("%s not found!", var_stdname(air_pressure));
  if ( zaxisIDh != -1 && psID == -1 && dpressID )
    cdoWarning("Surface pressure not found - set to upper level of %s!", var_stdname(air_pressure));

  for ( varID = 0; varID < nvars; ++varID )
    {
      if ( varinterp[varID] == TRUE && vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT )
	vlistDefVarTsteptype(vlistID2, varID, TSTEP_INSTANT);
    }

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      for ( varID = 0; varID < nvars; ++varID )
	{
          vars[varID] = FALSE;
	  nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    varnmiss[varID][levelID] = 0;
	}

      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  /*
	  zaxisID  = vlistInqVarZaxis(vlistID1, varID);
	  nlevel   = zaxisInqSize(zaxisID);
	  if ( sortlevels && zaxisIDh != -1 && zaxisID == zaxisIDh && nlevel == nhlev )
	    {
	      levelID = (int) (zaxisInqLevel(zaxisIDh, levelID)-1);
	      printf("levelID %d\n", levelID);
	    }
	  */
	  offset   = gridsize*levelID;
	  single1  = vardata1[varID] + offset;

	  streamReadRecord(streamID1, single1, &varnmiss[varID][levelID]);

	  vars[varID] = TRUE;
	}

      for ( varID = 0; varID < nvars; varID++ )
	if ( varinterp[varID] == TRUE ) vars[varID] = TRUE;

      if ( zaxisIDh != -1 )
	{
	  if ( psID != -1 )
	    {
	      memcpy(ps_prog, vardata1[psID], gridsize*sizeof(double)); 
	    }
          else if ( dpressID != -1 )
	    {
	      memcpy(dpress, vardata1[dpressID], gridsize*nhlevf*sizeof(double)); 
	      for ( i = 0; i < gridsize; i++ )  ps_prog[i] = 0;
	      for ( k = 0; k < nhlevf; ++k )
		for ( i = 0; i < gridsize; i++ )
		  ps_prog[i] += dpress[k*gridsize+i];
	    }
	  else
	    {
	      memcpy(ps_prog, vardata1[apressID]+gridsize*(nhlevf-1), gridsize*sizeof(double)); 
	      //for ( i = 0; i < gridsize; i++ )  ps_prog[i] = 110000;
	    }

	  /* check range of ps_prog */
	  minmaxval(gridsize, ps_prog, NULL, &minval, &maxval);
	  if ( minval < MIN_PS || maxval > MAX_PS )
	    cdoWarning("Surface pressure out of range (min=%g max=%g)!", minval, maxval);

	  memcpy(full_press, vardata1[apressID], gridsize*nhlevf*sizeof(double)); 

	  if ( opertype == type_log )
	    {
	      for ( i = 0; i < gridsize; i++ ) ps_prog[i] = log(ps_prog[i]);

	      for ( k = 0; k < nhlevf; k++ )
		for ( i = 0; i < gridsize; i++ )
		  full_press[k*gridsize+i] = log(full_press[k*gridsize+i]);
	    }

	  genind(vert_index, plev, full_press, gridsize, nplev, nhlevf);

	  if ( Extrapolate == 0 )
	    genindmiss(vert_index, plev, gridsize, nplev, ps_prog, pnmiss);
	}

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vars[varID] )
	    {
	      zaxisID  = vlistInqVarZaxis(vlistID1, varID);
	      missval  = vlistInqVarMissval(vlistID1, varID);
	      nlevel   = zaxisInqSize(zaxisID);
	      if ( varinterp[varID] )
		{
		  if ( nlevel == nhlevf )
		    {
		      hyb_press = full_press;
		    }
		  else
		    {
		      param = vlistInqVarParam(vlistID1, varID);
		      cdiParamToString(param, paramstr, sizeof(paramstr));
		      cdoAbort("Number of generalized height level differ from full/half level (param=%s)!", paramstr);
		    }

		  for ( levelID = 0; levelID < nlevel; levelID++ )
		    {
		      if ( varnmiss[varID][levelID] )
			cdoAbort("Missing values unsupported for this operator!");
		    }

		  interp_X(vardata1[varID], vardata2[varID], hyb_press,
			   vert_index, plev, nplev, gridsize, nlevel, missval);
		  
		  if ( Extrapolate == 0 )
		    memcpy(varnmiss[varID], pnmiss, nplev*sizeof(int));
		}
	    }
	}

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vars[varID] )
	    {
	      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  offset   = gridsize*levelID;
		  single2  = vardata2[varID] + offset;
		  streamDefRecord(streamID2, varID, levelID);
		  streamWriteRecord(streamID2, single2, varnmiss[varID][levelID]);
		}
	    }
	}

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  for ( varID = 0; varID < nvars; varID++ )
    {
      free(varnmiss[varID]);
      free(vardata1[varID]);
      if ( varinterp[varID] ) free(vardata2[varID]);
    }

  if ( pnmiss     ) free(pnmiss);
  if ( ps_prog    ) free(ps_prog);
  if ( vert_index ) free(vert_index);
  if ( full_press ) free(full_press);
  if ( dpress     ) free(dpress);
  if ( vct        ) free(vct);

  listDelete(flist);

  cdoFinish();

  return (0);
}

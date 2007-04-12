/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2007 Uwe Schulzweida, schulzweida@dkrz.de
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

      Remapeta     remapeta          Model to model level interpolation
*/


#include <ctype.h>
#include <string.h>
#include <math.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "functs.h"
#include "vinterp.h"
#include "list.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void hetaeta(int ltq, int ngp,
	     int nlev1, double *ah1, double *bh1,
             double *fis1, double *ps1, 
             double *t1, double *q1,
             int nlev2, double *ah2, double *bh2, 
             double *fis2, double *ps2, 
             double *t2, double *q2,
	     int nvars, double **vars1, double **vars2,
	     double *tscor, double *pscor, double *secor);

#define  MAX_VARS3D  1024

void *Remapeta(void *argument)
{
  static char func[] = "Remapeta";
  int REMAPETA;
  int operatorID;
  int streamID1, streamID2;
  int vlistID1, vlistID2;
  int gridsize, ngp = 0, nfis2gp = 0;
  int recID, nrecs;
  int i, offset, iv;
  int tsID, varID, levelID;
  int nvars, nvars3D = 0;
  int zaxisID2, zaxisIDh = -1, nzaxis;
  int ngrids, gridID, zaxisID;
  int nlevel;
  int nvct1, nvct2 = 0;
  int geopID = -1, tempID = -1, sqID = -1, psID = -1, lnpsID = -1, presID = -1;
  int code;
  char varname[128];
  double *single2;
  int taxisID1, taxisID2;
  int lhavevct;
  int nlev1 = 0, nlev2 = 0, nlev2p1;
  double *lev2;
  double *vct1 = NULL, *vct2 = NULL;
  double *a1 = NULL, *b1 = NULL, *a2 = NULL, *b2 = NULL;
  double *fis1 = NULL, *ps1 = NULL, *t1 = NULL, *q1 = NULL;
  double *fis2 = NULL, *ps2 = NULL, *t2 = NULL, *q2 = NULL;
  double *tscor = NULL, *pscor = NULL, *secor = NULL;
  int nmiss;
  int ltq = FALSE;
  int lfis2 = FALSE;
  int varids[MAX_VARS3D];
  double *array = NULL;
  double **vars1 = NULL, **vars2 = NULL;

  cdoInitialize(argument);

  REMAPETA = cdoOperatorAdd("remapeta", 0, 0, "VCT file name");

  operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  if ( operatorID == REMAPETA )
    {
      const char *fname = operatorArgv()[0];
      char line[1024], *pline;
      int num, i = 0;
      int maxvct = 8192;
      
      FILE *fp;

      fp = fopen(fname, "r");
      if ( fp == NULL ) { perror(fname); exit(EXIT_FAILURE); }

      vct2 = (double *) malloc(maxvct*sizeof(double));

      while ( readline(fp, line, 1024) )
	{
          if ( line[0] == '#' ) continue;
          if ( line[0] == '\0' ) continue;

	  pline = line;
	  num = (int) strtod(pline, &pline);
	  if ( pline == NULL ) cdoAbort("Format error in VCT file %s!", fname);
	  if ( num != i ) cdoWarning("Inconsistent VCT file, entry %d is %d.", i, num);

	  vct2[i] = strtod(pline, &pline);
	  if ( pline == NULL ) cdoAbort("Format error in VCT file %s!", fname);

	  vct2[i+maxvct/2] = strtod(pline, &pline);

	  i++;
	}

      fclose(fp);

      a2 = vct2;
      b2 = vct2 + i;
      nvct2 = 2*i;
      nlev2 = i - 1;
      nlev2p1 = nlev2 + 1;

      for ( i = 0; i < nlev2+1; ++i )
	vct2[i+nvct2/2] = vct2[i+maxvct/2];

      vct2 = (double *) realloc(vct2, 2*i*sizeof(double));

      if ( cdoVerbose )
	for ( i = 0; i < nlev2+1; ++i )
	  fprintf(stdout, "vct2: %5d %25.17f %25.17f\n", i, vct2[i], vct2[nvct2/2+i]);

      if ( operatorArgc() == 2 )
	{
	  lfis2 = TRUE;
	  fname = operatorArgv()[1];

	  streamID1 = streamOpenRead(fname);
	  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

	  vlistID1 = streamInqVlist(streamID1);

	  streamInqRecord(streamID1, &varID, &levelID);

	  gridID = vlistInqVarGrid(vlistID1, varID);
	  nfis2gp = gridInqSize(gridID);

	  fis2  = (double *) malloc(nfis2gp*sizeof(double));

	  streamReadRecord(streamID1, fis2, &nmiss);

	  /* check range of geop */
	  {
	    double minval =  DBL_MAX;
	    double maxval = -DBL_MAX;
	    for ( i = 0; i < nfis2gp; i++ )
	      {
		if      ( fis2[i] > maxval ) maxval = fis2[i];
		else if ( fis2[i] < minval ) minval = fis2[i];
	      }

	    if ( minval < -100000 || maxval > 100000 )
	      cdoWarning("Orography out of range (min=%g max=%g)!\n", minval, maxval);

	    if ( minval < -1.e10 || maxval > 1.e10 )
	      cdoAbort("Orography out of range!\n");
	  }

	  streamClose(streamID1); 
	}
    }

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

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

  zaxisID2 = zaxisCreate(ZAXIS_HYBRID, nlev2);
  lev2 = (double *) malloc(nlev2*sizeof(double));
  for ( i = 0; i < nlev2; ++i ) lev2[i] = i+1;
  zaxisDefLevels(zaxisID2, lev2);
  free(lev2);

  if ( nvct2 == 0 ) cdoAbort("Internal problem, vct2 undefined!");
  zaxisDefVct(zaxisID2, nvct2, vct2);

  nzaxis  = vlistNzaxis(vlistID1);
  lhavevct = FALSE;
  for ( i = 0; i < nzaxis; i++ )
    {
      zaxisID = vlistZaxis(vlistID1, i);
      nlevel  = zaxisInqSize(zaxisID);
      if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID && nlevel > 1 )
	{
	  nvct1 = zaxisInqVctSize(zaxisID);
	  if ( nlevel == (nvct1/2 - 1) )
	    {
	      if ( lhavevct == FALSE )
		{
		  lhavevct = TRUE;
		  zaxisIDh = zaxisID;
		  nlev1    = nlevel;
	      
		  vct1 = (double *) malloc(nvct1*sizeof(double));
		  memcpy(vct1, zaxisInqVctPtr(zaxisID), nvct1*sizeof(double));

		  vlistChangeZaxisIndex(vlistID2, i, zaxisID2);

		  a1 = vct1;
		  b1 = vct1 + nvct1/2;
		  if ( cdoVerbose )
		    for ( i = 0; i < nvct1/2; ++i )
		      fprintf(stdout, "vct1: %5d %25.17f %25.17f\n", i, vct1[i], vct1[nvct1/2+i]);
		}
	      else
		{
		  if ( memcmp(vct1, zaxisInqVctPtr(zaxisID), nvct1*sizeof(double)) == 0 )
		    vlistChangeZaxisIndex(vlistID2, i, zaxisID2);
		}
	    }
	}
    }

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);


  if ( zaxisIDh == -1 )
    cdoWarning("No data on hybrid model level found!");

  nvars = vlistNvars(vlistID1);

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID  = vlistInqVarGrid(vlistID1, varID);
      zaxisID = vlistInqVarZaxis(vlistID1, varID);

      code = vlistInqVarCode(vlistID1, varID);
      if ( code <= 0 )
	{
	  vlistInqVarName(vlistID1, varID, varname);

	  strtolower(varname);

	  if      ( strcmp(varname, "geosp")   == 0 ) code = 129;
	  else if ( strcmp(varname, "st")      == 0 ) code = 130;
	  else if ( strcmp(varname, "sq")      == 0 ) code = 133;
	  else if ( strcmp(varname, "aps")     == 0 ) code = 134;
	  else if ( strcmp(varname, "lsp")     == 0 ) code = 152;
	}

      if      ( code == 129 ) geopID    = varID;
      else if ( code == 130 ) tempID    = varID;
      else if ( code == 133 ) sqID      = varID;
      else if ( code == 134 ) psID      = varID;
      else if ( code == 152 ) lnpsID    = varID;

      if ( gridInqType(gridID) == GRID_SPECTRAL && zaxisInqType(zaxisID) == ZAXIS_HYBRID )
	cdoAbort("Spectral data on model level unsupported!");

      if ( gridInqType(gridID) == GRID_SPECTRAL )
	cdoAbort("Spectral data unsupported!");

      nlevel   = zaxisInqSize(zaxisID);

      if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID && zaxisIDh != -1 && nlevel == nlev1 )
	{
	  if ( ! (code == 130 || code == 133) )
	    varids[nvars3D++] = varID;
	}
      else
	{
	  if ( code == 130 ) tempID = -1;
	  if ( code == 133 ) sqID   = -1;
	}
    }

  if ( tempID != -1 && sqID != -1 )
    {
      ltq = TRUE;
    }
  else
    {
      if ( tempID != -1 ) cdoAbort("Temperature without humidity unsupported!");
      if ( sqID   != -1 ) cdoAbort("Humidity without temperature unsupported!");
    }

  array = (double *) malloc(ngp*sizeof(double));

  fis1  = (double *) malloc(ngp*sizeof(double));
  ps1   = (double *) malloc(ngp*sizeof(double));

  if ( lfis2 == FALSE ) fis2  = (double *) malloc(ngp*sizeof(double));
  if ( lfis2 == TRUE && ngp != nfis2gp ) cdoAbort("Orographies have different grid size!");

  ps2   = (double *) malloc(ngp*sizeof(double));

  if ( ltq )
    {
      tscor = (double *) malloc(ngp*sizeof(double));
      pscor = (double *) malloc(ngp*sizeof(double));
      secor = (double *) malloc(ngp*sizeof(double));

      t1    = (double *) malloc(ngp*nlev1*sizeof(double));
      q1    = (double *) malloc(ngp*nlev1*sizeof(double));

      t2    = (double *) malloc(ngp*nlev2*sizeof(double));
      q2    = (double *) malloc(ngp*nlev2*sizeof(double));
    }

  if ( nvars3D )
    {
      vars1  = (double **) malloc(nvars*sizeof(double));
      vars2  = (double **) malloc(nvars*sizeof(double));

      for ( varID = 0; varID < nvars3D; ++varID )
	{
	  vars1[varID] = (double *) malloc(ngp*nlev1*sizeof(double));
	  vars2[varID] = (double *) malloc(ngp*nlev2*sizeof(double));
	}
    }

  if ( zaxisIDh != -1 && geopID == -1 )
    {
      cdoWarning("Orography not found - using zero orography!");
      memset(fis1, 0, ngp*sizeof(double));
    }

  presID = lnpsID;
  if ( zaxisIDh != -1 && lnpsID == -1 )
    {
      presID = psID;
      if ( psID != -1 )
	cdoWarning("LOG surface pressure (code 152) not found - using surface pressure (code 134)!");
      else
	cdoAbort("Surface pressure not found!");
    }

  if ( cdoVerbose ) cdoPrint("nvars3D = %d   ltq = %d\n", nvars3D, ltq);

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
		memcpy(fis1, array, ngp*sizeof(double));
	      else if ( varID == presID )
		{
		  if ( lnpsID != -1 )
		    for ( i = 0; i < ngp; ++i ) ps1[i] = exp(array[i]);
		  else if ( psID != -1 )
		    memcpy(ps1, array, ngp*sizeof(double));
		}
	      else if ( ltq && varID == tempID )
		memcpy(t1+offset, array, ngp*sizeof(double));
	      else if ( ltq && varID == sqID )
		memcpy(q1+offset, array, ngp*sizeof(double));
	      else if ( zaxisID == zaxisIDh )
		{
		  for ( i = 0; i < nvars3D; ++i )
		    if ( varID == varids[i] ) break;

		  if ( i == nvars3D ) cdoAbort("Internal error, 3D variable not found!");

		  memcpy(vars1[i]+offset, array, ngp*sizeof(double));
		}
	      else
		{
		  streamDefRecord(streamID2, varID, levelID);
		  streamWriteRecord(streamID2, array, nmiss);
		}
	    }
	  else
	    {
	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, array, nmiss);
	    }
	}

      if ( zaxisIDh != -1 )
	{
	  /* check range of ps_prog */
	  {
	    double minval =  DBL_MAX;
	    double maxval = -DBL_MAX;
	    for ( i = 0; i < ngp; i++ )
	      {
		if      ( ps1[i] > maxval ) maxval = ps1[i];
		else if ( ps1[i] < minval ) minval = ps1[i];
	      }

	    if ( minval < 20000 || maxval > 150000 )
	      cdoWarning("Surface pressure out of range (min=%g max=%g)!\n", minval, maxval);

	    if ( minval < -1.e10 || maxval > 1.e10 )
	      cdoAbort("Surface pressure out of range!\n");
	  }
	  /* check range of geop */
	  {
	    double minval =  DBL_MAX;
	    double maxval = -DBL_MAX;
	    for ( i = 0; i < ngp; i++ )
	      {
		if      ( fis1[i] > maxval ) maxval = fis1[i];
		else if ( fis1[i] < minval ) minval = fis1[i];
	      }

	    if ( minval < -100000 || maxval > 100000 )
	      cdoWarning("Orography out of range (min=%g max=%g)!\n", minval, maxval);

	    if ( minval < -1.e10 || maxval > 1.e10 )
	      cdoAbort("Orography out of range!\n");
	  }
	}

      if ( lfis2 == FALSE )
	for ( i = 0; i < ngp; i++ ) fis2[i] = fis1[i];

      if ( nvars3D || ltq )
	hetaeta(ltq, ngp,
		nlev1, a1, b1,
		fis1, ps1,
		t1, q1,
		nlev2, a2, b2,
		fis2, ps2,
		t2, q2,
		nvars3D, vars1, vars2,
		tscor, pscor, secor);


      varID   = geopID;
      levelID = 0;
      nmiss   = 0;
      streamDefRecord(streamID2, varID, levelID);
      streamWriteRecord(streamID2, fis2, nmiss);

      if ( lnpsID != -1 )
	for ( i = 0; i < ngp; ++i ) ps2[i] = log(ps2[i]);

      varID   = presID;
      levelID = 0;
      nmiss   = 0;
      streamDefRecord(streamID2, varID, levelID);
      streamWriteRecord(streamID2, ps2, nmiss);

      if ( ltq )
	{
	  varID = tempID;
	  nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	      offset   = gridsize*levelID;
	      single2  = t2 + offset;
	      nmiss    = 0;

	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, single2, nmiss);
	    }

	  varID = sqID;
	  nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	      offset   = gridsize*levelID;
	      single2  = q2 + offset;
	      nmiss    = 0;

	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, single2, nmiss);
	    }
	}

      for ( iv = 0; iv < nvars3D; ++iv )
	{
	  varID = varids[iv];
	  nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	      offset   = gridsize*levelID;
	      single2  = vars2[iv] + offset;
	      nmiss    = 0;

	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, single2, nmiss);
	    }
	}

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( nvars3D )
    {
      for ( varID = 0; varID < nvars3D; varID++ )
	{
	  free(vars2[varID]);
	  free(vars1[varID]);
	}
      free(vars2);
      free(vars1);
    }

  if ( ltq )
    {
      free(q2);
      free(t2);
      free(q1);
      free(t1);
      free(secor);
      free(pscor);
      free(tscor);
    }

  free(ps2);
  free(fis2);
  free(ps1);
  free(fis1);

  free(array);

  cdoFinish();

  return (0);
}

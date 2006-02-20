/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2006 Uwe Schulzweida, schulzweida@dkrz.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include <string.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "specspace.h"
#include "list.h"

/*
@BeginDoc

@BeginModule

@Name      = Wind
@Title     = Spectral transformation
@Section   = Spectral transformation
@Class     = Transformation
@Arguments = ifile ofile
@Operators = uv2dv dv2uv

@EndModule


@BeginOperator_uv2dv

@Title     = u and v wind to divergence and vorticity

@BeginDesciption
@EndDesciption

@EndOperator


@BeginOperator_dv2uv

@Title     = Divergence and vorticity to u and v wind

@BeginDesciption
@EndDesciption

@EndOperator


@EndDoc
*/

#define  MAX_TRUNC  9999

void *Wind(void *argument)
{
  static char func[] = "Wind";
  int UV2DV, DV2UV;
  int operatorID;
  int streamID1, streamID2;
  int nrecs, nvars;
  int tsID, recID, varID, levelID;
  int gridsize;
  int index, ngrids;
  int vlistID1, vlistID2;
  int gridIDsp = -1, gridIDgp = -1;
  int gridID1 = -1, gridID2 = -1;
  int gridID;
  int nmiss;
  int ncut = 0;
  int *wnums = NULL, waves[MAX_TRUNC];
  int *vars;
  int lcopy = FALSE;
  double *array1 = NULL, *array2 = NULL;
  int taxisID1, taxisID2;
  int nlon, nlat, trunc;
  SPTRANS *sptrans = NULL;
  LIST *ilist = listNew(INT_LIST);

  cdoInitialize(argument);

  UV2DV = cdoOperatorAdd("uv2dv", 0, 0, NULL);
  DV2UV = cdoOperatorAdd("dv2uv", 0, 0, NULL);

  operatorID = cdoOperatorID();

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  ngrids = vlistNgrids(vlistID1);
  /* find first spectral grid */
  /*
  for ( index = 0; index < ngrids; index++ )
    {
      gridID = vlistGrid(vlistID1, index);
      if ( gridInqType(gridID) == GRID_SPECTRAL )
	{
	  gridIDsp = gridID;
	  break;
	}
    }
  */
  /* find first gaussian grid */
  /*
  for ( index = 0; index < ngrids; index++ )
    {
      gridID = vlistGrid(vlistID1, index);
      if ( gridInqType(gridID) == GRID_GAUSSIAN )
	{
	  gridIDgp = gridID;
	  break;
	}
    }
  */
  /* define output grid */
  /*
  if ( operatorID == GP2SP || operatorID == GP2SPL )
    {
      gridID1 = gridIDgp;

      if ( gridID1 != -1 )
	{
	  if ( operatorID == GP2SP )
	    trunc = nlat2trunc(gridInqYsize(gridID1));
	  else
	    trunc = nlat2trunc2(gridInqYsize(gridID1));

	  if ( gridIDsp != -1 )
	    if ( trunc != gridInqTrunc(gridIDsp) ) gridIDsp = -1;

	  if ( gridIDsp == -1 )
	    {
	      gridIDsp = gridNew(GRID_SPECTRAL, (trunc+1)*(trunc+2));
	      gridDefTrunc(gridIDsp, trunc);
	    }
	}

      if ( gridIDsp == -1 && gridInqType(vlistGrid(vlistID1, 0)) == GRID_GAUSSIAN_REDUCED )
	cdoAbort("Gaussian reduced grid found. Use option -R to convert it to a regular grid!");

      if ( gridIDsp == -1 ) cdoAbort("No Gaussian grid data found!");

      gridID2 = gridIDsp;

      nlon  = gridInqXsize(gridID1);
      nlat  = gridInqYsize(gridID1);
      trunc = gridInqTrunc(gridID2);

      sptrans = sptrans_new(nlon, nlat, trunc);
    }
  else if ( operatorID == SP2GP || operatorID == SP2GPL )
    {   
      if ( gridIDsp == -1 ) cdoWarning("No spectral data found!");

      gridID1 = gridIDsp;

      if ( gridID1 != -1 )
	{
	  if ( gridIDgp != -1 )
	    {
	      if ( operatorID == SP2GP )
		trunc = nlat2trunc(gridInqYsize(gridIDgp));
	      else
		trunc = nlat2trunc2(gridInqYsize(gridIDgp));

	      if ( gridInqTrunc(gridIDsp) != trunc ) gridIDgp = -1;
	    }

	  if ( gridIDgp == -1 )
	    {
	      char gridname[20];
	      if ( operatorID == SP2GP )
		sprintf(gridname, "t%dgrid", gridInqTrunc(gridIDsp));
	      else
		sprintf(gridname, "tl%dgrid", gridInqTrunc(gridIDsp));

	      gridIDgp = gridFromName(gridname);
	    }

	  gridID2 = gridIDgp;

	  trunc = gridInqTrunc(gridID1);
	  nlon  = gridInqXsize(gridID2);
	  nlat  = gridInqYsize(gridID2);
      
	  sptrans = sptrans_new(nlon, nlat, trunc);
	}
    }
  else if ( operatorID == SP2SP )
    {
      gridID1 = gridIDsp;

      operatorInputArg("truncation");
      if ( gridID1 != -1 )
	{
	  int trunc = atoi(operatorArgv()[0]);
	  int nsp = (trunc+1)*(trunc+2);
	  gridIDsp = gridNew(GRID_SPECTRAL, nsp);
	  gridDefTrunc(gridIDsp, trunc);
	}
      else
	cdoAbort("No spectral data found!");

      gridID2 = gridIDsp;
    }
  else if ( operatorID == SPCUT )
    {
      int i, j;
      gridID1 = gridIDsp;

      operatorInputArg("wave numbers");
      if ( gridID1 != -1 )
	{
	  ncut = args2intlist(operatorArgc(), operatorArgv(), ilist);
	  wnums = (int *) listArrayPtr(ilist);
	  for ( i = 0; i < MAX_TRUNC; i++ ) waves[i] = 1;
	  for ( i = 0; i < ncut; i++ )
	    {
	      j = wnums[i] - 1;
	      if ( j < 0 || j >= MAX_TRUNC )
		cdoAbort("wave number %d out of range!", wnums[i]);
	      waves[j] = 0;
	    }
	}
      else
	cdoAbort("No spectral data found!");

      gridID2 = gridIDsp;
    }
  */
  nvars = vlistNvars(vlistID2);
  vars  = (int *) malloc(nvars*sizeof(int));
  for ( varID = 0; varID < nvars; varID++ )
    {
      if ( gridID1 == vlistInqVarGrid(vlistID1, varID) )
	vars[varID] = TRUE;
      else
	vars[varID] = FALSE;
    }

  if ( gridID1 != -1 ) vlistChangeGrid(vlistID2, gridID1, gridID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);
  array1 = (double *) malloc(gridsize*sizeof(double));

  if ( gridID2 != -1 )
    {
      gridsize = gridInqSize(gridID2);
      array2 = (double *) malloc(gridsize*sizeof(double));
    }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  if ( vars[varID] )
	    {
	      streamReadRecord(streamID1, array1, &nmiss);
	      if ( nmiss ) cdoAbort("missing values unsupported for spectral data!");

	      gridID1 = vlistInqVarGrid(vlistID1, varID);
	      /*
	      if ( operatorID == GP2SP || operatorID == GP2SPL )
		grid2spec(sptrans, gridID1, array1, gridID2, array2);	      
	      else if ( operatorID == SP2GP || operatorID == SP2GPL )
		spec2grid(sptrans, gridID1, array1, gridID2, array2);
	      else if ( operatorID == SP2SP )
		spec2spec(gridID1, array1, gridID2, array2);
	      else if ( operatorID == SPCUT )
		speccut(gridID1, array1, array2, waves);
	      */

	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, array2, nmiss);  
	    }   
	  else
	    {
	      streamDefRecord(streamID2, varID, levelID);
	      if ( lcopy )
		{
		  streamCopyRecord(streamID2, streamID1);
		}
	      else
		{
		  streamReadRecord(streamID1, array1, &nmiss);
		  streamWriteRecord(streamID2, array1, nmiss);
		}
	    }    
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( array2 ) free(array2);
  if ( array1 ) free(array1);

  listDelete(ilist);

  sptrans_delete(sptrans);

  cdoFinish();

  return (0);
}

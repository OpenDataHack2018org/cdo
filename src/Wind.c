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


void *Wind(void *argument)
{
  static char func[] = "Wind";
  int UV2DV, DV2UV;
  int operatorID;
  int streamID1, streamID2;
  int nrecs, nvars;
  int tsID, recID, varID, levelID;
  int nlev, gridsize;
  int index, ngrids;
  int vlistID1, vlistID2;
  int gridIDsp = -1, gridIDgp = -1;
  int gridID1 = -1, gridID2 = -1;
  int gridID;
  int nmiss;
  int lcopy = FALSE;
  int taxisID1, taxisID2;
  int nlon, nlat, trunc;
  int code;
  int varID1 = -1, varID2 = -1;
  int offset;
  SPTRANS *sptrans = NULL;
  DVTRANS *dvtrans = NULL;
  char varname[128];
  double *array1 = NULL, *array2 = NULL;
  double *ivar1 = NULL, *ivar2 = NULL, *ovar1 = NULL, *ovar2 = NULL;

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

  /* find variables */
  nvars = vlistNvars(vlistID2);
  for ( varID = 0; varID < nvars; varID++ )
    {
      if ( operatorID == UV2DV ) 
	{
	  /* search for u and v wind */
	}
      else
	{
	  /* search for divergence and vorticity */
	  code = vlistInqVarCode(vlistID1, varID);
	  if ( code <= 0 )
	    {
	      vlistInqVarName(vlistID1, varID, varname);

	      strtolower(varname);

	      if      ( strcmp(varname, "sd")  == 0 ) code = 155;
	      else if ( strcmp(varname, "svo") == 0 ) code = 138;
	    }

	  if      ( code == 155 ) varID1 = varID;
	  else if ( code == 138 ) varID2 = varID;
	}
    }

  ngrids = vlistNgrids(vlistID1);

  /* find first spectral grid */
  for ( index = 0; index < ngrids; index++ )
    {
      gridID = vlistGrid(vlistID1, index);
      if ( gridInqType(gridID) == GRID_SPECTRAL )
	{
	  gridIDsp = gridID;
	  break;
	}
    }

  /* find first gaussian grid */
  for ( index = 0; index < ngrids; index++ )
    {
      gridID = vlistGrid(vlistID1, index);
      if ( gridInqType(gridID) == GRID_GAUSSIAN )
	{
	  gridIDgp = gridID;
	  break;
	}
    }

  /* define output grid */
  if ( operatorID == UV2DV )
    {
      gridID1 = gridIDgp;

      if ( gridID1 != -1 )
	{
	  trunc = nlat2trunc(gridInqYsize(gridID1));

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
  else
    {   
      if ( varID1 == -1 ) cdoWarning("Divergence not found!");
      if ( varID2 == -1 ) cdoWarning("Vorticity not found!");

      gridID1 = vlistInqVarGrid(vlistID1, varID2);

      if ( gridInqType(gridID1) != GRID_SPECTRAL )
	cdoAbort("Vorticity is not on spectral grid!");

      if ( gridID1 != vlistInqVarGrid(vlistID1, varID1) )
	cdoAbort("Divergence and vorticity must have the same grid represention!");

      if ( gridIDgp != -1 )
	{
	  trunc = nlat2trunc(gridInqYsize(gridIDgp));
	      
	  if ( gridInqTrunc(gridIDsp) != trunc ) gridIDgp = -1;
	}

      if ( gridIDgp == -1 )
	{
	  char gridname[20];

	  sprintf(gridname, "t%dgrid", gridInqTrunc(gridID1));

	  gridIDgp = gridFromName(gridname);
	}

      gridID2 = gridIDgp;

      vlistChangeVarGrid(vlistID2, varID1, gridID2);
      vlistChangeVarGrid(vlistID2, varID2, gridID2);
      vlistDefVarCode(vlistID2, varID1, 131);
      vlistDefVarCode(vlistID2, varID2, 132);
      /* define varname aso. !!! */

      trunc = gridInqTrunc(gridID1);
      nlon  = gridInqXsize(gridID2);
      nlat  = gridInqYsize(gridID2);
      
      sptrans = sptrans_new(nlon, nlat, trunc);
      dvtrans = dvtrans_new(trunc);
    }

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

  nlev     = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID1));

  gridsize = gridInqSize(gridID1);
  ivar1 = (double *) malloc(nlev*gridsize*sizeof(double));
  ivar2 = (double *) malloc(nlev*gridsize*sizeof(double));

  gridsize = gridInqSize(gridID2);
  ovar1 = (double *) malloc(nlev*gridsize*sizeof(double));
  ovar2 = (double *) malloc(nlev*gridsize*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  if ( varID == varID1 || varID == varID2 )
	    {
	      streamReadRecord(streamID1, array1, &nmiss);
	      if ( nmiss ) cdoAbort("missing values unsupported for spectral data!");

	      gridsize = gridInqSize(gridID1);
	      offset = gridsize*levelID;

	      if      ( varID == varID1 )
		memcpy(ivar1+offset, array1, gridsize*sizeof(double));
	      else if ( varID == varID2 )
	        memcpy(ivar2+offset, array1, gridsize*sizeof(double));
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
      /*
	if ( operatorID == UV2DV )
	uv2dv(sptrans, gridID1, ivar1, ivar2, gridID2, ovar1, ovar2);	      
	else if ( operatorID == DV2UV )
	dv2uv(sptrans, dvtrans, gridID1, ivar1, ivar2, gridID2, ovar1, ovar2);
      */
      if ( operatorID == DV2UV )
	trans_dv2uv(sptrans, dvtrans, nlev, gridID1, ivar1, ivar2, gridID2, ovar1, ovar2);

      gridsize = gridInqSize(gridID2);
      for ( levelID = 0; levelID < nlev; levelID++ )
	{
	  offset = gridsize*levelID;
	  streamDefRecord(streamID2, varID1, levelID);
	  streamWriteRecord(streamID2, ovar1+offset, 0);
	}
      for ( levelID = 0; levelID < nlev; levelID++ )
	{
	  offset = gridsize*levelID;
	  streamDefRecord(streamID2, varID2, levelID);
	  streamWriteRecord(streamID2, ovar2+offset, 0);
	}

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( array2 ) free(array2);
  if ( array1 ) free(array1);

  if ( ivar1 ) free(ivar1);
  if ( ivar2 ) free(ivar2);
  if ( ovar1 ) free(ovar1);
  if ( ovar2 ) free(ovar2);

  sptrans_delete(sptrans);
  if ( operatorID == DV2UV ) dvtrans_delete(dvtrans);

  cdoFinish();

  return (0);
}

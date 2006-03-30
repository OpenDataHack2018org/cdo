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

/*
   This module contains the following operators:

      Fldstat    fldmin          Field minimum
      Fldstat    fldmax          Field maximum
      Fldstat    fldsum          Field sum
      Fldstat    fldmean         Field mean
      Fldstat    fldavg          Field average
      Fldstat    fldstd          Field standard deviation
      Fldstat    fldvar          Field variance
*/


#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "functs.h"


void *Fldstat(void *argument)
{
  static char func[] = "Fldstat";
  int operatorID;
  int operfunc;
  int streamID1, streamID2;
  int vlistID1, vlistID2;
  int gridID2, lastgrid = -1;
  int index, ngrids;
  int recID, nrecs;
  int tsID, varID, levelID;
  int lim;
  int needWeights = FALSE;
  int nmiss;
  int ndiffgrids;
  double slon, slat;
  double sglval;
  FIELD field;
  int taxisID1, taxisID2;

  cdoInitialize(argument);

  cdoOperatorAdd("fldmin",  func_min,  0, NULL);
  cdoOperatorAdd("fldmax",  func_max,  0, NULL);
  cdoOperatorAdd("fldsum",  func_sum,  0, NULL);
  cdoOperatorAdd("fldmean", func_mean, 0, NULL);
  cdoOperatorAdd("fldavg",  func_avg,  0, NULL);
  cdoOperatorAdd("fldvar",  func_var,  0, NULL);
  cdoOperatorAdd("fldstd",  func_std,  0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

  if ( operfunc == func_mean || operfunc == func_avg ||
       operfunc == func_var  || operfunc == func_std )
    needWeights = TRUE;

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  slon = 0;
  slat = 0;
  gridID2 = gridNew(GRID_LONLAT, 1);
  gridDefXsize(gridID2, 1);
  gridDefYsize(gridID2, 1);
  gridDefXvals(gridID2, &slon);
  gridDefYvals(gridID2, &slat);

  ngrids = vlistNgrids(vlistID1);
  ndiffgrids = 0;
  for ( index = 1; index < ngrids; index++ )
    if ( vlistGrid(vlistID1, 0) != vlistGrid(vlistID1, index))
      ndiffgrids++;

  if ( ndiffgrids > 0 ) cdoAbort("Too many different grids!");

  for ( index = 0; index < ngrids; index++ )
    vlistChangeGridIndex(vlistID2, index, gridID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  lim = vlistGridsizeMax(vlistID1);
  field.ptr    = (double *) malloc(lim*sizeof(double));
  field.weight = NULL;
  if ( needWeights )
    field.weight = (double *) malloc(lim*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, field.ptr, &field.nmiss);

	  field.grid = vlistInqVarGrid(vlistID1, varID);
	  field.size = gridInqSize(field.grid);
	  if ( field.grid != lastgrid )
	    {
	      lastgrid = field.grid;
	      if ( needWeights ) gridWeights(field.grid, field.weight);
	    }
	  field.missval = vlistInqVarMissval(vlistID1, varID);

	  sglval = fldfun(field, operfunc);

	  if ( DBL_IS_EQUAL(sglval, field.missval) )
	    nmiss = 1;
	  else
	    nmiss = 0;

	  streamDefRecord(streamID2, varID,  levelID);
	  streamWriteRecord(streamID2, &sglval, nmiss);
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( field.ptr )    free(field.ptr);
  if ( field.weight ) free(field.weight);

  cdoFinish();

  return (0);
}

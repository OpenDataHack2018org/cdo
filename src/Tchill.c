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
   This module contains the following operators:

      Tchill     tchill          Compute the windchill temperature (K)
*/


#include <string.h>
#include <math.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

#define TO_KELVIN(x) ((x) + 273.15)


static double windchillTemperature(double t, double v, double missval)
{
  static const double tmax = TO_KELVIN(33.0);
  
  return v < 0.0 || t > tmax ? missval : tmax + (0.478 + 0.237 * sqrt(v) - 0.0124 * v) * (t - tmax);
}


static void farexpr(FIELD *field1, FIELD field2, double (*expression)(double, double, double))
{
  static char func[] = "farexpr";
  int   i, len;
  const int     grid1    = field1->grid;
  const int     nmiss1   = field1->nmiss;
  const double  missval1 = field1->missval;
  double       *array1   = field1->ptr;
  const int     grid2    = field2.grid;
  const int     nmiss2   = field2.nmiss;
  const double  missval2 = field2.missval;
  const double *array2   = field2.ptr;

  len = gridInqSize(grid1);

  if ( len != gridInqSize(grid2) )
    cdoAbort("Fields have different gridsize (%s)", func);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ )
        if ( DBL_IS_EQUAL(array1[i], missval1) || DBL_IS_EQUAL(array2[i], missval2) )  
	  array1[i] = missval1;
	else
	  array1[i] = expression(array1[i], array2[i], missval1);
    }
  else
    {
      for ( i = 0; i < len; i++ )
        array1[i] = expression(array1[i], array2[i], missval1);  
    }

  field1->nmiss = 0;
  for ( i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
}

   
void *Tchill(void *argument)
{
  static char func[] = "Tchill";
  enum {FILL_NONE, FILL_REC, FILL_TS};
  int streamIDm, streamIDs, streamID1, streamID2, streamID3;
  int gridsize;
  int nrecs, nrecs2, nvars = 0, nlev, recID;
  int tsID;
  int varID, levelID;
  int offset;
  int vlistIDm, vlistIDs, vlistID1, vlistID2, vlistID3;
  int taxisIDm, taxisID1, taxisID2, taxisID3;
  int filltype = FILL_NONE;
  FIELD *fieldm, *fields, fieldtmp, field1, field2;
  int **varnmiss = NULL;
  double **vardata = NULL;

  cdoInitialize(argument);
  cdoOperatorAdd("tchill", 0, 0, NULL);

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));
  streamID2 = streamOpenRead(cdoStreamName(1));
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamIDm = streamID1;
  streamIDs = streamID2;
  fieldm = &field1;
  fields = &field2;

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = streamInqVlist(streamID2);
  vlistIDm = vlistID1;
  vlistIDs = vlistID2;

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = vlistInqTaxis(vlistID2);
  taxisIDm = taxisID1;

  if ( vlistNrecs(vlistID1) != 1 && vlistNrecs(vlistID2) == 1 )
    {
      filltype = FILL_TS;
      cdoPrint("Filling up stream2 >%s< by copying the first record.", cdoStreamName(1));
    }
  else if ( vlistNrecs(vlistID1) == 1 && vlistNrecs(vlistID2) != 1 )
    {
      filltype = FILL_TS;
      cdoPrint("Filling up stream1 >%s< by copying the first record.", cdoStreamName(0));
      streamIDm = streamID2;
      streamIDs = streamID1;
      vlistIDm = vlistID2;
      vlistIDs = vlistID1;
      taxisIDm = taxisID2;
      fieldm = &field2;
      fields = &field1;
    }

  if ( filltype == FILL_NONE )
    vlistCompare(vlistID1, vlistID2, func_sft);
  
  nospec(vlistID1);

  gridsize = vlistGridsizeMax(vlistIDm);

  field1.ptr = (double *) malloc(gridsize*sizeof(double));
  field2.ptr = (double *) malloc(gridsize*sizeof(double));
  fieldtmp.ptr = NULL;
  fieldtmp.nmiss = 0;
  if ( filltype == FILL_TS )
    fieldtmp.ptr = (double *) malloc(gridsize*sizeof(double));

  if ( cdoVerbose )
    cdoPrint("Number of timesteps: file1 %d, file2 %d",
	     vlistNtsteps(vlistID1), vlistNtsteps(vlistID2));

  if ( filltype == FILL_NONE )
    {
      if ( vlistNtsteps(vlistID1) != 1 && vlistNtsteps(vlistID1) != 0 &&
	  (vlistNtsteps(vlistID2) == 1 || vlistNtsteps(vlistID2) == 0) )
	{
	  filltype = FILL_REC;
	  cdoPrint("Filling up stream2 >%s< by copying the first timestep.", cdoStreamName(1));
	}
      else if ( (vlistNtsteps(vlistID1) == 1 || vlistNtsteps(vlistID1) == 0) &&
		 vlistNtsteps(vlistID2) != 1 && vlistNtsteps(vlistID2) != 0 )
	{
	  filltype = FILL_REC;
	  cdoPrint("Filling up stream1 >%s< by copying the first timestep.", cdoStreamName(0));
	  streamIDm = streamID2;
          streamIDs = streamID1;
	  vlistIDm = vlistID2;
	  vlistIDs = vlistID1;
	  taxisIDm = taxisID2;
	  fieldm = &field2;
	  fields = &field1;
	}

      if ( filltype == FILL_REC )
	{
	  nvars  = vlistNvars(vlistIDs);
	  vardata  = (double **) malloc(nvars*sizeof(double *));
	  varnmiss = (int **) malloc(nvars*sizeof(int *));
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      gridsize = gridInqSize(vlistInqVarGrid(vlistIDs, varID));
	      nlev     = zaxisInqSize(vlistInqVarZaxis(vlistIDs, varID));
	      vardata[varID]  = (double *) malloc(nlev*gridsize*sizeof(double));
	      varnmiss[varID] = (int *) malloc(nlev*sizeof(int));
	    }
	}
    }

  vlistID3 = vlistDuplicate(vlistIDm);

  taxisID3 = taxisDuplicate(taxisIDm);
  vlistDefTaxis(vlistID3, taxisID3);

  streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());
  if ( streamID3 < 0 ) cdiError(streamID3, "Open failed on %s", cdoStreamName(2));

  streamDefVlist(streamID3, vlistID3);

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamIDm, tsID)) )
    {
      if ( tsID == 0 || filltype == FILL_NONE )
	{
	  nrecs2 = streamInqTimestep(streamIDs, tsID);
	  if ( nrecs2 == 0 )
	    cdoAbort("Input streams have different number of timesteps!");
	}

      taxisCopyTimestep(taxisID3, taxisIDm);

      streamDefTimestep(streamID3, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamIDm, &varID, &levelID);
	  streamReadRecord(streamIDm, fieldm->ptr, &fieldm->nmiss);

	  if ( tsID == 0 || filltype == FILL_NONE )
	    {
	      if ( recID == 0 || filltype != FILL_TS )
		{
		  streamInqRecord(streamIDs, &varID, &levelID);
		  streamReadRecord(streamIDs, fields->ptr, &fields->nmiss);
		}

	      if ( filltype == FILL_REC )
		{
		  gridsize = gridInqSize(vlistInqVarGrid(vlistIDs, varID));
		  offset   = gridsize*levelID;
		  memcpy(vardata[varID]+offset, fields->ptr, gridsize*sizeof(double));
		  varnmiss[varID][levelID] = fields->nmiss;
		}
	      else if ( recID == 0 && filltype == FILL_TS )
		{
		  gridsize = gridInqSize(vlistInqVarGrid(vlistIDs, 0));
		  memcpy(fieldtmp.ptr, fields->ptr, gridsize*sizeof(double));
		  fieldtmp.nmiss = fields->nmiss;
		}
	    }
	  else if ( filltype == FILL_REC )
	    {
	      gridsize = gridInqSize(vlistInqVarGrid(vlistIDs, varID));
	      offset   = gridsize*levelID;
	      memcpy(fields->ptr, vardata[varID]+offset, gridsize*sizeof(double));
	      fields->nmiss = varnmiss[varID][levelID];
	    }

	  fieldm->grid    = vlistInqVarGrid(vlistIDm, varID);
	  fieldm->missval = vlistInqVarMissval(vlistIDm, varID);

	  if ( filltype == FILL_TS )
	    {
	      gridsize = gridInqSize(vlistInqVarGrid(vlistIDs, 0));
	      memcpy(fields->ptr, fieldtmp.ptr, gridsize*sizeof(double));
	      fields->nmiss   = fieldtmp.nmiss;
	      fields->grid    = vlistInqVarGrid(vlistIDs, 0);
	      fields->missval = vlistInqVarMissval(vlistIDs, 0);
	    }
	  else
	    {
	      fields->grid    = vlistInqVarGrid(vlistIDs, varID);
	      fields->missval = vlistInqVarMissval(vlistIDs, varID);
	    }

	  farexpr(&field1, field2, windchillTemperature);

	  streamDefRecord(streamID3, varID, levelID);
	  streamWriteRecord(streamID3, field1.ptr, field1.nmiss);
	}

      tsID++;
    }

  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);

  if ( vardata )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  free(vardata[varID]);
	  free(varnmiss[varID]);
	}

      free(vardata);
      free(varnmiss);
    }

  if ( field1.ptr ) free(field1.ptr);
  if ( field2.ptr ) free(field2.ptr);
  if ( filltype == FILL_TS )
  if ( fieldtmp.ptr ) free(fieldtmp.ptr);

  cdoFinish();

  return (0);
}

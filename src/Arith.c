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

      Arith      add             Add two fields
      Arith      sub             Subtract two fields
      Arith      mul             Multiply two fields
      Arith      div             Divide two fields
      Arith      min             Minimum of two fields
      Arith      max             Maximum of two fields
*/


#include <string.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Arith(void *argument)
{
  static char func[] = "Arith";
  int operatorID;
  int operfunc;
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

  cdoOperatorAdd("add", func_add, 0, NULL);
  cdoOperatorAdd("sub", func_sub, 0, NULL);
  cdoOperatorAdd("mul", func_mul, 0, NULL);
  cdoOperatorAdd("div", func_div, 0, NULL);
  cdoOperatorAdd("min", func_min, 0, NULL);
  cdoOperatorAdd("max", func_max, 0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

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
  
  if ( operfunc == func_mul || operfunc == func_div )
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

	  farfun(&field1, field2, operfunc);

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

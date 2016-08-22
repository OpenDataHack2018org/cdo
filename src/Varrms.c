/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2016 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
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

*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"


void *Varrms(void *argument)
{
  int lastgrid = -1;
  int wstatus = FALSE;
  int oldcode = 0;
  int nrecs;
  int gridsize;
  int nmiss;
  int varID, levelID;
  long offset;
  double *single;
  double sglval;

  cdoInitialize(argument);

  bool needWeights = true;

  int streamID1 = streamOpenRead(cdoStreamName(0));
  int streamID2 = streamOpenRead(cdoStreamName(1));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = streamInqVlist(streamID2);

  double slon = 0;
  double slat = 0;
  int gridID3 = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID3, 1);
  gridDefYsize(gridID3, 1);
  gridDefXvals(gridID3, &slon);
  gridDefYvals(gridID3, &slat);

  vlistClearFlag(vlistID1);
  int nvars    = vlistNvars(vlistID1);
  for ( varID = 0; varID < nvars; varID++ )
    vlistDefFlag(vlistID1, varID, 0, TRUE);

  int vlistID3 = vlistCreate();
  vlistCopyFlag(vlistID3, vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID3);

  int ngrids = vlistNgrids(vlistID1);
  int index = 0;
  int gridID1 = vlistGrid(vlistID1, index);
  
  if ( needWeights &&
       gridInqType(gridID1) != GRID_LONLAT &&
       gridInqType(gridID1) != GRID_GAUSSIAN )
    cdoAbort("Unsupported gridtype: %s", gridNamePtr(gridInqType(gridID1)));

  vlistChangeGridIndex(vlistID3, index, gridID3);
  if ( ngrids > 1 ) cdoAbort("Too many different grids!");

  int streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());

  streamDefVlist(streamID3, vlistID3);

  double **vardata1 = (double**) Malloc(nvars*sizeof(double*));
  double **vardata2 = (double**) Malloc(nvars*sizeof(double*));

  for ( varID = 0; varID < nvars; varID++ )
    {
      int gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
      int nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      vardata1[varID] = (double*) Malloc(gridsize*nlevel*sizeof(double));
      vardata2[varID] = (double*) Malloc(gridsize*nlevel*sizeof(double));
    }

  field_t field1, field2, field3;
  field_init(&field1);
  field_init(&field2);
  field_init(&field2);

  int lim = vlistGridsizeMax(vlistID1);
  field1.weight = NULL;
  if ( needWeights )
    field1.weight = (double*) Malloc(lim*sizeof(double));

  field2.weight = NULL;

  field3.ptr  = &sglval;
  field3.grid = gridID3;

  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      nrecs = streamInqTimestep(streamID2, tsID);

      taxisCopyTimestep(taxisID3, taxisID1);

      streamDefTimestep(streamID3, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  offset   = gridsize*levelID;
	  single   = vardata1[varID] + offset;

	  streamReadRecord(streamID1, single, &nmiss);
	  if ( nmiss ) cdoAbort("Missing values unsupported for this operator!");

	  streamInqRecord(streamID2, &varID, &levelID);

	  single   = vardata2[varID] + offset;

	  streamReadRecord(streamID2, single, &nmiss);
	  if ( nmiss ) cdoAbort("this operator does not work with missing values!");
	}

      for ( varID = 0; varID < nvars; varID++ )
	{
	  field1.ptr = vardata1[varID];
	  field2.ptr = vardata2[varID];

	  field1.zaxis   = vlistInqVarZaxis(vlistID1, varID);
	  field1.grid    = vlistInqVarGrid(vlistID1, varID);
	  field2.grid    = vlistInqVarGrid(vlistID2, varID);
	  if ( needWeights && field1.grid != lastgrid )
	    {
	      lastgrid = field1.grid;
	      wstatus = gridWeights(field1.grid, field1.weight);
	    }
	  int code = vlistInqVarCode(vlistID1, varID);
	  if ( wstatus != 0 && tsID == 0 && code != oldcode )
	    cdoWarning("Using constant area weights for code %d!", oldcode=code);

	  field1.missval = vlistInqVarMissval(vlistID1, varID);
	  field2.missval = vlistInqVarMissval(vlistID1, varID);
	  field3.missval = vlistInqVarMissval(vlistID1, varID);

	  varrms(field1, field2, &field3);

	  streamDefRecord(streamID3, varID, 0);
	  streamWriteRecord(streamID3, &sglval, field3.nmiss);
	}
      tsID++;
    }

  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);

  vlistDestroy(vlistID3);

  if ( field1.weight ) Free(field1.weight);

  for ( varID = 0; varID < nvars; varID++ )
    {
      Free(vardata1[varID]);
      Free(vardata2[varID]);
    }

  Free(vardata1);
  Free(vardata2);

  cdoFinish();

  return 0;
}

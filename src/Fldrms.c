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

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "functs.h"


void *Fldrms(void *argument)
{
  static char func[] = "Fldrms";
  int streamID1, streamID2, streamID3;
  int vlistID1, vlistID2, vlistID3;
  int gridID1, gridID2, gridID3, lastgrid = -1;
  int index, ngrids;
  int recID, nrecs;
  int tsID, varID, levelID;
  int lim;
  int needWeights = FALSE;
  double slon, slat;
  double sglval;
  FIELD field1, field2, field3;
  int taxisID1, taxisID3;

  cdoInitialize(argument);

  needWeights = TRUE;

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  streamID2 = streamOpenRead(cdoStreamName(1));
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = streamInqVlist(streamID2);
  vlistID3 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID3);

  slon = 0;
  slat = 0;
  gridID3 = gridNew(GRID_LONLAT, 1);
  gridDefXsize(gridID3, 1);
  gridDefYsize(gridID3, 1);
  gridDefXvals(gridID3, &slon);
  gridDefYvals(gridID3, &slat);

  ngrids = vlistNgrids(vlistID1);
  index = 0;
  gridID1 = vlistGrid(vlistID1, index);
  gridID2 = vlistGrid(vlistID2, index);

  if ( gridInqSize(gridID1) != gridInqSize(gridID2) )
    cdoAbort("Fields have different grid size");
  
  if ( needWeights &&
       gridInqType(gridID1) != GRID_LONLAT &&
       gridInqType(gridID1) != GRID_GAUSSIAN )
    cdoAbort("Unsupported gridtype: %s", gridNamePtr(gridInqType(gridID1)));

  vlistChangeGridIndex(vlistID3, index, gridID3);
  if ( ngrids > 1 ) cdoAbort("Too many different grids!");

  streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());
  if ( streamID3 < 0 ) cdiError(streamID3, "Open failed on %s", cdoStreamName(2));

  streamDefVlist(streamID3, vlistID3);

  lim = vlistGridsizeMax(vlistID1);
  field1.ptr    = (double *) malloc(lim*sizeof(double));
  field1.weight = NULL;
  if ( needWeights )
    field1.weight = (double *) malloc(lim*sizeof(double));

  field2.ptr    = (double *) malloc(lim*sizeof(double));
  field2.weight = NULL;

  field3.ptr  = &sglval;
  field3.grid = gridID3;

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      nrecs = streamInqTimestep(streamID2, tsID);

      taxisCopyTimestep(taxisID3, taxisID1);

      streamDefTimestep(streamID3, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, field1.ptr, &field1.nmiss);
	  streamInqRecord(streamID2, &varID, &levelID);
	  streamReadRecord(streamID2, field2.ptr, &field2.nmiss);

	  field1.grid    = vlistInqVarGrid(vlistID1, varID);
	  field2.grid    = vlistInqVarGrid(vlistID2, varID);
	  if ( field1.grid != lastgrid )
	    {
	      lastgrid = field1.grid;
	      if ( needWeights ) gridWeights(field1.grid, field1.weight);
	    }
	  field1.missval = vlistInqVarMissval(vlistID1, varID);
	  field2.missval = vlistInqVarMissval(vlistID1, varID);
	  field3.missval = vlistInqVarMissval(vlistID1, varID);

	  fldrms(field1, field2, &field3);

	  streamDefRecord(streamID3, varID,  levelID);
	  streamWriteRecord(streamID3, &sglval, field3.nmiss);
	}
      tsID++;
    }

  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);

  if ( field1.ptr )    free(field1.ptr);
  if ( field1.weight ) free(field1.weight);
  if ( field2.ptr )    free(field2.ptr);

  cdoFinish();

  return (0);
}

/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2008 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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


#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "functs.h"


void *Gengriddes(void *argument)
{
  static char func[] = "Gengriddes";
  int streamID1, streamID2, streamID3;
  int vlistID1, vlistID2, vlistID3;
  int gridID1, gridID2, gridID3, lastgrid = -1;
  int wstatus = FALSE;
  int code = 0, oldcode = 0;
  int index, ngrids;
  int recID, nrecs;
  int tsID, varID, levelID;
  int lim;
  int ndiffgrids;
  int gridsize;
  double *array1, *array2, *array3;
  int taxisID1;

  cdoInitialize(argument);

  needWeights = TRUE;

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  streamID2 = streamOpenRead(cdoStreamName(1));
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = streamInqVlist(streamID2);

  taxisID1 = vlistInqTaxis(vlistID1);

  index = 0;
  gridID1 = vlistGrid(vlistID1, index);
  gridID2 = vlistGrid(vlistID2, index);

  if ( gridInqSize(gridID1) != gridInqSize(gridID2) )
    cdoAbort("Arrays have different grid size!");

  gridsize = gridInqSize(gridID1);
  xsize = gridInqXsize(gridID1);
  ysize = gridInqYsize(gridID1);

  array1 = (double *) malloc(gridsize*sizeof(double));
  array2 = (double *) malloc(gridsize*sizeof(double));
  array3 = (double *) malloc(gridsize*sizeof(double));

  tsID = 0;
  nrecs = streamInqTimestep(streamID1, tsID);
  nrecs = streamInqTimestep(streamID2, tsID);

  streamInqRecord(streamID1, &varID, &levelID);
  streamReadRecord(streamID1, array1, &nmiss1);
  streamInqRecord(streamID2, &varID, &levelID);
  streamReadRecord(streamID2, array2, &nmiss2);

  missval = vlistInqVarMissval(vlistID1, varID);

  streamClose(streamID2);
  streamClose(streamID1);

  if ( nmiss1 || nmiss2 ) cdoAbort("Missing values unsupported!");

  gridID3 = gridCreate(GRID_CURVILINEAR, gridsize);
  gridDefXsize(gridID3, xsize);
  gridDefYsize(gridID3, ysize);
  gridDefXvals(gridID3, array1);
  gridDefYvals(gridID3, array2);


  streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());
  if ( streamID3 < 0 ) cdiError(streamID3, "Open failed on %s", cdoStreamName(2));

  streamDefVlist(streamID3, vlistID3);

  tsID = 0;
  streamDefTimestep(streamID3, tsID);

  streamDefRecord(streamID3, 0, 0);
  streamWriteRecord(streamID3, array3, gridsize);

  streamClose(streamID3);

  if ( array1 ) free(array1);
  if ( array2 ) free(array2);
  if ( array3 ) free(array3);

  cdoFinish();

  return (0);
}

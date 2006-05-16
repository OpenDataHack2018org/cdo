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

      Intgrid    interpolate     Interpolate
*/


#include <string.h>
#include <math.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "interpol.h"


void *Intgrid(void *argument)
{
  static char func[] = "Intgrid";
  int INTGRID, INTPOINT, INTERPOLATE;
  int operatorID;
  int streamID1, streamID2;
  int nrecs, ngrids;
  int index;
  int tsID, recID, varID, levelID;
  int gridsize;
  int vlistID1, vlistID2;
  int gridID1 = -1, gridID2;
  int nmiss;
  double missval;
  double slon, slat;
  double *array1 = NULL, *array2 = NULL;
  FIELD field1, field2;
  int taxisID1, taxisID2;

  cdoInitialize(argument);

  INTGRID     = cdoOperatorAdd("intgridbil", 0, 0, NULL);
  INTPOINT    = cdoOperatorAdd("intpoint", 0, 0, NULL);
  INTERPOLATE = cdoOperatorAdd("interpolate", 0, 0, NULL);

  operatorID = cdoOperatorID();

  if ( operatorID == INTGRID || operatorID == INTERPOLATE )
    {
      operatorInputArg("grid description file or name");
      gridID2 = cdoDefineGrid(operatorArgv()[0]);
    }
  else
    {
      operatorInputArg("longitude and latitude");
      operatorCheckArgc(2);
      slon = atof(operatorArgv()[0]);
      slat = atof(operatorArgv()[1]);
      gridID2 = gridCreate(GRID_LONLAT, 1);
      gridDefXsize(gridID2, 1);
      gridDefYsize(gridID2, 1);
      gridDefXvals(gridID2, &slon);
      gridDefYvals(gridID2, &slat);
    }

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  ngrids = vlistNgrids(vlistID1);
  for ( index = 0; index < ngrids; index++ )
    {
      gridID1 = vlistGrid(vlistID1, index);

      if ( gridInqType(gridID1) != GRID_LONLAT &&
	   gridInqType(gridID1) != GRID_GAUSSIAN )
	cdoAbort("Interpolation of %s data unsupported!", gridNamePtr(gridInqType(gridID1)) );

      vlistChangeGridIndex(vlistID2, index, gridID2);
    }

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);
  array1   = (double *) malloc(gridsize*sizeof(double));

  gridsize = gridInqSize(gridID2);
  array2   = (double *) malloc(gridsize*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array1, &nmiss);

	  gridID1 = vlistInqVarGrid(vlistID1, varID);
	  missval = vlistInqVarMissval(vlistID1, varID);

	  field1.grid    = gridID1;
	  field1.nmiss   = nmiss;
	  field1.missval = missval;
	  field1.ptr     = array1;
	  field2.grid    = gridID2;
	  field2.ptr     = array2;
	  field2.nmiss   = 0;

	  if ( operatorID == INTGRID || operatorID == INTPOINT )
	    intgrid(&field1, &field2);
	  else if ( operatorID == INTERPOLATE )
	    interpolate(&field1, &field2);

	  nmiss = field2.nmiss;

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, array2, nmiss);
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( array2 ) free(array2);
  if ( array1 ) free(array1);

  cdoFinish();

  return (0);
}

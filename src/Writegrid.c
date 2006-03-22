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

*/


#include <string.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

/*
@BeginDoc

@BeginModule

@Name      = Setgrid
@Title     = Set grid type
@Section   = Information
@Class     = Information
@Arguments = ifile ofile
@Operators = writegrid

@EndModule


@BeginOperator_writegrid

@Title     = Write grid

@BeginDescription
Write the grid information to a file (SCRIP netCDF).
@EndDescription

@EndOperator

@EndDoc
*/

void *Writegrid(void *argument)
{
  static char func[] = "Writegrid";
  int streamID;
  int vlistID1;
  int gridID1;
  int i;
  int gridtype, gridsize;
  int varID, levelID;
  int nrecs;
  int nmiss;
  int *imask;
  double missval;
  double *array;

  cdoInitialize(argument);

  streamID = streamOpenRead(cdoStreamName(0));
  if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID);

  nrecs = streamInqTimestep(streamID, 0);

  streamInqRecord(streamID, &varID, &levelID);

  gridID1  = vlistInqVarGrid(vlistID1, varID);
  gridtype = gridInqType(gridID1);
  gridsize = gridInqSize(gridID1);

  if ( gridtype != GRID_CURVILINEAR && gridtype != GRID_CELL )
    {
      if ( gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT )
	gridID1 = gridToCurvilinear(gridID1);
      else if ( gridtype == GRID_GME )
	gridID1 = gridToCell(gridID1);
      else
	cdoAbort("%s grid unsupported!", gridNamePtr(gridtype));
    }

  array = (double *) malloc(gridsize*sizeof(double));
  imask = (int *) malloc(gridsize*sizeof(int));
  streamReadRecord(streamID, array, &nmiss);

  missval = vlistInqVarMissval(vlistID1, varID);
  for ( i = 0; i < gridsize; i++ )
    {
      if ( DBL_IS_EQUAL(array[i], missval) )
	imask[i] = 0;
      else
	imask[i] = 1;
    }
      
  writeNCgrid(cdoStreamName(1), gridID1, imask);

  streamClose(streamID);

  free(array);
  free(imask);

  cdoFinish();

  return (0);
}

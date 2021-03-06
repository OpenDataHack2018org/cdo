/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2018 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
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

#include "cdo_int.h"
#include "pstream_int.h"

void *
Gengrid(void *process)
{
  int varID, levelID;
  size_t nmiss1, nmiss2;
  double missval = 0;

  cdoInitialize(process);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));
  int streamID2 = cdoStreamOpenRead(cdoStreamName(1));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = cdoStreamInqVlist(streamID2);

  int gridID1 = vlistGrid(vlistID1, 0);
  int gridID2 = vlistGrid(vlistID2, 0);

  if (gridInqSize(gridID1) != gridInqSize(gridID2)) cdoAbort("Arrays have different grid size!");

  size_t gridsize = gridInqSize(gridID1);
  size_t xsize = gridInqXsize(gridID1);
  size_t ysize = gridInqYsize(gridID1);

  std::vector<double> array1(gridsize);
  std::vector<double> array2(gridsize);
  std::vector<double> array3(gridsize);

  cdoStreamInqTimestep(streamID1, 0);
  cdoStreamInqTimestep(streamID2, 0);

  pstreamInqRecord(streamID1, &varID, &levelID);
  pstreamReadRecord(streamID1, array1.data(), &nmiss1);
  pstreamInqRecord(streamID2, &varID, &levelID);
  pstreamReadRecord(streamID2, array2.data(), &nmiss2);

  int datatype = vlistInqVarDatatype(vlistID1, 0);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if (nmiss1 || nmiss2) cdoAbort("Missing values unsupported!");

  int gridID3 = gridCreate(GRID_CURVILINEAR, gridsize);

  if (cdoVerbose) cdoPrint("xsize %zu  ysize %zu", xsize, ysize);
  if (xsize * ysize != gridsize) cdoAbort("xsize*ysize != gridsize");

  gridDefXsize(gridID3, xsize);
  gridDefYsize(gridID3, ysize);
  gridDefXvals(gridID3, array1.data());
  gridDefYvals(gridID3, array2.data());

  if (datatype == CDI_DATATYPE_FLT64)
    gridDefDatatype(gridID3, CDI_DATATYPE_FLT64);
  else
    gridDefDatatype(gridID3, CDI_DATATYPE_FLT32);

  double xminval, xmaxval;
  double yminval, ymaxval;
  arrayMinMax(gridsize, array1.data(), &xminval, &xmaxval);
  arrayMinMax(gridsize, array2.data(), &yminval, &ymaxval);

  if (cdoVerbose) cdoPrint("xminval = %g, xmaxval = %g, yminval = %g, ymaxval = %g", xminval, xmaxval, yminval, ymaxval);

  /* check units */
  if (xminval > -4 && xmaxval < 8 && yminval > -2 && ymaxval < 2)
    {
      gridDefXunits(gridID3, "radians");
      gridDefYunits(gridID3, "radians");
    }
  else if (xminval > -181 && xmaxval < 361 && yminval > -91 && ymaxval < 91)
    {
      /* default is degrees */
    }
  else
    {
      cdoAbort("Units undefined!");
    }

  int zaxisID3 = zaxisCreate(ZAXIS_SURFACE, 1);

  int vlistID3 = vlistCreate();
  vlistDefVar(vlistID3, gridID3, zaxisID3, TIME_CONSTANT);
  vlistDefVarMissval(vlistID3, 0, missval);
  vlistDefVarName(vlistID3, 0, "dummy");
  vlistDefVarDatatype(vlistID3, 0, CDI_DATATYPE_INT8);

  int taxisID3 = taxisCreate(TAXIS_ABSOLUTE);

  vlistDefTaxis(vlistID3, taxisID3);

  int streamID3 = cdoStreamOpenWrite(cdoStreamName(2), cdoFiletype());

  pstreamDefVlist(streamID3, vlistID3);

  int tsID = 0;
  pstreamDefTimestep(streamID3, tsID);

  for (size_t i = 0; i < gridsize; ++i) array3[i] = missval;

  pstreamDefRecord(streamID3, 0, 0);
  pstreamWriteRecord(streamID3, array3.data(), gridsize);

  pstreamClose(streamID3);

  cdoFinish();

  return 0;
}

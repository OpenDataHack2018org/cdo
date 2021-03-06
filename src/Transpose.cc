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

      Transpose  transxy         Transpose X/Y
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

void
transxy(int gridID, double *array1, double *array2)
{
  size_t nx = gridInqXsize(gridID);
  size_t ny = gridInqYsize(gridID);
  size_t gridsize = nx * ny;

  if (gridsize > 0)
    {
      double **a2D1 = (double **) Malloc(ny * sizeof(double *));
      double **a2D2 = (double **) Malloc(nx * sizeof(double *));

      for (size_t j = 0; j < ny; ++j) a2D1[j] = array1 + j * nx;
      for (size_t i = 0; i < nx; ++i) a2D2[i] = array2 + i * ny;

      for (size_t j = 0; j < ny; ++j)
        for (size_t i = 0; i < nx; ++i) a2D2[i][j] = a2D1[j][i];

      Free(a2D1);
      Free(a2D2);
    }
  else
    {
      gridsize = gridInqSize(gridID);
      for (size_t i = 0; i < gridsize; ++i) array2[i] = array1[i];
    }
}

void *
Transpose(void *process)
{
  int nrecs;
  int varID, levelID;
  size_t nmiss;

  cdoInitialize(process);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int ngrids = vlistNgrids(vlistID1);
  for (int index = 0; index < ngrids; index++)
    {
      int gridID1 = vlistGrid(vlistID1, index);
      size_t nx = gridInqXsize(gridID1);
      int ny = gridInqYsize(gridID1);
      size_t gridsize = nx * ny;
      if (gridsize > 0)
        {
          int gridID2 = gridCreate(GRID_GENERIC, gridsize);
          gridDefXsize(gridID2, ny);
          gridDefYsize(gridID2, nx);
          vlistChangeGridIndex(vlistID2, index, gridID2);
        }
    }

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  size_t gridsize = vlistGridsizeMax(vlistID1);

  std::vector<double> array1(gridsize);
  std::vector<double> array2(gridsize);

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamReadRecord(streamID1, array1.data(), &nmiss);

          int gridID = vlistInqVarGrid(vlistID1, varID);
          transxy(gridID, array1.data(), array2.data());

          pstreamDefRecord(streamID2, varID, levelID);
          pstreamWriteRecord(streamID2, array2.data(), nmiss);
        }

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

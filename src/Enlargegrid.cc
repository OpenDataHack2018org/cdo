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
#include "grid.h"


void
genGridIndex(int gridID1, int gridID2, int *index)
{
  int i1;

  int gridtype1 = gridInqType(gridID1);
  int gridtype2 = gridInqType(gridID2);

  size_t gridsize2 = gridInqSize(gridID2);

  if (gridtype1 != gridtype2) cdoAbort("Input streams have different grid types!");

  if (index == NULL) cdoAbort("Internal problem, index not allocated!");

  for (size_t i = 0; i < gridsize2; i++) index[i] = -1;

  if (gridtype1 == GRID_LONLAT || gridtype1 == GRID_GAUSSIAN)
    {
      int nlon1 = gridInqXsize(gridID1);
      int nlat1 = gridInqYsize(gridID1);

      int nlon2 = gridInqXsize(gridID2);
      int nlat2 = gridInqYsize(gridID2);

      if (!(gridInqXvals(gridID1, NULL) && gridInqYvals(gridID1, NULL))) cdoAbort("Grid 1 has no values!");

      if (!(gridInqXvals(gridID2, NULL) && gridInqYvals(gridID2, NULL))) cdoAbort("Grid 2 has no values!");

      std::vector<double> xvals1(nlon1);
      std::vector<double> yvals1(nlat1);
      std::vector<double> xvals2(nlon2);
      std::vector<double> yvals2(nlat2);

      std::vector<int> xindex(nlon2);
      std::vector<int> yindex(nlat2);

      gridInqXvals(gridID1, xvals1.data());
      gridInqYvals(gridID1, yvals1.data());

      /* Convert lat/lon units if required */
      {
        char units[CDI_MAX_NAME];
        gridInqXunits(gridID1, units);
        grid_to_degree(units, nlon1, xvals1.data(), "grid1 center lon");
        gridInqYunits(gridID1, units);
        grid_to_degree(units, nlat1, yvals1.data(), "grid1 center lat");
      }

      gridInqXvals(gridID2, xvals2.data());
      gridInqYvals(gridID2, yvals2.data());

      /* Convert lat/lon units if required */
      {
        char units[CDI_MAX_NAME];
        gridInqXunits(gridID2, units);
        grid_to_degree(units, nlon2, xvals2.data(), "grid2 center lon");
        gridInqYunits(gridID2, units);
        grid_to_degree(units, nlat2, yvals2.data(), "grid2 center lat");
      }

      for (int i2 = 0; i2 < nlat2; i2++)
        {
          for (i1 = 0; i1 < nlat1; i1++)
            if (fabs(yvals2[i2] - yvals1[i1]) < 0.001) break;

          yindex[i2] = (i1 == nlat1) ? -1 : i1;
        }

      for (int i2 = 0; i2 < nlon2; i2++)
        {
          for (i1 = 0; i1 < nlon1; i1++)
            if (fabs(xvals2[i2] - xvals1[i1]) < 0.001) break;

          if (i1 == nlon1)
            {
              if (xvals2[i2] < 0)
                {
                  for (i1 = 0; i1 < nlon1; i1++)
                    if (fabs(xvals2[i2] + 360 - xvals1[i1]) < 0.001) break;
                }
              else if (xvals2[i2] > 180)
                {
                  for (i1 = 0; i1 < nlon1; i1++)
                    if (fabs(xvals2[i2] - 360 - xvals1[i1]) < 0.001) break;
                }
            }

          xindex[i2] = (i1 == nlon1) ? -1 : i1;
        }
      /*
      for ( i2 = 0; i2 < nlon2; i2++ )
        printf("x %d %d\n", i2, xindex[i2]);

      for ( i2 = 0; i2 < nlat2; i2++ )
        printf("y %d %d\n", i2, yindex[i2]);
      */
      int k = 0;
      for (int j = 0; j < nlat2; j++)
        for (int i = 0; i < nlon2; i++)
          {
            if (xindex[i] == -1 || yindex[j] == -1)
              index[k++] = -1;
            else
              index[k++] = yindex[j] * nlon1 + xindex[i];
          }
    }
  else
    cdoAbort("Unsupported grid type: %s", gridNamePtr(gridtype1));
}

void *
Enlargegrid(void *process)
{
  int nrecs = 0;

  cdoInitialize(process);

  operatorInputArg("grid description file or name");
  if (operatorArgc() < 1) cdoAbort("Too few arguments!");
  if (operatorArgc() > 2) cdoAbort("Too many arguments!");

  int gridID2 = cdoDefineGrid(operatorArgv()[0]);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);

  int ndiffgrids = 0;
  for (int index = 1; index < vlistNgrids(vlistID1); index++)
    if (vlistGrid(vlistID1, 0) != vlistGrid(vlistID1, index)) ndiffgrids++;

  if (ndiffgrids > 0) cdoAbort("Too many different grids in %s!", cdoGetStreamName(0).c_str());

  int gridID1 = vlistGrid(vlistID1, 0);

  size_t gridsize1 = gridInqSize(gridID1);
  size_t gridsize2 = gridInqSize(gridID2);

  std::vector<double> array1(gridsize1);
  std::vector<double> array2(gridsize2);
  std::vector<int> gindex(gridsize1);

  genGridIndex(gridID2, gridID1, gindex.data());

  int vlistID2 = vlistDuplicate(vlistID1);

  int ngrids = vlistNgrids(vlistID1);
  for (int index = 0; index < ngrids; index++)
    {
      vlistChangeGridIndex(vlistID2, index, gridID2);
    }

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());

  vlistDefTaxis(vlistID2, taxisID2);
  pstreamDefVlist(streamID2, vlistID2);

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          int varID, levelID;
          size_t nmiss;
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamReadRecord(streamID1, array1.data(), &nmiss);

          double missval1 = vlistInqVarMissval(vlistID1, varID);

          for (size_t i = 0; i < gridsize2; i++) array2[i] = missval1;
          for (size_t i = 0; i < gridsize1; i++)
            if (gindex[i] >= 0) array2[gindex[i]] = array1[i];

          nmiss = arrayNumMV(gridsize2, array2.data(), missval1);
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

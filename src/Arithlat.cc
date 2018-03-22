/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2018 Uwe Schulzweida, <uwe.schulzweida AT mpimet.m1pg.de>
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

      Arithlat   mulcoslat       Multiply with cos(lat)
      Arithlat   divcoslat       Divide by cos(lat)
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "grid.h"

void *
Arithlat(void *process)
{
  int gridID0 = -1;
  int nrecs;
  int varID, levelID;
  size_t nmiss;
  char units[CDI_MAX_NAME];
  std::vector<double> scale;

  cdoInitialize(process);

  cdoOperatorAdd("mulcoslat", func_mul, 0, NULL);
  cdoOperatorAdd("divcoslat", func_div, 0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  size_t gridsizemax = vlistGridsizeMax(vlistID1);
  std::vector<double> array(gridsizemax);

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      pstreamDefTimestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamReadRecord(streamID1, &array[0], &nmiss);

          int gridID = vlistInqVarGrid(vlistID1, varID);
          size_t gridsize = 0;

          if (gridID != gridID0)
            {
              gridID0 = gridID;

              int gridtype = gridInqType(gridID);
              int projtype = (gridtype == GRID_PROJECTION) ? gridInqProjType(gridID) : -1;
              if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || projtype == CDI_PROJ_LCC)
                {
                  gridID = gridToCurvilinear(gridID, 0);
                }
              else if (gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED)
                {
                  /* No conversion necessary */
                }
              else if (gridtype == GRID_GME)
                {
                  gridID = gridToUnstructured(gridID, 0);
                }
              else
                {
                  if (gridtype == GRID_GAUSSIAN_REDUCED)
                    cdoAbort("Unsupported grid type: %s, use CDO option -R to convert reduced to regular grid!",
                             gridNamePtr(gridtype));
                  else
                    cdoAbort("Unsupported grid type: %s", gridNamePtr(gridtype));
                }

              gridsize = gridInqSize(gridID);

              scale.resize(gridsize);
              gridInqYvals(gridID, &scale[0]);

              /* Convert lat/lon units if required */

              gridInqXunits(gridID, units);

              grid_to_radian(units, gridsize, &scale[0], "grid latitudes");

              if (operfunc == func_mul)
                for (size_t i = 0; i < gridsize; ++i) scale[i] = cos(scale[i]);
              else
                for (size_t i = 0; i < gridsize; ++i) scale[i] = 1. / cos(scale[i]);

              if (cdoVerbose)
                for (unsigned i = 0; i < 10; ++i) cdoPrint("coslat  %3d  %g", i + 1, scale[i]);
            }

          for (size_t i = 0; i < gridsize; ++i) array[i] *= scale[i];

          pstreamDefRecord(streamID2, varID, levelID);
          pstreamWriteRecord(streamID2, &array[0], nmiss);
        }

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

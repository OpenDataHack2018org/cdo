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

      Enlarge    enlarge         Enlarge fields
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

void *
Enlarge(void *process)
{
  int nrecs;
  bool linfo = true;

  cdoInitialize(process);

  operatorCheckArgc(1);

  int gridID2 = cdoDefineGrid(operatorArgv()[0]);
  size_t xsize2 = gridInqXsize(gridID2);
  size_t ysize2 = gridInqYsize(gridID2);

  if (cdoVerbose) fprintf(stderr, "gridID2 %d, xsize2 %zu, ysize2 %zu\n", gridID2, xsize2, ysize2);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  size_t gridsize2 = gridInqSize(gridID2);
  if (gridsize2 < vlistGridsizeMax(vlistID1)) cdoAbort("Gridsize of input stream is greater than new gridsize!");

  int ngrids = vlistNgrids(vlistID1);
  for (int index = 0; index < ngrids; index++)
    {
      vlistChangeGridIndex(vlistID2, index, gridID2);
    }

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());

  pstreamDefVlist(streamID2, vlistID2);

  double *array1 = (double *) Malloc(gridsize2 * sizeof(double));
  double *array2 = (double *) Malloc(gridsize2 * sizeof(double));

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
          pstreamReadRecord(streamID1, array1, &nmiss);

          double missval = vlistInqVarMissval(vlistID1, varID);
          int gridID1 = vlistInqVarGrid(vlistID1, varID);
          size_t xsize1 = gridInqXsize(gridID1);
          size_t ysize1 = gridInqYsize(gridID1);
          size_t gridsize1 = gridInqSize(gridID1);

          if (xsize1 == 0) xsize1 = 1;
          if (ysize1 == 0) ysize1 = 1;

          /* printf("%d %d %d %d\n", xsize1, ysize1, xsize2, ysize2); */
          if (xsize1 == 1 && ysize1 == ysize2 && xsize1 * ysize1 == gridsize1)
            {
              if (linfo)
                {
                  cdoPrint("Enlarge zonal");
                  linfo = false;
                }

              for (size_t iy = 0; iy < ysize2; iy++)
                for (size_t ix = 0; ix < xsize2; ix++) array2[ix + iy * xsize2] = array1[iy];

              if (nmiss) nmiss *= xsize2;
            }
          else if (ysize1 == 1 && xsize1 == xsize2 && xsize1 * ysize1 == gridsize1)
            {
              if (linfo)
                {
                  cdoPrint("Enlarge meridional");
                  linfo = false;
                }

              for (size_t iy = 0; iy < ysize2; iy++)
                for (size_t ix = 0; ix < xsize2; ix++) array2[ix + iy * xsize2] = array1[ix];

              if (nmiss) nmiss *= ysize2;
            }
          else
            {
              arrayCopy(gridsize1, array1, array2);
              for (size_t i = gridsize1; i < gridsize2; i++) array2[i] = array1[gridsize1 - 1];

              if (nmiss && DBL_IS_EQUAL(array1[gridsize1 - 1], missval)) nmiss += (gridsize2 - gridsize1);
            }

          pstreamDefRecord(streamID2, varID, levelID);
          pstreamWriteRecord(streamID2, array2, nmiss);
        }

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if (array1) Free(array1);
  if (array2) Free(array2);

  cdoFinish();

  return 0;
}

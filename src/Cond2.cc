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

      Cond2      ifthenelse      If then else
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

void *
Cond2(void *process)
{
  enum
  {
    FILL_NONE,
    FILL_TS,
    FILL_REC
  };
  int filltype = FILL_NONE;
  int nrecs, nrecs2, nvars = 0, nlev;
  int varID, levelID;
  double missval1 = -9.E33;
  size_t nmiss1;
  std::vector<std::vector<size_t>> varnmiss1;
  std::vector<std::vector<double>> vardata1;

  cdoInitialize(process);

  cdoOperatorAdd("ifthenelse", 0, 0, NULL);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));
  int streamID2 = cdoStreamOpenRead(cdoStreamName(1));
  int streamID3 = cdoStreamOpenRead(cdoStreamName(2));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = cdoStreamInqVlist(streamID2);
  int vlistID3 = cdoStreamInqVlist(streamID3);
  int vlistID4 = vlistDuplicate(vlistID2);

  int taxisID2 = vlistInqTaxis(vlistID2);
  int taxisID4 = taxisDuplicate(taxisID2);
  vlistDefTaxis(vlistID4, taxisID4);

  int ntsteps1 = vlistNtsteps(vlistID1);
  int ntsteps2 = vlistNtsteps(vlistID2);
  if (ntsteps1 == 0) ntsteps1 = 1;
  if (ntsteps2 == 0) ntsteps2 = 1;

  if (vlistNrecs(vlistID1) == 1 && vlistNrecs(vlistID2) != 1)
    {
      filltype = FILL_REC;
      cdoPrint("Filling up stream1 >%s< by copying the first record.", cdoGetStreamName(0).c_str());
    }

  if (filltype == FILL_NONE) vlistCompare(vlistID1, vlistID2, CMP_DIM);

  vlistCompare(vlistID2, vlistID3, CMP_DIM);

  nospec(vlistID1);
  nospec(vlistID2);
  nospec(vlistID3);

  int streamID4 = cdoStreamOpenWrite(cdoStreamName(3), cdoFiletype());
  pstreamDefVlist(streamID4, vlistID4);

  size_t gridsize = vlistGridsizeMax(vlistID1);

  if (filltype == FILL_REC && gridsize != gridInqSize(vlistGrid(vlistID1, 0)))
    cdoAbort("Stream1 >%s< has wrong gridsize!", cdoGetStreamName(0).c_str());

  std::vector<double> array1(gridsize);
  std::vector<double> array2(gridsize);
  std::vector<double> array3(gridsize);
  std::vector<double> array4(gridsize);

  if (cdoVerbose)
    cdoPrint("Number of timesteps: file1 %d, file2 %d, file3 %d", ntsteps1, ntsteps2, vlistNtsteps(vlistID3));

  if (filltype == FILL_NONE)
    {
      if (ntsteps1 == 1 && ntsteps2 != 1)
        {
          filltype = FILL_TS;
          cdoPrint("Filling up stream1 >%s< by copying the first timestep.", cdoGetStreamName(0).c_str());

          nvars = vlistNvars(vlistID1);
          vardata1.resize(nvars);
          varnmiss1.resize(nvars);
          for (varID = 0; varID < nvars; varID++)
            {
              gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
              nlev = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
              vardata1[varID].resize(nlev * gridsize);
              varnmiss1[varID].resize(nlev);
            }
        }
    }

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID2, tsID)))
    {
      nrecs = cdoStreamInqTimestep(streamID3, tsID);
      if (nrecs == 0) cdoAbort("Input streams have different number of timesteps!");

      if (tsID == 0 || filltype == FILL_NONE)
        {
          nrecs2 = cdoStreamInqTimestep(streamID1, tsID);
          if (nrecs2 == 0) cdoAbort("Input streams have different number of timesteps!");
        }

      taxisCopyTimestep(taxisID4, taxisID2);
      pstreamDefTimestep(streamID4, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          size_t nmiss;
          pstreamInqRecord(streamID2, &varID, &levelID);
          pstreamReadRecord(streamID2, &array2[0], &nmiss);

          pstreamInqRecord(streamID3, &varID, &levelID);
          pstreamReadRecord(streamID3, &array3[0], &nmiss);

          if (tsID == 0 || filltype == FILL_NONE)
            {
              if (recID == 0 || filltype != FILL_REC)
                {
                  pstreamInqRecord(streamID1, &varID, &levelID);
                  pstreamReadRecord(streamID1, &array1[0], &nmiss1);
                }

              if (filltype == FILL_TS)
                {
                  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
                  size_t offset = gridsize * levelID;
                  arrayCopy(gridsize, &array1[0], &vardata1[varID][offset]);
                  varnmiss1[varID][levelID] = nmiss1;
                }
            }
          else if (filltype == FILL_TS)
            {
              gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
              size_t offset = gridsize * levelID;
              arrayCopy(gridsize, &vardata1[varID][offset], &array1[0]);
              nmiss1 = varnmiss1[varID][levelID];
            }

          gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
          double missval2 = vlistInqVarMissval(vlistID2, varID);
          if (recID == 0 || filltype != FILL_REC) missval1 = vlistInqVarMissval(vlistID1, varID);

          for (size_t i = 0; i < gridsize; i++)
            array4[i]
                = DBL_IS_EQUAL(array1[i], missval1) ? missval2 : !DBL_IS_EQUAL(array1[i], 0.) ? array2[i] : array3[i];

          nmiss = arrayNumMV(gridsize, &array4[0], missval2);
          pstreamDefRecord(streamID4, varID, levelID);
          pstreamWriteRecord(streamID4, &array4[0], nmiss);
        }

      tsID++;
    }

  pstreamClose(streamID4);
  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

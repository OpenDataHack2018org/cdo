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

void genGridIndex(int gridID1, int gridID2, int *index);

void *
Mergegrid(void *process)
{
  int varID, levelID;
  int nrecs = 0;
  size_t nmiss1, nmiss2;
  int index;

  cdoInitialize(process);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID3 = taxisDuplicate(taxisID1);

  int streamID2 = cdoStreamOpenRead(cdoStreamName(1));
  int vlistID2 = cdoStreamInqVlist(streamID2);

  vlistCompare(vlistID1, vlistID2, CMP_NAME | CMP_NLEVEL);

  int ndiffgrids = 0;
  for (index = 1; index < vlistNgrids(vlistID1); index++)
    if (vlistGrid(vlistID1, 0) != vlistGrid(vlistID1, index)) ndiffgrids++;

  if (ndiffgrids > 0) cdoAbort("Too many different grids in %s!", cdoGetStreamName(0).c_str());

  ndiffgrids = 0;
  for (index = 1; index < vlistNgrids(vlistID2); index++)
    if (vlistGrid(vlistID2, 0) != vlistGrid(vlistID2, index)) ndiffgrids++;

  if (ndiffgrids > 0) cdoAbort("Too many different grids in %s!", cdoGetStreamName(1).c_str());

  int gridID1 = vlistGrid(vlistID1, 0);
  int gridID2 = vlistGrid(vlistID2, 0);

  size_t gridsize1 = gridInqSize(gridID1);
  size_t gridsize2 = gridInqSize(gridID2);

  std::vector<double> array1(gridsize1);
  std::vector<double> array2(gridsize2);
  std::vector<int> gindex(gridsize2);

  genGridIndex(gridID1, gridID2, gindex.data());

  int vlistID3 = vlistDuplicate(vlistID1);
  int streamID3 = cdoStreamOpenWrite(cdoStreamName(2), cdoFiletype());

  vlistDefTaxis(vlistID3, taxisID3);
  pstreamDefVlist(streamID3, vlistID3);

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID3, taxisID1);

      int nrecs2 = cdoStreamInqTimestep(streamID2, tsID);
      if (nrecs2 == 0) cdoAbort("Input streams have different number of timesteps!");

      if (nrecs != nrecs2) cdoAbort("Input streams have different number of records!");

      pstreamDefTimestep(streamID3, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID2, &varID, &levelID);
          pstreamReadRecord(streamID2, array2.data(), &nmiss2);

          double missval2 = vlistInqVarMissval(vlistID2, varID);

          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamReadRecord(streamID1, array1.data(), &nmiss1);

          double missval1 = vlistInqVarMissval(vlistID1, varID);

          for (size_t i = 0; i < gridsize2; i++)
            {
              if (gindex[i] >= 0 && !DBL_IS_EQUAL(array2[i], missval2))
                {
                  array1[gindex[i]] = array2[i];
                }
            }

          if (nmiss1)
            {
              nmiss1 = 0;
              for (size_t i = 0; i < gridsize1; i++)
                if (DBL_IS_EQUAL(array1[i], missval1)) nmiss1++;
            }

          pstreamDefRecord(streamID3, varID, levelID);
          pstreamWriteRecord(streamID3, array1.data(), nmiss1);
        }

      tsID++;
    }

  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

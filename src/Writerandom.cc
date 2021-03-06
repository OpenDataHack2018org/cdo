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

      Writerandom writerandom
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

void *
Writerandom(void *process)
{
  size_t gridsize;
  int nrecs;
  int varID, levelID;
  int rindex;

  cdoInitialize(process);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());

  pstreamDefVlist(streamID2, vlistID2);

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      std::vector<std::vector<double>> recdata(nrecs);
      std::vector<int> recvarID(nrecs);
      std::vector<int> reclevelID(nrecs);
      std::vector<size_t> recnmiss(nrecs);
      std::vector<int> recindex(nrecs);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
          recvarID[recID] = varID;
          reclevelID[recID] = levelID;
          recdata[recID].resize(gridsize);
          pstreamReadRecord(streamID1, recdata[recID].data(), &recnmiss[recID]);
        }

      for (int recID = 0; recID < nrecs; recID++) recindex[recID] = -1;

      for (rindex = nrecs - 1; rindex >= 0; rindex--)
        {
          int index = (int) (rindex * ((double) rand()) / ((double) RAND_MAX));
          /*	printf("rindex %d %d\n", rindex, index); */
          int ipos = -1;
          for (int recID = 0; recID < nrecs; recID++)
            {
              if (recindex[recID] == -1) ipos++;
              if (recindex[recID] == -1 && ipos == index)
                {
                  recindex[recID] = rindex;
                  break;
                }
            }
        }

      /*
      for ( int recID = 0; recID < nrecs; recID++ )
        printf("recID %d %d\n", recID, recindex[recID]);
      */
      for (int recID = 0; recID < nrecs; recID++)
        if (recindex[recID] == -1) cdoAbort("Internal problem! Random initialize.");

      for (int recID = 0; recID < nrecs; recID++)
        {
          rindex = recindex[recID];
          varID = recvarID[rindex];
          levelID = reclevelID[rindex];
          gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
          pstreamDefRecord(streamID2, varID, levelID);
          pstreamWriteRecord(streamID2, recdata[rindex].data(), recnmiss[rindex]);
        }

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

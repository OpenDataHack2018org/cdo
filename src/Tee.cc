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

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

void *
Tee(void *process)
{
  int nrecs;
  int varID, levelID;
  size_t nmiss;

  cdoInitialize(process);

  bool lcopy = UNCHANGED_RECORD;

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int taxisID1 = vlistInqTaxis(vlistID1);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  int streamID3 = cdoStreamOpenWrite(cdoStreamName(2), cdoFiletype());

  int vlistID2 = vlistDuplicate(vlistID1);
  int vlistID3 = vlistDuplicate(vlistID1);

  int taxisID2 = taxisDuplicate(taxisID1);
  int taxisID3 = taxisDuplicate(taxisID1);

  vlistDefTaxis(vlistID2, taxisID2);
  vlistDefTaxis(vlistID3, taxisID3);

  pstreamDefVlist(streamID2, vlistID2);
  pstreamDefVlist(streamID3, vlistID3);

  size_t gridsizemax = vlistGridsizeMax(vlistID1);
  std::vector<double> array(gridsizemax);

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      taxisCopyTimestep(taxisID3, taxisID1);

      pstreamDefTimestep(streamID2, tsID);
      pstreamDefTimestep(streamID3, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          if (lcopy)
            {
              pstreamInqRecord(streamID1, &varID, &levelID);

              pstreamDefRecord(streamID2, varID, levelID);
              pstreamCopyRecord(streamID2, streamID1);

              pstreamDefRecord(streamID3, varID, levelID);
              pstreamCopyRecord(streamID3, streamID1);
            }
          else
            {
              pstreamInqRecord(streamID1, &varID, &levelID);
              pstreamReadRecord(streamID1, array.data(), &nmiss);

              pstreamDefRecord(streamID2, varID, levelID);
              pstreamWriteRecord(streamID2, array.data(), nmiss);

              pstreamDefRecord(streamID3, varID, levelID);
              pstreamWriteRecord(streamID3, array.data(), nmiss);
            }
        }

      tsID++;
    }

  pstreamClose(streamID1);
  pstreamClose(streamID2);
  pstreamClose(streamID3);

  cdoFinish();

  return 0;
}

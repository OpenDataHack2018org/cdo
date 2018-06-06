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
Recttocomplex(void *process)
{
  int nrecs;
  int varID, levelID;
  size_t nmiss;

  cdoInitialize(process);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));
  int vlistID1 = cdoStreamInqVlist(streamID1);

  int streamID2 = cdoStreamOpenRead(cdoStreamName(1));
  int vlistID2 = cdoStreamInqVlist(streamID2);

  int vlistID3 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID3);


  int nvars = vlistNvars(vlistID3);
  for (varID = 0; varID < nvars; ++varID)
    {
      int datatype = vlistInqVarDatatype(vlistID2, varID);
      datatype =  (datatype == CDI_DATATYPE_FLT64) ? CDI_DATATYPE_CPX64 : CDI_DATATYPE_CPX32;
      vlistDefVarDatatype(vlistID3, varID, datatype);
    }

  int streamID3 = cdoStreamOpenWrite(cdoStreamName(2), cdoFiletype());
  pstreamDefVlist(streamID3, vlistID3);

  size_t gridsizemax = vlistGridsizeMax(vlistID1);
  std::vector<double> array1(gridsizemax);
  std::vector<double> array2(gridsizemax);
  std::vector<double> array3(2 *gridsizemax);

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID3, taxisID1);
      pstreamDefTimestep(streamID3, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamDefRecord(streamID3, varID, levelID);

          pstreamInqRecord(streamID2, &varID, &levelID);

          size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

          pstreamReadRecord(streamID1, array1.data(), &nmiss);
          pstreamReadRecord(streamID2, array2.data(), &nmiss);

          for (size_t i = 0; i < gridsize; ++i)
            {
              array3[2 * i]     = array1[i];
              array3[2 * i + 1] = array2[i];
            }

          pstreamWriteRecord(streamID3, array3.data(), nmiss);
        }

      tsID++;
    }

  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  vlistDestroy(vlistID3);

  cdoFinish();

  return 0;
}

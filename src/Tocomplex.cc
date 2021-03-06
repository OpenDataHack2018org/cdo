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
Tocomplex(void *process)
{
  int nrecs;
  int varID, levelID;
  size_t nmiss;

  cdoInitialize(process);

  int RETOCOMPLEX = cdoOperatorAdd("retocomplex", 0, 0, NULL);
  int IMTOCOMPLEX = cdoOperatorAdd("imtocomplex", 0, 0, NULL);

  int operatorID = cdoOperatorID();

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int nvars = vlistNvars(vlistID2);
  for (varID = 0; varID < nvars; ++varID)
    {
      int datatype = vlistInqVarDatatype(vlistID2, varID);
      datatype = (datatype == CDI_DATATYPE_FLT64) ? CDI_DATATYPE_CPX64 : CDI_DATATYPE_CPX32;
      vlistDefVarDatatype(vlistID2, varID, datatype);
    }

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  // if (cdoFiletype() != CDI_FILETYPE_EXT) cdoAbort("Complex numbers need EXTRA format; used CDO option -f ext!");
  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  size_t gridsize = vlistGridsizeMax(vlistID1);
  std::vector<double> array1(gridsize);
  std::vector<double> array2(2 * gridsize);

  int tsID = 0;
  int tsID2 = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID2++);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamDefRecord(streamID2, varID, levelID);

          gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

          pstreamReadRecord(streamID1, array1.data(), &nmiss);

          if (operatorID == RETOCOMPLEX)
            {
              for (size_t i = 0; i < gridsize; ++i)
                {
                  array2[2 * i] = array1[i];
                  array2[2 * i + 1] = 0;
                }
            }
          else if (operatorID == IMTOCOMPLEX)
            {
              for (size_t i = 0; i < gridsize; ++i)
                {
                  array2[2 * i] = 0;
                  array2[2 * i + 1] = array1[i];
                }
            }

          pstreamWriteRecord(streamID2, array2.data(), nmiss);
        }

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  vlistDestroy(vlistID2);

  cdoFinish();

  return 0;
}

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
Complextorect(void *process)
{
  int nrecs;
  int varID, levelID;
  size_t nmiss;

  cdoInitialize(process);

  int COMPLEXTORECT = cdoOperatorAdd("complextorect", 0, 0, NULL);
  int COMPLEXTOPOL = cdoOperatorAdd("complextopol", 0, 0, NULL);

  int operatorID = cdoOperatorID();

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);
  int vlistID3 = vlistDuplicate(vlistID1);

  int nvars = vlistNvars(vlistID2);
  for (varID = 0; varID < nvars; ++varID)
    {
      int datatype = vlistInqVarDatatype(vlistID2, varID);
      datatype =  (datatype == CDI_DATATYPE_CPX64) ? CDI_DATATYPE_FLT64 : CDI_DATATYPE_FLT32;
      vlistDefVarDatatype(vlistID2, varID, datatype);
      vlistDefVarDatatype(vlistID3, varID, datatype);
    }

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  int taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);
  vlistDefTaxis(vlistID3, taxisID3);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  int streamID3 = cdoStreamOpenWrite(cdoStreamName(2), cdoFiletype());

  pstreamDefVlist(streamID2, vlistID2);
  pstreamDefVlist(streamID3, vlistID3);

  size_t gridsizemax = vlistGridsizeMax(vlistID1);
  std::vector<double> array1(2 * gridsizemax);
  std::vector<double> array2(gridsizemax);
  std::vector<double> array3(gridsizemax);

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      taxisCopyTimestep(taxisID3, taxisID1);

      pstreamDefTimestep(streamID2, tsID);
      pstreamDefTimestep(streamID3, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamDefRecord(streamID2, varID, levelID);
          pstreamDefRecord(streamID3, varID, levelID);

          size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

          pstreamReadRecord(streamID1, array1.data(), &nmiss);

          double missval1 = vlistInqVarMissval(vlistID1, varID);
          double missval2 = missval1;

          if (operatorID == COMPLEXTORECT)
            {
              for (size_t i = 0; i < gridsize; ++i)
                {
                  array2[i] = array1[2 * i];
                  array3[i] = array1[2 * i + 1];
                }
            }
          else if (operatorID == COMPLEXTOPOL)
            {
              for (size_t i = 0; i < gridsize; ++i)
                {
                  array2[i] = SQRTMN(ADDMN(MULMN(array1[2 * i], array1[2 * i]), MULMN(array1[2 * i + 1], array1[2 * i + 1])));
                  array3[i] = (DBL_IS_EQUAL(array1[2 * i], missval1) || DBL_IS_EQUAL(array1[2 * i + 1], missval1))
                                  ? missval1
                                  : atan2(array1[2 * i + 1], array1[2 * i]);
                }
            }

          pstreamWriteRecord(streamID2, array2.data(), nmiss);
          pstreamWriteRecord(streamID3, array3.data(), nmiss);
        }

      tsID++;
    }

  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  vlistDestroy(vlistID2);
  vlistDestroy(vlistID3);

  cdoFinish();

  return 0;
}

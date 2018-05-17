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

      Selyearidx    selyearidx         Select index of year
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

void *
Selyearidx(void *process)
{
  int varID, levelID;
  int vtime;
  int year, month, day;
  int hour, minute, second;
  size_t nmiss;

  cdoInitialize(process);

  cdoOperatorAdd("selyearidx", 0, 0, NULL);
  cdoOperatorAdd("seltimeidx", 1, 0, NULL);

  bool ltime = cdoOperatorF1(cdoOperatorID());

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));
  int vlistID1 = cdoStreamInqVlist(streamID1);

  int taxisID1 = vlistInqTaxis(vlistID1);

  int streamID2 = cdoStreamOpenRead(cdoStreamName(1));
  int vlistID2 = cdoStreamInqVlist(streamID2);
  int taxisID2 = vlistInqTaxis(vlistID2);

  vlistCompare(vlistID1, vlistID2, CMP_ALL);

  int vlistID3 = vlistDuplicate(vlistID2);
  int taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID3);

  int streamID3 = cdoStreamOpenWrite(cdoStreamName(2), cdoFiletype());
  pstreamDefVlist(streamID3, vlistID3);

  int maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recinfo(maxrecs);

  size_t gridsizemax = vlistGridsizeMax(vlistID1);

  std::vector<double> array(gridsizemax);

  Field **vars1 = field_malloc(vlistID1, FIELD_PTR);
  Field **vars2 = field_malloc(vlistID1, FIELD_PTR);

  int nvars = vlistNvars(vlistID1);
  for (varID = 0; varID < nvars; ++varID)
    {
      size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
      int nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      double missval = vlistInqVarMissval(vlistID2, varID);
      for (levelID = 0; levelID < nlevs; ++levelID)
        {
          for (size_t i = 0; i < gridsize; ++i) vars2[varID][levelID].ptr[i] = missval;
        }
    }

  int calendar = taxisInqCalendar(taxisID2);

  int tsID = 0;
  int tsID2 = 0;
  int tsID3 = 0;
  while (TRUE)
    {
      int nrecs = cdoStreamInqTimestep(streamID1, tsID);
      if (nrecs == 0) break;

      int64_t vdate = taxisInqVdate(taxisID1);
      // int vtime = taxisInqVtime(taxisID1);
      cdiDecodeDate(vdate, &year, &month, &day);
      int year1 = year;

      bool lexactdate = gridsizemax == 1 && nrecs == 1;

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamReadRecord(streamID1, vars1[varID][levelID].ptr, &nmiss);
          vars1[varID][levelID].nmiss = nmiss;

          if (tsID == 0)
            {
              recinfo[recID].varID = varID;
              recinfo[recID].levelID = levelID;
              recinfo[recID].lconst = vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT;
            }
        }

      int nsets = 0;
      while ((nrecs = cdoStreamInqTimestep(streamID2, tsID2)))
        {
          vdate = taxisInqVdate(taxisID2);
          vtime = taxisInqVtime(taxisID2);
          cdiDecodeDate(vdate, &year, &month, &day);
          cdiDecodeTime(vtime, &hour, &minute, &second);

          if (year1 != year) break;

          for (int recID = 0; recID < nrecs; recID++)
            {
              pstreamInqRecord(streamID2, &varID, &levelID);
              pstreamReadRecord(streamID2, array.data(), &nmiss);

              size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
              for (size_t i = 0; i < gridsize; ++i)
                if (nsets == (int) lround(vars1[varID][levelID].ptr[i]))
                  {
                    if (lexactdate) taxisCopyTimestep(taxisID3, taxisID2);
                    if (ltime)
                      vars2[varID][levelID].ptr[i]
                          = date_to_julday(calendar, vdate) + (hour * 3600 + minute * 60 + second) / 86400.;
                    else
                      vars2[varID][levelID].ptr[i] = array[i];
                  }
            }

          nsets++;
          tsID2++;
        }

      if (nsets)
        {
          if (!lexactdate) taxisCopyTimestep(taxisID3, taxisID1);
          pstreamDefTimestep(streamID3, tsID3);

          for (int recID = 0; recID < maxrecs; recID++)
            {
              if (tsID && recinfo[recID].lconst) continue;

              varID = recinfo[recID].varID;
              levelID = recinfo[recID].levelID;
              pstreamDefRecord(streamID3, varID, levelID);
              size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
              double missval = vlistInqVarMissval(vlistID2, varID);
              nmiss = arrayNumMV(gridsize, vars2[varID][levelID].ptr, missval);
              pstreamWriteRecord(streamID3, vars2[varID][levelID].ptr, nmiss);
            }

          tsID3++;
        }

      if (nrecs == 0) break;
      tsID++;
    }

  field_free(vars2, vlistID1);
  field_free(vars1, vlistID1);

  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

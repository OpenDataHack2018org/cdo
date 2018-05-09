/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2007 Brockmann Consult
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

      Seascount   seascount         Seasonal counts
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

void *
Seascount(void *process)
{
  int vdate0 = 0, vtime0 = 0;
  int nrecs;
  int varID, levelID;
  size_t nmiss;
  int year, month, day, seas0 = 0;
  int oldmon = 0;

  cdoInitialize(process);

  cdoOperatorAdd("seascount", 0, 0, NULL);

  int season_start = get_season_start();

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recinfo(maxrecs);

  size_t gridsize = vlistGridsizeMax(vlistID1);
  if (vlistNumber(vlistID1) != CDI_REAL) gridsize *= 2;

  Field field;
  field_init(&field);
  field.ptr = (double *) Malloc(gridsize * sizeof(double));

  Field **vars1 = field_malloc(vlistID1, FIELD_PTR);

  int tsID = 0;
  int otsID = 0;
  while (TRUE)
    {
      int nsets = 0;
      bool newseas = false;
      while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
        {
          int vdate = taxisInqVdate(taxisID1);
          int vtime = taxisInqVtime(taxisID1);
          cdiDecodeDate(vdate, &year, &month, &day);

          int newmon = month;
          if (season_start == START_DEC && newmon == 12) newmon = 0;

          int seas = month_to_season(month);

          if (nsets == 0)
            {
              seas0 = seas;
              oldmon = newmon;
            }

          if (newmon < oldmon) newseas = true;

          if ((seas != seas0) || newseas) break;

          oldmon = newmon;

          for (int recID = 0; recID < nrecs; recID++)
            {
              pstreamInqRecord(streamID1, &varID, &levelID);

              if (tsID == 0)
                {
                  recinfo[recID].varID = varID;
                  recinfo[recID].levelID = levelID;
                  recinfo[recID].lconst = vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT;
                }
              // number of words per value; real:1  complex:2
              int nwpv = vars1[varID][levelID].nwpv;
              gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

              if (nsets == 0)
                {
                  for (size_t i = 0; i < nwpv * gridsize; i++) vars1[varID][levelID].ptr[i] = vars1[varID][levelID].missval;
                  vars1[varID][levelID].nmiss = gridsize;
                }

              pstreamReadRecord(streamID1, field.ptr, &nmiss);
              field.nmiss = nmiss;
              field.grid = vars1[varID][levelID].grid;
              field.missval = vars1[varID][levelID].missval;

              farcount(&vars1[varID][levelID], field);
            }

          vdate0 = vdate;
          vtime0 = vtime;
          nsets++;
          tsID++;
        }

      if (nrecs == 0 && nsets == 0) break;

      taxisDefVdate(taxisID2, vdate0);
      taxisDefVtime(taxisID2, vtime0);
      pstreamDefTimestep(streamID2, otsID);

      for (int recID = 0; recID < maxrecs; recID++)
        {
          if (otsID && recinfo[recID].lconst) continue;

          int varID = recinfo[recID].varID;
          int levelID = recinfo[recID].levelID;
          pstreamDefRecord(streamID2, varID, levelID);
          pstreamWriteRecord(streamID2, vars1[varID][levelID].ptr, vars1[varID][levelID].nmiss);
        }

      if (nrecs == 0) break;
      otsID++;
    }

  field_free(vars1, vlistID1);

  if (field.ptr) Free(field.ptr);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult
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

      Seaspctl   seaspctl        Seasonal percentiles
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "percentiles_hist.h"
#include "percentiles.h"
#include "datetime.h"

void *
Seaspctl(void *process)
{
  TimeStat timestat_date = TimeStat::MEAN;
  int nrecs;
  int gridID, varID, levelID;
  int year, month, day, seas0 = 0;
  size_t nmiss;
  int nlevels;
  int oldmon = 0;
  int season_start;
  double missval;

  cdoInitialize(process);

  cdoOperatorAdd("seaspctl", func_pctl, 0, NULL);

  operatorInputArg("percentile number");
  double pn = parameter2double(operatorArgv()[0]);
  percentile_check_number(pn);

  season_start = get_season_start();

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));
  int streamID2 = cdoStreamOpenRead(cdoStreamName(1));
  int streamID3 = cdoStreamOpenRead(cdoStreamName(2));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = cdoStreamInqVlist(streamID2);
  int vlistID3 = cdoStreamInqVlist(streamID3);
  int vlistID4 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, CMP_ALL);
  vlistCompare(vlistID1, vlistID3, CMP_ALL);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = vlistInqTaxis(vlistID2);
  int taxisID3 = vlistInqTaxis(vlistID3);
  /* TODO - check that time axes 2 and 3 are equal */

  int taxisID4 = taxisDuplicate(taxisID1);
  taxisWithBounds(taxisID4);
  vlistDefTaxis(vlistID4, taxisID4);

  int streamID4 = cdoStreamOpenWrite(cdoStreamName(3), cdoFiletype());
  pstreamDefVlist(streamID4, vlistID4);

  int nvars = vlistNvars(vlistID1);

  int maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recinfo(maxrecs);

  DateTimeList dtlist;
  dtlist.setStat(timestat_date);
  dtlist.setCalendar(taxisInqCalendar(taxisID1));

  size_t gridsize = vlistGridsizeMax(vlistID1);

  Field field;
  field_init(&field);
  field.ptr = (double *) Malloc(gridsize * sizeof(double));

  Field **vars1 = (Field **) Malloc(nvars * sizeof(Field *));
  HISTOGRAM_SET *hset = hsetCreate(nvars);

  for (varID = 0; varID < nvars; varID++)
    {
      gridID = vlistInqVarGrid(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      missval = vlistInqVarMissval(vlistID1, varID);

      vars1[varID] = (Field *) Malloc(nlevels * sizeof(Field));
      hsetCreateVarLevels(hset, varID, nlevels, gridID);

      for (levelID = 0; levelID < nlevels; levelID++)
        {
          vars1[varID][levelID].grid = gridID;
          vars1[varID][levelID].nmiss = 0;
          vars1[varID][levelID].missval = missval;
          vars1[varID][levelID].ptr = (double *) Malloc(gridsize * sizeof(double));
        }
    }

  int tsID = 0;
  int otsID = 0;
  while (TRUE)
    {
      nrecs = cdoStreamInqTimestep(streamID2, otsID);
      if (nrecs != cdoStreamInqTimestep(streamID3, otsID))
        cdoAbort("Number of records at time step %d of %s and %s differ!", otsID + 1, cdoGetStreamName(1).c_str(),
                 cdoGetStreamName(2).c_str());

      int64_t vdate2 = taxisInqVdate(taxisID2);
      int vtime2 = taxisInqVtime(taxisID2);
      int64_t vdate3 = taxisInqVdate(taxisID3);
      int vtime3 = taxisInqVtime(taxisID3);
      if (vdate2 != vdate3 || vtime2 != vtime3)
        cdoAbort("Verification dates at time step %d of %s and %s differ!", otsID + 1, cdoGetStreamName(1).c_str(),
                 cdoGetStreamName(2).c_str());

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID2, &varID, &levelID);
          pstreamReadRecord(streamID2, vars1[varID][levelID].ptr, &nmiss);
          vars1[varID][levelID].nmiss = nmiss;
        }

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID3, &varID, &levelID);
          pstreamReadRecord(streamID3, field.ptr, &nmiss);
          field.nmiss = nmiss;
          field.grid = vars1[varID][levelID].grid;
          field.missval = vars1[varID][levelID].missval;

          hsetDefVarLevelBounds(hset, varID, levelID, &vars1[varID][levelID], &field);
        }

      int nsets = 0;
      bool newseas = false;
      while (nrecs && (nrecs = cdoStreamInqTimestep(streamID1, tsID)))
        {
          dtlist.taxisInqTimestep(taxisID1, nsets);
          int64_t vdate1 = dtlist.getVdate(nsets);

          cdiDecodeDate(vdate1, &year, &month, &day);

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
              pstreamReadRecord(streamID1, vars1[varID][levelID].ptr, &nmiss);
              vars1[varID][levelID].nmiss = nmiss;

              hsetAddVarLevelValues(hset, varID, levelID, &vars1[varID][levelID]);
            }

          nsets++;
          tsID++;
        }

      if (nrecs == 0 && nsets == 0) break;

      for (varID = 0; varID < nvars; varID++)
        {
          if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT) continue;
          nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));

          for (levelID = 0; levelID < nlevels; levelID++)
            hsetGetVarLevelPercentiles(&vars1[varID][levelID], hset, varID, levelID, pn);
        }

      dtlist.statTaxisDefTimestep(taxisID4, nsets);
      pstreamDefTimestep(streamID4, otsID);

      for (int recID = 0; recID < maxrecs; recID++)
        {
          if (otsID && recinfo[recID].lconst) continue;

          int varID = recinfo[recID].varID;
          int levelID = recinfo[recID].levelID;
          pstreamDefRecord(streamID4, varID, levelID);
          pstreamWriteRecord(streamID4, vars1[varID][levelID].ptr, vars1[varID][levelID].nmiss);
        }

      if (nrecs == 0) break;
      otsID++;
    }

  for (varID = 0; varID < nvars; varID++)
    {
      nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      for (levelID = 0; levelID < nlevels; levelID++) Free(vars1[varID][levelID].ptr);
      Free(vars1[varID]);
    }

  Free(vars1);
  hsetDestroy(hset);

  if (field.ptr) Free(field.ptr);

  pstreamClose(streamID4);
  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

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

      Timpctl    timpctl         Time percentiles
      Hourpctl   hourpctl        Hourly percentiles
      Daypctl    daypctl         Daily percentiles
      Monpctl    monpctl         Monthly percentiles
      Yearpctl   yearpctl        Yearly percentiles
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "percentiles_hist.h"
#include "percentiles.h"
#include "datetime.h"

static void
timpctl(int operatorID)
{
  TimeStat timestat_date = TimeStat::MEAN;
  char indate1[DATE_LEN + 1], indate2[DATE_LEN + 1];
  int nrecs;
  int gridID, varID, levelID;
  size_t nmiss;
  int nlevels;

  operatorInputArg("percentile number");
  double pn = parameter2double(operatorArgv()[0]);
  percentile_check_number(pn);

  int cmplen = DATE_LEN - cdoOperatorF2(operatorID);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));
  int streamID2 = cdoStreamOpenRead(cdoStreamName(1));
  int streamID3 = cdoStreamOpenRead(cdoStreamName(2));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = cdoStreamInqVlist(streamID2);
  int vlistID3 = cdoStreamInqVlist(streamID3);
  int vlistID4 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, CMP_ALL);
  vlistCompare(vlistID1, vlistID3, CMP_ALL);

  if (cdoOperatorF2(operatorID) == 16) vlistDefNtsteps(vlistID4, 1);

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

  Field **vars1 = field_malloc(vlistID1, FIELD_PTR);
  HISTOGRAM_SET *hset = hsetCreate(nvars);

  for (varID = 0; varID < nvars; varID++)
    {
      gridID = vlistInqVarGrid(vlistID1, varID);
      nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));

      hsetCreateVarLevels(hset, varID, nlevels, gridID);
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
      while (nrecs && (nrecs = cdoStreamInqTimestep(streamID1, tsID)))
        {
          dtlist.taxisInqTimestep(taxisID1, nsets);
          int64_t vdate1 = dtlist.getVdate(nsets);
          int vtime1 = dtlist.getVtime(nsets);

          if (nsets == 0) SET_DATE(indate2, vdate1, vtime1);
          SET_DATE(indate1, vdate1, vtime1);

          if (DATE_IS_NEQ(indate1, indate2, cmplen)) break;

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

  field_free(vars1, vlistID1);
  hsetDestroy(hset);

  if (field.ptr) Free(field.ptr);

  pstreamClose(streamID4);
  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);
}

void *
Timpctl(void *process)
{
  int operatorID;

  cdoInitialize(process);

  // clang-format off
  cdoOperatorAdd("timpctl",  func_pctl, 31, NULL);
  cdoOperatorAdd("yearpctl", func_pctl, 10, NULL);
  cdoOperatorAdd("monpctl",  func_pctl,  8, NULL);
  cdoOperatorAdd("daypctl",  func_pctl,  6, NULL);
  cdoOperatorAdd("hourpctl", func_pctl,  4, NULL);
  // clang-format on

  operatorID = cdoOperatorID();

  timpctl(operatorID);

  cdoFinish();

  return 0;
}

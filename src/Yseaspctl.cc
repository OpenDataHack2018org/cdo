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

      Yseaspctl  yseaspctl       Multi-year seasonal percentiles
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "percentiles_hist.h"
#include "percentiles.h"

#define NSEAS 4

struct date_time_t
{
  int64_t vdate;
  int vtime;
};

void set_date(int64_t vdate_new, int vtime_new, date_time_t *datetime);

int getmonthday(int64_t date);

void *
Yseaspctl(void *process)
{
  int varID;
  int gridID;
  int64_t vdate;
  int vtime;
  int year, month, day, seas;
  int nrecs;
  int levelID;
  size_t nmiss;
  int nlevels;
  long nsets[NSEAS];
  date_time_t datetime1[NSEAS], datetime2[NSEAS];
  Field **vars1[NSEAS];
  HISTOGRAM_SET *hsets[NSEAS];

  cdoInitialize(process);
  cdoOperatorAdd("yseaspctl", func_pctl, 0, NULL);

  operatorInputArg("percentile number");
  double pn = parameter2double(operatorArgv()[0]);
  percentile_check_number(pn);

  for (seas = 0; seas < NSEAS; seas++)
    {
      vars1[seas] = NULL;
      hsets[seas] = NULL;
      nsets[seas] = 0;
      datetime1[seas].vdate = 0;
      datetime1[seas].vtime = 0;
      datetime2[seas].vdate = 0;
      datetime2[seas].vtime = 0;
    }

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
  if (taxisHasBounds(taxisID4)) taxisDeleteBounds(taxisID4);
  vlistDefTaxis(vlistID4, taxisID4);

  int streamID4 = cdoStreamOpenWrite(cdoStreamName(3), cdoFiletype());
  pstreamDefVlist(streamID4, vlistID4);

  int nvars = vlistNvars(vlistID1);

  int maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recinfo(maxrecs);

  size_t gridsize = vlistGridsizeMax(vlistID1);

  Field field;
  field_init(&field);
  field.ptr = (double *) Malloc(gridsize * sizeof(double));

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID2, tsID)))
    {
      if (nrecs != cdoStreamInqTimestep(streamID3, tsID))
        cdoAbort("Number of records at time step %d of %s and %s differ!", tsID + 1, cdoGetStreamName(1).c_str(),
                 cdoGetStreamName(2).c_str());

      vdate = taxisInqVdate(taxisID2);
      vtime = taxisInqVtime(taxisID2);

      if (vdate != taxisInqVdate(taxisID3))
        cdoAbort("Verification dates at time step %d of %s and %s differ!", tsID + 1, cdoGetStreamName(1).c_str(),
                 cdoGetStreamName(2).c_str());

      if (cdoVerbose) cdoPrint("process timestep: %d %d %d", tsID + 1, vdate, vtime);

      cdiDecodeDate(vdate, &year, &month, &day);

      seas = month_to_season(month);

      set_date(vdate, vtime, &datetime2[seas]);

      if (vars1[seas] == NULL)
        {
          vars1[seas] = field_malloc(vlistID1, FIELD_PTR);
          hsets[seas] = hsetCreate(nvars);

          for (varID = 0; varID < nvars; varID++)
            {
              gridID = vlistInqVarGrid(vlistID1, varID);
              nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));

              hsetCreateVarLevels(hsets[seas], varID, nlevels, gridID);
            }
        }

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID2, &varID, &levelID);
          pstreamReadRecord(streamID2, vars1[seas][varID][levelID].ptr, &nmiss);
          vars1[seas][varID][levelID].nmiss = nmiss;
        }
      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID3, &varID, &levelID);
          pstreamReadRecord(streamID3, field.ptr, &nmiss);
          field.nmiss = nmiss;
          field.grid = vars1[seas][varID][levelID].grid;
          field.missval = vars1[seas][varID][levelID].missval;

          hsetDefVarLevelBounds(hsets[seas], varID, levelID, &vars1[seas][varID][levelID], &field);
        }

      tsID++;
    }

  tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      cdiDecodeDate(vdate, &year, &month, &day);
      if (month < 0 || month > 16) cdoAbort("Month %d out of range!", month);

      if (month <= 12)
        seas = (month % 12) / 3;
      else
        seas = month - 13;
      if (seas < 0 || seas > 3) cdoAbort("Season %d out of range!", seas + 1);

      set_date(vdate, vtime, &datetime1[seas]);

      if (vars1[seas] == NULL)
        cdoAbort("No data for season %d in %s and %s", seas, cdoGetStreamName(1).c_str(), cdoGetStreamName(2).c_str());

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);

          if (tsID == 0)
            {
              recinfo[recID].varID = varID;
              recinfo[recID].levelID = levelID;
              recinfo[recID].lconst = vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT;
            }

          pstreamReadRecord(streamID1, vars1[seas][varID][levelID].ptr, &nmiss);
          vars1[seas][varID][levelID].nmiss = nmiss;

          hsetAddVarLevelValues(hsets[seas], varID, levelID, &vars1[seas][varID][levelID]);
        }

      nsets[seas]++;
      tsID++;
    }

  int otsID = 0;
  for (seas = 0; seas < NSEAS; seas++)
    if (nsets[seas])
      {
        if (getmonthday(datetime1[seas].vdate) != getmonthday(datetime2[seas].vdate))
          cdoAbort("Verification dates for the season %d of %s and %s are different!",
                   seas, cdoGetStreamName(0).c_str(), cdoGetStreamName(1).c_str());

        for (varID = 0; varID < nvars; varID++)
          {
            if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT) continue;
            nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));

            for (levelID = 0; levelID < nlevels; levelID++)
              hsetGetVarLevelPercentiles(&vars1[seas][varID][levelID], hsets[seas], varID, levelID, pn);
          }

        taxisDefVdate(taxisID4, datetime1[seas].vdate);
        taxisDefVtime(taxisID4, datetime1[seas].vtime);
        pstreamDefTimestep(streamID4, otsID);

        for (int recID = 0; recID < maxrecs; recID++)
          {
            if (otsID && recinfo[recID].lconst) continue;
            int varID = recinfo[recID].varID;
            int levelID = recinfo[recID].levelID;

            pstreamDefRecord(streamID4, varID, levelID);
            pstreamWriteRecord(streamID4, vars1[seas][varID][levelID].ptr, vars1[seas][varID][levelID].nmiss);
          }

        otsID++;
      }

  for (seas = 0; seas < NSEAS; seas++)
    {
      if (vars1[seas] != NULL)
        {
          field_free(vars1[seas], vlistID1);
          hsetDestroy(hsets[seas]);
        }
    }

  if (field.ptr) Free(field.ptr);

  pstreamClose(streamID4);
  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

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

      Ydrunpctl    ydrunpctl         Multi-year daily running percentiles
*/

#include <cdi.h>

#include "cdo_int.h"
#include "calendar.h"
#include "pstream_int.h"
#include "percentiles_hist.h"
#include "percentiles.h"
#include "datetime.h"

#define NDAY 373

int
getmonthday(int64_t date)
{
  int year, month, day;
  cdiDecodeDate(date, &year, &month, &day);
  return month * 100 + day;
}

void *
Ydrunpctl(void *process)
{
  int varID;
  int gridID;
  int nrecs;
  int levelID;
  int inp, its;
  size_t nmiss;
  int nlevels;
  int year, month, day, dayoy;
  int64_t vdate;
  int vtime;
  int64_t vdates1[NDAY], vdates2[NDAY];
  int vtimes1[NDAY] /*, vtimes2[NDAY]*/;
  int nsets[NDAY];
  Field **vars2[NDAY];
  HISTOGRAM_SET *hsets[NDAY];

  cdoInitialize(process);
  cdoOperatorAdd("ydrunpctl", func_pctl, 0, NULL);

  operatorInputArg("percentile number, number of timesteps");
  operatorCheckArgc(2);
  double pn = parameter2double(operatorArgv()[0]);
  percentile_check_number(pn);
  int ndates = parameter2int(operatorArgv()[1]);

  for (dayoy = 0; dayoy < NDAY; dayoy++)
    {
      vars2[dayoy] = NULL;
      hsets[dayoy] = NULL;
      nsets[dayoy] = 0;
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

  int calendar = taxisInqCalendar(taxisID1);
  int dpy = calendar_dpy(calendar);

  int streamID4 = cdoStreamOpenWrite(cdoStreamName(3), cdoFiletype());
  pstreamDefVlist(streamID4, vlistID4);

  int nvars = vlistNvars(vlistID1);

  int maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recinfo(maxrecs);

  size_t gridsize = vlistGridsizeMax(vlistID1);

  Field field;
  field_init(&field);
  field.ptr = (double *) Malloc(gridsize * sizeof(double));

  std::vector<CdoDateTime> datetime(ndates + 1);

  Field ***vars1 = (Field ***) Malloc((ndates + 1) * sizeof(Field **));

  for (its = 0; its < ndates; its++)
    {
      vars1[its] = field_malloc(vlistID1, FIELD_PTR);
    }

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

      if (month >= 1 && month <= 12)
        dayoy = (month - 1) * 31 + day;
      else
        dayoy = 0;

      if (dayoy < 0 || dayoy >= NDAY) cdoAbort("Day %d out of range!", dayoy);

      vdates2[dayoy] = vdate;
      // vtimes2[dayoy] = vtime;

      if (vars2[dayoy] == NULL)
        {
          vars2[dayoy] = field_malloc(vlistID2, FIELD_PTR);
          hsets[dayoy] = hsetCreate(nvars);

          for (varID = 0; varID < nvars; varID++)
            {
              gridID = vlistInqVarGrid(vlistID2, varID);
              nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));

              hsetCreateVarLevels(hsets[dayoy], varID, nlevels, gridID);
            }
        }

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID2, &varID, &levelID);
          pstreamReadRecord(streamID2, vars2[dayoy][varID][levelID].ptr, &nmiss);
          vars2[dayoy][varID][levelID].nmiss = nmiss;
        }
      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID3, &varID, &levelID);
          pstreamReadRecord(streamID3, field.ptr, &nmiss);
          field.nmiss = nmiss;
          field.grid = vars2[dayoy][varID][levelID].grid;
          field.missval = vars2[dayoy][varID][levelID].missval;

          hsetDefVarLevelBounds(hsets[dayoy], varID, levelID, &vars2[dayoy][varID][levelID], &field);
        }

      tsID++;
    }

  for (tsID = 0; tsID < ndates; tsID++)
    {
      nrecs = cdoStreamInqTimestep(streamID1, tsID);
      if (nrecs == 0) cdoAbort("File has less then %d timesteps!", ndates);

      datetime[tsID].date = taxisInqVdate(taxisID1);
      datetime[tsID].time = taxisInqVtime(taxisID1);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);

          if (tsID == 0)
            {
              recinfo[recID].varID = varID;
              recinfo[recID].levelID = levelID;
              recinfo[recID].lconst = vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT;
            }

          pstreamReadRecord(streamID1, vars1[tsID][varID][levelID].ptr, &nmiss);
          vars1[tsID][varID][levelID].nmiss = nmiss;
        }
    }

  while (TRUE)
    {
      datetime_avg(dpy, ndates, datetime.data());

      vdate = datetime[ndates].date;
      vtime = datetime[ndates].time;

      cdiDecodeDate(vdate, &year, &month, &day);

      if (month >= 1 && month <= 12)
        dayoy = (month - 1) * 31 + day;
      else
        dayoy = 0;

      if (dayoy < 0 || dayoy >= NDAY) cdoAbort("Day %d out of range!", dayoy);

      vdates1[dayoy] = vdate;
      vtimes1[dayoy] = vtime;

      if (vars2[dayoy] == NULL)
        cdoAbort("No data for day %d in %s and %s", dayoy, cdoGetStreamName(1).c_str(), cdoGetStreamName(2).c_str());

      for (varID = 0; varID < nvars; varID++)
        {
          if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT) continue;
          nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));

          for (levelID = 0; levelID < nlevels; levelID++)
            for (inp = 0; inp < ndates; inp++) hsetAddVarLevelValues(hsets[dayoy], varID, levelID, &vars1[inp][varID][levelID]);
        }

      datetime[ndates] = datetime[0];
      vars1[ndates] = vars1[0];

      for (inp = 0; inp < ndates; inp++)
        {
          datetime[inp] = datetime[inp + 1];
          vars1[inp] = vars1[inp + 1];
        }

      nrecs = cdoStreamInqTimestep(streamID1, tsID);
      if (nrecs == 0) break;

      datetime[ndates - 1].date = taxisInqVdate(taxisID1);
      datetime[ndates - 1].time = taxisInqVtime(taxisID1);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);

          pstreamReadRecord(streamID1, vars1[ndates - 1][varID][levelID].ptr, &nmiss);
          vars1[ndates - 1][varID][levelID].nmiss = nmiss;
        }

      nsets[dayoy] += ndates;
      tsID++;
    }
  /*
  int outyear = 1e9;
  for ( dayoy = 0; dayoy < NDAY; dayoy++ )
    if ( nsets[dayoy] )
      {
        int year, month, day;
        cdiDecodeDate(vdates1[dayoy], &year, &month, &day);
        if ( year < outyear ) outyear = year;
      }

  for ( dayoy = 0; dayoy < NDAY; dayoy++ )
    if ( nsets[dayoy] )
      {
        int year, month, day;
        cdiDecodeDate(vdates1[dayoy], &year, &month, &day);
        vdates1[dayoy] = cdiEncodeDate(outyear, month, day);
      }
  */
  int otsID = 0;
  for (dayoy = 0; dayoy < NDAY; dayoy++)
    if (nsets[dayoy])
      {
        if (getmonthday(vdates1[dayoy]) != getmonthday(vdates2[dayoy]))
          cdoAbort("Verification dates for day %d of %s, %s and %s are different!", dayoy, cdoGetStreamName(0).c_str(),
                   cdoGetStreamName(1).c_str());

        for (varID = 0; varID < nvars; varID++)
          {
            if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT) continue;
            nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));

            for (levelID = 0; levelID < nlevels; levelID++)
              hsetGetVarLevelPercentiles(&vars2[dayoy][varID][levelID], hsets[dayoy], varID, levelID, pn);
          }

        taxisDefVdate(taxisID4, vdates1[dayoy]);
        taxisDefVtime(taxisID4, vtimes1[dayoy]);
        pstreamDefTimestep(streamID4, otsID);

        for (int recID = 0; recID < maxrecs; recID++)
          {
            if (otsID && recinfo[recID].lconst) continue;

            int varID = recinfo[recID].varID;
            int levelID = recinfo[recID].levelID;
            pstreamDefRecord(streamID4, varID, levelID);
            pstreamWriteRecord(streamID4, vars2[dayoy][varID][levelID].ptr, vars2[dayoy][varID][levelID].nmiss);
          }

        otsID++;
      }

  for (its = 0; its < ndates; its++)
    {
      for (varID = 0; varID < nvars; varID++)
        {
          nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
          for (levelID = 0; levelID < nlevels; levelID++) Free(vars1[its][varID][levelID].ptr);
          Free(vars1[its][varID]);
        }
      Free(vars1[its]);
    }
  Free(vars1);

  for (dayoy = 0; dayoy < NDAY; dayoy++)
    {
      if (vars2[dayoy] != NULL)
        {
          field_free(vars2[dayoy], vlistID2);
          hsetDestroy(hsets[dayoy]);
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

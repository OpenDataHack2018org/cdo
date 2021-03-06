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

      Ymonpctl   ymonpctl        Multi-year monthly percentiles
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "percentiles_hist.h"
#include "percentiles.h"

#define NMONTH 17

static int
getmonth(int64_t date)
{
  int year, month, day;
  cdiDecodeDate(date, &year, &month, &day);
  return month;
}

void *
Ymonpctl(void *process)
{
  int varID;
  int gridID;
  int64_t vdate;
  int vtime;
  int year, month, day;
  int levelID;
  size_t nmiss;
  int nrecs, nlevels;
  int64_t vdates1[NMONTH], vdates2[NMONTH];
  int vtimes1[NMONTH];
  long nsets[NMONTH];
  Field **vars1[NMONTH];
  HISTOGRAM_SET *hsets[NMONTH];

  cdoInitialize(process);
  cdoOperatorAdd("ymonpctl", func_pctl, 0, NULL);

  operatorInputArg("percentile number");
  double pn = parameter2double(operatorArgv()[0]);
  percentile_check_number(pn);

  for (month = 0; month < NMONTH; month++)
    {
      vars1[month] = NULL;
      hsets[month] = NULL;
      nsets[month] = 0;
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
      if (month < 0 || month >= NMONTH) cdoAbort("Month %d out of range!", month);

      vdates2[month] = vdate;

      if (vars1[month] == NULL)
        {
          vars1[month] = field_malloc(vlistID1, FIELD_PTR);
          hsets[month] = hsetCreate(nvars);

          for (varID = 0; varID < nvars; varID++)
            {
              gridID = vlistInqVarGrid(vlistID1, varID);
              nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));

              hsetCreateVarLevels(hsets[month], varID, nlevels, gridID);
            }
        }

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID2, &varID, &levelID);
          pstreamReadRecord(streamID2, vars1[month][varID][levelID].ptr, &nmiss);
          vars1[month][varID][levelID].nmiss = nmiss;
        }
      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID3, &varID, &levelID);
          pstreamReadRecord(streamID3, field.ptr, &nmiss);
          field.nmiss = nmiss;
          field.grid = vars1[month][varID][levelID].grid;
          field.missval = vars1[month][varID][levelID].missval;

          hsetDefVarLevelBounds(hsets[month], varID, levelID, &vars1[month][varID][levelID], &field);
        }

      tsID++;
    }

  tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      if (cdoVerbose) cdoPrint("process timestep: %d %d %d", tsID + 1, vdate, vtime);

      cdiDecodeDate(vdate, &year, &month, &day);
      if (month < 0 || month >= NMONTH) cdoAbort("Month %d out of range!", month);

      vdates1[month] = vdate;
      vtimes1[month] = vtime;

      if (vars1[month] == NULL)
        cdoAbort("No data for month %d in %s and %s", month, cdoGetStreamName(1).c_str(), cdoGetStreamName(2).c_str());

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);

          if (tsID == 0)
            {
              recinfo[recID].varID = varID;
              recinfo[recID].levelID = levelID;
              recinfo[recID].lconst = vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT;
            }

          pstreamReadRecord(streamID1, vars1[month][varID][levelID].ptr, &nmiss);
          vars1[month][varID][levelID].nmiss = nmiss;

          hsetAddVarLevelValues(hsets[month], varID, levelID, &vars1[month][varID][levelID]);
        }

      nsets[month]++;
      tsID++;
    }

  int otsID = 0;
  for (month = 0; month < NMONTH; month++)
    if (nsets[month])
      {
        if (getmonth(vdates1[month]) != getmonth(vdates2[month]))
          cdoAbort("Verification dates for the month %d of %s and %s are different!", month, cdoGetStreamName(0).c_str(),
                   cdoGetStreamName(1).c_str());

        for (varID = 0; varID < nvars; varID++)
          {
            if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT) continue;
            nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));

            for (levelID = 0; levelID < nlevels; levelID++)
              hsetGetVarLevelPercentiles(&vars1[month][varID][levelID], hsets[month], varID, levelID, pn);
          }

        taxisDefVdate(taxisID4, vdates1[month]);
        taxisDefVtime(taxisID4, vtimes1[month]);
        pstreamDefTimestep(streamID4, otsID);

        for (int recID = 0; recID < maxrecs; recID++)
          {
            if (otsID && recinfo[recID].lconst) continue;

            int varID = recinfo[recID].varID;
            int levelID = recinfo[recID].levelID;
            pstreamDefRecord(streamID4, varID, levelID);
            pstreamWriteRecord(streamID4, vars1[month][varID][levelID].ptr, vars1[month][varID][levelID].nmiss);
          }

        otsID++;
      }

  for (month = 0; month < NMONTH; month++)
    {
      if (vars1[month] != NULL)
        {
          field_free(vars1[month], vlistID1);
          hsetDestroy(hsets[month]);
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

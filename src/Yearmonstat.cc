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

      Yearmonstat   yearmonmean        Yearly mean from monthly data
      Yearmonstat   yearmonavg         Yearly average from monthly data
*/

#include <cdi.h>

#include "cdo_int.h"
#include "calendar.h"
#include "pstream_int.h"
#include "datetime.h"

void *
Yearmonstat(void *process)
{
  TimeStat timestat_date = TimeStat::MEAN;
  int64_t vdate = 0, vdate0 = 0;
  int vtime = 0, vtime0 = 0;
  int nrecs;
  int varID, levelID;
  int dpm;
  int year0 = 0, month0 = 0;
  int year, month, day;
  size_t nmiss;
  int nlevel;
  long nsets;
  double dsets;
  char vdatestr[32], vtimestr[32];

  cdoInitialize(process);

  // clang-format off
  cdoOperatorAdd("yearmonmean",  func_mean, 0, NULL);
  cdoOperatorAdd("yearmonavg",   func_avg,  0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  taxisWithBounds(taxisID2);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int nvars = vlistNvars(vlistID1);

  int maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recinfo(maxrecs);

  int calendar = taxisInqCalendar(taxisID1);
  DateTimeList dtlist;
  dtlist.setStat(timestat_date);
  dtlist.setCalendar(calendar);

  size_t gridsize = vlistGridsizeMax(vlistID1);

  Field field;
  field_init(&field);
  field.ptr = (double *) Malloc(gridsize * sizeof(double));

  Field **vars1 = field_malloc(vlistID1, FIELD_PTR);
  Field **samp1 = field_malloc(vlistID1, FIELD_NONE);

  int tsID = 0;
  int otsID = 0;
  while (TRUE)
    {
      nsets = 0;
      dsets = 0;
      while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
        {
          dtlist.taxisInqTimestep(taxisID1, nsets);
          vdate = dtlist.getVdate(nsets);
          vtime = dtlist.getVtime(nsets);
          cdiDecodeDate(vdate, &year, &month, &day);

          if (nsets == 0) year0 = year;

          if (year != year0) break;

          if (nsets > 0 && month == month0)
            {
              date2str(vdate0, vdatestr, sizeof(vdatestr));
              time2str(vtime0, vtimestr, sizeof(vtimestr));
              cdoWarning("   last timestep: %s %s", vdatestr, vtimestr);
              date2str(vdate, vdatestr, sizeof(vdatestr));
              time2str(vtime, vtimestr, sizeof(vtimestr));
              cdoWarning("current timestep: %s %s", vdatestr, vtimestr);
              cdoAbort("Month does not change!");
            }

          dpm = days_per_month(calendar, year, month);

          for (int recID = 0; recID < nrecs; recID++)
            {
              pstreamInqRecord(streamID1, &varID, &levelID);

              if (tsID == 0)
                {
                  recinfo[recID].varID = varID;
                  recinfo[recID].levelID = levelID;
                  recinfo[recID].lconst = vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT;
                }

              gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

              if (nsets == 0)
                {
                  pstreamReadRecord(streamID1, vars1[varID][levelID].ptr, &nmiss);
                  vars1[varID][levelID].nmiss = nmiss;

                  farcmul(&vars1[varID][levelID], dpm);

                  if (nmiss > 0 || samp1[varID][levelID].ptr)
                    {
                      if (samp1[varID][levelID].ptr == NULL)
                        samp1[varID][levelID].ptr = (double *) Malloc(gridsize * sizeof(double));

                      for (size_t i = 0; i < gridsize; i++)
                        if (DBL_IS_EQUAL(vars1[varID][levelID].ptr[i], vars1[varID][levelID].missval))
                          samp1[varID][levelID].ptr[i] = 0;
                        else
                          samp1[varID][levelID].ptr[i] = dpm;
                    }
                }
              else
                {
                  pstreamReadRecord(streamID1, field.ptr, &nmiss);
                  field.nmiss = nmiss;
                  field.grid = vars1[varID][levelID].grid;
                  field.missval = vars1[varID][levelID].missval;

                  farcmul(&field, dpm);

                  if (field.nmiss > 0 || samp1[varID][levelID].ptr)
                    {
                      if (samp1[varID][levelID].ptr == NULL)
                        {
                          samp1[varID][levelID].ptr = (double *) Malloc(gridsize * sizeof(double));
                          for (size_t i = 0; i < gridsize; i++) samp1[varID][levelID].ptr[i] = dsets;
                        }

                      for (size_t i = 0; i < gridsize; i++)
                        if (!DBL_IS_EQUAL(field.ptr[i], vars1[varID][levelID].missval)) samp1[varID][levelID].ptr[i] += dpm;
                    }

                  farfun(&vars1[varID][levelID], field, operfunc);
                }
            }

          month0 = month;
          vdate0 = vdate;
          vtime0 = vtime;
          nsets++;
          dsets += dpm;
          tsID++;
        }

      if (nrecs == 0 && nsets == 0) break;

      for (varID = 0; varID < nvars; varID++)
        {
          if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT) continue;
          nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
          for (levelID = 0; levelID < nlevel; levelID++)
            {
              if (samp1[varID][levelID].ptr == NULL)
                farcdiv(&vars1[varID][levelID], dsets);
              else
                fardiv(&vars1[varID][levelID], samp1[varID][levelID]);
            }
        }

      if (cdoVerbose)
        {
          date2str(vdate0, vdatestr, sizeof(vdatestr));
          time2str(vtime0, vtimestr, sizeof(vtimestr));
          cdoPrint("%s %s  nsets = %ld", vdatestr, vtimestr, nsets);
        }

      dtlist.statTaxisDefTimestep(taxisID2, nsets);
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
  field_free(samp1, vlistID1);

  if (field.ptr) Free(field.ptr);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

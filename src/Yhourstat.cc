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

      Yhourstat   yhourrange       Multi-year hourly range
      Yhourstat   yhourmin         Multi-year hourly minimum
      Yhourstat   yhourmax         Multi-year hourly maximum
      Yhourstat   yhoursum         Multi-year hourly sum
      Yhourstat   yhourmean        Multi-year hourly mean
      Yhourstat   yhouravg         Multi-year hourly average
      Yhourstat   yhourvar         Multi-year hourly variance
      Yhourstat   yhourvar1        Multi-year hourly variance [Normalize by (n-1)]
      Yhourstat   yhourstd         Multi-year hourly standard deviation
      Yhourstat   yhourstd1        Multi-year hourly standard deviation [Normalize by (n-1)]
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"


static int
hour_of_year(int64_t vdate, int vtime)
{
  int year, month, day, houroy;
  int hour, minute, second;

  cdiDecodeDate(vdate, &year, &month, &day);
  cdiDecodeTime(vtime, &hour, &minute, &second);

  if (month >= 1 && month <= 12 && day >= 1 && day <= 31 && hour >= 0 && hour < 24)
    houroy = ((month - 1) * 31 + day - 1) * 25 + hour + 1;
  else
    houroy = 0;

  return houroy;
}

static int
hour_of_day(int64_t vdate, int vtime)
{
  int year, month, day, houroy;
  int hour, minute, second;

  cdiDecodeDate(vdate, &year, &month, &day);
  cdiDecodeTime(vtime, &hour, &minute, &second);

  if (month >= 1 && month <= 12 && day >= 1 && day <= 31 && hour >= 0 && hour < 24)
    houroy = hour + 1;
  else
    houroy = 0;

  return houroy;
}

void *
Yhourstat(void *process)
{
  int varID;
  int nrecs;
  int levelID;
  size_t nmiss;

  cdoInitialize(process);

  // clang-format off
  cdoOperatorAdd("yhourrange", func_range, 0, NULL);
  cdoOperatorAdd("yhourmin",   func_min,   0, NULL);
  cdoOperatorAdd("yhourmax",   func_max,   0, NULL);
  cdoOperatorAdd("yhoursum",   func_sum,   0, NULL);
  cdoOperatorAdd("yhourmean",  func_mean,  0, NULL);
  cdoOperatorAdd("yhouravg",   func_avg,   0, NULL);
  cdoOperatorAdd("yhourvar",   func_var,   0, NULL);
  cdoOperatorAdd("yhourvar1",  func_var1,  0, NULL);
  cdoOperatorAdd("yhourstd",   func_std,   0, NULL);
  cdoOperatorAdd("yhourstd1",  func_std1,  0, NULL);
  cdoOperatorAdd("dhourrange", func_range, 0, NULL);
  cdoOperatorAdd("dhourmin",   func_min,   1, NULL);
  cdoOperatorAdd("dhourmax",   func_max,   1, NULL);
  cdoOperatorAdd("dhoursum",   func_sum,   1, NULL);
  cdoOperatorAdd("dhourmean",  func_mean,  1, NULL);
  cdoOperatorAdd("dhouravg",   func_avg,   1, NULL);
  cdoOperatorAdd("dhourvar",   func_var,   1, NULL);
  cdoOperatorAdd("dhourvar1",  func_var1,  1, NULL);
  cdoOperatorAdd("dhourstd",   func_std,   1, NULL);
  cdoOperatorAdd("dhourstd1",  func_std1,  1, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  bool ldaily  = cdoOperatorF2(operatorID) == 1;

  bool lrange  = operfunc == func_range;
  bool lmean   = operfunc == func_mean || operfunc == func_avg;
  bool lstd    = operfunc == func_std || operfunc == func_std1;
  bool lvarstd = operfunc == func_std || operfunc == func_var || operfunc == func_std1 || operfunc == func_var1;
  int  divisor = operfunc == func_std1 || operfunc == func_var1;
  // clang-format on

  int maxHours = ldaily ? 25 : 9301; // 31*12*25 + 1
  std::vector<int> hourot_nsets(maxHours); // hour of time
  std::vector<int> vtimes(maxHours);
  std::vector<int64_t> vdates(maxHours);
  std::vector<Field **> vars1(maxHours), vars2(maxHours), samp1(maxHours);

  for (int hourot = 0; hourot < maxHours; ++hourot)
    {
      vars1[hourot] = NULL;
      vars2[hourot] = NULL;
      samp1[hourot] = NULL;
      hourot_nsets[hourot] = 0;
    }

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  if (taxisHasBounds(taxisID2)) taxisDeleteBounds(taxisID2);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recinfo(maxrecs);

  size_t gridsizemax = vlistGridsizeMax(vlistID1);

  Field field;
  field_init(&field);
  field.ptr = (double *) Malloc(gridsizemax * sizeof(double));

  int tsID = 0;
  int otsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      int64_t vdate = taxisInqVdate(taxisID1);
      int vtime = taxisInqVtime(taxisID1);

      if (cdoVerbose) cdoPrint("process timestep: %d %d %d", tsID + 1, vdate, vtime);

      int hourot = ldaily ? hour_of_day(vdate, vtime) : hour_of_year(vdate, vtime);
      if (hourot < 0 || hourot >= maxHours)
        {
          char vdatestr[32], vtimestr[32];
          date2str(vdate, vdatestr, sizeof(vdatestr));
          time2str(vtime, vtimestr, sizeof(vtimestr));
          cdoAbort("Hour of year %d out of range (%s %s)!", hourot, vdatestr, vtimestr);
        }

      vdates[hourot] = vdate;
      vtimes[hourot] = vtime;

      if (vars1[hourot] == NULL)
        {
          vars1[hourot] = field_malloc(vlistID1, FIELD_PTR);
          samp1[hourot] = field_malloc(vlistID1, FIELD_NONE);
          if (lvarstd || lrange) vars2[hourot] = field_malloc(vlistID1, FIELD_PTR);
        }

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);

          if (tsID == 0)
            {
              recinfo[recID].varID = varID;
              recinfo[recID].levelID = levelID;
              recinfo[recID].lconst = vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT;
            }

          Field *psamp1 = &samp1[hourot][varID][levelID];
          Field *pvars1 = &vars1[hourot][varID][levelID];
          Field *pvars2 = vars2[hourot] ? &vars2[hourot][varID][levelID] : NULL;
          int nsets = hourot_nsets[hourot];

          size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

          if (nsets == 0)
            {
              pstreamReadRecord(streamID1, pvars1->ptr, &nmiss);
              pvars1->nmiss = nmiss;
              if (lrange)
                {
                  pvars2->nmiss = pvars1->nmiss;
                  for (size_t i = 0; i < gridsize; i++) pvars2->ptr[i] = pvars1->ptr[i];
                }

              if (nmiss > 0 || psamp1->ptr)
                {
                  if (psamp1->ptr == NULL) psamp1->ptr = (double *) Malloc(gridsize * sizeof(double));

                  for (size_t i = 0; i < gridsize; i++) psamp1->ptr[i] = !DBL_IS_EQUAL(pvars1->ptr[i], pvars1->missval);
                }
            }
          else
            {
              pstreamReadRecord(streamID1, field.ptr, &nmiss);
              field.nmiss = nmiss;
              field.grid = pvars1->grid;
              field.missval = pvars1->missval;

              if (field.nmiss > 0 || psamp1->ptr)
                {
                  if (psamp1->ptr == NULL)
                    {
                      psamp1->ptr = (double *) Malloc(gridsize * sizeof(double));
                      for (size_t i = 0; i < gridsize; i++) psamp1->ptr[i] = nsets;
                    }

                  for (size_t i = 0; i < gridsize; i++)
                    if (!DBL_IS_EQUAL(field.ptr[i], pvars1->missval)) psamp1->ptr[i]++;
                }

              if (lvarstd)
                {
                  farsumq(pvars2, field);
                  farsum(pvars1, field);
                }
              else if (lrange)
                {
                  farmin(pvars2, field);
                  farmax(pvars1, field);
                }
              else
                {
                  farfun(pvars1, field, operfunc);
                }
            }
        }

      if (hourot_nsets[hourot] == 0 && lvarstd)
        for (int recID = 0; recID < maxrecs; recID++)
          {
            if (recinfo[recID].lconst) continue;

            int varID = recinfo[recID].varID;
            int levelID = recinfo[recID].levelID;
            Field *pvars1 = &vars1[hourot][varID][levelID];
            Field *pvars2 = &vars2[hourot][varID][levelID];

            farmoq(pvars2, *pvars1);
          }

      hourot_nsets[hourot]++;
      tsID++;
    }

  for (int hourot = 0; hourot < maxHours; ++hourot)
    if (hourot_nsets[hourot])
      {
        int nsets = hourot_nsets[hourot];
        for (int recID = 0; recID < maxrecs; recID++)
          {
            if (recinfo[recID].lconst) continue;

            int varID = recinfo[recID].varID;
            int levelID = recinfo[recID].levelID;
            Field *psamp1 = &samp1[hourot][varID][levelID];
            Field *pvars1 = &vars1[hourot][varID][levelID];
            Field *pvars2 = vars2[hourot] ? &vars2[hourot][varID][levelID] : NULL;

            if (lmean)
              {
                if (psamp1->ptr)
                  fardiv(pvars1, *psamp1);
                else
                  farcdiv(pvars1, (double) nsets);
              }
            else if (lvarstd)
              {
                if (psamp1->ptr)
                  {
                    if (lstd)
                      farstd(pvars1, *pvars2, *psamp1, divisor);
                    else
                      farvar(pvars1, *pvars2, *psamp1, divisor);
                  }
                else
                  {
                    if (lstd)
                      farcstd(pvars1, *pvars2, nsets, divisor);
                    else
                      farcvar(pvars1, *pvars2, nsets, divisor);
                  }
              }
            else if (lrange)
              {
                farsub(pvars1, *pvars2);
              }
          }

        taxisDefVdate(taxisID2, vdates[hourot]);
        taxisDefVtime(taxisID2, vtimes[hourot]);
        pstreamDefTimestep(streamID2, otsID);

        for (int recID = 0; recID < maxrecs; recID++)
          {
            if (otsID && recinfo[recID].lconst) continue;

            int varID = recinfo[recID].varID;
            int levelID = recinfo[recID].levelID;
            Field *pvars1 = &vars1[hourot][varID][levelID];

            pstreamDefRecord(streamID2, varID, levelID);
            pstreamWriteRecord(streamID2, pvars1->ptr, pvars1->nmiss);
          }

        otsID++;
      }

  for (int hourot = 0; hourot < maxHours; ++hourot)
    {
      if (vars1[hourot] != NULL)
        {
          field_free(samp1[hourot], vlistID1);
          field_free(vars1[hourot], vlistID1);
          if (lvarstd) field_free(vars2[hourot], vlistID1);
        }
    }

  if (field.ptr) Free(field.ptr);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

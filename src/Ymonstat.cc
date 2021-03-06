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

      Ymonstat   ymonrange       Multi-year monthly range
      Ymonstat   ymonmin         Multi-year monthly minimum
      Ymonstat   ymonmax         Multi-year monthly maximum
      Ymonstat   ymonsum         Multi-year monthly sum
      Ymonstat   ymonmean        Multi-year monthly mean
      Ymonstat   ymonavg         Multi-year monthly average
      Ymonstat   ymonvar         Multi-year monthly variance
      Ymonstat   ymonvar1        Multi-year monthly variance [Normalize by (n-1)]
      Ymonstat   ymonstd         Multi-year monthly standard deviation
      Ymonstat   ymonstd1        Multi-year monthly standard deviation [Normalize by (n-1)]
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

#define NMONTH 17

/*
static
int compareInt(const void *a, const void *b)
{
  const int *x = (const int*) a;
  const int *y = (const int*) b;
  return ((*x > *y) - (*x < *y)) * 2 + (*x > *y) - (*x < *y);
}
*/

void *
Ymonstat(void *process)
{
  int varID;
  int year, month, day;
  int nrecs;
  int levelID;
  int month_nsets[NMONTH];
  size_t nmiss;
  int64_t vdates[NMONTH];
  int vtimes[NMONTH];
  int mon[NMONTH];
  int nmon = 0;
  Field **vars1[NMONTH], **vars2[NMONTH], **samp1[NMONTH];
  Field field;

  cdoInitialize(process);

  // clang-format off
  cdoOperatorAdd("ymonrange", func_range, 0, NULL);
  cdoOperatorAdd("ymonmin",   func_min,   0, NULL);
  cdoOperatorAdd("ymonmax",   func_max,   0, NULL);
  cdoOperatorAdd("ymonsum",   func_sum,   0, NULL);
  cdoOperatorAdd("ymonmean",  func_mean,  0, NULL);
  cdoOperatorAdd("ymonavg",   func_avg,   0, NULL);
  cdoOperatorAdd("ymonvar",   func_var,   0, NULL);
  cdoOperatorAdd("ymonvar1",  func_var1,  0, NULL);
  cdoOperatorAdd("ymonstd",   func_std,   0, NULL);
  cdoOperatorAdd("ymonstd1",  func_std1,  0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  bool lrange  = operfunc == func_range;
  bool lmean   = operfunc == func_mean || operfunc == func_avg;
  bool lstd    = operfunc == func_std || operfunc == func_std1;
  bool lvarstd = operfunc == func_std || operfunc == func_var || operfunc == func_std1 || operfunc == func_var1;
  int  divisor = operfunc == func_std1 || operfunc == func_var1;
  // clang-format on

  for (month = 0; month < NMONTH; month++)
    {
      vars1[month] = NULL;
      vars2[month] = NULL;
      samp1[month] = NULL;
      month_nsets[month] = 0;
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
  field_init(&field);
  field.ptr = (double *) Malloc(gridsizemax * sizeof(double));

  int tsID = 0;
  int otsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      int64_t vdate = taxisInqVdate(taxisID1);
      int vtime = taxisInqVtime(taxisID1);

      if (cdoVerbose) cdoPrint("process timestep: %d %d %d", tsID + 1, vdate, vtime);

      cdiDecodeDate(vdate, &year, &month, &day);
      if (month < 0 || month >= NMONTH) cdoAbort("month %d out of range!", month);

      vdates[month] = vdate;
      vtimes[month] = vtime;
      // mon[month] = vdate;

      if (vars1[month] == NULL)
        {
          mon[nmon++] = month;
          vars1[month] = field_malloc(vlistID1, FIELD_PTR);
          samp1[month] = field_malloc(vlistID1, FIELD_NONE);
          if (lvarstd || lrange) vars2[month] = field_malloc(vlistID1, FIELD_PTR);
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

          Field *psamp1 = &samp1[month][varID][levelID];
          Field *pvars1 = &vars1[month][varID][levelID];
          Field *pvars2 = vars2[month] ? &vars2[month][varID][levelID] : NULL;
          int nsets = month_nsets[month];

          size_t gridsize = pvars1->size;

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

      if (month_nsets[month] == 0 && lvarstd)
        for (int recID = 0; recID < maxrecs; recID++)
          {
            if (recinfo[recID].lconst) continue;

            int varID = recinfo[recID].varID;
            int levelID = recinfo[recID].levelID;
            Field *pvars1 = &vars1[month][varID][levelID];
            Field *pvars2 = &vars2[month][varID][levelID];

            farmoq(pvars2, *pvars1);
          }

      month_nsets[month]++;
      tsID++;
    }

  if (nmon == 12)
    {
      int smon = 0;
      for (month = 1; month <= 12; month++)
        if (month_nsets[month]) smon++;
      if (smon == 12)
        for (month = 1; month <= 12; month++) mon[month - 1] = month;
    }

  /* sort output time steps */
  /*
  nmon = 0;
  for ( month = 0; month < NMONTH; month++ )
    {
      if ( month_nsets[month] == 0 )
        for ( i = month+1; i < NMONTH; i++ ) mon[i-1] = mon[i];
      else
        nmon++;
    }

  qsort(mon, nmon, sizeof(int), compareInt);

  for ( i = 0; i < nmon; i++ )
    {
      cdiDecodeDate(mon[i], &year, &month, &day);
      mon[i] = month;
    }
  */
  for (int i = 0; i < nmon; i++)
    {
      month = mon[i];
      int nsets = month_nsets[month];
      if (nsets == 0) cdoAbort("Internal problem, nsets[%d] not defined!", month);

      for (int recID = 0; recID < maxrecs; recID++)
        {
          if (recinfo[recID].lconst) continue;

          int varID = recinfo[recID].varID;
          int levelID = recinfo[recID].levelID;
          Field *psamp1 = &samp1[month][varID][levelID];
          Field *pvars1 = &vars1[month][varID][levelID];
          Field *pvars2 = vars2[month] ? &vars2[month][varID][levelID] : NULL;

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

      taxisDefVdate(taxisID2, vdates[month]);
      taxisDefVtime(taxisID2, vtimes[month]);
      pstreamDefTimestep(streamID2, otsID);

      for (int recID = 0; recID < maxrecs; recID++)
        {
          if (otsID && recinfo[recID].lconst) continue;

          int varID = recinfo[recID].varID;
          int levelID = recinfo[recID].levelID;
          Field *pvars1 = &vars1[month][varID][levelID];

          pstreamDefRecord(streamID2, varID, levelID);
          pstreamWriteRecord(streamID2, pvars1->ptr, pvars1->nmiss);
        }

      otsID++;
    }

  for (month = 0; month < NMONTH; month++)
    {
      if (vars1[month] != NULL)
        {
          field_free(vars1[month], vlistID1);
          field_free(samp1[month], vlistID1);
          if (lvarstd) field_free(vars2[month], vlistID1);
        }
    }

  if (field.ptr) Free(field.ptr);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

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

      Timstat    timrange        Time range
      Timstat    timmin          Time minimum
      Timstat    timmax          Time maximum
      Timstat    timsum          Time sum
      Timstat    timmean         Time mean
      Timstat    timavg          Time average
      Timstat    timvar          Time variance
      Timstat    timvar1         Time variance [Normalize by (n-1)]
      Timstat    timstd          Time standard deviation
      Timstat    timstd1         Time standard deviation [Normalize by (n-1)]
      Hourstat   hourrange       Hourly range
      Hourstat   hourmin         Hourly minimum
      Hourstat   hourmax         Hourly maximum
      Hourstat   hoursum         Hourly sum
      Hourstat   hourmean        Hourly mean
      Hourstat   houravg         Hourly average
      Hourstat   hourvar         Hourly variance
      Hourstat   hourvar1        Hourly variance [Normalize by (n-1)]
      Hourstat   hourstd         Hourly standard deviation
      Hourstat   hourstd1        Hourly standard deviation [Normalize by (n-1)]
      Daystat    dayrange        Daily range
      Daystat    daymin          Daily minimum
      Daystat    daymax          Daily maximum
      Daystat    daysum          Daily sum
      Daystat    daymean         Daily mean
      Daystat    dayavg          Daily average
      Daystat    dayvar          Daily variance
      Daystat    dayvar1         Daily variance [Normalize by (n-1)]
      Daystat    daystd          Daily standard deviation
      Daystat    daystd1         Daily standard deviation [Normalize by (n-1)]
      Monstat    monrange        Monthly range
      Monstat    monmin          Monthly minimum
      Monstat    monmax          Monthly maximum
      Monstat    monsum          Monthly sum
      Monstat    monmean         Monthly mean
      Monstat    monavg          Monthly average
      Monstat    monvar          Monthly variance
      Monstat    monvar1         Monthly variance [Normalize by (n-1)]
      Monstat    monstd          Monthly standard deviation
      Monstat    monstd1         Monthly standard deviation [Normalize by (n-1)]
      Yearstat   yearrange       Yearly range
      Yearstat   yearmin         Yearly minimum
      Yearstat   yearmax         Yearly maximum
      Yearstat   yearsum         Yearly sum
      Yearstat   yearmean        Yearly mean
      Yearstat   yearavg         Yearly average
      Yearstat   yearvar         Yearly variance
      Yearstat   yearvar1        Yearly variance [Normalize by (n-1)]
      Yearstat   yearstd         Yearly standard deviation
      Yearstat   yearstd1        Yearly standard deviation [Normalize by (n-1)]
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "datetime.h"

enum
{
  HOUR_LEN = 4,
  DAY_LEN = 6,
  MON_LEN = 8,
  YEAR_LEN = 10
};

static void
timstatAddOperators(void)
{
  // clang-format off
  cdoOperatorAdd("timrange",   func_range,  DATE_LEN, NULL);
  cdoOperatorAdd("timmin",     func_min,    DATE_LEN, NULL);
  cdoOperatorAdd("timmax",     func_max,    DATE_LEN, NULL);
  cdoOperatorAdd("timsum",     func_sum,    DATE_LEN, NULL);
  cdoOperatorAdd("timmean",    func_mean,   DATE_LEN, NULL);
  cdoOperatorAdd("timavg",     func_avg,    DATE_LEN, NULL);
  cdoOperatorAdd("timvar",     func_var,    DATE_LEN, NULL);
  cdoOperatorAdd("timvar1",    func_var1,   DATE_LEN, NULL);
  cdoOperatorAdd("timstd",     func_std,    DATE_LEN, NULL);
  cdoOperatorAdd("timstd1",    func_std1,   DATE_LEN, NULL);
  cdoOperatorAdd("yearrange",  func_range,  YEAR_LEN, NULL);
  cdoOperatorAdd("yearmin",    func_min,    YEAR_LEN, NULL);
  cdoOperatorAdd("yearmax",    func_max,    YEAR_LEN, NULL);
  cdoOperatorAdd("yearminidx", func_minidx, YEAR_LEN, NULL);
  cdoOperatorAdd("yearmaxidx", func_maxidx, YEAR_LEN, NULL);
  cdoOperatorAdd("yearsum",    func_sum,    YEAR_LEN, NULL);
  cdoOperatorAdd("yearmean",   func_mean,   YEAR_LEN, NULL);
  cdoOperatorAdd("yearavg",    func_avg,    YEAR_LEN, NULL);
  cdoOperatorAdd("yearvar",    func_var,    YEAR_LEN, NULL);
  cdoOperatorAdd("yearvar1",   func_var1,   YEAR_LEN, NULL);
  cdoOperatorAdd("yearstd",    func_std,    YEAR_LEN, NULL);
  cdoOperatorAdd("yearstd1",   func_std1,   YEAR_LEN, NULL);
  cdoOperatorAdd("monrange",   func_range,  MON_LEN, NULL);
  cdoOperatorAdd("monmin",     func_min,    MON_LEN, NULL);
  cdoOperatorAdd("monmax",     func_max,    MON_LEN, NULL);
  cdoOperatorAdd("monsum",     func_sum,    MON_LEN, NULL);
  cdoOperatorAdd("monmean",    func_mean,   MON_LEN, NULL);
  cdoOperatorAdd("monavg",     func_avg,    MON_LEN, NULL);
  cdoOperatorAdd("monvar",     func_var,    MON_LEN, NULL);
  cdoOperatorAdd("monvar1",    func_var1,   MON_LEN, NULL);
  cdoOperatorAdd("monstd",     func_std,    MON_LEN, NULL);
  cdoOperatorAdd("monstd1",    func_std1,   MON_LEN, NULL);
  cdoOperatorAdd("dayrange",   func_range,  DAY_LEN, NULL);
  cdoOperatorAdd("daymin",     func_min,    DAY_LEN, NULL);
  cdoOperatorAdd("daymax",     func_max,    DAY_LEN, NULL);
  cdoOperatorAdd("daysum",     func_sum,    DAY_LEN, NULL);
  cdoOperatorAdd("daymean",    func_mean,   DAY_LEN, NULL);
  cdoOperatorAdd("dayavg",     func_avg,    DAY_LEN, NULL);
  cdoOperatorAdd("dayvar",     func_var,    DAY_LEN, NULL);
  cdoOperatorAdd("dayvar1",    func_var1,   DAY_LEN, NULL);
  cdoOperatorAdd("daystd",     func_std,    DAY_LEN, NULL);
  cdoOperatorAdd("daystd1",    func_std1,   DAY_LEN, NULL);
  cdoOperatorAdd("hourrange",  func_range,  HOUR_LEN, NULL);
  cdoOperatorAdd("hourmin",    func_min,    HOUR_LEN, NULL);
  cdoOperatorAdd("hourmax",    func_max,    HOUR_LEN, NULL);
  cdoOperatorAdd("hoursum",    func_sum,    HOUR_LEN, NULL);
  cdoOperatorAdd("hourmean",   func_mean,   HOUR_LEN, NULL);
  cdoOperatorAdd("houravg",    func_avg,    HOUR_LEN, NULL);
  cdoOperatorAdd("hourvar",    func_var,    HOUR_LEN, NULL);
  cdoOperatorAdd("hourvar1",   func_var1,   HOUR_LEN, NULL);
  cdoOperatorAdd("hourstd",    func_std,    HOUR_LEN, NULL);
  cdoOperatorAdd("hourstd1",   func_std1,   HOUR_LEN, NULL);
  // clang-format on
}


void *Timstat(void *argument)
{
  TimeStat timestat_date = TimeStat::MEAN;
  int64_t vdate0 = 0;
  int vtime0 = 0;
  int nrecs;
  int varID, levelID;
  int streamID3 = -1;
  int vlistID3, taxisID3 = -1;
  size_t nmiss;
  bool lvfrac = false;
  int nwpv; // number of words per value; real:1  complex:2
  char indate1[DATE_LEN+1], indate2[DATE_LEN+1];
  double vfrac = 1;

  cdoInitialize(argument);

  timstatAddOperators();

  int operatorID = cdoOperatorID();
  int operfunc   = cdoOperatorF1(operatorID);
  int comparelen = cdoOperatorF2(operatorID);

  bool lrange  = operfunc == func_range;
  bool lminidx = operfunc == func_minidx;
  bool lmaxidx = operfunc == func_maxidx;
  bool lmean   = operfunc == func_mean || operfunc == func_avg;
  bool lstd    = operfunc == func_std || operfunc == func_std1;
  bool lvarstd = operfunc == func_std || operfunc == func_var || operfunc == func_std1 || operfunc == func_var1;
  int  divisor = operfunc == func_std1 || operfunc == func_var1;

  if (operfunc == func_mean)
    {
      int oargc = operatorArgc();
      char **oargv = operatorArgv();

      if (oargc == 1)
        {
          lvfrac = true;
          vfrac = atof(oargv[0]);
          if (cdoVerbose) cdoPrint("Set vfrac to %g", vfrac);
          if (vfrac < 0 || vfrac > 1) cdoAbort("vfrac out of range!");
        }
      else if (oargc > 1)
        cdoAbort("Too many arguments!");
    }

  int cmplen = DATE_LEN - comparelen;

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  if (cmplen == 0) vlistDefNtsteps(vlistID2, 1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  taxisWithBounds(taxisID2);
  if (taxisInqType(taxisID2) == TAXIS_FORECAST) taxisDefType(taxisID2, TAXIS_RELATIVE);
  vlistDefTaxis(vlistID2, taxisID2);

  int nvars = vlistNvars(vlistID1);

  const char *freq = NULL;
  if (comparelen == DAY_LEN)
    freq = "day";
  else if (comparelen == MON_LEN)
    freq = "mon";
  else if (comparelen == YEAR_LEN)
    freq = "year";
  if (freq) cdiDefAttTxt(vlistID2, CDI_GLOBAL, "frequency", (int) strlen(freq), freq);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  if (cdoDiag)
    {
      char filename[8192];

      strcpy(filename, cdoOperatorName(operatorID));
      strcat(filename, "_");
      strcat(filename, cdoGetStreamName(1).c_str());
      streamID3 = cdoStreamOpenWrite(filename, cdoFiletype());

      vlistID3 = vlistDuplicate(vlistID1);

      for (varID = 0; varID < nvars; ++varID)
        {
          vlistDefVarDatatype(vlistID3, varID, CDI_DATATYPE_INT32);
          vlistDefVarMissval(vlistID3, varID, -1);
          vlistDefVarUnits(vlistID3, varID, "");
          vlistDefVarAddoffset(vlistID3, varID, 0);
          vlistDefVarScalefactor(vlistID3, varID, 1);
        }

      taxisID3 = taxisDuplicate(taxisID1);
      taxisWithBounds(taxisID3);
      vlistDefTaxis(vlistID3, taxisID3);

      pstreamDefVlist(streamID3, vlistID3);
    }

  int maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recinfo(maxrecs);

  DateTimeList *dtlist = dtlist_new();
  dtlist_set_stat(dtlist, timestat_date);
  dtlist_set_calendar(dtlist, taxisInqCalendar(taxisID1));

  size_t gridsizemax = vlistGridsizeMax(vlistID1);
  if (vlistNumber(vlistID1) != CDI_REAL) gridsizemax *= 2;

  int FIELD_MEMTYPE = 0;
  if (CDO_Memtype == MEMTYPE_FLOAT) FIELD_MEMTYPE = MEMTYPE_FLOAT;

  Field field;
  field_init(&field);
  field.memtype = FIELD_MEMTYPE;
  if (FIELD_MEMTYPE == MEMTYPE_FLOAT)
    field.ptrf = (float *) Malloc(gridsizemax * sizeof(float));
  else
    field.ptr = (double *) Malloc(gridsizemax * sizeof(double));

  Field **samp1 = field_malloc(vlistID1, FIELD_NONE);
  Field **vars1 = field_malloc(vlistID1, FIELD_PTR);
  Field **vars2 = NULL;
  if (lvarstd || lrange || lminidx || lmaxidx) vars2 = field_malloc(vlistID1, FIELD_PTR);

  int tsID = 0;
  int otsID = 0;
  while (TRUE)
    {
      int nsets = 0;
      while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
        {
          dtlist_taxisInqTimestep(dtlist, taxisID1, nsets);
          int64_t vdate = dtlist_get_vdate(dtlist, nsets);
          int vtime = dtlist_get_vtime(dtlist, nsets);

          if (nsets == 0) SET_DATE(indate2, vdate, vtime);
          SET_DATE(indate1, vdate, vtime);

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

              Field *psamp1 = &samp1[varID][levelID];
              Field *pvars1 = &vars1[varID][levelID];
              Field *pvars2 = vars2 ? &vars2[varID][levelID] : NULL;

              nwpv = pvars1->nwpv;
              size_t gridsize = pvars1->size;

              if (nsets == 0)
                {
                  pstreamReadRecord(streamID1, pvars1->ptr, &nmiss);
                  pvars1->nmiss = nmiss;
                  if (lrange)
                    {
                      pvars2->nmiss = nmiss;
                      for (size_t i = 0; i < nwpv * gridsize; i++) pvars2->ptr[i] = pvars1->ptr[i];
                    }
                  else if (lminidx || lmaxidx)
                    {
                      pvars2->nmiss = nmiss;
                      for (size_t i = 0; i < nwpv * gridsize; i++) pvars2->ptr[i] = pvars1->ptr[i];
                      for (size_t i = 0; i < nwpv * gridsize; i++) pvars1->ptr[i] = 0;
                    }

                  if (nmiss > 0 || psamp1->ptr)
                    {
                      if (psamp1->ptr == NULL) psamp1->ptr = (double *) Malloc(nwpv * gridsize * sizeof(double));

                      for (size_t i = 0; i < nwpv * gridsize; i++)
                        psamp1->ptr[i] = !DBL_IS_EQUAL(pvars1->ptr[i], pvars1->missval);
                    }
                }
              else
                {
                  if (CDO_Memtype == MEMTYPE_FLOAT)
                    pstreamReadRecordF(streamID1, field.ptrf, &nmiss);
                  else
                    pstreamReadRecord(streamID1, field.ptr, &nmiss);
                  field.nmiss = nmiss;
                  field.size = gridsize;
                  field.grid = pvars1->grid;
                  field.missval = pvars1->missval;
                  if (field.nmiss > 0 || psamp1->ptr)
                    {
                      if (psamp1->ptr == NULL)
                        {
                          psamp1->ptr = (double *) Malloc(nwpv * gridsize * sizeof(double));
                          for (size_t i = 0; i < nwpv * gridsize; i++) psamp1->ptr[i] = nsets;
                        }

                      for (size_t i = 0; i < nwpv * gridsize; i++)
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
                  else if (lminidx)
                    {
                      farminidx(pvars1, pvars2, field, nsets);
                    }
                  else if (lmaxidx)
                    {
                      farmaxidx(pvars1, pvars2, field, nsets);
                    }
                  else
                    {
                      farfun(pvars1, field, operfunc);
                    }
                }
            }

          if (nsets == 0 && lvarstd)
            for (int recID = 0; recID < maxrecs; recID++)
              {
                if (recinfo[recID].lconst) continue;

                int varID = recinfo[recID].varID;
                int levelID = recinfo[recID].levelID;
                Field *pvars1 = &vars1[varID][levelID];
                Field *pvars2 = &vars2[varID][levelID];

                farmoq(pvars2, *pvars1);
              }

          vdate0 = vdate;
          vtime0 = vtime;
          nsets++;
          tsID++;
        }

      if (nrecs == 0 && nsets == 0) break;

      for (int recID = 0; recID < maxrecs; recID++)
        {
          if (recinfo[recID].lconst) continue;

          int varID = recinfo[recID].varID;
          int levelID = recinfo[recID].levelID;
          Field *psamp1 = &samp1[varID][levelID];
          Field *pvars1 = &vars1[varID][levelID];
          Field *pvars2 = vars2 ? &vars2[varID][levelID] : NULL;

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

      if (cdoVerbose)
        {
          char vdatestr[32], vtimestr[32];
          date2str(vdate0, vdatestr, sizeof(vdatestr));
          time2str(vtime0, vtimestr, sizeof(vtimestr));
          cdoPrint("%s %s  vfrac = %g, nsets = %d", vdatestr, vtimestr, vfrac, nsets);
        }

      if (lvfrac && operfunc == func_mean)
        for (int recID = 0; recID < maxrecs; recID++)
          {
            if (recinfo[recID].lconst) continue;

            int varID = recinfo[recID].varID;
            int levelID = recinfo[recID].levelID;
            Field *psamp1 = &samp1[varID][levelID];
            Field *pvars1 = &vars1[varID][levelID];

            int nwpv = pvars1->nwpv;
            size_t gridsize = gridInqSize(pvars1->grid);
            double missval = pvars1->missval;
            if (psamp1->ptr)
              {
                int irun = 0;
                for (size_t i = 0; i < nwpv * gridsize; ++i)
                  {
                    if ((psamp1->ptr[i] / nsets) < vfrac)
                      {
                        pvars1->ptr[i] = missval;
                        irun++;
                      }
                  }

                if (irun) pvars1->nmiss = arrayNumMV(nwpv * gridsize, pvars1->ptr, missval);
              }
          }

      dtlist_stat_taxisDefTimestep(dtlist, taxisID2, nsets);
      pstreamDefTimestep(streamID2, otsID);

      if (cdoDiag)
        {
          dtlist_stat_taxisDefTimestep(dtlist, taxisID3, nsets);
          pstreamDefTimestep(streamID3, otsID);
        }

      for (int recID = 0; recID < maxrecs; recID++)
        {
          if (otsID && recinfo[recID].lconst) continue;

          int varID = recinfo[recID].varID;
          int levelID = recinfo[recID].levelID;
          Field *psamp1 = &samp1[varID][levelID];
          Field *pvars1 = &vars1[varID][levelID];

          pstreamDefRecord(streamID2, varID, levelID);
          pstreamWriteRecord(streamID2, pvars1->ptr, pvars1->nmiss);

          if (cdoDiag)
            {
              double *sampptr = field.ptr;
              if (psamp1->ptr)
                sampptr = psamp1->ptr;
              else
                {
                  size_t gridsize = pvars1->size;
                  for (size_t i = 0; i < gridsize; ++i) sampptr[i] = nsets;
                }

              pstreamDefRecord(streamID3, varID, levelID);
              pstreamWriteRecord(streamID3, sampptr, 0);
            }
        }

      if (nrecs == 0) break;
      otsID++;
    }

  field_free(vars1, vlistID1);
  field_free(samp1, vlistID1);
  if (lvarstd || lrange || lminidx || lmaxidx) field_free(vars2, vlistID1);

  dtlist_delete(dtlist);

  if (cdoDiag) pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if (field.ptr) Free(field.ptr);

  cdoFinish();

  return 0;
}

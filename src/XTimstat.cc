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

      Timstat    timmin          Time minimum
      Timstat    timmax          Time maximum
      Timstat    timsum          Time sum
      Timstat    timmean         Time mean
      Timstat    timavg          Time average
      Timstat    timvar          Time variance
      Timstat    timvar1         Time variance [Normalize by (n-1)]
      Timstat    timstd          Time standard deviation
      Timstat    timstd1         Time standard deviation [Normalize by (n-1)]
      Hourstat   hourmin         Hourly minimum
      Hourstat   hourmax         Hourly maximum
      Hourstat   hoursum         Hourly sum
      Hourstat   hourmean        Hourly mean
      Hourstat   houravg         Hourly average
      Hourstat   hourvar         Hourly variance
      Hourstat   hourvar1        Hourly variance [Normalize by (n-1)]
      Hourstat   hourstd         Hourly standard deviation
      Hourstat   hourstd1        Hourly standard deviation [Normalize by (n-1)]
      Daystat    daymin          Daily minimum
      Daystat    daymax          Daily maximum
      Daystat    daysum          Daily sum
      Daystat    daymean         Daily mean
      Daystat    dayavg          Daily average
      Daystat    dayvar          Daily variance
      Daystat    dayvar1         Daily variance [Normalize by (n-1)]
      Daystat    daystd          Daily standard deviation
      Daystat    daystd1         Daily standard deviation [Normalize by (n-1)]
      Monstat    monmin          Monthly minimum
      Monstat    monmax          Monthly maximum
      Monstat    monsum          Monthly sum
      Monstat    monmean         Monthly mean
      Monstat    monavg          Monthly average
      Monstat    monvar          Monthly variance
      Monstat    monvar1         Monthly variance [Normalize by (n-1)]
      Monstat    monstd          Monthly standard deviation
      Monstat    monstd1         Monthly standard deviation [Normalize by (n-1)]
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
#include "cdo_task.h"
#include "pstream_int.h"
#include "timer.h"
#include "datetime.h"
//#include "pstream_write.h"

struct readarg_t
{
  int tsIDnext;
  int streamID, nrecs;
  RecordInfo *recinfo;
  Field **vars;
};

static int num_recs = 0;

static void *
cdoReadTimestep(void *rarg)
{
  int varID, levelID;
  size_t nmiss;
  readarg_t *readarg = (readarg_t *) rarg;
  Field **input_vars = readarg->vars;
  RecordInfo *recinfo = readarg->recinfo;
  int streamID = readarg->streamID;
  int tsIDnext = readarg->tsIDnext;
  int nrecs = readarg->nrecs;

  // timer_start(timer_read);

  for (int recID = 0; recID < nrecs; ++recID)
    {
      pstreamInqRecord(streamID, &varID, &levelID);

      if (tsIDnext == 1 && recinfo)
        {
          recinfo[recID].varID = varID;
          recinfo[recID].levelID = levelID;
        }

      if (CDO_Memtype == MEMTYPE_FLOAT)
        pstreamReadRecordF(streamID, input_vars[varID][levelID].ptr2f, &nmiss);
      else
        pstreamReadRecord(streamID, input_vars[varID][levelID].ptr2, &nmiss);

      input_vars[varID][levelID].nmiss2 = nmiss;
    }

  // timer_stop(timer_read);

  num_recs = cdoStreamInqTimestep(streamID, tsIDnext);

  return ((void *) &num_recs);
}

static void
cdoUpdateVars(int nvars, int vlistID, Field **vars)
{
  for (int varID = 0; varID < nvars; varID++)
    {
      int nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID, varID));
      for (int levelID = 0; levelID < nlevels; levelID++)
        {
          if (CDO_Memtype == MEMTYPE_FLOAT)
            {
              float *tmp = vars[varID][levelID].ptrf;
              vars[varID][levelID].ptrf = vars[varID][levelID].ptr2f;
              vars[varID][levelID].ptr2f = tmp;
            }
          else
            {
              double *tmp = vars[varID][levelID].ptr;
              vars[varID][levelID].ptr = vars[varID][levelID].ptr2;
              vars[varID][levelID].ptr2 = tmp;
            }
          vars[varID][levelID].nmiss = vars[varID][levelID].nmiss2;
        }
    }
}

void *
XTimstat(void *process)
{
  enum
  {
    HOUR_LEN = 4,
    DAY_LEN = 6,
    MON_LEN = 8,
    YEAR_LEN = 10
  };
  TimeStat timestat_date = TimeStat::MEAN;
  int64_t vdate = 0, vdate0 = 0;
  int vtime = 0, vtime0 = 0;
  int varID;
  int streamID3 = -1;
  int vlistID3, taxisID3 = -1;
  bool lvfrac = false;
  int nwpv;  // number of words per value; real:1  complex:2
  char indate1[DATE_LEN + 1], indate2[DATE_LEN + 1];
  double vfrac = 1;

  cdoInitialize(process);

  // clang-format off
  cdoOperatorAdd("xtimmin",    func_min,   DATE_LEN, NULL);
  cdoOperatorAdd("xtimmax",    func_max,   DATE_LEN, NULL);
  cdoOperatorAdd("xtimsum",    func_sum,   DATE_LEN, NULL);
  cdoOperatorAdd("xtimmean",   func_mean,  DATE_LEN, NULL);
  cdoOperatorAdd("xtimavg",    func_avg,   DATE_LEN, NULL);
  cdoOperatorAdd("xtimvar",    func_var,   DATE_LEN, NULL);
  cdoOperatorAdd("xtimvar1",   func_var1,  DATE_LEN, NULL);
  cdoOperatorAdd("xtimstd",    func_std,   DATE_LEN, NULL);
  cdoOperatorAdd("xtimstd1",   func_std1,  DATE_LEN, NULL);
  cdoOperatorAdd("xyearmin",   func_min,   YEAR_LEN, NULL);
  cdoOperatorAdd("xyearmax",   func_max,   YEAR_LEN, NULL);
  cdoOperatorAdd("xyearsum",   func_sum,   YEAR_LEN, NULL);
  cdoOperatorAdd("xyearmean",  func_mean,  YEAR_LEN, NULL);
  cdoOperatorAdd("xyearavg",   func_avg,   YEAR_LEN, NULL);
  cdoOperatorAdd("xyearvar",   func_var,   YEAR_LEN, NULL);
  cdoOperatorAdd("xyearvar1",  func_var1,  YEAR_LEN, NULL);
  cdoOperatorAdd("xyearstd",   func_std,   YEAR_LEN, NULL);
  cdoOperatorAdd("xyearstd1",  func_std1,  YEAR_LEN, NULL);
  cdoOperatorAdd("xmonmin",    func_min,   MON_LEN, NULL);
  cdoOperatorAdd("xmonmax",    func_max,   MON_LEN, NULL);
  cdoOperatorAdd("xmonsum",    func_sum,   MON_LEN, NULL);
  cdoOperatorAdd("xmonmean",   func_mean,  MON_LEN, NULL);
  cdoOperatorAdd("xmonavg",    func_avg,   MON_LEN, NULL);
  cdoOperatorAdd("xmonvar",    func_var,   MON_LEN, NULL);
  cdoOperatorAdd("xmonvar1",   func_var1,  MON_LEN, NULL);
  cdoOperatorAdd("xmonstd",    func_std,   MON_LEN, NULL);
  cdoOperatorAdd("xmonstd1",   func_std1,  MON_LEN, NULL);

  int operatorID = cdoOperatorID();
  int operfunc   = cdoOperatorF1(operatorID);
  int comparelen = cdoOperatorF2(operatorID);

  bool lmean   = operfunc == func_mean || operfunc == func_avg;
  bool lstd    = operfunc == func_std || operfunc == func_std1;
  bool lvarstd = operfunc == func_std || operfunc == func_var || operfunc == func_std1 || operfunc == func_var1;
  int divisor  = operfunc == func_std1 || operfunc == func_var1;
  // clang-format on

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
  // int streamID1 = cdoStreamOpenRead(cdoGetStreamName(0).c_str());

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

  DateTimeList dtlist;
  dtlist.setStat(timestat_date);
  dtlist.setCalendar(taxisInqCalendar(taxisID1));

  size_t gridsizemax = vlistGridsizeMax(vlistID1);
  if (vlistNumber(vlistID1) != CDI_REAL) gridsizemax *= 2;

  int FIELD_MEMTYPE = 0;
  if (CDO_Memtype == MEMTYPE_FLOAT) FIELD_MEMTYPE = FIELD_FLT;
  Field **input_vars = field_malloc(vlistID1, FIELD_PTR | FIELD_PTR2 | FIELD_MEMTYPE);
  Field **vars1 = field_malloc(vlistID1, FIELD_PTR);
  Field **samp1 = field_malloc(vlistID1, FIELD_NONE);
  Field **vars2 = NULL;
  if (lvarstd) vars2 = field_malloc(vlistID1, FIELD_PTR);

  readarg_t readarg;
  readarg.streamID = streamID1;
  readarg.vars = input_vars;

  bool lparallelread = CDO_Parallel_Read > 0;
  bool ltsfirst = true;
  void *read_task = NULL;
  void *readresult = NULL;

  if (lparallelread)
    {
      read_task = cdo_task_new();
      if (read_task == NULL)
        {
          lparallelread = false;
          cdoWarning("CDO tasks not available!");
        }
    }

  int tsID = 0;
  int otsID = 0;
  int nrecs = cdoStreamInqTimestep(streamID1, tsID);
  int maxrecs = nrecs;
  RecordInfo *recinfo = (RecordInfo *) Malloc(maxrecs * sizeof(RecordInfo));

  tsID++;
  while (TRUE)
    {
      int nsets = 0;
      while (nrecs > 0)
        {
          dtlist.taxisInqTimestep(taxisID1, nsets);
          vdate = dtlist.getVdate(nsets);
          vtime = dtlist.getVtime(nsets);

          if (nsets == 0) SET_DATE(indate2, vdate, vtime);
          SET_DATE(indate1, vdate, vtime);

          if (DATE_IS_NEQ(indate1, indate2, cmplen)) break;

          readarg.tsIDnext = tsID;
          readarg.nrecs = nrecs;
          readarg.recinfo = recinfo;

          if (ltsfirst || lparallelread == false)
            {
              ltsfirst = false;
              readresult = cdoReadTimestep(&readarg);
            }
          else
            {
              readresult = cdo_task_wait(read_task);
            }

          nrecs = *(int *) readresult;

          cdoUpdateVars(nvars, vlistID1, input_vars);

          if (nrecs && lparallelread)
            {
              readarg.tsIDnext = tsID + 1;
              cdo_task_start(read_task, cdoReadTimestep, &readarg);
            }

          if (nsets == 0)
            {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(maxrecs, recinfo, input_vars, vars1, samp1) if (maxrecs > 1)
#endif
              for (int recID = 0; recID < maxrecs; recID++)
                {
                  int varID = recinfo[recID].varID;
                  int levelID = recinfo[recID].levelID;

                  Field *pvars1 = &vars1[varID][levelID];
                  Field *pinput_var = &input_vars[varID][levelID];

                  int nwpv = pvars1->nwpv;
                  size_t gridsize = pvars1->size;
                  size_t nmiss = pinput_var->nmiss;

                  farcpy(pvars1, *pinput_var);
                  pvars1->nmiss = nmiss;
                  if (nmiss > 0 || samp1[varID][levelID].ptr)
                    {
                      if (samp1[varID][levelID].ptr == NULL)
                        samp1[varID][levelID].ptr = (double *) malloc(nwpv * gridsize * sizeof(double));

                      for (size_t i = 0; i < nwpv * gridsize; i++)
                        samp1[varID][levelID].ptr[i] = !DBL_IS_EQUAL(pvars1->ptr[i], pvars1->missval);
                    }
                }
            }
          else
            {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(lvarstd, nsets, maxrecs, recinfo, input_vars, vars1, samp1, vars2, \
                                              operfunc) if (maxrecs > 1)
#endif
              for (int recID = 0; recID < maxrecs; recID++)
                {
                  int varID = recinfo[recID].varID;
                  int levelID = recinfo[recID].levelID;

                  Field *pvars1 = &vars1[varID][levelID];
                  Field *pinput_var = &input_vars[varID][levelID];

                  int nwpv = pvars1->nwpv;
                  size_t gridsize = pvars1->size;
                  size_t nmiss = pinput_var->nmiss;

                  if (nmiss > 0 || samp1[varID][levelID].ptr)
                    {
                      if (samp1[varID][levelID].ptr == NULL)
                        {
                          samp1[varID][levelID].ptr = (double *) malloc(nwpv * gridsize * sizeof(double));
                          for (size_t i = 0; i < nwpv * gridsize; i++) samp1[varID][levelID].ptr[i] = nsets;
                        }

                      for (size_t i = 0; i < nwpv * gridsize; i++)
                        if (!DBL_IS_EQUAL(pinput_var->ptr[i], pvars1->missval)) samp1[varID][levelID].ptr[i]++;
                    }

                  if (lvarstd)
                    {
                      Field *pvars2 = &vars2[varID][levelID];
                      farsumq(pvars2, *pinput_var);
                      farsum(pvars1, *pinput_var);
                    }
                  else
                    {
                      farfun(pvars1, *pinput_var, operfunc);
                    }
                }
            }

          if (nsets == 0 && lvarstd)
            for (int recID = 0; recID < maxrecs; recID++)
              {
                int varID = recinfo[recID].varID;
                int levelID = recinfo[recID].levelID;
                Field *pvars1 = &vars1[varID][levelID];
                Field *pvars2 = &vars2[varID][levelID];

                if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT) continue;

                farmoq(pvars2, *pvars1);
              }

          vdate0 = vdate;
          vtime0 = vtime;
          nsets++;
          tsID++;
        }

      if (nrecs == 0 && nsets == 0) break;

      if (lmean)
        {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(vlistID1, nsets, maxrecs, recinfo, vars1, samp1) if (maxrecs > 1)
#endif
          for (int recID = 0; recID < maxrecs; recID++)
            {
              int varID = recinfo[recID].varID;
              int levelID = recinfo[recID].levelID;
              Field *pvars1 = &vars1[varID][levelID];

              if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT) continue;

              if (samp1[varID][levelID].ptr == NULL)
                farcdiv(pvars1, (double) nsets);
              else
                fardiv(pvars1, samp1[varID][levelID]);
            }
        }
      else if (lvarstd)
        {
          for (int recID = 0; recID < maxrecs; recID++)
            {
              int varID = recinfo[recID].varID;
              int levelID = recinfo[recID].levelID;
              Field *pvars1 = &vars1[varID][levelID];
              Field *pvars2 = &vars2[varID][levelID];

              if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT) continue;

              if (samp1[varID][levelID].ptr == NULL)
                {
                  if (lstd)
                    farcstd(pvars1, *pvars2, nsets, divisor);
                  else
                    farcvar(pvars1, *pvars2, nsets, divisor);
                }
              else
                {
                  if (lstd)
                    farstd(pvars1, *pvars2, samp1[varID][levelID], divisor);
                  else
                    farvar(pvars1, *pvars2, samp1[varID][levelID], divisor);
                }
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
            int varID = recinfo[recID].varID;
            int levelID = recinfo[recID].levelID;
            Field *pvars1 = &vars1[varID][levelID];

            if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT) continue;

            nwpv = pvars1->nwpv;
            size_t gridsize = pvars1->size;
            double missval = pvars1->missval;
            if (samp1[varID][levelID].ptr)
              {
                size_t irun = 0;
                for (size_t i = 0; i < nwpv * gridsize; ++i)
                  {
                    if ((samp1[varID][levelID].ptr[i] / nsets) < vfrac)
                      {
                        pvars1->ptr[i] = missval;
                        irun++;
                      }
                  }

                if (irun) pvars1->nmiss = arrayNumMV(nwpv * gridsize, pvars1->ptr, missval);
              }
          }

      dtlist.statTaxisDefTimestep(taxisID2, nsets);
      pstreamDefTimestep(streamID2, otsID);

      if (cdoDiag)
        {
          dtlist.statTaxisDefTimestep(taxisID3, nsets);
          pstreamDefTimestep(streamID3, otsID);
        }

      for (int recID = 0; recID < maxrecs; recID++)
        {
          int varID = recinfo[recID].varID;
          int levelID = recinfo[recID].levelID;
          Field *pvars1 = &vars1[varID][levelID];

          if (otsID && vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT) continue;

          pstreamDefRecord(streamID2, varID, levelID);
          pstreamWriteRecord(streamID2, pvars1->ptr, pvars1->nmiss);

          if (cdoDiag)
            {
              if (samp1[varID][levelID].ptr)
                {
                  pstreamDefRecord(streamID3, varID, levelID);
                  pstreamWriteRecord(streamID3, samp1[varID][levelID].ptr, 0);
                }
            }
        }

      if (nrecs == 0) break;
      otsID++;
    }

  field_free(input_vars, vlistID1);
  field_free(vars1, vlistID1);
  field_free(samp1, vlistID1);
  if (lvarstd) field_free(vars2, vlistID1);

  Free(recinfo);

  if (cdoDiag) pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

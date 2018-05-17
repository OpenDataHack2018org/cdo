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

      Runstat    runrange        Running range
      Runstat    runmin          Running minimum
      Runstat    runmax          Running maximum
      Runstat    runsum          Running sum
      Runstat    runmean         Running mean
      Runstat    runavg          Running average
      Runstat    runvar          Running variance
      Runstat    runvar1         Running variance [Normalize by (n-1)]
      Runstat    runstd          Running standard deviation
      Runstat    runstd1         Running standard deviation [Normalize by (n-1)]
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "datetime.h"

void *
Runstat(void *process)
{
  TimeStat timestat_date = TimeStat::MEAN;
  int varID;
  int levelID;
  size_t nmiss;
  bool runstat_nomiss = false;

  cdoInitialize(process);

  char *envstr = getenv("RUNSTAT_NOMISS");
  if (envstr)
    {
      char *endptr;
      int envval = (int) strtol(envstr, &endptr, 10);
      if (envval == 1) runstat_nomiss = true;
    }

  // clang-format off
  cdoOperatorAdd("runrange", func_range, 0, NULL);
  cdoOperatorAdd("runmin",   func_min,   0, NULL);
  cdoOperatorAdd("runmax",   func_max,   0, NULL);
  cdoOperatorAdd("runsum",   func_sum,   0, NULL);
  cdoOperatorAdd("runmean",  func_mean,  0, NULL);
  cdoOperatorAdd("runavg",   func_avg,   0, NULL);
  cdoOperatorAdd("runvar",   func_var,   0, NULL);
  cdoOperatorAdd("runvar1",  func_var1,  0, NULL);
  cdoOperatorAdd("runstd",   func_std,   0, NULL);
  cdoOperatorAdd("runstd1",  func_std1,  0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  bool lrange  = operfunc == func_range;
  bool lmean   = operfunc == func_mean || operfunc == func_avg;
  bool lstd    = operfunc == func_std || operfunc == func_std1;
  bool lvarstd = operfunc == func_std || operfunc == func_var || operfunc == func_std1 || operfunc == func_var1;
  int  divisor = operfunc == func_std1 || operfunc == func_var1;
  // clang-format on

  operatorInputArg("number of timesteps");
  operatorCheckArgc(1);
  int ndates = parameter2int(operatorArgv()[0]);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  taxisWithBounds(taxisID2);
  vlistDefTaxis(vlistID2, taxisID2);
  /*  Number of timestep will be reduced compared to the input
   *  error handling in case of not enough timesteps is done per record */
  int nsteps = vlistNtsteps(vlistID1);
  if (nsteps != -1)
    {
      nsteps -= ndates - 1;
      if (nsteps > 0) vlistDefNtsteps(vlistID2, nsteps);
    }

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recinfo(maxrecs);

  DateTimeList dtlist;
  dtlist.setStat(timestat_date);
  dtlist.setCalendar(taxisInqCalendar(taxisID1));

  Field ***vars1 = (Field ***) Malloc((ndates + 1) * sizeof(Field **));
  Field ***vars2 = NULL, ***samp1 = NULL;
  if (!runstat_nomiss) samp1 = (Field ***) Malloc((ndates + 1) * sizeof(Field **));
  if (lvarstd || lrange) vars2 = (Field ***) Malloc((ndates + 1) * sizeof(Field **));

  for (int its = 0; its < ndates; its++)
    {
      vars1[its] = field_malloc(vlistID1, FIELD_PTR);
      if (!runstat_nomiss) samp1[its] = field_malloc(vlistID1, FIELD_PTR);
      if (lvarstd || lrange) vars2[its] = field_malloc(vlistID1, FIELD_PTR);
    }

  size_t gridsizemax = vlistGridsizeMax(vlistID1);
  bool *imask = (bool *) Malloc(gridsizemax * sizeof(bool));

  int tsID = 0;
  for (tsID = 0; tsID < ndates; tsID++)
    {
      int nrecs = cdoStreamInqTimestep(streamID1, tsID);
      if (nrecs == 0) cdoAbort("File has less then %d timesteps!", ndates);

      dtlist.taxisInqTimestep(taxisID1, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);

          if (tsID == 0)
            {
              recinfo[recID].varID = varID;
              recinfo[recID].levelID = levelID;
              recinfo[recID].lconst = vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT;
            }

          Field *psamp1 = samp1 ? &samp1[tsID][varID][levelID] : NULL;
          Field *pvars1 = &vars1[tsID][varID][levelID];
          Field *pvars2 = vars2 ? &vars2[tsID][varID][levelID] : NULL;

          size_t gridsize = pvars1->size;

          pstreamReadRecord(streamID1, pvars1->ptr, &nmiss);
          pvars1->nmiss = nmiss;
          if (lrange)
            {
              pvars2->nmiss = pvars1->nmiss;
              for (size_t i = 0; i < gridsize; i++) pvars2->ptr[i] = pvars1->ptr[i];
            }

          if (runstat_nomiss && nmiss > 0)
            cdoAbort("Missing values supported was swichted off by env. "
                     "RUNSTAT_NOMISS!");

          if (!runstat_nomiss)
            {
              double missval = pvars1->missval;

              for (size_t i = 0; i < gridsize; i++) imask[i] = !DBL_IS_EQUAL(pvars1->ptr[i], missval);

              for (size_t i = 0; i < gridsize; i++) psamp1->ptr[i] = (double) imask[i];

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(tsID, gridsize, imask, samp1, varID, levelID)
#endif
              for (int inp = 0; inp < tsID; inp++)
                {
                  double *ptr = samp1[inp][varID][levelID].ptr;
                  for (size_t i = 0; i < gridsize; i++)
                    if (imask[i]) ptr[i]++;
                }
            }

          if (lvarstd)
            {
              farmoq(pvars2, *pvars1);
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(tsID, vars1, vars2, varID, levelID, pvars1)
#endif
              for (int inp = 0; inp < tsID; inp++)
                {
                  farsumq(&vars2[inp][varID][levelID], *pvars1);
                  farsum(&vars1[inp][varID][levelID], *pvars1);
                }
            }
          else if (lrange)
            {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(tsID, vars1, vars2, varID, levelID, pvars1)
#endif
              for (int inp = 0; inp < tsID; inp++)
                {
                  farmin(&vars2[inp][varID][levelID], *pvars1);
                  farmax(&vars1[inp][varID][levelID], *pvars1);
                }
            }
          else
            {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(tsID, vars1, operfunc, varID, levelID, pvars1)
#endif
              for (int inp = 0; inp < tsID; inp++)
                {
                  farfun(&vars1[inp][varID][levelID], *pvars1, operfunc);
                }
            }
        }
    }

  int otsID = 0;
  while (TRUE)
    {
      for (int recID = 0; recID < maxrecs; recID++)
        {
          if (recinfo[recID].lconst) continue;

          int varID = recinfo[recID].varID;
          int levelID = recinfo[recID].levelID;
          Field *psamp1 = samp1 ? &samp1[0][varID][levelID] : NULL;
          Field *pvars1 = &vars1[0][varID][levelID];
          Field *pvars2 = vars2 ? &vars2[0][varID][levelID] : NULL;
          int nsets = ndates;

          if (lmean)
            {
              if (!runstat_nomiss)
                fardiv(pvars1, *psamp1);
              else
                farcdiv(pvars1, (double) nsets);
            }
          else if (lvarstd)
            {
              if (!runstat_nomiss)
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

      dtlist.statTaxisDefTimestep(taxisID2, ndates);
      pstreamDefTimestep(streamID2, otsID);

      for (int recID = 0; recID < maxrecs; recID++)
        {
          if (otsID && recinfo[recID].lconst) continue;

          int varID = recinfo[recID].varID;
          int levelID = recinfo[recID].levelID;
          Field *pvars1 = &vars1[0][varID][levelID];

          pstreamDefRecord(streamID2, varID, levelID);
          pstreamWriteRecord(streamID2, pvars1->ptr, pvars1->nmiss);
        }

      otsID++;

      dtlist.shift();

      vars1[ndates] = vars1[0];
      if (!runstat_nomiss) samp1[ndates] = samp1[0];
      if (lvarstd || lrange) vars2[ndates] = vars2[0];

      for (int inp = 0; inp < ndates; inp++)
        {
          vars1[inp] = vars1[inp + 1];
          if (!runstat_nomiss) samp1[inp] = samp1[inp + 1];
          if (lvarstd || lrange) vars2[inp] = vars2[inp + 1];
        }

      int nrecs = cdoStreamInqTimestep(streamID1, tsID);
      if (nrecs == 0) break;

      dtlist.taxisInqTimestep(taxisID1, ndates - 1);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);

          Field *psamp1 = samp1 ? &samp1[ndates - 1][varID][levelID] : NULL;
          Field *pvars1 = &vars1[ndates - 1][varID][levelID];
          Field *pvars2 = vars2 ? &vars2[ndates - 1][varID][levelID] : NULL;

          size_t gridsize = pvars1->size;

          pstreamReadRecord(streamID1, pvars1->ptr, &nmiss);
          pvars1->nmiss = nmiss;
          if (lrange)
            {
              pvars2->nmiss = pvars1->nmiss;
              for (size_t i = 0; i < gridsize; i++) pvars2->ptr[i] = pvars1->ptr[i];
            }

          if (runstat_nomiss && nmiss > 0) cdoAbort("Missing values supported swichted off!");

          if (!runstat_nomiss)
            {
              double missval = pvars1->missval;

              for (size_t i = 0; i < gridsize; i++) imask[i] = !DBL_IS_EQUAL(pvars1->ptr[i], missval);

              for (size_t i = 0; i < gridsize; i++) psamp1->ptr[i] = (double) imask[i];

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(ndates, imask, gridsize, samp1, varID, levelID)
#endif
              for (int inp = 0; inp < ndates - 1; inp++)
                {
                  double *ptr = samp1[inp][varID][levelID].ptr;
                  for (size_t i = 0; i < gridsize; i++)
                    if (imask[i]) ptr[i]++;
                }
            }

          if (lvarstd)
            {
              farmoq(pvars2, *pvars1);
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(ndates, vars1, vars2, varID, levelID, pvars1)
#endif
              for (int inp = 0; inp < ndates - 1; inp++)
                {
                  farsumq(&vars2[inp][varID][levelID], *pvars1);
                  farsum(&vars1[inp][varID][levelID], *pvars1);
                }
            }
          else if (lrange)
            {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(ndates, vars1, vars2, varID, levelID, pvars1)
#endif
              for (int inp = 0; inp < ndates - 1; inp++)
                {
                  farmin(&vars2[inp][varID][levelID], *pvars1);
                  farmax(&vars1[inp][varID][levelID], *pvars1);
                }
            }
          else
            {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(ndates, vars1, varID, levelID, operfunc, pvars1)
#endif
              for (int inp = 0; inp < ndates - 1; inp++)
                {
                  farfun(&vars1[inp][varID][levelID], *pvars1, operfunc);
                }
            }
        }

      tsID++;
    }

  for (int its = 0; its < ndates; its++)
    {
      field_free(vars1[its], vlistID1);
      if (!runstat_nomiss) field_free(samp1[its], vlistID1);
      if (lvarstd || lrange) field_free(vars2[its], vlistID1);
    }

  Free(vars1);
  if (!runstat_nomiss) Free(samp1);
  if (lvarstd || lrange) Free(vars2);

  if (imask) Free(imask);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

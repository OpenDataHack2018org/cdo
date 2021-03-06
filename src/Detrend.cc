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

      Detrend    detrend         Detrend
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "cdoOptions.h"
#include "datetime.h"

#define NALLOC_INC 1024

static void
detrend(long nts, double missval1, double *array1, double *array2)
{
  double zj;
  double sumj, sumjj;
  double sumx, sumjx;
  double missval2 = missval1;

  sumx = sumjx = 0;
  sumj = sumjj = 0;
  long n = 0;
  for (long j = 0; j < nts; j++)
    if (!DBL_IS_EQUAL(array1[j], missval1))
      {
        zj = j;
        sumx += array1[j];
        sumjx += zj * array1[j];
        sumj += zj;
        sumjj += zj * zj;
        n++;
      }

  double work1 = DIVMN(SUBMN(sumjx, DIVMN(MULMN(sumx, sumj), n)), SUBMN(sumjj, DIVMN(MULMN(sumj, sumj), n)));
  double work2 = SUBMN(DIVMN(sumx, n), MULMN(work1, DIVMN(sumj, n)));

  for (long j = 0; j < nts; j++) array2[j] = SUBMN(array1[j], ADDMN(work2, MULMN(j, work1)));
}

void *
Detrend(void *process)
{
  size_t gridsize;
  int nrecs;
  int gridID, varID, levelID;
  size_t i;
  int nalloc = 0;
  size_t nmiss;
  int nlevel;
  double missval;
  DateTimeList dtlist;

  cdoInitialize(process);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int nvars = vlistNvars(vlistID1);
  std::vector<Field **> vars;

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      if (tsID >= nalloc)
        {
          nalloc += NALLOC_INC;
          vars.resize(nalloc);
        }

      dtlist.taxisInqTimestep(taxisID1, tsID);

      vars[tsID] = field_malloc(vlistID1, FIELD_NONE);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          gridID = vlistInqVarGrid(vlistID1, varID);
          gridsize = gridInqSize(gridID);
          vars[tsID][varID][levelID].ptr = (double *) Malloc(gridsize * sizeof(double));
          pstreamReadRecord(streamID1, vars[tsID][varID][levelID].ptr, &nmiss);
          vars[tsID][varID][levelID].nmiss = nmiss;
        }

      tsID++;
    }

  int nts = tsID;

  VECTOR_2D(double, array1, Threading::ompNumThreads, nts);
  VECTOR_2D(double, array2, Threading::ompNumThreads, nts);

  for (varID = 0; varID < nvars; varID++)
    {
      gridID = vlistInqVarGrid(vlistID1, varID);
      missval = vlistInqVarMissval(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      for (levelID = 0; levelID < nlevel; levelID++)
        {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(array1, array2, vars, varID, levelID, gridsize, nts, missval)
#endif
          for (i = 0; i < gridsize; i++)
            {
              int ompthID = cdo_omp_get_thread_num();

              for (int tsID = 0; tsID < nts; tsID++) array1[ompthID][tsID] = vars[tsID][varID][levelID].ptr[i];

              detrend(nts, missval, &array1[ompthID][0], &array2[ompthID][0]);

              for (int tsID = 0; tsID < nts; tsID++) vars[tsID][varID][levelID].ptr[i] = array2[ompthID][tsID];
            }
        }
    }

  for (int tsID = 0; tsID < nts; tsID++)
    {
      dtlist.taxisDefTimestep(taxisID2, tsID);
      pstreamDefTimestep(streamID2, tsID);

      for (varID = 0; varID < nvars; varID++)
        {
          nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
          for (levelID = 0; levelID < nlevel; levelID++)
            {
              if (vars[tsID][varID][levelID].ptr)
                {
                  nmiss = vars[tsID][varID][levelID].nmiss;
                  pstreamDefRecord(streamID2, varID, levelID);
                  pstreamWriteRecord(streamID2, vars[tsID][varID][levelID].ptr, nmiss);
                  Free(vars[tsID][varID][levelID].ptr);
                  vars[tsID][varID][levelID].ptr = NULL;
                }
            }
        }

      field_free(vars[tsID], vlistID1);
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

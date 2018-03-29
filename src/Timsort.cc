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

     Timsort    timsort         Sort over the time
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "cdoOptions.h"

#define NALLOC_INC 1024

static int
compareDouble(const void *a, const void *b)
{
  const double *x = (const double *) a;
  const double *y = (const double *) b;
  return ((*x > *y) - (*x < *y)) * 2 + (*x > *y) - (*x < *y);
}

void *
Timsort(void *process)
{
  int nrecs;
  int gridID, varID, levelID;
  int nalloc = 0;
  size_t nmiss;

  cdoInitialize(process);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisCreate(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int nvars = vlistNvars(vlistID1);
  std::vector<field_type **> vars;
  std::vector<int> vdate, vtime;

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      if (tsID >= nalloc)
        {
          nalloc += NALLOC_INC;
          vdate.resize(nalloc);
          vtime.resize(nalloc);
          vars.resize(nalloc);
        }

      vdate[tsID] = taxisInqVdate(taxisID1);
      vtime[tsID] = taxisInqVtime(taxisID1);

      vars[tsID] = field_malloc(vlistID1, FIELD_NONE);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          gridID = vlistInqVarGrid(vlistID1, varID);
          size_t gridsize = gridInqSize(gridID);
          vars[tsID][varID][levelID].ptr = (double *) Malloc(gridsize * sizeof(double));
          pstreamReadRecord(streamID1, vars[tsID][varID][levelID].ptr, &nmiss);
          vars[tsID][varID][levelID].nmiss = nmiss;
        }

      tsID++;
    }

  int nts = tsID;

  double **sarray = (double **) Malloc(Threading::ompNumThreads * sizeof(double *));
  for (int i = 0; i < Threading::ompNumThreads; i++) sarray[i] = (double *) Malloc(nts * sizeof(double));

  for (varID = 0; varID < nvars; varID++)
    {
      if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT) continue;

      gridID = vlistInqVarGrid(vlistID1, varID);
      size_t gridsize = gridInqSize(gridID);
      int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      for (levelID = 0; levelID < nlevel; levelID++)
        {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(gridsize, nts, sarray, vars, varID, levelID)
#endif
          for (size_t i = 0; i < gridsize; i++)
            {
              int ompthID = cdo_omp_get_thread_num();

              for (int tsID = 0; tsID < nts; tsID++) sarray[ompthID][tsID] = vars[tsID][varID][levelID].ptr[i];

              qsort(sarray[ompthID], nts, sizeof(double), compareDouble);

              for (int tsID = 0; tsID < nts; tsID++) vars[tsID][varID][levelID].ptr[i] = sarray[ompthID][tsID];
            }
        }
    }

  for (int i = 0; i < Threading::ompNumThreads; i++)
    if (sarray[i]) Free(sarray[i]);

  if (sarray) Free(sarray);

  for (tsID = 0; tsID < nts; tsID++)
    {
      taxisDefVdate(taxisID2, vdate[tsID]);
      taxisDefVtime(taxisID2, vtime[tsID]);
      pstreamDefTimestep(streamID2, tsID);

      for (varID = 0; varID < nvars; varID++)
        {
          int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
          for (levelID = 0; levelID < nlevel; levelID++)
            {
              if (vars[tsID][varID][levelID].ptr)
                {
                  nmiss = vars[tsID][varID][levelID].nmiss;
                  pstreamDefRecord(streamID2, varID, levelID);
                  pstreamWriteRecord(streamID2, vars[tsID][varID][levelID].ptr, nmiss);
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

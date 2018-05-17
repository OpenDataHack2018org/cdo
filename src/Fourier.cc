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

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "statistic.h"
#include "cdoOptions.h"

#define NALLOC_INC 1024

void *
Fourier(void *process)
{
  int nrecs;
  int varID, levelID;
  int nalloc = 0;
  size_t nmiss;
  struct FourierMemory
  {
    std::vector<double> real;
    std::vector<double> imag;
    std::vector<double> work_r;
    std::vector<double> work_i;
  };

  cdoInitialize(process);

  operatorInputArg("the sign of the exponent (-1 for normal or 1 for reverse transformation)!");
  int sign = parameter2int(operatorArgv()[0]);

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
  std::vector<int64_t> vdate;
  std::vector<int> vtime;

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
          int gridID = vlistInqVarGrid(vlistID1, varID);
          size_t gridsize = gridInqSize(gridID);
          vars[tsID][varID][levelID].ptr = (double *) Malloc(2 * gridsize * sizeof(double));
          pstreamReadRecord(streamID1, vars[tsID][varID][levelID].ptr, &nmiss);
          vars[tsID][varID][levelID].nmiss = nmiss;
        }

      tsID++;
    }

  int nts = tsID;

  bool lpower2 = ((nts & (nts - 1)) == 0);

  std::vector<FourierMemory> ompmem(Threading::ompNumThreads);
  for (int i = 0; i < Threading::ompNumThreads; i++)
    {
      ompmem[i].real.resize(nts);
      ompmem[i].imag.resize(nts);
      if (!lpower2) ompmem[i].work_r.resize(nts);
      if (!lpower2) ompmem[i].work_i.resize(nts);
    }

  for (varID = 0; varID < nvars; varID++)
    {
      int gridID = vlistInqVarGrid(vlistID1, varID);
      double missval = vlistInqVarMissval(vlistID1, varID);
      size_t gridsize = gridInqSize(gridID);
      int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      for (levelID = 0; levelID < nlevel; levelID++)
        {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(tsID)
#endif
          for (size_t i = 0; i < gridsize; i++)
            {
              int lmiss = 0;
              int ompthID = cdo_omp_get_thread_num();

              for (tsID = 0; tsID < nts; tsID++)
                {
                  ompmem[ompthID].real[tsID] = vars[tsID][varID][levelID].ptr[2 * i];
                  ompmem[ompthID].imag[tsID] = vars[tsID][varID][levelID].ptr[2 * i + 1];
                  if (DBL_IS_EQUAL(ompmem[ompthID].real[tsID], missval) || DBL_IS_EQUAL(ompmem[ompthID].imag[tsID], missval))
                    lmiss = 1;
                }

              if (lmiss == 0)
                {
                  if (lpower2) /* nts is a power of 2 */
                    fft(ompmem[ompthID].real.data(), ompmem[ompthID].imag.data(), nts, sign);
                  else
                    ft_r(ompmem[ompthID].real.data(), ompmem[ompthID].imag.data(), nts, sign, ompmem[ompthID].work_r.data(),
                         ompmem[ompthID].work_i.data());

                  for (tsID = 0; tsID < nts; tsID++)
                    {
                      vars[tsID][varID][levelID].ptr[2 * i] = ompmem[ompthID].real[tsID];
                      vars[tsID][varID][levelID].ptr[2 * i + 1] = ompmem[ompthID].imag[tsID];
                    }
                }
              else
                {
                  for (tsID = 0; tsID < nts; tsID++)
                    {
                      vars[tsID][varID][levelID].ptr[2 * i] = missval;
                      vars[tsID][varID][levelID].ptr[2 * i + 1] = missval;
                    }
                }
            }
        }
    }

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

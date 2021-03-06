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

      Spectrum    spectrum         Spectrum
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

#include "statistic.h"

#define NALLOC_INC 1024

static void
spectrum(int nrec, double *data, double *spectrum, double *real, double *imag, double *window, double wssum, int detrend, int seg_n,
         int seg_l)
{
  int k;
  double sumx, sumkx;
  double a, b;
  double *work_r = NULL, *work_i = NULL;
  int seg_i, offset;
  int bit;

  for (bit = seg_l; !(bit & 1); bit >>= 1)
    ;

  if (detrend == 1)
    {
      sumx = 0;
      for (k = 0; k < nrec; k++) sumx += data[k];
      sumx /= nrec;
      for (k = 0; k < nrec; k++) data[k] -= sumx;
    }
  else if (detrend == 2)
    {
      sumx = sumkx = 0;
      for (k = 0; k < nrec; k++)
        {
          sumx += data[k];
          sumkx += k * data[k];
        }
      b = (sumkx - sumx * (nrec - 1) / 2.) / ((nrec + 1) * nrec * (nrec - 1) / 12.);
      a = sumx / nrec - b * (nrec - 1) / 2.;
      for (k = 0; k < nrec; k++) data[k] -= a + b * k;
    }

  if (bit != 1)
    {
      work_r = (double *) Malloc(seg_l * sizeof(double));
      work_i = (double *) Malloc(seg_l * sizeof(double));
    }

  for (seg_i = 0; seg_i < seg_n; seg_i += 2)
    {
      offset = seg_n == 1 ? 0 : (int) ((double) (nrec - seg_l) / (seg_n - 1) * seg_i);

      for (k = 0; k < seg_l; k++) real[k] = data[offset + k];

      if (detrend == 3)
        {
          sumx = sumkx = 0;
          for (k = 0; k < seg_l; k++)
            {
              sumx += real[k];
              sumkx += k * real[k];
            }

          b = (sumkx - sumx * (seg_l - 1) / 2.) / ((seg_l + 1) * seg_l * (seg_l - 1) / 12.);
          a = sumx / seg_l - b * (seg_l - 1) / 2.;

          for (k = 0; k < seg_l; k++) real[k] -= a + b * k;
        }

      for (k = 0; k < seg_l; k++) real[k] *= window[k];

      if (seg_i + 1 < seg_n)
        {
          offset = seg_n == 1 ? 0 : (int) ((double) (nrec - seg_l) / (seg_n - 1) * (seg_i + 1));
          for (k = 0; k < seg_l; k++) imag[k] = data[offset + k];
          if (detrend == 3)
            {
              sumx = sumkx = 0;
              for (k = 0; k < seg_l; k++)
                {
                  sumx += imag[k];
                  sumkx += k * imag[k];
                }

              b = (sumkx - sumx * (seg_l - 1) / 2.) / ((seg_l + 1) * seg_l * (seg_l - 1) / 12.);
              a = sumx / seg_l - b * (seg_l - 1) / 2.;

              for (k = 0; k < seg_l; k++) imag[k] -= a + b * k;
            }

          for (k = 0; k < seg_l; k++) imag[k] *= window[k];
        }
      else
        for (k = 0; k < seg_l; k++) imag[k] = 0;

      if (bit == 1) /* seg_l is a power of 2 */
        fft(real, imag, seg_l, 1);
      else
        ft_r(real, imag, seg_l, 1, work_r, work_i);

      spectrum[0] += real[0] * real[0] + imag[0] * imag[0];

      for (k = 1; k < (seg_l + 1) / 2; k++)
        spectrum[k]
            += real[k] * real[k] + imag[k] * imag[k] + real[seg_l - k] * real[seg_l - k] + imag[seg_l - k] * imag[seg_l - k];

      if (!(seg_l & 1)) spectrum[seg_l / 2] += real[seg_l / 2] * real[seg_l / 2] + imag[seg_l / 2] * imag[seg_l / 2];
    }

  if (bit != 1)
    {
      Free(work_r);
      Free(work_i);
    }

  for (k = 0; k <= seg_l / 2; k++) spectrum[k] *= seg_l / (seg_n * wssum);

  spectrum[0] *= 2;

  if (!(seg_l & 1)) spectrum[seg_l / 2] *= 2;
}

void *
Spectrum(void *process)
{
  int nrecs;
  int gridID, varID, levelID;
  int k;
  int nalloc = 0;
  size_t nmiss;
  int freq;

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
          gridID = vlistInqVarGrid(vlistID1, varID);
          size_t gridsize = gridInqSize(gridID);
          vars[tsID][varID][levelID].ptr = (double *) Malloc(gridsize * sizeof(double));
          pstreamReadRecord(streamID1, vars[tsID][varID][levelID].ptr, &nmiss);
          vars[tsID][varID][levelID].nmiss = nmiss;

          if (nmiss > 0) cdoAbort("Missing values are not allowed!");
        }

      tsID++;
    }

  int nts = tsID;

  operatorInputArg("detrend type, length of segments, number of segments, window type\n\n"
                   "       detrend type: 0 - data should be used unchanged\n"
                   "                     1 - the mean of the whole time series should be subtracted\n"
                   "                     2 - the whole time series should be detrended\n"
                   "                     3 - every segment should be detrended\n\n"
                   "        window type: 0 - no data windowing\n"
                   "                     1 - Hann window\n"
                   "                     2 - Bartlett window\n"
                   "                     3 - Welch window\n");

  operatorCheckArgc(4);

  int detrend = parameter2int(operatorArgv()[0]);
  int seg_l = parameter2int(operatorArgv()[1]);
  int seg_n = parameter2int(operatorArgv()[2]);
  int which_window = parameter2int(operatorArgv()[3]);

  if (detrend < 0 || detrend > 3) cdoAbort("Illegal value for detrend (=%d)!", detrend);

  if (seg_l <= 2 || seg_l > nts) cdoAbort("Length must be at least 3 and at most the number of timesteps (=%d)", nts);

  if (seg_n <= 0 || seg_n > nts - seg_l + 1)
    cdoAbort("Number of segments must be positiv and not greater than %d!", nts - seg_l + 1);

  int nfreq = seg_l / 2 + 1;

  std::vector<Field **> vars2(nfreq);
  for (freq = 0; freq < nfreq; freq++) vars2[freq] = field_malloc(vlistID1, FIELD_PTR);

  double *array1 = (double *) Malloc(nts * sizeof(double));
  double *array2 = (double *) Malloc(nfreq * sizeof(double));
  double *real = (double *) Malloc(seg_l * sizeof(double));
  double *imag = (double *) Malloc(seg_l * sizeof(double));
  double *window = (double *) Malloc(seg_l * sizeof(double));

  switch (which_window)
    {
    case 0:
      for (k = 0; k < seg_l; k++) window[k] = 1;
      break;
    case 1:
      for (k = 0; k < seg_l; k++) window[k] = 1 - cos(2 * M_PI * (k + 1) / (seg_l + 1));
      break;
    case 2:
      for (k = 0; k < seg_l / 2; k++) window[k] = window[seg_l - 1 - k] = k;
      break;
    case 3:
      for (k = 0; k < seg_l; k++)
        {
          double temp;
          temp = ((k + 1.) - 0.5 * (seg_l + 1.)) / (0.5 * (seg_l + 1.));
          window[k] = 1 - temp * temp;
        }
      break;
    default: cdoAbort("Invalid window type %d!", which_window); break;
    }

  double wssum = 0;
  for (k = 0; k < seg_l; k++) wssum += window[k] * window[k];

  for (varID = 0; varID < nvars; varID++)
    {
      gridID = vlistInqVarGrid(vlistID1, varID);
      size_t gridsize = gridInqSize(gridID);
      int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      for (levelID = 0; levelID < nlevel; levelID++)
        {
          for (size_t i = 0; i < gridsize; i++)
            {
              for (tsID = 0; tsID < nts; tsID++) array1[tsID] = vars[tsID][varID][levelID].ptr[i];

              for (freq = 0; freq < nfreq; freq++) array2[freq] = 0;

              spectrum(nts, array1, array2, real, imag, window, wssum, detrend, seg_n, seg_l);

              for (freq = 0; freq < nfreq; freq++) vars2[freq][varID][levelID].ptr[i] = array2[freq];
            }
        }
    }

  if (array1) Free(array1);
  if (array2) Free(array2);

  for (tsID = 0; tsID < nts; tsID++) field_free(vars[tsID], vlistID1);

  for (tsID = 0; tsID < nfreq; tsID++)
    {
      taxisDefVdate(taxisID2, vdate[0]);
      taxisDefVtime(taxisID2, vtime[0]);
      pstreamDefTimestep(streamID2, tsID);

      for (varID = 0; varID < nvars; varID++)
        {
          int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
          for (levelID = 0; levelID < nlevel; levelID++)
            {
              if (vars2[tsID][varID][levelID].ptr)
                {
                  nmiss = vars2[tsID][varID][levelID].nmiss;
                  pstreamDefRecord(streamID2, varID, levelID);
                  pstreamWriteRecord(streamID2, vars2[tsID][varID][levelID].ptr, 0);
                }
            }
        }

      field_free(vars2[tsID], vlistID1);
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

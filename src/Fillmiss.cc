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

*/
#include <time.h>  // clock()

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "grid.h"
#include "grid_point_search.h"
#include "cdoOptions.h"

void
fillmiss(Field *field1, Field *field2, int nfill)
{
  int nx, ny, i, j;
  size_t nmiss2 = 0;
  int kr, ku, kl, ko;
  int ir, iu, il, io;
  int kh, kv, k1, k2, kk;
  int globgrid = FALSE, gridtype;
  double s1, s2;
  double xr, xu, xl, xo;
  double **matrix1, **matrix2;

  int gridID = field1->grid;
  size_t nmiss1 = field1->nmiss;
  double missval = field1->missval;
  double *array1 = field1->ptr;
  double *array2 = field2->ptr;

  nx = gridInqXsize(gridID);
  ny = gridInqYsize(gridID);
  globgrid = gridIsCircular(gridID);

  gridtype = gridInqType(gridID);
  if (!(gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN)) cdoAbort("Unsupported grid type: %s!", gridNamePtr(gridtype));

  matrix1 = (double **) Malloc(ny * sizeof(double *));
  matrix2 = (double **) Malloc(ny * sizeof(double *));

  for (j = 0; j < ny; j++)
    {
      matrix1[j] = array1 + j * nx;
      matrix2[j] = array2 + j * nx;
    }

  for (j = 0; j < ny; j++)
    for (i = 0; i < nx; i++)
      {
        if (DBL_IS_EQUAL(matrix1[j][i], missval))
          {
            nmiss2++;

            kr = ku = kl = ko = 0;
            xr = xu = xl = xo = 0.;

            for (ir = i + 1; ir < nx; ir++)
              if (!DBL_IS_EQUAL(matrix1[j][ir], missval))
                {
                  kr = ir - i;
                  xr = matrix1[j][ir];
                  break;
                }

            if (globgrid && ir == nx)
              {
                for (ir = 0; ir < i; ir++)
                  if (!DBL_IS_EQUAL(matrix1[j][ir], missval))
                    {
                      kr = nx + ir - i;
                      xr = matrix1[j][ir];
                      break;
                    }
              }

            for (il = i - 1; il >= 0; il--)
              if (!DBL_IS_EQUAL(matrix1[j][il], missval))
                {
                  kl = i - il;
                  xl = matrix1[j][il];
                  break;
                }

            if (globgrid && il == -1)
              {
                for (il = nx - 1; il > i; il--)
                  if (!DBL_IS_EQUAL(matrix1[j][il], missval))
                    {
                      kl = nx + i - il;
                      xl = matrix1[j][il];
                      break;
                    }
              }

            for (iu = j + 1; iu < ny; iu++)
              if (!DBL_IS_EQUAL(matrix1[iu][i], missval))
                {
                  ku = iu - j;
                  xu = matrix1[iu][i];
                  break;
                }

            for (io = j - 1; io >= 0; io--)
              if (!DBL_IS_EQUAL(matrix1[io][i], missval))
                {
                  ko = j - io;
                  xo = matrix1[io][i];
                  break;
                }

            /*  printf("%d %d %d %d %d %d %g %g %g %g\n",
             * j,i,kr,kl,ku,ko,xr,xl,xu,xo);*/

            kh = kl + kr;
            kv = ko + ku;
            if (kh == 0)
              {
                s1 = 0.;
                k1 = 0;
              }
            else if (kl == 0)
              {
                s1 = xr;
                k1 = 1;
              }
            else if (kr == 0)
              {
                s1 = xl;
                k1 = 1;
              }
            else
              {
                s1 = xr * kl / kh + xl * kr / kh;
                k1 = 2;
              }

            if (kv == 0)
              {
                s2 = 0.;
                k2 = 0;
              }
            else if (ku == 0)
              {
                s2 = xo;
                k2 = 1;
              }
            else if (ko == 0)
              {
                s2 = xu;
                k2 = 1;
              }
            else
              {
                s2 = xu * ko / kv + xo * ku / kv;
                k2 = 2;
              }

            kk = k1 + k2;
            if (kk >= nfill)
              {
                if (kk == 0)
                  cdoAbort("no point found!");
                else if (k1 == 0)
                  matrix2[j][i] = s2;
                else if (k2 == 0)
                  matrix2[j][i] = s1;
                else
                  matrix2[j][i] = s1 * k2 / kk + s2 * k1 / kk;
              }
            else
              matrix2[j][i] = matrix1[j][i];

            /* matrix1[j][i] = matrix2[j][i]; */
          }
        else
          {
            matrix2[j][i] = matrix1[j][i];
          }
      }

  if (nmiss1 != nmiss2) cdoAbort("found only %zu of %zu missing values!", nmiss2, nmiss1);

  Free(matrix2);
  Free(matrix1);
}

void
fillmiss_one_step(Field *field1, Field *field2, int maxfill)
{
  int gridID, nx, ny, i, j;
  size_t nmiss2 = 0;
  int kr, ku, kl, ko;
  int ir, iu, il, io;
  int kh, kv, k1, k2, kk;
  double s1, s2;
  double xr, xu, xl, xo;
  double missval;
  double *array1, *array2;
  double **matrix1, **matrix2;

  gridID = field1->grid;
  missval = field1->missval;
  array1 = field1->ptr;
  array2 = field2->ptr;

  nx = gridInqXsize(gridID);
  ny = gridInqYsize(gridID);

  matrix1 = (double **) Malloc(ny * sizeof(double *));
  matrix2 = (double **) Malloc(ny * sizeof(double *));

  for (j = 0; j < ny; j++)
    {
      matrix1[j] = array1 + j * nx;
      matrix2[j] = array2 + j * nx;
    }

  for (int fill_iterations = 0; fill_iterations < maxfill; fill_iterations++)
    {
      for (j = 0; j < ny; j++)
        for (i = 0; i < nx; i++)
          {
            if (DBL_IS_EQUAL(matrix1[j][i], missval))
              {
                nmiss2++;

                kr = ku = kl = ko = 0;
                xr = xu = xl = xo = 0.;

                for (ir = i + 1; ir < nx; ir++)
                  if (!DBL_IS_EQUAL(matrix1[j][ir], missval))
                    {
                      kr = ir - i;
                      xr = matrix1[j][ir];
                      break;
                    }

                for (il = i - 1; il >= 0; il--)
                  if (!DBL_IS_EQUAL(matrix1[j][il], missval))
                    {
                      kl = i - il;
                      xl = matrix1[j][il];
                      break;
                    }

                for (iu = j + 1; iu < ny; iu++)
                  if (!DBL_IS_EQUAL(matrix1[iu][i], missval))
                    {
                      ku = iu - j;
                      xu = matrix1[iu][i];
                      break;
                    }

                for (io = j - 1; io >= 0; io--)
                  if (!DBL_IS_EQUAL(matrix1[io][i], missval))
                    {
                      ko = j - io;
                      xo = matrix1[io][i];
                      break;
                    }

                kh = kl + kr;
                kv = ko + ku;
                if (kh == 0)
                  {
                    s1 = 0.;
                    k1 = 0;
                  }
                else if (kl == 0)
                  {
                    s1 = xr;
                    k1 = kr;
                  }
                else if (kr == 0)
                  {
                    s1 = xl;
                    k1 = kl;
                  }
                else
                  {
                    if (kl < kr)
                      {
                        s1 = xl;
                        k1 = kl;
                      }
                    else
                      {
                        s1 = xr;
                        k1 = kr;
                      }
                  }

                if (kv == 0)
                  {
                    s2 = 0.;
                    k2 = 0;
                  }
                else if (ku == 0)
                  {
                    s2 = xo;
                    k2 = ko;
                  }
                else if (ko == 0)
                  {
                    s2 = xu;
                    k2 = ku;
                  }
                else
                  {
                    if (ku < ko)
                      {
                        s2 = xu;
                        k2 = ku;
                      }
                    else
                      {
                        s2 = xo;
                        k2 = ko;
                      }
                  }

                kk = k1 + k2;
                if (kk == 0)
                  matrix2[j][i] = matrix1[j][i];
                else if (k1 == 0)
                  matrix2[j][i] = s2;
                else if (k2 == 0)
                  matrix2[j][i] = s1;
                else
                  {
                    if (k1 <= k2)
                      {
                        matrix2[j][i] = s1;
                      }
                    else
                      {
                        matrix2[j][i] = s2;
                      }
                  }

                // printf("%d %d %2d %2d %2d %2d %2g %2g %2g %2g %2g %2g %2g\n",
                // j,i,kr,kl,ku,ko,xr,xl,xu,xo,s1,s2,matrix2[j][i]);
                /* matrix1[j][i] = matrix2[j][i]; */
              }
            else
              {
                matrix2[j][i] = matrix1[j][i];
              }
          }
      for (j = 0; j < ny; j++)
        for (i = 0; i < nx; i++) matrix1[j][i] = matrix2[j][i];
    }

  Free(matrix2);
  Free(matrix1);
}

static void
setmisstodis(Field *field1, Field *field2, int numNeighbors)
{
  int gridID = field1->grid;
  int gridID0 = gridID;
  double missval = field1->missval;
  double *array1 = field1->ptr;
  double *array2 = field2->ptr;

  size_t gridsize = gridInqSize(gridID);

  size_t nmiss = field1->nmiss;
  size_t nvals = gridsize - nmiss;

  if (gridInqType(gridID) == GRID_GME) gridID = gridToUnstructured(gridID, 0);

  if (gridInqType(gridID) != GRID_UNSTRUCTURED && gridInqType(gridID) != GRID_CURVILINEAR) gridID = gridToCurvilinear(gridID, 0);

  std::vector<double> xvals(gridsize);
  std::vector<double> yvals(gridsize);

  gridInqXvals(gridID, &xvals[0]);
  gridInqYvals(gridID, &yvals[0]);

  /* Convert lat/lon units if required */
  char units[CDI_MAX_NAME];
  gridInqXunits(gridID, units);
  grid_to_radian(units, gridsize, &xvals[0], "grid center lon");
  gridInqYunits(gridID, units);
  grid_to_radian(units, gridsize, &yvals[0], "grid center lat");

  std::vector<size_t> mindex(nmiss, 1);
  std::vector<size_t> vindex(nvals, 1);
  std::vector<double> lons(nvals);
  std::vector<double> lats(nvals);

  size_t nv = 0, nm = 0;
  for (size_t i = 0; i < gridsize; ++i)
    {
      array2[i] = array1[i];
      if (DBL_IS_EQUAL(array1[i], missval))
        {
          mindex[nm] = i;
          nm++;
        }
      else
        {
          lons[nv] = xvals[i];
          lats[nv] = yvals[i];
          vindex[nv] = i;
          nv++;
        }
    }

  if (nv != nvals) cdoAbort("Internal problem, number of valid values differ!");

  std::vector<knnWeightsType> knnWeights;
  for (int i = 0; i < Threading::ompNumThreads; ++i) knnWeights.push_back(knnWeightsType(numNeighbors));

  clock_t start, finish;
  start = clock();

  GridPointSearch *gps = NULL;

  if (nmiss)
    {
      bool xIsCyclic = false;
      size_t dims[2] = { nvals, 0 };
      gps = gridPointSearchCreate(xIsCyclic, dims, nvals, &lons[0], &lats[0]);
      gridPointSearchExtrapolate(gps);
    }

  finish = clock();

  if (cdoVerbose) cdoPrint("Point search created: %.2f seconds", ((double) (finish - start)) / CLOCKS_PER_SEC);

  progressInit();

  start = clock();

  double findex = 0;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(knnWeights, findex, mindex, vindex, array1, array2, xvals, yvals, gps, nmiss, \
                                              numNeighbors)
#endif
  for (size_t i = 0; i < nmiss; ++i)
    {
#ifdef _OPENMP
#pragma omp atomic
#endif
      findex++;
      if (cdo_omp_get_thread_num() == 0) progressStatus(0, 1, findex / nmiss);

      int ompthID = cdo_omp_get_thread_num();

      gridSearchPoint(gps, xvals[mindex[i]], yvals[mindex[i]], knnWeights[ompthID]);

      // Compute weights based on inverse distance if mask is false, eliminate those points
      size_t nadds = knnWeights[ompthID].compute_weights();
      if (nadds)
        {
          double result = 0;
          for (size_t n = 0; n < nadds; ++n)
            result += array1[vindex[knnWeights[ompthID].m_addr[n]]] * knnWeights[ompthID].m_dist[n];
          array2[mindex[i]] = result;
        }
    }

  progressStatus(0, 1, 1);

  finish = clock();

  if (cdoVerbose) cdoPrint("Point search nearest: %.2f seconds", ((double) (finish - start)) / CLOCKS_PER_SEC);

  if (gps) gridPointSearchDelete(gps);

  if (gridID0 != gridID) gridDestroy(gridID);
}

void *
Fillmiss(void *process)
{
  size_t nmiss;
  int nrecs, varID, levelID;
  void (*fill_method)(Field * fin, Field * fout, int) = NULL;

  cdoInitialize(process);

  // clang-format off
  int FILLMISS        = cdoOperatorAdd("fillmiss"   ,   0, 0, "nfill");
  int FILLMISSONESTEP = cdoOperatorAdd("fillmiss2"  ,   0, 0, "nfill");
  int SETMISSTONN     = cdoOperatorAdd("setmisstonn" ,  0, 0, "");
  int SETMISSTODIS    = cdoOperatorAdd("setmisstodis" , 0, 0, "number of neighbors");
  // clang-format on

  int operatorID = cdoOperatorID();

  int nfill = 1;
  if (operatorID == FILLMISS)
    {
      fill_method = &fillmiss;
    }
  else if (operatorID == FILLMISSONESTEP)
    {
      fill_method = &fillmiss_one_step;
    }
  else if (operatorID == SETMISSTONN)
    {
      fill_method = &setmisstodis;
    }
  else if (operatorID == SETMISSTODIS)
    {
      nfill = 4;
      fill_method = &setmisstodis;
    }

  /* Argument handling */
  {
    int oargc = operatorArgc();
    char **oargv = operatorArgv();

    if (oargc == 1)
      {
        nfill = parameter2int(oargv[0]);
        if (operatorID == FILLMISS)
          {
            if (nfill < 1 || nfill > 4) cdoAbort("nfill out of range!");
          }
      }
    else if (oargc > 1)
      cdoAbort("Too many arguments!");
  }

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());

  pstreamDefVlist(streamID2, vlistID2);

  size_t gridsize = vlistGridsizeMax(vlistID1);

  Field field1, field2;
  field_init(&field1);
  field_init(&field2);
  field1.ptr = (double *) Malloc(gridsize * sizeof(double));
  field2.ptr = (double *) Malloc(gridsize * sizeof(double));

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      pstreamDefTimestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamReadRecord(streamID1, field1.ptr, &nmiss);
          field1.nmiss = nmiss;

          pstreamDefRecord(streamID2, varID, levelID);

          if (field1.nmiss == 0)
            {
              pstreamWriteRecord(streamID2, field1.ptr, 0);
            }
          else
            {
              int gridID = vlistInqVarGrid(vlistID1, varID);

              if ((operatorID == FILLMISS || operatorID == FILLMISSONESTEP)
                  && (gridInqType(gridID) == GRID_GME || gridInqType(gridID) == GRID_UNSTRUCTURED))
                cdoAbort("%s data unsupported!", gridNamePtr(gridInqType(gridID)));

              field1.grid = gridID;
              field1.missval = vlistInqVarMissval(vlistID1, varID);

              field2.grid = field1.grid;
              field2.nmiss = 0;
              field2.missval = field1.missval;

              fill_method(&field1, &field2, nfill);

              size_t gridsize = gridInqSize(field2.grid);
              size_t nmiss = arrayNumMV(gridsize, field2.ptr, field2.missval);
              pstreamWriteRecord(streamID2, field2.ptr, nmiss);
            }
        }
      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if (field2.ptr) Free(field2.ptr);
  if (field1.ptr) Free(field1.ptr);

  cdoFinish();

  return 0;
}

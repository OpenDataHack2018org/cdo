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

#include "cdo_int.h"
#include <cdi.h>
#include "percentiles.h"

void
zonfun(Field field1, Field *field2, int function)
{
  // clang-format off
  switch (function)
    {
    case func_min:   zonmin(field1, field2);    break;
    case func_max:   zonmax(field1, field2);    break;
    case func_range: zonrange(field1, field2);  break;
    case func_sum:   zonsum(field1, field2);    break;
    case func_mean:  zonmean(field1, field2);   break;
    case func_avg:   zonavg(field1, field2);    break;
    case func_std:   zonstd(field1, field2);    break;
    case func_std1:  zonstd1(field1, field2);   break;
    case func_var:   zonvar(field1, field2);    break;
    case func_var1:  zonvar1(field1, field2);   break;
    default: cdoAbort("function %d not implemented!", function);
    }
  // clang-format on
}

void
zonmin(Field field1, Field *field2)
{
  size_t rnmiss = 0;
  double rmin;

  size_t nx = gridInqXsize(field1.grid);
  size_t ny = gridInqYsize(field1.grid);

  for (size_t j = 0; j < ny; ++j)
    {
      if (field1.nmiss)
        {
          rmin = arrayMinMV(nx, &field1.ptr[j * nx], field1.missval);
          if (DBL_IS_EQUAL(rmin, field1.missval)) rnmiss++;
        }
      else
        {
          rmin = arrayMin(nx, &field1.ptr[j * nx]);
        }

      field2->ptr[j] = rmin;
    }

  field2->nmiss = rnmiss;
}

void
zonmax(Field field1, Field *field2)
{
  size_t rnmiss = 0;
  double rmax;

  size_t nx = gridInqXsize(field1.grid);
  size_t ny = gridInqYsize(field1.grid);

  for (size_t j = 0; j < ny; ++j)
    {
      if (field1.nmiss)
        {
          rmax = arrayMaxMV(nx, &field1.ptr[j * nx], field1.missval);
          if (DBL_IS_EQUAL(rmax, field1.missval)) rnmiss++;
        }
      else
        {
          rmax = arrayMax(nx, &field1.ptr[j * nx]);
        }

      field2->ptr[j] = rmax;
    }

  field2->nmiss = rnmiss;
}

void
zonrange(Field field1, Field *field2)
{
  size_t rnmiss = 0;
  double range;

  size_t nx = gridInqXsize(field1.grid);
  size_t ny = gridInqYsize(field1.grid);

  for (size_t j = 0; j < ny; ++j)
    {
      if (field1.nmiss)
        {
          range = arrayRangeMV(nx, &field1.ptr[j * nx], field1.missval);
          if (DBL_IS_EQUAL(range, field1.missval)) rnmiss++;
        }
      else
        {
          range = arrayRange(nx, &field1.ptr[j * nx]);
        }

      field2->ptr[j] = range;
    }

  field2->nmiss = rnmiss;
}

void
zonsum(Field field1, Field *field2)
{
  size_t rnmiss = 0;
  double rsum = 0;

  size_t nx = gridInqXsize(field1.grid);
  size_t ny = gridInqYsize(field1.grid);

  for (size_t j = 0; j < ny; ++j)
    {
      if (field1.nmiss)
        {
          rsum = arraySumMV(nx, &field1.ptr[j * nx], field1.missval);
          if (DBL_IS_EQUAL(rsum, field1.missval)) rnmiss++;
        }
      else
        {
          rsum = arraySum(nx, &field1.ptr[j * nx]);
        }

      field2->ptr[j] = rsum;
    }

  field2->nmiss = rnmiss;
}

void
zonmean(Field field1, Field *field2)
{
  size_t rnmiss = 0;
  double rmean = 0;

  size_t nx = gridInqXsize(field1.grid);
  size_t ny = gridInqYsize(field1.grid);

  for (size_t j = 0; j < ny; ++j)
    {
      if (field1.nmiss)
        {
          rmean = arrayMeanMV(nx, &field1.ptr[j * nx], field1.missval);
        }
      else
        {
          rmean = arrayMean(nx, &field1.ptr[j * nx]);
        }

      if (DBL_IS_EQUAL(rmean, field1.missval)) rnmiss++;

      field2->ptr[j] = rmean;
    }

  field2->nmiss = rnmiss;
}

void
zonavg(Field field1, Field *field2)
{
  size_t rnmiss = 0;
  double ravg = 0;

  size_t nx = gridInqXsize(field1.grid);
  size_t ny = gridInqYsize(field1.grid);

  for (size_t j = 0; j < ny; ++j)
    {
      if (field1.nmiss > 0)
        {
          ravg = arrayAvgMV(nx, &field1.ptr[j * nx], field1.missval);
        }
      else
        {
          ravg = arrayMean(nx, &field1.ptr[j * nx]);
        }

      if (DBL_IS_EQUAL(ravg, field1.missval)) rnmiss++;

      field2->ptr[j] = ravg;
    }

  field2->nmiss = rnmiss;
}

static void
prevarsum_zon(const double *restrict array, size_t nx, size_t nmiss, double missval, double *rsum, double *rsumw, double *rsumq,
              double *rsumwq)
{
  double w = 1. / nx;

  *rsum = 0;
  *rsumq = 0;
  *rsumw = 0;
  *rsumwq = 0;

  if (nmiss > 0)
    {
      for (size_t i = 0; i < nx; i++)
        if (!DBL_IS_EQUAL(array[i], missval))
          {
            *rsum += w * array[i];
            *rsumq += w * array[i] * array[i];
            *rsumw += w;
            *rsumwq += w * w;
          }
    }
  else
    {
      for (size_t i = 0; i < nx; i++)
        {
          *rsum += w * array[i];
          *rsumq += w * array[i] * array[i];
          *rsumw += w;
          *rsumwq += w * w;
        }
    }
}

void
zonvar(Field field1, Field *field2)
{
  size_t rnmiss = 0;
  int grid = field1.grid;
  size_t nmiss = field1.nmiss;
  double missval1 = field1.missval;
  double *array = field1.ptr;
  double rsum = 0, rsumw = 0, rvar = 0;
  double rsumq = 0, rsumwq = 0;

  size_t nx = gridInqXsize(grid);
  size_t ny = gridInqYsize(grid);

  for (size_t j = 0; j < ny; j++)
    {
      prevarsum_zon(array + j * nx, nx, nmiss, missval1, &rsum, &rsumw, &rsumq, &rsumwq);

      rvar = IS_NOT_EQUAL(rsumw, 0) ? (rsumq * rsumw - rsum * rsum) / (rsumw * rsumw) : missval1;
      if (rvar < 0 && rvar > -1.e-5) rvar = 0;

      if (DBL_IS_EQUAL(rvar, missval1)) rnmiss++;

      field2->ptr[j] = rvar;
    }

  field2->nmiss = rnmiss;
}

void
zonvar1(Field field1, Field *field2)
{
  size_t rnmiss = 0;
  int grid = field1.grid;
  size_t nmiss = field1.nmiss;
  double missval1 = field1.missval;
  double *array = field1.ptr;
  double rsum = 0, rsumw = 0, rvar = 0;
  double rsumq = 0, rsumwq = 0;

  size_t nx = gridInqXsize(grid);
  size_t ny = gridInqYsize(grid);

  for (size_t j = 0; j < ny; j++)
    {
      prevarsum_zon(array + j * nx, nx, nmiss, missval1, &rsum, &rsumw, &rsumq, &rsumwq);

      rvar = (rsumw * rsumw > rsumwq) ? (rsumq * rsumw - rsum * rsum) / (rsumw * rsumw - rsumwq) : missval1;
      if (rvar < 0 && rvar > -1.e-5) rvar = 0;

      if (DBL_IS_EQUAL(rvar, missval1)) rnmiss++;

      field2->ptr[j] = rvar;
    }

  field2->nmiss = rnmiss;
}

void
zonstd(Field field1, Field *field2)
{
  size_t rnmiss = 0;
  int grid = field1.grid;
  double missval = field1.missval;
  double rstd;

  size_t ny = gridInqYsize(grid);

  zonvar(field1, field2);

  for (size_t j = 0; j < ny; j++)
    {
      rstd = varToStd(field2->ptr[j], missval);

      if (DBL_IS_EQUAL(rstd, missval)) rnmiss++;

      field2->ptr[j] = rstd;
    }

  field2->nmiss = rnmiss;
}

void
zonstd1(Field field1, Field *field2)
{
  size_t rnmiss = 0;
  int grid = field1.grid;
  double missval = field1.missval;
  double rstd;

  size_t ny = gridInqYsize(grid);

  zonvar1(field1, field2);

  for (size_t j = 0; j < ny; j++)
    {
      rstd = varToStd(field2->ptr[j], missval);

      if (DBL_IS_EQUAL(rstd, missval)) rnmiss++;

      field2->ptr[j] = rstd;
    }

  field2->nmiss = rnmiss;
}

/* RQ */
void
zonpctl(Field field1, Field *field2, int p)
{
  size_t rnmiss = 0;
  int grid = field1.grid;
  size_t nmiss = field1.nmiss;
  double missval = field1.missval;
  double *array = field1.ptr;

  size_t nx = gridInqXsize(grid);
  size_t ny = gridInqYsize(grid);

  if (nmiss > 0)
    {
      double *array2 = (double *) Malloc(nx * sizeof(double));

      for (size_t j = 0; j < ny; j++)
        {
          size_t l = 0;
          for (size_t i = 0; i < nx; i++)
            if (!DBL_IS_EQUAL(array[j * nx + i], missval)) array2[l++] = array[j * nx + i];

          if (l > 0)
            {
              field2->ptr[j] = percentile(array2, l, p);
            }
          else
            {
              field2->ptr[j] = missval;
              rnmiss++;
            }
        }

      Free(array2);
    }
  else
    {
      for (size_t j = 0; j < ny; j++)
        {
          if (nx > 0)
            {
              field2->ptr[j] = percentile(&array[j * nx], nx, p);
            }
          else
            {
              field2->ptr[j] = missval;
              rnmiss++;
            }
        }
    }

  field2->nmiss = rnmiss;
}
/* QR */

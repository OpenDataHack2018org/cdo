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
#include <stdio.h>
#include <float.h>
#include <fenv.h>
#include <assert.h>

#include "compare.h"
#include "array.h"

//#pragma STDC FENV_ACCESS ON

const char *
fpe_errstr(int fpeRaised)
{
  const char *errstr = NULL;

  if (fpeRaised & FE_DIVBYZERO)
    errstr = "division by zero";
  else if (fpeRaised & FE_INEXACT)
    errstr = "inexact result";
  else if (fpeRaised & FE_INVALID)
    errstr = "invalid result";
  else if (fpeRaised & FE_OVERFLOW)
    errstr = "overflow";
  else if (fpeRaised & FE_UNDERFLOW)
    errstr = "underflow";

  return errstr;
}

void
arrayMinMax(size_t len, const double *array, double *rmin, double *rmax)
{
  double min = DBL_MAX;
  double max = -DBL_MAX;

  // #pragma omp parallel for default(none) shared(min, max, array, gridsize)
  // reduction(+:mean) #pragma omp simd reduction(+:mean) reduction(min:min)
  // reduction(max:max) aligned(array:16)
  for (size_t i = 0; i < len; ++i)
    {
      if (array[i] < min) min = array[i];
      if (array[i] > max) max = array[i];
    }

  if (rmin) *rmin = min;
  if (rmax) *rmax = max;
}

size_t
arrayMinMaxMV(size_t len, const double *array, double missval, double *rmin, double *rmax)
{
  double min = DBL_MAX;
  double max = -DBL_MAX;

  size_t nvals = 0;
  for (size_t i = 0; i < len; ++i)
    {
      if (!DBL_IS_EQUAL(array[i], missval))
        {
          if (array[i] < min) min = array[i];
          if (array[i] > max) max = array[i];
          nvals++;
        }
    }

  if (rmin) *rmin = min;
  if (rmax) *rmax = max;

  return nvals;
}

void
arrayMinMaxSum(size_t len, const double *array, double *rmin, double *rmax, double *rsum)
{
  double min = DBL_MAX;
  double max = -DBL_MAX;
  double sum = 0;

  // #pragma omp parallel for default(none) shared(min, max, array, gridsize)
  // reduction(+:mean) #pragma omp simd reduction(+:mean) reduction(min:min)
  // reduction(max:max) aligned(array:16)
  for (size_t i = 0; i < len; ++i)
    {
      if (array[i] < min) min = array[i];
      if (array[i] > max) max = array[i];
      sum += array[i];
    }

  if (rmin) *rmin = min;
  if (rmax) *rmax = max;
  if (rsum) *rsum = sum;
}

size_t
arrayMinMaxSumMV(size_t len, const double *array, double missval, double *rmin, double *rmax, double *rsum)
{
  double min = DBL_MAX;
  double max = -DBL_MAX;
  double sum = 0;

  size_t nvals = 0;
  for (size_t i = 0; i < len; ++i)
    {
      if (!DBL_IS_EQUAL(array[i], missval))
        {
          if (array[i] < min) min = array[i];
          if (array[i] > max) max = array[i];
          sum += array[i];
          nvals++;
        }
    }

  if (nvals == 0) min = missval;
  if (nvals == 0) max = missval;

  if (rmin) *rmin = min;
  if (rmax) *rmax = max;
  if (rsum) *rsum = sum;

  return nvals;
}

void
arrayMinMaxMean(size_t len, const double *array, double *rmin, double *rmax, double *rmean)
{
  double sum;
  arrayMinMaxSum(len, array, rmin, rmax, &sum);

  if (rmean) *rmean = len ? sum / (double) len : 0;
}

size_t
arrayMinMaxMeanMV(size_t len, const double *array, double missval, double *rmin, double *rmax, double *rmean)
{
  double sum;
  size_t nvals = arrayMinMaxSumMV(len, array, missval, rmin, rmax, &sum);

  if (rmean) *rmean = nvals ? sum / (double) nvals : missval;

  return nvals;
}

void
arrayMinMaxMask(size_t len, const double *array, int *mask, double *rmin, double *rmax)
{
  double xmin = DBL_MAX;
  double xmax = -DBL_MAX;

  if (mask)
    {
      for (size_t i = 0; i < len; ++i)
        {
          if (!mask[i])
            {
              if (array[i] > xmax)
                xmax = array[i];
              else if (array[i] < xmin)
                xmin = array[i];
            }
        }
    }
  else
    {
      for (size_t i = 0; i < len; ++i)
        {
          if (array[i] > xmax)
            xmax = array[i];
          else if (array[i] < xmin)
            xmin = array[i];
        }
    }

  if (rmin) *rmin = xmin;
  if (rmax) *rmax = xmax;
}

void
arrayAddArray(size_t len, double *restrict array1, const double *restrict array2)
{
  //#ifdef  _OPENMP
  //#pragma omp parallel for default(none) shared(array1,array2)
  //#endif
  for (size_t i = 0; i < len; ++i)
    array1[i] += array2[i];
}

void
arrayAddArrayMV(size_t len, double *restrict array1, const double *restrict array2, double missval)
{
  if (DBL_IS_NAN(missval))
    {
      for (size_t i = 0; i < len; i++)
        if (!DBL_IS_EQUAL(array2[i], missval))
          {
            if (!DBL_IS_EQUAL(array1[i], missval))
              array1[i] += array2[i];
            else
              array1[i] = array2[i];
          }
    }
  else
    {
      for (size_t i = 0; i < len; i++)
        if (IS_NOT_EQUAL(array2[i], missval))
          {
            if (IS_NOT_EQUAL(array1[i], missval))
              array1[i] += array2[i];
            else
              array1[i] = array2[i];
          }
    }
}

size_t
arrayNumMV(size_t len, const double *restrict array, double missval)
{
  size_t nmiss = 0;

  if (DBL_IS_NAN(missval))
    {
      for (size_t i = 0; i < len; ++i)
        if (DBL_IS_EQUAL(array[i], missval)) nmiss++;
    }
  else
    {
      for (size_t i = 0; i < len; ++i)
        if (IS_EQUAL(array[i], missval)) nmiss++;
    }

  return nmiss;
}

double
arrayMin(size_t len, const double *restrict array)
{
  assert(array != NULL);

  double min = array[0];

  for (size_t i = 0; i < len; ++i)
    if (array[i] < min) min = array[i];

  return min;
}

double
arrayMax(size_t len, const double *restrict array)
{
  assert(array != NULL);

  double max = array[0];

  for (size_t i = 0; i < len; ++i)
    if (array[i] > max) max = array[i];

  return max;
}

double
arrayRange(size_t len, const double *restrict array)
{
  assert(array != NULL);

  double min = array[0];
  double max = array[0];

  for (size_t i = 0; i < len; ++i)
    {
      if (array[i] < min) min = array[i];
      if (array[i] > max) max = array[i];
    }

  double range = max - min;

  return range;
}

double
arrayMinMV(size_t len, const double *restrict array, double missval)
{
  assert(array != NULL);

  double min = DBL_MAX;

  for (size_t i = 0; i < len; ++i)
    if (!DBL_IS_EQUAL(array[i], missval))
      if (array[i] < min) min = array[i];

  if (IS_EQUAL(min, DBL_MAX)) min = missval;

  return min;
}

double
arrayMaxMV(size_t len, const double *restrict array, double missval)
{
  assert(array != NULL);

  double max = -DBL_MAX;

  for (size_t i = 0; i < len; ++i)
    if (!DBL_IS_EQUAL(array[i], missval))
      if (array[i] > max) max = array[i];

  if (IS_EQUAL(max, -DBL_MAX)) max = missval;

  return max;
}

double
arrayRangeMV(size_t len, const double *restrict array, double missval)
{
  assert(array != NULL);

  double min = DBL_MAX;
  double max = -DBL_MAX;

  for (size_t i = 0; i < len; ++i)
    {
      if (!DBL_IS_EQUAL(array[i], missval))
        {
          if (array[i] < min) min = array[i];
          if (array[i] > max) max = array[i];
        }
    }

  double range;
  if (IS_EQUAL(min, DBL_MAX) && IS_EQUAL(max, -DBL_MAX))
    range = missval;
  else
    range = max - min;

  return range;
}

double
arraySum(size_t len, const double *restrict array)
{
  assert(array != NULL);

  double sum = 0;

  for (size_t i = 0; i < len; ++i)
    sum += array[i];

  return sum;
}

double
arraySumMV(size_t len, const double *restrict array, double missval)
{
  assert(array != NULL);

  double sum = 0;
  size_t nvals = 0;

  if (DBL_IS_NAN(missval))
    {
      for (size_t i = 0; i < len; ++i)
        if (!DBL_IS_EQUAL(array[i], missval))
          {
            sum += array[i];
            nvals++;
          }
    }
  else
    {
      for (size_t i = 0; i < len; ++i)
        if (IS_NOT_EQUAL(array[i], missval))
          {
            sum += array[i];
            nvals++;
          }
    }

  if (!nvals) sum = missval;

  return sum;
}

double
arrayMean(size_t len, const double *restrict array)
{
  assert(array != NULL);

  double sum = arraySum(len, array);

  return sum / len;
}

double
arrayMeanMV(size_t len, const double *restrict array, double missval)
{
  assert(array != NULL);

  double sum = 0, sumw = 0;

  for (size_t i = 0; i < len; ++i)
    if (!DBL_IS_EQUAL(array[i], missval))
      {
        sum += array[i];
        sumw += 1;
      }

  double missval1 = missval, missval2 = missval;
  return DIVMN(sum, sumw);
}

double
arrayWeightedMean(size_t len, const double *restrict array, const double *restrict w, double missval)
{
  assert(array != NULL);
  assert(w != NULL);

  double sum = 0, sumw = 0;

  for (size_t i = 0; i < len; ++i)
    {
      sum += w[i] * array[i];
      sumw += w[i];
    }

  return IS_EQUAL(sumw, 0.) ? missval : sum / sumw;
}

double
arrayWeightedMeanMV(size_t len, const double *restrict array, const double *restrict w, double missval)
{
  assert(array != NULL);
  assert(w != NULL);

  double missval1 = missval, missval2 = missval;
  double sum = 0, sumw = 0;

  for (size_t i = 0; i < len; ++i)
    if (!DBL_IS_EQUAL(array[i], missval1) && !DBL_IS_EQUAL(w[i], missval1))
      {
        sum += w[i] * array[i];
        sumw += w[i];
      }

  return DIVMN(sum, sumw);
}

double
arrayAvgMV(size_t len, const double *restrict array, double missval)
{
  assert(array != NULL);

  double missval1 = missval, missval2 = missval;
  double sum = 0, sumw = 0;

  for (size_t i = 0; i < len; ++i)
    {
      sum = ADDMN(sum, array[i]);
      sumw += 1;
    }

  return DIVMN(sum, sumw);
}

double
arrayWeightedAvgMV(size_t len, const double *restrict array, const double *restrict w, double missval)
{
  assert(array != NULL);
  assert(w != NULL);

  double missval1 = missval, missval2 = missval;
  double sum = 0, sumw = 0;

  for (size_t i = 0; i < len; ++i)
    if (!DBL_IS_EQUAL(w[i], missval))
      {
        sum = ADDMN(sum, MULMN(w[i], array[i]));
        sumw = ADDMN(sumw, w[i]);
      }

  return DIVMN(sum, sumw);
}

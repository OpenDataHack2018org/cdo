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

const char *fpe_errstr(int fpeRaised)
{
  const char *errstr = NULL;

  if      ( fpeRaised & FE_DIVBYZERO ) errstr = "division by zero";
  else if ( fpeRaised & FE_INEXACT   ) errstr = "inexact result";
  else if ( fpeRaised & FE_INVALID   ) errstr = "invalid result";
  else if ( fpeRaised & FE_OVERFLOW  ) errstr = "overflow";
  else if ( fpeRaised & FE_UNDERFLOW ) errstr = "underflow";

  return errstr;
}


int array_minmaxsum_val(size_t len, const double *array, double *rmin, double *rmax, double *rsum)
{
  double min = *rmin;
  double max = *rmax;
  double sum = *rsum;

  // #pragma omp parallel for default(none) shared(min, max, array, gridsize) reduction(+:mean)
  // #pragma omp simd reduction(+:mean) reduction(min:min) reduction(max:max) aligned(array:16)
  for ( size_t i = 0; i < len; ++i )
    {
      if ( array[i] < min ) min = array[i];
      if ( array[i] > max ) max = array[i];
      sum += array[i];
    }
    
  if ( rmin ) *rmin = min;
  if ( rmax ) *rmax = max;
  if ( rsum ) *rsum = sum;

  return 0;
}


int array_minmaxmean_val(size_t len, const double *array, double *rmin, double *rmax, double *rmean)
{
  // int excepts = FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW;
  // feclearexcept(FE_ALL_EXCEPT); // expensive !!!!

  double min =  DBL_MAX;
  double max = -DBL_MAX;
  double mean = 0;

  // #pragma omp parallel for default(none) shared(min, max, array, gridsize) reduction(+:mean)
  // #pragma omp simd reduction(+:mean) reduction(min:min) reduction(max:max) aligned(array:16)
  for ( size_t i = 0; i < len; ++i )
    {
      if ( array[i] < min ) min = array[i];
      if ( array[i] > max ) max = array[i];
      mean += array[i];
    }

  if ( len ) mean /= (double)len;
    
  if ( rmin ) *rmin = min;
  if ( rmax ) *rmax = max;
  if ( rmean ) *rmean = mean;

  // return fetestexcept(excepts);
  return 0;
}


int array_add_array(size_t len, double *restrict array1, const double *restrict array2)
{
  // int excepts = FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW;
  // feclearexcept(FE_ALL_EXCEPT); // expensive !!!!

  //#ifdef  _OPENMP
  //#pragma omp parallel for default(none) shared(array1,array2)
  //#endif
  for ( size_t i = 0; i < len; ++i ) array1[i] += array2[i];
  
  // return fetestexcept(excepts);
  return 0;
}


double arrayMin(size_t len, const double *restrict array)
{
  assert(array!=NULL);

  double min = array[0];

  for ( size_t i = 0; i < len; ++i )
    if ( array[i] < min ) min = array[i];

  return min;
}


double arrayMax(size_t len, const double *restrict array)
{
  assert(array!=NULL);

  double max = array[0];

  for ( size_t i = 0; i < len; ++i )
    if ( array[i] > max ) max = array[i];

  return max;
}


double arrayRange(size_t len, const double *restrict array)
{
  assert(array!=NULL);

  double min = array[0];
  double max = array[0];

  for ( size_t i = 0; i < len; ++i )
    {
      if ( array[i] < min ) min = array[i];
      if ( array[i] > max ) max = array[i];
    }

  double range = max - min;
  
  return range;
}


double arrayMinMV(size_t len, const double *restrict array, double missval)
{
  assert(array!=NULL);

  double min = DBL_MAX;

  for ( size_t i = 0; i < len; ++i )
    if ( !DBL_IS_EQUAL(array[i], missval) )
      if ( array[i] < min ) min = array[i];

  if ( IS_EQUAL(min, DBL_MAX) ) min = missval;

  return min;
}


double arrayMaxMV(size_t len, const double *restrict array, double missval)
{
  assert(array!=NULL);

  double max = -DBL_MAX;

  for ( size_t i = 0; i < len; ++i )
    if ( !DBL_IS_EQUAL(array[i], missval) )
      if ( array[i] > max ) max = array[i];

  if ( IS_EQUAL(max, -DBL_MAX) ) max = missval;

  return max;
}


double arrayRangeMV(size_t len, const double *restrict array, double missval)
{
  assert(array!=NULL);

  double min =  DBL_MAX;
  double max = -DBL_MAX;

  for ( size_t i = 0; i < len; ++i )
    {
      if ( !DBL_IS_EQUAL(array[i], missval) )
        {
          if ( array[i] < min ) min = array[i];
          if ( array[i] > max ) max = array[i];
        }
    }

  double range;
  if ( IS_EQUAL(min, DBL_MAX) && IS_EQUAL(max, -DBL_MAX) )
    range = missval;
  else
    range = max-min;
  
  return range;
}


double arraySum(size_t len, const double *restrict array)
{
  assert(array!=NULL);

  double sum = 0;

  for ( size_t i = 0; i < len; ++i ) sum += array[i];

  return sum;
}


double arraySumMV(size_t len, const double *restrict array, double missval)
{
  assert(array!=NULL);

  double sum = 0;
  size_t nvals = 0;

  if ( DBL_IS_NAN(missval) )
    {
      for ( size_t i = 0; i < len; ++i )
        if ( !DBL_IS_EQUAL(array[i], missval) )
          {
            sum += array[i];
            nvals++;
          }
    }
  else
    {
      for ( size_t i = 0; i < len; ++i )
        if ( IS_NOT_EQUAL(array[i], missval) )
          {
            sum += array[i];
            nvals++;
          }
    }

  if ( !nvals ) sum = missval;

  return sum;
}


double arrayMean(size_t len, const double *restrict array)
{
  assert(array!=NULL);

  double sum = arraySum(len, array);

  return sum/len;
}


double arrayMeanMV(size_t len, const double *restrict array, double missval)
{
  assert(array!=NULL);

  double sum = 0, sumw = 0;

  for ( size_t i = 0; i < len; ++i )
    if ( !DBL_IS_EQUAL(array[i], missval) )
      {
        sum  += array[i];
        sumw += 1;
      }

  double missval1 = missval, missval2 = missval;
  return DIVMN(sum, sumw);
}


double arrayWeightedMean(size_t len, const double *restrict array, const double *restrict w, double missval)
{
  assert(array!=NULL);
  assert(w!=NULL);

  double sum = 0, sumw = 0;

  for ( size_t i = 0; i < len; ++i ) 
    {
      sum  += w[i] * array[i];
      sumw += w[i];
    }

  return IS_EQUAL(sumw, 0.) ? missval : sum/sumw;
}


double arrayWeightedMeanMV(size_t len, const double *restrict array, const double *restrict w, double missval)
{
  assert(array!=NULL);
  assert(w!=NULL);

  double missval1 = missval, missval2 = missval;
  double sum = 0, sumw = 0;

  for ( size_t i = 0; i < len; ++i ) 
    if ( !DBL_IS_EQUAL(array[i], missval1) && !DBL_IS_EQUAL(w[i], missval1) )
      {
        sum  += w[i] * array[i];
        sumw += w[i];
      }

  return DIVMN(sum, sumw);
}


double arrayAvgMV(size_t len, const double *restrict array, double missval)
{
  assert(array!=NULL);

  double missval1 = missval, missval2 = missval;
  double sum = 0, sumw = 0;

  for ( size_t i = 0; i < len; ++i )
    {
      sum  = ADDMN(sum, array[i]);
      sumw += 1;
    }

  return DIVMN(sum, sumw);
}


double arrayWeightedAvgMV(size_t len, const double *restrict array, const double *restrict w, double missval)
{
  assert(array!=NULL);
  assert(w!=NULL);

  double missval1 = missval, missval2 = missval;
  double sum = 0, sumw = 0;

  for ( size_t i = 0; i < len; ++i ) 
    if ( !DBL_IS_EQUAL(w[i], missval) )
      {
        sum  = ADDMN(sum, MULMN(w[i], array[i]));
        sumw = ADDMN(sumw, w[i]);
      }

  return DIVMN(sum, sumw);
}

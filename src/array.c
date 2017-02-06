#include <stdio.h>
#include <float.h>
#include <fenv.h>

#include "compare.h"

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


int array_minmaxmean_val(size_t len, const double *array, double *rmin, double *rmax, double *rmean)
{
  int excepts = FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW;
  feclearexcept(FE_ALL_EXCEPT);

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

  return fetestexcept(excepts);
}


int array_mean_val_weighted(size_t len, const double *restrict array, const double *restrict w, double missval, double *rmean)
{
  int excepts = FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW;
  feclearexcept(FE_ALL_EXCEPT);

  double rsum = 0, rsumw = 0;

  for ( size_t i = 0; i < len; i++ ) 
    {
      rsum  += w[i] * array[i];
      rsumw += w[i];
    }

  *rmean = DBL_IS_EQUAL(rsumw, 0.) ? missval : rsum/rsumw;

  return fetestexcept(excepts);
}


int array_add_array(size_t len, double *restrict array1, const double *restrict array2)
{
  int excepts = FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW;
  feclearexcept(FE_ALL_EXCEPT);

  //#if defined(_OPENMP)
  //#pragma omp parallel for default(none) shared(array1,array2)
  //#endif
  for ( size_t i = 0; i < len; i++ ) array1[i] += array2[i];

  return fetestexcept(excepts);
}

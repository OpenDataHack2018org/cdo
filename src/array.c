#include <stdio.h>
#include <float.h>
#include <fenv.h>

#pragma STDC FENV_ACCESS ON

const char *fpe_errstr(int fpeRaised)
{
  const char *errstr = NULL;

  if      ( fpeRaised & FE_INEXACT   ) errstr = "inexact result";
  else if ( fpeRaised & FE_INVALID   ) errstr = "invalid result";
  else if ( fpeRaised & FE_DIVBYZERO ) errstr = "division by zero";
  else if ( fpeRaised & FE_OVERFLOW  ) errstr = "overflow";
  else if ( fpeRaised & FE_UNDERFLOW ) errstr = "underflow";

  return errstr;
}


int array_minmaxmean_val(const double *array, size_t len, double *rmin, double *rmax, double *rmean)
{
  double min =  DBL_MAX;
  double max = -DBL_MAX;
  double mean = 0;

  feclearexcept(FE_ALL_EXCEPT);

  // #pragma omp parallel for default(none) shared(min, max, array, gridsize) reduction(+:mean)
  // #pragma omp simd reduction(+:mean) reduction(min:min) reduction(max:max) aligned(array:16)
  for ( size_t i = 0; i < len; ++i )
    {
      if ( array[i] < min ) min = array[i];
      if ( array[i] > max ) max = array[i];
      mean += array[i];
    }

  int fpeRaised = fetestexcept(FE_ALL_EXCEPT);
    
  *rmin = min;
  *rmax = max;
  *rmean = mean;

  return fpeRaised;
}

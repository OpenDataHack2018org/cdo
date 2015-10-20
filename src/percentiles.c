#include <math.h>
#include <string.h>
#include "util.h"
#include "percentiles.h"
#include "nth_element.h"

enum percentile_methods {NRANK=1, NIST, NUMPY};
enum interpolation_methods {LINEAR=1, LOWER, HIGHER, NEAREST};

static int percentile_method = NRANK;
static int interpolation_method = LINEAR;

static
double percentile_nrank(double *array, int len, double pn)
{
  int irank = (int)ceil(len*(pn/100.0));
  if ( irank <   1 ) irank = 1;
  if ( irank > len ) irank = len;
  return nth_element(array, len, irank-1);
}

static
double percentile_nist(double *array, int len, double pn)
{
  double rank = (len+1)*(pn/100.0);
  int k = (int) rank;
  double d = rank - k;
  double percentil = 0;
  if      ( k ==   0 ) percentil = nth_element(array, len, 0);
  else if ( k >= len ) percentil = nth_element(array, len, len-1);
  else
    {
      double vk1 = nth_element(array, len, k);
      double vk  = array[k-1];
      percentil = vk + d*(vk1 - vk);
    }

  return percentil;
}

static
double percentile_numpy(double *array, int len, double pn)
{
  double rank = (len-1)*(pn/100.0) + 1;
  int k = (int) rank;
  double d = rank - k;
  double percentil = 0;
  if      ( k ==   0 ) percentil = nth_element(array, len, 0);
  else if ( k >= len ) percentil = nth_element(array, len, len-1);
  else
    {
      if ( interpolation_method == LINEAR )
        {
          double vk1 = nth_element(array, len, k);
          double vk  = array[k-1];
          percentil = vk + d*(vk1 - vk);
        }
      else
        {
          int irank = 0;
          if      ( interpolation_method == LOWER   ) irank = (int) rank;
          else if ( interpolation_method == HIGHER  ) irank = (int) rank + 1;
          else if ( interpolation_method == NEAREST ) irank = (int) lround(rank);

          if ( irank <   1 ) irank = 1;
          if ( irank > len ) irank = len;

          percentil = nth_element(array, len, irank-1);
        }
    }

  return percentil;
}


double percentile(double *array, int len, double pn)
{
  double percentil = 0;
  
  if      ( percentile_method == NRANK ) percentil = percentile_nrank(array, len, pn);
  else if ( percentile_method == NIST  ) percentil = percentile_nist(array, len, pn);
  else if ( percentile_method == NUMPY ) percentil = percentile_numpy(array, len, pn);
  else cdoAbort("Internal error: percentile method %d not implemented!", percentile_method);

  return percentil;
}


void percentile_set_method(const char *methodstr)
{
  char *methodname = strdup(methodstr);
  strtolower(methodname);

  if      ( strcmp("nrank", methodname) == 0 ) percentile_method = NRANK;
  else if ( strcmp("nist",  methodname) == 0 ) percentile_method = NIST;
  else if ( strcmp("numpy", methodname) == 0 ) percentile_method = NUMPY;
  else if ( strcmp("numpy_linear",  methodname) == 0 ) {percentile_method = NUMPY; interpolation_method = LINEAR;}
  else if ( strcmp("numpy_lower",   methodname) == 0 ) {percentile_method = NUMPY; interpolation_method = LOWER;}
  else if ( strcmp("numpy_higher",  methodname) == 0 ) {percentile_method = NUMPY; interpolation_method = HIGHER;}
  else if ( strcmp("numpy_nearest", methodname) == 0 ) {percentile_method = NUMPY; interpolation_method = NEAREST;}
  else cdoAbort("Percentile method %s not available!", methodstr);
}


void percentile_check_number(double pn)
{
  if ( pn < 0 || pn > 100 )
    cdoAbort("Illegal argument: percentile number %g is not in the range 0..100!", pn);
}

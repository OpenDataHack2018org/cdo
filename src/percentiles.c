#include <math.h>
#include <string.h>
#include "util.h"
#include "percentiles.h"
#include "nth_element.h"

enum percentile_methods {NEAREST_RANK, NIST, EXCEL};

static int percentile_method = NEAREST_RANK;


double percentile_nearest_rank(double *array, int len, double pn)
{
  int rank = (int)ceil(len*(pn/100.0));
  if ( rank <   1 ) rank = 1;
  if ( rank > len ) rank = len;
  double percentil = nth_element(array, len, rank-1);
  return percentil;
}


double percentile(double *array, int len, double pn)
{
  double percentil = 0;
  
  if ( percentile_method ) percentil = percentile_nearest_rank(array, len, pn);
  else cdoAbort("Internal error: percentile method %d not implemented", percentile_method);

  return percentil;
}


void percentile_set_method(const char *methodstr)
{
  char *methodname = strdup(methodstr);
  strtolower(methodname);

  if      ( strcmp("nearest_rank", methodname) == 0 ) percentile_method = NEAREST_RANK;
  else if ( strcmp("nist",         methodname) == 0 ) percentile_method = NIST;
  else if ( strcmp("excel",        methodname) == 0 ) percentile_method = EXCEL;
  else cdoAbort("Percentile method %s not available");
}


void percentile_check_number(double pn)
{
  if ( pn < 0 || pn > 100 )
    cdoAbort("Illegal argument: percentile number %g is not in the range 0..100!", pn);
}

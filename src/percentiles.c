#include <math.h>
#include "util.h"
#include "percentiles.h"
#include "nth_element.h"


double percentile(double *array, unsigned len, double pn)
{
  int rank = (int)ceil(len*(pn/100.0));
  if ( rank <   1 ) rank = 1;
  if ( rank > len ) rank = len;
  double percentil = nth_element(array, len, rank-1);
  return percentil;
}


void percentile_check_number(double pn)
{
  if ( pn < 0 || pn > 100 )
    cdoAbort("Illegal argument: percentile number %g is not in the range 0..100!", pn);
}

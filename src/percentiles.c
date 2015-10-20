#include <math.h>
#include "util.h"
#include "percentiles.h"
#include "nth_element.h"


double percentile(double *array, unsigned len, double pn)
{
  int element = (int)ceil(len*(pn/100.0))-1;
  if ( element <    0 ) element = 0;
  if ( element >= len ) element = len-1;
  double percentil = nth_element(array, len, element);
  return percentil;
}


void percentile_check_number(double pn)
{
  if ( pn < 0 || pn > 100 )
    cdoAbort("Illegal argument: percentile number %g is not in the range 0..100!", pn);
}

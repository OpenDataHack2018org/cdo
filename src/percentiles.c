#include <math.h>
#include "percentiles.h"
#include "nth_element.h"

double percentile(double *array, unsigned len, double p)
{
  double percentil;
  percentil = nth_element(array, len, (int)ceil(len*(p/100.0))-1);
  return percentil;
}

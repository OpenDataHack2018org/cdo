#ifdef _OPENMP
#include <omp.h>  // omp_get_wtime
#endif

inline double
cdo_get_wtime()
{
#ifdef _OPENMP
  return omp_get_wtime();
#else
  return 0;
#endif
}

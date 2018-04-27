#define USE_OMP_WTIME

#ifdef USE_OMP_WTIME

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

#else

#include <chrono>

inline double
cdo_get_wtime()
{
  return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count() / 1000.;
}

#endif

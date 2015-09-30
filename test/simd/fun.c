void fun1(const unsigned nelem, double *restrict array1, const double *restrict array2)
{
#if defined(_OPENMP)
#pragma  omp simd aligned(array1:64) aligned(array2:64)
#endif
  for ( unsigned i = 0; i < nelem; ++i )
    array1[i] += array2[i];
}

void fun2(const unsigned nelem, double *restrict array1, const double *restrict array2, const double *restrict array3)
{
#if defined(_OPENMP)
#pragma  omp simd aligned(array1:64) aligned(array2:64) aligned(array3:64)
#endif
  for ( unsigned i = 0; i < nelem; ++i )
    array1[i] += array2[i]*array3[i];
}

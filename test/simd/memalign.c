// icc -std=c99 -O3 -qopt-report=5 -march=native -openmp memalign.c fun.c

#include <stdio.h>
#include <stdlib.h>

void fun1(int nelem, double *restrict array1, double *restrict array2);

int get_alignment(double *ptr)
{
  int64_t iptr = (int64_t) ptr;
  int mk[4] = {64, 32, 16, 8};
  int malign = 0;

  for ( int i = 0; i < 4; ++i )
    if ( iptr%mk[i] == 0 )
      {
	malign = mk[i];
	break;
      }
  
  return malign;
}

int main(void)
{
  int nelem = 97;

  double *array1 = (double *) malloc(nelem*sizeof(double));
  double *array2 = (double *) malloc(nelem*sizeof(double));

  printf("mem alignment: %d %d\n", get_alignment(array1), get_alignment(array2));

  for ( int i = 0; i < nelem; ++i ) array1[i] = 0;
  for ( int i = 0; i < nelem; ++i ) array2[i] = 1;

  for ( int i = 0; i < 200000000; ++i )
    fun1(nelem, array1, array2);

  free(array1);
  free(array2);
  
  return 0;
}

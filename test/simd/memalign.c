// icc -std=c99 -O3 -march=native -qopt-report=5 -openmp memalign.c fun.c
// gcc -std=c99 -O3 -march=native -ftree-vectorize -fdump-tree-vect-blocks -fopt-info-optimized -fopenmp memalign.c fun.c

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600
#endif
/*
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200112L
#endif
*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

void fun1(int nelem, double *restrict array1, double *restrict array2);

void print_opt(void)
{
  fprintf(stderr, "opt: ");
#if defined(__AVX2__)
  fprintf(stderr, " AVX2");
#elif defined(__AVX__)
  fprintf(stderr, " AVX");
#elif defined(__SSE4_2__)
  fprintf(stderr, " SSE4_2");
#elif defined(__SSE4_1__)
  fprintf(stderr, " SSE4_1");
#elif defined(__SSE3__)
  fprintf(stderr, " SSE3");
#elif defined(__SSE2__)
  fprintf(stderr, " SSE2");
#endif 
#if defined(_OPENMP)
  fprintf(stderr, " OPENMP");  
#endif
  fprintf(stderr, "\n");
}

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
  double *array1, *array2;

  print_opt();

  malloc(nelem*sizeof(double));
  //array1 = (double *) malloc(nelem*sizeof(double));
  posix_memalign((void **)&array1, 64, nelem*sizeof(double));
  malloc(nelem*sizeof(double));
  //array2 = (double *) malloc(nelem*sizeof(double));
  posix_memalign((void **)&array2, 64, nelem*sizeof(double));

  printf("mem alignment: %d %d\n", get_alignment(array1), get_alignment(array2));

  for ( int i = 0; i < nelem; ++i ) array1[i] = 0;
  for ( int i = 0; i < nelem; ++i ) array2[i] = 1;

  for ( int i = 0; i < 200000000; ++i )
    fun1(nelem, array1, array2);

  free(array1);
  free(array2);
  
  return 0;
}

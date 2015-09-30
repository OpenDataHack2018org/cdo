// aligned access gives 3-5% speedup
// icc -g -std=c99 -O2 -xCORE-AVX2 -qopt-report=5 -openmp memalign.c fun.c
// gcc -g -std=c99 -O3 -march=native -ftree-vectorize -fdump-tree-vect-blocks -fopt-info-optimized -fopenmp memalign.c fun.c

/*
ICC16/hama2:
        fun1    fun2
SSE3    12.7    16.8  unaligned
SSE3    12.4    13.7    aligned
AVX2    10.4    10.4  unaligned
AVX2    10.4    10.4    aligned
 */
/*
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600
#endif
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200112L
#endif
*/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

//#define NITER 200000000
#define NITER 2000000

void fun1(const unsigned nelem, double *restrict array1, double *restrict array2);
void fun2(const unsigned nelem, double *restrict array1, double *restrict array2, double *restrict array3);

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
  double start_time;
  //  int nelem = 97;
  int nelem = 4*4096;
  double *array1, *array2;

  print_opt();

  malloc(nelem*sizeof(double));
  array1 = (double *) malloc(nelem*sizeof(double));
  //posix_memalign((void **)&array1, 64, nelem*sizeof(double));
  malloc(nelem*sizeof(double));
  array2 = (double *) malloc(nelem*sizeof(double));
  //posix_memalign((void **)&array2, 64, nelem*sizeof(double));

  printf("mem alignment: %d %d\n", get_alignment(array1), get_alignment(array2));

  for ( int i = 0; i < nelem; ++i ) array1[i] = 0;
  for ( int i = 0; i < nelem; ++i ) array2[i] = 1;

  start_time = omp_get_wtime();
  for ( int i = 0; i < NITER; ++i )
    fun1(nelem, array1, array2);
  printf("\n fun1 in %lf seconds\n ",omp_get_wtime() - start_time);

  start_time = omp_get_wtime();
  for ( int i = 0; i < NITER; ++i )
    fun2(nelem, array1, array2, array2);
  printf("\n fun2 in %lf seconds\n ",omp_get_wtime() - start_time);

  free(array1);
  free(array2);
  
  return 0;
}

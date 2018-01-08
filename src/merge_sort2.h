#ifndef _MERGE_SORT2_H_
#define _MERGE_SORT2_H_

/* MERGE SORT DEFINES */

#include "cdo_int.h"

#include <time.h>

#ifdef  _OPENMP
#include <omp.h>
#endif

void sort_iter_single(long num_links, double *restrict add1, int parent);

#endif

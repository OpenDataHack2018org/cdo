#ifndef _GRID_SEARCH_H_
#define _GRID_SEARCH_H_

#include "kdtreelib/kdtree.h"

struct gridsearch;

struct gridsearch *gridsearch_create(unsigned n, const double *restrict lons, const double *restrict lats);
struct gridsearch *gridsearch_index_create(unsigned n, const double *restrict lons, const double *restrict lats, const unsigned *restrict index);
void gridsearch_delete(struct gridsearch *gs);
void *gridsearch_nearest(struct gridsearch *gs, double lon, double lat, double *range);
unsigned gridsearch_item(void *gs_result);
struct pqueue *gridsearch_qnearest(struct gridsearch *gs, double lon, double lat, double *prange, unsigned nnn);

#endif

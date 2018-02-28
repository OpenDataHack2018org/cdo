#ifndef NBR_WEIGHTS_H
#define NBR_WEIGHTS_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <vector>

struct nbrWeightsType
{
  std::vector<uint8_t> mask;
  std::vector<size_t> add;
  std::vector<double> dist;
};

void
nbr_store_distance(size_t nadd, double distance, size_t numNeighbors, size_t *restrict nbr_add, double *restrict nbr_dist);
void
nbr_check_distance(size_t numNeighbors, const size_t *restrict nbr_add, double *restrict nbr_dist);
double
nbr_compute_weights(size_t numNeighbors, const int *restrict src_grid_mask, uint8_t *restrict nbr_mask, const size_t *restrict nbr_add, double *restrict nbr_dist);
size_t
nbr_normalize_weights(size_t numNeighbors, double dist_tot, const uint8_t *restrict nbr_mask, size_t *restrict nbr_add, double *restrict nbr_dist);

#endif

/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2018 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/
#ifdef _OPENMP
#include <omp.h>  // omp_get_wtime
#endif

#include "cdo_int.h"
#include "remap.h"
#include "remap_store_link.h"
#include "cdoOptions.h"

// Interpolation using a distance-weighted average

// This routine computes the inverse-distance weights for a nearest-neighbor interpolation
void
remapDistwgtWeights(size_t numNeighbors, RemapSearch &rsearch, RemapVars &rv)
{
  RemapGrid *src_grid = rsearch.srcGrid;
  RemapGrid *tgt_grid = rsearch.tgtGrid;

  if (cdoVerbose) cdoPrint("Called %s()", __func__);

  progressInit();

  // Compute mappings from source to target grid

  size_t tgt_grid_size = tgt_grid->size;

  std::vector<WeightLinks> weightLinks(tgt_grid_size);
  weightLinksAlloc(numNeighbors, tgt_grid_size, weightLinks);

  std::vector<knnWeightsType> knnWeights;
  for (int i = 0; i < Threading::ompNumThreads; ++i) knnWeights.push_back(knnWeightsType(numNeighbors));

#ifdef _OPENMP
  double start = cdoVerbose ? omp_get_wtime() : 0;
#endif

  // Loop over destination grid

  double findex = 0;

#ifdef HAVE_OPENMP4
#pragma omp parallel for default(none) reduction(+ : findex) shared(rsearch, weightLinks, numNeighbors, src_grid, tgt_grid, \
                                                                    tgt_grid_size, knnWeights)
#endif
  for (size_t tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add)
    {
      int ompthID = cdo_omp_get_thread_num();

      findex++;
      if (ompthID == 0) progressStatus(0, 1, findex / tgt_grid_size);

      weightLinks[tgt_cell_add].nlinks = 0;

      if (!tgt_grid->mask[tgt_cell_add]) continue;

      double plon = 0, plat = 0;
      remapgrid_get_lonlat(tgt_grid, tgt_cell_add, &plon, &plat);

      // Find nearest grid points on source grid and distances to each point
      remapSearchPoints(rsearch, plon, plat, knnWeights[ompthID]);

      // Compute weights based on inverse distance if mask is false, eliminate those points
      size_t nadds = knnWeights[ompthID].compute_weights(&src_grid->mask[0]);

      for (size_t n = 0; n < nadds; ++n)
        if (knnWeights[ompthID].m_mask[n]) tgt_grid->cell_frac[tgt_cell_add] = 1.0;

      // Store the link
      storeWeightLinks(0, nadds, &knnWeights[ompthID].m_addr[0], &knnWeights[ompthID].m_dist[0], tgt_cell_add, weightLinks);
    }

  progressStatus(0, 1, 1);

  if (rsearch.gs) gridsearch_delete(rsearch.gs);
  rsearch.gs = NULL;

  weightLinksToRemapLinks(0, tgt_grid_size, weightLinks, rv);

#ifdef _OPENMP
  if (cdoVerbose) cdoPrint("Point search nearest: %.2f seconds", omp_get_wtime() - start);
#endif
}  // remapDistwgtWeights

void
remapDistwgt(size_t numNeighbors, RemapSearch &rsearch, const double *restrict src_array, double *restrict tgt_array,
             double missval)
{
  RemapGrid *src_grid = rsearch.srcGrid;
  RemapGrid *tgt_grid = rsearch.tgtGrid;

  if (cdoVerbose) cdoPrint("Called %s()", __func__);

  progressInit();

  // Compute mappings from source to target grid

  size_t tgt_grid_size = tgt_grid->size;

  std::vector<knnWeightsType> knnWeights;
  for (int i = 0; i < Threading::ompNumThreads; ++i) knnWeights.push_back(knnWeightsType(numNeighbors));

#ifdef _OPENMP
  double start = cdoVerbose ? omp_get_wtime() : 0;
#endif

  // Loop over destination grid

  double findex = 0;

#ifdef HAVE_OPENMP4
#pragma omp parallel for default(none) reduction(+ : findex) shared(rsearch, numNeighbors, src_grid, tgt_grid, tgt_grid_size, \
                                                                    src_array, tgt_array, missval, knnWeights)
#endif
  for (size_t tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add)
    {
      int ompthID = cdo_omp_get_thread_num();

      findex++;
      if (ompthID == 0) progressStatus(0, 1, findex / tgt_grid_size);

      tgt_array[tgt_cell_add] = missval;

      if (!tgt_grid->mask[tgt_cell_add]) continue;

      double plon = 0, plat = 0;
      remapgrid_get_lonlat(tgt_grid, tgt_cell_add, &plon, &plat);

      // Find nearest grid points on source grid and distances to each point
      remapSearchPoints(rsearch, plon, plat, knnWeights[ompthID]);

      // Compute weights based on inverse distance if mask is false, eliminate those points
      size_t nadds = knnWeights[ompthID].compute_weights(&src_grid->mask[0]);
      if (nadds) tgt_array[tgt_cell_add] = knnWeights[ompthID].array_weights_sum(src_array);
    }

  progressStatus(0, 1, 1);

#ifdef _OPENMP
  if (cdoVerbose) cdoPrint("Point search nearest: %.2f seconds", omp_get_wtime() - start);
#endif
}  // remapDistwgt

void remapInit(remapType &remap);

void
intgriddis(field_type *field1, field_type *field2, size_t numNeighbors)
{
  RemapMethod mapType = RemapMethod::DISTWGT;
  int gridID1 = field1->grid;
  int gridID2 = field2->grid;
  double src_missval = field1->missval;
  double tgt_missval = field2->missval;
  double *src_array = field1->ptr;
  double *tgt_array = field2->ptr;

  if (cdoVerbose) cdoPrint("Called %s()", __func__);

  progressInit();

  // Interpolate from source to target grid

  remapType remap;
  remapInit(remap);

  bool remap_extrapolate = false;
  remapInitGrids(mapType, remap_extrapolate, gridID1, remap.src_grid, gridID2, remap.tgt_grid);

  size_t src_grid_size = remap.src_grid.size;
  size_t tgt_grid_size = remap.tgt_grid.size;

  std::vector<int> src_mask(src_grid_size);
  for (size_t i = 0; i < src_grid_size; ++i) src_mask[i] = !DBL_IS_EQUAL(src_array[i], src_missval);
  std::vector<int> tgt_mask(tgt_grid_size);
  for (size_t i = 0; i < tgt_grid_size; ++i) tgt_mask[i] = 1;

  std::vector<knnWeightsType> knnWeights;
  for (int i = 0; i < Threading::ompNumThreads; ++i) knnWeights.push_back(knnWeightsType(numNeighbors));

  remapSearchInit(mapType, remap.search, remap.src_grid, remap.tgt_grid);

#ifdef _OPENMP
  double start = cdoVerbose ? omp_get_wtime() : 0;
#endif

  // Loop over destination grid

  size_t nmiss = 0;
  double findex = 0;

  for (size_t tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add)
    {
      int ompthID = cdo_omp_get_thread_num();

      findex++;
      if (ompthID == 0) progressStatus(0, 1, findex / tgt_grid_size);

      tgt_array[tgt_cell_add] = tgt_missval;

      if (!tgt_mask[tgt_cell_add]) continue;

      double plon = 0, plat = 0;
      remapgrid_get_lonlat(&remap.tgt_grid, tgt_cell_add, &plon, &plat);

      // Find nearest grid points on source grid and distances to each point
      remapSearchPoints(remap.search, plon, plat, knnWeights[ompthID]);

      // Compute weights based on inverse distance if mask is false, eliminate those points
      size_t nadds = knnWeights[ompthID].compute_weights(&src_mask[0]);
      if (nadds)
        tgt_array[tgt_cell_add] = knnWeights[ompthID].array_weights_sum(src_array);
      else
        nmiss++;
    }

  progressStatus(0, 1, 1);

  field2->nmiss = nmiss;

  remapGridFree(remap.src_grid);
  remapGridFree(remap.tgt_grid);
  remapSearchFree(remap.search);

#ifdef _OPENMP
  if (cdoVerbose) cdoPrint("Point search nearest: %.2f seconds", omp_get_wtime() - start);
#endif
}  // intgriddis

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
#include "grid.h"
#include "remap.h"
#include "remap_store_link.h"
#include "grid_search.h"
#include "cdoOptions.h"

// Interpolation using a distance-weighted average

// This routine computes the inverse-distance weights for a nearest-neighbor interpolation.
void
remap_distwgt_weights(size_t numNeighbors, RemapSearch &rsearch, RemapGrid *src_grid, RemapGrid *tgt_grid, RemapVars &rv)
{
  if (cdoVerbose) cdoPrint("Called %s()", __func__);

  progressInit();

  // Compute mappings from source to target grid

  size_t tgt_grid_size = tgt_grid->size;

  std::vector<weightlinks_t> weightlinks(tgt_grid_size);
  weightlinks[0].addweights = (addweight_t *) Malloc(numNeighbors * tgt_grid_size * sizeof(addweight_t));
  for (size_t tgt_cell_add = 1; tgt_cell_add < tgt_grid_size; ++tgt_cell_add)
    weightlinks[tgt_cell_add].addweights = weightlinks[0].addweights + numNeighbors * tgt_cell_add;

  std::vector<knnWeightsType> knnWeights;
  for (int i = 0; i < Threading::ompNumThreads; ++i) knnWeights.push_back(knnWeightsType(numNeighbors));

#ifdef _OPENMP
  double start = cdoVerbose ? omp_get_wtime() : 0;
#endif

  // Loop over destination grid

  double findex = 0;

#ifdef HAVE_OPENMP4
#pragma omp parallel for default(none) reduction(+ : findex) shared(rsearch, weightlinks, numNeighbors, \
                                                                    src_grid, tgt_grid, tgt_grid_size, knnWeights)
#endif
  for (size_t tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add)
    {
      int ompthID = cdo_omp_get_thread_num();

      findex++;
      if (ompthID == 0) progressStatus(0, 1, findex / tgt_grid_size);

      weightlinks[tgt_cell_add].nlinks = 0;

      if (!tgt_grid->mask[tgt_cell_add]) continue;

      double plon = 0, plat = 0;
      remapgrid_get_lonlat(tgt_grid, tgt_cell_add, &plon, &plat);

      // Find nearest grid points on source grid and distances to each point
      remapSearchPoints(rsearch, plon, plat, knnWeights[ompthID]);

      // Compute weights based on inverse distance if mask is false, eliminate those points
      size_t nadds = knnWeights[ompthID].compute_weights(src_grid->mask);

      for (size_t n = 0; n < nadds; ++n)
        if (knnWeights[ompthID].m_mask[n]) tgt_grid->cell_frac[tgt_cell_add] = ONE;

      // Store the link
      store_weightlinks(0, nadds, &knnWeights[ompthID].m_addr[0], &knnWeights[ompthID].m_dist[0], tgt_cell_add,
                        &weightlinks[0]);
    }

  progressStatus(0, 1, 1);

  if (rsearch.gs) gridsearch_delete(rsearch.gs);
  rsearch.gs = NULL;

  weightlinks2remaplinks(0, tgt_grid_size, &weightlinks[0], rv);

#ifdef _OPENMP
  if (cdoVerbose) printf("gridsearch nearest: %.2f seconds\n", omp_get_wtime() - start);
#endif
}  // remap_distwgt_weights

void
remap_distwgt(size_t numNeighbors, RemapSearch &rsearch, RemapGrid *src_grid, RemapGrid *tgt_grid, const double *restrict src_array,
              double *restrict tgt_array, double missval)
{
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
#pragma omp parallel for default(none)                                                      \
  reduction(+ : findex) shared(rsearch, numNeighbors, src_grid, tgt_grid, \
                               tgt_grid_size) shared(src_array, tgt_array, missval, knnWeights)
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
      size_t nadds = knnWeights[ompthID].compute_weights(src_grid->mask);
      if (nadds) tgt_array[tgt_cell_add] = knnWeights[ompthID].array_weights_sum(src_array);
    }

  progressStatus(0, 1, 1);

#ifdef _OPENMP
  if (cdoVerbose) printf("gridsearch nearest: %.2f seconds\n", omp_get_wtime() - start);
#endif
}  // remap_distwgt

#include <cdi.h>

void
intgriddis(field_type *field1, field_type *field2, size_t numNeighbors)
{
  int gridID1 = field1->grid;
  int gridID2 = field2->grid;
  double src_missval = field1->missval;
  double tgt_missval = field2->missval;
  double *src_array = field1->ptr;
  double *tgt_array = field2->ptr;

  if (cdoVerbose) cdoPrint("Called %s()", __func__);

  progressInit();

  // Compute mappings from source to target grid

  int src_gridtype = gridInqType(gridID1);
  int tgt_gridtype = gridInqType(gridID2);
  if (src_gridtype != GRID_CURVILINEAR && src_gridtype != GRID_UNSTRUCTURED)
    cdoAbort("Source grid must be curvilinear or unstructured!");
  if (tgt_gridtype != GRID_CURVILINEAR && tgt_gridtype != GRID_UNSTRUCTURED)
    cdoAbort("Target grid must be curvilinear or unstructured!");

  size_t src_grid_size = gridInqSize(gridID1);
  size_t tgt_grid_size = gridInqSize(gridID2);

  int *src_mask = (int *) Malloc(src_grid_size * sizeof(int));
  for (size_t i = 0; i < src_grid_size; ++i) src_mask[i] = !DBL_IS_EQUAL(src_array[i], src_missval);
  int *tgt_mask = (int *) Malloc(tgt_grid_size * sizeof(int));
  for (size_t i = 0; i < tgt_grid_size; ++i) tgt_mask[i] = 1;

  double *src_cell_center_lon = (double *) Malloc(src_grid_size * sizeof(double));
  double *src_cell_center_lat = (double *) Malloc(src_grid_size * sizeof(double));
  gridInqXvals(gridID1, src_cell_center_lon);
  gridInqYvals(gridID1, src_cell_center_lat);

  double *tgt_cell_center_lon = (double *) Malloc(tgt_grid_size * sizeof(double));
  double *tgt_cell_center_lat = (double *) Malloc(tgt_grid_size * sizeof(double));
  gridInqXvals(gridID2, tgt_cell_center_lon);
  gridInqYvals(gridID2, tgt_cell_center_lat);

  char xunits[CDI_MAX_NAME];
  xunits[0] = 0;
  char yunits[CDI_MAX_NAME];
  yunits[0] = 0;
  cdiGridInqKeyStr(gridID1, CDI_KEY_XUNITS, CDI_MAX_NAME, xunits);
  cdiGridInqKeyStr(gridID1, CDI_KEY_YUNITS, CDI_MAX_NAME, yunits);
  grid_to_radian(xunits, src_grid_size, src_cell_center_lon, "src cell center lon");
  grid_to_radian(yunits, src_grid_size, src_cell_center_lat, "src cell center lat");
  xunits[0] = 0;
  yunits[0] = 0;
  cdiGridInqKeyStr(gridID2, CDI_KEY_XUNITS, CDI_MAX_NAME, xunits);
  cdiGridInqKeyStr(gridID2, CDI_KEY_YUNITS, CDI_MAX_NAME, yunits);
  grid_to_radian(xunits, tgt_grid_size, tgt_cell_center_lon, "tgt cell center lon");
  grid_to_radian(yunits, tgt_grid_size, tgt_cell_center_lat, "tgt cell center lat");

  std::vector<knnWeightsType> knnWeights;
  for (int i = 0; i < Threading::ompNumThreads; ++i) knnWeights.push_back(knnWeightsType(numNeighbors));

#ifdef _OPENMP
  double start = cdoVerbose ? omp_get_wtime() : 0;
#endif

  bool xIsCyclic = gridIsCircular(gridID1);
  size_t dims[2] = { src_grid_size, 0 };
  GridSearch *gs = NULL;
  // if ( src_remap_grid_type == REMAP_GRID_TYPE_REG2D )
  //  gs = gridsearch_create_reg2d(xIsCyclic, dims, src_grid->reg2d_center_lon,
  //  src_grid->reg2d_center_lat);
  gs = gridsearch_create(xIsCyclic, dims, src_grid_size, src_cell_center_lon, src_cell_center_lat);

// if ( src_grid->lextrapolate ) gridsearch_extrapolate(gs);
// gridsearch_extrapolate(gs);

#ifdef _OPENMP
  if (cdoVerbose) printf("gridsearch created: %.2f seconds\n", omp_get_wtime() - start);
  if (cdoVerbose) start = omp_get_wtime();
#endif

  // Loop over destination grid

  size_t nmiss = 0;
  double findex = 0;

#ifdef HAVE_OPENMP4
/*
#pragma omp parallel for default(none)  reduction(+:findex) \
shared(gs, numNeighbors, src_grid, tgt_grid, tgt_grid_size)  \
shared(src_array, tgt_array, missval, knnWeights)
*/
#endif
  for (size_t tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add)
    {
      int ompthID = cdo_omp_get_thread_num();

      findex++;
      if (ompthID == 0) progressStatus(0, 1, findex / tgt_grid_size);

      tgt_array[tgt_cell_add] = tgt_missval;

      if (!tgt_mask[tgt_cell_add]) continue;

      double plon = 0, plat = 0;
      // remapgrid_get_lonlat(tgt_grid, tgt_cell_add, &plon, &plat);
      plat = tgt_cell_center_lat[tgt_cell_add];
      plon = tgt_cell_center_lon[tgt_cell_add];

      // Find nearest grid points on source grid and distances to each point
      // if ( src_remap_grid_type == REMAP_GRID_TYPE_REG2D )
      //   grid_search_nbr_reg2d(gs, numNeighbors, nbr_add[ompthID],
      //   nbr_dist[ompthID], plon, plat);
      // else
      grid_search_nbr(gs, plon, plat, knnWeights[ompthID]);

      // Compute weights based on inverse distance if mask is false, eliminate those points
      size_t nadds = knnWeights[ompthID].compute_weights(src_mask);
      if (nadds)
        tgt_array[tgt_cell_add] = knnWeights[ompthID].array_weights_sum(src_array);
      else
        nmiss++;
    }

  progressStatus(0, 1, 1);

  field2->nmiss = nmiss;

  if (gs) gridsearch_delete(gs);

  Free(src_mask);
  Free(tgt_mask);
  Free(src_cell_center_lon);
  Free(src_cell_center_lat);
  Free(tgt_cell_center_lon);
  Free(tgt_cell_center_lat);

#ifdef _OPENMP
  if (cdoVerbose) printf("gridsearch nearest: %.2f seconds\n", omp_get_wtime() - start);
#endif
}  // intgriddis

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

//  Interpolation using a distance-weighted average

//  This routine finds the closest num_neighbor points to a search point and
//  computes a distance to each of the neighbors.

static inline void
LLtoXYZ(double lon, double lat, double *restrict xyz)
{
  double cos_lat = cos(lat);
  xyz[0] = cos_lat * cos(lon);
  xyz[1] = cos_lat * sin(lon);
  xyz[2] = sin(lat);
}

static constexpr double
square(const double x) noexcept
{
  return x * x;
}

static constexpr double
distance(const double *restrict a, const double *restrict b) noexcept
{
  return square(a[0] - b[0]) + square(a[1] - b[1]) + square(a[2] - b[2]);
}

#define MAX_SEARCH_CELLS 25
static void
grid_search_nbr_reg2d(GridSearch *gs, knnWeightsType &knnWeights, double plon, double plat)
{
  /*
    Output variables:

    int nbr_add[numNeighbors]     ! address of each of the closest points
    double nbr_dist[numNeighbors] ! distance to each of the closest points

    Input variables:

    double plat,         ! latitude  of the search point
    double plon,         ! longitude of the search point
  */

  size_t numNeighbors = knnWeights.maxNeighbors();
  size_t *restrict nbr_add = &knnWeights.m_addr[0];
  double *restrict nbr_dist = &knnWeights.m_dist[0];

  size_t src_add[MAX_SEARCH_CELLS];
  std::vector<size_t> src_add_tmp;
  size_t *psrc_add = src_add;
  size_t num_add = 0;
  double *restrict src_center_lon = gs->reg2d_center_lon;
  double *restrict src_center_lat = gs->reg2d_center_lat;

  long nx = gs->dims[0];
  long ny = gs->dims[1];

  size_t nxm = gs->is_cyclic ? nx + 1 : nx;

  if (plon < src_center_lon[0]) plon += PI2;
  if (plon > src_center_lon[nxm - 1]) plon -= PI2;

  size_t ii, jj;
  int lfound = rect_grid_search(&ii, &jj, plon, plat, nxm, ny, src_center_lon, src_center_lat);
  if (lfound)
    {
      if (gs->is_cyclic && ii == (nxm - 1)) ii = 0;

      long k;
      for (k = 3; k < 10000; k += 2)
        if (numNeighbors <= (k - 2) * (k - 2)) break;

      if ((k * k) > MAX_SEARCH_CELLS)
        {
          src_add_tmp.resize(k * k);
          psrc_add = &src_add_tmp[0];
        }

      k /= 2;

      long j0 = jj - k;
      long jn = jj + k;
      long i0 = ii - k;
      long in = ii + k;
      if (j0 < 0) j0 = 0;
      if (jn >= ny) jn = ny - 1;
      if ((in - i0) > nx)
        {
          i0 = 0;
          in = nx - 1;
        }

      for (long j = j0; j <= jn; ++j)
        for (long i = i0; i <= in; ++i)
          {
            long ix = i;
            if (gs->is_cyclic)
              {
                if (ix < 0) ix += nx;
                if (ix >= nx) ix -= nx;
              }

            if (ix >= 0 && ix < nx && j >= 0 && j < ny) psrc_add[num_add++] = j * nx + ix;
          }
    }

  // Initialize distance and address arrays
  knnWeights.init_addr();
  knnWeights.init_dist();

  if (lfound)
    {
      double *restrict coslon = gs->coslon;
      double *restrict sinlon = gs->sinlon;
      double *restrict coslat = gs->coslat;
      double *restrict sinlat = gs->sinlat;

      double xyz[3];
      double query_pt[3];
      LLtoXYZ(plon, plat, query_pt);
      double search_radius = SQR(gs->search_radius);

      for (size_t na = 0; na < num_add; ++na)
        {
          size_t nadd = psrc_add[na];
          size_t iy = nadd / nx;
          size_t ix = nadd - iy * nx;

          xyz[0] = coslat[iy] * coslon[ix];
          xyz[1] = coslat[iy] * sinlon[ix];
          xyz[2] = sinlat[iy];
          // Find distance to this point
          double dist = (float) distance(query_pt, xyz);
          if (dist <= search_radius)
            {
              // Store the address and distance if this is one of the smallest so far
              knnWeights.store_distance(nadd, sqrt(dist), numNeighbors);
            }
        }

      knnWeights.check_distance();
    }
  else if (gs->extrapolate)
    {
      int search_result;

      if (numNeighbors < 4)
        {
          size_t nbr_add4[4];
          double nbr_dist4[4];
          for (size_t n = 0; n < numNeighbors; ++n) nbr_add4[n] = SIZE_MAX;
          search_result = grid_search_reg2d_nn(nx, ny, nbr_add4, nbr_dist4, plat, plon, src_center_lat, src_center_lon);
          if (search_result < 0)
            {
              for (size_t n = 0; n < numNeighbors; ++n) nbr_add[n] = nbr_add4[n];
              for (size_t n = 0; n < numNeighbors; ++n) nbr_dist[n] = nbr_dist4[n];
            }
        }
      else
        {
          search_result = grid_search_reg2d_nn(nx, ny, nbr_add, nbr_dist, plat, plon, src_center_lat, src_center_lon);
        }

      if (search_result >= 0)
        for (size_t n = 0; n < numNeighbors; ++n) nbr_add[n] = SIZE_MAX;
    }
}  // grid_search_nbr_reg2d

void
grid_search_nbr(GridSearch *gs, knnWeightsType &knnWeights, double plon, double plat)
{
  /*
    Output variables:

    int nbr_add[numNeighbors]     ! address of each of the closest points
    double nbr_dist[numNeighbors] ! distance to each of the closest points

    Input variables:

    double plat,         ! latitude  of the search point
    double plon,         ! longitude of the search point
  */

  size_t numNeighbors = knnWeights.maxNeighbors();

  // Initialize distance and address arrays
  knnWeights.init_addr();
  knnWeights.init_dist();

  size_t ndist = numNeighbors;
  // check some more points if distance is the same use the smaller index (nadd)
  if (ndist > 8)
    ndist += 8;
  else
    ndist *= 2;
  if (ndist > gs->n) ndist = gs->n;

  if (knnWeights.m_tmpaddr.size() == 0) knnWeights.m_tmpaddr.resize(ndist);
  if (knnWeights.m_tmpdist.size() == 0) knnWeights.m_tmpdist.resize(ndist);
  size_t *adds = &knnWeights.m_tmpaddr[0];
  double *dist = &knnWeights.m_tmpdist[0];

  const double range0 = SQR(gs->search_radius);
  double range = range0;

  size_t nadds = 0;

  if (numNeighbors == 1)
    {
      size_t add = gridsearch_nearest(gs, plon, plat, &range);
      if (add != GS_NOT_FOUND)
        {
          // if ( range < range0 )
          {
            dist[nadds] = sqrt(range);
            adds[nadds] = add;
            nadds++;
          }
        }
    }
  else
    {
      nadds = gridsearch_qnearest(gs, plon, plat, &range, ndist, adds, dist);
      for (size_t i = 0; i < nadds; ++i) dist[i] = sqrt(dist[i]);
    }

  ndist = nadds;
  if (ndist < numNeighbors) numNeighbors = ndist;

  for (size_t i = 0; i < ndist; ++i) knnWeights.store_distance(adds[i], dist[i], numNeighbors);

  knnWeights.check_distance();
}  // grid_search_nbr

//  This routine computes the inverse-distance weights for a nearest-neighbor
//  interpolation.

void
remap_distwgt_weights(size_t numNeighbors, RemapSearch &rsearch, RemapGrid *src_grid, RemapGrid *tgt_grid, RemapVars &rv)
{
  int remap_grid_type = src_grid->remap_grid_type;

  if (cdoVerbose) cdoPrint("Called %s()", __func__);

  progressInit();

  // Compute mappings from source to target grid

  size_t src_grid_size = src_grid->size;
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

  bool xIsCyclic = src_grid->is_cyclic;
  size_t *dims = src_grid->dims;
  GridSearch *gs = NULL;
  if (remap_grid_type == REMAP_GRID_TYPE_REG2D)
    gs = gridsearch_create_reg2d(xIsCyclic, dims, src_grid->reg2d_center_lon, src_grid->reg2d_center_lat);
  else
    gs = gridsearch_create(xIsCyclic, dims, src_grid_size, src_grid->cell_center_lon, src_grid->cell_center_lat);

  if (src_grid->lextrapolate) gridsearch_extrapolate(gs);

#ifdef _OPENMP
  if (cdoVerbose) printf("gridsearch created: %.2f seconds\n", omp_get_wtime() - start);
  if (cdoVerbose) start = omp_get_wtime();
#endif

  // Loop over destination grid

  double findex = 0;

#ifdef HAVE_OPENMP4
#pragma omp parallel for default(none) reduction(+ : findex) shared(gs, weightlinks, numNeighbors, remap_grid_type, \
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
      if (remap_grid_type == REMAP_GRID_TYPE_REG2D)
        grid_search_nbr_reg2d(gs, knnWeights[ompthID], plon, plat);
      else
        grid_search_nbr(gs, knnWeights[ompthID], plon, plat);

      // Compute weights based on inverse distance if mask is false, eliminate those points
      size_t nadds = knnWeights[ompthID].compute_weights(src_grid->mask);

      for (size_t n = 0; n < nadds; ++n)
        if (knnWeights[ompthID].m_mask[n]) tgt_grid->cell_frac[tgt_cell_add] = ONE;

      // Store the link
      store_weightlinks(0, nadds, &knnWeights[ompthID].m_addr[0], &knnWeights[ompthID].m_dist[0], tgt_cell_add,
                        &weightlinks[0]);
    }

  progressStatus(0, 1, 1);

  if (gs) gridsearch_delete(gs);

  weightlinks2remaplinks(0, tgt_grid_size, &weightlinks[0], rv);

#ifdef _OPENMP
  if (cdoVerbose) printf("gridsearch nearest: %.2f seconds\n", omp_get_wtime() - start);
#endif
}  // remap_distwgt_weights

void
remap_distwgt(size_t numNeighbors, RemapSearch &rsearch, RemapGrid *src_grid, RemapGrid *tgt_grid, const double *restrict src_array,
              double *restrict tgt_array, double missval)
{
  int src_remap_grid_type = src_grid->remap_grid_type;

  if (cdoVerbose) cdoPrint("Called %s()", __func__);

  progressInit();

  // Compute mappings from source to target grid

  size_t src_grid_size = src_grid->size;
  size_t tgt_grid_size = tgt_grid->size;

  std::vector<knnWeightsType> knnWeights;
  for (int i = 0; i < Threading::ompNumThreads; ++i) knnWeights.push_back(knnWeightsType(numNeighbors));

#ifdef _OPENMP
  double start = cdoVerbose ? omp_get_wtime() : 0;
#endif

  bool xIsCyclic = src_grid->is_cyclic;
  size_t *dims = src_grid->dims;
  GridSearch *gs = NULL;
  if (src_remap_grid_type == REMAP_GRID_TYPE_REG2D)
    gs = gridsearch_create_reg2d(xIsCyclic, dims, src_grid->reg2d_center_lon, src_grid->reg2d_center_lat);
  else
    gs = gridsearch_create(xIsCyclic, dims, src_grid_size, src_grid->cell_center_lon, src_grid->cell_center_lat);

  if (src_grid->lextrapolate) gridsearch_extrapolate(gs);

#ifdef _OPENMP
  if (cdoVerbose) printf("gridsearch created: %.2f seconds\n", omp_get_wtime() - start);
  if (cdoVerbose) start = omp_get_wtime();
#endif

  // Loop over destination grid

  double findex = 0;

#ifdef HAVE_OPENMP4
#pragma omp parallel for default(none)                                                      \
    reduction(+ : findex) shared(gs, numNeighbors, src_remap_grid_type, src_grid, tgt_grid, \
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
      if (src_remap_grid_type == REMAP_GRID_TYPE_REG2D)
        grid_search_nbr_reg2d(gs, knnWeights[ompthID], plon, plat);
      else
        grid_search_nbr(gs, knnWeights[ompthID], plon, plat);

      // Compute weights based on inverse distance if mask is false, eliminate those points
      size_t nadds = knnWeights[ompthID].compute_weights(src_grid->mask);
      if (nadds) tgt_array[tgt_cell_add] = knnWeights[ompthID].array_weights_sum(src_array);
    }

  progressStatus(0, 1, 1);

  if (gs) gridsearch_delete(gs);

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
  // int src_remap_grid_type = src_grid->remap_grid_type;

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
      grid_search_nbr(gs, knnWeights[ompthID], plon, plat);

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

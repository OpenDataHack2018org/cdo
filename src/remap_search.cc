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

#include "cdo_int.h"
#include "remap.h"
#include "grid_search.h"

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

// This routine finds the closest num_neighbor points to a search point and computes a distance to each of the neighbors
#define MAX_SEARCH_CELLS 25
static void
grid_search_nbr_reg2d(GridSearch *gs, double plon, double plat, knnWeightsType &knnWeights)
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
grid_search_nbr(GridSearch *gs, double plon, double plat, knnWeightsType &knnWeights)
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


void remapSearchPoints(RemapSearch &rsearch, double plon, double plat, knnWeightsType &knnWeights)
{
  int remap_grid_type = rsearch.srcGrid->remap_grid_type;

  if (remap_grid_type == REMAP_GRID_TYPE_REG2D)
    grid_search_nbr_reg2d(rsearch.gs, plon, plat, knnWeights);
  else
    grid_search_nbr(rsearch.gs, plon, plat, knnWeights);
}


int remapSearchSquare(RemapSearch &rsearch, double plon, double plat, size_t src_add[4], double src_lats[4], double src_lons[4])
{
  RemapGrid *src_grid = rsearch.srcGrid;
  int remap_grid_type = src_grid->remap_grid_type;

  int search_result;
  if (remap_grid_type == REMAP_GRID_TYPE_REG2D)
    search_result = grid_search_reg2d(src_grid, src_add, src_lats, src_lons, plat, plon);
  else
    search_result
      = grid_search(src_grid, src_add, src_lats, src_lons, plat, plon, rsearch.srcBins);

  return search_result;
}

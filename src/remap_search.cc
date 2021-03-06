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

// This routine finds the closest numNeighbor points to a search point and computes a distance to each of the neighbors
#define MAX_SEARCH_CELLS 25
static void
gridSearchPointReg2d(GridPointSearch *gps, double plon, double plat, knnWeightsType &knnWeights)
{
  /*
    Input variables:

      double plat                   ! latitude  of the search point
      double plon                   ! longitude of the search point

    Output variables:

      int knnWeights.m_addr[numNeighbors]    ! address of each of the closest points
      double knnWeights.m_dist[numNeighbors] ! distance to each of the closest points
  */
  size_t numNeighbors = knnWeights.maxNeighbors();
  size_t *restrict nbr_add = &knnWeights.m_addr[0];
  double *restrict nbr_dist = &knnWeights.m_dist[0];

  size_t src_add[MAX_SEARCH_CELLS];
  std::vector<size_t> src_add_tmp;
  size_t *psrc_add = src_add;
  size_t num_add = 0;
  const double *restrict src_center_lon = gps->reg2d_center_lon;
  const double *restrict src_center_lat = gps->reg2d_center_lat;

  long nx = gps->dims[0];
  long ny = gps->dims[1];

  size_t nxm = gps->is_cyclic ? nx + 1 : nx;

  if (plon < src_center_lon[0]) plon += PI2;
  if (plon > src_center_lon[nxm - 1]) plon -= PI2;

  size_t ii, jj;
  int lfound = rect_grid_search(&ii, &jj, plon, plat, nxm, ny, src_center_lon, src_center_lat);
  if (lfound)
    {
      if (gps->is_cyclic && ii == (nxm - 1)) ii = 0;

      long k;
      for (k = 3; k < 10000; k += 2)
        if (numNeighbors <= (size_t)(k - 2) * (k - 2)) break;

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
            if (gps->is_cyclic)
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
      const double *restrict coslon = gps->coslon;
      const double *restrict sinlon = gps->sinlon;
      const double *restrict coslat = gps->coslat;
      const double *restrict sinlat = gps->sinlat;

      double xyz[3];
      double query_pt[3];
      cdoLLtoXYZ(plon, plat, query_pt);
      double sqrSearchRadius = SQR(gps->searchRadius);

      for (size_t na = 0; na < num_add; ++na)
        {
          size_t nadd = psrc_add[na];
          size_t iy = nadd / nx;
          size_t ix = nadd - iy * nx;

          xyz[0] = coslat[iy] * coslon[ix];
          xyz[1] = coslat[iy] * sinlon[ix];
          xyz[2] = sinlat[iy];
          // Find distance to this point
          double sqrDist = (float) squareDistance(query_pt, xyz);
          if (sqrDist <= sqrSearchRadius)
            {
              // Store the address and distance if this is one of the smallest so far
              knnWeights.store_distance(nadd, sqrt(sqrDist), numNeighbors);
            }
        }

      knnWeights.check_distance();
    }
  else if (gps->extrapolate)
    {
      int search_result;

      if (numNeighbors < 4)
        {
          size_t nbr_add4[4];
          double nbr_dist4[4];
          for (size_t n = 0; n < numNeighbors; ++n) nbr_add4[n] = SIZE_MAX;
          search_result = gridSearchSquareReg2dNN(nx, ny, nbr_add4, nbr_dist4, plat, plon, src_center_lat, src_center_lon);
          if (search_result < 0)
            {
              for (size_t n = 0; n < numNeighbors; ++n) nbr_add[n] = nbr_add4[n];
              for (size_t n = 0; n < numNeighbors; ++n) nbr_dist[n] = nbr_dist4[n];
            }
        }
      else
        {
          search_result = gridSearchSquareReg2dNN(nx, ny, nbr_add, nbr_dist, plat, plon, src_center_lat, src_center_lon);
        }

      if (search_result >= 0)
        for (size_t n = 0; n < numNeighbors; ++n) nbr_add[n] = SIZE_MAX;
    }
}

void
gridSearchPoint(GridPointSearch *gps, double plon, double plat, knnWeightsType &knnWeights)
{
  /*
    Input variables:

      plat: latitude  of the search point
      plon: longitude of the search point

    Output variables:

      knnWeights.m_addr[numNeighbors]: address of each of the closest points
      knnWeights.m_dist[numNeighbors]: distance to each of the closest points
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
  if (ndist > gps->n) ndist = gps->n;

  if (knnWeights.m_tmpaddr.size() == 0) knnWeights.m_tmpaddr.resize(ndist);
  if (knnWeights.m_tmpdist.size() == 0) knnWeights.m_tmpdist.resize(ndist);
  size_t *adds = &knnWeights.m_tmpaddr[0];
  double *dist = &knnWeights.m_tmpdist[0];

  size_t nadds = 0;
  if (numNeighbors == 1)
    nadds = gridPointSearchNearest(gps, plon, plat, adds, dist);
  else
    nadds = gridPointSearchQnearest(gps, plon, plat, ndist, adds, dist);

  ndist = nadds;
  if (ndist < numNeighbors) numNeighbors = ndist;

  for (size_t i = 0; i < ndist; ++i) knnWeights.store_distance(adds[i], dist[i], numNeighbors);

  knnWeights.check_distance();
}

void
remapSearchPoints(RemapSearch &rsearch, double plon, double plat, knnWeightsType &knnWeights)
{
  if (rsearch.srcGrid->remap_grid_type == REMAP_GRID_TYPE_REG2D)
    gridSearchPointReg2d(rsearch.gps, plon, plat, knnWeights);
  else
    gridSearchPoint(rsearch.gps, plon, plat, knnWeights);
}

static int
gridSearchSquareCurv2d(GridPointSearch *gps, RemapGrid *rgrid, size_t *restrict src_add, double *restrict src_lats,
                       double *restrict src_lons, double plat, double plon)
{
  /*
    Input variables:

      plat: latitude  of the search point
      plon: longitude of the search point

    Output variables:

      src_add[4]:   address of each corner point enclosing P
      src_lats[4]:  latitudes  of the four corner points
      src_lons[4]:  longitudes of the four corner points
  */
  int search_result = 0;

  for (unsigned n = 0; n < 4; ++n) src_add[n] = 0;

  size_t nx = rgrid->dims[0];
  size_t ny = rgrid->dims[1];

  double dist;
  size_t addr;
  size_t nadds = gridPointSearchNearest(gps, plon, plat, &addr, &dist);
  if (nadds > 0)
    {
      for (unsigned k = 0; k < 4; ++k)
        {
          // Determine neighbor addresses
          size_t j = addr / nx;
          size_t i = addr - j * nx;
          if (k == 0 || k == 2) i = (i > 0) ? i - 1 : rgrid->is_cyclic ? nx - 1 : 0;
          if (k == 0 || k == 1) j = (j > 0) ? j - 1 : 0;
          if (pointInQuad(rgrid->is_cyclic, nx, ny, i, j, src_add, src_lons, src_lats, plon, plat, rgrid->cell_center_lon,
                          rgrid->cell_center_lat))
            {
              search_result = 1;
              return search_result;
            }
        }
    }

  /*
    If no cell found, point is likely either in a box that straddles either pole or is outside the grid.
    Fall back to a distance-weighted average of the four closest points. Go ahead and compute weights here,
    but store in src_lats and return -add to prevent the parent routine from computing bilinear weights.
  */
  if (!rgrid->lextrapolate) return search_result;

  size_t ndist = 4;
  nadds = gridPointSearchQnearest(gps, plon, plat, ndist, src_add, src_lats);
  if (nadds == 4)
    {
      for (unsigned n = 0; n < 4; ++n) src_lats[n] = 1.0 / (src_lats[n] + TINY);
      double distance = 0.0;
      for (unsigned n = 0; n < 4; ++n) distance += src_lats[n];
      for (unsigned n = 0; n < 4; ++n) src_lats[n] /= distance;
      search_result = -1;
    }

  return search_result;
}

int
remapSearchSquare(RemapSearch &rsearch, double plon, double plat, size_t *src_add, double *src_lats, double *src_lons)
{
  int searchResult;
  if (rsearch.srcGrid->remap_grid_type == REMAP_GRID_TYPE_REG2D)
    searchResult = gridSearchSquareReg2d(rsearch.srcGrid, src_add, src_lats, src_lons, plat, plon);
  else if (rsearch.gps)
    searchResult = gridSearchSquareCurv2d(rsearch.gps, rsearch.srcGrid, src_add, src_lats, src_lons, plat, plon);
  else
    searchResult = gridSearchSquareCurv2dScrip(rsearch.srcGrid, src_add, src_lats, src_lons, plat, plon, rsearch.srcBins);

  return searchResult;
}

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
#include "grid.h"
#include "remap.h"

static void
calcBinAddr(GridSearchBins &searchBins)
{
  const size_t ncells = searchBins.ncells;
  const size_t nbins = searchBins.nbins;
  size_t *restrict bin_addr = &searchBins.bin_addr[0];
  const float *restrict bin_lats = &searchBins.bin_lats[0];
  const float *restrict cell_bound_box = &searchBins.cell_bound_box[0];

  for (size_t n = 0; n < nbins; ++n)
    {
      bin_addr[n * 2] = ncells;
      bin_addr[n * 2 + 1] = 0;
    }

  for (size_t nele = 0; nele < ncells; ++nele)
    {
      const size_t nele4 = nele << 2;
      const float cell_bound_box_lat1 = cell_bound_box[nele4];
      const float cell_bound_box_lat2 = cell_bound_box[nele4 + 1];
      for (size_t n = 0; n < nbins; ++n)
        {
          const size_t n2 = n << 1;
          if (cell_bound_box_lat1 <= bin_lats[n2 + 1] && cell_bound_box_lat2 >= bin_lats[n2])
            {
              bin_addr[n2] = MIN(nele, bin_addr[n2]);
              bin_addr[n2 + 1] = MAX(nele, bin_addr[n2 + 1]);
            }
        }
    }
}

void
calcLatBins(GridSearchBins &searchBins)
{
  const size_t nbins = searchBins.nbins;
  const double dlat = PI / nbins;  // lat interval for search bins

  if (cdoVerbose) cdoPrint("Using %zu latitude bins to restrict search.", nbins);

  if (nbins > 0)
    {
      searchBins.bin_lats.resize(2 * nbins);
      for (size_t n = 0; n < nbins; ++n)
        {
          searchBins.bin_lats[n * 2] = (n) *dlat - PIH;
          searchBins.bin_lats[n * 2 + 1] = (n + 1) * dlat - PIH;
        }

      searchBins.bin_addr.resize(2 * nbins);
      calcBinAddr(searchBins);
    }
}

size_t
get_srch_cells(size_t tgt_cell_addr, GridSearchBins &tgtBins, GridSearchBins &srcBins, float *tgt_cell_bound_box,
               size_t *srch_add)
{
  size_t nbins = srcBins.nbins;
  size_t src_grid_size = srcBins.ncells;
  const size_t *restrict bin_addr1 = &tgtBins.bin_addr[0];
  const size_t *restrict bin_addr2 = &srcBins.bin_addr[0];
  const float *restrict src_cell_bound_box = &srcBins.cell_bound_box[0];

  size_t src_cell_addm4;

  // Restrict searches first using search bins

  size_t min_add = src_grid_size - 1;
  size_t max_add = 0;

  for (size_t n = 0; n < nbins; ++n)
    {
      size_t n2 = n << 1;
      if (tgt_cell_addr >= bin_addr1[n2] && tgt_cell_addr <= bin_addr1[n2 + 1])
        {
          if (bin_addr2[n2] < min_add) min_add = bin_addr2[n2];
          if (bin_addr2[n2 + 1] > max_add) max_add = bin_addr2[n2 + 1];
        }
    }

  // Further restrict searches using bounding boxes

  float bound_box_lat1 = tgt_cell_bound_box[0];
  float bound_box_lat2 = tgt_cell_bound_box[1];
  float bound_box_lon1 = tgt_cell_bound_box[2];
  float bound_box_lon2 = tgt_cell_bound_box[3];

  size_t num_srch_cells = 0;
  for (size_t src_cell_add = min_add; src_cell_add <= max_add; ++src_cell_add)
    {
      src_cell_addm4 = src_cell_add << 2;
      if ((src_cell_bound_box[src_cell_addm4 + 2] <= bound_box_lon2)
          && (src_cell_bound_box[src_cell_addm4 + 3] >= bound_box_lon1))
        {
          if ((src_cell_bound_box[src_cell_addm4] <= bound_box_lat2)
              && (src_cell_bound_box[src_cell_addm4 + 1] >= bound_box_lat1))
            {
              srch_add[num_srch_cells] = src_cell_add;
              num_srch_cells++;
            }
        }
    }

  if (bound_box_lon1 < 0.0f || bound_box_lon2 > PI2_f)
    {
      if (bound_box_lon1 < 0.0f)
        {
          bound_box_lon1 += PI2_f;
          bound_box_lon2 += PI2_f;
        }
      else
        {
          bound_box_lon1 -= PI2_f;
          bound_box_lon2 -= PI2_f;
        }

      for (size_t src_cell_add = min_add; src_cell_add <= max_add; ++src_cell_add)
        {
          src_cell_addm4 = src_cell_add << 2;
          if ((src_cell_bound_box[src_cell_addm4 + 2] <= bound_box_lon2)
              && (src_cell_bound_box[src_cell_addm4 + 3] >= bound_box_lon1))
            {
              if ((src_cell_bound_box[src_cell_addm4] <= bound_box_lat2)
                  && (src_cell_bound_box[src_cell_addm4 + 1] >= bound_box_lat1))
                {
                  size_t ii;
                  for (ii = 0; ii < num_srch_cells; ++ii)
                    if (srch_add[ii] == src_cell_add) break;

                  if (ii == num_srch_cells)
                    {
                      srch_add[num_srch_cells] = src_cell_add;
                      num_srch_cells++;
                    }
                }
            }
        }
    }

  return num_srch_cells;
}

static int
gridSearchSquareCurv2dNNScrip(size_t min_add, size_t max_add, size_t *restrict nbr_add, double *restrict nbr_dist, double plat,
                              double plon, const double *restrict src_center_lat, const double *restrict src_center_lon)
{
  int search_result = 0;
  double distance;
  double coslat_dst = cos(plat);
  double sinlat_dst = sin(plat);
  double coslon_dst = cos(plon);
  double sinlon_dst = sin(plon);

  double dist_min = DBL_MAX;
  for (unsigned n = 0; n < 4; ++n) nbr_dist[n] = DBL_MAX;

  for (size_t srch_add = min_add; srch_add <= max_add; ++srch_add)
    {
      distance = acos(coslat_dst * cos(src_center_lat[srch_add])
                      * (coslon_dst * cos(src_center_lon[srch_add]) + sinlon_dst * sin(src_center_lon[srch_add]))
                      + sinlat_dst * sin(src_center_lat[srch_add]));

      if (distance < dist_min)
        {
          for (unsigned n = 0; n < 4; ++n)
            {
              if (distance < nbr_dist[n])
                {
                  for (unsigned i = 3; i > n; --i)
                    {
                      nbr_add[i] = nbr_add[i - 1];
                      nbr_dist[i] = nbr_dist[i - 1];
                    }
                  search_result = -1;
                  nbr_add[n] = srch_add;
                  nbr_dist[n] = distance;
                  dist_min = nbr_dist[3];
                  break;
                }
            }
        }
    }

  for (unsigned n = 0; n < 4; ++n) nbr_dist[n] = 1.0 / (nbr_dist[n] + TINY);
  distance = 0.0;
  for (unsigned n = 0; n < 4; ++n) distance += nbr_dist[n];
  for (unsigned n = 0; n < 4; ++n) nbr_dist[n] /= distance;

  return search_result;
}

static unsigned
quadCrossProducts(double plon, double plat, double *restrict lons, double *restrict lats)
{
  unsigned n;
  int scross[4], scross_last = 0;
  // Vectors for cross-product check
  double vec1_lat, vec1_lon;
  double vec2_lat, vec2_lon;

  // For consistency, we must make sure all lons are in same 2pi interval
  vec1_lon = lons[0] - plon;
  if (vec1_lon > PI)
    lons[0] -= PI2;
  else if (vec1_lon < -PI)
    lons[0] += PI2;

  for (n = 1; n < 4; ++n)
    {
      vec1_lon = lons[n] - lons[0];
      if (vec1_lon > PI)
        lons[n] -= PI2;
      else if (vec1_lon < -PI)
        lons[n] += PI2;
    }

  /* corner_loop */
  for (n = 0; n < 4; ++n)
    {
      unsigned next_n = (n + 1) % 4;
      /*
        Here we take the cross product of the vector making up each box side
        with the vector formed by the vertex and search point.
        If all the cross products are positive, the point is contained in the box.
      */
      vec1_lat = lats[next_n] - lats[n];
      vec1_lon = lons[next_n] - lons[n];
      vec2_lat = plat - lats[n];
      vec2_lon = plon - lons[n];

      // Check for 0,2pi crossings
      if (vec1_lon > 3.0 * PIH)
        vec1_lon -= PI2;
      else if (vec1_lon < -3.0 * PIH)
        vec1_lon += PI2;

      if (vec2_lon > 3.0 * PIH)
        vec2_lon -= PI2;
      else if (vec2_lon < -3.0 * PIH)
        vec2_lon += PI2;

      double cross_product = vec1_lon * vec2_lat - vec2_lon * vec1_lat;

      /* If cross product is less than ZERO, this cell doesn't work    */
      /* 2008-10-16 Uwe Schulzweida: bug fix for cross_product eq zero */
      scross[n] = cross_product < 0 ? -1 : cross_product > 0 ? 1 : 0;
      if (n == 0) scross_last = scross[n];
      if ((scross[n] < 0 && scross_last > 0) || (scross[n] > 0 && scross_last < 0)) break;
      scross_last = scross[n];
    }

  if (n >= 4)
    {
      n = 0;
      if (scross[0] >= 0 && scross[1] >= 0 && scross[2] >= 0 && scross[3] >= 0)
        n = 4;
      else if (scross[0] <= 0 && scross[1] <= 0 && scross[2] <= 0 && scross[3] <= 0)
        n = 4;
    }

  return n;
}

bool
pointInQuad(bool isCyclic, size_t nx, size_t ny, size_t i, size_t j, size_t adds[4], double lons[4], double lats[4],
            double plon, double plat, const double *restrict centerLon, const double *restrict centerLat)
{
  bool search_result = false;
  size_t ip1 = (i < (nx - 1)) ? i + 1 : isCyclic ? 0 : i;
  size_t jp1 = (j < (ny - 1)) ? j + 1 : j;

  if (i == ip1 || j == jp1) return search_result;

  size_t idx[4];
  idx[0] = j * nx + i;
  idx[1] = j * nx + ip1;    // east
  idx[2] = jp1 * nx + ip1;  // north-east
  idx[3] = jp1 * nx + i;    // north

  for (unsigned j = 0; j < 4; ++j) lons[j] = centerLon[idx[j]];
  for (unsigned j = 0; j < 4; ++j) lats[j] = centerLat[idx[j]];

  unsigned n = quadCrossProducts(plon, plat, lons, lats);

  // If cross products all same sign, we found the location
  if (n >= 4)
    {
      for (unsigned j = 0; j < 4; ++j) adds[j] = idx[j];
      search_result = true;
    }

  return search_result;
}

int
gridSearchSquareCurv2dScrip(RemapGrid *src_grid, size_t *restrict src_add, double *restrict src_lats, double *restrict src_lons,
                            double plat, double plon, GridSearchBins &srcBins)
{
  /*
    Output variables:

    int    src_add[4]              ! address of each corner point enclosing P
    double src_lats[4]             ! latitudes  of the four corner points
    double src_lons[4]             ! longitudes of the four corner points

    Input variables:

    double plat                    ! latitude  of the search point
    double plon                    ! longitude of the search point

    int src_grid_dims[2]           ! size of each src grid dimension

    double src_center_lat[]        ! latitude  of each src grid center
    double src_center_lon[]        ! longitude of each src grid center

    float src_grid_bound_box[][4] ! bound box for source grid

    int src_bin_addr[][2]           ! latitude bins for restricting
  */
  int search_result = 0;

  size_t nbins = srcBins.nbins;
  const size_t *restrict src_bin_addr = &srcBins.bin_addr[0];
  const float *restrict bin_lats = &srcBins.bin_lats[0];
  const float *restrict src_grid_bound_box = &srcBins.cell_bound_box[0];

  size_t n2, srch_add, srch_add4;

  float rlat = plat;
  float rlon = plon;

  // restrict search first using bins

  for (unsigned n = 0; n < 4; ++n) src_add[n] = 0;

  // addresses for restricting search
  size_t min_add = src_grid->size - 1;
  size_t max_add = 0;

  for (size_t n = 0; n < nbins; ++n)
    {
      n2 = n << 1;
      if (rlat >= bin_lats[n2] && rlat <= bin_lats[n2 + 1])
        {
          if (src_bin_addr[n2] < min_add) min_add = src_bin_addr[n2];
          if (src_bin_addr[n2 + 1] > max_add) max_add = src_bin_addr[n2 + 1];
        }
    }

  /* Now perform a more detailed search */

  size_t nx = src_grid->dims[0];
  size_t ny = src_grid->dims[1];

  for (srch_add = min_add; srch_add <= max_add; ++srch_add)
    {
      srch_add4 = srch_add << 2;
      /* First check bounding box */
      if (rlon >= src_grid_bound_box[srch_add4 + 2] && rlon <= src_grid_bound_box[srch_add4 + 3]
          && rlat >= src_grid_bound_box[srch_add4] && rlat <= src_grid_bound_box[srch_add4 + 1])
        {
          /* We are within bounding box so get really serious */

          /* Determine neighbor addresses */
          size_t j = srch_add / nx;
          size_t i = srch_add - j * nx;

          if (pointInQuad(src_grid->is_cyclic, nx, ny, i, j, src_add, src_lons, src_lats, plon, plat,
                          src_grid->cell_center_lon, src_grid->cell_center_lat))
            {
              search_result = 1;
              return search_result;
            }

          /* Otherwise move on to next cell */

        } /* Bounding box check */
    }

  /*
    If no cell found, point is likely either in a box that straddles either pole or is outside the grid.
    Fall back to a distance-weighted average of the four closest points. Go ahead and compute weights here,
    but store in src_lats and return -add to prevent the parent routine from computing bilinear weights.
  */
  if (!src_grid->lextrapolate) return search_result;

  /*
    printf("Could not find location for %g %g\n", plat*RAD2DEG, plon*RAD2DEG);
    printf("Using nearest-neighbor average for this point\n");
  */
  search_result = gridSearchSquareCurv2dNNScrip(min_add, max_add, src_add, src_lats, plat, plon,
                                                src_grid->cell_center_lat, src_grid->cell_center_lon);

  return search_result;
}

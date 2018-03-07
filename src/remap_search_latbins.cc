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
calc_bin_addr(size_t gridsize, size_t nbins, const float *restrict bin_lats, const float *restrict cell_bound_box,
              size_t *restrict bin_addr)
{
  for (size_t n = 0; n < nbins; ++n)
    {
      bin_addr[n*2] = gridsize;
      bin_addr[n*2 + 1] = 0;
    }

  for (size_t nele = 0; nele < gridsize; ++nele)
    {
      size_t n2;
      size_t nele4 = nele << 2;
      float cell_bound_box_lat1 = cell_bound_box[nele4];
      float cell_bound_box_lat2 = cell_bound_box[nele4 + 1];
      for (size_t n = 0; n < nbins; ++n)
        {
          n2 = n << 1;
          if (cell_bound_box_lat1 <= bin_lats[n2 + 1] && cell_bound_box_lat2 >= bin_lats[n2])
            {
              bin_addr[n2] = MIN(nele, bin_addr[n2]);
              bin_addr[n2 + 1] = MAX(nele, bin_addr[n2 + 1]);
            }
        }
    }
}

void
calc_lat_bins(remapGridType *src_grid, remapGridType *tgt_grid, RemapType mapType)
{
  size_t n2;
  size_t nbins = src_grid->num_srch_bins;
  double dlat = PI / nbins;  // lat/lon intervals for search bins

  if (cdoVerbose) cdoPrint("Using %zu latitude bins to restrict search.", nbins);

  if (nbins > 0)
    {
      float *bin_lats = src_grid->bin_lats = (float *) Realloc(src_grid->bin_lats, 2 * nbins * sizeof(float));

      for (size_t n = 0; n < nbins; ++n)
        {
          n2 = n << 1;
          bin_lats[n2] = (n) *dlat - PIH;
          bin_lats[n2 + 1] = (n + 1) * dlat - PIH;
        }

      src_grid->bin_addr = (size_t *) Malloc(2 * nbins * sizeof(size_t));
      calc_bin_addr(src_grid->size, nbins, bin_lats, src_grid->cell_bound_box, src_grid->bin_addr);

      if (mapType == RemapType::CONSERV || mapType == RemapType::CONSERV_YAC)
        {
          tgt_grid->bin_addr = (size_t *) Malloc(2 * nbins * sizeof(size_t));
          calc_bin_addr(tgt_grid->size, nbins, bin_lats, tgt_grid->cell_bound_box, tgt_grid->bin_addr);

          Free(src_grid->bin_lats);
          src_grid->bin_lats = NULL;
        }
    }

  if (mapType == RemapType::CONSERV_YAC)
    {
      Free(tgt_grid->cell_bound_box);
      tgt_grid->cell_bound_box = NULL;
    }

  if (mapType == RemapType::DISTWGT)
    {
      Free(src_grid->cell_bound_box);
      src_grid->cell_bound_box = NULL;
    }
}

size_t
get_srch_cells(size_t tgt_cell_add, size_t nbins, size_t *bin_addr1, size_t *bin_addr2, float *tgt_cell_bound_box,
               float *src_cell_bound_box, size_t src_grid_size, size_t *srch_add)
{
  size_t n2;
  size_t src_cell_addm4;

  /* Restrict searches first using search bins */

  size_t min_add = src_grid_size - 1;
  size_t max_add = 0;

  for (size_t n = 0; n < nbins; ++n)
    {
      n2 = n << 1;
      if (tgt_cell_add >= bin_addr1[n2] && tgt_cell_add <= bin_addr1[n2 + 1])
        {
          if (bin_addr2[n2] < min_add) min_add = bin_addr2[n2];
          if (bin_addr2[n2 + 1] > max_add) max_add = bin_addr2[n2 + 1];
        }
    }

  /* Further restrict searches using bounding boxes */

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
grid_search_nn(size_t min_add, size_t max_add, size_t *restrict nbr_add, double *restrict nbr_dist, double plat,
               double plon, const double *restrict src_center_lat, const double *restrict src_center_lon)
{
  int search_result = 0;
  double distance; /* For computing dist-weighted avg */
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

  for (unsigned n = 0; n < 4; ++n) nbr_dist[n] = ONE / (nbr_dist[n] + TINY);
  distance = 0.0;
  for (unsigned n = 0; n < 4; ++n) distance += nbr_dist[n];
  for (unsigned n = 0; n < 4; ++n) nbr_dist[n] /= distance;

  return search_result;
}

static unsigned
quad_cross_products(double plon, double plat, double lons[4], double lats[4])
{
  unsigned n;
  int scross[4], scross_last = 0;
  // Vectors for cross-product check
  double vec1_lat, vec1_lon;
  double vec2_lat, vec2_lon;

  /* For consistency, we must make sure all lons are in same 2pi interval */
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
        If all the cross products are positive, the point is contained in the
        box.
      */
      vec1_lat = lats[next_n] - lats[n];
      vec1_lon = lons[next_n] - lons[n];
      vec2_lat = plat - lats[n];
      vec2_lon = plon - lons[n];

      /* Check for 0,2pi crossings */
      if (vec1_lon > THREE * PIH)
        vec1_lon -= PI2;
      else if (vec1_lon < -THREE * PIH)
        vec1_lon += PI2;

      if (vec2_lon > THREE * PIH)
        vec2_lon -= PI2;
      else if (vec2_lon < -THREE * PIH)
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
point_in_quad(bool is_cyclic, size_t nx, size_t ny, size_t i, size_t j, size_t adds[4], double lons[4], double lats[4],
              double plon, double plat, const double *restrict center_lon, const double *restrict center_lat)
{
  bool search_result = false;
  size_t ip1 = (i < (nx - 1)) ? i + 1 : is_cyclic ? 0 : i;
  size_t jp1 = (j < (ny - 1)) ? j + 1 : j;

  if (i == ip1 || j == jp1) return search_result;

  size_t idx[4];
  idx[0] = j * nx + i;
  idx[1] = j * nx + ip1;    // east
  idx[2] = jp1 * nx + ip1;  // north-east
  idx[3] = jp1 * nx + i;    // north

  for (unsigned j = 0; j < 4; ++j) lons[j] = center_lon[idx[j]];
  for (unsigned j = 0; j < 4; ++j) lats[j] = center_lat[idx[j]];

  unsigned n = quad_cross_products(plon, plat, lons, lats);

  /* If cross products all same sign, we found the location */
  if (n >= 4)
    {
      for (unsigned j = 0; j < 4; ++j) adds[j] = idx[j];
      search_result = true;
    }

  return search_result;
}

int
grid_search(remapGridType *src_grid, size_t *restrict src_add, double *restrict src_lats, double *restrict src_lons,
            double plat, double plon, const size_t *restrict src_grid_dims, const double *restrict src_center_lat,
            const double *restrict src_center_lon, const float *restrict src_grid_bound_box,
            const size_t *restrict src_bin_add)
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

    int src_bin_add[][2]           ! latitude bins for restricting
  */
  /*  Local variables */
  size_t n2, srch_add, srch_add4; /* dummy indices                    */
  int search_result = 0;
  float *bin_lats = src_grid->bin_lats;

  size_t nbins = src_grid->num_srch_bins;

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
          if (src_bin_add[n2] < min_add) min_add = src_bin_add[n2];
          if (src_bin_add[n2 + 1] > max_add) max_add = src_bin_add[n2 + 1];
        }
    }

  /* Now perform a more detailed search */

  size_t nx = src_grid_dims[0];
  size_t ny = src_grid_dims[1];

  /* srch_loop */
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

          if (point_in_quad(src_grid->is_cyclic, nx, ny, i, j, src_add, src_lons, src_lats, plon, plat, src_center_lon,
                            src_center_lat))
            {
              search_result = 1;
              return search_result;
            }

          /* Otherwise move on to next cell */

        } /* Bounding box check */
    }     /* srch_loop */

  /*
    If no cell found, point is likely either in a box that straddles either pole
    or is outside the grid. Fall back to a distance-weighted average of the four
    closest points. Go ahead and compute weights here, but store in src_lats and
    return -add to prevent the parent routine from computing bilinear weights.
  */
  if (!src_grid->lextrapolate) return search_result;

  /*
    printf("Could not find location for %g %g\n", plat*RAD2DEG, plon*RAD2DEG);
    printf("Using nearest-neighbor average for this point\n");
  */
  search_result = grid_search_nn(min_add, max_add, src_add, src_lats, plat, plon, src_center_lat, src_center_lon);

  return search_result;
} /* grid_search */

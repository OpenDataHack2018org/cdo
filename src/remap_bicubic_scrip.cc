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
#include "remap_store_link.h"
#include "timer.h"

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/*                                                                         */
/*      BICUBIC INTERPOLATION                                              */
/*                                                                         */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

static void
set_bicubic_weights(double iw, double jw, double wgts[4][4])
{
  // clang-format off
  wgts[0][0] = (1.-jw*jw*(3.-2.*jw))  * (1.-iw*iw*(3.-2.*iw));
  wgts[1][0] = (1.-jw*jw*(3.-2.*jw))  *     iw*iw*(3.-2.*iw);
  wgts[2][0] =     jw*jw*(3.-2.*jw)   *     iw*iw*(3.-2.*iw);
  wgts[3][0] =     jw*jw*(3.-2.*jw)   * (1.-iw*iw*(3.-2.*iw));
  wgts[0][1] = (1.-jw*jw*(3.-2.*jw))  *     iw*(iw-1.)*(iw-1.);
  wgts[1][1] = (1.-jw*jw*(3.-2.*jw))  *     iw*iw*(iw-1.);
  wgts[2][1] =     jw*jw*(3.-2.*jw)   *     iw*iw*(iw-1.);
  wgts[3][1] =     jw*jw*(3.-2.*jw)   *     iw*(iw-1.)*(iw-1.);
  wgts[0][2] =     jw*(jw-1.)*(jw-1.) * (1.-iw*iw*(3.-2.*iw));
  wgts[1][2] =     jw*(jw-1.)*(jw-1.) *     iw*iw*(3.-2.*iw);
  wgts[2][2] =     jw*jw*(jw-1.)      *     iw*iw*(3.-2.*iw);
  wgts[3][2] =     jw*jw*(jw-1.)      * (1.-iw*iw*(3.-2.*iw));
  wgts[0][3] =     iw*(iw-1.)*(iw-1.) *     jw*(jw-1.)*(jw-1.);
  wgts[1][3] =     iw*iw*(iw-1.)      *     jw*(jw-1.)*(jw-1.);
  wgts[2][3] =     iw*iw*(iw-1.)      *     jw*jw*(jw-1.);
  wgts[3][3] =     iw*(iw-1.)*(iw-1.) *     jw*jw*(jw-1.);
  // clang-format on
}

unsigned num_src_points(const int *restrict mask, const size_t src_add[4], double src_lats[4]);

static void
renormalizeWeights(const double src_lats[4], double wgts[4][4])
{
  double sum_wgts = 0.0; // sum of weights for normalization
  for (unsigned n = 0; n < 4; ++n) sum_wgts += fabs(src_lats[n]);
  for (unsigned n = 0; n < 4; ++n) wgts[n][0] = fabs(src_lats[n]) / sum_wgts;
  for (unsigned n = 0; n < 4; ++n) wgts[n][1] = 0.;
  for (unsigned n = 0; n < 4; ++n) wgts[n][2] = 0.;
  for (unsigned n = 0; n < 4; ++n) wgts[n][3] = 0.;
}

static void
bicubic_warning(void)
{
  static bool lwarn = true;

  if (cdoVerbose || lwarn)
    {
      lwarn = false;
      // cdoWarning("Iteration for iw,jw exceed max iteration count of %d!",
      // remap_max_iter);
      cdoWarning("Bicubic interpolation failed for some grid points - used a "
                 "distance-weighted average instead!");
    }
}

static void
bicubic_remap(double *restrict tgt_point, const double *restrict src_array, double wgts[4][4], const size_t src_add[4],
              gradientsType &gradients)
{
  const double *restrict glat = &gradients.grad_lat[0];
  const double *restrict glon = &gradients.grad_lon[0];
  const double *restrict glatlon = &gradients.grad_latlon[0];

  *tgt_point = 0.;
  for (unsigned n = 0; n < 4; ++n)
    *tgt_point += src_array[src_add[n]] * wgts[n][0] + glat[src_add[n]] * wgts[n][1] + glon[src_add[n]] * wgts[n][2]
                  + glatlon[src_add[n]] * wgts[n][3];
}

/*
  -----------------------------------------------------------------------

  This routine computes the weights for a bicubic interpolation.

  -----------------------------------------------------------------------
*/
void
scrip_remap_bicubic_weights(RemapSearch &rsearch, RemapGrid *src_grid, RemapGrid *tgt_grid, RemapVars &rv)
{
  extern int timer_remap_bic;
  int remap_grid_type = src_grid->remap_grid_type;

  if (cdoVerbose) cdoPrint("Called %s()", __func__);

  if (cdoTimer) timer_start(timer_remap_bic);

  progressInit();

  // Compute mappings from source to target grid

  if (src_grid->rank != 2) cdoAbort("Can't do bicubic interpolation when source grid rank != 2");

  size_t tgt_grid_size = tgt_grid->size;

  std::vector<weightlinks4_t> weightlinks(tgt_grid_size);
  weightlinks[0].addweights = (addweight4_t *) Malloc(4 * tgt_grid_size * sizeof(addweight4_t));
  for (unsigned tgt_cell_add = 1; tgt_cell_add < tgt_grid_size; ++tgt_cell_add)
    weightlinks[tgt_cell_add].addweights = weightlinks[0].addweights + 4 * tgt_cell_add;

  // Loop over destination grid

  double findex = 0;

#ifdef HAVE_OPENMP4
#pragma omp parallel for default(none) reduction(+ : findex) shared(rsearch, weightlinks, remap_grid_type, tgt_grid_size, \
                                                                    src_grid, tgt_grid, rv)
#endif
  for (size_t tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add)
    {
      findex++;
      if (cdo_omp_get_thread_num() == 0) progressStatus(0, 1, findex / tgt_grid_size);

      weightlinks[tgt_cell_add].nlinks = 0;

      if (!tgt_grid->mask[tgt_cell_add]) continue;

      double plat = tgt_grid->cell_center_lat[tgt_cell_add];
      double plon = tgt_grid->cell_center_lon[tgt_cell_add];

      size_t src_add[4];   //  address for the four source points
      double src_lats[4];  //  latitudes  of four bilinear corners
      double src_lons[4];  //  longitudes of four bilinear corners
      double wgts[4][4];   //  bicubic weights for four corners

      // Find nearest square of grid points on source grid
      int search_result;
      if (remap_grid_type == REMAP_GRID_TYPE_REG2D)
        search_result = grid_search_reg2d(src_grid, src_add, src_lats, src_lons, plat, plon, src_grid->dims,
                                          src_grid->reg2d_center_lat, src_grid->reg2d_center_lon);
      else
        search_result
            = grid_search(src_grid, src_add, src_lats, src_lons, plat, plon, src_grid->dims, src_grid->cell_center_lat,
                          src_grid->cell_center_lon, rsearch.srcBins);

      // Check to see if points are land points
      if (search_result > 0)
        {
          for (unsigned n = 0; n < 4; ++n)
            if (!src_grid->mask[src_add[n]]) search_result = 0;
        }

      // If point found, find local iw,jw coordinates for weights
      if (search_result > 0)
        {
          tgt_grid->cell_frac[tgt_cell_add] = 1.;

          double iw, jw; // current guess for bilinear coordinate
          if (find_ij_weights(plon, plat, src_lons, src_lats, &iw, &jw))
            {
              // Successfully found iw,jw - compute weights
              set_bicubic_weights(iw, jw, wgts);
              store_weightlinks4(4, src_add, wgts, tgt_cell_add, &weightlinks[0]);
            }
          else
            {
              bicubic_warning();
              search_result = -1;
            }
        }

      /*
        Search for bicubic failed - use a distance-weighted average instead
        (this is typically near the pole) Distance was stored in src_lats!
      */
      if (search_result < 0)
        {
          if (num_src_points(src_grid->mask, src_add, src_lats) > 0)
            {
              tgt_grid->cell_frac[tgt_cell_add] = 1.;
              renormalizeWeights(src_lats, wgts);
              store_weightlinks4(4, src_add, wgts, tgt_cell_add, &weightlinks[0]);
            }
        }
    }

  if (cdoTimer) timer_stop(timer_remap_bic);

  weightlinks2remaplinks4(tgt_grid_size, &weightlinks[0], rv);

} /* scrip_remap_weights_bicubic */

/*
  -----------------------------------------------------------------------

  This routine computes ans apply the weights for a bicubic interpolation.

  -----------------------------------------------------------------------
*/

//#define TEST_KDTREE
#ifdef TEST_KDTREE
#include "grid_search.h"
int
grid_search_test(GridSearch *gs, size_t *restrict src_add, double *restrict src_lats, double *restrict src_lons,
                 double plat, double plon, const size_t *restrict src_grid_dims, double *restrict src_center_lat,
                 double *restrict src_center_lon)
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

  */
  bool is_cyclic = true;
  int search_result = 0;

  for (unsigned n = 0; n < 4; ++n) src_add[n] = 0;

  /* Now perform a more detailed search */

  size_t nx = src_grid_dims[0];
  size_t ny = src_grid_dims[1];

  double search_radius = gs->search_radius;
  const double range0 = SQR(search_radius);
  double range = range0;
  size_t add = gridsearch_nearest(gs, plon, plat, &range);
  // printf("plon, plat, add, range %g %g %g %g %zu %g\n", plon*RAD2DEG,
  // plat*RAD2DEG,
  //     src_center_lon[add]*RAD2DEG, src_center_lat[add]*RAD2DEG,add, range);
  if (add != GS_NOT_FOUND)
    {
      size_t idx[4];
      for (unsigned k = 0; k < 4; ++k)
        {
          /* Determine neighbor addresses */
          size_t j = add / nx;
          size_t i = add - j * nx;
          if (k == 1 || k == 3) i = (i > 0) ? i - 1 : (is_cyclic) ? nx - 1 : 0;
          if (k == 2 || k == 3) j = (j > 0) ? j - 1 : 0;

          if (point_in_quad(is_cyclic, nx, ny, i, j, src_add, src_lons, src_lats, plon, plat, src_center_lon,
                            src_center_lat))
            {
              search_result = 1;
              return search_result;
            }
          /* Otherwise move on to next cell */
        }
      /*
        If no cell found, point is likely either in a box that straddles either
        pole or is outside the grid. Fall back to a distance-weighted average of
        the four closest points. Go ahead and compute weights here, but store in
        src_lats and return -add to prevent the parent routine from computing
        bilinear weights.
      */
      // if ( !src_grid->lextrapolate ) return search_result;

      /*
        printf("Could not find location for %g %g\n", plat*RAD2DEG,
        plon*RAD2DEG); printf("Using nearest-neighbor for this point\n");
      */
      search_result = add;
    }

  return search_result;
} /* grid_search_test */
#endif

void
scrip_remap_bicubic(RemapSearch &rsearch, RemapGrid *src_grid, RemapGrid *tgt_grid, const double *restrict src_array,
                    double *restrict tgt_array, double missval)
{
  int remap_grid_type = src_grid->remap_grid_type;

  if (cdoVerbose) cdoPrint("Called %s()", __func__);

  progressInit();

#ifdef TEST_KDTREE
  bool xIsCyclic = false;
  size_t dims[2] = { src_grid->size, 0 };
  GridSearch *gs = NULL;
  if (remap_grid_type != REMAP_GRID_TYPE_REG2D)
    gs = gridsearch_create(xIsCyclic, dims, src_grid->size, src_grid->cell_center_lon, src_grid->cell_center_lat);
#endif

  size_t tgt_grid_size = tgt_grid->size;

  // Compute mappings from source to target grid

  if (src_grid->rank != 2) cdoAbort("Can't do bicubic interpolation when source grid rank != 2");

  gradientsType gradients(src_grid->size);
  remapGradients(*src_grid, src_array, gradients);

  // Loop over destination grid

  double findex = 0;

#ifdef HAVE_OPENMP4
#pragma omp parallel for default(none) reduction(+ : findex) shared(rsearch, remap_grid_type, tgt_grid_size, src_grid, \
                                                                    tgt_grid, src_array, tgt_array, missval,  \
                                                                    gradients)
#endif
  for (size_t tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add)
    {
      findex++;
      if (cdo_omp_get_thread_num() == 0) progressStatus(0, 1, findex / tgt_grid_size);

      tgt_array[tgt_cell_add] = missval;

      if (!tgt_grid->mask[tgt_cell_add]) continue;

      double plat = tgt_grid->cell_center_lat[tgt_cell_add];
      double plon = tgt_grid->cell_center_lon[tgt_cell_add];

      size_t src_add[4];   //  address for the four source points
      double src_lats[4];  //  latitudes  of four bilinear corners
      double src_lons[4];  //  longitudes of four bilinear corners
      double wgts[4][4];   //  bicubic weights for four corners

      // Find nearest square of grid points on source grid
      int search_result;
      if (remap_grid_type == REMAP_GRID_TYPE_REG2D)
        search_result = grid_search_reg2d(src_grid, src_add, src_lats, src_lons, plat, plon, src_grid->dims,
                                          src_grid->reg2d_center_lat, src_grid->reg2d_center_lon);
      else
#ifdef TEST_KDTREE
        search_result = grid_search_test(gs, src_add, src_lats, src_lons, plat, plon, src_grid->dims,
                                         src_grid->cell_center_lat, src_grid->cell_center_lon);
#else
        search_result
            = grid_search(src_grid, src_add, src_lats, src_lons, plat, plon, src_grid->dims, src_grid->cell_center_lat,
                          src_grid->cell_center_lon, rsearch.srcBins);
#endif

      // Check to see if points are land points
      if (search_result > 0)
        {
          for (unsigned n = 0; n < 4; ++n)
            if (!src_grid->mask[src_add[n]]) search_result = 0;
        }

      // If point found, find local iw,jw coordinates for weights
      if (search_result > 0)
        {
          tgt_grid->cell_frac[tgt_cell_add] = 1.;

          double iw, jw; // current guess for bilinear coordinate
          if (find_ij_weights(plon, plat, src_lons, src_lats, &iw, &jw))
            {
              // Successfully found iw,jw - compute weights
              set_bicubic_weights(iw, jw, wgts);
              sort_add_and_wgts4(4, src_add, wgts);
              bicubic_remap(&tgt_array[tgt_cell_add], src_array, wgts, src_add, gradients);
            }
          else
            {
              bicubic_warning();
              search_result = -1;
            }
        }

      /*
        Search for bicubic failed - use a distance-weighted average instead
        (this is typically near the pole) Distance was stored in src_lats!
      */
      if (search_result < 0)
        {
          if (num_src_points(src_grid->mask, src_add, src_lats) > 0)
            {
              tgt_grid->cell_frac[tgt_cell_add] = 1.;
              renormalizeWeights(src_lats, wgts);
              sort_add_and_wgts4(4, src_add, wgts);
              bicubic_remap(&tgt_array[tgt_cell_add], src_array, wgts, src_add, gradients);
            }
        }
    }

#ifdef TEST_KDTREE
  if (gs) gridsearch_delete(gs);
#endif

}  // scrip_remap_bicubic

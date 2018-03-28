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

// bicubic interpolation

static void
bicubicSetWeights(double iw, double jw, double wgts[4][4])
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
  double sum_wgts = 0.0;  // sum of weights for normalization
  for (unsigned n = 0; n < 4; ++n) sum_wgts += fabs(src_lats[n]);
  for (unsigned n = 0; n < 4; ++n) wgts[n][0] = fabs(src_lats[n]) / sum_wgts;
  for (unsigned n = 0; n < 4; ++n) wgts[n][1] = 0.;
  for (unsigned n = 0; n < 4; ++n) wgts[n][2] = 0.;
  for (unsigned n = 0; n < 4; ++n) wgts[n][3] = 0.;
}

static void
bicubicWarning(void)
{
  static bool lwarn = true;

  if (cdoVerbose || lwarn)
    {
      lwarn = false;
      cdoWarning("Bicubic interpolation failed for some grid points - used a distance-weighted average instead!");
    }
}

static void
bicubicRemap(double *restrict tgt_point, const double *restrict src_array, double wgts[4][4], const size_t src_add[4],
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
remapBicubicWeights(RemapSearch &rsearch, RemapVars &rv)
{
  RemapGrid *src_grid = rsearch.srcGrid;
  RemapGrid *tgt_grid = rsearch.tgtGrid;

  extern int timer_remap_bic;

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

  double findex = 0;

  // Loop over destination grid

#ifdef HAVE_OPENMP4
#pragma omp parallel for default(none) reduction(+ : findex) shared(rsearch, weightlinks, tgt_grid_size, src_grid, tgt_grid, rv)
#endif
  for (size_t tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add)
    {
      findex++;
      if (cdo_omp_get_thread_num() == 0) progressStatus(0, 1, findex / tgt_grid_size);

      weightlinks[tgt_cell_add].nlinks = 0;

      if (!tgt_grid->mask[tgt_cell_add]) continue;

      double plon = 0, plat = 0;
      remapgrid_get_lonlat(tgt_grid, tgt_cell_add, &plon, &plat);

      size_t src_add[4];   //  address for the four source points
      double src_lats[4];  //  latitudes  of four bilinear corners
      double src_lons[4];  //  longitudes of four bilinear corners
      double wgts[4][4];   //  bicubic weights for four corners

      // Find nearest square of grid points on source grid
      int search_result = remapSearchSquare(rsearch, plon, plat, src_add, src_lats, src_lons);

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

          double iw, jw;  // current guess for bilinear coordinate
          if (remapFindWeights(plon, plat, src_lons, src_lats, &iw, &jw))
            {
              // Successfully found iw,jw - compute weights
              bicubicSetWeights(iw, jw, wgts);
              storeWeightlinks4(4, src_add, wgts, tgt_cell_add, weightlinks);
            }
          else
            {
              bicubicWarning();
              search_result = -1;
            }
        }

      /*
        Search for bicubic failed - use a distance-weighted average instead
        (this is typically near the pole) Distance was stored in src_lats!
      */
      if (search_result < 0)
        {
          if (num_src_points(&src_grid->mask[0], src_add, src_lats) > 0)
            {
              tgt_grid->cell_frac[tgt_cell_add] = 1.;
              renormalizeWeights(src_lats, wgts);
              storeWeightlinks4(4, src_add, wgts, tgt_cell_add, weightlinks);
            }
        }
    }

  weightlinks2remaplinks4(tgt_grid_size, weightlinks, rv);

  if (cdoTimer) timer_stop(timer_remap_bic);
}  // scrip_remap_weights_bicubic

/*
  -----------------------------------------------------------------------

  This routine computes and apply the weights for a bicubic interpolation.

  -----------------------------------------------------------------------
*/

void
remapBicubic(RemapSearch &rsearch, const double *restrict src_array, double *restrict tgt_array, double missval)
{
  RemapGrid *src_grid = rsearch.srcGrid;
  RemapGrid *tgt_grid = rsearch.tgtGrid;

  if (cdoVerbose) cdoPrint("Called %s()", __func__);

  extern int timer_remap_bic;
  if (cdoTimer) timer_start(timer_remap_bic);

  progressInit();

  size_t tgt_grid_size = tgt_grid->size;

  // Compute mappings from source to target grid

  if (src_grid->rank != 2) cdoAbort("Can't do bicubic interpolation when source grid rank != 2");

  gradientsType gradients(src_grid->size);
  remapGradients(*src_grid, src_array, gradients);

  double findex = 0;

  // Loop over destination grid

#ifdef HAVE_OPENMP4
#pragma omp parallel for default(none) reduction(+ : findex) \
  shared(rsearch, tgt_grid_size, src_grid, tgt_grid, src_array, tgt_array, missval, gradients)
#endif
  for (size_t tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add)
    {
      findex++;
      if (cdo_omp_get_thread_num() == 0) progressStatus(0, 1, findex / tgt_grid_size);

      tgt_array[tgt_cell_add] = missval;

      if (!tgt_grid->mask[tgt_cell_add]) continue;

      double plon = 0, plat = 0;
      remapgrid_get_lonlat(tgt_grid, tgt_cell_add, &plon, &plat);

      size_t src_add[4];   //  address for the four source points
      double src_lats[4];  //  latitudes  of four bilinear corners
      double src_lons[4];  //  longitudes of four bilinear corners
      double wgts[4][4];   //  bicubic weights for four corners

      // Find nearest square of grid points on source grid
      int search_result = remapSearchSquare(rsearch, plon, plat, src_add, src_lats, src_lons);

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

          double iw, jw;  // current guess for bilinear coordinate
          if (remapFindWeights(plon, plat, src_lons, src_lats, &iw, &jw))
            {
              // Successfully found iw,jw - compute weights
              bicubicSetWeights(iw, jw, wgts);
              sort_add_and_wgts4(4, src_add, wgts);
              bicubicRemap(&tgt_array[tgt_cell_add], src_array, wgts, src_add, gradients);
            }
          else
            {
              bicubicWarning();
              search_result = -1;
            }
        }

      /*
        Search for bicubic failed - use a distance-weighted average instead
        (this is typically near the pole) Distance was stored in src_lats!
      */
      if (search_result < 0)
        {
          if (num_src_points(&src_grid->mask[0], src_add, src_lats) > 0)
            {
              tgt_grid->cell_frac[tgt_cell_add] = 1.;
              renormalizeWeights(src_lats, wgts);
              sort_add_and_wgts4(4, src_add, wgts);
              bicubicRemap(&tgt_array[tgt_cell_add], src_array, wgts, src_add, gradients);
            }
        }
    }

  if (cdoTimer) timer_stop(timer_remap_bic);
}  // remapBicubic

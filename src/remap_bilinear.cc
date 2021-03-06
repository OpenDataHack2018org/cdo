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
#include "cdo_wtime.h"
#include "grid.h"
#include "remap.h"
#include "remap_store_link.h"

// bilinear interpolation

bool
remapFindWeights(double plon, double plat, double *restrict src_lons, double *restrict src_lats, double *ig, double *jg)
{
  long iter;                           // iteration counters
  double deli, delj;                   // corrections to iw,jw
  double dthp, dphp;                   // difference between point and sw corner
  double mat1, mat2, mat3, mat4;       // matrix elements
  double determinant;                  // matrix determinant
  constexpr double converge = 1.e-10;  // Convergence criterion
  extern long remap_max_iter;

  /* Iterate to find iw,jw for bilinear approximation  */

  // some latitude  differences
  double dth1 = src_lats[1] - src_lats[0];
  double dth2 = src_lats[3] - src_lats[0];
  double dth3 = src_lats[2] - src_lats[1] - dth2;

  // some longitude differences
  double dph1 = src_lons[1] - src_lons[0];
  double dph2 = src_lons[3] - src_lons[0];
  double dph3 = src_lons[2] - src_lons[1];

  if (dph1 > 3.0 * PIH) dph1 -= PI2;
  if (dph2 > 3.0 * PIH) dph2 -= PI2;
  if (dph3 > 3.0 * PIH) dph3 -= PI2;
  if (dph1 < -3.0 * PIH) dph1 += PI2;
  if (dph2 < -3.0 * PIH) dph2 += PI2;
  if (dph3 < -3.0 * PIH) dph3 += PI2;

  dph3 = dph3 - dph2;

  // current guess for bilinear coordinate
  double iguess = 0.5;
  double jguess = 0.5;

  for (iter = 0; iter < remap_max_iter; ++iter)
    {
      dthp = plat - src_lats[0] - dth1 * iguess - dth2 * jguess - dth3 * iguess * jguess;
      dphp = plon - src_lons[0];

      if (dphp > 3.0 * PIH) dphp -= PI2;
      if (dphp < -3.0 * PIH) dphp += PI2;

      dphp = dphp - dph1 * iguess - dph2 * jguess - dph3 * iguess * jguess;

      mat1 = dth1 + dth3 * jguess;
      mat2 = dth2 + dth3 * iguess;
      mat3 = dph1 + dph3 * jguess;
      mat4 = dph2 + dph3 * iguess;

      determinant = mat1 * mat4 - mat2 * mat3;

      deli = (dthp * mat4 - dphp * mat2) / determinant;
      delj = (dphp * mat1 - dthp * mat3) / determinant;

      if (fabs(deli) < converge && fabs(delj) < converge) break;

      iguess += deli;
      jguess += delj;
    }

  *ig = iguess;
  *jg = jguess;

  return (iter < remap_max_iter);
}

static void
bilinearSetWeights(double iw, double jw, double wgts[4])
{
  // clang-format off
  wgts[0] = (1.-iw) * (1.-jw);
  wgts[1] =     iw  * (1.-jw);
  wgts[2] =     iw  *     jw;
  wgts[3] = (1.-iw) *     jw;
  // clang-format on
}

unsigned
num_src_points(const int *restrict mask, const size_t src_add[4], double src_lats[4])
{
  unsigned icount = 0;

  for (unsigned n = 0; n < 4; ++n)
    {
      if (mask[src_add[n]])
        icount++;
      else
        src_lats[n] = 0.;
    }

  return icount;
}

static void
renormalizeWeights(const double src_lats[4], double wgts[4])
{
  double sum_wgts = 0.0;  // sum of weights for normalization
  for (unsigned n = 0; n < 4; ++n) sum_wgts += fabs(src_lats[n]);
  for (unsigned n = 0; n < 4; ++n) wgts[n] = fabs(src_lats[n]) / sum_wgts;
}

static void
bilinearWarning(void)
{
  static bool lwarn = true;

  if (cdoVerbose || lwarn)
    {
      lwarn = false;
      cdoWarning("Bilinear interpolation failed for some grid points - used a distance-weighted average instead!");
    }
}

static void
bilinearRemap(double *restrict tgt_point, const double *restrict src_array, const double wgts[4], const size_t src_add[4])
{
  // *tgt_point = 0.;
  // for ( unsigned n = 0; n < 4; ++n ) *tgt_point += src_array[src_add[n]]*wgts[n];
  *tgt_point = src_array[src_add[0]] * wgts[0] + src_array[src_add[1]] * wgts[1] + src_array[src_add[2]] * wgts[2]
               + src_array[src_add[3]] * wgts[3];
}

/*
  -----------------------------------------------------------------------

  This routine computes the weights for a bilinear interpolation.

  -----------------------------------------------------------------------
*/
void
remapBilinearWeights(RemapSearch &rsearch, RemapVars &rv)
{
  RemapGrid *src_grid = rsearch.srcGrid;
  RemapGrid *tgt_grid = rsearch.tgtGrid;

  if (cdoVerbose) cdoPrint("Called %s()", __func__);

  if (src_grid->rank != 2) cdoAbort("Can't do bilinear interpolation when source grid rank != 2");

  double start = cdoVerbose ? cdo_get_wtime() : 0;

  progressInit();

  // Compute mappings from source to target grid

  size_t tgt_grid_size = tgt_grid->size;

  std::vector<WeightLinks> weightLinks(tgt_grid_size);
  weightLinksAlloc(4, tgt_grid_size, weightLinks);

  double findex = 0;

// Loop over destination grid

#ifdef _OPENMP
#pragma omp parallel for default(none) schedule(static) shared(findex, rsearch, weightLinks, tgt_grid_size, src_grid, tgt_grid, rv)
#endif
  for (size_t tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add)
    {
#ifdef _OPENMP
#pragma omp atomic
#endif
      findex++;
      if (cdo_omp_get_thread_num() == 0) progressStatus(0, 1, findex / tgt_grid_size);

      weightLinks[tgt_cell_add].nlinks = 0;

      if (!tgt_grid->mask[tgt_cell_add]) continue;

      double plon = 0, plat = 0;
      remapgrid_get_lonlat(tgt_grid, tgt_cell_add, &plon, &plat);

      size_t src_add[4];   //  address for the four source points
      double src_lats[4];  //  latitudes  of four bilinear corners
      double src_lons[4];  //  longitudes of four bilinear corners
      double wgts[4];      //  bilinear weights for four corners

      // Find nearest square of grid points on source grid
      int search_result = remapSearchSquare(rsearch, plon, plat, src_add, src_lats, src_lons);

      // Check to see if points are mask points
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
              bilinearSetWeights(iw, jw, wgts);
              storeWeightLinks(0, 4, src_add, wgts, tgt_cell_add, weightLinks);
            }
          else
            {
              bilinearWarning();
              search_result = -1;
            }
        }

      /*
        Search for bilinear failed - use a distance-weighted average instead
        (this is typically near the pole) Distance was stored in src_lats!
      */
      if (search_result < 0)
        {
          if (num_src_points(&src_grid->mask[0], src_add, src_lats) > 0)
            {
              tgt_grid->cell_frac[tgt_cell_add] = 1.;
              renormalizeWeights(src_lats, wgts);
              storeWeightLinks(0, 4, src_add, wgts, tgt_cell_add, weightLinks);
            }
        }
    }

  progressStatus(0, 1, 1);

  weightLinksToRemapLinks(0, tgt_grid_size, weightLinks, rv);

  if (cdoVerbose) cdoPrint("Square search nearest: %.2f seconds", cdo_get_wtime() - start);
}  // scrip_remap_weights_bilinear

/*
  -----------------------------------------------------------------------

  This routine computes and apply the weights for a bilinear interpolation.

  -----------------------------------------------------------------------
*/

void
remapBilinear(RemapSearch &rsearch, const double *restrict src_array, double *restrict tgt_array, double missval)
{
  RemapGrid *src_grid = rsearch.srcGrid;
  RemapGrid *tgt_grid = rsearch.tgtGrid;

  if (cdoVerbose) cdoPrint("Called %s()", __func__);

  if (src_grid->rank != 2) cdoAbort("Can't do bilinear interpolation when source grid rank != 2");

  double start = cdoVerbose ? cdo_get_wtime() : 0;

  progressInit();

  size_t tgt_grid_size = tgt_grid->size;

  // Compute mappings from source to target grid

  double findex = 0;

// Loop over destination grid

#ifdef _OPENMP
#pragma omp parallel for default(none) schedule(static) shared(findex, rsearch, tgt_grid_size, src_grid, tgt_grid, src_array, \
                                                               tgt_array, missval)
#endif
  for (size_t tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add)
    {
#ifdef _OPENMP
#pragma omp atomic
#endif
      findex++;
      if (cdo_omp_get_thread_num() == 0) progressStatus(0, 1, findex / tgt_grid_size);

      tgt_array[tgt_cell_add] = missval;

      if (!tgt_grid->mask[tgt_cell_add]) continue;

      double plon = 0, plat = 0;
      remapgrid_get_lonlat(tgt_grid, tgt_cell_add, &plon, &plat);

      size_t src_add[4];   //  address for the four source points
      double src_lats[4];  //  latitudes  of four bilinear corners
      double src_lons[4];  //  longitudes of four bilinear corners
      double wgts[4];      //  bilinear weights for four corners

      // Find nearest square of grid points on source grid
      int search_result = remapSearchSquare(rsearch, plon, plat, src_add, src_lats, src_lons);

      // Check to see if points are mask points
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
              bilinearSetWeights(iw, jw, wgts);
              sort_add_and_wgts(4, src_add, wgts);
              bilinearRemap(&tgt_array[tgt_cell_add], src_array, wgts, src_add);
            }
          else
            {
              bilinearWarning();
              search_result = -1;
            }
        }

      /*
        Search for bilinear failed - use a distance-weighted average instead
        (this is typically near the pole) Distance was stored in src_lats!
      */
      if (search_result < 0)
        {
          if (num_src_points(&src_grid->mask[0], src_add, src_lats) > 0)
            {
              tgt_grid->cell_frac[tgt_cell_add] = 1.;
              renormalizeWeights(src_lats, wgts);
              sort_add_and_wgts(4, src_add, wgts);
              bilinearRemap(&tgt_array[tgt_cell_add], src_array, wgts, src_add);
            }
        }
    }

  progressStatus(0, 1, 1);

  if (cdoVerbose) cdoPrint("Square search nearest: %.2f seconds", cdo_get_wtime() - start);
}  // remapBilinear

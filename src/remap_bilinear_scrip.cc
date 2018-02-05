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
#include "cdoOptions.h"


/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/*                                                                         */
/*      BILINEAR INTERPOLATION                                             */
/*                                                                         */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */


bool find_ij_weights(double plon, double plat, double *restrict src_lons, double *restrict src_lats, double *ig, double *jg)
{
  bool lfound = false;
  long iter;                     /*  iteration counters   */
  double deli, delj;             /*  corrections to iw,jw                   */
  double dthp, dphp;             /*  difference between point and sw corner */
  double mat1, mat2, mat3, mat4; /*  matrix elements                        */
  double determinant;            /*  matrix determinant                     */
  double converge = 1.e-10;      /* Convergence criterion                   */
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

  if ( dph1 >  THREE*PIH ) dph1 -= PI2;
  if ( dph2 >  THREE*PIH ) dph2 -= PI2;
  if ( dph3 >  THREE*PIH ) dph3 -= PI2;
  if ( dph1 < -THREE*PIH ) dph1 += PI2;
  if ( dph2 < -THREE*PIH ) dph2 += PI2;
  if ( dph3 < -THREE*PIH ) dph3 += PI2;

  dph3 = dph3 - dph2;

  // current guess for bilinear coordinate
  double iguess = HALF;
  double jguess = HALF;

  for ( iter = 0; iter < remap_max_iter; ++iter )
    {
      dthp = plat - src_lats[0] - dth1*iguess - dth2*jguess - dth3*iguess*jguess;
      dphp = plon - src_lons[0];
      
      if ( dphp >  THREE*PIH ) dphp -= PI2;
      if ( dphp < -THREE*PIH ) dphp += PI2;

      dphp = dphp - dph1*iguess - dph2*jguess - dph3*iguess*jguess;

      mat1 = dth1 + dth3*jguess;
      mat2 = dth2 + dth3*iguess;
      mat3 = dph1 + dph3*jguess;
      mat4 = dph2 + dph3*iguess;

      determinant = mat1*mat4 - mat2*mat3;

      deli = (dthp*mat4 - dphp*mat2)/determinant;
      delj = (dphp*mat1 - dthp*mat3)/determinant;

      if ( fabs(deli) < converge && fabs(delj) < converge ) break;

      iguess += deli;
      jguess += delj;
    }

  *ig = iguess;
  *jg = jguess;

  if ( iter < remap_max_iter ) lfound = true;

  return lfound;
}

static
void set_bilinear_weights(double iw, double jw, double wgts[4])
{
  // clang-format off
  wgts[0] = (1.-iw) * (1.-jw);
  wgts[1] =     iw  * (1.-jw);
  wgts[2] =     iw  *     jw;
  wgts[3] = (1.-iw) *     jw;
  // clang-format on
}


unsigned num_src_points(const int* restrict mask, const size_t src_add[4], double src_lats[4])
{
  unsigned icount = 0;

  for ( unsigned n = 0; n < 4; ++n )
    {
      if ( mask[src_add[n]] )
	icount++;
      else
	src_lats[n] = 0.;
    }

  return icount;
}

static
void renormalize_weights(const double src_lats[4], double wgts[4])
{
  double sum_wgts = 0.0; /* sum of weights for normalization */
  /* 2012-05-08 Uwe Schulzweida: using absolute value of src_lats (bug fix) */
  for ( unsigned n = 0; n < 4; ++n ) sum_wgts += fabs(src_lats[n]);
  for ( unsigned n = 0; n < 4; ++n ) wgts[n] = fabs(src_lats[n])/sum_wgts;
}

static
void bilinear_warning(double plon, double plat, double iw, double jw, size_t* src_add, double* src_lons, double* src_lats, remapgrid_t* src_grid)
{
  static bool lwarn = true;

  if ( cdoVerbose )
    {
      cdoPrint("Point coords: %g %g", plat*RAD2DEG, plon*RAD2DEG);
      cdoPrint("Src grid lats: %g %g %g %g", src_lats[0]*RAD2DEG, src_lats[1]*RAD2DEG, src_lats[2]*RAD2DEG, src_lats[3]*RAD2DEG);
      cdoPrint("Src grid lons: %g %g %g %g", src_lons[0]*RAD2DEG, src_lons[1]*RAD2DEG, src_lons[2]*RAD2DEG, src_lons[3]*RAD2DEG);
      cdoPrint("Src grid addresses: %zu %zu %zu %zu", src_add[0], src_add[1], src_add[2], src_add[3]);
      cdoPrint("Src grid lats: %g %g %g %g",
	       src_grid->cell_center_lat[src_add[0]]*RAD2DEG, src_grid->cell_center_lat[src_add[1]]*RAD2DEG,
	       src_grid->cell_center_lat[src_add[2]]*RAD2DEG, src_grid->cell_center_lat[src_add[3]]*RAD2DEG);
      cdoPrint("Src grid lons: %g %g %g %g",
	       src_grid->cell_center_lon[src_add[0]]*RAD2DEG, src_grid->cell_center_lon[src_add[1]]*RAD2DEG,
	       src_grid->cell_center_lon[src_add[2]]*RAD2DEG, src_grid->cell_center_lon[src_add[3]]*RAD2DEG);
      cdoPrint("Current iw,jw : %g %g", iw, jw);
    }

  if ( cdoVerbose || lwarn )
    {
      lwarn = false;
      //  cdoWarning("Iteration for iw,jw exceed max iteration count of %d!", remap_max_iter);
      cdoWarning("Bilinear interpolation failed for some grid points - used a distance-weighted average instead!");
    }
}

static
void bilinear_remap(double* restrict tgt_point, const double *restrict src_array, const double wgts[4], const size_t src_add[4])
{
  // *tgt_point = 0.;
  // for ( unsigned n = 0; n < 4; ++n ) *tgt_point += src_array[src_add[n]]*wgts[n];
  *tgt_point = src_array[src_add[0]]*wgts[0] + src_array[src_add[1]]*wgts[1]
             + src_array[src_add[2]]*wgts[2] + src_array[src_add[3]]*wgts[3];
}

/*
  -----------------------------------------------------------------------

  This routine computes the weights for a bilinear interpolation.

  -----------------------------------------------------------------------
*/
void scrip_remap_bilinear_weights(remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv)
{
  extern int timer_remap_bil;
  int remap_grid_type = src_grid->remap_grid_type;

  if ( cdoVerbose ) cdoPrint("Called %s()", __func__);

  if ( cdoTimer ) timer_start(timer_remap_bil);

  progressInit();

  /* Compute mappings from source to target grid */

  if ( src_grid->rank != 2 )
    cdoAbort("Can not do bilinear interpolation when source grid rank != 2"); 

  size_t tgt_grid_size = tgt_grid->size;

  weightlinks_t *weightlinks = (weightlinks_t *) Malloc(tgt_grid_size*sizeof(weightlinks_t));
  weightlinks[0].addweights = (addweight_t *) Malloc(4*tgt_grid_size*sizeof(addweight_t));
  for ( size_t tgt_cell_add = 1; tgt_cell_add < tgt_grid_size; ++tgt_cell_add )
    weightlinks[tgt_cell_add].addweights = weightlinks[0].addweights + 4*tgt_cell_add;

  double findex = 0;

  /* Loop over destination grid */

#ifdef  HAVE_OPENMP4
#pragma omp parallel for default(none)  schedule(static)  reduction(+:findex) \
  shared(weightlinks, remap_grid_type, tgt_grid_size, src_grid, tgt_grid, rv)
#endif
  for ( size_t tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add )
    {
      findex++;
      if ( cdo_omp_get_thread_num() == 0 ) progressStatus(0, 1, findex/tgt_grid_size);

      weightlinks[tgt_cell_add].nlinks = 0;	

      if ( ! tgt_grid->mask[tgt_cell_add] ) continue;

      double plon = 0, plat = 0;
      remapgrid_get_lonlat(tgt_grid, tgt_cell_add, &plon, &plat);

      size_t src_add[4];    //  address for the four source points
      double src_lats[4];   //  latitudes  of four bilinear corners
      double src_lons[4];   //  longitudes of four bilinear corners
      double wgts[4];       //  bilinear weights for four corners

      // Find nearest square of grid points on source grid
      int search_result;
      if ( remap_grid_type == REMAP_GRID_TYPE_REG2D )
	search_result = grid_search_reg2d(src_grid, src_add, src_lats, src_lons, 
					  plat, plon, src_grid->dims,
					  src_grid->reg2d_center_lat, src_grid->reg2d_center_lon);
      else
	search_result = grid_search(src_grid, src_add, src_lats, src_lons, 
				    plat, plon, src_grid->dims,
				    src_grid->cell_center_lat, src_grid->cell_center_lon,
				    src_grid->cell_bound_box, src_grid->bin_addr);

      // Check to see if points are mask points
      if ( search_result > 0 )
	{
	  for ( unsigned n = 0; n < 4; ++n )
	    if ( ! src_grid->mask[src_add[n]] ) search_result = 0;
	}

      // If point found, find local iw,jw coordinates for weights
      if ( search_result > 0 )
	{
          tgt_grid->cell_frac[tgt_cell_add] = 1.;

	  double iw, jw;  // current guess for bilinear coordinate 
          if ( find_ij_weights(plon, plat, src_lons, src_lats, &iw, &jw) )
	    {
	      // Successfully found iw,jw - compute weights
	      set_bilinear_weights(iw, jw, wgts);

	      store_weightlinks(0, 4, src_add, wgts, tgt_cell_add, weightlinks);
	    }
          else
	    {
	      bilinear_warning(plon, plat, iw, jw, src_add, src_lons, src_lats, src_grid);

	      search_result = -1;
	    }
	}

      /*
	Search for bilinear failed - use a distance-weighted average instead (this is typically near the pole)
	Distance was stored in src_lats!
      */
      if ( search_result < 0 )
	{
          if ( num_src_points(src_grid->mask, src_add, src_lats) > 0 )
	    {
	      renormalize_weights(src_lats, wgts);

	      tgt_grid->cell_frac[tgt_cell_add] = 1.;

	      store_weightlinks(0, 4, src_add, wgts, tgt_cell_add, weightlinks);
	    }
        }
    }

  weightlinks2remaplinks(0, tgt_grid_size, weightlinks, rv);

  if ( weightlinks ) Free(weightlinks);

  if ( cdoTimer ) timer_stop(timer_remap_bil);
} // scrip_remap_weights_bilinear

/*
  -----------------------------------------------------------------------

  This routine computes and apply the weights for a bilinear interpolation.

  -----------------------------------------------------------------------
*/

//#define TEST_KDTREE
#ifdef  TEST_KDTREE
#include "grid_search.h"
int grid_search_test(struct gridsearch *gs, size_t *restrict src_add, double *restrict src_lats, 
                     double *restrict src_lons,  double plat, double plon, const size_t *restrict src_grid_dims,
                     double *restrict src_center_lat, double *restrict src_center_lon)
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

  for ( unsigned n = 0; n < 4; ++n ) src_add[n] = 0;
 
  /* Now perform a more detailed search */

  size_t nx = src_grid_dims[0];
  size_t ny = src_grid_dims[1];

  double search_radius = gs->search_radius;
  const double range0 = SQR(search_radius);
  double range = range0;
  size_t add = gridsearch_nearest(gs, plon, plat, &range);
  // printf("plon, plat, add, range %g %g %g %g %zu %g\n", plon*RAD2DEG, plat*RAD2DEG,
  //     src_center_lon[add]*RAD2DEG, src_center_lat[add]*RAD2DEG,add, range);
  if ( add != GS_NOT_FOUND )
    {
      size_t idx[4];
      for ( unsigned k = 0; k < 4; ++k )
        {
          /* Determine neighbor addresses */
          size_t j = add/nx;
          size_t i = add - j*nx;
          if ( k == 1 || k == 3 )  i = (i > 0) ? i - 1 : (is_cyclic) ? nx-1 : 0;
          if ( k == 2 || k == 3 )  j = (j > 0) ? j - 1 : 0;

          if ( point_in_quad(is_cyclic, nx, ny, i, j, src_add, src_lons, src_lats,
                             plon, plat, src_center_lon, src_center_lat) )
	    {
	      search_result = 1;
	      return search_result;
	    }
	  /* Otherwise move on to next cell */
        }
      /*
        If no cell found, point is likely either in a box that straddles either pole or is outside 
        the grid. Fall back to a distance-weighted average of the four closest points.
        Go ahead and compute weights here, but store in src_lats and return -add to prevent the 
        parent routine from computing bilinear weights.
      */
      // if ( !src_grid->lextrapolate ) return search_result;

      /*
        printf("Could not find location for %g %g\n", plat*RAD2DEG, plon*RAD2DEG);
        printf("Using nearest-neighbor for this point\n");
      */
      search_result = add;
    }

  return search_result;
}  /* grid_search_test */
#endif

void scrip_remap_bilinear(remapgrid_t *src_grid, remapgrid_t *tgt_grid, const double *restrict src_array, double *restrict tgt_array, double missval)
{
  extern int timer_remap_bil;
  int remap_grid_type = src_grid->remap_grid_type;

  if ( cdoVerbose ) cdoPrint("Called %s()", __func__);

  if ( cdoTimer ) timer_start(timer_remap_bil);

  progressInit();

#ifdef  TEST_KDTREE
  bool xIsCyclic = false;
  size_t dims[2] = {src_grid->size, 0};
  struct gridsearch *gs = NULL;
  if ( remap_grid_type != REMAP_GRID_TYPE_REG2D )
    gs = gridsearch_create(xIsCyclic, dims, src_grid->size, src_grid->cell_center_lon, src_grid->cell_center_lat);
#else
  void *gs;
#endif

  size_t tgt_grid_size = tgt_grid->size;

  /* Compute mappings from source to target grid */

  if ( src_grid->rank != 2 )
    cdoAbort("Can not do bilinear interpolation when source grid rank != 2"); 

  double findex = 0;

  /* Loop over destination grid */

#ifdef  HAVE_OPENMP4
#pragma omp parallel for default(none)  schedule(static)  reduction(+:findex) \
  shared(gs) \
  shared(Options::silentMode, remap_grid_type, tgt_grid_size, src_grid, tgt_grid, src_array, tgt_array, missval)
#endif
  for ( size_t tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add )
    {
      findex++;
      if ( cdo_omp_get_thread_num() == 0 ) progressStatus(0, 1, findex/tgt_grid_size);

      tgt_array[tgt_cell_add] = missval;

      if ( ! tgt_grid->mask[tgt_cell_add] ) continue;

      double plon = 0, plat = 0;
      remapgrid_get_lonlat(tgt_grid, tgt_cell_add, &plon, &plat);

      size_t src_add[4];    //  address for the four source points
      double src_lats[4];   //  latitudes  of four bilinear corners
      double src_lons[4];   //  longitudes of four bilinear corners
      double wgts[4];       //  bilinear weights for four corners

      // Find nearest square of grid points on source grid
      int search_result;
      if ( remap_grid_type == REMAP_GRID_TYPE_REG2D )
	search_result = grid_search_reg2d(src_grid, src_add, src_lats, src_lons, 
					  plat, plon, src_grid->dims,
					  src_grid->reg2d_center_lat, src_grid->reg2d_center_lon);
      else
#ifdef  TEST_KDTREE
        search_result = grid_search_test(gs, src_add, src_lats, src_lons, plat, plon, src_grid->dims,
                                         src_grid->cell_center_lat, src_grid->cell_center_lon);
#else
	search_result = grid_search(src_grid, src_add, src_lats, src_lons, 
				    plat, plon, src_grid->dims,
				    src_grid->cell_center_lat, src_grid->cell_center_lon,
				    src_grid->cell_bound_box, src_grid->bin_addr);
#endif

      // Check to see if points are mask points
      if ( search_result > 0 )
	{
	  for ( unsigned n = 0; n < 4; ++n )
	    if ( ! src_grid->mask[src_add[n]] ) search_result = 0;
	}

      // If point found, find local iw,jw coordinates for weights
      if ( search_result > 0 )
	{
          tgt_grid->cell_frac[tgt_cell_add] = 1.;

	  double iw, jw;  // current guess for bilinear coordinate
          if ( find_ij_weights(plon, plat, src_lons, src_lats, &iw, &jw) )
	    {
	      // Successfully found iw,jw - compute weights
	      set_bilinear_weights(iw, jw, wgts);

	      sort_add_and_wgts(4, src_add, wgts);

	      bilinear_remap(&tgt_array[tgt_cell_add], src_array, wgts, src_add);
	    }
          else
	    {
	      bilinear_warning(plon, plat, iw, jw, src_add, src_lons, src_lats, src_grid);

	      search_result = -1;
	    }
	}

      /*
	Search for bilinear failed - use a distance-weighted average instead (this is typically near the pole)
	Distance was stored in src_lats!
      */
      if ( search_result < 0 )
	{
          if ( num_src_points(src_grid->mask, src_add, src_lats) > 0 )
	    {
	      renormalize_weights(src_lats, wgts);

	      tgt_grid->cell_frac[tgt_cell_add] = 1.;

	      sort_add_and_wgts(4, src_add, wgts);

	      bilinear_remap(&tgt_array[tgt_cell_add], src_array, wgts, src_add);
	    }
        }
    }

#ifdef  TEST_KDTREE
  if ( gs ) gridsearch_delete(gs);
#endif

  if ( cdoTimer ) timer_stop(timer_remap_bil);
} // scrip_remap_bilinear

#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "remap.h"


/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/*                                                                         */
/*      BICUBIC INTERPOLATION                                              */
/*                                                                         */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

/*
  This routine stores the address and weight for four links associated with one destination 
  point in the appropriate address and weight arrays and resizes those arrays if necessary.
*/
static
void store_link_bicub(remapvars_t *rv, int dst_add, int *restrict src_add, double weights[4][4])
{
  /*
    Input variables:
    int dst_add          ! address on destination grid
    int src_add[4]       ! addresses on source grid
    double weights[4][4] ! array of remapping weights for these links
  */
  long n, k;
  long nlink;

  /*
     Increment number of links and check to see if remap arrays need
     to be increased to accomodate the new link. Then store the link.
  */
  nlink = rv->num_links;
  rv->num_links += 4;

  if ( rv->num_links >= rv->max_links ) 
    resize_remap_vars(rv, rv->resize_increment);

  if ( rv->sort_add == FALSE )
    {
      for ( n = 0; n < 3; ++n )
	if ( src_add[n] > src_add[n+1] ) break;

      if ( n < 2 ) rv->sort_add = TRUE;
      else if ( n == 2 ) // swap 3rd and 4th elements
	{
	  {
	    int itmp;
	    itmp = src_add[2];
	    src_add[2] = src_add[3];
	    src_add[3] = itmp;
	  }
	  {
	    double dtmp[4];
	    for ( k = 0; k < 4; ++k ) dtmp[k] = weights[k][2];
	    for ( k = 0; k < 4; ++k ) weights[k][2] = weights[k][3];
	    for ( k = 0; k < 4; ++k ) weights[k][3] = dtmp[k];
	  }
	}
    }

  for ( n = 0; n < 4; ++n )
    {
      rv->src_grid_add[nlink+n] = src_add[n];
      rv->tgt_grid_add[nlink+n] = dst_add;
      for ( k = 0; k < 4; ++k )
	rv->wts[4*(nlink+n)+k] = weights[k][n];
    }

} /* store_link_bicub */

/*
  -----------------------------------------------------------------------

  This routine computes the weights for a bicubic interpolation.

  -----------------------------------------------------------------------
*/
void remap_bicubic(remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv)
{
  /*   Local variables */
  int  lwarn = TRUE;
  int  search_result;
  long tgt_grid_size;
  long n, icount;
  long dst_add;        /*  destination addresss */
  long iter;           /*  iteration counters   */

  int src_add[4];     /* address for the four source points */

  double src_lats[4]; /*  latitudes  of four bilinear corners */
  double src_lons[4]; /*  longitudes of four bilinear corners */
  double wgts[4][4];  /*  bicubic weights for four corners    */

  double plat, plon;             /*  lat/lon coords of destination point    */
  double findex = 0;
  extern long remap_max_iter;
  int remap_grid_type = src_grid->remap_grid_type;

  if ( cdoVerbose ) cdoPrint("Called %s()", __func__);

  progressInit();

  tgt_grid_size = tgt_grid->size;

  /* Compute mappings from source to target grid */

  if ( src_grid->rank != 2 )
    cdoAbort("Can not do bicubic interpolation when source grid rank != 2"); 

  /* Loop over destination grid */

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(ompNumThreads, cdoTimer, cdoVerbose, remap_grid_type, tgt_grid_size, src_grid, tgt_grid, rv, remap_max_iter, lwarn, findex) \
  private(dst_add, n, icount, iter, src_add, src_lats, src_lons, wgts, plat, plon, search_result) \
  schedule(dynamic,1)
#endif
  /* grid_loop1 */
  for ( dst_add = 0; dst_add < tgt_grid_size; ++dst_add )
    {
      int lprogress = 1;
#if defined(_OPENMP)
      if ( omp_get_thread_num() != 0 ) lprogress = 0;
#endif
#if defined(_OPENMP)
#pragma omp atomic
#endif
      findex++;
      if ( lprogress ) progressStatus(0, 1, findex/tgt_grid_size);

      if ( ! tgt_grid->mask[dst_add] ) continue;

      plat = tgt_grid->cell_center_lat[dst_add];
      plon = tgt_grid->cell_center_lon[dst_add];

      /* Find nearest square of grid points on source grid  */
      if ( remap_grid_type == REMAP_GRID_TYPE_REG2D )
	search_result = grid_search_reg2d(src_grid, src_add, src_lats, src_lons, 
					  plat, plon, src_grid->dims,
					  src_grid->reg2d_center_lat, src_grid->reg2d_center_lon);
      else
	search_result = grid_search(src_grid, src_add, src_lats, src_lons, 
				    plat, plon, src_grid->dims,
				    src_grid->cell_center_lat, src_grid->cell_center_lon,
				    src_grid->cell_bound_box, src_grid->bin_addr);

      /* Check to see if points are land points */
      if ( search_result > 0 )
	{
	  for ( n = 0; n < 4; ++n )
	    if ( ! src_grid->mask[src_add[n]] ) search_result = 0;
	}

      /* If point found, find local iw,jw coordinates for weights  */
      if ( search_result > 0 )
	{
	  double iw, jw;  /*  current guess for bilinear coordinate  */

          tgt_grid->cell_frac[dst_add] = ONE;

	  iter = find_ij_weights(plon, plat, src_lats, src_lons, &iw, &jw);

          if ( iter < remap_max_iter )
	    {
	      /* Successfully found iw,jw - compute weights */

	      wgts[0][0] = (ONE-jw*jw*(THREE-TWO*jw)) * (ONE-iw*iw*(THREE-TWO*iw));
	      wgts[0][1] = (ONE-jw*jw*(THREE-TWO*jw)) *      iw*iw*(THREE-TWO*iw);
	      wgts[0][2] =      jw*jw*(THREE-TWO*jw)  *      iw*iw*(THREE-TWO*iw);
	      wgts[0][3] =      jw*jw*(THREE-TWO*jw)  * (ONE-iw*iw*(THREE-TWO*iw));
	      wgts[1][0] = (ONE-jw*jw*(THREE-TWO*jw)) *      iw*(iw-ONE)*(iw-ONE);
	      wgts[1][1] = (ONE-jw*jw*(THREE-TWO*jw)) *      iw*iw*(iw-ONE);
	      wgts[1][2] =      jw*jw*(THREE-TWO*jw)  *      iw*iw*(iw-ONE);
	      wgts[1][3] =      jw*jw*(THREE-TWO*jw)  *      iw*(iw-ONE)*(iw-ONE);
	      wgts[2][0] =      jw*(jw-ONE)*(jw-ONE)  * (ONE-iw*iw*(THREE-TWO*iw));
	      wgts[2][1] =      jw*(jw-ONE)*(jw-ONE)  *      iw*iw*(THREE-TWO*iw);
	      wgts[2][2] =      jw*jw*(jw-ONE)        *      iw*iw*(THREE-TWO*iw);
	      wgts[2][3] =      jw*jw*(jw-ONE)        * (ONE-iw*iw*(THREE-TWO*iw));
	      wgts[3][0] =      iw*(iw-ONE)*(iw-ONE)  *      jw*(jw-ONE)*(jw-ONE);
              wgts[3][1] =      iw*iw*(iw-ONE)        *      jw*(jw-ONE)*(jw-ONE);
	      wgts[3][2] =      iw*iw*(iw-ONE)        *      jw*jw*(jw-ONE);
	      wgts[3][3] =      iw*(iw-ONE)*(iw-ONE)  *      jw*jw*(jw-ONE);

#if defined(_OPENMP)
#pragma omp critical
#endif
	      store_link_bicub(rv, dst_add, src_add, wgts);
	    }
          else
	    {
	      if ( cdoVerbose || lwarn )
		{
		  lwarn = FALSE;
		  // cdoWarning("Iteration for iw,jw exceed max iteration count of %d!", remap_max_iter);
		  cdoWarning("Bicubic interpolation failed for some grid points - used a distance-weighted average instead!");
		}

	      search_result = -1;
	    }
	}
	  
      /*
	Search for bicubic failed - use a distance-weighted average instead (this is typically near the pole)
      */
      if ( search_result < 0 )
	{
          icount = 0;
          for ( n = 0; n < 4; ++n )
	    {
	      if ( src_grid->mask[src_add[n]] )
		icount++;
	      else
		src_lats[n] = ZERO;
	    }

          if ( icount > 0 )
	    {
	      /* Renormalize weights */
	      double sum_wgts = 0.0; /* sum of weights for normalization */
	      /* 2012-05-08 Uwe Schulzweida: using absolute value of src_lats (bug fix) */
	      for ( n = 0; n < 4; ++n ) sum_wgts += fabs(src_lats[n]);
	      for ( n = 0; n < 4; ++n ) wgts[0][n] = fabs(src_lats[n])/sum_wgts;
	      for ( n = 0; n < 4; ++n ) wgts[1][n] = ZERO;
	      for ( n = 0; n < 4; ++n ) wgts[2][n] = ZERO;
	      for ( n = 0; n < 4; ++n ) wgts[3][n] = ZERO;

	      tgt_grid->cell_frac[dst_add] = ONE;

#if defined(_OPENMP)
#pragma omp critical
#endif
	      store_link_bicub(rv, dst_add, src_add, wgts);
	    }
        }
    } /* grid_loop1 */

} /* remap_bicub */

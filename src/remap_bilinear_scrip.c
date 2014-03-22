#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "remap.h"


/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/*                                                                         */
/*      BILINEAR INTERPOLATION                                             */
/*                                                                         */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

long find_element(double x, long nelem, const double *restrict array);
int rect_grid_search(long *ii, long *jj, double x, double y, long nxm, long nym, const double *restrict xm, const double *restrict ym);


int grid_search_reg2d_nn(long nx, long ny, int *restrict nbr_add, double *restrict nbr_dist, double plat, double plon,
			 const double *restrict src_center_lat, const double *restrict src_center_lon)
{
  int search_result = 0;
  long n, srch_add;
  long i;
  long ii, jj;
  long jjskip;
  double coslat, sinlat;
  double dist_min, distance; /* For computing dist-weighted avg */
  double *sincoslon;
  double coslat_dst = cos(plat);
  double sinlat_dst = sin(plat);
  double coslon_dst = cos(plon);
  double sinlon_dst = sin(plon);

  dist_min = BIGNUM;
  for ( n = 0; n < 4; ++n ) nbr_dist[n] = BIGNUM;  

  long jjf = 0, jjl = ny-1;
  if ( plon >= src_center_lon[0] && plon <= src_center_lon[nx-1] )
    {
      if ( src_center_lat[0] < src_center_lat[ny-1] )
	{
	  if ( plat <= src_center_lat[0] )
	    { jjf = 0; jjl = 1; }
	  else
	    { jjf = ny-2; jjl = ny-1; }
	}
      else
	{
	  if ( plat >= src_center_lat[0] )
	    { jjf = 0; jjl = 1; }
	  else
	    { jjf = ny-2; jjl = ny-1; }
	}
    }

  sincoslon = malloc(nx*sizeof(double));

  for ( ii = 0; ii < nx; ++ii )
    sincoslon[ii] = coslon_dst*cos(src_center_lon[ii]) + sinlon_dst*sin(src_center_lon[ii]);

  for ( jj = jjf; jj <= jjl; ++jj )
    {
      coslat = coslat_dst*cos(src_center_lat[jj]);
      sinlat = sinlat_dst*sin(src_center_lat[jj]);

      jjskip = jj > 1 && jj < (ny-2);

      for ( ii = 0; ii < nx; ++ii )
	{
	  if ( jjskip && ii > 1 && ii < (nx-2) ) continue;

	  srch_add = jj*nx + ii;

	  distance = acos(coslat*sincoslon[ii] + sinlat);

	  if ( distance < dist_min )
	    {
	      for ( n = 0; n < 4; ++n )
		{
		  if ( distance < nbr_dist[n] )
		    {
		      for ( i = 3; i > n; --i )
			{
			  nbr_add [i] = nbr_add [i-1];
			  nbr_dist[i] = nbr_dist[i-1];
			}
		      search_result = -1;
		      nbr_add [n] = srch_add;
		      nbr_dist[n] = distance;
		      dist_min = nbr_dist[3];
		      break;
		    }
		}
	    }
	}
    }

  free(sincoslon);

  for ( n = 0; n < 4; ++n ) nbr_dist[n] = ONE/(nbr_dist[n] + TINY);
  distance = 0.0;
  for ( n = 0; n < 4; ++n ) distance += nbr_dist[n];
  for ( n = 0; n < 4; ++n ) nbr_dist[n] /= distance;

  return (search_result);
}


int grid_search_reg2d(remapgrid_t *src_grid, int *restrict src_add, double *restrict src_lats, 
		      double *restrict src_lons,  double plat, double plon, const int *restrict src_grid_dims,
		      const double *restrict src_center_lat, const double *restrict src_center_lon)
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
  */
  /*  Local variables */
  int search_result = 0;
  int lfound;
  long n;
  long nx, nxm, ny;
  long ii, iix, jj;

  for ( n = 0; n < 4; ++n ) src_add[n] = 0;

  nx = src_grid_dims[0];
  ny = src_grid_dims[1];

  nxm = nx;
  if ( src_grid->is_cyclic ) nxm++;

  if ( /*plon < 0   &&*/ plon < src_center_lon[0]     ) plon += PI2;
  if ( /*plon > PI2 &&*/ plon > src_center_lon[nxm-1] ) plon -= PI2;

  lfound = rect_grid_search(&ii, &jj, plon, plat, nxm, ny, src_center_lon, src_center_lat);

  if ( lfound )
    {
      iix = ii;
      if ( src_grid->is_cyclic && iix == (nxm-1) ) iix = 0;
      src_add[0] = (jj-1)*nx+(ii-1);
      src_add[1] = (jj-1)*nx+(iix);
      src_add[2] = (jj)*nx+(iix);
      src_add[3] = (jj)*nx+(ii-1);

      src_lons[0] = src_center_lon[ii-1];
      src_lons[1] = src_center_lon[iix];
      /* For consistency, we must make sure all lons are in same 2pi interval */
      if ( src_lons[0] > PI2 ) src_lons[0] -= PI2;
      if ( src_lons[0] < 0   ) src_lons[0] += PI2;
      if ( src_lons[1] > PI2 ) src_lons[1] -= PI2;
      if ( src_lons[1] < 0   ) src_lons[1] += PI2;
      src_lons[2] = src_lons[1];
      src_lons[3] = src_lons[0];

      src_lats[0] = src_center_lat[jj-1];
      src_lats[1] = src_lats[0];
      src_lats[2] = src_center_lat[jj];
      src_lats[3] = src_lats[2];

      search_result = 1;
      
      return (search_result);
    }

  /*
    If no cell found, point is likely either in a box that straddles either pole or is outside 
    the grid. Fall back to a distance-weighted average of the four closest points.
    Go ahead and compute weights here, but store in src_lats and return -add to prevent the 
    parent routine from computing bilinear weights.
  */
  if ( !src_grid->lextrapolate ) return (search_result);

  /*
    printf("Could not find location for %g %g\n", plat*RAD2DEG, plon*RAD2DEG);
    printf("Using nearest-neighbor average for this point\n");
  */
  search_result = grid_search_reg2d_nn(nx, ny, src_add, src_lats, plat, plon, src_center_lat, src_center_lon);

  return (search_result);
}  /* grid_search_reg2d */


static
int grid_search_nn(long min_add, long max_add, int *restrict nbr_add, double *restrict nbr_dist, 
		   double plat, double plon,
		   const double *restrict src_center_lat, const double *restrict src_center_lon)
{
  int search_result = 0;
  long n, srch_add;
  long i;
  double dist_min, distance; /* For computing dist-weighted avg */
  double coslat_dst = cos(plat);
  double sinlat_dst = sin(plat);
  double coslon_dst = cos(plon);
  double sinlon_dst = sin(plon);

  dist_min = BIGNUM;
  for ( n = 0; n < 4; ++n ) nbr_dist[n] = BIGNUM;
  for ( srch_add = min_add; srch_add <= max_add; ++srch_add )
    {
      distance = acos(coslat_dst*cos(src_center_lat[srch_add])*
		     (coslon_dst*cos(src_center_lon[srch_add]) +
                      sinlon_dst*sin(src_center_lon[srch_add]))+
		      sinlat_dst*sin(src_center_lat[srch_add]));

      if ( distance < dist_min )
	{
          for ( n = 0; n < 4; ++n )
	    {
	      if ( distance < nbr_dist[n] )
		{
		  for ( i = 3; i > n; --i )
		    {
		      nbr_add [i] = nbr_add [i-1];
		      nbr_dist[i] = nbr_dist[i-1];
		    }
		  search_result = -1;
		  nbr_add [n] = srch_add;
		  nbr_dist[n] = distance;
		  dist_min = nbr_dist[3];
		  break;
		}
	    }
        }
    }

  for ( n = 0; n < 4; ++n ) nbr_dist[n] = ONE/(nbr_dist[n] + TINY);
  distance = 0.0;
  for ( n = 0; n < 4; ++n ) distance += nbr_dist[n];
  for ( n = 0; n < 4; ++n ) nbr_dist[n] /= distance;

  return (search_result);
}


int grid_search(remapgrid_t *src_grid, int *restrict src_add, double *restrict src_lats, 
		double *restrict src_lons,  double plat, double plon, const int *restrict src_grid_dims,
		const double *restrict src_center_lat, const double *restrict src_center_lon,
		const restr_t *restrict src_grid_bound_box, const int *restrict src_bin_add)
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

    restr_t src_grid_bound_box[][4] ! bound box for source grid

    int src_bin_add[][2]           ! latitude bins for restricting
  */
  /*  Local variables */
  long n, n2, next_n, srch_add, srch_add4;    /* dummy indices                    */
  long nx, ny;                                /* dimensions of src grid           */
  long min_add, max_add;                      /* addresses for restricting search */
  long i, j, jp1, ip1, n_add, e_add, ne_add;  /* addresses                        */
  long nbins;
  /* Vectors for cross-product check */
  double vec1_lat, vec1_lon;
  double vec2_lat, vec2_lon, cross_product;
  int scross[4], scross_last = 0;
  int search_result = 0;
  restr_t rlat, rlon;
  restr_t *bin_lats = src_grid->bin_lats;

  nbins = src_grid->num_srch_bins;

  rlat = RESTR_SCALE(plat);
  rlon = RESTR_SCALE(plon);

  /* restrict search first using bins */

  for ( n = 0; n < 4; ++n ) src_add[n] = 0;

  min_add = src_grid->size-1;
  max_add = 0;

  for ( n = 0; n < nbins; ++n )
    {
      n2 = n<<1;
      if ( rlat >= bin_lats[n2] && rlat <= bin_lats[n2+1] )
	{
	  if ( src_bin_add[n2  ] < min_add ) min_add = src_bin_add[n2  ];
	  if ( src_bin_add[n2+1] > max_add ) max_add = src_bin_add[n2+1];
	}
    }
 
  /* Now perform a more detailed search */

  nx = src_grid_dims[0];
  ny = src_grid_dims[1];

  /* srch_loop */
  for ( srch_add = min_add; srch_add <= max_add; ++srch_add )
    {
      srch_add4 = srch_add<<2;
      /* First check bounding box */
      if ( rlon >= src_grid_bound_box[srch_add4+2] &&
	   rlon <= src_grid_bound_box[srch_add4+3] &&
	   rlat >= src_grid_bound_box[srch_add4  ] &&
	   rlat <= src_grid_bound_box[srch_add4+1])
	{
	  /* We are within bounding box so get really serious */

          /* Determine neighbor addresses */
          j = srch_add/nx;
          i = srch_add - j*nx;

          if ( i < (nx-1) )
            ip1 = i + 1;
          else
	    {
	      /* 2009-01-09 Uwe Schulzweida: bug fix */
	      if ( src_grid->is_cyclic )
		ip1 = 0;
	      else
		ip1 = i;
	    }

          if ( j < (ny-1) )
            jp1 = j + 1;
          else
	    {
	      /* 2008-12-17 Uwe Schulzweida: latitute cyclic ??? (bug fix) */
	      jp1 = j;
	    }

          n_add  = jp1*nx + i;
          e_add  = j  *nx + ip1;
	  ne_add = jp1*nx + ip1;

          src_lons[0] = src_center_lon[srch_add];
          src_lons[1] = src_center_lon[e_add];
          src_lons[2] = src_center_lon[ne_add];
          src_lons[3] = src_center_lon[n_add];

          src_lats[0] = src_center_lat[srch_add];
          src_lats[1] = src_center_lat[e_add];
          src_lats[2] = src_center_lat[ne_add];
          src_lats[3] = src_center_lat[n_add];

	  /* For consistency, we must make sure all lons are in same 2pi interval */

          vec1_lon = src_lons[0] - plon;
          if      ( vec1_lon >  PI ) src_lons[0] -= PI2;
          else if ( vec1_lon < -PI ) src_lons[0] += PI2;

          for ( n = 1; n < 4; ++n )
	    {
	      vec1_lon = src_lons[n] - src_lons[0];
	      if      ( vec1_lon >  PI ) src_lons[n] -= PI2;
	      else if ( vec1_lon < -PI ) src_lons[n] += PI2;
	    }

          /* corner_loop */
          for ( n = 0; n < 4; ++n )
	    {
	      next_n = (n+1)%4;

	      /*
		Here we take the cross product of the vector making 
		up each box side with the vector formed by the vertex
		and search point.  If all the cross products are 
		positive, the point is contained in the box.
	      */
	      vec1_lat = src_lats[next_n] - src_lats[n];
	      vec1_lon = src_lons[next_n] - src_lons[n];
	      vec2_lat = plat - src_lats[n];
	      vec2_lon = plon - src_lons[n];

	      /* Check for 0,2pi crossings */

	      if      ( vec1_lon >  THREE*PIH ) vec1_lon -= PI2;
	      else if ( vec1_lon < -THREE*PIH ) vec1_lon += PI2;

	      if      ( vec2_lon >  THREE*PIH ) vec2_lon -= PI2;
	      else if ( vec2_lon < -THREE*PIH ) vec2_lon += PI2;

	      cross_product = vec1_lon*vec2_lat - vec2_lon*vec1_lat;

	      /* If cross product is less than ZERO, this cell doesn't work    */
	      /* 2008-10-16 Uwe Schulzweida: bug fix for cross_product eq zero */

	      scross[n] = cross_product < 0 ? -1 : cross_product > 0 ? 1 : 0;

	      if ( n == 0 ) scross_last = scross[n];

	      if ( (scross[n] < 0 && scross_last > 0) || (scross[n] > 0 && scross_last < 0) ) break;

	      scross_last = scross[n];
	    } /* corner_loop */

	  if ( n >= 4 )
	    {
	      n = 0;
	      if      ( scross[0]>=0 && scross[1]>=0 && scross[2]>=0 && scross[3]>=0 ) n = 4;
	      else if ( scross[0]<=0 && scross[1]<=0 && scross[2]<=0 && scross[3]<=0 ) n = 4;
	    }

	  /* If cross products all same sign, we found the location */
          if ( n >= 4 )
	    {
	      src_add[0] = srch_add;
	      src_add[1] = e_add;
	      src_add[2] = ne_add;
	      src_add[3] = n_add;

	      search_result = 1;

	      return (search_result);
	    }

	  /* Otherwise move on to next cell */

        } /* Bounding box check */
    } /* srch_loop */

  /*
    If no cell found, point is likely either in a box that straddles either pole or is outside 
    the grid. Fall back to a distance-weighted average of the four closest points.
    Go ahead and compute weights here, but store in src_lats and return -add to prevent the 
    parent routine from computing bilinear weights.
  */
  if ( !src_grid->lextrapolate ) return (search_result);

  /*
    printf("Could not find location for %g %g\n", plat*RAD2DEG, plon*RAD2DEG);
    printf("Using nearest-neighbor average for this point\n");
  */
  search_result = grid_search_nn(min_add, max_add, src_add, src_lats, plat, plon, src_center_lat, src_center_lon);

  return (search_result);
}  /* grid_search */

/*
  This routine stores the address and weight for four links associated with one destination
  point in the appropriate address and weight arrays and resizes those arrays if necessary.
*/
void store_link_bilin(remapvars_t *rv, int dst_add, int *restrict src_add, double *restrict weights)
{
  /*
    Input variables:
    int dst_add       ! address on destination grid
    int src_add[4]    ! addresses on source grid
    double weights[4] ! array of remapping weights for these links
  */
  long n;
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
	    double dtmp;
	    dtmp = weights[2];
	    weights[2] = weights[3];
	    weights[3] = dtmp;
	  }
	}
    }

  for ( n = 0; n < 4; ++n )
    {
      rv->src_grid_add[nlink+n] = src_add[n];
      rv->tgt_grid_add[nlink+n] = dst_add;
      rv->wts         [nlink+n] = weights[n];
    }

} /* store_link_bilin */


long find_ij_weights(double plon, double plat, double* restrict src_lats, double* restrict src_lons, double *ig, double *jg)
{
  long iter;                     /*  iteration counters   */
  double iguess, jguess;         /*  current guess for bilinear coordinate  */
  double deli, delj;             /*  corrections to iw,jw                   */
  double dth1, dth2, dth3;       /*  some latitude  differences             */
  double dph1, dph2, dph3;       /*  some longitude differences             */
  double dthp, dphp;             /*  difference between point and sw corner */
  double mat1, mat2, mat3, mat4; /*  matrix elements                        */
  double determinant;            /*  matrix determinant                     */
  double converge = 1.e-10;      /* Convergence criterion                   */
  extern long remap_max_iter;

  /* Iterate to find iw,jw for bilinear approximation  */

  dth1 = src_lats[1] - src_lats[0];
  dth2 = src_lats[3] - src_lats[0];
  dth3 = src_lats[2] - src_lats[1] - dth2;

  dph1 = src_lons[1] - src_lons[0];
  dph2 = src_lons[3] - src_lons[0];
  dph3 = src_lons[2] - src_lons[1];

  if ( dph1 >  THREE*PIH ) dph1 -= PI2;
  if ( dph2 >  THREE*PIH ) dph2 -= PI2;
  if ( dph3 >  THREE*PIH ) dph3 -= PI2;
  if ( dph1 < -THREE*PIH ) dph1 += PI2;
  if ( dph2 < -THREE*PIH ) dph2 += PI2;
  if ( dph3 < -THREE*PIH ) dph3 += PI2;

  dph3 = dph3 - dph2;

  iguess = HALF;
  jguess = HALF;

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

  return (iter);
}

/*
  -----------------------------------------------------------------------

  This routine computes the weights for a bilinear interpolation.

  -----------------------------------------------------------------------
*/
void remap_bilinear(remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv)
{
  /*   Local variables */
  int  lwarn = TRUE;
  int  search_result;
  long tgt_grid_size;
  long dst_add;                  /*  destination addresss */
  long n, icount;
  long iter;                     /*  iteration counters   */

  int src_add[4];                /*  address for the four source points     */

  double src_lats[4];            /*  latitudes  of four bilinear corners    */
  double src_lons[4];            /*  longitudes of four bilinear corners    */
  double wgts[4];                /*  bilinear weights for four corners      */

  double plat, plon;             /*  lat/lon coords of destination point    */
  double findex = 0;
  extern int timer_remap_bil;
  extern long remap_max_iter;
  int remap_grid_type = src_grid->remap_grid_type;

  if ( cdoVerbose ) cdoPrint("Called %s()", __func__);

  if ( cdoTimer ) timer_start(timer_remap_bil);

  progressInit();

  tgt_grid_size = tgt_grid->size;

  /* Compute mappings from source to target grid */

  if ( src_grid->rank != 2 )
    cdoAbort("Can not do bilinear interpolation when source grid rank != 2"); 

  /* Loop over destination grid */

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(ompNumThreads, cdoTimer, cdoVerbose, remap_grid_type, tgt_grid_size, src_grid, tgt_grid, rv, remap_max_iter, lwarn, findex) \
  private(dst_add, n, icount, iter, src_add, src_lats, src_lons, wgts, plat, plon, search_result)    \
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

	      wgts[0] = (ONE-iw) * (ONE-jw);
	      wgts[1] =      iw  * (ONE-jw);
	      wgts[2] =      iw  *      jw;
	      wgts[3] = (ONE-iw) *      jw;

#if defined(_OPENMP)
#pragma omp critical
#endif
	      store_link_bilin(rv, dst_add, src_add, wgts);
	    }
          else
	    {
	      if ( cdoVerbose )
		{
		  cdoPrint("Point coords: %g %g", plat, plon);
		  cdoPrint("Src grid lats: %g %g %g %g", src_lats[0], src_lats[1], src_lats[2], src_lats[3]);
		  cdoPrint("Src grid lons: %g %g %g %g", src_lons[0], src_lons[1], src_lons[2], src_lons[3]);
		  cdoPrint("Src grid addresses: %d %d %d %d", src_add[0], src_add[1], src_add[2], src_add[3]);
		  cdoPrint("Src grid lats: %g %g %g %g",
			   src_grid->cell_center_lat[src_add[0]], src_grid->cell_center_lat[src_add[1]],
			   src_grid->cell_center_lat[src_add[2]], src_grid->cell_center_lat[src_add[3]]);
		  cdoPrint("Src grid lons: %g %g %g %g",
			   src_grid->cell_center_lon[src_add[0]], src_grid->cell_center_lon[src_add[1]],
			   src_grid->cell_center_lon[src_add[2]], src_grid->cell_center_lon[src_add[3]]);
		  cdoPrint("Current iw,jw : %g %g", iw, jw);
		}

	      if ( cdoVerbose || lwarn )
		{
		  lwarn = FALSE;
		  //  cdoWarning("Iteration for iw,jw exceed max iteration count of %d!", remap_max_iter);
		  cdoWarning("Bilinear interpolation failed for some grid points - used a distance-weighted average instead!");
		}

	      search_result = -1;
	    }
	}

      /*
	Search for bilinear failed - use a distance-weighted average instead (this is typically near the pole)
	Distance was stored in src_lats!
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
	      for ( n = 0; n < 4; ++n ) wgts[n] = fabs(src_lats[n])/sum_wgts;

	      tgt_grid->cell_frac[dst_add] = ONE;

#if defined(_OPENMP)
#pragma omp critical
#endif
	      store_link_bilin(rv, dst_add, src_add, wgts);
	    }
        }
    } /* grid_loop1 */

  if ( cdoTimer ) timer_stop(timer_remap_bil);
} /* remap_bilinear */

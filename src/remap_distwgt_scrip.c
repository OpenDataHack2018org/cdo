#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "remap.h"

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/*                                                                         */
/*      INTERPOLATION USING A DISTANCE-WEIGHTED AVERAGE                    */
/*                                                                         */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

static
void get_restrict_add(remapgrid_t *src_grid, double plat, const int *restrict src_bin_add, long *minadd, long *maxadd)
{
  long n, n2;
  long min_add = 0, max_add = 0, nm1, np1;
  long nbins;
  restr_t rlat;
  restr_t *bin_lats = src_grid->bin_lats;

  nbins = src_grid->num_srch_bins;

  rlat = RESTR_SCALE(plat);

  for ( n = 0; n < nbins; ++n )
    {
      n2 = n<<1;
      if ( rlat >= bin_lats[n2  ] && rlat <= bin_lats[n2+1] )
	{
	  min_add = src_bin_add[n2  ];
	  max_add = src_bin_add[n2+1];

	  nm1 = MAX(n-1, 0);
	  np1 = MIN(n+1, nbins-1);

	  min_add = MIN(min_add, src_bin_add[2*nm1  ]);
	  max_add = MAX(max_add, src_bin_add[2*nm1+1]);
	  min_add = MIN(min_add, src_bin_add[2*np1  ]);
	  max_add = MAX(max_add, src_bin_add[2*np1+1]);
	}
    }

  *minadd = min_add;
  *maxadd = max_add;
  /*
  if ( cdoVerbose )
    printf("plon %g plat %g min_add %ld max_add %ld diff %ld\n",
	   plon, plat, min_add, max_add, max_add-min_add);
  */
}

/*
   This routine finds the closest num_neighbor points to a search 
   point and computes a distance to each of the neighbors.
*/
static
void grid_search_nbr_reg2d(int num_neighbors, remapgrid_t *src_grid, int *restrict nbr_add, double *restrict nbr_dist, 
			   double plat, double plon, const int *restrict src_grid_dims,
			   double coslat_dst, double coslon_dst, double sinlat_dst, double sinlon_dst,
			   const double *restrict sinlat, const double *restrict coslat,
			   const double *restrict sinlon, const double *restrict coslon,
			   const double *restrict src_center_lat, const double *restrict src_center_lon)
{
  /*
    Output variables:

    int nbr_add[num_neighbors]     ! address of each of the closest points
    double nbr_dist[num_neighbors] ! distance to each of the closest points

    Input variables:

    double plat,         ! latitude  of the search point
    double plon,         ! longitude of the search point
  */
  /*  Local variables */
  int lfound;
  long n, nadd, nchk;
  long nx, nxm, ny;
  long ii, jj;
  long i, j, ix;
  int src_add[25];
  long num_add = 0;
  double distance;   //  Angular distance
  /*
  double coslat_dst = cos(plat);  // cos(lat)  of the search point
  double coslon_dst = cos(plon);  // cos(lon)  of the search point
  double sinlat_dst = sin(plat);  // sin(lat)  of the search point
  double sinlon_dst = sin(plon);  // sin(lon)  of the search point
  */
  nx = src_grid_dims[0];
  ny = src_grid_dims[1];

  nxm = nx;
  if ( src_grid->is_cyclic ) nxm++;

  if ( plon < src_center_lon[0]     ) plon += PI2;
  if ( plon > src_center_lon[nxm-1] ) plon -= PI2;

  lfound = rect_grid_search(&ii, &jj, plon, plat, nxm, ny, src_center_lon, src_center_lat);

  if ( lfound )
    {
      if ( src_grid->is_cyclic && ii == (nxm-1) ) ii = 0;

      for ( j = (jj-2); j <= (jj+2); ++j )
	for ( i = (ii-2); i <= (ii+2); ++i )
	  {
	    ix = i;
	    
	    if ( src_grid->is_cyclic )
	      {
		if ( ix <   0 ) ix += nx;
		if ( ix >= nx ) ix -= nx;
	      }

	    if ( ix >= 0 && ix < nx && j >= 0 && j < ny )
	      src_add[num_add++] = j*nx+ix;
	  }
      /*
      num_add = 0;

      for ( j = (jj-1); j <= jj; ++j )
	for ( i = (ii-1); i <= ii; ++i )
	  {
	    ix = i;
	    if ( src_grid->is_cyclic && ix == (nxm-1) ) ix = 0;

	    src_add[num_add++] = j*nx+ix;
	  }
      */
    }

  /* Initialize distance and address arrays */
  for ( n = 0; n < num_neighbors; ++n )
    {
      nbr_add[n]  = -1;
      nbr_dist[n] = BIGNUM;
    }

  if ( lfound )
    {
      long ix, iy;

      for ( long na = 0; na < num_add; ++na )
	{
	  nadd = src_add[na];

	  iy = nadd/nx;
	  ix = nadd - iy*nx;

	  /* Find distance to this point */
	  distance =  sinlat_dst*sinlat[iy] + coslat_dst*coslat[iy]*
	             (coslon_dst*coslon[ix] + sinlon_dst*sinlon[ix]);
	  /*
	  distance =  sinlat_dst*sinlat[nadd] + coslat_dst*coslat[nadd]*
	             (coslon_dst*coslon[nadd] + sinlon_dst*sinlon[nadd]);
	  */
	  /* 2008-07-30 Uwe Schulzweida: check that distance is inside the range of -1 to 1,
	                                 otherwise the result of acos(distance) is NaN */
	  if ( distance >  1 ) distance =  1;
	  if ( distance < -1 ) distance = -1;
	  distance = acos(distance);

	  /* Uwe Schulzweida: if distance is zero, set to small number */
	  if ( IS_EQUAL(distance, 0) ) distance = TINY;

	  /* Store the address and distance if this is one of the smallest four so far */
	  for ( nchk = 0; nchk < num_neighbors; ++nchk )
	    {
	      if ( distance < nbr_dist[nchk] )
		{
		  for ( n = num_neighbors-1; n > nchk; --n )
		    {
		      nbr_add[n]  = nbr_add[n-1];
		      nbr_dist[n] = nbr_dist[n-1];
		    }
		  nbr_add[nchk]  = nadd;
		  nbr_dist[nchk] = distance;
		  break;
		}
	      else if ( num_neighbors == 1 && distance <= nbr_dist[0] && nadd < nbr_add[0] )
		{
		  nbr_add[0]  = nadd;
		  nbr_dist[0] = distance;
		}
	    }
	}
    }
  else if ( src_grid->lextrapolate )
    {
      int search_result;
      search_result = grid_search_reg2d_nn(nx, ny, nbr_add, nbr_dist, plat, plon, src_center_lat, src_center_lon);
      
      if ( search_result >= 0 )
	for ( n = 0; n < 4; ++n ) nbr_add[n] = -1;
    }
}  /*  grid_search_nbr_reg2d  */

static
void grid_search_nbr(int num_neighbors, remapgrid_t *src_grid, int *restrict nbr_add, double *restrict nbr_dist, 
		     double plat, double plon, const int *restrict src_bin_add,
		     double coslat_dst, double coslon_dst, double sinlat_dst, double sinlon_dst,
		     const double *restrict sinlat, const double *restrict coslat,
		     const double *restrict sinlon, const double *restrict coslon)
{
  /*
    Output variables:

    int nbr_add[num_neighbors]     ! address of each of the closest points
    double nbr_dist[num_neighbors] ! distance to each of the closest points

    Input variables:

    int src_bin_add[][2]  ! search bins for restricting search

    double plat,         ! latitude  of the search point
    double plon,         ! longitude of the search point
  */
  /*  Local variables */
  long n, nadd, nchk;
  long min_add, max_add;
  double distance;     /* Angular distance */
  /* result changed a little on a few points with high resolution grid
  double xcoslat_dst = cos(plat);  // cos(lat)  of the search point
  double xcoslon_dst = cos(plon);  // cos(lon)  of the search point
  double xsinlat_dst = sin(plat);  // sin(lat)  of the search point
  double xsinlon_dst = sin(plon);  // sin(lon)  of the search point
  */
  /* Loop over source grid and find nearest neighbors                         */
  /* restrict the search using search bins expand the bins to catch neighbors */

  get_restrict_add(src_grid, plat, src_bin_add, &min_add, &max_add);

  /* Initialize distance and address arrays */
  for ( n = 0; n < num_neighbors; ++n )
    {
      nbr_add[n]  = -1;
      nbr_dist[n] = BIGNUM;
    }

  for ( nadd = min_add; nadd <= max_add; ++nadd )
    {
      /* Find distance to this point */
      distance =  sinlat_dst*sinlat[nadd] + coslat_dst*coslat[nadd]*
	         (coslon_dst*coslon[nadd] + sinlon_dst*sinlon[nadd]);
      /* 2008-07-30 Uwe Schulzweida: check that distance is inside the range of -1 to 1,
                                     otherwise the result of acos(distance) is NaN */
      if ( distance >  1 ) distance =  1;
      if ( distance < -1 ) distance = -1;
      distance = acos(distance);

      /* Uwe Schulzweida: if distance is zero, set to small number */
      if ( IS_EQUAL(distance, 0) ) distance = TINY;

      /* Store the address and distance if this is one of the smallest four so far */
      for ( nchk = 0; nchk < num_neighbors; ++nchk )
	{
          if ( distance < nbr_dist[nchk] )
	    {
	      for ( n = num_neighbors-1; n > nchk; --n )
		{
		  nbr_add[n]  = nbr_add[n-1];
		  nbr_dist[n] = nbr_dist[n-1];
		}
	      nbr_add[nchk]  = nadd;
	      nbr_dist[nchk] = distance;
	      break;
	    }
        }
    }

}  /*  grid_search_nbr  */

/*
  This routine stores the address and weight for this link in the appropriate 
  address and weight arrays and resizes those arrays if necessary.
*/
static
void store_link_nbr(remapvars_t *rv, int add1, int add2, double weights)
{
  /*
    Input variables:
    int  add1         ! address on source grid
    int  add2         ! address on target grid
    double weights    ! remapping weight for this link
  */
  long nlink;

  /*
     Increment number of links and check to see if remap arrays need
     to be increased to accomodate the new link. Then store the link.
  */
  nlink = rv->num_links;
  rv->num_links++;

  if ( rv->num_links >= rv->max_links ) 
    resize_remap_vars(rv, rv->resize_increment);

  rv->src_grid_add[nlink] = add1;
  rv->tgt_grid_add[nlink] = add2;
  rv->wts[nlink]          = weights;

} /* store_link_nbr */

/*
  -----------------------------------------------------------------------

   This routine computes the inverse-distance weights for a
   nearest-neighbor interpolation.

  -----------------------------------------------------------------------
*/
void scrip_remap_weights_distwgt(int num_neighbors, remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv)
{
  /*  Local variables */

  long src_grid_size;
  long tgt_grid_size;
  long n;
  long dst_add;                   /* destination address                     */
  int nbr_mask[num_neighbors];    /* mask at nearest neighbors               */
  int nbr_add[num_neighbors];     /* source address at nearest neighbors     */
  double nbr_dist[num_neighbors]; /* angular distance four nearest neighbors */
  double dist_tot;         /* sum of neighbor distances (for normalizing) */
  double coslat_dst;       /* cos(lat) of destination grid point */
  double coslon_dst;       /* cos(lon) of destination grid point */
  double sinlat_dst;       /* sin(lat) of destination grid point */
  double sinlon_dst;       /* sin(lon) of destination grid point */
  double *coslat, *sinlat; /* cosine, sine of grid lats (for distance)    */
  double *coslon, *sinlon; /* cosine, sine of grid lons (for distance)    */
  double wgtstmp;          /* hold the link weight                        */
  double plat, plon;             /*  lat/lon coords of destination point    */
  double findex = 0;
  int remap_grid_type = src_grid->remap_grid_type;

  if ( cdoVerbose ) cdoPrint("Called %s()", __func__);

  progressInit();

  /* Compute mappings from source to target grid */

  src_grid_size = src_grid->size;
  tgt_grid_size = tgt_grid->size;

  /* Compute cos, sin of lat/lon on source grid for distance calculations */

  if ( remap_grid_type == REMAP_GRID_TYPE_REG2D )
    {
      long nx = src_grid->dims[0];
      long ny = src_grid->dims[1];

      coslat = malloc(ny*sizeof(double));
      coslon = malloc(nx*sizeof(double));
      sinlat = malloc(ny*sizeof(double));
      sinlon = malloc(nx*sizeof(double));

      for ( n = 0; n < nx; ++n )
	{
	  double rlon = src_grid->reg2d_center_lon[n];
	  if ( rlon > PI2  ) rlon -= PI2;
	  if ( rlon < ZERO ) rlon += PI2;
	  coslon[n] = cos(rlon);
	  sinlon[n] = sin(rlon);
	}
      for ( n = 0; n < ny; ++n )
	{
	  coslat[n] = cos(src_grid->reg2d_center_lat[n]);
	  sinlat[n] = sin(src_grid->reg2d_center_lat[n]);
	}
    }
  else
    {
      coslat = malloc(src_grid_size*sizeof(double));
      coslon = malloc(src_grid_size*sizeof(double));
      sinlat = malloc(src_grid_size*sizeof(double));
      sinlon = malloc(src_grid_size*sizeof(double));

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(src_grid, src_grid_size, coslat, coslon, sinlat, sinlon)
#endif
      for ( n = 0; n < src_grid_size; ++n )
	{
	  coslat[n] = cos(src_grid->cell_center_lat[n]);
	  coslon[n] = cos(src_grid->cell_center_lon[n]);
	  sinlat[n] = sin(src_grid->cell_center_lat[n]);
	  sinlon[n] = sin(src_grid->cell_center_lon[n]);
	}
    }

  /* Loop over destination grid  */
#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(ompNumThreads, cdoTimer, num_neighbors, remap_grid_type, src_grid, tgt_grid, rv, tgt_grid_size, coslat, coslon, sinlat, sinlon, findex) \
  private(dst_add, n, coslat_dst, coslon_dst, sinlat_dst, sinlon_dst, dist_tot, nbr_add, nbr_dist, nbr_mask, wgtstmp, plat, plon) \
  schedule(dynamic,1)
#endif
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

      coslat_dst = cos(plat);
      coslon_dst = cos(plon);
      sinlat_dst = sin(plat);
      sinlon_dst = sin(plon);

      /* Find nearest grid points on source grid and distances to each point */
      if ( remap_grid_type == REMAP_GRID_TYPE_REG2D )
	grid_search_nbr_reg2d(num_neighbors, src_grid, nbr_add, nbr_dist, 
			      plat, plon, src_grid->dims,
			      coslat_dst, coslon_dst, sinlat_dst, sinlon_dst,
			      sinlat, coslat, sinlon, coslon,
			      src_grid->reg2d_center_lat, src_grid->reg2d_center_lon);
      else
	grid_search_nbr(num_neighbors, src_grid, nbr_add, nbr_dist, 
			plat, plon, src_grid->bin_addr,
			coslat_dst, coslon_dst, sinlat_dst, sinlon_dst,
			sinlat, coslat, sinlon, coslon);

      /* Compute weights based on inverse distance if mask is false, eliminate those points */

      dist_tot = ZERO;
      for ( n = 0; n < num_neighbors; ++n )
	{
	  // printf("dst_add %ld %ld %d %g\n", dst_add, n, nbr_add[n], nbr_dist[n]);
	  nbr_mask[n] = FALSE;

	  /* Uwe Schulzweida: check if nbr_add is valid */
	  if ( nbr_add[n] >= 0 )
	    if ( src_grid->mask[nbr_add[n]] )
	      {
		nbr_dist[n] = ONE/nbr_dist[n];
		dist_tot = dist_tot + nbr_dist[n];
		nbr_mask[n] = TRUE;
	      }
	}

      /* Normalize weights and store the link */

      for ( n = 0; n < num_neighbors; ++n )
	{
          if ( nbr_mask[n] )
	    {
	      wgtstmp = nbr_dist[n]/dist_tot;

	      tgt_grid->cell_frac[dst_add] = ONE;
#if defined(_OPENMP)
#pragma omp critical
#endif
	      store_link_nbr(rv, nbr_add[n], dst_add, wgtstmp);
	    }
	}
    } /* for ( dst_add = 0; dst_add < tgt_grid_size; ++dst_add ) */

  free(coslat);
  free(coslon);
  free(sinlat);
  free(sinlon);

}  /* remap_distwgt */
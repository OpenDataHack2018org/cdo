#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "remap.h"
#include "remap_store_link.h"
#include "grid_search.h"


/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/*                                                                         */
/*      INTERPOLATION USING A DISTANCE-WEIGHTED AVERAGE                    */
/*                                                                         */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

static
void nbr_store_distance(int nadd, double distance, int num_neighbors, int *restrict nbr_add, double *restrict nbr_dist)
{
  if ( num_neighbors == 1 )
    {
      // if ( (distance+1.e-10) < nbr_dist[0] || ((fabs(distance-nbr_dist[0]) < 1.e-10) && nadd < nbr_add[0]) )
      if ( distance < nbr_dist[0] || (distance <= nbr_dist[0] && nadd < nbr_add[0]) )
	{
	  nbr_add[0]  = nadd;
	  nbr_dist[0] = distance;
	}
    }
  else
    {
      int n, nchk;
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
}

static
void nbr_check_distance(unsigned num_neighbors, const int *restrict nbr_add, double *restrict nbr_dist)
{
  double distance;

  /* Uwe Schulzweida: if distance is zero, set to small number */
  for ( unsigned nchk = 0; nchk < num_neighbors; ++nchk )
    {
      if ( nbr_add[nchk] >= 0 )
	{
	  distance = nbr_dist[nchk];
	  if ( IS_EQUAL(distance, 0.) ) distance = TINY;
	  nbr_dist[nchk] = distance;
	}
    }
}

static
double nbr_compute_weights(unsigned num_neighbors, const int *restrict src_grid_mask, int *restrict nbr_mask, const int *restrict nbr_add, double *restrict nbr_dist)
{
  /* Compute weights based on inverse distance if mask is false, eliminate those points */

  double dist_tot = 0.; /* sum of neighbor distances (for normalizing) */

  for ( unsigned n = 0; n < num_neighbors; ++n )
    {
      // printf("tgt_cell_add %ld %ld %d %g\n", tgt_cell_add, n, nbr_add[n], nbr_dist[n]);
      nbr_mask[n] = FALSE;

      /* Uwe Schulzweida: check if nbr_add is valid */
      if ( nbr_add[n] >= 0 )
        if ( src_grid_mask[nbr_add[n]] )
          {
            nbr_dist[n] = ONE/nbr_dist[n];
            dist_tot = dist_tot + nbr_dist[n];
            nbr_mask[n] = TRUE;
          }
    }

  return dist_tot;
}

static
unsigned nbr_normalize_weights(unsigned num_neighbors, double dist_tot, const int *restrict nbr_mask, int *restrict nbr_add, double *restrict nbr_dist)
{
  /* Normalize weights and store the link */

  unsigned nadds = 0;

  for ( unsigned n = 0; n < num_neighbors; ++n )
    {
      if ( nbr_mask[n] )
        {
          nbr_dist[nadds] = nbr_dist[n]/dist_tot;
          nbr_add[nadds]  = nbr_add[n];
          nadds++;
        }
    }

  return nadds;
}

static
double get_search_radius(void)
{
  extern double remap_search_radius;

  double search_radius = remap_search_radius;

  if ( search_radius <    0. ) search_radius = 0.;
  if ( search_radius >  180. ) search_radius = 180.;

  search_radius = search_radius*DEG2RAD;

  return search_radius;
}

/*
   This routine finds the closest num_neighbor points to a search 
   point and computes a distance to each of the neighbors.
*/
static
void grid_search_nbr_reg2d(struct gridsearch *gs, int num_neighbors, remapgrid_t *src_grid, int *restrict nbr_add, double *restrict nbr_dist, 
			   double plat, double plon, const int *restrict src_grid_dims)
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
  int n, nadd;
  long ii, jj;
  int i, j, ix;
  int src_add[25];
  int num_add = 0;
  double distance;   //  Angular distance
  double cos_search_radius = cos(get_search_radius());
  double coslat_dst = cos(plat);  // cos(lat)  of the search point
  double coslon_dst = cos(plon);  // cos(lon)  of the search point
  double sinlat_dst = sin(plat);  // sin(lat)  of the search point
  double sinlon_dst = sin(plon);  // sin(lon)  of the search point
  double *restrict coslon = gs->coslon;
  double *restrict sinlon = gs->sinlon;
  double *restrict coslat = gs->coslat;
  double *restrict sinlat = gs->sinlat;
  double *restrict src_center_lon = gs->reg2d_center_lon;
  double *restrict src_center_lat = gs->reg2d_center_lat;

  long nx = src_grid_dims[0];
  long ny = src_grid_dims[1];

  long nxm = nx;
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
      int ix, iy;

      for ( int na = 0; na < num_add; ++na )
	{
	  nadd = src_add[na];

	  iy = nadd/nx;
	  ix = nadd - iy*nx;

	  /* Find distance to this point */
	  distance =  sinlat_dst*sinlat[iy] + coslat_dst*coslat[iy]*
	             (coslon_dst*coslon[ix] + sinlon_dst*sinlon[ix]);

	  /* 2008-07-30 Uwe Schulzweida: check that distance is inside the range of -1 to 1,
	                                 otherwise the result of acos(distance) is NaN */
	  if ( distance > 1. ) distance =  1.;

	  if ( distance >= cos_search_radius )
	    {
	      distance = acos(distance);

	      /* Store the address and distance if this is one of the smallest four so far */
	      nbr_store_distance(nadd, distance, num_neighbors, nbr_add, nbr_dist);
	    }
	}

      nbr_check_distance(num_neighbors, nbr_add, nbr_dist);
    }
  else if ( src_grid->lextrapolate )
    {
      int search_result;

      if ( num_neighbors < 4 )
	{
	  int nbr4_add[4];
	  double nbr4_dist[4];
	  for ( n = 0; n < num_neighbors; ++n ) nbr4_add[n] = -1;
	  search_result = grid_search_reg2d_nn(nx, ny, nbr4_add, nbr4_dist, plat, plon, src_center_lat, src_center_lon);
	  if ( search_result < 0 )
	    {
	      for ( n = 0; n < num_neighbors; ++n ) nbr_add[n]  = nbr4_add[n];
	      for ( n = 0; n < num_neighbors; ++n ) nbr_dist[n] = nbr4_dist[n];
	    }
	}
      else
	{
	  search_result = grid_search_reg2d_nn(nx, ny, nbr_add, nbr_dist, plat, plon, src_center_lat, src_center_lon);
	}

      if ( search_result >= 0 )
	for ( n = 0; n < num_neighbors; ++n ) nbr_add[n] = -1;
    }
}  /*  grid_search_nbr_reg2d  */

static
void grid_search_nbr(struct gridsearch *gs, int num_neighbors, int *restrict nbr_add, double *restrict nbr_dist, double plat, double plon)
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
  int n;
  double search_radius = get_search_radius();

  /* Initialize distance and address arrays */
  for ( n = 0; n < num_neighbors; ++n ) nbr_add[n]  = -1;
  for ( n = 0; n < num_neighbors; ++n ) nbr_dist[n] = BIGNUM;

  int ndist = num_neighbors;
  ndist = ndist*2; // check some more points if distance is the same use the smaller index (nadd)
  double dist[ndist];
  int    adds[ndist];

  const double range0 = SQR(2*search_radius);
  double range = range0;

  int j = 0;

  if ( num_neighbors == 1 )
    {
      unsigned nadd = gridsearch_nearest(gs, plon, plat, &range);
      if ( nadd != GS_NOT_FOUND )
        {
          //if ( range < range0 )
            {
              dist[j] = sqrt(range);
              adds[j] = nadd;
              j++;
            }
        }
    }
  else
    {
      struct pqueue *gs_result = gridsearch_qnearest(gs, plon, plat, &range, ndist);
      if ( gs_result )
        {
          unsigned nadd;
          struct resItem *p;
          while ( pqremove_min(gs_result, &p) )
            {
              nadd  = p->node->index;
              range = p->dist_sq;
              Free(p); // Free the result node taken from the heap
              
              if ( range < range0 )
                {
                  dist[j] = sqrt(range);
                  adds[j] = nadd;
                  j++;
                }
            }
          Free(gs_result->d); // free the heap
          Free(gs_result);    // and free the heap information structure
        }
    }
  
  ndist = j;

  for ( j = 0; j < ndist; ++j )
    nbr_store_distance(adds[j], dist[j], num_neighbors, nbr_add, nbr_dist);

  nbr_check_distance(num_neighbors, nbr_add, nbr_dist);

}  /*  grid_search_nbr  */

/*
  -----------------------------------------------------------------------------------------

   This routine computes the inverse-distance weights for a nearest-neighbor interpolation.

  -----------------------------------------------------------------------------------------
*/
void remap_distwgt_weights(unsigned num_neighbors, remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv)
{
  /*  Local variables */
  int remap_grid_type = src_grid->remap_grid_type;

  if ( cdoVerbose ) cdoPrint("Called %s()", __func__);

  progressInit();

  /* Compute mappings from source to target grid */

  unsigned src_grid_size = src_grid->size;
  unsigned tgt_grid_size = tgt_grid->size;
  unsigned nx = src_grid->dims[0];
  unsigned ny = src_grid->dims[1];

  weightlinks_t *weightlinks = (weightlinks_t *) Malloc(tgt_grid_size*sizeof(weightlinks_t));

  int nbr_mask[num_neighbors];    /* mask at nearest neighbors                   */
  int nbr_add[num_neighbors];     /* source address at nearest neighbors         */
  double nbr_dist[num_neighbors]; /* angular distance four nearest neighbors     */

  struct gridsearch *gs = NULL;
  if ( remap_grid_type == REMAP_GRID_TYPE_REG2D )
    gs = gridsearch_create_reg2d(nx, ny, src_grid->reg2d_center_lon, src_grid->reg2d_center_lat);
  else if ( num_neighbors == 1 )
    gs = gridsearch_create_nn(src_grid_size, src_grid->cell_center_lon, src_grid->cell_center_lat);
  else
    gs = gridsearch_create(src_grid_size, src_grid->cell_center_lon, src_grid->cell_center_lat);

  /* Loop over destination grid  */

  double findex = 0;

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(gs, weightlinks, num_neighbors, remap_grid_type, src_grid, tgt_grid, tgt_grid_size, findex) \
  private(nbr_mask, nbr_add, nbr_dist)
#endif
  for ( unsigned tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add )
    {
#if defined(_OPENMP)
#include "pragma_omp_atomic_update.h"
#endif
      findex++;
      if ( cdo_omp_get_thread_num() == 0 ) progressStatus(0, 1, findex/tgt_grid_size);
      
      weightlinks[tgt_cell_add].nlinks = 0;	

      if ( ! tgt_grid->mask[tgt_cell_add] ) continue;
	
      double plat = tgt_grid->cell_center_lat[tgt_cell_add];
      double plon = tgt_grid->cell_center_lon[tgt_cell_add];

      /* Find nearest grid points on source grid and distances to each point */
      if ( remap_grid_type == REMAP_GRID_TYPE_REG2D )
	grid_search_nbr_reg2d(gs, num_neighbors, src_grid, nbr_add, nbr_dist, 
			      plat, plon, src_grid->dims);
      else
        grid_search_nbr(gs, num_neighbors, nbr_add, nbr_dist, plat, plon);

      /* Compute weights based on inverse distance if mask is false, eliminate those points */

      double dist_tot = nbr_compute_weights(num_neighbors, src_grid->mask, nbr_mask, nbr_add, nbr_dist);

      /* Normalize weights and store the link */

      unsigned nadds = nbr_normalize_weights(num_neighbors, dist_tot, nbr_mask, nbr_add, nbr_dist);

      for ( unsigned n = 0; n < nadds; ++n )
        if ( nbr_mask[n] ) tgt_grid->cell_frac[tgt_cell_add] = ONE;

      store_weightlinks(nadds, nbr_add, nbr_dist, tgt_cell_add, weightlinks);
    }

  progressStatus(0, 1, 1);

  if ( gs ) gridsearch_delete(gs);

  weightlinks2remaplinks(tgt_grid_size, weightlinks, rv);

  if ( weightlinks ) Free(weightlinks);

}  /* scrip_remap_weights_distwgt */

static
void distwgt_remap(double* restrict tgt_point, const double* restrict src_array, long nadds, const double wgts[4], const int src_add[4])
{
  *tgt_point = 0.;
  for ( int n = 0; n < nadds; ++n ) *tgt_point += src_array[src_add[n]]*wgts[n];
}


void remap_distwgt(unsigned num_neighbors, remapgrid_t *src_grid, remapgrid_t *tgt_grid, const double* restrict src_array, double* restrict tgt_array, double missval)
{
  /*  Local variables */
  int src_remap_grid_type = src_grid->remap_grid_type;

  if ( cdoVerbose ) cdoPrint("Called %s()", __func__);

  progressInit();

  /* Compute mappings from source to target grid */

  unsigned src_grid_size = src_grid->size;
  unsigned tgt_grid_size = tgt_grid->size;
  unsigned nx = src_grid->dims[0];
  unsigned ny = src_grid->dims[1];

  int nbr_mask[num_neighbors];    /* mask at nearest neighbors                   */
  int nbr_add[num_neighbors];     /* source address at nearest neighbors         */
  double nbr_dist[num_neighbors]; /* angular distance four nearest neighbors     */

  clock_t start;

  start = clock();

  struct gridsearch *gs = NULL;
  if ( src_remap_grid_type == REMAP_GRID_TYPE_REG2D )
    gs = gridsearch_create_reg2d(nx, ny, src_grid->reg2d_center_lon, src_grid->reg2d_center_lat);
  else if ( num_neighbors == 1 )
    gs = gridsearch_create_nn(src_grid_size, src_grid->cell_center_lon, src_grid->cell_center_lat);
  else
    gs = gridsearch_create(src_grid_size, src_grid->cell_center_lon, src_grid->cell_center_lat);

  if ( cdoVerbose ) printf("gridsearch created: %.2f seconds\n", ((double)(clock()-start))/CLOCKS_PER_SEC);

  start = clock();

  /* Loop over destination grid  */

  double findex = 0;

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  shared(gs, num_neighbors, src_remap_grid_type, src_grid, tgt_grid, tgt_grid_size, findex) \
  shared(src_array, tgt_array, missval) \
  private(nbr_mask, nbr_add, nbr_dist)
#endif
  for ( unsigned tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add )
    {
#if defined(_OPENMP)
#include "pragma_omp_atomic_update.h"
#endif
      findex++;
      if ( cdo_omp_get_thread_num() == 0 ) progressStatus(0, 1, findex/tgt_grid_size);
      
      tgt_array[tgt_cell_add] = missval;

      if ( ! tgt_grid->mask[tgt_cell_add] ) continue;

      double plon = 0, plat = 0;
      remapgrid_get_lonlat(tgt_grid, tgt_cell_add, &plon, &plat);

      /* Find nearest grid points on source grid and distances to each point */
      if ( src_remap_grid_type == REMAP_GRID_TYPE_REG2D )
	grid_search_nbr_reg2d(gs, num_neighbors, src_grid, nbr_add, nbr_dist, 
			      plat, plon, src_grid->dims);
      else
        grid_search_nbr(gs, num_neighbors, nbr_add, nbr_dist, plat, plon);
      
      /* Compute weights based on inverse distance if mask is false, eliminate those points */

      double dist_tot = nbr_compute_weights(num_neighbors, src_grid->mask, nbr_mask, nbr_add, nbr_dist);

      /* Normalize weights and store the link */

      unsigned nadds = nbr_normalize_weights(num_neighbors, dist_tot, nbr_mask, nbr_add, nbr_dist);

      for ( unsigned n = 0; n < nadds; ++n )
        if ( nbr_mask[n] ) tgt_grid->cell_frac[tgt_cell_add] = ONE;

      if ( nadds > 1 ) sort_add_and_wgts(nadds, nbr_add, nbr_dist);

      if ( nadds ) distwgt_remap(&tgt_array[tgt_cell_add], src_array, nadds, nbr_dist, nbr_add);
    }

  progressStatus(0, 1, 1);

  if ( gs ) gridsearch_delete(gs);

  if ( cdoVerbose ) printf("gridsearch nearest: %.2f seconds\n", ((double)(clock()-start))/CLOCKS_PER_SEC);

}  /* scrip_remap_distwgt */

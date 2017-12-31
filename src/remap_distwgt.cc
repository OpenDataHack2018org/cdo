#ifdef  _OPENMP
#include <omp.h> // omp_get_wtime
#endif

#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "remap.h"
#include "remap_store_link.h"
#include "grid_search.h"

//  Interpolation using a distance-weighted average

static
void nbr_store_distance(size_t nadd, double distance, size_t numNeighbors, size_t *restrict nbr_add, double *restrict nbr_dist)
{
  if ( numNeighbors == 1 )
    {
      if ( distance < nbr_dist[0] || (distance <= nbr_dist[0] && nadd < nbr_add[0]) )
	{
	  nbr_add[0]  = nadd;
	  nbr_dist[0] = distance;
	}
    }
  else
    {
      for ( size_t nchk = 0; nchk < numNeighbors; ++nchk )
	{
	  if ( distance < nbr_dist[nchk] || (distance <= nbr_dist[nchk] && nadd < nbr_add[nchk]) )
	    {
	      for ( size_t n = numNeighbors-1; n > nchk; --n )
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
void nbr_check_distance(size_t numNeighbors, const size_t *restrict nbr_add, double *restrict nbr_dist)
{
  // If distance is zero, set to small number
  for ( size_t nchk = 0; nchk < numNeighbors; ++nchk )
    if ( nbr_add[nchk] < SIZE_MAX && nbr_dist[nchk] <= 0. ) nbr_dist[nchk] = TINY;
}


double nbr_compute_weights(size_t numNeighbors, const int *restrict src_grid_mask, bool *restrict nbr_mask, const size_t *restrict nbr_add, double *restrict nbr_dist)
{
  // Compute weights based on inverse distance if mask is false, eliminate those points

  double dist_tot = 0.; // sum of neighbor distances (for normalizing)

  if ( src_grid_mask )
    {
      for ( size_t n = 0; n < numNeighbors; ++n )
        {
          nbr_mask[n] = false;
          if ( nbr_add[n] < SIZE_MAX )
            if ( src_grid_mask[nbr_add[n]] )
              {
                nbr_dist[n] = 1./nbr_dist[n];
                dist_tot += nbr_dist[n];
                nbr_mask[n] = true;
              }
        }
    }
  else
    {
      for ( size_t n = 0; n < numNeighbors; ++n )
        {
          nbr_mask[n] = false;
          if ( nbr_add[n] < SIZE_MAX )
            {
              nbr_dist[n] = 1./nbr_dist[n];
              dist_tot += nbr_dist[n];
              nbr_mask[n] = true;
            }
        }
    }
  

  return dist_tot;
}


size_t nbr_normalize_weights(size_t numNeighbors, double dist_tot, const bool *restrict nbr_mask, size_t *restrict nbr_add, double *restrict nbr_dist)
{
  // Normalize weights and store the link

  size_t nadds = 0;

  for ( size_t n = 0; n < numNeighbors; ++n )
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


//  This routine finds the closest num_neighbor points to a search point and computes a distance to each of the neighbors.

template <typename T>
static inline
void LLtoXYZ(double lon, double lat, T *restrict xyz)
{
   double cos_lat = cos(lat);
   xyz[0] = cos_lat * cos(lon);
   xyz[1] = cos_lat * sin(lon);
   xyz[2] = sin(lat);
}

template <typename T>
static constexpr
T square(const T x)
{
  return x*x;
}

template <typename T>
static constexpr
T distance(const T *restrict a, const T *restrict b)
{
  return square(a[0]-b[0]) + square(a[1]-b[1]) + square(a[2]-b[2]);
}


#define MAX_SEARCH_CELLS 25
static
void grid_search_nbr_reg2d(struct gridsearch *gs, size_t numNeighbors, size_t *restrict nbr_add, double *restrict nbr_dist, 
			   double plon, double plat)
{
  /*
    Output variables:

    int nbr_add[numNeighbors]     ! address of each of the closest points
    double nbr_dist[numNeighbors] ! distance to each of the closest points

    Input variables:

    double plat,         ! latitude  of the search point
    double plon,         ! longitude of the search point
  */
  size_t src_add[MAX_SEARCH_CELLS];
  size_t *src_add_tmp = NULL;
  size_t *psrc_add = src_add;
  size_t num_add = 0;
  double *restrict src_center_lon = gs->reg2d_center_lon;
  double *restrict src_center_lat = gs->reg2d_center_lat;

  long nx = gs->dims[0];
  long ny = gs->dims[1];

  size_t nxm = gs->is_cyclic ? nx+1 : nx;

  if ( plon < src_center_lon[0]     ) plon += PI2;
  if ( plon > src_center_lon[nxm-1] ) plon -= PI2;

  size_t ii, jj;
  int lfound = rect_grid_search(&ii, &jj, plon, plat, nxm, ny, src_center_lon, src_center_lat);
  if ( lfound )
    {
      if ( gs->is_cyclic && ii == (nxm-1) ) ii = 0;

      long k;
      for ( k = 3; k < 10000; k+=2 )
        if ( numNeighbors <= (k-2)*(k-2) ) break;

      if ( (k*k) > MAX_SEARCH_CELLS ) psrc_add = src_add_tmp = (size_t *) Malloc(k*k*sizeof(size_t));

      k /= 2;

      long j0 = jj-k;
      long jn = jj+k;
      long i0 = ii-k;
      long in = ii+k;
      if ( j0 < 0 ) j0 = 0;
      if ( jn >= ny ) jn = ny-1;
      if ( (in-i0) > nx ) { i0 = 0; in = nx-1; }
      
      for ( long j = j0; j <= jn; ++j )
	for ( long i = i0; i <= in; ++i )
	  {
	    long ix = i;
	    if ( gs->is_cyclic )
	      {
		if ( ix <   0 ) ix += nx;
		if ( ix >= nx ) ix -= nx;
	      }

	    if ( ix >= 0 && ix < nx && j >= 0 && j < ny )
	      psrc_add[num_add++] = j*nx+ix;
	  }
    }

  // Initialize distance and address arrays
  for ( size_t n = 0; n < numNeighbors; ++n )
    {
      nbr_add[n]  = SIZE_MAX;
      nbr_dist[n] = BIGNUM;
    }

  if ( lfound )
    {
      double *restrict coslon = gs->coslon;
      double *restrict sinlon = gs->sinlon;
      double *restrict coslat = gs->coslat;
      double *restrict sinlat = gs->sinlat;

      double xyz[3];
      double query_pt[3];
      LLtoXYZ<double>(plon, plat, query_pt);
      double search_radius = SQR(gs->search_radius);

      for ( size_t na = 0; na < num_add; ++na )
	{
	  size_t nadd = psrc_add[na];
	  size_t iy = nadd/nx;
	  size_t ix = nadd - iy*nx;

          xyz[0] = coslat[iy] * coslon[ix];
          xyz[1] = coslat[iy] * sinlon[ix];
          xyz[2] = sinlat[iy];
	  // Find distance to this point
          double dist = (float) distance<double>(query_pt, xyz);
	  if ( dist <= search_radius )
	    {
	      // Store the address and distance if this is one of the smallest so far
	      nbr_store_distance(nadd, sqrt(dist), numNeighbors, nbr_add, nbr_dist);
	    }
	}

      nbr_check_distance(numNeighbors, nbr_add, nbr_dist);
    }
  else if ( gs->extrapolate )
    {
      int search_result;

      if ( numNeighbors < 4 )
	{
	  size_t nbr_add4[4];
	  double nbr_dist4[4];
	  for ( size_t n = 0; n < numNeighbors; ++n ) nbr_add4[n] = SIZE_MAX;
	  search_result = grid_search_reg2d_nn(nx, ny, nbr_add4, nbr_dist4, plat, plon, src_center_lat, src_center_lon);
	  if ( search_result < 0 )
	    {
	      for ( size_t n = 0; n < numNeighbors; ++n ) nbr_add[n]  = nbr_add4[n];
	      for ( size_t n = 0; n < numNeighbors; ++n ) nbr_dist[n] = nbr_dist4[n];
	    }
	}
      else
	{
	  search_result = grid_search_reg2d_nn(nx, ny, nbr_add, nbr_dist, plat, plon, src_center_lat, src_center_lon);
	}

      if ( search_result >= 0 )
	for ( size_t n = 0; n < numNeighbors; ++n ) nbr_add[n] = SIZE_MAX;
    }

  if ( src_add_tmp ) Free(src_add_tmp);
} // grid_search_nbr_reg2d


int grid_search_nbr(struct gridsearch *gs, size_t numNeighbors, size_t *restrict nbr_add, double *restrict nbr_dist, double plon, double plat)
{
  /*
    Output variables:

    int nbr_add[numNeighbors]     ! address of each of the closest points
    double nbr_dist[numNeighbors] ! distance to each of the closest points

    Input variables:

    double plat,         ! latitude  of the search point
    double plon,         ! longitude of the search point
  */

  // Initialize distance and address arrays
  for ( size_t n = 0; n < numNeighbors; ++n ) nbr_add[n]  = SIZE_MAX;
  for ( size_t n = 0; n < numNeighbors; ++n ) nbr_dist[n] = BIGNUM;

  size_t ndist = numNeighbors;
  // check some more points if distance is the same use the smaller index (nadd)
  if ( ndist > 8 ) ndist += 8;
  else             ndist *= 2; 
  if ( ndist > gs->n ) ndist = gs->n;

  double zdist[32];
  size_t zadds[32];
  double *dist = numNeighbors > 16 ? (double*) Malloc(ndist*sizeof(double)) : zdist;
  size_t *adds = numNeighbors > 16 ? (size_t*) Malloc(ndist*sizeof(size_t)) : zadds;
  
  const double range0 = SQR(gs->search_radius);
  double range = range0;

  size_t nadds = 0;

  if ( numNeighbors == 1 )
    {
      size_t add = gridsearch_nearest(gs, plon, plat, &range);
      if ( add != GS_NOT_FOUND )
        {
          //if ( range < range0 )
            {
              dist[nadds] = sqrt(range);
              adds[nadds] = add;
              nadds++;
            }
        }
    }
  else
    {
      nadds = gridsearch_qnearest(gs, plon, plat, &range, ndist, adds, dist);
      for ( size_t i = 0; i < nadds; ++i ) dist[i] = sqrt(dist[i]);
    }

  ndist = nadds;
  size_t max_neighbors = (ndist < numNeighbors) ? ndist : numNeighbors;

  for ( size_t i = 0; i < ndist; ++i )
    nbr_store_distance(adds[i], dist[i], max_neighbors, nbr_add, nbr_dist);

  nbr_check_distance(max_neighbors, nbr_add, nbr_dist);

  if ( numNeighbors > 16 )
    {
      Free(dist);
      Free(adds);
    }

  if ( ndist > numNeighbors ) ndist = numNeighbors;

  return ndist;
} // grid_search_nbr

//  This routine computes the inverse-distance weights for a nearest-neighbor interpolation.

void remap_distwgt_weights(size_t numNeighbors, remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv)
{
  int remap_grid_type = src_grid->remap_grid_type;

  if ( cdoVerbose ) cdoPrint("Called %s()", __func__);

  progressInit();

  // Compute mappings from source to target grid

  size_t src_grid_size = src_grid->size;
  size_t tgt_grid_size = tgt_grid->size;

  weightlinks_t *weightlinks = (weightlinks_t *) Malloc(tgt_grid_size*sizeof(weightlinks_t));
  weightlinks[0].addweights = (addweight_t *) Malloc(numNeighbors*tgt_grid_size*sizeof(addweight_t));
  for ( size_t tgt_cell_add = 1; tgt_cell_add < tgt_grid_size; ++tgt_cell_add )
    weightlinks[tgt_cell_add].addweights = weightlinks[0].addweights + numNeighbors*tgt_cell_add;

  NEW_2D(bool, nbr_mask, ompNumThreads, numNeighbors);   // mask at nearest neighbors
  NEW_2D(size_t, nbr_add, ompNumThreads, numNeighbors);  // source address at nearest neighbors
  NEW_2D(double, nbr_dist, ompNumThreads, numNeighbors); // angular distance four nearest neighbors

#ifdef  _OPENMP
  double start = cdoVerbose ? omp_get_wtime() : 0;
#endif

  bool xIsCyclic = src_grid->is_cyclic;
  size_t *dims = src_grid->dims;
  struct gridsearch *gs = NULL;
  if ( remap_grid_type == REMAP_GRID_TYPE_REG2D )
    gs = gridsearch_create_reg2d(xIsCyclic, dims, src_grid->reg2d_center_lon, src_grid->reg2d_center_lat);
  else
    gs = gridsearch_create(xIsCyclic, dims, src_grid_size, src_grid->cell_center_lon, src_grid->cell_center_lat);

  if ( src_grid->lextrapolate ) gridsearch_extrapolate(gs);
  //else
  /*
   if ( remap_grid_type != REMAP_GRID_TYPE_REG2D )
    gridsearch_bound_poly(gs, src_grid->dims,  src_grid_size, src_grid->cell_center_lon, src_grid->cell_center_lat);
  */
#ifdef  _OPENMP
  if ( cdoVerbose ) printf("gridsearch created: %.2f seconds\n", omp_get_wtime()-start);
  if ( cdoVerbose ) start = omp_get_wtime();
#endif

  // Loop over destination grid

  double findex = 0;

#ifdef  HAVE_OPENMP4
#pragma omp parallel for default(none)  reduction(+:findex) \
  shared(gs, weightlinks, numNeighbors, remap_grid_type, src_grid, tgt_grid, tgt_grid_size) \
  shared(nbr_mask, nbr_add, nbr_dist)
#endif
  for ( size_t tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add )
    {
      findex++;
      if ( cdo_omp_get_thread_num() == 0 ) progressStatus(0, 1, findex/tgt_grid_size);
      
      int ompthID = cdo_omp_get_thread_num();

      weightlinks[tgt_cell_add].nlinks = 0;	

      if ( ! tgt_grid->mask[tgt_cell_add] ) continue;
	
      double plon = 0, plat = 0;
      remapgrid_get_lonlat(tgt_grid, tgt_cell_add, &plon, &plat);

      // Find nearest grid points on source grid and distances to each point
      if ( remap_grid_type == REMAP_GRID_TYPE_REG2D )
	grid_search_nbr_reg2d(gs, numNeighbors, nbr_add[ompthID], nbr_dist[ompthID], plon, plat);
      else
        grid_search_nbr(gs, numNeighbors, nbr_add[ompthID], nbr_dist[ompthID], plon, plat);

      // Compute weights based on inverse distance if mask is false, eliminate those points
      double dist_tot = nbr_compute_weights(numNeighbors, src_grid->mask, nbr_mask[ompthID], nbr_add[ompthID], nbr_dist[ompthID]);

      // Normalize weights and store the link
      size_t nadds = nbr_normalize_weights(numNeighbors, dist_tot, nbr_mask[ompthID], nbr_add[ompthID], nbr_dist[ompthID]);

      for ( size_t n = 0; n < nadds; ++n )
        if ( nbr_mask[ompthID][n] ) tgt_grid->cell_frac[tgt_cell_add] = ONE;

      store_weightlinks(0, nadds, nbr_add[ompthID], nbr_dist[ompthID], tgt_cell_add, weightlinks);
    }

  progressStatus(0, 1, 1);

  DELETE_2D(nbr_mask);
  DELETE_2D(nbr_add);
  DELETE_2D(nbr_dist);

  if ( gs ) gridsearch_delete(gs);

  weightlinks2remaplinks(0, tgt_grid_size, weightlinks, rv);

  if ( weightlinks ) Free(weightlinks);

#ifdef  _OPENMP
  if ( cdoVerbose ) printf("gridsearch nearest: %.2f seconds\n", omp_get_wtime()-start);
#endif
} // remap_distwgt_weights

static
void distwgt_remap(double *restrict tgt_point, const double* restrict src_array, size_t nadds, const double wgts[4], const size_t src_add[4])
{
  if ( nadds ) *tgt_point = src_array[src_add[0]]*wgts[0];
  for ( size_t n = 1; n < nadds; ++n ) *tgt_point += src_array[src_add[n]]*wgts[n];
}


void remap_distwgt(size_t numNeighbors, remapgrid_t *src_grid, remapgrid_t *tgt_grid, const double* restrict src_array, double* restrict tgt_array, double missval)
{
  int src_remap_grid_type = src_grid->remap_grid_type;

  if ( cdoVerbose ) cdoPrint("Called %s()", __func__);

  progressInit();

  // Compute mappings from source to target grid 

  size_t src_grid_size = src_grid->size;
  size_t tgt_grid_size = tgt_grid->size;

  NEW_2D(bool, nbr_mask, ompNumThreads, numNeighbors);   // mask at nearest neighbors
  NEW_2D(size_t, nbr_add, ompNumThreads, numNeighbors);  // source address at nearest neighbors
  NEW_2D(double, nbr_dist, ompNumThreads, numNeighbors); // angular distance four nearest neighbors

#ifdef  _OPENMP
  double start = cdoVerbose ? omp_get_wtime() : 0;
#endif

  bool xIsCyclic = src_grid->is_cyclic;
  size_t *dims = src_grid->dims;
  struct gridsearch *gs = NULL;
  if ( src_remap_grid_type == REMAP_GRID_TYPE_REG2D )
    gs = gridsearch_create_reg2d(xIsCyclic, dims, src_grid->reg2d_center_lon, src_grid->reg2d_center_lat);
  else
    gs = gridsearch_create(xIsCyclic, dims, src_grid_size, src_grid->cell_center_lon, src_grid->cell_center_lat);

  if ( src_grid->lextrapolate ) gridsearch_extrapolate(gs);

#ifdef  _OPENMP
  if ( cdoVerbose ) printf("gridsearch created: %.2f seconds\n", omp_get_wtime()-start);
  if ( cdoVerbose ) start = omp_get_wtime();
#endif

  // Loop over destination grid

  double findex = 0;

#ifdef  HAVE_OPENMP4
#pragma omp parallel for default(none)  reduction(+:findex) \
  shared(gs, numNeighbors, src_remap_grid_type, src_grid, tgt_grid, tgt_grid_size) \
  shared(src_array, tgt_array, missval, nbr_mask, nbr_add, nbr_dist)
#endif
  for ( size_t tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add )
    {
      findex++;
      if ( cdo_omp_get_thread_num() == 0 ) progressStatus(0, 1, findex/tgt_grid_size);
      
      int ompthID = cdo_omp_get_thread_num();

      tgt_array[tgt_cell_add] = missval;

      if ( ! tgt_grid->mask[tgt_cell_add] ) continue;

      double plon = 0, plat = 0;
      remapgrid_get_lonlat(tgt_grid, tgt_cell_add, &plon, &plat);

      // Find nearest grid points on source grid and distances to each point
      if ( src_remap_grid_type == REMAP_GRID_TYPE_REG2D )
	grid_search_nbr_reg2d(gs, numNeighbors, nbr_add[ompthID], nbr_dist[ompthID], plon, plat);
      else
        grid_search_nbr(gs, numNeighbors, nbr_add[ompthID], nbr_dist[ompthID], plon, plat);
      
      // Compute weights based on inverse distance if mask is false, eliminate those points
      double dist_tot = nbr_compute_weights(numNeighbors, src_grid->mask, nbr_mask[ompthID], nbr_add[ompthID], nbr_dist[ompthID]);

      // Normalize weights and store the link
      size_t nadds = nbr_normalize_weights(numNeighbors, dist_tot, nbr_mask[ompthID], nbr_add[ompthID], nbr_dist[ompthID]);

      if ( nadds > 1 ) sort_add_and_wgts(nadds, nbr_add[ompthID], nbr_dist[ompthID]);

      if ( nadds ) distwgt_remap(&tgt_array[tgt_cell_add], src_array, nadds, nbr_dist[ompthID], nbr_add[ompthID]);
    }

  progressStatus(0, 1, 1);

  DELETE_2D(nbr_mask);
  DELETE_2D(nbr_add);
  DELETE_2D(nbr_dist);

  if ( gs ) gridsearch_delete(gs);

#ifdef  _OPENMP
  if ( cdoVerbose ) printf("gridsearch nearest: %.2f seconds\n", omp_get_wtime()-start);
#endif
} // remap_distwgt


#include <cdi.h>

void intgriddis(field_type *field1, field_type *field2, size_t numNeighbors)
{
  int gridID1 = field1->grid;
  int gridID2 = field2->grid;
  double src_missval = field1->missval;
  double tgt_missval = field2->missval;
  double *src_array = field1->ptr;
  double *tgt_array = field2->ptr;
  // int src_remap_grid_type = src_grid->remap_grid_type;

  if ( cdoVerbose ) cdoPrint("Called %s()", __func__);

  progressInit();

  // Compute mappings from source to target grid

  int src_gridtype = gridInqType(gridID1);
  int tgt_gridtype = gridInqType(gridID2);
  if ( src_gridtype != GRID_CURVILINEAR && src_gridtype != GRID_UNSTRUCTURED )
    cdoAbort("Source grid must be curvilinear or unstructured!");
  if ( tgt_gridtype != GRID_CURVILINEAR && tgt_gridtype != GRID_UNSTRUCTURED )
    cdoAbort("Target grid must be curvilinear or unstructured!");

  size_t src_grid_size = gridInqSize(gridID1);
  size_t tgt_grid_size = gridInqSize(gridID2);

  int *src_mask = (int*) Malloc(src_grid_size*sizeof(int));
  for ( size_t i = 0; i < src_grid_size; ++i ) src_mask[i] = !DBL_IS_EQUAL(src_array[i], src_missval);
  int *tgt_mask = (int*) Malloc(tgt_grid_size*sizeof(int));
  for ( size_t i = 0; i < tgt_grid_size; ++i ) tgt_mask[i] = 1;

  double *src_cell_center_lon = (double*) Malloc(src_grid_size*sizeof(double));
  double *src_cell_center_lat = (double*) Malloc(src_grid_size*sizeof(double));
  gridInqXvals(gridID1, src_cell_center_lon);
  gridInqYvals(gridID1, src_cell_center_lat);

  double *tgt_cell_center_lon = (double*) Malloc(tgt_grid_size*sizeof(double));
  double *tgt_cell_center_lat = (double*) Malloc(tgt_grid_size*sizeof(double));
  gridInqXvals(gridID2, tgt_cell_center_lon);
  gridInqYvals(gridID2, tgt_cell_center_lat);

  char xunits[CDI_MAX_NAME]; xunits[0] = 0;
  char yunits[CDI_MAX_NAME]; yunits[0] = 0;
  cdiGridInqKeyStr(gridID1, CDI_KEY_XUNITS, CDI_MAX_NAME, xunits);
  cdiGridInqKeyStr(gridID1, CDI_KEY_YUNITS, CDI_MAX_NAME, yunits);
  grid_to_radian(xunits, src_grid_size, src_cell_center_lon, "src cell center lon");
  grid_to_radian(yunits, src_grid_size, src_cell_center_lat, "src cell center lat");
  xunits[0] = 0; yunits[0] = 0;
  cdiGridInqKeyStr(gridID2, CDI_KEY_XUNITS, CDI_MAX_NAME, xunits);
  cdiGridInqKeyStr(gridID2, CDI_KEY_YUNITS, CDI_MAX_NAME, yunits);
  grid_to_radian(xunits, tgt_grid_size, tgt_cell_center_lon, "tgt cell center lon");
  grid_to_radian(yunits, tgt_grid_size, tgt_cell_center_lat, "tgt cell center lat");


  NEW_2D(bool, nbr_mask, ompNumThreads, numNeighbors);   // mask at nearest neighbors
  NEW_2D(size_t, nbr_add, ompNumThreads, numNeighbors);  // source address at nearest neighbors
  NEW_2D(double, nbr_dist, ompNumThreads, numNeighbors); // angular distance four nearest neighbors

#ifdef  _OPENMP
  double start = cdoVerbose ? omp_get_wtime() : 0;
#endif

  bool xIsCyclic = gridIsCircular(gridID1);
  size_t dims[2] = {src_grid_size, 0};
  struct gridsearch *gs = NULL;
  // if ( src_remap_grid_type == REMAP_GRID_TYPE_REG2D )
  //  gs = gridsearch_create_reg2d(xIsCyclic, dims, src_grid->reg2d_center_lon, src_grid->reg2d_center_lat);
  gs = gridsearch_create(xIsCyclic, dims, src_grid_size, src_cell_center_lon, src_cell_center_lat);

  // if ( src_grid->lextrapolate ) gridsearch_extrapolate(gs);
  // gridsearch_extrapolate(gs);

#ifdef  _OPENMP
  if ( cdoVerbose ) printf("gridsearch created: %.2f seconds\n", omp_get_wtime()-start);
  if ( cdoVerbose ) start = omp_get_wtime();
#endif

  // Loop over destination grid

  size_t nmiss = 0;
  double findex = 0;

#ifdef  HAVE_OPENMP4
  /*
#pragma omp parallel for default(none)  reduction(+:findex) \
  shared(gs, numNeighbors, src_grid, tgt_grid, tgt_grid_size)  \
  shared(src_array, tgt_array, missval, nbr_mask, nbr_add, nbr_dist)
  */
#endif
  for ( size_t tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add )
    {
      findex++;
      if ( cdo_omp_get_thread_num() == 0 ) progressStatus(0, 1, findex/tgt_grid_size);
      
      int ompthID = cdo_omp_get_thread_num();

      tgt_array[tgt_cell_add] = tgt_missval;

      if ( ! tgt_mask[tgt_cell_add] ) continue;

      double plon = 0, plat = 0;
      // remapgrid_get_lonlat(tgt_grid, tgt_cell_add, &plon, &plat);
      plat = tgt_cell_center_lat[tgt_cell_add];
      plon = tgt_cell_center_lon[tgt_cell_add];

      // Find nearest grid points on source grid and distances to each point
      // if ( src_remap_grid_type == REMAP_GRID_TYPE_REG2D )
      //   grid_search_nbr_reg2d(gs, numNeighbors, nbr_add[ompthID], nbr_dist[ompthID], plon, plat);
      // else
        grid_search_nbr(gs, numNeighbors, nbr_add[ompthID], nbr_dist[ompthID], plon, plat);
      
      // Compute weights based on inverse distance if mask is false, eliminate those points
      double dist_tot = nbr_compute_weights(numNeighbors, src_mask, nbr_mask[ompthID], nbr_add[ompthID], nbr_dist[ompthID]);

      // Normalize weights and store the link
      size_t nadds = nbr_normalize_weights(numNeighbors, dist_tot, nbr_mask[ompthID], nbr_add[ompthID], nbr_dist[ompthID]);

      if ( nadds > 1 ) sort_add_and_wgts(nadds, nbr_add[ompthID], nbr_dist[ompthID]);

      if ( nadds ) distwgt_remap(&tgt_array[tgt_cell_add], src_array, nadds, nbr_dist[ompthID], nbr_add[ompthID]);
      else nmiss++;
    }

  progressStatus(0, 1, 1);
  
  field2->nmiss = nmiss;

  DELETE_2D(nbr_mask);
  DELETE_2D(nbr_add);
  DELETE_2D(nbr_dist);

  if ( gs ) gridsearch_delete(gs);

  Free(src_mask);
  Free(tgt_mask);
  Free(src_cell_center_lon);
  Free(src_cell_center_lat);
  Free(tgt_cell_center_lon);
  Free(tgt_cell_center_lat);

#ifdef  _OPENMP
  if ( cdoVerbose ) printf("gridsearch nearest: %.2f seconds\n", omp_get_wtime()-start);
#endif
} // intgriddis

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "dmemory.h"
#include "util.h"
#include "grid_search.h"


#ifndef  M_PI
#define  M_PI        3.14159265358979323846264338327950288  /* pi */
#endif

#define  PI       M_PI
#define  PI2      (2.0*PI)


static int gridsearch_method_nn = GS_KDTREE;

static inline void LLtoXYZ_f(double lon, double lat, float *restrict xyz)
{
   double cos_lat = cos(lat);
   xyz[0] = cos_lat * cos(lon);
   xyz[1] = cos_lat * sin(lon);
   xyz[2] = sin(lat);
}


void gridsearch_set_method(int method)
{
  gridsearch_method_nn = method;
}


struct gridsearch *gridsearch_create_reg2d(unsigned nx, unsigned ny, const double *restrict lons, const double *restrict lats)
{
  struct gridsearch *gs = (struct gridsearch *) calloc(1, sizeof(struct gridsearch));

  gs->nx = nx;
  gs->ny = ny;

  double *reg2d_center_lon = (double *) malloc((nx+1)*sizeof(double));
  double *reg2d_center_lat = (double *) malloc(ny*sizeof(double));

  memcpy(reg2d_center_lon, lons, (nx+1)*sizeof(double));
  memcpy(reg2d_center_lat, lats, ny*sizeof(double));

  double *coslon = (double *) malloc(nx*sizeof(double));
  double *sinlon = (double *) malloc(nx*sizeof(double));
  double *coslat = (double *) malloc(ny*sizeof(double));
  double *sinlat = (double *) malloc(ny*sizeof(double));

  for ( unsigned n = 0; n < nx; ++n )
    {
      double rlon = lons[n];
      if ( rlon > PI2 ) rlon -= PI2;
      if ( rlon < 0   ) rlon += PI2;
      coslon[n] = cos(rlon);
      sinlon[n] = sin(rlon);
    }
  for ( unsigned n = 0; n < ny; ++n )
    {
      coslat[n] = cos(lats[n]);
      sinlat[n] = sin(lats[n]);
    }

  gs->reg2d_center_lon = reg2d_center_lon;
  gs->reg2d_center_lat = reg2d_center_lat;

  gs->coslon = coslon;
  gs->sinlon = sinlon;
  gs->coslat = coslat;
  gs->sinlat = sinlat;

  return gs;
}


struct kdNode *kdtree_index_create(unsigned n, const double *restrict lons, const double *restrict lats, const unsigned *restrict index)
{
  struct kd_point *pointlist = (struct kd_point *) malloc(n * sizeof(struct kd_point));  
  // see  example_cartesian.c
  if ( cdoVerbose) printf("kdtree lib init 3D: n=%d  nthreads=%d\n", n, ompNumThreads);
  float min[3], max[3];
  min[0] = min[1] = min[2] =  1e9;
  max[0] = max[1] = max[2] = -1e9;
  float *restrict point;
  for ( unsigned i = 0; i < n; i++ ) 
    {
      point = pointlist[i].point;
      LLtoXYZ_f(lons[index[i]], lats[index[i]], point);
      for ( unsigned j = 0; j < 3; ++j )
        {
          min[j] = point[j] < min[j] ? point[j] : min[j];
          max[j] = point[j] > max[j] ? point[j] : max[j];
        }
      pointlist[i].index = index[i];
    }

  struct kdNode *kdt = kd_buildTree(pointlist, n, min, max, 3, ompNumThreads);
  if ( pointlist ) free(pointlist);

   return kdt;
}


struct gridsearch *gridsearch_index_create(unsigned n, const double *restrict lons, const double *restrict lats, const unsigned *restrict index)
{
  struct gridsearch *gs = (struct gridsearch *) calloc(1, sizeof(struct gridsearch));

  gs->n = n;
  if ( n == 0 ) return gs;

  gs->kdt = kdtree_index_create(n, lons, lats, index);

  return gs;
}


struct gridsearch *gridsearch_index_create_nn(unsigned n, const double *restrict lons, const double *restrict lats, const unsigned *restrict index)
{
  struct gridsearch *gs = (struct gridsearch *) calloc(1, sizeof(struct gridsearch));

  gs->method_nn = gridsearch_method_nn;
  gs->n = n;
  if ( n == 0 ) return gs;

  if ( gs->method_nn == GS_KDTREE )
    gs->kdt = kdtree_index_create(n, lons, lats, index);
  // else if ( gridsearch_method_nn == GS_NEARPT3 )
   
  return gs;
}


struct gridsearch *gridsearch_create(unsigned n, const double *restrict lons, const double *restrict lats)
{
  struct gridsearch *gs = (struct gridsearch *) calloc(1, sizeof(struct gridsearch));

  gs->n = n;

  if ( n == 0 ) return gs;

  struct kd_point *pointlist = (struct kd_point *) malloc(n * sizeof(struct kd_point));  
  // see  example_cartesian.c
  if ( cdoVerbose) printf("kdtree lib init 3D: n=%d  nthreads=%d\n", n, ompNumThreads);
  float min[3], max[3];
  min[0] = min[1] = min[2] =  1e9;
  max[0] = max[1] = max[2] = -1e9;
  float *restrict point;
  for ( unsigned i = 0; i < n; i++ ) 
    {
      point = pointlist[i].point;
      LLtoXYZ_f(lons[i], lats[i], point);
      for ( unsigned j = 0; j < 3; ++j )
        {
          min[j] = point[j] < min[j] ? point[j] : min[j];
          max[j] = point[j] > max[j] ? point[j] : max[j];
        }
      pointlist[i].index = i;
    }

  gs->kdt = kd_buildTree(pointlist, n, min, max, 3, ompNumThreads);
  if ( pointlist ) free(pointlist);
 
  return gs;
}


void gridsearch_delete(struct gridsearch *gs)
{
  if ( gs )
    {
      if ( gs->kdt ) kd_destroyTree(gs->kdt);
      
      if ( gs->reg2d_center_lon ) free(gs->reg2d_center_lon);
      if ( gs->reg2d_center_lat ) free(gs->reg2d_center_lat);

      if ( gs->coslat ) free(gs->coslat);
      if ( gs->coslon ) free(gs->coslon);
      if ( gs->sinlat ) free(gs->sinlat);
      if ( gs->sinlon ) free(gs->sinlon);

      free(gs);
    }
}


kdNode *kdtree_nearest(kdNode *kdt, double lon, double lat, double *prange)
{
  if ( kdt == NULL ) return NULL;
  
  float point[3];
  float range0;
  if ( prange )
    range0 = *prange;
  else
    range0 = SQR(2 * M_PI);     /* This has to be bigger than the presumed
                                 * maximum distance to the NN but smaller
                                 * than once around the sphere. The content
                                 * of this variable is replaced with the
                                 * distance to the NN squared. */
  float range = range0;

  LLtoXYZ_f(lon, lat, point);

  kdNode *node = kd_nearest(kdt, point, &range, 3);
  // printf("range %g %g %g %p\n", lon, lat, range, node);

  if ( !(range < range0) ) node = NULL;
  if ( prange ) *prange = range;

  return node;
}


unsigned gridsearch_nearest(struct gridsearch *gs, double lon, double lat, double *prange)
{
  unsigned index = GS_NOT_FOUND;

  if ( gs->kdt )
    {
      kdNode *node = kdtree_nearest(gs->kdt, lon, lat, prange);
      if ( node ) index = (int) node->index;
    }

  return index;
}


struct pqueue *gridsearch_qnearest(struct gridsearch *gs, double lon, double lat, double *prange, unsigned nnn)
{
  if ( gs->kdt == NULL ) return NULL;
  
  float point[3];
  float range0;
  if ( prange )
    range0 = *prange;
  else
    range0 = SQR(2 * M_PI);     /* This has to be bigger than the presumed
                                 * maximum distance to the NN but smaller
                                 * than once around the sphere. The content
                                 * of this variable is replaced with the
                                 * distance to the NN squared. */
  float range = range0;

  LLtoXYZ_f(lon, lat, point);

  struct pqueue *result = kd_qnearest(gs->kdt, point, &range, nnn, 3);
  // printf("range %g %g %g %p\n", lon, lat, range, node);

  if ( !(range < range0) ) result = NULL;
  if ( prange ) *prange = range;

  return result;
}


unsigned gridsearch_item(void *gs_result)
{
  unsigned index = 0;

  if ( gs_result == NULL ) return 0;

  index = ((kdNode *) gs_result)->index;
  
  return index;
}

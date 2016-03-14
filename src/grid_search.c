#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "cdo_int.h"
#include "dmemory.h"
#include "util.h"
#include "grid.h"
#include "grid_search.h"


#define  PI       M_PI
#define  PI2      (2.0*PI)


static int gridsearch_method_nn = GS_KDTREE;

static
double cdo_default_search_radius(void)
{
  extern double gridsearch_radius;

  double search_radius = gridsearch_radius;

  if ( search_radius <    0. ) search_radius = 0.;
  if ( search_radius >  180. ) search_radius = 180.;

  search_radius = search_radius*DEG2RAD;

  return search_radius;
}

static inline void LLtoXYZ_f(double lon, double lat, float *restrict xyz)
{
   double cos_lat = cos(lat);
   xyz[0] = cos_lat * cos(lon);
   xyz[1] = cos_lat * sin(lon);
   xyz[2] = sin(lat);
}

static inline void LLtoXYZ_kd(double lon, double lat, kdata_t *restrict xyz)
{
   double cos_lat = cos(lat);
   xyz[0] = KDATA_SCALE(cos_lat * cos(lon));
   xyz[1] = KDATA_SCALE(cos_lat * sin(lon));
   xyz[2] = KDATA_SCALE(sin(lat));
}

static
float square(const float x)
{
  return x*x;
}

static
float distance(const float *restrict a, const float *restrict b)
{
  return (square((a[0]-b[0]))+square((a[1]-b[1]))+square((a[2]-b[2])));
}


void gridsearch_set_method(const char *methodstr)
{
  if      ( strcmp(methodstr, "kdtree")  == 0 ) gridsearch_method_nn = GS_KDTREE;
  else if ( strcmp(methodstr, "nearpt3") == 0 ) gridsearch_method_nn = GS_NEARPT3;
  else if ( strcmp(methodstr, "full")    == 0 ) gridsearch_method_nn = GS_FULL;
  else
    cdoAbort("gridsearch method %s not available!\n", methodstr);
}


struct gridsearch *gridsearch_create_reg2d(unsigned lcyclic, unsigned nx, unsigned ny, const double *restrict lons, const double *restrict lats)
{
  struct gridsearch *gs = (struct gridsearch *) Calloc(1, sizeof(struct gridsearch));

  gs->nx = nx;
  gs->ny = ny;

  unsigned nxm = nx;
  if ( lcyclic ) nxm++;

  double *reg2d_center_lon = (double *) Malloc(nxm*sizeof(double));
  double *reg2d_center_lat = (double *) Malloc(ny*sizeof(double));

  memcpy(reg2d_center_lon, lons, nxm*sizeof(double));
  memcpy(reg2d_center_lat, lats, ny*sizeof(double));

  double *coslon = (double *) Malloc(nx*sizeof(double));
  double *sinlon = (double *) Malloc(nx*sizeof(double));
  double *coslat = (double *) Malloc(ny*sizeof(double));
  double *sinlat = (double *) Malloc(ny*sizeof(double));

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

  gs->search_radius = cdo_default_search_radius();

  return gs;
}


struct kdNode *gs_create_kdtree(unsigned n, const double *restrict lons, const double *restrict lats)
{
  struct kd_point *pointlist = (struct kd_point *) Malloc(n * sizeof(struct kd_point));  
  // see  example_cartesian.c
  if ( cdoVerbose ) printf("kdtree lib init 3D: n=%d  nthreads=%d\n", n, ompNumThreads);
  kdata_t min[3], max[3];
  min[0] = min[1] = min[2] =  1e9;
  max[0] = max[1] = max[2] = -1e9;
  kdata_t *restrict point;
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
  for ( unsigned i = 0; i < n; i++ ) 
    {
      point = pointlist[i].point;
      LLtoXYZ_kd(lons[i], lats[i], point);
      for ( unsigned j = 0; j < 3; ++j )
        {
          min[j] = point[j] < min[j] ? point[j] : min[j];
          max[j] = point[j] > max[j] ? point[j] : max[j];
        }
      pointlist[i].index = i;
    }

  struct kdNode *kdt = kd_buildTree(pointlist, n, min, max, 3, ompNumThreads);
  if ( pointlist ) Free(pointlist);

  return kdt;
}


void gs_destroy_nearpt3(struct gsNear *near)
{
  if ( near )
    {
#if defined(ENABLE_NEARPT3)
      if ( near->nearpt3 ) nearpt3_destroy(near->nearpt3);
#endif
      if ( near->pts )
        {
          Free(near->pts[0]);
          Free(near->pts);
        }

      Free(near);
    }
}


struct gsNear *gs_create_nearpt3(unsigned n, const double *restrict lons, const double *restrict lats)
{
  struct gsNear *near = (struct gsNear *) Calloc(1, sizeof(struct gsNear));

  Coord_T **p = (Coord_T **) Malloc(n*sizeof(Coord_T *));
  p[0] = (Coord_T *) Malloc(3*n*sizeof(Coord_T));
  for ( unsigned i = 1; i < n; i++ ) p[i] = p[0] + i*3;

  float point[3];

#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
  for ( unsigned i = 0; i < n; i++ )
    {
      LLtoXYZ_f(lons[i], lats[i], point);

      p[i][0] = NPT3SCALE(point[0]);
      p[i][1] = NPT3SCALE(point[1]);
      p[i][2] = NPT3SCALE(point[2]);
    }

  near->n = n;
  near->plons = lons;
  near->plats = lats;
  near->pts = p;
#if defined(ENABLE_NEARPT3)
  near->nearpt3 = nearpt3_preprocess(n, p);
#else
  cdoAbort("nearpt3 support not compiled in!");
#endif
  
  return near;
}


void gs_destroy_full(struct gsFull *full)
{
  if ( full )
    {
      if ( full->pts )
        {
          Free(full->pts[0]);
          Free(full->pts);
        }

      Free(full);
    }
}


struct gsFull *gs_create_full(unsigned n, const double *restrict lons, const double *restrict lats)
{
  struct gsFull *full = (struct gsFull *) Calloc(1, sizeof(struct gsFull));

  float **p = (float **) Malloc(n*sizeof(float *));
  p[0] = (float *) Malloc(3*n*sizeof(float));
  for ( unsigned i = 1; i < n; i++ ) p[i] = p[0] + i*3;

#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
  for ( unsigned i = 0; i < n; i++ )
    {
      LLtoXYZ_f(lons[i], lats[i], p[i]);
    }
  
  full->n = n;
  full->plons = lons;
  full->plats = lats;
  full->pts = p;

  return full;
}


struct gridsearch *gridsearch_create(unsigned n, const double *restrict lons, const double *restrict lats)
{
  struct gridsearch *gs = (struct gridsearch *) Calloc(1, sizeof(struct gridsearch));

  gs->n = n;
  if ( n == 0 ) return gs;

  gs->kdt  = gs_create_kdtree(n, lons, lats);

  gs->search_radius = cdo_default_search_radius();

  return gs;
}


struct gridsearch *gridsearch_create_nn(unsigned n, const double *restrict lons, const double *restrict lats)
{
  struct gridsearch *gs = (struct gridsearch *) Calloc(1, sizeof(struct gridsearch));

  gs->method_nn = gridsearch_method_nn;
  gs->n = n;
  if ( n == 0 ) return gs;

  if      ( gs->method_nn == GS_KDTREE  ) gs->kdt  = gs_create_kdtree(n, lons, lats);
  else if ( gs->method_nn == GS_NEARPT3 ) gs->near = gs_create_nearpt3(n, lons, lats);
  else if ( gs->method_nn == GS_FULL    ) gs->full = gs_create_full(n, lons, lats);

  gs->search_radius = cdo_default_search_radius();

  return gs;
}


void gridsearch_delete(struct gridsearch *gs)
{
  if ( gs )
    {
      if ( gs->kdt ) kd_destroyTree(gs->kdt);
      
      if ( gs->reg2d_center_lon ) Free(gs->reg2d_center_lon);
      if ( gs->reg2d_center_lat ) Free(gs->reg2d_center_lat);

      if ( gs->coslat ) Free(gs->coslat);
      if ( gs->coslon ) Free(gs->coslon);
      if ( gs->sinlat ) Free(gs->sinlat);
      if ( gs->sinlon ) Free(gs->sinlon);

      if ( gs->near ) gs_destroy_nearpt3(gs->near);
      if ( gs->full ) gs_destroy_full(gs->full);

      Free(gs);
    }
}

static
double gs_set_range(double *prange)
{
  double range;

  if ( prange )
    range = *prange;
  else
    range = SQR(2 * M_PI);     /* This has to be bigger than the presumed
                                * maximum distance to the NN but smaller
                                * than once around the sphere. The content
                                * of this variable is replaced with the
                                * distance to the NN squared. */
  return range;
}


kdNode *gs_nearest_kdtree(kdNode *kdt, double lon, double lat, double *prange)
{
  if ( kdt == NULL ) return NULL;
  
  float range0 = gs_set_range(prange);
  kdata_t range = KDATA_SCALE(range0);

  kdata_t point[3];
  LLtoXYZ_kd(lon, lat, point);

  kdNode *node = kd_nearest(kdt, point, &range, 3);

  float frange = KDATA_INVSCALE(range);
  if ( !(frange < range0) ) node = NULL;
  if ( prange ) *prange = frange;

  return node;
}


unsigned gs_nearest_nearpt3(struct gsNear *near, double lon, double lat, double *prange)
{
  unsigned index = GS_NOT_FOUND;
  if ( near == NULL ) return index;
  
#if defined(ENABLE_NEARPT3)
  float range0 = gs_set_range(prange);

  float point[3];
  LLtoXYZ_f(lon, lat, point);

  Coord_T q[3];
  q[0] = NPT3SCALE(point[0]);
  q[1] = NPT3SCALE(point[1]);
  q[2] = NPT3SCALE(point[2]);

  int closestpt = nearpt3_query(near->nearpt3, q);

  if ( closestpt >= 0 )
    {
      float point0[3];
      LLtoXYZ_f(near->plons[closestpt], near->plats[closestpt], point0);
      
      float range = distance(point, point0);
      if ( range < range0 )
        {
           index = (unsigned) closestpt;
           *prange = range;
        }
    }
#else
  UNUSED(lon);
  UNUSED(lat);
  UNUSED(prange);
#endif

  return index;
}


unsigned gs_nearest_full(struct  gsFull *full, double lon, double lat, double *prange)
{
  unsigned index = GS_NOT_FOUND;
  if ( full == NULL ) return index;
  
  float range0 = gs_set_range(prange);

  float point[3];
  LLtoXYZ_f(lon, lat, point);

  int n = full->n;
  float **pts = full->pts;
  int closestpt = -1;
  float dist = FLT_MAX;
  for ( int i = 0; i < n; i++ )
    {
      float d = distance(point, pts[i]);
      if ( closestpt < 0 || d < dist || (d<=dist && i < closestpt) )
        {
          dist = d;
          closestpt = i;
        }
    }

  if ( closestpt >= 0 )
    {
      if ( dist < range0 )
        {
          *prange = dist;
          index = (unsigned) closestpt;
        }
    }
  
  return index;
}


unsigned gridsearch_nearest(struct gridsearch *gs, double lon, double lat, double *prange)
{
  unsigned index = GS_NOT_FOUND;

  if ( gs )
    {
      if ( gs->method_nn == GS_KDTREE )
        {
          kdNode *node = gs_nearest_kdtree(gs->kdt, lon, lat, prange);
          if ( node ) index = (int) node->index;
        }
      else if ( gs->method_nn == GS_NEARPT3 )
        {
          index = gs_nearest_nearpt3(gs->near, lon, lat, prange);
        }
      else if ( gs->method_nn == GS_FULL )
        {
          index = gs_nearest_full(gs->full, lon, lat, prange);
        }
      else
        {
          cdoAbort("gridsearch_nearest::method_nn undefined!");
        }
    }

  return index;
}


struct pqueue *gridsearch_qnearest(struct gridsearch *gs, double lon, double lat, double *prange, unsigned nnn)
{
  if ( gs->kdt == NULL ) return NULL;
  
  kdata_t point[3];
  float range0 = gs_set_range(prange);
  kdata_t range = KDATA_SCALE(range0);
  struct pqueue *result = NULL;

  LLtoXYZ_kd(lon, lat, point);

  if ( gs )
    {
      result = kd_qnearest(gs->kdt, point, &range, nnn, 3);
      // printf("range %g %g %g %p\n", lon, lat, range, node);

      float frange = KDATA_INVSCALE(range);
      if ( !(frange < range0) )
        {
          if ( result )
            {
              struct resItem *p;
              while ( pqremove_min(result, &p) ) Free(p); // Free the result node taken from the heap
              Free(result->d); // free the heap
              Free(result);    // and free the heap information structure
            }
          result = NULL;
        }
      if ( prange ) *prange = frange;
    }
  
  return result;
}

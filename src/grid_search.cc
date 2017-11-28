#ifdef  HAVE_CONFIG_H
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
#include "kdtreelib/kdtree.h"
#include "nanoflann.hpp"
#include "nearpt3c.h"


#define  PI       M_PI
#define  PI2      (2.0*PI)


static int gridsearch_method_nn = GS_KDTREE;


struct gsFull {
  size_t n;
  const double *plons;
  const double *plats;
  float **pts;
};

struct gsNear {
  size_t n;
  const double *plons;
  const double *plats;
  float **pts;
  void *nearpt3;
};

template <typename T>
struct PointCloud
{
  struct Point { T  x,y,z; };
  std::vector<Point>  pts;
  T min[3], max[3];

  // Must return the number of data points
  inline size_t kdtree_get_point_count() const { return pts.size(); }

  // Returns the dim'th component of the idx'th point in the class:
  // Since this is inlined and the "dim" argument is typically an immediate value, the
  //  "if/else's" are actually solved at compile time.
  inline T kdtree_get_pt(const size_t idx, int dim) const
  {
    if (dim == 0) return pts[idx].x;
    else if (dim == 1) return pts[idx].y;
    else return pts[idx].z;
  }

  // Optional bounding-box computation: return false to default to a standard bbox computation loop.
  //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
  //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
  template <class BBOX>
  bool kdtree_get_bbox(BBOX& bb) const
  {
    for ( unsigned j = 0; j < 3; ++j )
      {
        bb[j].low  = min[j];
        bb[j].high = max[j];
      }
    return true;
  }
  // bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
};

typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<float, PointCloud<float> > ,
    PointCloud<float>,
    3 /* dim */
    > nfTree_t;


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

static inline
void LLtoXYZ_f(double lon, double lat, float *restrict xyz)
{
   double cos_lat = cos(lat);
   xyz[0] = cos_lat * cos(lon);
   xyz[1] = cos_lat * sin(lon);
   xyz[2] = sin(lat);
}

static inline
void LLtoXYZ_kd(double lon, double lat, kdata_t *restrict xyz)
{
   double cos_lat = cos(lat);
   xyz[0] = KDATA_SCALE(cos_lat * cos(lon));
   xyz[1] = KDATA_SCALE(cos_lat * sin(lon));
   xyz[2] = KDATA_SCALE(sin(lat));
}

static constexpr
float square(const float x)
{
  return x*x;
}

static constexpr
float distance(const float *restrict a, const float *restrict b)
{
  return (square((a[0]-b[0]))+square((a[1]-b[1]))+square((a[2]-b[2])));
}


void gridsearch_set_method(const char *methodstr)
{
  if      ( strcmp(methodstr, "kdtree")    == 0 ) gridsearch_method_nn = GS_KDTREE;
  else if ( strcmp(methodstr, "nanoflann") == 0 ) gridsearch_method_nn = GS_NANOFLANN;
  else if ( strcmp(methodstr, "kdsph")     == 0 ) gridsearch_method_nn = GS_KDSPH;
  else if ( strcmp(methodstr, "nearpt3")   == 0 ) gridsearch_method_nn = GS_NEARPT3;
  else if ( strcmp(methodstr, "full")      == 0 ) gridsearch_method_nn = GS_FULL;
  else
    cdoAbort("gridsearch method %s not available!", methodstr);
}


void gridsearch_extrapolate(struct gridsearch *gs)
{
  gs->extrapolate = true;
}


struct gridsearch *gridsearch_create_reg2d(bool is_cyclic, size_t dims[2], const double *restrict lons, const double *restrict lats)
{
  struct gridsearch *gs = (struct gridsearch *) Calloc(1, sizeof(struct gridsearch));

  gs->is_cyclic = is_cyclic;
  gs->is_reg2d = true;
  gs->dims[0] = dims[0];
  gs->dims[1] = dims[1];
  size_t nx = dims[0];
  size_t ny = dims[0];

  size_t nxm = nx;
  if ( is_cyclic ) nxm++;

  double *reg2d_center_lon = (double *) Malloc(nxm*sizeof(double));
  double *reg2d_center_lat = (double *) Malloc(ny*sizeof(double));

  memcpy(reg2d_center_lon, lons, nxm*sizeof(double));
  memcpy(reg2d_center_lat, lats, ny*sizeof(double));

  double *coslon = (double *) Malloc(nx*sizeof(double));
  double *sinlon = (double *) Malloc(nx*sizeof(double));
  double *coslat = (double *) Malloc(ny*sizeof(double));
  double *sinlat = (double *) Malloc(ny*sizeof(double));

  for ( size_t n = 0; n < nx; ++n )
    {
      double rlon = lons[n];
      if ( rlon > PI2 ) rlon -= PI2;
      if ( rlon < 0   ) rlon += PI2;
      coslon[n] = cos(rlon);
      sinlon[n] = sin(rlon);
    }
  for ( size_t n = 0; n < ny; ++n )
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

static
void *gs_create_kdtree(size_t n, const double *restrict lons, const double *restrict lats, struct gridsearch *gs)
{
  struct kd_point *pointlist = (struct kd_point *) Malloc(n*sizeof(struct kd_point));  
  // see  example_cartesian.c
  if ( cdoVerbose ) printf("kdtree lib init 3D: n=%zu  nthreads=%d\n", n, ompNumThreads);
  kdata_t min[3], max[3];
  min[0] = min[1] = min[2] =  1e9;
  max[0] = max[1] = max[2] = -1e9;
#if defined(HAVE_OPENMP4)
  // failed with INTEL CC: error: 'min' has invalid type for 'reduction'
  // #pragma omp parallel for reduction(min: min) reduction(max: max)
#endif
  for ( size_t i = 0; i < n; i++ ) 
    {
      kdata_t *restrict point = pointlist[i].point;
      LLtoXYZ_kd(lons[i], lats[i], point);
      for ( unsigned j = 0; j < 3; ++j )
        {
          min[j] = point[j] < min[j] ? point[j] : min[j];
          max[j] = point[j] > max[j] ? point[j] : max[j];
        }
      pointlist[i].index = i;
    }

  for ( unsigned j = 0; j < 3; ++j )
    {
      gs->min[j] = min[j];
      gs->max[j] = max[j];
    }

  kdTree_t *kdt = kd_buildTree(pointlist, n, min, max, 3, ompNumThreads);
  if ( pointlist ) Free(pointlist);
  if ( kdt == NULL ) cdoAbort("kd_buildTree failed!");

  return (void*) kdt;
}

static
void *gs_create_nanoflann(size_t n, const double *restrict lons, const double *restrict lats, struct gridsearch *gs)
{
  PointCloud<float> *pointcloud = new PointCloud<float>();
  if ( cdoVerbose ) printf("nanoflann init 3D: n=%zu  nthreads=%d\n", n, ompNumThreads);

  float min[3], max[3];
  min[0] = min[1] = min[2] =  1e9;
  max[0] = max[1] = max[2] = -1e9;
  // Generating  Point Cloud
  pointcloud->pts.resize(n);
#if defined(HAVE_OPENMP4)
  // failed with INTEL CC: error: 'min' has invalid type for 'reduction'
  // #pragma omp parallel for reduction(min: min) reduction(max: max)
#endif
  for ( size_t i = 0; i < n; i++ ) 
    {
      float point[3];
      LLtoXYZ_f(lons[i], lats[i], point);
      pointcloud->pts[i].x = point[0];
      pointcloud->pts[i].y = point[1];
      pointcloud->pts[i].z = point[2];
      for ( unsigned j = 0; j < 3; ++j )
        {
          min[j] = point[j] < min[j] ? point[j] : min[j];
          max[j] = point[j] > max[j] ? point[j] : max[j];
        }
    }

  for ( unsigned j = 0; j < 3; ++j )
    {
      gs->min[j] = min[j];
      gs->max[j] = max[j];
    }

  for ( unsigned j = 0; j < 3; ++j )
    {
      pointcloud->min[j] = min[j];
      pointcloud->max[j] = max[j];
    }

  // construct a kd-tree index:
  nfTree_t *nft = new nfTree_t(3 /*dim*/, *pointcloud, nanoflann::KDTreeSingleIndexAdaptorParams(50 /* max leaf */) );
  nft->buildIndex();

  return (void*)nft;
}

static
void *gs_create_kdsph(size_t n, const double *restrict lons, const double *restrict lats)
{
  struct kd_point *pointlist = (struct kd_point *) Malloc(n*sizeof(struct kd_point)); // kd_point contains 3d point
  // see  example_cartesian.c
  if ( cdoVerbose ) printf("kdtree lib spherical init: n=%zu  nthreads=%d\n", n, ompNumThreads);
  kdata_t min[2], max[2];
  min[0] = min[1] =  1e9;
  max[0] = max[1] = -1e9;
  kdata_t *restrict point;
#if defined(HAVE_OPENMP4)
  //#pragma omp simd
#endif
  for ( size_t i = 0; i < n; i++ ) 
    {
      point = pointlist[i].point;
      point[0] = lons[i];
      point[1] = lats[i];
      point[2] = 0; // dummy
      for ( unsigned j = 0; j < 2; ++j )
        {
          min[j] = point[j] < min[j] ? point[j] : min[j];
          max[j] = point[j] > max[j] ? point[j] : max[j];
        }
      pointlist[i].index = i;
    }

  kdTree_t *kdt = kd_sph_buildTree(pointlist, n, min, max, ompNumThreads);
  if ( pointlist ) Free(pointlist);

  return (void*)kdt;
}

static
void gs_destroy_kdtree(void *search_container)
{
  kdTree_t *kdt = (kdTree_t *) search_container;
  if ( kdt ) kd_destroyTree(kdt);
}

static
void gs_destroy_nearpt3(void *search_container)
{
  struct gsNear *near = (struct gsNear*) search_container;
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

static
void *gs_create_nearpt3(size_t n, const double *restrict lons, const double *restrict lats)
{
  struct gsNear *near = (struct gsNear *) Calloc(1, sizeof(struct gsNear));

  Coord_T **p = (Coord_T **) Malloc(n*sizeof(Coord_T *));
  p[0] = (Coord_T *) Malloc(3*n*sizeof(Coord_T));
  for ( size_t i = 1; i < n; i++ ) p[i] = p[0] + i*3;

  float point[3];

#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
  for ( size_t i = 0; i < n; i++ )
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
  
  return (void*)near;
}

static
void gs_destroy_full(void *search_container)
{
  struct gsFull *full = (struct gsFull *) search_container;
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

static
void *gs_create_full(size_t n, const double *restrict lons, const double *restrict lats)
{
  struct gsFull *full = (struct gsFull *) Calloc(1, sizeof(struct gsFull));

  float **p = (float **) Malloc(n*sizeof(float *));
  p[0] = (float *) Malloc(3*n*sizeof(float));
  for ( size_t i = 1; i < n; i++ ) p[i] = p[0] + i*3;

#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
  for ( size_t i = 0; i < n; i++ )
    {
      LLtoXYZ_f(lons[i], lats[i], p[i]);
    }
  
  full->n = n;
  full->plons = lons;
  full->plats = lats;
  full->pts = p;

  return (void*) full;
}

static
void cal_bound_box(struct gridsearch *gs, size_t n, const double *restrict lons, const double *restrict lats)
{
  enum {mlon=180, mlat=90};
  const size_t masksize = mlon*mlat;
  bool *mask1 = (bool*) Malloc(masksize*sizeof(bool));
  for ( size_t i = 0; i < masksize; ++i ) mask1[i] = 0;

  for ( size_t i = 0; i < n; ++i )
    {
      double fi = lons[i];
      double fj = lats[i]+0.5*M_PI;
      if ( fi < 0 ) fi += 2*M_PI;
      fi *= mlon/(2*M_PI);
      fj *= mlat/M_PI;
      int64_t ii = (int64_t) llround(fi);
      int64_t jj = (int64_t) llround(fj);
      if ( i < 1000 )
        printf("lon=%g/lat=%g  ii=%ld  jj=%ld  fi=%g fj=%g\n", lons[i]*RAD2DEG, lats[i]*RAD2DEG, (long)ii, (long)jj, fi, fj);
      mask1[jj*mlon+ii] = 1;
    }

  double xlonmin = 999., xlonmax = 999.;
  double xlatmin = 999., xlatmax = 999.;
  for ( size_t j = 0; j < mlat; ++j )
    {
      for ( size_t i = 0; i < mlon; ++i )
        if ( mask1[j*mlon+i] )
          {
            xlatmin = (j+0.)*M_PI/mlat - 0.5*M_PI;
            break;
          }
      if ( IS_NOT_EQUAL(xlatmin, 999.) ) break;
    }
  for ( size_t jj = 0; jj < mlat; ++jj )
    {
      size_t j = mlat-jj-1;
      for ( size_t i = 0; i < mlon; ++i )
        if ( mask1[j*mlon+i] )
          {
            xlatmax = (j+1.)*M_PI/mlat - 0.5*M_PI;
            break;
          }
      if ( IS_NOT_EQUAL(xlatmax, 999.) ) break;
    }

  size_t ii = mlon;
  for ( size_t i = 0; i < mlon; ++i )
    {
      for ( size_t j = 0; j < mlat; ++j )
        if ( mask1[j*mlon+i] )
          {
            ii = i;
            xlonmin = (i+.0)*2*M_PI/mlon;
            break;
          }
      if ( IS_NOT_EQUAL(xlonmin, 999.) ) break;
    }
  if ( cdoVerbose ) printf("ii= %zu\n", ii);
  for ( size_t i = ii; i < mlon; ++i )
    {
      ii = mlon;
      size_t imask = 0;
      for ( size_t j = 0; j < mlat; ++j )
        {
          if ( mask1[j*mlon+i] )
            {
              imask++;
              xlonmax = (i+1.)*2*M_PI/mlon;
            }
        }
      if ( imask == 0 )
        {
          ii = i+1;
          break;
        }
    }
  if ( cdoVerbose ) cdoPrint("boundbox: lonmin=%g lonmax=%g latmin=%g latmax=%g",
                             xlonmin*RAD2DEG, xlonmax*RAD2DEG, xlatmin*RAD2DEG, xlatmax*RAD2DEG);
  if ( cdoVerbose ) printf("ii= %zu\n", ii);
  for ( size_t i = ii; i < mlon; ++i )
    {
      bool lfound = false;
      for ( size_t j = 0; j < mlat; ++j )
        if ( mask1[j*mlon+i] )
          {
            lfound = true;
            xlonmin = (i+.0)*2*M_PI/mlon;
            break;
          }
      if ( lfound ) break;
    }

  if ( xlonmin > xlonmax ) xlonmin -= 2*M_PI;

  if ( cdoVerbose ) cdoPrint("boundbox: lonmin=%g lonmax=%g latmin=%g latmax=%g",
                             xlonmin*RAD2DEG, xlonmax*RAD2DEG, xlatmin*RAD2DEG, xlatmax*RAD2DEG);

  for ( size_t j = 0; j < mlat; ++j )
    {
      for ( size_t i = 0; i < mlon; ++i )
        {
           printf("%1d", mask1[j*mlon+i]);
        }
        printf("\n");
    }

  double lonmin = 1.e33, lonmax = -1.e33;
  double latmin = 1.e33, latmax = -1.e33;
  for ( size_t i = 0; i < n; ++i )
    {
      if ( lons[i] < lonmin ) lonmin = lons[i];
      if ( lons[i] > lonmax ) lonmax = lons[i];
      if ( lats[i] < latmin ) latmin = lats[i];
      if ( lats[i] > latmax ) latmax = lats[i];
    }

  gs->lonmin = xlonmin;
  gs->lonmax = xlonmax;
  gs->latmin = xlatmin;
  gs->latmax = xlatmax;
  if ( cdoVerbose ) cdoPrint("boundbox: lonmin=%g lonmax=%g latmin=%g latmax=%g",
                             lonmin*RAD2DEG, lonmax*RAD2DEG, latmin*RAD2DEG, latmax*RAD2DEG);

  if ( mask1 ) free(mask1);
}

static
void cal_mask(struct gridsearch *gs)
{
  enum {mmin = 100};
  const double fact = 1;
  double dlon = gs->lonmax - gs->lonmin;
  double dlat = gs->latmax - gs->latmin;
  if ( cdoVerbose ) cdoPrint("mask: dlon=%g, dlat=%g", dlon*RAD2DEG, dlat*RAD2DEG);
  size_t sqrtn = (size_t) sqrt((double)gs->n);
  if ( cdoVerbose ) cdoPrint("n=%zu  sqrt(n)=%zu", gs->n, sqrtn);
  size_t mlat = fact*(sqrtn/2.)*(1+dlat/dlon);
  size_t mlon = mlat*dlon/dlat;
  if ( cdoVerbose ) cdoPrint("mlon=%zu mlat=%zu mn=%zu", mlon, mlat, mlon*mlat);
}


struct gridsearch *gridsearch_create(size_t n, const double *restrict lons, const double *restrict lats)
{
  struct gridsearch *gs = (struct gridsearch *) Calloc(1, sizeof(struct gridsearch));

  gs->method_nn = gridsearch_method_nn;
  gs->n = n;
  if ( n == 0 ) return gs;

  if      ( gs->method_nn == GS_KDTREE    ) gs->search_container = gs_create_kdtree(n, lons, lats, gs);
  else if ( gs->method_nn == GS_NANOFLANN ) gs->search_container = gs_create_nanoflann(n, lons, lats, gs);

  gs->search_radius = cdo_default_search_radius();

  return gs;
}


struct gridsearch *gridsearch_create_nn(size_t n, const double *restrict lons, const double *restrict lats)
{
  struct gridsearch *gs = (struct gridsearch *) Calloc(1, sizeof(struct gridsearch));

  gs->method_nn = gridsearch_method_nn;
  gs->n = n;
  if ( n == 0 ) return gs;

  if      ( gs->method_nn == GS_KDTREE    ) gs->search_container = gs_create_kdtree(n, lons, lats, gs);
  else if ( gs->method_nn == GS_NANOFLANN ) gs->search_container = gs_create_nanoflann(n, lons, lats, gs);
  else if ( gs->method_nn == GS_KDSPH     ) gs->search_container = gs_create_kdsph(n, lons, lats);
  else if ( gs->method_nn == GS_NEARPT3   ) gs->search_container = gs_create_nearpt3(n, lons, lats);
  else if ( gs->method_nn == GS_FULL      ) gs->search_container = gs_create_full(n, lons, lats);

  gs->search_radius = cdo_default_search_radius();

  // cal_bound_box(gs, n, lons, lats);
  // cal_mask(gs);

  return gs;
}


void gridsearch_delete(struct gridsearch *gs)
{
  if ( gs )
    {      
      if ( gs->reg2d_center_lon ) Free(gs->reg2d_center_lon);
      if ( gs->reg2d_center_lat ) Free(gs->reg2d_center_lat);

      if ( gs->coslat ) Free(gs->coslat);
      if ( gs->coslon ) Free(gs->coslon);
      if ( gs->sinlat ) Free(gs->sinlat);
      if ( gs->sinlon ) Free(gs->sinlon);

      if      ( gs->method_nn == GS_KDTREE    ) gs_destroy_kdtree(gs->search_container );
      else if ( gs->method_nn == GS_NANOFLANN ) ;
      else if ( gs->method_nn == GS_KDSPH     ) gs_destroy_kdtree(gs->search_container);
      else if ( gs->method_nn == GS_NEARPT3   ) gs_destroy_nearpt3(gs->search_container);
      else if ( gs->method_nn == GS_FULL      ) gs_destroy_full(gs->search_container);

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

static
size_t gs_nearest_kdtree(void *search_container, double lon, double lat, double *prange, struct gridsearch *gs)
{
  size_t index = GS_NOT_FOUND;
  kdTree_t *kdt = (kdTree_t *) search_container;
  if ( kdt == NULL ) return index;
  
  float range0 = gs_set_range(prange);
  kdata_t range = KDATA_SCALE(range0);

  kdata_t query_pt[3];
  LLtoXYZ_kd(lon, lat, query_pt);

  if ( !gs->extrapolate )
    for ( unsigned j = 0; j < 3; ++j )
      if ( query_pt[j] < gs->min[j] || query_pt[j] > gs->max[j] ) return index;

  kdNode *node = kd_nearest(kdt->node, query_pt, &range, 3);

  float frange = KDATA_INVSCALE(range);
  if ( !(frange < range0) ) node = NULL;
  if ( prange ) *prange = frange;

  if ( node ) index = node->index;
  //printf("%zu %g\n", index, range);

  return index;
}

static
size_t gs_nearest_nanoflann(void *search_container, double lon, double lat, double *prange, struct gridsearch *gs)
{
  size_t index = GS_NOT_FOUND;
  nfTree_t *nft = (nfTree_t *) search_container;
  if ( nft == NULL ) return index;
  
  float range0 = gs_set_range(prange);
  float range = range0;

  float query_pt[3];
  LLtoXYZ_f(lon, lat, query_pt);

  if ( !gs->extrapolate )
    for ( unsigned j = 0; j < 3; ++j )
      if ( query_pt[j] < gs->min[j] || query_pt[j] > gs->max[j] ) return index;

  const size_t num_results = 1;
  size_t ret_index;
  float out_dist_sqr;
  nanoflann::KNNResultSet<float> resultSet(range, num_results);
  resultSet.init(&ret_index, &out_dist_sqr);
  nft->findNeighbors(resultSet, query_pt, nanoflann::SearchParams(10));
  //printf("%zu %g\n", ret_index, out_dist_sqr);

  index = ret_index;
  *prange = out_dist_sqr;
  //float frange = range;
  //if ( !(frange < range0) ) node = NULL;
  //if ( prange ) *prange = frange;

  //if ( node ) index = node->index;

  return index;
}

static
size_t gs_nearest_kdsph(void *search_container, double lon, double lat, double *prange)
{
  size_t index = GS_NOT_FOUND;
  kdTree_t *kdt = (kdTree_t *) search_container;
  if ( kdt == NULL ) return index;
  
  float range0 = gs_set_range(prange);
  kdata_t range = KDATA_SCALE(range0);

  kdata_t query_pt[2];
  query_pt[0] = lon;
  query_pt[1] = lat;

  kdNode *node = kd_nearest(kdt->node, query_pt, &range, 3);

  float frange = KDATA_INVSCALE(range);
  if ( !(frange < range0) ) node = NULL;
  if ( prange ) *prange = frange;

  if ( node ) index = node->index;

  return index;
}

static
size_t gs_nearest_nearpt3(void *search_container, double lon, double lat, double *prange)
{
  size_t index = GS_NOT_FOUND;
  struct gsNear *near = (struct gsNear *) search_container;
  if ( near == NULL ) return index;
  
#if defined(ENABLE_NEARPT3)
  float range0 = gs_set_range(prange);

  float query_pt[3];
  LLtoXYZ_f(lon, lat, query_pt);

  Coord_T q[3];
  q[0] = NPT3SCALE(query_pt[0]);
  q[1] = NPT3SCALE(query_pt[1]);
  q[2] = NPT3SCALE(query_pt[2]);

  int closestpt = nearpt3_query(near->nearpt3, q);

  if ( closestpt >= 0 )
    {
      float query_pt0[3];
      LLtoXYZ_f(near->plons[closestpt], near->plats[closestpt], query_pt0);
      
      float range = distance(query_pt, query_pt0);
      if ( range < range0 )
        {
          index = (size_t) closestpt;
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

static
size_t gs_nearest_full(void *search_container, double lon, double lat, double *prange)
{
  size_t index = GS_NOT_FOUND;
  struct  gsFull *full = (struct  gsFull *) search_container;
  if ( full == NULL ) return index;
  
  float range0 = gs_set_range(prange);

  float query_pt[3];
  LLtoXYZ_f(lon, lat, query_pt);

  size_t n = full->n;
  float **pts = full->pts;
  size_t closestpt = n;
  float dist = FLT_MAX;
  for ( size_t i = 0; i < n; i++ )
    {
      float d = distance(query_pt, pts[i]);
      if ( closestpt >=n || d < dist || (d<=dist && i < closestpt) )
        {
          dist = d;
          closestpt = i;
        }
    }

  if ( closestpt < n )
    {
      if ( dist < range0 )
        {
          *prange = dist;
          index = closestpt;
        }
    }
  
  return index;
}


size_t gridsearch_nearest(struct gridsearch *gs, double lon, double lat, double *prange)
{
  size_t index = GS_NOT_FOUND;

  if ( gs )
    {
      void *sc = gs->search_container;
      // clang-format off
      if      ( gs->method_nn == GS_KDTREE )    index = gs_nearest_kdtree(sc, lon, lat, prange, gs);
      else if ( gs->method_nn == GS_NANOFLANN ) index = gs_nearest_nanoflann(sc, lon, lat, prange, gs);
      else if ( gs->method_nn == GS_KDSPH )     index = gs_nearest_kdsph(sc, lon, lat, prange);
      else if ( gs->method_nn == GS_NEARPT3 )   index = gs_nearest_nearpt3(sc, lon, lat, prange);
      else if ( gs->method_nn == GS_FULL )      index = gs_nearest_full(sc, lon, lat, prange);
      else cdoAbort("%s::method_nn undefined!", __func__);
      // clang-format on
    }

  return index;
}

static
size_t gs_qnearest_kdtree(struct gridsearch *gs, double lon, double lat, double *prange, size_t nnn, size_t *adds, double *dist)
{
  size_t nadds = 0;

  kdTree_t *kdt = (kdTree_t *) gs->search_container;
  if ( kdt == NULL ) return nadds;
  
  kdata_t query_pt[3];
  float range0 = gs_set_range(prange);
  kdata_t range = KDATA_SCALE(range0);
  struct pqueue *result = NULL;

  LLtoXYZ_kd(lon, lat, query_pt);

  if ( gs )
    {
      result = kd_qnearest(kdt->node, query_pt, &range, nnn, 3);
      // printf("range %g %g %g %p\n", lon, lat, range, node);

      float frange = KDATA_INVSCALE(range);

      if ( result )
        {
          size_t index;
          struct resItem *p;
          while ( pqremove_min(result, &p) )
            {
              index = p->node->index;
              range = p->dist_sq;
              Free(p); // Free the result node taken from the heap

              if ( range < range0 )
                {
                  dist[nadds] = range;
                  adds[nadds] = index;
                  nadds++;
                }
            }
          Free(result->d); // free the heap
          Free(result);    // and free the heap information structure
        }

      if ( prange ) *prange = frange;
    }
  
  return nadds;
}

static
size_t gs_qnearest_nanoflann(struct gridsearch *gs, double lon, double lat, double *prange, size_t nnn, size_t *adds, double *dist)
{
  size_t nadds = 0;

  nfTree_t *nft = (nfTree_t*) gs->search_container;
  if ( nft == NULL ) return nadds;
  
  float range0 = gs_set_range(prange);
  float range = range0;

  float query_pt[3];
  LLtoXYZ_f(lon, lat, query_pt);

  if ( gs )
    {
      std::vector<float> out_dist_sqr(nnn);
      nadds = nft->knnRangeSearch(&query_pt[0], range, nnn, &adds[0], &out_dist_sqr[0]);

      for ( size_t i = 0; i < nadds; ++i ) dist[i] = out_dist_sqr[i];

      float frange = range;
      if ( prange ) *prange = frange;
    }
  
  return nadds;
}


size_t gridsearch_qnearest(struct gridsearch *gs, double lon, double lat, double *prange, size_t nnn, size_t *adds, double *dist)
{
  size_t nadds = 0;

  if ( gs )
    {
      // clang-format off
      if      ( gs->method_nn == GS_KDTREE )    nadds = gs_qnearest_kdtree(gs, lon, lat, prange, nnn, adds, dist);
      else if ( gs->method_nn == GS_NANOFLANN ) nadds = gs_qnearest_nanoflann(gs, lon, lat, prange, nnn, adds, dist);
      else cdoAbort("%s::method_nn undefined!", __func__);
      // clang-format on
    }

  return nadds;
}

#define  BIGNUM   1.e+20
#define  TINY     1.e-14

static
void knn_store_distance(size_t nadd, double distance, size_t num_neighbors, size_t *restrict nbr_add, double *restrict nbr_dist)
{
  if ( num_neighbors == 1 )
    {
      if ( distance < nbr_dist[0] || (distance <= nbr_dist[0] && nadd < nbr_add[0]) )
	{
	  nbr_add[0]  = nadd;
	  nbr_dist[0] = distance;
	}
    }
  else
    {
      for ( size_t nchk = 0; nchk < num_neighbors; ++nchk )
	{
	  if ( distance < nbr_dist[nchk] || (distance <= nbr_dist[nchk] && nadd < nbr_add[nchk]) )
	    {
	      for ( size_t n = num_neighbors-1; n > nchk; --n )
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
void knn_check_distance(size_t num_neighbors, const size_t *restrict nbr_add, double *restrict nbr_dist)
{
  // If distance is zero, set to small number
  for ( size_t nchk = 0; nchk < num_neighbors; ++nchk )
    if ( nbr_add[nchk] != GS_NOT_FOUND && nbr_dist[nchk] <= 0. ) nbr_dist[nchk] = TINY;
}


void gridsearch_knn_init(struct gsknn *knn)
{
  size_t ndist = knn->ndist;
  size_t *restrict add = knn->add;
  double *restrict dist = knn->dist;

  for ( size_t i = 0; i < ndist; ++i )
    {
      add[i]  = GS_NOT_FOUND;
      dist[i] = BIGNUM;
    }
}


struct gsknn *gridsearch_knn_new(size_t size)
{
  struct gsknn *knn = (struct gsknn *) Malloc(sizeof(struct gsknn));
  
  knn->ndist   = size;
  knn->size    = size;
  knn->mask    = (bool*) Malloc(size*sizeof(bool));     // mask at nearest neighbors
  knn->add     = (size_t*) Malloc(size*sizeof(size_t)); // source address at nearest neighbors
  knn->dist    = (double*) Malloc(size*sizeof(double)); // angular distance of the nearest neighbors
  knn->tmpadd  = NULL;
  knn->tmpdist = NULL;

  gridsearch_knn_init(knn);

  return knn;
}


void gridsearch_knn_delete(struct gsknn *knn)
{
  if ( knn )
    {
      knn->size = 0;
      if ( knn->dist    ) Free(knn->dist);
      if ( knn->add     ) Free(knn->add);
      if ( knn->tmpdist ) Free(knn->tmpdist);
      if ( knn->tmpadd  ) Free(knn->tmpadd);
      if ( knn->mask    ) Free(knn->mask);
      Free(knn);
    }
}


size_t gridsearch_knn(struct gridsearch *gs, struct gsknn *knn, double plon, double plat)
{
  /*
    Output variables:

    int nbr_add[num_neighbors]     ! address of each of the closest points
    double nbr_dist[num_neighbors] ! distance to each of the closest points

    Input variables:

    double plat,         ! latitude  of the search point
    double plon,         ! longitude of the search point
  */

  double search_radius = gs->search_radius;

  // Initialize distance and address arrays
  gridsearch_knn_init(knn);

  size_t num_neighbors = knn->size;
  size_t *restrict nbr_add = knn->add;
  double *restrict nbr_dist = knn->dist;

  size_t ndist = num_neighbors;
  // check some more points if distance is the same use the smaller index (nadd)
  if ( ndist > 8 ) ndist += 8;
  else             ndist *= 2;
  if ( ndist > gs->n ) ndist = gs->n;

  if ( knn->tmpadd  == NULL ) knn->tmpadd  = (size_t*) Malloc(ndist*sizeof(size_t));
  if ( knn->tmpdist == NULL ) knn->tmpdist = (double*) Malloc(ndist*sizeof(double));

  size_t *adds = knn->tmpadd;
  double *dist = knn->tmpdist;
  
  const double range0 = SQR(search_radius);
  double range = range0;

  size_t nadds = 0;

  if ( num_neighbors == 1 )
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
  size_t max_neighbors = (ndist < num_neighbors) ? ndist : num_neighbors;

  for ( size_t i = 0; i < ndist; ++i )
    knn_store_distance(adds[i], dist[i], max_neighbors, nbr_add, nbr_dist);

  knn_check_distance(max_neighbors, nbr_add, nbr_dist);

  if ( ndist > num_neighbors ) ndist = num_neighbors;

  knn->ndist = ndist;

  return ndist;
} // gridsearch_knn

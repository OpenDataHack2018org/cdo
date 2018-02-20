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
#ifdef HAVE_CONFIG_H
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
#include "cdoOptions.h"

#define PI M_PI
#define PI2 (2.0 * PI)

static GridsearchMethod gridsearch_method_nn(GridsearchMethod::kdtree);

struct gsFull
{
  size_t n;
  const double *plons;
  const double *plats;
  double **pts;
};

template <typename T>
struct PointCloud
{
  struct Point
  {
    T x, y, z;
  };
  std::vector<Point> pts;
  T min[3], max[3];

  // Must return the number of data points
  inline size_t
  kdtree_get_point_count() const
  {
    return pts.size();
  }

  // Returns the dim'th component of the idx'th point in the class:
  // Since this is inlined and the "dim" argument is typically an immediate
  // value, the
  //  "if/else's" are actually solved at compile time.
  inline T
  kdtree_get_pt(const size_t idx, int dim) const
  {
    if (dim == 0)
      return pts[idx].x;
    else if (dim == 1)
      return pts[idx].y;
    else
      return pts[idx].z;
  }

  // Optional bounding-box computation: return false to default to a standard
  // bbox computation loop.
  //   Return true if the BBOX was already computed by the class and returned in
  //   "bb" so it can be avoided to redo it again. Look at bb.size() to find out
  //   the expected dimensionality (e.g. 2 or 3 for point clouds)
  template <class BBOX>
  bool
  kdtree_get_bbox(BBOX &bb) const
  {
    for (unsigned j = 0; j < 3; ++j)
      {
        bb[j].low = min[j];
        bb[j].high = max[j];
      }
    return true;
  }
  // bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
};

typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>,
    PointCloud<double>, 3 /* dim */
    >
    nfTree_t;

static double
cdo_default_search_radius(void)
{
  extern double gridsearch_radius;

  double search_radius = gridsearch_radius;

  if (search_radius < 0.) search_radius = 0.;
  if (search_radius > 180.) search_radius = 180.;

  search_radius = search_radius * DEG2RAD;

  return search_radius;
}

static inline void
LLtoXYZ(double lon, double lat, double *restrict xyz)
{
  double cos_lat = cos(lat);
  xyz[0] = cos_lat * cos(lon);
  xyz[1] = cos_lat * sin(lon);
  xyz[2] = sin(lat);
}

static constexpr double
square(const double x) noexcept
{
  return x * x;
}

static constexpr double
distance(const double *restrict a, const double *restrict b) noexcept
{
  return square(a[0] - b[0]) + square(a[1] - b[1]) + square(a[2] - b[2]);
}

void
gridsearch_set_method(const char *methodstr)
{
  if (strcmp(methodstr, "kdtree") == 0)
    gridsearch_method_nn = GridsearchMethod::kdtree;
  else if (strcmp(methodstr, "nanoflann") == 0)
    gridsearch_method_nn = GridsearchMethod::nanoflann;
  else if (strcmp(methodstr, "full") == 0)
    gridsearch_method_nn = GridsearchMethod::full;
  else
    cdoAbort("gridsearch method %s not available!", methodstr);
}

void
gridsearch_extrapolate(struct gridsearch *gs)
{
  gs->extrapolate = true;
}

struct gridsearch *
gridsearch_create_reg2d(bool xIsCyclic, size_t dims[2],
                        const double *restrict lons,
                        const double *restrict lats)
{
  struct gridsearch *gs
      = (struct gridsearch *) Calloc(1, sizeof(struct gridsearch));

  gs->is_cyclic = xIsCyclic;
  gs->is_reg2d = true;
  gs->dims[0] = dims[0];
  gs->dims[1] = dims[1];
  size_t nx = dims[0];
  size_t ny = dims[0];

  size_t nxm = xIsCyclic ? nx + 1 : nx;

  double *reg2d_center_lon = (double *) Malloc(nxm * sizeof(double));
  double *reg2d_center_lat = (double *) Malloc(ny * sizeof(double));

  arrayCopy(nxm, lons, reg2d_center_lon);
  arrayCopy(ny, lats, reg2d_center_lat);

  double *coslon = (double *) Malloc(nx * sizeof(double));
  double *sinlon = (double *) Malloc(nx * sizeof(double));
  double *coslat = (double *) Malloc(ny * sizeof(double));
  double *sinlat = (double *) Malloc(ny * sizeof(double));

  for (size_t n = 0; n < nx; ++n)
    {
      double rlon = lons[n];
      if (rlon > PI2) rlon -= PI2;
      if (rlon < 0) rlon += PI2;
      coslon[n] = cos(rlon);
      sinlon[n] = sin(rlon);
    }

  for (size_t n = 0; n < ny; ++n)
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

static void *
gs_create_kdtree(size_t n, const double *restrict lons,
                 const double *restrict lats, struct gridsearch *gs)
{
  struct kd_point *pointlist
      = (struct kd_point *) Malloc(n * sizeof(struct kd_point));
  // see  example_cartesian.c
  if (cdoVerbose)
    printf("kdtree lib init 3D: n=%zu  nthreads=%d\n", n,
           Threading::ompNumThreads);

  kdata_t min[3] = { 1.e9, 1.e9, 1.e9 };
  kdata_t max[3] = { -1.e9, -1.e9, -1.e9 };

#ifdef HAVE_OPENMP45
#pragma omp parallel for reduction(min : min[:3]) reduction(max : max[:3])
#endif
  for (size_t i = 0; i < n; i++)
    {
      kdata_t *restrict point = pointlist[i].point;
      LLtoXYZ(lons[i], lats[i], point);
      for (unsigned j = 0; j < 3; ++j)
        {
          min[j] = point[j] < min[j] ? point[j] : min[j];
          max[j] = point[j] > max[j] ? point[j] : max[j];
        }
      pointlist[i].index = i;
    }

  for (unsigned j = 0; j < 3; ++j)
    {
      min[j] = min[j] < 0 ? min[j] * 1.01 : min[j] * 0.99;
      max[j] = max[j] < 0 ? max[j] * 0.99 : max[j] * 1.01;
      gs->min[j] = min[j];
      gs->max[j] = max[j];
    }

  if (cdoVerbose)
    printf("BBOX: min=%g/%g/%g  max=%g/%g/%g\n", min[0], min[1], min[2], max[0],
           max[1], max[2]);

  kdTree_t *kdt
      = kd_buildTree(pointlist, n, min, max, 3, Threading::ompNumThreads);
  if (pointlist) Free(pointlist);
  if (kdt == NULL) cdoAbort("kd_buildTree failed!");

  return (void *) kdt;
}

static void *
gs_create_nanoflann(size_t n, const double *restrict lons,
                    const double *restrict lats, struct gridsearch *gs)
{
  PointCloud<double> *pointcloud = new PointCloud<double>();
  if (cdoVerbose)
    printf("nanoflann init 3D: n=%zu  nthreads=%d\n", n,
           Threading::ompNumThreads);

  double min[3] = { 1.e9, 1.e9, 1.e9 };
  double max[3] = { -1.e9, -1.e9, -1.e9 };

  // Generating  Point Cloud
  pointcloud->pts.resize(n);
#ifdef HAVE_OPENMP45
#pragma omp parallel for reduction(min : min[:3]) reduction(max : max[:3])
#endif
  for (size_t i = 0; i < n; i++)
    {
      double point[3];
      LLtoXYZ(lons[i], lats[i], point);
      pointcloud->pts[i].x = point[0];
      pointcloud->pts[i].y = point[1];
      pointcloud->pts[i].z = point[2];
      for (unsigned j = 0; j < 3; ++j)
        {
          min[j] = point[j] < min[j] ? point[j] : min[j];
          max[j] = point[j] > max[j] ? point[j] : max[j];
        }
    }

  gs->pointcloud = (void *) pointcloud;

  for (unsigned j = 0; j < 3; ++j)
    {
      min[j] = min[j] < 0 ? min[j] * 1.01 : min[j] * 0.99;
      max[j] = max[j] < 0 ? max[j] * 0.99 : max[j] * 1.01;
      gs->min[j] = min[j];
      gs->max[j] = max[j];
    }

  for (unsigned j = 0; j < 3; ++j)
    {
      pointcloud->min[j] = min[j];
      pointcloud->max[j] = max[j];
    }

  // construct a kd-tree index:
  nfTree_t *nft = new nfTree_t(
      3 /*dim*/, *pointcloud,
      nanoflann::KDTreeSingleIndexAdaptorParams(50 /* max leaf */));
  nft->buildIndex();

  return (void *) nft;
}

static void
gs_destroy_kdtree(void *search_container)
{
  kdTree_t *kdt = (kdTree_t *) search_container;
  if (kdt) kd_destroyTree(kdt);
}

static void
gs_destroy_full(void *search_container)
{
  struct gsFull *full = (struct gsFull *) search_container;
  if (full)
    {
      if (full->pts)
        {
          Free(full->pts[0]);
          Free(full->pts);
        }

      Free(full);
    }
}

static void *
gs_create_full(size_t n, const double *restrict lons,
               const double *restrict lats)
{
  struct gsFull *full = (struct gsFull *) Calloc(1, sizeof(struct gsFull));

  double **p = (double **) Malloc(n * sizeof(double *));
  p[0] = (double *) Malloc(3 * n * sizeof(double));
  for (size_t i = 1; i < n; i++)
    p[i] = p[0] + i * 3;

#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
  for (size_t i = 0; i < n; i++)
    {
      LLtoXYZ(lons[i], lats[i], p[i]);
    }

  full->n = n;
  full->plons = lons;
  full->plats = lats;
  full->pts = p;

  return (void *) full;
}

struct gridsearch *
gridsearch_create(bool xIsCyclic, size_t dims[2], size_t n,
                  const double *restrict lons, const double *restrict lats)
{
  struct gridsearch *gs
      = (struct gridsearch *) Calloc(1, sizeof(struct gridsearch));

  gs->is_cyclic = xIsCyclic;
  gs->is_curve = n != 1 && n == dims[0] * dims[1];
  gs->dims[0] = dims[0];
  gs->dims[1] = dims[1];

  gs->method_nn = gridsearch_method_nn;
  gs->n = n;
  if (n == 0) return gs;

  gs->plons = lons;
  gs->plats = lats;

  if (gs->method_nn == GridsearchMethod::kdtree)
    gs->search_container = gs_create_kdtree(n, lons, lats, gs);
  else if (gs->method_nn == GridsearchMethod::nanoflann)
    gs->search_container = gs_create_nanoflann(n, lons, lats, gs);
  else if (gs->method_nn == GridsearchMethod::full)
    gs->search_container = gs_create_full(n, lons, lats);

  gs->search_radius = cdo_default_search_radius();

  return gs;
}

void
gridsearch_delete(struct gridsearch *gs)
{
  if (gs)
    {
      if (gs->reg2d_center_lon) Free(gs->reg2d_center_lon);
      if (gs->reg2d_center_lat) Free(gs->reg2d_center_lat);

      if (gs->coslat) Free(gs->coslat);
      if (gs->coslon) Free(gs->coslon);
      if (gs->sinlat) Free(gs->sinlat);
      if (gs->sinlon) Free(gs->sinlon);

      if (gs->method_nn == GridsearchMethod::kdtree)
        gs_destroy_kdtree(gs->search_container);
      else if (gs->method_nn == GridsearchMethod::nanoflann)
        delete ((PointCloud<double> *) gs->pointcloud);
      else if (gs->method_nn == GridsearchMethod::full)
        gs_destroy_full(gs->search_container);

      Free(gs);
    }
}

static double
gs_set_range(double *prange)
{
  double range;

  if (prange)
    range = *prange;
  else
    range = SQR(2 * M_PI); /* This has to be bigger than the presumed
                            * maximum distance to the NN but smaller
                            * than once around the sphere. The content
                            * of this variable is replaced with the
                            * distance to the NN squared. */
  return range;
}

static size_t
gs_nearest_kdtree(void *search_container, double lon, double lat,
                  double *prange, struct gridsearch *gs)
{
  size_t index = GS_NOT_FOUND;
  kdTree_t *kdt = (kdTree_t *) search_container;
  if (kdt == NULL) return index;

  kdata_t range0 = gs_set_range(prange);
  kdata_t range = range0;

  kdata_t query_pt[3];
  LLtoXYZ(lon, lat, query_pt);
  /*
  if ( lon*RAD2DEG > -27 && lon*RAD2DEG < 60 &&  lat*RAD2DEG > 30 &&
  lat*RAD2DEG < 35 )
    {
      printf("lon %g lat %g\n", lon*RAD2DEG, lat*RAD2DEG);
      for ( unsigned j = 0; j < 3; ++j )
        if ( query_pt[j] < gs->min[j] || query_pt[j] > gs->max[j] )
          printf("  %g %g\n", query_pt[j] - gs->min[j],  query_pt[j] -
  gs->max[j]);
    }
  */
  if (!gs->extrapolate)
    for (unsigned j = 0; j < 3; ++j)
      if (query_pt[j] < gs->min[j] || query_pt[j] > gs->max[j]) return index;

  kdNode *node = kd_nearest(kdt->node, query_pt, &range, 3);

  kdata_t frange = range;
  if (!(frange < range0)) node = NULL;
  if (prange) *prange = frange;

  if (node) index = node->index;
  // printf("%zu %g\n", index, range);

  return index;
}

bool point_in_quad(bool is_cyclic, size_t nx, size_t ny, size_t i, size_t j,
                   size_t adds[4], double lons[4], double lats[4], double plon,
                   double plat, const double *restrict center_lon,
                   const double *restrict center_lat);

static size_t
llindex_in_quad(struct gridsearch *gs, size_t index, double lon, double lat)
{
  size_t ret_index = index;
  if (ret_index != GS_NOT_FOUND)
    {
      ret_index = GS_NOT_FOUND;
      size_t nx = gs->dims[0];
      size_t ny = gs->dims[1];
      size_t adds[4];
      double lons[4];
      double lats[4];
      bool is_cyclic = gs->is_cyclic;
      for (unsigned k = 0; k < 4; ++k)
        {
          /* Determine neighbor addresses */
          size_t j = index / nx;
          size_t i = index - j * nx;
          if (k == 1 || k == 3) i = (i > 0) ? i - 1 : (is_cyclic) ? nx - 1 : 0;
          if (k == 2 || k == 3) j = (j > 0) ? j - 1 : 0;

          if (point_in_quad(is_cyclic, nx, ny, i, j, adds, lons, lats, lon, lat,
                            gs->plons, gs->plats))
            {
              ret_index = index;
              break;
            }
        }
    }

  return ret_index;
}

static size_t
gs_nearest_nanoflann(void *search_container, double lon, double lat,
                     double *prange, struct gridsearch *gs)
{
  size_t index = GS_NOT_FOUND;
  nfTree_t *nft = (nfTree_t *) search_container;
  if (nft == NULL) return index;

  double range0 = gs_set_range(prange);
  double range = range0;

  double query_pt[3];
  LLtoXYZ(lon, lat, query_pt);

  if (!gs->extrapolate)
    for (unsigned j = 0; j < 3; ++j)
      if (query_pt[j] < gs->min[j] || query_pt[j] > gs->max[j]) return index;

  const size_t num_results = 1;
  size_t ret_index;
  double out_dist_sqr;
  nanoflann::KNNResultSet<double> resultSet(range, num_results);
  resultSet.init(&ret_index, &out_dist_sqr);
  nft->findNeighbors(resultSet, query_pt, nanoflann::SearchParams(10));
  // printf("%zu %g\n", ret_index, out_dist_sqr);

  index = ret_index;
  *prange = out_dist_sqr;

  // float frange = range;
  // if ( !(frange < range0) ) node = NULL;
  // if ( prange ) *prange = frange;

  // if ( node ) index = node->index;

  return index;
}

static size_t
gs_nearest_full(void *search_container, double lon, double lat, double *prange)
{
  size_t index = GS_NOT_FOUND;
  struct gsFull *full = (struct gsFull *) search_container;
  if (full == NULL) return index;

  double range0 = gs_set_range(prange);

  double query_pt[3];
  LLtoXYZ(lon, lat, query_pt);

  size_t n = full->n;
  size_t closestpt = n;
  double **pts = full->pts;
  double dist = FLT_MAX;
  for (size_t i = 0; i < n; i++)
    {
      double d = (float) distance(query_pt, pts[i]);
      if (closestpt >= n || d < dist || (d <= dist && i < closestpt))
        {
          dist = d;
          closestpt = i;
        }
    }

  if (closestpt < n)
    {
      if (dist < range0)
        {
          *prange = dist;
          index = closestpt;
        }
    }

  return index;
}

size_t
gridsearch_nearest(struct gridsearch *gs, double lon, double lat,
                   double *prange)
{
  size_t index = GS_NOT_FOUND;

  if (gs)
    {
      void *sc = gs->search_container;
      // clang-format off
      if      ( gs->method_nn == GridsearchMethod::kdtree )    index = gs_nearest_kdtree(sc, lon, lat, prange, gs);
      else if ( gs->method_nn == GridsearchMethod::nanoflann ) index = gs_nearest_nanoflann(sc, lon, lat, prange, gs);
      else if ( gs->method_nn == GridsearchMethod::full )      index = gs_nearest_full(sc, lon, lat, prange);
      else cdoAbort("%s::method_nn undefined!", __func__);
      // clang-format on

      if (!gs->extrapolate && gs->is_curve)
        index = llindex_in_quad(gs, index, lon, lat);
    }

  return index;
}

static size_t
gs_qnearest_kdtree(struct gridsearch *gs, double lon, double lat,
                   double *prange, size_t nnn, size_t *adds, double *dist)
{
  size_t nadds = 0;

  kdTree_t *kdt = (kdTree_t *) gs->search_container;
  if (kdt == NULL) return nadds;

  kdata_t range0 = gs_set_range(prange);
  kdata_t range = range0;
  struct pqueue *result = NULL;

  kdata_t query_pt[3];
  LLtoXYZ(lon, lat, query_pt);

  if (!gs->extrapolate)
    for (unsigned j = 0; j < 3; ++j)
      if (query_pt[j] < gs->min[j] || query_pt[j] > gs->max[j]) return nadds;

  if (gs)
    {
      result = kd_qnearest(kdt->node, query_pt, &range, nnn, 3);
      // printf("range %g %g %g %p\n", lon, lat, range, node);

      kdata_t frange = range;

      if (result)
        {
          size_t index;
          struct resItem *p;
          while (pqremove_min(result, &p))
            {
              index = p->node->index;
              range = p->dist_sq;
              Free(p);  // Free the result node taken from the heap

              if (range < range0)
                {
                  dist[nadds] = range;
                  adds[nadds] = index;
                  nadds++;
                }
            }
          Free(result->d);  // free the heap
          Free(result);     // and free the heap information structure
        }

      if (prange) *prange = frange;
    }

  return nadds;
}

static size_t
gs_qnearest_nanoflann(struct gridsearch *gs, double lon, double lat,
                      double *prange, size_t nnn, size_t *adds, double *dist)
{
  size_t nadds = 0;

  nfTree_t *nft = (nfTree_t *) gs->search_container;
  if (nft == NULL) return nadds;

  double range0 = gs_set_range(prange);
  double range = range0;

  double query_pt[3];
  LLtoXYZ(lon, lat, query_pt);

  if (!gs->extrapolate)
    for (unsigned j = 0; j < 3; ++j)
      if (query_pt[j] < gs->min[j] || query_pt[j] > gs->max[j]) return nadds;

  if (gs)
    {
      std::vector<double> out_dist_sqr(nnn);
      nadds = nft->knnRangeSearch(&query_pt[0], range, nnn, &adds[0],
                                  &out_dist_sqr[0]);

      for (size_t i = 0; i < nadds; ++i)
        dist[i] = out_dist_sqr[i];

      double frange = range;
      if (prange) *prange = frange;
    }

  return nadds;
}

size_t
gridsearch_qnearest(struct gridsearch *gs, double lon, double lat,
                    double *prange, size_t nnn, size_t *adds, double *dist)
{
  size_t nadds = 0;

  if (gs)
    {
      // clang-format off
      if      ( gs->method_nn == GridsearchMethod::kdtree )    nadds = gs_qnearest_kdtree(gs, lon, lat, prange, nnn, adds, dist);
      else if ( gs->method_nn == GridsearchMethod::nanoflann ) nadds = gs_qnearest_nanoflann(gs, lon, lat, prange, nnn, adds, dist);
      else cdoAbort("%s::method_nn undefined!", __func__);
      // clang-format on

      if (!gs->extrapolate && gs->is_curve)
        {
          size_t nadds_old = nadds;
          nadds = 0;
          for (size_t i = 0; i < nadds_old; ++i)
            {
              size_t index = adds[i];
              index = llindex_in_quad(gs, index, lon, lat);
              if (index != GS_NOT_FOUND)
                {
                  dist[nadds] = dist[i];
                  adds[nadds] = adds[i];
                  nadds++;
                }
            }
        }
    }

  return nadds;
}

#define BIGNUM 1.e+20
#define TINY 1.e-14

static void
knn_store_distance(size_t nadd, double distance, size_t numNeighbors,
                   size_t *restrict nbr_add, double *restrict nbr_dist)
{
  if (numNeighbors == 1)
    {
      if (distance < nbr_dist[0]
          || (distance <= nbr_dist[0] && nadd < nbr_add[0]))
        {
          nbr_add[0] = nadd;
          nbr_dist[0] = distance;
        }
    }
  else
    {
      for (size_t nchk = 0; nchk < numNeighbors; ++nchk)
        {
          if (distance < nbr_dist[nchk]
              || (distance <= nbr_dist[nchk] && nadd < nbr_add[nchk]))
            {
              for (size_t n = numNeighbors - 1; n > nchk; --n)
                {
                  nbr_add[n] = nbr_add[n - 1];
                  nbr_dist[n] = nbr_dist[n - 1];
                }
              nbr_add[nchk] = nadd;
              nbr_dist[nchk] = distance;
              break;
            }
        }
    }
}

static void
knn_check_distance(size_t numNeighbors, const size_t *restrict nbr_add,
                   double *restrict nbr_dist)
{
  // If distance is zero, set to small number
  for (size_t nchk = 0; nchk < numNeighbors; ++nchk)
    if (nbr_add[nchk] != GS_NOT_FOUND && nbr_dist[nchk] <= 0.)
      nbr_dist[nchk] = TINY;
}

void
gridsearch_knn_init(struct gsknn *knn)
{
  size_t ndist = knn->ndist;
  size_t *restrict add = knn->add;
  double *restrict dist = knn->dist;

  for (size_t i = 0; i < ndist; ++i)
    {
      add[i] = GS_NOT_FOUND;
      dist[i] = BIGNUM;
    }
}

struct gsknn *
gridsearch_knn_new(size_t size)
{
  struct gsknn *knn = (struct gsknn *) Malloc(sizeof(struct gsknn));

  knn->ndist = size;
  knn->size = size;
  knn->mask
      = (bool *) Malloc(size * sizeof(bool));  // mask at nearest neighbors
  knn->add = (size_t *) Malloc(
      size * sizeof(size_t));  // source address at nearest neighbors
  knn->dist = (double *) Malloc(
      size * sizeof(double));  // angular distance of the nearest neighbors
  knn->tmpadd = NULL;
  knn->tmpdist = NULL;

  gridsearch_knn_init(knn);

  return knn;
}

void
gridsearch_knn_delete(struct gsknn *knn)
{
  if (knn)
    {
      knn->size = 0;
      if (knn->dist) Free(knn->dist);
      if (knn->add) Free(knn->add);
      if (knn->tmpdist) Free(knn->tmpdist);
      if (knn->tmpadd) Free(knn->tmpadd);
      if (knn->mask) Free(knn->mask);
      Free(knn);
    }
}

size_t
gridsearch_knn(struct gridsearch *gs, struct gsknn *knn, double plon,
               double plat)
{
  /*
    Output variables:

    int nbr_add[numNeighbors]     ! address of each of the closest points
    double nbr_dist[numNeighbors] ! distance to each of the closest points

    Input variables:

    double plat,         ! latitude  of the search point
    double plon,         ! longitude of the search point
  */

  double search_radius = gs->search_radius;

  // Initialize distance and address arrays
  gridsearch_knn_init(knn);

  size_t numNeighbors = knn->size;
  size_t *restrict nbr_add = knn->add;
  double *restrict nbr_dist = knn->dist;

  size_t ndist = numNeighbors;
  // check some more points if distance is the same use the smaller index (nadd)
  if (ndist > 8)
    ndist += 8;
  else
    ndist *= 2;
  if (ndist > gs->n) ndist = gs->n;

  if (knn->tmpadd == NULL)
    knn->tmpadd = (size_t *) Malloc(ndist * sizeof(size_t));
  if (knn->tmpdist == NULL)
    knn->tmpdist = (double *) Malloc(ndist * sizeof(double));

  size_t *adds = knn->tmpadd;
  double *dist = knn->tmpdist;

  const double range0 = SQR(search_radius);
  double range = range0;

  size_t nadds = 0;

  if (numNeighbors == 1)
    {
      size_t add = gridsearch_nearest(gs, plon, plat, &range);
      if (add != GS_NOT_FOUND)
        {
          // if ( range < range0 )
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
      for (size_t i = 0; i < nadds; ++i)
        dist[i] = sqrt(dist[i]);
    }

  ndist = nadds;
  size_t max_neighbors = (ndist < numNeighbors) ? ndist : numNeighbors;

  for (size_t i = 0; i < ndist; ++i)
    knn_store_distance(adds[i], dist[i], max_neighbors, nbr_add, nbr_dist);

  knn_check_distance(max_neighbors, nbr_add, nbr_dist);

  if (ndist > numNeighbors) ndist = numNeighbors;

  knn->ndist = ndist;

  return ndist;
}  // gridsearch_knn

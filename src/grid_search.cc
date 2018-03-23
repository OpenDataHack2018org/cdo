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
#include "cdoOptions.h"
#include "dmemory.h"
#include "util.h"
#include "grid.h"
#include "grid_search.h"
#include "kdtreelib/kdtree.h"
#include "nanoflann.hpp"
extern "C" {
#include "lib/yac/sphere_part.h"
}

#define PI M_PI
#define PI2 (2.0 * PI)

PointSearchMethod pointSearchMethod(PointSearchMethod::nanoflann);

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

typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>,
                                            PointCloud<double>, 3 /* dim */
                                            >
    nfTree_t;

static double
cdoDefaultSearchRadius(void)
{
  extern double pointSearchRadius;

  double searchRadius = pointSearchRadius;

  if (searchRadius < 0.) searchRadius = 0.;
  if (searchRadius > 180.) searchRadius = 180.;

  searchRadius *= DEG2RAD;

  return searchRadius;
}

void
setPointSearchMethod(const char *methodstr)
{
  // clang-format off
  if      (STR_IS_EQ(methodstr, "kdtree"))     pointSearchMethod = PointSearchMethod::kdtree;
  else if (STR_IS_EQ(methodstr, "nanoflann"))  pointSearchMethod = PointSearchMethod::nanoflann;
  else if (STR_IS_EQ(methodstr, "spherepart")) pointSearchMethod = PointSearchMethod::spherepart;
  else if (STR_IS_EQ(methodstr, "full"))       pointSearchMethod = PointSearchMethod::full;
  else if (STR_IS_EQ(methodstr, "latbins"))    pointSearchMethod = PointSearchMethod::latbins;
  else cdoAbort("gridsearch method %s not available!", methodstr);
  // clang-format on
}

void
gridsearch_extrapolate(GridSearch *gs)
{
  gs->extrapolate = true;
}

GridSearch *
gridsearch_create_reg2d(bool xIsCyclic, size_t dims[2], const double *restrict lons, const double *restrict lats)
{
  GridSearch *gs = (GridSearch *) Calloc(1, sizeof(GridSearch));

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

  gs->searchRadius = cdoDefaultSearchRadius();

  return gs;
}

static void *
gs_create_kdtree(size_t n, const double *restrict lons, const double *restrict lats, GridSearch *gs)
{
  struct kd_point *pointlist = (struct kd_point *) Malloc(n * sizeof(struct kd_point));
  // see  example_cartesian.c
  if (cdoVerbose) cdoPrint("Init kdtree lib 3D: n=%zu  nthreads=%d", n, Threading::ompNumThreads);

  kdata_t min[3] = { 1.e9, 1.e9, 1.e9 };
  kdata_t max[3] = { -1.e9, -1.e9, -1.e9 };

#ifdef HAVE_OPENMP45
#pragma omp parallel for reduction(min : min[ : 3]) reduction(max : max[ : 3])
#endif
  for (size_t i = 0; i < n; i++)
    {
      kdata_t *restrict point = pointlist[i].point;
      cdoLLtoXYZ(lons[i], lats[i], point);
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

  if (cdoVerbose) cdoPrint("BBOX: min=%g/%g/%g  max=%g/%g/%g", min[0], min[1], min[2], max[0], max[1], max[2]);

  kdTree_t *kdt = kd_buildTree(pointlist, n, min, max, 3, Threading::ompNumThreads);
  if (pointlist) Free(pointlist);
  if (kdt == NULL) cdoAbort("kd_buildTree failed!");

  return (void *) kdt;
}

static void *
gs_create_nanoflann(size_t n, const double *restrict lons, const double *restrict lats, GridSearch *gs)
{
  PointCloud<double> *pointcloud = new PointCloud<double>();
  if (cdoVerbose) cdoPrint("Init nanoflann 3D: n=%zu  nthreads=%d", n, Threading::ompNumThreads);

  double min[3] = { 1.e9, 1.e9, 1.e9 };
  double max[3] = { -1.e9, -1.e9, -1.e9 };

  // Generating  Point Cloud
  pointcloud->pts.resize(n);
#ifdef HAVE_OPENMP45
#pragma omp parallel for reduction(min : min[ : 3]) reduction(max : max[ : 3])
#endif
  for (size_t i = 0; i < n; i++)
    {
      double point[3];
      cdoLLtoXYZ(lons[i], lats[i], point);
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
  nfTree_t *nft = new nfTree_t(3 /*dim*/, *pointcloud, nanoflann::KDTreeSingleIndexAdaptorParams(50 /* max leaf */));
  nft->buildIndex();

  return (void *) nft;
}

static void *
gs_create_spherepart(size_t n, const double *restrict lons, const double *restrict lats, GridSearch *gs)
{
  if (cdoVerbose) cdoPrint("Init spherepart 3D: n=%zu  nthreads=%d", n, Threading::ompNumThreads);

  double *coordinates_xyz = (double*) malloc(3 * n * sizeof(*coordinates_xyz));
  gs->coordinates_xyz = coordinates_xyz;

  double min[3] = { 1.e9, 1.e9, 1.e9 };
  double max[3] = { -1.e9, -1.e9, -1.e9 };

#ifdef HAVE_OPENMP45
#pragma omp parallel for reduction(min : min[ : 3]) reduction(max : max[ : 3])
#endif
  for (size_t i = 0; i < n; i++)
    {
      double *restrict point = coordinates_xyz + i * 3;
      cdoLLtoXYZ(lons[i], lats[i], point);
      for (unsigned j = 0; j < 3; ++j)
        {
          min[j] = point[j] < min[j] ? point[j] : min[j];
          max[j] = point[j] > max[j] ? point[j] : max[j];
        }
    }

  for (unsigned j = 0; j < 3; ++j)
    {
      min[j] = min[j] < 0 ? min[j] * 1.01 : min[j] * 0.99;
      max[j] = max[j] < 0 ? max[j] * 0.99 : max[j] * 1.01;
      gs->min[j] = min[j];
      gs->max[j] = max[j];
    }

  if (cdoVerbose) cdoPrint("BBOX: min=%g/%g/%g  max=%g/%g/%g", min[0], min[1], min[2], max[0], max[1], max[2]);

  return (void *) yac_point_sphere_part_search_new(n, coordinates_xyz);;
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

static void
gs_destroy_spherepart(void *search_container)
{
  yac_delete_point_sphere_part_search((struct point_sphere_part_search *)search_container);
}

static void *
gs_create_full(size_t n, const double *restrict lons, const double *restrict lats)
{
  if (cdoVerbose) cdoPrint("Init full grid search: n=%zu", n);

  struct gsFull *full = (struct gsFull *) Calloc(1, sizeof(struct gsFull));

  double **p = (double **) Malloc(n * sizeof(double *));
  p[0] = (double *) Malloc(3 * n * sizeof(double));
  for (size_t i = 1; i < n; i++) p[i] = p[0] + i * 3;

#ifdef HAVE_OPENMP4
#pragma omp simd
#endif
  for (size_t i = 0; i < n; i++)
    {
      cdoLLtoXYZ(lons[i], lats[i], p[i]);
    }

  full->n = n;
  full->plons = lons;
  full->plats = lats;
  full->pts = p;

  return (void *) full;
}

GridSearch *
gridsearch_create(bool xIsCyclic, size_t dims[2], size_t n, const double *restrict lons, const double *restrict lats)
{
  GridSearch *gs = (GridSearch *) Calloc(1, sizeof(GridSearch));

  gs->is_cyclic = xIsCyclic;
  gs->is_curve = n != 1 && n == dims[0] * dims[1];
  gs->dims[0] = dims[0];
  gs->dims[1] = dims[1];

  gs->method = pointSearchMethod;
  if ( gs->method == PointSearchMethod::latbins ) gs->method = PointSearchMethod::nanoflann;

  gs->n = n;
  if (n == 0) return gs;

  gs->plons = lons;
  gs->plats = lats;

  // clang-format off
  if      (gs->method == PointSearchMethod::kdtree)     gs->search_container = gs_create_kdtree(n, lons, lats, gs);
  else if (gs->method == PointSearchMethod::full)       gs->search_container = gs_create_full(n, lons, lats);
  else if (gs->method == PointSearchMethod::nanoflann)  gs->search_container = gs_create_nanoflann(n, lons, lats, gs);
  else if (gs->method == PointSearchMethod::spherepart) gs->search_container = gs_create_spherepart(n, lons, lats, gs);
  else cdoAbort("%s::method undefined!", __func__);
  // clang-format on
  
  gs->searchRadius = cdoDefaultSearchRadius();

  return gs;
}

void
gridsearch_delete(GridSearch *gs)
{
  if (gs)
    {
      if (gs->reg2d_center_lon) Free(gs->reg2d_center_lon);
      if (gs->reg2d_center_lat) Free(gs->reg2d_center_lat);

      if (gs->coslat) Free(gs->coslat);
      if (gs->coslon) Free(gs->coslon);
      if (gs->sinlat) Free(gs->sinlat);
      if (gs->sinlon) Free(gs->sinlon);

      // clang-format off
      if      (gs->method == PointSearchMethod::kdtree)     gs_destroy_kdtree(gs->search_container);
      else if (gs->method == PointSearchMethod::nanoflann)  delete ((PointCloud<double> *) gs->pointcloud);
      else if (gs->method == PointSearchMethod::spherepart)
        {
          free(gs->coordinates_xyz);
          gs_destroy_spherepart(gs->search_container);
        }
      else if (gs->method == PointSearchMethod::full)       gs_destroy_full(gs->search_container);
      // clang-format on

      Free(gs);
    }
}

static size_t
gs_nearest_kdtree(void *search_container, double lon, double lat, double searchRadius, size_t *addr, double *dist, GridSearch *gs)
{
  kdTree_t *kdt = (kdTree_t *) search_container;
  if (kdt == NULL) return 0;

  double sqrDistMax = SQR(searchRadius);
  double sqrDist = sqrDistMax;

  double query_pt[3];
  cdoLLtoXYZ(lon, lat, query_pt);

  if (!gs->extrapolate)
    for (unsigned j = 0; j < 3; ++j)
      if (query_pt[j] < gs->min[j] || query_pt[j] > gs->max[j]) return 0;

  kdNode *node = kd_nearest(kdt->node, query_pt, &sqrDist, 3);


  if (node && sqrDist < sqrDistMax)
    {
      *addr = node->index;
      *dist = sqrt(sqrDist);
      return 1;
    }

  return 0;
}

static size_t
gs_nearest_nanoflann(void *search_container, double lon, double lat, double searchRadius, size_t *addr, double *dist, GridSearch *gs)
{
  nfTree_t *nft = (nfTree_t *) search_container;
  if (nft == NULL) return 0;

  double sqrDistMax = SQR(searchRadius);

  double query_pt[3];
  cdoLLtoXYZ(lon, lat, query_pt);

  if (!gs->extrapolate)
    for (unsigned j = 0; j < 3; ++j)
      if (query_pt[j] < gs->min[j] || query_pt[j] > gs->max[j]) return 0;

  const size_t num_results = 1;
  size_t retIndex;
  double sqrDist;
  nanoflann::KNNResultSet<double> resultSet(sqrDistMax, num_results);
  resultSet.init(&retIndex, &sqrDist);
  nft->findNeighbors(resultSet, query_pt, nanoflann::SearchParams(10));

  if (retIndex != GS_NOT_FOUND)
    {
      *addr = retIndex;
      *dist = sqrt(sqrDist);
      return 1;
    }

  return 0;
}

static size_t
gs_nearest_spherepart(void *search_container, double lon, double lat, double searchRadius, size_t *addr, double *dist, GridSearch *gs)
{
  double query_pt[3];
  cdoLLtoXYZ(lon, lat, query_pt);

  if (!gs->extrapolate)
    for (unsigned j = 0; j < 3; ++j)
      if (query_pt[j] < gs->min[j] || query_pt[j] > gs->max[j]) return 0;

  size_t local_point_ids_array_size = 0;
  size_t num_local_point_ids;
  unsigned *local_point_ids = NULL;
  double cos_angle;

  yac_point_sphere_part_search_NN((struct point_sphere_part_search *)search_container, 1, query_pt, &cos_angle, NULL, NULL,
                                  &local_point_ids, &local_point_ids_array_size, &num_local_point_ids);

  size_t nadd = 0;
  if ( num_local_point_ids > 0 )
    {
      if ( cos_angle < -1 ) cos_angle = -1;
      if ( cos_angle >  1 ) cos_angle =  1;
      *dist = acos(cos_angle);
      if ( *dist <= searchRadius )
        {
          *addr = local_point_ids[0];
          nadd = 1;
        }
    }

  if (local_point_ids) free(local_point_ids);

  return nadd;
}

static size_t
gs_nearest_full(void *search_container, double lon, double lat, double searchRadius, size_t *addr, double *dist)
{
  struct gsFull *full = (struct gsFull *) search_container;
  if (full == NULL) return 0;

  double sqrDistMax = SQR(searchRadius);

  double query_pt[3];
  cdoLLtoXYZ(lon, lat, query_pt);

  size_t n = full->n;
  size_t closestpt = n;
  double **pts = full->pts;
  double sqrDist = FLT_MAX;
  for (size_t i = 0; i < n; i++)
    {
      double d = (float) squareDistance(query_pt, pts[i]);
      if (closestpt >= n || d < sqrDist || (d <= sqrDist && i < closestpt))
        {
          sqrDist = d;
          closestpt = i;
        }
    }

  if (closestpt < n && sqrDist < sqrDistMax)
    {
      *addr = closestpt;
      *dist = sqrt(sqrDist);
      return 1;
    }

  return 0;
}

bool point_in_quad(bool is_cyclic, size_t nx, size_t ny, size_t i, size_t j, size_t adds[4], double lons[4],
                   double lats[4], double plon, double plat, const double *restrict center_lon,
                   const double *restrict center_lat);

static size_t
llindex_in_quad(GridSearch *gs, size_t index, double lon, double lat)
{
  if (index != GS_NOT_FOUND)
    {
      size_t nx = gs->dims[0];
      size_t ny = gs->dims[1];
      size_t adds[4];
      double lons[4];
      double lats[4];
      bool is_cyclic = gs->is_cyclic;
      for (unsigned k = 0; k < 4; ++k)
        {
          // Determine neighbor addresses
          size_t j = index / nx;
          size_t i = index - j * nx;
          if (k == 1 || k == 3) i = (i > 0) ? i - 1 : (is_cyclic) ? nx - 1 : 0;
          if (k == 2 || k == 3) j = (j > 0) ? j - 1 : 0;

          if (point_in_quad(is_cyclic, nx, ny, i, j, adds, lons, lats, lon, lat, gs->plons, gs->plats))
            return index;
        }
    }

  return GS_NOT_FOUND;
}

size_t
gridsearch_nearest(GridSearch *gs, double lon, double lat, size_t *addr, double *dist)
{
  if (gs)
    {
      size_t nadds = 0;
      double searchRadius = gs->searchRadius;
      void *sc = gs->search_container;
      // clang-format off
      if      ( gs->method == PointSearchMethod::kdtree )     nadds = gs_nearest_kdtree(sc, lon, lat, searchRadius, addr, dist, gs);
      else if ( gs->method == PointSearchMethod::nanoflann )  nadds = gs_nearest_nanoflann(sc, lon, lat, searchRadius, addr, dist, gs);
      else if ( gs->method == PointSearchMethod::spherepart ) nadds = gs_nearest_spherepart(sc, lon, lat, searchRadius, addr, dist, gs);
      else if ( gs->method == PointSearchMethod::full )       nadds = gs_nearest_full(sc, lon, lat, searchRadius, addr, dist);
      else cdoAbort("%s::method undefined!", __func__);
      // clang-format on

      if (nadds > 0)
        {
          size_t index = *addr;
          if (!gs->extrapolate && gs->is_curve) index = llindex_in_quad(gs, *addr, lon, lat);
          if (index != GS_NOT_FOUND) return 1;
        }
    }

  return 0;
}

static size_t
gs_qnearest_kdtree(GridSearch *gs, double lon, double lat, double searchRadius, size_t nnn, size_t *adds, double *dist)
{
  size_t nadds = 0;

  kdTree_t *kdt = (kdTree_t *) gs->search_container;
  if (kdt == NULL) return nadds;

  double sqrDistMax = SQR(searchRadius);
  kdata_t sqrDist = sqrDistMax;
  struct pqueue *result = NULL;

  kdata_t query_pt[3];
  cdoLLtoXYZ(lon, lat, query_pt);

  if (!gs->extrapolate)
    for (unsigned j = 0; j < 3; ++j)
      if (query_pt[j] < gs->min[j] || query_pt[j] > gs->max[j]) return nadds;

  if (gs)
    {
      result = kd_qnearest(kdt->node, query_pt, &sqrDist, nnn, 3);
      if (result)
        {
          size_t index;
          struct resItem *p;
          while (pqremove_min(result, &p))
            {
              index = p->node->index;
              sqrDist = p->dist_sq;
              Free(p);  // Free the result node taken from the heap

              if (sqrDist < sqrDistMax)
                {
                  adds[nadds] = index;
                  dist[nadds] = sqrt(sqrDist);
                  nadds++;
                }
            }
          Free(result->d);  // free the heap
          Free(result);     // and free the heap information structure
        }
    }

  return nadds;
}

static size_t
gs_qnearest_nanoflann(GridSearch *gs, double lon, double lat, double searchRadius, size_t nnn, size_t *adds, double *dist)
{
  size_t nadds = 0;

  nfTree_t *nft = (nfTree_t *) gs->search_container;
  if (nft == NULL) return nadds;

  double sqrDistMax = SQR(searchRadius);

  double query_pt[3];
  cdoLLtoXYZ(lon, lat, query_pt);

  if (!gs->extrapolate)
    for (unsigned j = 0; j < 3; ++j)
      if (query_pt[j] < gs->min[j] || query_pt[j] > gs->max[j]) return nadds;

  nadds = nft->knnRangeSearch(&query_pt[0], sqrDistMax, nnn, &adds[0], &dist[0]);
  for ( size_t i = 0; i < nadds; ++i ) dist[i] = sqrt(dist[i]);

  return nadds;
}

static size_t
gs_qnearest_spherepart(GridSearch *gs, double lon, double lat, double searchRadius, size_t nnn, size_t *adds, double *dist)
{
  size_t nadds = 0;

  double query_pt[3];
  cdoLLtoXYZ(lon, lat, query_pt);

  if (!gs->extrapolate)
    for (unsigned j = 0; j < 3; ++j)
      if (query_pt[j] < gs->min[j] || query_pt[j] > gs->max[j]) return nadds;

  if (gs)
    {
      size_t local_point_ids_array_size = 0;
      size_t num_local_point_ids;
      unsigned *local_point_ids = NULL;

      size_t cos_angles_array_size = 0;
      double *cos_angles = NULL;

      yac_point_sphere_part_search_NNN((struct point_sphere_part_search *)gs->search_container, 1, query_pt, nnn, &cos_angles, &cos_angles_array_size,
                                       NULL, NULL, &local_point_ids, &local_point_ids_array_size, &num_local_point_ids);
      nadds = num_local_point_ids;
      if ( nadds )
        {
          size_t naddsmax = nadds;
          nadds = 0;
          for (size_t i = 0; i < naddsmax; ++i)
            {
              if ( cos_angles[i] < -1 ) cos_angles[i] = -1;
              if ( cos_angles[i] >  1 ) cos_angles[i] =  1;
              double angle = acos(cos_angles[i]);
              if ( angle < searchRadius )
                {
                  adds[nadds] = local_point_ids[i];
                  dist[nadds] = angle;
                  nadds++;
                }       
            }
        }

      free(cos_angles);
      free(local_point_ids);
    }

  return nadds;
}

size_t
gridsearch_qnearest(GridSearch *gs, double lon, double lat, size_t nnn, size_t *adds, double *dist)
{
  size_t nadds = 0;

  if (gs)
    {
      double searchRadius = gs->searchRadius;
      // clang-format off
      if      ( gs->method == PointSearchMethod::kdtree )     nadds = gs_qnearest_kdtree(gs, lon, lat, searchRadius, nnn, adds, dist);
      else if ( gs->method == PointSearchMethod::nanoflann )  nadds = gs_qnearest_nanoflann(gs, lon, lat, searchRadius, nnn, adds, dist);
      else if ( gs->method == PointSearchMethod::spherepart ) nadds = gs_qnearest_spherepart(gs, lon, lat, searchRadius, nnn, adds, dist);
      else cdoAbort("%s::method undefined!", __func__);
      // clang-format on

      if (!gs->extrapolate && gs->is_curve)
        {
          size_t naddsmax = nadds;
          nadds = 0;
          for (size_t i = 0; i < naddsmax; ++i)
            {
              size_t index = llindex_in_quad(gs, adds[i], lon, lat);
              if (index != GS_NOT_FOUND)
                {
                  adds[nadds] = adds[i];
                  dist[nadds] = dist[i];
                  nadds++;
                }
            }
        }
    }

  return nadds;
}

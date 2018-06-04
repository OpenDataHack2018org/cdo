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
#include "grid_point_search.h"
#include "kdtreelib/kdtree.h"
#include "nanoflann.hpp"
extern "C" {
#include "lib/yac/sphere_part.h"
}

#define PI M_PI
#define PI2 (2.0 * PI)

PointSearchMethod pointSearchMethod(PointSearchMethod::nanoflann);

struct gpsFull
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

typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>, PointCloud<double>, 3>
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
  else cdoAbort("Grid point search method %s not available!", methodstr);
  // clang-format on
}

void
gridSearchExtrapolate(GridPointSearch *gps)
{
  gps->extrapolate = true;
}

GridPointSearch *
gridPointSearchCreateReg2d(bool xIsCyclic, size_t dims[2], const double *restrict lons, const double *restrict lats)
{
  GridPointSearch *gps = (GridPointSearch *) Calloc(1, sizeof(GridPointSearch));

  gps->is_cyclic = xIsCyclic;
  gps->is_reg2d = true;
  gps->dims[0] = dims[0];
  gps->dims[1] = dims[1];
  size_t nx = dims[0];
  size_t ny = dims[1];

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

  gps->reg2d_center_lon = reg2d_center_lon;
  gps->reg2d_center_lat = reg2d_center_lat;

  gps->coslon = coslon;
  gps->sinlon = sinlon;
  gps->coslat = coslat;
  gps->sinlat = sinlat;

  gps->searchRadius = cdoDefaultSearchRadius();

  return gps;
}

static void *
gps_create_kdtree(size_t n, const double *restrict lons, const double *restrict lats, GridPointSearch *gps)
{
  struct kd_point *pointlist = (struct kd_point *) Malloc(n * sizeof(struct kd_point));
  // see  example_cartesian.c

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
      gps->min[j] = min[j];
      gps->max[j] = max[j];
    }

  if (cdoVerbose) cdoPrint("BBOX: min=%g/%g/%g  max=%g/%g/%g", min[0], min[1], min[2], max[0], max[1], max[2]);

  kdTree_t *kdt = kd_buildTree(pointlist, n, min, max, 3, Threading::ompNumThreads);
  if (pointlist) Free(pointlist);
  if (kdt == NULL) cdoAbort("kd_buildTree failed!");

  return (void *) kdt;
}

static void *
gps_create_nanoflann(size_t n, const double *restrict lons, const double *restrict lats, GridPointSearch *gps)
{
  PointCloud<double> *pointcloud = new PointCloud<double>();

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

  gps->pointcloud = (void *) pointcloud;

  for (unsigned j = 0; j < 3; ++j)
    {
      min[j] = min[j] < 0 ? min[j] * 1.01 : min[j] * 0.99;
      max[j] = max[j] < 0 ? max[j] * 0.99 : max[j] * 1.01;
      gps->min[j] = min[j];
      gps->max[j] = max[j];
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
gps_create_spherepart(size_t n, const double *restrict lons, const double *restrict lats, GridPointSearch *gps)
{
  double (*coordinates_xyz)[3];
  coordinates_xyz = (double (*)[3]) malloc(n * sizeof(*coordinates_xyz));
  gps->coordinates_xyz = coordinates_xyz;

  double min[3] = { 1.e9, 1.e9, 1.e9 };
  double max[3] = { -1.e9, -1.e9, -1.e9 };

#ifdef HAVE_OPENMP45
#pragma omp parallel for reduction(min : min[ : 3]) reduction(max : max[ : 3])
#endif
  for (size_t i = 0; i < n; i++)
    {
      double *restrict point = coordinates_xyz[i];
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
      gps->min[j] = min[j];
      gps->max[j] = max[j];
    }

  if (cdoVerbose) cdoPrint("BBOX: min=%g/%g/%g  max=%g/%g/%g", min[0], min[1], min[2], max[0], max[1], max[2]);

  return (void *) yac_point_sphere_part_search_new(n, gps->coordinates_xyz);
}

static void
gps_destroy_kdtree(void *search_container)
{
  kdTree_t *kdt = (kdTree_t *) search_container;
  if (kdt) kd_destroyTree(kdt);
}

static void
gps_destroy_full(void *search_container)
{
  struct gpsFull *full = (struct gpsFull *) search_container;
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
gps_destroy_spherepart(void *search_container)
{
  yac_delete_point_sphere_part_search((struct point_sphere_part_search *) search_container);
}

static void *
gps_create_full(size_t n, const double *restrict lons, const double *restrict lats)
{
  if (cdoVerbose) cdoPrint("Init full grid search: n=%zu", n);

  struct gpsFull *full = (struct gpsFull *) Calloc(1, sizeof(struct gpsFull));

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

GridPointSearch *
gridPointSearchCreate(bool xIsCyclic, size_t dims[2], size_t n, const double *restrict lons, const double *restrict lats)
{
  GridPointSearch *gps = (GridPointSearch *) Calloc(1, sizeof(GridPointSearch));

  gps->is_cyclic = xIsCyclic;
  gps->is_curve = n != 1 && n == dims[0] * dims[1];
  gps->dims[0] = dims[0];
  gps->dims[1] = dims[1];

  gps->method = pointSearchMethod;
  if (gps->method == PointSearchMethod::latbins) gps->method = PointSearchMethod::nanoflann;

  gps->n = n;
  if (n == 0) return gps;

  gps->plons = lons;
  gps->plats = lats;

  // clang-format off
  if (cdoVerbose)
    {
      if      (gps->method == PointSearchMethod::kdtree)     cdoPrint("Point search method: kdtree");
      else if (gps->method == PointSearchMethod::full)       cdoPrint("Point search method: full");
      else if (gps->method == PointSearchMethod::nanoflann)  cdoPrint("Point search method: nanoflann");
      else if (gps->method == PointSearchMethod::spherepart) cdoPrint("Point search method: spherepart");
    }

  if      (gps->method == PointSearchMethod::kdtree)     gps->search_container = gps_create_kdtree(n, lons, lats, gps);
  else if (gps->method == PointSearchMethod::full)       gps->search_container = gps_create_full(n, lons, lats);
  else if (gps->method == PointSearchMethod::nanoflann)  gps->search_container = gps_create_nanoflann(n, lons, lats, gps);
  else if (gps->method == PointSearchMethod::spherepart) gps->search_container = gps_create_spherepart(n, lons, lats, gps);
  else cdoAbort("%s::method undefined!", __func__);
  // clang-format on

  gps->searchRadius = cdoDefaultSearchRadius();

  return gps;
}

void
gridSearchDelete(GridPointSearch *gps)
{
  if (gps)
    {
      if (gps->reg2d_center_lon) Free(gps->reg2d_center_lon);
      if (gps->reg2d_center_lat) Free(gps->reg2d_center_lat);

      if (gps->coslat) Free(gps->coslat);
      if (gps->coslon) Free(gps->coslon);
      if (gps->sinlat) Free(gps->sinlat);
      if (gps->sinlon) Free(gps->sinlon);

      // clang-format off
      if      (gps->method == PointSearchMethod::kdtree)     gps_destroy_kdtree(gps->search_container);
      else if (gps->method == PointSearchMethod::nanoflann)
        {
          delete ((PointCloud<double> *) gps->pointcloud);
          delete ((nfTree_t *) gps->search_container);
        }
      else if (gps->method == PointSearchMethod::spherepart)
        {
          free(gps->coordinates_xyz);
          gps_destroy_spherepart(gps->search_container);
        }
      else if (gps->method == PointSearchMethod::full)       gps_destroy_full(gps->search_container);
      // clang-format on

      Free(gps);
    }
}

static size_t
gps_nearest_kdtree(void *search_container, double lon, double lat, double searchRadius, size_t *addr, double *dist, GridPointSearch *gps)
{
  kdTree_t *kdt = (kdTree_t *) search_container;
  if (kdt == NULL) return 0;

  double sqrDistMax = SQR(searchRadius);
  double sqrDist = sqrDistMax;

  double query_pt[3];
  cdoLLtoXYZ(lon, lat, query_pt);

  if (!gps->extrapolate)
    for (unsigned j = 0; j < 3; ++j)
      if (query_pt[j] < gps->min[j] || query_pt[j] > gps->max[j]) return 0;

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
gps_nearest_nanoflann(void *search_container, double lon, double lat, double searchRadius, size_t *addr, double *dist,
                     GridPointSearch *gps)
{
  nfTree_t *nft = (nfTree_t *) search_container;
  if (nft == NULL) return 0;

  double sqrDistMax = SQR(searchRadius);

  double query_pt[3];
  cdoLLtoXYZ(lon, lat, query_pt);

  if (!gps->extrapolate)
    for (unsigned j = 0; j < 3; ++j)
      if (query_pt[j] < gps->min[j] || query_pt[j] > gps->max[j]) return 0;

  const size_t num_results = 1;
  size_t retIndex;
  double sqrDist;
  nanoflann::KNNResultSet<double> resultSet(sqrDistMax, num_results);
  resultSet.init(&retIndex, &sqrDist);
  nft->findNeighbors(resultSet, query_pt, nanoflann::SearchParams(10));

  if (retIndex != GPS_NOT_FOUND)
    {
      *addr = retIndex;
      *dist = sqrt(sqrDist);
      return 1;
    }

  return 0;
}

static size_t
gps_nearest_spherepart(void *search_container, double lon, double lat, double searchRadius, size_t *addr, double *dist,
                      GridPointSearch *gps)
{
  double query_pt[1][3];
  cdoLLtoXYZ(lon, lat, query_pt[0]);

  if (!gps->extrapolate)
    for (unsigned j = 0; j < 3; ++j)
      if (query_pt[0][j] < gps->min[j] || query_pt[0][j] > gps->max[j]) return 0;

  size_t local_point_ids_array_size = 0;
  size_t num_local_point_ids;
  unsigned *local_point_ids = NULL;
  double cos_angle;

  yac_point_sphere_part_search_NN((struct point_sphere_part_search *) search_container, 1, query_pt, &cos_angle, NULL, NULL,
                                  &local_point_ids, &local_point_ids_array_size, &num_local_point_ids);

  size_t nadd = 0;
  if (num_local_point_ids > 0)
    {
      *dist = acos(cos_angle);
      if (*dist <= searchRadius)
        {
          nadd = 1;
          *addr = local_point_ids[0];
          for (size_t i = 1; i < num_local_point_ids; ++i)
            if (local_point_ids[i] < *addr) *addr = local_point_ids[i];
        }
    }

  if (local_point_ids) free(local_point_ids);

  return nadd;
}

static size_t
gps_nearest_full(void *search_container, double lon, double lat, double searchRadius, size_t *addr, double *dist)
{
  struct gpsFull *full = (struct gpsFull *) search_container;
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

bool pointInQuad(bool is_cyclic, size_t nx, size_t ny, size_t i, size_t j, size_t adds[4], double lons[4], double lats[4],
                 double plon, double plat, const double *restrict centerLon, const double *restrict centerLat);

static size_t
llindex_in_quad(GridPointSearch *gps, size_t index, double lon, double lat)
{
  if (index != GPS_NOT_FOUND)
    {
      size_t nx = gps->dims[0];
      size_t ny = gps->dims[1];
      size_t adds[4];
      double lons[4];
      double lats[4];
      bool isCyclic = gps->is_cyclic;
      for (unsigned k = 0; k < 4; ++k)
        {
          // Determine neighbor addresses
          size_t j = index / nx;
          size_t i = index - j * nx;
          if (k == 1 || k == 3) i = (i > 0) ? i - 1 : (isCyclic) ? nx - 1 : 0;
          if (k == 2 || k == 3) j = (j > 0) ? j - 1 : 0;

          if (pointInQuad(isCyclic, nx, ny, i, j, adds, lons, lats, lon, lat, gps->plons, gps->plats)) return index;
        }
    }

  return GPS_NOT_FOUND;
}

size_t
gridSearchNearest(GridPointSearch *gps, double lon, double lat, size_t *addr, double *dist)
{
  if (gps)
    {
      size_t nadds = 0;
      double searchRadius = gps->searchRadius;
      void *sc = gps->search_container;
      // clang-format off
      if      ( gps->method == PointSearchMethod::kdtree )     nadds = gps_nearest_kdtree(sc, lon, lat, searchRadius, addr, dist, gps);
      else if ( gps->method == PointSearchMethod::nanoflann )  nadds = gps_nearest_nanoflann(sc, lon, lat, searchRadius, addr, dist, gps);
      else if ( gps->method == PointSearchMethod::spherepart ) nadds = gps_nearest_spherepart(sc, lon, lat, searchRadius, addr, dist, gps);
      else if ( gps->method == PointSearchMethod::full )       nadds = gps_nearest_full(sc, lon, lat, searchRadius, addr, dist);
      else cdoAbort("%s::method undefined!", __func__);
      // clang-format on

      if (nadds > 0)
        {
          size_t index = *addr;
          if (!gps->extrapolate && gps->is_curve) index = llindex_in_quad(gps, *addr, lon, lat);
          if (index != GPS_NOT_FOUND) return 1;
        }
    }

  return 0;
}

static size_t
gps_qnearest_kdtree(GridPointSearch *gps, double lon, double lat, double searchRadius, size_t nnn, size_t *adds, double *dist)
{
  size_t nadds = 0;

  kdTree_t *kdt = (kdTree_t *) gps->search_container;
  if (kdt == NULL) return nadds;

  double sqrDistMax = SQR(searchRadius);
  kdata_t sqrDist = sqrDistMax;
  struct pqueue *result = NULL;

  kdata_t query_pt[3];
  cdoLLtoXYZ(lon, lat, query_pt);

  if (!gps->extrapolate)
    for (unsigned j = 0; j < 3; ++j)
      if (query_pt[j] < gps->min[j] || query_pt[j] > gps->max[j]) return nadds;

  if (gps)
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
gps_qnearest_nanoflann(GridPointSearch *gps, double lon, double lat, double searchRadius, size_t nnn, size_t *adds, double *dist)
{
  size_t nadds = 0;

  nfTree_t *nft = (nfTree_t *) gps->search_container;
  if (nft == NULL) return nadds;

  double sqrDistMax = SQR(searchRadius);

  double query_pt[3];
  cdoLLtoXYZ(lon, lat, query_pt);

  if (!gps->extrapolate)
    for (unsigned j = 0; j < 3; ++j)
      if (query_pt[j] < gps->min[j] || query_pt[j] > gps->max[j]) return nadds;

  nadds = nft->knnRangeSearch(&query_pt[0], sqrDistMax, nnn, &adds[0], &dist[0]);
  for (size_t i = 0; i < nadds; ++i) dist[i] = sqrt(dist[i]);

  return nadds;
}

static size_t
gps_qnearest_spherepart(GridPointSearch *gps, double lon, double lat, double searchRadius, size_t nnn, size_t *adds, double *dist)
{
  size_t nadds = 0;

  if (gps)
    {
      double query_pt[1][3];
      cdoLLtoXYZ(lon, lat, query_pt[0]);

      if (!gps->extrapolate)
        for (unsigned j = 0; j < 3; ++j)
          if (query_pt[0][j] < gps->min[j] || query_pt[0][j] > gps->max[j]) return nadds;

      size_t local_point_ids_array_size = 0;
      size_t num_local_point_ids;
      unsigned *local_point_ids = NULL;

      size_t cos_angles_array_size = 0;
      double *cos_angles = NULL;

      yac_point_sphere_part_search_NNN((struct point_sphere_part_search *) gps->search_container, 1, query_pt, nnn, &cos_angles,
                                       &cos_angles_array_size, NULL, NULL, &local_point_ids, &local_point_ids_array_size,
                                       &num_local_point_ids);

      if (num_local_point_ids > 0)
        {
          size_t maxadds = (num_local_point_ids < nnn) ? num_local_point_ids : nnn;
          nadds = 0;
          for (size_t i = 0; i < maxadds; ++i)
            {
              double angle = acos(cos_angles[i]);
              if (angle < searchRadius)
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
gridSearchQnearest(GridPointSearch *gps, double lon, double lat, size_t nnn, size_t *adds, double *dist)
{
  size_t nadds = 0;

  if (gps)
    {
      double searchRadius = gps->searchRadius;
      // clang-format off
      if      ( gps->method == PointSearchMethod::kdtree )     nadds = gps_qnearest_kdtree(gps, lon, lat, searchRadius, nnn, adds, dist);
      else if ( gps->method == PointSearchMethod::nanoflann )  nadds = gps_qnearest_nanoflann(gps, lon, lat, searchRadius, nnn, adds, dist);
      else if ( gps->method == PointSearchMethod::spherepart ) nadds = gps_qnearest_spherepart(gps, lon, lat, searchRadius, nnn, adds, dist);
      else cdoAbort("%s::method undefined!", __func__);
      // clang-format on

      if (!gps->extrapolate && gps->is_curve)
        {
          size_t naddsmax = nadds;
          nadds = 0;
          for (size_t i = 0; i < naddsmax; ++i)
            {
              size_t index = llindex_in_quad(gps, adds[i], lon, lat);
              if (index != GPS_NOT_FOUND)
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

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
#ifndef GRID_SEARCH_H_
#define GRID_SEARCH_H_

#include <stdbool.h>
#include "nbr_weights.h"

#define GS_NOT_FOUND SIZE_MAX

enum struct GridsearchMethod
{
  full,
  nanoflann,
  kdtree
};

struct gridsearch
{
  bool extrapolate;
  bool is_cyclic;
  bool is_reg2d;
  bool is_curve;
  GridsearchMethod method_nn;
  size_t n;
  size_t dims[2];

  void *search_container;
  double search_radius;

  // reg2d search
  double *reg2d_center_lon, *reg2d_center_lat;
  double *coslat, *sinlat;  // cosine, sine of grid lats (for distance)
  double *coslon, *sinlon;  // cosine, sine of grid lons (for distance)

  const double *plons, *plats;

  double lonmin, lonmax, latmin, latmax;
  float min[3], max[3];
  void *pointcloud;
};

struct gsknn
{
  size_t ndist;
  size_t *tmpadd;
  double *tmpdist;

  //nbrWeightsType nbrWeights;
  size_t size;
  bool *mask;
  size_t *add;
  double *dist;
};

struct gsknn *gridsearch_knn_new(size_t size);
void gridsearch_knn_delete(struct gsknn *knn);
size_t gridsearch_knn(struct gridsearch *gs, struct gsknn *knn, double plon, double plat);

struct gridsearch *gridsearch_create_reg2d(bool xIsCyclic, size_t dims[2], const double *restrict lons,
                                           const double *restrict lats);
struct gridsearch *gridsearch_create(bool xIsCyclic, size_t dims[2], size_t n, const double *restrict lons,
                                     const double *restrict lats);
void gridsearch_delete(struct gridsearch *gs);
size_t gridsearch_nearest(struct gridsearch *gs, double lon, double lat, double *range);
size_t gridsearch_qnearest(struct gridsearch *gs, double lon, double lat, double *prange, size_t nnn, size_t *adds,
                           double *dist);
void gridsearch_extrapolate(struct gridsearch *gs);
void gridsearch_bound_poly(struct gridsearch *gs, size_t dims[2], size_t n, const double *restrict lons,
                           const double *restrict lats);

#endif

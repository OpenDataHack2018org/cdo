#ifndef  GRID_SEARCH_H_
#define  GRID_SEARCH_H_

#include <stdbool.h>

#define  GS_NOT_FOUND  SIZE_MAX

enum struct GridsearchMethod {full, nanoflann, kdtree};


struct gridsearch {
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
  double *coslat, *sinlat;   // cosine, sine of grid lats (for distance)
  double *coslon, *sinlon;   // cosine, sine of grid lons (for distance)

  const double *plons, *plats;

  double lonmin, lonmax, latmin, latmax;
  float min[3], max[3];
  void *pointcloud;
};

struct gsknn {
  size_t   ndist;
  size_t   size;
  bool    *mask;
  size_t  *add;
  size_t  *tmpadd;
  double  *dist;
  double  *tmpdist;
};

struct gsknn *gridsearch_knn_new(size_t size);
void gridsearch_knn_delete(struct gsknn *knn);
size_t gridsearch_knn(struct gridsearch *gs, struct gsknn *knn, double plon, double plat);

struct gridsearch *gridsearch_create_reg2d(bool xIsCyclic, size_t dims[2], const double *restrict lons, const double *restrict lats);
struct gridsearch *gridsearch_create(bool xIsCyclic, size_t dims[2], size_t n, const double *restrict lons, const double *restrict lats);
void gridsearch_delete(struct gridsearch *gs);
size_t gridsearch_nearest(struct gridsearch *gs, double lon, double lat, double *range);
size_t gridsearch_qnearest(struct gridsearch *gs, double lon, double lat, double *prange, size_t nnn, size_t *adds, double *dist);
void gridsearch_extrapolate(struct gridsearch *gs);
void gridsearch_bound_poly(struct gridsearch *gs, size_t dims[2],  size_t n, const double *restrict lons, const double *restrict lats);

#endif

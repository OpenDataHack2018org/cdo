#ifndef  GRID_SEARCH_H_
#define  GRID_SEARCH_H_

#include <stdbool.h>
#include "kdtreelib/kdtree.h"
#include "nearpt3c.h"
#include "nanoflann.hpp"

#define GS_NOT_FOUND  SIZE_MAX

enum T_GRIDSEARCH_METHOD_NN  {GS_FULL=1, GS_NANOFLANN, GS_KDTREE, GS_KDSPH, GS_NEARPT3};

template <typename T>
struct PointCloud
{
  struct Point
  {
    T  x,y,z;
  };

  std::vector<Point>  pts;

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
  bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
};

typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<float, PointCloud<float> > ,
    PointCloud<float>,
    3 /* dim */
    > nfTree_t;

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
  Coord_T **pts;
  void *nearpt3;
};

struct gridsearch {
  bool extrapolate;
  bool is_cyclic;
  bool is_reg2d;
  int method_nn;
  size_t n;
  size_t dims[2];

  struct gsNear *near;
  struct kdTree *kdt;
  nfTree_t *nft;
  struct gsFull *full;

  double search_radius;

  // reg2d search
  double *reg2d_center_lon, *reg2d_center_lat;
  double *coslat, *sinlat;   // cosine, sine of grid lats (for distance)
  double *coslon, *sinlon;   // cosine, sine of grid lons (for distance)

  double lonmin, lonmax, latmin, latmax;
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

struct gridsearch *gridsearch_create_reg2d(bool is_cyclic, size_t dims[2], const double *restrict lons, const double *restrict lats);
struct gridsearch *gridsearch_create(size_t n, const double *restrict lons, const double *restrict lats);
struct gridsearch *gridsearch_create_nn(size_t n, const double *restrict lons, const double *restrict lats);
void gridsearch_delete(struct gridsearch *gs);
size_t gridsearch_nearest(struct gridsearch *gs, double lon, double lat, double *range);
struct pqueue *gridsearch_qnearest(struct gridsearch *gs, double lon, double lat, double *prange, size_t nnn);
void gridsearch_extrapolate(struct gridsearch *gs);

#endif

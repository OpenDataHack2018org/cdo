#ifndef _GRID_SEARCH_H_
#define _GRID_SEARCH_H_

#include <stdbool.h>
#include <limits.h>
#include "kdtreelib/kdtree.h"
#include "nearpt3c.h"

#define GS_NOT_FOUND  INT_MAX


enum T_GRIDSEARCH_METHOD_NN  {GS_FULL=1, GS_KDTREE, GS_NEARPT3};

struct gsFull {
  unsigned n;
  const double *plons;
  const double *plats;
  float **pts;
};

struct gsNear {
  unsigned n;
  const double *plons;
  const double *plats;
  Coord_T **pts;
  void *nearpt3;
};

struct gridsearch {
  int method_nn;
  unsigned n;
  unsigned nx, ny;

  struct gsNear *near;
  struct kdNode *kdt;
  struct gsFull *full;

  double search_radius;

  // reg2d search
  double *reg2d_center_lon, *reg2d_center_lat;
  double *coslat, *sinlat;   // cosine, sine of grid lats (for distance)
  double *coslon, *sinlon;   // cosine, sine of grid lons (for distance)
};

struct gsknn {
  unsigned ndist;
  unsigned size;
  int     *mask;
  int     *add;
  int     *tmpadd;
  double  *dist;
  double  *tmpdist;
};

struct gsknn *gridsearch_knn_new(unsigned size);
void gridsearch_knn_delete(struct gsknn *knn);
int gridsearch_knn(struct gridsearch *gs, struct gsknn *knn, double plon, double plat);

struct gridsearch *gridsearch_create_reg2d(bool lcyclic, unsigned nx, unsigned ny, const double *restrict lons, const double *restrict lats);
struct gridsearch *gridsearch_create(unsigned n, const double *restrict lons, const double *restrict lats);
struct gridsearch *gridsearch_create_nn(unsigned n, const double *restrict lons, const double *restrict lats);
void gridsearch_delete(struct gridsearch *gs);
unsigned gridsearch_nearest(struct gridsearch *gs, double lon, double lat, double *range);
struct pqueue *gridsearch_qnearest(struct gridsearch *gs, double lon, double lat, double *prange, unsigned nnn);

#endif

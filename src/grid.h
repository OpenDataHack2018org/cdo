#ifndef _GRID_H
#define _GRID_H

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <math.h>
#include <stdbool.h>

#ifndef  M_PI
#define  M_PI        3.14159265358979323846264338327950288  /* pi */
#endif


#ifndef  RAD2DEG
#define  RAD2DEG  (180./M_PI)   /* conversion for rad to deg */
#endif

#ifndef  DEG2RAD
#define  DEG2RAD  (M_PI/180.)   /* conversion for deg to rad */
#endif


int nfc_to_nlat(int nfc, int ntr);
int nlat_to_ntr(int nlat);
int nlat_to_ntr_linear(int nlat);
int ntr_to_nlat(int ntr);
int ntr_to_nlat_linear(int ntr);
int nlat_to_nlon(int nlat);

void grid_copy_attributes(int gridID1, int gridID2);
void grid_copy_mapping(int gridID1, int gridID2);

void grid_def_param_laea(int gridID, double a, double lon0, double lat0);
void grid_def_param_sinu(int gridID);

bool grid_is_distance_generic(int gridID);

void grid_to_radian(const char *units, long nvals, double *restrict values, const char *description);
void grid_to_degree(const char *units, long nvals, double *restrict values, const char *description);

void grid_gen_corners(int n, const double* restrict vals, double* restrict corners);
void grid_gen_bounds(int n, const double *restrict vals, double *restrict bounds);
void grid_check_lat_borders(int n, double *ybounds);

void grid_gen_xbounds2D(int nx, int ny, const double *restrict xbounds, double *restrict xbounds2D);
void grid_gen_ybounds2D(int nx, int ny, const double *restrict ybounds, double *restrict ybounds2D);

void grid_cell_center_to_bounds_X2D(const char* xunitstr, long xsize, long ysize,
				    const double* restrict grid_center_lon, double* restrict grid_corner_lon, double dlon);
void grid_cell_center_to_bounds_Y2D(const char* yunitstr, long xsize, long ysize,
				    const double* restrict grid_center_lat, double* restrict grid_corner_lat);

int  gridWeights(int gridID, double *weights);
int  gridGenArea(int gridID, double *area);

int referenceToGrid(int gridID);
int gridToZonal(int gridID);
int gridToMeridional(int gridID);
int gridToUnstructured(int gridID, int lbounds);
int gridToUnstructuredSelecton(int gridID1, int selectionSize, int *selectionIndexList, int nocoords ,int nobounds);
int gridToCurvilinear(int gridID, int lbounds);
int gridCurvilinearToRegular(int gridID);
int gridToRegular(int gridID);
void field2regular(int gridID1, int gridID2, double missval, double *array, int nmiss, int lnearest);

/* GME grid */
struct cart {
  double x[3];
};

struct geo {
  double lon;
  double lat;
};

void correct_sinxvals(int xsize, int ysize, double *xvals);

struct cart gc2cc(struct geo *position);
void factorni(int kni, int *kni2, int *kni3);
void gme_grid_restore(double *p, int ni, int nd);
void gme_grid(int lbounds, int gridsize, double *rlon, double *rlat,
	      double *blon, double *blat, int *imask,
              int ni, int nd, int ni2, int ni3);

/* Rotated grid */
double lamrot_to_lam(double phis, double rlas, double polphi, double pollam, double polgam);
double phirot_to_phi(double phis, double rlas, double polphi, double polgam);
void usvs_to_uv(double us, double vs, double phi, double rla,
		double polphi, double pollam, double *u, double *v);


// Projection codes for proj_info structure:
#define PROJ_LATLON  0
#define PROJ_MERC    1
#define PROJ_LC      3
#define PROJ_PS      5


typedef struct {
  int       code;     // Integer code for projection type
  double    lat1;     // SW latitude (1,1) in degrees (-90->90N)
  double    lon1;     // SW longitude (1,1) in degrees (-180->180E)
  double    dx;       // Grid spacing in meters at truelats, used
                      // only for ps, lc, and merc projections
  double    dlat;     // Lat increment for lat/lon grids
  double    dlon;     // Lon increment for lat/lon grids
  double    stdlon;   // Longitude parallel to y-axis (-180->180E)
  double    truelat1; // First true latitude (all projections)
  double    truelat2; // Second true lat (LC only)
  double    cone;     // Cone factor for LC projections
  double    polei;    // Computed i-location of pole point
  double    polej;    // Computed j-location of pole point
  double    rsw;      // Computed radius to SW corner
  double    rebydx;   // Earth radius divided by dx
  int       hemi;     // 1 for NH, -1 for SH
  int       init;     // Flag to indicate if this struct is ready for use
} proj_info_t;

/* Lambert Conformal grid (new version) */
void map_set(int proj_code, double lat1, double lon1, double dx, double stdlon,
	     double truelat1, double truelat2, proj_info_t *proj);
void ijll_lc(double i, double j, proj_info_t proj, double *lat, double *lon);

/* Lambert Conformal grid (old version) */
/*
int W3FB12(double xi, double xj, double alat1, double elon1, double dx,
	   double elonv, double alatan, double *alat, double *elon);
*/

// Define a sampled grid of another grid
int cdo_define_sample_grid(int gridID, int sampleFactor);

// Define a sub-grid of another grid
int cdo_define_subgrid_grid(int gridSrcID, int gridIDcurvl, int subI0, int subI1, int subJ0, int subJ1);

#endif  /* _GRID_H */

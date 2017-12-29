#ifndef  REMAP_H
#define  REMAP_H

#include <stdint.h>
#include <math.h>

#ifndef  M_PI
#define  M_PI        3.14159265358979323846264338327950288  /* pi */
#endif

#define  PI       M_PI
#define  PI2      (2.0*PI)
#define  PIH      (0.5*PI)

#define  ZERO     0.0
#define  ONE      1.0
#define  TWO      2.0
#define  THREE    3.0
#define  HALF     0.5
#define  QUART    0.25
#define  BIGNUM   1.e+20
#define  TINY     1.e-14


#define  REMAP_GRID_TYPE_REG2D     1
#define  REMAP_GRID_TYPE_CURVE2D   2
#define  REMAP_GRID_TYPE_UNSTRUCT  3

#define  REMAP_GRID_BASIS_SRC      1
#define  REMAP_GRID_BASIS_TGT      2

#define  RESTR_TYPE  int  /* restrict data types: 0 -> double, float; 1 -> int */

typedef RESTR_TYPE restr_t;
/* short
#  define RESTR_SFAC     4000
#  define RESTR_SCALE(x) ((short) (0.5+RESTR_SFAC*(x)))
#  define RESTR_ABS(x)   abs(x)
*/
/* int */
#  define RESTR_SFAC     100000000
#  define RESTR_SCALE(x) ((int) (0.5+RESTR_SFAC*(x)))
#  define RESTR_ABS(x)   abs(x)
/*
#  define RESTR_SFAC     1.
#  define RESTR_SCALE(x) (x)
#  define RESTR_ABS(x)   fabs(x)
*/

#define  TINY_FRAC     1.e-10

enum struct RemapType   {UNDEF, BILINEAR, BICUBIC, DISTWGT, CONSERV, CONSERV_YAC};
enum struct SubmapType  {NONE, LAF, SUM};
enum struct NormOpt     {NONE, DESTAREA, FRACAREA};


typedef struct {
  int      gridID;
  int      remap_grid_type;
  int      rank;                  /* rank of the grid */
  size_t   size;                  /* total points on the grid */
  size_t   num_cell_corners;      /* number of corners for each grid cell */

  bool     lneed_cell_corners;
  bool     luse_cell_corners;     /* use corners for bounding boxes  */

  bool     lextrapolate;
  bool     non_global;
  bool     is_cyclic;

  size_t   dims[2];               /* size of grid dimension */

  int      nvgp;                  /* size of vgpm           */
  int*     vgpm;                  /* flag which cells are valid   */

  int*     mask;                  /* flag which cells participate */

  double*  reg2d_center_lon;      /* reg2d lon/lat coordinates for */
  double*  reg2d_center_lat;      /* each grid center in radians   */
  double*  reg2d_corner_lon;      /* reg2d lon/lat coordinates for */
  double*  reg2d_corner_lat;      /* each grid corner in radians   */

  double*  cell_center_lon;       /* lon/lat coordinates for       */
  double*  cell_center_lat;       /* each grid center in radians   */
  double*  cell_corner_lon;       /* lon/lat coordinates for       */
  double*  cell_corner_lat;       /* each grid corner in radians   */

  double*  cell_area;             /* tot area of each grid cell    */
  double*  cell_frac;             /* fractional area of grid cells participating in remapping  */

  restr_t *cell_bound_box;        /* lon/lat bounding box for use    */
  int      num_srch_bins;         /* num of bins for restricted srch */
  size_t  *bin_addr;              /* min,max adds for grid cells in this lat bin  */
  restr_t* bin_lats;              /* min,max latitude for each search bin   */
}
remapgrid_t;

typedef struct {
  bool     option;
  size_t   max_links;
  size_t   num_blks;
  size_t  *num_links;
  size_t **src_add;
  size_t **dst_add;
  size_t **w_index;
}
remaplink_t;

typedef struct {
  long      links_per_value;
  bool      sort_add;
  bool      pinit;            /* true: if the pointers are initialized    */
  size_t    max_links;        /* current size of link arrays              */
  size_t    num_links;        /* actual number of links for remapping     */
  size_t    num_wts;          /* num of weights used in remapping         */
  RemapType mapType;          /* identifier for remapping method          */
  NormOpt   normOpt;          /* option for normalization (conserv only)  */
  size_t    resize_increment; /* default amount to increase array size    */

  size_t   *src_cell_add;     /* source grid address for each link        */
  size_t   *tgt_cell_add;     /* target grid address for each link        */

  double   *wts;              /* map weights for each link [max_links*num_wts] */

  remaplink_t links;
}
remapvars_t;

typedef struct {
  int      nused;
  int      gridID;
  size_t   gridsize;
  size_t   nmiss;
  remapgrid_t src_grid;
  remapgrid_t tgt_grid;
  remapvars_t vars;
}
remap_t;

#define  REMAP_STORE_LINK_FAST  1
#define  REMAP_WRITE_REMAP      2
#define  REMAP_MAX_ITER         3
#define  REMAP_NUM_SRCH_BINS    4
#define  REMAP_GENWEIGHTS       5

void remap_set_threshhold(double threshhold);
void remap_set_int(int remapvar, int value);


void remap_grids_init(RemapType mapType, bool lextrapolate, int gridID1, remapgrid_t *src_grid, int gridID2, remapgrid_t *tgt_grid);
void remap_vars_init(RemapType mapType, size_t src_grid_size, size_t tgt_grid_size, remapvars_t *rv);

void remapVarsFree(remapvars_t *rv);
void remapGridFree(remapgrid_t *grid);

void remap(double *restrict dst_array, double missval, size_t dst_size, size_t num_links, double *restrict map_wts, 
	   size_t num_wts, const size_t *restrict dst_add, const size_t *restrict src_add, const double *restrict src_array, 
	   const double *restrict src_grad1, const double *restrict src_grad2, const double *restrict src_grad3,
	   remaplink_t links, long links_per_value);

void remap_laf(double *restrict dst_array, double missval, size_t dst_size, size_t num_links, double *restrict map_wts,
	       size_t num_wts, const size_t *restrict dst_add, const size_t *restrict src_add, const double *restrict src_array);

void remap_sum(double *restrict dst_array, double missval, size_t dst_size, size_t num_links, double *restrict map_wts,
	       size_t num_wts, const size_t *restrict dst_add, const size_t *restrict src_add, const double *restrict src_array);

void scrip_remap_bilinear_weights(remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv);
void scrip_remap_bicubic_weights(remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv);
void remap_distwgt_weights(size_t num_neighbors, remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv);
void scrip_remap_conserv_weights(remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv);
void remap_conserv_weights(remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv);

void scrip_remap_bilinear(remapgrid_t *src_grid, remapgrid_t *tgt_grid, const double *restrict src_array, double *restrict tgt_array, double missval);
void scrip_remap_bicubic(remapgrid_t *src_grid, remapgrid_t *tgt_grid, const double *restrict src_array, double *restrict tgt_array, double missval);
void remap_distwgt(size_t num_neighbors, remapgrid_t *src_grid, remapgrid_t *tgt_grid, const double *restrict src_array, double *restrict tgt_array, double missval);
void remap_conserv(remapgrid_t *src_grid, remapgrid_t *tgt_grid, const double *restrict src_array, double *restrict tgt_array, double missval);


void resize_remap_vars(remapvars_t *rv, int64_t increment);

void remap_stat(int remap_order, remapgrid_t src_grid, remapgrid_t tgt_grid, remapvars_t rv, const double *restrict array1, 
		const double *restrict array2, double missval);
void remap_gradients(remapgrid_t grid, const double *restrict array, double *restrict grad_lat,
		     double *restrict grad_lon, double *restrict grad_latlon);

void reorder_links(remapvars_t *rv);

void sort_add(size_t num_links, size_t num_wts, size_t *restrict add1, size_t *restrict add2, double *restrict weights);
void sort_iter(size_t num_links, size_t num_wts, size_t *restrict add1, size_t *restrict add2, double *restrict weights, int parent);

void write_remap_scrip(const char *interp_file, RemapType mapType, SubmapType submapType, int num_neighbors,
		       int remap_order, remapgrid_t src_grid, remapgrid_t tgt_grid, remapvars_t rv);
void read_remap_scrip(const char *interp_file, int gridID1, int gridID2, RemapType *mapType, SubmapType *submapType, int *num_neighbors,
		      int *remap_order, remapgrid_t *src_grid, remapgrid_t *tgt_grid, remapvars_t *rv);

void calc_lat_bins(remapgrid_t* src_grid, remapgrid_t* tgt_grid, RemapType mapType);
size_t get_srch_cells(size_t tgt_cell_add, size_t nbins, size_t *bin_addr1, size_t *bin_addr2,
                      restr_t *tgt_cell_bound_box, restr_t *src_cell_bound_box, size_t src_grid_size, size_t *srch_add);

int grid_search_reg2d_nn(size_t nx, size_t ny, size_t *restrict nbr_add, double *restrict nbr_dist, double plat, double plon,
                         const double *restrict src_center_lat, const double *restrict src_center_lon);

int grid_search_reg2d(remapgrid_t *src_grid, size_t *restrict src_add, double *restrict src_lats, 
                      double *restrict src_lons,  double plat, double plon, const size_t *restrict src_grid_dims,
                      const double *restrict src_center_lat, const double *restrict src_center_lon);

bool point_in_quad(bool is_cyclic, size_t nx, size_t ny, size_t i, size_t j, size_t adds[4], double lons[4], double lats[4],
                   double plon, double plat, const double *restrict center_lon, const double *restrict center_lat);

int grid_search(remapgrid_t *src_grid, size_t *restrict src_add, double *restrict src_lats, 
		double *restrict src_lons,  double plat, double plon, const size_t *restrict src_grid_dims,
		const double *restrict src_center_lat, const double *restrict src_center_lon,
		const restr_t *restrict src_grid_bound_box, const size_t *restrict src_bin_add);

bool find_ij_weights(double plon, double plat, double* restrict src_lons, double* restrict src_lats, double *ig, double *jg);
int rect_grid_search(size_t *ii, size_t *jj, double x, double y, size_t nxm, size_t nym, const double *restrict xm, const double *restrict ym);

void remapgrid_get_lonlat(remapgrid_t *grid, size_t cell_add, double *plon, double *plat);

void remapCheckArea(size_t grid_size, double *restrict cell_area, const char *name);
void remapCheckWeights(size_t num_links, size_t num_wts, NormOpt normOpt, size_t *src_cell_add, size_t *tgt_cell_add, double *wts);

#endif  /* REMAP_H */

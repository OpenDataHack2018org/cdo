

#define  NORM_OPT_NONE      1
#define  NORM_OPT_DESTAREA  2
#define  NORM_OPT_FRACAREA  3

#define  MAP_TYPE_CONSERV   1
#define  MAP_TYPE_BILINEAR  2
#define  MAP_TYPE_BICUBIC   3
#define  MAP_TYPE_DISTWGT   4
#define  MAP_TYPE_DISTWGT1  5

#define  SUBMAP_TYPE_NONE   0
#define  SUBMAP_TYPE_LAF    1

#define  RESTRICT_LATITUDE  1
#define  RESTRICT_LATLON    2


typedef struct {
  int      pinit;            /* TRUE if the pointers are initialized     */
  int      gridID1;
  int      gridID2;
  int      lextrapolate;
  int      non_global;
  int      grid1_is_cyclic, grid2_is_cyclic;
  int      grid1_size, grid2_size; /* total points on each grid */
  int      grid1_rank, grid2_rank; /* rank of each grid */
  int      grid1_corners, grid2_corners; /* number of corners for each grid cell */

  int      grid1_dims[2], grid2_dims[2]; /* size of each grid dimension */

  int      grid1_nvgp;         /* size of grid1_vgpm           */
  int     *grid1_vgpm;         /* flag which cells are valid   */

  int      grid2_nvgp;         /* size of grid2_vgpm           */
  int     *grid2_vgpm;         /* flag which cells are valid   */

  int     *grid1_mask;         /* flag which cells participate */
  int     *grid2_mask;         /* flag which cells participate */

  double  *grid1_center_lat;   /* lat/lon coordinates for      */
  double  *grid1_center_lon;   /* each grid center in radians  */
  double  *grid2_center_lat; 
  double  *grid2_center_lon;
  double  *grid1_area;         /* tot area of each grid1 cell     */
  double  *grid2_area;         /* tot area of each grid2 cell     */
  /* double  *grid1_area_in; */     /* area of grid1 cell from file    */
  /* double  *grid2_area_in; */     /* area of grid2 cell from file    */
  double  *grid1_frac;         /* fractional area of grid cells   */
  double  *grid2_frac;         /* participating in remapping      */

  double  *grid1_corner_lat;   /* lat/lon coordinates for         */
  double  *grid1_corner_lon;   /* each grid corner in radians     */
  double  *grid2_corner_lat; 
  double  *grid2_corner_lon;

  int      lneed_grid1_corners;
  int      lneed_grid2_corners;
  int      luse_grid1_corners;  /* use corners for bounding boxes  */
  int      luse_grid2_corners;  /* use corners for bounding boxes  */
  /* int      luse_grid1_area;   */ /* use area from grid file         */
  /* int      luse_grid2_area;   */ /* use area from grid file         */

  double  *grid1_bound_box;    /* lat/lon bounding box for use    */
  double  *grid2_bound_box;    /* in restricting grid searches    */

  int      restrict_type;
  int      num_srch_bins;      /* num of bins for restricted srch */

  int     *bin_addr1;       /* min,max adds for grid1 cells in this lat bin  */
  int     *bin_addr2;       /* min,max adds for grid2 cells in this lat bin  */

  double  *bin_lats;        /* min,max latitude for each search bin   */
  double  *bin_lons;        /* min,max longitude for each search bin  */
}
remapgrid_t;

typedef struct {
  int   option;
  int   max_links;
  int   num_blks;
  int  *num_links;
  int **src_add;
  int **dst_add;
  int **w_index;
}
remaplink_t;

typedef struct {
  int   pinit;            /* TRUE if the pointers are initialized     */
  int   max_links;        /* current size of link arrays              */
  int   num_links;        /* actual number of links for remapping     */
  int   num_wts;          /* num of weights used in remapping         */
  int   map_type;         /* identifier for remapping method          */
  int   norm_opt;         /* option for normalization (conserv only)  */
  int   resize_increment; /* default amount to increase array size    */

  int  *grid1_add;        /* grid1 address for each link              */
  int  *grid2_add;        /* grid2 address for each link              */

  double *wts[4];         /* map weights for each link [max_links][num_wts] */

  remaplink_t  links;
}
remapvars_t;

typedef struct {
  int gridID;
  int gridsize;
  int nmiss;
  remapgrid_t grid;
  remapvars_t vars;
}
remap_t;

void remap_set_max_iter(int max_iter);

void remapGridInit(int map_type, int lextrapolate, int gridID1, int gridID2, remapgrid_t *rg);
void remapVarsInit(int map_type, remapgrid_t *rg, remapvars_t *rv);

void remapVarsFree(remapvars_t *rv);
void remapGridFree(remapgrid_t *rg);

void remap(double *dst_array, double missval, int dst_size, int dst_array_dim,
	   double **map_wts, int map_wts_dim,
	   int *dst_add, int *src_add, 
	   double *src_array, double *src_grad1, double *src_grad2, double *src_grad3,
	   remaplink_t links);

void remap_laf(double *dst_array, double missval, int dst_size, int num_links, double **map_wts,
	       int *dst_add, int *src_add, double *src_array);

void remap_bilin(remapgrid_t *rg, remapvars_t *rv);
void remap_bicub(remapgrid_t *rg, remapvars_t *rv);
void remap_conserv(remapgrid_t *rg, remapvars_t *rv);
void remap_distwgt(remapgrid_t *rg, remapvars_t *rv);
void remap_distwgt1(remapgrid_t *rg, remapvars_t *rv);

void resize_remap_vars(remapvars_t *rv, int increment);

void remap_stat(int remap_order, remapgrid_t rg, remapvars_t rv, double *array1, double *array2, double missval);
void remap_gradients(remapgrid_t rg, double *array, double *grad1_lat,
		     double *grad1_lon, double *grad1_latlon);

void reorder_links(remapvars_t *rv);

void sort_add(int num_links, int num_wts, int *add1, int *add2, double **weights);

void write_remap_scrip(const char *interp_file, int map_type, int submap_type, 
		       int remap_order, remapgrid_t rg, remapvars_t rv);
void read_remap_scrip(const char *interp_file, int gridID1, int gridID2, int *map_type, int *submap_type,
		      int *remap_order, remapgrid_t *rg, remapvars_t *rv);

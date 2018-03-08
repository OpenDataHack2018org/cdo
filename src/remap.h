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
#ifndef REMAP_H
#define REMAP_H

#include <stdint.h>
#include <cmath>
#include "remap_vars.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288 /* pi */
#endif

constexpr double PI = M_PI;
constexpr double PI2 = (2.0 * PI);
constexpr double PIH = (0.5 * PI);
constexpr float  PI_f = PI;
constexpr float  PI2_f = PI2;
constexpr float  PIH_f = PIH;

#define ZERO 0.0
#define ONE 1.0
#define TWO 2.0
#define THREE 3.0
#define HALF 0.5
#define QUART 0.25
#define TINY 1.e-14

#define REMAP_GRID_TYPE_REG2D 1
#define REMAP_GRID_TYPE_CURVE2D 2
#define REMAP_GRID_TYPE_UNSTRUCT 3

#define REMAP_GRID_BASIS_SRC 1
#define REMAP_GRID_BASIS_TGT 2

enum struct SubmapType
{
  NONE,
  LAF,
  SUM
};

struct RemapGridType
{
  int gridID;
  int remap_grid_type;
  int rank;                 // rank of the grid
  size_t size;              // total points on the grid
  size_t num_cell_corners;  // number of corners for each grid cell

  bool lneed_cell_corners;
  bool luse_cell_corners;   // use corners for bounding boxes

  bool lextrapolate;
  bool non_global;
  bool is_cyclic;

  size_t dims[2];           // size of grid dimension

  int nvgp;                 // size of vgpm
  int *vgpm;                // flag which cells are valid

  int *mask;                // flag which cells participate

  double *reg2d_center_lon; // reg2d lon/lat coordinates for
  double *reg2d_center_lat; // each grid center in radians
  double *reg2d_corner_lon; // reg2d lon/lat coordinates for
  double *reg2d_corner_lat; // each grid corner in radians

  double *cell_center_lon;  // lon/lat coordinates for
  double *cell_center_lat;  // each grid center in radians
  double *cell_corner_lon;  // lon/lat coordinates for
  double *cell_corner_lat;  // each grid corner in radians

  double *cell_area;        // tot area of each grid cell
  double *cell_frac;        // fractional area of grid cells participating in remapping

  int num_srch_bins;        // num of bins for restricted srch
  size_t *bin_addr;         // min,max adds for grid cells in this lat bin
  float *bin_lats;          // min,max latitude for each search bin
  float *cell_bound_box;    // lon/lat bounding box for use
};


struct remapType
{
  int nused;
  int gridID;
  size_t gridsize;
  size_t nmiss;
  RemapGridType src_grid;
  RemapGridType tgt_grid;
  RemapVarsType vars;
};

#define REMAP_WRITE_REMAP 2
#define REMAP_MAX_ITER 3
#define REMAP_NUM_SRCH_BINS 4
#define REMAP_GENWEIGHTS 5

void remap_set_threshhold(double threshhold);
void remap_set_int(int remapvar, int value);

void remapInitGrids(RemapType mapType, bool lextrapolate, int gridID1, RemapGridType &src_grid, int gridID2,
                    RemapGridType &tgt_grid);

void remapGridInit(RemapGridType &grid);
void remapGridFree(RemapGridType &grid);
void remapGridAlloc(RemapType mapType, RemapGridType &grid);

void scrip_remap_bilinear_weights(RemapGridType *src_grid, RemapGridType *tgt_grid, RemapVarsType &rv);
void scrip_remap_bicubic_weights(RemapGridType *src_grid, RemapGridType *tgt_grid, RemapVarsType &rv);
void remap_distwgt_weights(size_t numNeighbors, RemapGridType *src_grid, RemapGridType *tgt_grid, RemapVarsType &rv);
void scrip_remap_conserv_weights(RemapGridType *src_grid, RemapGridType *tgt_grid, RemapVarsType &rv);
void remap_conserv_weights(RemapGridType *src_grid, RemapGridType *tgt_grid, RemapVarsType &rv);

void scrip_remap_bilinear(RemapGridType *src_grid, RemapGridType *tgt_grid, const double *restrict src_array,
                          double *restrict tgt_array, double missval);
void scrip_remap_bicubic(RemapGridType *src_grid, RemapGridType *tgt_grid, const double *restrict src_array,
                         double *restrict tgt_array, double missval);
void remap_distwgt(size_t numNeighbors, RemapGridType *src_grid, RemapGridType *tgt_grid, const double *restrict src_array,
                   double *restrict tgt_array, double missval);
void remap_conserv(RemapGridType *src_grid, RemapGridType *tgt_grid, const double *restrict src_array,
                   double *restrict tgt_array, double missval);

void remap_stat(int remap_order, RemapGridType &src_grid, RemapGridType &tgt_grid, RemapVarsType &rv,
                const double *restrict array1, const double *restrict array2, double missval);
void remap_gradients(RemapGridType &grid, const double *restrict array, gradientsType &gradients);

void sort_add(size_t num_links, size_t num_wts, size_t *restrict add1, size_t *restrict add2, double *restrict weights);
void sort_iter(size_t num_links, size_t num_wts, size_t *restrict add1, size_t *restrict add2, double *restrict weights,
               int parent);

void write_remap_scrip(const char *interp_file, RemapType mapType, SubmapType submapType, int numNeighbors,
                       int remap_order, RemapGridType &src_grid, RemapGridType &tgt_grid, RemapVarsType &rv);
void read_remap_scrip(const char *interp_file, int gridID1, int gridID2, RemapType *mapType, SubmapType *submapType,
                      int *numNeighbors, int *remap_order, RemapGridType &src_grid, RemapGridType &tgt_grid,
                      RemapVarsType &rv);

void calc_lat_bins(RemapGridType &src_grid, RemapGridType &tgt_grid, RemapType mapType);
size_t get_srch_cells(size_t tgt_cell_add, size_t nbins, size_t *bin_addr1, size_t *bin_addr2,
                      float *tgt_cell_bound_box, float *src_cell_bound_box, size_t src_grid_size, size_t *srch_add);

int grid_search_reg2d_nn(size_t nx, size_t ny, size_t *restrict nbr_add, double *restrict nbr_dist, double plat,
                         double plon, const double *restrict src_center_lat, const double *restrict src_center_lon);

int grid_search_reg2d(RemapGridType *src_grid, size_t *restrict src_add, double *restrict src_lats,
                      double *restrict src_lons, double plat, double plon, const size_t *restrict src_grid_dims,
                      const double *restrict src_center_lat, const double *restrict src_center_lon);

bool point_in_quad(bool is_cyclic, size_t nx, size_t ny, size_t i, size_t j, size_t adds[4], double lons[4],
                   double lats[4], double plon, double plat, const double *restrict center_lon,
                   const double *restrict center_lat);

int grid_search(RemapGridType *src_grid, size_t *restrict src_add, double *restrict src_lats, double *restrict src_lons,
                double plat, double plon, const size_t *restrict src_grid_dims, const double *restrict src_center_lat,
                const double *restrict src_center_lon, const float *restrict src_grid_bound_box,
                const size_t *restrict src_bin_add);

bool find_ij_weights(double plon, double plat, double *restrict src_lons, double *restrict src_lats, double *ig,
                     double *jg);
int rect_grid_search(size_t *ii, size_t *jj, double x, double y, size_t nxm, size_t nym, const double *restrict xm,
                     const double *restrict ym);

void remapgrid_get_lonlat(RemapGridType *grid, size_t cell_add, double *plon, double *plat);

void remapCheckArea(size_t grid_size, double *restrict cell_area, const char *name);

#endif /* REMAP_H */

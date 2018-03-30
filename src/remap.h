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
#include "grid_search.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

constexpr double PI = M_PI;
constexpr double PI2 = (2.0 * PI);
constexpr double PIH = (0.5 * PI);
constexpr float PI_f = PI;
constexpr float PI2_f = PI2;
constexpr float PIH_f = PIH;

#define TINY 1.e-14

#define REMAP_GRID_TYPE_REG2D 1
#define REMAP_GRID_TYPE_CURVE2D 2
#define REMAP_GRID_TYPE_UNSTRUCT 3
// enum struct RemapGridType {REG2D, CURVE2D, UNSTRUCT}

#define REMAP_GRID_BASIS_SRC 1
#define REMAP_GRID_BASIS_TGT 2

enum struct SubmapType
{
  NONE,
  LAF,
  SUM
};

struct RemapGrid
{
  int gridID;
  int remap_grid_type;
  int rank;                 // rank of the grid
  size_t size;              // total points on the grid
  size_t num_cell_corners;  // number of corners for each grid cell

  bool lneed_cell_corners;
  bool luse_cell_corners;  // use corners for bounding boxes

  bool lextrapolate;
  bool non_global;
  bool is_cyclic;

  size_t dims[2];  // size of grid dimension

  int nvgp;               // size of vgpm
  std::vector<int> vgpm;  // flag which cells are valid

  std::vector<int> mask;  // flag which cells participate

  double *reg2d_center_lon;  // reg2d lon/lat coordinates for
  double *reg2d_center_lat;  // each grid center in radians
  double *reg2d_corner_lon;  // reg2d lon/lat coordinates for
  double *reg2d_corner_lat;  // each grid corner in radians

  double *cell_center_lon;  // lon/lat coordinates for
  double *cell_center_lat;  // each grid center in radians
  double *cell_corner_lon;  // lon/lat coordinates for
  double *cell_corner_lat;  // each grid corner in radians

  std::vector<double> cell_area;  // total area of each grid cell
  std::vector<double> cell_frac;  // fractional area of grid cells participating in remapping
};

struct GridSearchBins
{
  unsigned nbins;                     // num of bins for restricted search
  size_t ncells;                      // total number of grid cells (cell_bound_box)
  std::vector<size_t> bin_addr;       // min,max adds for grid cells in this lat bin
  std::vector<float> bin_lats;        // min,max latitude for each search bin
  std::vector<float> cell_bound_box;  // lon/lat bounding box for use
};

struct RemapSearch
{
  RemapGrid *srcGrid;
  RemapGrid *tgtGrid;

  GridSearchBins srcBins;
  GridSearchBins tgtBins;

  GridSearch *gs;
};

struct remapType
{
  int nused;
  int gridID;
  size_t gridsize;
  size_t nmiss;
  RemapGrid src_grid;
  RemapGrid tgt_grid;
  RemapVars vars;
  RemapSearch search;
};

#define REMAP_WRITE_REMAP 2
#define REMAP_MAX_ITER 3
#define REMAP_NUM_SRCH_BINS 4
#define REMAP_GENWEIGHTS 5

void remap_set_threshhold(double threshhold);
void remap_set_int(int remapvar, int value);

void remapInitGrids(RemapMethod mapType, bool lextrapolate, int gridID1, RemapGrid &src_grid, int gridID2, RemapGrid &tgt_grid);

void remapGridInit(RemapGrid &grid);
void remapGridFree(RemapGrid &grid);
void remapGridAlloc(RemapMethod mapType, RemapGrid &grid);
void remapSearchInit(RemapMethod mapType, RemapSearch &search, RemapGrid &src_grid, RemapGrid &tgt_grid);
void remapSearchFree(RemapSearch &search);

void remapSearchPoints(RemapSearch &rsearch, double plon, double plat, knnWeightsType &knnWeights);
int remapSearchSquare(RemapSearch &rsearch, double plon, double plat, size_t *src_add, double *src_lats, double *src_lons);

void remapBilinearWeights(RemapSearch &rsearch, RemapVars &rv);
void remapBicubicWeights(RemapSearch &rsearch, RemapVars &rv);
void remapDistwgtWeights(size_t numNeighbors, RemapSearch &rsearch, RemapVars &rv);
void remapConservWeights(RemapSearch &rsearch, RemapVars &rv);
void remapConservWeightsScrip(RemapSearch &rsearch, RemapVars &rv);

void remapBilinear(RemapSearch &rsearch, const double *restrict src_array, double *restrict tgt_array, double missval);
void remapBicubic(RemapSearch &rsearch, const double *restrict src_array, double *restrict tgt_array, double missval);
void remapDistwgt(size_t numNeighbors, RemapSearch &rsearch, const double *restrict src_array, double *restrict tgt_array,
                  double missval);
void remapConserv(RemapSearch &rsearch, const double *restrict src_array, double *restrict tgt_array, double missval);

void remapStat(int remapOrder, RemapGrid &src_grid, RemapGrid &tgt_grid, RemapVars &rv, const double *restrict array1,
               const double *restrict array2, double missval);
void remapGradients(RemapGrid &grid, const double *restrict array, gradientsType &gradients);

void sort_add(size_t num_links, size_t num_wts, size_t *restrict add1, size_t *restrict add2, double *restrict weights);
void sort_iter(size_t num_links, size_t num_wts, size_t *restrict add1, size_t *restrict add2, double *restrict weights,
               int parent);

void remapWriteDataScrip(const char *interp_file, RemapMethod mapType, SubmapType submapType, int numNeighbors, int remapOrder,
                         RemapGrid &src_grid, RemapGrid &tgt_grid, RemapVars &rv);
void remapReadDataScrip(const char *interp_file, int gridID1, int gridID2, RemapMethod *mapType, SubmapType *submapType,
                        int *numNeighbors, int *remapOrder, RemapGrid &src_grid, RemapGrid &tgt_grid, RemapVars &rv);

void calc_lat_bins(GridSearchBins &searchBins);
size_t get_srch_cells(size_t tgt_cell_addr, GridSearchBins &tgtBins, GridSearchBins &srcBins, float *tgt_cell_bound_box,
                      size_t *srch_add);

int gridSearchSquareReg2dNN(size_t nx, size_t ny, size_t *restrict nbr_add, double *restrict nbr_dist, double plat, double plon,
                            const double *restrict src_center_lat, const double *restrict src_center_lon);

int gridSearchSquareReg2d(RemapGrid *src_grid, size_t *restrict src_add, double *restrict src_lats, double *restrict src_lons,
                          double plat, double plon);

bool pointInQuad(bool isCyclic, size_t nx, size_t ny, size_t i, size_t j, size_t adds[4], double lons[4], double lats[4],
                 double plon, double plat, const double *restrict centerLon, const double *restrict centerLat);

int gridSearchSquareCurv2dScrip(RemapGrid *src_grid, size_t *restrict src_add, double *restrict src_lats,
                                double *restrict src_lons, double plat, double plon, GridSearchBins &srcBins);

bool remapFindWeights(double plon, double plat, double *restrict src_lons, double *restrict src_lats, double *ig, double *jg);
int rect_grid_search(size_t *ii, size_t *jj, double x, double y, size_t nxm, size_t nym, const double *restrict xm,
                     const double *restrict ym);

void remapgrid_get_lonlat(RemapGrid *grid, size_t cell_add, double *plon, double *plat);

void remapCheckArea(size_t grid_size, double *restrict cell_area, const char *name);

#endif /* REMAP_H */

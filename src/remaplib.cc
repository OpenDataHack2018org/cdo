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
/*
  This is a C library of the Fortran SCRIP version 1.4

  ===>>> Please send bug reports to <http://mpimet.mpg.de/cdo> <<<===

  Spherical Coordinate Remapping and Interpolation Package (SCRIP)
  ================================================================

  SCRIP is a software package which computes addresses and weights for
  remapping and interpolating fields between grids in spherical coordinates.
  It was written originally for remapping fields to other grids in a coupled
  climate model, but is sufficiently general that it can be used in other
  applications as well. The package should work for any grid on the surface
  of a sphere. SCRIP currently supports four remapping options:

  Conservative remapping
  ----------------------
  First- and second-order conservative remapping as described in
  Jones (1999, Monthly Weather Review, 127, 2204-2210).

  Bilinear interpolation
  ----------------------
  Slightly generalized to use a local bilinear approximation
  (only logically-rectangular grids).

  Bicubic interpolation
  ----------------------
  Similarly generalized (only logically-rectangular grids).

  Distance-weighted averaging
  ---------------------------
  Distance-weighted average of a user-specified number of nearest neighbor
  values.

  Documentation
  =============

  http://climate.lanl.gov/Software/SCRIP/SCRIPusers.pdf

*/
/*
  2013-11-08 Uwe Schulzweida: split remapgrid class to src_grid and tgt_grid
  2012-01-16 Uwe Schulzweida: alloc grid2_bound_box only for conservative
  remapping 2011-01-07 Uwe Schulzweida: Changed remap weights from 2D to 1D
  array 2009-05-25 Uwe Schulzweida: Changed restrict data type from double to
  int 2009-01-11 Uwe Schulzweida: OpenMP parallelization
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cdi.h>

#include "cdo_int.h"
#include "grid.h"
#include "remap.h"

#define IS_REG2D_GRID(gridID) (gridInqType(gridID) == GRID_LONLAT || gridInqType(gridID) == GRID_GAUSSIAN)

static bool remap_gen_weights = true;
static bool remap_write_remap = false;
static int remap_num_srch_bins = 180;
#define DEFAULT_MAX_ITER 100
size_t remap_max_iter = DEFAULT_MAX_ITER;  // Max iteration count for i, j iteration

void
remap_set_int(int remapvar, int value)
{
  if (remapvar == REMAP_WRITE_REMAP)
    remap_write_remap = value > 0;
  else if (remapvar == REMAP_MAX_ITER)
    remap_max_iter = value;
  else if (remapvar == REMAP_NUM_SRCH_BINS)
    remap_num_srch_bins = value;
  else if (remapvar == REMAP_GENWEIGHTS)
    remap_gen_weights = value > 0;
  else
    cdoAbort("Unsupported remap variable (%d)!", remapvar);
}

/*****************************************************************************/

void
remapGridAlloc(RemapType mapType, RemapGridType &grid)
{
  if (grid.nvgp) grid.vgpm = (int *) Malloc(grid.nvgp * sizeof(int));

  grid.mask = (int *) Malloc(grid.size * sizeof(int));

  if (remap_write_remap || grid.remap_grid_type != REMAP_GRID_TYPE_REG2D)
    {
      grid.cell_center_lon = (double *) Malloc(grid.size * sizeof(double));
      grid.cell_center_lat = (double *) Malloc(grid.size * sizeof(double));
    }

  if (mapType == RemapType::CONSERV || mapType == RemapType::CONSERV_YAC)
    {
      grid.cell_area = (double *) Calloc(grid.size, sizeof(double));
    }

  grid.cell_frac = (double *) Calloc(grid.size, sizeof(double));

  if (grid.lneed_cell_corners)
    {
      if (grid.num_cell_corners > 0)
        {
          size_t nalloc = grid.num_cell_corners * grid.size;
          grid.cell_corner_lon = (double *) Calloc(nalloc, sizeof(double));
          grid.cell_corner_lat = (double *) Calloc(nalloc, sizeof(double));
        }
    }
}

/*****************************************************************************/
static void
boundbox_from_corners(size_t size, size_t nc, const double *restrict corner_lon, const double *restrict corner_lat,
                      float *restrict bound_box)
{
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(bound_box, corner_lat, corner_lon, nc, size)
#endif
  for (size_t i = 0; i < size; ++i)
    {
      size_t i4 = i << 2;  // *4
      size_t inc = i * nc;
      float clat = corner_lat[inc];
      float clon = corner_lon[inc];
      bound_box[i4] = clat;
      bound_box[i4 + 1] = clat;
      bound_box[i4 + 2] = clon;
      bound_box[i4 + 3] = clon;
      for (size_t j = 1; j < nc; ++j)
        {
          clat = corner_lat[inc + j];
          clon = corner_lon[inc + j];
          if (clat < bound_box[i4]) bound_box[i4] = clat;
          if (clat > bound_box[i4 + 1]) bound_box[i4 + 1] = clat;
          if (clon < bound_box[i4 + 2]) bound_box[i4 + 2] = clon;
          if (clon > bound_box[i4 + 3]) bound_box[i4 + 3] = clon;
        }
    }
}

static void
boundbox_from_center(bool lonIsCyclic, size_t size, size_t nx, size_t ny, const double *restrict center_lon,
                     const double *restrict center_lat, float *restrict bound_box)
{
  size_t idx[4];
  float tmp_lats[4], tmp_lons[4];

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(lonIsCyclic, size, nx, ny, center_lon, center_lat, \
                                              bound_box) private(idx, tmp_lats, tmp_lons)
#endif
  for (size_t n = 0; n < size; n++)
    {
      size_t n4 = n << 2;

      /* Find N,S and NE points to this grid point */

      size_t j = n / nx;
      size_t i = n - j * nx;

      size_t ip1 = (i < (nx - 1)) ? i + 1 : lonIsCyclic ? 0 : i;
      size_t jp1 = (j < (ny - 1)) ? j + 1 : j;

      idx[0] = n;
      idx[1] = j * nx + ip1;    // east
      idx[2] = jp1 * nx + ip1;  // north-east
      idx[3] = jp1 * nx + i;    // north

      /* Find N,S and NE lat/lon coords and check bounding box */

      for (unsigned j = 0; j < 4; ++j) tmp_lons[j] = center_lon[idx[j]];
      for (unsigned j = 0; j < 4; ++j) tmp_lats[j] = center_lat[idx[j]];

      bound_box[n4] = tmp_lats[0];
      bound_box[n4 + 1] = tmp_lats[0];
      bound_box[n4 + 2] = tmp_lons[0];
      bound_box[n4 + 3] = tmp_lons[0];

      for (unsigned k = 1; k < 4; k++)
        {
          if (tmp_lats[k] < bound_box[n4]) bound_box[n4] = tmp_lats[k];
          if (tmp_lats[k] > bound_box[n4 + 1]) bound_box[n4 + 1] = tmp_lats[k];
          if (tmp_lons[k] < bound_box[n4 + 2]) bound_box[n4 + 2] = tmp_lons[k];
          if (tmp_lons[k] > bound_box[n4 + 3]) bound_box[n4 + 3] = tmp_lons[k];
        }
    }
}

void
remapgrid_get_lonlat(RemapGridType *grid, size_t cell_add, double *plon, double *plat)
{
  if (grid->remap_grid_type == REMAP_GRID_TYPE_REG2D)
    {
      size_t nx = grid->dims[0];
      size_t iy = cell_add / nx;
      size_t ix = cell_add - iy * nx;
      *plat = grid->reg2d_center_lat[iy];
      *plon = grid->reg2d_center_lon[ix];
      if (*plon < 0) *plon += PI2;
    }
  else
    {
      *plat = grid->cell_center_lat[cell_add];
      *plon = grid->cell_center_lon[cell_add];
    }
}

void
check_lon_range(size_t nlons, double *lons)
{
  assert(lons != NULL);

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(nlons, lons)
#endif
  for (size_t n = 0; n < nlons; ++n)
    {
      // remove missing values
      if (lons[n] < -PI2) lons[n] = 0;
      if (lons[n] > 2 * PI2) lons[n] = PI2;

      if (lons[n] > PI2) lons[n] -= PI2;
      if (lons[n] < ZERO) lons[n] += PI2;
    }
}

void
check_lat_range(size_t nlats, double *lats)
{
  assert(lats != NULL);

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(nlats, lats)
#endif
  for (size_t n = 0; n < nlats; ++n)
    {
      if (lats[n] > PIH) lats[n] = PIH;
      if (lats[n] < -PIH) lats[n] = -PIH;
    }
}

static void
check_lon_boundbox_range(size_t nlons, float *bound_box)
{
  size_t n4;

  assert(bound_box != NULL);

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(nlons, bound_box) private(n4)
#endif
  for (size_t n = 0; n < nlons; ++n)
    {
      n4 = n << 2;
      if (fabsf(bound_box[n4 + 3] - bound_box[n4 + 2]) > PI_f)
        {
          bound_box[n4 + 2] = 0.0f;
          bound_box[n4 + 3] = PI2_f;
        }
    }
}

static void
check_lat_boundbox_range(size_t nlats, float *restrict bound_box, double *restrict lats)
{
  size_t n4;

  assert(bound_box != NULL);

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(nlats, bound_box, lats) private(n4)
#endif
  for (size_t n = 0; n < nlats; ++n)
    {
      n4 = n << 2;
      if ((float)lats[n] < bound_box[n4]) bound_box[n4] = -PIH_f;
      if ((float)lats[n] > bound_box[n4 + 1]) bound_box[n4 + 1] = PIH_f;
    }
}

/*****************************************************************************/

static void
grid_check_lat_borders_rad(size_t n, double *ybounds)
{
  constexpr double YLIM = 88 * DEG2RAD;
  if (ybounds[0] > ybounds[n - 1])
    {
      if (ybounds[0] > YLIM) ybounds[0] = PIH;
      if (ybounds[n - 1] < -YLIM) ybounds[n - 1] = -PIH;
    }
  else
    {
      if (ybounds[0] < -YLIM) ybounds[0] = -PIH;
      if (ybounds[n - 1] > YLIM) ybounds[n - 1] = PIH;
    }
}

static void
remap_define_reg2d(int gridID, RemapGridType &grid)
{
  size_t nx = grid.dims[0];
  size_t ny = grid.dims[1];
  size_t nxp1 = nx + 1;
  size_t nyp1 = ny + 1;

  size_t nxm = nx;
  if (grid.is_cyclic) nxm++;

  if (grid.size != nx * ny) cdoAbort("Internal error, wrong dimensions!");

  grid.reg2d_center_lon = (double *) Malloc(nxm * sizeof(double));
  grid.reg2d_center_lat = (double *) Malloc(ny * sizeof(double));

  grid.reg2d_center_lon[0] = 0;
  grid.reg2d_center_lat[0] = 0;
  gridInqXvals(gridID, grid.reg2d_center_lon);
  gridInqYvals(gridID, grid.reg2d_center_lat);

  // Convert lat/lon units if required
  char yunits[CDI_MAX_NAME];
  yunits[0] = 0;
  cdiGridInqKeyStr(gridID, CDI_KEY_YUNITS, CDI_MAX_NAME, yunits);

  grid_to_radian(yunits, nx, grid.reg2d_center_lon, "grid reg2d center lon");
  grid_to_radian(yunits, ny, grid.reg2d_center_lat, "grid reg2d center lat");

  if (grid.reg2d_center_lon[nx - 1] < grid.reg2d_center_lon[0])
    for (size_t i = 1; i < nx; ++i)
      if (grid.reg2d_center_lon[i] < grid.reg2d_center_lon[i - 1]) grid.reg2d_center_lon[i] += PI2;

  if (grid.is_cyclic) grid.reg2d_center_lon[nx] = grid.reg2d_center_lon[0] + PI2;

  grid.reg2d_corner_lon = (double *) Malloc(nxp1 * sizeof(double));
  grid.reg2d_corner_lat = (double *) Malloc(nyp1 * sizeof(double));

  grid_gen_corners(nx, grid.reg2d_center_lon, grid.reg2d_corner_lon);
  grid_gen_corners(ny, grid.reg2d_center_lat, grid.reg2d_corner_lat);
  grid_check_lat_borders_rad(ny + 1, grid.reg2d_corner_lat);
}

static void
remapDefineGrid(RemapType mapType, int gridID, RemapGridType &grid, const char *txt)
{
  bool lgrid_destroy = false;
  bool lgrid_gen_bounds = false;
  int gridID_gme = -1;

  if (gridInqType(grid.gridID) != GRID_UNSTRUCTURED && gridInqType(grid.gridID) != GRID_CURVILINEAR)
    {
      if (gridInqType(grid.gridID) == GRID_GME)
        {
          gridID_gme = gridToUnstructured(grid.gridID, 1);
          grid.nvgp = gridInqSize(gridID_gme);
          gridID = gridDuplicate(gridID_gme);
          gridCompress(gridID);
          grid.luse_cell_corners = true;
        }
      else if (remap_write_remap || grid.remap_grid_type != REMAP_GRID_TYPE_REG2D)
        {
          lgrid_destroy = true;
          gridID = gridToCurvilinear(grid.gridID, 1);
          lgrid_gen_bounds = true;
        }
    }

  size_t gridsize = grid.size = gridInqSize(gridID);

  grid.dims[0] = gridInqXsize(gridID);
  grid.dims[1] = gridInqYsize(gridID);
  if (gridInqType(grid.gridID) != GRID_UNSTRUCTURED)
    {
      if (grid.dims[0] == 0)
        cdoAbort("%s grid without longitude coordinates!", gridNamePtr(gridInqType(grid.gridID)));
      if (grid.dims[1] == 0) cdoAbort("%s grid without latitude coordinates!", gridNamePtr(gridInqType(grid.gridID)));
    }

  grid.is_cyclic = (gridIsCircular(gridID) > 0);

  grid.rank = (gridInqType(gridID) == GRID_UNSTRUCTURED) ? 1 : 2;

  grid.num_cell_corners = (gridInqType(gridID) == GRID_UNSTRUCTURED) ? gridInqNvertex(gridID) : 4;

  remapGridAlloc(mapType, grid);

/* Initialize logical mask */

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(gridsize, grid)
#endif
  for (size_t i = 0; i < gridsize; ++i) grid.mask[i] = TRUE;

  if (gridInqMask(gridID, NULL))
    {
      std::vector<int> mask(gridsize);
      gridInqMask(gridID, &mask[0]);
      for (size_t i = 0; i < gridsize; ++i)
        if (mask[i] == 0) grid.mask[i] = FALSE;
    }

  if (!remap_write_remap && grid.remap_grid_type == REMAP_GRID_TYPE_REG2D) return;

  if (!(gridInqXvals(gridID, NULL) && gridInqYvals(gridID, NULL)))
    cdoAbort("%s grid cell center coordinates missing!", txt);

  gridInqXvals(gridID, grid.cell_center_lon);
  gridInqYvals(gridID, grid.cell_center_lat);

  char xunits[CDI_MAX_NAME];
  xunits[0] = 0;
  char yunits[CDI_MAX_NAME];
  yunits[0] = 0;
  cdiGridInqKeyStr(gridID, CDI_KEY_XUNITS, CDI_MAX_NAME, xunits);
  cdiGridInqKeyStr(gridID, CDI_KEY_YUNITS, CDI_MAX_NAME, yunits);

  if (grid.lneed_cell_corners)
    {
      if (gridInqXbounds(gridID, NULL) && gridInqYbounds(gridID, NULL))
        {
          gridInqXbounds(gridID, grid.cell_corner_lon);
          gridInqYbounds(gridID, grid.cell_corner_lat);
        }
      else if (lgrid_gen_bounds)
        {
          grid_cell_center_to_bounds_X2D(xunits, grid.dims[0], grid.dims[1], grid.cell_center_lon,
                                         grid.cell_corner_lon, 0);
          grid_cell_center_to_bounds_Y2D(yunits, grid.dims[0], grid.dims[1], grid.cell_center_lat,
                                         grid.cell_corner_lat);
        }
      else
        {
          cdoAbort("%s grid cell corner coordinates missing!", txt);
        }
    }

  if (gridInqType(grid.gridID) == GRID_GME) gridInqMaskGME(gridID_gme, grid.vgpm);

  /* Convert lat/lon units if required */

  grid_to_radian(xunits, grid.size, grid.cell_center_lon, "grid center lon");
  grid_to_radian(yunits, grid.size, grid.cell_center_lat, "grid center lat");
  /* Note: using units from cell center instead from bounds */
  if (grid.num_cell_corners && grid.lneed_cell_corners)
    {
      grid_to_radian(xunits, grid.num_cell_corners * grid.size, grid.cell_corner_lon, "grid corner lon");
      grid_to_radian(yunits, grid.num_cell_corners * grid.size, grid.cell_corner_lat, "grid corner lat");
    }

  if (lgrid_destroy) gridDestroy(gridID);

  /* Convert longitudes to 0,2pi interval */

  check_lon_range(grid.size, grid.cell_center_lon);

  if (grid.num_cell_corners && grid.lneed_cell_corners)
    check_lon_range(grid.num_cell_corners * grid.size, grid.cell_corner_lon);

  /*  Make sure input latitude range is within the machine values for +/- pi/2
   */

  check_lat_range(grid.size, grid.cell_center_lat);

  if (grid.num_cell_corners && grid.lneed_cell_corners)
    check_lat_range(grid.num_cell_corners * grid.size, grid.cell_corner_lat);
}

/*  Compute bounding boxes for restricting future grid searches */
static void
cell_bounding_boxes(RemapGridType &grid, float *cell_bound_box, int remap_grid_basis)
{
  if (grid.luse_cell_corners)
    {
      if (grid.lneed_cell_corners)
        {
          if (cdoVerbose) cdoPrint("Grid: boundbox_from_corners");
          boundbox_from_corners(grid.size, grid.num_cell_corners, grid.cell_corner_lon, grid.cell_corner_lat,
                                cell_bound_box);
        }
      else // full grid search
        {
          if (cdoVerbose) cdoPrint("Grid: bounds missing -> full grid search!");

          size_t gridsize = grid.size;
          for (size_t i = 0; i < gridsize; ++i)
            {
              cell_bound_box[i*4] = -PIH_f;
              cell_bound_box[i*4 + 1] = PIH_f;
              cell_bound_box[i*4 + 2] = 0.0f;
              cell_bound_box[i*4 + 3] = PI2_f;
            }
        }
    }
  else if (remap_grid_basis == REMAP_GRID_BASIS_SRC)
    {
      if (cdoVerbose) cdoPrint("Grid: boundbox_from_center");
      if (grid.rank != 2) cdoAbort("Internal problem, grid rank = %d!", grid.rank);

      size_t nx = grid.dims[0];
      size_t ny = grid.dims[1];
      boundbox_from_center(grid.is_cyclic, grid.size, nx, ny, grid.cell_center_lon, grid.cell_center_lat,
                           cell_bound_box);
    }

  if (remap_grid_basis == REMAP_GRID_BASIS_SRC || grid.lneed_cell_corners)
    check_lon_boundbox_range(grid.size, cell_bound_box);

  // Try to check for cells that overlap poles
  if (remap_grid_basis == REMAP_GRID_BASIS_SRC || grid.lneed_cell_corners)
    check_lat_boundbox_range(grid.size, cell_bound_box, grid.cell_center_lat);
}

void
remapGridInit(RemapGridType &grid)
{
  grid.remap_grid_type = -1;

  grid.num_cell_corners = 0;
  grid.luse_cell_corners = false;
  grid.lneed_cell_corners = false;

  grid.nvgp = 0;
  grid.vgpm = NULL;

  grid.mask = NULL;

  grid.reg2d_center_lon = NULL;
  grid.reg2d_center_lat = NULL;
  grid.reg2d_corner_lon = NULL;
  grid.reg2d_corner_lat = NULL;

  grid.cell_center_lon = NULL;
  grid.cell_center_lat = NULL;
  grid.cell_corner_lon = NULL;
  grid.cell_corner_lat = NULL;

  grid.cell_area = NULL;
  grid.cell_frac = NULL;
}

void
remapGridFree(RemapGridType &grid)
{
  if (grid.vgpm) Free(grid.vgpm);
  if (grid.mask) Free(grid.mask);

  if (grid.reg2d_center_lat) Free(grid.reg2d_center_lat);
  if (grid.reg2d_center_lon) Free(grid.reg2d_center_lon);
  if (grid.reg2d_corner_lat) Free(grid.reg2d_corner_lat);
  if (grid.reg2d_corner_lon) Free(grid.reg2d_corner_lon);

  if (grid.cell_center_lat) Free(grid.cell_center_lat);
  if (grid.cell_center_lon) Free(grid.cell_center_lon);
  if (grid.cell_corner_lat) Free(grid.cell_corner_lat);
  if (grid.cell_corner_lon) Free(grid.cell_corner_lon);

  if (grid.cell_area) Free(grid.cell_area);
  if (grid.cell_frac) Free(grid.cell_frac);
}

void
remapSearchInit(RemapType mapType, RemapSearch &search, RemapGridType &src_grid, RemapGridType &tgt_grid)
{
  search.src_bins.ncells = src_grid.size;
  search.tgt_bins.ncells = tgt_grid.size;

  search.src_bins.nbins = remap_num_srch_bins;
  search.tgt_bins.nbins = remap_num_srch_bins;

  search.src_bins.cell_bound_box = NULL;
  search.src_bins.bin_addr = NULL;
  search.src_bins.bin_lats = NULL;

  search.tgt_bins.cell_bound_box = NULL;
  search.tgt_bins.bin_addr = NULL;
  search.tgt_bins.bin_lats = NULL;

  if (!(src_grid.remap_grid_type == REMAP_GRID_TYPE_REG2D || tgt_grid.remap_grid_type == REMAP_GRID_TYPE_REG2D))
    {
      if (mapType != RemapType::DISTWGT
          //            && mapType != RemapType::BILINEAR
          )
        {
          search.src_bins.cell_bound_box = (float *) Malloc(4 * src_grid.size * sizeof(float));
          if ( tgt_grid.luse_cell_corners )
            search.tgt_bins.cell_bound_box = (float *) Malloc(4 * tgt_grid.size * sizeof(float));

          cell_bounding_boxes(src_grid, search.src_bins.cell_bound_box, REMAP_GRID_BASIS_SRC);
          cell_bounding_boxes(tgt_grid, search.tgt_bins.cell_bound_box, REMAP_GRID_BASIS_TGT);
          // Set up and assign address ranges to search bins in order to further restrict later searches
          calc_lat_bins(search.src_bins);
          if (mapType == RemapType::CONSERV || mapType == RemapType::CONSERV_YAC)
            {
              calc_lat_bins(search.tgt_bins);
              Free(search.src_bins.bin_lats);
              search.src_bins.bin_lats = NULL;
              Free(search.tgt_bins.bin_lats);
              search.tgt_bins.bin_lats = NULL;
              if (mapType == RemapType::CONSERV_YAC)
                {
                  Free(search.tgt_bins.cell_bound_box);
                  search.tgt_bins.cell_bound_box = NULL;
                }
            }
        }
    }
}

void
remapSearchFree(RemapSearch &search)
{
  if (search.src_bins.cell_bound_box) Free(search.src_bins.cell_bound_box);
  if (search.src_bins.bin_addr) Free(search.src_bins.bin_addr);
  if (search.src_bins.bin_lats) Free(search.src_bins.bin_lats);

  if (search.tgt_bins.cell_bound_box) Free(search.tgt_bins.cell_bound_box);
  if (search.tgt_bins.bin_addr) Free(search.tgt_bins.bin_addr);
  if (search.tgt_bins.bin_lats) Free(search.tgt_bins.bin_lats);
}

void
remapInitGrids(RemapType mapType, bool lextrapolate, int gridID1, RemapGridType &src_grid, int gridID2,
               RemapGridType &tgt_grid)
{
  int reg2d_src_gridID = gridID1;
  int reg2d_tgt_gridID = gridID2;

  remapGridInit(src_grid);
  remapGridInit(tgt_grid);

  if (mapType == RemapType::BILINEAR || mapType == RemapType::BICUBIC || mapType == RemapType::DISTWGT
      || mapType == RemapType::CONSERV_YAC)
    {
      if (IS_REG2D_GRID(gridID1)) src_grid.remap_grid_type = REMAP_GRID_TYPE_REG2D;
      // src_grid.remap_grid_type = 0;
    }

  if (src_grid.remap_grid_type == REMAP_GRID_TYPE_REG2D)
    {
      if (IS_REG2D_GRID(gridID2) && mapType == RemapType::CONSERV_YAC)
        tgt_grid.remap_grid_type = REMAP_GRID_TYPE_REG2D;
      // else src_grid.remap_grid_type = -1;
    }

  if (!remap_gen_weights && IS_REG2D_GRID(gridID2) && tgt_grid.remap_grid_type != REMAP_GRID_TYPE_REG2D)
    {
      if (mapType == RemapType::DISTWGT) tgt_grid.remap_grid_type = REMAP_GRID_TYPE_REG2D;
      if (mapType == RemapType::BILINEAR && src_grid.remap_grid_type == REMAP_GRID_TYPE_REG2D)
        tgt_grid.remap_grid_type = REMAP_GRID_TYPE_REG2D;
    }

  src_grid.lextrapolate = lextrapolate;

  if (mapType == RemapType::CONSERV || mapType == RemapType::CONSERV_YAC)
    {
      if (src_grid.remap_grid_type != REMAP_GRID_TYPE_REG2D)
        {
          src_grid.luse_cell_corners = true;
          src_grid.lneed_cell_corners = true;
        }

      if (tgt_grid.remap_grid_type != REMAP_GRID_TYPE_REG2D)
        {
          tgt_grid.luse_cell_corners = true;
          tgt_grid.lneed_cell_corners = true;
        }
    }

  src_grid.gridID = gridID1;
  tgt_grid.gridID = gridID2;

  if (gridInqType(gridID1) == GRID_UNSTRUCTURED)
    {
      if (gridInqYvals(gridID1, NULL) == 0 || gridInqXvals(gridID1, NULL) == 0)
        {
          if (gridInqNumber(gridID1) > 0)
            {
              src_grid.gridID = gridID1 = referenceToGrid(gridID1);
              if (gridID1 == -1) cdoAbort("Reference to source grid not found!");
            }
        }
    }

  if (gridInqType(gridID2) == GRID_UNSTRUCTURED)
    {
      if (gridInqYvals(gridID2, NULL) == 0 || gridInqXvals(gridID2, NULL) == 0)
        {
          if (gridInqNumber(gridID2) > 0)
            {
              tgt_grid.gridID = gridID2 = referenceToGrid(gridID2);
              if (gridID2 == -1) cdoAbort("Reference to target grid not found!");
            }
        }
    }

  int sgridID = src_grid.gridID;
  if (gridInqSize(sgridID) > 1
      && ((gridInqType(sgridID) == GRID_PROJECTION && gridInqProjType(sgridID) == CDI_PROJ_LCC)
          || (gridInqType(sgridID) == GRID_PROJECTION && gridInqProjType(sgridID) == CDI_PROJ_RLL)
          || (gridInqType(sgridID) == GRID_PROJECTION && gridInqProjType(sgridID) == CDI_PROJ_LAEA)
          || (gridInqType(sgridID) == GRID_PROJECTION && gridInqProjType(sgridID) == CDI_PROJ_SINU)))
    {
      int lbounds = TRUE;
      src_grid.gridID = gridID1 = gridToCurvilinear(src_grid.gridID, lbounds);
    }

  // if ( src_grid.remap_grid_type != REMAP_GRID_TYPE_REG2D )
  remapDefineGrid(mapType, gridID1, src_grid, "Source");
  remapDefineGrid(mapType, gridID2, tgt_grid, "Target");

  if (src_grid.remap_grid_type == REMAP_GRID_TYPE_REG2D || tgt_grid.remap_grid_type == REMAP_GRID_TYPE_REG2D)
    {
      if (src_grid.remap_grid_type == REMAP_GRID_TYPE_REG2D) remap_define_reg2d(reg2d_src_gridID, src_grid);
      if (tgt_grid.remap_grid_type == REMAP_GRID_TYPE_REG2D) remap_define_reg2d(reg2d_tgt_gridID, tgt_grid);
    }
}

/*****************************************************************************/

void
remapStat(int remapOrder, RemapGridType &src_grid, RemapGridType &tgt_grid, RemapVarsType &rv, const double *restrict array1,
          const double *restrict array2, double missval)
{
  if (remapOrder == 2)
    cdoPrint("Second order mapping from grid1 to grid2:");
  else
    cdoPrint("First order mapping from grid1 to grid2:");
  cdoPrint("----------------------------------------");

  double mean, minval, maxval;
  arrayMinMaxMeanMV(src_grid.size, array1, missval, &minval, &maxval, &mean);
  cdoPrint("Grid1 min,mean,max: %g %g %g", minval, mean, maxval);

  arrayMinMaxMeanMV(tgt_grid.size, array2, missval, &minval, &maxval, &mean);
  cdoPrint("Grid2 min,mean,max: %g %g %g", minval, mean, maxval);

  /* Conservation Test */

  if (src_grid.cell_area)
    {
      cdoPrint("Conservation:");
      double sum = 0;
      for (size_t n = 0; n < src_grid.size; ++n)
        if (!DBL_IS_EQUAL(array1[n], missval)) sum += array1[n] * src_grid.cell_area[n] * src_grid.cell_frac[n];
      cdoPrint("Grid1 Integral = %g", sum);

      sum = 0;
      for (size_t n = 0; n < tgt_grid.size; ++n)
        if (!DBL_IS_EQUAL(array2[n], missval)) sum += array2[n] * tgt_grid.cell_area[n] * tgt_grid.cell_frac[n];
      cdoPrint("Grid2 Integral = %g", sum);
      /*
      for ( n = 0; n < src_grid.size; n++ )
       fprintf(stderr, "1 %d %g %g %g\n", n, array1[n], src_grid.cell_area[n],
      src_grid.cell_frac[n]); for ( n = 0; n < tgt_grid.size; n++ )
        fprintf(stderr, "2 %d %g %g %g\n", n, array2[n], tgt_grid.cell_area[n],
      tgt_grid.cell_frac[n]);
      */
    }

  cdoPrint("Number of weights %zu", rv.num_wts);
  cdoPrint("Number of sparse matrix entries %zu", rv.num_links);
  cdoPrint("Total number of dest cells %zu", tgt_grid.size);

  std::vector<size_t> tgt_count(tgt_grid.size, 0);

#if defined(SX)
#pragma vdir nodep
#endif
  for (size_t n = 0; n < rv.num_links; ++n) tgt_count[rv.tgt_cell_add[n]]++;

  size_t imin = SIZE_MAX;
  size_t imax = 0;
  for (size_t n = 0; n < tgt_grid.size; ++n)
    {
      if (tgt_count[n] > 0)
        {
          if (tgt_count[n] < imin) imin = tgt_count[n];
          if (tgt_count[n] > imax) imax = tgt_count[n];
        }
    }

  size_t idiff = (imax - imin) / 10 + 1;
  size_t icount = 0;
  for (size_t i = 0; i < tgt_grid.size; ++i)
    if (tgt_count[i] > 0) icount++;

  cdoPrint("Number of cells participating in remap %zu", icount);

  if (icount)
    {
      cdoPrint("Min no of entries/row = %zu", imin);
      cdoPrint("Max no of entries/row = %zu", imax);

      imax = imin + idiff;
      for (size_t n = 0; n < 10; ++n)
        {
          icount = 0;
          for (size_t i = 0; i < tgt_grid.size; ++i)
            if (tgt_count[i] >= imin && tgt_count[i] < imax) icount++;

          if (icount) cdoPrint("Num of rows with entries between %zu - %zu  %zu", imin, imax - 1, icount);

          imin = imin + idiff;
          imax = imax + idiff;
        }
    }

  if (rv.sort_add) cdoPrint("Sparse matrix entries are explicitly sorted.");
}

/*****************************************************************************/

void
remapGradients(RemapGridType &grid, const double *restrict array, gradientsType &gradients)
{
  if (grid.rank != 2) cdoAbort("Internal problem (remapGradients), grid rank = %d!", grid.rank);

  double *restrict grad_lat = &gradients.grad_lat[0];
  double *restrict grad_lon = &gradients.grad_lon[0];
  double *restrict grad_latlon = &gradients.grad_latlon[0];
      
  size_t grid_size = grid.size;
  size_t nx = grid.dims[0];
  size_t ny = grid.dims[1];

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(grid_size, grad_lat, grad_lon, grad_latlon, grid, nx, ny, array)
#endif
  for (size_t n = 0; n < grid_size; ++n)
    {
      grad_lat[n] = ZERO;
      grad_lon[n] = ZERO;
      grad_latlon[n] = ZERO;

      if (grid.mask[n])
        {
          double delew = HALF;
          double delns = HALF;

          size_t j = n / nx + 1;
          size_t i = n - (j - 1) * nx + 1;

          size_t ip1 = i + 1;
          size_t im1 = i - 1;
          size_t jp1 = j + 1;
          size_t jm1 = j - 1;

          if (ip1 > nx) ip1 = ip1 - nx;
          if (im1 < 1) im1 = nx;
          if (jp1 > ny)
            {
              jp1 = j;
              delns = ONE;
            }
          if (jm1 < 1)
            {
              jm1 = j;
              delns = ONE;
            }

          size_t in = (jp1 - 1) * nx + i - 1;
          size_t is = (jm1 - 1) * nx + i - 1;
          size_t ie = (j - 1) * nx + ip1 - 1;
          size_t iw = (j - 1) * nx + im1 - 1;

          size_t ine = (jp1 - 1) * nx + ip1 - 1;
          size_t inw = (jp1 - 1) * nx + im1 - 1;
          size_t ise = (jm1 - 1) * nx + ip1 - 1;
          size_t isw = (jm1 - 1) * nx + im1 - 1;

          // Compute i-gradient
          if (!grid.mask[ie])
            {
              ie = n;
              delew = ONE;
            }
          if (!grid.mask[iw])
            {
              iw = n;
              delew = ONE;
            }

          grad_lat[n] = delew * (array[ie] - array[iw]);

          // Compute j-gradient
          if (!grid.mask[in])
            {
              in = n;
              delns = ONE;
            }
          if (!grid.mask[is])
            {
              is = n;
              delns = ONE;
            }

          grad_lon[n] = delns * (array[in] - array[is]);

          // Compute ij-gradient
          delew = HALF;
          delns = (jp1 == j || jm1 == j) ? ONE : HALF;

          if (!grid.mask[ine])
            {
              if (in != n)
                {
                  ine = in;
                  delew = ONE;
                }
              else if (ie != n)
                {
                  ine = ie;
                  inw = iw;
                  if (inw == n) delew = ONE;
                  delns = ONE;
                }
              else
                {
                  ine = n;
                  inw = iw;
                  delew = ONE;
                  delns = ONE;
                }
            }

          if (!grid.mask[inw])
            {
              if (in != n)
                {
                  inw = in;
                  delew = ONE;
                }
              else if (iw != n)
                {
                  inw = iw;
                  ine = ie;
                  if (ie == n) delew = ONE;
                  delns = ONE;
                }
              else
                {
                  inw = n;
                  ine = ie;
                  delew = ONE;
                  delns = ONE;
                }
            }

          double grad_lat_zero = delew * (array[ine] - array[inw]);

          if (!grid.mask[ise])
            {
              if (is != n)
                {
                  ise = is;
                  delew = ONE;
                }
              else if (ie != n)
                {
                  ise = ie;
                  isw = iw;
                  if (isw == n) delew = ONE;
                  delns = ONE;
                }
              else
                {
                  ise = n;
                  isw = iw;
                  delew = ONE;
                  delns = ONE;
                }
            }

          if (!grid.mask[isw])
            {
              if (is != n)
                {
                  isw = is;
                  delew = ONE;
                }
              else if (iw != n)
                {
                  isw = iw;
                  ise = ie;
                  if (ie == n) delew = ONE;
                  delns = ONE;
                }
              else
                {
                  isw = n;
                  ise = ie;
                  delew = ONE;
                  delns = ONE;
                }
            }

          double grad_lon_zero = delew * (array[ise] - array[isw]);

          grad_latlon[n] = delns * (grad_lat_zero - grad_lon_zero);
        }
    }
} /* remapGradients */

/*****************************************************************************/

void
remapCheckArea(size_t grid_size, double *restrict cell_area, const char *name)
{
  for (size_t n = 0; n < grid_size; ++n)
    {
      if (cell_area[n] < -.01) cdoPrint("%s grid area error: %zu %g", name, n, cell_area[n]);
    }
}

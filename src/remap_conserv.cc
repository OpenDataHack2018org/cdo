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

#include "cdo_int.h"
#include "cdo_wtime.h"
#include "grid.h"
#include "remap.h"
#include "remap_store_link.h"
#include "cdoOptions.h"

extern "C" {
#include "lib/yac/clipping.h"
#include "lib/yac/area.h"
#include "lib/yac/geometry.h"
#include "lib/yac/grid_search.h"
}

struct search_t
{
  enum yac_edge_type *src_edge_type;
  size_t srch_corners;
  size_t max_srch_cells;
  std::vector<double> partial_areas;
  std::vector<double> partial_weights;
  std::vector<struct grid_cell> src_grid_cells;
  std::vector<struct grid_cell> overlap_buffer;
};

static void
cellSearchRealloc(size_t num_srch_cells, search_t &search)
{
  if (num_srch_cells > search.max_srch_cells)
    {
      search.partial_areas.resize(num_srch_cells);
      search.partial_weights.resize(num_srch_cells);
      search.overlap_buffer.resize(num_srch_cells);
      search.src_grid_cells.resize(num_srch_cells);

      for (size_t n = search.max_srch_cells; n < num_srch_cells; ++n)
        {
          search.overlap_buffer[n].array_size = 0;
          search.overlap_buffer[n].num_corners = 0;
          search.overlap_buffer[n].edge_type = NULL;
          search.overlap_buffer[n].coordinates_x = NULL;
          search.overlap_buffer[n].coordinates_y = NULL;
          search.overlap_buffer[n].coordinates_xyz = NULL;

          search.src_grid_cells[n].array_size = search.srch_corners;
          search.src_grid_cells[n].num_corners = search.srch_corners;
          search.src_grid_cells[n].edge_type = search.src_edge_type;
          search.src_grid_cells[n].coordinates_x = new double[search.srch_corners];
          search.src_grid_cells[n].coordinates_y = new double[search.srch_corners];
          search.src_grid_cells[n].coordinates_xyz = new double[3 * search.srch_corners];
        }

      search.max_srch_cells = num_srch_cells;
    }
}

static void
cellSearchFree(search_t &search)
{
  size_t max_srch_cells = search.max_srch_cells;
  for (size_t n = 0; n < max_srch_cells; n++)
    {
      if (search.overlap_buffer[n].array_size > 0)
        {
          Free(search.overlap_buffer[n].coordinates_x);
          Free(search.overlap_buffer[n].coordinates_y);
          if (search.overlap_buffer[n].coordinates_xyz) Free(search.overlap_buffer[n].coordinates_xyz);
          if (search.overlap_buffer[n].edge_type) Free(search.overlap_buffer[n].edge_type);
        }

      delete[] search.src_grid_cells[n].coordinates_x;
      delete[] search.src_grid_cells[n].coordinates_y;
      delete[] search.src_grid_cells[n].coordinates_xyz;
    }

  vectorFree(search.partial_areas);
  vectorFree(search.partial_weights);
  vectorFree(search.overlap_buffer);
  vectorFree(search.src_grid_cells);
}

int rect_grid_search2(long *imin, long *imax, double xmin, double xmax, long nxm, const double *restrict xm);

static size_t
get_srch_cells_reg2d(const size_t *restrict src_grid_dims, const double *restrict src_corner_lat,
                     const double *restrict src_corner_lon, const double *restrict tgt_cell_bound_box, size_t *srch_add)
{
  bool debug = false;
  long nx = src_grid_dims[0];
  long ny = src_grid_dims[1];
  size_t num_srch_cells = 0;  // num cells in restricted search arrays

  long nxp1 = nx + 1;
  long nyp1 = ny + 1;

  double src_lon_min = src_corner_lon[0];
  double src_lon_max = src_corner_lon[nx];

  long imin = nxp1, imax = -1, jmin = nyp1, jmax = -1;

  int lfound = rect_grid_search2(&jmin, &jmax, tgt_cell_bound_box[0], tgt_cell_bound_box[1], nyp1, src_corner_lat);
  if (!lfound) return 0;
  // printf("lfound, jmin, jmax %d %ld %ld\n", lfound, jmin, jmax);
  // if ( jmin > 0 ) jmin--;
  // if ( jmax < (ny-2) ) jmax++;

  double bound_lon1 = tgt_cell_bound_box[2];
  double bound_lon2 = tgt_cell_bound_box[3];
  if (bound_lon1 <= src_lon_max && bound_lon2 >= src_lon_min)
    {
      if (debug) printf("  b1 %g %g\n", bound_lon1 * RAD2DEG, bound_lon2 * RAD2DEG);
      if (bound_lon1 < src_lon_min && bound_lon2 > src_lon_min) bound_lon1 = src_lon_min;
      if (bound_lon2 > src_lon_max && bound_lon1 < src_lon_max) bound_lon2 = src_lon_max;
      lfound = rect_grid_search2(&imin, &imax, bound_lon1, bound_lon2, nxp1, src_corner_lon);
      if (lfound)
        {
          if (debug)
            printf("   %g %g imin %ld  imax %ld  jmin %ld jmax %ld\n", RAD2DEG * src_corner_lon[imin],
                   RAD2DEG * src_corner_lon[imax + 1], imin, imax, jmin, jmax);
          for (long jm = jmin; jm <= jmax; ++jm)
            for (long im = imin; im <= imax; ++im) srch_add[num_srch_cells++] = jm * nx + im;
        }
    }

  bound_lon1 = tgt_cell_bound_box[2];
  bound_lon2 = tgt_cell_bound_box[3];
  if (bound_lon1 <= src_lon_min && bound_lon2 >= src_lon_min)
    {
      long imin2 = nxp1, imax2 = -1;
      bound_lon1 += 2 * M_PI;
      bound_lon2 += 2 * M_PI;
      if (debug) printf("  b2 %g %g\n", bound_lon1 * RAD2DEG, bound_lon2 * RAD2DEG);
      if (bound_lon1 < src_lon_min && bound_lon2 > src_lon_min) bound_lon1 = src_lon_min;
      if (bound_lon2 > src_lon_max && bound_lon1 < src_lon_max) bound_lon2 = src_lon_max;
      lfound = rect_grid_search2(&imin2, &imax2, bound_lon1, bound_lon2, nxp1, src_corner_lon);
      if (lfound)
        {
          if (imax != -1 && imin2 <= imax) imin2 = imax + 1;
          if (imax != -1 && imax2 <= imax) imax2 = imax + 1;
          if (imin2 >= 0 && imax2 < nxp1)
            {
              if (debug)
                printf("   %g %g imin %ld  imax %ld  jmin %ld jmax %ld\n", RAD2DEG * src_corner_lon[imin2],
                       RAD2DEG * src_corner_lon[imax2 + 1], imin2, imax2, jmin, jmax);
              for (long jm = jmin; jm <= jmax; ++jm)
                for (long im = imin2; im <= imax2; ++im) srch_add[num_srch_cells++] = jm * nx + im;
            }
        }
    }

  bound_lon1 = tgt_cell_bound_box[2];
  bound_lon2 = tgt_cell_bound_box[3];
  if (bound_lon1 <= src_lon_max && bound_lon2 >= src_lon_max)
    {
      long imin3 = nxp1, imax3 = -1;
      bound_lon1 -= 2 * M_PI;
      bound_lon2 -= 2 * M_PI;
      if (debug) printf("  b3 %g %g\n", bound_lon1 * RAD2DEG, bound_lon2 * RAD2DEG);
      if (bound_lon1 < src_lon_min && bound_lon2 > src_lon_min) bound_lon1 = src_lon_min;
      if (bound_lon2 > src_lon_max && bound_lon1 < src_lon_max) bound_lon2 = src_lon_max;
      lfound = rect_grid_search2(&imin3, &imax3, bound_lon1, bound_lon2, nxp1, src_corner_lon);
      if (lfound)
        {
          if (imin != nxp1 && imin3 >= imin) imin3 = imin - 1;
          if (imax != nxp1 && imax3 >= imin) imax3 = imin - 1;
          if (imin3 >= 0 && imin3 < nxp1)
            {
              if (debug)
                printf("   %g %g imin %ld  imax %ld  jmin %ld jmax %ld\n", RAD2DEG * src_corner_lon[imin3],
                       RAD2DEG * src_corner_lon[imax3 + 1], imin3, imax3, jmin, jmax);
              for (long jm = jmin; jm <= jmax; ++jm)
                for (long im = imin3; im <= imax3; ++im) srch_add[num_srch_cells++] = jm * nx + im;
            }
        }
    }

  if (debug) printf(" . num_srch_cells: %zu\n", num_srch_cells);

  return num_srch_cells;
}

static void
restrict_boundbox(const double *restrict grid_bound_box, double *restrict bound_box)
{
  if (bound_box[0] < grid_bound_box[0] && bound_box[1] > grid_bound_box[0]) bound_box[0] = grid_bound_box[0];
  if (bound_box[1] > grid_bound_box[1] && bound_box[0] < grid_bound_box[1]) bound_box[1] = grid_bound_box[1];

  if (bound_box[2] >= grid_bound_box[3] && (bound_box[3] - 2 * M_PI) > grid_bound_box[2])
    {
      bound_box[2] -= 2 * M_PI;
      bound_box[3] -= 2 * M_PI;
    }
  if (bound_box[3] <= grid_bound_box[2] && (bound_box[2] - 2 * M_PI) < grid_bound_box[3])
    {
      bound_box[2] += 2 * M_PI;
      bound_box[3] += 2 * M_PI;
    }
  //  if ( bound_box[2] < grid_bound_box[2] && bound_box[3] > grid_bound_box[2]
  //  ) bound_box[2] = grid_bound_box[2]; if ( bound_box[3] > grid_bound_box[3]
  //  && bound_box[2] < grid_bound_box[3] ) bound_box[3] = grid_bound_box[3];
}

static void
boundbox_from_corners_reg2d(size_t grid_add, const size_t *restrict grid_dims, const double *restrict corner_lon,
                            const double *restrict corner_lat, double *restrict bound_box)
{
  size_t nx = grid_dims[0];
  size_t iy = grid_add / nx;
  size_t ix = grid_add - iy * nx;

  double clat1 = corner_lat[iy];
  double clat2 = corner_lat[iy + 1];

  if (clat2 > clat1)
    {
      bound_box[0] = clat1;
      bound_box[1] = clat2;
    }
  else
    {
      bound_box[0] = clat2;
      bound_box[1] = clat1;
    }

  bound_box[2] = corner_lon[ix];
  bound_box[3] = corner_lon[ix + 1];
}

static void
boundbox_from_corners1(size_t ic, size_t nc, const double *restrict corner_lon, const double *restrict corner_lat,
                       double *restrict bound_box)
{
  size_t inc = ic * nc;

  double clat = corner_lat[inc];
  double clon = corner_lon[inc];

  bound_box[0] = clat;
  bound_box[1] = clat;
  bound_box[2] = clon;
  bound_box[3] = clon;

  for (size_t j = 1; j < nc; ++j)
    {
      clat = corner_lat[inc + j];
      clon = corner_lon[inc + j];

      if (clat < bound_box[0]) bound_box[0] = clat;
      if (clat > bound_box[1]) bound_box[1] = clat;
      if (clon < bound_box[2]) bound_box[2] = clon;
      if (clon > bound_box[3]) bound_box[3] = clon;
    }

  if (fabs(bound_box[3] - bound_box[2]) > PI)
    {
      bound_box[2] = 0;
      bound_box[3] = PI2;
    }

  /*
  double dlon = fabsf(bound_box[3] - bound_box[2]);

  if ( dlon > PI )
    {
      if ( bound_box[3] > bound_box[2] && (bound_box[3]-PI2) < 0. )
        {
          double tmp = bound_box[2];
          bound_box[2] = bound_box[3] - PI2;
          bound_box[3] = tmp;
        }
      else
        {
          bound_box[2] = 0;
          bound_box[3] = PI2;
        }
    }
  */
}

static void
boundbox_from_corners1r(size_t ic, size_t nc, const double *restrict corner_lon, const double *restrict corner_lat,
                        float *restrict bound_box)
{
  size_t inc = ic * nc;

  float clat = corner_lat[inc];
  float clon = corner_lon[inc];

  bound_box[0] = clat;
  bound_box[1] = clat;
  bound_box[2] = clon;
  bound_box[3] = clon;

  for (size_t j = 1; j < nc; ++j)
    {
      clat = corner_lat[inc + j];
      clon = corner_lon[inc + j];

      if (clat < bound_box[0]) bound_box[0] = clat;
      if (clat > bound_box[1]) bound_box[1] = clat;
      if (clon < bound_box[2]) bound_box[2] = clon;
      if (clon > bound_box[3]) bound_box[3] = clon;
    }

  if (fabsf(bound_box[3] - bound_box[2]) > PI_f)
    {
      bound_box[2] = 0;
      bound_box[3] = PI2_f;
    }
  /*
  if ( fabsf(bound_box[3] - bound_box[2]) > PI_f )
    {
      if ( bound_box[3] > bound_box[2] && (bound_box[3]-PI2_f) <
  0.0f )
        {
          float tmp = bound_box[2];
          bound_box[2] = bound_box[3] - PI2_f;
          bound_box[3] = tmp;
        }
    }
  */
}

static double
gridcell_area(struct grid_cell cell)
{
  return yac_huiliers_area(cell);
}

static void
cdo_compute_overlap_areas(size_t N, search_t *search, struct grid_cell target_cell)
{
  struct grid_cell *overlap_buffer = &search->overlap_buffer[0];
  struct grid_cell *source_cells = &search->src_grid_cells[0];

  // Do the clipping and get the cell for the overlapping area
  yac_cell_clipping(N, source_cells, target_cell, overlap_buffer);

  // Get the partial areas for the overlapping regions
  for (size_t n = 0; n < N; n++)
    {
      search->partial_areas[n] = gridcell_area(overlap_buffer[n]);
    }

#ifdef VERBOSE
  for (size_t n = 0; n < N; n++) printf("overlap area : %lf\n", search->partial_areas[n]);
#endif
}

static double const tol = 1.0e-12;

enum cell_type
{
  UNDEF_CELL,
  LON_LAT_CELL,
  LAT_CELL,
  GREAT_CIRCLE_CELL,
  MIXED_CELL
};
/*
static enum cell_type get_cell_type(struct grid_cell target_cell) {

  int count_lat_edges = 0, count_great_circle_edges = 0;

   if ((target_cell.num_corners == 4) &&
       ((target_cell.edge_type[0] == LAT_CIRCLE &&
         target_cell.edge_type[1] == LON_CIRCLE &&
         target_cell.edge_type[2] == LAT_CIRCLE &&
         target_cell.edge_type[3] == LON_CIRCLE) ||
        (target_cell.edge_type[0] == LON_CIRCLE &&
         target_cell.edge_type[1] == LAT_CIRCLE &&
         target_cell.edge_type[2] == LON_CIRCLE &&
         target_cell.edge_type[3] == LAT_CIRCLE)))
      return LON_LAT_CELL;
   else
      for (unsigned i = 0; i < target_cell.num_corners; ++i)
         if (target_cell.edge_type[i] == LON_CIRCLE ||
             target_cell.edge_type[i] == GREAT_CIRCLE)
            count_great_circle_edges++;
         else
            count_lat_edges++;

   if (count_lat_edges && count_great_circle_edges)
      return MIXED_CELL;
   else if (count_lat_edges)
      return LAT_CELL;
   else
      return GREAT_CIRCLE_CELL;
}
*/
static void
cdo_compute_concave_overlap_areas(size_t N, search_t *search, struct grid_cell target_cell, double target_node_x,
                                  double target_node_y)
{
  double *partial_areas = &search->partial_areas[0];
  struct grid_cell *overlap_buffer = &search->overlap_buffer[0];
  struct grid_cell *source_cell = &search->src_grid_cells[0];

  /*
  enum cell_type target_cell_type = UNDEF_CELL;

  if ( target_cell.num_corners > 3 )
    target_cell_type = get_cell_type(target_cell);

  if ( target_cell.num_corners < 4 || target_cell_type == LON_LAT_CELL )
    {
      cdo_compute_overlap_areas(N, overlap_buffer, source_cell, target_cell,
  partial_areas); return;
    }

  if ( target_node_x == NULL || target_node_y == NULL )
    cdoAbort("Internal problem (cdo_compute_concave_overlap_areas): missing
  target point coordinates!");
  */
  /*
  struct grid_cell target_partial_cell =
    {.coordinates_x   = (double[3]){-1, -1, -1},
     .coordinates_y   = (double[3]){-1, -1, -1},
     .coordinates_xyz = (double[3*3]){-1, -1, -1},
     .edge_type       = (enum yac_edge_type[3]) {GREAT_CIRCLE, GREAT_CIRCLE,
  GREAT_CIRCLE}, .num_corners     = 3};
  */
  double coordinates_x[3] = { -1, -1, -1 };
  double coordinates_y[3] = { -1, -1, -1 };
  double coordinates_xyz[9] = { -1, -1, -1 };
  enum yac_edge_type edge_types[3] = { GREAT_CIRCLE, GREAT_CIRCLE, GREAT_CIRCLE };
  struct grid_cell target_partial_cell;
  target_partial_cell.coordinates_x = coordinates_x;
  target_partial_cell.coordinates_y = coordinates_y;
  target_partial_cell.coordinates_xyz = coordinates_xyz;
  target_partial_cell.edge_type = edge_types;
  target_partial_cell.num_corners = 3;

  /* Do the clipping and get the cell for the overlapping area */

  for (size_t n = 0; n < N; n++) partial_areas[n] = 0.0;

  // common node point to all partial target cells
  target_partial_cell.coordinates_x[0] = target_node_x;
  target_partial_cell.coordinates_y[0] = target_node_y;

  LLtoXYZ(target_node_x, target_node_y, target_partial_cell.coordinates_xyz);

  for (unsigned num_corners = 0; num_corners < target_cell.num_corners; ++num_corners)
    {
      unsigned corner_a = num_corners;
      unsigned corner_b = (num_corners + 1) % target_cell.num_corners;

      // skip clipping and area calculation for degenerated triangles
      //
      // If this is not sufficient, instead we can try something like:
      //
      //     struct point_list target_list
      //     init_point_list(&target_list);
      //     generate_point_list(&target_list, target_cell);
      //     struct grid_cell temp_target_cell;
      //     generate_overlap_cell(target_list, temp_target_cell);
      //     free_point_list(&target_list);
      //
      // and use temp_target_cell for triangulation.
      //
      // Compared to the if statement below the alternative seems
      // to be quite costly.

      if (((fabs(target_cell.coordinates_xyz[0 + 3 * corner_a] - target_cell.coordinates_xyz[0 + 3 * corner_b]) < tol)
           && (fabs(target_cell.coordinates_xyz[1 + 3 * corner_a] - target_cell.coordinates_xyz[1 + 3 * corner_b]) < tol)
           && (fabs(target_cell.coordinates_xyz[2 + 3 * corner_a] - target_cell.coordinates_xyz[2 + 3 * corner_b]) < tol))
          || ((fabs(target_cell.coordinates_xyz[0 + 3 * corner_a] - target_partial_cell.coordinates_xyz[0]) < tol)
              && (fabs(target_cell.coordinates_xyz[1 + 3 * corner_a] - target_partial_cell.coordinates_xyz[1]) < tol)
              && (fabs(target_cell.coordinates_xyz[2 + 3 * corner_a] - target_partial_cell.coordinates_xyz[2]) < tol))
          || ((fabs(target_cell.coordinates_xyz[0 + 3 * corner_b] - target_partial_cell.coordinates_xyz[0]) < tol)
              && (fabs(target_cell.coordinates_xyz[1 + 3 * corner_b] - target_partial_cell.coordinates_xyz[1]) < tol)
              && (fabs(target_cell.coordinates_xyz[2 + 3 * corner_b] - target_partial_cell.coordinates_xyz[2]) < tol)))
        continue;

      target_partial_cell.coordinates_x[1] = target_cell.coordinates_x[corner_a];
      target_partial_cell.coordinates_y[1] = target_cell.coordinates_y[corner_a];
      target_partial_cell.coordinates_x[2] = target_cell.coordinates_x[corner_b];
      target_partial_cell.coordinates_y[2] = target_cell.coordinates_y[corner_b];

      target_partial_cell.coordinates_xyz[0 + 3 * 1] = target_cell.coordinates_xyz[0 + 3 * corner_a];
      target_partial_cell.coordinates_xyz[1 + 3 * 1] = target_cell.coordinates_xyz[1 + 3 * corner_a];
      target_partial_cell.coordinates_xyz[2 + 3 * 1] = target_cell.coordinates_xyz[2 + 3 * corner_a];
      target_partial_cell.coordinates_xyz[0 + 3 * 2] = target_cell.coordinates_xyz[0 + 3 * corner_b];
      target_partial_cell.coordinates_xyz[1 + 3 * 2] = target_cell.coordinates_xyz[1 + 3 * corner_b];
      target_partial_cell.coordinates_xyz[2 + 3 * 2] = target_cell.coordinates_xyz[2 + 3 * corner_b];

      yac_cell_clipping((unsigned) N, source_cell, target_partial_cell, overlap_buffer);

      // Get the partial areas for the overlapping regions as sum over the partial target cells.
      for (size_t n = 0; n < N; n++)
        {
          partial_areas[n] += gridcell_area(overlap_buffer[n]);
          // we cannot use pole_area because it is rather inaccurate for great
          // circle edges that are nearly circles of longitude
          // partial_areas[n] = pole_area (overlap_buffer[n]);
        }
    }

#ifdef VERBOSE
  for (size_t n = 0; n < N; n++) printf("overlap area %zu: %lf\n", n, partial_areas[n]);
#endif
}

//#endif

static int
get_lonlat_circle_index(RemapGrid *remap_grid)
{
  int lonlat_circle_index = -1;

  if (remap_grid->num_cell_corners == 4)
    {
      if (remap_grid->remap_grid_type == REMAP_GRID_TYPE_REG2D)
        {
          lonlat_circle_index = 1;
        }
      else
        {
          const double *cell_corner_lon = remap_grid->cell_corner_lon;
          const double *cell_corner_lat = remap_grid->cell_corner_lat;
          size_t gridsize = remap_grid->size;
          size_t num_i = 0, num_eq0 = 0, num_eq1 = 0;
          long iadd = gridsize / 3 - 1;

          if (iadd == 0) iadd++;

          for (size_t i = 0; i < gridsize; i += iadd)
            {
              num_i++;

              if (IS_EQUAL(cell_corner_lon[i * 4 + 1], cell_corner_lon[i * 4 + 2])
                  && IS_EQUAL(cell_corner_lon[i * 4 + 3], cell_corner_lon[i * 4 + 0])
                  && IS_EQUAL(cell_corner_lat[i * 4 + 0], cell_corner_lat[i * 4 + 1])
                  && IS_EQUAL(cell_corner_lat[i * 4 + 2], cell_corner_lat[i * 4 + 3]))
                {
                  num_eq1++;
                }
              else if (IS_EQUAL(cell_corner_lon[i * 4 + 0], cell_corner_lon[i * 4 + 1])
                       && IS_EQUAL(cell_corner_lon[i * 4 + 2], cell_corner_lon[i * 4 + 3])
                       && IS_EQUAL(cell_corner_lat[i * 4 + 1], cell_corner_lat[i * 4 + 2])
                       && IS_EQUAL(cell_corner_lat[i * 4 + 3], cell_corner_lat[i * 4 + 0]))
                {
                  num_eq0++;
                }
            }

          if (num_i == num_eq1) lonlat_circle_index = 1;
          if (num_i == num_eq0) lonlat_circle_index = 0;
        }
    }

  // printf("lonlat_circle_index %d\n", lonlat_circle_index);

  return lonlat_circle_index;
}

static void
remapNormalizeWeights(RemapGrid *tgt_grid, RemapVars &rv)
{
  // Include centroids in weights and normalize using destination area if requested
  size_t num_links = rv.num_links;
  size_t num_wts = rv.num_wts;
  size_t tgt_cell_add;     // current linear address for target grid cell
  double norm_factor = 0;  // factor for normalizing wts

  if (rv.normOpt == NormOpt::DESTAREA)
    {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(num_wts, num_links, rv, tgt_grid) private(tgt_cell_add, norm_factor)
#endif
      for (size_t n = 0; n < num_links; ++n)
        {
          tgt_cell_add = rv.tgt_cell_add[n];
          norm_factor = IS_NOT_EQUAL(tgt_grid->cell_area[tgt_cell_add], 0) ? 1. / tgt_grid->cell_area[tgt_cell_add] : 0.;
          rv.wts[n * num_wts] *= norm_factor;
        }
    }
  else if (rv.normOpt == NormOpt::FRACAREA)
    {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(num_wts, num_links, rv, tgt_grid) private(tgt_cell_add, norm_factor)
#endif
      for (size_t n = 0; n < num_links; ++n)
        {
          tgt_cell_add = rv.tgt_cell_add[n];
          norm_factor = IS_NOT_EQUAL(tgt_grid->cell_frac[tgt_cell_add], 0) ? 1. / tgt_grid->cell_frac[tgt_cell_add] : 0.;
          rv.wts[n * num_wts] *= norm_factor;
        }
    }
  else if (rv.normOpt == NormOpt::NONE)
    {
    }
}

static void
set_yac_coordinates(int remap_grid_type, size_t cell_add, size_t num_cell_corners, RemapGrid *remap_grid,
                    struct grid_cell *yac_grid_cell)
{
  double *coordinates_x = yac_grid_cell->coordinates_x;
  double *coordinates_y = yac_grid_cell->coordinates_y;

  if (remap_grid_type == REMAP_GRID_TYPE_REG2D)
    {
      size_t nx = remap_grid->dims[0];
      size_t iy = cell_add / nx;
      size_t ix = cell_add - iy * nx;
      const double *reg2d_corner_lon = remap_grid->reg2d_corner_lon;
      const double *reg2d_corner_lat = remap_grid->reg2d_corner_lat;
      coordinates_x[0] = reg2d_corner_lon[ix];
      coordinates_y[0] = reg2d_corner_lat[iy];
      coordinates_x[1] = reg2d_corner_lon[ix + 1];
      coordinates_y[1] = reg2d_corner_lat[iy];
      coordinates_x[2] = reg2d_corner_lon[ix + 1];
      coordinates_y[2] = reg2d_corner_lat[iy + 1];
      coordinates_x[3] = reg2d_corner_lon[ix];
      coordinates_y[3] = reg2d_corner_lat[iy + 1];
    }
  else
    {
      const double *cell_corner_lon = remap_grid->cell_corner_lon;
      const double *cell_corner_lat = remap_grid->cell_corner_lat;
      for (size_t ic = 0; ic < num_cell_corners; ++ic)
        {
          coordinates_x[ic] = cell_corner_lon[cell_add * num_cell_corners + ic];
          coordinates_y[ic] = cell_corner_lat[cell_add * num_cell_corners + ic];
        }
    }

  double *coordinates_xyz = yac_grid_cell->coordinates_xyz;
  for (size_t ic = 0; ic < num_cell_corners; ++ic) LLtoXYZ(coordinates_x[ic], coordinates_y[ic], coordinates_xyz + ic * 3);
}

static void
reg2d_bound_box(RemapGrid *remap_grid, double *grid_bound_box)
{
  size_t nx = remap_grid->dims[0];
  size_t ny = remap_grid->dims[1];
  const double *reg2d_corner_lon = remap_grid->reg2d_corner_lon;
  const double *reg2d_corner_lat = remap_grid->reg2d_corner_lat;

  grid_bound_box[0] = reg2d_corner_lat[0];
  grid_bound_box[1] = reg2d_corner_lat[ny];
  if (grid_bound_box[0] > grid_bound_box[1])
    {
      grid_bound_box[0] = reg2d_corner_lat[ny];
      grid_bound_box[1] = reg2d_corner_lat[0];
    }
  grid_bound_box[2] = reg2d_corner_lon[0];
  grid_bound_box[3] = reg2d_corner_lon[nx];
}

void
remapConservWeights(RemapSearch &rsearch, RemapVars &rv)
{
  RemapGrid *src_grid = rsearch.srcGrid;
  RemapGrid *tgt_grid = rsearch.tgtGrid;

  bool lcheck = true;

  // Variables necessary if segment manages to hit pole
  int src_remap_grid_type = src_grid->remap_grid_type;
  int tgt_remap_grid_type = tgt_grid->remap_grid_type;

  if (cdoVerbose) cdoPrint("Called %s()", __func__);

  progressInit();

  double start = cdoVerbose ? cdo_get_wtime() : 0;

  size_t src_grid_size = src_grid->size;
  size_t tgt_grid_size = tgt_grid->size;

  size_t src_num_cell_corners = src_grid->num_cell_corners;
  size_t tgt_num_cell_corners = tgt_grid->num_cell_corners;

  size_t max_num_cell_corners = src_num_cell_corners;
  if (tgt_num_cell_corners > max_num_cell_corners) max_num_cell_corners = tgt_num_cell_corners;

  std::vector<enum yac_edge_type> great_circle_type(max_num_cell_corners);
  for (size_t i = 0; i < max_num_cell_corners; ++i) great_circle_type[i] = GREAT_CIRCLE;

  enum yac_edge_type lonlat_circle_type[] = { LON_CIRCLE, LAT_CIRCLE, LON_CIRCLE, LAT_CIRCLE, LON_CIRCLE };

  enum yac_edge_type *src_edge_type = &great_circle_type[0];
  enum yac_edge_type *tgt_edge_type = &great_circle_type[0];

  enum cell_type target_cell_type = UNDEF_CELL;

  if (src_num_cell_corners == 4)
    {
      int lonlat_circle_index = get_lonlat_circle_index(src_grid);
      if (lonlat_circle_index >= 0) src_edge_type = &lonlat_circle_type[lonlat_circle_index];
    }

  if (tgt_num_cell_corners == 4)
    {
      int lonlat_circle_index = get_lonlat_circle_index(tgt_grid);
      if (lonlat_circle_index >= 0)
        {
          target_cell_type = LON_LAT_CELL;
          tgt_edge_type = &lonlat_circle_type[lonlat_circle_index];
        }
    }

  if (!(tgt_num_cell_corners < 4 || target_cell_type == LON_LAT_CELL))
    {
      if (tgt_grid->cell_center_lon == NULL || tgt_grid->cell_center_lat == NULL)
        cdoAbort("Internal problem (%s): missing target point coordinates!", __func__);
    }

  std::vector<struct grid_cell> tgt_grid_cell(Threading::ompNumThreads);
  for (int i = 0; i < Threading::ompNumThreads; ++i)
    {
      tgt_grid_cell[i].array_size = tgt_num_cell_corners;
      tgt_grid_cell[i].num_corners = tgt_num_cell_corners;
      tgt_grid_cell[i].edge_type = tgt_edge_type;
      tgt_grid_cell[i].coordinates_x = new double[tgt_num_cell_corners];
      tgt_grid_cell[i].coordinates_y = new double[tgt_num_cell_corners];
      tgt_grid_cell[i].coordinates_xyz = new double[3 * tgt_num_cell_corners];
    }

  std::vector<search_t> search(Threading::ompNumThreads);
  for (int i = 0; i < Threading::ompNumThreads; ++i)
    {
      search[i].srch_corners = src_num_cell_corners;
      search[i].src_edge_type = src_edge_type;
      search[i].max_srch_cells = 0;
    }

  std::vector<size_t *> srch_add(Threading::ompNumThreads);
  for (int i = 0; i < Threading::ompNumThreads; ++i) srch_add[i] = new size_t[src_grid_size];

  size_t srch_corners = src_num_cell_corners;  // num of corners of srch cells

  double src_grid_bound_box[4];
  if (src_remap_grid_type == REMAP_GRID_TYPE_REG2D) reg2d_bound_box(src_grid, src_grid_bound_box);

  std::vector<WeightLinks> weightLinks(tgt_grid_size);

  double findex = 0;

  size_t num_srch_cells_stat[3];
  num_srch_cells_stat[0] = 0;
  num_srch_cells_stat[1] = 100000;
  num_srch_cells_stat[2] = 0;

  extern CellSearchMethod cellSearchMethod;

// Loop over destination grid

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(none)                                                                        \
    shared(findex, rsearch, src_remap_grid_type, tgt_remap_grid_type, src_grid_bound_box, rv, cdoVerbose, tgt_num_cell_corners, \
           target_cell_type, weightLinks, srch_corners, src_grid, tgt_grid, tgt_grid_size, src_grid_size, search, srch_add,     \
           tgt_grid_cell, num_srch_cells_stat, cellSearchMethod)
#endif
  for (size_t tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add)
    {
      double partial_weight;
      size_t src_cell_add;  // current linear address for source grid cell
      size_t num_srch_cells;
      size_t n, num_weights, num_weights_old;
      int ompthID = cdo_omp_get_thread_num();

#ifdef _OPENMP
#pragma omp atomic
#endif
      findex++;
      if (ompthID == 0) progressStatus(0, 1, findex / tgt_grid_size);

      weightLinks[tgt_cell_add].nlinks = 0;

      // Get search cells

      if ( cellSearchMethod == CellSearchMethod::spherepart )
        {
          set_yac_coordinates(tgt_remap_grid_type, tgt_cell_add, tgt_num_cell_corners, tgt_grid, &tgt_grid_cell[ompthID]);

          struct bounding_circle bnd_circle;
          yac_get_cell_bounding_circle(tgt_grid_cell[ompthID], &bnd_circle);

          struct dep_list result_list;
          yac_do_bnd_circle_search((struct grid_search *)rsearch.yacSearch, &bnd_circle, 1, &result_list);

          // unsigned num_matching_cells = yac_get_total_num_dependencies (result_list);
          // printf("num_matching_cells %u\n", num_matching_cells);

          num_srch_cells = result_list.num_deps_per_element[0];
          unsigned const *curr_neighs = yac_get_dependencies_of_element(result_list,0);
          for (size_t i = 0; i < num_srch_cells; ++i)
            srch_add[ompthID][i] = curr_neighs[i];

          yac_free_dep_list(&result_list);
        }
      else
        {
          if (src_remap_grid_type == REMAP_GRID_TYPE_REG2D && tgt_remap_grid_type == REMAP_GRID_TYPE_REG2D)
            {
              double tgt_cell_bound_box[4];
              boundbox_from_corners_reg2d(tgt_cell_add, tgt_grid->dims, tgt_grid->reg2d_corner_lon, tgt_grid->reg2d_corner_lat,
                                          tgt_cell_bound_box);
              restrict_boundbox(src_grid_bound_box, tgt_cell_bound_box);

              num_srch_cells = get_srch_cells_reg2d(src_grid->dims, src_grid->reg2d_corner_lat, src_grid->reg2d_corner_lon,
                                                    tgt_cell_bound_box, srch_add[ompthID]);

              if (num_srch_cells == 1 && src_grid->dims[0] == 1 && src_grid->dims[1] == 1
                  && IS_EQUAL(src_grid->reg2d_corner_lat[0], src_grid->reg2d_corner_lat[1])
                  && IS_EQUAL(src_grid->reg2d_corner_lon[0], src_grid->reg2d_corner_lon[1]))
                num_srch_cells = 0;
            }
          else if (src_remap_grid_type == REMAP_GRID_TYPE_REG2D)
            {
              double tgt_cell_bound_box[4];
              boundbox_from_corners1(tgt_cell_add, tgt_num_cell_corners, tgt_grid->cell_corner_lon, tgt_grid->cell_corner_lat,
                                     tgt_cell_bound_box);
              restrict_boundbox(src_grid_bound_box, tgt_cell_bound_box);

              num_srch_cells = get_srch_cells_reg2d(src_grid->dims, src_grid->reg2d_corner_lat, src_grid->reg2d_corner_lon,
                                                    tgt_cell_bound_box, srch_add[ompthID]);

              if (num_srch_cells == 1 && src_grid->dims[0] == 1 && src_grid->dims[1] == 1
                  && IS_EQUAL(src_grid->reg2d_corner_lat[0], src_grid->reg2d_corner_lat[1])
                  && IS_EQUAL(src_grid->reg2d_corner_lon[0], src_grid->reg2d_corner_lon[1]))
                num_srch_cells = 0;
            }
          else
            {
              float tgt_cell_bound_box_r[4];
              boundbox_from_corners1r(tgt_cell_add, tgt_num_cell_corners, tgt_grid->cell_corner_lon, tgt_grid->cell_corner_lat,
                                      tgt_cell_bound_box_r);

              num_srch_cells = get_srch_cells(tgt_cell_add, rsearch.tgtBins, rsearch.srcBins, tgt_cell_bound_box_r, srch_add[ompthID]);
            }
        }

      if (1 && cdoVerbose)
        {
          num_srch_cells_stat[0] += num_srch_cells;
          if (num_srch_cells < num_srch_cells_stat[1]) num_srch_cells_stat[1] = num_srch_cells;
          if (num_srch_cells > num_srch_cells_stat[2]) num_srch_cells_stat[2] = num_srch_cells;
        }

      if (0 && cdoVerbose) printf("tgt_cell_add %zu  num_srch_cells %zu\n", tgt_cell_add, num_srch_cells);

      if (num_srch_cells == 0) continue;

      if ( cellSearchMethod != CellSearchMethod::spherepart )
        set_yac_coordinates(tgt_remap_grid_type, tgt_cell_add, tgt_num_cell_corners, tgt_grid, &tgt_grid_cell[ompthID]);

      // Create search arrays

      cellSearchRealloc(num_srch_cells, search[ompthID]);

      double *partial_areas = &search[ompthID].partial_areas[0];
      double *partial_weights = &search[ompthID].partial_weights[0];
      struct grid_cell *src_grid_cells = &search[ompthID].src_grid_cells[0];

      for (n = 0; n < num_srch_cells; ++n)
        {
          src_cell_add = srch_add[ompthID][n];
          set_yac_coordinates(src_remap_grid_type, src_cell_add, srch_corners, src_grid, &src_grid_cells[n]);
        }

      if (tgt_num_cell_corners < 4 || target_cell_type == LON_LAT_CELL)
        {
          cdo_compute_overlap_areas(num_srch_cells, &search[ompthID], tgt_grid_cell[ompthID]);
        }
      else
        {
          cdo_compute_concave_overlap_areas(num_srch_cells, &search[ompthID], tgt_grid_cell[ompthID],
                                            tgt_grid->cell_center_lon[tgt_cell_add],
                                            tgt_grid->cell_center_lat[tgt_cell_add]);
        }

      double tgt_area = gridcell_area(tgt_grid_cell[ompthID]);
      // double tgt_area = cell_area(&tgt_grid_cell[ompthID]);

      for (num_weights = 0, n = 0; n < num_srch_cells; ++n)
        {
          if (partial_areas[n] > 0)
            {
              partial_areas[num_weights] = partial_areas[n];
              srch_add[ompthID][num_weights] = srch_add[ompthID][n];
              num_weights++;
            }
        }

      for (n = 0; n < num_weights; ++n) partial_weights[n] = partial_areas[n] / tgt_area;

      if (rv.normOpt == NormOpt::FRACAREA) yac_correct_weights((unsigned) num_weights, partial_weights);

      for (n = 0; n < num_weights; ++n) partial_weights[n] *= tgt_area;

      num_weights_old = num_weights;
      num_weights = 0;
      for (n = 0; n < num_weights_old; ++n)
        {
          src_cell_add = srch_add[ompthID][n];
          if (partial_weights[n] <= 0.) src_cell_add = src_grid_size;
          if (src_cell_add != src_grid_size)
            {
              partial_weights[num_weights] = partial_weights[n];
              srch_add[ompthID][num_weights] = src_cell_add;
              num_weights++;
            }
        }

      for (n = 0; n < num_weights; ++n)
        {
          partial_weight = partial_weights[n];
          src_cell_add = srch_add[ompthID][n];

#ifdef _OPENMP
#pragma omp atomic
#endif
          src_grid->cell_area[src_cell_add] += partial_weight;
        }

      num_weights_old = num_weights;
      for (num_weights = 0, n = 0; n < num_weights_old; ++n)
        {
          src_cell_add = srch_add[ompthID][n];

          /*
            Store the appropriate addresses and weights.
            Also add contributions to cell areas.
            The source grid mask is the master mask
          */
          if (src_grid->mask[src_cell_add])
            {
              partial_weights[num_weights] = partial_weights[n];
              srch_add[ompthID][num_weights] = src_cell_add;
              num_weights++;
            }
        }

      for (n = 0; n < num_weights; ++n)
        {
          partial_weight = partial_weights[n];
          src_cell_add = srch_add[ompthID][n];

#ifdef _OPENMP
#pragma omp atomic
#endif
          src_grid->cell_frac[src_cell_add] += partial_weight;

          tgt_grid->cell_frac[tgt_cell_add] += partial_weight;
        }

      storeWeightLinks(1, num_weights, srch_add[ompthID], partial_weights, tgt_cell_add, weightLinks);

      tgt_grid->cell_area[tgt_cell_add] = tgt_area;
      // printf("area %d %g %g\n", tgt_cell_add, tgt_grid->cell_area[tgt_cell_add], tgt_area);
    }

  progressStatus(0, 1, 1);

  if (1 && cdoVerbose)
    {
      cdoPrint("Num search cells min,mean,max :  %zu  %3.1f  %zu", num_srch_cells_stat[1],
               num_srch_cells_stat[0] / (double) tgt_grid_size, num_srch_cells_stat[2]);
    }

  // Finished with all cells: deallocate search arrays

  for (int ompthID = 0; ompthID < Threading::ompNumThreads; ++ompthID)
    {
      cellSearchFree(search[ompthID]);

      delete[] tgt_grid_cell[ompthID].coordinates_x;
      delete[] tgt_grid_cell[ompthID].coordinates_y;
      delete[] tgt_grid_cell[ompthID].coordinates_xyz;
      delete[] srch_add[ompthID];
    }

  weightLinksToRemapLinks(1, tgt_grid_size, weightLinks, rv);

  // Normalize weights using destination area if requested
  remapNormalizeWeights(tgt_grid, rv);

  if (cdoVerbose) cdoPrint("Total number of links = %zu", rv.num_links);

  for (size_t n = 0; n < src_grid_size; ++n)
    if (IS_NOT_EQUAL(src_grid->cell_area[n], 0)) src_grid->cell_frac[n] /= src_grid->cell_area[n];

  for (size_t n = 0; n < tgt_grid_size; ++n)
    if (IS_NOT_EQUAL(tgt_grid->cell_area[n], 0)) tgt_grid->cell_frac[n] /= tgt_grid->cell_area[n];

  // Perform some error checking on final weights
  if (lcheck)
    {
      remapCheckArea(src_grid_size, &src_grid->cell_area[0], "Source");
      remapCheckArea(tgt_grid_size, &tgt_grid->cell_area[0], "Target");

      remapVarsCheckWeights(rv);
    }

  if (cdoVerbose) cdoPrint("Cells search: %.2f seconds", cdo_get_wtime() - start);
}  // remapConservWeights

void
remapConserv(RemapSearch &rsearch, const double *restrict, double *restrict, double)
{
}  // remapConserv

/**
 * @file grid_search_utils.c
 *
 * @copyright Copyright  (C)  2014 Moritz Hanke <hanke@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 * URL: https://doc.redmine.dkrz.de/YAC/html/index.html
 *
 * This file is part of YAC.
 *
 * YAC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * YAC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with YAC.  If not, see <http://www.gnu.org/licenses/gpl.txt>.
 */

#include <stdlib.h>

#include "grid_search_utils.h"
#include "ensure_array_size.h"
#include "utils.h"
#include "string.h"

void yac_grid_search_utils_do_point_search_p(struct grid_search * search,
                                             struct grid * search_grid_data,
                                             struct grid * grid_data,
                                             struct dep_list * target_to_src_points) {

   unsigned num_target_corners;
   // determine the number of target points
   num_target_corners = yac_get_num_grid_corners(grid_data);

   struct dep_list tgt_corner_to_src_cell;

   // do the actual point search
   yac_do_point_search_c(search, grid_data, &tgt_corner_to_src_cell);

   unsigned i;
   unsigned * num_src_per_tgt;
   unsigned num_total_src_deps;
   unsigned curr_src_cell;

   // compute the number of src corners per cell and the total number of
   // "target corner to source corner"-dependencies
   num_src_per_tgt = calloc (num_target_corners, sizeof(num_src_per_tgt[0]));
   num_total_src_deps = 0;

   for (i = 0; i < num_target_corners; ++i) {
      if (tgt_corner_to_src_cell.num_deps_per_element[i] == 1) {

         curr_src_cell = yac_get_dependencies_of_element(tgt_corner_to_src_cell, i)[0];
         num_src_per_tgt[i] = yac_get_num_cell_corners(search_grid_data, curr_src_cell);
         num_total_src_deps += num_src_per_tgt[i];
      } else {
         num_src_per_tgt[i] = 0;
      }
   }

   unsigned j;
   unsigned * tgt_to_src;
   unsigned const * curr_src_corners;

   // generate dependency array
   tgt_to_src = malloc (num_total_src_deps * sizeof (tgt_to_src[0]));
   num_total_src_deps = 0;

   for (i = 0; i < num_target_corners; ++i) {

      if (tgt_corner_to_src_cell.num_deps_per_element[i] == 1) {

         curr_src_cell = yac_get_dependencies_of_element(tgt_corner_to_src_cell, i)[0];
         curr_src_corners = yac_get_cell_corner_indices (search_grid_data,
                                                     curr_src_cell);
         for (j = 0; j < num_src_per_tgt[i]; ++j) {

            tgt_to_src[num_total_src_deps++] = curr_src_corners[j];
         }
      }
   }

   yac_free_dep_list(&tgt_corner_to_src_cell);

   yac_set_dependencies(target_to_src_points, num_target_corners, num_src_per_tgt, tgt_to_src);
}

void yac_grid_search_utils_do_point_search_p2 (struct grid_search * search,
                                               struct grid * search_grid_data,
                                               double (*coordinates_xyz)[3],
                                               unsigned num_points,
                                               struct dep_list * target_to_src_points) {

   struct dep_list tgt_point_to_src_cell;

   // do the actual point search
   yac_do_point_search_c2(
      search, coordinates_xyz, num_points, &tgt_point_to_src_cell);

   unsigned i;
   unsigned * num_src_per_tgt;
   unsigned num_total_src_deps;
   unsigned * src_cells;

   // compute the number of src corners per cell and the total number of
   // "target corner to source corner"-dependencies
   num_src_per_tgt = calloc (num_points, sizeof(num_src_per_tgt[0]));
   num_total_src_deps = 0;

   src_cells = malloc (num_points * sizeof(*src_cells));

   for (i = 0; i < num_points; ++i) {

      unsigned num_cells = tgt_point_to_src_cell.num_deps_per_element[i];

      if (num_cells > 0) {

         unsigned const * dependencies =
            yac_get_dependencies_of_element(tgt_point_to_src_cell, i);

         unsigned src_cell = dependencies[0];

         for (unsigned i = 1; i < num_cells; ++i)
            if (src_cell > dependencies[i])
               src_cell = dependencies[i];

         src_cells[i] = src_cell;

         num_src_per_tgt[i] = yac_get_num_cell_corners(search_grid_data, src_cell);
         num_total_src_deps += num_src_per_tgt[i];
      } else {
         src_cells[i] = -1;
         num_src_per_tgt[i] = 0;
      }
   }

   unsigned j;
   unsigned * tgt_to_src;
   unsigned const * curr_src_corners;

   // generate dependency array
   tgt_to_src = malloc (num_total_src_deps * sizeof (tgt_to_src[0]));
   num_total_src_deps = 0;

   for (i = 0; i < num_points; ++i) {

      if (src_cells[i] != -1) {

         curr_src_corners = yac_get_cell_corner_indices (search_grid_data,
                                                     src_cells[i]);
         for (j = 0; j < num_src_per_tgt[i]; ++j) {

            tgt_to_src[num_total_src_deps++] = curr_src_corners[j];
         }
      }
   }

   free(src_cells);

   yac_free_dep_list(&tgt_point_to_src_cell);

   yac_set_dependencies(target_to_src_points, num_points, num_src_per_tgt, tgt_to_src);
}

static unsigned check_grid_cells(
  double point_coords[3], struct grid * grid,
  unsigned const * cell_indices, size_t num_cells,
  struct grid_cell * tmp_cell) {

  struct grid_cell cell = *tmp_cell;

  for (size_t i = 0; i < num_cells; ++i) {

    struct bounding_circle bnd_circle;
    yac_get_grid_cell2(grid, cell_indices[i], &cell, &bnd_circle);

    if (yac_point_in_cell2(point_coords, cell, bnd_circle)) {

      *tmp_cell = cell;
      return cell_indices[i];
    }
  }

  *tmp_cell = cell;
  return -1;
}

static unsigned check_point_grid_cells_CELL_(
  double point_coords[3], unsigned base_cell_index, struct grid * base_grid,
  struct grid * point_grid, struct grid_cell * tmp_cell) {

  // the base cell index is a corner of the point grid
  size_t num_point_cells =
    (size_t)yac_get_num_corner_cells(point_grid, base_cell_index);
  unsigned const * point_cell_indices =
    yac_get_corner_cell_indices(point_grid, base_cell_index);

  return check_grid_cells(point_coords, point_grid, point_cell_indices,
                          num_point_cells, tmp_cell);
}

static unsigned check_point_grid_cells_CORNER_(
  double point_coords[3], unsigned base_cell_index, struct grid * base_grid,
  struct grid * point_grid, struct grid_cell * tmp_cell) {

  struct bounding_circle bnd_circle;
  yac_get_grid_cell2(point_grid, base_cell_index, tmp_cell, &bnd_circle);

  if (yac_point_in_cell2(point_coords, *tmp_cell, bnd_circle))
    return base_cell_index;

  // check the neighbour cells
  struct dep_list cell_neigh_deps = yac_get_cell_neigh_dep_list(base_grid);
  size_t num_neigh_cells =
    (size_t)cell_neigh_deps.num_deps_per_element[base_cell_index];
  unsigned const * neigh_cells =
    yac_get_dependencies_of_element(cell_neigh_deps, base_cell_index);

  return check_grid_cells(point_coords, point_grid, neigh_cells,
                          num_neigh_cells, tmp_cell);
}

static unsigned check_point_grid_cells_EDGE_(
  double point_coords[3], unsigned base_cell_index, struct grid * base_grid,
  struct grid * point_grid, struct grid_cell * tmp_cell) {

  unsigned const * base_cell_edges =
    yac_get_cell_edge_indices(base_grid, base_cell_index);
  size_t num_base_cell_edges =
    (size_t)yac_get_num_cell_edges(base_grid, base_cell_index);

  // with this algorithm many cell are checked twice...

  for (size_t i = 0; i < num_base_cell_edges; ++i) {

    unsigned const * point_cell_indices =
      yac_get_corner_cell_indices(point_grid, base_cell_edges[i]);
    size_t num_point_cells =
      (size_t)yac_get_num_corner_cells(point_grid, base_cell_edges[i]);

    unsigned ret_val =
      check_grid_cells(point_coords, point_grid, point_cell_indices,
                       num_point_cells, tmp_cell);

    if (ret_val != -1) return ret_val;
  }

  return -1;
}

static int check_point_grid_cells(
  double coordinate_xyz[3], unsigned base_cell_index,
  struct grid * base_grid, struct grid * point_grid,
  unsigned (*check_point_grid_cells_)(
    double point_coords[3], unsigned base_cell_index,
    struct grid * base_grid, struct grid * point_grid,
    struct grid_cell * tmp_cell), struct grid_cell * tmp_cell,
  unsigned * num_matching_cell_corners, unsigned ** matching_cell_corners,
  unsigned * matching_cell_corners_array_size, size_t offset) {

  unsigned point_cell_index =
    check_point_grid_cells_(
      coordinate_xyz, base_cell_index, base_grid, point_grid, tmp_cell);

  // if a matching point cell was found
  if (point_cell_index != -1) {

    size_t num_point_cell_corners =
      (size_t)yac_get_num_cell_corners(point_grid, point_cell_index);

    unsigned const * corner_indices =
      yac_get_cell_corner_indices(point_grid, point_cell_index);

    ENSURE_ARRAY_SIZE(*matching_cell_corners, *matching_cell_corners_array_size,
                      num_point_cell_corners + offset);
    *num_matching_cell_corners = (unsigned)num_point_cell_corners;

    for (size_t i = 0; i < num_point_cell_corners; ++i)
      (*matching_cell_corners)[offset+i] = corner_indices[i];
  } else return 0;

  return 1;
}

void yac_grid_search_utils_do_point_search_p3(struct grid_search * search,
                                              double (*coordinates_xyz)[3],
                                              unsigned num_points,
                                              struct dep_list * target_to_src_points,
                                              struct points * points) {

  struct dep_list tgt_point_to_src_cell;

  // start with a cell search
  yac_do_point_search_c2(
    search, coordinates_xyz, num_points, &tgt_point_to_src_cell);

  unsigned * tgt_point_to_src_points_dep = NULL;
  unsigned tgt_point_to_src_points_dep_array_size = 0;
  unsigned * num_dep_per_tgt_point =
    calloc(num_points, sizeof(*num_dep_per_tgt_point));
  unsigned total_num_deps = 0;

  unsigned (*check_point_grid_cells_)(
    double point_coords[3], unsigned base_cell_index, struct grid * base_grid,
    struct grid * point_grid, struct grid_cell * tmp_cell);
  struct grid_cell tmp_cell;
  struct grid * base_grid = yac_get_base_grid(points);
  struct grid * point_grid = yac_get_point_grid(points);

  switch (points->location) {
    case (CORNER):
      check_point_grid_cells_ = check_point_grid_cells_CORNER_;
      break;
    case (CELL):
      check_point_grid_cells_ = check_point_grid_cells_CELL_;
      break;
    case (EDGE):
      check_point_grid_cells_ = check_point_grid_cells_EDGE_;
      break;
    default:
      yac_internal_abort_message(
        "ERROR: unknown point location type in do_point_search_p3",
        __FILE__, __LINE__);
  }

  yac_init_grid_cell(&tmp_cell);

  // for all points to be search
  for (unsigned i = 0; i < num_points; ++i) {

    unsigned num_matching_cells =
      tgt_point_to_src_cell.num_deps_per_element[i];
    unsigned const * matching_cells =
      yac_get_dependencies_of_element(tgt_point_to_src_cell,i);

    for (unsigned j = 0; j < num_matching_cells; ++j) {

      // check method point cells of current cell
      if (check_point_grid_cells(
            coordinates_xyz[i], matching_cells[j], base_grid, point_grid,
            check_point_grid_cells_, &tmp_cell, num_dep_per_tgt_point+i,
            &tgt_point_to_src_points_dep,
            &tgt_point_to_src_points_dep_array_size, total_num_deps)) {

        total_num_deps += num_dep_per_tgt_point[i];
        break;
      }
    }
  }

  yac_free_grid_cell(&tmp_cell);

  tgt_point_to_src_points_dep =
    realloc(tgt_point_to_src_points_dep,
            total_num_deps * sizeof(*tgt_point_to_src_points_dep));

  yac_set_dependencies(
    target_to_src_points, num_points, num_dep_per_tgt_point,
    tgt_point_to_src_points_dep);

  yac_free_dep_list(&tgt_point_to_src_cell);
}

void yac_grid_search_utils_do_point_search_p4 (struct grid_search * search,
                                               struct grid * search_grid_data,
                                               double coordinate_xyz[3],
                                               unsigned * n_points,
                                               unsigned * points_size,
                                               unsigned ** points) {

   struct dep_list tgt_point_to_src_cell;

   // do the actual point search
   yac_do_point_search_c2(search, (double(*)[3])&(coordinate_xyz[0]), 1, &tgt_point_to_src_cell);

   unsigned num_cells = tgt_point_to_src_cell.num_deps_per_element[0];

   if (num_cells > 0) {

      unsigned const * dependencies =
         yac_get_dependencies_of_element(tgt_point_to_src_cell, 0);

      unsigned src_cell = dependencies[0];

      for (unsigned i = 1; i < num_cells; ++i)
         if (src_cell > dependencies[i])
            src_cell = dependencies[i];

      *n_points = yac_get_num_cell_corners(search_grid_data, src_cell);
      ENSURE_ARRAY_SIZE(*points, *points_size, *n_points);

      unsigned const * src_corners =
         yac_get_cell_corner_indices(search_grid_data, src_cell);
      for (unsigned i = 0; i < *n_points; ++i)
         (*points)[i] = src_corners[i];
   } else {
      *n_points = 0;
   }

   yac_free_dep_list(&tgt_point_to_src_cell);
}

static unsigned get_num_point_grid_cells_CELL(
  struct grid * point_grid, struct grid * base_grid, unsigned num_cells,
  unsigned const * base_grid_cell_indices) {

  unsigned num_point_grid_cells = 0;

  // the base cell index is a corner of the point grid
  for (size_t i = 0; i < (size_t)num_cells; ++i)
    num_point_grid_cells +=
      yac_get_num_corner_cells(point_grid, base_grid_cell_indices[i]);

  return num_point_grid_cells;
}

static unsigned get_num_point_grid_cells_CORNER(
  struct grid * point_grid, struct grid * base_grid, unsigned num_cells,
  unsigned const * base_grid_cell_indices) {

  unsigned num_point_grid_cells = num_cells;

  // the neighbour cells
  struct dep_list cell_neigh_deps = yac_get_cell_neigh_dep_list(base_grid);
  for (size_t i = 0; i < (size_t)num_cells; ++i)
    num_point_grid_cells +=
      cell_neigh_deps.num_deps_per_element[base_grid_cell_indices[i]];

  return num_point_grid_cells;
}

static unsigned get_num_point_grid_cells_EDGE(
  struct grid * point_grid, struct grid * base_grid, unsigned num_cells,
  unsigned const * base_grid_cell_indices) {

  unsigned num_point_grid_cells = 0;

  for (size_t i = 0; i < (size_t)num_cells; ++i) {

    unsigned const * base_cell_edges =
      yac_get_cell_edge_indices(base_grid, base_grid_cell_indices[i]);
    size_t num_base_cell_edges =
      (size_t)yac_get_num_cell_edges(base_grid, base_grid_cell_indices[i]);

    // with this algorithm many cell are checked twice...

    for (size_t j = 0; j < num_base_cell_edges; ++j)
      num_point_grid_cells +=
        yac_get_num_corner_cells(point_grid, base_cell_edges[j]);
  }

  return num_point_grid_cells;
}

void get_point_grid_cell_indices_CELL(
  struct grid * point_grid, struct grid * base_grid, unsigned num_cells,
  unsigned const * base_grid_cell_indices, unsigned * point_grid_cell_indices) {

  size_t num_point_grid_cells = 0;

  // the base cell index is a corner of the point grid
  for (size_t i = 0; i < (size_t)num_cells; ++i) {

    // the base cell index is a corner of the point grid
    size_t curr_num_point_cells =
      (size_t)yac_get_num_corner_cells(point_grid, base_grid_cell_indices[i]);

    memcpy(
      point_grid_cell_indices + num_point_grid_cells,
      yac_get_corner_cell_indices(point_grid, base_grid_cell_indices[i]),
      curr_num_point_cells * sizeof(*point_grid_cell_indices));

    num_point_grid_cells += curr_num_point_cells;
  }
}

void get_point_grid_cell_indices_CORNER(
  struct grid * point_grid, struct grid * base_grid, unsigned num_cells,
  unsigned const * base_grid_cell_indices, unsigned * point_grid_cell_indices) {

  size_t num_point_grid_cells;

  for (num_point_grid_cells = 0; num_point_grid_cells < (size_t)num_cells;
       ++num_point_grid_cells)
    point_grid_cell_indices[num_point_grid_cells] =
      base_grid_cell_indices[num_point_grid_cells];

  // the neighbour cells
  struct dep_list cell_neigh_deps = yac_get_cell_neigh_dep_list(base_grid);
  for (size_t i = 0; i < (size_t)num_cells; ++i) {

    size_t curr_num_point_cells =
      (size_t)(cell_neigh_deps.num_deps_per_element[base_grid_cell_indices[i]]);

    memcpy(
      point_grid_cell_indices + num_point_grid_cells,
      yac_get_dependencies_of_element(
        cell_neigh_deps, base_grid_cell_indices[i]),
      curr_num_point_cells * sizeof(*point_grid_cell_indices));

    num_point_grid_cells += curr_num_point_cells;
  }
}

void get_point_grid_cell_indices_EDGE(
  struct grid * point_grid, struct grid * base_grid, unsigned num_cells,
  unsigned const * base_grid_cell_indices, unsigned * point_grid_cell_indices) {

  size_t num_point_grid_cells = 0;

  for (size_t i = 0; i < (size_t)num_cells; ++i) {

    unsigned const * base_cell_edges =
      yac_get_cell_edge_indices(base_grid, base_grid_cell_indices[i]);
    size_t num_base_cell_edges =
      (size_t)yac_get_num_cell_edges(base_grid, base_grid_cell_indices[i]);

    // with this algorithm many cell are checked twice...

    for (size_t j = 0; j < num_base_cell_edges; ++j) {

      size_t curr_num_point_cells =
        (size_t)yac_get_num_corner_cells(point_grid, base_cell_edges[j]);

      memcpy(
        point_grid_cell_indices + num_point_grid_cells,
        yac_get_corner_cell_indices(point_grid, base_cell_edges[j]),
        curr_num_point_cells * sizeof(*point_grid_cell_indices));

      num_point_grid_cells += curr_num_point_cells;
    }
  }
}

static int compare_uint (const void * a, const void * b) {
  return ( *(unsigned*)a - *(unsigned*)b );
}

void yac_grid_search_utils_do_point_search_c3 (struct grid_search * search,
                                               double (*coordinates_xyz)[3],
                                               unsigned num_points,
                                               struct dep_list * tgt_to_src_cells,
                                               struct points * points) {

  struct dep_list tgt_point_to_src_cell;

  // start with a cell search
  yac_do_point_search_c2(
    search, coordinates_xyz, num_points, &tgt_point_to_src_cell);

  unsigned (*get_num_point_grid_cells)(
    struct grid * point_grid, struct grid * base_grid, unsigned num_cells,
    unsigned const * base_grid_cell_indices);
  void (*get_point_grid_cell_indices)(
    struct grid * point_grid, struct grid * base_grid, unsigned num_cells,
    unsigned const * base_grid_cell_indices, unsigned * point_grid_cell_indices);
  struct grid * point_grid = yac_get_point_grid(points);
  struct grid * base_grid = yac_get_base_grid(points);

  switch (points->location) {
    case (CORNER):
      get_num_point_grid_cells = get_num_point_grid_cells_CORNER;
      get_point_grid_cell_indices = get_point_grid_cell_indices_CORNER;
      break;
    case (CELL):
      get_num_point_grid_cells = get_num_point_grid_cells_CELL;
      get_point_grid_cell_indices = get_point_grid_cell_indices_CELL;
      break;
    case (EDGE):
      get_num_point_grid_cells = get_num_point_grid_cells_EDGE;
      get_point_grid_cell_indices = get_point_grid_cell_indices_EDGE;
      break;
    default:
      yac_internal_abort_message(
        "ERROR: unknown point location type in do_point_search_p3",
        __FILE__, __LINE__);
  }

  unsigned * tgt_point_to_src_cell_dep = NULL;
  size_t tgt_point_to_src_points_dep_array_size = 0;
  unsigned * num_dep_per_tgt_point =
    malloc(num_points * sizeof(*num_dep_per_tgt_point));
  unsigned total_num_deps = 0;

  unsigned max_num_point_grid_cells = 0;
  {
    // for all points to be search
    for (unsigned i = 0; i < num_points; ++i) {

      // get number of potentially matching point grid cells
      unsigned num_point_grid_cells =
        get_num_point_grid_cells(
          point_grid, base_grid, tgt_point_to_src_cell.num_deps_per_element[i],
          yac_get_dependencies_of_element(tgt_point_to_src_cell,i));

      if (num_point_grid_cells > max_num_point_grid_cells)
        max_num_point_grid_cells = num_point_grid_cells;

      num_dep_per_tgt_point[i] = num_point_grid_cells;
    }
  }

  unsigned * tmp_point_grid_cell_indices =
    malloc(max_num_point_grid_cells * sizeof(*tmp_point_grid_cell_indices));
  struct grid_cell tmp_cell;
  yac_init_grid_cell(&tmp_cell);

  // for all points to be search
  for (unsigned i = 0; i < num_points; ++i) {

    // get point grid cell indices
    get_point_grid_cell_indices(
      point_grid, base_grid, tgt_point_to_src_cell.num_deps_per_element[i],
      yac_get_dependencies_of_element(tgt_point_to_src_cell,i),
      tmp_point_grid_cell_indices);

    unsigned num_point_grid_cells = num_dep_per_tgt_point[i];

    // sort point grid cell indices and remove duplicates
    qsort(tmp_point_grid_cell_indices, num_point_grid_cells,
          sizeof(*tmp_point_grid_cell_indices), compare_uint);
    yac_remove_duplicates_uint(
      tmp_point_grid_cell_indices, &num_point_grid_cells);

    ENSURE_ARRAY_SIZE(
      tgt_point_to_src_cell_dep, tgt_point_to_src_points_dep_array_size,
      total_num_deps + num_point_grid_cells);

    unsigned matching_point_grid_cells = 0;

    // check point grid cells
    for (unsigned j = 0; j < num_point_grid_cells; ++j) {

      struct bounding_circle bnd_circle;
      yac_get_grid_cell2(
        point_grid, tmp_point_grid_cell_indices[j], &tmp_cell, &bnd_circle);

      if (yac_point_in_cell2(coordinates_xyz[i], tmp_cell, bnd_circle)) {

        tgt_point_to_src_cell_dep[total_num_deps++] =
          tmp_point_grid_cell_indices[j];
        ++matching_point_grid_cells;
      }
    }

    num_dep_per_tgt_point[i] = matching_point_grid_cells;
  }

  yac_free_grid_cell(&tmp_cell);
  free(tmp_point_grid_cell_indices);

  tgt_point_to_src_cell_dep =
    realloc(tgt_point_to_src_cell_dep,
            total_num_deps * sizeof(*tgt_point_to_src_cell_dep));

  yac_set_dependencies(
    tgt_to_src_cells, num_points, num_dep_per_tgt_point,
    tgt_point_to_src_cell_dep);

  yac_free_dep_list(&tgt_point_to_src_cell);
}

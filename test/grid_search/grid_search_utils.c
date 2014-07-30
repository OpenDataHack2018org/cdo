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
 * URL: https://redmine.dkrz.de/doc/YAC/html/index.html
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

void grid_search_utils_do_point_search_p(struct grid_search * search,
                                         struct grid * search_grid_data,
                                         struct grid * grid_data,
                                         struct dep_list * target_to_src_points) {

   unsigned num_target_corners;
   // determine the number of target points
   num_target_corners = get_num_grid_corners(grid_data);

   struct dep_list tgt_corner_to_src_cell;
   
   // do the actual point search
   do_point_search_c(search, grid_data, &tgt_corner_to_src_cell);

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

         curr_src_cell = get_dependencies_of_element(tgt_corner_to_src_cell, i)[0];
         num_src_per_tgt[i] = get_num_cell_corners(search_grid_data, curr_src_cell);
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

         curr_src_cell = get_dependencies_of_element(tgt_corner_to_src_cell, i)[0];
         curr_src_corners = get_cell_corner_indices (search_grid_data,
                                                     curr_src_cell);
         for (j = 0; j < num_src_per_tgt[i]; ++j) {

            tgt_to_src[num_total_src_deps++] = curr_src_corners[j];
         }
      }
   }

   free_dep_list(&tgt_corner_to_src_cell);

   set_dependencies(target_to_src_points, num_target_corners, num_src_per_tgt, tgt_to_src);
}

void grid_search_utils_do_point_search_p2 (struct grid_search * search,
                                           struct grid * search_grid_data,
                                           double * x_coordinates,
                                           double * y_coordinates,
                                           unsigned num_points,
                                           struct dep_list * target_to_src_points) {

   struct dep_list tgt_point_to_src_cell;
   
   // do the actual point search
   do_point_search_c2(search, x_coordinates, y_coordinates, num_points,
                      &tgt_point_to_src_cell);

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
            get_dependencies_of_element(tgt_point_to_src_cell, i);

         unsigned src_cell = dependencies[0];

         for (unsigned i = 1; i < num_cells; ++i)
            if (src_cell > dependencies[i])
               src_cell = dependencies[i];

         src_cells[i] = src_cell;

         num_src_per_tgt[i] = get_num_cell_corners(search_grid_data, src_cell);
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

         curr_src_corners = get_cell_corner_indices (search_grid_data,
                                                     src_cells[i]);
         for (j = 0; j < num_src_per_tgt[i]; ++j) {

            tgt_to_src[num_total_src_deps++] = curr_src_corners[j];
         }
      }
   }

   free(src_cells);

   free_dep_list(&tgt_point_to_src_cell);

   set_dependencies(target_to_src_points, num_points, num_src_per_tgt, tgt_to_src);
}

static int check_point_grid_cells_(struct point point, double point_coords[3],
                                   unsigned cell_index, struct points * points,
                                   unsigned * matching_point_cell_index) {

   struct grid * point_grid;

   point_grid = get_point_grid(points);

   unsigned num_cells;
   unsigned const * cell_indices;

   // cell indices of the base grid are the corner indices of the point grid
   num_cells = get_num_corner_cells(point_grid, cell_index);
   cell_indices = get_corner_cell_indices(point_grid, cell_index);

   unsigned i;

   struct grid_cell cell;

   init_grid_cell(&cell);

   for (i = 0; i < num_cells; ++i) {

      struct bounding_circle bnd_circle;

      get_grid_cell2(point_grid, cell_indices[i], &cell, &bnd_circle);

      if (point_in_cell2(point, point_coords, cell, bnd_circle)) {

         *matching_point_cell_index = cell_indices[i];

         free_grid_cell(&cell);

         return 1 == 1;
      }
   }
   free_grid_cell(&cell);

   return 1 == 0;
}

static int check_point_grid_cells(double x_coordinate, double y_coordinate,
                                  unsigned cell_index, struct points * points,
                                  unsigned * num_matching_cell_corners,
                                  unsigned ** matching_cell_corners,
                                  unsigned * matching_cell_corners_array_size,
                                  unsigned offset) {

   unsigned i;

   struct point point;
   double point_coords[3];

   point.lon = x_coordinate;
   point.lat = y_coordinate;
   LLtoXYZ(x_coordinate, y_coordinate, point_coords);

   unsigned point_cell_index;
   int found;

   found = check_point_grid_cells_(point, point_coords, cell_index, points,
                                   &point_cell_index);

   if (!found) {

      struct grid * base_grid;

      base_grid = get_base_grid(points);

      unsigned num_neighbour_cells;

      struct dep_list cell_neigh_deps = get_cell_neigh_dep_list(base_grid);

      num_neighbour_cells = cell_neigh_deps.num_deps_per_element[cell_index];

      unsigned const * neigh_cells;

      neigh_cells = get_dependencies_of_element(cell_neigh_deps, cell_index);

      for (i = 0; i < num_neighbour_cells; ++i) {
         found = check_point_grid_cells_(point, point_coords, neigh_cells[i],
                                         points, &point_cell_index);
         if (found) break;
      }
   }

   if (found) {

      struct grid * point_grid;

      point_grid = get_point_grid(points);

      *num_matching_cell_corners = get_num_cell_corners(point_grid, point_cell_index);

      unsigned const * corner_indices;

      corner_indices = get_cell_corner_indices(point_grid, point_cell_index);

      ENSURE_ARRAY_SIZE(*matching_cell_corners, *matching_cell_corners_array_size,
                        *num_matching_cell_corners + offset);

      for (i = 0; i < *num_matching_cell_corners; ++i)
         (*matching_cell_corners)[offset+i] = corner_indices[i];
   }

   return found;
}

void grid_search_utils_do_point_search_p3(struct grid_search * search,
                                          double * x_coordinates,
                                          double * y_coordinates,
                                          unsigned num_points,
                                          struct dep_list * target_to_src_points,
                                          struct points * points) {

   // I do not know how to generate auxiliary grid cells from points that are on the
   // edges
   if (points->location == EDGE)
      abort_message ( "ERROR: unsupported point location in do_point_search_p3 (EDGE).",
                     __FILE__, __LINE__ );

   else if (points->location == CORNER) {
      do_point_search_p2(search, x_coordinates, y_coordinates, num_points,
                         target_to_src_points);
      return;

   } else if (points->location != CELL)
      abort_message("ERROR: unknown point location type in do_point_search_p3",
                    __FILE__, __LINE__);

   struct dep_list tgt_point_to_src_cell;
   
   // do the actual point search
   do_point_search_c2(search, x_coordinates, y_coordinates, num_points,
                      &tgt_point_to_src_cell);

   unsigned * tgt_point_to_src_points_dep;
   unsigned tgt_point_to_src_points_dep_array_size;
   unsigned * num_dep_per_tgt_point;
   unsigned total_num_deps;

   tgt_point_to_src_points_dep = NULL;
   tgt_point_to_src_points_dep_array_size = 0;
   num_dep_per_tgt_point = calloc(num_points, sizeof(*num_dep_per_tgt_point));
   total_num_deps = 0;

   unsigned i, j;

   // for all points to be search
   for (i = 0; i < num_points; ++i) {

      unsigned num_matching_cells;
      unsigned const * matching_cells;

      num_matching_cells = tgt_point_to_src_cell.num_deps_per_element[i];
      matching_cells = get_dependencies_of_element(tgt_point_to_src_cell,i);

      for (j = 0; j < num_matching_cells; ++j) {

         // check method point cells of current cell
         if (check_point_grid_cells(x_coordinates[i], y_coordinates[i],
                                    matching_cells[j], points,
                                    num_dep_per_tgt_point+i,
                                    &tgt_point_to_src_points_dep,
                                    &tgt_point_to_src_points_dep_array_size,
                                    total_num_deps)) {

            total_num_deps += num_dep_per_tgt_point[i];
            break;
         }
      }
   }

   tgt_point_to_src_points_dep = realloc(tgt_point_to_src_points_dep,
                                         total_num_deps * sizeof(*tgt_point_to_src_points_dep));

   set_dependencies(target_to_src_points, num_points, num_dep_per_tgt_point,
                    tgt_point_to_src_points_dep);

   free_dep_list(&tgt_point_to_src_cell);
}

void grid_search_utils_do_point_search_p4 (struct grid_search * search,
                                           struct grid * search_grid_data,
                                           double x_coordinate,
                                           double y_coordinate,
                                           unsigned * n_points,
                                           unsigned * points_size,
                                           unsigned ** points) {

   struct dep_list tgt_point_to_src_cell;
   
   // do the actual point search
   do_point_search_c2(search, &x_coordinate, &y_coordinate, 1,
                      &tgt_point_to_src_cell);

   unsigned num_cells = tgt_point_to_src_cell.num_deps_per_element[0];

   if (num_cells > 0) {

      unsigned const * dependencies =
         get_dependencies_of_element(tgt_point_to_src_cell, 0);

      unsigned src_cell = dependencies[0];

      for (unsigned i = 1; i < num_cells; ++i)
         if (src_cell > dependencies[i])
            src_cell = dependencies[i];

      *n_points = get_num_cell_corners(search_grid_data, src_cell);                                
      ENSURE_ARRAY_SIZE(*points, *points_size, *n_points);

      unsigned const * src_corners =
         get_cell_corner_indices(search_grid_data, src_cell);
      for (unsigned i = 0; i < *n_points; ++i)
         (*points)[i] = src_corners[i];
   } else {
      *n_points = 0;
   }

   free_dep_list(&tgt_point_to_src_cell);
}

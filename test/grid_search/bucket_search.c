/**
 * @file bucket_search.c
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *                                 Rene Redler <rene.redler@mpimet.mpg.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler <rene.redler@mpimet.mpg.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
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
#include <math.h>

#include "grid_search.h"
#include "search.h"
#include "grid_reg2d.h"
#include "utils.h"
#include "ensure_array_size.h"
#include "grid_search_utils.h"

static double const tol = 1.0e-12;

static
void bucket_search_do_cell_search (struct grid_search * search,
                                   struct grid * grid_data,
                                   struct dep_list * tgt_to_src_cells);
static
void bucket_search_do_cell_search_single (struct grid_search * search,
                                          struct grid_cell grid_cell,
                                          unsigned * n_cells, unsigned * cells_size,
                                          unsigned ** cells);
static
void bucket_search_do_point_search_c (struct grid_search * search,
                                      struct grid * grid_data,
                                      struct dep_list * tgt_to_src_cells);
static
void bucket_search_do_point_search_c2 (struct grid_search * search,
                                       double * x_coordinates,
                                       double * y_coordinates, unsigned num_points,
                                       struct dep_list * tgt_to_src_cells);
static
void bucket_search_do_point_search_p (struct grid_search * search,
                                      struct grid * grid_data,
                                      struct dep_list * target_to_src_points);
static
void bucket_search_do_point_search_p2 (struct grid_search * search,
                                       double * x_coordinates,
                                       double * y_coordinates, unsigned num_points,
                                       struct dep_list * target_to_src_points);
static
void bucket_search_do_point_search_p3 (struct grid_search * search,
                                       double * x_coordinates,
                                       double * y_coordinates, unsigned num_points,
                                       struct dep_list * target_to_src_points,
                                       struct points * points);
static
void bucket_search_do_point_search_p4 (struct grid_search * search,
                                       double x_coordinate, double y_coordinate,
                                       unsigned * n_points, unsigned * points_size,
                                       unsigned ** points);

static void delete_bucket_search(struct grid_search * search);

struct bucket_search {
   struct grid_search_vtable * vtable;
   struct grid * bucket_grid;
   double * bucket_grid_x_coords;
   double * bucket_grid_y_coords;
   struct grid * grid_data;
   unsigned num_buckets[2];
   struct dep_list bucket_to_cell;
};

static struct grid_search_vtable bucket_search_vtable =
{
   .do_cell_search        = bucket_search_do_cell_search,
   .do_cell_search_single = bucket_search_do_cell_search_single,
   .do_point_search_c     = bucket_search_do_point_search_c,
   .do_point_search_c2    = bucket_search_do_point_search_c2,
   .do_point_search_p     = bucket_search_do_point_search_p,
   .do_point_search_p2    = bucket_search_do_point_search_p2,
   .do_point_search_p3    = bucket_search_do_point_search_p3,
   .do_point_search_p4    = bucket_search_do_point_search_p4,
   .delete_grid_search    = delete_bucket_search
};

static void delete_bucket_search(struct grid_search * search) {

   struct bucket_search * bucket_search = (struct bucket_search *)search;

   free(bucket_search->bucket_grid_x_coords);
   free(bucket_search->bucket_grid_y_coords);
   delete_grid(bucket_search->bucket_grid);
   free_dep_list(&(bucket_search->bucket_to_cell));
   free(search);
}

// ==============================================================================
// Creation of the bucket grid
static
void init_buckets(struct bucket_search * bucket_search, double grid_extent[2][2],
                  unsigned num_cells) {

   double dx, dy;                            // grid spacing of the bucket grid
   double * coordinates_x, * coordinates_y;  // coordinates of the bucket grid
   unsigned bucket_cells[2];                 // number of bucket cells in x and y direction
   unsigned cyclic[2] = {0,0};               // marker to indicate cyclic coordinates
   unsigned i;

   bucket_cells[0] = MAX (sqrt((double)num_cells * 
                          (grid_extent[0][1] - grid_extent[0][0]) / 
                          (grid_extent[1][1] - grid_extent[1][0])),1);

   bucket_cells[1] = num_cells / bucket_cells[0];


   // if the grid is cyclic
   if (fabs((2*M_PI) - (grid_extent[0][1] - grid_extent[0][0])) < 1e-10)
      cyclic[0] = 1;

   // even though we do not need "bucket_cells[0]+1" coordinates in case of a
   // cyclic grid we still allocate them for the bisection search done later
   // MoHa: I do not like this
   coordinates_x = malloc((bucket_cells[0]+1)*sizeof(double));
   coordinates_y = malloc((bucket_cells[1]+1)*sizeof(double));

   dx = (grid_extent[0][1] - grid_extent[0][0]) / (double)(bucket_cells[0]);
   dy = (grid_extent[1][1] - grid_extent[1][0]) / (double)(bucket_cells[1]);

   for (i = 0; i < bucket_cells[0]; ++i) {
      coordinates_x[i] = grid_extent[0][0] + dx * i;
   }

   coordinates_x[bucket_cells[0]] = grid_extent[0][1];

   for (i = 0; i < bucket_cells[1]; ++i) {
      coordinates_y[i] = grid_extent[1][0] + dy * i;
   }
   coordinates_y[bucket_cells[1]] = grid_extent[1][1];

   bucket_search->num_buckets[0] = bucket_cells[0];
   bucket_search->num_buckets[1] = bucket_cells[1];

   bucket_search->bucket_grid_x_coords = coordinates_x;
   bucket_search->bucket_grid_y_coords = coordinates_y;

   // register the bucket grid as a regular grid in lon and lat

   bucket_search->bucket_grid = reg2d_grid_new(coordinates_x, coordinates_y,
                                               bucket_cells, cyclic);

   init_dep_list(&bucket_search->bucket_to_cell);
}

// ==============================================================================
// Search for the grid cells on the bucket grid
//------------------------------------------------------------------
// Because of the properties of the bucket grid this search is
// being optimised by usage of the function bisection_search.
//------------------------------------------------------------------
static
void search_on_bucket_grid (struct bucket_search * bucket_search, struct grid * grid_data,
                            struct dep_list * cell_to_bucket) {

   unsigned num_x_coords, num_y_coords;
   int * found_x, * position_x, * found_y, * position_y;
   double const * x_coords;
   double const * y_coords;

   x_coords = get_x_coords(grid_data);
   y_coords = get_y_coords(grid_data);

   num_x_coords = get_size_x_coords(grid_data);
   num_y_coords = get_size_y_coords(grid_data);

   found_x    = malloc (num_x_coords * sizeof(*found_x));
   found_y    = malloc (num_y_coords * sizeof(*found_y));
   position_x = malloc (num_x_coords * sizeof(*position_x));
   position_y = malloc (num_y_coords * sizeof(*position_y));

   bisection_search(x_coords, num_x_coords, 
            bucket_search->bucket_grid_x_coords, 
            bucket_search->num_buckets[0]+1, 
            position_x, found_x, 2*M_PI);

   bisection_search(y_coords, num_y_coords, 
            bucket_search->bucket_grid_y_coords, 
            bucket_search->num_buckets[1]+1, 
            position_y, found_y, 0);

   //----------------------------------------------
   // generate initial cell->bucket dependency list
   //----------------------------------------------

   unsigned curr_cell_x_corner, curr_cell_y_corner;
   unsigned num_cells;
   unsigned * cell_to_bucket_dependencies;
   unsigned * num_buckets_per_cell;
   struct dep_list initial_cell_to_bucket;
   unsigned cell_to_bucket_dependencies_size;
   unsigned num_total_dependencies;

   unsigned const * curr_x_coords, * curr_y_coords;

   unsigned i, j, k, num_corners;

   num_cells = get_num_grid_cells(grid_data);
   cell_to_bucket_dependencies_size = num_cells * 4;
   cell_to_bucket_dependencies = malloc (cell_to_bucket_dependencies_size *
                                         sizeof (cell_to_bucket_dependencies[0]));
   num_buckets_per_cell = calloc (num_cells, sizeof (num_buckets_per_cell[0]));

   num_total_dependencies = 0;

   // for all cells
   for (i = 0; i < num_cells; ++i) {

      num_corners = get_num_cell_corners(grid_data, i);

      curr_x_coords = get_cell_x_coord_indices(grid_data, i);
      curr_y_coords = get_cell_y_coord_indices(grid_data, i);

      for (j = 0; j < num_corners; ++j) {

         curr_cell_x_corner = curr_x_coords[j];
         curr_cell_y_corner = curr_y_coords[j];

         if (!found_x[curr_cell_x_corner] || !found_y[curr_cell_y_corner])
            continue;

         unsigned num_matching_buckets = 1;
         unsigned bucket_flags = 0;
         unsigned matching_buckets[4] = {(unsigned)position_x[curr_cell_x_corner] +
                                         (unsigned)position_y[curr_cell_y_corner] *
                                         bucket_search->num_buckets[0]};

         // if the current point falls onto the left edge of the bucket
         if ((position_x[curr_cell_x_corner] != 0) &&
             (fabs(x_coords[curr_cell_x_corner] -
                   bucket_search->bucket_grid_x_coords[position_x[curr_cell_x_corner]]) < tol)) {

            matching_buckets[num_matching_buckets++] =
               (unsigned)position_x[curr_cell_x_corner] - 1 +
               (unsigned)position_y[curr_cell_y_corner] *
               bucket_search->num_buckets[0];

            bucket_flags |= 1;

         // if the current point falls onto the right edge of the bucket
         } else if ((position_x[curr_cell_x_corner] != bucket_search->num_buckets[0] - 1) &&
             (fabs(x_coords[curr_cell_x_corner] -
                   bucket_search->bucket_grid_x_coords[position_x[curr_cell_x_corner] + 1]) < tol)) {

            matching_buckets[num_matching_buckets++] =
               (unsigned)position_x[curr_cell_x_corner] + 1 +
               (unsigned)position_y[curr_cell_y_corner] *
               bucket_search->num_buckets[0];

            bucket_flags |= 2;
         }

         // if the current point falls onto the lower edge of the bucket
         if ((position_y[curr_cell_y_corner] != 0) &&
             (fabs(y_coords[curr_cell_y_corner] -
                   bucket_search->bucket_grid_y_coords[position_y[curr_cell_y_corner]]) < tol)) {

            matching_buckets[num_matching_buckets++] =
               (unsigned)position_x[curr_cell_x_corner] +
               ((unsigned)position_y[curr_cell_y_corner] - 1) *
               bucket_search->num_buckets[0];

            bucket_flags |= 4;

         // if the current point falls onto the upper edge of the bucket
         } else if ((position_y[curr_cell_y_corner] != bucket_search->num_buckets[1] - 1) &&
             (fabs(y_coords[curr_cell_y_corner] -
                   bucket_search->bucket_grid_y_coords[position_y[curr_cell_y_corner]+1]) < tol)) {

            matching_buckets[num_matching_buckets++] =
               (unsigned)position_x[curr_cell_x_corner] +
               ((unsigned)position_y[curr_cell_y_corner] + 1) *
               bucket_search->num_buckets[0];

            bucket_flags |= 8;
         }

         if (bucket_flags == 5)
            matching_buckets[num_matching_buckets++] =
               (unsigned)position_x[curr_cell_x_corner] - 1 +
               ((unsigned)position_y[curr_cell_y_corner] - 1) *
               bucket_search->num_buckets[0];
         else if (bucket_flags == 9)
            matching_buckets[num_matching_buckets++] =
               (unsigned)position_x[curr_cell_x_corner] - 1 +
               ((unsigned)position_y[curr_cell_y_corner] + 1) *
               bucket_search->num_buckets[0];
         else if (bucket_flags == 6)
            matching_buckets[num_matching_buckets++] =
               (unsigned)position_x[curr_cell_x_corner] + 1 +
               ((unsigned)position_y[curr_cell_y_corner] - 1) *
               bucket_search->num_buckets[0];
         else if (bucket_flags == 10)
            matching_buckets[num_matching_buckets++] =
               (unsigned)position_x[curr_cell_x_corner] + 1 +
               ((unsigned)position_y[curr_cell_y_corner] + 1) *
               bucket_search->num_buckets[0];

         if (num_matching_buckets > 4)
            abort_message("ERROR: too many matching buckets\n", __FILE__, __LINE__);

         ENSURE_ARRAY_SIZE(cell_to_bucket_dependencies,
                           cell_to_bucket_dependencies_size,
                           num_total_dependencies+num_matching_buckets);

         for (k = 0; k < num_matching_buckets; ++k) {

            cell_to_bucket_dependencies[num_total_dependencies++] = matching_buckets[k];
            num_buckets_per_cell[i]++;
         }
      }
   }

   if (cell_to_bucket_dependencies_size != num_total_dependencies)
      cell_to_bucket_dependencies = realloc(cell_to_bucket_dependencies,
                                            num_total_dependencies *
                                            sizeof(*cell_to_bucket_dependencies));

   init_dep_list(&initial_cell_to_bucket);
   set_dependencies(&initial_cell_to_bucket, num_cells, num_buckets_per_cell,
                    cell_to_bucket_dependencies);

   //----------------------------------------------
   // generate cell->bucket dependency list
   //----------------------------------------------

   init_dep_list(cell_to_bucket);
   find_overlapping_cells (grid_data, bucket_search->bucket_grid,
                           initial_cell_to_bucket, cell_to_bucket);
   //-------------
   // free memory
   //-------------

   free_dep_list(&initial_cell_to_bucket);
   free(found_x); free(found_y); free(position_x); free(position_y);
}

static
void search_on_bucket_grid_single (struct bucket_search * bucket_search,
                                   struct grid_cell grid_cell,
                                   unsigned ** deps, unsigned * deps_size,
                                   unsigned * num_deps, unsigned src_index,
                                   unsigned * tgts_already_touched,
                                   unsigned ** stack, unsigned * stack_size) {

   int found_x[grid_cell.num_corners];
   int position_x[grid_cell.num_corners];
   int found_y[grid_cell.num_corners];
   int position_y[grid_cell.num_corners];

   bisection_search(grid_cell.coordinates_x, grid_cell.num_corners, 
                    bucket_search->bucket_grid_x_coords, 
                    bucket_search->num_buckets[0]+1, 
                    position_x, found_x, 2*M_PI);

   bisection_search(grid_cell.coordinates_y, grid_cell.num_corners, 
                    bucket_search->bucket_grid_y_coords, 
                    bucket_search->num_buckets[1]+1, 
                    position_y, found_y, 0);

   //----------------------------------------------
   // generate initial cell->bucket dependency list
   //----------------------------------------------

   unsigned j, k;
   unsigned * buckets = NULL;
   unsigned buckets_size = 0;
   unsigned num_buckets = 0;

   for (j = 0; j < grid_cell.num_corners; ++j) {

      if (!found_x[j] || !found_y[j])
         continue;

      unsigned num_matching_buckets = 1;
      unsigned bucket_flags = 0;
      unsigned matching_buckets[4] = {(unsigned)position_x[j] +
                                      (unsigned)position_y[j] *
                                      bucket_search->num_buckets[0]};

      // if the current point falls onto the left edge of the bucket
      if ((position_x[j] != 0) &&
          (fabs(grid_cell.coordinates_x[j] -
                bucket_search->bucket_grid_x_coords[position_x[j]]) < tol)) {

         matching_buckets[num_matching_buckets++] =
            (unsigned)position_x[j] - 1 +
            (unsigned)position_y[j] *
            bucket_search->num_buckets[0];

         bucket_flags |= 1;

      // if the current point falls onto the right edge of the bucket
      } else if ((position_x[j] != bucket_search->num_buckets[0] - 1) &&
          (fabs(grid_cell.coordinates_x[j] -
                bucket_search->bucket_grid_x_coords[position_x[j] + 1]) < tol)) {

         matching_buckets[num_matching_buckets++] =
            (unsigned)position_x[j] + 1 +
            (unsigned)position_y[j] *
            bucket_search->num_buckets[0];

         bucket_flags |= 2;
      }

      // if the current point falls onto the lower edge of the bucket
      if ((position_y[j] != 0) &&
          (fabs(grid_cell.coordinates_y[j] -
                bucket_search->bucket_grid_y_coords[position_y[j]]) < tol)) {

         matching_buckets[num_matching_buckets++] =
            (unsigned)position_x[j] +
            ((unsigned)position_y[j] - 1) *
            bucket_search->num_buckets[0];

         bucket_flags |= 4;

      // if the current point falls onto the upper edge of the bucket
      } else if ((position_y[j] != bucket_search->num_buckets[1] - 1) &&
          (fabs(grid_cell.coordinates_y[j] -
                bucket_search->bucket_grid_y_coords[position_y[j]+1]) < tol)) {

         matching_buckets[num_matching_buckets++] =
            (unsigned)position_x[j] +
            ((unsigned)position_y[j] + 1) *
            bucket_search->num_buckets[0];

         bucket_flags |= 8;
      }

      if (bucket_flags == 5)
         matching_buckets[num_matching_buckets++] =
            (unsigned)position_x[j] - 1 +
            ((unsigned)position_y[j] - 1) *
            bucket_search->num_buckets[0];
      else if (bucket_flags == 9)
         matching_buckets[num_matching_buckets++] =
            (unsigned)position_x[j] - 1 +
            ((unsigned)position_y[j] + 1) *
            bucket_search->num_buckets[0];
      else if (bucket_flags == 6)
         matching_buckets[num_matching_buckets++] =
            (unsigned)position_x[j] + 1 +
            ((unsigned)position_y[j] - 1) *
            bucket_search->num_buckets[0];
      else if (bucket_flags == 10)
         matching_buckets[num_matching_buckets++] =
            (unsigned)position_x[j] + 1 +
            ((unsigned)position_y[j] + 1) *
            bucket_search->num_buckets[0];

      if (num_matching_buckets > 4)
         abort_message("ERROR: too many matching buckets\n", __FILE__, __LINE__);

      ENSURE_ARRAY_SIZE(buckets, buckets_size,
                        num_buckets+num_matching_buckets);

      for (k = 0; k < num_matching_buckets; ++k)
         buckets[num_buckets++] = matching_buckets[k];
   }

   //----------------------------------------------
   // generate cell->bucket dependency list
   //----------------------------------------------

   struct bounding_circle bnd_circle;

   get_cell_bounding_circle(grid_cell, &bnd_circle);
   find_overlapping_cells_s (grid_cell, bnd_circle, bucket_search->bucket_grid,
                             buckets, num_buckets, deps, deps_size, num_deps,
                             src_index, tgts_already_touched, stack, stack_size);

   free(buckets);
}


// ==============================================================================
// Generation of the grid_search object
struct grid_search * bucket_search_new (struct grid * grid_data) {

   struct bucket_search * search;

   search = malloc(1 * sizeof(*search));

   search->vtable = &bucket_search_vtable;

   double grid_extent[2][2];

   if (get_num_grid_corners(grid_data) == 0)
      abort_message("ERROR: bucket_search_new empty input grid\n",
                    __FILE__, __LINE__);

   //-------------------------
   // initiate the search data
   //-------------------------

   get_2d_grid_extent(grid_data, grid_extent);

   search->grid_data = grid_data;

   init_buckets (search, grid_extent, get_num_grid_cells(grid_data));

   //----------------------------------------------
   // generate cell->bucket dependency list
   //----------------------------------------------

   struct dep_list cell_to_bucket;

   search_on_bucket_grid (search, grid_data, &cell_to_bucket);

   //----------------------------------------------
   // generate bucket->cell dependency list
   //----------------------------------------------

   invert_dep_list(cell_to_bucket, &search->bucket_to_cell);

   //-------------
   // free memory
   //-------------
   
   free_dep_list(&cell_to_bucket);

   return (struct grid_search *)search;
}

// ==============================================================================
// This routine searches the given grid data in the data associated with
//  the search handle and returns a list of dependencies between the grids.
static
void bucket_search_do_cell_search (struct grid_search * search,
                                   struct grid * grid_data,
                                   struct dep_list * tgt_to_src_cells) {

   struct bucket_search * bucket_search = (struct bucket_search *)search;

   //--------------------------------
   // find grid cells on bucket grid
   //--------------------------------

   struct dep_list tgt_to_bucket;

   search_on_bucket_grid (bucket_search, grid_data, &tgt_to_bucket);

   //----------------------------------------------------------
   // use initial cell-to-bucket-results to fill initial result
   // for src to tgt grid dependencies
   //----------------------------------------------------------

   struct dep_list initial_tgt_to_src;
   unsigned * initial_num_srcs_per_tgt;
   unsigned * initial_tgt_to_src_dependencies;
   unsigned curr_size_tgt_to_src_dependencies;
   unsigned num_total_initial_tgt_to_src_deps;
   unsigned * src_cell_mask; // array that helps double entries of the same tgt-to-src
                             // dependency in initial_tgt_to_src_dependencies
   unsigned const * curr_buckets;
   unsigned curr_bucket;
   unsigned const * curr_src_cells;
   unsigned num_tgt_grid_cells;
   unsigned num_src_grid_cells;

   unsigned i, j, k;

   num_tgt_grid_cells = get_num_grid_cells(grid_data);
   initial_num_srcs_per_tgt = calloc(num_tgt_grid_cells, sizeof(initial_num_srcs_per_tgt[0]));
   num_src_grid_cells = get_num_grid_cells(bucket_search->grid_data);
   src_cell_mask = calloc(num_src_grid_cells, sizeof(src_cell_mask[0]));

   num_total_initial_tgt_to_src_deps = 0;
   curr_size_tgt_to_src_dependencies = 0;
   initial_tgt_to_src_dependencies = NULL;

   // for all target grid cells
   for (i = 0; i < num_tgt_grid_cells; ++i) {

      curr_buckets = get_dependencies_of_element(tgt_to_bucket, i);

      // for all matching bucket grid cells of the current tgt grid cell
      for (j = 0; j < tgt_to_bucket.num_deps_per_element[i]; ++j) {

         curr_bucket = curr_buckets[j];

         if (curr_bucket >= bucket_search->bucket_to_cell.num_elements)
            continue;

         curr_src_cells = get_dependencies_of_element(
            (bucket_search->bucket_to_cell), curr_bucket);

         // for all src grid cells that overlap with the current bucket cell
         for (k = 0; k < bucket_search->bucket_to_cell.num_deps_per_element[curr_bucket]; ++k) {

            // check for double src entries
            if (src_cell_mask[curr_src_cells[k]] != i + 1) {

               ENSURE_ARRAY_SIZE(initial_tgt_to_src_dependencies,
                                 curr_size_tgt_to_src_dependencies,
                                 num_total_initial_tgt_to_src_deps+1);
               initial_tgt_to_src_dependencies[num_total_initial_tgt_to_src_deps++] =
                  curr_src_cells[k];
               initial_num_srcs_per_tgt[i]++;
               src_cell_mask[curr_src_cells[k]] = i + 1;
            }
         } // (k = 0; k < bucket_search->bucket_to_cell.num_deps_per_element[curr_bucket], ++k)
      } // (j = 0; j < tgt_to_bucket.num_deps_per_element[i]; ++j)
   } // (i = 0; i < num_tgt_grid_cells; ++i)

   init_dep_list(&initial_tgt_to_src);
   set_dependencies(&initial_tgt_to_src, num_tgt_grid_cells, initial_num_srcs_per_tgt,
                    initial_tgt_to_src_dependencies);

   //----------------------------------------------
   // generate actual src to tgt grid dependencies
   //----------------------------------------------

   find_overlapping_cells(grid_data, bucket_search->grid_data,
                          initial_tgt_to_src, tgt_to_src_cells);
   //-------------
   // free memory
   //-------------

   free(src_cell_mask);

   free_dep_list(&tgt_to_bucket);
   free_dep_list(&initial_tgt_to_src);
}

static
void bucket_search_do_cell_search_single (struct grid_search * search,
                                          struct grid_cell cell,
                                          unsigned * n_cells, unsigned * cells_size,
                                          unsigned ** cells) {

   struct bucket_search * bucket_search = (struct bucket_search *)search;

   //--------------------------------
   // find grid cells on bucket grid
   //--------------------------------

   unsigned * cell_to_bucket = NULL, cell_to_bucket_size = 0, num_buckets = 0;
   unsigned * helper_array =
      calloc(MAX(get_num_grid_cells(bucket_search->bucket_grid),
                 get_num_grid_cells(bucket_search->grid_data)),
             sizeof(*helper_array));
   unsigned * stack = NULL;
   unsigned stack_size = 0;

   search_on_bucket_grid_single (bucket_search, cell, &cell_to_bucket,
                                 &cell_to_bucket_size, &num_buckets, 0,
                                 helper_array, &stack, &stack_size);

   //----------------------------------------------------------
   // use initial cell-to-bucket-results to fill initial result
   // for src to tgt grid dependencies
   //----------------------------------------------------------

   unsigned * initial_tgt_to_src_cells = NULL;
   unsigned initial_tgt_to_src_cells_size = 0;
   unsigned initial_num_cells = 0;

   unsigned curr_bucket;
   unsigned const * curr_src_cells;

   unsigned j, k;

   // for all matching bucket grid cells of the current tgt grid cell
   for (j = 0; j < num_buckets; ++j) {

      curr_bucket = cell_to_bucket[j];

      if (curr_bucket >= bucket_search->bucket_to_cell.num_elements)
         continue;

      curr_src_cells = get_dependencies_of_element(
         (bucket_search->bucket_to_cell), curr_bucket);

      // for all src grid cells that overlap with the current bucket cell
      for (k = 0; k < bucket_search->bucket_to_cell.num_deps_per_element[curr_bucket]; ++k) {

         // check for double src entries
         if (helper_array[curr_src_cells[k]] != 2) {

            ENSURE_ARRAY_SIZE(initial_tgt_to_src_cells,
                              initial_tgt_to_src_cells_size,
                              initial_num_cells+1);
            initial_tgt_to_src_cells[initial_num_cells++] = curr_src_cells[k];
            helper_array[curr_src_cells[k]] = 2;
         }
      } // (k = 0; k < bucket_search->bucket_to_cell.num_deps_per_element[curr_bucket], ++k)
   } // (j = 0; j < tgt_to_bucket.num_deps_per_element[i]; ++j)

   //----------------------------------------------
   // generate actual src to tgt grid dependencies
   //----------------------------------------------

   struct bounding_circle bnd_circle;

   get_cell_bounding_circle(cell, &bnd_circle);
   find_overlapping_cells_s (cell, bnd_circle, bucket_search->grid_data,
                             initial_tgt_to_src_cells, initial_num_cells,
                             cells, cells_size, n_cells, 3, helper_array,
                             &stack, &stack_size);
                             
   //-------------
   // free memory
   //-------------

   free(initial_tgt_to_src_cells);
   free(helper_array);
   free(stack);
   free(cell_to_bucket);
}

// ==============================================================================
// This routine searches the given grid data in the data associated with
// the search handle and returns a list of dependencies between the grids.
static
void bucket_search_do_point_search_c (struct grid_search * search,
                                      struct grid * grid_data,
                                      struct dep_list * tgt_to_src_cells) {

   struct bucket_search * bucket_search = (struct bucket_search *)search;

   //--------------------------------
   // find grid cells on bucket grid
   //--------------------------------

   unsigned num_x_coords, num_y_coords;
   int * found_x, * position_x, * found_y, * position_y;
   double const * x_coords;
   double const * y_coords;

   x_coords = get_x_coords(grid_data);
   y_coords = get_y_coords(grid_data);

   num_x_coords = get_size_x_coords(grid_data);
   num_y_coords = get_size_y_coords(grid_data);

   found_x    = malloc (num_x_coords * sizeof(*found_x));
   found_y    = malloc (num_y_coords * sizeof(*found_y));
   position_x = malloc (num_x_coords * sizeof(*position_x));
   position_y = malloc (num_y_coords * sizeof(*position_y));

   bisection_search(x_coords, num_x_coords, 
            bucket_search->bucket_grid_x_coords, 
            bucket_search->num_buckets[0]+1, 
            position_x, found_x, 2*M_PI);

   bisection_search(y_coords, num_y_coords, 
            bucket_search->bucket_grid_y_coords, 
            bucket_search->num_buckets[1]+1, 
            position_y, found_y, 0);

   //----------------------------------------------------------
   // use inital cell-to-bucket-results to find the find the
   // actual tgt point to src cell dependency
   //----------------------------------------------------------

   unsigned i, j;
   unsigned curr_bucket;
   unsigned curr_corner_index_x, curr_corner_index_y;
   unsigned const * curr_src_cells;
   struct grid_cell src_cell;
   unsigned num_tgt_corners;
   unsigned * tgt_to_src_cell;
   unsigned * num_src_per_tgt;
   unsigned num_total_dependencies;

   num_tgt_corners = get_num_grid_corners(grid_data);
   num_total_dependencies = 0;
   num_src_per_tgt = calloc (num_tgt_corners, sizeof (num_src_per_tgt[0]));
   tgt_to_src_cell = malloc (num_tgt_corners * sizeof (tgt_to_src_cell[0]));

   init_grid_cell(&src_cell);

   // for all corners of the target grid
   for (i = 0; i < num_tgt_corners; ++i) {

      curr_corner_index_x = get_corner_x_coord_index(grid_data, i);
      curr_corner_index_y = get_corner_y_coord_index(grid_data, i);

      // if the current target corner was not found on the source grid
      if (!found_x[curr_corner_index_x] || !found_y[curr_corner_index_y]) continue; // continue with next corner

      // compute the current bucket cell index
      curr_bucket = (unsigned)position_x[curr_corner_index_x] +
         (unsigned)position_y[curr_corner_index_y] *
         bucket_search->num_buckets[0];

      if (curr_bucket >= bucket_search->bucket_to_cell.num_elements)
         continue;

      curr_src_cells = get_dependencies_of_element(
         (bucket_search->bucket_to_cell), curr_bucket);
      // for all source cells associated with the current bucket cell
      for (j = 0; j < bucket_search->bucket_to_cell.num_deps_per_element[curr_bucket]; ++j) {

         get_grid_cell(bucket_search->grid_data, curr_src_cells[j],
                       &src_cell);

         struct point tgt_point;
         double tgt_coords[3];

         tgt_point.lon = x_coords[curr_corner_index_x];
         tgt_point.lat = y_coords[curr_corner_index_y];

         LLtoXYZ(tgt_point.lon, tgt_point.lat, tgt_coords);

         // if the target point is within the current source cell
         if (point_in_cell(tgt_point, tgt_coords, src_cell)) {

            num_src_per_tgt[i] = 1;
            tgt_to_src_cell[num_total_dependencies++] = curr_src_cells[j];
            break; // continue with next target point
         }
      } // for j
   } // for i

   tgt_to_src_cell = realloc (tgt_to_src_cell, num_total_dependencies *
                              sizeof(tgt_to_src_cell[0]));

   set_dependencies(tgt_to_src_cells, num_tgt_corners, num_src_per_tgt, tgt_to_src_cell);
   //-------------
   // free memory
   //-------------

   free_grid_cell(&src_cell);
   free(found_x); free(found_y); free(position_x); free(position_y);
}

// this routine searches the given point data in the data associated with the search
// handle and returns a dependency list containing the source cell for each point
static
void bucket_search_do_point_search_c2 (struct grid_search * search,
                                       double * x_coordinates,
                                       double * y_coordinates, unsigned num_points,
                                       struct dep_list * tgt_to_src_cells) {

   struct bucket_search * bucket_search = (struct bucket_search *)search;

   //--------------------------------
   // find grid cells on bucket grid
   //--------------------------------

   int * found_x, * position_x, * found_y, * position_y;

   found_x    = malloc (num_points * sizeof(*found_x));
   found_y    = malloc (num_points * sizeof(*found_y));
   position_x = malloc (num_points * sizeof(*position_x));
   position_y = malloc (num_points * sizeof(*position_y));

   bisection_search(x_coordinates, num_points, 
            bucket_search->bucket_grid_x_coords, 
            bucket_search->num_buckets[0]+1, 
            position_x, found_x, 2*M_PI);

   bisection_search(y_coordinates, num_points, 
            bucket_search->bucket_grid_y_coords, 
            bucket_search->num_buckets[1]+1, 
            position_y, found_y, 0);

   //----------------------------------------------------------
   // use initial cell-to-bucket-results to find the find the
   // actual tgt point to src cell dependency
   //----------------------------------------------------------

   unsigned i, j;
   unsigned curr_bucket;
   unsigned const * curr_src_cells;
   struct grid_cell src_cell;
   unsigned * tgt_to_src_cell;
   unsigned * num_src_per_tgt;
   unsigned num_total_dependencies;

   num_total_dependencies = 0;
   num_src_per_tgt = calloc (num_points, sizeof (num_src_per_tgt[0]));
   tgt_to_src_cell = malloc (num_points * sizeof (tgt_to_src_cell[0]));

   init_grid_cell(&src_cell);

   // for all corners of the target grid
   for (i = 0; i < num_points; ++i) {

      // if the current target corner was not found on the source grid
      if (!found_x[i] || !found_y[i]) continue; // continue with next corner

      // compute the current bucket cell index
      curr_bucket = (unsigned)position_x[i] + (unsigned)position_y[i] *
         bucket_search->num_buckets[0];

      if (curr_bucket >= bucket_search->bucket_to_cell.num_elements)
         continue;

      curr_src_cells = get_dependencies_of_element(
         bucket_search->bucket_to_cell, curr_bucket);
      // for all source cells associated with the current bucket cell
      for (j = 0; j < bucket_search->bucket_to_cell.num_deps_per_element[curr_bucket]; ++j) {

         get_grid_cell(bucket_search->grid_data, curr_src_cells[j],
                       &src_cell);

         struct point tgt_point;
         double tgt_coords[3];

         tgt_point.lon = x_coordinates[i];
         tgt_point.lat = y_coordinates[i];

         LLtoXYZ(tgt_point.lon, tgt_point.lat, tgt_coords);

         // if the target point is within the current source cell
         if (point_in_cell(tgt_point, tgt_coords, src_cell)) {

            num_src_per_tgt[i] = 1;
            tgt_to_src_cell[num_total_dependencies++] = curr_src_cells[j];
            break; // continue with next target point
         }
      } // for j
   } // for i

   tgt_to_src_cell = realloc (tgt_to_src_cell, num_total_dependencies *
                              sizeof(tgt_to_src_cell[0]));

   set_dependencies(tgt_to_src_cells, num_points, num_src_per_tgt, tgt_to_src_cell);
   //-------------
   // free memory
   //-------------

   free_grid_cell(&src_cell);
   free(found_x); free(found_y); free(position_x); free(position_y);
}

static
void bucket_search_do_point_search_p (struct grid_search * search,
                                      struct grid * grid_data,
                                      struct dep_list * target_to_src_points) {


   struct bucket_search * bucket_search = (struct bucket_search *)search;

   grid_search_utils_do_point_search_p(search, bucket_search->grid_data,
                                       grid_data, target_to_src_points);
}

static
void bucket_search_do_point_search_p2 (struct grid_search * search,
                                       double * x_coordinates,
                                       double * y_coordinates, unsigned num_points,
                                       struct dep_list * target_to_src_points) {

   struct bucket_search * bucket_search = (struct bucket_search *)search;

   grid_search_utils_do_point_search_p2(search, bucket_search->grid_data,
                                        x_coordinates, y_coordinates,
                                        num_points, target_to_src_points);
}

static
void bucket_search_do_point_search_p3 (struct grid_search * search,
                                       double * x_coordinates,
                                       double * y_coordinates, unsigned num_points,
                                       struct dep_list * target_to_src_points,
                                       struct points * points) {

   grid_search_utils_do_point_search_p3(search, x_coordinates, y_coordinates,
                                        num_points, target_to_src_points,
                                        points);
}

static
void bucket_search_do_point_search_p4 (struct grid_search * search,
                                       double x_coordinate, double y_coordinate,
                                       unsigned * n_points, unsigned * points_size,
                                       unsigned ** points) {

   struct bucket_search * bucket_search = (struct bucket_search *)search;

   grid_search_utils_do_point_search_p4(search, bucket_search->grid_data,
                                        x_coordinate, y_coordinate, n_points,
                                        points_size, points);
}

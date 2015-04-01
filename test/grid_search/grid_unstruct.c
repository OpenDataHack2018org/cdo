/**
 * @file grid_unstruct.c
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
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
#include <string.h>
#include <math.h>

#include "grid.h"
#include "grid_unstruct.h"
#include "utils.h"
#include "dep_list.h"
#include "geometry.h"
#include "utils.h"
#include "ensure_array_size.h"
#include "sphere_part.h"
#include "bucket_search.h"

// routine declarations

static struct grid * copy_grid_unstruct(struct grid * grid);
static void get_2d_grid_extent_unstruct(struct grid * grid, double (* extent)[2]);
static void get_grid_cell_unstruct(struct grid * grid, unsigned cell_index,
                                   struct grid_cell * cell);
static void get_grid_cell2_unstruct(struct grid * grid, unsigned cell_index,
                                    struct grid_cell * cell,
                                    struct bounding_circle * bnd_circle);
static unsigned get_size_x_coords_unstruct(struct grid * grid);
static unsigned get_size_y_coords_unstruct(struct grid * grid);
static double const * get_x_coords_unstruct(struct grid * grid);
static double const * get_y_coords_unstruct(struct grid * grid);
static void set_x_coords_unstruct(struct grid * grid, double * x_coords);
static void set_y_coords_unstruct(struct grid * grid, double * y_coords);
static unsigned get_size_cell_grid_x_coords_unstruct(struct grid * grid);
static unsigned get_size_cell_grid_y_coords_unstruct(struct grid * grid);
static unsigned get_num_grid_cells_unstruct(struct grid * grid);
static unsigned get_num_grid_corners_unstruct(struct grid * grid);
static unsigned get_num_cell_corners_unstruct(struct grid * grid, unsigned cell_index);
static unsigned get_num_corner_cells_unstruct(struct grid * grid, unsigned corner_index);
static unsigned get_num_grid_edges_unstruct(struct grid * grid);
static unsigned get_num_corner_edges_unstruct(struct grid * grid, unsigned corner_index);
static unsigned get_num_cell_edges_unstruct(struct grid * grid, unsigned cell_index);
static unsigned const * get_corner_edges_unstruct(struct grid * grid, unsigned corner_index);
static unsigned const * get_cell_edge_indices_unstruct(struct grid * grid,
                                                       unsigned cell_index);
static enum yac_edge_type get_edge_type_unstruct(struct grid * grid, unsigned edge_index);
static unsigned const * get_cell_corner_indices_unstruct(struct grid * grid,
                                                         unsigned cell_index);
static unsigned const * get_corner_cell_indices_unstruct(struct grid * grid,
                                                         unsigned corner_index);
static unsigned const * get_cell_x_coord_indices_unstruct(struct grid * grid,
                                                          unsigned cell_index);
static unsigned const * get_cell_y_coord_indices_unstruct(struct grid * grid,
                                                          unsigned cell_index);
static unsigned get_corner_x_coord_index_unstruct(struct grid * grid, unsigned corner_index);
static unsigned get_corner_y_coord_index_unstruct(struct grid * grid, unsigned corner_index);
static double get_corner_x_coord_unstruct(struct grid * grid, unsigned corner_index);
static double get_corner_y_coord_unstruct(struct grid * grid, unsigned corner_index);
static int get_aux_grid_cell_unstruct(struct grid * grid, unsigned corner_index,
                                      unsigned * cell_indices, enum yac_edge_type * edge_type);
static struct dep_list get_cell_neigh_dep_list_unstruct(struct grid * grid);
static void get_boundary_corners_unstruct(struct grid * grid, unsigned * bnd_corners,
                                          unsigned * num_bnd_corners);
static struct grid * generate_cell_grid_unstruct(struct grid * grid,
                                                 double * coordinates_x, 
                                                 double * coordinates_y);
static struct grid * generate_subgrid_unstruct(struct grid * grid,
                                               unsigned * selected_local_cell_ids,
                                               unsigned num_local_cells,
                                               unsigned ** local_cell_ids,
                                               unsigned ** local_corner_ids,
                                               unsigned ** local_edge_ids);
static void pack_grid_unstruct(struct grid * grid, double ** dble_buf,
                               unsigned dble_buf_offset, unsigned * dble_buf_data_size,
                               unsigned * dble_buf_size, unsigned ** uint_buf,
                               unsigned uint_buf_offset, unsigned * uint_buf_data_size,
                               unsigned * uint_buf_size);
static struct grid_search * get_grid_search_unstruct(struct grid * grid);
static void delete_grid_unstruct(struct grid * grid);

static struct grid_vtable unstruct_grid_vtable = {

   .copy                        = copy_grid_unstruct,
   .get_2d_extent               = get_2d_grid_extent_unstruct,
   .get_grid_cell               = get_grid_cell_unstruct,
   .get_grid_cell2              = get_grid_cell2_unstruct,
   .get_size_x_coords           = get_size_x_coords_unstruct,
   .get_size_y_coords           = get_size_y_coords_unstruct,
   .get_x_coords                = get_x_coords_unstruct,
   .get_y_coords                = get_y_coords_unstruct,
   .set_x_coords                = set_x_coords_unstruct,
   .set_y_coords                = set_y_coords_unstruct,
   .get_size_cell_grid_x_coords = get_size_cell_grid_x_coords_unstruct,
   .get_size_cell_grid_y_coords = get_size_cell_grid_y_coords_unstruct,
   .get_num_grid_cells          = get_num_grid_cells_unstruct,
   .get_num_grid_corners        = get_num_grid_corners_unstruct,
   .get_num_cell_corners        = get_num_cell_corners_unstruct,
   .get_num_corner_cells        = get_num_corner_cells_unstruct,
   .get_num_grid_edges          = get_num_grid_edges_unstruct,
   .get_num_corner_edges        = get_num_corner_edges_unstruct,
   .get_num_cell_edges          = get_num_cell_edges_unstruct,
   .get_corner_edges            = get_corner_edges_unstruct,
   .get_cell_edge_indices       = get_cell_edge_indices_unstruct,
   .get_edge_type               = get_edge_type_unstruct,
   .get_cell_corner_indices     = get_cell_corner_indices_unstruct,
   .get_corner_cell_indices     = get_corner_cell_indices_unstruct,
   .get_cell_x_coord_indices    = get_cell_x_coord_indices_unstruct,
   .get_cell_y_coord_indices    = get_cell_y_coord_indices_unstruct,
   .get_corner_x_coord          = get_corner_x_coord_unstruct,
   .get_corner_y_coord          = get_corner_y_coord_unstruct,
   .get_corner_x_coord_index    = get_corner_x_coord_index_unstruct,
   .get_corner_y_coord_index    = get_corner_y_coord_index_unstruct,
   .get_aux_grid_cell           = get_aux_grid_cell_unstruct,
   .get_cell_neigh_dep_list     = get_cell_neigh_dep_list_unstruct,
   .get_boundary_corners        = get_boundary_corners_unstruct,
   .generate_cell_grid          = generate_cell_grid_unstruct,
   .generate_subgrid            = generate_subgrid_unstruct,
   .pack_grid                   = pack_grid_unstruct,
   .get_grid_search             = get_grid_search_unstruct,
   .delete                      = delete_grid_unstruct
};

struct unstruct_grid {

   struct grid_vtable * vtable;

   struct dep_list cell_to_neigh; //!< dependency list containing cell neighbourhood
                                  //!< information (automatically generated once it is
                                  //!< needed)

   double * cell_corners_x,       //!< latitude data
          * cell_corners_y,       //!< longitude data
          * cell_corners_xyz;     //!< 3d coordinates

   unsigned cell_corners_x_by_user; //!< indicates whether cell_corners_x was provided by
                                    //!< the user or automatically generated by a grid
                                    //!< routine
   unsigned cell_corners_y_by_user; //!< indicates whether cell_corners_y was provided by
                                    //!< the user or automatically generated by a grid
                                    //!< routine

   struct dep_list cell_to_vertex;
   struct dep_list vertex_to_cell;
   unsigned cell_to_vertex_by_user; //!< indicates whether cell_to_vertex was provided by
                                    //!< the user or automatically generated by a grid
                                    //!< routine (for example generate_cell_grid)
   unsigned num_vertices;

   struct dep_list corner_to_corner; //!< dependency list that contains for each corner
                                     //!< the other corners that are connected via an edge
                                     //!< it only stores each edge once (for the point with
                                     //!< the lowest id)
   struct dep_list inv_corner_to_corner; //!< inverted version of corner_to_corner

   struct grid_search * grid_search;
};

static unsigned get_num_cell_corners_unstruct (struct grid * grid, unsigned cell_index) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   return unstruct_grid->cell_to_vertex.num_deps_per_element[cell_index];
}

static unsigned get_num_cell_edges_unstruct (struct grid * grid, unsigned cell_index) {

   return get_num_cell_edges_unstruct(grid, cell_index);
}

static unsigned get_num_corner_cells_unstruct (struct grid * grid, unsigned corner_index) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   return unstruct_grid->vertex_to_cell.num_deps_per_element[corner_index];
}

static unsigned const * get_cell_x_coord_indices_unstruct (struct grid * grid,
                                                           unsigned cell_index) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   return yac_get_dependencies_of_element(
      unstruct_grid->cell_to_vertex, cell_index);
}

static unsigned const * get_cell_y_coord_indices_unstruct (struct grid * grid,
                                                           unsigned cell_index) {

   return get_cell_x_coord_indices_unstruct(grid, cell_index);
}

static double get_corner_x_coord_unstruct (struct grid * grid, unsigned corner_index) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   return unstruct_grid->cell_corners_x[corner_index];
}

static double get_corner_y_coord_unstruct (struct grid * grid, unsigned corner_index) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   return unstruct_grid->cell_corners_y[corner_index];
}

static unsigned get_corner_x_coord_index_unstruct (struct grid * grid, unsigned corner_index) {

   return corner_index;
}

static unsigned get_corner_y_coord_index_unstruct (struct grid * grid, unsigned corner_index) {

   return corner_index;
}

static unsigned get_size_x_coords_unstruct(struct grid * grid) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   return unstruct_grid->num_vertices;
}

static unsigned get_size_y_coords_unstruct(struct grid * grid) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   return unstruct_grid->num_vertices;
}

static double const * get_x_coords_unstruct(struct grid * grid) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   return unstruct_grid->cell_corners_x;
}

static double const * get_y_coords_unstruct(struct grid * grid) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   return unstruct_grid->cell_corners_y;
}

static void set_x_coords_unstruct(struct grid * grid, double * x_coords) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   if (!unstruct_grid->cell_corners_x_by_user)
      free(unstruct_grid->cell_corners_x);

   unstruct_grid->cell_corners_x_by_user = 1 == 1;
   unstruct_grid->cell_corners_x = x_coords;
}

static void set_y_coords_unstruct(struct grid * grid, double * y_coords) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   if (!unstruct_grid->cell_corners_y_by_user)
      free(unstruct_grid->cell_corners_y);

   unstruct_grid->cell_corners_y_by_user = 1 == 1;
   unstruct_grid->cell_corners_y = y_coords;
}

static unsigned get_size_cell_grid_x_coords_unstruct(struct grid * grid) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   return unstruct_grid->cell_to_vertex.num_elements;
}

static unsigned get_size_cell_grid_y_coords_unstruct(struct grid * grid) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   return unstruct_grid->cell_to_vertex.num_elements;
}

static unsigned get_num_grid_cells_unstruct(struct grid * grid) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   return unstruct_grid->cell_to_vertex.num_elements;
}

static unsigned get_num_grid_corners_unstruct(struct grid * grid) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   return unstruct_grid->num_vertices;
}

static void generate_unstruct_edges(struct grid * grid) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   unsigned num_corners;
   unsigned num_edges;
   unsigned * num_edges_per_corner;
   unsigned * edge_list;

   yac_init_dep_list(&(unstruct_grid->corner_to_corner));

   num_corners = get_num_grid_corners_unstruct(grid);
   num_edges_per_corner = calloc (num_corners, sizeof (num_edges_per_corner[0]));

   unsigned i, j;
   unsigned const * curr_cell_corners;
   unsigned * temp_num_edges_per_corner;
   unsigned max_num_edges_per_corner;

   temp_num_edges_per_corner = calloc (2 * num_corners, sizeof (temp_num_edges_per_corner[0]));

   //----------------------------------------------
   // count the maximum amount of edges per corner
   //----------------------------------------------

   // for all cells
   for (i = 0; i < unstruct_grid->cell_to_vertex.num_elements; ++i) {

      // get the corners of the current cell
      curr_cell_corners =
         yac_get_dependencies_of_element(unstruct_grid->cell_to_vertex, i);

      // for all corners of the current cell
      for (j = 0; j < unstruct_grid->cell_to_vertex.num_deps_per_element[i] - 1; ++j) {

         if (curr_cell_corners[j] > curr_cell_corners[j+1]) {
            ++temp_num_edges_per_corner[curr_cell_corners[j]];
            ++temp_num_edges_per_corner[curr_cell_corners[j+1] + num_corners];
         } else {
            ++temp_num_edges_per_corner[curr_cell_corners[j+1]];
            ++temp_num_edges_per_corner[curr_cell_corners[j] + num_corners];
         }
      }

      if (curr_cell_corners[0] > curr_cell_corners[j]) {
            ++temp_num_edges_per_corner[curr_cell_corners[0]];
            ++temp_num_edges_per_corner[curr_cell_corners[j] + num_corners];
         } else {
            ++temp_num_edges_per_corner[curr_cell_corners[j]];
            ++temp_num_edges_per_corner[curr_cell_corners[0] + num_corners];
         }
   }

   max_num_edges_per_corner = 0;
   for (i = 0; i < num_corners; ++i) {
   
      if (temp_num_edges_per_corner[i] > max_num_edges_per_corner)
         max_num_edges_per_corner = temp_num_edges_per_corner[i];
      if (temp_num_edges_per_corner[i + num_corners] > max_num_edges_per_corner)
         max_num_edges_per_corner = temp_num_edges_per_corner[i + num_corners];
   }

   //--------------------------------------------------------------
   // get all reverse edges (array is too big, because we do not
   // know yet the actual number of edges per corner)
   //--------------------------------------------------------------
   
   unsigned * temp_edge_list;
   
   temp_edge_list = malloc (num_corners * max_num_edges_per_corner *
                            sizeof (temp_edge_list[0]));
   for (i = 0; i < num_corners; ++i)
      temp_num_edges_per_corner[i] = 0;

   // for all cells
   for (i = 0; i < unstruct_grid->cell_to_vertex.num_elements; ++i) {

      // get the corners of the current cell
      curr_cell_corners =
         yac_get_dependencies_of_element(unstruct_grid->cell_to_vertex, i);

      // for all corners of the current cell
      for (j = 0; j < unstruct_grid->cell_to_vertex.num_deps_per_element[i] - 1; ++j) {

         if (curr_cell_corners[j] < curr_cell_corners[j+1])
            temp_edge_list[curr_cell_corners[j+1] * max_num_edges_per_corner +
                           temp_num_edges_per_corner[curr_cell_corners[j+1]]++] =
               curr_cell_corners[j];
         else
            temp_edge_list[curr_cell_corners[j] * max_num_edges_per_corner +
                           temp_num_edges_per_corner[curr_cell_corners[j]]++] =
               curr_cell_corners[j+1];
      }

      if (curr_cell_corners[0] < curr_cell_corners[j])
         temp_edge_list[curr_cell_corners[j] * max_num_edges_per_corner +
                        temp_num_edges_per_corner[curr_cell_corners[j]]++] =
            curr_cell_corners[0];
      else
         temp_edge_list[curr_cell_corners[0] * max_num_edges_per_corner +
                        temp_num_edges_per_corner[curr_cell_corners[0]]++] =
            curr_cell_corners[j];
   }

   //----------------------
   // get the actual edges
   //----------------------

   num_edges = 0;

   // for all corners (we can skip the first corner)
   for (i = 1; i < num_corners; ++i) {

      // for all reverse edges of the current corner
      for (j = 0; j < temp_num_edges_per_corner[i]; ++j) {

         unsigned curr_corner;

         curr_corner = temp_edge_list[i * max_num_edges_per_corner + j];

         // if the edge is already existing -> get the next one
         if ((curr_corner == i) ||
             (num_edges_per_corner[curr_corner] > 0 &&
              temp_edge_list[curr_corner * max_num_edges_per_corner +
                             num_edges_per_corner[curr_corner] - 1] == i)) continue;

         // store edge
         temp_edge_list[curr_corner * max_num_edges_per_corner +
                        num_edges_per_corner[curr_corner]++] = i;
         ++num_edges;
      }
   }

   free(temp_num_edges_per_corner);

   //-------------------------------------------------------------
   // copy the edges from the temporary edge list to the final one
   //-------------------------------------------------------------

   unsigned * curr_edge, * curr_temp_edge;

   curr_edge = temp_edge_list;
   curr_temp_edge = temp_edge_list;

   for (i = 0; i < num_corners; ++i, curr_temp_edge += max_num_edges_per_corner) {

      for (j = 0; j < num_edges_per_corner[i]; ++j) {

         *curr_edge = curr_temp_edge[j];
         ++curr_edge;
      }
   }

   edge_list = realloc(temp_edge_list, num_edges * sizeof (*edge_list));

   yac_set_dependencies(&(unstruct_grid->corner_to_corner), num_corners,
                        num_edges_per_corner, edge_list);
   yac_invert_dep_list(unstruct_grid->corner_to_corner,
                     &(unstruct_grid->inv_corner_to_corner));
}

// inserts an element into an array and increases the corresponding size
// if the element already exists in the array nothing is done
static void insertion_sort(unsigned element, unsigned * a, unsigned * curr_length) {

   unsigned i, j;

   if (*curr_length > 0) {

      for (i = 0; i < *curr_length; ++i)
         if (a[i] >= element) break;

      // if the element is already in the array
      if ((i != *curr_length) && (a[i] == element)) return;

      // copy new element into array and move bigger elements one position up
      for (j = *curr_length; j > i; --j) {

         a[j] = a[j-1];
      }

      a[i] = element;

      // increase array length indicator
      ++(*curr_length);

   } else {

      a[0] = element;
      *curr_length = 1;
   }
}

static void generate_cell_neigh_dep_unstruct(struct grid * grid) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   struct dep_list vertex_to_cell;

   // generate a vertex to cell mapping
   yac_invert_dep_list(unstruct_grid->cell_to_vertex, &vertex_to_cell);

   unsigned total_num_neighs;
   unsigned num_cells;

   total_num_neighs = 0;
   num_cells = get_num_grid_cells_unstruct(grid);

   unsigned * num_neigh_per_cell;
   unsigned * cell_neigh_dependencies;
   unsigned cell_neigh_dependencies_array_size;

   num_neigh_per_cell = calloc (num_cells, sizeof (num_neigh_per_cell[0]));
   cell_neigh_dependencies = NULL;
   cell_neigh_dependencies_array_size = 0;

   unsigned i, j, k;

   // for all cells
   for (i = 0; i < num_cells; ++i) {

      unsigned curr_num_corners;
      unsigned const * curr_corners;

      curr_num_corners = unstruct_grid->cell_to_vertex.num_deps_per_element[i];
      curr_corners = yac_get_dependencies_of_element(unstruct_grid->cell_to_vertex, i);

      // for all corners of the current cell
      for (j = 0; j < curr_num_corners; ++j) {

         unsigned curr_num_cells;
         unsigned const * curr_cells;

         curr_num_cells = vertex_to_cell.num_deps_per_element[curr_corners[j]];
         curr_cells = yac_get_dependencies_of_element(vertex_to_cell, curr_corners[j]);

         // for all cells associated to this corner
         for (k = 0; k < curr_num_cells; ++k) {

            if (curr_cells[k] == i) continue;

            ENSURE_ARRAY_SIZE(cell_neigh_dependencies,
                              cell_neigh_dependencies_array_size,
                              total_num_neighs + num_neigh_per_cell[i] + 1);

            insertion_sort(curr_cells[k], cell_neigh_dependencies + total_num_neighs,
                           num_neigh_per_cell+i);
         }
      }
      total_num_neighs += num_neigh_per_cell[i];
   }

   cell_neigh_dependencies = realloc (cell_neigh_dependencies, total_num_neighs *
                                      sizeof(cell_neigh_dependencies[0]));

   yac_set_dependencies (&(unstruct_grid->cell_to_neigh), num_cells, num_neigh_per_cell,
                         cell_neigh_dependencies);

   yac_free_dep_list(&vertex_to_cell);
}

static unsigned get_num_grid_edges_unstruct(struct grid * grid) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   if (unstruct_grid->corner_to_corner.num_elements == 0)
      generate_unstruct_edges (grid);

   return yac_get_total_num_dependencies(unstruct_grid->corner_to_corner);
}

static unsigned get_num_corner_edges_unstruct(struct grid * grid, unsigned corner_index) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   if (unstruct_grid->corner_to_corner.num_elements == 0)
      generate_unstruct_edges (grid);

   return unstruct_grid->corner_to_corner.num_deps_per_element[corner_index] +
          unstruct_grid->inv_corner_to_corner.num_deps_per_element[corner_index];
}

static struct dep_list get_cell_neigh_dep_list_unstruct(struct grid * grid) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   if (unstruct_grid->cell_to_neigh.num_elements == 0) generate_cell_neigh_dep_unstruct(grid);
   return unstruct_grid->cell_to_neigh;
}

static void get_boundary_corners_unstruct(struct grid * grid, unsigned * bnd_corners,
                                          unsigned * num_bnd_corners) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   unsigned i, j, k;

   unsigned * num_cells_per_edge; // contains the number cell associated with each edge

   num_cells_per_edge = calloc (yac_get_num_grid_edges(grid),
                                sizeof (num_cells_per_edge[0]));

   //---------------------------------------------
   // each inner edges belongs to two edges
   // outer edges belong to only on local cell ->
   // corners of outer edges are boundary corners
   //---------------------------------------------

   unsigned const * curr_cell_corners;
   
   // for all cells
   for (i = 0; i < yac_get_num_grid_cells(grid); ++i) {

      curr_cell_corners = yac_get_dependencies_of_element(unstruct_grid->cell_to_vertex, i);

      k = unstruct_grid->cell_to_vertex.num_deps_per_element[i] - 1;
      // for all corners of the current cell
      for (j = 0; j < unstruct_grid->cell_to_vertex.num_deps_per_element[i]; k = j++) {

         if (curr_cell_corners[k] < curr_cell_corners[j]) {

            ++num_cells_per_edge[yac_get_dependency_index (unstruct_grid->corner_to_corner,
                                                           curr_cell_corners[k], curr_cell_corners[j])];
         } else {

            ++num_cells_per_edge[yac_get_dependency_index (unstruct_grid->corner_to_corner,
                                                           curr_cell_corners[j], curr_cell_corners[k])];
         }
      }
   }

   unsigned temp_num_bnd_corners, corner_a, corner_b;

   temp_num_bnd_corners = 0;
   // for all edges
   for (i = 0; i < yac_get_num_grid_edges(grid); ++i) {

      if (num_cells_per_edge[i] == 1) {

         yac_get_dependency(unstruct_grid->corner_to_corner, i,
                            &corner_a, &corner_b);

         insertion_sort(corner_a, bnd_corners, &temp_num_bnd_corners);
         insertion_sort(corner_b, bnd_corners, &temp_num_bnd_corners);
      }
   }

   *num_bnd_corners = temp_num_bnd_corners;

   free (num_cells_per_edge);
}

static void get_grid_cell_unstruct(struct grid * grid, unsigned cell_index,
                                   struct grid_cell * cell) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   unsigned num_corners;

   num_corners = get_num_cell_corners_unstruct(grid, cell_index);

   // if the memory for the coordinates needs to be reallocated
   if (num_corners > cell->array_size) {
      cell->coordinates_x = realloc (cell->coordinates_x, num_corners * sizeof(cell->coordinates_x[0]));
      cell->coordinates_y = realloc (cell->coordinates_y, num_corners * sizeof(cell->coordinates_y[0]));
      cell->coordinates_xyz = realloc (cell->coordinates_xyz, 3 * num_corners * sizeof(cell->coordinates_xyz[0]));
      cell->edge_type = realloc (cell->edge_type, num_corners * sizeof(cell->edge_type[0]));
      cell->array_size = num_corners;
   }

   cell->num_corners = num_corners;

   // set the corners for of the cell
   unsigned i;
   unsigned const * curr_vertex;

   curr_vertex = yac_get_dependencies_of_element (
      unstruct_grid->cell_to_vertex, cell_index);

   for (i = 0; i < cell->num_corners; ++i) {

      cell->coordinates_x[i] = unstruct_grid->cell_corners_x[curr_vertex[i]];
      cell->coordinates_y[i] = unstruct_grid->cell_corners_y[curr_vertex[i]];
      cell->coordinates_xyz[0+3*i] = unstruct_grid->cell_corners_xyz[0+3*curr_vertex[i]];
      cell->coordinates_xyz[1+3*i] = unstruct_grid->cell_corners_xyz[1+3*curr_vertex[i]];
      cell->coordinates_xyz[2+3*i] = unstruct_grid->cell_corners_xyz[2+3*curr_vertex[i]];
   }

   // set the edge type for the cell
   for (i = 0; i < cell->num_corners; ++i)
         cell->edge_type[i] = GREAT_CIRCLE;
}

static void get_grid_cell2_unstruct(struct grid * grid, unsigned cell_index,
                                    struct grid_cell * cell,
                                    struct bounding_circle * bnd_circle) {

   get_grid_cell_unstruct(grid, cell_index, cell);
   if (cell->num_corners == 3) {
      yac_get_cell_bounding_circle_unstruct_triangle(cell->coordinates_xyz + 0*3,
                                                     cell->coordinates_xyz + 1*3,
                                                     cell->coordinates_xyz + 2*3,
                                                     bnd_circle);
   } else
      yac_get_cell_bounding_circle(*cell, bnd_circle);
}

static unsigned const * get_corner_edges_unstruct(struct grid * grid, unsigned corner_index) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   static unsigned * corners = NULL;
   static unsigned size_corners = 0;

   if (unstruct_grid->corner_to_corner.num_elements == 0)
      generate_unstruct_edges (grid);

   unsigned num_edges = unstruct_grid->corner_to_corner.num_deps_per_element[corner_index] +
                        unstruct_grid->inv_corner_to_corner.num_deps_per_element[corner_index];

   if (num_edges > size_corners) {

      corners = realloc (corners, num_edges * sizeof (corners[0]));
      size_corners = num_edges;
   }

   unsigned i, j;
   unsigned const * curr_edges;
   unsigned curr_num_edges;

   curr_num_edges = unstruct_grid->inv_corner_to_corner.num_deps_per_element[corner_index];
   curr_edges = yac_get_dependencies_of_element(unstruct_grid->inv_corner_to_corner, corner_index);

   for (i = 0; i < curr_num_edges; ++i)
      corners[i] = curr_edges[i];

   curr_num_edges = unstruct_grid->corner_to_corner.num_deps_per_element[corner_index];
   curr_edges = yac_get_dependencies_of_element(unstruct_grid->corner_to_corner, corner_index);

   for (j = 0; j < curr_num_edges; ++j)
      corners[i+j] = curr_edges[j];

   return corners;
}

static unsigned const * get_cell_corner_indices_unstruct(struct grid * grid,
                                                         unsigned cell_index) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   return yac_get_dependencies_of_element(
      unstruct_grid->cell_to_vertex, cell_index);
}

static unsigned const * get_corner_cell_indices_unstruct(struct grid * grid, unsigned corner_index) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   return yac_get_dependencies_of_element(
      unstruct_grid->vertex_to_cell, corner_index);
}

static unsigned const * get_cell_edge_indices_unstruct(struct grid * grid, unsigned cell_index) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   static unsigned * edges = NULL;
   static unsigned size_edges = 0;

   if (unstruct_grid->corner_to_corner.num_elements == 0)
      generate_unstruct_edges (grid);

   unsigned num_edges;

   num_edges = get_num_cell_corners_unstruct(grid, cell_index);

   if (size_edges < num_edges) {

      edges = realloc (edges, num_edges * sizeof (edges[0]));
      size_edges = num_edges;
   }

   unsigned const * curr_cell_corners;

   curr_cell_corners = get_cell_corner_indices_unstruct(grid, cell_index);

   unsigned i;

   for (i = 0; i < num_edges-1; ++i) {

      if (curr_cell_corners[i] > curr_cell_corners[i+1])
         edges[i] = yac_get_dependency_index(unstruct_grid->corner_to_corner,
                                             curr_cell_corners[i+1],
                                             curr_cell_corners[i]);
      else
         edges[i] = yac_get_dependency_index(unstruct_grid->corner_to_corner,
                                             curr_cell_corners[i],
                                             curr_cell_corners[i+1]);
   }

   if (curr_cell_corners[num_edges-1] > curr_cell_corners[0])
      edges[num_edges-1] = yac_get_dependency_index(unstruct_grid->corner_to_corner,
                                                    curr_cell_corners[0],
                                                    curr_cell_corners[num_edges-1]);
   else
      edges[num_edges-1] = yac_get_dependency_index(unstruct_grid->corner_to_corner,
                                                    curr_cell_corners[num_edges-1],
                                                    curr_cell_corners[0]);

   return edges;
}

static enum yac_edge_type get_edge_type_unstruct(struct grid * grid, unsigned edge_index) {

   return GREAT_CIRCLE;
}

static void get_2d_grid_extent_unstruct(struct grid * grid, double (* extent)[2]) {

   double const tol = 1.0e-12;
   struct bounding_circle circle;

   yac_get_grid_bounding_circle(grid, &circle);

   // check if the grid covers the whole sphere
   if (circle.inc_angle - tol >= M_PI) {

      extent[0][0] =  -M_PI_2;
      extent[0][1] =   M_PI_2;
      extent[1][0] = -M_PI;
      extent[1][1] =  M_PI;

      return;
   }

   double rad_inc_angle;
   double base_point[2];

   rad_inc_angle = circle.inc_angle + tol;
   base_point[0] = circle.base_point[0] + tol;
   base_point[1] = circle.base_point[1] + tol;

   extent[1][0] = MAX(-M_PI_2, base_point[1] - rad_inc_angle);
   extent[1][1] = MIN( M_PI_2, base_point[1] + rad_inc_angle);

   // check if the circle covers a pole
   if (fabs(extent[1][0] + M_PI_2) <= tol ||
       fabs(extent[1][1] - M_PI_2) <= tol) {

      extent[0][0] = -M_PI;
      extent[0][1] =  M_PI;
   } else {
      extent[0][0] = base_point[0] - rad_inc_angle;
      extent[0][1] = base_point[0] + rad_inc_angle;
   }
}

static struct grid * generate_cell_grid_unstruct(struct grid * grid,
                                                 double * coordinates_x, 
                                                 double * coordinates_y) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   // construct corner to cell dependency list
   struct dep_list corner_to_cell;

   yac_init_dep_list(&corner_to_cell);
   yac_invert_dep_list(unstruct_grid->cell_to_vertex, &corner_to_cell);

   // generate neighbour dependency list if it does no yet exist
   // struct dep_list * cell_to_neigh;
   if (unstruct_grid->cell_to_neigh.num_elements == 0) generate_cell_neigh_dep_unstruct(grid);
   // cell_to_neigh = &(unstruct_grid->cell_to_neigh);

   unsigned i, j, k;

   unsigned * num_corners_per_cell;
   unsigned num_total_cells;
   unsigned * cell_corner_dependencies;
   unsigned size_cell_corner_dependencies;
   unsigned num_total_cell_corner_dependencies;
   unsigned const * curr_adjacent_cells;
   unsigned num_curr_adjacent_cells;
   unsigned const * adjacent_corners;
   unsigned num_adjacent_corners;
   unsigned * edge_to_cell;
   unsigned * num_cell_per_edge;
   unsigned size_num_cell_per_edge;

   size_num_cell_per_edge = 0;
   edge_to_cell = NULL;
   num_cell_per_edge = NULL;

   num_corners_per_cell = malloc (yac_get_num_grid_corners(grid) *
                                  sizeof (num_corners_per_cell[0]));
   num_total_cells = 0;
   cell_corner_dependencies = NULL;
   size_cell_corner_dependencies = 0;
   num_total_cell_corner_dependencies = 0;

   // for all corners of the original grid (possible centres of cells of the cell grid)
   for (i = 0; i < get_num_grid_corners_unstruct(grid); ++i) {

      // check whether the current corner is the centre of a cell of the cell grid

      // if the current corner has not enough adjacent cell to form the centre
      // of a new cell of the cell grid
      if (corner_to_cell.num_deps_per_element[i] <= 2)
         continue;

      // gets the corners that are directly connected to the current corner (i) by
      // an edge
      num_adjacent_corners = get_num_corner_edges_unstruct(grid, i);
      adjacent_corners = get_corner_edges_unstruct(grid, i);

      // check whether num_cell_per_edge is big enough
      if (size_num_cell_per_edge < num_adjacent_corners) {

         size_num_cell_per_edge = num_adjacent_corners;
         free(num_cell_per_edge);
         free(edge_to_cell);
         num_cell_per_edge = malloc (size_num_cell_per_edge * sizeof (num_cell_per_edge[0]));
         edge_to_cell = malloc (2 * size_num_cell_per_edge * sizeof (edge_to_cell[0]));
      }

      // initialise num_cell_per_edge
      for (j = 0; j < num_adjacent_corners; ++j)
         num_cell_per_edge[j] = 0;

      // gets all cells that have the current corner (i) on their boundary definition
      num_curr_adjacent_cells = corner_to_cell.num_deps_per_element[i];
      curr_adjacent_cells = yac_get_dependencies_of_element(corner_to_cell, i);
         
      // for all cell adjacent to the current corner
      for (j = 0; j < num_curr_adjacent_cells; ++j) {

         unsigned const * curr_corners_of_cell;
         unsigned num_curr_corners_of_cell;

         // gets the corners of the current cell (curr_adjacent_cells[j])
         num_curr_corners_of_cell = unstruct_grid->cell_to_vertex.num_deps_per_element[curr_adjacent_cells[j]];
         curr_corners_of_cell = yac_get_dependencies_of_element (
            unstruct_grid->cell_to_vertex, curr_adjacent_cells[j]);

         // for all corners of the current cell
         // find the current corner (i) in the boundary definition of the
         // current cell (curr_adjacent_cells[j])
         for (k = 0; k < num_curr_corners_of_cell; ++k)
            if (curr_corners_of_cell[k] == i) break;

         unsigned curr_adjacent_corners[2];

         // get the two associated edges
         if (k == 0) {
            curr_adjacent_corners[0] = curr_corners_of_cell[num_curr_corners_of_cell-1];
            curr_adjacent_corners[1] = curr_corners_of_cell[k+1];
         } else if (k == num_curr_corners_of_cell - 1){
            curr_adjacent_corners[0] = curr_corners_of_cell[k-1];
            curr_adjacent_corners[1] = curr_corners_of_cell[0];
         } else {
            curr_adjacent_corners[0] = curr_corners_of_cell[k-1];
            curr_adjacent_corners[1] = curr_corners_of_cell[k+1];
         }

         // for all adjacent corners of the current corner (i)
         for (k = 0; k < num_adjacent_corners; ++k) {

            if (adjacent_corners[k] == curr_adjacent_corners[0]) {

               edge_to_cell[k + num_cell_per_edge[k] * num_adjacent_corners] =
                  curr_adjacent_cells[j];
               num_cell_per_edge[k]++;

            } else if (adjacent_corners[k] == curr_adjacent_corners[1]) {

               edge_to_cell[k + num_cell_per_edge[k] * num_adjacent_corners] =
                  curr_adjacent_cells[j];
               num_cell_per_edge[k]++;
            }
         }
      } // (j = 0; j < num_curr_adjacent_cells; ++j)

      // check whether we found two cells for each edge
      for (j = 0; j < num_adjacent_corners; ++j)
         if (num_cell_per_edge[j] != 2) break;

      // if the current corner is the centre of a cell grid cell
      if (j == num_adjacent_corners) {

         // check whether we have not enough space to store the new cell
         if (size_cell_corner_dependencies < num_total_cell_corner_dependencies +
             num_adjacent_corners) {

            size_cell_corner_dependencies = num_total_cell_corner_dependencies +
                                            50 * num_adjacent_corners;
            cell_corner_dependencies = realloc (cell_corner_dependencies,
                                                size_cell_corner_dependencies *
                                                sizeof (cell_corner_dependencies[0]));
         }

         unsigned curr_corner;
         unsigned curr_cell;

         curr_corner = adjacent_corners[0];
         curr_cell = edge_to_cell[0];

         cell_corner_dependencies[num_total_cell_corner_dependencies] = curr_cell;

         // bring the corners into the correct order to form a cell
         for (j = 1; j < num_adjacent_corners; ++ j) {

            for (k = 0; k < num_adjacent_corners; ++k)
               if ((edge_to_cell[k] == curr_cell) && (adjacent_corners[k] != curr_corner)) {

                  curr_corner = adjacent_corners[k];
                  curr_cell = edge_to_cell[k + num_adjacent_corners];
                  break;
               }

            if (k == num_adjacent_corners) {
               for (; k < 2*num_adjacent_corners; ++k)
                  if ((edge_to_cell[k] == curr_cell) && (k != curr_corner +
                                                         num_adjacent_corners)) {
                        k -= num_adjacent_corners;
                        curr_corner = adjacent_corners[k];
                        curr_cell = edge_to_cell[k];
                        break;
                     }
            }

            cell_corner_dependencies[num_total_cell_corner_dependencies + j] = curr_cell;
         }

         num_corners_per_cell[num_total_cells++] = num_adjacent_corners;
         num_total_cell_corner_dependencies += num_adjacent_corners;
      }
   } // (i = 0; i < get_num_grid_corners(*grid); ++i)

   free(num_cell_per_edge);
   free(edge_to_cell);

   struct dep_list cell_to_vertex;

   num_corners_per_cell = realloc(num_corners_per_cell, num_total_cells *
                                  sizeof (num_corners_per_cell[0]));
   cell_corner_dependencies = realloc(cell_corner_dependencies,
                                      num_total_cell_corner_dependencies *
                                      sizeof (cell_corner_dependencies[0]));

   yac_init_dep_list(&cell_to_vertex);
   yac_set_dependencies(&cell_to_vertex, num_total_cells, num_corners_per_cell,
                        cell_corner_dependencies);

   struct unstruct_grid * cell_grid = (struct unstruct_grid *)yac_unstruct_grid_new(
      coordinates_x, coordinates_y, get_num_grid_cells_unstruct(grid), cell_to_vertex);

   cell_grid->cell_to_vertex_by_user = 1 == 0;

   yac_free_dep_list(&corner_to_cell);

   return (struct grid *)cell_grid;
}

static int get_aux_grid_cell_unstruct(struct grid * grid, unsigned corner_index,
                                      unsigned * cell_indices, enum yac_edge_type * edge_type) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   unsigned j, k;

   unsigned const * adjacent_cells;
   unsigned num_adjacent_cells;
   unsigned const * adjacent_corners;
   unsigned num_adjacent_corners;
   unsigned * edge_to_cell;
   unsigned * num_cell_per_edge;

   // if the current corner has not enough adjacent cell to form the centre
   // of a new cell of the cell grid
   if (unstruct_grid->vertex_to_cell.num_deps_per_element[corner_index] <= 2)
      return 1 == 0;

   // gets the corners that are directly connected to the given corner by
   // an edge
   num_adjacent_corners = get_num_corner_edges_unstruct(grid, corner_index);
   adjacent_corners = get_corner_edges_unstruct(grid, corner_index);

   num_cell_per_edge = calloc (num_adjacent_corners, sizeof (num_cell_per_edge[0]));
   edge_to_cell = malloc (2 * num_adjacent_corners * sizeof (edge_to_cell[0]));

   // gets all cells that have the given corner on their boundary definition
   num_adjacent_cells = unstruct_grid->vertex_to_cell.num_deps_per_element[corner_index];
   adjacent_cells = yac_get_dependencies_of_element(unstruct_grid->vertex_to_cell, corner_index);
      
   // for all cell adjacent to the given corner
   for (j = 0; j < num_adjacent_cells; ++j) {

      unsigned const * curr_corners_of_cell;
      unsigned num_curr_corners_of_cell;

      // gets the corners of the current cell (adjacent_cells[j])
      num_curr_corners_of_cell = unstruct_grid->cell_to_vertex.num_deps_per_element[adjacent_cells[j]];
      curr_corners_of_cell = yac_get_dependencies_of_element (
         unstruct_grid->cell_to_vertex, adjacent_cells[j]);

      // for all corners of the current cell
      // find the given corner in the boundary definition of the
      // current cell (adjacent_cells[j])
      for (k = 0; k < num_curr_corners_of_cell; ++k)
         if (curr_corners_of_cell[k] == corner_index) break;

      unsigned curr_adjacent_corners[2];

      // get the two associated edges
      if (k == 0) {
         curr_adjacent_corners[0] = curr_corners_of_cell[num_curr_corners_of_cell-1];
         curr_adjacent_corners[1] = curr_corners_of_cell[k+1];
      } else if (k == num_curr_corners_of_cell - 1){
         curr_adjacent_corners[0] = curr_corners_of_cell[k-1];
         curr_adjacent_corners[1] = curr_corners_of_cell[0];
      } else {
         curr_adjacent_corners[0] = curr_corners_of_cell[k-1];
         curr_adjacent_corners[1] = curr_corners_of_cell[k+1];
      }

      // for all adjacent corners of the given corner
      for (k = 0; k < num_adjacent_corners; ++k) {

         if (adjacent_corners[k] == curr_adjacent_corners[0]) {

            edge_to_cell[k + num_cell_per_edge[k] * num_adjacent_corners] =
               adjacent_cells[j];
            num_cell_per_edge[k]++;

         } else if (adjacent_corners[k] == curr_adjacent_corners[1]) {

            edge_to_cell[k + num_cell_per_edge[k] * num_adjacent_corners] =
               adjacent_cells[j];
            num_cell_per_edge[k]++;
         }
      }
   } // (j = 0; j < num_curr_adjacent_cells; ++j)

   // check whether we found two cells for each edge
   for (j = 0; j < num_adjacent_corners; ++j)
      if (num_cell_per_edge[j] != 2) break;

   free(num_cell_per_edge);

   // if the current corner is the centre of a cell grid cell
   if (j == num_adjacent_corners) {

      unsigned curr_corner;
      unsigned curr_cell;

      curr_corner = adjacent_corners[0];
      curr_cell = edge_to_cell[0];

      cell_indices[0] = curr_cell;

      // bring the corners into the correct order to form a cell
      for (j = 1; j < num_adjacent_corners; ++ j) {

         for (k = 0; k < num_adjacent_corners; ++k)
            if ((edge_to_cell[k] == curr_cell) && (adjacent_corners[k] != curr_corner)) {

               curr_corner = adjacent_corners[k];
               curr_cell = edge_to_cell[k + num_adjacent_corners];
               break;
            }

         if (k == num_adjacent_corners) {
            for (; k < 2*num_adjacent_corners; ++k)
               if ((edge_to_cell[k] == curr_cell) && (k != curr_corner +
                                                      num_adjacent_corners)) {
                     k -= num_adjacent_corners;
                     curr_corner = adjacent_corners[k];
                     curr_cell = edge_to_cell[k];
                     break;
                  }
         }

         cell_indices[j] = curr_cell;
      }

      for (j = 0; j < num_adjacent_corners; ++j)
         edge_type[j] = GREAT_CIRCLE;
   }

   free(edge_to_cell);

   return j == num_adjacent_corners;
}

static struct grid * copy_grid_unstruct(struct grid * grid) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   struct unstruct_grid * copy;
   double * coordinates_x, * coordinates_y;
   unsigned x_coords_size = get_size_x_coords_unstruct(grid) * sizeof(*coordinates_x);
   unsigned y_coords_size = get_size_y_coords_unstruct(grid) * sizeof(*coordinates_y);
   unsigned num_vertices;
   struct dep_list cell_to_vertex;

   coordinates_x = malloc(x_coords_size);
   coordinates_y = malloc(y_coords_size);
   memcpy(coordinates_x, unstruct_grid->cell_corners_x, x_coords_size);
   memcpy(coordinates_y, unstruct_grid->cell_corners_y, y_coords_size);

   num_vertices = get_num_grid_corners_unstruct(grid);

   yac_init_dep_list(&cell_to_vertex);
   yac_copy_dep_list(unstruct_grid->cell_to_vertex, &cell_to_vertex);

   copy = (struct unstruct_grid *)yac_unstruct_grid_new(coordinates_x, coordinates_y,
                                                        num_vertices, cell_to_vertex);

   copy->cell_to_vertex_by_user = 1 == 0;
   copy->cell_corners_x_by_user = 1 == 0;
   copy->cell_corners_y_by_user = 1 == 0;

   return (struct grid *)copy;
}

static void pack_grid_unstruct(struct grid * grid, double ** dble_buf,
                               unsigned dble_buf_offset, unsigned * dble_buf_data_size,
                               unsigned * dble_buf_size, unsigned ** uint_buf,
                               unsigned uint_buf_offset, unsigned * uint_buf_data_size,
                               unsigned * uint_buf_size) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   *dble_buf_data_size = 2 * unstruct_grid->num_vertices;

   ENSURE_ARRAY_SIZE(*dble_buf, *dble_buf_size, dble_buf_offset+*dble_buf_data_size);

   memcpy((*dble_buf)+dble_buf_offset, unstruct_grid->cell_corners_x,
          unstruct_grid->num_vertices * sizeof(**dble_buf));
   memcpy((*dble_buf)+dble_buf_offset+unstruct_grid->num_vertices, unstruct_grid->cell_corners_y,
          unstruct_grid->num_vertices * sizeof(**dble_buf));

   ENSURE_ARRAY_SIZE(*uint_buf, *uint_buf_size, uint_buf_offset+2);

   (*uint_buf)[uint_buf_offset+0] = yac_hash("UNSTRUCT");
   (*uint_buf)[uint_buf_offset+1] = unstruct_grid->num_vertices;

   yac_pack_dep_list(unstruct_grid->cell_to_vertex, uint_buf,
                     uint_buf_offset+2, uint_buf_data_size, uint_buf_size);

   *uint_buf_data_size += 2;
}

struct grid * yac_unpack_unstruct_grid(double * dble_buf, unsigned * dble_buf_data_size,
                                       unsigned * uint_buf, unsigned * uint_buf_data_size) {

   unsigned num_vertices = uint_buf[1];

   double * cell_corners_x, * cell_corners_y;

   cell_corners_x = malloc(num_vertices * sizeof(*cell_corners_x));
   cell_corners_y = malloc(num_vertices * sizeof(*cell_corners_y));

   memcpy(cell_corners_x, dble_buf, num_vertices * sizeof(*cell_corners_x));
   memcpy(cell_corners_y, dble_buf+num_vertices, num_vertices * sizeof(*cell_corners_y));

   *dble_buf_data_size = 2 * num_vertices;

   struct dep_list cell_to_vertex;

   yac_unpack_dep_list(&cell_to_vertex, uint_buf+2, uint_buf_data_size);

   *uint_buf_data_size += 2;

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)yac_unstruct_grid_new(cell_corners_x,
                                                                 cell_corners_y,
                                                                 num_vertices,
                                                                 cell_to_vertex);

   unstruct_grid->cell_to_vertex_by_user = 1 == 0;
   unstruct_grid->cell_corners_x_by_user = 1 == 0;
   unstruct_grid->cell_corners_y_by_user = 1 == 0;

   return (struct grid *)unstruct_grid;
}

static void delete_grid_unstruct(struct grid * grid) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   yac_free_dep_list(&(unstruct_grid->cell_to_neigh));
   if (!unstruct_grid->cell_to_vertex_by_user)
      yac_free_dep_list(&(unstruct_grid->cell_to_vertex));
   yac_free_dep_list(&(unstruct_grid->vertex_to_cell));
   yac_free_dep_list(&(unstruct_grid->corner_to_corner));
   yac_free_dep_list(&(unstruct_grid->inv_corner_to_corner));

   if (!unstruct_grid->cell_corners_x_by_user)
      free(unstruct_grid->cell_corners_x);
   if (!unstruct_grid->cell_corners_y_by_user)
      free(unstruct_grid->cell_corners_y);
   free(unstruct_grid->cell_corners_xyz);

   if (unstruct_grid->grid_search != NULL)
      yac_delete_grid_search(unstruct_grid->grid_search);

   free(unstruct_grid);
}

static struct grid * generate_subgrid_unstruct(struct grid * grid,
                                               unsigned * selected_local_cell_ids,
                                               unsigned num_local_cells,
                                               unsigned ** local_cell_ids,
                                               unsigned ** local_corner_ids,
                                               unsigned ** local_edge_ids) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   unsigned i, j, k;

   unsigned * subgrid_corners;
   unsigned num_subgrid_corners;

   unsigned * num_corners_per_cell;
   unsigned total_num_cell_corners;

   //--------------------------------------------------------
   // generate a list of all corners required for the subgrid
   //--------------------------------------------------------

// NBITS is the number of bits in an unsigned
#define NBITS (sizeof(unsigned)*8)

   // the mask contains a bit (initially set to zero) for each corner of the grid
   unsigned * corner_mask;
   unsigned num_grid_corners;

   num_grid_corners = yac_get_num_grid_corners(grid);

   // an unsigned normally has 32 bit -> for a grid with 128 corners we need a
   // corner_mask with an array size of 4
   corner_mask = calloc(((num_grid_corners + NBITS - 1) / NBITS ), sizeof(*corner_mask));

   num_corners_per_cell = malloc(num_local_cells * sizeof(*num_corners_per_cell));
   total_num_cell_corners = 0;
   num_subgrid_corners = 0;

   // for all selected cells
   for (i = 0; i < num_local_cells; ++i) {

      unsigned const * corner_ids;
      unsigned num_corners;

      // get the ids of corners of the current cell
      corner_ids = get_cell_corner_indices_unstruct(grid, selected_local_cell_ids[i]);
      num_corners = get_num_cell_corners_unstruct(grid, selected_local_cell_ids[i]);

      num_corners_per_cell[i] = num_corners;
      total_num_cell_corners += num_corners;

      // for all the corners of the current cell
      for (j = 0; j < num_corners; ++j) {

         // if the respective bit in corner_mask for the current corner is not yet set
         if (!(corner_mask[corner_ids[j] / NBITS ] & (1 << (corner_ids[j] & ( NBITS - 1))))) {

            // set the respective bit
            corner_mask[corner_ids[j] / NBITS ] |= 1 << (corner_ids[j] & ( NBITS - 1));
            num_subgrid_corners++;
         }
      }
   }

   subgrid_corners = malloc(num_subgrid_corners * sizeof(*subgrid_corners));

   num_subgrid_corners = 0;

   // read the data from corner_mask and convert the set bits into corner ids
   for (i = 0; i < (num_grid_corners + NBITS - 1) / NBITS; ++i)
      for (j = 0; j < NBITS; ++j)
         if (corner_mask[i] & 1 << j)
            subgrid_corners[num_subgrid_corners++] = i * NBITS + j;

   free(corner_mask);

#undef NBITS

   //----------------------------------------------------------
   // generate a cell to vertex dependency list for the subgrid
   //----------------------------------------------------------

   unsigned curr_index;

   curr_index = 0;

   unsigned * cell_to_corner_dep;

   cell_to_corner_dep = malloc(total_num_cell_corners * sizeof(*cell_to_corner_dep));

   total_num_cell_corners = 0;

   for (i = 0; i < num_local_cells; ++i) {

      unsigned const * corner_ids;
      unsigned num_corners;

      corner_ids = get_cell_corner_indices_unstruct(grid, selected_local_cell_ids[i]);
      num_corners = get_num_cell_corners_unstruct(grid, selected_local_cell_ids[i]);

      for (j = 0; j < num_corners; ++j) {

         while (subgrid_corners[curr_index] > corner_ids[j]) --curr_index;
         while (subgrid_corners[curr_index] < corner_ids[j]) ++curr_index;

         cell_to_corner_dep[total_num_cell_corners+j] = curr_index;
      }
      total_num_cell_corners += num_corners;
   }

   struct dep_list cell_to_corner;

   yac_set_dependencies(&cell_to_corner, num_local_cells, num_corners_per_cell,
                        cell_to_corner_dep);

   // get the corner coordinates of the subgrid

   double * cell_corners_x, * cell_corners_y;

   cell_corners_x = malloc(num_subgrid_corners * sizeof(*cell_corners_x));
   cell_corners_y = malloc(num_subgrid_corners * sizeof(*cell_corners_y));

   for (i = 0; i < num_subgrid_corners; ++i) {

      cell_corners_x[i] = unstruct_grid->cell_corners_x[subgrid_corners[i]];
      cell_corners_y[i] = unstruct_grid->cell_corners_y[subgrid_corners[i]];
   }

   // generate subgrid

   struct unstruct_grid * subgrid;

   subgrid = (struct unstruct_grid *)yac_unstruct_grid_new(cell_corners_x, cell_corners_y,
                                                           num_subgrid_corners,
                                                           cell_to_corner);

   subgrid->cell_corners_x_by_user = 1 == 0;
   subgrid->cell_corners_y_by_user = 1 == 0;
   subgrid->cell_to_vertex_by_user = 1 == 0;

   // set the local_cell_ids array

   *local_cell_ids = malloc(num_local_cells * sizeof(**local_cell_ids));

   memcpy(*local_cell_ids, selected_local_cell_ids, num_local_cells * sizeof(**local_cell_ids));

   // set the local_corner_ids array

   *local_corner_ids = subgrid_corners;

   // set the local_edge_ids array

   unsigned num_subgrid_edges;

   num_subgrid_edges = get_num_grid_edges_unstruct((struct grid*)subgrid);

   *local_edge_ids = malloc(num_subgrid_edges * sizeof(**local_edge_ids));

   k = 0;

   if (unstruct_grid->corner_to_corner.num_elements == 0)
      generate_unstruct_edges (grid);

   for (i = 0; i < num_subgrid_corners; ++i) {

      unsigned const * subgrid_edges;
      unsigned num_curr_subgrid_edges;

      subgrid_edges = yac_get_dependencies_of_element(subgrid->corner_to_corner, i);
      num_curr_subgrid_edges = subgrid->corner_to_corner.num_deps_per_element[i];

      for (j = 0; j < num_curr_subgrid_edges; ++j) {

         (*local_edge_ids)[k++] = yac_get_dependency_index(unstruct_grid->corner_to_corner,
                                                           (*local_corner_ids)[i],
                                                           (*local_corner_ids)[subgrid_edges[j]]);
      }
   }

   return (struct grid *)subgrid;
}

struct grid * yac_unstruct_grid_new(double * coordinates_x, double * coordinates_y,
                                    unsigned num_vertices, struct dep_list cell_to_vertex) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = malloc(1 * sizeof(*unstruct_grid));

   unstruct_grid->vtable = &unstruct_grid_vtable;

   yac_init_dep_list(&(unstruct_grid->cell_to_neigh));

   unstruct_grid->cell_corners_x = coordinates_x;
   unstruct_grid->cell_corners_y = coordinates_y;
   unstruct_grid->cell_corners_xyz =
      malloc(3 * num_vertices * sizeof(*(unstruct_grid->cell_corners_xyz)));

   for (unsigned i = 0; i < num_vertices; ++i)
      LLtoXYZ(unstruct_grid->cell_corners_x[i],
              unstruct_grid->cell_corners_y[i],
              unstruct_grid->cell_corners_xyz + i * 3);

   unstruct_grid->cell_to_vertex = cell_to_vertex;
   yac_invert_dep_list(cell_to_vertex, &(unstruct_grid->vertex_to_cell));
   unstruct_grid->num_vertices = num_vertices;

   yac_init_dep_list(&(unstruct_grid->corner_to_corner));
   yac_init_dep_list(&(unstruct_grid->inv_corner_to_corner));

   unstruct_grid->cell_to_vertex_by_user = 1 == 1;
   unstruct_grid->cell_corners_x_by_user = 1 == 1;
   unstruct_grid->cell_corners_y_by_user = 1 == 1;

   unstruct_grid->grid_search = NULL;

   return (struct grid *)unstruct_grid;
}

static struct grid_search * get_grid_search_unstruct(struct grid * grid) {

   struct unstruct_grid * unstruct_grid;

   unstruct_grid = (struct unstruct_grid *)grid;

   if (unstruct_grid->grid_search == NULL)
      unstruct_grid->grid_search = yac_sphere_part_search_new(grid);

   return unstruct_grid->grid_search;
}

/**
 * @file grid_reg2d.c
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
#include <string.h>
#include <math.h>

#include "grid.h"
#include "grid_reg2d.h"
#include "grid_reg2d_common.h"
#include "utils.h"
#include "dep_list.h"
#include "geometry.h"
#include "ensure_array_size.h"
#include "sphere_part.h"

// routine declarations

static struct grid * copy_grid_reg2d(struct grid * grid);
static void get_2d_grid_extent_reg2d(struct grid * grid, double (* extent)[2]);
static void get_grid_cell_reg2d(struct grid * grid, unsigned cell_index,
                                struct grid_cell * cell);
static void get_grid_cell2_reg2d(struct grid * grid, unsigned cell_index,
                                 struct grid_cell * cell,
                                 struct bounding_circle * bnd_circle);
static unsigned get_size_x_coords_reg2d(struct grid * grid);
static unsigned get_size_y_coords_reg2d(struct grid * grid);
static double const * get_x_coords_reg2d(struct grid * grid);
static double const * get_y_coords_reg2d(struct grid * grid);
static void set_x_coords_reg2d(struct grid * grid, double * x_coords);
static void set_y_coords_reg2d(struct grid * grid, double * y_coords);
static unsigned get_size_cell_grid_x_coords_reg2d(struct grid * grid);
static unsigned get_size_cell_grid_y_coords_reg2d(struct grid * grid);
static unsigned get_num_grid_cells_reg2d(struct grid * grid);
static unsigned get_num_grid_corners_reg2d(struct grid * grid);
static unsigned get_num_cell_corners_reg2d(struct grid * grid, unsigned cell_index);
static unsigned get_num_corner_cells_reg2d(struct grid * grid, unsigned corner_index);
static unsigned get_num_grid_edges_reg2d(struct grid * grid);
static unsigned get_num_corner_edges_reg2d(struct grid * grid, unsigned corner_index);
static unsigned get_num_cell_edges_reg2d(struct grid * grid, unsigned cell_index);
static unsigned const * get_corner_edges_reg2d(struct grid * grid, unsigned corner_index);
static unsigned const * get_cell_edge_indices_reg2d(struct grid * grid,
                                                    unsigned cell_index);
static enum yac_edge_type get_edge_type_reg2d(struct grid * grid, unsigned edge_index);
static unsigned const * get_cell_corner_indices_reg2d(struct grid * grid,
                                                      unsigned cell_index);
static unsigned const * get_corner_cell_indices_reg2d(struct grid * grid,
                                                      unsigned corner_index);
static unsigned const * get_cell_x_coord_indices_reg2d(struct grid * grid,
                                                       unsigned cell_index);
static unsigned const * get_cell_y_coord_indices_reg2d(struct grid * grid,
                                                       unsigned cell_index);
static unsigned get_corner_x_coord_index_reg2d(struct grid * grid, unsigned corner_index);
static unsigned get_corner_y_coord_index_reg2d(struct grid * grid, unsigned corner_index);
static double get_corner_x_coord_reg2d(struct grid * grid, unsigned corner_index);
static double get_corner_y_coord_reg2d(struct grid * grid, unsigned corner_index);
static int get_aux_grid_cell_reg2d(struct grid * grid, unsigned corner_index,
                                   unsigned * cell_indices, enum yac_edge_type * edge_type);
static struct dep_list get_cell_neigh_dep_list_reg2d(struct grid * grid);
static void get_boundary_corners_reg2d(struct grid * grid, unsigned * bnd_corners,
                                       unsigned * num_bnd_corners);
static struct grid * generate_cell_grid_reg2d(struct grid * grid, double * coordinates_x, 
                                              double * coordinates_y);
static struct grid * generate_edge_grid_reg2d(struct grid * grid, double * coordinates_x, 
                                              double * coordinates_y);
static struct grid * generate_subgrid_reg2d(struct grid * grid,
                                            unsigned * selected_local_cell_ids,
                                            unsigned num_local_cells,
                                            unsigned ** local_cell_ids,
                                            unsigned ** local_corner_ids,
                                            unsigned ** local_edge_ids);
static void pack_grid_reg2d(struct grid * grid, double ** dble_buf,
                            unsigned dble_buf_offset, unsigned * dble_buf_data_size,
                            unsigned * dble_buf_size, unsigned ** uint_buf,
                            unsigned uint_buf_offset, unsigned * uint_buf_data_size,
                            unsigned * uint_buf_size);
static struct grid_search * get_grid_search_reg2d(struct grid * grid);
static int cell_contains_non_gc_edge_reg2d(struct grid * grid, unsigned local_cell_id);
static void delete_grid_reg2d(struct grid * grid);

static struct grid_vtable reg2d_grid_vtable = {

   .copy                        = copy_grid_reg2d,
   .get_2d_extent               = get_2d_grid_extent_reg2d,
   .get_grid_cell               = get_grid_cell_reg2d,
   .get_grid_cell2              = get_grid_cell2_reg2d,
   .get_size_x_coords           = get_size_x_coords_reg2d,
   .get_size_y_coords           = get_size_y_coords_reg2d,
   .get_x_coords                = get_x_coords_reg2d,
   .get_y_coords                = get_y_coords_reg2d,
   .set_x_coords                = set_x_coords_reg2d,
   .set_y_coords                = set_y_coords_reg2d,
   .get_size_cell_grid_x_coords = get_size_cell_grid_x_coords_reg2d,
   .get_size_cell_grid_y_coords = get_size_cell_grid_y_coords_reg2d,
   .get_num_grid_cells          = get_num_grid_cells_reg2d,
   .get_num_grid_corners        = get_num_grid_corners_reg2d,
   .get_num_cell_corners        = get_num_cell_corners_reg2d,
   .get_num_corner_cells        = get_num_corner_cells_reg2d,
   .get_num_grid_edges          = get_num_grid_edges_reg2d,
   .get_num_corner_edges        = get_num_corner_edges_reg2d,
   .get_num_cell_edges          = get_num_cell_edges_reg2d,
   .get_corner_edges            = get_corner_edges_reg2d,
   .get_cell_edge_indices       = get_cell_edge_indices_reg2d,
   .get_edge_type               = get_edge_type_reg2d,
   .get_cell_corner_indices     = get_cell_corner_indices_reg2d,
   .get_corner_cell_indices     = get_corner_cell_indices_reg2d,
   .get_cell_x_coord_indices    = get_cell_x_coord_indices_reg2d,
   .get_cell_y_coord_indices    = get_cell_y_coord_indices_reg2d,
   .get_corner_x_coord          = get_corner_x_coord_reg2d,
   .get_corner_y_coord          = get_corner_y_coord_reg2d,
   .get_corner_x_coord_index    = get_corner_x_coord_index_reg2d,
   .get_corner_y_coord_index    = get_corner_y_coord_index_reg2d,
   .get_aux_grid_cell           = get_aux_grid_cell_reg2d,
   .get_cell_neigh_dep_list     = get_cell_neigh_dep_list_reg2d,
   .get_boundary_corners        = get_boundary_corners_reg2d,
   .generate_cell_grid          = generate_cell_grid_reg2d,
   .generate_edge_grid          = generate_edge_grid_reg2d,
   .generate_subgrid            = generate_subgrid_reg2d,
   .pack_grid                   = pack_grid_reg2d,
   .get_grid_search             = get_grid_search_reg2d,
   .cell_contains_non_gc_edge   = cell_contains_non_gc_edge_reg2d,
   .grid_delete                 = delete_grid_reg2d
};

struct reg2d_grid {

   struct grid_vtable * vtable;

   struct dep_list cell_to_neigh; //!< dependency list containing cell neighbourhood
                                  //!< information (automatically generated once it is
                                  //!< needed)

   double * cell_corners_x,       //!< latitude data
          * cell_corners_y,       //!< longitude data
          * cos_cell_corners_x,
          * sin_cell_corners_x,
          * cos_cell_corners_y,
          * sin_cell_corners_y;

   unsigned cell_corners_x_by_user; //!< indicates whether cell_corners_x was provided by
                                    //!< the user or automatically generated by a grid
                                    //!< routine
   unsigned cell_corners_y_by_user; //!< indicates whether cell_corners_y was provided by
                                    //!< the user or automatically generated by a grid
                                    //!< routine

   unsigned num_cells[2];
   unsigned cyclic[2];

   struct grid_search * grid_search;
};

static unsigned get_num_cell_edges_reg2d(struct grid * grid, unsigned cell_index) {

   return 4;
}

static unsigned get_num_cell_corners_reg2d(struct grid * grid, unsigned cell_index) {

   return 4;
}

static unsigned const * get_cell_x_coord_indices_reg2d(struct grid * grid,
                                                       unsigned cell_index) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   static unsigned cell_corners[4];
   unsigned x_index, y_index;

   y_index = cell_index / reg2d_grid->num_cells[0];
   x_index = cell_index - y_index * reg2d_grid->num_cells[0];

   cell_corners[0] = x_index;
   
   if (reg2d_grid->cyclic[0] &&
       x_index+1 == reg2d_grid->num_cells[0]) {
       
      cell_corners[1] = 0;
      cell_corners[2] = 0;
   } else {
      cell_corners[1] = x_index + 1;
      cell_corners[2] = x_index + 1;
   }
   cell_corners[3] = x_index;

   return cell_corners;
}

static unsigned const * get_cell_y_coord_indices_reg2d(struct grid * grid,
                                                       unsigned cell_index) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   static unsigned cell_corners[4];
   unsigned y_index;

   y_index = cell_index / reg2d_grid->num_cells[0];

   cell_corners[0] = y_index;
   cell_corners[1] = y_index;
   cell_corners[2] = y_index + 1;
   cell_corners[3] = y_index + 1;

   return cell_corners;
}

static double get_corner_x_coord_reg2d(struct grid * grid, unsigned corner_index) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   unsigned y_index;

   if (!reg2d_grid->cyclic[0]) {
      y_index = corner_index / (reg2d_grid->num_cells[0] + 1);
      return reg2d_grid->cell_corners_x[corner_index - y_index * (reg2d_grid->num_cells[0] + 1)];
   } else {
      y_index = corner_index / reg2d_grid->num_cells[0];
      return reg2d_grid->cell_corners_x[corner_index - y_index * reg2d_grid->num_cells[0]];
   }
}

static double get_corner_y_coord_reg2d(struct grid * grid, unsigned corner_index) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   if (!reg2d_grid->cyclic[0]) {
      return reg2d_grid->cell_corners_y[corner_index / (reg2d_grid->num_cells[0] + 1)];
   } else {
      return reg2d_grid->cell_corners_y[corner_index / reg2d_grid->num_cells[0]];
   }
}

static unsigned get_corner_x_coord_index_reg2d(struct grid * grid, unsigned corner_index) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   unsigned y_index;

   if (!reg2d_grid->cyclic[0]) {
      y_index = corner_index / (reg2d_grid->num_cells[0] + 1);
      return corner_index - y_index * (reg2d_grid->num_cells[0] + 1);
   } else {
      y_index = corner_index / reg2d_grid->num_cells[0];
      return corner_index - y_index * reg2d_grid->num_cells[0];
   }
}

static unsigned get_corner_y_coord_index_reg2d(struct grid * grid, unsigned corner_index) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   if (!reg2d_grid->cyclic[0]) {
      return corner_index / (reg2d_grid->num_cells[0] + 1);
   } else {
      return corner_index / reg2d_grid->num_cells[0];
   }
}

static int get_aux_grid_cell_reg2d(struct grid * grid, unsigned corner_index,
                                   unsigned * cell_indices, enum yac_edge_type * edge_type) {

   struct reg2d_grid * reg2d_grid = (struct reg2d_grid *)grid;

  edge_type[0] = LAT_CIRCLE;
  edge_type[1] = LON_CIRCLE;
  edge_type[2] = LAT_CIRCLE;
  edge_type[3] = LON_CIRCLE;

  return
    yac_get_aux_grid_cell_reg2d_common(
      reg2d_grid->num_cells, reg2d_grid->cyclic[0], corner_index, cell_indices);
}

static unsigned get_size_x_coords_reg2d(struct grid * grid) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   return reg2d_grid->num_cells[0] + ((reg2d_grid->cyclic[0])?0:1); 
}

static unsigned get_size_y_coords_reg2d(struct grid * grid) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   return reg2d_grid->num_cells[1] + 1;
}

static double const * get_x_coords_reg2d(struct grid * grid) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   return reg2d_grid->cell_corners_x;
}

static double const * get_y_coords_reg2d(struct grid * grid) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   return reg2d_grid->cell_corners_y;
}

static void set_x_coords_reg2d(struct grid * grid, double * x_coords) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   if (!reg2d_grid->cell_corners_x_by_user)
      free(reg2d_grid->cell_corners_x);

   reg2d_grid->cell_corners_x = x_coords;
   reg2d_grid->cell_corners_x_by_user = 1 == 1;

   if (x_coords != NULL) {

      unsigned num_x_coords = get_size_x_coords_reg2d(grid);

      if (reg2d_grid->cos_cell_corners_x == NULL)
         reg2d_grid->cos_cell_corners_x =
            malloc(num_x_coords * sizeof(*(reg2d_grid->cos_cell_corners_x)));
      if (reg2d_grid->sin_cell_corners_x == NULL)
         reg2d_grid->sin_cell_corners_x =
            malloc(num_x_coords * sizeof(*(reg2d_grid->sin_cell_corners_x)));

      for (unsigned i = 0; i < num_x_coords; ++i) {

         reg2d_grid->cos_cell_corners_x[i] = cos(reg2d_grid->cell_corners_x[i]);
         reg2d_grid->sin_cell_corners_x[i] = sin(reg2d_grid->cell_corners_x[i]);
      }
   } else {

      free(reg2d_grid->cos_cell_corners_x);
      reg2d_grid->cos_cell_corners_x = NULL;
      free(reg2d_grid->sin_cell_corners_x);
      reg2d_grid->sin_cell_corners_x = NULL;
   }
}

static void set_y_coords_reg2d(struct grid * grid, double * y_coords) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   if (!reg2d_grid->cell_corners_y_by_user)
      free(reg2d_grid->cell_corners_y);

   reg2d_grid->cell_corners_y = y_coords;
   reg2d_grid->cell_corners_y_by_user = 1 == 1;

   if (y_coords != NULL) {

      unsigned num_y_coords = get_size_y_coords_reg2d(grid);

      if (reg2d_grid->cos_cell_corners_y == NULL)
         reg2d_grid->cos_cell_corners_y =
            malloc(num_y_coords * sizeof(*(reg2d_grid->cos_cell_corners_y)));
      if (reg2d_grid->sin_cell_corners_y == NULL)
         reg2d_grid->sin_cell_corners_y =
            malloc(num_y_coords * sizeof(*(reg2d_grid->sin_cell_corners_y)));

      for (unsigned i = 0; i < num_y_coords; ++i) {

         reg2d_grid->cos_cell_corners_y[i] = cos(reg2d_grid->cell_corners_y[i]);
         reg2d_grid->sin_cell_corners_y[i] = sin(reg2d_grid->cell_corners_y[i]);
      }
   } else {

      free(reg2d_grid->cos_cell_corners_y);
      reg2d_grid->cos_cell_corners_y = NULL;
      free(reg2d_grid->sin_cell_corners_y);
      reg2d_grid->sin_cell_corners_y = NULL;
   }
}

static unsigned get_size_cell_grid_x_coords_reg2d(struct grid * grid) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   return reg2d_grid->num_cells[0]; 
}

static unsigned get_size_cell_grid_y_coords_reg2d(struct grid * grid) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   return reg2d_grid->num_cells[1];
}

static unsigned get_num_grid_cells_reg2d(struct grid * grid) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   return reg2d_grid->num_cells[0] * reg2d_grid->num_cells[1];
}

static unsigned get_num_grid_corners_reg2d(struct grid * grid) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   return (reg2d_grid->num_cells[0]+ ((reg2d_grid->cyclic[0])?0:1)) *
          (reg2d_grid->num_cells[1]+1);
}

static unsigned get_num_grid_edges_reg2d(struct grid * grid) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   return (reg2d_grid->num_cells[0] + ((reg2d_grid->cyclic[0])?0:1)) *
          reg2d_grid->num_cells[1] + reg2d_grid->num_cells[0] *
          (reg2d_grid->num_cells[1] + 1);
}

static unsigned get_num_corner_edges_reg2d(struct grid * grid, unsigned corner_index) {

  struct reg2d_grid * reg2d_grid = (struct reg2d_grid *)grid;

  return
    yac_get_num_corner_edges_reg2d_common(
      reg2d_grid->num_cells, reg2d_grid->cyclic[0], corner_index);
}

static unsigned get_num_corner_cells_reg2d(struct grid * grid, unsigned corner_index) {

   // struct reg2d_grid * reg2d_grid;

   // reg2d_grid = (struct reg2d_grid *)grid;

   return 1 << (get_num_corner_edges_reg2d(grid, corner_index) - 2);
}

static void generate_cell_neigh_dep_reg2d(struct grid * grid) {

  struct reg2d_grid * reg2d_grid = (struct reg2d_grid *)grid;

  yac_generate_cell_neigh_dep_reg2d_common(
    &(reg2d_grid->cell_to_neigh), reg2d_grid->num_cells, reg2d_grid->cyclic);
}

static struct dep_list get_cell_neigh_dep_list_reg2d(struct grid * grid) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   if (reg2d_grid->cell_to_neigh.num_elements == 0)
      generate_cell_neigh_dep_reg2d(grid);
   return reg2d_grid->cell_to_neigh;
}

static void get_boundary_corners_reg2d(struct grid * grid, unsigned * bnd_corners,
                                       unsigned * num_bnd_corners) {

  struct reg2d_grid * reg2d_grid = (struct reg2d_grid *)grid;

  yac_get_boundary_corners_reg2d_common(
    reg2d_grid->num_cells, reg2d_grid->cyclic[0], bnd_corners, num_bnd_corners);
}

static void get_grid_cell_reg2d(struct grid * grid, unsigned cell_index,
                                struct grid_cell * cell) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   unsigned num_corners;

   num_corners = get_num_cell_corners_reg2d(grid, cell_index);

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
   unsigned const * x_indices, * y_indices;
   unsigned i;

   x_indices = get_cell_x_coord_indices_reg2d(grid, cell_index);
   y_indices = get_cell_y_coord_indices_reg2d(grid, cell_index);

   for (i = 0; i < 4; ++i) {
      cell->coordinates_x[i] = reg2d_grid->cell_corners_x[x_indices[i]];
      cell->coordinates_y[i] = reg2d_grid->cell_corners_y[y_indices[i]];
      cell->coordinates_xyz[0+i*3] = reg2d_grid->cos_cell_corners_x[x_indices[i]] *
                                     reg2d_grid->cos_cell_corners_y[y_indices[i]];
      cell->coordinates_xyz[1+i*3] = reg2d_grid->sin_cell_corners_x[x_indices[i]] *
                                     reg2d_grid->cos_cell_corners_y[y_indices[i]];
      cell->coordinates_xyz[2+i*3] = reg2d_grid->sin_cell_corners_y[y_indices[i]];
   }


   // set the edge type for the cell
   cell->edge_type[0] = LAT_CIRCLE;
   cell->edge_type[1] = LON_CIRCLE;
   cell->edge_type[2] = LAT_CIRCLE;
   cell->edge_type[3] = LON_CIRCLE;
}

static void get_grid_cell2_reg2d(struct grid * grid, unsigned cell_index,
                                 struct grid_cell * cell,
                                 struct bounding_circle * bnd_circle) {

   get_grid_cell_reg2d(grid, cell_index, cell);

   yac_get_cell_bounding_circle_reg_quad(cell->coordinates_xyz+0, cell->coordinates_xyz+3, cell->coordinates_xyz+6,
                                         cell->coordinates_xyz+9, bnd_circle);
}

static unsigned const * get_corner_edges_reg2d(struct grid * grid, unsigned corner_index) {

  struct reg2d_grid * reg2d_grid = (struct reg2d_grid *)grid;

  return
    yac_get_corner_edges_reg2d_common(
      reg2d_grid->num_cells, reg2d_grid->cyclic[0], corner_index);
}

static unsigned const * get_cell_edge_indices_reg2d(struct grid * grid,
                                                    unsigned cell_index) {

   struct reg2d_grid * reg2d_grid = (struct reg2d_grid *)grid;

   if (!reg2d_grid->cyclic[0]) {

      unsigned row = cell_index / reg2d_grid->num_cells[0];
      unsigned column = cell_index - reg2d_grid->num_cells[0] * row;

      return
        yac_get_cell_edge_indices_reg2d_common_non_cyclic(
          reg2d_grid->num_cells, row, column);

   } else {

      unsigned row = cell_index / reg2d_grid->num_cells[0];
      unsigned column = cell_index - reg2d_grid->num_cells[0] * row;

      return
        yac_get_cell_edge_indices_reg2d_common_cyclic(
          reg2d_grid->num_cells, row, column);

   } // if (!reg2d_grid->cyclic[0])
}

static enum yac_edge_type get_edge_type_reg2d(struct grid * grid, unsigned edge_index) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   if (!reg2d_grid->cyclic[0]) {

      unsigned row_index;

      row_index = edge_index / (2*reg2d_grid->num_cells[0]+1);

      if (row_index == reg2d_grid->num_cells[1])
         return LAT_CIRCLE;

      unsigned temp;

      temp = edge_index - row_index * (2*reg2d_grid->num_cells[0]+1);

      if (temp == 2*reg2d_grid->num_cells[0])
         return LON_CIRCLE;

      return (temp&1)?LON_CIRCLE:LAT_CIRCLE;

   } else {

      unsigned row_index;

      row_index = edge_index / (2*reg2d_grid->num_cells[0]);

      if (row_index == reg2d_grid->num_cells[1])
         return LAT_CIRCLE;

      unsigned temp;

      temp = edge_index - row_index * 2*reg2d_grid->num_cells[0];

      if (temp == 0)
         return LAT_CIRCLE;

      if (temp == 2*reg2d_grid->num_cells[0] - 1)
         return LON_CIRCLE;

      return (edge_index&1)?LAT_CIRCLE:LON_CIRCLE;
   }
}

static unsigned const * get_cell_corner_indices_reg2d(
  struct grid * grid, unsigned cell_index) {

  struct reg2d_grid * reg2d_grid = (struct reg2d_grid *)grid;

  return
    yac_get_cell_corner_indices_reg2d_common(
      reg2d_grid->num_cells, reg2d_grid->cyclic[0], cell_index);
}

static unsigned const * get_corner_cell_indices_reg2d(struct grid * grid,
                                                      unsigned corner_index) {

  struct reg2d_grid * reg2d_grid = (struct reg2d_grid *)grid;

  return
    yac_get_corner_cell_indices_reg2d_common(
      reg2d_grid->num_cells, reg2d_grid->cyclic[0], corner_index);
}

static void get_2d_grid_extent_reg2d(struct grid * grid, double (* extent)[2]) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   unsigned i;

   double const tol = 1e-12;

   //--------------------------------
   // compute the extent of the grid
   //--------------------------------

   if (get_num_grid_corners_reg2d(grid) == 0) return;
   if (extent == NULL) return;

   extent[1][0] = extent[1][1] = reg2d_grid->cell_corners_y[0];

   for (i = 1; i < yac_get_size_y_coords(grid); ++i) {

      if (extent[1][0] > reg2d_grid->cell_corners_y[i])
          extent[1][0] = reg2d_grid->cell_corners_y[i];
      if (extent[1][1] < reg2d_grid->cell_corners_y[i])
          extent[1][1] = reg2d_grid->cell_corners_y[i];
   }

   if (reg2d_grid->cyclic[0]) {

      extent[0][0] = 0;
      extent[0][1] = 2*M_PI;
   } else {

      unsigned * checked_corners;

      checked_corners = calloc(get_num_grid_corners_reg2d(grid), sizeof (checked_corners[0]));

      unsigned * corner_stack;
      unsigned corner_stack_size;
      
      corner_stack_size = 0;
      corner_stack = malloc(yac_get_num_grid_corners(grid) * sizeof(corner_stack[0]));

      corner_stack[corner_stack_size++] = 0;
      extent[0][0] = extent[0][1] = reg2d_grid->cell_corners_x[yac_get_corner_x_coord_index(grid, 0)];

      unsigned curr_corner;
      unsigned const * curr_neigh_corners;
      unsigned curr_num_neigh_corners;

      while (corner_stack_size > 0) {

         curr_corner = corner_stack[--corner_stack_size];
         curr_num_neigh_corners = get_num_corner_edges_reg2d(grid, curr_corner);
         curr_neigh_corners = get_corner_edges_reg2d(grid, curr_corner);

         double edge_angle;
         double edge_x_coord[2];

         edge_x_coord[0] =reg2d_grid->cell_corners_x[yac_get_corner_x_coord_index(grid, curr_corner)];

         for (i = 0; i < curr_num_neigh_corners; ++i) {

            // process every edge only once
            if (curr_corner > curr_neigh_corners[i]) continue;
            else if (!checked_corners[curr_neigh_corners[i]]) {
               corner_stack[corner_stack_size++] = curr_neigh_corners[i];
               checked_corners[curr_neigh_corners[i]] = 1 == 1;
            }

            edge_x_coord[1] = reg2d_grid->cell_corners_x[get_corner_x_coord_index_reg2d(grid, curr_neigh_corners[i])];
            edge_angle = fabs(get_angle(edge_x_coord[0], edge_x_coord[1]));

            // if the current edge crosses the up boundary of the current extent
            if (fabs(get_angle(edge_x_coord[0], extent[0][1])) <= edge_angle + tol &&
                fabs(get_angle(edge_x_coord[1], extent[0][1])) <= edge_angle + tol &&
                get_angle(extent[0][1], edge_x_coord[1]) <= 0)
               extent[0][1] = edge_x_coord[1];

            // if the current edge crosses the up boundary of the current extent
            if (fabs(get_angle(edge_x_coord[0], extent[0][0])) <= edge_angle + tol &&
                fabs(get_angle(edge_x_coord[1], extent[0][0])) <= edge_angle + tol &&
                get_angle(extent[0][0], edge_x_coord[1]) >= 0)
               extent[0][0] = edge_x_coord[1];

            // if we already have 360 degree
            if (get_angle(extent[0][1], extent[0][0]) <= tol) {
               extent[0][0] = 0;
               extent[0][1] = 2*M_PI;
               corner_stack_size = 0;
               break;
            }
         }
      }

      free(corner_stack);
      free(checked_corners);
   }

   //----------------------
   // check for pole cells
   //----------------------

   struct grid_cell cell;

   yac_init_grid_cell(&cell);

   // test whether the grid data includes the pole
   for (i = 0; i < get_num_grid_cells_reg2d(grid); ++i) {
        
      // get the current cell
      get_grid_cell_reg2d(grid, i, &cell);

      // check the cell
      if (yac_cell_covers_pole(cell.num_corners, cell.coordinates_x, cell.coordinates_y)) {

         if (cell.coordinates_y[0] > 0)
            extent[1][1] = M_PI_2;
         else
            extent[1][0] = -M_PI_2;
      }
   } // i

   yac_free_grid_cell(&cell);
}

static struct grid * generate_cell_grid_reg2d(struct grid * grid, double * coordinates_x, 
                                              double * coordinates_y) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   unsigned num_cells[2] = {reg2d_grid->num_cells[0]-1,
                            reg2d_grid->num_cells[1]-1};

   if (reg2d_grid->cyclic[0])
      num_cells[0]++;

   return yac_reg2d_grid_new(coordinates_x, coordinates_y, num_cells, reg2d_grid->cyclic);
}

static struct grid * generate_edge_grid_reg2d(
  struct grid * grid, double * coordinates_x, double * coordinates_y) {

  struct reg2d_grid * reg2d_grid = (struct reg2d_grid *)grid;

  return
    yac_generate_edge_grid_reg2d_common(
      reg2d_grid->num_cells, reg2d_grid->cyclic, coordinates_x, coordinates_y);
}

static struct grid * copy_grid_reg2d(struct grid * grid) {

   struct reg2d_grid * reg2d_grid = (struct reg2d_grid *)grid;

   struct reg2d_grid * copy = malloc(1 * sizeof(*copy));

   copy->vtable = &reg2d_grid_vtable;
   yac_init_dep_list(&copy->cell_to_neigh);

   unsigned x_coords_size = get_size_x_coords_reg2d(grid) * sizeof(*copy->cell_corners_x);
   unsigned y_coords_size = get_size_y_coords_reg2d(grid) * sizeof(*copy->cell_corners_y);

   copy->cell_corners_x = malloc(x_coords_size);
   copy->cell_corners_y = malloc(y_coords_size);
   copy->sin_cell_corners_x = malloc(x_coords_size);
   copy->cos_cell_corners_x = malloc(x_coords_size);
   copy->sin_cell_corners_y = malloc(y_coords_size);
   copy->cos_cell_corners_y = malloc(y_coords_size);
   memcpy(copy->cell_corners_x, reg2d_grid->cell_corners_x, x_coords_size);
   memcpy(copy->cell_corners_y, reg2d_grid->cell_corners_y, y_coords_size);
   memcpy(copy->sin_cell_corners_x, reg2d_grid->sin_cell_corners_x, x_coords_size);
   memcpy(copy->cos_cell_corners_x, reg2d_grid->cos_cell_corners_x, x_coords_size);
   memcpy(copy->sin_cell_corners_y, reg2d_grid->sin_cell_corners_y, y_coords_size);
   memcpy(copy->cos_cell_corners_y, reg2d_grid->cos_cell_corners_y, y_coords_size);

   copy->cell_corners_x_by_user = 1 == 0;
   copy->cell_corners_y_by_user = 1 == 0;

   copy->num_cells[0] = reg2d_grid->num_cells[0];
   copy->num_cells[1] = reg2d_grid->num_cells[1];
   copy->cyclic[0] = reg2d_grid->cyclic[0];
   copy->cyclic[1] = reg2d_grid->cyclic[1];

   copy->grid_search = NULL;

   return (struct grid *)copy;
}

static void pack_grid_reg2d(struct grid * grid, double ** dble_buf,
                            unsigned dble_buf_offset, unsigned * dble_buf_data_size,
                            unsigned * dble_buf_size, unsigned ** uint_buf,
                            unsigned uint_buf_offset, unsigned * uint_buf_data_size,
                            unsigned * uint_buf_size) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   unsigned size_x_coords, size_y_coords;

   size_x_coords = get_size_x_coords_reg2d(grid);
   size_y_coords = get_size_y_coords_reg2d(grid);

   unsigned required_dble_buf_size, required_uint_buf_size;

   required_dble_buf_size = size_x_coords + size_y_coords;
   required_uint_buf_size = 5;

   ENSURE_ARRAY_SIZE(*dble_buf, *dble_buf_size, dble_buf_offset+required_dble_buf_size);
   ENSURE_ARRAY_SIZE(*uint_buf, *uint_buf_size, uint_buf_offset+required_uint_buf_size);

   memcpy((*dble_buf)+dble_buf_offset, reg2d_grid->cell_corners_x, size_x_coords * sizeof(**dble_buf));
   dble_buf_offset += size_x_coords;
   memcpy((*dble_buf)+dble_buf_offset, reg2d_grid->cell_corners_y, size_y_coords * sizeof(**dble_buf));

   (*uint_buf)[uint_buf_offset+0] = yac_hash("REG2D");
   (*uint_buf)[uint_buf_offset+1] = reg2d_grid->num_cells[0];
   (*uint_buf)[uint_buf_offset+2] = reg2d_grid->num_cells[1];
   (*uint_buf)[uint_buf_offset+3] = reg2d_grid->cyclic[0];
   (*uint_buf)[uint_buf_offset+4] = reg2d_grid->cyclic[1];

   *dble_buf_data_size = required_dble_buf_size;
   *uint_buf_data_size = required_uint_buf_size;
}

struct grid * yac_unpack_reg2d_grid(double * dble_buf, unsigned * dble_buf_data_size,
                                    unsigned * uint_buf, unsigned * uint_buf_data_size) {

   unsigned num_cells[2], cyclic[2];

   num_cells[0] = uint_buf[1];
   num_cells[1] = uint_buf[2];
   cyclic[0] = uint_buf[3];
   cyclic[1] = uint_buf[4];

   unsigned size_x_coords, size_y_coords;

   size_x_coords = num_cells[0] + ((cyclic[0])?0:1);
   size_y_coords = num_cells[1] + 1;

   double * cell_corners_x, * cell_corners_y;

   cell_corners_x = malloc(size_x_coords * sizeof(*cell_corners_x));
   cell_corners_y = malloc(size_y_coords * sizeof(*cell_corners_y));

   memcpy(cell_corners_x, dble_buf, size_x_coords * sizeof(*cell_corners_x));
   memcpy(cell_corners_y, dble_buf+size_x_coords, size_y_coords * sizeof(*cell_corners_y));

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)yac_reg2d_grid_new(cell_corners_x, cell_corners_y, num_cells, cyclic);

   *dble_buf_data_size = size_x_coords + size_y_coords;
   *uint_buf_data_size = 5;

   reg2d_grid->cell_corners_x_by_user = 1 == 0;
   reg2d_grid->cell_corners_y_by_user = 1 == 0;

   return (struct grid *)reg2d_grid;
}

static void delete_grid_reg2d(struct grid * grid) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   yac_free_dep_list(&(reg2d_grid->cell_to_neigh));
   if (!reg2d_grid->cell_corners_x_by_user)
      free(reg2d_grid->cell_corners_x);
   if (!reg2d_grid->cell_corners_y_by_user)
      free(reg2d_grid->cell_corners_y);
   free(reg2d_grid->sin_cell_corners_x);
   free(reg2d_grid->cos_cell_corners_x);
   free(reg2d_grid->sin_cell_corners_y);
   free(reg2d_grid->cos_cell_corners_y);

   if (reg2d_grid->grid_search != NULL)
      yac_delete_grid_search(reg2d_grid->grid_search);

   free(reg2d_grid);
}

static struct grid * generate_subgrid_reg2d(struct grid * grid,
                                            unsigned * selected_local_cell_ids,
                                            unsigned num_local_cells,
                                            unsigned ** local_cell_ids,
                                            unsigned ** local_corner_ids,
                                            unsigned ** local_edge_ids) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   if (num_local_cells == 0)
      yac_internal_abort_message ( "ERROR: cells need to be defined.", __FILE__, __LINE__ );

   // compute the index range of the subgrid

   unsigned index_range[2][2];

   index_range[0][0] = index_range[0][1] = get_cell_x_coord_indices_reg2d(grid, selected_local_cell_ids[0])[0];
   index_range[1][0] = index_range[1][1] = get_cell_y_coord_indices_reg2d(grid, selected_local_cell_ids[0])[0];

   unsigned i, j, k, l; 

   for (i = 0; i < num_local_cells; ++i) {

      unsigned const * cell_x_coords, * cell_y_coords;

      cell_x_coords = get_cell_x_coord_indices_reg2d(grid, selected_local_cell_ids[i]);
      cell_y_coords = get_cell_y_coord_indices_reg2d(grid, selected_local_cell_ids[i]);

      for (j = 0; j < 4; ++j) {

         if (index_range[0][0] > cell_x_coords[j]) index_range[0][0] = cell_x_coords[j];
         if (index_range[0][1] < cell_x_coords[j]) index_range[0][1] = cell_x_coords[j];
         if (index_range[1][0] > cell_y_coords[j]) index_range[1][0] = cell_y_coords[j];
         if (index_range[1][1] < cell_y_coords[j]) index_range[1][1] = cell_y_coords[j];
      }
   }

   unsigned index_range_extent[2];

   index_range_extent[0] = index_range[0][1] - index_range[0][0] + 1;
   index_range_extent[1] = index_range[1][1] - index_range[1][0] + 1;

   // generate the subgrid

   unsigned num_cells[2], cyclic[2];
   double * cell_x_coords, * cell_y_coords;

   cell_x_coords = malloc(index_range_extent[0] * sizeof(*cell_x_coords));
   cell_y_coords = malloc(index_range_extent[1] * sizeof(*cell_y_coords));

   memcpy(cell_x_coords, reg2d_grid->cell_corners_x + index_range[0][0], index_range_extent[0] * sizeof(*cell_x_coords));
   memcpy(cell_y_coords, reg2d_grid->cell_corners_y + index_range[1][0], index_range_extent[1] * sizeof(*cell_y_coords));

   cyclic[0] = (index_range_extent[0] == get_size_x_coords_reg2d(grid)) &&
               reg2d_grid->cyclic[0];
   cyclic[1] = 0;

   num_cells[0] = index_range_extent[0] - ((cyclic[0])?0:1);
   num_cells[1] = index_range_extent[1] - 1;

   struct grid * subgrid;

   subgrid = yac_reg2d_grid_new(cell_x_coords, cell_y_coords, num_cells, cyclic);

   ((struct reg2d_grid *)subgrid)->cell_corners_x_by_user = 1 == 0;
   ((struct reg2d_grid *)subgrid)->cell_corners_y_by_user = 1 == 0;

   //generate the local id mapping between grid and subgrid
   *local_cell_ids = malloc(get_num_grid_cells_reg2d(subgrid) * sizeof(**local_cell_ids));
   *local_corner_ids = malloc(get_num_grid_corners_reg2d(subgrid) * sizeof(**local_corner_ids));
   *local_edge_ids = malloc(get_num_grid_edges_reg2d(subgrid) * sizeof(**local_edge_ids));

   unsigned index_range_stride;

   index_range_stride = get_size_x_coords_reg2d(grid);

   unsigned num_grid_edges;
   unsigned * edge_mask;
   unsigned const * cell_edges;

   num_grid_edges = get_num_grid_edges_reg2d(grid);
   edge_mask = calloc(num_grid_edges, sizeof(*edge_mask));

   k = 0;

   for (i = index_range[1][0]; i < index_range[1][1]; ++i) {

      for (j = index_range[0][0]; j < index_range[0][1] + ((cyclic[0])?1:0); ++j) {

         (*local_cell_ids)[k] = j + i * (index_range_stride - ((reg2d_grid->cyclic[0])?0:1));

         cell_edges = get_cell_edge_indices_reg2d(grid, (*local_cell_ids)[k++]);

         for (l = 0; l < 4; ++l)
            edge_mask[cell_edges[l]] = 1;
      }
   }

   k = 0;

   for (i = 0; i < num_grid_edges; ++i)
      if (edge_mask[i]) (*local_edge_ids)[k++] = i;

   free(edge_mask);

   if (k != get_num_grid_edges_reg2d(subgrid))
      yac_internal_abort_message ( "ERROR: Number of edges does not match.", __FILE__, __LINE__ );

   k = 0;

   for (i = index_range[1][0]; i <= index_range[1][1]; ++i) {

      for (j = index_range[0][0]; j <= index_range[0][1]; ++j) {

         (*local_corner_ids)[k++] = j + i * index_range_stride;
      }
   }

   return subgrid;
}

struct grid * yac_reg2d_grid_new(double * coordinates_x, double * coordinates_y,
                                 unsigned const num_cells[2], unsigned const cyclic[2]) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = malloc(1 * sizeof(*reg2d_grid));

   reg2d_grid->vtable = &reg2d_grid_vtable;
   yac_init_dep_list(&reg2d_grid->cell_to_neigh);

   reg2d_grid->cell_corners_x = coordinates_x;
   reg2d_grid->cell_corners_y = coordinates_y;
   reg2d_grid->num_cells[0] = num_cells[0];
   reg2d_grid->num_cells[1] = num_cells[1];
   reg2d_grid->cyclic[0] = cyclic[0];
   reg2d_grid->cyclic[1] = cyclic[1];
   reg2d_grid->cell_corners_x_by_user = 1 == 1;
   reg2d_grid->cell_corners_y_by_user = 1 == 1;
   reg2d_grid->grid_search = NULL;

   if (coordinates_x != NULL) {
      unsigned num_x_coords = get_size_x_coords_reg2d((struct grid *) reg2d_grid);
      reg2d_grid->sin_cell_corners_x =
         malloc(num_x_coords * sizeof(*(reg2d_grid->sin_cell_corners_x)));
      reg2d_grid->cos_cell_corners_x = 
         malloc(num_x_coords * sizeof(*(reg2d_grid->cos_cell_corners_x)));

      for (unsigned i = 0; i < num_x_coords; ++i) {

         reg2d_grid->sin_cell_corners_x[i] = sin(reg2d_grid->cell_corners_x[i]);
         reg2d_grid->cos_cell_corners_x[i] = cos(reg2d_grid->cell_corners_x[i]);
      }
   } else {
      reg2d_grid->sin_cell_corners_x = NULL;
      reg2d_grid->cos_cell_corners_x = NULL;
   }

   if (coordinates_y != NULL) {
      unsigned num_y_coords = get_size_y_coords_reg2d((struct grid *) reg2d_grid);
      reg2d_grid->sin_cell_corners_y =
         malloc(num_y_coords * sizeof(*(reg2d_grid->sin_cell_corners_y)));
      reg2d_grid->cos_cell_corners_y = 
         malloc(num_y_coords * sizeof(*(reg2d_grid->cos_cell_corners_y)));

      for (unsigned i = 0; i < num_y_coords; ++i) {

         reg2d_grid->sin_cell_corners_y[i] = sin(reg2d_grid->cell_corners_y[i]);
         reg2d_grid->cos_cell_corners_y[i] = cos(reg2d_grid->cell_corners_y[i]);
      }
   } else {
      reg2d_grid->sin_cell_corners_y = NULL;
      reg2d_grid->cos_cell_corners_y = NULL;
   }

   return (struct grid *) reg2d_grid;
}

static int cell_contains_non_gc_edge_reg2d(
  struct grid * grid, unsigned local_cell_id) {

  return 1;
}

static struct grid_search * get_grid_search_reg2d(struct grid * grid) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   if (reg2d_grid->grid_search == NULL)
      reg2d_grid->grid_search = yac_sphere_part_search_new(grid);

   return reg2d_grid->grid_search;
}

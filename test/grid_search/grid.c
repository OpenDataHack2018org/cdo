/**
 * @file grid.c
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *                                 Rene Redler <rene.redler@mpimet.mpg.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler <rene.redler@mpimet.mpg.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
 *             Thomas Jahns <jahns@dkrz.de>
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
#include <stdio.h>

#include "grid.h"
#include "dep_list.h"
#include "math.h"
#include "geometry.h"
#include "utils.h"
#include "ensure_array_size.h"
#include "grid_cell.h"

struct grid * copy_grid(struct grid * grid) {

   if (grid->vtable->copy != NULL)
      return grid->vtable->copy(grid);
   else
      abort_message ( "ERROR: missing implementation for copy_grid",
                      __FILE__, __LINE__ );

   return NULL;
}

void get_grid_cell (struct grid * grid, unsigned cell_index,
                    struct grid_cell * cell) {

   if (grid->vtable->get_grid_cell != NULL)
       grid->vtable->get_grid_cell(grid, cell_index, cell);
   else
      abort_message ( "ERROR: missing implementation for get_grid_cell",
                      __FILE__, __LINE__ );
}

void get_grid_cell2 (struct grid * grid, unsigned cell_index,
                     struct grid_cell * cell,
                     struct bounding_circle * bnd_circle) {

   if (grid->vtable->get_grid_cell2 != NULL)
       grid->vtable->get_grid_cell2(grid, cell_index, cell, bnd_circle);
   else
      abort_message ( "ERROR: missing implementation for get_grid_cell2",
                      __FILE__, __LINE__ );
}

void get_2d_grid_extent(struct grid * grid, double (* extent)[2]) {

   if (grid->vtable->get_2d_extent != NULL)
      grid->vtable->get_2d_extent(grid, extent);
   else
      abort_message ( "ERROR: missing implementation vor get_2d_grid_extent",
                      __FILE__, __LINE__ );
}

unsigned get_size_x_coords(struct grid * grid) {

   if (grid->vtable->get_size_x_coords != NULL)
      return grid->vtable->get_size_x_coords(grid);
   else
      abort_message ( "ERROR: missing implementation for get_size_x_coords",
                       __FILE__, __LINE__ );

   // routine will never reach this point
   return -1;
}

unsigned get_size_y_coords(struct grid * grid) {

   if (grid->vtable->get_size_y_coords != NULL)
      return grid->vtable->get_size_y_coords(grid);
   else
      abort_message ( "ERROR: missing implementation for get_size_y_coords",
                      __FILE__, __LINE__ );

   // routine will never reach this point
   return -1;
}

double const * get_x_coords(struct grid * grid) {

   if (grid->vtable->get_x_coords != NULL)
      return grid->vtable->get_x_coords(grid);
   else
      abort_message ( "ERROR: missing implementation for get_x_coords",
                      __FILE__, __LINE__ );

   // routine will never reach this point
   return NULL;
}

double const * get_y_coords(struct grid * grid) {

   if (grid->vtable->get_y_coords != NULL)
      return grid->vtable->get_y_coords(grid);
   else
      abort_message ( "ERROR: missing implementation for get_y_coords",
                      __FILE__, __LINE__ );

   // routine will never reach this point
   return NULL;
}

void set_x_coords(struct grid * grid, double * x_coords) {

   if (grid->vtable->set_x_coords != NULL)
      grid->vtable->set_x_coords(grid, x_coords);
   else
      abort_message ( "ERROR: missing implementation for set_x_coords",
                      __FILE__, __LINE__ );
}

void set_y_coords(struct grid * grid, double * y_coords) {

   if (grid->vtable->set_y_coords != NULL)
      grid->vtable->set_y_coords(grid, y_coords);
   else
      abort_message ( "ERROR: missing implementation for set_y_coords",
                      __FILE__, __LINE__ );
}

unsigned get_size_cell_grid_x_coords(struct grid * grid) {

   if (grid->vtable->get_size_cell_grid_x_coords != NULL)
      return grid->vtable->get_size_cell_grid_x_coords(grid);
   else
      abort_message ( "ERROR: missing implementation for get_size_cell_grid_x_coords",
                      __FILE__, __LINE__ );

   // routine will never reach this point
   return -1;
}

unsigned get_size_cell_grid_y_coords(struct grid * grid) {

   if (grid->vtable->get_size_cell_grid_y_coords != NULL)
      return grid->vtable->get_size_cell_grid_y_coords(grid);
   else
      abort_message ( "ERROR: missing implementation for get_size_cell_grid_y_coords",
                      __FILE__, __LINE__ );

   // routine will never reach this point
   return -1;
}

unsigned get_num_grid_cells (struct grid * grid) {

   if (grid->vtable->get_num_grid_cells != NULL)
      return grid->vtable->get_num_grid_cells(grid);
   else
      abort_message ( "ERROR: missing implementation for get_num_grid_cells",
                      __FILE__, __LINE__ );

   // routine will never reach this point
   return -1;
}

unsigned get_num_grid_corners (struct grid * grid) {

   if (grid->vtable->get_num_grid_corners != NULL)
      return grid->vtable->get_num_grid_corners(grid);
   else
      abort_message ( "ERROR: missing implementation for get_num_grid_corners",
                      __FILE__, __LINE__ );

   // routine will never reach this point
   return -1;
}

unsigned get_num_cell_corners (struct grid * grid, unsigned cell_index) {

   if (grid->vtable->get_num_cell_corners != NULL)
      return grid->vtable->get_num_cell_corners(grid, cell_index);
   else
      abort_message ( "ERROR: missing implementation for get_num_cell_corners",
                      __FILE__, __LINE__ );

   // routine will never reach this point
   return -1;
}

unsigned get_num_corner_cells (struct grid * grid, unsigned corner_index) {

   if (grid->vtable->get_num_corner_cells != NULL)
      return grid->vtable->get_num_corner_cells(grid, corner_index);
   else
      abort_message ( "ERROR: missing implementation for get_num_corner_cells",
                      __FILE__, __LINE__ );

   // routine will never reach this point
   return -1;
}

unsigned get_num_cell_edges (struct grid * grid, unsigned cell_index) {

   if (grid->vtable->get_num_cell_edges != NULL)
      return get_num_cell_corners(grid, cell_index);
   else
      abort_message ( "ERROR: missing implementation for get_num_cell_edges",
                      __FILE__, __LINE__);

   return -1;
}

unsigned get_num_grid_edges (struct grid * grid) {

   if (grid->vtable->get_num_grid_edges != NULL)
      return grid->vtable->get_num_grid_edges(grid);
   else
      abort_message ( "ERROR: missing for get_num_grid_edges", __FILE__, __LINE__ );

   // routine will never reach this point
   return -1;
}

unsigned get_num_corner_edges (struct grid * grid, unsigned corner_index) {

   if (grid->vtable->get_num_corner_edges != NULL)
      return grid->vtable->get_num_corner_edges(grid, corner_index);
   else
      abort_message ( "ERROR: missing implementation for get_num_corner_edges",
                      __FILE__, __LINE__ );

   // routine will never reach this point
   return -1;
}

unsigned const * get_corner_edges (struct grid * grid, unsigned corner_index) {

   if (grid->vtable->get_corner_edges != NULL)
      return grid->vtable->get_corner_edges(grid, corner_index);
   else
      abort_message ( "ERROR: missing implementation for get_corner_edges",
                      __FILE__, __LINE__ );

   // routine will never reach this point
   return NULL;
}

unsigned const * get_cell_edge_indices (struct grid * grid, unsigned cell_index) {

   if (grid->vtable->get_cell_edge_indices != NULL)
      return grid->vtable->get_cell_edge_indices(grid, cell_index);
   else
      abort_message ( "ERROR: missing implementation for get_cell_edge_indices",
                      __FILE__, __LINE__ );

   // routine will never reach this point
   return NULL;
}

enum edge_type get_edge_type(struct grid * grid, unsigned edge_index) {

   if (grid->vtable->get_edge_type != NULL)
      return grid->vtable->get_edge_type(grid, edge_index);
   else
      abort_message ( "ERROR: missing implementation for get_edge_type",
                      __FILE__, __LINE__ );

   // routine will never reach this point
   exit(EXIT_FAILURE);
}

unsigned const * get_cell_corner_indices (struct grid * grid, unsigned cell_index) {

   if (grid->vtable->get_cell_corner_indices != NULL)
      return grid->vtable->get_cell_corner_indices(grid, cell_index);
   else
      abort_message ( "ERROR: missing implementation for get_cell_corner_indices",
                      __FILE__, __LINE__ );

   // routine will never reach this point
   return NULL;
}

unsigned const * get_corner_cell_indices (struct grid * grid, unsigned corner_index) {

   if (grid->vtable->get_corner_cell_indices != NULL)
      return grid->vtable->get_corner_cell_indices(grid, corner_index);
   else
      abort_message ( "ERROR: missing implementation for get_corner_cell_indices",
                      __FILE__, __LINE__ );

   // routine will never reach this point
   return NULL;
}

unsigned const * get_cell_x_coord_indices (struct grid * grid, unsigned cell_index) {

   if (grid->vtable->get_cell_x_coord_indices != NULL)
      return grid->vtable->get_cell_x_coord_indices(grid, cell_index);
   else
      abort_message ( "ERROR: missing implementation for get_cell_x_coord_indices",
                      __FILE__, __LINE__ );

   // routine will never reach this point
   return NULL;
}

unsigned const * get_cell_y_coord_indices (struct grid * grid, unsigned cell_index) {

   if (grid->vtable->get_cell_y_coord_indices != NULL)
      return grid->vtable->get_cell_y_coord_indices(grid, cell_index);
   else
      abort_message ( "ERROR: missing implementation for get_cell_y_coord_indices",
                      __FILE__, __LINE__ );

   // routine will never reach this point
   return NULL;
}

double get_corner_x_coord (struct grid * grid, unsigned corner_index) {

   if (grid->vtable->get_corner_x_coord != NULL)
      return grid->vtable->get_corner_x_coord(grid, corner_index);
   else
      abort_message ( "ERROR: missing implementation for get_corner_x_coord",
                      __FILE__, __LINE__ );

   // routine will never reach this point
   return -1;
}

double get_corner_y_coord (struct grid * grid, unsigned corner_index) {

   if (grid->vtable->get_corner_y_coord != NULL)
      return grid->vtable->get_corner_y_coord(grid, corner_index);
   else
      abort_message ( "ERROR: missing implementation for get_corner_y_coord",
                      __FILE__, __LINE__ );

   // routine will never reach this point
   return -1;
}

unsigned get_corner_x_coord_index (struct grid * grid, unsigned corner_index) {

   if (grid->vtable->get_corner_x_coord_index != NULL)
      return grid->vtable->get_corner_x_coord_index(grid, corner_index);
   else
      abort_message ( "ERROR: missing implementation for get_corner_x_coord_index",
                      __FILE__, __LINE__ );

   // routine will never reach this point
   return -1;
}

unsigned get_corner_y_coord_index (struct grid * grid, unsigned corner_index) {

   if (grid->vtable->get_corner_y_coord_index != NULL)
      return grid->vtable->get_corner_y_coord_index(grid, corner_index);
   else
      abort_message ( "ERROR: missing implementation for get_corner_y_coord_index",
                      __FILE__, __LINE__ );

   // routine will never reach this point
   return -1;
}

int get_aux_grid_cell(struct grid * grid, unsigned corner_index,
                      unsigned * cell_indices, enum edge_type * edge_type) {

   if (grid->vtable->get_aux_grid_cell != NULL)
      return grid->vtable->get_aux_grid_cell(grid, corner_index, cell_indices, edge_type);
   else
      abort_message ( "ERROR: missing implementation for get_aux_grid_cell",
                      __FILE__, __LINE__ );

   // routine will never reach this point
   return -1;
}

struct dep_list get_cell_neigh_dep_list(struct grid * grid) {

   if (grid->vtable->get_cell_neigh_dep_list != NULL)
      return grid->vtable->get_cell_neigh_dep_list(grid);
   else
      abort_message ( "ERROR: missing implementation for get_cell_neigh_dep_list",
                      __FILE__, __LINE__ );

   // routine will never reach this point
   exit(EXIT_FAILURE);
}

unsigned cell_covers_pole (unsigned num_corners, double * const corners_lon, 
                                                 double * const corners_lat) {

   double const tol = 1e-10;

   unsigned i;

   int sense;

   // test if any corner is directly on the pole
   for (i = 0; i < num_corners; ++i)
      if (fabs(M_PI_2 - fabs(corners_lat[i])) < tol ) return 1;
   sense = 0;

   // for all edges

   for (i = 0; i < num_corners-1; ++i) {

      if (get_angle(corners_lon[i+1], corners_lon[i]) >= 0.0)
         ++sense;
      else 
         --sense;
   }

   if (get_angle(corners_lon[0], corners_lon[num_corners-1]) >= 0.0)
      ++sense;
   else 
      --sense;

   return abs(sense) == (int)num_corners;
}

void get_boundary_corners (struct grid * grid, unsigned * bnd_corners,
                           unsigned * num_bnd_corners) {

   if (grid->vtable->get_boundary_corners != NULL)
      grid->vtable->get_boundary_corners(grid, bnd_corners, num_bnd_corners);
   else
      abort_message ( "ERROR: missing implementation get_boundary_corners",
                      __FILE__, __LINE__ );
}

struct grid * generate_cell_grid(struct grid * grid, double * coordinates_x, 
                                 double * coordinates_y) {

   if (grid->vtable->generate_cell_grid != NULL)
      return grid->vtable->generate_cell_grid(grid, coordinates_x, coordinates_y);
   else
      abort_message ( "ERROR: missing implementation for generate_cell_grid",
                      __FILE__, __LINE__ );

   return NULL;
}

void pack_grid(struct grid * grid, double ** dble_buf,
               unsigned dble_buf_offset, unsigned * dble_buf_data_size,
               unsigned * dble_buf_size, unsigned ** uint_buf,
               unsigned uint_buf_offset, unsigned * uint_buf_data_size,
               unsigned * uint_buf_size) {

   if (grid->vtable->pack_grid != NULL)
      grid->vtable->pack_grid(grid, dble_buf, dble_buf_offset, dble_buf_data_size,
                             dble_buf_size, uint_buf, uint_buf_offset,
                             uint_buf_data_size, uint_buf_size);
   else
      abort_message ( "ERROR: missing implementation for pack_grid",
                      __FILE__, __LINE__ );
}

struct grid * generate_subgrid(struct grid * grid, unsigned * selected_local_cell_ids,
                               unsigned num_local_cells, unsigned ** local_cell_ids,
                               unsigned ** local_corner_ids, unsigned ** local_edge_ids) {

   if (grid->vtable->generate_subgrid != NULL)
      return grid->vtable->generate_subgrid( grid, selected_local_cell_ids,
                                             num_local_cells, local_cell_ids,
                                             local_corner_ids, local_edge_ids);
   else
      abort_message ( "ERROR: missing implementation for generate_subgrid",
                      __FILE__, __LINE__ );

   return NULL;
}

struct grid_search * get_grid_search(struct grid * grid) {

   if (grid->vtable->get_grid_search != NULL)
      return grid->vtable->get_grid_search(grid);
   else
      abort_message("ERROR: missing implementation for get_grid_search",
                    __FILE__, __LINE__);

   return NULL;
}

void delete_grid(struct grid * grid) {

   if (grid == NULL)
      return;

   if (grid->vtable->delete != NULL)
      grid->vtable->delete(grid);
   else
      abort_message("ERROR: missing implementation for delete_grid", __FILE__, __LINE__);
}

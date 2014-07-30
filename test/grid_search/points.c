/**
 * @file points.c
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

#include "points.h"
#include "utils.h"

void init_points(struct points * points, struct grid * base_grid, enum location location,
                 double * coordinates_x, double * coordinates_y) {

   points->location = location;
   points->coordinates_x = coordinates_x;
   points->coordinates_y = coordinates_y;
   points->base_grid = base_grid;
   points->point_grid = NULL;
}

struct grid * get_base_grid(struct points * points) {

   return points->base_grid;
}

struct grid * get_point_grid(struct points * points) {

   if (points->point_grid == NULL) {

      switch (points->location) {
      
         case (CELL):
         {
            points->point_grid = generate_cell_grid(points->base_grid,
                                                    points->coordinates_x,
                                                    points->coordinates_y);
            break;
         }
         case (CORNER):
         {
            points->point_grid = copy_grid(points->base_grid);
            set_x_coords(points->point_grid, points->coordinates_x);
            set_y_coords(points->point_grid, points->coordinates_y);

            break;;
         }
         case (EDGE):
            abort_message ( "ERROR: get_point_grid: cannot generate point grid for location EDGE.", __FILE__, __LINE__ );
         default:
            abort_message ( "ERROR: get_point_grid: location must be one of CORNER/EDGE/CELL.", __FILE__, __LINE__ );
      };
   }

   return points->point_grid;
}

void get_coordinate_array_sizes (struct grid * grid, enum location location,
                                 unsigned * sizes) {

   switch (location) {

      case (CELL):
      {
         sizes[0] = get_size_cell_grid_x_coords(grid);
         sizes[1] = get_size_cell_grid_y_coords(grid);

         break;
      }
      case (CORNER):
      {
         sizes[0] = get_size_x_coords(grid);
         sizes[1] = get_size_y_coords(grid);

         break;
      }
      case (EDGE):
      {
         sizes[0] = get_num_grid_edges(grid);
         sizes[1] = sizes[0];

         break;
      }
   };
}

unsigned get_data_size(struct points points) {

   switch (points.location) {

      case (CELL):
         return get_num_grid_cells(points.base_grid);
      case (CORNER):
         return get_num_grid_corners(points.base_grid);
      case (EDGE):
         return get_num_grid_edges(points.base_grid);
      default:
         abort_message ( "ERROR: get_data_size: location must be one of CORNER/EDGE/CELL.", __FILE__, __LINE__ );
   };

   // routine will never reach this points
   return -1;
}

void free_points(struct points * points) {

   delete_grid(points->point_grid);
}

enum location get_location(int const location) {

  switch (location) {
  case (CELL):
  case (CORNER):
  case (EDGE):
    return (enum location) location;
  default:
    abort_message ( "ERROR: get_location: location must be one of CORNER/EDGE/CELL.", __FILE__, __LINE__ );
  };

  // routine will never reach this point
  return CORNER;
}

void get_point_coordinates (struct points * points, unsigned local_point_id,
                            double * coordinates) {

   switch (points->location) {

      case (CELL):
      case (CORNER):
      {
         struct grid * grid = get_point_grid(points);

         coordinates[0] = points->coordinates_x[get_corner_x_coord_index(grid, local_point_id)];
         coordinates[1] = points->coordinates_y[get_corner_y_coord_index(grid, local_point_id)];

         break;
      }
      case (EDGE):
      {
         coordinates[0] = points->coordinates_x[local_point_id];
         coordinates[1] = points->coordinates_y[local_point_id];

         break;
      }
      default:
         abort_message ("ERROR: get_point_coordinates: location must be one of CORNER/EDGE/CELL.",
                        __FILE__, __LINE__);
   };
}

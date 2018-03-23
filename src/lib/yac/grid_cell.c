/**
 * @file grid_cell.c
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
#include <stdio.h>

#include "grid_cell.h"
#include "utils.h"
#include "ensure_array_size.h"
#include "geometry.h"

void yac_init_grid_cell(struct grid_cell * cell) {

   cell->coordinates_x = NULL;
   cell->coordinates_y = NULL;
   cell->coordinates_xyz = NULL;
   cell->edge_type = NULL;
   cell->num_corners = 0;
   cell->array_size = 0;
}

static void ensure_grid_cell_size(
  struct grid_cell * cell, unsigned num_corners) {

  if (num_corners > cell->array_size) {

    cell->coordinates_x = realloc(cell->coordinates_x, num_corners *
                                  sizeof(*(cell->coordinates_x)));
    cell->coordinates_y = realloc(cell->coordinates_y, num_corners *
                                  sizeof(*(cell->coordinates_y)));
    cell->coordinates_xyz = realloc(cell->coordinates_xyz, 3 * num_corners *
                                    sizeof(*(cell->coordinates_xyz)));
    cell->edge_type = realloc(cell->edge_type, num_corners *
                              sizeof(*(cell->edge_type)));
    cell->array_size = num_corners;
  }
}

void yac_copy_grid_cell(struct grid_cell in_cell, struct grid_cell * out_cell) {

   ensure_grid_cell_size(out_cell, in_cell.num_corners);

   memcpy(out_cell->coordinates_x, in_cell.coordinates_x,
          in_cell.num_corners * sizeof(*(out_cell->coordinates_x)));
   memcpy(out_cell->coordinates_y, in_cell.coordinates_y,
          in_cell.num_corners * sizeof(*(out_cell->coordinates_y)));
   memcpy(out_cell->coordinates_xyz, in_cell.coordinates_xyz,
          3 * in_cell.num_corners * sizeof(*(out_cell->coordinates_xyz)));
   memcpy(out_cell->edge_type, in_cell.edge_type,
          in_cell.num_corners * sizeof(*(out_cell->edge_type)));
   out_cell->num_corners = in_cell.num_corners;
}

void yac_free_grid_cell(struct grid_cell * cell) {

   if (cell->coordinates_x != NULL) free(cell->coordinates_x);
   if (cell->coordinates_y != NULL) free(cell->coordinates_y);
   if (cell->coordinates_xyz != NULL) free(cell->coordinates_xyz);
   if (cell->edge_type != NULL) free(cell->edge_type);

   yac_init_grid_cell(cell);
}

static void set_triangle(
  struct grid_cell cell, struct grid_cell * triangle, size_t idx[3]) {

  ensure_grid_cell_size(triangle, 3);
  triangle->coordinates_x[0] = cell.coordinates_x[idx[0]];
  triangle->coordinates_y[0] = cell.coordinates_y[idx[0]];
  triangle->coordinates_xyz[0+0] = cell.coordinates_xyz[3*idx[0]+0];
  triangle->coordinates_xyz[0+1] = cell.coordinates_xyz[3*idx[0]+1];
  triangle->coordinates_xyz[0+2] = cell.coordinates_xyz[3*idx[0]+2];
  triangle->edge_type[0] = GREAT_CIRCLE;
  triangle->coordinates_x[1] = cell.coordinates_x[idx[1]];
  triangle->coordinates_y[1] = cell.coordinates_y[idx[1]];
  triangle->coordinates_xyz[3+0] = cell.coordinates_xyz[3*idx[1]+0];
  triangle->coordinates_xyz[3+1] = cell.coordinates_xyz[3*idx[1]+1];
  triangle->coordinates_xyz[3+2] = cell.coordinates_xyz[3*idx[1]+2];
  triangle->edge_type[1] = GREAT_CIRCLE;
  triangle->coordinates_x[2] = cell.coordinates_x[idx[2]];
  triangle->coordinates_y[2] = cell.coordinates_y[idx[2]];
  triangle->coordinates_xyz[6+0] = cell.coordinates_xyz[3*idx[2]+0];
  triangle->coordinates_xyz[6+1] = cell.coordinates_xyz[3*idx[2]+1];
  triangle->coordinates_xyz[6+2] = cell.coordinates_xyz[3*idx[2]+2];
  triangle->edge_type[2] = GREAT_CIRCLE;
  triangle->num_corners = 3;
}

void yac_triangulate_cell(
  struct grid_cell cell, unsigned start_corner, struct grid_cell * triangles) {

  switch (cell.num_corners) {
    case(0):
    case(1):
    case(2):
      yac_internal_abort_message(
        "ERROR(yac_triangulate_cell): number < 3", __FILE__, __LINE__ );
      break;
    case(3):

      if (start_corner == 0) {
        yac_copy_grid_cell(cell, triangles);
        return;
      }

      set_triangle(
        cell, triangles,
        (size_t [3]){(size_t)((start_corner + 0)%3),
                     (size_t)((start_corner + 1)%3),
                     (size_t)((start_corner + 2)%3)});
      return;
    case(4): {
      size_t idx[5] = {(size_t)((start_corner + 0)%cell.num_corners),
                       (size_t)((start_corner + 1)%cell.num_corners),
                       (size_t)((start_corner + 2)%cell.num_corners),
                       (size_t)((start_corner + 3)%cell.num_corners),
                       (size_t)((start_corner + 0)%cell.num_corners)};
      set_triangle(cell, triangles + 0, &(idx[0]));
      set_triangle(cell, triangles + 1, &(idx[2]));
      return;
    }
    case(6): {
      size_t idx[7] = {(size_t)((start_corner + 0)%cell.num_corners),
                       (size_t)((start_corner + 1)%cell.num_corners),
                       (size_t)((start_corner + 2)%cell.num_corners),
                       (size_t)((start_corner + 3)%cell.num_corners),
                       (size_t)((start_corner + 4)%cell.num_corners),
                       (size_t)((start_corner + 5)%cell.num_corners),
                       (size_t)((start_corner + 0)%cell.num_corners)};
      set_triangle(cell, triangles + 0, &(idx[0]));
      set_triangle(cell, triangles + 1, &(idx[2]));
      set_triangle(cell, triangles + 2, &(idx[4]));
      set_triangle(cell, triangles + 3, (size_t [3]){idx[0], idx[2], idx[4]});
      return;
    }
    case(5):
    default:

      for (size_t i = 0; i < cell.num_corners -2; ++i)
        set_triangle(
          cell, triangles + i,
          (size_t [3]){(size_t)(start_corner),
                       (size_t)((start_corner + i + 1)%cell.num_corners),
                       (size_t)((start_corner + i + 2)%cell.num_corners)});
      return;
  }
}

void yac_triangulate_cell_indices(
  unsigned const * cell_indices, unsigned num_corners, unsigned start_corner,
  unsigned * triangle_indices) {

  switch (num_corners) {
    case(0):
    case(1):
    case(2):
      yac_internal_abort_message(
        "ERROR(yac_triangulate_cell_indices): number < 3", __FILE__, __LINE__ );
      break;
    case(3):

      if (start_corner == 0) {
        memcpy(triangle_indices, cell_indices, 3 * sizeof(*cell_indices));
        return;
      }

      triangle_indices[0] = cell_indices[(start_corner + 0)%3];
      triangle_indices[1] = cell_indices[(start_corner + 1)%3];
      triangle_indices[2] = cell_indices[(start_corner + 2)%3];
      return;
    case(4): {
      triangle_indices[0] = cell_indices[start_corner];
      triangle_indices[1] = cell_indices[(start_corner + 1)%4];
      triangle_indices[2] = cell_indices[(start_corner + 2)%4];
      triangle_indices[3] = cell_indices[(start_corner + 2)%4];
      triangle_indices[4] = cell_indices[(start_corner + 3)%4];
      triangle_indices[5] = cell_indices[start_corner];
      return;
    }
    case(6): {
      triangle_indices[ 0] = cell_indices[start_corner];
      triangle_indices[ 1] = cell_indices[(start_corner + 1)%6];
      triangle_indices[ 2] = cell_indices[(start_corner + 2)%6];
      triangle_indices[ 3] = cell_indices[(start_corner + 2)%6];
      triangle_indices[ 4] = cell_indices[(start_corner + 3)%6];
      triangle_indices[ 5] = cell_indices[(start_corner + 4)%6];
      triangle_indices[ 6] = cell_indices[(start_corner + 4)%6];
      triangle_indices[ 7] = cell_indices[(start_corner + 5)%6];
      triangle_indices[ 8] = cell_indices[start_corner];
      triangle_indices[ 9] = cell_indices[start_corner];
      triangle_indices[10] = cell_indices[(start_corner + 2)%6];
      triangle_indices[11] = cell_indices[(start_corner + 4)%6];
      return;
    }
    case(5):
    default:

      for (size_t i = 0; i < num_corners-2; ++i) {
        triangle_indices[3*i+0] = cell_indices[start_corner];
        triangle_indices[3*i+1] = cell_indices[(start_corner + i + 1)%num_corners];
        triangle_indices[3*i+2] = cell_indices[(start_corner + i + 2)%num_corners];
      }
      return;
  }
}

#ifdef YAC_DEBUG_GRID_CELL
void print_grid_cell(FILE * stream, struct grid_cell cell, char * name) {

  char * out = NULL;
  unsigned out_array_size = 0;
  unsigned out_size = 0;

  if (name != NULL) {

      out_size = strlen(name) + 1 + 1 + 1;
      ENSURE_ARRAY_SIZE(out, out_array_size, out_size);

      strcpy(out, name);
      strcat(out, ":\n");
  }

  for (unsigned i = 0; i < cell.num_corners; ++i) {

      char buffer[1024];

      sprintf(buffer, "%d x %.16f y %.16f %s\n", i, cell.coordinates_x[i],
              cell.coordinates_y[i],
             (cell.edge_type[i] == LAT_CIRCLE)?("LAT_CIRCLE"):
            ((cell.edge_type[i] == LON_CIRCLE)?("LON_CIRCLE"):
             ("GREAT_CIRCLE")));

      out_size += strlen(buffer);

      ENSURE_ARRAY_SIZE(out, out_array_size, out_size);

      strcat(out, buffer);
  }

  if (out != NULL)
    fputs(out, stream);

  free(out);
}
#endif

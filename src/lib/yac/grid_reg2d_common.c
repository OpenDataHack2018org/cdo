/**
 * @file grid_reg2d_common.c
 *
 * @copyright Copyright  (C)  2017 Moritz Hanke <hanke@dkrz.de>
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

#include "dep_list.h"
#include "grid_unstruct.h"
#include "utils.h"

int yac_get_aux_grid_cell_reg2d_common(
  unsigned num_cells[2], unsigned cyclic, unsigned corner_index,
  unsigned * cell_indices) {

  if (!cyclic) {

    unsigned row = corner_index / (num_cells[0] + 1);
    unsigned column = corner_index - (num_cells[0] + 1) * row;
  
    if (row == 0 || row == num_cells[1] ||
        column == 0 || column == num_cells[0])
      return 1 == 0;
  
    cell_indices[0] = (row - 1) * num_cells[0] + column - 1;
    cell_indices[1] = (row - 1) * num_cells[0] + column;
    cell_indices[2] = row * num_cells[0] + column;
    cell_indices[3] = row * num_cells[0] + column - 1;

  } else {

    unsigned row = corner_index / num_cells[0];
    unsigned column = corner_index - num_cells[0] * row;

    if (row == 0 || row == num_cells[1])
      return 1 == 0;

    if (column == 0)
      cell_indices[0] = corner_index - 1;
    else
      cell_indices[0] = corner_index - num_cells[0] - 1;
    cell_indices[1] = corner_index - num_cells[0];
    cell_indices[2] = corner_index;
    if (column == 0)
      cell_indices[3] = corner_index + num_cells[0] - 1;
    else
      cell_indices[3] = corner_index - 1;
  } // if (!cyclic)

  return 1 == 1;
}

unsigned yac_get_num_corner_edges_reg2d_common(
  unsigned num_cells[2], unsigned cyclic, unsigned corner_index) {

  unsigned num_edges = 4;

  unsigned temp;

  if (!cyclic) {

    temp = corner_index / (num_cells[0] + 1);
    if (temp == 0 || temp == num_cells[1]) --num_edges;

    temp = corner_index - (num_cells[0] + 1) * temp;
    if (temp == 0 || temp == num_cells[0]) --num_edges;
  } else {

    temp = corner_index / num_cells[0];
    if (temp == 0 || temp == num_cells[1]) --num_edges;
  }

  return num_edges;
}

void yac_generate_cell_neigh_dep_reg2d_common(
  struct dep_list * cell_to_neigh, unsigned num_cells[2], unsigned cyclic[2]) {

   unsigned total_num_neighs;
   unsigned num_grid_cells;
   unsigned * num_neigh_per_cell;
   unsigned * cell_neigh_dependencies;
   unsigned * curr_cell_neigh_dep;
   unsigned curr_cell_index;

   unsigned i, j, k;

   total_num_neighs = 0;

   num_grid_cells = num_cells[0] * num_cells[1];
   num_neigh_per_cell = calloc (num_grid_cells, sizeof (num_neigh_per_cell[0]));

   // compute the number of total neighbours
   total_num_neighs = 0;
   // if there are inner cells (which always have 8 neighbours)
   if (num_cells[0] > 2 &&
       num_cells[1] > 2)
      total_num_neighs = (num_cells[0]-2) *
                         (num_cells[1]-2) * 8;
   // if the right and left boundary contains cells
   if (num_cells[1] > 2) {

      if (num_cells[0] > 1) {
         // if this boundary is cyclic
         if (cyclic[0])
            total_num_neighs += 2 * 8 * (num_cells[1]-2);
         else
            total_num_neighs += 2 * 5 * (num_cells[1]-2);
      } else
         total_num_neighs += 2 * (num_cells[1]-2);
   }
   // if the upper and lower boundary contains cells
   if (num_cells[0] > 2) {

      if (num_cells[1] > 1) {
         // if this boundary is cyclic
         if (cyclic[1])
            total_num_neighs += 2 * 8 * (num_cells[0]-2);
         else
            total_num_neighs += 2 * 5 * (num_cells[0]-2);
      } else
         total_num_neighs += 2 * (num_cells[0]-2);
   }
   // add the number of neighbours for the four corners
   if ((num_cells[0] > 1) && (num_cells[0] > 1)) {
      total_num_neighs += 4 * 3;
      if (cyclic[0]) total_num_neighs += 4 * 2;
      if (cyclic[1]) total_num_neighs += 4 * 2;
      if (cyclic[0] && cyclic[1])
         total_num_neighs += 4 * 1;
   } else if (!((num_cells[0] == 1) && (num_cells[1] == 1))) {
      total_num_neighs += 2;
      if ((num_cells[0] == 1) && cyclic[1])
         total_num_neighs += 2;
      if ((num_cells[1] == 1) && cyclic[0])
         total_num_neighs += 2;
   }

   //as an upper estimate we assume that there are 4 neighbours per cell
   cell_neigh_dependencies = malloc (total_num_neighs * sizeof(cell_neigh_dependencies[0]));

   enum position {
      LOW   = 1 << 0,
      LORI  = 1 << 1, // lower right
      RIGHT = 1 << 2,
      UPRI  = 1 << 3, // upper right
      UP    = 1 << 4,
      UPLE  = 1 << 5, // upper left
      LEFT  = 1 << 6, 
      LOLE  = 1 << 7, // lower left
      FULL  = 255
   };

   unsigned valid_neigh, temp_vaid_neigh;
   unsigned neigh_offset[8];

   // set the neighbours
   curr_cell_neigh_dep = cell_neigh_dependencies;
   for (i = 0; i < num_cells[1]; ++i) {

      neigh_offset[0] =     - num_cells[0];
      neigh_offset[1] =   1 - num_cells[0];
      neigh_offset[2] =   1;
      neigh_offset[3] =   1 + num_cells[0];
      neigh_offset[4] =       num_cells[0];
      neigh_offset[5] = - 1 + num_cells[0];
      neigh_offset[6] = - 1;
      neigh_offset[7] = - 1 - num_cells[0];

      temp_vaid_neigh = FULL;

      // if the grid is not cyclic in vertical direction
      if (!cyclic[1]) {
         // if we are at the lower boundary of the grid
         if (i == 0) temp_vaid_neigh &= ~(LOW + LOLE + LORI); // deactivates the 3 lower neighs
         // if we are at the upper boundary of the grid
         else if (i == num_cells[1]-1)
            temp_vaid_neigh &= ~(UP + UPRI + UPLE); // deactivates the 3 upper neighs
      } else {
         // if we are at the lower boundary of the grid
         if (i == 0) {
            neigh_offset[0] += num_grid_cells;
            neigh_offset[1] += num_grid_cells;
            neigh_offset[7] += num_grid_cells;
                // if we are at the upper boundary of the grid
         } else if (i == num_cells[0]-1) {
            neigh_offset[1] -= num_grid_cells;
            neigh_offset[2] -= num_grid_cells;
            neigh_offset[3] -= num_grid_cells;
         }
      }

      // if we only have one column of cells
      if (num_cells[0] == 1)
         temp_vaid_neigh &= ~(LORI + RIGHT + UPRI + LOLE + LEFT + UPLE);
            // deactivates left and right neighbours

      // if we only have one row of cells
      if (num_cells[1] == 1)
         temp_vaid_neigh &= ~(UPLE + UP + UPRI + LOLE + LOW + LORI);
            // deactivates lower and upper neighbours

      for (j = 0; j < num_cells[0]; ++j) {

         valid_neigh = temp_vaid_neigh;
         // if the grid is not cyclic in horizontal direction
         if (!cyclic[0]) {
            // if we are at the left boundary of the grid
            if (j == 0) valid_neigh &= ~(LEFT + LOLE + UPLE); // deactivates the 3 left neighs
            // if we are at the right boundary of the grid
            else if (j == num_cells[0]-1)
               valid_neigh &= ~(RIGHT + UPRI + LORI); // deactivates the 3 right neighs
         } else {
            // if we are at the left boundary of the grid
            if (j == 0) {
               neigh_offset[5] += num_cells[0];
               neigh_offset[6] += num_cells[0];
               neigh_offset[7] += num_cells[0];
             // if we are at the right boundary of the grid
            } else if (j == num_cells[0]-1) {
               neigh_offset[1] -= num_cells[0];
               neigh_offset[2] -= num_cells[0];
               neigh_offset[3] -= num_cells[0];
            }
         }

         curr_cell_index = j + i * num_cells[0];

         // for all neighbours
         for (k = 0; k < 8; ++k) {

            if (valid_neigh & (1 << k)) {

               ++num_neigh_per_cell[curr_cell_index];
               *(curr_cell_neigh_dep++) = curr_cell_index + neigh_offset[k];
            }
         }

         // reset neigh_offset array if necessary
         if (cyclic[0]) {
            if (j == 0) {
               neigh_offset[5] -= num_cells[0];
               neigh_offset[6] -= num_cells[0];
               neigh_offset[7] -= num_cells[0];
            } else if (j == num_cells[0]-1) {
               neigh_offset[5] += num_cells[0];
               neigh_offset[6] += num_cells[0];
               neigh_offset[7] += num_cells[0];
            }
         }
      }
   }

   yac_set_dependencies(
      cell_to_neigh, num_grid_cells, num_neigh_per_cell,
      cell_neigh_dependencies);
}

void yac_get_boundary_corners_reg2d_common(
  unsigned num_cells[2], unsigned cyclic, unsigned * bnd_corners,
  unsigned * num_bnd_corners) {

   unsigned i;

   if (!cyclic) {

      *num_bnd_corners = 2 * num_cells[0] +
                         2 * num_cells[1];

      for (i = 0; i <= num_cells[0]; ++i) {

         bnd_corners[2*i]   = i;
         bnd_corners[2*i+1] = i + (num_cells[0]+1) *
                              num_cells[1];
      }

      for (i = 1; i < num_cells[1]; ++i) {

         bnd_corners[2*(i+num_cells[0])]   = i * (num_cells[0]+1);
         bnd_corners[2*(i+num_cells[0])+1] =
            i * (num_cells[0]+1) + num_cells[0];
      }
   } else {

      *num_bnd_corners = 2 * num_cells[0];

      for (i = 0; i < num_cells[0]; ++i) {

         bnd_corners[i] = i;
         bnd_corners[num_cells[0]+i] =
            i + num_cells[0] * num_cells[1];
      }
   }
}

unsigned const * yac_get_corner_edges_reg2d_common(
  unsigned num_cells[2], unsigned cyclic, unsigned corner_index) {

   static unsigned corners[4];

   if (!cyclic) {

      unsigned row = corner_index / (num_cells[0] + 1);
      unsigned column = corner_index - (num_cells[0] + 1) * row;

      unsigned edge_idx = 0;

      if (row != 0)
         corners[edge_idx++] = corner_index - num_cells[0] - 1;

      if (column != 0)
         corners[edge_idx++] = corner_index - 1;

      if (column != num_cells[0])
         corners[edge_idx++] = corner_index + 1;

      if (row != num_cells[1])
         corners[edge_idx] = corner_index + num_cells[0] + 1;
   } else {

      unsigned row = corner_index / num_cells[0];
      unsigned column = corner_index - num_cells[0] * row;

      unsigned edge_idx = 0;

      if (row != 0)
         corners[edge_idx++] = corner_index - num_cells[0];

      if (column == 0)
         corners[edge_idx++] = corner_index - 1 + num_cells[0];
      else
         corners[edge_idx++] = corner_index - 1;

      if (column+1 == num_cells[0])
         corners[edge_idx++] = corner_index + 1 - num_cells[0];
      else
         corners[edge_idx++] = corner_index + 1;

      if (row != num_cells[1])
         corners[edge_idx] = corner_index + num_cells[0];
   } // if (!cyclic)

   return corners;
}

unsigned const * yac_get_corner_edge_indices_reg2d_common_cyclic(
  unsigned num_cells[2], unsigned row, unsigned column) {

  static unsigned edges[4];

  unsigned temp = 2 * (num_cells[0] * row + column);

  // upper edge
  if (row == num_cells[1]) edges[0] = -1;
  else {
    edges[0] = temp + 2;
    if (column + 1 == num_cells[0]) edges[0] -= 1;
  }

  // left edge
  if (column == 0) edges[1] = 1 + 2 * num_cells[0] * row;
  else {
    edges[1] = temp - 1;
    if (column == 1) edges[1] -= 1;
  }

  // lower edge
  if (row == 0) edges[2] = -1;
  else {
    edges[2] = temp + 2 - 2 * num_cells[0];
    if (column + 1 == num_cells[0]) edges[2] -= 1;
  }

  // right edge
  if (column + 1 == num_cells[0]) edges[3] = 1 + 2 * num_cells[0] * row;
  else {
    edges[3] = temp + 1;
    if (column == 0) edges[3] -= 1;
  }

  return edges;
}

unsigned const * yac_get_corner_edge_indices_reg2d_common_non_cyclic(
  unsigned num_cells[2], unsigned row, unsigned column) {

  unsigned temp = 2 * num_cells[0] + 1;

  static unsigned edges[4];

  // upper edge
  if (row == num_cells[1]) edges[0] = -1;
  else {
    edges[0] = temp * row + 2 * column + 1;
    if (column == num_cells[0]) edges[0] -= 1;
  }

  // left edge
  if (column == 0) edges[1] = -1;
  else edges[1] = temp * row + 2 * column - 2;

  // lower edge
  if (row == 0) edges[2] = -1;
  else {
    edges[2] = temp * (row - 1) + 2 * column + 1;
    if (column == num_cells[0]) edges[2] -= 1;
  }

  // right edge
  if (row == num_cells[0]) edges[3] = -1;
  else edges[3] = temp * row + 2 * column;

  return edges;
}

unsigned const * yac_get_cell_edge_indices_reg2d_common_cyclic(
   unsigned num_cells[2], unsigned row, unsigned column) {

   static unsigned edges[4];

   // left edge
   if (column+1 == num_cells[0])
      edges[0] = 1 + 2 * num_cells[0] * row + 2 * column;
   else
      edges[0] = 2 + 2 * num_cells[0] * row + 2 * column;

   // lower edge
   if (column+1 == num_cells[0])
      edges[1] = 1 + 2 * num_cells[0] * row;
   else if (column == 0)
      edges[1] = edges[0] - 2;
   else
      edges[1] = edges[0] - 1;

   // right edge
   if (column+1 == num_cells[0])
      edges[2] = 2 + 2 * num_cells[0] * row;
   else if (column+2 == num_cells[0])
      edges[2] = edges[0] + 1;
   else
      edges[2] = edges[0] + 2;

   // upper edge
   if (row+1 == num_cells[1]) {
      if (column == 0 || column+1 == num_cells[0])
         edges[3] = edges[1] + 2 * num_cells[0];
      else
         edges[3] = 2 * num_cells[0] * num_cells[1] + 1 + column;
   } else
      edges[3] = edges[1] + 2 * num_cells[0];

   return edges;
}

unsigned const * yac_get_cell_edge_indices_reg2d_common_non_cyclic(
   unsigned num_cells[2], unsigned row, unsigned column) {

   static unsigned edges[4];

   edges[0] = (2 * num_cells[0] + 1) * row + 2 * column;
   if (column+1 == num_cells[0])
      edges[1] = edges[0] + 2;
   else
      edges[1] = edges[0] + 3;
   if (row+1 == num_cells[1])
      edges[2] = (2 * num_cells[0] + 1) * (row + 1) + column;
   else
      edges[2] = edges[0] + 2 * num_cells[0] + 1;
   edges[3] = edges[0] + 1;

   return edges;
}

unsigned const * yac_get_corner_cell_indices_reg2d_common(
  unsigned num_cells[2], unsigned cyclic, unsigned corner_index) {

   static unsigned corner_cells[4];

   unsigned cell_idx = 0;

   if (!cyclic) {

      unsigned row = corner_index / (num_cells[0] + 1);
      unsigned column = corner_index - (num_cells[0] + 1) * row;

      if (row != 0) {

         if (column != 0)
            corner_cells[cell_idx++] = (row - 1) * num_cells[0] + column - 1;

         if (column != num_cells[0])
            corner_cells[cell_idx++] = (row - 1) * num_cells[0] + column;
      }

      if (row != num_cells[1]) {

         if (column != 0)
            corner_cells[cell_idx++] = row * num_cells[0] + column - 1;

         if (column != num_cells[0])
            corner_cells[cell_idx] = row * num_cells[0] + column;
      }
   } else {

      unsigned row = corner_index / num_cells[0];
      unsigned column = corner_index - num_cells[0] * row;

      if (row != 0) {

         if (column == 0)
            corner_cells[cell_idx++] = corner_index - 1;
         else
            corner_cells[cell_idx++] = corner_index - num_cells[0] - 1;

         corner_cells[cell_idx++] = corner_index - num_cells[0];
      }

      if (row != num_cells[1]) {

         if (column == 0)
            corner_cells[cell_idx++] = corner_index + num_cells[0] - 1;
         else
            corner_cells[cell_idx++] = corner_index - 1;

         corner_cells[cell_idx] = corner_index;
      }

   } // if (!cyclic)

   return corner_cells;
}

unsigned const * yac_get_cell_corner_indices_reg2d_common(
  unsigned num_cells[2], unsigned cyclic, unsigned cell_index) {

   static unsigned cell_corners[4];

   unsigned y_index = cell_index / num_cells[0];
   unsigned x_index = cell_index - y_index * num_cells[0];

   if (!cyclic) {

      cell_corners[0] =  y_index      * (num_cells[0] + 1) + x_index;
      cell_corners[1] =  y_index      * (num_cells[0] + 1) + x_index + 1;
      cell_corners[2] = (y_index + 1) * (num_cells[0] + 1) + x_index + 1;
      cell_corners[3] = (y_index + 1) * (num_cells[0] + 1) + x_index;
   } else {

      cell_corners[0] =  y_index      * num_cells[0] + x_index;
      if (x_index + 1 != num_cells[0]) {
         cell_corners[1] =  y_index      * num_cells[0] + x_index + 1;
         cell_corners[2] = (y_index + 1) * num_cells[0] + x_index + 1;
      } else {
         cell_corners[1] =  y_index      * num_cells[0];
         cell_corners[2] = (y_index + 1) * num_cells[0];
      }
      cell_corners[3] = (y_index + 1) * num_cells[0] + x_index;
   }

   return cell_corners;
}

struct grid * yac_generate_edge_grid_reg2d_common(
  unsigned num_cells[2], unsigned cyclic[2],
  double * coordinates_x, double * coordinates_y) {

  if (cyclic[1])
    yac_internal_abort_message(
      "ERROR(generate_edge_grid_reg2d): internal error\n", __FILE__, __LINE__);

  unsigned num_edge_grid_cells =
    // edge cells within reg2d cells
    num_cells[0] * num_cells[1] +
    // edge cell around every inner vertex
    (num_cells[0] - ((cyclic[0])?(0):(1))) * (num_cells[1] - 1);

  unsigned total_num_corners_per_edge_grid_cell =
    num_edge_grid_cells * 4;

  unsigned * num_corners_per_edge_grid_cell =
    malloc(num_edge_grid_cells * sizeof(*num_corners_per_edge_grid_cell));
  unsigned * corners_of_edge_grid_cells =
    malloc(total_num_corners_per_edge_grid_cell *
           sizeof(*corners_of_edge_grid_cells));

  total_num_corners_per_edge_grid_cell = 0;

  for (unsigned i = 0; i < num_edge_grid_cells; ++i)
    num_corners_per_edge_grid_cell[i] = 4;

  if (cyclic[0]) {
    // edge cells in each reg2d cell
    for (unsigned j = 0; j < num_cells[1]; ++j)
      for (unsigned i = 0; i < num_cells[0];
           ++i, total_num_corners_per_edge_grid_cell += 4)
        memcpy(
          corners_of_edge_grid_cells + total_num_corners_per_edge_grid_cell,
          yac_get_cell_edge_indices_reg2d_common_cyclic(num_cells, j, i),
          4 * sizeof(*corners_of_edge_grid_cells));
    // edge cells around each inner vertex
    for (unsigned j = 1; j < num_cells[1]; ++j)
      for (unsigned i = 0; i < num_cells[0];
           ++i, total_num_corners_per_edge_grid_cell += 4)
        memcpy(
          corners_of_edge_grid_cells + total_num_corners_per_edge_grid_cell,
          yac_get_corner_edge_indices_reg2d_common_cyclic(num_cells, j, i),
          4 * sizeof(*corners_of_edge_grid_cells));
  } else {
    // edge cells in each reg2d cell
    for (unsigned j = 0; j < num_cells[1]; ++j)
      for (unsigned i = 0; i < num_cells[0];
           ++i, total_num_corners_per_edge_grid_cell += 4)
        memcpy(
          corners_of_edge_grid_cells + total_num_corners_per_edge_grid_cell,
          yac_get_cell_edge_indices_reg2d_common_non_cyclic(num_cells, j, i),
          4 * sizeof(*corners_of_edge_grid_cells));
    // edge cells around each inner vertex
    for (unsigned j = 1; j < num_cells[1]; ++j)
      for (unsigned i = 1; i < num_cells[0];
           ++i, total_num_corners_per_edge_grid_cell += 4)
        memcpy(
          corners_of_edge_grid_cells + total_num_corners_per_edge_grid_cell,
          yac_get_corner_edge_indices_reg2d_common_non_cyclic(num_cells, j, i),
          4 * sizeof(*corners_of_edge_grid_cells));
  }

  struct dep_list edge_grid_cell_to_vertex;
  yac_set_dependencies(
    &edge_grid_cell_to_vertex, num_edge_grid_cells,
    num_corners_per_edge_grid_cell, corners_of_edge_grid_cells);

  unsigned num_edge_grid_edges =
    (num_cells[0] + ((cyclic[0])?0:1)) * num_cells[1] +
    num_cells[0] * (num_cells[1] + 1);

  struct grid * grid =
    yac_unstruct_grid_new(
      coordinates_x, coordinates_y, num_edge_grid_edges,
      edge_grid_cell_to_vertex);
  yac_unstruct_grid_set_cell_to_vertex_by_user(grid, 0);

  return grid;
}

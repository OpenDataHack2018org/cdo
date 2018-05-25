/**
 * @file grid_reg2d_common.h
 * @brief routines shared between grid_reg2d and grid_curve2d
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

#ifndef GRID_REG2D_COMMON_H
#define GRID_REG2D_COMMON_H

#include "grid.h"

int yac_get_aux_grid_cell_reg2d_common(
  unsigned num_cells[2], unsigned cyclic, unsigned corner_index,
  unsigned * cell_indices);

unsigned yac_get_num_corner_edges_reg2d_common(
  unsigned num_cells[2], unsigned cyclic, unsigned corner_index);

void yac_generate_cell_neigh_dep_reg2d_common(
  struct dep_list * cell_to_neigh, unsigned num_cells[2], unsigned cyclic[2]);

void yac_get_boundary_corners_reg2d_common(
  unsigned num_cells[2], unsigned cyclic, unsigned * bnd_corners,
  unsigned * num_bnd_corners);

unsigned const * yac_get_corner_edges_reg2d_common(
  unsigned num_cells[2], unsigned cyclic, unsigned corner_index);

unsigned const * yac_get_corner_edge_indices_reg2d_common_cyclic(
  unsigned num_cells[2], unsigned row, unsigned column);

unsigned const * yac_get_corner_edge_indices_reg2d_common_non_cyclic(
  unsigned num_cells[2], unsigned row, unsigned column);

unsigned const * yac_get_cell_edge_indices_reg2d_common_cyclic(
   unsigned num_cells[2], unsigned row, unsigned column);

unsigned const * yac_get_cell_edge_indices_reg2d_common_non_cyclic(
   unsigned num_cells[2], unsigned row, unsigned column);

unsigned const * yac_get_cell_corner_indices_reg2d_common(
  unsigned num_cells[2], unsigned cyclic, unsigned cell_index);

unsigned const * yac_get_corner_cell_indices_reg2d_common(
  unsigned num_cells[2], unsigned cyclic, unsigned corner_index);

struct grid * yac_generate_edge_grid_reg2d_common(
  unsigned num_cells[2], unsigned cyclic[2],
  double * coordinates_x, double * coordinates_y);

#endif // GRID_REG2D_COMMON_H
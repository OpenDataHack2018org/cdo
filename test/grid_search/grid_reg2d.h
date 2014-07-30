/**
 * @file grid_reg2d.h
 * @brief Initialisation for 2-dimensional fully regular grids
 *
 * grids with longitudinal coordinates aligned along constant latitudes
 * and latitudinal coordinates aligned along constant longitudes. 
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

#ifndef GRID_REG2D_H
#define GRID_REG2D_H

#include "grid.h"

struct grid * reg2d_grid_new(double * coordinates_x, double * coordinates_y,
                             unsigned const num_cells[2],
                             unsigned const cyclic[2]);

struct grid * unpack_reg2d_grid(double * dble_buf, unsigned * dble_buf_data_size,
                                unsigned * uint_buf, unsigned * uint_buf_data_size);

#endif // GRID_REG2D_H

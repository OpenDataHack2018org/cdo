/**
 * @file grid_search_utils.h
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

#ifndef GRID_SEARCH_UTILS_H
#define GRID_SEARCH_UTILS_H

#include "dep_list.h"
#include "grid_search.h"
#include "grid.h"
#include "geometry.h"

void grid_search_utils_do_point_search_p(struct grid_search * search,
                                         struct grid * search_grid_data,
                                         struct grid * grid_data,
                                         struct dep_list * target_to_src_points);

void grid_search_utils_do_point_search_p2(struct grid_search * search,
                                          struct grid * search_grid_data,
                                          double * x_coordinates,
                                          double * y_coordinates,
                                          unsigned num_points,
                                          struct dep_list * target_to_src_points);

void grid_search_utils_do_point_search_p3(struct grid_search * search,
                                          double * x_coordinates,
                                          double * y_coordinates,
                                          unsigned num_points,
                                          struct dep_list * target_to_src_points,
                                          struct points * points);

void grid_search_utils_do_point_search_p4 (struct grid_search * search,
                                           struct grid * search_grid_data,
                                           double x_coordinate,
                                           double y_coordinate,
                                           unsigned * n_points,
                                           unsigned * points_size,
                                           unsigned ** points);

#endif // GRID_SEARCH_UTILS_H

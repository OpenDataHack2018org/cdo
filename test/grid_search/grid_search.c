/**
 * @file grid_search.c
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
#include "grid_search.h"
#include "utils.h"

void do_cell_search (struct grid_search * search, struct grid * grid_data,
                     struct dep_list * tgt_to_src_cells) {

   if (search->vtable->do_cell_search == NULL)
      abort_message("ERROR: routine not implemented: do_cell_search",
                    __FILE__, __LINE__);

   search->vtable->do_cell_search(search, grid_data, tgt_to_src_cells);
}

void do_cell_search_single (struct grid_search * search,
                            struct grid_cell grid_cell,
                            unsigned * n_cells, unsigned * cells_size,
                            unsigned ** cells) {

   if (search->vtable->do_cell_search_single == NULL)
      abort_message("ERROR: routine not implemented: do_cell_search_single",
                    __FILE__, __LINE__);

   search->vtable->do_cell_search_single(search, grid_cell, n_cells,
                                         cells_size, cells);
}

void do_point_search_c (struct grid_search * search, struct grid * grid_data,
                        struct dep_list * tgt_to_src_cells) {

   if (search->vtable->do_point_search_c == NULL)
      abort_message("ERROR: routine not implemented: do_point_search_c",
                    __FILE__, __LINE__);

   search->vtable->do_point_search_c(search, grid_data, tgt_to_src_cells);
}

void do_point_search_c2 (struct grid_search * search, double * x_coordinates,
                         double * y_coordinates, unsigned num_points,
                         struct dep_list * tgt_to_src_cells) {

   if (search->vtable->do_point_search_c2 == NULL)
      abort_message("ERROR: routine not implemented: do_point_search_c2",
                    __FILE__, __LINE__);

   search->vtable->do_point_search_c2(search, x_coordinates, y_coordinates,
                                     num_points, tgt_to_src_cells);
}

void do_point_search_p (struct grid_search * search, struct grid * grid_data,
                        struct dep_list * target_to_src_points) {

   if (search->vtable->do_point_search_p == NULL)
      abort_message("ERROR: routine not implemented: do_point_search_p",
                    __FILE__, __LINE__);

   search->vtable->do_point_search_p(search, grid_data, target_to_src_points);
}

void do_point_search_p2 (struct grid_search * search, double * x_coordinates,
                         double * y_coordinates, unsigned num_points,
                         struct dep_list * target_to_src_points) {

   if (search->vtable->do_point_search_p2 == NULL)
      abort_message("ERROR: routine not implemented: do_point_search_p2",
                    __FILE__, __LINE__);

   search->vtable->do_point_search_p2(search, x_coordinates, y_coordinates, 
                                     num_points, target_to_src_points);
}

void do_point_search_p3 (struct grid_search * search, double * x_coordinates,
                         double * y_coordinates, unsigned num_points,
                         struct dep_list * target_to_src_points,
                         struct points * points) {

   if (search->vtable->do_point_search_p3 == NULL)
      abort_message("ERROR: routine not implemented: do_point_search_p3",
                    __FILE__, __LINE__);

   search->vtable->do_point_search_p3(search, x_coordinates, y_coordinates, 
                                      num_points, target_to_src_points, points);
}

void do_point_search_p4 (struct grid_search * search, double x_coordinate,
                         double y_coordinate, unsigned * n_points,
                         unsigned * points_size, unsigned ** points) {

   if (search->vtable->do_point_search_p4 == NULL)
      abort_message("ERROR: routine not implemented: do_point_search_p4",
                    __FILE__, __LINE__);

   search->vtable->do_point_search_p4(search, x_coordinate, y_coordinate,
                                      n_points, points_size, points);
}

void delete_grid_search (struct grid_search * search) {

   if (search->vtable->delete_grid_search == NULL)
      abort_message("ERROR: routine not implemented: delete_grid_search",
                    __FILE__, __LINE__);

   search->vtable->delete_grid_search(search);
}

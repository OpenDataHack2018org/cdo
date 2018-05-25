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
#include "grid_search.h"
#include "utils.h"

void yac_do_cell_search (struct grid_search * search, struct grid * grid_data,
                         struct dep_list * tgt_to_src_cells) {

   if (search->vtable->do_cell_search == NULL)
      yac_internal_abort_message("ERROR: routine not implemented: do_cell_search",
                                 __FILE__, __LINE__);

   search->vtable->do_cell_search(search, grid_data, tgt_to_src_cells);
}

void yac_do_cell_search_single (struct grid_search * search,
                                struct grid_cell grid_cell,
                                unsigned * n_cells, unsigned * cells_size,
                                unsigned ** cells) {

   if (search->vtable->do_cell_search_single == NULL)
      yac_internal_abort_message("ERROR: routine not implemented: do_cell_search_single",
                                 __FILE__, __LINE__);

   search->vtable->do_cell_search_single(search, grid_cell, n_cells,
                                         cells_size, cells);
}

void yac_do_point_search_c (struct grid_search * search, struct grid * grid_data,
                            struct dep_list * tgt_to_src_cells) {

   if (search->vtable->do_point_search_c == NULL)
      yac_internal_abort_message("ERROR: routine not implemented: do_point_search_c",
                                 __FILE__, __LINE__);

   search->vtable->do_point_search_c(search, grid_data, tgt_to_src_cells);
}

void yac_do_point_search_c2 (struct grid_search * search,
                             double (*coordinates_xyz)[3], unsigned num_points,
                             struct dep_list * tgt_to_src_cells) {

   if (search->vtable->do_point_search_c2 == NULL)
      yac_internal_abort_message("ERROR: routine not implemented: do_point_search_c2",
                                 __FILE__, __LINE__);

   search->vtable->do_point_search_c2(
      search, coordinates_xyz, num_points, tgt_to_src_cells);
}

void yac_do_point_search_c3 (struct grid_search * search,
                             double (*coordinates_xyz)[3], unsigned num_points,
                             struct dep_list * tgt_to_src_cells,
                             struct points * points) {

   if (search->vtable->do_point_search_c3 == NULL)
      yac_internal_abort_message("ERROR: routine not implemented: do_point_search_c3",
                                 __FILE__, __LINE__);

   search->vtable->do_point_search_c3(
      search, coordinates_xyz, num_points, tgt_to_src_cells, points);
}

void yac_do_point_search_p (struct grid_search * search, struct grid * grid_data,
                            struct dep_list * target_to_src_points) {

   if (search->vtable->do_point_search_p == NULL)
      yac_internal_abort_message("ERROR: routine not implemented: do_point_search_p",
                                 __FILE__, __LINE__);

   search->vtable->do_point_search_p(search, grid_data, target_to_src_points);
}

void yac_do_point_search_p2 (struct grid_search * search,
                             double (*coordinates_xyz)[3], unsigned num_points,
                             struct dep_list * target_to_src_points) {

   if (search->vtable->do_point_search_p2 == NULL)
      yac_internal_abort_message("ERROR: routine not implemented: do_point_search_p2",
                                 __FILE__, __LINE__);

   search->vtable->do_point_search_p2(
      search, coordinates_xyz, num_points, target_to_src_points);
}

void yac_do_point_search_p3 (struct grid_search * search,
                             double (*coordinates_xyz)[3], unsigned num_points,
                             struct dep_list * target_to_src_points,
                             struct points * points) {

   if (search->vtable->do_point_search_p3 == NULL)
      yac_internal_abort_message("ERROR: routine not implemented: do_point_search_p3",
                                 __FILE__, __LINE__);

   search->vtable->do_point_search_p3(
      search, coordinates_xyz, num_points, target_to_src_points, points);
}

void yac_do_point_search_p4 (struct grid_search * search,
                             double coordinate_xyz[3], unsigned * n_points,
                             unsigned * points_size, unsigned ** points) {

   if (search->vtable->do_point_search_p4 == NULL)
      yac_internal_abort_message("ERROR: routine not implemented: do_point_search_p4",
                                 __FILE__, __LINE__);

   search->vtable->do_point_search_p4(
      search, coordinate_xyz, n_points, points_size, points);
}

void yac_do_bnd_circle_search (struct grid_search * search,
                               struct bounding_circle * bnd_circles,
                               unsigned n, struct dep_list * bnd_to_cells) {

   if (search->vtable->do_bnd_circle_search == NULL)
      yac_internal_abort_message("ERROR: routine not implemented: do_bnd_circle_search",
                                 __FILE__, __LINE__);

   search->vtable->do_bnd_circle_search(search, bnd_circles, n, bnd_to_cells);
}

void yac_delete_grid_search (struct grid_search * search) {

   if (search->vtable->delete_grid_search == NULL)
      yac_internal_abort_message("ERROR: routine not implemented: delete_grid_search",
                                 __FILE__, __LINE__);

   search->vtable->delete_grid_search(search);
}

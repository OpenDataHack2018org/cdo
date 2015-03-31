/**
 * @file find_overlapping_cells.c
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

#include "search.h"
#include <stdlib.h>
#include "ensure_array_size.h"

void yac_find_overlapping_cells_s (struct grid_cell src_cell,
                                   struct bounding_circle src_bnd_circle,
                                   struct grid * tgt_grid,
                                   unsigned const * initial_dep,
                                   unsigned num_initial_deps, unsigned ** deps,
                                   unsigned * deps_size, unsigned * num_deps,
                                   unsigned src_index,
                                   unsigned * tgts_already_touched,
                                   unsigned ** stack, unsigned * stack_size) {

   struct bounding_circle tgt_bnd_circle;

   unsigned tgt_index;
   struct grid_cell tgt_grid_cell;

   unsigned num_elements_on_stack;
   unsigned const * curr_tgt_neighs;

   // allocates memory and initialises it with 0
   *num_deps = 0;

   yac_init_grid_cell (&tgt_grid_cell);

   struct dep_list tgt_cell_neigh_dep = yac_get_cell_neigh_dep_list(tgt_grid);

   //check whether the current stack can hold the initial number of tgt cells
   ENSURE_ARRAY_SIZE(*stack, *stack_size, num_initial_deps);
   // add initial results to stack
   num_elements_on_stack = 0;
   for (int i = 0; i < num_initial_deps; ++i) {

      if (tgts_already_touched[initial_dep[i]] != src_index + 1) {
         (*stack)[num_elements_on_stack++] = initial_dep[i];
         tgts_already_touched[initial_dep[i]] = src_index + 1;
      }
   }

   // while there are still elements on the stack
   while (num_elements_on_stack > 0) {

      // pop top element of stack
      tgt_index = (*stack)[--num_elements_on_stack];

      yac_get_grid_cell2 (tgt_grid, tgt_index, &tgt_grid_cell, &tgt_bnd_circle);

      // check overlap between src and tgt cell
      // TODO: this check needs to be revisited because it has problems at the
      //       poles and does not take different properties of regular and
      //       unstructured grids into account
      if (yac_check_overlap_cells2(src_cell, src_bnd_circle,
                                   tgt_grid_cell, tgt_bnd_circle)) {

         //check whether the current dependency list can hold another dependency
         ENSURE_ARRAY_SIZE(*deps, *deps_size, *num_deps+1);

         // set dependency
         (*deps)[*num_deps] = tgt_index;
         ++*num_deps;

         // check size of stack
         ENSURE_ARRAY_SIZE(*stack, *stack_size, num_elements_on_stack + 
                           tgt_cell_neigh_dep.num_deps_per_element[tgt_index]);
         // add neighbours to stack
         curr_tgt_neighs = yac_get_dependencies_of_element(tgt_cell_neigh_dep, tgt_index);
         for (int i = 0; i < tgt_cell_neigh_dep.num_deps_per_element[tgt_index]; ++i) {

            if (tgts_already_touched[curr_tgt_neighs[i]] != src_index + 1) {
               (*stack)[num_elements_on_stack++] = curr_tgt_neighs[i];
               tgts_already_touched[curr_tgt_neighs[i]] = src_index + 1;
            }
         }
      } // if (check_overlap_cells(src_grid_cell, tgt_grid_cell))
   } // while (num_elements_on_stack > 0)

   *deps = realloc (*deps, *num_deps * sizeof (**deps));
   *deps_size = *num_deps;

   yac_free_grid_cell(&tgt_grid_cell);
}

void yac_find_overlapping_cells (struct grid * src_grid, struct grid * tgt_grid,
                                 struct dep_list initial_src_to_tgt_dep,
                                 struct dep_list * src_to_tgt_dep) {

   unsigned src_index;
   unsigned num_tgt_cells, num_src_cells;
   struct grid_cell src_grid_cell;
   unsigned num_total_deps;

   unsigned * stack = NULL;
   unsigned stack_size = 0;

   unsigned * curr_src_tgt_dependencies = NULL;
   unsigned curr_src_tgt_dependencies_size = 0;

   unsigned * num_tgt_per_src_cell;
   unsigned * src_tgt_dependencies = NULL;
   unsigned curr_dependency_size = 0;

   unsigned * tgts_already_touched;

   num_tgt_cells = yac_get_num_grid_cells (tgt_grid);
   num_src_cells = yac_get_num_grid_cells (src_grid);

   // allocates memory and initialises it with 0
   num_tgt_per_src_cell = calloc (num_src_cells, sizeof (num_tgt_per_src_cell[0]));

   tgts_already_touched = calloc (num_tgt_cells, sizeof (tgts_already_touched[0]));

   num_total_deps = 0;

   yac_init_grid_cell (&src_grid_cell);

   // for all src cells
   for (src_index = 0; src_index < num_src_cells; ++src_index) {

      struct bounding_circle src_bnd_circle;
      yac_get_grid_cell2 (src_grid, src_index, &src_grid_cell, &src_bnd_circle);

      // search for overlapping cells
      yac_find_overlapping_cells_s(src_grid_cell, src_bnd_circle, tgt_grid,
         yac_get_dependencies_of_element(initial_src_to_tgt_dep, src_index),
         initial_src_to_tgt_dep.num_deps_per_element[src_index],
         &curr_src_tgt_dependencies, &curr_src_tgt_dependencies_size,
         num_tgt_per_src_cell + src_index, src_index, tgts_already_touched,
         &stack, &stack_size);

      //check whether the current dependency list can hold the dependencies of
      // the current source cell
      ENSURE_ARRAY_SIZE(src_tgt_dependencies, curr_dependency_size,
                        num_total_deps + num_tgt_per_src_cell[src_index]);

      for (int i = 0; i < num_tgt_per_src_cell[src_index]; ++i)
         src_tgt_dependencies[num_total_deps++] = curr_src_tgt_dependencies[i];

   } // for (src_index = 0; src_index < num_src_cells; ++src_index)

   src_tgt_dependencies = realloc (src_tgt_dependencies,
      num_total_deps * sizeof (src_tgt_dependencies[0]));

   yac_init_dep_list(src_to_tgt_dep);
   yac_set_dependencies(src_to_tgt_dep, num_src_cells, num_tgt_per_src_cell,
                        src_tgt_dependencies);

   yac_free_grid_cell(&src_grid_cell);
   free(curr_src_tgt_dependencies);
   free(stack);
   free(tgts_already_touched);
}

/**
 * @file test_point_search.c
 *
 * @copyright Copyright  (C)  2013 DKRZ, MPI-M
 *
 * @version 1.0
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler  <rene.redler@mpimet.mpg.de>
 *
 */
/*
 * Keywords:
 *
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler  <rene.redler@mpimet.mpg.de>
 *
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
 *
 */

#include "tests.h"
#include "grid_search.h"
#include "bucket_search.h"
#include "sphere_part.h"
#include "grid.h"
#include "grid_reg2d.h"
#include "grid_unstruct.h"

typedef struct grid_search * (*grid_search_constructor)(struct grid *);

int main (void) {

   grid_search_constructor grid_search_construct[] =
      {bucket_search_new, sphere_part_search_new};
   unsigned num_grid_search_construct = sizeof(grid_search_construct) /
                                        sizeof(grid_search_construct[0]);

   //---------------
   // setup
   //---------------
   unsigned num_corners_per_cell[4] = {4,4,4,4};
   unsigned cell_vertices[16] = {0,1,4,3, 1,2,5,4, 3,4,7,6, 4,5,8,7}; 
   double vertex_coordinates_x[2][9] = {{0,1,2, 0,1,2, 0,1,2},
                                        {0.5,1.5,2.5, 0.5,1.5,2.5, 0.5,1.5,2.5}};
   double vertex_coordinates_y[2][9] = {{0,0,0, 1,1,1, 2,2,2},
                                        {0.5,0.5,0.5, 1.5,1.5,1.5, 2.5,2.5,2.5}};

   for ( unsigned i = 0; i < 2; ++i ) {
      for ( unsigned j = 0; j < 9; ++j ) {
         vertex_coordinates_x[i][j] *= rad;
         vertex_coordinates_y[i][j] *= rad;
      }
   }

   struct grid * grid_data_a, * grid_data_b;
   struct dep_list cell_to_vertex_a, cell_to_vertex_b;

   { // set up grid data
   init_dep_list(&cell_to_vertex_a);
   set_dependencies(&cell_to_vertex_a, 4, num_corners_per_cell, cell_vertices);

   grid_data_a = unstruct_grid_new(vertex_coordinates_x[0], vertex_coordinates_y[0], 
                                   9, cell_to_vertex_a);

   init_dep_list(&cell_to_vertex_b);
   set_dependencies(&cell_to_vertex_b, 4, num_corners_per_cell, cell_vertices);

   grid_data_b = unstruct_grid_new(vertex_coordinates_x[1], vertex_coordinates_y[1], 
                                   9, cell_to_vertex_b);
   }

   //---------------
   // testing
   //---------------
   { // test do_point_search

      for (int i = 0; i < num_grid_search_construct; ++i) {

         struct grid_search * search;
         struct dep_list deps;

         search = grid_search_construct[i](grid_data_a);

         do_point_search_c(search, grid_data_b, &deps);

         {
            unsigned ref_num_deps[9] = {1,1,0,1,1,0,0,0,0};
            unsigned ref_deps[9] = {0,1,-1,2,3,-1,-1,-1,-1};
            unsigned i;

            if (get_total_num_dependencies(deps) != 4)
               PUT_ERR("error in do_point_search_c(get_total_num_dependencies(deps))\n");

            for (i = 0; i < 9; ++i)
               if (deps.num_deps_per_element[i] != ref_num_deps[i])
                  PUT_ERR("error in do_point_search_c(deps.num_deps_per_element[i])\n");

            for (i = 0; i < 9; ++i) {

               if (deps.num_deps_per_element[i] == 0) continue;

               if (get_dependencies_of_element(deps,i)[0] !=
                   ref_deps[i])
                  PUT_ERR("error in do_point_search_c(get_dependencies_of_element)\n");
            }
         }

         free_dep_list(&deps);
         delete_grid_search(search);
      }
   }

   { // test do_point_search_p

      for (int i = 0; i < num_grid_search_construct; ++i) {

         struct grid_search * search;
         struct dep_list deps;

         search = grid_search_construct[i](grid_data_a);

         do_point_search_p(search, grid_data_b, &deps);

         {
            unsigned ref_num_deps[9] = {4,4,0,4,4,0,0,0,0};
            unsigned ref_deps[9][4] = {{0,1,4,3}, {1,2,5,4}, {0,0,0,0},
                                       {3,4,7,6}, {4,5,8,7}, {0,0,0,0},
                                       {0,0,0,0}, {0,0,0,0}, {0,0,0,0}};
            unsigned i, j, k;

            if (get_total_num_dependencies(deps) != 16)
               PUT_ERR("error in do_point_search_p(get_total_num_dependencies(deps))\n");

            for (i = 0; i < 9; ++i)
               if (deps.num_deps_per_element[i] != ref_num_deps[i])
                  PUT_ERR("error in do_point_search_p(deps.num_deps_per_element[i])\n");

            unsigned const * curr_src_corners;
            unsigned num_matches;
            for (i = 0; i < 9; ++i) {

               curr_src_corners = get_dependencies_of_element(deps,i);

               num_matches = 0;

               for (j = 0; j < ref_num_deps[i]; ++j)
                  for (k = 0; k < ref_num_deps[i]; ++k)
                     if (ref_deps[i][j] == curr_src_corners[k])
                        ++num_matches;

               if (num_matches != ref_num_deps[i])
                  PUT_ERR("error in do_point_search_p(get_dependencies_of_element)\n");
            }
         }

         free_dep_list(&deps);
         delete_grid_search(search);
      }
   }

   delete_grid(grid_data_a);
   delete_grid(grid_data_b);

   { // test do_point_search_p3

      //-------------------
      // set up
      //-------------------

      // set up points structure

      double grid_coord_x[] = {0,1,2,3,4,5};
      double grid_coord_y[] = {0,1,2,3,4,5};
      double points_coord_x[] = {0.5, 1.5, 2.5, 3.5, 4.5};
      double points_coord_y[] = {0.5, 1.5, 2.5, 3.5, 4.5};

      for ( unsigned i = 0; i < 6; ++i ) {
         grid_coord_x[i] *= rad;
         grid_coord_y[i] *= rad;
      }

      for ( unsigned i = 0; i < 5; ++i ) {
         points_coord_x[i] *= rad;
         points_coord_y[i] *= rad;
       }

      unsigned num_cells[2] = {5,5};
      unsigned cyclic[2] = {0,0};
      struct grid * reg_grid;
      struct points points;

      reg_grid = reg2d_grid_new(grid_coord_x, grid_coord_y, num_cells, cyclic);
      init_points(&points, reg_grid, CELL, points_coord_x, points_coord_y);

      // define points that are to be seached for

      double search_points_coord_x[] = {4.75,1.25,3.5,0.25,2.25,0.75,3.25,5.25};
      double search_points_coord_y[] = {0.75,1.25,1.5, 2.5,2.75,4.25,4.25,4.25};

      for ( unsigned i = 0; i < 8; ++i ) {
         search_points_coord_x[i] *= rad;
         search_points_coord_y[i] *= rad;
      }

      unsigned num_search_points = sizeof(search_points_coord_x) /
                                   sizeof(search_points_coord_x[0]);

      // init search data structure

      for (int i = 0; i < 2; ++i) {

         struct grid_search * search = grid_search_construct[i](reg_grid);

         //--------------
         // testing
         //--------------

         struct dep_list search_results;

         // do search
         do_point_search_p3 (search, search_points_coord_x, search_points_coord_y,
                             num_search_points, &search_results, &points);

         {
            unsigned ref_num_deps[8] = {0,4,4,0,4,4,4,0};
            unsigned ref_deps[8][4] = {{-1}, {0,1,6,5}, {2,3,8,7},{-1},
                                       {11,12,17,16}, {15,16,21,20},{17,18,23,22},{-1}};
            unsigned i, j, k;

            if (get_total_num_dependencies(search_results) != 20)
               PUT_ERR("error in do_point_search_p3(get_total_num_dependencies(search_results))\n");

            for (i = 0; i < 7; ++i)
               if (search_results.num_deps_per_element[i] != ref_num_deps[i])
                  PUT_ERR("error in do_point_search_p3(search_results.num_deps_per_element[i])\n");

            unsigned const * curr_src_corners;
            unsigned num_matches;
            for (i = 0; i < 7; ++i) {

               curr_src_corners = get_dependencies_of_element(search_results,i);

               num_matches = 0;

               for (j = 0; j < ref_num_deps[i]; ++j)
                  for (k = 0; k < ref_num_deps[i]; ++k)
                     if (ref_deps[i][j] == curr_src_corners[k])
                        ++num_matches;

               if (num_matches != ref_num_deps[i])
                  PUT_ERR("error in do_point_search_p3(get_dependencies_of_element)\n");
            }
         }

         free_dep_list(&search_results);
         delete_grid_search(search);
      }

      delete_grid(reg_grid);
      free_points(&points);
   }

   { // test do_point_search_p4

      //-------------------
      // set up
      //-------------------

      // set up points structure

      double grid_coord_x[] = {0.5*rad, 1.5*rad, 2.5*rad, 3.5*rad, 4.5*rad};
      double grid_coord_y[] = {0.5*rad, 1.5*rad, 2.5*rad, 3.5*rad, 4.5*rad};

      unsigned num_cells[2] = {4,4};
      unsigned cyclic[2] = {0,0};
      struct grid * reg_grid;

      reg_grid = reg2d_grid_new(grid_coord_x, grid_coord_y, num_cells, cyclic);

      // define points that are to be seached for

      double search_points_coord_x[] = {4.75*rad,1.25*rad,3.5*rad,0.25*rad,
                                        2.25*rad,0.75*rad,3.25*rad,5.25*rad};
      double search_points_coord_y[] = {0.75*rad,1.25*rad,1.5*rad, 2.5*rad,
                                        2.75*rad,4.25*rad,4.25*rad,4.25*rad};

      unsigned num_search_points = sizeof(search_points_coord_x) /
                                   sizeof(search_points_coord_x[0]);

      for (int i = 0; i < num_grid_search_construct; ++i) {

         // init search data structure
         struct grid_search * search = grid_search_construct[i](reg_grid);

         //--------------
         // testing
         //--------------

         unsigned * deps = NULL;
         unsigned deps_size = 0;

         for (int i = 0; i < num_search_points; ++i) {

            unsigned num_deps;

            // do search
            do_point_search_p4 (search, search_points_coord_x[i],
                                search_points_coord_y[i], &num_deps, &deps_size,
                                &deps);

            {
               unsigned ref_num_deps[8] = {0,4,4,0,4,4,4,0};
               unsigned ref_deps[8][4] = {{-1}, {0,1,6,5}, {2,3,8,7},{-1},
                                          {11,12,17,16}, {15,16,21,20},{17,18,23,22},{-1}};

               if (num_deps != ref_num_deps[i])
                  PUT_ERR("error in do_point_search_p3(search_results.num_deps_per_element[i])\n");

               unsigned num_matches = 0;

               for (unsigned j = 0; j < ref_num_deps[i]; ++j)
                  for (unsigned k = 0; k < ref_num_deps[i]; ++k)
                     if (ref_deps[i][j] == deps[k])
                        ++num_matches;

               if (num_matches != ref_num_deps[i])
                  PUT_ERR("error in do_point_search_p3(get_dependencies_of_element)\n");
            }
         }

         free(deps);
         delete_grid_search(search);
      }

      delete_grid(reg_grid);
   }

   return TEST_EXIT_CODE;
}

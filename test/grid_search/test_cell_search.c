/**
 * @file test_cell_search.c
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

#include <string.h>
#include <time.h>

#include "tests.h"
#include "grid_search.h"
#include "reg2d_search.h"
#include "bucket_search.h"
#include "sphere_part.h"
#include "grid.h"
#include "grid_reg2d.h"
#include "grid_unstruct.h"
#include "ensure_array_size.h"


typedef struct grid_search * (*grid_search_constructor)(struct grid *);

void check_dep_lists(struct dep_list * lists);

int main (void) {

   grid_search_constructor grid_search_construct[] =
      {reg2d_search_new, yac_bucket_search_new, yac_sphere_part_search_new};
   unsigned num_grid_search_construct = sizeof(grid_search_construct) /
                                        sizeof(grid_search_construct[0]);

   {
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

      struct dep_list deps;

      // set up grid data
      yac_init_dep_list(&cell_to_vertex_a);
      yac_set_dependencies(&cell_to_vertex_a, 4, num_corners_per_cell, cell_vertices);

      grid_data_a = yac_unstruct_grid_new(vertex_coordinates_x[0], vertex_coordinates_y[0], 
					  9, cell_to_vertex_a);

      yac_init_dep_list(&cell_to_vertex_b);
      yac_set_dependencies(&cell_to_vertex_b, 4, num_corners_per_cell, cell_vertices);

      grid_data_b = yac_unstruct_grid_new(vertex_coordinates_x[1], vertex_coordinates_y[1], 
					  9, cell_to_vertex_b);

      //---------------
      // testing
      //---------------

      for (int i = 0; i < num_grid_search_construct; ++i) {

         struct grid_search * search = grid_search_construct[i](grid_data_a);
      
         yac_do_cell_search(search, grid_data_b, &deps);

         // test whether the search returned something in deps
         if (deps.num_elements == 0 || deps.num_deps_per_element == NULL || 
             deps.dependencies == NULL)
            PUT_ERR("do_cell_search returned no result\n\n")

         // test whether the number of found cells is correct
         if (deps.num_elements != 4)
            PUT_ERR("wrong number of elements in deps\n")

         { // test whether the number of overlaps per cell is correct
            unsigned ref_num_deps_per_element[4] = {4,2,2,1};

            for (unsigned i = 0; i < 4; ++i)
               if (deps.num_deps_per_element[i] != ref_num_deps_per_element[i])
                  PUT_ERR("wrong number of matching cells in deps.num_deps_per_element[i]\n")
         }
         { // test whether the right cells were found
            unsigned ref_dependencies[4][4] = {{0,1,2,3}, {1,3,-1,-1},
                                             {2,3,-1,-1}, {3,-1,-1,-1}};
            unsigned curr_num_deps;
            unsigned const * curr_deps;

            for (unsigned i = 0; i < 4; ++i) {

               curr_num_deps = 0;
               curr_deps = yac_get_dependencies_of_element(deps, i);

               for (unsigned j = 0; j < deps.num_deps_per_element[i]; ++j) {
                  for (unsigned k = 0; k < deps.num_deps_per_element[i]; ++k) {

                     if (curr_deps[j] == ref_dependencies[i][k]) ++curr_num_deps;
                  }
               }

               if (curr_num_deps != deps.num_deps_per_element[i])
                  PUT_ERR("wrong dependencies in deps\n")
            }
         }

      //---------------
      // cleanup
      //---------------

         yac_free_dep_list(&deps);
         yac_delete_grid_search(search);
      }

      yac_delete_grid(grid_data_a);
      yac_delete_grid(grid_data_b);
   }

   {
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

      unsigned j, k;

      // set up grid data
      yac_init_dep_list(&cell_to_vertex_a);
      yac_set_dependencies(&cell_to_vertex_a, 4, num_corners_per_cell, cell_vertices);

      grid_data_a = yac_unstruct_grid_new(vertex_coordinates_x[0], vertex_coordinates_y[0], 
					  9, cell_to_vertex_a);

      yac_init_dep_list(&cell_to_vertex_b);
      yac_set_dependencies(&cell_to_vertex_b, 4, num_corners_per_cell, cell_vertices);

      grid_data_b = yac_unstruct_grid_new(vertex_coordinates_x[1], vertex_coordinates_y[1], 
					  9, cell_to_vertex_b);

      //---------------
      // testing
      //---------------

      for (int i = 0; i < num_grid_search_construct; ++i) {

         struct grid_search * search = grid_search_construct[i](grid_data_a);

         struct grid_cell cell;
         unsigned n_cells = 0, cells_size = 0, *cells = NULL;

         yac_init_grid_cell(&cell);

         for (unsigned i = 0; i < yac_get_num_grid_cells(grid_data_b); ++i) {

            yac_get_grid_cell(grid_data_b, i, &cell);

            yac_do_cell_search_single(search, cell, &n_cells, &cells_size, &cells);

            { // test whether the number of overlaps is correct
               unsigned ref_num_deps_per_element[4] = {4,2,2,1};

               if (n_cells != ref_num_deps_per_element[i])
                  PUT_ERR("wrong number of matching cells in n_cells\n")
            }
            { // test whether the right cells were found
               unsigned ref_dependencies[4][4] = {{0,1,2,3}, {1,3,-1,-1},
                                                {2,3,-1,-1}, {3,-1,-1,-1}};
               unsigned curr_num_deps = 0;

               for (j = 0; j < n_cells; ++j)
                  for (k = 0; k < n_cells; ++k)
                     if (cells[j] == ref_dependencies[i][k]) ++curr_num_deps;

               if (curr_num_deps != n_cells)
                  PUT_ERR("wrong dependencies in deps\n")
            }

         }

      //---------------
      // cleanup
      //---------------

         yac_delete_grid_search(search);
         free(cells);
         yac_free_grid_cell(&cell);
      }

      yac_delete_grid(grid_data_a);
      yac_delete_grid(grid_data_b);
   }

   { // test cell search for cyclic grids
      //---------------
      // setup
      //---------------
      struct grid * grid_a;
      double coords_x_a[12] = {-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150};
      double coords_y_a[7] = {-90, -60, -30, 0, 30, 60, 90};

      for ( unsigned i = 0; i < 12; ++i ) coords_x_a[i] *= rad;
      for ( unsigned i = 0; i <  7; ++i ) coords_y_a[i] *= rad;

      unsigned num_cells_a[2] = {12,6};
      unsigned cyclic_a[2] = {1,0};
   
      grid_a = yac_reg2d_grid_new(coords_x_a, coords_y_a, num_cells_a, cyclic_a);

      struct grid * grid_b;
      double coords_x_b[11] = {-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50};
      double coords_y_b[7] = {-30, -20, -10, 0, 10, 20, 30};

      for ( unsigned i = 0; i < 11; ++i ) coords_x_b[i] *= rad;
      for ( unsigned i = 0; i <  7; ++i ) coords_y_b[i] *= rad;

      unsigned num_cells_b[2] = {10,6};
      unsigned cyclic_b[2] = {0,0};
   
      grid_b = yac_reg2d_grid_new(coords_x_b, coords_y_b, num_cells_b, cyclic_b);

      //---------------
      // testing
      //---------------

      struct dep_list deps[2];

      for (int i = 0; i < num_grid_search_construct; ++i) {

         struct grid_search * search = grid_search_construct[i](grid_a);

         yac_do_cell_search(search, grid_b, deps + i);

         // test whether the number of found cells is correct
         if (deps[i].num_elements != 60)
            PUT_ERR("error in get_num_elements(deps)\n");

         yac_delete_grid_search(search);
      }

      check_dep_lists(deps);

      yac_free_dep_list(deps+0);
      yac_free_dep_list(deps+1);

      //---------------
      // cleanup
      //---------------
      yac_delete_grid(grid_a);
      yac_delete_grid(grid_b);
   }

   { // test cell search for cyclic grids ---> large grids

     printf("test cell search for cyclic grids\n");
      //---------------
      // setup
      //---------------
      struct grid * grid_a;
      const int nxa = 720;
      const int nya = 361;
      const int nxb = 360;
      const int nyb = 181;
      double coords_x_a[nxa];
      double coords_y_a[nya];

      for ( unsigned i = 0; i < nxa; ++i ) coords_x_a[i] = -180 + i*360./nxa;
      for ( unsigned i = 0; i < nya; ++i ) coords_y_a[i] =  -90 + i*360./nxa;
      // for ( unsigned i = 0; i < nxa; ++i ) printf("xa %d %g\n", i, coords_x_a[i]);
      // for ( unsigned i = 0; i < nya; ++i ) printf("ya %d %g\n", i, coords_y_a[i]);

      printf("a: %g %g ... %g %g\n", coords_x_a[0], coords_x_a[1], coords_x_a[nxa-2], coords_x_a[nxa-1]);

      for ( unsigned i = 0; i < nxa; ++i ) coords_x_a[i] *= rad;
      for ( unsigned i = 0; i < nya; ++i ) coords_y_a[i] *= rad;

      unsigned num_cells_a[2] = {nxa,nya-1};
      unsigned cyclic_a[2] = {1,0};
   
      grid_a = yac_reg2d_grid_new(coords_x_a, coords_y_a, num_cells_a, cyclic_a);

      struct grid * grid_b;
      double coords_x_b[nxb];
      double coords_y_b[nyb];

      for ( unsigned i = 0; i < nxb; ++i ) coords_x_b[i] = -180 + i*360./nxb;
      for ( unsigned i = 0; i < nyb; ++i ) coords_y_b[i] =  -90 + i*360./nxb;

      printf("b: %g %g ... %g %g\n", coords_x_b[0], coords_x_b[1], coords_x_b[nxb-2], coords_x_b[nxb-1]);

      for ( unsigned i = 0; i < nxb; ++i ) coords_x_b[i] *= rad;
      for ( unsigned i = 0; i < nyb; ++i ) coords_y_b[i] *= rad;

      unsigned num_cells_b[2] = {nxb,nyb-1};
      unsigned cyclic_b[2] = {1,0};
   
      grid_b = yac_reg2d_grid_new(coords_x_b, coords_y_b, num_cells_b, cyclic_b);

      //---------------
      // testing
      //---------------

      struct dep_list deps[2];

      clock_t start, finish;

      for (int i = 0; i < num_grid_search_construct; ++i) {
  	 start = clock();

         struct grid_search * search = grid_search_construct[i](grid_a);

         yac_do_cell_search(search, grid_b, deps + i);

         // test whether the number of found cells is correct
	 printf("deps[i].num_elements %d  ngp %d\n", deps[i].num_elements, nxb*(nyb-1));
         if (deps[i].num_elements != nxb*(nyb-1))
            PUT_ERR("error in get_num_elements(deps)\n");

         yac_delete_grid_search(search);

	 finish = clock();
	 printf("search %d: %.2f seconds\n", i+1, ((double)(finish-start))/CLOCKS_PER_SEC);
      }


      check_dep_lists(deps);

      for (int i = 0; i < num_grid_search_construct; ++i)
	yac_free_dep_list(deps+i);

      //---------------
      // cleanup
      //---------------
      yac_delete_grid(grid_a);
      yac_delete_grid(grid_b);
   }

   return TEST_EXIT_CODE;
}

static int compare_uint(const void * a, const void * b) {

  return ((*(unsigned*)a) > (*(unsigned*)b)) -
         ((*(unsigned*)a) < (*(unsigned*)b));
}

void check_dep_lists(struct dep_list * lists) {

   unsigned *results[2] = {NULL, NULL};
   unsigned results_size[2] = {0,0};

   if (lists[0].num_elements != lists[1].num_elements)
      PUT_ERR("num_elements does not match\n");

   for (unsigned i = 0; i < lists[0].num_elements; ++i) {

      if (lists[0].num_deps_per_element[i] != lists[1].num_deps_per_element[i])
         PUT_ERR("num_deps_per_element does not match\n");

      for (unsigned j = 0; j < 2; ++j) {

         ENSURE_ARRAY_SIZE(results[j], results_size[j],
                           lists[j].num_deps_per_element[i]);
         memcpy(results[j], yac_get_dependencies_of_element(lists[j],i),
                lists[j].num_deps_per_element[i] * sizeof(**results));
         qsort(results[j], lists[j].num_deps_per_element[i], sizeof(**results),
               compare_uint);
      }

      for (unsigned j = 0; j < lists[0].num_deps_per_element[i]; ++j)
	{
	  //	  printf("%d %d %d\n", lists[0].num_deps_per_element[i], j , results[0][j]);
	  if (results[0][j] != results[1][j])
            PUT_ERR("dependencies do not match\n");
	}
   }

   free(results[0]), free(results[1]);
}

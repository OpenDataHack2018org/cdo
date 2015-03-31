/**
 * @file sphere_part.c
 *
 * @copyright Copyright  (C)  2014 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
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

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "sphere_part.h"
#include "geometry.h"
#include "grid.h"
#include "interval_tree.h"
#include "utils.h"
#include "ensure_array_size.h"
#include "grid_search.h"
#include "grid_search_utils.h"


enum {
   I_list_tree_min_size = 2, //!< make I list into tree when list is
                             //!<  larger than this
};

struct sphere_part_search {

   struct grid_search_vtable * vtable;
   struct sphere_part_node base_node;
   struct grid * grid_data;
};

static void sphere_part_do_cell_search(struct grid_search * search,
                                       struct grid * grid_data,
                                       struct dep_list * tgt_to_src_cells);
static void sphere_part_do_cell_search_single(struct grid_search * search,
                                              struct grid_cell cell,
                                              unsigned * n_cells,
                                              unsigned * cells_size,
                                              unsigned ** cells);
static void sphere_part_do_point_search_c(struct grid_search * search,
                                          struct grid * grid_data,
                                          struct dep_list * tgt_to_src_cells);
static void sphere_part_do_point_search_c2(struct grid_search * search,
                                           double * x_coordinates,
                                           double * y_coordinates,
                                           unsigned num_points,
                                           struct dep_list * tgt_to_src_cells);
static void sphere_part_do_point_search_p(struct grid_search * search,
                                          struct grid * grid_data,
                                          struct dep_list * target_to_src_points);
static void sphere_part_do_point_search_p2(struct grid_search * search,
                                           double * x_coordinates,
                                           double * y_coordinates,
                                           unsigned num_points,
                                           struct dep_list * target_to_src_points);
static void sphere_part_do_point_search_p3(struct grid_search * search,
                                           double * x_coordinates,
                                           double * y_coordinates,
                                           unsigned num_points,
                                           struct dep_list * target_to_src_points,
                                           struct points * points);
static void sphere_part_do_point_search_p4 (struct grid_search * search,
                                            double x_coordinate,
                                            double y_coordinate,
                                            unsigned * n_points,
                                            unsigned * points_size,
                                            unsigned ** points);
static void delete_sphere_part_search(struct grid_search * search);

static struct grid_search_vtable sphere_part_search_vtable =
{
   .do_cell_search        = sphere_part_do_cell_search,
   .do_cell_search_single = sphere_part_do_cell_search_single,
   .do_point_search_c     = sphere_part_do_point_search_c,
   .do_point_search_c2    = sphere_part_do_point_search_c2,
   .do_point_search_p     = sphere_part_do_point_search_p,
   .do_point_search_p2    = sphere_part_do_point_search_p2,
   .do_point_search_p3    = sphere_part_do_point_search_p3,
   .do_point_search_p4    = sphere_part_do_point_search_p4,
   .delete_grid_search    = delete_sphere_part_search
};

static void normalise_vector(double v[]) {

   double norm;

   norm = 1.0 / sqrt(v[0]*v[0] +
                     v[1]*v[1] +
                     v[2]*v[2]);

   v[0] *= norm;
   v[1] *= norm;
   v[2] *= norm;
}

static void init_sphere_part_node(struct sphere_part_node * node) {

   node->flags = 0;
   node->I_size = 0;
   node->U_size = 0;
   node->T_size = 0;
   node->gc_norm_vector[0] = 0;
   node->gc_norm_vector[1] = 0;
   node->gc_norm_vector[2] = 1;
}

static struct sphere_part_node * get_sphere_part_node() {

   struct sphere_part_node * node;

   node = malloc(1 * sizeof(*node));

   init_sphere_part_node(node);

   return node;
}

static void partition_data (struct grid * grid, unsigned * local_cell_ids,
                            unsigned num_cell_ids, unsigned threshold,
                            struct sphere_part_node * parent_node,
                            double prev_gc_norm_vector[]) {

   unsigned i;

   double balance_point[3] = {0.0,0.0,0.0};

   // compute balance point

   struct grid_cell cell;
   struct bounding_circle bnd_circle;

   yac_init_grid_cell(&cell);

   for (i = 0; i < num_cell_ids; ++i) {

      yac_get_grid_cell2(grid, local_cell_ids[i], &cell, &bnd_circle);

      balance_point[0] += bnd_circle.base_vector[0];
      balance_point[1] += bnd_circle.base_vector[1];
      balance_point[2] += bnd_circle.base_vector[2];
   }

   normalise_vector(balance_point);

   // compute the great circle that partitions the data in half (more or less)

   crossproduct_ld(balance_point, prev_gc_norm_vector, parent_node->gc_norm_vector);
   normalise_vector(parent_node->gc_norm_vector);

   // partition data into cells that overlap with the great circle and cells that are
   // on side of the circle

   unsigned * I = NULL; // list of cells on the great circle
   unsigned * U = NULL; // list of cells on one side of the great circle
   unsigned * T = NULL; // list of cells on that other side of the great circle

   unsigned I_size = 0, I_array_size = 0;
   unsigned U_size = 0, U_array_size = 0;
   unsigned T_size = 0, T_array_size = 0;

   double max_inc_angle = 0.0;

   for (i = 0; i < num_cell_ids; ++i) {

      yac_get_grid_cell2(grid, local_cell_ids[i], &cell, &bnd_circle);

      // get angle between the norm vector of the great circle and the base
      // point of the bounding circle
      double angle;

      angle = get_vector_angle(bnd_circle.base_vector, parent_node->gc_norm_vector);

      // if the point intersects with the great circle
      if ((angle >= M_PI_2 && angle - bnd_circle.inc_angle <= M_PI_2) ||
          (angle <= M_PI_2 && angle + bnd_circle.inc_angle >= M_PI_2)) {

         angle = fabs(angle - M_PI_2);

         // add to list I
         ENSURE_ARRAY_SIZE(I, I_array_size, I_size + 1);
         I[I_size++] = local_cell_ids[i];
         max_inc_angle = MAX(max_inc_angle, angle + bnd_circle.inc_angle);

      } else if (angle > M_PI_2) {

         // add to list U
         ENSURE_ARRAY_SIZE(U, U_array_size, U_size + 1);
         U[U_size++] = local_cell_ids[i];

      } else {

         // add to list T
         ENSURE_ARRAY_SIZE(T, T_array_size, T_size + 1);
         T[T_size++] = local_cell_ids[i];
      }
   }

   parent_node->I_size = I_size;
   parent_node->U_size = U_size;
   parent_node->T_size = T_size;
   parent_node->I_angle = max_inc_angle;

   if (I_size > 0) {

      if (I_size > I_list_tree_min_size) {

         assert(sizeof(struct interval_node) > sizeof(unsigned));
         parent_node->I.ivt.head_node
            = (void *)(I = realloc(I, I_size * sizeof(struct interval_node)));
         parent_node->I.ivt.num_nodes = I_size;

         for (i = I_size - 1 ; i < I_size; --i) {

            double GCp[3], bVp[3], base_angle;
            struct interval iv;
            unsigned cell_idx = I[i];

            yac_get_grid_cell2(grid, I[i], &cell, &bnd_circle);
            crossproduct_ld(parent_node->gc_norm_vector,
                         bnd_circle.base_vector, GCp);
            crossproduct_ld(GCp, parent_node->gc_norm_vector, bVp);
            normalise_vector(bVp);
            base_angle = get_vector_angle(bVp, prev_gc_norm_vector);

            iv.left = base_angle - bnd_circle.inc_angle;
            iv.right = base_angle + bnd_circle.inc_angle;

            parent_node->I.ivt.head_node[i].range = iv;
            parent_node->I.ivt.head_node[i].value = cell_idx;
         }

         parent_node->I.ivt.num_nodes = I_size;
         parent_node->I.ivt.head_node = realloc(parent_node->I.ivt.head_node,
            I_size * sizeof(*(parent_node->I.ivt.head_node)));

         yac_generate_interval_tree(parent_node->I.ivt.head_node,
                                    parent_node->I.ivt.num_nodes);
         parent_node->flags |= I_IS_INTERVAL_TREE;
      } else
         parent_node->I.list = realloc(I, I_size * sizeof(unsigned));
   } else
      parent_node->I.list = NULL;

   yac_free_grid_cell(&cell);

   // check whether the lists are small enough (if not -> partition again)
   if (U_size <= threshold) {

      parent_node->U = realloc(U, U_size * sizeof(unsigned));
      parent_node->flags |= U_IS_LEAF;

   } else {

      parent_node->U = get_sphere_part_node();
      partition_data(grid, U, U_size, threshold, parent_node->U,
                     parent_node->gc_norm_vector);
      free(U);
   }

   if (T_size <= threshold) {

      parent_node->T = realloc(T, T_size * sizeof(unsigned));
      parent_node->flags |= T_IS_LEAF;

   } else {

      parent_node->T = get_sphere_part_node();
      partition_data(grid, T, T_size, threshold, parent_node->T,
                     parent_node->gc_norm_vector);
      free(T);
   }
}

struct grid_search * yac_sphere_part_search_new (struct grid * grid) {

   struct sphere_part_search * search = malloc(1 * sizeof(*search));

   search->vtable = &sphere_part_search_vtable;
   search->grid_data = grid;

   double gc_norm_vector[3] = {0.0,0.0,1.0};
   
   init_sphere_part_node(&(search->base_node));

   unsigned * local_cell_ids;
   unsigned num_grid_cells;

   num_grid_cells = yac_get_num_grid_cells(grid);
   local_cell_ids = malloc(num_grid_cells * sizeof(*local_cell_ids));

   unsigned i;

   for (i = 0; i < num_grid_cells; ++i)
      local_cell_ids[i] = i;

   partition_data(grid, local_cell_ids, num_grid_cells, I_list_tree_min_size,
                  &(search->base_node), gc_norm_vector);

   free(local_cell_ids);

   return (struct grid_search *)search;
}

static void search_bnd_circle(struct sphere_part_node * node,
                              struct bounding_circle bnd_circle,
                              unsigned ** overlap_cells,
                              unsigned * overlap_cells_array_size,
                              unsigned * num_overlap_cells,
                              struct overlaps * search_interval_tree_buffer,
                              double prev_gc_norm_vector[]) {

   double angle = get_vector_angle(node->gc_norm_vector,
                                   bnd_circle.base_vector);

   if (angle - bnd_circle.inc_angle < M_PI_2) {

      if (node->flags & T_IS_LEAF) {

         ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                           *num_overlap_cells + node->T_size);
         memcpy((*overlap_cells) + (*num_overlap_cells),
                node->T, node->T_size * sizeof(**overlap_cells));
         *num_overlap_cells += node->T_size;

      } else {
         search_bnd_circle(node->T, bnd_circle, overlap_cells,
                           overlap_cells_array_size, num_overlap_cells,
                           search_interval_tree_buffer, node->gc_norm_vector);
      }
   }

   if (angle + bnd_circle.inc_angle > M_PI_2) {

      if (node->flags & U_IS_LEAF) {

         ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                           *num_overlap_cells + node->U_size);
         memcpy((*overlap_cells) + (*num_overlap_cells),
                node->U, node->U_size * sizeof(**overlap_cells));
         *num_overlap_cells += node->U_size;

      } else {
         search_bnd_circle(node->U, bnd_circle, overlap_cells,
                           overlap_cells_array_size, num_overlap_cells,
                           search_interval_tree_buffer, node->gc_norm_vector);
      }
   }

   if (fabs(angle - M_PI_2) <= (bnd_circle.inc_angle + node->I_angle)) {

      if (node->flags & I_IS_INTERVAL_TREE) {

         double GCp[3], bVp[3], base_angle;
         crossproduct_ld(node->gc_norm_vector,
                      bnd_circle.base_vector, GCp);
         crossproduct_ld(GCp, node->gc_norm_vector, bVp);
         normalise_vector(bVp);
         base_angle = get_vector_angle(bVp, prev_gc_norm_vector);

         struct interval search_interval =
            {.left = base_angle - bnd_circle.inc_angle,
             .right = base_angle + bnd_circle.inc_angle};

         search_interval_tree_buffer->num_overlaps = 0;

         yac_search_interval_tree(node->I.ivt.head_node, node->I.ivt.num_nodes,
                                  search_interval, search_interval_tree_buffer);

         ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                           *num_overlap_cells +
                           search_interval_tree_buffer->num_overlaps);

         for (size_t i = 0; i < search_interval_tree_buffer->num_overlaps;
              ++i) {
            (*overlap_cells)[(*num_overlap_cells)+i] =
               node->I.ivt.head_node[search_interval_tree_buffer->overlap_iv[i]].value;
         }

         *num_overlap_cells += search_interval_tree_buffer->num_overlaps;

      } else {

         ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                           *num_overlap_cells + node->I_size);
         memcpy((*overlap_cells) + (*num_overlap_cells),
                node->I.list, node->I_size * sizeof(**overlap_cells));
         *num_overlap_cells += node->I_size;
      }
   }
}

static void search_point(struct sphere_part_node * node,
                         double point[],
                         unsigned ** overlap_cells,
                         unsigned * overlap_cells_array_size,
                         unsigned * num_overlap_cells,
                         struct overlaps * search_interval_tree_buffer,
                         double prev_gc_norm_vector[]) {

   double angle = get_vector_angle(node->gc_norm_vector, point);

   if (angle < M_PI_2) {

      if (node->flags & T_IS_LEAF) {

         ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                           *num_overlap_cells + node->T_size);
         memcpy((*overlap_cells) + (*num_overlap_cells),
                node->T, node->T_size * sizeof(**overlap_cells));
         *num_overlap_cells += node->T_size;

      } else {
         search_point(node->T, point, overlap_cells,
                      overlap_cells_array_size, num_overlap_cells,
                      search_interval_tree_buffer, node->gc_norm_vector);
      }
   }

   if (angle > M_PI_2) {

      if (node->flags & U_IS_LEAF) {

         ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                           *num_overlap_cells + node->U_size);
         memcpy((*overlap_cells) + (*num_overlap_cells),
                node->U, node->U_size * sizeof(**overlap_cells));
         *num_overlap_cells += node->U_size;

      } else {
         search_point(node->U, point, overlap_cells,
                      overlap_cells_array_size, num_overlap_cells,
                      search_interval_tree_buffer, node->gc_norm_vector);
      }
   }

   if (fabs(angle - M_PI_2) <= node->I_angle) {

      if (node->flags & I_IS_INTERVAL_TREE) {

         double GCp[3], bVp[3], base_angle;
         crossproduct_ld(node->gc_norm_vector, point, GCp);
         crossproduct_ld(GCp, node->gc_norm_vector, bVp);
         normalise_vector(bVp);
         base_angle = get_vector_angle(bVp, prev_gc_norm_vector);

         struct interval search_interval =
            {.left = base_angle, .right = base_angle};

         search_interval_tree_buffer->num_overlaps = 0;

         yac_search_interval_tree(node->I.ivt.head_node, node->I.ivt.num_nodes,
                                  search_interval, search_interval_tree_buffer);

         ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                           *num_overlap_cells +
                           search_interval_tree_buffer->num_overlaps);

         for (size_t i = 0; i < search_interval_tree_buffer->num_overlaps;
              ++i) {
            (*overlap_cells)[(*num_overlap_cells)+i] =
               node->I.ivt.head_node[search_interval_tree_buffer->overlap_iv[i]].value;
         }

         *num_overlap_cells += search_interval_tree_buffer->num_overlaps;

      } else {

         ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                           *num_overlap_cells + node->I_size);
         memcpy((*overlap_cells) + (*num_overlap_cells),
                node->I.list, node->I_size * sizeof(**overlap_cells));
         *num_overlap_cells += node->I_size;
      }
   }
}

static void sphere_part_do_cell_search(struct grid_search * search,
                                       struct grid * grid_data,
                                       struct dep_list * tgt_to_src_cells) {

   struct sphere_part_search * sp_search = (struct sphere_part_search *)search;
   struct sphere_part_node * base_node = &(sp_search->base_node);

   unsigned num_cells = yac_get_num_grid_cells(grid_data);

   unsigned * temp_search_results = NULL;
   unsigned temp_search_results_array_size = 0;
   unsigned num_temp_search_results = 0;

   unsigned i, j;
   struct grid_cell cell_a, cell_b;
   struct bounding_circle circle_a, circle_b;

   yac_init_grid_cell(&cell_a);
   yac_init_grid_cell(&cell_b);

   struct overlaps search_interval_tree_buffer = {0, 0, NULL};

   unsigned * num_src_per_tgt_cell = calloc(num_cells,
                                            sizeof(*num_src_per_tgt_cell));
   unsigned * tgt_src_dependencies = NULL;
   unsigned tgt_src_dependencies_size = 0;
   unsigned total_num_dependencies = 0;

   for (i = 0; i < num_cells; ++i) {

      yac_get_grid_cell2(grid_data, i, &cell_a, &circle_a);

      num_temp_search_results = 0;

      double gc_norm_vector[3] = {0.0,0.0,1.0};

      search_bnd_circle(base_node, circle_a, &temp_search_results,
                        &temp_search_results_array_size,
                        &num_temp_search_results, &search_interval_tree_buffer,
                        gc_norm_vector);

      ENSURE_ARRAY_SIZE(tgt_src_dependencies, tgt_src_dependencies_size,
                        total_num_dependencies + num_temp_search_results);

      for (j = 0; j < num_temp_search_results; ++j) {

         yac_get_grid_cell2(sp_search->grid_data, temp_search_results[j],
                            &cell_b, &circle_b);

         if (yac_check_overlap_cells2(cell_a, circle_a, cell_b, circle_b)) {

            tgt_src_dependencies[total_num_dependencies++] =
               temp_search_results[j];
            num_src_per_tgt_cell[i]++;
         }
      }
   }

   yac_free_grid_cell(&cell_a);
   yac_free_grid_cell(&cell_b);
   free(temp_search_results);
   free(search_interval_tree_buffer.overlap_iv);

   tgt_src_dependencies = realloc(tgt_src_dependencies, total_num_dependencies *
                                  sizeof(tgt_src_dependencies));
   yac_init_dep_list(tgt_to_src_cells);
   yac_set_dependencies(tgt_to_src_cells, num_cells, num_src_per_tgt_cell,
                        tgt_src_dependencies);
}

static void sphere_part_do_cell_search_single(struct grid_search * search,
                                              struct grid_cell cell,
                                              unsigned * n_cells,
                                              unsigned * cells_size,
                                              unsigned ** cells) {

   struct sphere_part_search * sp_search = (struct sphere_part_search *)search;
   struct sphere_part_node * base_node = &(sp_search->base_node);

   struct bounding_circle circle_a;

   yac_get_cell_bounding_circle(cell, &circle_a);

   struct overlaps search_interval_tree_buffer = {0, 0, NULL};

   *n_cells = 0;

   double gc_norm_vector[3] = {0.0,0.0,1.0};

   search_bnd_circle(base_node, circle_a, cells, cells_size, n_cells,
                     &search_interval_tree_buffer, gc_norm_vector);

   unsigned n_cells_final = 0;

   struct grid_cell cell_b;
   struct bounding_circle circle_b;

   yac_init_grid_cell(&cell_b);

   for (unsigned j = 0; j < *n_cells; ++j) {

      yac_get_grid_cell2(sp_search->grid_data, (*cells)[j],
                         &cell_b, &circle_b);

      if (yac_check_overlap_cells2(cell, circle_a, cell_b, circle_b))
         (*cells)[n_cells_final++] = (*cells)[j];
   }

   *n_cells = n_cells_final;

   yac_free_grid_cell(&cell_b);
   free(search_interval_tree_buffer.overlap_iv);
}

static void sphere_part_do_point_search_c(struct grid_search * search,
                                          struct grid * grid_data,
                                          struct dep_list * tgt_to_src_cells) {

   struct sphere_part_search * sp_search = (struct sphere_part_search *)search;
   struct sphere_part_node * base_node = &(sp_search->base_node);

   unsigned num_corners = yac_get_num_grid_corners(grid_data);

   unsigned * temp_search_results = NULL;
   unsigned temp_search_results_array_size = 0;
   unsigned num_temp_search_results = 0;

   unsigned i, j;
   struct grid_cell cell;
   struct bounding_circle bnd_circle;

   yac_init_grid_cell(&cell);

   struct overlaps search_interval_tree_buffer = {0, 0, NULL};

   unsigned * num_src_per_tgt_corner = calloc(num_corners,
                                              sizeof(*num_src_per_tgt_corner));
   unsigned * tgt_src_dependencies = NULL;
   unsigned tgt_src_dependencies_size = 0;
   unsigned total_num_dependencies = 0;

   for (i = 0; i < num_corners; ++i) {

      struct point point = {.lon = yac_get_corner_x_coord(grid_data, i),
                            .lat = yac_get_corner_y_coord(grid_data, i)};
      double point_3d[3];

      LLtoXYZ(point.lon, point.lat, point_3d);

      num_temp_search_results = 0;

      double gc_norm_vector[3] = {0.0,0.0,1.0};

      search_point(base_node, point_3d, &temp_search_results,
                   &temp_search_results_array_size, &num_temp_search_results,
                   &search_interval_tree_buffer, gc_norm_vector);

      ENSURE_ARRAY_SIZE(tgt_src_dependencies, tgt_src_dependencies_size,
                        total_num_dependencies + num_temp_search_results);

      for (j = 0; j < num_temp_search_results; ++j) {

         yac_get_grid_cell2(sp_search->grid_data, temp_search_results[j],
                            &cell, &bnd_circle);

         if (yac_point_in_cell2(point, point_3d, cell, bnd_circle)) {

            tgt_src_dependencies[total_num_dependencies++] =
               temp_search_results[j];
            num_src_per_tgt_corner[i]++;
         }
      }
   }

   yac_free_grid_cell(&cell);
   free(temp_search_results);
   free(search_interval_tree_buffer.overlap_iv);

   tgt_src_dependencies = realloc(tgt_src_dependencies, total_num_dependencies *
                                  sizeof(tgt_src_dependencies));
   yac_init_dep_list(tgt_to_src_cells);
   yac_set_dependencies(tgt_to_src_cells, num_corners, num_src_per_tgt_corner,
                        tgt_src_dependencies);
}

static void sphere_part_do_point_search_c2(struct grid_search * search,
                                           double * x_coordinates,
                                           double * y_coordinates,
                                           unsigned num_points,
                                           struct dep_list * tgt_to_src_cells) {

   struct sphere_part_search * sp_search = (struct sphere_part_search *)search;
   struct sphere_part_node * base_node = &(sp_search->base_node);

   unsigned * temp_search_results = NULL;
   unsigned temp_search_results_array_size = 0;
   unsigned num_temp_search_results = 0;

   unsigned i, j;
   struct grid_cell cell;
   struct bounding_circle bnd_circle;

   yac_init_grid_cell(&cell);

   struct overlaps search_interval_tree_buffer = {0, 0, NULL};

   unsigned * num_src_per_tgt_corner = calloc(num_points,
                                              sizeof(*num_src_per_tgt_corner));
   unsigned * tgt_src_dependencies = NULL;
   unsigned tgt_src_dependencies_size = 0;
   unsigned total_num_dependencies = 0;

   for (i = 0; i < num_points; ++i) {

      struct point point = {.lon = x_coordinates[i],
                            .lat = y_coordinates[i]};
      double point_3d[3];

      LLtoXYZ(point.lon, point.lat, point_3d);

      num_temp_search_results = 0;

      double gc_norm_vector[3] = {0.0,0.0,1.0};

      search_point(base_node, point_3d, &temp_search_results,
                   &temp_search_results_array_size, &num_temp_search_results,
                   &search_interval_tree_buffer, gc_norm_vector);

      ENSURE_ARRAY_SIZE(tgt_src_dependencies, tgt_src_dependencies_size,
                        total_num_dependencies + num_temp_search_results);

      for (j = 0; j < num_temp_search_results; ++j) {

         yac_get_grid_cell2(sp_search->grid_data, temp_search_results[j],
                            &cell, &bnd_circle);

         if (yac_point_in_cell2(point, point_3d, cell, bnd_circle)) {

            tgt_src_dependencies[total_num_dependencies++] =
               temp_search_results[j];
            num_src_per_tgt_corner[i]++;
         }
      }
   }

   yac_free_grid_cell(&cell);
   free(temp_search_results);
   free(search_interval_tree_buffer.overlap_iv);

   tgt_src_dependencies = realloc(tgt_src_dependencies, total_num_dependencies *
                                  sizeof(tgt_src_dependencies));
   yac_init_dep_list(tgt_to_src_cells);
   yac_set_dependencies(tgt_to_src_cells, num_points, num_src_per_tgt_corner,
                        tgt_src_dependencies);
}

static void sphere_part_do_point_search_p(struct grid_search * search,
                                          struct grid * grid_data,
                                          struct dep_list * target_to_src_points) {

   struct sphere_part_search * sp_search = (struct sphere_part_search *)search;

   yac_grid_search_utils_do_point_search_p(search, sp_search->grid_data, grid_data,
                                           target_to_src_points);
}

static void sphere_part_do_point_search_p2(struct grid_search * search,
                                           double * x_coordinates,
                                           double * y_coordinates,
                                           unsigned num_points,
                                           struct dep_list * target_to_src_points) {

   struct sphere_part_search * sp_search = (struct sphere_part_search *)search;

   yac_grid_search_utils_do_point_search_p2(search, sp_search->grid_data,
                                            x_coordinates, y_coordinates,
                                            num_points, target_to_src_points);
}

static void sphere_part_do_point_search_p3(struct grid_search * search,
                                           double * x_coordinates,
                                           double * y_coordinates,
                                           unsigned num_points,
                                           struct dep_list * target_to_src_points,
                                           struct points * points) {

   yac_grid_search_utils_do_point_search_p3(search, x_coordinates, y_coordinates,
                                            num_points, target_to_src_points,
                                            points);
}


static void sphere_part_do_point_search_p4 (struct grid_search * search,
                                            double x_coordinate,
                                            double y_coordinate,
                                            unsigned * n_points,
                                            unsigned * points_size,
                                            unsigned ** points) {

   struct sphere_part_search * sp_search = (struct sphere_part_search *)search;

   yac_grid_search_utils_do_point_search_p4(search, sp_search->grid_data,
                                            x_coordinate, y_coordinate, n_points,
                                            points_size, points);
}

static void free_sphere_part_tree (struct sphere_part_node * tree) {


   // free I_list
   if (tree->flags & I_IS_INTERVAL_TREE)
      free(tree->I.ivt.head_node);
   else
      free(tree->I.list);

   if ((tree->flags & U_IS_LEAF) == 0)
      free_sphere_part_tree(tree->U);
   free(tree->U);

   if ((tree->flags & T_IS_LEAF) == 0)
      free_sphere_part_tree(tree->T);
   free(tree->T);
}

static void delete_sphere_part_search(struct grid_search * search) {

   struct sphere_part_search * sp_search = (struct sphere_part_search *)search;

   free_sphere_part_tree(&(sp_search->base_node));

   free(sp_search);
}

#ifdef SPHERE_PART_DEBUG
struct sphere_part_node * yac_get_sphere_part_tree(struct grid_search * search) {

   return &(((struct sphere_part_search *)search)->base_node);
}
#endif

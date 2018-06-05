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

union I_list {
  struct {
    struct interval_node * head_node;
    size_t num_nodes;
  } ivt;
  unsigned *list;
};

enum {
   I_list_tree_min_size = 2, //!< make I list into tree when list is
                             //!<  larger than this
};

enum yac_node_flags {
   U_IS_LEAF = 1,
   T_IS_LEAF = 2,
   I_IS_INTERVAL_TREE = 4,
};

struct sphere_part_node {

   unsigned flags;

   union I_list I;
   void * U, * T;

   size_t I_size, U_size, T_size;

   struct sin_cos_angle I_angle;

   double gc_norm_vector[3];
};

struct sphere_part_search {

   struct grid_search_vtable * vtable;
   struct sphere_part_node base_node;
   unsigned * local_cell_ids;
   struct grid * grid_data;
};

enum node_type {
  I_NODE = 0,
  U_NODE = 1,
  T_NODE = 2,
};

struct temp_partition_data {
  unsigned local_id;
  struct bounding_circle bnd_circle;
  int node_type;
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
                                           double (*coordinates_xyz)[3],
                                           unsigned num_points,
                                           struct dep_list * tgt_to_src_cells);
static void sphere_part_do_point_search_c3(struct grid_search * search,
                                           double (*coordinates_xyz)[3],
                                           unsigned num_points,
                                           struct dep_list * tgt_to_src_cells,
                                           struct points * points);
static void sphere_part_do_point_search_p(struct grid_search * search,
                                          struct grid * grid_data,
                                          struct dep_list * target_to_src_points);
static void sphere_part_do_point_search_p2(struct grid_search * search,
                                           double (*coordinates_xyz)[3],
                                           unsigned num_points,
                                           struct dep_list * target_to_src_points);
static void sphere_part_do_point_search_p3(struct grid_search * search,
                                           double (*coordinates_xyz)[3],
                                           unsigned num_points,
                                           struct dep_list * target_to_src_points,
                                           struct points * points);
static void sphere_part_do_point_search_p4 (struct grid_search * search,
                                            double coordinate_xyz[3],
                                            unsigned * n_points,
                                            unsigned * points_size,
                                            unsigned ** points);
static void sphere_part_do_bnd_circle_search (struct grid_search * search,
                                              struct bounding_circle * bnd_circles,
                                              unsigned n,
                                              struct dep_list * bnd_to_cells);
static void delete_sphere_part_search(struct grid_search * search);

static struct grid_search_vtable sphere_part_search_vtable =
{
   .do_cell_search        = sphere_part_do_cell_search,
   .do_cell_search_single = sphere_part_do_cell_search_single,
   .do_point_search_c     = sphere_part_do_point_search_c,
   .do_point_search_c2    = sphere_part_do_point_search_c2,
   .do_point_search_c3    = sphere_part_do_point_search_c3,
   .do_point_search_p     = sphere_part_do_point_search_p,
   .do_point_search_p2    = sphere_part_do_point_search_p2,
   .do_point_search_p3    = sphere_part_do_point_search_p3,
   .do_point_search_p4    = sphere_part_do_point_search_p4,
   .do_bnd_circle_search  = sphere_part_do_bnd_circle_search,
   .delete_grid_search    = delete_sphere_part_search
};

struct point_id_xyz {
  size_t idx; // index of the point in the coordinates array passed to
              // the constructor
  double coordinates_xyz[3];
};

struct point_sphere_part_node {

   unsigned flags;

   void * U, * T;

   size_t U_size, T_size;

   double gc_norm_vector[3];
};

struct point_sphere_part_search {

   struct point_sphere_part_node base_node;
   struct point_id_xyz * points;
};

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

   struct sphere_part_node * node = malloc(1 * sizeof(*node));

   init_sphere_part_node(node);

   return node;
}

static int compare_temp_partition_data(const void * a, const void * b) {

  return ((const struct temp_partition_data *)a)->node_type -
         ((const struct temp_partition_data *)b)->node_type;
}

static void partition_data (
   struct grid * grid, unsigned * local_cell_ids,
   struct temp_partition_data * part_data, size_t num_cell_ids,
   size_t threshold, struct sphere_part_node * parent_node,
   double prev_gc_norm_vector[]) {

   double balance_point[3] = {0.0,0.0,0.0};

   // compute balance point
   for (size_t i = 0; i < num_cell_ids; ++i) {

     balance_point[0] += part_data[i].bnd_circle.base_vector[0];
     balance_point[1] += part_data[i].bnd_circle.base_vector[1];
     balance_point[2] += part_data[i].bnd_circle.base_vector[2];
   }

   normalise_vector(balance_point);

   // compute the great circle that partitions the data in half (more or less)

   crossproduct_ld(balance_point, prev_gc_norm_vector,
                   parent_node->gc_norm_vector);
   normalise_vector(parent_node->gc_norm_vector);

   // partition data into cells that overlap with the great circle and cells
   // that are on side of the circle

   size_t I_size = 0;
   size_t U_size = 0;
   size_t T_size = 0;

   struct sin_cos_angle max_inc_angle = SIN_COS_ZERO;

   for (size_t i = 0; i < num_cell_ids; ++i) {

      struct bounding_circle curr_bnd_circle = part_data[i].bnd_circle;

      // get angle between the norm vector of the great circle and the base
      // point of the bounding circle
      struct sin_cos_angle angle =
        get_vector_angle_2(
          curr_bnd_circle.base_vector, parent_node->gc_norm_vector);

      // get the angle between between the plane of the great circle and base
      // point of the bounding circle
      struct sin_cos_angle diff_angle_gc =
        sin_cos_angle_new(fabs(angle.cos), angle.sin);

      // if the point intersects with the great circle
      if (compare_angles(diff_angle_gc, curr_bnd_circle.inc_angle) <= 0) {

         // set node type for current cell
         part_data[i].node_type = I_NODE;
         I_size++;

         struct sin_cos_angle inc_angle =
            sum_angles_no_check(diff_angle_gc, curr_bnd_circle.inc_angle);

         if (compare_angles(max_inc_angle, inc_angle) < 0)
            max_inc_angle = inc_angle;

      // angle > M_PI_2
      } else if (angle.cos < 0.0) {

         // set node type for current cell
         part_data[i].node_type = U_NODE;
         U_size++;

      } else {

         // set node type for current cell
         part_data[i].node_type = T_NODE;
         T_size++;
      }
   }

   qsort(
      part_data, num_cell_ids, sizeof(*part_data), compare_temp_partition_data);

   parent_node->I_size = I_size;
   parent_node->U_size = U_size;
   parent_node->T_size = T_size;

   // if max_inc_angle > PI/2
   if (compare_angles(max_inc_angle, SIN_COS_M_PI_2) >= 0) {
      parent_node->I_angle = SIN_COS_M_PI_2;
   } else {
      parent_node->I_angle = max_inc_angle;
   }

   if (I_size > 0) {

      if (I_size > I_list_tree_min_size) {

         assert(sizeof(struct interval_node) > sizeof(unsigned));
         struct interval_node * head_node = malloc(I_size * sizeof(*head_node));
         parent_node->I.ivt.head_node = head_node;
         parent_node->I.ivt.num_nodes = I_size;

         for (size_t i = 0; i < I_size; ++i) {

            struct sin_cos_angle base_angle, corrected_inc_angle;
            int big_sum, neg;
            double GCp[3], bVp[3];
            struct bounding_circle curr_bnd_circle = part_data[i].bnd_circle;

            crossproduct_ld(parent_node->gc_norm_vector,
                            curr_bnd_circle.base_vector, GCp);
            crossproduct_ld(GCp, parent_node->gc_norm_vector, bVp);
            normalise_vector(bVp);
            base_angle = get_vector_angle_2(bVp, prev_gc_norm_vector);
            big_sum =
              sum_angles(curr_bnd_circle.inc_angle,
                         get_vector_angle_2(bVp, curr_bnd_circle.base_vector),
                         &corrected_inc_angle);

            // if the angle is bigger then PI
            if ((big_sum) || (corrected_inc_angle.sin < 0.0))
              corrected_inc_angle = SIN_COS_M_PI;

            struct sin_cos_angle left, right;
            // base_angle - corrected_inc_angle
            neg = sub_angles(base_angle, corrected_inc_angle, &left);
            // base_angle + corrected_inc_angle
            big_sum = sum_angles(base_angle, corrected_inc_angle, &right);

            head_node[i].range.left = compute_angle(left) - (neg?2.0*M_PI:0.0);
            head_node[i].range.right = compute_angle(right) +
                                       (big_sum?2.0*M_PI:0.0);
            head_node[i].value = part_data[i].local_id;
         }

         yac_generate_interval_tree(head_node, I_size);
         parent_node->flags |= I_IS_INTERVAL_TREE;
      } else {
         for (size_t i = 0; i < I_size; ++i)
            local_cell_ids[i] = part_data[i].local_id;
         parent_node->I.list = (void*)local_cell_ids;
      }
   } else
      parent_node->I.list = NULL;

   part_data += I_size;
   local_cell_ids += I_size;

   // check whether the lists are small enough (if not -> partition again)
   if (U_size <= threshold) {

      for (size_t i = 0; i < U_size; ++i)
         local_cell_ids[i] = part_data[i].local_id;
      parent_node->U = (void*)local_cell_ids;
      parent_node->flags |= U_IS_LEAF;

   } else {

      parent_node->U = get_sphere_part_node();
      partition_data(grid, local_cell_ids, part_data, U_size, threshold,
                     parent_node->U, parent_node->gc_norm_vector);
   }
   local_cell_ids += U_size;
   part_data += U_size;

   if (T_size <= threshold) {

      for (size_t i = 0; i < T_size; ++i)
         local_cell_ids[i] = part_data[i].local_id;
      parent_node->T = (void*)local_cell_ids;
      parent_node->flags |= T_IS_LEAF;
      local_cell_ids += T_size;

   } else {

      parent_node->T = get_sphere_part_node();
      partition_data(grid, local_cell_ids, part_data, T_size, threshold,
                     parent_node->T, parent_node->gc_norm_vector);
   }
}

static int compare_point_idx_xyz(void const * a, void const * b) {
  return (((struct point_id_xyz *)a)->idx > ((struct point_id_xyz *)b)->idx) -
         (((struct point_id_xyz *)a)->idx < ((struct point_id_xyz *)b)->idx);
}

static struct point_sphere_part_node * partition_point_data (
  struct point_id_xyz * points, size_t num_points, size_t threshold,
  double prev_gc_norm_vector[]) {

  struct point_sphere_part_node * node = malloc(1 * sizeof(*node));

  double balance_point[3] = {0.0,0.0,0.0};

  // compute balance point

  for (size_t i = 0; i < num_points; ++i) {

    double * curr_coordinates_xyz = &(points[i].coordinates_xyz[0]);
    balance_point[0] += curr_coordinates_xyz[0];
    balance_point[1] += curr_coordinates_xyz[1];
    balance_point[2] += curr_coordinates_xyz[2];
  }

  normalise_vector(balance_point);

  // compute the great circle that partitions the data in half (more or less)

  double * gc_norm_vector = &(node->gc_norm_vector[0]);
  crossproduct_ld(balance_point, prev_gc_norm_vector, gc_norm_vector);
  normalise_vector(gc_norm_vector);

  // angle between a point and the great circle plane
  // acos(dot(gc_norm_vector, point_xyz)) = angle(gc_norm_vector, point_xyz)
  // acos(dot(gc_norm_vector, point_xyz)) - PI/2 = angle(gc_plane, point_xyz)
  // dot <= 0.0    -> U list
  // dot >  0.0    -> T list

  struct point_id_xyz * left = points, * right = points + num_points - 1;

  // sort such that all points for the U list come first, followed by the
  // elements of the T list
  while (1) {
    // find element that does not belong into U-list
    while (left <= right) {
      double * curr_coordinates_xyz = &(left->coordinates_xyz[0]);
      double dot = curr_coordinates_xyz[0] * gc_norm_vector[0] +
                   curr_coordinates_xyz[1] * gc_norm_vector[1] +
                   curr_coordinates_xyz[2] * gc_norm_vector[2];

      // if (angle < M_PI_2)
      if (dot > 0.0) break;
      ++left;
    };

    // find element that does not belong into T-list
    while (left < right) {
      double * curr_coordinates_xyz = &(right->coordinates_xyz[0]);
      double dot = curr_coordinates_xyz[0] * gc_norm_vector[0] +
                   curr_coordinates_xyz[1] * gc_norm_vector[1] +
                   curr_coordinates_xyz[2] * gc_norm_vector[2];

      // if (angle >= M_PI_2)
      if (dot <= 0.0) break;
      --right;
    }

    if (left < right) {
      struct point_id_xyz tmp_point = *left;
      *left = *right;
      *right = tmp_point;
      ++left;
      --right;
    } else {
      break;
    }
  }

  size_t U_size = left - points;
  size_t T_size = num_points - U_size;

  node->U_size = U_size;
  node->T_size = T_size;
  node->flags = 0;

  // check whether the lists are small enough (if not -> partition again)
  if (U_size <= threshold) {

    node->U = points;
    node->flags |= U_IS_LEAF;
    qsort(points, U_size, sizeof(*points), compare_point_idx_xyz);

  } else {

    node->U = partition_point_data(points, U_size, threshold, gc_norm_vector);
  }

  if (T_size <= threshold) {

    node->T = points + U_size;
    node->flags |= T_IS_LEAF;
    qsort(points + U_size, T_size, sizeof(*points), compare_point_idx_xyz);

  } else {

    node->T =
      partition_point_data(points + U_size, T_size, threshold, gc_norm_vector);
  }

  return node;
}

struct grid_search * yac_sphere_part_search_new (struct grid * grid) {

   struct sphere_part_search * search = malloc(1 * sizeof(*search));

   search->vtable = &sphere_part_search_vtable;
   search->grid_data = grid;

   double gc_norm_vector[3] = {0.0,0.0,1.0};
   
   init_sphere_part_node(&(search->base_node));

   size_t num_grid_cells = (size_t)yac_get_num_grid_cells(grid);
   unsigned * local_cell_ids = malloc(num_grid_cells * sizeof(*local_cell_ids));
   search->local_cell_ids = local_cell_ids;

   struct grid_cell dummy_cell;
   struct temp_partition_data * part_data =
      malloc(num_grid_cells * sizeof(*part_data));
   yac_init_grid_cell(&dummy_cell);
   for (size_t i = 0; i < num_grid_cells; ++i) {
      yac_get_grid_cell2(grid, i, &dummy_cell, &(part_data[i].bnd_circle));
      part_data[i].local_id = i;
   }
   yac_free_grid_cell(&dummy_cell);

   partition_data(grid, local_cell_ids, part_data, num_grid_cells,
                  I_list_tree_min_size, &(search->base_node), gc_norm_vector);

   free(part_data);

   return (struct grid_search *)search;
}

struct point_sphere_part_search * yac_point_sphere_part_search_new (
  size_t num_points, double (*coordinates_xyz)[3]) {

  if (num_points == 0) return NULL;

  struct point_sphere_part_search * search = malloc(1 * sizeof(*search));
  struct point_id_xyz * points = malloc(num_points * sizeof(*points));
  search->points = points;

  for (size_t i = 0; i < num_points; ++i) {
    points[i].idx = i;
    points[i].coordinates_xyz[0] = coordinates_xyz[i][0];
    points[i].coordinates_xyz[1] = coordinates_xyz[i][1];
    points[i].coordinates_xyz[2] = coordinates_xyz[i][2];
  }

  struct point_sphere_part_node * tmp_node =
    partition_point_data(
      points, num_points, I_list_tree_min_size, (double[3]){0.0,0.0,1.0});

  search->base_node = *tmp_node;
  free(tmp_node);

  return search;
}

static void search_bnd_circle_I_node(
  struct sphere_part_node * node, struct bounding_circle bnd_circle,
  unsigned ** restrict overlap_cells, size_t * overlap_cells_array_size,
  size_t * restrict num_overlap_cells,
  struct overlaps * search_interval_tree_buffer, double prev_gc_norm_vector[]) {

  if (node->flags & I_IS_INTERVAL_TREE) {

     struct sin_cos_angle base_angle, corrected_inc_angle;
     int big_sum, neg;
     double GCp[3], bVp[3];
     crossproduct_ld(node->gc_norm_vector,
                     bnd_circle.base_vector, GCp);
     crossproduct_ld(GCp, node->gc_norm_vector, bVp);
     normalise_vector(bVp);
     base_angle = get_vector_angle_2(bVp, prev_gc_norm_vector);
     big_sum =
        sum_angles(
          bnd_circle.inc_angle, get_vector_angle_2(bVp, bnd_circle.base_vector),
          &corrected_inc_angle);

     // if the angle is bigger then PI
     if ((big_sum) || (corrected_inc_angle.sin < 0.0))
       corrected_inc_angle = SIN_COS_M_PI;

     struct sin_cos_angle left, right;
     // base_angle - corrected_inc_angle
     neg = sub_angles(base_angle, corrected_inc_angle, &left);
     // base_angle + corrected_inc_angle
     big_sum = sum_angles(base_angle, corrected_inc_angle, &right);

     search_interval_tree_buffer->num_overlaps = 0;

     yac_search_interval_tree(
        node->I.ivt.head_node, node->I.ivt.num_nodes,
        (struct interval){
           .left = compute_angle(left) - (neg?2.0*M_PI:0.0),
           .right = compute_angle(right) + (big_sum?2.0*M_PI:0.0)},
        search_interval_tree_buffer);

     ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                       *num_overlap_cells +
                       search_interval_tree_buffer->num_overlaps);

     for (size_t i = 0; i < search_interval_tree_buffer->num_overlaps;
          ++i) {
        (*overlap_cells)[(*num_overlap_cells)+i] =
           node->I.ivt.head_node[
              search_interval_tree_buffer->overlap_iv[i]].value;
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

//! TODO change to iterative implementation and allocate overlap_cells first
static void search_big_bnd_circle(
  struct sphere_part_node * node, struct bounding_circle bnd_circle,
  unsigned ** restrict overlap_cells, size_t * overlap_cells_array_size,
  size_t * restrict num_overlap_cells,
  struct overlaps * search_interval_tree_buffer, double prev_gc_norm_vector[]) {

  if (node->flags & T_IS_LEAF) {

    ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                      *num_overlap_cells + node->T_size);
    for (size_t i = 0; i < node->T_size; ++i)
      (*overlap_cells)[(*num_overlap_cells) + i] =
         ((unsigned*)(node->T))[i];
    *num_overlap_cells += node->T_size;

  } else {
    search_big_bnd_circle(
       node->T, bnd_circle, overlap_cells, overlap_cells_array_size,
       num_overlap_cells, search_interval_tree_buffer, node->gc_norm_vector);
  }

  if (node->flags & U_IS_LEAF) {

    ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                      *num_overlap_cells + node->U_size);
    for (size_t i = 0; i < node->U_size; ++i)
      (*overlap_cells)[(*num_overlap_cells) + i] =
         ((unsigned*)(node->U))[i];
    *num_overlap_cells += node->U_size;

  } else {
    search_big_bnd_circle(
       node->U, bnd_circle, overlap_cells, overlap_cells_array_size,
       num_overlap_cells, search_interval_tree_buffer, node->gc_norm_vector);
  }

  search_bnd_circle_I_node(
     node, bnd_circle, overlap_cells, overlap_cells_array_size,
     num_overlap_cells, search_interval_tree_buffer, prev_gc_norm_vector);
}

static void search_small_bnd_circle(
  struct sphere_part_node * node, struct bounding_circle bnd_circle,
  unsigned ** restrict overlap_cells, size_t * overlap_cells_array_size,
  size_t * restrict num_overlap_cells,
  struct overlaps * search_interval_tree_buffer, double prev_gc_norm_vector[]) {

   double dot = bnd_circle.base_vector[0] * node->gc_norm_vector[0] +
                bnd_circle.base_vector[1] * node->gc_norm_vector[1] +
                bnd_circle.base_vector[2] * node->gc_norm_vector[2];

   // angle < M_PI_2 + bnd_circle.inc_angle
   if (dot > - bnd_circle.inc_angle.sin) {

      if (node->flags & T_IS_LEAF) {

         ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                           *num_overlap_cells + node->T_size);
         for (size_t i = 0; i < node->T_size; ++i)
           (*overlap_cells)[(*num_overlap_cells) + i] =
              ((unsigned*)(node->T))[i];
         *num_overlap_cells += node->T_size;

      } else {
         search_small_bnd_circle(
            node->T, bnd_circle, overlap_cells, overlap_cells_array_size,
            num_overlap_cells, search_interval_tree_buffer, node->gc_norm_vector);
      }
   }

   // angle > M_PI_2 - bnd_circle.inc_angle
   if (dot < bnd_circle.inc_angle.sin) {

      if (node->flags & U_IS_LEAF) {

         ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                           *num_overlap_cells + node->U_size);
         for (size_t i = 0; i < node->U_size; ++i)
           (*overlap_cells)[(*num_overlap_cells) + i] =
              ((unsigned*)(node->U))[i];
         *num_overlap_cells += node->U_size;

      } else {
         search_small_bnd_circle(
            node->U, bnd_circle, overlap_cells, overlap_cells_array_size,
            num_overlap_cells, search_interval_tree_buffer, node->gc_norm_vector);
      }
   }

   struct sin_cos_angle angle_sum =
      sum_angles_no_check(node->I_angle, bnd_circle.inc_angle);

   // if (I_angle + inc_angle > PI/2) ||
   //    (fabs(angle - M_PI_2) <= (I_angle + inc_angle))
   //
   // assumtion:
   //   I_angle >= 0 && I_angle <= PI/2
   //   inc_angle >= 0 && inc_angle <= PI/2
   //   angle >= 0 && angle <= PI
   //
   //   => I_angle + inc_angle >= 0 && I_angle + inc_angle <= PI
   //
   // I_angle + inc_angle >= PI/2
   //
   // fabs(angle - M_PI_2) <= (I_angle + inc_angle)
   // => sin(fabs(angle - M_PI_2)) <= sin(I_angle + inc_angle)
   //    this is wrong for (I_angle + inc_angle) > PI/2, however that case is
   //    already covered by the first condition
   // => fabs(cos(angle)) <= sin(I_angle + inc_angle)
   if ((compare_angles(angle_sum, SIN_COS_M_PI_2) >= 0) ||
       (fabs(dot) <= angle_sum.sin)) {
      search_bnd_circle_I_node(
         node, bnd_circle, overlap_cells, overlap_cells_array_size,
         num_overlap_cells, search_interval_tree_buffer, prev_gc_norm_vector);
   }
}

static void search_bnd_circle(struct sphere_part_node * node,
                              struct bounding_circle bnd_circle,
                              unsigned ** restrict overlap_cells,
                              size_t * overlap_cells_array_size,
                              size_t * restrict num_overlap_cells,
                              struct overlaps * search_interval_tree_buffer,
                              double prev_gc_norm_vector[]) {

  // if the bounding circle has an angle in the range of [0;PI/2[
  if (compare_angles(bnd_circle.inc_angle, SIN_COS_M_PI_2) < 0)
    search_small_bnd_circle(
      node, bnd_circle, overlap_cells, overlap_cells_array_size,
      num_overlap_cells, search_interval_tree_buffer, prev_gc_norm_vector);
  else
    search_big_bnd_circle(
      node, bnd_circle, overlap_cells, overlap_cells_array_size,
      num_overlap_cells, search_interval_tree_buffer, prev_gc_norm_vector);
}

static void point_search_small_bnd_circle(
   struct point_sphere_part_node * node, struct bounding_circle bnd_circle,
   struct point_id_xyz ** restrict overlap_points,
   size_t * overlap_points_array_size, size_t * restrict num_overlap_points) {

   double dot = node->gc_norm_vector[0]*bnd_circle.base_vector[0] +
                node->gc_norm_vector[1]*bnd_circle.base_vector[1] +
                node->gc_norm_vector[2]*bnd_circle.base_vector[2];

   // angle + inc_angle >= M_PI_2
   if (dot < bnd_circle.inc_angle.sin) {

      if (node->flags & U_IS_LEAF) {

         ENSURE_ARRAY_SIZE(*overlap_points, *overlap_points_array_size,
                           *num_overlap_points + node->U_size);
         memcpy(*overlap_points + *num_overlap_points, node->U,
                node->U_size * sizeof(**overlap_points));
         *num_overlap_points += node->U_size;

      } else {
         point_search_small_bnd_circle(
            node->U, bnd_circle, overlap_points, overlap_points_array_size,
            num_overlap_points);
      }
   }

   // angle - inc_angle < M_PI_2
   if (dot > - bnd_circle.inc_angle.sin) {

      if (node->flags & T_IS_LEAF) {

         ENSURE_ARRAY_SIZE(*overlap_points, *overlap_points_array_size,
                           *num_overlap_points + node->T_size);
         memcpy(*overlap_points + *num_overlap_points, node->T,
                node->T_size * sizeof(**overlap_points));
         *num_overlap_points += node->T_size;

      } else {
         point_search_small_bnd_circle(
            node->T, bnd_circle, overlap_points, overlap_points_array_size,
            num_overlap_points);
      }
   }
}

static void point_search_bnd_circle(
   struct point_sphere_part_search * search, struct bounding_circle bnd_circle,
   struct point_id_xyz ** restrict overlap_points,
   size_t * overlap_points_array_size, size_t * restrict num_overlap_points) {

  struct sin_cos_angle tmp_inc_angle;
  int big_sum = sum_angles(bnd_circle.inc_angle, SIN_COS_TOL, &tmp_inc_angle);
  bnd_circle.inc_angle = tmp_inc_angle;

  // if the inc_angle of the bounding circle is small than M_PI_2
  if (!big_sum && (compare_angles(bnd_circle.inc_angle, SIN_COS_M_PI_2) < 0)) {
    point_search_small_bnd_circle(
      &(search->base_node), bnd_circle, overlap_points,
      overlap_points_array_size, num_overlap_points);
  } else {
    size_t total_size = search->base_node.U_size + search->base_node.T_size;
    ENSURE_ARRAY_SIZE(*overlap_points, *overlap_points_array_size, total_size);
    memcpy(*overlap_points, search->points,
           total_size * sizeof(**overlap_points));
    *num_overlap_points = total_size;
  }
}

static int point_check_bnd_circle(
   struct point_sphere_part_node * node, struct bounding_circle bnd_circle) {

   double dot = node->gc_norm_vector[0]*bnd_circle.base_vector[0] +
                node->gc_norm_vector[1]*bnd_circle.base_vector[1] +
                node->gc_norm_vector[2]*bnd_circle.base_vector[2];

   int ret = 0;

   // angle + inc_angle >= M_PI_2
   if (dot <= bnd_circle.inc_angle.sin) {

      if (node->flags & U_IS_LEAF) {

         struct point_id_xyz * U = (struct point_id_xyz *)(node->U);
         size_t U_size = node->U_size;
         for (size_t i = 0; i < U_size; ++i) {
            double cos_angle =
               U[i].coordinates_xyz[0] * bnd_circle.base_vector[0] +
               U[i].coordinates_xyz[1] * bnd_circle.base_vector[1] +
               U[i].coordinates_xyz[2] * bnd_circle.base_vector[2];
            if (cos_angle > bnd_circle.inc_angle.cos) return 1;
         }

      } else {
         ret = point_check_bnd_circle(node->U, bnd_circle);
      }
   }

   // angle - inc_angle < M_PI_2
   if ((!ret) && (dot > - bnd_circle.inc_angle.sin)) {

      if (node->flags & T_IS_LEAF) {

         struct point_id_xyz * T = (struct point_id_xyz *)(node->T);
         size_t T_size = node->T_size;
         for (size_t i = 0; i < T_size; ++i) {
            double cos_angle =
               T[i].coordinates_xyz[0] * bnd_circle.base_vector[0] +
               T[i].coordinates_xyz[1] * bnd_circle.base_vector[1] +
               T[i].coordinates_xyz[2] * bnd_circle.base_vector[2];
            if (cos_angle > bnd_circle.inc_angle.cos) return 1;
         }

      } else {
         ret = point_check_bnd_circle(node->T, bnd_circle);
      }
   }

   return ret;
}

static inline struct sin_cos_angle get_min_angle(
  struct point_id_xyz * points, size_t num_points, double coordinate_xyz[3]) {

  double max_cos = -1.0;
  size_t min_angle_idx = 0;
  double const tol = 1e-10;

  for (size_t i = 0; i < num_points; ++i) {

    // if the points are nearly identical
    if ((fabs(points[i].coordinates_xyz[0] - coordinate_xyz[0]) < tol) &&
        (fabs(points[i].coordinates_xyz[1] - coordinate_xyz[1]) < tol) &&
        (fabs(points[i].coordinates_xyz[2] - coordinate_xyz[2]) < tol))
      return SIN_COS_ZERO;

    double curr_cos = points[i].coordinates_xyz[0] * coordinate_xyz[0] +
                      points[i].coordinates_xyz[1] * coordinate_xyz[1] +
                      points[i].coordinates_xyz[2] * coordinate_xyz[2];
    if (curr_cos > max_cos) {
      max_cos = curr_cos;
      min_angle_idx = i;
    }
  }

  if (max_cos != -1.0) {

    double cross_ab[3];
    crossproduct_ld(
      points[min_angle_idx].coordinates_xyz, coordinate_xyz, cross_ab);

    return sin_cos_angle_new(sqrt(cross_ab[0]*cross_ab[0] +
                                  cross_ab[1]*cross_ab[1] +
                                  cross_ab[2]*cross_ab[2]), max_cos);
  } else {
    return SIN_COS_M_PI;
  }
}

void yac_point_sphere_part_search_NN(struct point_sphere_part_search * search,
                                     size_t num_points,
                                     double (*coordinates_xyz)[3],
                                     double * cos_angles,
                                     double (**result_coordinates_xyz)[3],
                                     size_t * result_coordinates_xyz_array_size,
                                     unsigned ** local_point_ids,
                                     size_t * local_point_ids_array_size,
                                     size_t * num_local_point_ids) {

  for (size_t i = 0; i < num_points; ++i) num_local_point_ids[i] = 0;
  if (cos_angles != NULL)
    for (size_t i = 0; i < num_points; ++i) cos_angles[i] = -1.0;

  if (search == NULL) return;

  struct point_sphere_part_node * base_node = &(search->base_node);

  struct point_id_xyz * overlap_leaf_points = NULL;
  size_t overlap_leaf_points_array_size = 0;
  size_t num_overlap_leaf_points;

  size_t total_num_local_point_ids = 0;

  double * cos_distances = NULL;
  size_t cos_distances_array_size = 0;

  // coarse search
  for (size_t i = 0; i < num_points; ++i) {

    struct point_sphere_part_node * curr_node = base_node;

    double * curr_coordinates_xyz = coordinates_xyz[i];

    struct bounding_circle bnd_circle;
    bnd_circle.base_vector[0] = curr_coordinates_xyz[0];
    bnd_circle.base_vector[1] = curr_coordinates_xyz[1];
    bnd_circle.base_vector[2] = curr_coordinates_xyz[2];
    bnd_circle.inc_angle.sin = -2.0;

    // get the matching leaf for the current point
    do {
      double dot = curr_node->gc_norm_vector[0]*curr_coordinates_xyz[0] +
                   curr_node->gc_norm_vector[1]*curr_coordinates_xyz[1] +
                   curr_node->gc_norm_vector[2]*curr_coordinates_xyz[2];

      // angle > M_PI_2
      if (dot < yac_angle_tol) {

        if (curr_node->flags & U_IS_LEAF) {
          if (curr_node->U_size > 0) {
            bnd_circle.inc_angle = get_min_angle(
                (struct point_id_xyz*)(curr_node->U),
                (size_t)(curr_node->U_size), curr_coordinates_xyz);
            break;
          } else if (curr_node->flags & T_IS_LEAF) {
            bnd_circle.inc_angle = get_min_angle(
                (struct point_id_xyz*)(curr_node->T),
                (size_t)(curr_node->T_size), curr_coordinates_xyz);
            break;
          } else curr_node = curr_node->T;
        } else curr_node = curr_node->U;
      }

      // angle < M_PI_2
      if (dot > -yac_angle_tol) {

        if (curr_node->flags & T_IS_LEAF) {
          if (curr_node->T_size > 0) {
            bnd_circle.inc_angle = get_min_angle(
                (struct point_id_xyz*)(curr_node->T),
                (size_t)(curr_node->T_size), curr_coordinates_xyz);
            break;
          } else if (curr_node->flags & U_IS_LEAF) {
            bnd_circle.inc_angle = get_min_angle(
                (struct point_id_xyz*)(curr_node->U),
                (size_t)(curr_node->U_size), curr_coordinates_xyz);
            break;
          } else curr_node = curr_node->U;
        } else curr_node = curr_node->T;
      }
    } while (1);

    assert(bnd_circle.inc_angle.sin != -2.0);

    num_overlap_leaf_points = 0;

    // get all leaf points that match the bounding circle of the current point
    point_search_bnd_circle(
      search, bnd_circle, &overlap_leaf_points,
      &overlap_leaf_points_array_size, &num_overlap_leaf_points);

    assert(num_overlap_leaf_points > 0);

    ENSURE_ARRAY_SIZE(
      cos_distances, cos_distances_array_size, num_overlap_leaf_points);

    for (size_t j = 0; j < num_overlap_leaf_points; ++j)
      cos_distances[j] =
        overlap_leaf_points[j].coordinates_xyz[0] * curr_coordinates_xyz[0] +
        overlap_leaf_points[j].coordinates_xyz[1] * curr_coordinates_xyz[1] +
        overlap_leaf_points[j].coordinates_xyz[2] * curr_coordinates_xyz[2];
    double max_cos_distance = cos_distances[0];
    for (size_t j = 1; j < num_overlap_leaf_points; ++j)
      if (cos_distances[j] > max_cos_distance)
        max_cos_distance = cos_distances[j];
    if (cos_angles != NULL) cos_angles[i] = max_cos_distance;

    ENSURE_ARRAY_SIZE(*local_point_ids, *local_point_ids_array_size,
                      total_num_local_point_ids + num_overlap_leaf_points);
    if (result_coordinates_xyz != NULL)
      ENSURE_ARRAY_SIZE(
        *result_coordinates_xyz, *result_coordinates_xyz_array_size,
        total_num_local_point_ids + num_overlap_leaf_points);

    for (size_t j = 0; j < num_overlap_leaf_points; ++j) {
      if (cos_distances[j] == max_cos_distance) {
        (*local_point_ids)[total_num_local_point_ids] =
          (unsigned)(overlap_leaf_points[j].idx);
        if (result_coordinates_xyz != NULL) {
          (*result_coordinates_xyz)[total_num_local_point_ids][0] =
            overlap_leaf_points[j].coordinates_xyz[0];
          (*result_coordinates_xyz)[total_num_local_point_ids][1] =
            overlap_leaf_points[j].coordinates_xyz[1];
          (*result_coordinates_xyz)[total_num_local_point_ids][2] =
            overlap_leaf_points[j].coordinates_xyz[2];
        }
        total_num_local_point_ids++;
        num_local_point_ids[i]++;
      }
    }
  }

  free(cos_distances);
  free(overlap_leaf_points);
}

// this insertion sort is slightly adjusted to the usage in
// yac_point_sphere_part_search_NNN
static void
inv_insertion_sort_dble(double element, double * a, size_t * curr_length) {

  size_t i;
  for (i = 0; i < *curr_length; ++i)
    if (a[i] <= element) break;

  // copy new element into array and move bigger elements one position up
  for (size_t j = *curr_length; j > i; --j) a[j] = a[j-1];
  a[i] = element;

  // increase array length indicator
  ++(*curr_length);
}

static size_t get_leaf_points(
  struct point_id_xyz * leaf_points, struct point_sphere_part_node * node) {

  size_t size;

  if (node->flags & U_IS_LEAF) {
    memcpy(leaf_points, node->U, node->U_size * sizeof(*leaf_points));
    size = node->U_size;
  } else
    size = get_leaf_points(leaf_points, node->U);

  if (node->flags & T_IS_LEAF) {
    memcpy(leaf_points + size, node->T, node->T_size * sizeof(*leaf_points));
    return size + node->T_size;
  } else {
    return size + get_leaf_points(leaf_points + size, node->T);
  }
}

void yac_point_sphere_part_search_NNN(struct point_sphere_part_search * search,
                                      size_t num_points,
                                      double (*coordinates_xyz)[3], unsigned n,
                                      double ** cos_angles,
                                      size_t * cos_angles_array_size,
                                      double (**result_coordinates_xyz)[3],
                                      size_t * result_coordinates_xyz_array_size,
                                      unsigned ** local_point_ids,
                                      size_t * local_point_ids_array_size,
                                      size_t * num_local_point_ids) {

  if (n == 1) {
    if (cos_angles != NULL)
      ENSURE_ARRAY_SIZE(*cos_angles, *cos_angles_array_size, num_points);
    yac_point_sphere_part_search_NN(
      search, num_points, coordinates_xyz, (cos_angles!=NULL)?*cos_angles:NULL,
      result_coordinates_xyz, result_coordinates_xyz_array_size,
      local_point_ids, local_point_ids_array_size, num_local_point_ids);

    size_t total_num_local_points = 0;
    for (size_t i = 0; i < num_points; ++i)
      total_num_local_points += num_local_point_ids[i];

    if (total_num_local_points > num_points) {

      ENSURE_ARRAY_SIZE(*cos_angles, *cos_angles_array_size,
                        total_num_local_points);

      for (size_t i = num_points - 1, offset = total_num_local_points - 1;
           i < num_points; i--) {

        for (size_t j = 0; j < num_local_point_ids[i]; ++j, --offset)
          (*cos_angles)[offset] = (*cos_angles)[i];
      }
    }
    return;
  }

  for (size_t i = 0; i < num_points; ++i) num_local_point_ids[i] = 0;

  if (search == NULL) return;

  struct point_sphere_part_node * base_node = &(search->base_node);

  struct point_id_xyz * overlap_leaf_points = NULL;
  size_t overlap_leaf_points_array_size = 0;
  size_t num_overlap_leaf_points;

  size_t total_num_local_point_ids = 0;

  double * cos_distances = NULL;
  size_t cos_distances_array_size = 0;

  double * n_cos_distances = malloc((n + 1) * sizeof(*n_cos_distances));

  struct point_id_xyz * leaf_points = NULL;
  size_t leaf_points_array_size = 0;

  // coarse search
  for (size_t i = 0; i < num_points; ++i) {

    double * curr_coordinates_xyz = coordinates_xyz[i];

    struct point_sphere_part_node * curr_node = base_node;
    struct point_sphere_part_node * next_node;
    size_t next_node_size;
    unsigned next_node_is_leaf;

    size_t num_leaf_points = 0;

    // get the matching leaf for the current point
    do {

      double dot = curr_node->gc_norm_vector[0]*curr_coordinates_xyz[0] +
                   curr_node->gc_norm_vector[1]*curr_coordinates_xyz[1] +
                   curr_node->gc_norm_vector[2]*curr_coordinates_xyz[2];

      // angle >= M_PI_2
      if (dot <= 0.0) {
        next_node = curr_node->U;
        next_node_size = curr_node->U_size;
        next_node_is_leaf = curr_node->flags & U_IS_LEAF;
      } else {
        next_node = curr_node->T;
        next_node_size = curr_node->T_size;
        next_node_is_leaf = curr_node->flags & T_IS_LEAF;
      }

      // if the next node is too small
      if (next_node_size < n) {
        break;
      // if the next node is a leaf
      } else if (next_node_is_leaf) {
        ENSURE_ARRAY_SIZE(leaf_points, leaf_points_array_size, next_node_size);
        memcpy(leaf_points, next_node, next_node_size * sizeof(*leaf_points));
        num_leaf_points = next_node_size;
        break;
      // go down one level in the tree
      } else {
        curr_node = next_node;
      }
    } while (1);

    if (num_leaf_points == 0) {
      num_leaf_points = curr_node->U_size + curr_node->T_size;
      ENSURE_ARRAY_SIZE(leaf_points, leaf_points_array_size, num_leaf_points);
      num_leaf_points = get_leaf_points(leaf_points, curr_node);
    }

    // compute the n best points of the found leaf points
    // since the angle is in the interval [-Pi;Pi], we can use the dot product
    // instead of computing the exact angle (the higher the dot product the
    // smaller the angle)
    n_cos_distances[0] =
      leaf_points[0].coordinates_xyz[0]*curr_coordinates_xyz[0] +
      leaf_points[0].coordinates_xyz[1]*curr_coordinates_xyz[1] +
      leaf_points[0].coordinates_xyz[2]*curr_coordinates_xyz[2];
    size_t num_distances = 1;
    for (size_t j = 1; j < num_leaf_points; ++j) {
      double cos_distance =
        leaf_points[j].coordinates_xyz[0]*curr_coordinates_xyz[0] +
        leaf_points[j].coordinates_xyz[1]*curr_coordinates_xyz[1] +
        leaf_points[j].coordinates_xyz[2]*curr_coordinates_xyz[2];
      if ((num_distances < n) ||
          (n_cos_distances[num_distances-1] < cos_distance)) {
        inv_insertion_sort_dble(cos_distance, n_cos_distances, &num_distances);
        // we only need n points
        if (num_distances > n) num_distances = n;
      }
    }
    double max_cos_distance = n_cos_distances[num_distances-1];

    num_overlap_leaf_points = 0;

    // get all leaf points that match the bounding circle of the current point
    // cos(acos(cos_distance)+PI/2) = sin(acos(cos_distance))
    //                              = sqrt(1-cos_distance*cos_distance)
    struct bounding_circle bnd_circle;
    bnd_circle.base_vector[0] = curr_coordinates_xyz[0];
    bnd_circle.base_vector[1] = curr_coordinates_xyz[1];
    bnd_circle.base_vector[2] = curr_coordinates_xyz[2];
    bnd_circle.inc_angle = sin_cos_angle_new(
      sqrt(1.0 - MIN(max_cos_distance * max_cos_distance, 1.0)),
      max_cos_distance);

    point_search_bnd_circle(
      search, bnd_circle, &overlap_leaf_points,
      &overlap_leaf_points_array_size, &num_overlap_leaf_points);

    assert(num_overlap_leaf_points > 0);

    ENSURE_ARRAY_SIZE(cos_distances, cos_distances_array_size,
                      num_overlap_leaf_points);

    for (size_t j = 0; j < num_overlap_leaf_points; ++j)
      cos_distances[j] =
        overlap_leaf_points[j].coordinates_xyz[0] * curr_coordinates_xyz[0] +
        overlap_leaf_points[j].coordinates_xyz[1] * curr_coordinates_xyz[1] +
        overlap_leaf_points[j].coordinates_xyz[2] * curr_coordinates_xyz[2];
    n_cos_distances[0] = cos_distances[0];
    num_distances = 1;
    for (size_t j = 1; j < num_overlap_leaf_points; ++j) {
      if ((num_distances < n) ||
          (n_cos_distances[num_distances-1] < cos_distances[j])) {
        inv_insertion_sort_dble(cos_distances[j], n_cos_distances,
                                &num_distances);
        // we only need n points
        if (num_distances > n) num_distances = n;
      }
    }
    max_cos_distance = n_cos_distances[num_distances-1];

    ENSURE_ARRAY_SIZE(*local_point_ids, *local_point_ids_array_size,
                      total_num_local_point_ids + num_overlap_leaf_points);
    if (cos_angles != NULL)
      ENSURE_ARRAY_SIZE(*cos_angles, *cos_angles_array_size,
                        total_num_local_point_ids + num_overlap_leaf_points);
    if (result_coordinates_xyz != NULL)
      ENSURE_ARRAY_SIZE(
        *result_coordinates_xyz, *result_coordinates_xyz_array_size,
        total_num_local_point_ids + num_overlap_leaf_points);

    for (size_t j = 0; j < num_overlap_leaf_points; ++j) {
      if (cos_distances[j] >= max_cos_distance) {
        if (cos_angles != NULL)
          (*cos_angles)[total_num_local_point_ids] = cos_distances[j];
        (*local_point_ids)[total_num_local_point_ids] =
          (unsigned)(overlap_leaf_points[j].idx);
        if (result_coordinates_xyz != NULL) {
          (*result_coordinates_xyz)[total_num_local_point_ids][0] =
            overlap_leaf_points[j].coordinates_xyz[0];
          (*result_coordinates_xyz)[total_num_local_point_ids][1] =
            overlap_leaf_points[j].coordinates_xyz[1];
          (*result_coordinates_xyz)[total_num_local_point_ids][2] =
            overlap_leaf_points[j].coordinates_xyz[2];
        }
        total_num_local_point_ids++;
        num_local_point_ids[i]++;
      }
    }
  }

  free(leaf_points);
  free(cos_distances);
  free(n_cos_distances);
  free(overlap_leaf_points);
}

int yac_point_sphere_part_search_bnd_circle_contains_points(
  struct point_sphere_part_search * search, struct bounding_circle circle) {

  if (search == NULL) return 0;

  return point_check_bnd_circle(&(search->base_node), circle);
}

static void search_point(struct sphere_part_node * node,
                         double point[],
                         unsigned ** overlap_cells,
                         size_t * overlap_cells_array_size,
                         size_t * num_overlap_cells,
                         struct overlaps * search_interval_tree_buffer,
                         double prev_gc_norm_vector[]) {

   double dot = point[0] * node->gc_norm_vector[0] +
                point[1] * node->gc_norm_vector[1] +
                point[2] * node->gc_norm_vector[2];

   // angle < M_PI_2
   if (dot > -yac_angle_tol) {

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

   // angle > M_PI_2
   if (dot < yac_angle_tol) {

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

   // fabs(angle - M_PI_2) <= (node->I_angle)
   // fabs(cos(angle)) <= sin(node->I_angle)
   if (fabs(dot) <= node->I_angle.sin) {

      if (node->flags & I_IS_INTERVAL_TREE) {

         double GCp[3], bVp[3];
         crossproduct_ld(node->gc_norm_vector, point, GCp);
         crossproduct_ld(GCp, node->gc_norm_vector, bVp);
         normalise_vector(bVp);
         double base_angle = get_vector_angle(bVp, prev_gc_norm_vector);

         struct interval search_interval =
          {.left=base_angle, .right=base_angle};

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

   size_t num_cells = (size_t)yac_get_num_grid_cells(grid_data);

   unsigned * temp_search_results = NULL;
   size_t temp_search_results_array_size = 0;
   size_t num_temp_search_results = 0;

   struct grid_cell cell_a, cell_b;
   struct bounding_circle circle_a, circle_b;

   yac_init_grid_cell(&cell_a);
   yac_init_grid_cell(&cell_b);

   struct overlaps search_interval_tree_buffer = {0, 0, NULL};

   unsigned * num_src_per_tgt_cell =
      calloc(num_cells, sizeof(*num_src_per_tgt_cell));
   unsigned * tgt_src_dependencies = NULL;
   size_t tgt_src_dependencies_size = 0;
   size_t total_num_dependencies = 0;

   for (size_t i = 0; i < num_cells; ++i) {

      yac_get_grid_cell2(grid_data, i, &cell_a, &circle_a);

      num_temp_search_results = 0;

      double gc_norm_vector[3] = {0.0,0.0,1.0};

      search_bnd_circle(
         base_node, circle_a, &temp_search_results,
         &temp_search_results_array_size, &num_temp_search_results,
         &search_interval_tree_buffer, gc_norm_vector);

      ENSURE_ARRAY_SIZE(tgt_src_dependencies, tgt_src_dependencies_size,
                        total_num_dependencies + num_temp_search_results);

      for (size_t j = 0; j < num_temp_search_results; ++j) {

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

   double gc_norm_vector[3] = {0.0,0.0,1.0};

   size_t cells_size_ = *cells_size;
   size_t n_cells_ = 0;

   search_bnd_circle(base_node, circle_a, cells, &cells_size_, &n_cells_,
                     &search_interval_tree_buffer, gc_norm_vector);

   size_t n_cells_final = 0;

   struct grid_cell cell_b;
   struct bounding_circle circle_b;

   yac_init_grid_cell(&cell_b);

   for (size_t j = 0; j < n_cells_; ++j) {

      yac_get_grid_cell2(sp_search->grid_data, (*cells)[j],
                         &cell_b, &circle_b);

      if (yac_check_overlap_cells2(cell, circle_a, cell_b, circle_b))
         (*cells)[n_cells_final++] = (*cells)[j];
   }

   *n_cells = (unsigned)n_cells_final;
   *cells_size = (unsigned)cells_size_;

   yac_free_grid_cell(&cell_b);
   free(search_interval_tree_buffer.overlap_iv);
}

static void sphere_part_do_point_search_c(struct grid_search * search,
                                          struct grid * grid_data,
                                          struct dep_list * tgt_to_src_cells) {

   struct sphere_part_search * sp_search = (struct sphere_part_search *)search;
   struct sphere_part_node * base_node = &(sp_search->base_node);

   size_t num_corners = (size_t)yac_get_num_grid_corners(grid_data);

   unsigned * temp_search_results = NULL;
   size_t temp_search_results_array_size = 0;
   size_t num_temp_search_results = 0;

   struct grid_cell cell;
   struct bounding_circle bnd_circle;

   yac_init_grid_cell(&cell);

   struct overlaps search_interval_tree_buffer = {0, 0, NULL};

   unsigned * num_src_per_tgt_corner =
      calloc(num_corners, sizeof(*num_src_per_tgt_corner));
   unsigned * tgt_src_dependencies = NULL;
   size_t tgt_src_dependencies_size = 0;
   size_t total_num_dependencies = 0;

   for (size_t i = 0; i < num_corners; ++i) {
      double point_3d[3];
      LLtoXYZ(yac_get_corner_x_coord(grid_data, i),
              yac_get_corner_y_coord(grid_data, i), point_3d);

      num_temp_search_results = 0;

      double gc_norm_vector[3] = {0.0,0.0,1.0};

      search_point(base_node, point_3d, &temp_search_results,
                   &temp_search_results_array_size, &num_temp_search_results,
                   &search_interval_tree_buffer, gc_norm_vector);

      ENSURE_ARRAY_SIZE(tgt_src_dependencies, tgt_src_dependencies_size,
                        total_num_dependencies + num_temp_search_results);

      for (size_t j = 0; j < num_temp_search_results; ++j) {

         yac_get_grid_cell2(sp_search->grid_data, temp_search_results[j],
                            &cell, &bnd_circle);

         if (yac_point_in_cell2(point_3d, cell, bnd_circle)) {

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

static int compare_uint (const void * a, const void * b) {
  return ( *(unsigned*)a - *(unsigned*)b );
}

static void sphere_part_do_point_search_c2(struct grid_search * search,
                                           double (*coordinates_xyz)[3],
                                           unsigned num_points,
                                           struct dep_list * tgt_to_src_cells) {

   struct sphere_part_search * sp_search = (struct sphere_part_search *)search;
   struct sphere_part_node * base_node = &(sp_search->base_node);

   unsigned * temp_search_results = NULL;
   size_t temp_search_results_array_size = 0;
   size_t num_temp_search_results = 0;

   struct grid_cell cell;
   struct bounding_circle bnd_circle;

   yac_init_grid_cell(&cell);

   struct overlaps search_interval_tree_buffer = {0, 0, NULL};

   unsigned * num_src_per_tgt_corner =
      calloc(num_points, sizeof(*num_src_per_tgt_corner));
   unsigned * tgt_src_dependencies = NULL;
   size_t tgt_src_dependencies_size = 0;
   size_t total_num_dependencies = 0;

   for (size_t i = 0; i < num_points; ++i) {

      double * curr_coordinates_xyz = coordinates_xyz[i];

      num_temp_search_results = 0;

      double gc_norm_vector[3] = {0.0,0.0,1.0};

      search_point(base_node, curr_coordinates_xyz, &temp_search_results,
                   &temp_search_results_array_size, &num_temp_search_results,
                   &search_interval_tree_buffer, gc_norm_vector);

      ENSURE_ARRAY_SIZE(tgt_src_dependencies, tgt_src_dependencies_size,
                        total_num_dependencies + num_temp_search_results);

      size_t num_matches = 0;

      for (size_t j = 0; j < num_temp_search_results; ++j) {

         yac_get_grid_cell2(sp_search->grid_data, temp_search_results[j],
                            &cell, &bnd_circle);

         if (yac_point_in_cell2(curr_coordinates_xyz, cell, bnd_circle)) {

            tgt_src_dependencies[total_num_dependencies + num_matches] =
               temp_search_results[j];
            num_matches++;
         }
      }
      qsort(
        tgt_src_dependencies + total_num_dependencies, num_matches,
        sizeof(*tgt_src_dependencies), compare_uint);
      num_src_per_tgt_corner[i] = num_matches;
      total_num_dependencies += num_matches;
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

static void sphere_part_do_point_search_c3(struct grid_search * search,
                                           double (*coordinates_xyz)[3],
                                           unsigned num_points,
                                           struct dep_list * tgt_to_src_cells,
                                           struct points * points) {

  yac_grid_search_utils_do_point_search_c3(
    search, coordinates_xyz, num_points, tgt_to_src_cells, points);
}

static void sphere_part_do_point_search_p(struct grid_search * search,
                                          struct grid * grid_data,
                                          struct dep_list * target_to_src_points) {

   struct sphere_part_search * sp_search = (struct sphere_part_search *)search;

   yac_grid_search_utils_do_point_search_p(search, sp_search->grid_data, grid_data,
                                           target_to_src_points);
}

static void sphere_part_do_point_search_p2(struct grid_search * search,
                                           double (*coordinates_xyz)[3],
                                           unsigned num_points,
                                           struct dep_list * target_to_src_points) {

   struct sphere_part_search * sp_search = (struct sphere_part_search *)search;

   yac_grid_search_utils_do_point_search_p2(
      search, sp_search->grid_data, coordinates_xyz, num_points,
      target_to_src_points);
}

static void sphere_part_do_point_search_p3(struct grid_search * search,
                                           double (*coordinates_xyz)[3],
                                           unsigned num_points,
                                           struct dep_list * target_to_src_points,
                                           struct points * points) {

   yac_grid_search_utils_do_point_search_p3(
      search, coordinates_xyz, num_points, target_to_src_points, points);
}


static void sphere_part_do_point_search_p4 (struct grid_search * search,
                                            double coordinate_xyz[3],
                                            unsigned * n_points,
                                            unsigned * points_size,
                                            unsigned ** points) {

   struct sphere_part_search * sp_search = (struct sphere_part_search *)search;

   yac_grid_search_utils_do_point_search_p4(
      search, sp_search->grid_data, coordinate_xyz, n_points, points_size,
      points);
}

static void sphere_part_do_bnd_circle_search (struct grid_search * search,
                                              struct bounding_circle * bnd_circles,
                                              unsigned num_bnd_circles,
                                              struct dep_list * bnd_to_cells) {

   struct sphere_part_search * sp_search = (struct sphere_part_search *)search;
   struct sphere_part_node * base_node = &(sp_search->base_node);

   unsigned * temp_search_results = NULL;
   size_t temp_search_results_array_size = 0;
   size_t num_temp_search_results = 0;

   struct grid_cell cell;

   yac_init_grid_cell(&cell);

   struct overlaps search_interval_tree_buffer = {0, 0, NULL};

   unsigned * num_cells_per_bnd =
      calloc(num_bnd_circles, sizeof(*num_cells_per_bnd));
   unsigned * bnd_to_cells_dependencies = NULL;
   size_t bnd_to_cells_dependencies_size = 0;
   size_t total_num_dependencies = 0;

   for (unsigned i = 0; i < num_bnd_circles; ++i) {

      num_temp_search_results = 0;

      double gc_norm_vector[3] = {0.0,0.0,1.0};

      search_bnd_circle(
         base_node, bnd_circles[i], &temp_search_results,
         &temp_search_results_array_size, &num_temp_search_results,
         &search_interval_tree_buffer, gc_norm_vector);

      ENSURE_ARRAY_SIZE(bnd_to_cells_dependencies,
                        bnd_to_cells_dependencies_size,
                        total_num_dependencies + num_temp_search_results);

      for (size_t j = 0; j < num_temp_search_results; ++j) {

         struct bounding_circle cell_bnd_circle;
         yac_get_grid_cell2(sp_search->grid_data, temp_search_results[j], &cell,
                            &cell_bnd_circle);

         // if the bounding circle of the current cell overlaps with the current
         // bounding circle
         if (yac_extents_overlap(&cell_bnd_circle, bnd_circles + i)) {
            bnd_to_cells_dependencies[total_num_dependencies++] =
               temp_search_results[j];
            num_cells_per_bnd[i]++;
         }
      }
   }

   yac_free_grid_cell(&cell);
   free(temp_search_results);
   free(search_interval_tree_buffer.overlap_iv);

   bnd_to_cells_dependencies =
    realloc(bnd_to_cells_dependencies, total_num_dependencies *
            sizeof(bnd_to_cells_dependencies));
   yac_init_dep_list(bnd_to_cells);
   yac_set_dependencies(bnd_to_cells, num_bnd_circles, num_cells_per_bnd,
                        bnd_to_cells_dependencies);
}

static void free_sphere_part_tree (struct sphere_part_node tree) {

   // free I_list
   if (tree.flags & I_IS_INTERVAL_TREE)
      free(tree.I.ivt.head_node);

   if ((tree.flags & U_IS_LEAF) == 0) {
      free_sphere_part_tree(*(struct sphere_part_node*)(tree.U));
      free(tree.U);
   }

   if ((tree.flags & T_IS_LEAF) == 0) {
      free_sphere_part_tree(*(struct sphere_part_node*)(tree.T));
      free(tree.T);
   }
}

static void free_point_sphere_part_tree (struct point_sphere_part_node * tree) {

   if ((tree->flags & U_IS_LEAF) == 0) {
      free_point_sphere_part_tree(tree->U);
      free(tree->U);
   }

   if ((tree->flags & T_IS_LEAF) == 0) {
      free_point_sphere_part_tree(tree->T);
      free(tree->T);
   }
}

static void delete_sphere_part_search(struct grid_search * search) {

   struct sphere_part_search * sp_search = (struct sphere_part_search *)search;

   free_sphere_part_tree(sp_search->base_node);

   free(sp_search->local_cell_ids);
   free(sp_search);
}

void yac_delete_point_sphere_part_search(
   struct point_sphere_part_search * search) {

   if (search == NULL) return;

   free_point_sphere_part_tree(&(search->base_node));
   free(search->points);
   free(search);
}

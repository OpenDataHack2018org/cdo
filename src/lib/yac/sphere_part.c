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
//#include "grid_search.h"
//#include "grid_search_utils.h"


enum {
   I_list_tree_min_size = 2, //!< make I list into tree when list is
                             //!<  larger than this
};

struct sphere_part_search {

  //struct grid_search_vtable * vtable;
   struct sphere_part_node base_node;
   struct grid * grid_data;
};

/*
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
   .do_point_search_p     = sphere_part_do_point_search_p,
   .do_point_search_p2    = sphere_part_do_point_search_p2,
   .do_point_search_p3    = sphere_part_do_point_search_p3,
   .do_point_search_p4    = sphere_part_do_point_search_p4,
   .do_bnd_circle_search  = sphere_part_do_bnd_circle_search,
   .delete_grid_search    = delete_sphere_part_search
};
*/

struct point_sphere_part_search {

   struct point_sphere_part_node base_node;
   double * coordinates_xyz;
   unsigned * inclusive_node_size;
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

static void init_point_sphere_part_node(struct point_sphere_part_node * node) {

   node->flags = 0;
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

      yac_get_grid_cell(grid, local_cell_ids[i], &cell);

      for (unsigned j = 0; j < cell.num_corners; ++j) {

        balance_point[0] += cell.coordinates_xyz[0+j*3];
        balance_point[1] += cell.coordinates_xyz[1+j*3];
        balance_point[2] += cell.coordinates_xyz[2+j*3];
      }
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
      double dot = bnd_circle.base_vector[0] * parent_node->gc_norm_vector[0] +
                   bnd_circle.base_vector[1] * parent_node->gc_norm_vector[1] +
                   bnd_circle.base_vector[2] * parent_node->gc_norm_vector[2];
      double sin_inc_angle =
        (bnd_circle.inc_angle >= M_PI_2)?1.0:sin(bnd_circle.inc_angle);

      // if the point intersects with the great circle
      // (angle >= M_PI_2 && angle <= M_PI_2 + bnd_circle.inc_angle) ||
      // (angle <= M_PI_2 && angle >= M_PI_2 - bnd_circle.inc_angle)
      if ((dot <= 0.0 && dot >= - sin_inc_angle) ||
          (dot >= 0.0 && dot <=   sin_inc_angle)) {

         // angle = fabs(angle - M_PI_2);
         double angle = fabs(asin(dot));

         // add to list I
         ENSURE_ARRAY_SIZE(I, I_array_size, I_size + 1);
         I[I_size++] = local_cell_ids[i];
         max_inc_angle = MAX(max_inc_angle, angle + bnd_circle.inc_angle);

      // angle > M_PI_2
      } else if (dot < 0.0) {

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

   if (max_inc_angle >= M_PI_2) {
     parent_node->I_angle = M_PI_2;
     parent_node->sin_I_angle = 1.0;
     parent_node->cos_I_angle = 0.0;
   } else {
     parent_node->I_angle = max_inc_angle;
     parent_node->sin_I_angle = sin(max_inc_angle);
     parent_node->cos_I_angle = cos(max_inc_angle);
   }

   if (I_size > 0) {

      if (I_size > I_list_tree_min_size) {

         assert(sizeof(struct interval_node) > sizeof(unsigned));
         parent_node->I.ivt.head_node
            = (void *)(I = realloc(I, I_size * sizeof(struct interval_node)));
         parent_node->I.ivt.num_nodes = I_size;

         for (i = I_size - 1 ; i < I_size; --i) {

            double GCp[3], bVp[3], base_angle, corrected_inc_angle;
            unsigned cell_idx = I[i];

            yac_get_grid_cell2(grid, I[i], &cell, &bnd_circle);
            crossproduct_ld(parent_node->gc_norm_vector,
                            bnd_circle.base_vector, GCp);
            crossproduct_ld(GCp, parent_node->gc_norm_vector, bVp);
            normalise_vector(bVp);
            base_angle = get_vector_angle(bVp, prev_gc_norm_vector);
            corrected_inc_angle = bnd_circle.inc_angle +
                                  get_vector_angle(bVp, bnd_circle.base_vector);

            parent_node->I.ivt.head_node[i].range.left  =
              base_angle - corrected_inc_angle;
            parent_node->I.ivt.head_node[i].range.right =
              base_angle + corrected_inc_angle;
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

static void partition_point_data (
  double * coordinates_xyz, unsigned * local_point_ids, unsigned num_points,
  unsigned threshold, struct point_sphere_part_node * parent_node,
  double prev_gc_norm_vector[]) {

  double balance_point[3] = {0.0,0.0,0.0};

  // compute balance point

  for (unsigned i = 0; i < num_points; ++i) {

    balance_point[0] += coordinates_xyz[0+local_point_ids[i]*3];
    balance_point[1] += coordinates_xyz[1+local_point_ids[i]*3];
    balance_point[2] += coordinates_xyz[2+local_point_ids[i]*3];
  }

  normalise_vector(balance_point);

  // compute the great circle that partitions the data in half (more or less)

  crossproduct_ld(balance_point, prev_gc_norm_vector,
                  parent_node->gc_norm_vector);
  normalise_vector(parent_node->gc_norm_vector);

  // partition data into points that are on the great circle and points that
  // are on side of the circle

  unsigned * U = NULL; // list of points on one side of the great circle
  unsigned * T = NULL; // list of points on that other side of the great circle

  unsigned U_size = 0, U_array_size = 0;
  unsigned T_size = 0, T_array_size = 0;

  for (unsigned i = 0; i < num_points; ++i) {

    // angle between the point of the point coordinates and the great circle
    // plane
    // acos(dot(gc_norm_vector, point_xyz)) =
    //   angle(gc_norm_vector, point_xyz)
    // acos(dot(gc_norm_vector, point_xyz)) - PI/2 =
    //   angle(gc_plane, point_xyz)
    // acos(x) = PI/2 - asin(x)
    // - asin(dot(gc_norm_vector, point_xyz)) =
    //   angle(gc_plane, point_xyz)
    // since we just want to know whether (angle >= 0.0), we can test
    // (dot <= 0.0) instead
    double dot = coordinates_xyz[0+local_point_ids[i]*3]*
                 parent_node->gc_norm_vector[0] +
                 coordinates_xyz[1+local_point_ids[i]*3]*
                 parent_node->gc_norm_vector[1] +
                 coordinates_xyz[2+local_point_ids[i]*3]*
                 parent_node->gc_norm_vector[2];

    // (angle >= M_PI_2)
    if (dot <= 0.0) {

      // add to list U
      ENSURE_ARRAY_SIZE(U, U_array_size, U_size + 1);
      U[U_size++] = local_point_ids[i];

    } else {

      // add to list T
      ENSURE_ARRAY_SIZE(T, T_array_size, T_size + 1);
      T[T_size++] = local_point_ids[i];
    }
  }

  parent_node->U_size = U_size;
  parent_node->T_size = T_size;

  // check whether the lists are small enough (if not -> partition again)
  if (U_size <= threshold) {

    parent_node->U = realloc(U, U_size * sizeof(unsigned));
    parent_node->flags |= U_IS_LEAF;

  } else {

    parent_node->U = get_sphere_part_node();
    partition_point_data(coordinates_xyz, U, U_size, threshold,
                         parent_node->U, parent_node->gc_norm_vector);
    free(U);
  }

  if (T_size <= threshold) {

    parent_node->T = realloc(T, T_size * sizeof(unsigned));
    parent_node->flags |= T_IS_LEAF;

  } else {

    parent_node->T = get_sphere_part_node();
    partition_point_data(coordinates_xyz, T, T_size, threshold,
                         parent_node->T, parent_node->gc_norm_vector);
    free(T);
  }
}
/*
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

   for (unsigned i = 0; i < num_grid_cells; ++i)
      local_cell_ids[i] = i;

   partition_data(grid, local_cell_ids, num_grid_cells, I_list_tree_min_size,
                  &(search->base_node), gc_norm_vector);

   free(local_cell_ids);

   return (struct grid_search *)search;
}
*/
struct point_sphere_part_search * yac_point_sphere_part_search_new (
  unsigned num_points, double * x_coordinates, double * y_coordinates) {

  if (num_points == 0) return NULL;

  struct point_sphere_part_search * search = malloc(1 * sizeof(*search));

  double gc_norm_vector[3] = {0.0,0.0,1.0};

  unsigned * local_point_ids = malloc(num_points * sizeof(*local_point_ids));
  double * coordinates_xyz = malloc(3 * num_points * sizeof(*coordinates_xyz));

  for (unsigned i = 0; i < num_points; ++i) {
    local_point_ids[i] = i;
    LLtoXYZ(x_coordinates[i], y_coordinates[i], coordinates_xyz + i * 3);
  }

  init_point_sphere_part_node(&(search->base_node));

  partition_point_data(coordinates_xyz, local_point_ids, num_points,
                       I_list_tree_min_size, &(search->base_node),
                       gc_norm_vector);

  search->coordinates_xyz = coordinates_xyz;

  free(local_point_ids);

  search->inclusive_node_size = NULL;

  return search;
}

void *cdo_point_sphere_part_search_new(unsigned num_points, double *coordinates_xyz) {

  if (num_points == 0) return NULL;

  struct point_sphere_part_search * search = malloc(1 * sizeof(*search));

  double gc_norm_vector[3] = {0.0,0.0,1.0};

  unsigned * local_point_ids = malloc(num_points * sizeof(*local_point_ids));

  for (unsigned i = 0; i < num_points; ++i) local_point_ids[i] = i;

  init_point_sphere_part_node(&(search->base_node));

  partition_point_data(coordinates_xyz, local_point_ids, num_points,
                       I_list_tree_min_size, &(search->base_node),
                       gc_norm_vector);

  search->coordinates_xyz = coordinates_xyz;

  free(local_point_ids);

  search->inclusive_node_size = NULL;

  return (void*)search;
}

static void search_bnd_circle(struct sphere_part_node * node,
                              double * bnd_circle_base_vector, double inc_angle,
                              double sin_inc_angle, double cos_inc_angle,
                              unsigned ** restrict overlap_cells,
                              unsigned * overlap_cells_array_size,
                              unsigned * restrict num_overlap_cells,
                              struct overlaps * search_interval_tree_buffer,
                              double prev_gc_norm_vector[]) {


   double dot = bnd_circle_base_vector[0] * node->gc_norm_vector[0] +
                bnd_circle_base_vector[1] * node->gc_norm_vector[1] +
                bnd_circle_base_vector[2] * node->gc_norm_vector[2];

   // angle < M_PI_2 + bnd_circle.inc_angle
   if (dot > - sin_inc_angle) {

      if (node->flags & T_IS_LEAF) {

         ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                           *num_overlap_cells + node->T_size);
         for (unsigned i = 0; i < node->T_size; ++i)
           (*overlap_cells)[(*num_overlap_cells) + i] =
              ((unsigned*)(node->T))[i];
         *num_overlap_cells += node->T_size;

      } else {
         search_bnd_circle(node->T, bnd_circle_base_vector, inc_angle,
                           sin_inc_angle, cos_inc_angle, overlap_cells,
                           overlap_cells_array_size, num_overlap_cells,
                           search_interval_tree_buffer, node->gc_norm_vector);
      }
   }

   // angle > M_PI_2 - bnd_circle.inc_angle
   if (dot < sin_inc_angle) {

      if (node->flags & U_IS_LEAF) {

         ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                           *num_overlap_cells + node->U_size);
         for (unsigned i = 0; i < node->U_size; ++i)
           (*overlap_cells)[(*num_overlap_cells) + i] =
              ((unsigned*)(node->U))[i];
         *num_overlap_cells += node->U_size;

      } else {
         search_bnd_circle(node->U, bnd_circle_base_vector, inc_angle,
                           sin_inc_angle, cos_inc_angle, overlap_cells,
                           overlap_cells_array_size, num_overlap_cells,
                           search_interval_tree_buffer, node->gc_norm_vector);
      }
   }

   // fabs(angle - M_PI_2) <= (bnd_circle.inc_angle + node->I_angle)
   // sin(x+y) = sin(x)*cos(y)+cos(x)*sin(y)
   if ((node->I_angle + inc_angle >= M_PI_2) ||
       (fabs(dot) <= sin_inc_angle*node->cos_I_angle +
                     cos_inc_angle*node->sin_I_angle)) {

      if (node->flags & I_IS_INTERVAL_TREE) {

         double GCp[3], bVp[3], base_angle, corrected_inc_angle;
         crossproduct_ld(node->gc_norm_vector,
                         bnd_circle_base_vector, GCp);
         crossproduct_ld(GCp, node->gc_norm_vector, bVp);
         normalise_vector(bVp);
         base_angle = get_vector_angle(bVp, prev_gc_norm_vector);
         corrected_inc_angle =
            inc_angle + get_vector_angle(bVp, bnd_circle_base_vector);

         search_interval_tree_buffer->num_overlaps = 0;

         yac_search_interval_tree(
            node->I.ivt.head_node, node->I.ivt.num_nodes,
            (struct interval){.left = base_angle - corrected_inc_angle,
                              .right = base_angle + corrected_inc_angle},
            search_interval_tree_buffer);

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

static void point_search_bnd_circle(struct point_sphere_part_node * node,
                                    double * bnd_circle_base_vector,
                                    double sin_inc_angle,
                                    unsigned ** restrict overlap_cells,
                                    unsigned * overlap_cells_array_size,
                                    unsigned * restrict num_overlap_cells) {

   double dot = node->gc_norm_vector[0]*bnd_circle_base_vector[0] +
                node->gc_norm_vector[1]*bnd_circle_base_vector[1] +
                node->gc_norm_vector[2]*bnd_circle_base_vector[2];

   // angle + inc_angle >= M_PI_2
   if (dot <= sin_inc_angle) {

      if (node->flags & U_IS_LEAF) {

         ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                           *num_overlap_cells + node->U_size);
         for (unsigned i = 0; i < node->U_size; ++i)
            (*overlap_cells)[(*num_overlap_cells) + i] =
               ((unsigned*)(node->U))[i];
         *num_overlap_cells += node->U_size;

      } else {
         point_search_bnd_circle(node->U, bnd_circle_base_vector, sin_inc_angle,
                                 overlap_cells, overlap_cells_array_size,
                                 num_overlap_cells);
      }
   }

   // angle - inc_angle < M_PI_2
   if (dot > - sin_inc_angle) {

      if (node->flags & T_IS_LEAF) {

         ENSURE_ARRAY_SIZE(*overlap_cells, *overlap_cells_array_size,
                           *num_overlap_cells + node->T_size);
         for (unsigned i = 0; i < node->T_size; ++i)
            (*overlap_cells)[(*num_overlap_cells) + i] =
               ((unsigned*)(node->T))[i];
         *num_overlap_cells += node->T_size;

      } else {
         point_search_bnd_circle(node->T, bnd_circle_base_vector, sin_inc_angle,
                                 overlap_cells, overlap_cells_array_size,
                                 num_overlap_cells);
      }
   }
}

static int point_check_bnd_circle(struct point_sphere_part_node * node,
                                  double * bnd_circle_base_vector,
                                  double sin_inc_angle, double cos_inc_angle,
                                  struct point_sphere_part_search * search) {

   double dot = node->gc_norm_vector[0]*bnd_circle_base_vector[0] +
                node->gc_norm_vector[1]*bnd_circle_base_vector[1] +
                node->gc_norm_vector[2]*bnd_circle_base_vector[2];

   int ret = 0;

   // angle + inc_angle >= M_PI_2
   if (dot <= sin_inc_angle) {

      if (node->flags & U_IS_LEAF) {

         unsigned * node_point_indices = (unsigned*)(node->U);
         unsigned num_node_point_indices = node->U_size;
         double * coordinates_xyz = search->coordinates_xyz;
         for (unsigned i = 0; i < num_node_point_indices; ++i) {
            double cos_angle = coordinates_xyz[3*node_point_indices[i]+0] *
                               bnd_circle_base_vector[0] +
                               coordinates_xyz[3*node_point_indices[i]+1] *
                               bnd_circle_base_vector[1] +
                               coordinates_xyz[3*node_point_indices[i]+2] *
                               bnd_circle_base_vector[2];
            if (cos_angle > cos_inc_angle) return 1;
         }

      } else {
         ret = point_check_bnd_circle(node->U, bnd_circle_base_vector,
                                      sin_inc_angle, cos_inc_angle, search);
      }
   }

   // angle - inc_angle < M_PI_2
   if ((!ret) && (dot > - sin_inc_angle)) {

      if (node->flags & T_IS_LEAF) {

         unsigned * node_point_indices = (unsigned*)(node->T);
         unsigned num_node_point_indices = node->T_size;
         double * coordinates_xyz = search->coordinates_xyz;
         for (unsigned i = 0; i < num_node_point_indices; ++i) {
            double cos_angle = coordinates_xyz[3*node_point_indices[i]+0] *
                               bnd_circle_base_vector[0] +
                               coordinates_xyz[3*node_point_indices[i]+1] *
                               bnd_circle_base_vector[1] +
                               coordinates_xyz[3*node_point_indices[i]+2] *
                               bnd_circle_base_vector[2];
            if (cos_angle > cos_inc_angle) return 1;
         }

      } else {
         ret = point_check_bnd_circle(node->T, bnd_circle_base_vector,
                                      sin_inc_angle, cos_inc_angle, search);
      }
   }

   return ret;
}

void yac_point_sphere_part_search_NN(struct point_sphere_part_search * search,
                                     unsigned num_points, double * x_coordinates,
                                     double * y_coordinates,
                                     double * cos_angles,
                                     unsigned ** local_point_ids,
                                     unsigned * local_point_ids_array_size,
                                     unsigned * num_local_point_ids) {

  for (unsigned i = 0; i < num_points; ++i) num_local_point_ids[i] = 0;
  if (cos_angles != NULL)
    for (unsigned i = 0; i < num_points; ++i) cos_angles[i] = -1.0;

  if (search == NULL) return;

  struct point_sphere_part_node * base_node = &(search->base_node);

  unsigned * overlap_leaf_points = NULL;
  unsigned overlap_leaf_points_array_size = 0;
  unsigned num_overlap_leaf_points;

  unsigned total_num_local_point_ids = 0;

  double * cos_distances = NULL;
  unsigned cos_distances_array_size = 0;

  // coarse search
  for (unsigned i = 0; i < num_points; ++i) {

    struct point_sphere_part_node * curr_node = base_node;

    double coordinates_xyz[3];
    LLtoXYZ(x_coordinates[i], y_coordinates[i], coordinates_xyz);

    unsigned  * leaf_points = NULL;
    unsigned num_leaf_points = -1;

    // get the matching leaf for the current point
    do {

      double dot = curr_node->gc_norm_vector[0]*coordinates_xyz[0] +
                   curr_node->gc_norm_vector[1]*coordinates_xyz[1] +
                   curr_node->gc_norm_vector[2]*coordinates_xyz[2];

      // angle >= M_PI_2
      if (dot <= 0.0) {

        if (curr_node->flags & U_IS_LEAF) {
          if (curr_node->U_size > 0) {
            leaf_points = (unsigned*)(curr_node->U);
            num_leaf_points = curr_node->U_size;
          } else if (curr_node->flags & T_IS_LEAF) {
            leaf_points = (unsigned*)(curr_node->T);
            num_leaf_points = curr_node->T_size;
          } else curr_node = curr_node->T;
        } else curr_node = curr_node->U;

      } else {

        if (curr_node->flags & T_IS_LEAF) {
          if (curr_node->T_size > 0) {
            leaf_points = (unsigned*)(curr_node->T);
            num_leaf_points = curr_node->T_size;
          } else if (curr_node->flags & U_IS_LEAF) {
            leaf_points = (unsigned*)(curr_node->U);
            num_leaf_points = curr_node->U_size;
          } else curr_node = curr_node->U;
        } else curr_node = curr_node->T;
      }
    } while (num_leaf_points == -1);

    assert(num_leaf_points > 0);

    // compute the best point of the found leaf
    // since the angle is in the interval [-Pi;Pi], we can use the dot product
    // instead of computing the exact angle (the higher the dot product the
    // smaller the angle)
    double cos_distance =
      search->coordinates_xyz[0+leaf_points[0]*3]*coordinates_xyz[0] +
      search->coordinates_xyz[1+leaf_points[0]*3]*coordinates_xyz[1] +
      search->coordinates_xyz[2+leaf_points[0]*3]*coordinates_xyz[2];
    for (unsigned j = 1; j < num_leaf_points; ++j) {
      double tmp_cos_distance =
        search->coordinates_xyz[0+leaf_points[j]*3]*coordinates_xyz[0] +
        search->coordinates_xyz[1+leaf_points[j]*3]*coordinates_xyz[1] +
        search->coordinates_xyz[2+leaf_points[j]*3]*coordinates_xyz[2];
      if (tmp_cos_distance > cos_distance) cos_distance = tmp_cos_distance;
    }

    num_overlap_leaf_points = 0;

    // get all leaf points that match the bounding circle of the current point
    // cos(acos(cos_distance)+PI/2) = sin(acos(cos_distance))
    //                              = sqrt(1-cos_distance*cos_distance)
    point_search_bnd_circle(base_node, coordinates_xyz,
                            sqrt(1.0 - MIN(cos_distance * cos_distance, 1.0)),
                            &overlap_leaf_points,
                            &overlap_leaf_points_array_size,
                            &num_overlap_leaf_points);

    assert(num_overlap_leaf_points > 0);

    ENSURE_ARRAY_SIZE(cos_distances, cos_distances_array_size,
                      num_overlap_leaf_points);

    for (unsigned j = 0; j < num_overlap_leaf_points; ++j)
      cos_distances[j] =
        search->coordinates_xyz[0+overlap_leaf_points[j]*3]*coordinates_xyz[0] +
        search->coordinates_xyz[1+overlap_leaf_points[j]*3]*coordinates_xyz[1] +
        search->coordinates_xyz[2+overlap_leaf_points[j]*3]*coordinates_xyz[2];
    double max_cos_distance = cos_distances[0];
    for (unsigned j = 1; j < num_overlap_leaf_points; ++j)
      if (cos_distances[j] > max_cos_distance)
        max_cos_distance = cos_distances[j];
    if (cos_angles != NULL) cos_angles[i] = max_cos_distance;

    ENSURE_ARRAY_SIZE(*local_point_ids, *local_point_ids_array_size,
                      total_num_local_point_ids + num_overlap_leaf_points);

    for (unsigned j = 0; j < num_overlap_leaf_points; ++j) {
      if (cos_distances[j] == max_cos_distance) {
        (*local_point_ids)[total_num_local_point_ids++] =
          overlap_leaf_points[j];
        num_local_point_ids[i]++;
      }
    }
  }

  free(cos_distances);
  free(overlap_leaf_points);
}

void cdo_point_sphere_part_search_NN(void * search_container,
                                     unsigned num_points, double * x_coordinates,
                                     double * y_coordinates,
                                     double * cos_angles,
                                     unsigned ** local_point_ids,
                                     unsigned * local_point_ids_array_size,
                                     unsigned * num_local_point_ids)
{
  struct point_sphere_part_search * search = (struct point_sphere_part_search *) search_container;
  yac_point_sphere_part_search_NN(search, num_points, x_coordinates, y_coordinates, cos_angles,
                                  local_point_ids, local_point_ids_array_size, num_local_point_ids);
}

static unsigned
get_max_free_depth_mask(struct point_sphere_part_node * node, unsigned mask) {

  unsigned ret = mask;

  if (!(node->flags & U_IS_LEAF))
    ret |= get_max_free_depth_mask(node->U, mask << 1);
  if (!(node->flags & T_IS_LEAF))
    ret |= get_max_free_depth_mask(node->T, mask << 1);

  return ret;
}

// this insertion sort is slightly adjusted to the usage in
// yac_point_sphere_part_search_NNN
static void
inv_insertion_sort_dble(double element, double * a, unsigned * curr_length) {

  unsigned i;
  for (i = 0; i < *curr_length; ++i)
    if (a[i] <= element) break;

  // copy new element into array and move bigger elements one position up
  for (unsigned j = *curr_length; j > i; --j) a[j] = a[j-1];
  a[i] = element;

  // increase array length indicator
  ++(*curr_length);
}

static unsigned
get_inclusive_node_size(struct point_sphere_part_node * node,
                        unsigned * inclusive_node_size, unsigned access_mask,
                        unsigned depth_mask) {

  unsigned U_size = 0;
  unsigned T_size = 0;

  if (node->flags & U_IS_LEAF)
    U_size = node->U_size;
  else
    U_size = get_inclusive_node_size(
      node->U, inclusive_node_size, access_mask << 1, (depth_mask << 1) | 1);

  if (node->flags & T_IS_LEAF)
    T_size = node->T_size;
  else
    T_size = get_inclusive_node_size(
      node->T, inclusive_node_size, (access_mask | 1) << 1, (depth_mask << 1) | 1);

  inclusive_node_size[access_mask + (depth_mask << 1)] = U_size;
  inclusive_node_size[access_mask + 1 + (depth_mask << 1)] = T_size;

  return U_size + T_size;
}

static unsigned *
generate_inclusive_node_size(struct point_sphere_part_node * base_node) {

  unsigned * inclusive_node_size =
    calloc(2 * get_max_free_depth_mask(base_node, 1),
           sizeof(*inclusive_node_size));

  unsigned total_num_points =
    get_inclusive_node_size(base_node, inclusive_node_size, 0, 0);

  return inclusive_node_size;
}

static unsigned get_leaf_points(unsigned * leaf_points,
                                struct point_sphere_part_node * node) {

  unsigned size;

  if (node->flags & U_IS_LEAF) {
    for (unsigned i = 0; i < node->U_size; ++i)
      leaf_points[i] = ((unsigned*)(node->U))[i];
    size = node->U_size;
  } else
    size = get_leaf_points(leaf_points, node->U);

  if (node->flags & T_IS_LEAF) {
    for (unsigned i = 0; i < node->T_size; ++i)
      leaf_points[i+size] = ((unsigned*)(node->T))[i];
    return size + node->T_size;
  } else {
    return size + get_leaf_points(leaf_points + size, node->T);
  }
}

void yac_point_sphere_part_search_NNN(struct point_sphere_part_search * search,
                                      unsigned num_points, double * x_coordinates,
                                      double * y_coordinates, unsigned n,
                                      double ** cos_angles,
                                      unsigned * cos_angles_array_size,
                                      unsigned ** local_point_ids,
                                      unsigned * local_point_ids_array_size,
                                      unsigned * num_local_point_ids) {

  if (n == 1) {
    if (cos_angles != NULL)
      ENSURE_ARRAY_SIZE(*cos_angles, *cos_angles_array_size, num_points);
    yac_point_sphere_part_search_NN(search, num_points, x_coordinates,
                                    y_coordinates,
                                    (cos_angles!=NULL)?*cos_angles:NULL,
                                    local_point_ids, local_point_ids_array_size,
                                    num_local_point_ids);

    unsigned total_num_local_points = 0;
    for (unsigned i = 0; i < num_points; ++i)
      total_num_local_points += num_local_point_ids[i];

    if (total_num_local_points > num_points) {

      ENSURE_ARRAY_SIZE(*cos_angles, *cos_angles_array_size,
                        total_num_local_points);

      for (unsigned i = num_points - 1, offset = total_num_local_points - 1;
           i < num_points; i--) {

        for (unsigned j = 0; j < num_local_point_ids[i]; ++j, --offset)
          (*cos_angles)[offset] = (*cos_angles)[i];
      }
    }
    return;
  }

  for (unsigned i = 0; i < num_points; ++i)
    num_local_point_ids[i] = 0;

  if (search == NULL) return;

  if (search->inclusive_node_size == NULL)
    search->inclusive_node_size =
      generate_inclusive_node_size(&(search->base_node));

  struct point_sphere_part_node * base_node = &(search->base_node);

  unsigned * overlap_leaf_points = NULL;
  unsigned overlap_leaf_points_array_size = 0;
  unsigned num_overlap_leaf_points;

  unsigned total_num_local_point_ids = 0;

  double * cos_distances = NULL;
  unsigned cos_distances_array_size = 0;

  double * n_cos_distances = malloc((n + 1) * sizeof(*n_cos_distances));

  unsigned  * leaf_points = NULL;
  unsigned leaf_points_array_size = 0;

  // coarse search
  for (unsigned i = 0; i < num_points; ++i) {

    struct point_sphere_part_node * curr_node = base_node, * prev_node;

    double coordinates_xyz[3];
    LLtoXYZ(x_coordinates[i], y_coordinates[i], coordinates_xyz);

    unsigned access_mask = 0;
    unsigned depth_mask = 0;
    unsigned prev_access_mask;
    unsigned prev_depth_mask;
    unsigned num_leaf_points = 0;

    // get the matching leaf for the current point
    do {

      prev_node = curr_node;
      prev_access_mask = access_mask;
      prev_depth_mask = depth_mask;

      double dot = curr_node->gc_norm_vector[0]*coordinates_xyz[0] +
                   curr_node->gc_norm_vector[1]*coordinates_xyz[1] +
                   curr_node->gc_norm_vector[2]*coordinates_xyz[2];

      // angle >= M_PI_2
      if (dot <= 0.0) {
        if (!(curr_node->flags & U_IS_LEAF)) {
          access_mask = access_mask << 1;
          depth_mask = (depth_mask << 1) | 1;
          curr_node = curr_node->U;
        } else {
          if (curr_node->U_size >= n) {
            ENSURE_ARRAY_SIZE(leaf_points, leaf_points_array_size,
                              curr_node->U_size);
            memcpy(leaf_points, curr_node->U,
                   curr_node->U_size * sizeof(*leaf_points));
            num_leaf_points = curr_node->U_size;
          }
          break;
        }
      } else {
        if (!(curr_node->flags & T_IS_LEAF)) {
          access_mask = (access_mask | 1) << 1;
          depth_mask = (depth_mask << 1) | 1;
          curr_node = curr_node->T;
        } else {
          if (curr_node->T_size >= n) {
            ENSURE_ARRAY_SIZE(leaf_points, leaf_points_array_size,
                              curr_node->T_size);
            memcpy(leaf_points, curr_node->T,
                   curr_node->T_size * sizeof(*leaf_points));
            num_leaf_points = curr_node->T_size;
          }
          break;
        }
      }

      if (search->inclusive_node_size[access_mask + (depth_mask << 1)] < n) {
        curr_node = prev_node;
        access_mask = prev_access_mask;
        depth_mask = prev_depth_mask;
        break;
      }
    } while (1);

    if (num_leaf_points == 0) {
      num_leaf_points =
        search->inclusive_node_size[access_mask + (depth_mask << 1)];
      ENSURE_ARRAY_SIZE(leaf_points, leaf_points_array_size, num_leaf_points);
      num_leaf_points = get_leaf_points(leaf_points, curr_node);
    }

    // compute the n best points of the found leaf points
    // since the angle is in the interval [-Pi;Pi], we can use the dot product
    // instead of computing the exact angle (the higher the dot product the
    // smaller the angle)
    n_cos_distances[0] =
      search->coordinates_xyz[0+leaf_points[0]*3]*coordinates_xyz[0] +
      search->coordinates_xyz[1+leaf_points[0]*3]*coordinates_xyz[1] +
      search->coordinates_xyz[2+leaf_points[0]*3]*coordinates_xyz[2];
    unsigned num_distances = 1;
    for (unsigned j = 1; j < num_leaf_points; ++j) {
      double cos_distance =
        search->coordinates_xyz[0+leaf_points[j]*3]*coordinates_xyz[0] +
        search->coordinates_xyz[1+leaf_points[j]*3]*coordinates_xyz[1] +
        search->coordinates_xyz[2+leaf_points[j]*3]*coordinates_xyz[2];
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
    point_search_bnd_circle(base_node, coordinates_xyz,
                            sqrt(1.0 - MIN(max_cos_distance * max_cos_distance,
                                 1.0)), &overlap_leaf_points,
                            &overlap_leaf_points_array_size,
                            &num_overlap_leaf_points);

    assert(num_overlap_leaf_points > 0);

    ENSURE_ARRAY_SIZE(cos_distances, cos_distances_array_size,
                      num_overlap_leaf_points);

    for (unsigned j = 0; j < num_overlap_leaf_points; ++j)
      cos_distances[j] =
        search->coordinates_xyz[0+overlap_leaf_points[j]*3]*coordinates_xyz[0] +
        search->coordinates_xyz[1+overlap_leaf_points[j]*3]*coordinates_xyz[1] +
        search->coordinates_xyz[2+overlap_leaf_points[j]*3]*coordinates_xyz[2];
    n_cos_distances[0] = cos_distances[0];
    num_distances = 1;
    for (unsigned j = 1; j < num_overlap_leaf_points; ++j) {
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

    for (unsigned j = 0; j < num_overlap_leaf_points; ++j) {
      if (cos_distances[j] >= max_cos_distance) {
        if (cos_angles != NULL)
          (*cos_angles)[total_num_local_point_ids] = cos_distances[j];
        (*local_point_ids)[total_num_local_point_ids++] =
          overlap_leaf_points[j];
        num_local_point_ids[i]++;
      }
    }
  }

  free(leaf_points);
  free(cos_distances);
  free(n_cos_distances);
  free(overlap_leaf_points);
}

void cdo_point_sphere_part_search_NNN(void * search_container,
                                      unsigned num_points, double * x_coordinates,
                                      double * y_coordinates, unsigned n,
                                      double ** cos_angles,
                                      unsigned * cos_angles_array_size,
                                      unsigned ** local_point_ids,
                                      unsigned * local_point_ids_array_size,
                                      unsigned * num_local_point_ids)
{
  struct point_sphere_part_search * search = (struct point_sphere_part_search *) search_container;
  yac_point_sphere_part_search_NNN(search, num_points, x_coordinates, y_coordinates, n,
                                   cos_angles, cos_angles_array_size, local_point_ids, local_point_ids_array_size, num_local_point_ids);
}
 
/*
int yac_point_sphere_part_search_bnd_circle_contains_points(
  struct point_sphere_part_search * search,
  struct reduced_bounding_circle circle) {

  // cos(acos(cos_angle)+PI/2) = sin(acos(cos_angle))
  //                           = sqrt(1-cos_angle*cos_angle)
  return
    point_check_bnd_circle(
      &(search->base_node), circle.base_vector,
      sqrt(1.0 - MIN(circle.cos_inc_angle * circle.cos_inc_angle, 1.0)),
      circle.cos_inc_angle, search);
}
*/
static void search_point(struct sphere_part_node * node,
                         double point[],
                         unsigned ** overlap_cells,
                         unsigned * overlap_cells_array_size,
                         unsigned * num_overlap_cells,
                         struct overlaps * search_interval_tree_buffer,
                         double prev_gc_norm_vector[]) {

   double dot = point[0] * node->gc_norm_vector[0] +
                point[1] * node->gc_norm_vector[1] +
                point[2] * node->gc_norm_vector[2];

   // angle < M_PI_2
   if (dot > 0.0) {

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
   if (dot < 0.0) {

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
   if (fabs(dot) <= node->sin_I_angle) {

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
/*
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

      double sin_inc_angle =
        (circle_a.inc_angle >= M_PI_2)?1.0:sin(circle_a.inc_angle);

      search_bnd_circle(base_node, circle_a.base_vector, circle_a.inc_angle,
                        sin_inc_angle, cos(circle_a.inc_angle),
                        &temp_search_results, &temp_search_results_array_size,
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

   double sin_inc_angle =
      (circle_a.inc_angle >= M_PI_2)?1.0:sin(circle_a.inc_angle);

   search_bnd_circle(base_node, circle_a.base_vector, circle_a.inc_angle,
                     sin_inc_angle, cos(circle_a.inc_angle), cells, cells_size,
                     n_cells, &search_interval_tree_buffer, gc_norm_vector);

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

static void sphere_part_do_bnd_circle_search (struct grid_search * search,
                                              struct bounding_circle * bnd_circles,
                                              unsigned num_bnd_circles,
                                              struct dep_list * bnd_to_cells) {

   struct sphere_part_search * sp_search = (struct sphere_part_search *)search;
   struct sphere_part_node * base_node = &(sp_search->base_node);

   unsigned * temp_search_results = NULL;
   unsigned temp_search_results_array_size = 0;
   unsigned num_temp_search_results = 0;

   struct grid_cell cell;

   yac_init_grid_cell(&cell);

   struct overlaps search_interval_tree_buffer = {0, 0, NULL};

   unsigned * num_cells_per_bnd = calloc(num_bnd_circles,
                                            sizeof(*num_cells_per_bnd));
   unsigned * bnd_to_cells_dependencies = NULL;
   unsigned bnd_to_cells_dependencies_size = 0;
   unsigned total_num_dependencies = 0;

   for (unsigned i = 0; i < num_bnd_circles; ++i) {

      num_temp_search_results = 0;

      double gc_norm_vector[3] = {0.0,0.0,1.0};

      double sin_inc_angle =
        (bnd_circles[i].inc_angle >= M_PI_2)?1.0:sin(bnd_circles[i].inc_angle);

      search_bnd_circle(base_node, bnd_circles[i].base_vector,
                        bnd_circles[i].inc_angle, sin_inc_angle,
                        cos(bnd_circles[i].inc_angle), &temp_search_results,
                        &temp_search_results_array_size,
                        &num_temp_search_results, &search_interval_tree_buffer,
                        gc_norm_vector);

      ENSURE_ARRAY_SIZE(bnd_to_cells_dependencies,
                        bnd_to_cells_dependencies_size,
                        total_num_dependencies + num_temp_search_results);

      for (unsigned j = 0; j < num_temp_search_results; ++j) {

         yac_get_grid_cell(sp_search->grid_data, temp_search_results[j], &cell);

         unsigned k;
         // for all corners of the cell
         for (k = 0; k < cell.num_corners; ++k)
           // if the current corner is inside the bounding circle
           if (yac_point_in_bounding_circle_vec(cell.coordinates_xyz+k*3,
                                                bnd_circles + i)) break;
         if (k != cell.num_corners) {
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
*/
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

static void free_point_sphere_part_tree (struct point_sphere_part_node * tree) {

   if ((tree->flags & U_IS_LEAF) == 0)
      free_point_sphere_part_tree(tree->U);
   free(tree->U);

   if ((tree->flags & T_IS_LEAF) == 0)
      free_point_sphere_part_tree(tree->T);
   free(tree->T);
}
/*
static void delete_sphere_part_search(struct grid_search * search) {

   struct sphere_part_search * sp_search = (struct sphere_part_search *)search;

   free_sphere_part_tree(&(sp_search->base_node));

   free(sp_search);
}
*/
void yac_delete_point_sphere_part_search(
   struct point_sphere_part_search * search) {

   if (search == NULL) return;

   free_point_sphere_part_tree(&(search->base_node));
   free(search->coordinates_xyz);
   free(search->inclusive_node_size);
   free(search);
}

void cdo_delete_point_sphere_part_search(
   void * search_container) {

   struct point_sphere_part_search * search = (struct point_sphere_part_search *) search_container;
   if (search == NULL) return;

   free_point_sphere_part_tree(&(search->base_node));
   free(search->coordinates_xyz);
   free(search->inclusive_node_size);
   free(search);
}

#ifdef YAC_DEBUG_SPHERE_PART
struct sphere_part_node * yac_get_sphere_part_tree(struct grid_search * search) {

   return &(((struct sphere_part_search *)search)->base_node);
}
#endif

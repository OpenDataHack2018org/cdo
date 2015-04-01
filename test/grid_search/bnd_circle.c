/**
 * @file bnd_circle.c
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

#include <math.h>
#include <stdlib.h>

#include "geometry.h"
#include "utils.h"
#include "grid.h"
#include "ensure_array_size.h"

/** \file bnd_circle.c
 *  \brief Set of functions to calculate a bounding circle around a certain set of points.
 *
 **/

static double const tol = 1.0e-12;

double yac_get_point_angle(struct point * a, struct point * b) {

  /* method 1*/
  // return acos(sin(a->lat)*sin(b->lat)+cos(a->lat)*cos(b->lat)*cos(a->lon-b->lon));

  /* method 2 - haversine formula (supposed to be more accurate)*/
  double t_1 = sin((a->lat - b->lat)*0.5);
  double t_2 = sin((a->lon - b->lon)*0.5);
  return 2.0 * asin(sqrt(t_1*t_1 + cos(a->lat)*cos(b->lat)*t_2*t_2));
}

// this routine computes the distance between two points in rad
// be careful when using this routine for small angles, because it is very inaccurate
double yac_get_dist(double a_lon, double a_lat, double b_lon, double b_lat) {

   // this formula is very inaccurate for small distances
   double temp;

   temp = sin(a_lat)*sin(b_lat)+cos(a_lat)*cos(b_lat)*cos(a_lon - b_lon);

   if (temp >= 1) return 0;
   if (temp <= -1) return M_PI;

   return acos(temp);

   /*
   double delta_lon;

   delta_lon = a_lon - b_lon;

   double sin_lat_a, sin_lat_b, cos_lat_a, cos_lat_b, sin_delta_lon, cos_delta_lon;

   sin_lat_a = sin(a_lat);
   sin_lat_b = sin(b_lat);
   cos_lat_a = cos(a_lat);
   cos_lat_b = cos(b_lat);
   sin_delta_lon = sin(delta_lon);
   cos_delta_lon = cos(delta_lon);

   double temp_a, temp_b;

   temp_a = cos_lat_a * sin_delta_lon;
   temp_b = cos_lat_b * sin_lat_a - sin_lat_b * cos_lat_a * cos_delta_lon;

   return atan2(sqrt(temp_a * temp_a + temp_b * temp_b),
                sin_lat_b * sin_lat_a + cos_lat_b * cos_lat_a * cos_delta_lon);*/

   /*double temp_a, temp_b;

   temp_a = sin((a_lat - b_lat) / 2.0);
   temp_a *= temp_a;

   temp_b = sin((a_lon - b_lon) / 2.0);
   temp_b *= temp_b;

   return 2.0 * asin(sqrt(temp_a + cos(a_lat) * cos(b_lat) * temp_b));*/
}

static void normalise_vector(double v[]) {

   double norm;

   norm = 1.0 / sqrt(v[0]*v[0] +
                     v[1]*v[1] +
                     v[2]*v[2]);

   v[0] *= norm;
   v[1] *= norm;
   v[2] *= norm;
}

// computes the circumscribe circle of a quad on the sphere
// it is assumed that the edges are circles of longitude and latitude
void yac_get_cell_circumscribe_circle_reg_quad(
   double a[3], double b[3], double c[3], double d[3],
   struct bounding_circle * bnd_circle) {

   yac_get_cell_circumscribe_circle_unstruct_triangle(a, b, c, bnd_circle);
}

// computes the bounding circle of a quad on the sphere
// it is assumed that the edges are circles of longitude and latitude
void yac_get_cell_bounding_circle_reg_quad(
   double a[3], double b[3], double c[3], double d[3],
   struct bounding_circle * bnd_circle) {

   yac_get_cell_circumscribe_circle_unstruct_triangle(a, b, c, bnd_circle);
   bnd_circle->inc_angle += tol;
}

// computes the circumscribe circle of a triangle on the sphere
// it is assumed that all edges are great circles
void yac_get_cell_circumscribe_circle_unstruct_triangle(
   double a[3], double b[3], double c[3],
   struct bounding_circle * bnd_circle) {

   double ab[3] = {a[0]-b[0], a[1]-b[1], a[2]-b[2]},
          ac[3] = {b[0]-c[0], b[1]-c[1], b[2]-c[2]};

   // it is assumed that the angles of a triangle do not get too small...
   // crossproduct_ld(ab, ac, bnd_circle->base_vector);
   crossproduct_d(ab, ac, bnd_circle->base_vector);
   normalise_vector(bnd_circle->base_vector);

   int biggest_component_index = 0;
   // find biggest component of base_vector
   if (fabs(bnd_circle->base_vector[0]) <
       fabs(bnd_circle->base_vector[1]))
      biggest_component_index = 1;
   if (fabs(bnd_circle->base_vector[biggest_component_index]) <
       fabs(bnd_circle->base_vector[2]))
      biggest_component_index = 2;

   if ((bnd_circle->base_vector[biggest_component_index] > 0) ^
       (a[biggest_component_index] > 0)) {
      bnd_circle->base_vector[0] = -bnd_circle->base_vector[0];
      bnd_circle->base_vector[1] = -bnd_circle->base_vector[1];
      bnd_circle->base_vector[2] = -bnd_circle->base_vector[2];
   }

   XYZtoLL(bnd_circle->base_vector, bnd_circle->base_point+0,
           bnd_circle->base_point+1);
   bnd_circle->inc_angle = get_vector_angle(bnd_circle->base_vector, a);;
}

// computes the bounding circle of a triangle on the sphere
// it is assumed that all edges are great circles
void yac_get_cell_bounding_circle_unstruct_triangle(
   double a[3], double b[3], double c[3],
   struct bounding_circle * bnd_circle) {

   double ab[3] = {a[0]-b[0], a[1]-b[1], a[2]-b[2]},
          bc[3] = {b[0]-c[0], b[1]-c[1], b[2]-c[2]},
          ca[3] = {c[0]-a[0], c[1]-a[1], c[2]-a[2]};

   double length_ab = ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2],
          length_bc = bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2],
          length_ca = ca[0]*ca[0] + ca[1]*ca[1] + ca[2]*ca[2];

   int flag = 0;

   if (length_ab > length_bc) flag |= 1;
   if (length_ab > length_ca) flag |= 2;
   if (length_bc > length_ca) flag |= 4;

   double * longest_edge_start;
   double * longest_edge;
   double longest_edge_length;
   double * other_point;

   // find longest edge
   switch (flag) {
      case (7): // edge AB is the longest one
      case (3):

         longest_edge_start = b;
         longest_edge = ab;
         longest_edge_length = length_ab;
         other_point = c;
         break;

      case (6): // edge BC is the longest one
      case (4):

         longest_edge_start = c;
         longest_edge = bc;
         longest_edge_length = length_bc;
         other_point = a;
         break;

      case (1): // edge CA is the longest one
      case (0):

         longest_edge_start = a;
         longest_edge = ca;
         longest_edge_length = length_ca;
         other_point = b;
         break;

      default:
         longest_edge_start = NULL;
         longest_edge = NULL;
         longest_edge_length = 0.0;
         other_point = NULL;
         yac_internal_abort_message("internal error", __FILE__, __LINE__);
         // this function should never reach this point...
   };

   // compute middle point
   double longest_edge_middle_point[3] = {longest_edge_start[0] + longest_edge[0]*0.5,
                                          longest_edge_start[1] + longest_edge[1]*0.5,
                                          longest_edge_start[2] + longest_edge[2]*0.5};

   double other_point_distance = (longest_edge_middle_point[0] - other_point[0]) *
                                 (longest_edge_middle_point[0] - other_point[0]) +
                                 (longest_edge_middle_point[1] - other_point[1]) *
                                 (longest_edge_middle_point[1] - other_point[1]) +
                                 (longest_edge_middle_point[2] - other_point[2]) *
                                 (longest_edge_middle_point[2] - other_point[2]);

   // if the other point would be included
   if (longest_edge_length >= 4.0*other_point_distance) {

      bnd_circle->base_vector[0] = longest_edge_middle_point[0];
      bnd_circle->base_vector[1] = longest_edge_middle_point[1];
      bnd_circle->base_vector[2] = longest_edge_middle_point[2];
      XYZtoLL(bnd_circle->base_vector, bnd_circle->base_point+0,
              bnd_circle->base_point+1);
      bnd_circle->inc_angle = get_vector_angle(bnd_circle->base_vector,
                                               longest_edge_start);
      bnd_circle->inc_angle += tol;

   // else compute circumscribe circle for all three points
   }  else {
      yac_get_cell_circumscribe_circle_unstruct_triangle(a, b, c, bnd_circle);
      bnd_circle->inc_angle += tol;
   }
}

void yac_get_cell_bounding_circle(struct grid_cell cell,
                                  struct bounding_circle * bnd_circle) {

   unsigned i;
   double middle_point[3];

   middle_point[0] = 0;
   middle_point[1] = 0;
   middle_point[2] = 0;

   // compute the coordinates in rad and 3d
   for (i = 0; i < cell.num_corners; ++i) {

      middle_point[0] += (cell.coordinates_xyz + 3 * i)[0];
      middle_point[1] += (cell.coordinates_xyz + 3 * i)[1];
      middle_point[2] += (cell.coordinates_xyz + 3 * i)[2];
   }

   normalise_vector(middle_point);

   double max_angle;

   max_angle = 0;

   // compute the angle required for the bounding circle
   for (i = 0; i < cell.num_corners-1; ++i) {

      double edge_middle_point[3];

      edge_middle_point[0] = (cell.coordinates_xyz + 3 * i)[0] +
                             (cell.coordinates_xyz + 3 * i + 3)[0];
      edge_middle_point[1] = (cell.coordinates_xyz + 3 * i)[1] +
                             (cell.coordinates_xyz + 3 * i + 3)[1];
      edge_middle_point[2] = (cell.coordinates_xyz + 3 * i)[2] +
                             (cell.coordinates_xyz + 3 * i + 3)[2];

      normalise_vector(edge_middle_point);

      double corner_corner_angle;

      corner_corner_angle = get_vector_angle((cell.coordinates_xyz + 3 * i),
                                             (cell.coordinates_xyz + 3 * i + 3));

      double edge_middle_angle;

      edge_middle_angle = get_vector_angle(edge_middle_point, middle_point);

      max_angle = MAX(max_angle, edge_middle_angle + corner_corner_angle * 0.5);
   }

   bnd_circle->base_vector[0] = middle_point[0];
   bnd_circle->base_vector[1] = middle_point[1];
   bnd_circle->base_vector[2] = middle_point[2];

   bnd_circle->inc_angle = max_angle;

   XYZtoLL(bnd_circle->base_vector, bnd_circle->base_point+0,
           bnd_circle->base_point+1);
}

// based on http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/
static void rotate_vector(double rotated_vector[], double rotation_axis[],
                          double vector[], double angle) {

   double sin_angle, cos_angle;

   sin_angle = sin(angle);
   cos_angle = cos(angle);

   double u, v, w;

   u = rotation_axis[0];
   v = rotation_axis[1];
   w = rotation_axis[2];

   double x, y, z;

   x = vector[0];
   y = vector[1];
   z = vector[2];

   double temp;

   temp = (u*x+v*y+w*z)*(1-cos_angle);

   rotated_vector[0] = u*temp+x*cos_angle+(-w*y+v*z)*sin_angle;
   rotated_vector[1] = v*temp+y*cos_angle+( w*x-u*z)*sin_angle;
   rotated_vector[2] = w*temp+z*cos_angle+(-v*x+u*y)*sin_angle;

   normalise_vector(rotated_vector);
}

static void compute_norm_vector(double norm_vector[], double a[],
                                double b[]) {

   norm_vector[0] = a[1]*b[2] - a[2]*b[1];
   norm_vector[1] = a[2]*b[0] - a[0]*b[2];
   norm_vector[2] = a[0]*b[1] - a[1]*b[0];

   normalise_vector(norm_vector);
}

static void merge_bounding_circles(struct bounding_circle * dest_circle,
                                   struct bounding_circle * circle) {

   double base_vector_angle;

   // check whether one circle is within the other
   base_vector_angle = get_vector_angle(circle->base_vector,
                                            dest_circle->base_vector);

   if (dest_circle->inc_angle >= base_vector_angle + circle->inc_angle) {

      return;

   } else if (circle->inc_angle >= base_vector_angle + dest_circle->inc_angle) {

      *dest_circle = *circle;
      return;

   } else {

      double norm_vector[3];

      compute_norm_vector(norm_vector, dest_circle->base_vector, circle->base_vector);

      double angle;

      angle = (base_vector_angle + circle->inc_angle - dest_circle->inc_angle) * 0.5;

      double rotated_vector[3];

      rotate_vector(rotated_vector, norm_vector, dest_circle->base_vector, angle);

      dest_circle->inc_angle = MIN((base_vector_angle +
                                    circle->inc_angle +
                                    dest_circle->inc_angle) * 0.5,
                                    M_PI);

      dest_circle->base_vector[0] = rotated_vector[0];
      dest_circle->base_vector[1] = rotated_vector[1];
      dest_circle->base_vector[2] = rotated_vector[2];
   }
}

void yac_get_grid_bounding_circle(struct grid * grid,
                                  struct bounding_circle * bnd_circle) {

   unsigned num_grid_cells;

   num_grid_cells = yac_get_num_grid_cells(grid);

   if (num_grid_cells == 0) {

      bnd_circle->inc_angle = 0;
      return;
   }

   struct grid_cell cell;

   yac_init_grid_cell(&cell);

   yac_get_grid_cell(grid, 0, &cell);

   yac_get_cell_bounding_circle(cell, bnd_circle);

   unsigned i;
   struct bounding_circle curr_bnd_circle;
   double base_vector_angle;

   for (i = 1; i < num_grid_cells; ++i) {

      yac_get_grid_cell(grid, i, &cell);
      yac_get_cell_bounding_circle(cell, &curr_bnd_circle);

      // check whether the bounding circle of the cell is within the current
      // global bounding circle
      base_vector_angle = get_vector_angle(bnd_circle->base_vector,
                                           curr_bnd_circle.base_vector);

      if (bnd_circle->inc_angle <
          base_vector_angle + curr_bnd_circle.inc_angle) {

          merge_bounding_circles(bnd_circle, &curr_bnd_circle);
      }
   }

   yac_free_grid_cell(&cell);

   XYZtoLL(bnd_circle->base_vector, bnd_circle->base_point+0,
           bnd_circle->base_point+1);
}

unsigned yac_extents_overlap(struct bounding_circle * extent_a,
                             struct bounding_circle * extent_b) {

   double angle;

   angle = get_vector_angle(extent_a->base_vector, extent_b->base_vector);

   return angle - tol <= extent_a->inc_angle + extent_b->inc_angle;
}

unsigned yac_point_in_bounding_circle(struct point point,
                                      struct bounding_circle * bnd_circle) {

   double point_vector[3];

   LLtoXYZ(point.lon, point.lat, point_vector);

   return bnd_circle->inc_angle + tol >=
          get_vector_angle(bnd_circle->base_vector, point_vector);
}

unsigned yac_point_in_bounding_circle_vec(double point_vector[3],
                                          struct bounding_circle * bnd_circle) {

   return bnd_circle->inc_angle + tol >=
          get_vector_angle(bnd_circle->base_vector, point_vector);
}

void yac_get_matching_grid_cells(struct grid * grid, struct bounding_circle extent,
                                 struct grid_cell ** matching_cells,
                                 unsigned * curr_matching_cells_array_size,
                                 unsigned ** local_ids, unsigned * curr_local_ids_array_size,
                                 unsigned * num_matching_cells, unsigned offset) {

   unsigned num_grid_cells;

   num_grid_cells = yac_get_num_grid_cells(grid);

   unsigned i, j;
   struct grid_cell cell;
   double dist;

   *num_matching_cells = 0;

   yac_init_grid_cell(&cell);

   for (i = 0; i < num_grid_cells; ++i) {

      yac_get_grid_cell(grid, i, &cell);

      for (j = 0; j < cell.num_corners; ++j) {

         dist = yac_get_dist(cell.coordinates_x[j], cell.coordinates_y[j],
                             extent.base_point[0], extent.base_point[1]);

         if (dist - tol <= extent.inc_angle) {

            if (matching_cells != NULL) {
               ENSURE_ARRAY_SIZE(*matching_cells, *curr_matching_cells_array_size,
                                 *num_matching_cells + offset + 1);
               yac_copy_grid_cell(cell, (*matching_cells) + offset +
                                  *num_matching_cells);
            }

            if (local_ids != NULL) {
               ENSURE_ARRAY_SIZE(*local_ids, *curr_local_ids_array_size,
                                 *num_matching_cells + offset + 1);
               (*local_ids)[offset + *num_matching_cells] = i;
            }

            ++*num_matching_cells;

            break;
         }
      }
   }
   yac_free_grid_cell(&cell);
}

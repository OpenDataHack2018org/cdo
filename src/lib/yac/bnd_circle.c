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

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "utils.h"
#include "grid.h"
#include "grid_search.h"
#include "ensure_array_size.h"

/** \file bnd_circle.c
 *  \brief Set of functions to calculate a bounding circle around a certain set of points.
 *
 **/

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

   bnd_circle->inc_angle =
      sum_angles_no_check(
         get_vector_angle_2(
            bnd_circle->base_vector, a), SIN_COS_TOL);
}

static double inline get_sin_vector_angle(
  double a[3], double b[3]) {

  double cross_ab[3];
  crossproduct_ld(a, b, cross_ab);

  return sqrt(cross_ab[0]*cross_ab[0] +
              cross_ab[1]*cross_ab[1] +
              cross_ab[2]*cross_ab[2]);
}

// computes the bounding circle of a triangle on the sphere
// it is assumed that all edges are great circles
void yac_get_cell_bounding_circle_unstruct_triangle(
  double a[3], double b[3], double c[3],
  struct bounding_circle * bnd_circle) {

  double middle_point[3];

  middle_point[0] = a[0] + b[0] + c[0];
  middle_point[1] = a[1] + b[1] + c[1];
  middle_point[2] = a[2] + b[2] + c[2];

  normalise_vector(middle_point);

  double cos_angles[3] = {middle_point[0] * a[0] +
                          middle_point[1] * a[1] +
                          middle_point[2] * a[2],
                          middle_point[0] * b[0] +
                          middle_point[1] * b[1] +
                          middle_point[2] * b[2],
                          middle_point[0] * c[0] +
                          middle_point[1] * c[1] +
                          middle_point[2] * c[2]};

  struct sin_cos_angle inc_angle;

  // find the biggest angle

  if (cos_angles[0] < cos_angles[1]) {
    if (cos_angles[0] < cos_angles[2]) {
      inc_angle =
        sin_cos_angle_new(get_sin_vector_angle(middle_point, a), cos_angles[0]);
    } else {
      inc_angle =
        sin_cos_angle_new(get_sin_vector_angle(middle_point, c), cos_angles[2]);
    }
  } else {
    if (cos_angles[1] < cos_angles[2]) {
      inc_angle =
        sin_cos_angle_new(get_sin_vector_angle(middle_point, b), cos_angles[1]);
    } else {
      inc_angle =
        sin_cos_angle_new(get_sin_vector_angle(middle_point, c), cos_angles[2]);
    }
  }

  bnd_circle->base_vector[0] = middle_point[0];
  bnd_circle->base_vector[1] = middle_point[1];
  bnd_circle->base_vector[2] = middle_point[2];
  bnd_circle->inc_angle = sum_angles_no_check(inc_angle, SIN_COS_TOL);
}

void yac_get_cell_bounding_circle(struct grid_cell cell,
                                  struct bounding_circle * bnd_circle) {

   double middle_point[3];

   middle_point[0] = 0;
   middle_point[1] = 0;
   middle_point[2] = 0;

   // compute the coordinates in rad and 3d
   for (unsigned i = 0; i < cell.num_corners; ++i) {

      middle_point[0] += (cell.coordinates_xyz + 3 * i)[0];
      middle_point[1] += (cell.coordinates_xyz + 3 * i)[1];
      middle_point[2] += (cell.coordinates_xyz + 3 * i)[2];
   }

   normalise_vector(middle_point);

   // compute the angle required for the bounding circle
   double edge_middle_point[3] = {
      (cell.coordinates_xyz)[0] + (cell.coordinates_xyz + 3)[0],
      (cell.coordinates_xyz)[1] + (cell.coordinates_xyz + 3)[1],
      (cell.coordinates_xyz)[2] + (cell.coordinates_xyz + 3)[2]};

   normalise_vector(edge_middle_point);

   struct sin_cos_angle max_angle =
      sum_angles_no_check(
        get_vector_angle_2(edge_middle_point, cell.coordinates_xyz),
        get_vector_angle_2(edge_middle_point, middle_point));

   for (unsigned i = 1; i < cell.num_corners-1; ++i) {

      double edge_middle_point[3] = {
         (cell.coordinates_xyz + 3*i)[0] + (cell.coordinates_xyz + 3*i + 3)[0],
         (cell.coordinates_xyz + 3*i)[1] + (cell.coordinates_xyz + 3*i + 3)[1],
         (cell.coordinates_xyz + 3*i)[2] + (cell.coordinates_xyz + 3*i + 3)[2]};

      normalise_vector(edge_middle_point);

      struct sin_cos_angle angle =
        sum_angles_no_check(
          get_vector_angle_2(edge_middle_point, cell.coordinates_xyz + 3 * i),
          get_vector_angle_2(edge_middle_point, middle_point));

      if (compare_angles(max_angle, angle) < 0) max_angle = angle;
   }

   bnd_circle->base_vector[0] = middle_point[0];
   bnd_circle->base_vector[1] = middle_point[1];
   bnd_circle->base_vector[2] = middle_point[2];

   bnd_circle->inc_angle = sum_angles_no_check(max_angle, SIN_COS_TOL);
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

  // compute the angle between the two base vectors
  struct sin_cos_angle base_vector_angle =
    get_vector_angle_2(circle->base_vector, dest_circle->base_vector);

  int big_sum;
  struct sin_cos_angle tmp_angle;

  // if circle is already covered by dest_circle
  big_sum = sum_angles(base_vector_angle, circle->inc_angle, &tmp_angle);
  if (big_sum || (compare_angles(dest_circle->inc_angle, tmp_angle) >= 0))
    return;

  // if dest->inc_angle already covers the whole sphere
  if (compare_angles(dest_circle->inc_angle, SIN_COS_M_PI) >= 0)
    return;

  // if circle->inc_angle already covers the whole sphere
  if (compare_angles(circle->inc_angle, SIN_COS_M_PI) >= 0) {
    *dest_circle = *circle;
    return;
  }

  struct sin_cos_angle new_inc_angle;

  big_sum = sum_angles(tmp_angle, dest_circle->inc_angle, &new_inc_angle);
  new_inc_angle = half_angle(new_inc_angle);

  // if the new inc angle is >= PI
  if (big_sum || (compare_angles(new_inc_angle, SIN_COS_M_PI) >= 0)) {
    dest_circle->inc_angle = SIN_COS_M_PI;
    return;
  }

  // compute the middle point of the merged bounding circle

  double rotation_axis[3];
  compute_norm_vector(
    rotation_axis, dest_circle->base_vector, circle->base_vector);
  struct sin_cos_angle rotation_angle;
  // we already know:
  // (base_vector_angle + circle->inc_angle) > dest_circle->inc_angle
  // therefore the following should work without problems
  sub_angles(tmp_angle, dest_circle->inc_angle, &rotation_angle);
  rotation_angle = half_angle(rotation_angle);

  double rotated_vector[3];
  rotate_vector2(
    rotation_axis, rotation_angle, dest_circle->base_vector, rotated_vector);
  normalise_vector(rotated_vector);

  dest_circle->base_vector[0] = rotated_vector[0];
  dest_circle->base_vector[1] = rotated_vector[1];
  dest_circle->base_vector[2] = rotated_vector[2];
  dest_circle->inc_angle = new_inc_angle;
}

void yac_get_grid_bounding_circle(struct grid * grid,
                                  struct bounding_circle * bnd_circle) {

   size_t num_grid_cells = (size_t)yac_get_num_grid_cells(grid);

   if (num_grid_cells == 0) {

      bnd_circle->inc_angle = SIN_COS_ZERO;
      return;
   }

   struct grid_cell cell;

   yac_init_grid_cell(&cell);

   yac_get_grid_cell2(grid, 0, &cell, bnd_circle);

   struct bounding_circle curr_bnd_circle;

   for (size_t i = 1; i < num_grid_cells; ++i) {

      yac_get_grid_cell2(grid, (unsigned)i, &cell, &curr_bnd_circle);

      merge_bounding_circles(bnd_circle, &curr_bnd_circle);
   }

   yac_free_grid_cell(&cell);
}

unsigned yac_extents_overlap(struct bounding_circle * extent_a,
                             struct bounding_circle * extent_b) {

  struct sin_cos_angle base_vector_angle =
    get_vector_angle_2(extent_a->base_vector, extent_b->base_vector);

  struct sin_cos_angle tmp_angle, inc_angle_sum;
  int big_sum =
    sum_angles(extent_a->inc_angle, extent_b->inc_angle, &tmp_angle);

  if (big_sum ||
      sum_angles(tmp_angle, SIN_COS_TOL, &inc_angle_sum))
    return 1;

  return compare_angles(base_vector_angle, inc_angle_sum) <= 0;
}

unsigned yac_point_in_bounding_circle(
  struct point point, struct bounding_circle * bnd_circle) {

  double point_vector[3];

  LLtoXYZ(point.lon, point.lat, point_vector);

  return yac_point_in_bounding_circle_vec(point_vector, bnd_circle);
}

unsigned yac_point_in_bounding_circle_vec(
  double point_vector[3], struct bounding_circle * bnd_circle) {

  return
    compare_angles(
      get_vector_angle_2(bnd_circle->base_vector, point_vector),
      bnd_circle->inc_angle) <= 0;
}

void yac_get_matching_grid_cells(struct grid * grid, struct bounding_circle extent,
                                 struct grid_cell ** matching_cells,
                                 unsigned * curr_matching_cells_array_size,
                                 unsigned ** local_ids, unsigned * curr_local_ids_array_size,
                                 unsigned * num_matching_cells, unsigned offset) {

   struct grid_search * search = yac_get_grid_search(grid);
   struct dep_list bnd_to_cell;

   yac_do_bnd_circle_search(search, &extent, 1, &bnd_to_cell);

   *num_matching_cells = yac_get_total_num_dependencies (bnd_to_cell);

   unsigned const * cells = yac_get_dependencies_of_element (bnd_to_cell, 0);

   if (matching_cells != NULL) {

      ENSURE_ARRAY_SIZE(*matching_cells, *curr_matching_cells_array_size,
                        *num_matching_cells + offset);

      for (unsigned i = offset; i < *num_matching_cells; ++i) {
        yac_init_grid_cell(*matching_cells + offset + i);
        yac_get_grid_cell(grid, cells[i], *matching_cells + offset + i);
      }
   }

   if (local_ids != NULL) {
      ENSURE_ARRAY_SIZE(*local_ids, *curr_local_ids_array_size,
                        *num_matching_cells + offset);
      memcpy(*local_ids + offset, cells,
             *num_matching_cells * sizeof(**local_ids));
   }

   yac_free_dep_list(&bnd_to_cell);
}

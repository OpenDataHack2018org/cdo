/**
 * @file check_overlap.c
 * @brief Set of functions to determine an overlap between any two cells 
 *
 * Compared to \ref check_overlap.c these routine take care of how vertex
 * points of a cell are connected, either by great circles or along the loxodrome.
 *
 * very interesting literature:
 * - http://geospatialmethods.org/spheres/GCIntersect.html
 *
 * Note: Not all functions are documented by Doxygen. See the source code
 * and \ref geometry.h for further details.
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

#include <stdlib.h>
#include <math.h>
#include "utils.h"
#include "geometry.h"

static double const tol = 1.0e-12;

// returns the direction of the given edge along the equator
// switching the points a and b of the edge give the following result:
// a->b == -(b-a)
// (exeption: if the edge is a meridian the result in always 0)
static inline int edge_direction(double * a, double * b) {

   double cross_ab_2 = a[0] * b[1] - a[1] * b[0];

   return (cross_ab_2 > - tol) - (cross_ab_2 < tol);
}

static inline void get_edge(struct grid_cell const cell, size_t edge_index,
                            size_t num_corners, struct edge * edge,
                            double ** a, double ** b,
                            struct bounding_circle * bnd_circle) {

   edge->points[0].lon = cell.coordinates_x[edge_index];
   edge->points[0].lat = cell.coordinates_y[edge_index];
   double * a_ = cell.coordinates_xyz + edge_index * 3, * b_;
   if (edge_index < num_corners - 1) {
      edge->points[1].lon = cell.coordinates_x[edge_index+1];
      edge->points[1].lat = cell.coordinates_y[edge_index+1];
      b_ = cell.coordinates_xyz + edge_index * 3 + 3;
   } else {
      edge->points[1].lon = cell.coordinates_x[0];
      edge->points[1].lat = cell.coordinates_y[0];
      b_ = cell.coordinates_xyz;
   }
   *a = a_;
   *b = b_;
   edge->edge_type = cell.edge_type[edge_index];

   bnd_circle->base_vector[0] = a_[0] + b_[0];
   bnd_circle->base_vector[1] = a_[1] + b_[1];
   bnd_circle->base_vector[2] = a_[2] + b_[2];
   normalise_vector(bnd_circle->base_vector);
   bnd_circle->inc_angle = half_angle(get_vector_angle_2(a_, b_));
}

static int point_in_cell_generic(
  double point_coords[3], struct grid_cell cell) {

   double second_point[3];

   // if the point is on the pole
   if (fabs(fabs(point_coords[2]) - 1.0) < tol) {
      second_point[0] = 1;
      second_point[1] = 0;
      second_point[2] = 0;
   } else if (point_coords[2] > 0) {
      second_point[0] = 0;
      second_point[1] = 0;
      second_point[2] = -1;
   } else {
      second_point[0] = 0;
      second_point[1] = 0;
      second_point[2] = 1;
   }

   int through_cell_corner;
   int edge_crossings;
   through_cell_corner = 0;
   edge_crossings = 0;

   // for all edges of cell
   for (size_t i = 0; i < (size_t)(cell.num_corners); ++i) {

      double * a = cell.coordinates_xyz + i * 3;
      double * b = cell.coordinates_xyz + ((i+1)%cell.num_corners) * 3;

      int ret_value;
      double p[3], q[3];

      ret_value = yac_intersect_vec(cell.edge_type[i], a, b, LON_CIRCLE,
                                    point_coords, second_point, p, q);

      // if both edges do not intersect
      if ((ret_value == -1) || (ret_value == 0)) continue;

      // if p is not the intersection point of both edges
      if ((ret_value & ((1 << 0) + (1 << 2))) != (1 << 0) + (1 << 2)) {

         // if q is the intersection point of both edges
         if ((ret_value & ((1 << 1) + (1 << 3))) == (1 << 1) + (1 << 3))
            p[0] = q[0], p[1] = q[1], p[2] = q[2];
         else
            continue;
      }

      // if the intersection is the point itself
      if (points_are_identically(point_coords, p))
         return 1;

      //----------
      // in case the point edge goes through a endpoint of a cell edge
      // there are three special cases that need to be taken care of:
      // 1: the cell is concave an the point edge only touched an inner corner of the cell -> does not switch in/outside
      // 2: the point edge only touched an outer corner of the cell -> does not switch in/outside
      // 3: the point edge entered/left the cell through a corner -> switch in/outside
      //---------
      if ((points_are_identically(p, a)) ||
          (points_are_identically(p, b))) {

         // check the direction of the cell edge (edge can have a positive or
         // negative direction in terms of longitudes)
         // each cell endpoint is checked twice (for each adjacent edge)
         // If the direction for two adjacent edges sharing an endpoint that
         // is crossed by the point edge, then the related endpoint
         // intersection can be classified as case 3.
         // case 1 and 2 result in no change to through_cell_corner
         // case 3 result in an in/decrease of through_cell_corner by 2

         through_cell_corner += edge_direction(a, b);
      } else {

         //standard case -> we crossed an edge
         edge_crossings++;
      }
   } // (j = 0; j < cell.num_corners; ++j)

   edge_crossings += through_cell_corner / 2;

   // if we have an odd number of edge crossings -> point was inside cell
   // or the point was directly on an longitude edge
   return ( (edge_crossings & 1) || (through_cell_corner & 1) );
}

static int point_in_cell_unstruct_triangle(
  double point_coords[3], struct grid_cell cell) {

  double * a = cell.coordinates_xyz + 0;
  double * b = cell.coordinates_xyz + 3;
  double * c = cell.coordinates_xyz + 6;

  double cross_edges[3][3];

  crossproduct_ld(a, b, &(cross_edges[0][0]));
  crossproduct_ld(b, c, &(cross_edges[1][0]));
  crossproduct_ld(c, a, &(cross_edges[2][0]));

  double sq_sin_edge[3] = {
    cross_edges[0][0] * cross_edges[0][0] +
    cross_edges[0][1] * cross_edges[0][1] +
    cross_edges[0][2] * cross_edges[0][2],
    cross_edges[1][0] * cross_edges[1][0] +
    cross_edges[1][1] * cross_edges[1][1] +
    cross_edges[1][2] * cross_edges[1][2],
    cross_edges[2][0] * cross_edges[2][0] +
    cross_edges[2][1] * cross_edges[2][1] +
    cross_edges[2][2] * cross_edges[2][2]};

  normalise_vector(&(cross_edges[0][0]));
  normalise_vector(&(cross_edges[1][0]));
  normalise_vector(&(cross_edges[2][0]));

  double sin_point_edge[3] = {
    cross_edges[0][0] * point_coords[0] +
    cross_edges[0][1] * point_coords[1] +
    cross_edges[0][2] * point_coords[2],
    cross_edges[1][0] * point_coords[0] +
    cross_edges[1][1] * point_coords[1] +
    cross_edges[1][2] * point_coords[2],
    cross_edges[2][0] * point_coords[0] +
    cross_edges[2][1] * point_coords[1] +
    cross_edges[2][2] * point_coords[2]};

  // if all edges of the triangle have a length of at least yac_angle_tol
  if ((sq_sin_edge[0] > (yac_angle_tol * yac_angle_tol)) &&
      (sq_sin_edge[1] > (yac_angle_tol * yac_angle_tol)) &&
      (sq_sin_edge[2] > (yac_angle_tol * yac_angle_tol))) {

    int cell_corner_order =
      cross_edges[0][0] * c[0] + cross_edges[0][1] * c[1] + cross_edges[0][2] * c[2] >= 0;

    if (cell_corner_order)
      return (sin_point_edge[0] > - yac_angle_tol) &&
             (sin_point_edge[1] > - yac_angle_tol) &&
             (sin_point_edge[2] > - yac_angle_tol);
    else
      return (sin_point_edge[0] < yac_angle_tol) &&
             (sin_point_edge[1] < yac_angle_tol) &&
             (sin_point_edge[2] < yac_angle_tol);

  // if the length of the first edge is shorter than yac_angle_tol
  } else if (sq_sin_edge[0] <= (yac_angle_tol * yac_angle_tol)) {

    // if the length of the second edge is also shorter than yac_angle_tol
    // -> triangle is a point
    if (sq_sin_edge[1] <= (yac_angle_tol * yac_angle_tol))
      return
        compare_angles(get_vector_angle_2(a, point_coords), SIN_COS_TOL) <= 0;
    // -> triangle is an edge
    else {
      struct sin_cos_angle edge_angle =
        sum_angles_no_check(get_vector_angle_2(a, c), SIN_COS_TOL);
      struct sin_cos_angle ap_angle = get_vector_angle_2(a, point_coords);
      struct sin_cos_angle cp_angle = get_vector_angle_2(c, point_coords);
      // point is on the plane of the edge and between b and c
      return (fabs(sin_point_edge[1]) <= yac_angle_tol) &&
             (compare_angles(ap_angle, edge_angle) <= 0) &&
             (compare_angles(cp_angle, edge_angle) <= 0);
    }

  // -> triangle is an edge
  } else {
    struct sin_cos_angle edge_angle =
      sum_angles_no_check(get_vector_angle_2(a, b), SIN_COS_TOL);
    struct sin_cos_angle ap_angle = get_vector_angle_2(a, point_coords);
    struct sin_cos_angle bp_angle = get_vector_angle_2(b, point_coords);
    // point is on the plane of the edge and between b and c
    return (fabs(sin_point_edge[0]) <= yac_angle_tol) &&
           (compare_angles(ap_angle, edge_angle) <= 0) &&
           (compare_angles(bp_angle, edge_angle) <= 0);
  }
}

static int any_point_in_cell_unstruct_triangle(
  double * point_coords, size_t num_points, struct grid_cell cell) {

  double * a = cell.coordinates_xyz + 0;
  double * b = cell.coordinates_xyz + 3;
  double * c = cell.coordinates_xyz + 6;

  double cross_edges[3][3];

  crossproduct_ld(a, b, &(cross_edges[0][0]));
  crossproduct_ld(b, c, &(cross_edges[1][0]));
  crossproduct_ld(c, a, &(cross_edges[2][0]));

  double sq_sin_edge[3] = {
    cross_edges[0][0] * cross_edges[0][0] +
    cross_edges[0][1] * cross_edges[0][1] +
    cross_edges[0][2] * cross_edges[0][2],
    cross_edges[1][0] * cross_edges[1][0] +
    cross_edges[1][1] * cross_edges[1][1] +
    cross_edges[1][2] * cross_edges[1][2],
    cross_edges[2][0] * cross_edges[2][0] +
    cross_edges[2][1] * cross_edges[2][1] +
    cross_edges[2][2] * cross_edges[2][2]};

  normalise_vector(&(cross_edges[0][0]));
  normalise_vector(&(cross_edges[1][0]));
  normalise_vector(&(cross_edges[2][0]));

  // if all edges of the triangle have a length of at least yac_angle_tol
  if ((sq_sin_edge[0] > (yac_angle_tol * yac_angle_tol)) &&
      (sq_sin_edge[1] > (yac_angle_tol * yac_angle_tol)) &&
      (sq_sin_edge[2] > (yac_angle_tol * yac_angle_tol))) {

    int cell_corner_order =
      cross_edges[0][0] * c[0] + cross_edges[0][1] * c[1] + cross_edges[0][2] * c[2] >= 0;

    if (cell_corner_order) {
      for (size_t i = 0; i < num_points; ++i) {
        double sin_point_edge[3] = {
          cross_edges[0][0] * point_coords[0+3*i] +
          cross_edges[0][1] * point_coords[1+3*i] +
          cross_edges[0][2] * point_coords[2+3*i],
          cross_edges[1][0] * point_coords[0+3*i] +
          cross_edges[1][1] * point_coords[1+3*i] +
          cross_edges[1][2] * point_coords[2+3*i],
          cross_edges[2][0] * point_coords[0+3*i] +
          cross_edges[2][1] * point_coords[1+3*i] +
          cross_edges[2][2] * point_coords[2+3*i]};
        if ((sin_point_edge[0] > - yac_angle_tol) &&
            (sin_point_edge[1] > - yac_angle_tol) &&
            (sin_point_edge[2] > - yac_angle_tol)) return 1;
      }
    } else {
      for (size_t i = 0; i < num_points; ++i) {
        double sin_point_edge[3] = {
          cross_edges[0][0] * point_coords[0+3*i] +
          cross_edges[0][1] * point_coords[1+3*i] +
          cross_edges[0][2] * point_coords[2+3*i],
          cross_edges[1][0] * point_coords[0+3*i] +
          cross_edges[1][1] * point_coords[1+3*i] +
          cross_edges[1][2] * point_coords[2+3*i],
          cross_edges[2][0] * point_coords[0+3*i] +
          cross_edges[2][1] * point_coords[1+3*i] +
          cross_edges[2][2] * point_coords[2+3*i]};
        if ((sin_point_edge[0] < yac_angle_tol) &&
            (sin_point_edge[1] < yac_angle_tol) &&
            (sin_point_edge[2] < yac_angle_tol)) return 1;
      }
    }

  // if the length of the first edge is shorter than yac_angle_tol
  } else if (sq_sin_edge[0] <= (yac_angle_tol * yac_angle_tol)) {

    // if the length of the second edge is also shorter than yac_angle_tol
    // -> triangle is a point
    if (sq_sin_edge[1] <= (yac_angle_tol * yac_angle_tol))
      return
        compare_angles(get_vector_angle_2(a, point_coords), SIN_COS_TOL) <= 0;
    // -> triangle is an edge
    else {
      struct sin_cos_angle edge_angle =
        sum_angles_no_check(get_vector_angle_2(a, c), SIN_COS_TOL);

      for (size_t i = 0; i < num_points; ++i) {
        struct sin_cos_angle ap_angle = get_vector_angle_2(a, point_coords+3*i);
        struct sin_cos_angle cp_angle = get_vector_angle_2(c, point_coords+3*i);
        double sin_point_edge = cross_edges[1][0] * point_coords[0+3*i] +
                                cross_edges[1][1] * point_coords[1+3*i] +
                                cross_edges[1][2] * point_coords[2+3*i];
        // point is on the plane of the edge and between b and c
        if ((fabs(sin_point_edge) <= yac_angle_tol) &&
            (compare_angles(ap_angle, edge_angle) <= 0) &&
            (compare_angles(cp_angle, edge_angle) <= 0)) return 1;
      }
    }

  // -> triangle is an edge
  } else {
    struct sin_cos_angle edge_angle =
      sum_angles_no_check(get_vector_angle_2(a, b), SIN_COS_TOL);

    for (size_t i = 0; i < num_points; ++i) {
      struct sin_cos_angle ap_angle = get_vector_angle_2(a, point_coords+3*i);
      struct sin_cos_angle bp_angle = get_vector_angle_2(b, point_coords+3*i);
      double sin_point_edge = cross_edges[0][0] * point_coords[0+3*i] +
                              cross_edges[0][1] * point_coords[1+3*i] +
                              cross_edges[0][2] * point_coords[2+3*i];
      // point is on the plane of the edge and between b and c
      if  ((fabs(sin_point_edge) <= yac_angle_tol) &&
           (compare_angles(ap_angle, edge_angle) <= 0) &&
           (compare_angles(bp_angle, edge_angle) <= 0)) return 1;
    }
  }

  return 0;
}

static inline double dotproduct(double a[], double b[]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static int point_in_cell_unstruct_quad(
  double point_coords[3], struct grid_cell cell) {

  double * a = cell.coordinates_xyz + 0;
  double * b = cell.coordinates_xyz + 3;
  double * c = cell.coordinates_xyz + 6;
  double * d = cell.coordinates_xyz + 9;

  double cross[4][3];

  crossproduct_ld(a, b, &(cross[0][0]));
  crossproduct_ld(b, c, &(cross[1][0]));
  crossproduct_ld(c, d, &(cross[2][0]));
  crossproduct_ld(d, a, &(cross[3][0]));

  double sq_sin_edge[4] = {
    cross[0][0]*cross[0][0]+cross[0][1]*cross[0][1]+cross[0][2]*cross[0][2],
    cross[1][0]*cross[1][0]+cross[1][1]*cross[1][1]+cross[1][2]*cross[1][2],
    cross[2][0]*cross[2][0]+cross[2][1]*cross[2][1]+cross[2][2]*cross[2][2],
    cross[3][0]*cross[3][0]+cross[3][1]*cross[3][1]+cross[3][2]*cross[3][2]};

  normalise_vector(&(cross[0][0]));
  normalise_vector(&(cross[1][0]));
  normalise_vector(&(cross[2][0]));
  normalise_vector(&(cross[3][0]));

  // if all edges of the quad have a length of at least yac_angle_tol
  if ((sq_sin_edge[0] > (yac_angle_tol * yac_angle_tol)) &&
      (sq_sin_edge[1] > (yac_angle_tol * yac_angle_tol)) &&
      (sq_sin_edge[2] > (yac_angle_tol * yac_angle_tol)) &&
      (sq_sin_edge[3] > (yac_angle_tol * yac_angle_tol))) {

    double tmp = copysign(1.0, dotproduct(c, cross[0])) +
                 copysign(1.0, dotproduct(d, cross[0])) +
                 copysign(1.0, dotproduct(a, cross[1])) +
                 copysign(1.0, dotproduct(d, cross[1])) +
                 copysign(1.0, dotproduct(a, cross[2])) +
                 copysign(1.0, dotproduct(b, cross[2])) +
                 copysign(1.0, dotproduct(b, cross[3])) +
                 copysign(1.0, dotproduct(c, cross[3]));
    // check whether the quad is convex
    if (tmp == 8.0)
      return
        copysign(1.0, dotproduct(point_coords, cross[0]) + yac_angle_tol) +
        copysign(1.0, dotproduct(point_coords, cross[1]) + yac_angle_tol) +
        copysign(1.0, dotproduct(point_coords, cross[2]) + yac_angle_tol) +
        copysign(1.0, dotproduct(point_coords, cross[3]) + yac_angle_tol) ==
        4.0;
    else if (tmp == -8.0)
      return
        copysign(1.0, dotproduct(point_coords, cross[0]) - yac_angle_tol) +
        copysign(1.0, dotproduct(point_coords, cross[1]) - yac_angle_tol) +
        copysign(1.0, dotproduct(point_coords, cross[2]) - yac_angle_tol) +
        copysign(1.0, dotproduct(point_coords, cross[3]) - yac_angle_tol) ==
        -4.0;
  }

  return point_in_cell_generic(point_coords, cell);
}

static int any_point_in_cell_unstruct_quad(
  double * point_coords, size_t num_points, struct grid_cell cell) {

  double * a = cell.coordinates_xyz + 0;
  double * b = cell.coordinates_xyz + 3;
  double * c = cell.coordinates_xyz + 6;
  double * d = cell.coordinates_xyz + 9;

  double cross[4][3];

  crossproduct_ld(a, b, &(cross[0][0]));
  crossproduct_ld(b, c, &(cross[1][0]));
  crossproduct_ld(c, d, &(cross[2][0]));
  crossproduct_ld(d, a, &(cross[3][0]));

  double sq_sin_edge[4] = {
    cross[0][0]*cross[0][0]+cross[0][1]*cross[0][1]+cross[0][2]*cross[0][2],
    cross[1][0]*cross[1][0]+cross[1][1]*cross[1][1]+cross[1][2]*cross[1][2],
    cross[2][0]*cross[2][0]+cross[2][1]*cross[2][1]+cross[2][2]*cross[2][2],
    cross[3][0]*cross[3][0]+cross[3][1]*cross[3][1]+cross[3][2]*cross[3][2]};

  normalise_vector(&(cross[0][0]));
  normalise_vector(&(cross[1][0]));
  normalise_vector(&(cross[2][0]));
  normalise_vector(&(cross[3][0]));

  // if all edges of the quad have a length of at least yac_angle_tol
  if ((sq_sin_edge[0] > (yac_angle_tol * yac_angle_tol)) &&
      (sq_sin_edge[1] > (yac_angle_tol * yac_angle_tol)) &&
      (sq_sin_edge[2] > (yac_angle_tol * yac_angle_tol)) &&
      (sq_sin_edge[3] > (yac_angle_tol * yac_angle_tol))) {

    double tmp = copysign(1.0, dotproduct(c, cross[0])) +
                 copysign(1.0, dotproduct(d, cross[0])) +
                 copysign(1.0, dotproduct(a, cross[1])) +
                 copysign(1.0, dotproduct(d, cross[1])) +
                 copysign(1.0, dotproduct(a, cross[2])) +
                 copysign(1.0, dotproduct(b, cross[2])) +
                 copysign(1.0, dotproduct(b, cross[3])) +
                 copysign(1.0, dotproduct(c, cross[3]));
    // check whether the quad is convex
    if (tmp == 8.0) {
      for (size_t i = 0; i < num_points; ++i) {
        if (copysign(1.0, dotproduct(point_coords+3*i, cross[0]) + yac_angle_tol) +
            copysign(1.0, dotproduct(point_coords+3*i, cross[1]) + yac_angle_tol) +
            copysign(1.0, dotproduct(point_coords+3*i, cross[2]) + yac_angle_tol) +
            copysign(1.0, dotproduct(point_coords+3*i, cross[3]) + yac_angle_tol) ==
            4.0) return 1;
      }
    } else if (tmp == -8.0) {
      for (size_t i = 0; i < num_points; ++i) {
        if (copysign(1.0, dotproduct(point_coords+3*i, cross[0]) - yac_angle_tol) +
            copysign(1.0, dotproduct(point_coords+3*i, cross[1]) - yac_angle_tol) +
            copysign(1.0, dotproduct(point_coords+3*i, cross[2]) - yac_angle_tol) +
            copysign(1.0, dotproduct(point_coords+3*i, cross[3]) - yac_angle_tol) ==
            -4.0) return 1;
      }
    }
  } else {
    for (size_t i = 0; i < num_points; ++i)
      return point_in_cell_generic(point_coords+3*i, cell);
  }
  return 0;
}

// checks whether a point is within a cell
int yac_point_in_cell (double point_coords[3], struct grid_cell cell) {

   struct bounding_circle bnd_circle;
   yac_get_cell_bounding_circle(cell, &bnd_circle);

   return yac_point_in_cell2(point_coords, cell, bnd_circle);
}

// checks whether a point is within a cell
int yac_point_in_cell2 (double point_coords[3], struct grid_cell cell,
                        struct bounding_circle bnd_circle) {

   // check whether the point is within the bounding circle of the cell;

   if (!yac_point_in_bounding_circle_vec(point_coords, &bnd_circle))
      return 1 == 0;

   if ((cell.num_corners == 3) &&
       (cell.edge_type[0] == GREAT_CIRCLE) &&
       (cell.edge_type[1] == GREAT_CIRCLE) &&
       (cell.edge_type[2] == GREAT_CIRCLE))
      return point_in_cell_unstruct_triangle(point_coords, cell);
   else if ((cell.num_corners == 4) &&
       (cell.edge_type[0] == GREAT_CIRCLE) &&
       (cell.edge_type[1] == GREAT_CIRCLE) &&
       (cell.edge_type[2] == GREAT_CIRCLE) &&
       (cell.edge_type[3] == GREAT_CIRCLE))
      return point_in_cell_unstruct_quad(point_coords, cell);
   else
      return point_in_cell_generic(point_coords, cell);
}

static int any_point_in_cell2(
  double * point_coords, size_t num_points, struct grid_cell cell,
  struct bounding_circle bnd_circle) {

  double point_coords_[num_points][3];

  size_t count = 0;
  for (size_t i = 0; i < num_points; ++i) {
    if (yac_point_in_bounding_circle_vec(point_coords + 3 * i, &bnd_circle)) {
      point_coords_[count][0] = point_coords[0+3*i];
      point_coords_[count][1] = point_coords[1+3*i];
      point_coords_[count][2] = point_coords[2+3*i];
      ++count;
    }
  }

  if (count == 0) return 0;

  if ((cell.num_corners == 3) &&
      (cell.edge_type[0] == GREAT_CIRCLE) &&
      (cell.edge_type[1] == GREAT_CIRCLE) &&
      (cell.edge_type[2] == GREAT_CIRCLE))
    return
      any_point_in_cell_unstruct_triangle(&(point_coords_[0][0]), count, cell);
  else if ((cell.num_corners == 4) &&
           (cell.edge_type[0] == GREAT_CIRCLE) &&
           (cell.edge_type[1] == GREAT_CIRCLE) &&
           (cell.edge_type[2] == GREAT_CIRCLE) &&
           (cell.edge_type[3] == GREAT_CIRCLE))
    return any_point_in_cell_unstruct_quad(&(point_coords_[0][0]), count, cell);
  else {
    for (size_t i = 0; i < num_points; ++i)
     return point_in_cell_generic(&(point_coords_[i][0]), cell);
  }
  return 0;
}

// checks whether any point of cell b is inside of cell a
static inline int inside (struct grid_cell const cell_a,
                          struct bounding_circle circle_a,
                          struct grid_cell const cell_b,
                          struct bounding_circle circle_b) {

  return
    any_point_in_cell2(
      cell_b.coordinates_xyz, cell_b.num_corners, cell_a, circle_a);
}

static inline int exact (struct grid_cell const cell_a, 
                         struct grid_cell const cell_b) {

   size_t n = 0;

   for (size_t i = 0; i < (size_t)(cell_a.num_corners); ++i) {
      for (size_t j = 0; j < (size_t)(cell_b.num_corners); ++j) {

         if ((fabs(cell_a.coordinates_xyz[0+3*i] -
                   cell_b.coordinates_xyz[0+3*j]) < tol) &&
             (fabs(cell_a.coordinates_xyz[1+3*i] -
                   cell_b.coordinates_xyz[1+3*j]) < tol) &&
             (fabs(cell_a.coordinates_xyz[2+3*i] -
                   cell_b.coordinates_xyz[2+3*j]) < tol)) {

            if (++n > 2) return 1;
            break;
         }
      }
   }
   return 0;
}

/** \brief simple test to identify whether cells are given on a regular longitude-latitude grid 
 *
 **/
static inline int is_regular_cell(struct grid_cell const cell) {

   return (cell.num_corners == 4) &&
          (cell.edge_type[0] == LAT_CIRCLE) &&
          (cell.edge_type[1] == LON_CIRCLE) &&
          (cell.edge_type[2] == LAT_CIRCLE) &&
          (cell.edge_type[3] == LON_CIRCLE);
}

/** \brief checks overlap of two cells defined on a regular longitude-latitude grid 
 *
 **/
static int regular_cells_intersect(struct grid_cell const cell_a,
                                   struct grid_cell const cell_b) {

   double lon_diff_a, lon_diff_b;

   lon_diff_a = fabs(get_angle(cell_a.coordinates_x[0], cell_a.coordinates_x[1]));
   lon_diff_b = fabs(get_angle(cell_b.coordinates_x[0], cell_b.coordinates_x[1]));

   double lon_diff_ab[4];

   lon_diff_ab[0] = fabs(get_angle(cell_a.coordinates_x[0], cell_b.coordinates_x[0]));
   lon_diff_ab[1] = fabs(get_angle(cell_a.coordinates_x[0], cell_b.coordinates_x[1]));
   lon_diff_ab[2] = fabs(get_angle(cell_a.coordinates_x[1], cell_b.coordinates_x[0]));
   lon_diff_ab[3] = fabs(get_angle(cell_a.coordinates_x[1], cell_b.coordinates_x[1]));

   double min_lat_a, max_lat_a, min_lat_b, max_lat_b;

   if (cell_a.coordinates_y[0] > cell_a.coordinates_y[2]) {
      max_lat_a = cell_a.coordinates_y[0];
      min_lat_a = cell_a.coordinates_y[2];
   } else {
      max_lat_a = cell_a.coordinates_y[2];
      min_lat_a = cell_a.coordinates_y[0];
   }

   if (cell_b.coordinates_y[0] > cell_b.coordinates_y[2]) {
      max_lat_b = cell_b.coordinates_y[0];
      min_lat_b = cell_b.coordinates_y[2];
   } else {
      max_lat_b = cell_b.coordinates_y[2];
      min_lat_b = cell_b.coordinates_y[0];
   }

   return ((min_lat_a < max_lat_b && min_lat_b < max_lat_a) &&
           ((lon_diff_a + tol >= lon_diff_ab[0] + lon_diff_ab[2]) ||
            (lon_diff_a + tol >= lon_diff_ab[1] + lon_diff_ab[3]) ||
            (lon_diff_b + tol >= lon_diff_ab[0] + lon_diff_ab[1]) ||
            (lon_diff_b + tol >= lon_diff_ab[2] + lon_diff_ab[3])));
}

/** \brief checks whether two cells overlap
 *
 * the cells can be concave but the extents of the cells
 * should be below PI/2 in lon and lat direction
 *
 **/
int yac_check_overlap_cells (struct grid_cell const cell_a,
                             struct grid_cell const cell_b) {

   struct bounding_circle circle_a, circle_b;

   // get the bounding circles of both cells
   yac_get_cell_bounding_circle(cell_a, &circle_a);
   yac_get_cell_bounding_circle(cell_b, &circle_b);

   return yac_check_overlap_cells2(cell_a, circle_a, cell_b, circle_b);
}

/** \brief checks whether two cells overlap
 *
 * the cells can be concave but the extents of the cells
 * should be below PI/2 in lon and lat direction
 *
 **/
int yac_check_overlap_cells2 (struct grid_cell const cell_a,
                              struct bounding_circle circle_a,
                              struct grid_cell const cell_b,
                              struct bounding_circle circle_b) {

   // check whether the cells have the potential to overlap
   if (!yac_extents_overlap(&circle_a, &circle_b)) return 0;

   // Test for exact matches
   if (exact (cell_a, cell_b)) return 1;

   // check for the special case of two regular cells
   if (is_regular_cell(cell_a) && is_regular_cell(cell_b)) {

      return regular_cells_intersect(cell_a, cell_b);

   } else {

      // Test whether one point is inside the cell
      if (inside (cell_a, circle_a, cell_b, circle_b)) return 1;
      if (inside (cell_b, circle_b, cell_a, circle_a)) return 1;

      struct edge edge_a, edge_b;
      double * a, * b, * c, * d;
      struct bounding_circle edge_bnd_circle_ab;
      struct bounding_circle edge_bnd_circle_cd;
      size_t num_corners_a = (size_t)cell_a.num_corners;
      size_t num_corners_b = (size_t)cell_b.num_corners;

      // Test intersection of two edges
      for (size_t i = 0; i < num_corners_a; ++i) {

         get_edge(
            cell_a, i, num_corners_a, &edge_a, &a, &b, &edge_bnd_circle_ab);

         for (size_t j = 0; j < num_corners_b; ++j) {

            get_edge(
               cell_b, j, num_corners_b, &edge_b, &c, &d, &edge_bnd_circle_cd);

            if (yac_extents_overlap(&edge_bnd_circle_ab, &edge_bnd_circle_cd) &&
                yac_do_intersect (edge_a, a, b, edge_b, c, d)) return 1;
         }
      }
   }

   return 0;
}

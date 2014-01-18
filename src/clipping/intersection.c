/**
 * @file intersection.c
 * @brief Set of functions to determine the intersection between two edges
 *
 * very interesting literature:
 * - http://geospatialmethods.org/spheres/GCIntersect.html
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

#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "geometry.h"

static double const tol = 1.0e-12;

static void crossproduct (double a[], double b[], double cross[]) {

/* crossproduct in cartesian coordinates */

   cross[0] = a[1] * b[2] - a[2] * b[1];
   cross[1] = a[2] * b[0] - a[0] * b[2];
   cross[2] = a[0] * b[1] - a[1] * b[0];
}

static int vector_is_between (double a[], double b[], double p[], double e_ab[]) {

/* determines whether p is between a and b
   (a, b, p are in the same plane AB)
   e_ab is the crossproduct of a and b */

   double cross_ap, cross_pb; // we only need one element of the cross product
   int needed_index;

   needed_index = 0;
   if (fabs(e_ab[1]) > fabs(e_ab[0])) needed_index |= 1;
   if (fabs(e_ab[2]) > fabs(e_ab[0])) needed_index |= 2;
   if (fabs(e_ab[2]) > fabs(e_ab[1])) needed_index |= 4;

   switch (needed_index) {
      case (0): // index 0 is biggest value in e_ab
      case (4):
         cross_ap = a[1] * p[2] - a[2] * p[1];
         cross_pb = p[1] * b[2] - p[2] * b[1];

         return (e_ab[0] > - tol && cross_ap > - tol && cross_pb > - tol) ||
                (e_ab[0] < + tol && cross_ap < + tol && cross_pb < + tol);

      case (1): // index 1 is biggest value in e_ab
      case (3):
         cross_ap = a[2] * p[0] - a[0] * p[2];
         cross_pb = p[2] * b[0] - p[0] * b[2];

         return (e_ab[1] > - tol && cross_ap > - tol && cross_pb > - tol) ||
                (e_ab[1] < + tol && cross_ap < + tol && cross_pb < + tol);

      case (6): // index 2 is biggest value in e_ab
      case (7):
         cross_ap = a[0] * p[1] - a[1] * p[0];
         cross_pb = p[0] * b[1] - p[1] * b[0];


         return (e_ab[2] > - tol && cross_ap > - tol && cross_pb > - tol) ||
                (e_ab[2] < + tol && cross_ap < + tol && cross_pb < + tol);

      default:
         abort_message("internal error", __FILE__, __LINE__);
         // this function should never reach this point...
         return -1;
   };
}

/** \brief compute the intersection points of two great circles
  *
  * if p and q are != NULL they contain the intersection points
  *
  * the return value is :
  *    -  0 if the intersection points are neither between (a and b) or (c and d)
  *    - 1st bit will be set if p is between a and b
  *    - 2nd bit will be set if q is between a and b
  *    - 3rd bit will be set if p is between c and d
  *    - 4th bit will be set if q is between c and d
  *    - 5th bit will be set if both great circles are identically
  *
  * based on
  * - http://www.geoclub.de/viewtopic.php?f=54&t=29689
  **/

 int gcxgc (struct edge edge_a, struct edge edge_b,
            struct point * p, struct point * q) {

   double a_[3], b_[3], c_[3], d_[3], e_ab[3], e_cd[3], n;

   // TODO check for special case: meridian

   LLtoXYZ( edge_a.points[0].lon, edge_a.points[0].lat, a_);
   LLtoXYZ( edge_a.points[1].lon, edge_a.points[1].lat, b_);
   LLtoXYZ( edge_b.points[0].lon, edge_b.points[0].lat, c_);
   LLtoXYZ( edge_b.points[1].lon, edge_b.points[1].lat, d_);

   // compute unit vector of ab plane
   crossproduct(a_, b_, e_ab);
   n = 1.0 / sqrt(e_ab[0] * e_ab[0] + e_ab[1] * e_ab[1] + e_ab[2] * e_ab[2]);
   e_ab[0] *= n;
   e_ab[1] *= n;
   e_ab[2] *= n;

   // compute unit vector of cd plane
   crossproduct(c_, d_, e_cd);
   n = 1.0 / sqrt(e_cd[0] * e_cd[0] + e_cd[1] * e_cd[1] + e_cd[2] * e_cd[2]);
   e_cd[0] *= n;
   e_cd[1] *= n;
   e_cd[2] *= n;

   double f1, f2;

   // compute cos between e and c_/d_ times length of e
   f1 = e_ab[0] * c_[0] + e_ab[1] * c_[1] + e_ab[2] * c_[2];
   f2 = e_ab[0] * d_[0] + e_ab[1] * d_[1] + e_ab[2] * d_[2];

   // if both great circles are nearly identically
   if ((fabs(f1) < tol) && (fabs(f2) < tol)) {

      int ret_value = 1 << 4;
   
      // if c is between a and b
      if (vector_is_between(a_, b_, c_, e_ab) > 0) {

         if (p != NULL)
            *p = edge_b.points[0];

         if (q != NULL) {
            q->lon = edge_b.points[0].lon + M_PI;
            q->lat = - edge_b.points[0].lat;
         }

         ret_value |= (1 << 0) | (1 << 2);

      // if d is between a and b
      } else if (vector_is_between(a_, b_, d_, e_ab) > 0) {

         if (p != NULL)
            *p = edge_b.points[1];

         if (q != NULL) {
            q->lon = edge_b.points[1].lon + M_PI;
            q->lat = - edge_b.points[1].lat;
         }

         ret_value |= (1 << 0) | (1 << 2);

      // if a is between c and d
      } else if (vector_is_between(c_, d_, a_, e_cd)) {

         if (p != NULL)
            *p = edge_a.points[0];

         if (q != NULL) {
            q->lon = edge_a.points[0].lon + M_PI;
            q->lat = - edge_a.points[0].lat;
         }

         ret_value |= (1 << 0) | (1 << 2);

      } else {

         if (p != NULL)
            *p = edge_a.points[0];

         if (q != NULL) {
            q->lon = edge_a.points[0].lon + M_PI;
            q->lat = - edge_a.points[0].lat;
         }

         ret_value |= 1 << 0;
      }

      return ret_value;
  }

   double g[3];

   // compute line intersection of both planes defined by the two great circles
   if (fabs(f1) > fabs(f2)) {

      n = - (f2 / f1);

      g[0] = n * c_[0] + d_[0];
      g[1] = n * c_[1] + d_[1];
      g[2] = n * c_[2] + d_[2];

   } else {

      n = - (f1 / f2);

      g[0] = c_[0] + n * d_[0];
      g[1] = c_[1] + n * d_[1];
      g[2] = c_[2] + n * d_[2];
   }

   // normalise g
   n = 1.0 / sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);

   g[0]=g[0]*n;
   g[1]=g[1]*n;
   g[2]=g[2]*n;

   // determine p and q
   double p_[3], q_[3];

   p_[0]= g[0];
   p_[1]= g[1];
   p_[2]= g[2];

   q_[0]=-p_[0];
   q_[1]=-p_[1];
   q_[2]=-p_[2];

   // set p and q
   if (p != 0) XYZtoLL(p_, &p->lon, &p->lat);
   if (q != 0) XYZtoLL(q_, &q->lon, &q->lat);

   int result;

   result = 0;
   if (vector_is_between(a_, b_, p_, e_ab)) result |= 1 << 0;
   if (vector_is_between(a_, b_, q_, e_ab)) result |= 1 << 1;
   if (vector_is_between(c_, d_, p_, e_cd)) result |= 1 << 2;
   if (vector_is_between(c_, d_, q_, e_cd)) result |= 1 << 3;

   return result;
}

/** \brief compute the intersection points of two great circles
  *
  * if p and q are != NULL they contain the intersection points
  *
  * the return value is :
  *    -  0 if the intersection points are neither between (a and b) or (c and d)
  *    - 1st bit will be set if p is between a and b
  *    - 2nd bit will be set if q is between a and b
  *    - 3rd bit will be set if p is between c and d
  *    - 4th bit will be set if q is between c and d
  *    - 5th bit will be set if both great circles are identically
  *
  * based on
  * - http://www.geoclub.de/viewtopic.php?f=54&t=29689
  **/
 int gcxgc_vec (double a[3], double b[3], double c[3], double d[3],
                double p[3], double q[3]) {

   double e_ab[3], e_cd[3], n;

   // compute unit vector of ab plane
   crossproduct(a, b, e_ab);
   n = 1.0 / sqrt(e_ab[0] * e_ab[0] + e_ab[1] * e_ab[1] + e_ab[2] * e_ab[2]);
   e_ab[0] *= n;
   e_ab[1] *= n;
   e_ab[2] *= n;

   // compute unit vector of cd plane
   crossproduct(c, d, e_cd);
   n = 1.0 / sqrt(e_cd[0] * e_cd[0] + e_cd[1] * e_cd[1] + e_cd[2] * e_cd[2]);
   e_cd[0] *= n;
   e_cd[1] *= n;
   e_cd[2] *= n;

   double f1, f2;

   // compute cos between e and c_/d_ times length of e
   f1 = e_ab[0] * c[0] + e_ab[1] * c[1] + e_ab[2] * c[2];
   f2 = e_ab[0] * d[0] + e_ab[1] * d[1] + e_ab[2] * d[2];

   // if both great circles are nearly identically
   if ((fabs(f1) < tol) && (fabs(f2) < tol)) {

      int ret_value = 1 << 4;
   
      // if c is between a and b
      if (vector_is_between(a, b, c, e_ab) > 0) {

         if (p != NULL) {
            p[0] = c[0];
            p[1] = c[1];
            p[2] = c[2];
         }

         if (q != NULL) {
            q[0] = -c[0];
            q[1] = -c[1];
            q[2] = -c[2];
         }

         ret_value |= (1 << 0) | (1 << 2);

      // if d is between a and b
      } else if (vector_is_between(a, b, d, e_ab) > 0) {

         if (p != NULL) {
            p[0] = d[0];
            p[1] = d[1];
            p[2] = d[2];
         }

         if (q != NULL) {
            q[0] = -d[0];
            q[1] = -d[1];
            q[2] = -d[2];
         }

         ret_value |= (1 << 0) | (1 << 2);

      // if a is between c and d
      } else if (vector_is_between(c, d, a, e_cd)) {

         if (p != NULL) {
            p[0] = a[0];
            p[1] = a[1];
            p[2] = a[2];
         }

         if (q != NULL) {
            q[0] = -a[0];
            q[1] = -a[1];
            q[2] = -a[2];
         }

         ret_value |= (1 << 0) | (1 << 2);

      } else {

         if (p != NULL) {
            p[0] = a[0];
            p[1] = a[1];
            p[2] = a[2];
         }

         if (q != NULL) {
            q[0] = -a[0];
            q[1] = -a[1];
            q[2] = -a[2];
         }

         ret_value |= 1 << 0;
      }

      return ret_value;
  }

   double g[3];

   // compute line intersection of both planes defined by the two great circles
   if (fabs(f1) > fabs(f2)) {

      n = - (f2 / f1);

      g[0] = n * c[0] + d[0];
      g[1] = n * c[1] + d[1];
      g[2] = n * c[2] + d[2];

   } else {

      n = - (f1 / f2);

      g[0] = c[0] + n * d[0];
      g[1] = c[1] + n * d[1];
      g[2] = c[2] + n * d[2];
   }

   // normalise g
   n = 1.0 / sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);

   g[0]=g[0]*n;
   g[1]=g[1]*n;
   g[2]=g[2]*n;

   // determine p and q
   double p_[3], q_[3];

   p_[0]= g[0];
   p_[1]= g[1];
   p_[2]= g[2];

   q_[0]=-p_[0];
   q_[1]=-p_[1];
   q_[2]=-p_[2];

   // set p and q
   if (p != 0) {
      p[0] = p_[0];
      p[1] = p_[1];
      p[2] = p_[2];
   }
   if (q != 0) {
      q[0] = q_[0];
      q[1] = q_[1];
      q[2] = q_[2];
   }

   int result;

   result = 0;
   if (vector_is_between(a, b, p_, e_ab)) result |= 1 << 0;
   if (vector_is_between(a, b, q_, e_ab)) result |= 1 << 1;
   if (vector_is_between(c, d, p_, e_cd)) result |= 1 << 2;
   if (vector_is_between(c, d, q_, e_cd)) result |= 1 << 3;

   return result;
}

/** \brief compute the intersection point two circles of latitude
 *
 * compute the intersection points of two circle of latitude
 * if p and q are != NULL they contain the intersection points
 * the return value is:
 *      - 0 if the intersection points are neither between (a and b) or (c and d)
 *      - -1 if an error occurred
 *      - 1st bit will be set if p is between a and b
 *      - 2nd bit will be set if q is between a and b
 *      - 3rd bit will be set if p is between c and d
 *      - 4th bit will be set if q is between c and d
 *      - 5th bit will be set if both edges are on the same circle of latitude
 **/
int latcxlatc (struct edge edge_a, struct edge edge_b,
               struct point * p, struct point * q) {

   // two circles of latitude can only intersect if they are on the same latitude
   if (fabs(edge_a.points[0].lat - edge_b.points[0].lat) > tol)
      return 0;

   // check whether the two circles overlap
   
   double angle_ab;
   double angle_ac;
   double angle_ad;
   double angle_bc;
   double angle_bd;
   double angle_cd;

   angle_ab = fabs(get_angle(edge_a.points[0].lon, edge_a.points[1].lon));
   angle_ac = fabs(get_angle(edge_a.points[0].lon, edge_b.points[0].lon));
   angle_ad = fabs(get_angle(edge_a.points[0].lon, edge_b.points[1].lon));
   angle_bc = fabs(get_angle(edge_a.points[1].lon, edge_b.points[0].lon));
   angle_bd = fabs(get_angle(edge_a.points[1].lon, edge_b.points[1].lon));
   angle_cd = fabs(get_angle(edge_b.points[0].lon, edge_b.points[1].lon));

   int ret_value = 1 << 4;

   // if the both edges are on the north pole
   if (fabs(M_PI_2 - edge_a.points[0].lat) < tol) {

      if (p != NULL)
         p->lon = 0, p->lat = M_PI_2;

      if (q != NULL)
         q->lon = M_PI, q->lat = M_PI_2;

      ret_value |= (1 << 0) | (1 << 1) | (1 << 2) | (1 << 3);

   // if the both edges are on the south pole
   } else if (fabs(M_PI_2 + edge_a.points[0].lat) < tol) {

      if (p != NULL)
         p->lon = 0, p->lat = -M_PI_2;

      if (q != NULL)
         q->lon = M_PI, q->lat = -M_PI_2;

      ret_value |= (1 << 0) | (1 << 1) | (1 << 2) | (1 << 3);

   // if c is between a and b
   } else if (angle_ac + tol <= angle_ab && angle_bc + tol <= angle_ab) {

      if (p != NULL)
         *p = edge_b.points[0];

      if (q != NULL) {
         *q = edge_b.points[0];
         q->lon += M_PI;
      }

      ret_value |= (1 << 0) | (1 << 2);

   // if d is betwenn a and b
   } else if (angle_ad + tol <= angle_ab && angle_bd + tol <= angle_ab) {

      if (p != NULL)
         *p = edge_b.points[1];

      if (q != NULL) {
         *q = edge_b.points[1];
         q->lon += M_PI;
      }

      ret_value |= (1 << 0) | (1 << 2);

   // if a between c and d
   } else if (angle_cd + tol <= angle_ac && angle_cd + tol <= angle_ad) {

      if (p != NULL)
         *p = edge_a.points[0];

      if (q != NULL) {
         *q = edge_a.points[0];
         q->lon += M_PI;
      }

      ret_value |= (1 << 0) | (1 << 2);

   } else {

      if (p != NULL)
         *p = edge_a.points[0];

      if (q != NULL) {
         *q = edge_a.points[0];
         q->lon += M_PI;
      }

      ret_value |= 1 << 0;
   }

   return ret_value;
}

/** \brief compute the intersection point two circles of latitude
 *
 * compute the intersection points of two circle of latitude
 * if p and q are != NULL they contain the intersection points
 * the return value is:
 *      - 0 if the intersection points are neither between (a and b) or (c and d)
 *      - -1 if an error occurred
 *      - 1st bit will be set if p is between a and b
 *      - 2nd bit will be set if q is between a and b
 *      - 3rd bit will be set if p is between c and d
 *      - 4th bit will be set if q is between c and d
 *      - 5th bit will be set if both edges are on the same circle of latitude
 **/
int latcxlatc_vec (double a[3], double b[3], double c[3], double d[3],
                   double p[3], double q[3]) {

   // two circles of latitude can only intersect if they are on the same latitude
   if (fabs(a[2] - c[2]) > tol)
      return -1;

   int result = 16;

   double cross_ab[3], cross_cd[3];

   crossproduct(a, b, cross_ab);
   crossproduct(c, d, cross_cd);

   double a_[3] = {a[0], a[1], 0};
   double b_[3] = {b[0], b[1], 0};
   double c_[3] = {c[0], c[1], 0};
   double d_[3] = {d[0], d[1], 0};

   int a_between_cd, b_between_cd, c_between_ab, d_between_ab;

   a_between_cd = vector_is_between(c_, d_, a_, cross_cd);
   b_between_cd = vector_is_between(c_, d_, b_, cross_cd);
   c_between_ab = vector_is_between(a_, b_, c_, cross_ab);
   d_between_ab = vector_is_between(a_, b_, d_, cross_ab);

   if (a_between_cd && b_between_cd && c_between_ab && d_between_ab) {

      p[0] = a[0], p[1] = a[1], p[2] = a[2];
      q[0] = a[0], q[1] = a[1], q[2] = a[2];

      result |= 1 + 4;

   } else if (a_between_cd) {

      p[0] = a[0], p[1] = a[1], p[2] = a[2];

      result |= 1 + 2 + 4 + 8;

      if (b_between_cd) q[0] = b[0], q[1] = b[1], q[2] = b[2];
      else if (c_between_ab) q[0] = c[0], q[1] = c[1], q[2] = c[2];
      else if (d_between_ab) q[0] = d[0], q[1] = d[1], q[2] = d[2];
      else abort_message("internal error", __FILE__, __LINE__);

   } else if (b_between_cd) {

      p[0] = b[0], p[1] = b[1], p[2] = b[2];

      result |= 1 + 2 + 4 + 8;

      if (c_between_ab) q[0] = c[0], q[1] = c[1], q[2] = c[2];
      else if (d_between_ab) q[0] = d[0], q[1] = d[1], q[2] = d[2];
      else abort_message("internal error", __FILE__, __LINE__);

   } else if (c_between_ab && d_between_ab) {

      p[0] = c[0], p[1] = c[1], p[2] = c[2];
      q[0] = d[0], q[1] = d[1], q[2] = d[2];

      result |= 1 + 2 + 4 + 8;

   } else {

      p[0] = a[0], p[1] = a[1], p[2] = a[2];
      q[0] = b[0], q[1] = b[1], q[2] = b[2];

      result |= 1 + 2;
   }

   return result;
}

/** \brief compute the intersection point two circles of longitude
 *
 * compute the intersection points of two circle of longitude
 * if p and q are != NULL they contain the intersection points
 * the return value is:
 *      - 0 if the intersection points are neither between (a and b) or (c and d)
 *      - -1 if an error occurred
 *      - 1st bit will be set if p is between a and b
 *      - 2nd bit will be set if q is between a and b
 *      - 3rd bit will be set if p is between c and d
 *      - 4th bit will be set if q is between c and d
 *      - 5th bit will be set if both edges are on the same circle of longitude
 **/
int loncxlonc (struct edge edge_a, struct edge edge_b,
               struct point * p, struct point * q) {

   double angle_ac = fabs(get_angle(edge_a.points[0].lon, edge_b.points[0].lon));
   double angle_ab = fabs(get_angle(edge_a.points[0].lon, edge_a.points[1].lon));
   double angle_cd = fabs(get_angle(edge_b.points[0].lon, edge_b.points[1].lon));

   int ret_value = 0;

   // if both edges are on the same circle of longitude
   if ((angle_ac < tol) || (fabs(M_PI - angle_ac) < tol))
      ret_value |= 1 << 4;

   // check whether both edges cross a pole
   if (fabs(M_PI - angle_ab) < tol && fabs(M_PI - angle_cd) < tol) {

      if (p != NULL)
         p->lon = edge_a.points[0].lon, p->lat = M_PI_2;
      if (q != NULL)
         q->lon = edge_a.points[0].lon+M_PI, q->lat = -M_PI_2;

      //check if both edges cross the north pole
      if (edge_a.points[0].lat > 0 && edge_b.points[0].lat > 0)
         ret_value |= (1 << 0) | (1 << 2);

      // if both edges corr the south pole
       else if (edge_a.points[0].lat < 0 && edge_b.points[0].lat < 0)
         ret_value |= (1 << 1) | (1 << 3);

   // if both edges are on the same circle of longitude
   } else if (ret_value & (1 << 4) ) {

      // if the first edge crosses the pole
      if (fabs(M_PI - angle_ab) < tol) {

         if (angle_ac < tol) {
            edge_a.points[1].lon = edge_a.points[0].lon;
            edge_a.points[1].lat = (edge_a.points[1].lat > 0)?M_PI_2:-M_PI_2;
         } else {
            edge_a.points[0].lon = edge_a.points[1].lon;
            edge_a.points[0].lat = (edge_a.points[0].lat > 0)?M_PI_2:-M_PI_2;
         }

      // if the second edge crosse the pole
      } else if (fabs(M_PI - angle_cd) < tol) {

         if (angle_ac < tol) {
            edge_b.points[1].lon = edge_b.points[0].lon;
            edge_b.points[1].lat = (edge_b.points[1].lat > 0)?M_PI_2:-M_PI_2;
         } else {
            edge_b.points[0].lon = edge_b.points[1].lon;
            edge_b.points[0].lat = (edge_b.points[0].lat > 0)?M_PI_2:-M_PI_2;
         }
      }

      double lat_diff_ab = fabs(edge_a.points[0].lat - edge_a.points[1].lat);
      double lat_diff_ac = fabs(edge_a.points[0].lat - edge_b.points[0].lat);
      double lat_diff_ad = fabs(edge_a.points[0].lat - edge_b.points[1].lat);
      double lat_diff_bc = fabs(edge_a.points[1].lat - edge_b.points[0].lat);
      double lat_diff_bd = fabs(edge_a.points[1].lat - edge_b.points[1].lat);
      double lat_diff_cd = fabs(edge_b.points[0].lat - edge_b.points[1].lat);

      // if c is between a and b
      if (lat_diff_ab + tol > lat_diff_ac && lat_diff_ab + tol > lat_diff_bc) {

         if (p != NULL)
            *p = edge_b.points[0];

         if (q != NULL) {
            q->lon = edge_b.points[0].lon + M_PI;
            q->lat = - edge_b.points[0].lat;
         }

         ret_value |= (1 << 0) | (1 << 2);

      // if d is between a and b
      } else if (lat_diff_ab + tol > lat_diff_ad && lat_diff_ab + tol > lat_diff_bd) {

         if (p != NULL)
            *p = edge_b.points[1];

         if (q != NULL) {
            q->lon = edge_b.points[1].lon + M_PI;
            q->lat = - edge_b.points[1].lat;
         }

         ret_value |= (1 << 0) | (1 << 2);

      // if a is between c and d
      } else if (lat_diff_cd + tol > lat_diff_ac && lat_diff_cd + tol > lat_diff_ad) {

         if (p != NULL)
            *p = edge_a.points[0];

         if (q != NULL) {
            q->lon = edge_a.points[0].lon + M_PI;
            q->lat = - edge_a.points[0].lat;
         }

         ret_value |= (1 << 0) | (1 << 2);

      } else {

         if (p != NULL)
            *p = edge_a.points[0];

         if (q != NULL) {
            q->lon = edge_a.points[0].lon + M_PI;
            q->lat = - edge_a.points[0].lat;
         }

         ret_value |= 1 << 0;
      }
   }

   return ret_value;
}

/** \brief compute the intersection point two circles of longitude
 *
 * compute the intersection points of two circle of longitude
 * if p and q are != NULL they contain the intersection points
 * the return value is:
 *      - 0 if the intersection points are neither between (a and b) or (c and d)
 *      - -1 if an error occurred
 *      - 1st bit will be set if p is between a and b
 *      - 2nd bit will be set if q is between a and b
 *      - 3rd bit will be set if p is between c and d
 *      - 4th bit will be set if q is between c and d
 *      - 5th bit will be set if both edges are on the same circle of longitude
 **/
int loncxlonc_vec (double a[3], double b[3], double c[3], double d[3],
                   double p[3], double q[3]) {

   int ret_value = 0;

   double cross_ab[3], cross_cd[3];

   crossproduct(a, b, cross_ab);
   crossproduct(c, d, cross_cd);

   double norm_cross_ab[3], norm_cross_cd[3];

   norm_cross_ab[2] = 0;
   norm_cross_cd[2] = 0;

   if (fabs(cross_ab[0]) < tol) {
      norm_cross_ab[0] = 0;
      norm_cross_ab[1] = (cross_ab[1] > 0) ? 1 : -1;
   } else if (fabs(cross_ab[1]) < tol) {
      norm_cross_ab[0] = (cross_ab[0] > 0) ? 1 : -1;
      norm_cross_ab[1] = 0;
   } else {
      double scale = 1.0 / sqrt(cross_ab[0] * cross_ab[0] + cross_ab[1] * cross_ab[1]);
      norm_cross_ab[0] = cross_ab[0] * scale;
      norm_cross_ab[1] = cross_ab[1] * scale;
   }

   if (fabs(cross_cd[0]) < tol) {
      norm_cross_cd[0] = 0;
      norm_cross_cd[1] = (cross_cd[1] > 0) ? 1 : -1;
   } else if (fabs(cross_cd[1]) < tol) {
      norm_cross_cd[0] = (cross_cd[0] > 0) ? 1 : -1;
      norm_cross_cd[1] = 0;
   } else {
      double scale = 1.0 / sqrt(cross_cd[0] * cross_cd[0] + cross_cd[1] * cross_cd[1]);
      norm_cross_cd[0] = cross_cd[0] * scale;
      norm_cross_cd[1] = cross_cd[1] * scale;
   }

   // if both edges are on the same circle of longitude
   if ((fabs(norm_cross_ab[0] - norm_cross_cd[0]) < tol &&
        fabs(norm_cross_ab[1] - norm_cross_cd[1]) < tol) ||
       (fabs(norm_cross_ab[0] + norm_cross_cd[0]) < tol &&
        fabs(norm_cross_ab[1] + norm_cross_cd[1]) < tol)) {

      ret_value |= 16;

      int a_between_cd, b_between_cd, c_between_ab, d_between_ab;

      a_between_cd = vector_is_between(c, d, a, cross_cd);
      b_between_cd = vector_is_between(c, d, b, cross_cd);
      c_between_ab = vector_is_between(a, b, c, cross_ab);
      d_between_ab = vector_is_between(a, b, d, cross_ab);

      if (a_between_cd) {

         p[0] = a[0], p[1] = a[1], p[2] = a[2];

         ret_value |= 1 + 2 + 4 + 8;

         if (b_between_cd) q[0] = b[0], q[1] = b[1], q[2] = b[2];
         else if (c_between_ab) q[0] = c[0], q[1] = c[1], q[2] = c[2];
         else if (d_between_ab) q[0] = d[0], q[1] = d[1], q[2] = d[2];
         else abort_message("internal error", __FILE__, __LINE__);

      } else if (b_between_cd) {

         p[0] = b[0], p[1] = b[1], p[2] = b[2];

         ret_value |= 1 + 2 + 4 + 8;

         if (c_between_ab) q[0] = c[0], q[1] = c[1], q[2] = c[2];
         else if (d_between_ab) q[0] = d[0], q[1] = d[1], q[2] = d[2];
         else abort_message("internal error", __FILE__, __LINE__);

      } else if (c_between_ab && d_between_ab) {

         p[0] = c[0], p[1] = c[1], p[2] = c[2];
         q[0] = d[0], q[1] = d[1], q[2] = d[2];

         ret_value |= 1 + 2 + 4 + 8;

      } else {

         p[0] = 0, p[1] = 0, p[2] = 1;
         q[0] = 0, q[1] = 0, q[2] = -1;
      }

   } else {

      p[0] = 0, p[1] = 0; p[2] = 1;
      q[0] = 0, q[1] = 0; q[2] = -1;

      if (vector_is_between(a, b, p, cross_ab)) ret_value |= 1;
      if (vector_is_between(a, b, q, cross_ab)) ret_value |= 2;
      if (vector_is_between(c, d, p, cross_cd)) ret_value |= 4;
      if (vector_is_between(c, d, q, cross_cd)) ret_value |= 8;
   }

   return ret_value;
}

/** \brief compute the intersection point of a meridian and a parallel
 *
 * compute the intersection points of a circle of longitude (defined by a and b)
 * and a circle of latitude (defined by c and d)
 * if p and q are != NULL they contain the intersection points
 * the return value is:
 *      - 0 if the intersection points are neither between (a and b) or (c and d)
 *      - -1 if an error occurred
 *      - 1st bit will be set if p is between a and b
 *      - 2nd bit will be set if q is between a and b
 *      - 3rd bit will be set if p is between c and d
 *      - 4th bit will be set if q is between c and d
 **/
int loncxlatc (struct edge edge_a, struct edge edge_b,
               struct point * p, struct point * q) {

   double lon_a, lat_a[2];
   double lon_b[2], lat_b;

   unsigned ret_value;

   ret_value = 0;

   lon_a = edge_a.points[0].lon;

   if (edge_a.points[0].lat > edge_a.points[1].lat) {

      lat_a[0] = edge_a.points[1].lat;
      lat_a[1] = edge_a.points[0].lat;
   } else {
      lat_a[0] = edge_a.points[0].lat;
      lat_a[1] = edge_a.points[1].lat;
   }

   lat_b = edge_b.points[0].lat;

   if (edge_b.points[0].lon > edge_b.points[1].lon) {

      lon_b[0] = edge_b.points[1].lon;
      lon_b[1] = edge_b.points[0].lon;
   } else {
      lon_b[0] = edge_b.points[0].lon;
      lon_b[1] = edge_b.points[1].lon;
   }

   unsigned a_goes_across_pole;
   unsigned b_is_on_pole;

   a_goes_across_pole = fabs(get_angle(edge_a.points[0].lon,
                                           edge_a.points[1].lon)) > tol;
   b_is_on_pole = fabs(M_PI_2 - fabs(lat_b)) < tol;

   if (b_is_on_pole) {

      if (p != NULL) *p = edge_b.points[0];
      if (q != NULL) *q = edge_b.points[0], q->lon += M_PI;

      if (((a_goes_across_pole) && (fabs(lat_b - lat_a[0]) < M_PI_2)) ||
          ((fabs(lat_b - lat_a[0]) < tol) ||
           (fabs(lat_b - lat_a[1]) < tol))){

         ret_value |= 1 + 2;
      }

      ret_value |= 4 + 8;

   } else {

      if (p != NULL) p->lon = lon_a, p->lat = lat_b;
      if (q != NULL) q->lon = lon_a + M_PI, q->lat = lat_b;

      double angle_cd, angle_cp, angle_cq;

      angle_cd = get_angle(lon_b[0], lon_b[1]);
      angle_cp = get_angle(lon_b[0], lon_a);
      angle_cq = get_angle(lon_b[0], lon_a + M_PI);

      if (angle_cd > 0) {
         if (angle_cp <= angle_cd + tol && angle_cp > -tol) ret_value |= 1 << 2;
         if (angle_cq <= angle_cd + tol && angle_cq > -tol) ret_value |= 1 << 3;
      } else {
         if (angle_cp >= angle_cd - tol && angle_cp < tol) ret_value |= 1 << 2;
         if (angle_cq >= angle_cd - tol && angle_cq < tol) ret_value |= 1 << 3;
      }

      if (a_goes_across_pole) {

         if (lat_a[0] > 0.0) {

            if (lat_b > lat_a[0]) ret_value |= 1 << 0;
            if (lat_b > lat_a[1]) ret_value |= 1 << 1;

         } else {

            if (lat_b < lat_a[0]) ret_value |= 1 << 0;
            if (lat_b < lat_a[1]) ret_value |= 1 << 1;
         }
      } else if ((lat_b >= lat_a[0]) && (lat_b <= lat_a[1])) ret_value |= 1 << 0;
   }

   return ret_value;
}

/** \brief compute the intersection point of a meridian and a parallel
 *
 * compute the intersection points of a circle of longitude (defined by a and b)
 * and a circle of latitude (defined by c and d)
 * if p and q are != NULL they contain the intersection points
 * the return value is:
 *      - 0 if the intersection points are neither between (a and b) or (c and d)
 *      - -1 if an error occurred
 *      - 1st bit will be set if p is between a and b
 *      - 2nd bit will be set if q is between a and b
 *      - 3rd bit will be set if p is between c and d
 *      - 4th bit will be set if q is between c and d
 **/
int loncxlatc_vec (double a[3], double b[3], double c[3], double d[3],
                   double p[3], double q[3]) {

   unsigned ret_value;

   ret_value = 0;

   if (fabs(a[0] * b[1] - a[1] * b[0]) > tol)
      abort_message("edge is not a circle of longitude", __FILE__, __LINE__);

   unsigned ab_goes_across_pole;
   unsigned cd_is_on_pole;

   ab_goes_across_pole =
      ((fabs(1.0 - fabs(a[2])) < tol || fabs(1.0 - fabs(b[2])) < tol) ||
       (((a[0] > 0.0) ^ (b[0] > 0.0)) && ((a[1] > 0.0) ^ (b[1] > 0.0))));

   cd_is_on_pole = fabs(1.0 - fabs(c[2])) < tol;

   if (cd_is_on_pole) {

      if (((ab_goes_across_pole) && (fabs(a[2] - c[2]) < 1.0)) ||
          ((fabs(a[2] - c[2]) < tol)) || (fabs(b[2] - c[2]) < tol))
         ret_value |= 1;

         if (p != NULL)
            p[0] = c[0], p[1] = c[1], p[2] = c[2];
         if (q != NULL)
            q[0] = c[0], q[1] = c[1], q[2] = c[2];

      ret_value |= 4;

   } else {

      /*
      // the cos is too inaccurate close to the equator
      {
         if (fabs(a[2]) < fabs(b[2])) {

            double scale = cos(c[2]) / cos(a[2]);

            if (p != NULL)
               p[0] = a[0] * scale, p[1] = a[1] * scale, p[2] = c[2];
            if (q != NULL)
               q[0] = -a[0] * scale, q[1] = -a[1] * scale, q[2] = c[2];
         } else {d

            double scale = cos(c[2]) / cos(b[2]);

            if (p != NULL)
               p[0] = b[0] * scale, p[1] = b[1] * scale, p[2] = c[2];
            if (q != NULL)
               q[0] = -b[0] * scale, q[1] = -b[1] * scale, q[2] = c[2];
         }
      }
      */
      {

         if (fabs(a[2]) < fabs(b[2])) {

            double scale = sqrt((1.0 - c[2] * c[2])/
                                (a[0] * a[0] + a[1] * a[1]));

            if (p != NULL)
               p[0] = a[0] * scale, p[1] = a[1] * scale, p[2] = c[2];
            if (q != NULL)
               q[0] = -a[0] * scale, q[1] = -a[1] * scale, q[2] = c[2];

         } else {

            double scale = sqrt((1.0 - c[2] * c[2])/
                                (b[0] * b[0] + b[1] * b[1]));

            if (p != NULL)
               p[0] = b[0] * scale, p[1] = b[1] * scale, p[2] = c[2];
            if (q != NULL)
               q[0] = -b[0] * scale, q[1] = -b[1] * scale, q[2] = c[2];
         }
      }

      double angle_cd, angle_cp, angle_dp, angle_cq, angle_dq;

      angle_cd = get_vector_angle(c, d);
      angle_cp = get_vector_angle(c, p);
      angle_dp = get_vector_angle(d, p);
      angle_cq = get_vector_angle(c, q);
      angle_dq = get_vector_angle(d, q);

      if (angle_cp < tol || angle_dp < tol ||
          (angle_cp < angle_cd + tol && angle_dp < angle_cd + tol))
         ret_value |= 1 << 2;
      if (angle_cq < tol || angle_dq < tol ||
          (angle_cq < angle_cd + tol && angle_dq < angle_cd + tol))
         ret_value |= 1 << 3;

      double cross_ab[3];

      crossproduct(a, b, cross_ab);

      if (vector_is_between(a, b, p, cross_ab)) ret_value |= 1 << 0;
      if (vector_is_between(a, b, q, cross_ab)) ret_value |= 1 << 1;
   }

   return ret_value;
}

/** \brief compute the intersection of a great circle with the parallel
 *
 *  compute the intersection points of a great circle (defined by a and b)
 * and a circle of latitude (defined by c and d)
 * if p and q are != NULL they contain the intersection points
 * the return value is:
 *    - 0 if the intersection points are neither between (a and b) or (c and d)
 *    - -1 if the two circles do not intersect or an error occurred
 *    - 1st bit will be set if p is between a and b
 *    - 2nd bit will be set if q is between a and b
 *    - 3rd bit will be set if p is between c and d
 *    - 4th bit will be set if q is between c and d
 *    - 5th bit will be set if both circles are identically
 *
 *   based on
 *   - http://geospatialmethods.org/spheres/GCIntersect.html
 **/
int gcxlatc(struct edge edge_a, struct edge edge_b,
            struct point * p, struct point * q) {

   // if the great circle is nearly a lon circle, then the accuracy of the normal
   // computation gets messy, therefore we handle it as a lon circle...

   if (fabs(edge_a.points[0].lon - edge_a.points[1].lon) < tol) {

      edge_a.points[0].lon = edge_a.points[1].lon = (edge_a.points[0].lon + edge_a.points[1].lon) / 2.0;
      edge_a.edge_type = LON_CIRCLE;

      return loncxlatc(edge_a, edge_b, p, q);
   }

   // if the great circle is the equator, we handle the great circle as a circle of latitude
   if (fabs(edge_a.points[0].lat) < tol && fabs(edge_a.points[1].lat) < tol) {

      edge_a.points[0].lat = edge_a.points[1].lat = 0;
      edge_a.edge_type = LAT_CIRCLE;
      return latcxlatc(edge_a, edge_b, p, q);
   }

   double a[3], b[3], c[3], d[3], p_[3], q_[3];

   LLtoXYZ(edge_a.points[0].lon, edge_a.points[0].lat, a);
   LLtoXYZ(edge_a.points[1].lon, edge_a.points[1].lat, b);
   LLtoXYZ(edge_b.points[0].lon, edge_b.points[0].lat, c);
   LLtoXYZ(edge_b.points[1].lon, edge_b.points[1].lat, d);

   int ret_value = gcxlatc_vec(a, b, c, d, p_, q_);

   if (ret_value == -1) return -1;

   if (p != NULL) XYZtoLL(p_, &p->lon, &p->lat);
   if (q != NULL) XYZtoLL(q_, &q->lon, &q->lat);

   return ret_value;
}

/** \brief compute the intersection of a great circle with the parallel
 *
 *  compute the intersection points of a great circle (defined by a and b)
 * and a circle of latitude (defined by c and d)
 * if p and q are != NULL they contain the intersection points
 * the return value is:
 *    - 0 if the intersection points are neither between (a and b) or (c and d)
 *    - -1 if the two circles do not intersect or an error occurred
 *    - 1st bit will be set if p is between a and b
 *    - 2nd bit will be set if q is between a and b
 *    - 3rd bit will be set if p is between c and d
 *    - 4th bit will be set if q is between c and d
 *    - 5th bit will be set if both circles are identically
 * \remarks if -1 is returned neither p or q is set
 * \remarks if the two circles only have one intersection point,
 *          p and q will be identically, but only the p bits will be set
 **/

int gcxlatc_vec(double a[3], double b[3], double c[3], double d[3],
                double p[3], double q[3]) {

   unsigned result = 0;

   // if the great circle is the equator
   if (fabs(a[2]) < tol && fabs(b[2]) < tol) {

      // if the circle of latitude is also the equator
      if (fabs(c[2]) < tol) {

         result |= 16;

         double cross_ab[3], cross_cd[3];

         crossproduct(a, b, cross_ab);
         crossproduct(c, d, cross_cd);

         int a_between_cd, b_between_cd, c_between_ab, d_between_ab;

         a_between_cd = vector_is_between(c, d, a, cross_cd);
         b_between_cd = vector_is_between(c, d, b, cross_cd);
         c_between_ab = vector_is_between(a, b, c, cross_ab);
         d_between_ab = vector_is_between(a, b, d, cross_ab);

         if (a_between_cd) {

            p[0] = a[0], p[1] = a[1], p[2] = a[2];

            result |= 1 + 2 + 4 + 8;

            if (b_between_cd) q[0] = b[0], q[1] = b[1], q[2] = b[2];
            else if (c_between_ab) q[0] = c[0], q[1] = c[1], q[2] = c[2];
            else if (d_between_ab) q[0] = d[0], q[1] = d[1], q[2] = d[2];
            else abort_message("internal error", __FILE__, __LINE__);

         } else if (b_between_cd) {

            p[0] = b[0], p[1] = b[1], p[2] = b[2];

            result |= 1 + 2 + 4 + 8;

            if (c_between_ab) q[0] = c[0], q[1] = c[1], q[2] = c[2];
            else if (d_between_ab) q[0] = d[0], q[1] = d[1], q[2] = d[2];
            else abort_message("internal error", __FILE__, __LINE__);

         } else if (c_between_ab && d_between_ab) {

            p[0] = c[0], p[1] = c[1], p[2] = c[2];
            q[0] = d[0], q[1] = d[1], q[2] = d[2];

            result |= 1 + 2 + 4 + 8;

         } else {

            p[0] = c[0], p[1] = c[1], p[2] = c[2];
            q[0] = d[0], q[1] = d[1], q[2] = d[2];
         }
      }
   }

   double t[3], s[3];

   if (fabs(a[2]) > fabs(b[2])) {

      double scale = c[2] / a[2];

      t[0] = scale * a[0];
      t[1] = scale * a[1];
      
   } else {

      double scale = c[2] / b[2];

      t[0] = scale * b[0];
      t[1] = scale * b[1];
   }

   t[2] = c[2];

   s[2] = 0;

   if (fabs(a[2]) < tol)
      s[0] = a[0], s[1] = a[1];
   else if (fabs(b[2]) < tol)
      s[0] = b[0], s[1] = b[1];
   else if (fabs(a[2]) > fabs(b[2])) {
      double scale = b[2] / a[2];
      s[0] = b[0] - scale * a[0];
      s[1] = b[1] - scale * a[1];
   } else {
      double scale = a[2] / b[2];
      s[0] = a[0] - scale * b[0];
      s[1] = a[1] - scale * b[1];
   }

   if (fabs(s[0]) < tol && fabs(s[1]) < tol)
      abort_message("internal error", __FILE__, __LINE__);

   {
      // the intersection of the planes of both circles is defined by:
      // x = t + n * s

      // x_0^2 + x_1^2 + x_2^2 = 1
      // x_2 = c_2

      double a_ = s[0] * s[0] + s[1] * s[1];
      double b_ = 2.0 * (t[0] * s[0] + t[1] * s[1]);
      double c_ = t[0] * t[0] + t[1] * t[1] + c[2] * c[2] - 1.0;

      double temp = b_ * b_ - 4.0 * a_ * c_;

      // no intersection possible
      if (temp < 0.0) {
        if (temp < -tol)
          return -1;
        else
          temp = 0;
      }

      double n[2];

      if (fabs(temp) < tol) {

         n[0] = n[1] = - b_ / (2.0 * a_);

      } else {

         n[0] = - (b_ + sqrt(temp)) / (2.0 * a_);
         n[1] = - (b_ - sqrt(temp)) / (2.0 * a_);
      }

      p[0] = t[0] + n[0] * s[0];
      p[1] = t[1] + n[0] * s[1];
      p[2] = t[2] + n[0] * s[2];

      double cross_ab[3], cross_cd[3] = {0, 0, c[0]*d[1]-c[1]*d[0]};
      double temp_c[3] = {c[0], c[1], 0};
      double temp_d[3] = {d[0], d[1], 0};
      double temp_p[3] = {p[0], p[1], 0};

      crossproduct(a, b, cross_ab);

      if (vector_is_between(a, b, p, cross_ab)) result |= 1;
      if (vector_is_between(temp_c, temp_d, temp_p, cross_cd)) result |= 4;

      if (fabs(n[0] - n[1]) >= tol) {

         q[0] = t[0] + n[1] * s[0];
         q[1] = t[1] + n[1] * s[1];
         q[2] = t[2] + n[1] * s[2];

         double temp_q[3] = {q[0], q[1], 0};

         if (vector_is_between(a, b, q, cross_ab)) result |= 2;
         if (vector_is_between(temp_c, temp_d, temp_q, cross_cd)) result |= 8;
      } else
         q[0] = p[0], q[1] = p[1], q[2] = p[2];
   }

   return result;
}

int intersect (struct edge const edge_a, struct edge const edge_b,
               struct point * intersection) {

   int switch_edges;

   switch_edges = 0;

   int (*intersect_func)(struct edge, struct edge, struct point *, struct point *);

   // if both edges are on circles of latitude
   if (edge_a.edge_type == LAT_CIRCLE && edge_b.edge_type == LAT_CIRCLE) {

      intersect_func = latcxlatc;

   // if both edges are on circle of longitude
   } else if (edge_a.edge_type == LON_CIRCLE && edge_b.edge_type == LON_CIRCLE) {

      intersect_func = loncxlonc;

   // if both edges are on great circles
   } else if ((edge_a.edge_type == GREAT_CIRCLE &&
        edge_b.edge_type == GREAT_CIRCLE) ||
       (edge_a.edge_type == LON_CIRCLE   &&
        edge_b.edge_type == GREAT_CIRCLE) ||
       (edge_a.edge_type == GREAT_CIRCLE &&
        edge_b.edge_type == LON_CIRCLE)) {

      intersect_func = gcxgc;

   // if one edge a is on a great circle and edge b on a circle of latitude
   } else if (edge_a.edge_type == GREAT_CIRCLE &&
              edge_b.edge_type == LAT_CIRCLE) {

      intersect_func = gcxlatc;

   // if one edge a is on a circle of latitude and edge b on a great circle
   } else if (edge_a.edge_type == LAT_CIRCLE &&
              edge_b.edge_type == GREAT_CIRCLE ) {

      switch_edges = 1;
      intersect_func = gcxlatc;

   // if one edge a is on a circle of longitude and edge b on a circle of latitude
   } else if (edge_a.edge_type == LON_CIRCLE &&
              edge_b.edge_type == LAT_CIRCLE) {

      intersect_func = loncxlatc;

   // if one edge a is on a circle of latitude and edge b on a circle of longitude
   } else if (edge_a.edge_type == LAT_CIRCLE &&
              edge_b.edge_type == LON_CIRCLE ) {

      switch_edges = 1;
      intersect_func = loncxlatc;

   } else {

      abort_message ( "ERROR: unknown edge type.", __FILE__, __LINE__ );
      exit(EXIT_FAILURE);
   }

   unsigned const p_on_a = 1 << 0;
   unsigned const q_on_a = 1 << 1;
   unsigned const p_on_b = 1 << 2;
   unsigned const q_on_b = 1 << 3;

   int ret_value;

   struct point p, q;

   // compute the intersection between both circles
   if (switch_edges) ret_value = intersect_func(edge_b, edge_a, &p, &q);
   else              ret_value = intersect_func(edge_a, edge_b, &p, &q);

   // check whether the circles defined by the edges intersect
   if (ret_value > 0) {

      // check for intersection of the edges
      if ( (ret_value & p_on_a) && (ret_value & p_on_b) ) {

         if (intersection != NULL) *intersection = p;
         return 1;

      } else if ( (ret_value & q_on_a) && (ret_value & q_on_b) ) {

         if (intersection != NULL) *intersection = q;
         return 1;
      }
   }
   return 0;
}

int intersect_vec (enum edge_type edge_type_a, double a[3], double b[3],
                   enum edge_type edge_type_b, double c[3], double d[3],
                   double p[3], double q[3]) {

   int switch_edges;

   switch_edges = 0;

   int (*intersect_func)(double *, double *, double *, double *, double *, double *);

   // if both edges are on circles of latitude
   if (edge_type_a == LAT_CIRCLE &&
       edge_type_b == LAT_CIRCLE) {

      intersect_func = latcxlatc_vec;

   // if both edges are on circle of longitude
   } else if (edge_type_a == LON_CIRCLE &&
              edge_type_b == LON_CIRCLE) {

      intersect_func = loncxlonc_vec;

   // if both edges are on great circles
   } else if ((edge_type_a == GREAT_CIRCLE &&
               edge_type_b == GREAT_CIRCLE) ||
              (edge_type_a == LON_CIRCLE   &&
               edge_type_b == GREAT_CIRCLE) ||
              (edge_type_a == GREAT_CIRCLE &&
               edge_type_b == LON_CIRCLE)) {

      intersect_func = gcxgc_vec;

   // if one edge a is on a great circle and edge b on a circle of latitude
   } else if (edge_type_a == GREAT_CIRCLE &&
              edge_type_b == LAT_CIRCLE) {

      intersect_func = gcxlatc_vec;

   // if one edge a is on a circle of latitude and edge b on a great circle
   } else if (edge_type_a == LAT_CIRCLE &&
              edge_type_b == GREAT_CIRCLE ) {

      switch_edges = 1;
      intersect_func = gcxlatc_vec;

   // if one edge a is on a circle of longitude and edge b on a circle of latitude
   } else if (edge_type_a == LON_CIRCLE &&
              edge_type_b == LAT_CIRCLE) {

      intersect_func = loncxlatc_vec;

   // if one edge a is on a circle of latitude and edge b on a circle of longitude
   } else if (edge_type_a == LAT_CIRCLE &&
              edge_type_b == LON_CIRCLE ) {

      switch_edges = 1;
      intersect_func = loncxlatc_vec;

   } else {

      abort_message ( "ERROR: unknown edge type.", __FILE__, __LINE__ );
      exit(EXIT_FAILURE);
   }

   int ret_value;

   // compute the intersection between both circles
   if (switch_edges) ret_value = intersect_func(c, d, a, b, p, q);
   else              ret_value = intersect_func(a, b, c, d, p, q);

   if (switch_edges)
      ret_value = (ret_value & (~(1 + 2 + 4 + 8))) +
                  ((ret_value & (1 + 2)) << 2) +
                  ((ret_value & (4 + 8)) >> 2);

   return ret_value;
}

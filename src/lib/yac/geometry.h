/**
 * @file geometry.h
 * @brief Structs and interfaces for investigating the geometry and relations of cells
 *
 * The functions are used to determine relations between
 * source and target cells. Complete cells are constructed
 * on the sphere by connecting cell vertices along either
 * great circles (along the orthodrome) or directly (along
 * the loxodrome).
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

/** \example test_geometry.c
 * This contains some tests for basic routines of \ref geometry.h.
 */

#ifndef GEOMETRY_H
#define GEOMETRY_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "grid.h"
#include "utils.h"

#ifndef M_SQRT1_2
#define M_SQRT1_2 (0.707106781186547524401) /* sqrt(0.5) */
#endif

#ifndef M_SQRT3_4
#define M_SQRT3_4 (0.866025403784438646764) /* sqrt(0.75) */
#endif

#define SIN_COS_M_PI_2 ((struct sin_cos_angle){.sin = 1.0, .cos = 0.0}) /* PI/2 */
#define SIN_COS_M_PI ((struct sin_cos_angle){.sin = 0.0, .cos = -1.0}) /* PI */
#define SIN_COS_ZERO ((struct sin_cos_angle){.sin = 0.0, .cos = 1.0}) /* 0.0 */
#define SIN_COS_TOL \
  ((struct sin_cos_angle){.sin = yac_angle_tol, .cos = yac_cos_angle_tol}) /* yac_angle_tol */
#define SIN_COS_LOW_TOL \
  ((struct sin_cos_angle){.sin = yac_angle_low_tol, .cos = yac_cos_angle_low_tol}) /* yac_angle_low_tol */

extern const double yac_angle_tol;
extern const double yac_cos_angle_tol;
extern const double yac_angle_low_tol;
extern const double yac_cos_angle_low_tol;

struct line {
   struct {
      double x;
      double y;
   } p1, p2;
};

struct point {
   double lon, lat;
};

struct edge {
   struct point points[2];
   enum yac_edge_type edge_type;
};

struct sin_cos_angle {
  double sin, cos;
};

/**
 * defines a spherical cap, which is used as a convex hull to describe the extents of subdomain
 */
struct bounding_circle {

   //! the middle point of the spherical cap in cartesian coordinates (is a unit vector)
   double base_vector[3];
   //! angle between the middle point and the boundary of the spherical cap
   struct sin_cos_angle inc_angle;
};

/** \example test_check_overlap.c
 * This contains examples on how to use check_overlap_cells.
 */

/**
 * checks whether two cells overlap \n
 * @param[in] cell_a
 * @param[in] cell_b
 * @return 0 if both cells do not overlap
 * @see check_overlap_cells
 */
int yac_check_overlap_cells (struct grid_cell const cell_a, 
                             struct grid_cell const cell_b);

/**
 * checks whether two cells overlap \n
 * @param[in] cell_a
 * @param[in] circle_a bounding circle of cell a
 * @param[in] cell_b
 * @param[in] circle_b bounding circle of cell b
 * @return 0 if both cells do not overlap
 * @see check_overlap_cells
 */
int yac_check_overlap_cells2 (struct grid_cell const cell_a,
                              struct bounding_circle circle_a,
                              struct grid_cell const cell_b,
                              struct bounding_circle circle_b);

/** \example test_point_in_cell.c
 * This contains examples on how to use point_in_cell.
 */

/**
 * checks whether a given point is within a given cell \n
 * @param[in] point_coords
 * @param[in] cell
 * @return 0 if the point is not in the cell
 */
int yac_point_in_cell (double point_coords[3], struct grid_cell cell);

/**
 * checks whether a given point is within a given cell \n
 * @param[in] point_coords
 * @param[in] cell
 * @param[in] bnd_circle
 * @return 0 if the point is not in the cell
 */
int yac_point_in_cell2 (double point_coords[3], struct grid_cell cell,
                        struct bounding_circle bnd_circle);

/**
 * computes the angle between two longitude coordinates (in rad) \n
 * takes into account that longitude have a period of 2 pi
 * @param[in] a_lon
 * @param[in] b_lon
 * @return angle between both coordinates (in rad)
 */
static inline double get_angle (double a_lon, double b_lon) {
   double diff = a_lon - b_lon;
#if defined(CDO)
   return diff - lround(diff / (2.0 * M_PI)) * (2.0 * M_PI);
#else
   return diff - round(diff / (2.0 * M_PI)) * (2.0 * M_PI);
#endif
}

/** \example test_find_overlap.c
 * This contains an example on how to use find_overlapping_cells.
 */

/**
 * searches for overlapping cells between grids \n
 * the user needs to provide a dependency list containing an initial search result that is used as a starting point, which is then refined
 * @param[in] src_grid source %grid
 * @param[in] tgt_grid target %grid
 * @param[in] initial_src_to_tgt_dep dependency list containing for each source cell a list of target cell that might overlap
 * @param[out] src_to_tgt_dep dependency list containing the result of the search
 */
void yac_find_overlapping_cells (struct grid * src_grid, struct grid * tgt_grid,
                                 struct dep_list initial_src_to_tgt_dep,
                                 struct dep_list * src_to_tgt_dep);

/**
 * searches for all overlapping cells of a given grid that overlap with a
 * given cell\n
 * the user needs to provide initial search results that are used as a starting
 * point, which is then refined
 * @param[in]     src_cell             cell for which the search is to conducted
 * @param[in]     src_bnd_circle       bounding circle of source cell
 * @param[in]     tgt_grid             grid on which overlapping cells are to
 *                                     be searched
 * @param[in]     initial_dep          initial search results
 * @param[in]     num_initial_deps     number of dependencies in initial_dep
 * @param[in,out] deps                 NULL or pointer to pointer returned by a
 *                                     standard allocation routine
 * @param[in,out] deps_size            number of elements that fit into deps or
 *                                     0 if deps equals NULL
 * @param[out]    num_deps             number of overlapping cells
 * @param[in]     src_index            any number (no entry in
 *                                     tgts_already_touched is allowed to be
 *                                     bigger than src_index)
 * @param[in,out] tgts_already_touched array with the size of the target grid
 *                                     and no entry bigger than src_index
 * @param[in,out] stack                NULL or pointer to pointer returned by a
 *                                     standard allocation routine
 * @param[in,out] stack_size           number of elements that fit into stack
 *                                     or 0 if stack equals NULL
 * @remarks the user is responsible to free the memory associated with deps
 *          and stack
 */
void yac_find_overlapping_cells_s (struct grid_cell src_cell,
                                   struct bounding_circle src_bnd_circle,
                                   struct grid * tgt_grid,
                                   unsigned const * initial_dep,
                                   unsigned num_initial_deps, unsigned ** deps,
                                   unsigned * deps_size, unsigned * num_deps,
                                   unsigned src_index,
                                   unsigned * tgts_already_touched,
                                   unsigned ** stack, unsigned * stack_size);

/** \example test_gcxgc.c
 * This contains examples on how to use gcxgc and gcxgc_vec
 */

/**
 * computes the intersection points of two great circles \n
 * based on http://www.geoclub.de/viewtopic.php?f=54&t=29689
 * @param[in] edge_a edge defining a great circle
 * @param[in] edge_b edge defining a great circle
 * @param[out] p intersection point
 * @param[out] q intersection point
 * @return  0 if the intersection points are neither on edge a or b \n
 *          -1 if an error occurred \n
 *          1st bit will be set if p is on edge a \n
 *          2nd bit will be set if q is on edge a \n
 *          3rd bit will be set if p is on edge b \n
 *          4th bit will be set if q is on edge b
  *         5th bit will be set if both great circles are identically
 */
int yac_gcxgc (struct edge edge_a, struct edge edge_b,
               struct point * p, struct point * q);

/**
 * computes the intersection points of two great circles \n
 * based on http://www.geoclub.de/viewtopic.php?f=54&t=29689
 * @param[in] a first point of edge a the is on a great circle
 * @param[in] b second point of edge a the is on a great circle
 * @param[in] c first point of edge b the is on a great circle
 * @param[in] d second point of edge b the is on a great circle
 * @param[out] p intersection point
 * @param[out] q intersection point
 * @return  0 if the intersection points are neither on edge a or b \n
 *          -1 if an error occurred \n
 *          1st bit will be set if p is on edge a \n
 *          2nd bit will be set if q is on edge a \n
 *          3rd bit will be set if p is on edge b \n
 *          4th bit will be set if q is on edge b
  *         5th bit will be set if both great circles are identically
 */
 int yac_gcxgc_vec (double a[3], double b[3], double c[3], double d[3],
                    double p[3], double q[3]);

/** \example test_loncxlatc.c
 * This contains examples on loncxlatc and loncxlatc_vec
 */

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
int yac_loncxlatc_vec (double a[3], double b[3], double c[3], double d[3],
                       double p[3], double q[3]);

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
int yac_loncxlatc (struct edge edge_a, struct edge edge_b,
                   struct point * p, struct point * q);

/** \example test_latcxlatc.c
 * This contains examples on latcxlatc and latcxlatc_vec
 */

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
int yac_latcxlatc_vec (double a[3], double b[3], double c[3], double d[3],
                       double p[3], double q[3]);

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
int yac_latcxlatc (struct edge edge_a, struct edge edge_b,
                   struct point * p, struct point * q);

/** \example test_loncxlonc.c
 * This contains examples on loncxlonc loncxlonc_vec
 */

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
int yac_loncxlonc_vec (double a[3], double b[3], double c[3], double d[3],
                       double p[3], double q[3]);

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
int yac_loncxlonc (struct edge edge_a, struct edge edge_b,
                   struct point * p, struct point * q);

/** \example test_gcxlatc.c
 * This contains examples on gcxlatc and gcxlatc_vec.
 */

/**
 * compute the intersection points of a great circles
 * and a circle of latitude \n
 * based on http://geospatialmethods.org/spheres/GCIntersect.html
 * @param[in] edge_a edge defining a great circle
 * @param[in] edge_b edge defining a circle of latitude
 * @param[out] p intersection point
 * @param[out] q intersection point
 * @return  0 if the intersection points are neither on edge a or b \n
 *         -1 if the two circles do not intersect \n
 *          1st bit will be set if p is on edge a \n
 *          2nd bit will be set if q is on edge a \n
 *          3rd bit will be set if p is on edge b \n
 *          4th bit will be set if q is on edge b
 *
 * \remarks if -1 is returned neither p or q is set
 * \remarks if the two circles only have one intersection point,
 *          p and q will be identically, but only the p bits will be set
 */
int yac_gcxlatc (struct edge edge_a, struct edge edge_b,
                 struct point * p, struct point * q);

/**
 * compute the intersection points of a great circles
 * and a circle of latitude \n
 * based on http://geospatialmethods.org/spheres/GCIntersect.html
 * @param[in] a first point of edge on a great circle
 * @param[in] b second point of edge on a great circle
 * @param[in] c first point of edge on a circle of latitude
 * @param[in] d second point of edge on a circle of latitude
 * @param[out] p intersection point
 * @param[out] q intersection point
 * @return  0 if the intersection points are neither on edge a or b \n
 *         -1 if the two circles do not intersect \n
 *          1st bit will be set if p is on edge a \n
 *          2nd bit will be set if q is on edge a \n
 *          3rd bit will be set if p is on edge b \n
 *          4th bit will be set if q is on edge b
 *
 * \remarks if -1 is returned neither p or q is set
 * \remarks if the two circles only have one intersection point,
 *          p and q will be identically, but only the p bits will be set
 */
int yac_gcxlatc_vec(double a[3], double b[3], double c[3], double d[3],
                    double p[3], double q[3]);

/**
 * computes the intersection points of two edges
 * @param[in]  edge_a       edge a
 * @param[in]  edge_b       edge b
 * @param[out] intersection intersection point
 * @return  0 if the edges do not intersect \n
 *          1 if the edges intersect
 *
 * \remarks if intersection is not NULL and the edges intersect
 *          the intersection point is returned
 * \remarks only one intersection is returned even if the edges intersect twice
 */
int yac_intersect (struct edge const edge_a, struct edge const edge_b,
                   struct point * intersection);

/**
 * computes the intersection points of two edges
 * @param[in]  edge_type_a type of edge a
 * @param[in]  a           first point of edge a
 * @param[in]  b           second point of edge a
 * @param[in]  edge_type_b type of edge b
 * @param[in]  c           first point of edge b
 * @param[in]  d           second point of edge b
 * @param[out] p           first intersection point
 * @param[out] q           second intersection point
 * @return  0 if the intersection points are neither on edge a or b \n
 *         -1 if the two circles do not intersect \n
 *          1st bit will be set if p is on edge a \n
 *          2nd bit will be set if q is on edge a \n
 *          3rd bit will be set if p is on edge b \n
 *          4th bit will be set if q is on edge b \n
 *          5th bit will be set if both circles are identically
 *
 * \remarks if -1 is returned neither p or q is set
 * \remarks if the two circles only have one intersection point,
 *          p and q will be identically, but only the p bits will be set
 */
int yac_intersect_vec (enum yac_edge_type edge_type_a, double a[3], double b[3],
                       enum yac_edge_type edge_type_b, double c[3], double d[3],
                       double p[3], double q[3]);

/** \example test_cell_bnd_circle.c
 * These are some examples on how to use \ref yac_get_cell_bounding_circle.
 */

/**
 * gets the bounding circle for a grid cell
 * @param[in] cell grid cell (coordinates have to be in radian)
 * @param[out] bnd_circle bounding circle of the grid cell
 */
void yac_get_cell_bounding_circle(struct grid_cell cell,
                                  struct bounding_circle * bnd_circle);

/**
 * computes the circumscribe circle for a triangle on the sphere
 * @param[in]  a          coordinates of first point (xyz)
 * @param[in]  b          coordinates of second point (xyz)
 * @param[in]  c          coordinates of thrid point (xyz)
 * @param[out] bnd_circle circumscribe circle
 * @remark it is assumed that all three edges of the triangle are great circles
 */
void yac_get_cell_circumscribe_circle_unstruct_triangle(
   double a[3], double b[3], double c[3], struct bounding_circle * bnd_circle);

/**
 * computes the smallest bounding circle for a triangle on the sphere
 * @param[in]  a          coordinates of first point (xyz)
 * @param[in]  b          coordinates of second point (xyz)
 * @param[in]  c          coordinates of thrid point (xyz)
 * @param[out] bnd_circle bounding circle
 * @remark it is assumed that all three edges of the triangle are great circles
 */
void yac_get_cell_bounding_circle_unstruct_triangle(
   double a[3], double b[3], double c[3], struct bounding_circle * bnd_circle);

/**
 * computes the circumscribe circle for a quad on the sphere
 * @param[in]  a          coordinates of first point (xyz)
 * @param[in]  b          coordinates of second point (xyz)
 * @param[in]  c          coordinates of thrid point (xyz)
 * @param[in]  d          coordinates of fourth point (xyz)
 * @param[out] bnd_circle circumscribe circle
 * @remark it is assumed that all edges of the quad are either circles of
 *         longitude or latitude
 */
void yac_get_cell_circumscribe_circle_reg_quad(
   double a[3], double b[3], double c[3], double d[3],
   struct bounding_circle * bnd_circle);

/**
 * computes the smallest bounding circle for a triangle on the sphere
 * @param[in]  a          coordinates of first point (xyz)
 * @param[in]  b          coordinates of second point (xyz)
 * @param[in]  c          coordinates of thrid point (xyz)
 * @param[in]  d          coordinates of fourth point (xyz)
 * @param[out] bnd_circle bounding circle
 * @remark it is assumed that all edges of the quad are either circles of
 *         longitude or latitude
 */
void yac_get_cell_bounding_circle_reg_quad(
   double a[3], double b[3], double c[3], double d[3],
   struct bounding_circle * bnd_circle);

/** \example test_grid_bnd_circle.c
 * These are some examples on how to use \ref yac_get_grid_bounding_circle.
 */

/**
 * gets a bounding circle for a grid
 * @param[in] grid
 * @param[out] bnd_circle bounding circle of the grid
 */
void yac_get_grid_bounding_circle(struct grid * grid,
                                  struct bounding_circle * bnd_circle);

/**
 * gets all cells of a grid that have an overlap with the area defined by the given bounding circle
 * @param[in] grid
 * @param[in] extent bounding circle defining an area against which all cells of the grid are going to checked
 * @param[in,out] matching_cells array containing the matching cells (if required the array will be reallocated)
 * @param[in,out] curr_matching_cells_array_size size of matching_cells array (if the matching_cells array is reallocated the size will be updated as well)
 * @param[in,out] local_ids local ids of the matching cells (if required the array will be reallocated)
 * @param[in,out] curr_local_ids_array_size size of matching_cells array (if the local_ids array is reallocated the size will be updated as well)
 * @param[out] num_matching_cells number of matching cells
 * @param[in] offset number of cells/local ids already in the matching_cells and local_ids array (these will not be overwritten)
 *
 * \remark *matching_cells == NULL and *curr_matching_cells_array_size == 0 are valid input values
 * \remark *local_ids == NULL and *curr_local_ids_array_size == 0 are valid input values
 * \remark if matching_cells == NULL or local_ids == NULL the respective data is not returned by this routine
 */
void yac_get_matching_grid_cells(struct grid * grid, struct bounding_circle extent,
                                 struct grid_cell ** matching_cells,
                                 unsigned * curr_matching_cells_array_size,
                                 unsigned ** local_ids, unsigned * curr_local_ids_array_size,
                                 unsigned * num_matching_cells, unsigned offset);

/**
 * checks whether two extents overlap
 * @param[in] extent_a bounding circle
 * @param[in] extent_b bounding circle
 * @return 0 if the bounding circles do not overlap
 */
unsigned yac_extents_overlap(struct bounding_circle * extent_a,
                             struct bounding_circle * extent_b);

/**
 * checks whether a point is within a bounding circle
 * @param[in] point point to be checked (coordinates in radian)
 * @param[in] bnd_circle bounding circle
 * @return 0 if point is not within the bounding circle
 */
unsigned yac_point_in_bounding_circle(struct point point,
                                      struct bounding_circle * bnd_circle);

/**
 * checks whether a point is within a bounding circle
 * @param[in] point_vector point to be checked
 * @param[in] bnd_circle bounding circle
 * @return 0 if point is not within the bounding circle
 */
unsigned yac_point_in_bounding_circle_vec(double point_vector[3],
                                          struct bounding_circle * bnd_circle);

/**
 * converts lon-lat coordinates into xyz ones
 *
 * Further information:
 * http://en.wikipedia.org/wiki/List_of_common_coordinate_transformations
 *
 * @param[in]  lon   longitude coordinates in radian
 * @param[in]  lat   latitude coordinates in radian
 * @param[out] p_out xyz coordinates
 */
static inline void LLtoXYZ(double lon, double lat, double p_out[]) {

   while (lon < -M_PI) lon += 2.0 * M_PI;
   while (lon >= M_PI) lon -= 2.0 * M_PI;

   double cos_lat = cos(lat);
   p_out[0] = cos_lat * cos(lon);
   p_out[1] = cos_lat * sin(lon);
   p_out[2] = sin(lat);
}

/**
 * converts lon-lat coordinates into xyz ones
 * @param[in]  lon   longitude coordinates in deg
 * @param[in]  lat   latitude coordinates in deg
 * @param[out] p_out xyz coordinates
 */
static inline void LLtoXYZ_deg(double lon, double lat, double p_out[]) {
   LLtoXYZ(lon*rad, lat*rad, p_out);
}

/**
 * converts lon-lat coordinates into xyz ones
 *
 * Further information:
 * http://en.wikipedia.org/wiki/List_of_common_coordinate_transformations
 *
 * @param[in]  p_in xyz coordinates
 * @param[out] lon  longitude coordinate in radian
 * @param[out] lat  latitude coordinate in radian
 */
static inline void XYZtoLL (double p_in[], double * lon, double * lat) {

   *lon = atan2(p_in[1] , p_in[0]);
   *lat = M_PI_2 - acos(p_in[2]);
}

static inline void crossproduct_ld (double a[], double b[], double cross[]) {

/* crossproduct in Cartesian coordinates */

   long double a_[3] = {a[0], a[1], a[2]};
   long double b_[3] = {b[0], b[1], b[2]};

   cross[0] = a_[1] * b_[2] - a_[2] * b_[1];
   cross[1] = a_[2] * b_[0] - a_[0] * b_[2];
   cross[2] = a_[0] * b_[1] - a_[1] * b_[0];
}

/**
 * for small angles <= 1e-?8? the crossproduct is inaccurate\n
 * use \ref crossproduct_ld for these cases
 */
static inline void crossproduct_d (double a[], double b[], double cross[]) {

/* crossproduct in Cartesian coordinates */

   cross[0] = a[1] * b[2] - a[2] * b[1];
   cross[1] = a[2] * b[0] - a[0] * b[2];
   cross[2] = a[0] * b[1] - a[1] * b[0];
}

/**
 * computes the great circle distance in rad for two points given in xyz coordinates
 * taken from http://johnblackburne.blogspot.de/2012/05/angle-between-two-3d-vectors.html
 * @param[in] a point coordinates of point a
 * @param[in] b point coordinates of point b
 * @return great circle distance in rad between both points
 */
static inline double get_vector_angle(double a[3], double b[3]) {

   double dot_product = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];

   // the acos most accurate in the range [-M_SQRT1_2;M_SQRT1_2]
   if (fabs(dot_product) > M_SQRT1_2) {

      double cross_ab[3];

      crossproduct_ld(a, b, cross_ab);

      double asin_tmp = asin(sqrt(cross_ab[0]*cross_ab[0]+
                                  cross_ab[1]*cross_ab[1]+
                                  cross_ab[2]*cross_ab[2]));

      if (dot_product > 0.0) // if the angle is smaller than (PI / 2)
         return MAX(asin_tmp,0.0);
      else
         return MIN(M_PI - asin_tmp, M_PI);

   } else {

      return acos(dot_product);
   }

   /*
   // this solution is simpler, but has a much worse performance
   double cross[3], dot, cross_abs;

   crossproduct_ld(a, b, cross);
   dot = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
   cross_abs = sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);

   return fabs(atan2(cross_abs, dot));
   */
}

static inline struct sin_cos_angle sin_cos_angle_new(double sin, double cos) {

  struct sin_cos_angle angle;

  if (fabs(sin) <= 1.0) angle.sin = sin;
  else if (sin < 0.0) angle.sin = -1.0;
  else angle.sin = 1.0;

  if (fabs(cos) <= 1.0) angle.cos = cos;
  else if (cos < 0.0) angle.cos = -1.0;
  else angle.cos = 1.0;

  return angle;
}

static inline struct sin_cos_angle get_vector_angle_2(
  double a[3], double b[3]) {

  double cross_ab[3];
  crossproduct_ld(a, b, cross_ab);

  return sin_cos_angle_new(sqrt(cross_ab[0]*cross_ab[0] +
                                cross_ab[1]*cross_ab[1] +
                                cross_ab[2]*cross_ab[2]),
                           a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

// works for angles in the range [0;2*PI[
static inline int compare_angles(
  struct sin_cos_angle a, struct sin_cos_angle b) {

  // there are 0 section:
  // 0: 0      <= angle < PI/4
  // 1: PI/4   <= angle < 3*PI/4
  // 2: 3*PI/4 <= angle < 5*PI/4
  // 3: 5*PI/4 <= angle < 7*PI/4
  // 4: 7*PI/4 <= angle < 2*PI
  int t_a = fabs(a.cos) <= M_SQRT1_2;
  int t_b = fabs(b.cos) <= M_SQRT1_2;
  int a_section = t_a | ((((a.sin < 0.0) & t_a) |
                          ((a.cos < 0.0) & (fabs(a.sin) < M_SQRT1_2))) << 1);
  int b_section = t_b | ((((b.sin < 0.0) & t_b) |
                          ((b.cos < 0.0) & (fabs(b.sin) < M_SQRT1_2))) << 1);
  if (!a_section) a_section = (a.sin < 0.0) << 2;
  if (!b_section) b_section = (b.sin < 0.0) << 2;

  if (a_section != b_section)
    return (a_section > b_section) - (a_section < b_section);

  int ret;

  switch (a_section) {
    case(0):
    case(4):
    default:
      if ((a.sin < 0.0) == (b.sin < 0.0)) {
        ret = (a.sin > b.sin) - (a.sin < b.sin);
        if (ret) return ret;
        else {
          ret = (a.cos < b.cos) - (a.cos > b.cos);
          if (a.sin >= 0.0) return ret;
          else return -ret;
        }
      } else {
        return (a.sin < 0.0) - (b.sin < 0.0);
      }
    case(1):
      ret = (a.cos < b.cos) - (a.cos > b.cos);
      if (ret) return ret;
      else {
        ret = (a.sin > b.sin) - (a.sin < b.sin);
        if (a.cos >= 0.0) return ret;
        else return -ret;
      }
    case(2):
      ret = (a.sin < b.sin) - (a.sin > b.sin);
      if (ret) return ret;
      else {
        ret = (a.cos < b.cos) - (a.cos > b.cos);
        if (a.sin >= 0.0) return ret;
        else return -ret;
      }
    case(3):
      ret = (a.cos > b.cos) - (a.cos < b.cos);
      if (ret) return ret;
      else {
        ret = (a.sin < b.sin) - (a.sin > b.sin);
        if (a.cos <= 0.0) return ret;
        else return -ret;
      }
  }
}

//! computes (a+b), if (a+b) >= 2*PI, the result is (a+b-2*PI)
//! a and be have to be in the range [0;2*PI[
static inline struct sin_cos_angle sum_angles_no_check(
  struct sin_cos_angle a, struct sin_cos_angle b) {

  return sin_cos_angle_new(a.sin * b.cos + a.cos * b.sin,
                           a.cos * b.cos - a.sin * b.sin);
}

//! computes (a+b), if (a+b) >= 2*PI, the result is (a+b-2*PI)
//! a and be have to be in the range [0;2*PI[
//! \returns 1 if (a+b) is >= 2*PI
static inline int sum_angles(
  struct sin_cos_angle a, struct sin_cos_angle b,
  struct sin_cos_angle * restrict sum) {

  struct sin_cos_angle sum_ = sum_angles_no_check(a, b);

  *sum = sum_;

  // if a or b is smaller than the result
  return (compare_angles(sum_, a) < 0) || (compare_angles(sum_, b) < 0);
}

//! computes (a-b), if (a-b) < 0, the result is (a+b+2*PI)
//! a and be have to be in the range [0;2*PI[
static inline struct sin_cos_angle sub_angles_no_check(
  struct sin_cos_angle a, struct sin_cos_angle b) {

  return sin_cos_angle_new(a.sin * b.cos - a.cos * b.sin,
                           a.cos * b.cos + a.sin * b.sin);
}

//! computes (a-b), if (a-b) < 0, the result is (a+b+2*PI)
//! a and be have to be in the range [0;2*PI[
//! \returns 1 if (a-b) < 0
static inline int sub_angles(
  struct sin_cos_angle a, struct sin_cos_angle b,
  struct sin_cos_angle * restrict sub) {

  *sub = sub_angles_no_check(a, b);

  return compare_angles(a, b) < 0;
}

//! return angles in the range of [0;2PI[
static inline double compute_angle(struct sin_cos_angle angle) {

  // the acos and asin are most accurate in the range [-M_SQRT1_2;M_SQRT1_2]
  if (angle.cos > M_SQRT1_2) {

    double angle_ = asin(angle.sin);

    if (angle_ < 0.0) angle_ += 2.0 * M_PI;
    return angle_;

  } else if (angle.cos > - M_SQRT1_2) {

    double angle_ = acos(angle.cos);

    if (angle.sin > 0.0) return angle_;
    else return 2.0 * M_PI - angle_;

  } else {
    return M_PI - asin(angle.sin);
  }
}

//! computes angle / 2
//! The basic idea is to imagin angle being a vector v with (cos(x); sin(x)) and
//! w being a vector with (1.0; 0.0). The normalised sum of v and w gives us v_h
//! with (cos(x/2); sin(x/2)).
//! For quadrants 2, 3, and 4 we have to apply some little tricks for higher
//! accuracy.
static inline struct sin_cos_angle half_angle(struct sin_cos_angle angle) {

  double x = (1.0 + fabs(angle.cos));

  double scale = 1.0 / sqrt(x * x + angle.sin * angle.sin);

  // first or fourth quadrant
  if (angle.cos >= 0) {
    scale = copysign(scale, angle.sin);
    return sin_cos_angle_new(angle.sin * scale, x * scale);

  // second and third quadrant
  } else return sin_cos_angle_new(x * scale, angle.sin * scale);
}

//! computes angle / 4
//! I derived it based on \ref half_angle
static inline struct sin_cos_angle quarter_angle(struct sin_cos_angle angle) {

  double tan = fabs(angle.sin) / (1.0 + fabs(angle.cos));
  double one_plus_sq_tan = 1.0 + tan * tan;
  double sqrt_one_plus_sq_tan = sqrt(one_plus_sq_tan);

  double a = 1.0;
  double b = tan;

  // second and third quadrant
  if (angle.cos < 0.0) {
    a = tan;
    b = 1.0;
  }

  double scale = M_SQRT1_2 / sqrt(one_plus_sq_tan + a * sqrt_one_plus_sq_tan);
  double x = b * scale;
  double y = (a + sqrt_one_plus_sq_tan) * scale;

  // first and second quadrant
  if (angle.sin >= 0.0) return sin_cos_angle_new(x, y);
  // third and fourth quadrant
  else return sin_cos_angle_new(y, x);
}

/**
 * determines whether two given points are (nearly) identically
 * @param[in] a point coordinates of point a
 * @param[in] b point coordinates of point b
 * @return true if both points are (nearly) identically
 */
static inline int points_are_identically(double * a, double * b) {

   double dot_product = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];

   // if the angle is bigger than ~0.81 degree
   // (in this range the acos is still rather accurate)
   if (dot_product < 0.9999) // (acos(0.9999) = ~0.81 degree)

      return 0;

   else { // both points are close to each other -> use cross product for higher
          // accuracy

      double cross_ab[3];

      crossproduct_ld(a, b, cross_ab);

      // for very small angles: asin(alpha) = ~alpha   (alpha in rad)
      return sqrt(cross_ab[0]*cross_ab[0] +
                  cross_ab[1]*cross_ab[1] +
                  cross_ab[2]*cross_ab[2]) < yac_angle_tol;
   }
}

/** normalises a vector */
static inline void normalise_vector(double v[]) {

   double norm = 1.0 / sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

   v[0] *= norm;
   v[1] *= norm;
   v[2] *= norm;
}

/** rotate vector v_in around the given axis by a given angle
 * @param[in] axis axis around which v_in is to rotated
 * @param[in] angle rotation angle
 * @param[in] v_in vector to be rotated
 * @param[out] v_out rotated vector
 */
static inline void rotate_vector2(
  double axis[], struct sin_cos_angle angle, double v_in[], double v_out[]) {

  // using Rodrigues' rotation formula
  // v_out = v_in * cos(angle) +
  //         (axis x v_in) * sin(angle) +
  //         axis * (axis * v_in) * (1 - cos(angle))

  double cross_axis_v_in[3];
  crossproduct_d(axis, v_in, cross_axis_v_in);

  double dot_axis_v_in = axis[0]*v_in[0] + axis[1]*v_in[1] + axis[2]*v_in[2];
  double temp = dot_axis_v_in * (1.0 - angle.cos);

  v_out[0] =
    v_in[0] * angle.cos + cross_axis_v_in[0] * angle.sin + axis[0] * temp;
  v_out[1] =
    v_in[1] * angle.cos + cross_axis_v_in[1] * angle.sin + axis[1] * temp;
  v_out[2] =
    v_in[2] * angle.cos + cross_axis_v_in[2] * angle.sin + axis[2] * temp;
}

/** rotate vector v_in around the given axis by a given angle
 * @param[in] axis axis around which v_in is to rotated
 * @param[in] angle rotation angle
 * @param[in] v_in vector to be rotated
 * @param[out] v_out rotated vector
 */
static inline void rotate_vector(
  double axis[], double angle, double v_in[], double v_out[]) {

  rotate_vector2(
    axis, sin_cos_angle_new(sin(angle), cos(angle)), v_in, v_out);
}

/**
 * computes the great circle distance in rad for two points given in lon-lat coordinates
 * @param[in] a point coordinates of point a (in rad)
 * @param[in] b point coordinates of point b (in rad)
 * @return great circle distance in rad between both points
 */
double yac_get_point_angle(struct point * a, struct point * b);

/**
 * determines whether two edges intersect
 * @param[in] edge_a first edge
 * @param[in] a      3d coordinate of first point of first edge
 * @param[in] b      3d coordinate of second point of first edge
 * @param[in] edge_b second edge
 * @param[in] c      3d coordinate of first point of second edge
 * @param[in] d      3d coordinate of second point of second edge
 * @return 0 if edges do not intersect\n
 *         1 if edges intersect
 */
int yac_do_intersect (struct edge edge_a, double a[3], double b[3],
                      struct edge edge_b, double c[3], double d[3]);

/**
 * splits given cell into triangles
 * @param[in] cell cell to be triangulated
 * @param[in] start_corner start algorithm at corner with this index
 *                         (0 <= start_corner < n; with n being number of cell
 *                          corners)
 * @param[out] triangles triangles that are the result of the triangulation
 * @remark the user needs to provide (n-2) initialised grid_cells for the
 *         argument triangles (n being the number of cell corners of the cell)
 * @remark this routine currently only works for convex grid cells
 * @remark the algorithm ignores that actual edge types for cells with more than
 *         3 corners and assumes them to be great circle edges
 */
void yac_triangulate_cell(
  struct grid_cell cell, unsigned start_corner, struct grid_cell * triangles);

/**
 * splits given indices of a cell into triangles
 * @param[in] cell_indices cell indices to be triangulated
 * @param[in] num_corners  number of corners of the cell
 * @param[in] start_corner start algorithm at corner with this index
 *                         (0 <= start_corner < n; with n being number of cell
 *                          corners)
 * @param[out] triangle_indices triangle indices that are the result of the
 *                              triangulation
 * @remark the user needs to provide provide an array of size 3*(num_corners-2)
 *         for the argument triangle_indices
 */
void yac_triangulate_cell_indices(
  unsigned const * cell_indices, unsigned num_corners, unsigned start_corner,
  unsigned * triangle_indices);

/**
 * computes a spherical triangle that contains all given vertices,
 * the algorithm tries to minimise the size of the triangle
 * @param[in] vertices vertices that are supposed to be within the triangle
 * @param[in] num_vertices number of vertices in vertices
 * @param[out] triangle bounding triangle
 * @param[in] num_tests with increasing number of tests computation time and
 *                      chance of smaller bounding triangle increases
 * @remark num_tests has to be > 0
 */
void yac_compute_bnd_triangle(
  double * vertices, size_t num_vertices, double * triangle, size_t num_tests);

#endif // GEOMETRY_H

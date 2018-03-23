/**
 * @file sphere_part.h
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

/** \example test_sphere_part.c
* This contains a test of the sphere_part grid search algorithm.
*/


#ifndef SPHERE_PART_H
#define SPHERE_PART_H

#include "grid.h"
#include "interval_tree.h"
#ifdef YAC
#include "grid_search.h"
#endif

/**
 * \file sphere_part.h
 * \brief algorithm for searching cells and points on a grid
 *
 * \ref yac_sphere_part_search_new generates a tree structure, which makes it
 * easy to look for cells and points. A documentation of the respective
 * algorithm can be found at \ref sphere_part_docu.
 */
 
/**
 * \page sphere_part_docu Sphere Partitioning Algorithm
 *
 * The following describes how the Sphere Partitioning algorithm generates a tree data structure for a given set of polygons on a sphere.
 * This tree structure allows to easily search all cells in the given data set that overlaps with another given cell or point (set of cells or points).
 *
 * \section sphere_part_tree_gen Generation of tree structure
 * to partition a set of polygons on the sphere the following data structure can be constructed:
 * - c = center of sphere
 * - S = set of all points
 *
 * \subsection recurspart_sec Description for a routine that generates the tree
 * recurspart(P, v, t)
 * -# Set p to balance point of P
 * -# Compute great circle L with base plane collinear to |c p| and orthogonal to v
 * -# Partition data
 *    -# Set I to subset of S intersecting L
 *    -# Set T to subset of S/I with positive scalar product to normal vector of L
 *       - Ti = recurspart(T, norm(L)) if |T| > threshold t else Ti = list(T)
 *    -# Set U to subset of S/I with negative scalar product to normal vector of L
 *       - Ui = recurspart(U, norm(L)) if |U| > threshold t else Ui = list(U)
 * -# return node(list(I), Ti, Ui, L, alpha=angle over L containing I)
 *
 * \subsection recuspart_observ_sec Observations
 * - I, T, U and form a partition of S
 * - first L is orthogonal to equator, i.e. a longitude circle
 * - initial call: recurspart(S, (0,0,1), t)
 *
 * \section sphere_part_search Searching in the tree structure
 *
 * \subsection search_sec Description for a routine that searches for a list of cells
 * search(n, p)
 * -# if n is leaf
 *    -# P = search_list(n, p)
 * -# else
 *    -# P = {}
 *    -# if p in n.alpha
 *       - P = P united search_list(n.I, p)
 *    -# if p * norm(n.L) > 0
 *       - P = P united search(n.Ti, p)
 *    -# else if p * norm(n.L) < 0
 *       - P = P united search(n.Ui, p)
 * -# return P
 *
 * \subsection search_observ_sec Observations
 * - returns list of matching polygons
 */

#ifdef YAC
struct grid_search * yac_sphere_part_search_new (struct grid * grid_data);
#endif
struct point_sphere_part_search;

struct point_sphere_part_search * yac_point_sphere_part_search_new (
  size_t num_points, double * coordinates_xyz);

void yac_delete_point_sphere_part_search(
  struct point_sphere_part_search * search);

/**
 * This routine does a nearest neighbour search between the points provided to
 * this routine and the matching yac_point_sphere_part_search_new call.
 */
void yac_point_sphere_part_search_NN(struct point_sphere_part_search * search,
                                     size_t num_points,
                                     double * coordinates_xyz,
                                     double * cos_angles,
                                     double ** result_coordinates_xyz,
                                     size_t * result_coordinates_xyz_array_size,
                                     unsigned ** local_point_ids,
                                     size_t * local_point_ids_array_size,
                                     size_t * num_local_point_ids);

/**
 * This routine does a n nearest neighbour search between the points provided to
 * this routine and the matching yac_point_sphere_part_search_new call.
 */
void yac_point_sphere_part_search_NNN(struct point_sphere_part_search * search,
                                      size_t num_points,
                                      double * coordinates_xyz, unsigned n,
                                      double ** cos_angles,
                                      size_t * cos_angles_array_size,
                                      double ** result_coordinates_xyz,
                                      size_t * result_coordinates_xyz_array_size,
                                      unsigned ** local_point_ids,
                                      size_t * local_point_ids_array_size,
                                      size_t * num_local_point_ids);

/**
 * This routine returns true if the provided point_sphere_part_search contains
 * a point that is within the provided bounding circle.
 */
#ifdef YAC
int yac_point_sphere_part_search_bnd_circle_contains_points(
  struct point_sphere_part_search * search, struct bounding_circle circle);
#endif
#endif // SPHERE_PART_H

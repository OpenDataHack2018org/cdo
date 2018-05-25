/**
 * @file grid_search.h
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
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

#ifndef GRID_SEARCH_H
#define GRID_SEARCH_H

#include "grid.h"
#include "dep_list.h"
#include "points.h"
#include "geometry.h"

struct grid_search;

struct grid_search_vtable {

   void (*do_cell_search)(struct grid_search * search, struct grid * grid_data,
                          struct dep_list * tgt_to_src_cells);
   void (*do_cell_search_single)(struct grid_search * search,
                                 struct grid_cell grid_cell,
                                 unsigned * n_cells, unsigned * cells_size,
                                 unsigned ** cells);
   void (*do_point_search_c) (struct grid_search * search, struct grid * grid_data,
                              struct dep_list * tgt_to_src_cells);
   void (*do_point_search_c2) (struct grid_search * search, double (*coordinates_xyz)[3],
                               unsigned num_points, struct dep_list * tgt_to_src_cells);
   void (*do_point_search_c3) (struct grid_search * search, double (*coordinates_xyz)[3],
                               unsigned num_points, struct dep_list * tgt_to_src_cells,
                               struct points * points);
   void (*do_point_search_p) (struct grid_search * search, struct grid * grid_data,
                              struct dep_list * target_to_src_points);
   void (*do_point_search_p2) (struct grid_search * search, double (*coordinates_xyz)[3],
                               unsigned num_points, struct dep_list * target_to_src_points);
   void (*do_point_search_p3) (struct grid_search * search, double (*coordinates)[3],
                               unsigned num_points, struct dep_list * target_to_src_points,
                               struct points * points);
   void (*do_point_search_p4) (struct grid_search * search, double coordinate_xyz[3],
                               unsigned * n_points, unsigned * points_size,
                               unsigned ** points);
   void (*do_bnd_circle_search) (struct grid_search * search,
                                 struct bounding_circle * bnd_circles,
                                 unsigned n, struct dep_list * bnd_to_cells);
   void (*delete_grid_search) (struct grid_search * search);
};

struct grid_search {

   struct grid_search_vtable *vtable;
};

/** \example test_cell_search.c
 * These are some examples on how to use \ref yac_do_cell_search.
 */

/**
 * does a cell search \n determines for every target cell all overlapping source cells
 * @param[in]  search           grid search object
 * @param[in]  grid_data        grid whose cells are to be searched
 * @param[out] tgt_to_src_cells dependency containing the mapping of target to
 *                              source cells
 */
void yac_do_cell_search (struct grid_search * search, struct grid * grid_data,
                         struct dep_list * tgt_to_src_cells);

/**
 * does a cell search for a single cell \n determines for the target cell all
 * overlapping source cells
 * @param[in]        search     grid search object
 * @param[in]        grid_cell  grid cell for which the search is to conducted
 * @param[out]       n_cells    number of overlapping source cells
 * @param[in,out]    cells_size number of elements that can be stored in *cells
 * @param[in,out]    cells      indices of overlapping source cells
 * @remarks the user is responsible to free the memory associated to *cells
 * @remarks *cells has to be either NULL or a pointer that was returned by
 *          C-standard memory allocation routine \n in case *cells is not NULL
 *          cells_size must contain the number of elements that fit into *cells
 */
void yac_do_cell_search_single (struct grid_search * search,
                                struct grid_cell grid_cell,
                                unsigned * n_cells, unsigned * cells_size,
                                unsigned ** cells);

/** \example test_point_search.c
 * This contains examples on how to use the point search.
 */

/**
 * does a point search \n searches for source cells that match the target corners
 * @param[in]  search           grid search object
 * @param[in]  grid_data        grid whose corners are to be searched
 * @param[out] tgt_to_src_cells dependency list that contains for every target
 *                              grid corner the respective source cell index
 */
void yac_do_point_search_c (struct grid_search * search, struct grid * grid_data,
                            struct dep_list * tgt_to_src_cells);

/**
 * does a point search \n searches for every target point the source cell into
 * which the point falls\n returns a dependency list that contains for every
 * target point the respective source cell index
 * @param[in]  search           grid search object
 * @param[in]  coordinates_xyz  cartesian coordinates of the target points
 * @param[in]  num_points       number of target points
 * @param[out] tgt_to_src_cells dependency list that contains for every target
 *                              point the respective source cell index
 */
void yac_do_point_search_c2 (struct grid_search * search,
                             double (*coordinates_xyz)[3], unsigned num_points,
                             struct dep_list * tgt_to_src_cells);

/**
 * does a point search \n searches for every target point the source cell into
 * which the point falls\n returns a dependency list that contains for every
 * target point the respective source cell index
 * @param[in]  search           grid search object
 * @param[in]  coordinates_xyz  cartesian coordinates of the target points
 * @param[in]  num_points       number of target points
 * @param[out] tgt_to_src_cells dependency list that contains for every target
 *                              point the respective source cell index
 * @param[in]  points               point set to be used for the search
 */
void yac_do_point_search_c3 (struct grid_search * search,
                             double (*coordinates_xyz)[3], unsigned num_points,
                             struct dep_list * tgt_to_src_cells,
                             struct points * points);
/**
 * does a point search \n searches for source cells that matches the target corners
 * @param[in]  search               grid search object
 * @param[in]  grid_data            grid whose corners are to be searched
 * @param[out] target_to_src_points dependency list that contains for every target
 *                                  grid corner the source grid corner indices of
 *                                  the matching cell
 */
void yac_do_point_search_p (struct grid_search * search, struct grid * grid_data,
                            struct dep_list * target_to_src_points);

/**
 * does a point search\n searches for every target point the source cell into
 * which the point falls\n  returns a dependency list that contains for every
 * target point the indices of the corners of the respective source cell
 * @param[in]  search               grid search object
 * @param[in]  coordinates_xyz      cartesian coordinates of the target points
 * @param[in]  num_points           number of target points
 * @param[out] target_to_src_points dependency list that contains for every target point
 *                                  the source points, that build the source cell into
 *                                  which the respective target point falls
 */
void yac_do_point_search_p2 (struct grid_search * search,
                             double (*coordinates_xyz)[3], unsigned num_points,
                             struct dep_list * target_to_src_points);
/**
 * does a point search\n searches for every target point the source cell into
 * which the point falls\n (source cells are built from the provided point set)\n
 * the results are returned in form of a dependency list that contains for every
 * target point the indices of the corners of the respective source cell
 * @param[in]  search               grid search object
 * @param[in]  coordinates_xyz      cartesian coordinates of the target points
 * @param[in]  num_points           number of target points
 * @param[out] target_to_src_points dependency list that contains for every target point
 *                                  the source points, that build the source cell into
 *                                  which the respective target point falls
 * @param[in]  points               point set to be used for the search
 * @see global_search_new
 */
void yac_do_point_search_p3 (struct grid_search * search,
                             double (*coordinates_xyz)[3], unsigned num_points,
                             struct dep_list * target_to_src_points,
                             struct points * points);

/**
 * does a point search\n searches for the target point the matching source cell\n
 * returns a list that contains for the target point the indices of the corners
 * of the respective source cell
 * @param[in]        search         grid search object
 * @param[in]        coordinate_xyz cartesian coordinate of the target point
 * @param[out]       n_points       number of source points
 * @param[in,out]    points_size    number of elements that can be stored in *points
 * @param[in,out]    points         indices of corners of source cell into which the
 * @remarks the user is responsible to free the memory associated to *cells
 * @remarks *cells has to be either NULL or a pointer that was returned by
 *          C-standard memory allocation routine \n in case *cells is not NULL
 *          cells_size must contain the number of elements that fit into *cells
 */
void yac_do_point_search_p4 (struct grid_search * search,
                             double coordinate_xyz[3], unsigned * n_points,
                             unsigned * points_size, unsigned ** points);

/**
 * does a bounding circle search\n
 * searches for all cells that overlap with the provided bounding circles\n
 * returns a list that contains for each bounding circle a list of cells that
 * overlap with it
 * of the respective source cell
 * @param[in]        search       grid search object
 * @param[in]        bnd_circles  bounding circles
 * @param[in]        n            number of entries in bnd_circles
 * @param[in]        bnd_to_cells results
 */
void yac_do_bnd_circle_search (struct grid_search * search,
                               struct bounding_circle * bnd_circles,
                               unsigned n, struct dep_list * bnd_to_cells);

/**
 * frees all memory associated with a grid search object
 * @param[in,out] search grid search object
 */
void yac_delete_grid_search(struct grid_search * search);

#endif

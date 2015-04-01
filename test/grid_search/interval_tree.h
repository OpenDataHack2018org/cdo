/**
 * @file interval_tree.h
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
 *             Thomas Jahns <jahns@dkrz.de>
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

/** \example test_interval_tree.c
 * This contains a test of the interval_tree.
 */

#ifndef INTERVAL_TREE_H
#define INTERVAL_TREE_H

#include <stdlib.h>

struct interval
{
  double left, right;
};

static inline int
overlap_test(struct interval a, struct interval b)
{
  return (a.left <= b.left && a.right >= b.left) ||
    (a.left > b.left && a.left <= b.right);
}

struct interval_node
{
  struct interval range;
  double max;
  unsigned value;
};

void
yac_generate_interval_tree(struct interval_node intervals[], size_t num_nodes);

struct overlaps
{
  size_t num_overlaps, a_size;
  size_t *overlap_iv;
};

void
yac_search_interval_tree(struct interval_node tree[], size_t num_nodes,
                         struct interval query, struct overlaps *overlaps);

#endif

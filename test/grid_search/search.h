/**
 * @file search.h
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

#include "dep_list.h"
#include "grid.h"
//#include "component.h"
#include "geometry.h"
#include "points.h"

/** \example test_bisection.c
 * These are some examples on how to use \ref bisection_search.
 */

/**
 * does a 1D bisection search for double precision data
 * @param[in]  data        data to be searched
 * @param[in]  data_length number of elements in data
 * @param[in]  axis_data   1D data array with continuous values in which
 *                         the elements of data are to be searched
 * @param[in]  axis_length number of points in axis_data
 * @param[out] position    position of elements of data in axis_data
 * @param[out] found       0 for elements of data that were not found in axis_data
 * @param[in]  period      if > 0 it indicates the period for cyclic behaviour
 *                         of elements in data
 */
void bisection_search (double const * data, unsigned data_length,
                       double const * axis_data, unsigned axis_length,
                       int * position, int * found, double period);

/**
 * does a 1D bisection search for integer data
 * @param[in]  data        data to be searched
 * @param[in]  data_length number of elements in data
 * @param[in]  axis_data   1D data array with continuous values in
 *                         which the elements of data are to be searched
 * @param[in]  axis_length number of points in axis_data
 * @param[out] position    position of elements of data in axis_data
 * @param[out] found       0 for elements of data that were not found in axis_data
 */
void bisection_search_int (int const * data, unsigned data_length,
                           int const * axis_data, unsigned axis_length,
                           int * position, int * found);

/**
 * End of the definition phase, invocation of the search
 *
 * @param[in] nbr_comps  number of components defined on the calling process
 * @param[in] comp_ids   list of component IDs as they are returned by yac_cdef_comp
 * @param[in] nbr_fields number of fields defined on the calling process
 * @param[in] field_ids  list of field IDs as they are returned by yac_cdef_field
 */
int start_search (int nbr_comps, int * comp_ids, int nbr_fields, int * field_ids );

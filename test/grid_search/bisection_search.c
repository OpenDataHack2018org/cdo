/**
 * @file bisection_search.c
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *                                 Rene Redler <rene.redler@mpimet.mpg.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler <rene.redler@mpimet.mpg.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/** \file bisection_search.c
 *  \brief bisection search along 1d arrays
 **/

/** bisection search for double arrays
 *
 * @param[in]  data         double data array to be searched
 * @param[in]  data_length  length of data array
 * @param[in]  axis_data    double axis against which data is compared
 * @param[in]  axis_length  size of axis array
 * @param[out] position     location of data element on axis 
 * @param[out] found        1 if something was found, else 0
 * @param[out] period       ???
 *
 **/

void yac_bisection_search (double const * data, unsigned data_length, 
                           double const * axis_data, unsigned axis_length,
                           int * restrict position, int * restrict found,
                           double period) {

   int upper, lower, middle;
   unsigned i;

   double lower_bound, upper_bound, curr_item;

   int ascending_order;

   ascending_order = axis_data[0] < axis_data[axis_length-1];

   if (ascending_order) {
      lower_bound = axis_data[0];
      upper_bound = axis_data[axis_length-1];
   } else {
      upper_bound = axis_data[0];
      lower_bound = axis_data[axis_length-1];
   }

   for (i = 0; i < data_length; ++i) {

      lower = -1, upper = axis_length;

      curr_item = data[i];

      if (period > 0.0) {

         while (curr_item < lower_bound) curr_item += period;
         while (curr_item > upper_bound) curr_item -= period;
      }

      /*---------------------------------------------------------------------
       * test upper and lower bounds
       *---------------------------------------------------------------------*/

      if ( curr_item <= lower_bound || curr_item >= upper_bound ) {

         if ( curr_item == axis_data[0] ) {
            found[i]    = 1;
            position[i] = 0;

         } else if ( curr_item == axis_data[axis_length-1] ) {

            found[i]    = 1;
            position[i] = axis_length-2;

         } else if ( (curr_item < lower_bound) != ascending_order ) {

            found[i]    = 0;
            position[i] = upper;

         } else {

            found[i]    = 0;
            position[i] = lower;
         }

      } else {

         /*---------------------------------------------------------------------
          * bisectional search
          *---------------------------------------------------------------------*/

         while ( upper - lower > 1 ) {

            middle = ( upper + lower ) / 2;

            if ( ascending_order != (curr_item > axis_data[middle]) )
               upper = middle;
            else
               lower = middle;
         }

         found[i] = 1;
         position[i] = lower;

      } // ( curr_item <= lower_bound || curr_item >= upper_bound )

   } // i
}


/** exact location on integer array using bisectional search
 *
 * @param[in]  data         integer data array to be searched
 * @param[in]  data_length  length of data array
 * @param[in]  axis_data    integer axis against which data is compared
 * @param[in]  axis_length  size of axis array
 * @param[out] position     location of data element on axis 
 * @param[out] found        1 if something was found, else 0
 *
 **/

void yac_bisection_search_int (int const * data, unsigned data_length,
                               int const * axis_data, unsigned axis_length,
                               int * restrict position, int * restrict found) {
  int upper, lower, middle;
  unsigned i;

  double lower_bound, upper_bound, curr_item;

  int ascending_order;

  ascending_order = axis_data[0] < axis_data[axis_length-1];

  if (ascending_order) {
    lower_bound = axis_data[0];
    upper_bound = axis_data[axis_length-1];
  } else {
    upper_bound = axis_data[0];
    lower_bound = axis_data[axis_length-1];
  }

  for (i = 0; i < data_length; ++i) {

    lower = -1, upper = axis_length;

    curr_item = data[i];

    /* ---------------------------------------------------------------------
       test upper and lower bounds
       ---------------------------------------------------------------------*/

    if ( curr_item <= lower_bound || curr_item >= upper_bound ) {

      if ( curr_item == axis_data[0] ) {
	found[i]    = 1;
	position[i] = 0;

      } else if ( curr_item == axis_data[axis_length-1] ) {

	found[i]    = 1;
	position[i] = axis_length-1;

      } else if ( (curr_item < lower_bound) != ascending_order ) {

	found[i]    = 0;
	position[i] = upper;

      } else {

	found[i]    = 0;
	position[i] = lower;
      }

    } else {

      /* ---------------------------------------------------------------------
	 bisectional search
	 ---------------------------------------------------------------------*/

      while ( upper - lower > 1 ) {

	middle = ( upper + lower ) / 2;

	if ( ascending_order != (curr_item > axis_data[middle]) )
	  upper = middle;
	else
	  lower = middle;
      }

      found[i] = 1;
      if (ascending_order) {
	position[i] = upper;
      }
      else {
	position[i] = lower;
      }
    } // ( curr_item <= lower_bound || curr_item >= upper_bound )

  } // for (i = 0; i < data_length; ++i)
}

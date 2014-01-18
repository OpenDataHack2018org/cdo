/**
 * @file utils.h
 * @brief Utlity functions
 *
 * Small general utility functions:
 *  - pointer-id conversion
 *  - hash
 *  - sorting
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

#ifndef UTILS_H
#define UTILS_H

/**
 * gives a unique index for a given pointer
 * @param[in] pointer
 * @return unique value associated to pointer
 * @see unique_id_to_pointer
 */
unsigned pointer_to_unique_id(void * pointer);

/**
 * gives the pointer that is associated to the given id
 * @param[in] id unique index previously returned by pointer_to_unique_id
 * @return pointer the is associated to the given id \n NULL if the id is invalid
 */
void *   unique_id_to_pointer(unsigned id);

/**
 * frees all memory used for the pointer/unique_id conversion
 * \remarks this should only be called after the last call to 
 *          \ref pointer_to_unique_id and \ref unique_id_to_pointer, because afterwords
 *          \ref unique_id_to_pointer will not be able to return the respective pointers
 *          for previously valid unique ids
 */
void free_pointer_unique_lookup();

/**
 * prints a short error message and info from where it was called
 * followed by an exit. 
 */
void abort_message ( char * text, char * file, int line );

/** \example test_quicksort.c
 * This contains an example of how to use quicksort_index.
 */

void quicksort_index ( int * a, int n, int * idx);

/* --------------------------------------------------------------------

   Hash function

   This algorithm (k=33) was first reported by dan bernstein many
   years ago in comp.lang.c. another version of this algorithm (now
   favored by bernstein) uses xor: hash(i) = hash(i - 1) * 33 ^
   str[i]; the magic of number 33 (why it works better than many other
   constants, prime or not) has never been adequately explained.

   Source: http://www.cse.yorku.ca/~oz/hash.html 

   --------------------------------------------------------------------- */

unsigned int hash(const char *str);

/* =======================================================================
   Macros
   ======================================================================= */

#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define MIN(a,b) ((a) < (b) ? (a) : (b))

#define ASSERT(c) \
if (!(c)) {\
   fprintf(stderr, "### Assertion violation: %s in %s:%d\n",\
           #c, __FILE__, __LINE__);\
   abort ();\
}

#define ASSERT2(c, a, b) \
if (!(c)) {\
   fprintf(stderr, "### Assertion violation: %s (%s = %d, %s = %d) in %s:%d\n", #c, #a, a, #b, b, __FILE__, __LINE__);\
   abort();\
}

#endif // UTILS_H

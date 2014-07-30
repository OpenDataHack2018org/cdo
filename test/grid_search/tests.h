/**
 * @file tests.h
 *
 * @copyright Copyright  (C)  2013 DKRZ, MPI-M
 *
 * @version 1.0
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler  <rene.redler@mpimet.mpg.de>
 */
/*
 *
 * Keywords:
 *
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler  <rene.redler@mpimet.mpg.de>
 *
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
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>

extern unsigned err_count__;

#ifdef VERBOSE
#define PUT_ERR(string) err_count__++, fputs((string), stderr);
#else
#define PUT_ERR(string) err_count__++;
#endif

#define INC_ERR (err_count__++)

#define TEST_EXIT_CODE ((!err_count__)?EXIT_SUCCESS:EXIT_FAILURE)
#define EXIT_SKIP_TEST (77)


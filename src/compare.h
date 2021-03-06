/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2018 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#ifndef COMPARE_H
#define COMPARE_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cmath>
#include <string.h>

#ifdef __xlC__ /* performance problem on IBM */
#ifndef DBL_IS_NAN
#define DBL_IS_NAN(x) ((x) != (x))
#endif
#else
#ifndef DBL_IS_NAN
#if defined(HAVE_DECL_ISNAN) && !defined(DISABLE_ISNAN)
#define DBL_IS_NAN(x) (std::isnan(x))
#elif defined(FP_NAN)
#define DBL_IS_NAN(x) (std::fpclassify(x) == FP_NAN)
#else
#define DBL_IS_NAN(x) ((x) != (x))
#endif
#endif
#endif

#ifndef DBL_IS_EQUAL
/*#define DBL_IS_EQUAL(x,y) (!(x < y || y < x)) */
#define DBL_IS_EQUAL(x, y) (DBL_IS_NAN(x) || DBL_IS_NAN(y) ? (DBL_IS_NAN(x) && DBL_IS_NAN(y) ? 1 : 0) : !(x < y || y < x))
#endif

#ifndef IS_EQUAL
#define IS_NOT_EQUAL(x, y) (x < y || y < x)
#define IS_EQUAL(x, y) (!IS_NOT_EQUAL(x, y))
#endif

#define STR_IS_EQ(x, y) ((*x) == (*y) && strcmp(x, y) == 0)

#endif /* COMPARE_H */

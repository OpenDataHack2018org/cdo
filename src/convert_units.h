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
#ifndef __CONVERT_UNITS_H_
#define __CONVERT_UNITS_H_

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#if defined(HAVE_LIBUDUNITS2) && (defined(HAVE_UDUNITS2_H) || defined(HAVE_UDUNITS2_UDUNITS2_H))
#define HAVE_UDUNITS2
#endif

#ifdef HAVE_UDUNITS2
#ifdef HAVE_UDUNITS2_UDUNITS2_H
#include <udunits2/udunits2.h>
#else
#include <udunits2.h>
#endif

void cdoConvertFree(void *ut_converter);
void cdoConvertDestroy();
#endif

#include "cdo_int.h"

void cdoConvertUnits(void **ut_converter, bool *changeunits, char *units, char *units_old, const char *name);

#endif  // __CONVERT_UNITS_H_

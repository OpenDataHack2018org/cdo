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
#ifndef INTERPOL_H
#define INTERPOL_H

#include <stdio.h>

void intgridbil(field_type *field1, field_type *field2);
void intgridcon(field_type *field1, field_type *field2);
void interpolate(field_type *field1, field_type *field2);

void intgriddis(field_type *field1, field_type *field2, size_t nnn);

#endif

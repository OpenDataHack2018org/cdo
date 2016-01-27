/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2016 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include <stdio.h>
#include "field.h"

static
void fld_field_init(field_t *field, size_t nmiss, double missval, size_t ngp, double *restrict array)
{
  field_init(field);

  field->size    = ngp;
  field->nmiss   = nmiss;
  field->missval = missval;
  field->ptr     = array;
  field->weight  = NULL;
}


double fun_fldmin(size_t nmiss, double missval, size_t ngp, double *restrict array)
{
  field_t field;
  fld_field_init(&field, nmiss, missval, ngp, array);
  
  return fldmin(field);
}


double fun_fldmax(size_t nmiss, double missval, size_t ngp, double *restrict array)
{
  field_t field;
  fld_field_init(&field, nmiss, missval, ngp, array);
  
  return fldmax(field);
}

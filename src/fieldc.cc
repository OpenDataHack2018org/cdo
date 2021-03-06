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

#include "cdo_int.h"
#include <cdi.h>

void
farcfun(Field *field, double rconst, int function)
{
  switch (function)
    {
    case func_add: farcadd(field, rconst); break;
    case func_sub: farcsub(field, rconst); break;
    case func_mul: farcmul(field, rconst); break;
    case func_div: farcdiv(field, rconst); break;
    case func_mod: farmod(field, rconst); break;
    default: cdoAbort("%s: function %d not implemented!", __func__, function);
    }
}

void
farcmul(Field *field, double rconst)
{
  int nwpv = field->nwpv;
  int grid = field->grid;
  size_t nmiss = field->nmiss;
  double missval1 = field->missval;
  double missval2 = field->missval;
  double *array = field->ptr;

  if (nwpv != 2) nwpv = 1;

  size_t len = nwpv * gridInqSize(grid);

  if (nmiss > 0)
    {
      for (size_t i = 0; i < len; i++) array[i] = MULMN(array[i], rconst);
    }
  else
    {
      /*
#ifdef  _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
      */
      for (size_t i = 0; i < len; i++) array[i] *= rconst;
    }
}

void
farcdiv(Field *field, double rconst)
{
  int nwpv = field->nwpv;
  int grid = field->grid;
  size_t nmiss = field->nmiss;
  double missval1 = field->missval;
  double missval2 = field->missval;
  double *array = field->ptr;

  if (nwpv != 2) nwpv = 1;

  size_t len = nwpv * gridInqSize(grid);

  if (nmiss > 0 || IS_EQUAL(rconst, 0))
    {
      for (size_t i = 0; i < len; i++) array[i] = DIVMN(array[i], rconst);

      if (IS_EQUAL(rconst, 0)) field->nmiss = len;
    }
  else
    {
      for (size_t i = 0; i < len; i++) array[i] /= rconst;
    }
}

void
farcadd(Field *field, double rconst)
{
  int nwpv = field->nwpv;
  int grid = field->grid;
  size_t nmiss = field->nmiss;
  double missval1 = field->missval;
  double missval2 = field->missval;
  double *array = field->ptr;

  if (nwpv != 2) nwpv = 1;

  size_t len = nwpv * gridInqSize(grid);

  if (nmiss > 0)
    {
      for (size_t i = 0; i < len; i++) array[i] = ADDMN(array[i], rconst);
    }
  else
    {
      for (size_t i = 0; i < len; i++) array[i] += rconst;
    }
}

void
farcsub(Field *field, double rconst)
{
  farcadd(field, -rconst);
}

void
farinv(Field *field)
{
  int grid = field->grid;
  double missval1 = field->missval;
  double missval2 = field->missval;
  double *array = field->ptr;

  size_t len = gridInqSize(grid);

  for (size_t i = 0; i < len; i++) array[i] = DIVMN(1.0, array[i]);

  field->nmiss = arrayNumMV(len, array, missval1);
}

void
farround(Field *field)
{
  int grid = field->grid;
  double missval1 = field->missval;
  double *array = field->ptr;

  size_t len = gridInqSize(grid);

  for (size_t i = 0; i < len; i++)
    if (!DBL_IS_EQUAL(array[i], missval1)) array[i] = round(array[i]);
}

void
farmod(Field *field, double divisor)
{
  int grid = field->grid;
  double missval1 = field->missval;
  double *array = field->ptr;

  size_t len = gridInqSize(grid);

  for (size_t i = 0; i < len; i++)
    {
      array[i] = DBL_IS_EQUAL(array[i], missval1) ? missval1 : fmod(array[i], divisor);
    }
}

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

static void
farcmulcplx(Field *field, const double rconst[2])
{
  int nwpv = field->nwpv;
  int grid = field->grid;
  double missval1 = field->missval;
  double missval2 = field->missval;
  double *array = field->ptr;

  if (nwpv != 2) nwpv = 1;

  size_t gridsize = gridInqSize(grid);
  size_t len = nwpv * gridInqSize(grid);

  // z1 x z2 = (x1x2 - y1y2) + i(x1y2 + x2y1)
  for (size_t i = 0; i < gridsize; i++)
    {
      double a1r = array[2 * i];
      double a1i = array[2 * i + 1];
      array[2 * i] = SUBMN(MULMN(a1r, rconst[0]), MULMN(a1i, rconst[1]));
      array[2 * i + 1] = ADDMN(MULMN(a1r, rconst[1]), MULMN(a1i, rconst[0]));
    }
}

static void
farcdivcplx(Field *field, const double rconst[2])
{
  int nwpv = field->nwpv;
  int grid = field->grid;
  double missval1 = field->missval;
  double missval2 = field->missval;
  double *array = field->ptr;

  if (nwpv != 2) nwpv = 1;

  size_t gridsize = gridInqSize(grid);
  size_t len = nwpv * gridInqSize(grid);

  // z1 / z2 = (x1x2 + y1y2) / (x2x2 + y2y2) + i (y1x2 - x1y2) / (x2x2 + y2y2)
  for (size_t i = 0; i < gridsize; i++)
    {
      double a1r = array[2 * i];
      double a1i = array[2 * i + 1];
      double denominator = ADDMN(MULMN(rconst[0], rconst[0]), MULMN(rconst[1], rconst[1]));
      array[2 * i] = DIVMN(ADDMN(MULMN(a1r, rconst[0]), MULMN(a1i, rconst[1])), denominator);
      array[2 * i + 1] = DIVMN(SUBMN(MULMN(a1i, rconst[0]), MULMN(a1r, rconst[1])), denominator);
    }
}

static void
farcaddcplx(Field *field, const double rconst[2])
{
  int nwpv = field->nwpv;
  int grid = field->grid;
  double missval1 = field->missval;
  double missval2 = field->missval;
  double *array = field->ptr;

  if (nwpv != 2) nwpv = 1;

  size_t gridsize = gridInqSize(grid);
  size_t len = nwpv * gridInqSize(grid);

  for (size_t i = 0; i < gridsize; i++)
    {
      array[2 * i] = ADDMN(array[2 * i], rconst[0]);
      array[2 * i + 1] = ADDMN(array[2 * i + 1], rconst[1]);
    }
}

static void
farcsubcplx(Field *field, const double rconst[2])
{
  int nwpv = field->nwpv;
  int grid = field->grid;
  double missval1 = field->missval;
  double missval2 = field->missval;
  double *array = field->ptr;

  if (nwpv != 2) nwpv = 1;

  size_t gridsize = gridInqSize(grid);
  size_t len = nwpv * gridInqSize(grid);

  for (size_t i = 0; i < gridsize; i++)
    {
      array[2 * i] = SUBMN(array[2 * i], rconst[0]);
      array[2 * i + 1] = SUBMN(array[2 * i + 1], rconst[1]);
    }
}

void
farcfuncplx(Field *field, const double rconst[2], int function)
{
  switch (function)
    {
    case func_add: farcaddcplx(field, rconst); break;
    case func_sub: farcsubcplx(field, rconst); break;
    case func_mul: farcmulcplx(field, rconst); break;
    case func_div: farcdivcplx(field, rconst); break;
    default: cdoAbort("%s: function %d not implemented!", __func__, function);
    }
}

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

#include <cdi.h>

#include "cdo_int.h"
#include "array.h"

static void
faraddcplx(Field *field1, Field field2)
{
  int nwpv = field1->nwpv;
  int grid1 = field1->grid;
  int grid2 = field2.grid;
  double missval1 = field1->missval;
  double missval2 = field2.missval;
  double *restrict array1 = field1->ptr;
  const double *restrict array2 = field2.ptr;

  if (nwpv != 2) nwpv = 1;

  size_t gridsize = gridInqSize(grid1);
  size_t len = nwpv * gridsize;
  if (len != (nwpv * gridInqSize(grid2))) cdoAbort("Fields have different gridsize (%s)", __func__);

  for (size_t i = 0; i < gridsize; i++)
    {
      array1[2 * i] = ADDMN(array1[2 * i], array2[2 * i]);
      array1[2 * i + 1] = ADDMN(array1[2 * i + 1], array2[2 * i + 1]);
    }
}

static void
farsubcplx(Field *field1, Field field2)
{
  int nwpv = field1->nwpv;
  int grid1 = field1->grid;
  int grid2 = field2.grid;
  double missval1 = field1->missval;
  double missval2 = field2.missval;
  double *restrict array1 = field1->ptr;
  const double *restrict array2 = field2.ptr;

  if (nwpv != 2) nwpv = 1;

  size_t gridsize = gridInqSize(grid1);
  size_t len = nwpv * gridsize;
  if (len != (nwpv * gridInqSize(grid2))) cdoAbort("Fields have different gridsize (%s)", __func__);

  for (size_t i = 0; i < gridsize; i++)
    {
      array1[2 * i] = SUBMN(array1[2 * i], array2[2 * i]);
      array1[2 * i + 1] = SUBMN(array1[2 * i + 1], array2[2 * i + 1]);
    }
}

static void
farmulcplx(Field *field1, Field field2)
{
  int nwpv = field1->nwpv;
  int grid1 = field1->grid;
  int grid2 = field2.grid;
  double missval1 = field1->missval;
  double missval2 = field2.missval;
  double *restrict array1 = field1->ptr;
  const double *restrict array2 = field2.ptr;

  if (nwpv != 2) nwpv = 1;

  size_t gridsize = gridInqSize(grid1);
  size_t len = nwpv * gridsize;
  if (len != (nwpv * gridInqSize(grid2))) cdoAbort("Fields have different gridsize (%s)", __func__);

  // z1 x z2 = (x1x2 - y1y2) + i(x1y2 + x2y1)
  for (size_t i = 0; i < gridsize; i++)
    {
      double a1r = array1[2 * i];
      double a1i = array1[2 * i + 1];
      array1[2 * i] = SUBMN(MULMN(a1r, array2[2 * i]), MULMN(a1i, array2[2 * i + 1]));
      array1[2 * i + 1] = ADDMN(MULMN(a1r, array2[2 * i + 1]), MULMN(a1i, array2[2 * i]));
    }
}

static void
fardivcplx(Field *field1, Field field2)
{
  int nwpv = field1->nwpv;
  int grid1 = field1->grid;
  int grid2 = field2.grid;
  double missval1 = field1->missval;
  double missval2 = field2.missval;
  double *restrict array1 = field1->ptr;
  const double *restrict array2 = field2.ptr;

  if (nwpv != 2) nwpv = 1;

  size_t gridsize = gridInqSize(grid1);
  size_t len = nwpv * gridsize;
  if (len != (nwpv * gridInqSize(grid2))) cdoAbort("Fields have different gridsize (%s)", __func__);

  // z1 / z2 = (x1x2 + y1y2) / (x2x2 + y2y2) + i (y1x2 - x1y2) / (x2x2 + y2y2)
  for (size_t i = 0; i < gridsize; i++)
    {
      double a1r = array1[2 * i];
      double a1i = array1[2 * i + 1];
      double denominator = ADDMN(MULMN(array2[2 * i], array2[2 * i]), MULMN(array2[2 * i + 1], array2[2 * i + 1]));
      array1[2 * i] = DIVMN(ADDMN(MULMN(a1r, array2[2 * i]), MULMN(a1i, array2[2 * i + 1])), denominator);
      array1[2 * i + 1] = DIVMN(SUBMN(MULMN(a1i, array2[2 * i]), MULMN(a1r, array2[2 * i + 1])), denominator);
    }
}

void
farfuncplx(Field *field1, Field field2, int function)
{
  // clang-format off
  switch (function)
    {
    case func_add:     faraddcplx(field1, field2);   break;
    case func_sub:     farsubcplx(field1, field2);   break;
    case func_mul:     farmulcplx(field1, field2);   break;
    case func_div:     fardivcplx(field1, field2);   break;
    default: cdoAbort("%s: function %d not implemented!", __func__, function);
    }
  // clang-format on
}

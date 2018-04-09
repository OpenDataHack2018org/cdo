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
#include "grid_cell_search.h"

CellSearchMethod cellSearchMethod(CellSearchMethod::latbins);

void
setCellSearchMethod(const char *methodstr)
{
  // clang-format off
  if      (STR_IS_EQ(methodstr, "spherepart")) cellSearchMethod = CellSearchMethod::spherepart;
  else if (STR_IS_EQ(methodstr, "latbins"))    cellSearchMethod = CellSearchMethod::latbins;
  else cdoAbort("Grid cell search method %s not available!", methodstr);
  // clang-format on
}

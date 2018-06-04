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

#ifndef EXPR_FUN_H
#define EXPR_FUN_H

void fld_field_init(Field *field, size_t nmiss, double missval, size_t ngp, double *array, double *w);
double *fld_weights(int gridID, size_t ngp);
void vert_weights(int zaxisID, size_t nlev, std::vector<double> &weights);

#endif

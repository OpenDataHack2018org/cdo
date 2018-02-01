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
#ifndef _ARRAY_H
#define _ARRAY_H

const char *fpe_errstr(int fpeRaised);

int array_minmaxsum_val(size_t len, const double *array, double *rmin, double *rmax, double *rsum);
int array_minmaxmean_val(size_t len, const double *array, double *rmin, double *rmax, double *rmean);

int array_mean_val(size_t len, const double *restrict array, double *rmean);
int array_mean_val_weighted(size_t len, const double *restrict array, const double *restrict w, double missval, double *rmean);

int array_add_array(size_t len, double *restrict array1, const double *restrict array2);

#endif // _ARRAY_H


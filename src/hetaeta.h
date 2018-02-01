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
#ifndef _HETAETA_H
#define _HETAETA_H

#ifdef  HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdbool.h>

void hetaeta(bool ltq, int ngp, const int *imiss,
	     int nlev1, const double *ah1, const double *bh1,
             const double *fis1, const double *ps1, 
             const double *t1, const double *q1,
             int nlev2, const double *ah2, const double *bh2, 
             const double *fis2, double *ps2, 
             double *t2, double *q2,
	     int nvars, double **vars1, double **vars2,
	     double *tscor, double *pscor, double *secor);

#endif  /* _HETAETA_H */

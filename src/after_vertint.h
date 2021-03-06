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
#ifndef VINTERP_H
#define VINTERP_H

#include <stdbool.h>

void height2pressure(double *restrict phlev, const double *restrict hlev, long nphlev);

void presh(double *restrict fullp, double *halfp, const double *restrict vct, const double *restrict ps, long nhlev, long ngp);

void genind(int *nx, const double *restrict plev, const double *restrict fullp, long ngp, long nplev, long nhlev);
void genindmiss(int *nx, const double *restrict plev, int ngp, int nplev, const double *restrict ps_prog,
                size_t *restrict pnmiss);

void extra_P(double *restrict slp, const double *restrict halfp, const double *restrict fullp, const double *restrict geop,
             const double *restrict temp, long ngp);

void interp_T(const double *restrict geop, const double *restrict gt, double *pt, const double *restrict fullp,
              const double *restrict halfp, const int *nx, const double *restrict plev, long nplev, long ngp, long nhlev,
              double missval);
void interp_Z(const double *restrict geop, const double *restrict gz, double *pz, const double *restrict fullp,
              const double *restrict halfp, const int *nx, const double *restrict gt, const double *restrict plev, long nplev,
              long ngp, long nhlev, double missval);
void interp_X(const double *restrict gt, double *pt, const double *restrict hyb_press, const int *nx,
              const double *restrict plev, long nplev, long ngp, long nhlev, double missval);

#endif /* VINTERP_H */

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
#ifndef SPECSPACE_H
#define SPECSPACE_H

#include "afterburner.h"

struct SPTRANS
{
  long nlon;
  long nlat;
  long ntr;
  long poldim;
  long ifax[10];
  double *trig;
  double *poli;
  double *pold;
  double *pol2;    /* only for uv2dv  */
  double *pol3;    /* only for uv2dv  */
  double *coslat;  /* only for scaluv with uv2dv */
  double *rcoslat; /* only for scaluv with dv2uv */
};

struct DVTRANS
{
  int ntr;
  int fdim;
  double *f1;
  double *f2;
};

void dv2ps(const double *restrict div, double *restrict pot, long nlev, long ntr);

SPTRANS *sptrans_new(int nlon, int nlat, int ntr, int flag);
void sptrans_delete(SPTRANS *sptrans);

DVTRANS *dvtrans_new(int ntr);
void dvtrans_delete(DVTRANS *dvtrans);

void trans_uv2dv(SPTRANS *sptrans, int nlev, int gridID1, double *gu, double *gv, int gridID2, double *sd, double *svo);

void trans_dv2uv(SPTRANS *sptrans, DVTRANS *dvtrans, int nlev, int gridID1, double *sd, double *svo, int gridID2, double *gu,
                 double *gv);

void grid2spec(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut);
void spec2grid(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut);
void four2spec(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut);
void spec2four(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut);
void four2grid(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut);
void grid2four(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut);

void spec2spec(int gridIDin, double *arrayIn, int gridIDout, double *arrayOut);
void speccut(int gridIDin, double *arrayIn, double *arrayOut, int *waves);

void spcut(double *arrayIn, double *arrayOut, int ntr, int *waves);

#endif

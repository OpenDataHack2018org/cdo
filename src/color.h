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
#ifndef _COLOR_H
#define _COLOR_H

typedef struct{
  double z_low, z_high, i_dz;
  int rgb_low[3], rgb_high[3], rgb_diff[3];
  int annot;
  int skip;
}
LUT;


typedef struct {	/* For back-, fore-, and nan-colors */
  int rgb[3];
  int skip;
}
BFN_COLOR;


typedef struct {
  int       ncolors;
  LUT      *lut;
  BFN_COLOR bfn[3];
}
CPT;


int cptRead(FILE *fp, CPT *cpt);
int cptWrite(FILE *fp, CPT cpt);
int cptWriteC(FILE *fp, CPT cpt, const char *name);

#endif  /* _COLOR_H */

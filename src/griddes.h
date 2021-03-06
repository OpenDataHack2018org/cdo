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
#ifndef GRIDDES_H
#define GRIDDES_H

#include <stdbool.h>
#include <vector>

class GridDesciption
{
 public:
  std::vector<int> mask;
  std::vector<double> xvals;
  std::vector<double> yvals;
  std::vector<double> xbounds;
  std::vector<double> ybounds;
  std::vector<double> area;
  std::vector<int> rowlon;
  double xfirst, yfirst;
  double xlast, ylast;
  double xinc, yinc;
  double xpole, ypole, angle; // rotated north pole
  int scanningMode;
  /*
    scanningMode  = 128 * iScansNegatively + 64 * jScansPositively + 32 * jPointsAreConsecutive;
              64  = 128 * 0                + 64 *        1         + 32 * 0 
              00  = 128 * 0                + 64 *        0         + 32 * 0
              96  = 128 * 0                + 64 *        1         + 32 * 1
    Default  implicit scanning mode is 64: i and j scan positively, i points are consecutive (row-major)
  */
  bool uvRelativeToGrid;
  double a;
  int datatype;
  int isRotated; // TRUE for rotated grids
  int type;
  int ntr;
  bool genBounds;
  int nvertex;
  size_t size;
  size_t xsize;
  size_t ysize;
  int np;
  int lcomplex;
  bool def_xfirst;
  bool def_yfirst;
  bool def_xlast;
  bool def_ylast;
  bool def_xinc;
  bool def_yinc;
  int nd, ni, ni2, ni3;
  int number, position;
  unsigned char uuid[CDI_UUID_SIZE];
  char path[16384];
  char xname[CDI_MAX_NAME];
  char xlongname[CDI_MAX_NAME];
  char xunits[CDI_MAX_NAME];
  char xdimname[CDI_MAX_NAME];
  char yname[CDI_MAX_NAME];
  char ylongname[CDI_MAX_NAME];
  char yunits[CDI_MAX_NAME];
  char ydimname[CDI_MAX_NAME];
  char vdimname[CDI_MAX_NAME];

  void init();
  GridDesciption()
    {
      init();
    }
};

int gridDefine(GridDesciption &grid);

int gridFromNCfile(const char *gridfile);
int gridFromH5file(const char *gridfile);

#endif /* GRIDDES_H */

/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2006 Uwe Schulzweida, schulzweida@dkrz.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include <stdio.h>

#include "cdo.h"
#include "cdo_int.h"
#include "namelist.h"

#define  NPARAM  5

void *Nmltest(void *argument)
{
  int i1[5] = {-99, -99, -99, -99, -99};
  int i2    = -99;
  char lop[99] = "";
  double dm = 0;
  const int nparam = NPARAM;
  int pocc[NPARAM];
  int plen[NPARAM], ptyp[NPARAM], pdis[NPARAM];
  char *pnam[NPARAM];
  int i;
  char *var[3];

  cdoInitialize(argument);

  for ( i = 0; i < nparam; i++ ) pocc[i] = 0;

  for ( i = 0; i < 3; i++ ) var[i] = NULL;

  i = 0;
  pnam[i] = "i1";  ptyp[i] = NML_INT;    plen[i] = 5; pdis[i] = 0; i++;
  pnam[i] = "i2";  ptyp[i] = NML_INT;    plen[i] = 1; pdis[i] = 1; i++;
  pnam[i] = "lop"; ptyp[i] = NML_TEXT;   plen[i] = 99; pdis[i] = 2; i++;
  pnam[i] = "dm";  ptyp[i] = NML_DOUBLE; plen[i] = 1; pdis[i] = 1; i++;
  pnam[i] = "var"; ptyp[i] = NML_WORD;   plen[i] = 3; pdis[i] = 0; i++;

  for ( i = 0; i < nparam; i++)
    printf("%4s %4d %4d %4d %4d\n", pnam[i], ptyp[i], plen[i], pdis[i], pocc[i]);

  namelist(nparam, pnam, ptyp, plen, pocc, pdis, i1, &i2, lop, &dm, var);
  
  for ( i = 0; i < nparam; i++)
    printf("%4s %4d %4d %4d %4d\n", pnam[i], ptyp[i], plen[i], pdis[i], pocc[i]);
  
  printf("&NAMEL\n");
  printf(" I1 = %d %d %d %d %d,\n", i1[0], i1[1], i1[2], i1[3], i1[4]);
  printf(" I2 = %d,\n", i2);
  printf(" LOP = '%s',\n", lop);
  printf(" VAR = %s, %s, %s\n", var[0], var[1], var[2]);
  printf(" DM = %#g/\n", dm);

  cdoFinish();

  return (0);
}

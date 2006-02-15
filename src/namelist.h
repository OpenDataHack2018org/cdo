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

#ifndef _NAMELIST_H
#define _NAMELIST_H

#define  func_1      -1
#define  func_2      -2
#define  func_3      -3

#define  func_nix     1
#define  func_not     2
#define  func_all     8

#define  NML_INT     -1
#define  NML_DOUBLE  -2
#define  func_npr     0
#define  NML_WORD   143
#define  NML_TEXT   144
#define  func_ntu   145
#define  func_ntl   146
#define  func_nir   147
#define  func_nkw   148

struct _NAMELIST
{
  int size;
};

typedef struct _NAMELIST NAMELIST;

NAMELIST *namelistNew(void);
void      namelistDelete(NAMELIST *nml);
void      namelistDebug(int debug);
void      namelistAdd(NAMELIST *nml, char *name);
void      namelistPrint(NAMELIST *nml);

void namelist(int nparam, char *cn[], int nt[], int nl[], int nc[], int no[], ...);

#endif  /* _NAMELIST_H */

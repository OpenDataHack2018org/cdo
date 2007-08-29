/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2006 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

#include <stdio.h>


#define  NML_INT         1
#define  NML_DOUBLE      2
#define  NML_WORD        3
#define  NML_TEXT        4

#define MAX_NML_ENTRY  256

#define MAX_LINE_LEN  4096

struct _NML_LINE
{
  int nptype, namitf, namitl;
  char lineac[MAX_LINE_LEN], lineuc[MAX_LINE_LEN], linelc[MAX_LINE_LEN];
};

struct _NML_ENTRY
{
  char *name;
  void *ptr;
  int type;
  int occ;
  int dis;
  size_t size;
};

typedef struct _NML_LINE   NML_LINE;
typedef struct _NML_ENTRY  NML_ENTRY;

struct _NAMELIST
{
  int size;
  int dis;
  char *name;
  NML_LINE line;
  NML_ENTRY *entry[MAX_NML_ENTRY];
};

typedef struct _NAMELIST   NAMELIST;

NAMELIST *namelistNew(const char *name);
void namelistDelete(NAMELIST *nml);
void namelistClear(NAMELIST *nml);
void namelistDebug(int debug);
int  namelistAdd(NAMELIST *nml, const char *name, int type, int dis, void *ptr, size_t size);
void namelistPrint(NAMELIST *nml);
void namelistRead(FILE *nmlfp, NAMELIST *nml);
int  namelistNum(NAMELIST *nml, const char *name);

#endif  /* _NAMELIST_H */

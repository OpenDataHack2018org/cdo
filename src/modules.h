/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2017 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#ifndef _MODULES_H
#define _MODULES_H

#include <stdbool.h>

void *(*operatorModule(const char *operatorName))(void *);

const char **operatorHelp(const char *operatorName);

int operatorStreamInCnt(const char *operatorName);
int operatorStreamOutCnt(const char *operatorName);
int operatorStreamNumber(const char *operatorName);

void operatorPrintAll(void);
void operatorPrintList(bool print_no_output);

#endif  /* _MODULES_H */

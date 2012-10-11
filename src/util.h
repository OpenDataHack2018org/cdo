/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2012 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#ifndef _UTIL_H
#define _UTIL_H

char *getProgname(char *string);
char *getOperator(const char *argument);
char *getOperatorName(const char *xoperator);

char *makeArgument(int argc, char *argv[]);
char *getFileArg(char *argument);

enum {START_DEC, START_JAN};
int get_season_start(void);
void get_season_name(const char *seas_name[4]);

void init_is_tty(void);

void progressInit(void);
void progressStatus(double offset, double refval, double curval);

int fileExists(const char *filename);
int userFileOverwrite(const char *filename);

/* convert a CDI datatype to string */
int datatype2str(int datatype, char *datatypestr);
int str2datatype(const char *datatypestr);

#endif  /* _UTIL_H */

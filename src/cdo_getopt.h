/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2014 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#ifndef _CDO_GETOPT_H
#define _CDO_GETOPT_H

struct cdo_option {
  char *name;
  int   has_arg;
  int  *flag;
  int   val;
};

#define  no_argument        1   // no argument to the option is expect
#define  required_argument  2   // an argument to the option is required
#define  optional_argument  3   // an argument to the option may be presented.

int cdo_getopt(int argc, char * const *argv, const char *optstring);
int cdo_getopt_long(int argc, char * const *argv, const char *optstring, const struct cdo_option *longopts, int *longindex);

#endif  /* _CDO_GETOPT_H */

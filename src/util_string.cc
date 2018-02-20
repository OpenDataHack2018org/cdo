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
#include <string>
#include <stdlib.h>

std::string
string2lower(std::string str)
{
  std::string lower_case_string = str;
  for (char c : str)
    c = tolower(c);
  return lower_case_string;
}

void
strtolower(char *str)
{
  if (str)
    for (size_t i = 0; str[i]; ++i)
      str[i] = (char) tolower((int) str[i]);
}

void
strtoupper(char *str)
{
  if (str)
    for (size_t i = 0; str[i]; ++i)
      str[i] = (char) toupper((int) str[i]);
}

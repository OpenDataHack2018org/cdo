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
#ifndef STR_UTILITIES_H
#define STR_UTILITIES_H


int StringSplitWithSeperator(const char *source_string, const char *seperator, char*** ptr_split_string );

int IsNumeric (const char *s);

void StrToUpperCase ( char *sPtr );

void StrToLowerCase ( char *sPtr );

void StrReplaceChar( char *str_in, char orig_char, char rep_char );

#endif

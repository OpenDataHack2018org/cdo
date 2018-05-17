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
#ifndef NAMELIST_H_
#define NAMELIST_H_

enum namelisttype
{
  NAMELIST_UNDEFINED = 0,
  NAMELIST_OBJECT = 1,
  NAMELIST_KEY = 2,
  NAMELIST_STRING = 3,
  NAMELIST_WORD = 4
};

enum namelisterr
{
  NAMELIST_ERROR_INVAL = -1,  // Invalid character inside NAMELIST string/word
  NAMELIST_ERROR_PART = -2,   // The string is not a full NAMELIST packet, more bytes expected
  NAMELIST_ERROR_INKEY = -3,  // Invalid character inside NAMELIST key
  NAMELIST_ERROR_INTYP = -4,  // Invalid NAMELIST key type
  NAMELIST_ERROR_INOBJ = -5,  // Invalid NAMELIST object
  NAMELIST_ERROR_EMKEY = -6   // Empty key name
};

// NAMELIST token description.
struct NamelistToken
{
  int type;   // type (object, key, string word)
  int start;  // start position in NAMELIST buffer
  int end;    // end position in NAMELIST buffer
};

struct NamelistParser
{
  NamelistToken *tokens;
  unsigned int num_tokens;
  unsigned int toknext;
  unsigned int pos;
  unsigned int lineno;
};

NamelistParser *namelist_new(void);

void namelist_destroy(NamelistParser *parser);

int namelist_parse(NamelistParser *parser, const char *buf, size_t len);

void namelist_dump(NamelistParser *parser, const char *buf);

int namelist_verify(NamelistParser *parser, const char *buf);

#endif  // NAMELIST_H_

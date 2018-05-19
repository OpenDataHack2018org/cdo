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

#include <vector>

enum class NamelistType
{
  UNDEFINED = 0,
  OBJECT = 1,
  KEY = 2,
  STRING = 3,
  WORD = 4
};

enum class NamelistError
{
  UNDEFINED = 0,
  INVAL = -1,  // Invalid character inside NAMELIST string/word
  PART = -2,   // The string is not a full NAMELIST packet, more bytes expected
  INKEY = -3,  // Invalid character inside NAMELIST key
  INTYP = -4,  // Invalid NAMELIST key type
  INOBJ = -5,  // Invalid NAMELIST object
  EMKEY = -6   // Empty key name
};

// NAMELIST token description.
struct NamelistToken
{
  NamelistType type;   // type (object, key, string word)
  int start;           // start position in NAMELIST buffer
  int end;             // end position in NAMELIST buffer
};

class NamelistParser
{
 public:
  std::vector<NamelistToken> tokens;
  unsigned int num_tokens;
  unsigned int toknext;
  unsigned int pos;
  unsigned int lineno;

  void init()
    {
      num_tokens = 0;
      toknext = 0;
      pos = 0;
      lineno = 0;
    }

  NamelistParser()
    {
      init();
    }

  ~NamelistParser()
    {
    }
  
  NamelistError parse(const char *buf, size_t len);
  void dump(const char *buf);
  int verify();
};

#endif  // NAMELIST_H_

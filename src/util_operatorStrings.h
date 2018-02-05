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
#ifndef UTIL_OPERATORSTRINGS_H
#define UTIL_OPERATORSTRINGS_H

#include <array>
#include <vector>

typedef  std::vector<std::array<std::string,2>> CdoArgv;
typedef  std::array<std::string,2> cdoCommand;

cdoCommand split(std::string p_command);

CdoArgv preProcessArgv(std::vector<std::string> argv);

const char *getOperatorName(const char *operatorCommand);

char *getOperatorArg(const char *p_operatorCommand);

#endif

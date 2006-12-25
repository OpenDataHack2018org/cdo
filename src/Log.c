/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2007 Uwe Schulzweida, schulzweida@dkrz.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include "cdo.h"
#include "cdo_int.h"

void dumplogs(const char *logfilename);
void daylogs(const char *logfilename);
void monlogs(const char *logfilename);
void dumplogo(const char *logfilename);

void *Log(void *argument)
{
  int DUMPLOGS, DAYLOGS, MONLOGS, DUMPLOGO;
  int operatorID;

  cdoInitialize(argument);

  DUMPLOGS = cdoOperatorAdd("dumplogs",  0, 0, NULL);
  DAYLOGS  = cdoOperatorAdd("daylogs",   0, 0, NULL);
  MONLOGS  = cdoOperatorAdd("monlogs",   0, 0, NULL);
  DUMPLOGO = cdoOperatorAdd("dumplogo",  0, 0, NULL);

  operatorID = cdoOperatorID();

  if ( cdoStreamName(0)[0] == '-' )
    cdoAbort("This operator does not work with pipes!");

  if ( operatorID == DUMPLOGS )
    {
      dumplogs(cdoStreamName(0));
    }
  else if ( operatorID == DAYLOGS )
    {
      daylogs(cdoStreamName(0));
    }
  else if ( operatorID == MONLOGS )
    {
      monlogs(cdoStreamName(0));
    }
  else if ( operatorID == DUMPLOGO )
    {
      dumplogo(cdoStreamName(0));
    }

  cdoFinish();

  return (0);
}

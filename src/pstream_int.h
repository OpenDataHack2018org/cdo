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

#ifndef PSTREAM_INT_H
#define PSTREAM_INT_H

#include "process.h"
#include "pstream.h"
void
pstreamClose(int pstreamID);

int
pstreamInqVlist(int pstreamID);

void
pstreamDefVlist(int pstreamID, int vlistID);

int
pstreamInqTimestep(int pstreamID, int tsID);

void
pstreamDefTimestep(int pstreamID, int tsID);

int
pstreamInqFiletype(int pstreamID);
#endif 

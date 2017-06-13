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

#ifndef PSTREAM_H
#define PSTREAM_H

#include "pstream_write.h"

#include <sys/types.h> /* off_t */

int     pstreamOpenRead(const argument_t *argument);
int     pstreamOpenAppend(const argument_t *argument);
void    pstreamClose(int pstreamID);

int     pstreamInqFiletype(int pstreamID);
int     pstreamInqByteorder(int pstreamID);

int     pstreamInqVlist(int pstreamID);

int     pstreamInqTimestep(int pstreamID, int tsID);

int     pstreamInqRecord(int pstreamID, int *varID, int *levelID);

void    pstreamReadRecord(int pstreamID, double *data, int *nmiss);
void    pstreamReadRecordF(int pstreamID, float *data, int *nmiss);
void    pstreamCopyRecord(int pstreamIDdest, int pstreamIDsrc);

void    pstreamInqGRIBinfo(int pstreamID, int *intnum, float *fltnum, off_t *bignum);

int     pstreamFileID(int pstreamID);

void    cdoVlistCopyFlag(int vlistID2, int vlistID1);

#endif  /* PSTREAM_H */

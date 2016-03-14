/*
   This file is part of CDO. CDO is a collection of Operators to
   manipulate and analyse Climate model Data.

   Copyright (C) 2003-2016 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
   See COPYING file for copying and redistribution conditions.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; version 2 of the License.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   */

/*
   This module contains the following operators:

   Pack      reduce
   Pack    unreduce
   */

#if defined(_OPENMP)
#  include <omp.h>
#endif

#include <limits.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

void *MapReduce(void *argument)
{
  int gridsize;
  int nrecs;
  int gridID, varID, levelID, recID;
  int i;
  int nts;
  int nalloc = 0;
  int nmiss;
  int nlevel;
  int datatype = DATATYPE_INT16;
  dtlist_type *dtlist = dtlist_new();
  double missval1, missval2;
  field_t ***vars = NULL;

  cdoInitialize(argument);

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int nvars = vlistNvars(vlistID1);

  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
  {
    tsID++;
  }

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return 0;
}

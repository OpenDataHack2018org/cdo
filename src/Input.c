/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2006 Uwe Schulzweida, schulzweida@dkrz.de
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

*/


#include <stdio.h>
#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Input(void *argument)
{
  /*
  static char func[] = "Input";
  int i;
  int varID, recID;
  int nvars, nlevs;
  int gridsize = 0;
  int gridID, zaxisID, code, vdate, vtime;
  int nrecs;
  int levelID;
  int tsID, taxisID;
  int streamID = 0;
  int vlistID;
  int nmiss, nout;
  int nlon, nlat, level;
  int hour, minute;
  const char *format = NULL;
  double *array = NULL;
  double xdate;

  cdoInitialize(argument);

  streamID = streamOpenRead(cdoStreamName(0));
  if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", cdoStreamName(0));

  streamDefVlist(streamID, vlistID);

  gridsize = vlistGridsizeMax(vlistID);

  array = (double *) malloc(gridsize*sizeof(double));

  tsID = 0;
  taxisID = vlistInqTaxis(vlistID);
  while ( TRUE )
    {
      taxisDefVdate(taxisID, vdate);
      taxisDefVtime(taxisID, vtime);
      streamDefTimestep(streamID, tsID);

      for ( varID = 0; varID < nvars; varID++ )
	{
	  for ( levelID = 0; levelID < nlevs; levelID++ )
	    {
	      streamDefRecord(streamID, varID, levelID);

	      streamWriteRecord(streamID, array, nmiss);
	    }
	}
      tsID++;
    }

  streamClose(streamID);

  if ( array ) free(array);

  cdoFinish();
  */
  return (0);
}

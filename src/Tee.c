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

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Tee(void *argument)
{
  int nrecs;
  int varID, levelID;
  int nmiss;

  cdoInitialize(argument);

  bool lcopy = UNCHANGED_RECORD;

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int taxisID1 = vlistInqTaxis(vlistID1);

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  int streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());

  int vlistID2 = vlistDuplicate(vlistID1);
  int vlistID3 = vlistDuplicate(vlistID1);

  int taxisID2 = taxisDuplicate(taxisID1);
  int taxisID3 = taxisDuplicate(taxisID1);

  vlistDefTaxis(vlistID2, taxisID2);
  vlistDefTaxis(vlistID3, taxisID3);

  streamDefVlist(streamID2, vlistID2);
  streamDefVlist(streamID3, vlistID3);

  int gridsize = vlistGridsizeMax(vlistID1);
  double *array = (double*) Malloc(gridsize*sizeof(double));

  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      taxisCopyTimestep(taxisID3, taxisID1);

      streamDefTimestep(streamID2, tsID);
      streamDefTimestep(streamID3, tsID);
	       
      for ( int recID = 0; recID < nrecs; recID++ )
	{ 
	  if ( lcopy )
	    {
	      streamInqRecord(streamID1, &varID, &levelID);

	      streamDefRecord(streamID2,  varID,  levelID);
	      streamCopyRecord(streamID2, streamID1);

	      streamDefRecord(streamID3,  varID,  levelID);
	      streamCopyRecord(streamID3, streamID1);
	    }
	  else
	    {
	      streamInqRecord(streamID1, &varID, &levelID);
	      streamReadRecord(streamID1, array, &nmiss);

	      streamDefRecord(streamID2,  varID,  levelID);
	      streamWriteRecord(streamID2, array, nmiss);

	      streamDefRecord(streamID3,  varID,  levelID);
	      streamWriteRecord(streamID3, array, nmiss);
	    }
	}

      tsID++;
    }

  streamClose(streamID1);

  streamClose(streamID2);
  streamClose(streamID3);

  if ( array ) Free(array);

  cdoFinish();

  return 0;
}

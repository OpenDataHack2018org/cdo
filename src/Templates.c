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


void *Template1(void *argument)
{
  int nrecs;
  int varID, levelID;
  int gridsize, nmiss;

  cdoInitialize(argument);

  bool lcopy = UNCHANGED_RECORD;

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  double *array = NULL;
  if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(vlistID1);
      array = (double*) Malloc(gridsize*sizeof(double));
    }

  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamDefRecord(streamID2,  varID,  levelID);
	  
	  if ( lcopy )
	    {
	      streamCopyRecord(streamID2, streamID1);
	    }
	  else
	    {
	      streamReadRecord(streamID1, array, &nmiss);
	      streamWriteRecord(streamID2, array, nmiss);
	    }
	}

      tsID++;
    }

  streamClose(streamID1);
  streamClose(streamID2);

  vlistDestroy(vlistID2);

  if ( array ) Free(array);

  cdoFinish();

  return 0;
}


void *Template2(void *argument)
{
  int nrecs;
  int varID, levelID;
  int nmiss;

  cdoInitialize(argument);

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  int gridsize = vlistGridsizeMax(vlistID1);
  double *array = (double*) Malloc(gridsize*sizeof(double));

  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamDefRecord(streamID2,  varID,  levelID);
	  
	  streamReadRecord(streamID1, array, &nmiss);
	  streamWriteRecord(streamID2, array, nmiss);
	}

      tsID++;
    }

  streamClose(streamID1);
  streamClose(streamID2);

  vlistDestroy(vlistID2);

  if ( array ) Free(array);

  cdoFinish();

  return 0;
}

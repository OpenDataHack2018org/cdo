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

#include <string.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

/*
@BeginDoc

@BeginModule

@Name      = Copy
@Title     = Copy
@Section   = File operations
@Class     = File operation
@Arguments = ifiles ofile
@Operators = copy

@EndModule


@BeginOperator_copy

@Title     = Copy files

@BeginDesciption
Copies all input files to ofile. Each input file must have the same variables with complete timesteps.
@EndDesciption

@EndOperator

@EndDoc
*/

void *Copy(void *argument)
{
  static char func[] = "Copy";
  int streamID1, streamID2 = CDI_UNDEFID;
  int nrecs;
  int tsID1, tsID2, recID, varID, levelID;
  int gridsize;
  int vlistID1, vlistID2;
  int nmiss;
  int streamCnt, nfiles, indf;
  int taxisID1, taxisID2 = CDI_UNDEFID;
  double *array = NULL;

  cdoInitialize(argument);

  streamCnt = cdoStreamCnt();
  nfiles = streamCnt - 1;

  tsID2 = 0;
  for ( indf = 0; indf < nfiles; indf++ )
    {
      streamID1 = streamOpenRead(cdoStreamName(indf));
      if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(indf));

      vlistID1 = streamInqVlist(streamID1);
      taxisID1 = vlistInqTaxis(vlistID1);

      if ( indf == 0 )
	{
	  streamID2 = streamOpenWrite(cdoStreamName(nfiles), cdoFiletype());
	  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(nfiles));

	  vlistID2 = vlistDuplicate(vlistID1);
	  taxisID2 = taxisDuplicate(taxisID1);
	  vlistDefTaxis(vlistID2, taxisID2);

	  streamDefVlist(streamID2, vlistID2);

	  gridsize = vlistGridsizeMax(vlistID1);
	  array = (double *) malloc(gridsize*sizeof(double));
	}

      tsID1 = 0;
      while ( (nrecs = streamInqTimestep(streamID1, tsID1)) )
	{
	  taxisCopyTimestep(taxisID2, taxisID1);

	  streamDefTimestep(streamID2, tsID2);
	       
	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      streamInqRecord(streamID1, &varID, &levelID);
	      streamDefRecord(streamID2,  varID,  levelID);
	  
	      streamReadRecord(streamID1, array, &nmiss);
	      streamWriteRecord(streamID2, array, nmiss);
	    }
	  tsID1++;
	  tsID2++;
	}
      streamClose(streamID1);
    }

  streamClose(streamID2);

  if ( array ) free(array);

  cdoFinish();

  return (0);
}

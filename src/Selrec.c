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

      Selrec     selrec          Select records
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"  /* processSelf */
#include "pstream.h"
#include "error.h"
#include "util.h"
#include "list.h"


void *Selrec(void *argument)
{
  int nrecs;
  int varID, levelID;
  LIST *ilist = listNew(INT_LIST);

  cdoInitialize(argument);

  if ( processSelf() != 0 ) cdoAbort("This operator can't be combined with other operators!");

  operatorInputArg("records");

  int nsel = args2intlist(operatorArgc(), operatorArgv(), ilist);

  int *intarr = (int *) listArrayPtr(ilist);

  if ( cdoVerbose )
    {
      for ( int i = 0; i < nsel; i++ )
	cdoPrint("intarr entry: %d %d", i, intarr[i]);
    }

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int filetype = streamInqFiletype(streamID1);

  if ( filetype == CDI_FILETYPE_NC || filetype == CDI_FILETYPE_NC2 || filetype == CDI_FILETYPE_NC4 || filetype == CDI_FILETYPE_NC4C )
    cdoAbort("This operator does not work on NetCDF data!");

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  int recordID = 0;
  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
     
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  recordID++;
	  streamInqRecord(streamID1, &varID, &levelID);

	  for ( int i = 0; i < nsel; i++ )
	    {
	      if ( recordID == intarr[i] )
		{
		  streamDefRecord(streamID2, varID, levelID);
		  streamCopyRecord(streamID2, streamID1);

		  break;
		}
	    }
	}
       
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  listDelete(ilist);

  cdoFinish();

  return 0;
}

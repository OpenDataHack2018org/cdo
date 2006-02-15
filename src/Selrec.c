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
#include <math.h>   /* fabs */

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"  /* processSelf */
#include "pstream.h"
#include "error.h"
#include "util.h"
#include "functs.h"
#include "list.h"

/*
@BeginDoc

@BeginModule

@Name      = Selrec
@Title     = Selrec
@Section   = Selection
@Class     = Selection
@Arguments = ifile ofile
@Operators = selrec

@EndModule


@BeginOperator_selrec

@Title     = Select records
@Parameter = records

@BeginDesciption
Selects all fields with a record number in a user given list.
This operator does not work on netCDF data!
@EndDesciption

@BeginParameter
@Item = records
INTEGER  Comma separated list of records
@EndParameter

@EndOperator

@EndDoc
*/


void *Selrec(void *argument)
{
  int streamID1, streamID2;
  int tsID, nrecs;
  int recID, varID, levelID;
  int *intarr, nsel = 0;
  int vlistID1 = -1, vlistID2 = -1;
  int i;
  int recordID;
  int filetype;
  int taxisID1, taxisID2;
  LIST *ilist = listNew(INT_LIST);

  cdoInitialize(argument);

  if ( processSelf() != 0 && *(char *)argument == '-' )
    cdoAbort("This operator does not work with pipes!");

  operatorInputArg("records");

  nsel = args2intlist(operatorArgc(), operatorArgv(), ilist);

  intarr = (int *) listArrayPtr(ilist);

  if ( cdoVerbose )
    {
      for ( i = 0; i < nsel; i++ )
	cdoPrint("intarr entry: %d %d", i, intarr[i]);
    }

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  filetype = streamInqFiletype(streamID1);

  if ( filetype == FILETYPE_NC || filetype == FILETYPE_NC2 )
    cdoAbort("This operator does not work on netCDF data!");

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  recordID = 0;
  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
     
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  recordID++;
	  streamInqRecord(streamID1, &varID, &levelID);
	  for ( i = 0; i < nsel; i++ )
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

  return (NULL);
}

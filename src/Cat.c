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

      Copy       cat             Concatenate datasets
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"


void *Cat(void *argument)
{
  bool lconstvars = true;
  int nrecs;
  int tsID2 = 0, varID, levelID;
  int streamID2 = CDI_UNDEFID;
  int vlistID2 = CDI_UNDEFID;
  int taxisID2 = CDI_UNDEFID;
  int nmiss;
  double tw0 = 0, tw = 0;
  double *array = NULL;

  cdoInitialize(argument);

  bool lcopy = UNCHANGED_RECORD;

  int timer_cat = timer_new("cat");
  if ( cdoTimer ) timer_start(timer_cat);

  int streamCnt = cdoStreamCnt();
  int nfiles = streamCnt - 1;

  for ( int indf = 0; indf < nfiles; ++indf )
    {
      if ( cdoVerbose ) cdoPrint("Process file: %s", cdoStreamName(indf)->args);
      if ( cdoTimer ) tw0 = timer_val(timer_cat);

      int streamID1 = streamOpenRead(cdoStreamName(indf));

      int vlistID1 = streamInqVlist(streamID1);
      int taxisID1 = vlistInqTaxis(vlistID1);

      if ( indf == 0 )
	{
	  int file_exists = fileExists(cdoStreamName(nfiles)->args);
	  if ( cdoOverwriteMode ) file_exists = 0;

	  if ( file_exists )
	    {
	      streamID2 = streamOpenAppend(cdoStreamName(nfiles));

	      vlistID2 = streamInqVlist(streamID2);
	      taxisID2 = vlistInqTaxis(vlistID2);

	      vlistCompare(vlistID1, vlistID2, CMP_ALL);

	      tsID2 = vlistNtsteps(vlistID2);
	      if ( tsID2 == 0 ) tsID2 = 1; /* bug fix for time constant data only */
	    }
	  else
	    {
	      if ( cdoVerbose )
		cdoPrint("Output file doesn't exist, creating: %s", cdoStreamName(nfiles)->args);

	      streamID2 = streamOpenWrite(cdoStreamName(nfiles), cdoFiletype());

	      vlistID2 = vlistDuplicate(vlistID1);
	      taxisID2 = taxisDuplicate(taxisID1);
	      vlistDefTaxis(vlistID2, taxisID2);
	  
	      int ntsteps = vlistNtsteps(vlistID1);
	      int nvars   = vlistNvars(vlistID1);
	      
	      if ( ntsteps == 1 )
		{
		  for ( varID = 0; varID < nvars; ++varID )
		    if ( vlistInqVarTsteptype(vlistID1, varID) != TSTEP_CONSTANT ) break;
		  
		  if ( varID == nvars ) ntsteps = 0;
		}

	      if ( ntsteps == 0 && nfiles > 1 )
		{
                  lconstvars = false;
		  for ( varID = 0; varID < nvars; ++varID )
		    vlistDefVarTsteptype(vlistID2, varID, TSTEP_INSTANT);
		}

	      streamDefVlist(streamID2, vlistID2);
	    }

	  if ( ! lcopy )
	    {
	      size_t gridsize = (size_t) vlistGridsizeMax(vlistID1);
	      array = (double*) Malloc(gridsize*sizeof(double));
	    }
	}
      else
	{
	  vlistCompare(vlistID1, vlistID2, CMP_ALL);
	}

      int ntsteps = vlistNtsteps(vlistID1);

      int tsID1 = 0;
      while ( (nrecs = streamInqTimestep(streamID1, tsID1)) )
	{          
	  {
	    double fstatus = indf+1.;
	    if ( ntsteps > 1 ) fstatus = indf+(tsID1+1.)/ntsteps;
	    if ( !cdoVerbose ) progressStatus(0, 1, fstatus/nfiles);
	  }

	  taxisCopyTimestep(taxisID2, taxisID1);

	  streamDefTimestep(streamID2, tsID2);
	       
	  for ( int recID = 0; recID < nrecs; recID++ )
	    {
	      streamInqRecord(streamID1, &varID, &levelID);

              if ( lconstvars && tsID2 > 0 && tsID1 == 0 )
                if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT )
                  continue;

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

	  tsID1++;
	  tsID2++;
	}
      
      streamClose(streamID1);

      if ( cdoTimer ) tw = timer_val(timer_cat) - tw0;
      if ( cdoTimer ) cdoPrint("Processed file: %s   %.2f seconds", cdoStreamName(indf)->args, tw);
    }

  streamClose(streamID2);
 
  if ( cdoTimer ) timer_stop(timer_cat);

  if ( array ) Free(array);

  cdoFinish();

  return 0;
}

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

/*
   This module contains the following operators:

      Copy       cat             Concatenate datasets
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Cat(void *argument)
{
  bool lconstvars = true;
  int nrecs;
  int tsID2 = 0, varID, levelID;
  int streamID2 = CDI_UNDEFID;
  int vlistID2 = CDI_UNDEFID;
  int taxisID2 = CDI_UNDEFID;
  size_t nmiss;
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
      if ( cdoVerbose ) cdoPrint("Process file: %s", cdoGetStreamName(indf).c_str());
      if ( cdoTimer ) tw0 = timer_val(timer_cat);

      int streamID1 = cdoStreamOpenRead(cdoStreamName(indf));

      int vlistID1 = pstreamInqVlist(streamID1);
      int taxisID1 = vlistInqTaxis(vlistID1);

      if ( indf == 0 )
	{
          int ntsteps = vlistNtsteps(vlistID1);
          int nvars   = vlistNvars(vlistID1);
          if ( ntsteps == 1 )
            {
              for ( varID = 0; varID < nvars; ++varID )
                if ( vlistInqVarTimetype(vlistID1, varID) != TIME_CONSTANT ) break;
		  
              if ( varID == nvars ) ntsteps = 0;
            }
	      
	  bool file_exists = (!cdoOverwriteMode) ? fileExists(cdoGetStreamName(nfiles).c_str()) : false;
	  if ( file_exists )
	    {
      std::cout << "here we are 1 " << std::endl;
	      streamID2 = pstreamOpenAppend(cdoStreamName(nfiles));
      std::cout << "here we are 2" << std::endl;

	      vlistID2 = pstreamInqVlist(streamID2);
	      taxisID2 = vlistInqTaxis(vlistID2);
      std::cout << "here we are 3" << std::endl;

	      vlistCompare(vlistID1, vlistID2, CMP_ALL);

      std::cout << "here we are 4" << std::endl;
	      tsID2 = vlistNtsteps(vlistID2);
      std::cout << "here we are 5" << std::endl;
	      if ( tsID2 == 0 ) tsID2 = 1; /* bug fix for time constant data only */

              if ( ntsteps == 0 ) lconstvars = false;
            }
	  else
	    {
	      if ( cdoVerbose )
		cdoPrint("Output file doesn't exist, creating: %s", cdoGetStreamName(nfiles).c_str());

	      streamID2 = cdoStreamOpenWrite(cdoStreamName(nfiles), cdoFiletype());

	      vlistID2 = vlistDuplicate(vlistID1);
	      taxisID2 = taxisDuplicate(taxisID1);
	      vlistDefTaxis(vlistID2, taxisID2);

	      if ( ntsteps == 0 && nfiles > 1 )
		{
                  lconstvars = false;
		  for ( varID = 0; varID < nvars; ++varID )
		    vlistDefVarTimetype(vlistID2, varID, TIME_VARYING);
		}

	      pstreamDefVlist(streamID2, vlistID2);
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
      while ( (nrecs = pstreamInqTimestep(streamID1, tsID1)) )
	{          
          double fstatus = (ntsteps > 1) ? indf+(tsID1+1.)/ntsteps : indf+1.;
          if ( !cdoVerbose ) progressStatus(0, 1, fstatus/nfiles);

	  taxisCopyTimestep(taxisID2, taxisID1);

	  pstreamDefTimestep(streamID2, tsID2);
	       
	  for ( int recID = 0; recID < nrecs; recID++ )
	    {
	      pstreamInqRecord(streamID1, &varID, &levelID);

              if ( lconstvars && tsID2 > 0 && tsID1 == 0 )
                if ( vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT )
                  continue;

	      pstreamDefRecord(streamID2, varID, levelID);

	      if ( lcopy )
		{
		  pstreamCopyRecord(streamID2, streamID1); 
		}
	      else
		{
		  pstreamReadRecord(streamID1, array, &nmiss);
		  pstreamWriteRecord(streamID2, array, nmiss);
		}
	    }

	  tsID1++;
	  tsID2++;
	}
      
      pstreamClose(streamID1);

      if ( cdoTimer ) tw = timer_val(timer_cat) - tw0;
      if ( cdoTimer ) cdoPrint("Processed file: %s   %.2f seconds", cdoGetStreamName(indf).c_str(), tw);
    }

  pstreamClose(streamID2);
 
  if ( cdoTimer ) timer_stop(timer_cat);

  if ( array ) Free(array);

  cdoFinish();

  return 0;
}

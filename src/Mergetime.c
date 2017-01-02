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

/*
   This module contains the following operators:

      Merge      mergetime       Merge datasets sorted by date and time
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"


void *Mergetime(void *argument)
{
  int streamID1;
  int tsID2 = 0, varID, levelID;
  int vlistID1, vlistID2;
  int fileID;
  int taxisID1, taxisID2 = CDI_UNDEFID;
  int nmiss;
  int vdate, vtime;
  int last_vdate = -1, last_vtime = -1;
  int next_fileID;
  bool skip_same_time = false;
  double *array = NULL;
  typedef struct
  {
    int streamID;
    int vlistID;
    int taxisID;
    int tsID;
    int vdate;
    int vtime;
    int nrecs;
  } sfile_t;

  cdoInitialize(argument);

  {
    char *envstr = getenv("SKIP_SAME_TIME");
    if ( envstr )
      {
	int ival;
	ival = atoi(envstr);
	if ( ival == 1 )
	  {
	    skip_same_time = true;
	    if ( cdoVerbose )
	      cdoPrint("Set SKIP_SAME_TIME to %d", ival);
	  }
      }
  }

  bool lcopy = UNCHANGED_RECORD;

  int nfiles = cdoStreamCnt() - 1;

  sfile_t *sf = (sfile_t*) Malloc(nfiles*sizeof(sfile_t));

  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      if ( cdoVerbose ) cdoPrint("process: %s", cdoStreamName(fileID)->args);

      streamID1 = streamOpenRead(cdoStreamName(fileID));

      vlistID1 = streamInqVlist(streamID1);
      taxisID1 = vlistInqTaxis(vlistID1);

      sf[fileID].streamID = streamID1;
      sf[fileID].vlistID  = vlistID1;
      sf[fileID].taxisID  = taxisID1;
    }

  
  /* check that the contents is always the same */
  for ( fileID = 1; fileID < nfiles; fileID++ )
    vlistCompare(sf[0].vlistID, sf[fileID].vlistID, CMP_ALL);

  /* read the first time step */
  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      sf[fileID].tsID = 0;
      sf[fileID].nrecs = streamInqTimestep(sf[fileID].streamID, sf[fileID].tsID);
      if ( sf[fileID].nrecs == 0 )
	{
	  streamClose(sf[fileID].streamID);
	  sf[fileID].streamID = -1;
	}
      else
	{
	  sf[fileID].vdate = taxisInqVdate(sf[fileID].taxisID);
	  sf[fileID].vtime = taxisInqVtime(sf[fileID].taxisID);
	}
    }

  const char *ofilename = cdoStreamName(nfiles)->args;

  if ( !cdoOverwriteMode && fileExists(ofilename) && !userFileOverwrite(ofilename) )
    cdoAbort("Outputfile %s already exists!", ofilename);

  int streamID2 = streamOpenWrite(cdoStreamName(nfiles), cdoFiletype());

  if ( ! lcopy )
    {
      size_t gridsize = vlistGridsizeMax(sf[0].vlistID);
      array = (double*) Malloc(gridsize*sizeof(double));
    }

  while ( TRUE )
    {
      bool process_timestep = true;

      next_fileID = -1;
      vdate = 0;
      vtime = 0;
      for ( fileID = 0; fileID < nfiles; fileID++ )
	{
	  if ( sf[fileID].streamID != -1 )
	    if ( next_fileID == -1 || sf[fileID].vdate < vdate ||
		 (sf[fileID].vdate == vdate && sf[fileID].vtime < vtime) )
	      {
		next_fileID = fileID;
		vdate = sf[fileID].vdate;
		vtime = sf[fileID].vtime;
	      }
	}

      fileID = next_fileID;

      if ( cdoVerbose )
	cdoPrint("nextstep = %d  vdate = %d  vtime = %d", next_fileID, vdate, vtime);

      if ( next_fileID == -1 ) break;

      if ( skip_same_time )
	if ( vdate == last_vdate && vtime == last_vtime )
	  {
	    char vdatestr[32], vtimestr[32];
	    date2str(vdate, vdatestr, sizeof(vdatestr));
	    time2str(vtime, vtimestr, sizeof(vtimestr));
	    cdoPrint("Timestep %4d in stream %d (%s %s) already exists, skipped!",
		     sf[fileID].tsID+1, sf[fileID].streamID, vdatestr, vtimestr);
	    process_timestep = false;
	  }

      if ( process_timestep )
	{
	  if ( tsID2 == 0 )
	    {
	      vlistID1 = sf[0].vlistID;
	      vlistID2 = vlistDuplicate(vlistID1);
	      taxisID1 = vlistInqTaxis(vlistID1);
	      taxisID2 = taxisDuplicate(taxisID1);
	      vlistDefTaxis(vlistID2, taxisID2);
	      
	      streamDefVlist(streamID2, vlistID2);
	    }

	  last_vdate = vdate;
	  last_vtime = vtime;

	  taxisCopyTimestep(taxisID2, sf[fileID].taxisID);

	  streamDefTimestep(streamID2, tsID2);
	       
	  for ( int recID = 0; recID < sf[fileID].nrecs; recID++ )
	    {
	      streamInqRecord(sf[fileID].streamID, &varID, &levelID);

              if ( tsID2 > 0 && sf[fileID].tsID == 0 )
                if ( vlistInqVarTsteptype(sf[fileID].vlistID, varID) == TSTEP_CONSTANT )
                  continue;

              streamDefRecord(streamID2, varID, levelID);
	  
	      if ( lcopy )
		{
		  streamCopyRecord(streamID2, sf[fileID].streamID); 
		}
	      else
		{
		  streamReadRecord(sf[fileID].streamID, array, &nmiss);
		  streamWriteRecord(streamID2, array, nmiss);
		}
	    }

	  tsID2++;
	}

      sf[fileID].nrecs = streamInqTimestep(sf[fileID].streamID, ++sf[fileID].tsID);
      if ( sf[fileID].nrecs == 0 )
	{
	  streamClose(sf[fileID].streamID);
	  sf[fileID].streamID = -1;
	}
      else
	{
	  sf[fileID].vdate = taxisInqVdate(sf[fileID].taxisID);
	  sf[fileID].vtime = taxisInqVtime(sf[fileID].taxisID);
	}
    }

  streamClose(streamID2);

  if ( ! lcopy )
    if ( array ) Free(array);

  if ( sf ) Free(sf);

  cdoFinish();

  return 0;
}

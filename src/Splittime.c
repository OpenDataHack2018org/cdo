/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2010 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Splittime  splithour       Split hours
      Splittime  splitday        Split days
      Splittime  splitmon        Split months
      Splittime  splitseas       Split seasons
*/


#include <string.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"


#define  MAX_STREAMS 32

void *Splittime(void *argument)
{
  static char func[] = "Splittime";
  int SPLITHOUR, SPLITDAY, SPLITMON, SPLITSEAS;
  int operatorID;
  int operfunc, operintval;
  int nchars;
  int streamID1, streamID2;
  int varID;
  int nrecs;
  int tsID, recID, levelID;
  int vlistID1, vlistID2;
  int  streamIDs[MAX_STREAMS], tsIDs[MAX_STREAMS];
  char filesuffix[32];
  char filename[1024];
  int index = 0;
  int i;
  int taxisID1, taxisID2;
  int vdate, vtime;
  int lcopy = FALSE;
  int gridsize;
  int nmiss;
  double *array = NULL;
  int season_start;
  const char *seas_name[4];

  cdoInitialize(argument);

  SPLITHOUR = cdoOperatorAdd("splithour", func_time, 10000, NULL);
  SPLITDAY  = cdoOperatorAdd("splitday",  func_date,     1, NULL);
  SPLITMON  = cdoOperatorAdd("splitmon",  func_date,   100, NULL);
  SPLITSEAS = cdoOperatorAdd("splitseas", func_date,   100, NULL);

  operatorID = cdoOperatorID();
  operfunc   = cdoOperatorFunc(operatorID);
  operintval = cdoOperatorIntval(operatorID);

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  season_start = get_season_start();
  get_season_name(seas_name);

  for ( i = 0; i < MAX_STREAMS; i++ ) streamIDs[i] = -1;
  for ( i = 0; i < MAX_STREAMS; i++ ) tsIDs[i] = 0;

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  strcpy(filename, cdoStreamName(1));
  nchars = strlen(filename);

  filesuffix[0] = 0;
  if ( cdoDisableFilesuffix == FALSE )
    {
      strcat(filesuffix, streamFilesuffix(cdoDefaultFileType));
      if ( cdoDefaultFileType == FILETYPE_GRB )
	if ( vlistIsSzipped(vlistID1) || cdoZtype == COMPRESS_SZIP )
	  strcat(filesuffix, ".sz");
    }

  if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(vlistID1);
      array = (double *) malloc(gridsize*sizeof(double));
    }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      if ( operfunc == func_date )
	{
	  index = (vdate/operintval)%100;

	  if ( operatorID == SPLITSEAS )
	    {
	      if ( index < 0 || index > 16 )
		cdoAbort("Month %d out of range!", index);

	      if ( season_start == START_DEC )
		{
		  if ( index <= 12 )
		    index = (index % 12) / 3;
		  else
		    index = index - 13;
		}
	      else
		{
		  if ( index <= 12 )
		    index = (index - 1) / 3;
		  else
		    index = index - 13;
		}
	      
	      if ( index < 0 || index > 3 )
		cdoAbort("Season %d out of range!", index+1);
	    }
	}
      else if ( operfunc == func_time )
	{
	  index = (vtime/operintval)%100;
	}

      if ( index < 0 || index >= MAX_STREAMS )
	cdoAbort("Index out of range!");

      streamID2 = streamIDs[index];
      if ( streamID2 < 0 )
	{
	  if ( operatorID == SPLITSEAS )
	    {
	      sprintf(filename+nchars, "%3s", seas_name[index]);
	      if ( filesuffix[0] )
		sprintf(filename+nchars+3, "%s", filesuffix);
	    }
	  else
	    {
	      sprintf(filename+nchars, "%02d", index);
	      if ( filesuffix[0] )
		sprintf(filename+nchars+2, "%s", filesuffix);
	    }

	  if ( cdoVerbose ) cdoPrint("create file %s", filename);

	  streamID2 = streamOpenWrite(filename, cdoFiletype());
	  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", filename);

	  streamDefVlist(streamID2, vlistID2);

	  streamIDs[index] = streamID2;
	}

      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsIDs[index]++);

      for ( recID = 0; recID < nrecs; recID++ )
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

  for ( index = 0; index < MAX_STREAMS; index++ )
    {
      streamID2 = streamIDs[index];
      if ( streamID2 >= 0 )  streamClose(streamID2);
    }
 
  if ( ! lcopy )
    if ( array ) free(array);

  cdoFinish();

  return (0);
}

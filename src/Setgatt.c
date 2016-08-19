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

      Setgatt    setgatt         Set global attribute
      Setgatt    setgatts        Set global attributes
*/

#include <ctype.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Setgatt(void *argument)
{
  int nrecs;
  int varID, levelID;
  int gridsize;
  int nmiss;
  char *attname = NULL, *attstring = NULL, *attfile = NULL;

  cdoInitialize(argument);

  int SETGATT  = cdoOperatorAdd("setgatt",  0, 0, "attribute name and string");
                 cdoOperatorAdd("setgatts", 0, 0, NULL);

  int operatorID = cdoOperatorID();

  if ( operatorID == SETGATT )
    {
      operatorInputArg(cdoOperatorEnter(operatorID));
      attname   = operatorArgv()[0];
      attstring = operatorArgv()[1];
    }
  else
    {
      attfile   = operatorArgv()[0];
    }

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  if ( operatorID == SETGATT )
    {
      cdiDefAttTxt(vlistID2, CDI_GLOBAL, attname, (int)strlen(attstring), attstring);
    }
  else
    {
      char line[1024];
      int attlen = 0;

      FILE *fp = fopen(attfile, "r");
      if ( fp == 0 ) cdoAbort("Open failed on %s", attfile);

      while ( readline(fp, line, 1024) )
	{
	  attlen = 0;
	  if ( line[0] == '#' ) continue;
	  if ( line[0] == '\0' ) continue;
	  attname = line;
	  while ( isspace((int) *attname) ) attname++;
	  if ( attname[0] == '\0' ) continue;
	  attstring = attname;
	  while ( *attstring != ' ' && *attstring != '\0' &&
		  *attstring != '=' && *attstring != '"' ) attstring++;
	  if ( *attstring == '\0' )
	    attstring = NULL;
	  else
	    {
	      *attstring = '\0';
	      attstring++;
	      while ( isspace((int) *attstring) || (int) *attstring == '=' ||
		      (int) *attstring == '"' || (int) *attstring == '\'' ) attstring++;
	      attlen = strlen(attstring);
	      if ( attstring[attlen-1] == '"' || attstring[attlen-1] == '\'' )
		attstring[--attlen] = 0;
	    }

	  if ( attstring && attlen)
	    cdiDefAttTxt(vlistID2, CDI_GLOBAL, attname, attlen, attstring);
	}

      fclose(fp);
    }

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);
  if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
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

  if ( array ) Free(array);

  cdoFinish();

  return 0;
}

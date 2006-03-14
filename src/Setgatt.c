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
#include <ctype.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

/*
@BeginDoc

@BeginModule

@Name      = Setgatt
@Title     = Set global attribute
@Section   = Manipulating the header/field
@Class     = Manipulation
@Arguments = ifile ofile
@Operators = setgatt setgatts

@EndModule


@BeginOperator_setgatt

@Title     = Set global attribute
@Parameter = attname attstring

@BeginDescription
Sets one user defined global text attribute.
@EndDescription

@BeginParameter
@Item = attname,attstring
STRING  Name and text of the global attribute
@EndParameter

@EndOperator


@BeginOperator_setgatts

@Title     = Set global attributes
@Parameter = attfile

@BeginDescription
Sets user defined global text attributes. The name and text
of the global attrubutes are read from a file.
@EndDescription

@BeginParameter
@Item = attfile
STRING  File name which contains global attributes
@EndParameter

@EndOperator

@EndDoc
*/

void *Setgatt(void *argument)
{
  static char func[] = "Setgatt";
  int SETGATT, SETGATTS;
  int operatorID;
  int streamID1, streamID2 = CDI_UNDEFID;
  int nrecs;
  int tsID, recID, varID, levelID;
  int vlistID1, vlistID2;
  int gridsize;
  int nmiss;
  char *attname = NULL, *attstring = NULL, *attfile = NULL;
  double *array = NULL;
  int taxisID1, taxisID2;

  cdoInitialize(argument);

  SETGATT  = cdoOperatorAdd("setgatt",  0, 0, "attribute name and string");
  SETGATTS = cdoOperatorAdd("setgatts", 0, 0, NULL);

  operatorID = cdoOperatorID();

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

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  if ( operatorID == SETGATT )
    {
      vlistDefAttribute(vlistID2, attname, attstring);
    }
  else
    {
      char line[1024];
      FILE *fp;

      fp = fopen(attfile, "r");
      if ( fp == 0 ) cdoAbort("Open failed on %s", attfile);

      while ( readline(fp, line, 1024) )
	{
	  if ( line[0] == '#' ) continue;
	  if ( line[0] == '\0' ) continue;
	  attname = line;
	  while ( isspace((int) *attname) ) attname++;
	  if ( attname[0] == '\0' ) continue;
	  attstring = attname;
	  while ( *attstring != ' ' && *attstring != '\0' ) attstring++;
	  if ( *attstring == '\0' )
	    attstring = NULL;
	  else
	    {
	      *attstring = '\0';
	      attstring++;
	      while ( isspace((int) *attstring) ) attstring++;
	    }

	  vlistDefAttribute(vlistID2, attname, attstring);
	}

      fclose(fp);
    }

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);
  array = (double *) malloc(gridsize*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
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

  if ( array ) free(array);

  cdoFinish();

  return (0);
}

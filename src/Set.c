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

@Name      = Set
@Title     = Set
@Section   = Manipulating the header/field
@Class     = Manipulation
@Arguments = ifile ofile
@Operators = setpartab setcode setvar setlevel

@EndModule


@BeginOperator_setpartab

@Title     = Set parameter table
@Parameter = table

@BeginDesciption
Sets the parameter table for all variables.
@EndDesciption

@BeginParameter
@Item = table
STRING  Parameter table file or name
@EndParameter

@EndOperator


@BeginOperator_setcode

@Title     = Set code
@Parameter = code

@BeginDesciption
Sets the code for all variables to the same given value.
@EndDesciption

@BeginParameter
@Item = code
INTEGER Code number
@EndParameter

@EndOperator


@BeginOperator_setvar

@Title     = Set variable name
@Parameter = name

@BeginDesciption
Sets the name of the first variable.
@EndDesciption

@BeginParameter
@Item = name
STRING  Variable name
@EndParameter

@EndOperator


@BeginOperator_setlevel

@Title     = Set level
@Parameter = level

@BeginDesciption
Sets the first level of all variables.
@EndDesciption

@BeginParameter
@Item = level
FLOAT New level
@EndParameter

@EndOperator

@EndDoc
*/


void *Set(void *argument)
{
  static char func[] = "Set";
  int SETPARTAB, SETCODE, SETVAR, SETLEVEL;
  int operatorID;
  int streamID1, streamID2 = CDI_UNDEFID;
  int nrecs, nvars, newval = -1;
  int tsID1, recID, varID, levelID;
  int vlistID1, vlistID2;
  int taxisID1, taxisID2;
  int nmiss;
  int gridsize;
  int index, zaxisID1, zaxisID2, nzaxis, nlevs;
  int tableID = -1;
  char *newname = NULL, *partab = NULL;
  double newlevel = 0;
  double *levels = NULL;
  double *array = NULL;

  cdoInitialize(argument);

  SETPARTAB = cdoOperatorAdd("setpartab", 0, 0, "parameter table");
  SETCODE   = cdoOperatorAdd("setcode",   0, 0, "code");
  SETVAR    = cdoOperatorAdd("setvar",    0, 0, "variable name");
  SETLEVEL  = cdoOperatorAdd("setlevel",  0, 0, "level");

  operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));
  if ( operatorID == SETCODE )
    {
      newval = atoi(operatorArgv()[0]);
    }
  else if ( operatorID == SETVAR )
    {
      newname = operatorArgv()[0];
    }
  else if ( operatorID == SETPARTAB )
    {
      partab = operatorArgv()[0];
      tableID = defineTable(partab);
    }
  else if ( operatorID == SETLEVEL )
    {
      newlevel = atof(operatorArgv()[0]);
    }

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if ( operatorID == SETCODE )
    {
      nvars = vlistNvars(vlistID2);
      for ( varID = 0; varID < nvars; varID++ )
	vlistDefVarCode(vlistID2, varID, newval);
    }
  else if ( operatorID == SETVAR )
    {
      vlistDefVarName(vlistID2, 0, newname);
    }
  else if ( operatorID == SETPARTAB )
    {
      nvars = vlistNvars(vlistID2);
      for ( varID = 0; varID < nvars; varID++ )
	vlistDefVarTable(vlistID2, varID, tableID);
    }
  else if ( operatorID == SETLEVEL )
    {
      nzaxis = vlistNzaxis(vlistID2);
      for ( index = 0; index < nzaxis; index++ )
	{
	  zaxisID1 = vlistZaxis(vlistID2, index);
	  zaxisID2 = zaxisDuplicate(zaxisID1);
	  nlevs = zaxisInqSize(zaxisID2);
	  levels = (double *) malloc(nlevs*sizeof(double));
	  zaxisInqLevels(zaxisID2, levels);
	  levels[0] = newlevel;
	  zaxisDefLevels(zaxisID2, levels);
	  vlistChangeZaxis(vlistID2, index, zaxisID2);
	  free(levels);
	}
    }

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID2);
  array = (double *) malloc(gridsize*sizeof(double));

  tsID1 = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID1)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID1);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamDefRecord(streamID2,  varID,  levelID);
	  
	  streamReadRecord(streamID1, array, &nmiss);
	  streamWriteRecord(streamID2, array, nmiss);
	}
      tsID1++;
    }

  streamClose(streamID1);
  streamClose(streamID2);

  if ( array ) free(array);

  cdoFinish();

  return (0);
}

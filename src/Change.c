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
#include "cdo_int.h"
#include "pstream.h"

/*
@BeginDoc

@BeginModule

@Name      = Change
@Title     = Change
@Section   = Manipulating the header/field
@Class     = Manipulation
@Arguments = ifile ofile
@Operators = chcode chvar chlevel chlevelc chlevelv

@EndModule


@BeginOperator_chcode

@Title     = Change code
@Parameter = ocode ncode ...

@BeginDesciption
Changes some user given codes to new user given values.
@EndDesciption

@BeginParameter
@Item = ocode,ncode,...
INTEGER  Pairs of old and new code
@EndParameter

@EndOperator


@BeginOperator_chvar

@Title     = Change variable name
@Parameter = ovar nvar ...

@BeginDesciption
Changes some user given variable names to new user given names.
@EndDesciption

@BeginParameter
@Item = ovar,nvar,...
STRING  Pairs of old and new variable name
@EndParameter

@EndOperator


@BeginOperator_chlevel

@Title     = Change level
@Parameter = olevel nlevel ...

@BeginDesciption
Changes some user given levels to new user given values.
@EndDesciption

@BeginParameter
@Item = olevel,nlevel,...
FLOAT  Pairs of old and new level
@EndParameter

@EndOperator


@BeginOperator_chlevelc

@Title     = Change level of one code
@Parameter = code olevel nlevel

@BeginDesciption
Changes one level of a user given code number.
@EndDesciption

@BeginParameter olevel
@Item = code
INTEGER Code number
@Item = olevel
FLOAT   Old level
@Item = nlevel
FLOAT   New level
@EndParameter

@EndOperator


@BeginOperator_chlevelv

@Title     = Change level of one variable
@Parameter = var olevel nlevel

@BeginDesciption
Changes one level of a user given variable.
@EndDesciption

@BeginParameter olevel
@Item = var
STRING Variable name
@Item = olevel
FLOAT  Old level
@Item = nlevel
FLOAT  New level
@EndParameter

@EndOperator

@EndDoc
*/

#define  MAXARG     16384

void *Change(void *argument)
{
  static char func[] = "Change";
  int CHCODE, CHVAR, CHLEVEL, CHLEVELC, CHLEVELV;
  int operatorID;
  int streamID1, streamID2 = CDI_UNDEFID;
  int nrecs, nvars;
  int tsID1, recID, varID = 0, levelID;
  int vlistID1, vlistID2;
  int taxisID1, taxisID2;
  int chcodes[MAXARG], nch = 0;
  char *chvars[MAXARG];
  char varname[128];
  char *chvar = NULL;
  int chcode = 0;
  int code, i;
  int nmiss;
  int gridsize;
  int nfound;
  int nzaxis, zaxisID1, zaxisID2, k, nlevs, index;
  double chlevels[MAXARG];
  double *levels = NULL;
  double *array = NULL;

  cdoInitialize(argument);

  CHCODE   = cdoOperatorAdd("chcode",   0, 0, "pairs of old and new code");
  CHVAR    = cdoOperatorAdd("chvar",    0, 0, "pairs of old and new variable name");
  CHLEVEL  = cdoOperatorAdd("chlevel",  0, 0, "pairs of old and new level");
  CHLEVELC = cdoOperatorAdd("chlevelc", 0, 0, "code number, old and new level");
  CHLEVELV = cdoOperatorAdd("chlevelv", 0, 0, "variable name, old and new level");

  operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  nch = operatorArgc();

  if ( operatorID == CHCODE )
    {
      if ( nch%2 ) cdoAbort("Odd number of input arguments!");
      for ( i = 0; i < nch; i++ )
	chcodes[i] = atoi(operatorArgv()[i]);
    }
  else if ( operatorID == CHVAR )
    {
      if ( nch%2 ) cdoAbort("Odd number of input arguments!");
      for ( i = 0; i < nch; i++ )
	chvars[i] = operatorArgv()[i];
    }
  else if ( operatorID == CHLEVEL )
    {
      if ( nch%2 ) cdoAbort("Odd number of input arguments!");
      for ( i = 0; i < nch; i++ )
	chlevels[i] = atof(operatorArgv()[i]);
    }
  else if ( operatorID == CHLEVELC )
    {
      operatorCheckArgc(3);
      
      chcode = atoi(operatorArgv()[0]);
      chlevels[0] = atof(operatorArgv()[1]);
      chlevels[1] = atof(operatorArgv()[2]);
    }
  else if ( operatorID == CHLEVELV )
    {
      operatorCheckArgc(3);
      
      chvar = operatorArgv()[0];
      chlevels[0] = atof(operatorArgv()[1]);
      chlevels[1] = atof(operatorArgv()[2]);
    }

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if ( operatorID == CHCODE )
    {
      nvars = vlistNvars(vlistID2);
      for ( varID = 0; varID < nvars; varID++ )
	{
	  code = vlistInqVarCode(vlistID2, varID);
	  for ( i = 0; i < nch; i += 2 )
	    if ( code == chcodes[i] )
	      vlistDefVarCode(vlistID2, varID, chcodes[i+1]);
	}
    }
  else if ( operatorID == CHVAR )
    {
      nvars = vlistNvars(vlistID2);
      for ( varID = 0; varID < nvars; varID++ )
	{
	  vlistInqVarName(vlistID2, varID, varname);
	  for ( i = 0; i < nch; i += 2 )
	    if ( strcmp(varname, chvars[i]) == 0 )
	      vlistDefVarName(vlistID2, varID, chvars[i+1]);
	}
    }
  else if ( operatorID == CHLEVEL )
    {
      nzaxis = vlistNzaxis(vlistID2);
      for ( index = 0; index < nzaxis; index++ )
	{
	  zaxisID1 = vlistZaxis(vlistID2, index);
	  nlevs = zaxisInqSize(zaxisID1);
	  levels = (double *) malloc(nlevs*sizeof(double));
	  zaxisInqLevels(zaxisID1, levels);
	  nfound = 0;
	  for ( i = 0; i < nch; i += 2 )
	    for ( k = 0; k < nlevs; k++ )
	      if ( fabs(levels[k] - chlevels[i]) < 0.0001 ) nfound++;

	  if ( nfound )
	    {
	      zaxisID2 = zaxisDuplicate(zaxisID1);
	      for ( i = 0; i < nch; i += 2 )
		for ( k = 0; k < nlevs; k++ )
		  if ( fabs(levels[k] - chlevels[i]) < 0.001 )
		    levels[k] = chlevels[i+1];

	      zaxisDefLevels(zaxisID2, levels);
	      vlistChangeZaxis(vlistID2, index, zaxisID2);
	    }

	  free(levels);
	}
    }
  else if ( operatorID == CHLEVELC || operatorID == CHLEVELV )
    {
      nvars = vlistNvars(vlistID2);
      if ( operatorID == CHLEVELC )
	{
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      code = vlistInqVarCode(vlistID2, varID);
	      if ( code == chcode ) break;
	    }
	  if ( varID == nvars ) cdoAbort("Code %d not found!", chcode);
	}
      else
	{
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      vlistInqVarName(vlistID2, varID, varname);
	      if ( strcmp(varname, chvar) == 0 ) break;
	    }
	  if ( varID == nvars ) cdoAbort("Variable name %s not found!", chvar);
	}

      zaxisID1 = vlistInqVarZaxis(vlistID2, varID);
      nlevs = zaxisInqSize(zaxisID1);
      levels = (double *) malloc(nlevs*sizeof(double));
      zaxisInqLevels(zaxisID1, levels);
      nfound = 0;
      for ( k = 0; k < nlevs; k++ )
	if ( fabs(levels[k] - chlevels[0]) < 0.0001 ) nfound++;

      if ( nfound )
	{
	  zaxisID2 = zaxisDuplicate(zaxisID1);
	  for ( k = 0; k < nlevs; k++ )
	    if ( fabs(levels[k] - chlevels[0]) < 0.001 )
	      levels[k] = chlevels[1];

	  zaxisDefLevels(zaxisID2, levels);
	  vlistChangeVarZaxis(vlistID2, varID, zaxisID2);
	}
      else
	cdoAbort("Level %g not found!", chlevels[0]);

      free(levels);
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

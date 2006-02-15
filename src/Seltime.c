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
#include "cdo_int.h"  /* processSelf */
#include "pstream.h"
#include "error.h"
#include "util.h"
#include "functs.h"
#include "list.h"

/*
@BeginDoc

@BeginModule

@Name      = Seltime
@Title     = Select time
@Section   = Selection
@Class     = Selection
@Arguments = ifile ofile
@Operators = seltimestep seltime selhour selday selmon selseas selyear seldate

@EndModule


@BeginOperator_seltimestep

@Title     = Select timesteps
@Parameter = timesteps

@BeginDesciption
Selects all timesteps with a timestep in a user given list.
@EndDesciption

@BeginParameter
@Item = timesteps
INTEGER  Comma separated list of timesteps
@EndParameter

@EndOperator


@BeginOperator_seltime

@Title     = Select times
@Parameter = times

@BeginDesciption
Selects all timesteps with a time in a user given list.
@EndDesciption

@BeginParameter
@Item = times
STRING  Comma separated list of times (format HH:MM)
@EndParameter

@EndOperator


@BeginOperator_selhour

@Title     = Select hours
@Parameter = hours

@BeginDesciption
Selects all timesteps with a hour in a user given list.
@EndDesciption

@BeginParameter
@Item = hours
INTEGER  Comma separated list of hours
@EndParameter

@EndOperator


@BeginOperator_selday

@Title     = Select days
@Parameter = days

@BeginDesciption
Selects all timesteps with a day in a user given list.
@EndDesciption

@BeginParameter
@Item = days
INTEGER  Comma separated list of days
@EndParameter

@EndOperator


@BeginOperator_selmon

@Title     = Select months
@Parameter = months

@BeginDesciption
Selects all timesteps with a month in a user given list.
@EndDesciption

@BeginParameter
@Item = months
INTEGER  Comma separated list of months
@EndParameter

@EndOperator


@BeginOperator_selseas

@Title     = Select seasons
@Parameter = seasons

@BeginDesciption
Selects all timesteps with a month of a season in a user given list.
@EndDesciption

@BeginParameter
@Item = seasons
STRING   Comma separated list of seasons (DJF, MAM, JJA, SON)
@EndParameter

@EndOperator


@BeginOperator_selyear

@Title     = Select years
@Parameter = years

@BeginDesciption
Selects all timesteps with a year in a user given list.
@EndDesciption

@BeginParameter
@Item = years
INTEGER  Comma separated list of years
@EndParameter

@EndOperator


@BeginOperator_seldate

@Title     = Select dates
@Parameter = date1 [date2]

@BeginDesciption
Selects all timesteps with a date in a given range.
@EndDesciption

@BeginParameter
@Item = date1
STRING  Start date (format YYYY-MM-DD)
@Item = date2
STRING  End date (format YYYY-MM-DD)
@EndParameter

@EndOperator

@EndDoc
*/

#define  NOPERATORS   8

void *Seltime(void *argument)
{
  const char func[] = "Seltime";
  int SELTIMESTEP, SELDATE, SELTIME, SELHOUR, SELDAY, SELMON, SELYEAR, SELSEAS;
  int operatorID;
  int operfunc, intval;
  int moddat[NOPERATORS];
  int streamID1, streamID2;
  int tsID, tsID2, nrecs;
  int recID, varID, levelID;
  int *intarr, nsel = 0, selval;
  int vlistID1 = -1, vlistID2 = -1;
  int taxisID1, taxisID2;
  int vdate, vtime;
  int copytimestep;
  int i;
  int lcopy = FALSE;
  int gridsize;
  int nmiss;
  int year, month, day, hour, minute;
  double *array = NULL;
  LIST *ilist = listNew(INT_LIST);

  cdoInitialize(argument);

  SELTIMESTEP = cdoOperatorAdd("seltimestep", func_step,     1, "timesteps");
  SELDATE     = cdoOperatorAdd("seldate",     func_date,     1, "start date and end date (format YYYY-MM-DD)");
  SELTIME     = cdoOperatorAdd("seltime",     func_time,     1, "times (format HH:MM)");
  SELHOUR     = cdoOperatorAdd("selhour",     func_time,   100, "hours");
  SELDAY      = cdoOperatorAdd("selday",      func_date,     1, "days");
  SELMON      = cdoOperatorAdd("selmon",      func_date,   100, "months");
  SELYEAR     = cdoOperatorAdd("selyear",     func_date, 10000, "years");
  SELSEAS     = cdoOperatorAdd("selseas",     func_date,   100, "seasons");

  moddat[SELTIMESTEP] =          1;
  moddat[SELDATE]     = 1000000000;
  moddat[SELTIME]     =      10000;
  moddat[SELHOUR]     =        100;
  moddat[SELDAY]      =        100;
  moddat[SELMON]      =        100;
  moddat[SELYEAR]     =      10000;
  moddat[SELSEAS]     =        100;

  operatorID = cdoOperatorID();

  operfunc = cdoOperatorFunc(operatorID);
  intval   = cdoOperatorIntval(operatorID);

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  operatorInputArg(cdoOperatorEnter(operatorID));

  if ( operatorID == SELSEAS )
    {
      char Seas[3];
      int seas[4] = {FALSE, FALSE, FALSE, FALSE};
      int ival[4] = {0, 0, 0, 0};
      size_t len;

      nsel = operatorArgc();
      if ( isdigit(*operatorArgv()[0]))
	for ( i = 0; i < nsel; i++ )
	  {
	    if      ( atoi(operatorArgv()[i]) == 1 || atoi(operatorArgv()[i]) == 13 ) seas[0] = TRUE;
	    else if ( atoi(operatorArgv()[i]) == 2 || atoi(operatorArgv()[i]) == 14 ) seas[1] = TRUE;
	    else if ( atoi(operatorArgv()[i]) == 3 || atoi(operatorArgv()[i]) == 15 ) seas[2] = TRUE;
	    else if ( atoi(operatorArgv()[i]) == 4 || atoi(operatorArgv()[i]) == 16 ) seas[3] = TRUE;
	  }
      else
	for ( i = 0; i < nsel; i++ )
	  {
	    len = strlen(operatorArgv()[i]);
	    if ( len > 3 ) len = 3;
	    while ( len-- > 0 ) Seas[len] = toupper(operatorArgv()[i][len]);
	    if      ( Seas[0] == 'D' && Seas[1] == 'J'  && Seas[2] == 'F' ) seas[0] = TRUE;
	    else if ( Seas[0] == 'M' && Seas[1] == 'A'  && Seas[2] == 'M' ) seas[1] = TRUE;
	    else if ( Seas[0] == 'J' && Seas[1] == 'J'  && Seas[2] == 'A' ) seas[2] = TRUE;
	    else if ( Seas[0] == 'S' && Seas[1] == 'O'  && Seas[2] == 'N' ) seas[3] = TRUE;
	  }

      nsel = 4;

      if ( seas[0] ) { ival[0] = 12; ival[1] =  1; ival[2] =  2; ival[3] = 13; }
      if ( seas[1] ) { ival[0] =  3; ival[1] =  4; ival[2] =  5; ival[3] = 14; }
      if ( seas[2] ) { ival[0] =  6; ival[1] =  7; ival[2] =  8; ival[3] = 15; }
      if ( seas[3] ) { ival[0] =  9; ival[1] = 10; ival[2] = 11; ival[3] = 16; }

      listSetInt(ilist, 0, ival[0]);
      listSetInt(ilist, 1, ival[1]);
      listSetInt(ilist, 2, ival[2]);
      listSetInt(ilist, 3, ival[3]);
    }
  else if ( operatorID == SELDATE )
    {
      nsel = operatorArgc();
      if ( nsel < 1 ) cdoAbort("Not enough arguments!");
      for ( i = 0; i < nsel; i++)
	{
	  if ( strchr(operatorArgv()[i], '-') == NULL )
	    {
	      listSetInt(ilist, i, atoi(operatorArgv()[i]));
	    }
	  else
	    {
	      sscanf(operatorArgv()[i], "%d-%d-%d", &year, &month, &day);
	      listSetInt(ilist, i, year*10000 + month*100 + day);
	    }
	}
    }
  else if ( operatorID == SELTIME )
    {
      nsel = operatorArgc();
      if ( nsel < 1 ) cdoAbort("Not enough arguments!");
      for ( i = 0; i < nsel; i++ )
	{
	  if ( strchr(operatorArgv()[i], ':') == NULL )
	    {
	      listSetInt(ilist, i, atoi(operatorArgv()[i]));
	    }
	  else
	    {
	      sscanf(operatorArgv()[i], "%d:%d", &hour, &minute);
	      listSetInt(ilist, i, hour*100 + minute);
	    }
	}
    }
  else
    nsel = args2intlist(operatorArgc(), operatorArgv(), ilist);

  intarr = (int *) listArrayPtr(ilist);

  if ( cdoVerbose )
    {
      for ( i = 0; i < nsel; i++ )
	cdoPrint("intarr entry: %d %d", i, intarr[i]);
    }

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);
  if ( nsel == 1 && operfunc == func_step )  vlistDefNtsteps(vlistID2, 1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(vlistID1);
      array = (double *) malloc(gridsize*sizeof(double));
    }

  tsID  = 0;
  tsID2 = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      copytimestep = FALSE;
      selval = -1;

      if ( operfunc == func_step )
	{
	  selval = tsID + 1;
	  if ( selval > intarr[nsel-1] ) break;
	}
      else if ( operfunc == func_date )
	{
	  selval = (vdate/intval)%moddat[operatorID];
	}
      else if ( operfunc == func_time )
	{
	  selval = (vtime/intval)%moddat[operatorID];
	}

      if ( operatorID == SELDATE )
	{
	  if ( selval >= intarr[0] && selval <= intarr[nsel-1] ) copytimestep = TRUE;
	}
      else
	{
	  for ( i = 0; i < nsel; i++ )
	    if ( selval == intarr[i] )
	      {
		copytimestep = TRUE;
		break;
	      }
	}

      if ( copytimestep )
	{
	  taxisCopyTimestep(taxisID2, taxisID1);

	  streamDefTimestep(streamID2, tsID2++);

	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      streamInqRecord(streamID1, &varID, &levelID);
	      streamDefRecord(streamID2, varID, levelID);
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
	}
       
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);
 
  if ( ! lcopy )
    if ( array ) free(array);

  listDelete(ilist);

  cdoFinish();

  return (NULL);
}

/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2007 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Yhourstat   yhourmin         Multi-year hourly minimum
      Yhourstat   yhourmax         Multi-year hourly maximum
      Yhourstat   yhoursum         Multi-year hourly sum
      Yhourstat   yhourmean        Multi-year hourly mean
      Yhourstat   yhouravg         Multi-year hourly average
      Yhourstat   yhourvar         Multi-year hourly variance
      Yhourstat   yhourstd         Multi-year hourly standard deviation
*/


#include <stdio.h>
#include <math.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "field.h"


#define  NHOUR       8952  /* 31*12*24 */


void *Yhourstat(void *argument)
{
  static char func[] = "Yhourstat";
  int operatorID;
  int operfunc;
  int gridsize;
  int varID;
  int recID;
  int gridID;
  int vdate, vtime;
  int year, month, day, houroy;
  int hour, minute;
  int nrecs, nrecords;
  int levelID;
  int tsID;
  int otsID;
  long nsets[NHOUR];
  int streamID1, streamID2;
  int vlistID1, vlistID2, taxisID1, taxisID2;
  int nmiss;
  int nvars, nlevel;
  int *recVarID, *recLevelID;
  int vdates[NHOUR], vtimes[NHOUR];
  double missval;
  FIELD **vars1[NHOUR], **vars2[NHOUR];
  FIELD field;

  cdoInitialize(argument);

  cdoOperatorAdd("yhourmin",  func_min,  0, NULL);
  cdoOperatorAdd("yhourmax",  func_max,  0, NULL);
  cdoOperatorAdd("yhoursum",  func_sum,  0, NULL);
  cdoOperatorAdd("yhourmean", func_mean, 0, NULL);
  cdoOperatorAdd("yhouravg",  func_avg,  0, NULL);
  cdoOperatorAdd("yhourvar",  func_var,  0, NULL);
  cdoOperatorAdd("yhourstd",  func_std,  0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

  for ( houroy = 0; houroy < NHOUR; houroy++ )
    {
      vars1[houroy] = NULL;
      vars2[houroy] = NULL;
      nsets[houroy] = 0;
    }

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisCreate(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  nvars    = vlistNvars(vlistID1);
  nrecords = vlistNrecs(vlistID1);

  recVarID   = (int *) malloc(nrecords*sizeof(int));
  recLevelID = (int *) malloc(nrecords*sizeof(int));

  gridsize = vlistGridsizeMax(vlistID1);
  field.ptr = (double *) malloc(gridsize*sizeof(double));

  tsID = 0;
  otsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      if ( cdoVerbose ) cdoPrint("process timestep: %d %d %d", tsID+1, vdate, vtime);

      decode_date(vdate, &year, &month, &day);
      decode_time(vtime, &hour, &minute);

      if ( month >= 1 && month <= 12 && hour >= 0 && hour < 24 )
	houroy = ((month-1)*31 + day - 1)*24 + hour;
      else
	houroy = 0;

      if ( houroy < 0 || houroy >= NHOUR )
	cdoAbort("hour of year %d out of range (date=%d time=%d)!", houroy, vdate, vtime);

      vdates[houroy] = vdate;
      vtimes[houroy] = vtime;

      if ( vars1[houroy] == NULL )
	{
	  vars1[houroy] = (FIELD **) malloc(nvars*sizeof(FIELD *));
	  if ( operfunc == func_std || operfunc == func_var )
	    vars2[houroy] = (FIELD **) malloc(nvars*sizeof(FIELD *));

	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      gridID   = vlistInqVarGrid(vlistID1, varID);
	      gridsize = gridInqSize(gridID);
	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      missval  = vlistInqVarMissval(vlistID1, varID);

	      vars1[houroy][varID] = (FIELD *)  malloc(nlevel*sizeof(FIELD));
	      if ( operfunc == func_std || operfunc == func_var )
		vars2[houroy][varID] = (FIELD *)  malloc(nlevel*sizeof(FIELD));
	      
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  vars1[houroy][varID][levelID].grid    = gridID;
		  vars1[houroy][varID][levelID].nmiss   = 0;
		  vars1[houroy][varID][levelID].missval = missval;
		  vars1[houroy][varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
		  if ( operfunc == func_std || operfunc == func_var )
		    {
		      vars2[houroy][varID][levelID].grid    = gridID;
		      vars2[houroy][varID][levelID].nmiss   = 0;
		      vars2[houroy][varID][levelID].missval = missval;
		      vars2[houroy][varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
		    }
		}
	    }
	}

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  recVarID[recID]   = varID;
	  recLevelID[recID] = levelID;

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

	  if ( nsets[houroy] == 0 )
	    {
	      streamReadRecord(streamID1, vars1[houroy][varID][levelID].ptr, &nmiss);
	      vars1[houroy][varID][levelID].nmiss = nmiss;
	    }
	  else
	    {
	      streamReadRecord(streamID1, field.ptr, &field.nmiss);
	      field.grid    = vars1[houroy][varID][levelID].grid;
	      field.missval = vars1[houroy][varID][levelID].missval;

	      if ( operfunc == func_std || operfunc == func_var )
		{
		  farsumq(&vars2[houroy][varID][levelID], field);
		  farsum(&vars1[houroy][varID][levelID], field);
		}
	      else
		{
		  farfun(&vars1[houroy][varID][levelID], field, operfunc);
		}
	    }
	}

      if ( nsets[houroy] == 0 && (operfunc == func_std || operfunc == func_var) )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	    gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      farmoq(&vars2[houroy][varID][levelID], vars1[houroy][varID][levelID]);
	  }

      nsets[houroy]++;
      tsID++;
    }

  for ( houroy = 0; houroy < NHOUR; houroy++ )
    if ( nsets[houroy] )
      {
	if ( operfunc == func_mean || operfunc == func_avg )
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		farcmul(&vars1[houroy][varID][levelID], 1.0/nsets[houroy]);
	    }
	else if ( operfunc == func_std || operfunc == func_var )
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		if ( operfunc == func_std )
		  farcstd(&vars1[houroy][varID][levelID], vars2[houroy][varID][levelID], 1.0/nsets[houroy]);
		else
		  farcvar(&vars1[houroy][varID][levelID], vars2[houroy][varID][levelID], 1.0/nsets[houroy]);
	    }

	taxisDefVdate(taxisID2, vdates[houroy]);
	taxisDefVtime(taxisID2, vtimes[houroy]);
	streamDefTimestep(streamID2, otsID++);

	for ( recID = 0; recID < nrecords; recID++ )
	  {
	    varID    = recVarID[recID];
	    levelID  = recLevelID[recID];

	    if ( otsID == 1 || vlistInqVarTime(vlistID1, varID) == TIME_VARIABLE )
	      {
		streamDefRecord(streamID2, varID, levelID);
		streamWriteRecord(streamID2, vars1[houroy][varID][levelID].ptr,
				  vars1[houroy][varID][levelID].nmiss);
	      }
	  }
      }

  for ( houroy = 0; houroy < NHOUR; houroy++ )
    {
      if ( vars1[houroy] != NULL )
	{
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  free(vars1[houroy][varID][levelID].ptr);
		  if ( operfunc == func_std || operfunc == func_var ) free(vars2[houroy][varID][levelID].ptr);
		}
	      
	      free(vars1[houroy][varID]);
	      if ( operfunc == func_std || operfunc == func_var ) free(vars2[houroy][varID]);
	    }

	  free(vars1[houroy]);
	  if ( operfunc == func_std || operfunc == func_var ) free(vars2[houroy]);
	}
    }

  if ( field.ptr ) free(field.ptr);

  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}

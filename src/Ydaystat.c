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

/*
   This module contains the following operators:

      Ydaystat   ydaymin         Multi-year daily minimum
      Ydaystat   ydaymax         Multi-year daily maximum
      Ydaystat   ydaysum         Multi-year daily sum
      Ydaystat   ydaymean        Multi-year daily mean
      Ydaystat   ydayavg         Multi-year daily average
      Ydaystat   ydaystd         Multi-year daily standard deviation
*/


#include <stdio.h>
#include <math.h>

#include "cdi.h"
#include "cdo.h"
#include "pstream.h"
#include "functs.h"
#include "field.h"
#include "dmemory.h"


#define  NDAY       373


void *Ydaystat(void *argument)
{
  static char func[] = "Ydaystat";
  int operatorID;
  int operfunc;
  int gridsize;
  int varID;
  int recID;
  int gridID;
  int vdate, vtime;
  int year, month, day, dayoy;
  int nrecs, nrecords;
  int levelID;
  int tsID;
  int otsID;
  long nsets[NDAY];
  int streamID1, streamID2;
  int vlistID1, vlistID2, taxisID1, taxisID2;
  int nmiss;
  int nvars, nlevel;
  int *recVarID, *recLevelID;
  int vdates[NDAY], vtimes[NDAY];
  double missval;
  FIELD **vars1[NDAY], **vars2[NDAY];
  FIELD field;

  cdoInitialize(argument);

  cdoOperatorAdd("ydaymin",  func_min,  0, NULL);
  cdoOperatorAdd("ydaymax",  func_max,  0, NULL);
  cdoOperatorAdd("ydaysum",  func_sum,  0, NULL);
  cdoOperatorAdd("ydaymean", func_mean, 0, NULL);
  cdoOperatorAdd("ydayavg",  func_avg,  0, NULL);
  cdoOperatorAdd("ydaystd",  func_std,  0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

  for ( dayoy = 0; dayoy < NDAY; dayoy++ )
    {
      vars1[dayoy] = NULL;
      vars2[dayoy] = NULL;
      nsets[dayoy] = 0;
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

      year  =  vdate / 10000;
      month = (vdate - year*10000) / 100;
      day   =  vdate - year*10000 - month*100;

      if ( month >= 1 && month <= 12 )
	dayoy = (month-1)*31 + day;
      else
	dayoy = 0;

      if ( dayoy < 0 || dayoy >= NDAY )
	cdoAbort("day %d out of range!", dayoy);

      vdates[dayoy] = vdate;
      vtimes[dayoy] = vtime;

      if ( vars1[dayoy] == NULL )
	{
	  vars1[dayoy] = (FIELD **) malloc(nvars*sizeof(FIELD *));
	  if ( operfunc == func_std )
	    vars2[dayoy] = (FIELD **) malloc(nvars*sizeof(FIELD *));

	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      gridID   = vlistInqVarGrid(vlistID1, varID);
	      gridsize = gridInqSize(gridID);
	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      missval  = vlistInqVarMissval(vlistID1, varID);

	      vars1[dayoy][varID] = (FIELD *)  malloc(nlevel*sizeof(FIELD));
	      if ( operfunc == func_std )
		vars2[dayoy][varID] = (FIELD *)  malloc(nlevel*sizeof(FIELD));
	      
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  vars1[dayoy][varID][levelID].grid    = gridID;
		  vars1[dayoy][varID][levelID].nmiss   = 0;
		  vars1[dayoy][varID][levelID].missval = missval;
		  vars1[dayoy][varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
		  if ( operfunc == func_std )
		    {
		      vars2[dayoy][varID][levelID].grid    = gridID;
		      vars2[dayoy][varID][levelID].nmiss   = 0;
		      vars2[dayoy][varID][levelID].missval = missval;
		      vars2[dayoy][varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
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

	  if ( nsets[dayoy] == 0 )
	    {
	      streamReadRecord(streamID1, vars1[dayoy][varID][levelID].ptr, &nmiss);
	      vars1[dayoy][varID][levelID].nmiss = nmiss;
	    }
	  else
	    {
	      streamReadRecord(streamID1, field.ptr, &field.nmiss);
	      field.grid    = vars1[dayoy][varID][levelID].grid;
	      field.missval = vars1[dayoy][varID][levelID].missval;

	      if ( operfunc == func_std )
		{
		  farsumq(&vars2[dayoy][varID][levelID], field);
		  farsum(&vars1[dayoy][varID][levelID], field);
		}
	      else
		{
		  farfun(&vars1[dayoy][varID][levelID], field, operfunc);
		}
	    }
	}

      if ( nsets[dayoy] == 0 && operfunc == func_std )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	    gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      farmoq(&vars2[dayoy][varID][levelID], vars1[dayoy][varID][levelID]);
	  }

      nsets[dayoy]++;
      tsID++;
    }

  for ( dayoy = 0; dayoy < NDAY; dayoy++ )
    if ( nsets[dayoy] )
      {
	if ( operfunc == func_mean || operfunc == func_avg )
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		farcmul(&vars1[dayoy][varID][levelID], 1.0/nsets[dayoy]);
	    }
	else if ( operfunc == func_std )
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		farcstd(&vars1[dayoy][varID][levelID], vars2[dayoy][varID][levelID], 1.0/nsets[dayoy]);
	    }

	taxisDefVdate(taxisID2, vdates[dayoy]);
	taxisDefVtime(taxisID2, vtimes[dayoy]);
	streamDefTimestep(streamID2, otsID++);

	for ( recID = 0; recID < nrecords; recID++ )
	  {
	    varID    = recVarID[recID];
	    levelID  = recLevelID[recID];

	    if ( otsID == 1 || vlistInqVarTime(vlistID1, varID) == TIME_VARIABLE )
	      {
		streamDefRecord(streamID2, varID, levelID);
		streamWriteRecord(streamID2, vars1[dayoy][varID][levelID].ptr,
				  vars1[dayoy][varID][levelID].nmiss);
	      }
	  }
      }

  for ( dayoy = 0; dayoy < NDAY; dayoy++ )
    {
      if ( vars1[dayoy] != NULL )
	{
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  free(vars1[dayoy][varID][levelID].ptr);
		  if ( operfunc == func_std ) free(vars2[dayoy][varID][levelID].ptr);
		}
	      
	      free(vars1[dayoy][varID]);
	      if ( operfunc == func_std ) free(vars2[dayoy][varID]);
	    }
	  /* RQ */
	  free(vars1[dayoy]);
	  if ( operfunc == func_std ) free(vars2[dayoy]);
	  /* QR */ 
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

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

      Yhourstat   yhourmin         Multi-year hourly minimum
      Yhourstat   yhourmax         Multi-year hourly maximum
      Yhourstat   yhoursum         Multi-year hourly sum
      Yhourstat   yhourmean        Multi-year hourly mean
      Yhourstat   yhouravg         Multi-year hourly average
      Yhourstat   yhourvar         Multi-year hourly variance
      Yhourstat   yhourvar1        Multi-year hourly variance [Normalize by (n-1)]
      Yhourstat   yhourstd         Multi-year hourly standard deviation
      Yhourstat   yhourstd1        Multi-year hourly standard deviation [Normalize by (n-1)]
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  MAX_HOUR  9301  /* 31*12*25 + 1 */

static
int hour_of_year(int vdate, int vtime)
{
  int year, month, day, houroy;
  int hour, minute, second;

  cdiDecodeDate(vdate, &year, &month, &day);
  cdiDecodeTime(vtime, &hour, &minute, &second);
      
  if ( month >= 1 && month <= 12 && day >= 1 && day <=31 && hour >= 0 && hour < 24 )
    houroy = ((month-1)*31 + day - 1)*25 + hour + 1;
  else
    houroy = 0;

  if ( houroy < 0 || houroy >= MAX_HOUR )
    {
      char vdatestr[32], vtimestr[32];
      date2str(vdate, vdatestr, sizeof(vdatestr));
      time2str(vtime, vtimestr, sizeof(vtimestr));
      cdoAbort("Hour of year %d out of range (%s %s)!", houroy, vdatestr, vtimestr);
    }

  return houroy;
}


void *Yhourstat(void *argument)
{
  int i;
  int varID;
  int vdate, vtime;
  int houroy;
  int nrecs;
  int levelID;
  int nsets[MAX_HOUR];
  int nmiss;
  int nlevel;
  int vdates[MAX_HOUR], vtimes[MAX_HOUR];
  field_type **vars1[MAX_HOUR], **vars2[MAX_HOUR], **samp1[MAX_HOUR];

  cdoInitialize(argument);

  cdoOperatorAdd("yhourmin",  func_min,  0, NULL);
  cdoOperatorAdd("yhourmax",  func_max,  0, NULL);
  cdoOperatorAdd("yhoursum",  func_sum,  0, NULL);
  cdoOperatorAdd("yhourmean", func_mean, 0, NULL);
  cdoOperatorAdd("yhouravg",  func_avg,  0, NULL);
  cdoOperatorAdd("yhourvar",  func_var,  0, NULL);
  cdoOperatorAdd("yhourvar1", func_var1, 0, NULL);
  cdoOperatorAdd("yhourstd",  func_std,  0, NULL);
  cdoOperatorAdd("yhourstd1", func_std1, 0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  bool lmean   = operfunc == func_mean || operfunc == func_avg;
  bool lstd    = operfunc == func_std || operfunc == func_std1;
  bool lvarstd = operfunc == func_std || operfunc == func_var || operfunc == func_std1 || operfunc == func_var1;
  int divisor = operfunc == func_std1 || operfunc == func_var1;

  for ( houroy = 0; houroy < MAX_HOUR; ++houroy )
    {
      vars1[houroy] = NULL;
      vars2[houroy] = NULL;
      samp1[houroy] = NULL;
      nsets[houroy] = 0;
    }

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  if ( taxisHasBounds(taxisID2) ) taxisDeleteBounds(taxisID2);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  int nvars    = vlistNvars(vlistID1);
  int nrecords = vlistNrecs(vlistID1);

  int *recVarID   = (int*) Malloc(nrecords*sizeof(int));
  int *recLevelID = (int*) Malloc(nrecords*sizeof(int));

  int gridsize = vlistGridsizeMax(vlistID1);

  field_type field;
  field_init(&field);
  field.ptr = (double*) Malloc(gridsize*sizeof(double));

  int tsID = 0;
  int otsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      if ( cdoVerbose ) cdoPrint("process timestep: %d %d %d", tsID+1, vdate, vtime);

      houroy = hour_of_year(vdate, vtime);

      vdates[houroy] = vdate;
      vtimes[houroy] = vtime;

      if ( vars1[houroy] == NULL )
	{
	  vars1[houroy] = field_malloc(vlistID1, FIELD_PTR);
	  samp1[houroy] = field_malloc(vlistID1, FIELD_NONE);
	  if ( lvarstd )
	    vars2[houroy] = field_malloc(vlistID1, FIELD_PTR);
	}

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  if ( tsID == 0 )
	    {
	      recVarID[recID]   = varID;
	      recLevelID[recID] = levelID;
	    }

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

	  if ( nsets[houroy] == 0 )
	    {
	      streamReadRecord(streamID1, vars1[houroy][varID][levelID].ptr, &nmiss);
	      vars1[houroy][varID][levelID].nmiss = (size_t) nmiss;

	      if ( nmiss > 0 || samp1[houroy][varID][levelID].ptr )
		{
		  if ( samp1[houroy][varID][levelID].ptr == NULL )
		    samp1[houroy][varID][levelID].ptr = (double*) Malloc(gridsize*sizeof(double));

		  for ( i = 0; i < gridsize; i++ )
		    if ( DBL_IS_EQUAL(vars1[houroy][varID][levelID].ptr[i],
				      vars1[houroy][varID][levelID].missval) )
		      samp1[houroy][varID][levelID].ptr[i] = 0;
		    else
		      samp1[houroy][varID][levelID].ptr[i] = 1;
		}
	    }
	  else
	    {
	      streamReadRecord(streamID1, field.ptr, &nmiss);
              field.nmiss   = (size_t) nmiss;
	      field.grid    = vars1[houroy][varID][levelID].grid;
	      field.missval = vars1[houroy][varID][levelID].missval;

	      if ( field.nmiss > 0 || samp1[houroy][varID][levelID].ptr )
		{
		  if ( samp1[houroy][varID][levelID].ptr == NULL )
		    {
		      samp1[houroy][varID][levelID].ptr = (double*) Malloc(gridsize*sizeof(double));
		      for ( i = 0; i < gridsize; i++ )
			samp1[houroy][varID][levelID].ptr[i] = nsets[houroy];
		    }
		  
		  for ( i = 0; i < gridsize; i++ )
		    if ( !DBL_IS_EQUAL(field.ptr[i], vars1[houroy][varID][levelID].missval) )
		      samp1[houroy][varID][levelID].ptr[i]++;
		}

	      if ( lvarstd )
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

      if ( nsets[houroy] == 0 && lvarstd )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      farmoq(&vars2[houroy][varID][levelID], vars1[houroy][varID][levelID]);
	  }

      nsets[houroy]++;
      tsID++;
    }

  for ( houroy = 0; houroy < MAX_HOUR; ++houroy )
    if ( nsets[houroy] )
      {
	if ( lmean )
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  if ( samp1[houroy][varID][levelID].ptr == NULL )
		    farcdiv(&vars1[houroy][varID][levelID], (double)nsets[houroy]);
		  else
		    fardiv(&vars1[houroy][varID][levelID], samp1[houroy][varID][levelID]);
		}
	    }
	else if ( lvarstd )
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  if ( samp1[houroy][varID][levelID].ptr == NULL )
		    {
		      if ( lstd )
			farcstd(&vars1[houroy][varID][levelID], vars2[houroy][varID][levelID], nsets[houroy], divisor);
		      else
			farcvar(&vars1[houroy][varID][levelID], vars2[houroy][varID][levelID], nsets[houroy], divisor);
		    }
		  else
		    {
		      if ( lstd )
			farstd(&vars1[houroy][varID][levelID], vars2[houroy][varID][levelID], samp1[houroy][varID][levelID], divisor);
		      else
			farvar(&vars1[houroy][varID][levelID], vars2[houroy][varID][levelID], samp1[houroy][varID][levelID], divisor);
		    }
		}
	    }

	taxisDefVdate(taxisID2, vdates[houroy]);
	taxisDefVtime(taxisID2, vtimes[houroy]);
	streamDefTimestep(streamID2, otsID);

	for ( int recID = 0; recID < nrecords; recID++ )
	  {
	    varID    = recVarID[recID];
	    levelID  = recLevelID[recID];

	    if ( otsID && vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

	    streamDefRecord(streamID2, varID, levelID);
	    streamWriteRecord(streamID2, vars1[houroy][varID][levelID].ptr,
			      (int)vars1[houroy][varID][levelID].nmiss);
	  }

	otsID++;
      }

  for ( houroy = 0; houroy < MAX_HOUR; ++houroy )
    {
      if ( vars1[houroy] != NULL )
	{
	  field_free(samp1[houroy], vlistID1);
	  field_free(vars1[houroy], vlistID1);
	  if ( lvarstd ) field_free(vars2[houroy], vlistID1);
	}
    }

  if ( field.ptr ) Free(field.ptr);

  if ( recVarID   ) Free(recVarID);
  if ( recLevelID ) Free(recLevelID);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return 0;
}

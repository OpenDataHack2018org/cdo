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

      Ymonstat   ymonmin         Multi-year monthly minimum
      Ymonstat   ymonmax         Multi-year monthly maximum
      Ymonstat   ymonsum         Multi-year monthly sum
      Ymonstat   ymonmean        Multi-year monthly mean
      Ymonstat   ymonavg         Multi-year monthly average
      Ymonstat   ymonvar         Multi-year monthly variance
      Ymonstat   ymonvar1        Multi-year monthly variance [Normalize by (n-1)]
      Ymonstat   ymonstd         Multi-year monthly standard deviation
      Ymonstat   ymonstd1        Multi-year monthly standard deviation [Normalize by (n-1)]
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  NMONTH     17

/*
static
int cmpint(const void *s1, const void *s2)
{
  int cmp = 0;
  const int *x = s1;
  const int *y = s2;

  if      ( *x < *y ) cmp = -1;
  else if ( *x > *y ) cmp =  1;

  return cmp;
}
*/

void *Ymonstat(void *argument)
{
  int i;
  int varID;
  int vdate, vtime;
  int year, month, day;
  int nrecs;
  int levelID;
  int nsets[NMONTH];
  int nmiss;
  int nlevel;
  int vdates[NMONTH], vtimes[NMONTH];
  int mon[NMONTH];
  int nmon = 0;
  field_t **vars1[NMONTH], **vars2[NMONTH], **samp1[NMONTH];
  field_t field;

  cdoInitialize(argument);

  cdoOperatorAdd("ymonmin",  func_min,  0, NULL);
  cdoOperatorAdd("ymonmax",  func_max,  0, NULL);
  cdoOperatorAdd("ymonsum",  func_sum,  0, NULL);
  cdoOperatorAdd("ymonmean", func_mean, 0, NULL);
  cdoOperatorAdd("ymonavg",  func_avg,  0, NULL);
  cdoOperatorAdd("ymonvar",  func_var,  0, NULL);
  cdoOperatorAdd("ymonvar1", func_var1, 0, NULL);
  cdoOperatorAdd("ymonstd",  func_std,  0, NULL);
  cdoOperatorAdd("ymonstd1", func_std1, 0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  bool lmean   = operfunc == func_mean || operfunc == func_avg;
  bool lstd    = operfunc == func_std || operfunc == func_std1;
  bool lvarstd = operfunc == func_std || operfunc == func_var || operfunc == func_std1 || operfunc == func_var1;
  int divisor  = operfunc == func_std1 || operfunc == func_var1;

  for ( month = 0; month < NMONTH; month++ )
    {
      vars1[month] = NULL;
      vars2[month] = NULL;
      samp1[month] = NULL;
      nsets[month] = 0;
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
  field_init(&field);
  field.ptr = (double*) Malloc(gridsize*sizeof(double));

  int tsID = 0;
  int otsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      if ( cdoVerbose ) cdoPrint("process timestep: %d %d %d", tsID+1, vdate, vtime);

      cdiDecodeDate(vdate, &year, &month, &day);
      if ( month < 0 || month >= NMONTH )
	cdoAbort("month %d out of range!", month);

      vdates[month] = vdate;
      vtimes[month] = vtime;
      // mon[month] = vdate;

      if ( vars1[month] == NULL )
	{
	  mon[nmon++] = month;
	  vars1[month] = field_malloc(vlistID1, FIELD_PTR);
	  samp1[month] = field_malloc(vlistID1, FIELD_NONE);
	  if ( lvarstd )
	    vars2[month] = field_malloc(vlistID1, FIELD_PTR);
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

	  if ( nsets[month] == 0 )
	    {
	      streamReadRecord(streamID1, vars1[month][varID][levelID].ptr, &nmiss);
	      vars1[month][varID][levelID].nmiss = (size_t) nmiss;

	      if ( nmiss > 0 || samp1[month][varID][levelID].ptr )
		{
		  if ( samp1[month][varID][levelID].ptr == NULL )
		    samp1[month][varID][levelID].ptr = (double*) Malloc(gridsize*sizeof(double));

		  for ( i = 0; i < gridsize; i++ )
		    if ( DBL_IS_EQUAL(vars1[month][varID][levelID].ptr[i],
				      vars1[month][varID][levelID].missval) )
		      samp1[month][varID][levelID].ptr[i] = 0;
		    else
		      samp1[month][varID][levelID].ptr[i] = 1;
		}
	    }
	  else
	    {
	      streamReadRecord(streamID1, field.ptr, &nmiss);
              field.nmiss   = (size_t) nmiss;
	      field.grid    = vars1[month][varID][levelID].grid;
	      field.missval = vars1[month][varID][levelID].missval;

	      if ( field.nmiss > 0 || samp1[month][varID][levelID].ptr )
		{
		  if ( samp1[month][varID][levelID].ptr == NULL )
		    {
		      samp1[month][varID][levelID].ptr = (double*) Malloc(gridsize*sizeof(double));
		      for ( i = 0; i < gridsize; i++ )
			samp1[month][varID][levelID].ptr[i] = nsets[month];
		    }
		  
		  for ( i = 0; i < gridsize; i++ )
		    if ( !DBL_IS_EQUAL(field.ptr[i], vars1[month][varID][levelID].missval) )
		      samp1[month][varID][levelID].ptr[i]++;
		}

	      if ( lvarstd )
		{
		  farsumq(&vars2[month][varID][levelID], field);
		  farsum(&vars1[month][varID][levelID], field);
		}
	      else
		{
		  farfun(&vars1[month][varID][levelID], field, operfunc);
		}
	    }
	}

      if ( nsets[month] == 0 && lvarstd )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      farmoq(&vars2[month][varID][levelID], vars1[month][varID][levelID]);
	  }

      nsets[month]++;
      tsID++;
    }

  if ( nmon == 12 )
    {
      int smon = 0;
      for ( month = 1; month <= 12; month++ ) if ( nsets[month] ) smon++;
      if ( smon == 12 ) for ( month = 1; month <= 12; month++ ) mon[month-1] = month;
    }

  /* sort output time steps */
  /*
  nmon = 0;
  for ( month = 0; month < NMONTH; month++ )
    {
      if ( nsets[month] == 0 )
	for ( i = month+1; i < NMONTH; i++ ) mon[i-1] = mon[i];
      else
	nmon++;
    }

  qsort(mon, nmon, sizeof(int), cmpint);
	      
  for ( i = 0; i < nmon; i++ )
    {
      cdiDecodeDate(mon[i], &year, &month, &day);
      mon[i] = month;
    }
  */
  for ( i = 0; i < nmon; i++ )
    {
      month = mon[i];
      if ( nsets[month] == 0 ) cdoAbort("Internal problem, nsets[%d] not defined!", month);

      if ( lmean )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      {
		if ( samp1[month][varID][levelID].ptr == NULL )
		  farcdiv(&vars1[month][varID][levelID], (double)nsets[month]);
		else
		  fardiv(&vars1[month][varID][levelID], samp1[month][varID][levelID]);
	      }
	  }
      else if ( lvarstd )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      {
		if ( samp1[month][varID][levelID].ptr == NULL )
		  {
		    if ( lstd )
		      farcstd(&vars1[month][varID][levelID], vars2[month][varID][levelID], nsets[month], divisor);
		    else
		      farcvar(&vars1[month][varID][levelID], vars2[month][varID][levelID], nsets[month], divisor);
		  }
		else
		  {
		    if ( lstd )
		      farstd(&vars1[month][varID][levelID], vars2[month][varID][levelID], samp1[month][varID][levelID], divisor);
		    else
		      farvar(&vars1[month][varID][levelID], vars2[month][varID][levelID], samp1[month][varID][levelID], divisor);
		  }
	      }
	  }

      taxisDefVdate(taxisID2, vdates[month]);
      taxisDefVtime(taxisID2, vtimes[month]);
      streamDefTimestep(streamID2, otsID);

      for ( int recID = 0; recID < nrecords; recID++ )
	{
	  varID    = recVarID[recID];
	  levelID  = recLevelID[recID];
	  
	  if ( otsID && vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, vars1[month][varID][levelID].ptr,
			    (int)vars1[month][varID][levelID].nmiss);
	}

      otsID++;
    }

  for ( month = 0; month < NMONTH; month++ )
    {
      if ( vars1[month] != NULL )
	{
	  field_free(vars1[month], vlistID1);
	  field_free(samp1[month], vlistID1);
	  if ( lvarstd ) field_free(vars2[month], vlistID1);
	}
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( field.ptr ) Free(field.ptr);

  if ( recVarID   ) Free(recVarID);
  if ( recLevelID ) Free(recLevelID);

  cdoFinish();

  return 0;
}

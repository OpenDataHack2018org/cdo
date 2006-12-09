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

      Timstat    timmin          Time minimum
      Timstat    timmax          Time maximum
      Timstat    timsum          Time sum
      Timstat    timmean         Time mean
      Timstat    timavg          Time average
      Timstat    timvar          Time variance
      Timstat    timstd          Time standard deviation
      Hourstat   hourmin         Hourly minimum
      Hourstat   hourmax         Hourly maximum
      Hourstat   hoursum         Hourly sum
      Hourstat   hourmean        Hourly mean
      Hourstat   houravg         Hourly average
      Hourstat   hourvar         Hourly variance
      Hourstat   hourstd         Hourly standard deviation
      Daystat    daymin          Daily minimum
      Daystat    daymax          Daily maximum
      Daystat    daysum          Daily sum
      Daystat    daymean         Daily mean
      Daystat    dayavg          Daily average
      Daystat    dayvar          Daily variance
      Daystat    daystd          Daily standard deviation
      Monstat    monmin          Monthly minimum
      Monstat    monmax          Monthly maximum
      Monstat    monsum          Monthly sum
      Monstat    monmean         Monthly mean
      Monstat    monavg          Monthly average
      Monstat    monvar          Monthly variance
      Monstat    monstd          Monthly standard deviation
      Yearstat   yearmin         Yearly minimum
      Yearstat   yearmax         Yearly maximum
      Yearstat   yearsum         Yearly sum
      Yearstat   yearmean        Yearly mean
      Yearstat   yearavg         Yearly average
      Yearstat   yearvar         Yearly variance
      Yearstat   yearstd         Yearly standard deviation
*/


#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "dtypes.h"


void *Timstat(void *argument)
{
  static char func[] = "Timstat";
  int operatorID;
  int operfunc;
  INT64 intvdat;
  INT64 indate1, indate2 = 0;
  int gridsize;
  int vdate = 0, vtime = 0;
  int vdate0 = 0, vtime0 = 0;
  int nrecs, nrecords;
  int gridID, varID, levelID, recID;
  int tsID;
  int otsID;
  long nsets;
  int i;
  int streamID1, streamID2;
  int vlistID1, vlistID2, taxisID1, taxisID2;
  int nmiss;
  int nvars, nlevel;
  int *recVarID, *recLevelID;
  double missval;
  FIELD **vars1 = NULL, **vars2 = NULL, **samp1 = NULL;
  FIELD field;

  cdoInitialize(argument);

  cdoOperatorAdd("timmin",   func_min,  17, NULL);
  cdoOperatorAdd("timmax",   func_max,  17, NULL);
  cdoOperatorAdd("timsum",   func_sum,  17, NULL);
  cdoOperatorAdd("timmean",  func_mean, 17, NULL);
  cdoOperatorAdd("timavg",   func_avg,  17, NULL);
  cdoOperatorAdd("timvar",   func_var,  17, NULL);
  cdoOperatorAdd("timstd",   func_std,  17, NULL);
  cdoOperatorAdd("yearmin",  func_min,   8, NULL);
  cdoOperatorAdd("yearmax",  func_max,   8, NULL);
  cdoOperatorAdd("yearsum",  func_sum,   8, NULL);
  cdoOperatorAdd("yearmean", func_mean,  8, NULL);
  cdoOperatorAdd("yearavg",  func_avg,   8, NULL);
  cdoOperatorAdd("yearvar",  func_var,   8, NULL);
  cdoOperatorAdd("yearstd",  func_std,   8, NULL);
  cdoOperatorAdd("monmin",   func_min,   6, NULL);
  cdoOperatorAdd("monmax",   func_max,   6, NULL);
  cdoOperatorAdd("monsum",   func_sum,   6, NULL);
  cdoOperatorAdd("monmean",  func_mean,  6, NULL);
  cdoOperatorAdd("monavg",   func_avg,   6, NULL);
  cdoOperatorAdd("monvar",   func_var,   6, NULL);
  cdoOperatorAdd("monstd",   func_std,   6, NULL);
  cdoOperatorAdd("daymin",   func_min,   4, NULL);
  cdoOperatorAdd("daymax",   func_max,   4, NULL);
  cdoOperatorAdd("daysum",   func_sum,   4, NULL);
  cdoOperatorAdd("daymean",  func_mean,  4, NULL);
  cdoOperatorAdd("dayavg",   func_avg,   4, NULL);
  cdoOperatorAdd("dayvar",   func_var,   4, NULL);
  cdoOperatorAdd("daystd",   func_std,   4, NULL);
  cdoOperatorAdd("hourmin",  func_min,   2, NULL);
  cdoOperatorAdd("hourmax",  func_max,   2, NULL);
  cdoOperatorAdd("hoursum",  func_sum,   2, NULL);
  cdoOperatorAdd("hourmean", func_mean,  2, NULL);
  cdoOperatorAdd("houravg",  func_avg,   2, NULL);
  cdoOperatorAdd("hourvar",  func_var,   2, NULL);
  cdoOperatorAdd("hourstd",  func_std,   2, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

  intvdat = (INT64) pow(10.0, (double) cdoOperatorIntval(operatorID));

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  if ( cdoOperatorIntval(operatorID) == 17 ) vlistDefNtsteps(vlistID2, 1);

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

  vars1 = (FIELD **) malloc(nvars*sizeof(FIELD *));
  samp1 = (FIELD **) malloc(nvars*sizeof(FIELD *));
  if ( operfunc == func_std || operfunc == func_var )
    vars2 = (FIELD **) malloc(nvars*sizeof(FIELD *));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      missval  = vlistInqVarMissval(vlistID1, varID);

      vars1[varID] = (FIELD *)  malloc(nlevel*sizeof(FIELD));
      samp1[varID] = (FIELD *)  malloc(nlevel*sizeof(FIELD));
      if ( operfunc == func_std || operfunc == func_var )
	vars2[varID] = (FIELD *)  malloc(nlevel*sizeof(FIELD));

      for ( levelID = 0; levelID < nlevel; levelID++ )
	{
	  vars1[varID][levelID].grid    = gridID;
	  vars1[varID][levelID].nmiss   = 0;
	  vars1[varID][levelID].missval = missval;
	  vars1[varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
	  samp1[varID][levelID].grid    = gridID;
	  samp1[varID][levelID].nmiss   = 0;
	  samp1[varID][levelID].missval = missval;
	  samp1[varID][levelID].ptr     = NULL;
	  if ( operfunc == func_std || operfunc == func_var )
	    {
	      vars2[varID][levelID].grid    = gridID;
	      vars2[varID][levelID].nmiss   = 0;
	      vars2[varID][levelID].missval = missval;
	      vars2[varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
	    }
	}
    }

  indate1 = 0;
  tsID    = 0;
  otsID   = 0;
  while ( TRUE )
    {
      nsets = 0;
      while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
	{
	  vdate = taxisInqVdate(taxisID1);
	  vtime = taxisInqVtime(taxisID1);

	  if ( nsets == 0 ) indate2 = (INT64)vdate*10000 + vtime;
	  indate1 = (INT64)vdate*10000 + vtime;

	  if ( indate1/intvdat != indate2/intvdat ) break;

	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      streamInqRecord(streamID1, &varID, &levelID);

	      if ( tsID == 0 )
		{
		  recVarID[recID]   = varID;
		  recLevelID[recID] = levelID;
		}

	      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

	      if ( nsets == 0 )
		{
		  streamReadRecord(streamID1, vars1[varID][levelID].ptr, &nmiss);
		  vars1[varID][levelID].nmiss = nmiss;

		  if ( nmiss > 0 || samp1[varID][levelID].ptr )
		    {
		      if ( samp1[varID][levelID].ptr == NULL )
			samp1[varID][levelID].ptr = (double *) malloc(gridsize*sizeof(double));

		      for ( i = 0; i < gridsize; i++ )
			if ( DBL_IS_EQUAL(vars1[varID][levelID].ptr[i], vars1[varID][levelID].missval) )
			  samp1[varID][levelID].ptr[i] = 0;
			else
			  samp1[varID][levelID].ptr[i] = 1;
		    }
		}
	      else
		{
		  streamReadRecord(streamID1, field.ptr, &field.nmiss);
		  field.grid    = vars1[varID][levelID].grid;
		  field.missval = vars1[varID][levelID].missval;

		  if ( field.nmiss > 0 || samp1[varID][levelID].ptr )
		    {
		      if ( samp1[varID][levelID].ptr == NULL )
			{
			  samp1[varID][levelID].ptr = (double *) malloc(gridsize*sizeof(double));
			  for ( i = 0; i < gridsize; i++ )
			    samp1[varID][levelID].ptr[i] = nsets;
			}

		      for ( i = 0; i < gridsize; i++ )
			if ( !DBL_IS_EQUAL(field.ptr[i], vars1[varID][levelID].missval) )
			  samp1[varID][levelID].ptr[i]++;
		    }

		  if ( operfunc == func_std || operfunc == func_var )
		    {
		      farsumq(&vars2[varID][levelID], field);
		      farsum(&vars1[varID][levelID], field);
		    }
		  else
		    {
		      farfun(&vars1[varID][levelID], field, operfunc);
		    }
		}
	    }

	  if ( nsets == 0 && (operfunc == func_std || operfunc == func_var) )
	    for ( varID = 0; varID < nvars; varID++ )
	      {
		if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
		nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
		for ( levelID = 0; levelID < nlevel; levelID++ )
		  farmoq(&vars2[varID][levelID], vars1[varID][levelID]);
	      }

	  vdate0 = vdate;
	  vtime0 = vtime;
	  nsets++;
	  tsID++;
	}

      if ( nrecs == 0 && nsets == 0 ) break;

      if ( operfunc == func_mean || operfunc == func_avg )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      {
		if ( samp1[varID][levelID].ptr == NULL )
		  farcmul(&vars1[varID][levelID], 1.0/nsets);
		else
		  fardiv(&vars1[varID][levelID], samp1[varID][levelID]);
	      }
	  }
      else if ( operfunc == func_std || operfunc == func_var )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      {
		if ( samp1[varID][levelID].ptr == NULL )
		  {
		    if ( operfunc == func_std )
		      farcstd(&vars1[varID][levelID], vars2[varID][levelID], 1.0/nsets);
		    else
		      farcvar(&vars1[varID][levelID], vars2[varID][levelID], 1.0/nsets);
		  }
		else
		  {
		    farinv(&samp1[varID][levelID]);
		    if ( operfunc == func_std )
		      farstd(&vars1[varID][levelID], vars2[varID][levelID], samp1[varID][levelID]);
		    else
		      farvar(&vars1[varID][levelID], vars2[varID][levelID], samp1[varID][levelID]);
		  }
	      }
	  }

      taxisDefVdate(taxisID2, vdate0);
      taxisDefVtime(taxisID2, vtime0);
      streamDefTimestep(streamID2, otsID++);

      for ( recID = 0; recID < nrecords; recID++ )
	{
	  varID   = recVarID[recID];
	  levelID = recLevelID[recID];

	  if ( otsID == 1 || vlistInqVarTime(vlistID1, varID) == TIME_VARIABLE )
	    {
	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, vars1[varID][levelID].ptr,  vars1[varID][levelID].nmiss);
	    }
	}

      if ( nrecs == 0 ) break;
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      for ( levelID = 0; levelID < nlevel; levelID++ )
	{
	  free(vars1[varID][levelID].ptr);
	  if ( samp1[varID][levelID].ptr ) free(samp1[varID][levelID].ptr);
	  if ( operfunc == func_std || operfunc == func_var ) free(vars2[varID][levelID].ptr);
	}

      free(vars1[varID]);
      free(samp1[varID]);
      if ( operfunc == func_std || operfunc == func_var ) free(vars2[varID]);
    }

  free(vars1);
  free(samp1);
  if ( operfunc == func_std || operfunc == func_var ) free(vars2);

  if ( field.ptr ) free(field.ptr);

  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}

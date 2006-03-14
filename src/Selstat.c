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

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "dtypes.h"

/*
@BeginDoc

@BeginModule

@Name      = Selstat
@Title     = 
@Section   = Statistical description of the data
@Class     = Statistic
@Arguments = ifile ofile
@Operators = selmin selmax selsum selmean selavg selstd

@EndModule


@BeginOperator_selmin

@Title     = Time range minimum
@Parameter = nsets [noffset] [nskip]

@BeginDescription
@IfMan
For every adjacent sequence t1, ...., tn of timesteps of the same 
selected time range, it is

o(t,x) = min{i(t',x), t1 < t' <= tn}
@EndifMan
@IfDoc
For every adjacent sequence \begin{math}t_1, ...,t_n\end{math} of timesteps of the same 
selected time range, it is
@BeginMath
o(t,x) = \mbox{\bf min}\{i(t',x), t_1 < t' \le t_n\}
@EndMath
@EndifDoc
@EndDescription

@BeginParameter noffset
@Item = nsets
INTEGER  Number of input timesteps for each output timestep
@Item = noffset
INTEGER  Number of input timesteps skipped before the first timestep range (optional)
@Item = nskip
INTEGER  Number of input timesteps skipped between timestep ranges (optional)
@EndParameter

@EndOperator


@BeginOperator_selmax

@Title     = Time range maximum
@Parameter = nsets [noffset] [nskip]

@BeginDescription
@IfMan
For every adjacent sequence t1, ...., tn of timesteps of the same 
selected time range, it is

o(t,x) = max{i(t',x), t1 < t' <= tn}
@EndifMan
@IfDoc
For every adjacent sequence \begin{math}t_1, ...,t_n\end{math} of timesteps of the same 
selected time range, it is
@BeginMath
o(t,x) = \mbox{\bf max}\{i(t',x), t_1 < t' \le t_n\}
@EndMath
@EndifDoc
@EndDescription

@BeginParameter noffset
@Item = nsets
INTEGER  Number of input timesteps for each output timestep
@Item = noffset
INTEGER  Number of input timesteps skipped before the first timestep range (optional)
@Item = nskip
INTEGER  Number of input timesteps skipped between timestep ranges (optional)
@EndParameter

@EndOperator


@BeginOperator_selsum

@Title     = Time range sum
@Parameter = nsets [noffset] [nskip]

@BeginDescription
@IfMan
For every adjacent sequence t1, ...., tn of timesteps of the same 
selected time range, it is

o(t,x) = sum{i(t',x), t1 < t' <= tn}
@EndifMan
@IfDoc
For every adjacent sequence \begin{math}t_1, ...,t_n\end{math} of timesteps of the same 
selected time range, it is
@BeginMath
o(t,x) = \sum\limits_{t=1}^{n} i(t',x)
@EndMath
@EndifDoc
@EndDescription

@BeginParameter noffset
@Item = nsets
INTEGER  Number of input timesteps for each output timestep
@Item = noffset
INTEGER  Number of input timesteps skipped before the first timestep range (optional)
@Item = nskip
INTEGER  Number of input timesteps skipped between timestep ranges (optional)
@EndParameter

@EndOperator


@BeginOperator_selmean

@Title     = Time range mean
@Parameter = nsets [noffset] [nskip]

@BeginDescription
@IfMan
For every adjacent sequence t1, ...., tn of timesteps of the same 
selected time range, it is

o(t,x) = mean{i(t',x), t1 < t' <= tn}
@EndifMan
@IfDoc
For every adjacent sequence \begin{math}t_1, ...,t_n\end{math} of timesteps of the same 
selected time range, it is
@BeginMath
o(t,x) = \mbox{\bf mean}\{i(t',x), t_1 < t' \le t_n\}
@EndMath
@EndifDoc
@EndDescription

@BeginParameter noffset
@Item = nsets
INTEGER  Number of input timesteps for each output timestep
@Item = noffset
INTEGER  Number of input timesteps skipped before the first timestep range (optional)
@Item = nskip
INTEGER  Number of input timesteps skipped between timestep ranges (optional)
@EndParameter

@EndOperator


@BeginOperator_selavg

@Title     = Time range average
@Parameter = nsets [noffset] [nskip]

@BeginDescription
@IfMan
For every adjacent sequence t1, ...., tn of timesteps of the same 
selected time range, it is

o(t,x) = avg{i(t',x), t1 < t' <= tn}
@EndifMan
@IfDoc
For every adjacent sequence \begin{math}t_1, ...,t_n\end{math} of timesteps of the same 
selected time range, it is
@BeginMath
o(t,x) = \mbox{\bf avg}\{i(t',x), t_1 < t' \le t_n\}
@EndMath
@EndifDoc
@EndDescription

@BeginParameter noffset
@Item = nsets
INTEGER  Number of input timesteps for each output timestep
@Item = noffset
INTEGER  Number of input timesteps skipped before the first timestep range (optional)
@Item = nskip
INTEGER  Number of input timesteps skipped between timestep ranges (optional)
@EndParameter

@EndOperator


@BeginOperator_selstd

@Title     = Time range standard deviation
@Parameter = nsets [noffset] [nskip]

@BeginDescription
@IfMan
For every adjacent sequence t1, ...., tn of timesteps of the same 
selected time range, it is

o(t,x) = sqrt{var{i(t',x), t1 < t' <= tn}}
@EndifMan
@IfDoc
For every adjacent sequence \begin{math}t_1, ...,t_n\end{math} of timesteps of the same 
selected time range, it is
@BeginMath
o(t,x) = \sqrt{\mbox{\bf var}\{i(t',x'), t_1 < t' \le t_n\}}
@EndMath
@EndifDoc
@EndDescription

@BeginParameter noffset
@Item = nsets
INTEGER  Number of input timesteps for each output timestep
@Item = noffset
INTEGER  Number of input timesteps skipped before the first timestep range (optional)
@Item = nskip
INTEGER  Number of input timesteps skipped between timestep ranges (optional)
@EndParameter

@EndOperator

@EndDoc
*/


void *Selstat(void *argument)
{
  static char func[] = "Selstat";
  int operatorID;
  int operfunc;
  int gridsize;
  int vdate = 0, vtime = 0;
  int vdate0 = 0, vtime0 = 0;
  int nrecs = 0, nrecords;
  int gridID, varID, levelID, recID;
  int tsID;
  int otsID;
  int nsets;
  int i;
  int streamID1, streamID2;
  int vlistID1, vlistID2, taxisID1, taxisID2;
  int nmiss;
  int nvars, nlevel;
  int ndates = 0, noffset = 0, nskip = 0, nargc;
  int *recVarID, *recLevelID;
  double missval;
  FIELD **vars1 = NULL, **vars2 = NULL, **samp1 = NULL;
  FIELD field;

  cdoInitialize(argument);

  cdoOperatorAdd("selmin",  func_min,  0, NULL);
  cdoOperatorAdd("selmax",  func_max,  0, NULL);
  cdoOperatorAdd("selsum",  func_sum,  0, NULL);
  cdoOperatorAdd("selmean", func_mean, 0, NULL);
  cdoOperatorAdd("selavg",  func_avg,  0, NULL);
  cdoOperatorAdd("selstd",  func_std,  0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

  operatorInputArg("nsets <noffset <nskip>>");

  nargc = operatorArgc();

  ndates = atoi(operatorArgv()[0]);
  if ( nargc > 1 ) noffset = atoi(operatorArgv()[1]);
  if ( nargc > 2 ) nskip   = atoi(operatorArgv()[2]);

  if ( cdoVerbose ) cdoPrint("nsets = %d, noffset = %d, nskip = %d", ndates, noffset, nskip);

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisNew(TAXIS_ABSOLUTE);
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
  if ( operfunc == func_std )
    vars2 = (FIELD **) malloc(nvars*sizeof(FIELD *));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      missval  = vlistInqVarMissval(vlistID1, varID);

      vars1[varID] = (FIELD *)  malloc(nlevel*sizeof(FIELD));
      samp1[varID] = (FIELD *)  malloc(nlevel*sizeof(FIELD));
      if ( operfunc == func_std )
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
	  if ( operfunc == func_std )
	    {
	      vars2[varID][levelID].grid    = gridID;
	      vars2[varID][levelID].nmiss   = 0;
	      vars2[varID][levelID].missval = missval;
	      vars2[varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
	    }
	}
    }

  for ( tsID = 0; tsID < noffset; tsID++ )
    {
      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) break;
      tsID++;
    }

  if ( tsID < noffset )
    {
      cdoWarning("noffset larger than number of timesteps!");
      goto LABEL_END;
    }

  otsID   = 0;
  while ( TRUE )
    {
      for ( nsets = 0; nsets < ndates; nsets++ )
	{
	  nrecs = streamInqTimestep(streamID1, tsID);
	  if ( nrecs == 0 ) break;

	  vdate = taxisInqVdate(taxisID1);
	  vtime = taxisInqVtime(taxisID1);

	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      streamInqRecord(streamID1, &varID, &levelID);

	      recVarID[recID]   = varID;
	      recLevelID[recID] = levelID;

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

		  if ( operfunc == func_std )
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

	  if ( nsets == 0 && operfunc == func_std )
	    for ( varID = 0; varID < nvars; varID++ )
	      {
		if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
		nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
		for ( levelID = 0; levelID < nlevel; levelID++ )
		  farmoq(&vars2[varID][levelID], vars1[varID][levelID]);
	      }

	  vdate0 = vdate;
	  vtime0 = vtime;
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
      else if ( operfunc == func_std )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      {
		if ( samp1[varID][levelID].ptr == NULL )
		  farcstd(&vars1[varID][levelID], vars2[varID][levelID], 1.0/nsets);
		else
		  {
		    farinv(&samp1[varID][levelID]);
		    farstd(&vars1[varID][levelID], vars2[varID][levelID], samp1[varID][levelID]);
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

      for ( i = 0; i < nskip; i++ )
	{
	  nrecs = streamInqTimestep(streamID1, tsID);
	  if ( nrecs == 0 ) break;
	  tsID++;
	}

      if ( nrecs == 0 ) break;
    }

 LABEL_END:

  for ( varID = 0; varID < nvars; varID++ )
    {
      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      for ( levelID = 0; levelID < nlevel; levelID++ )
	{
	  free(vars1[varID][levelID].ptr);
	  if ( samp1[varID][levelID].ptr ) free(samp1[varID][levelID].ptr);
	  if ( operfunc == func_std ) free(vars2[varID][levelID].ptr);
	}

      free(vars1[varID]);
      free(samp1[varID]);
      if ( operfunc == func_std ) free(vars2[varID]);
    }

  free(vars1);
  free(samp1);
  if ( operfunc == func_std ) free(vars2);

  if ( field.ptr ) free(field.ptr);

  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}

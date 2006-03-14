/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2005 Uwe Schulzweida, schulzweida@dkrz.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include <stdio.h>
#include <math.h>

#include "cdi.h"
#include "cdo.h"
#include "pstream.h"
#include "functs.h"
#include "field.h"
#include "dmemory.h"

/*
@BeginDoc

@BeginModule

@Name      = Ymonstat
@Title     = Multi-year monthly statistics
@Section   = Statistical description of the data
@Class     = Statistic
@Arguments = ifile ofile
@Operators = ymonmin ymonmax ymonmean ymonavg ymonstd

@EndModule


@BeginOperator_ymonmin

@Title     = Multi-year monthly minimum

@BeginDescription
@IfMan
o(01,x) = min{i(t,x), month(i(t)) = 01}
                 ...
o(12,x) = min{i(t,x), month(i(t)) = 12}
@EndifMan
@IfDoc
@BeginMath
\begin{array}{c}
o(\mbox{01},x) = \mbox{\bf min}\{i(t,x), \mbox{month}(i(t)) = \mbox{01}\} \\
\vdots \\
o(\mbox{12},x) = \mbox{\bf min}\{i(t,x), \mbox{month}(i(t)) = \mbox{12}\} \\
\end{array}
@EndMath
@EndifDoc
@EndDescription

@EndOperator


@BeginOperator_ymonmax

@Title     = Multi-year monthly maximum

@BeginDescription
@IfMan
o(01,x) = max{i(t,x), month(i(t)) = 01}
                 ...
o(12,x) = max{i(t,x), month(i(t)) = 12}
@EndifMan
@IfDoc
@BeginMath
\begin{array}{c}
o(\mbox{01},x) = \mbox{\bf max}\{i(t,x), \mbox{month}(i(t)) = \mbox{01}\} \\
\vdots \\
o(\mbox{12},x) = \mbox{\bf max}\{i(t,x), \mbox{month}(i(t)) = \mbox{12}\} \\
\end{array}
@EndMath
@EndifDoc
@EndDescription

@EndOperator


@BeginOperator_ymonmean

@Title     = Multi-year monthly mean

@BeginDescription
@IfMan
o(01,x) = mean{i(t,x), month(i(t)) = 01}
                 ...
o(12,x) = mean{i(t,x), month(i(t)) = 12}
@EndifMan
@IfDoc
@BeginMath
\begin{array}{c}
o(\mbox{01},x) = \mbox{\bf mean}\{i(t,x), \mbox{month}(i(t)) = \mbox{01}\} \\
\vdots \\
o(\mbox{12},x) = \mbox{\bf mean}\{i(t,x), \mbox{month}(i(t)) = \mbox{12}\} \\
\end{array}
@EndMath
@EndifDoc
@EndDescription

@EndOperator


@BeginOperator_ymonavg

@Title     = Multi-year monthly average

@BeginDescription
@IfMan
o(01,x) = avg{i(t,x), month(i(t)) = 01}
                 ...
o(12,x) = avg{i(t,x), month(i(t)) = 12}
@EndifMan
@IfDoc
@BeginMath
\begin{array}{c}
o(\mbox{01},x) = \mbox{\bf avg}\{i(t,x), \mbox{month}(i(t)) = \mbox{01}\} \\
\vdots \\
o(\mbox{12},x) = \mbox{\bf avg}\{i(t,x), \mbox{month}(i(t)) = \mbox{12}\} \\
\end{array}
@EndMath
@EndifDoc
@EndDescription

@EndOperator


@BeginOperator_ymonstd

@Title     = Multi-year monthly standard deviation

@BeginDescription
@IfMan
o(01,x) = sqrt{var{i(t,x), month(i(t)) = 01}}
                 ...
o(12,x) = sqrt{var{i(t,x), month(i(t)) = 12}}
@EndifMan
@IfDoc
@BeginMath
\begin{array}{c}
o(\mbox{01},x) = \sqrt{\mbox{\bf var}\{i(t,x), \mbox{month}(i(t)) = \mbox{01}\}} \\
\vdots \\
o(\mbox{12},x) = \sqrt{\mbox{\bf var}\{i(t,x), \mbox{month}(i(t)) = \mbox{12}\}} \\
\end{array}
@EndMath
@EndifDoc
@EndDescription

@EndOperator

@EndDoc
*/

#define  NMONTH     17

void *Ymonstat(void *argument)
{
  static char func[] = "Ymonstat";
  int operatorID;
  int operfunc;
  int gridsize;
  int varID;
  int recID;
  int gridID;
  int vdate, vtime;
  int year, month;
  int nrecs, nrecords;
  int levelID;
  int tsID;
  int otsID;
  long nsets[NMONTH];
  int streamID1, streamID2;
  int vlistID1, vlistID2, taxisID1, taxisID2;
  int nmiss;
  int nvars, nlevel;
  int *recVarID, *recLevelID;
  int vdates[NMONTH], vtimes[NMONTH];
  double missval;
  FIELD **vars1[NMONTH], **vars2[NMONTH];
  FIELD field;

  cdoInitialize(argument);

  cdoOperatorAdd("ymonmin",  func_min,  0, NULL);
  cdoOperatorAdd("ymonmax",  func_max,  0, NULL);
  cdoOperatorAdd("ymonmean", func_mean, 0, NULL);
  cdoOperatorAdd("ymonavg",  func_avg,  0, NULL);
  cdoOperatorAdd("ymonstd",  func_std,  0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

  for ( month = 0; month < NMONTH; month++ )
    {
      vars1[month] = NULL;
      vars2[month] = NULL;
      nsets[month] = 0;
    }

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

  tsID = 0;
  otsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      if ( cdoVerbose ) cdoPrint("process timestep: %d %d %d", tsID+1, vdate, vtime);

      year  =  vdate / 10000;
      month = (vdate - year*10000) / 100;
      if ( month < 0 || month >= NMONTH )
	cdoAbort("month %d out of range!", month);

      vdates[month] = vdate;
      vtimes[month] = vtime;

      if ( vars1[month] == NULL )
	{
	  vars1[month] = (FIELD **) malloc(nvars*sizeof(FIELD *));
	  if ( operfunc == func_std )
	    vars2[month] = (FIELD **) malloc(nvars*sizeof(FIELD *));

	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      gridID   = vlistInqVarGrid(vlistID1, varID);
	      gridsize = gridInqSize(gridID);
	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      missval  = vlistInqVarMissval(vlistID1, varID);

	      vars1[month][varID] = (FIELD *)  malloc(nlevel*sizeof(FIELD));
	      if ( operfunc == func_std )
		vars2[month][varID] = (FIELD *)  malloc(nlevel*sizeof(FIELD));
	      
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  vars1[month][varID][levelID].grid    = gridID;
		  vars1[month][varID][levelID].nmiss   = 0;
		  vars1[month][varID][levelID].missval = missval;
		  vars1[month][varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
		  if ( operfunc == func_std )
		    {
		      vars2[month][varID][levelID].grid    = gridID;
		      vars2[month][varID][levelID].nmiss   = 0;
		      vars2[month][varID][levelID].missval = missval;
		      vars2[month][varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
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

	  if ( nsets[month] == 0 )
	    {
	      streamReadRecord(streamID1, vars1[month][varID][levelID].ptr, &nmiss);
	      vars1[month][varID][levelID].nmiss = nmiss;
	    }
	  else
	    {
	      streamReadRecord(streamID1, field.ptr, &field.nmiss);
	      field.grid    = vars1[month][varID][levelID].grid;
	      field.missval = vars1[month][varID][levelID].missval;

	      if ( operfunc == func_std )
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

      if ( nsets[month] == 0 && operfunc == func_std )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	    gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	    nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevel; levelID++ )
	      farmoq(&vars2[month][varID][levelID], vars1[month][varID][levelID]);
	  }

      nsets[month]++;
      tsID++;
    }

  for ( month = 0; month < NMONTH; month++ )
    if ( nsets[month] )
      {
	if ( operfunc == func_mean || operfunc == func_avg )
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		farcmul(&vars1[month][varID][levelID], 1.0/nsets[month]);
	    }
	else if ( operfunc == func_std )
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		farcstd(&vars1[month][varID][levelID], vars2[month][varID][levelID], 1.0/nsets[month]);
	    }

	taxisDefVdate(taxisID2, vdates[month]);
	taxisDefVtime(taxisID2, vtimes[month]);
	streamDefTimestep(streamID2, otsID++);

	for ( recID = 0; recID < nrecords; recID++ )
	  {
	    varID    = recVarID[recID];
	    levelID  = recLevelID[recID];

	    if ( otsID == 1 || vlistInqVarTime(vlistID1, varID) == TIME_VARIABLE )
	      {
		streamDefRecord(streamID2, varID, levelID);
		streamWriteRecord(streamID2, vars1[month][varID][levelID].ptr,
				  vars1[month][varID][levelID].nmiss);
	      }
	  }
      }

  for ( month = 0; month < NMONTH; month++ )
    {
      if ( vars1[month] != NULL )
	{
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  free(vars1[month][varID][levelID].ptr);
		  if ( operfunc == func_std ) free(vars2[month][varID][levelID].ptr);
		}
	      
	      free(vars1[month][varID]);
	      if ( operfunc == func_std ) free(vars2[month][varID]);
	    }
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

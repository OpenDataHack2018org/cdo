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

@Name      = Comp
@Title     = Comparison of two fields
@Section   = Comparisons
@Class     = Comparison
@Arguments = ifile1 ifile2 ofile
@Operators = eq  ne  le  lt  ge  gt

@EndModule


@BeginOperator_eq

@Title     = Equal

@BeginDescription
@IfMan
          /   1   if i_1(t,x) EQ i_2(t,x)  AND  i_1(t,x),i_2(t,x) NE miss
o(t,x) = <    0   if i_1(t,x) NE i_2(t,x)  AND  i_1(t,x),i_2(t,x) NE miss
          \  miss if i_1(t,x) EQ miss      OR   i_2(t,x) EQ miss
@EndifMan
@IfDoc
@BeginMath
o(t,x) = \left\{
\begin{array}{cll}
  1   & \mbox{if} \;\; i_1(t,x) = i_2(t,x)    & \wedge \;\; i_1(t,x), i_2(t,x) \neq \mbox{miss}  \\
  0   & \mbox{if} \;\; i_1(t,x) \neq i_2(t,x) & \wedge \;\; i_1(t,x), i_2(t,x) \neq \mbox{miss}  \\
 \mbox{miss} & \mbox{if} \;\; i_1(t,x) = \mbox{miss} & \vee   \;\; i_2(t,x) = \mbox{miss}        \\
\end{array}   \right.
@EndMath
@EndifDoc
@EndDescription

@EndOperator


@BeginOperator_ne

@Title     = Not equal

@BeginDescription
@IfMan
          /   1   if i_1(t,x) NE i_2(t,x)  AND  i_1(t,x),i_2(t,x) NE miss
o(t,x) = <    0   if i_1(t,x) EQ i_2(t,x)  AND  i_1(t,x),i_2(t,x) NE miss
          \  miss if i_1(t,x) EQ miss      OR   i_2(t,x) EQ miss
@EndifMan
@IfDoc
@BeginMath
o(t,x) = \left\{
\begin{array}{cll}
  1   & \mbox{if} \;\; i_1(t,x) \neq i_2(t,x) & \wedge \;\; i_1(t,x), i_2(t,x) \neq \mbox{miss}  \\
  0   & \mbox{if} \;\; i_1(t,x) = i_2(t,x)    & \wedge \;\; i_1(t,x), i_2(t,x) \neq \mbox{miss}  \\
 \mbox{miss} & \mbox{if} \;\; i_1(t,x) = \mbox{miss} & \vee   \;\; i_2(t,x) = \mbox{miss}        \\
\end{array}   \right.
@EndMath
@EndifDoc
@EndDescription

@EndOperator


@BeginOperator_le

@Title     = Less equal

@BeginDescription
@IfMan
          /   1   if i_1(t,x) LE i_2(t,x)  AND  i_1(t,x),i_2(t,x) NE miss
o(t,x) = <    0   if i_1(t,x) GT i_2(t,x)  AND  i_1(t,x),i_2(t,x) NE miss
          \  miss if i_1(t,x) EQ miss      OR   i_2(t,x) EQ miss
@EndifMan
@IfDoc
@BeginMath
o(t,x) = \left\{
\begin{array}{cll}
  1   & \mbox{if} \;\; i_1(t,x) \leq i_2(t,x) & \wedge \;\; i_1(t,x), i_2(t,x) \neq \mbox{miss}  \\
  0   & \mbox{if} \;\; i_1(t,x) > i_2(t,x)    & \wedge \;\; i_1(t,x), i_2(t,x) \neq \mbox{miss}  \\
 \mbox{miss} & \mbox{if} \;\; i_1(t,x) = \mbox{miss} & \vee   \;\; i_2(t,x) = \mbox{miss}        \\
\end{array}   \right.
@EndMath
@EndifDoc
@EndDescription

@EndOperator


@BeginOperator_lt

@Title     = Less then

@BeginDescription
@IfMan
          /   1   if i_1(t,x) LT i_2(t,x)  AND  i_1(t,x),i_2(t,x) NE miss
o(t,x) = <    0   if i_1(t,x) GE i_2(t,x)  AND  i_1(t,x),i_2(t,x) NE miss
          \  miss if i_1(t,x) EQ miss      OR   i_2(t,x) EQ miss
@EndifMan
@IfDoc
@BeginMath
o(t,x) = \left\{
\begin{array}{cll}
  1   & \mbox{if} \;\; i_1(t,x) < i_2(t,x)    & \wedge \;\; i_1(t,x), i_2(t,x) \neq \mbox{miss}  \\
  0   & \mbox{if} \;\; i_1(t,x) \geq i_2(t,x) & \wedge \;\; i_1(t,x), i_2(t,x) \neq \mbox{miss}  \\
 \mbox{miss} & \mbox{if} \;\; i_1(t,x) = \mbox{miss} & \vee   \;\; i_2(t,x) = \mbox{miss}        \\
\end{array}   \right.
@EndMath
@EndifDoc
@EndDescription

@EndOperator


@BeginOperator_ge

@Title     = Greater equal

@BeginDescription
@IfMan
          /   1   if i_1(t,x) GE i_2(t,x)  AND  i_1(t,x),i_2(t,x) NE miss
o(t,x) = <    0   if i_1(t,x) LT i_2(t,x)  AND  i_1(t,x),i_2(t,x) NE miss
          \  miss if i_1(t,x) EQ miss      OR   i_2(t,x) EQ miss
@EndifMan
@IfDoc
@BeginMath
o(t,x) = \left\{
\begin{array}{cll}
  1   & \mbox{if} \;\; i_1(t,x) \geq i_2(t,x) & \wedge \;\; i_1(t,x), i_2(t,x) \neq \mbox{miss}  \\
  0   & \mbox{if} \;\; i_1(t,x) < i_2(t,x)    & \wedge \;\; i_1(t,x), i_2(t,x) \neq \mbox{miss}  \\
 \mbox{miss} & \mbox{if} \;\; i_1(t,x) = \mbox{miss} & \vee   \;\; i_2(t,x) = \mbox{miss}        \\
\end{array}   \right.
@EndMath
@EndifDoc
@EndDescription

@EndOperator


@BeginOperator_gt

@Title     = Greater then

@BeginDescription
@IfMan
          /   1   if i_1(t,x) GT i_2(t,x)  AND  i_1(t,x),i_2(t,x) NE miss
o(t,x) = <    0   if i_1(t,x) LE i_2(t,x)  AND  i_1(t,x),i_2(t,x) NE miss
          \  miss if i_1(t,x) EQ miss      OR   i_2(t,x) EQ miss
@EndifMan
@IfDoc
@BeginMath
o(t,x) = \left\{
\begin{array}{cll}
  1   & \mbox{if} \;\; i_1(t,x) > i_2(t,x)    & \wedge \;\; i_1(t,x), i_2(t,x) \neq \mbox{miss}  \\
  0   & \mbox{if} \;\; i_1(t,x) \leq i_2(t,x) & \wedge \;\; i_1(t,x), i_2(t,x) \neq \mbox{miss}  \\
 \mbox{miss} & \mbox{if} \;\; i_1(t,x) = \mbox{miss} & \vee   \;\; i_2(t,x) = \mbox{miss}        \\
\end{array}   \right.
@EndMath
@EndifDoc
@EndDescription

@EndOperator

@EndDoc
*/


void *Comp(void *argument)
{
  static char func[] = "Comp";
  int EQ, NE, LE, LT, GE, GT;
  int operatorID;
  int streamID1, streamID2, streamID3;
  int gridsize;
  int nrecs, nrecs2, nvars = 0, nlev, recID;
  int tsID;
  int varID, levelID;
  int offset;
  int vlistID1, vlistID2, vlistID3;
  int nmiss, nmiss2, nmiss3;
  int i;
  int lptype = 0;
  double missval1, missval2 = 0;
  double *array1, *array2, *array3;
  int **varnmiss2 = NULL;
  double **vardata2 = NULL;
  int taxisID1, taxisID3;

  cdoInitialize(argument);

  EQ = cdoOperatorAdd("eq", 0, 0, NULL);
  NE = cdoOperatorAdd("ne", 0, 0, NULL);
  LE = cdoOperatorAdd("le", 0, 0, NULL);
  LT = cdoOperatorAdd("lt", 0, 0, NULL);
  GE = cdoOperatorAdd("ge", 0, 0, NULL);
  GT = cdoOperatorAdd("gt", 0, 0, NULL);

  operatorID = cdoOperatorID();

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));
  streamID2 = streamOpenRead(cdoStreamName(1));
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = streamInqVlist(streamID2);
  vlistID3 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID3);

  if ( vlistNrecs(vlistID2) == 1 && vlistNrecs(vlistID1) != 1 )
    lptype = 2;
  else
    vlistCompare(vlistID1, vlistID2, func_sft);

  nospec(vlistID1);

  streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());
  if ( streamID3 < 0 ) cdiError(streamID3, "Open failed on %s", cdoStreamName(2));

  streamDefVlist(streamID3, vlistID3);

  gridsize = vlistGridsizeMax(vlistID1);

  array1 = (double *) malloc(gridsize*sizeof(double));
  array2 = (double *) malloc(gridsize*sizeof(double));
  array3 = (double *) malloc(gridsize*sizeof(double));

  if ( cdoVerbose )
    cdoPrint("Number of timesteps: file1 %d, file2 %d", vlistNtsteps(vlistID1), vlistNtsteps(vlistID2));

  if ( vlistNtsteps(vlistID2) == 1 && vlistNtsteps(vlistID1) != 1 && lptype == 0 )
    {
      lptype = 1;
      nvars  = vlistNvars(vlistID2);
      vardata2  = (double **) malloc(nvars*sizeof(double *));
      varnmiss2 = (int **) malloc(nvars*sizeof(int *));
      for ( varID = 0; varID < nvars; varID++ )
	{
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  nlev     = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	  vardata2[varID]  = (double *) malloc(nlev*gridsize*sizeof(double));
	  varnmiss2[varID] = (int *) malloc(nlev*sizeof(int));
	}

      cdoPrint("Filling up stream >%s< by copying the first timestep.", cdoStreamName(1));
    }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      if ( tsID == 0 || lptype == 0 )
	{
	  nrecs2 = streamInqTimestep(streamID2, tsID);
	  if ( nrecs2 == 0 )
	    cdoAbort("Input streams have different number of timesteps!");
	}
	  
      taxisCopyTimestep(taxisID3, taxisID1);

      streamDefTimestep(streamID3, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array1, &nmiss);

	  if ( tsID == 0 || lptype == 0 )
	    {
	      if ( recID == 0 || lptype != 2 )
		{
		  streamInqRecord(streamID2, &varID, &levelID);
		  streamReadRecord(streamID2, array2, &nmiss2);
		}

	      if ( lptype == 1 )
		{
		  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
		  offset   = gridsize*levelID;
		  memcpy(vardata2[varID]+offset, array2, gridsize*sizeof(double));
		  varnmiss2[varID][levelID] = nmiss2;
		}
	    }
	  else if ( lptype == 1 )
	    {
	      gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	      offset   = gridsize*levelID;
	      memcpy(array2, vardata2[varID]+offset, gridsize*sizeof(double));
	      nmiss2 = varnmiss2[varID][levelID];
	    }

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  missval1 = vlistInqVarMissval(vlistID1, varID);

	  if ( recID == 0 || lptype != 2 )
	    {
	      /* gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID)); */
	      missval2 = vlistInqVarMissval(vlistID2, varID);
	    }

	  if ( operatorID == EQ )
	    {
	      for ( i = 0; i < gridsize; i++ )
		array3[i] = (DBL_IS_EQUAL(array1[i], missval1) || DBL_IS_EQUAL(array2[i], missval2) ?
			     missval1 : DBL_IS_EQUAL(array1[i], array2[i]));
	    }
	  else if ( operatorID == NE )
	    {
	      for ( i = 0; i < gridsize; i++ )
		array3[i] = (DBL_IS_EQUAL(array1[i], missval1) || DBL_IS_EQUAL(array2[i], missval2) ?
			     missval1 : !DBL_IS_EQUAL(array1[i], array2[i]));
	    }
	  else if ( operatorID == LE )
	    {
	      for ( i = 0; i < gridsize; i++ )
		array3[i] = (DBL_IS_EQUAL(array1[i], missval1) || DBL_IS_EQUAL(array2[i], missval2) ?
			     missval1 : array1[i] <= array2[i]);
	    }
	  else if ( operatorID == LT )
	    {
	      for ( i = 0; i < gridsize; i++ )
		array3[i] = (DBL_IS_EQUAL(array1[i], missval1) || DBL_IS_EQUAL(array2[i], missval2) ?
			     missval1 : array1[i] < array2[i]);
	    }
	  else if ( operatorID == GE )
	    {
	      for ( i = 0; i < gridsize; i++ )
		array3[i] = (DBL_IS_EQUAL(array1[i], missval1) || DBL_IS_EQUAL(array2[i], missval2) ?
			     missval1 : array1[i] >= array2[i]);
	    }
	  else if ( operatorID == GT )
	    {
	      for ( i = 0; i < gridsize; i++ )
		array3[i] = (DBL_IS_EQUAL(array1[i], missval1) || DBL_IS_EQUAL(array2[i], missval2) ?
			     missval1 : array1[i] > array2[i]);
	    }
	  else
	    {
	      cdoAbort("Operator not implemented!");
	    }

	  nmiss3 = 0;
	  for ( i = 0; i < gridsize; i++ )
	    if ( DBL_IS_EQUAL(array3[i], missval1) ) nmiss3++;

	  streamDefRecord(streamID3, varID, levelID);
	  streamWriteRecord(streamID3, array3, nmiss3);
	}

      tsID++;
    }

  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);

  if ( vardata2 )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  free(vardata2[varID]);
	  free(varnmiss2[varID]);
	}

      free(vardata2);
      free(varnmiss2);
    }

  if ( array3 ) free(array3);
  if ( array2 ) free(array2);
  if ( array1 ) free(array1);

  cdoFinish();

  return (0);
}

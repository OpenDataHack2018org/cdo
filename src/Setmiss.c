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

#if  defined  (HAVE_CONFIG_H)
#  include "config.h"
#endif

#if  defined  (HAVE_ISNAN) && ! defined(__cplusplus)
int isnan(const double x);
#endif

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

/*
@BeginDoc

@BeginModule

@Name      = Setmiss
@Title     = Set missing data
@Section   = Missing values
@Arguments = ifile ofile
@Operators = setmissval setctomiss setmisstoc setrtomiss

@EndModule


@BeginOperator_setmissval

@Title     = Set a new missing value
@Parameter = miss

@BeginDescription
@IfMan
         / miss   if i(t,x) EQ miss
o(t,x) = 
         \ i(t,x) if i(t,x) NE miss
@EndifMan
@IfDoc
@BeginMath
o(t,x) = \left\{
\begin{array}{cll}
 \mbox{miss} & \mbox{if} \;\; i(t,x) = \mbox{\sl miss}      \\
 i(t,x)      & \mbox{if} \;\; i(t,x) \neq \mbox{\sl miss}   \\
\end{array}   \right.
@EndMath
@EndifDoc
@EndDescription

@BeginParameter
@Item = miss
FLOAT  New missing value
@EndParameter

@EndOperator


@BeginOperator_setctomiss

@Title     = Set constant to missing value
@Parameter = c

@BeginDescription
@IfMan
         / miss   if i(t,x) EQ c
o(t,x) = 
         \ i(t,x) if i(t,x) NE c
@EndifMan
@IfDoc
@BeginMath
o(t,x) = \left\{
\begin{array}{cll}
 \mbox{miss} & \mbox{if} \;\; i(t,x) = \mbox{\sl c}      \\
 i(t,x)      & \mbox{if} \;\; i(t,x) \neq \mbox{\sl c}   \\
\end{array}   \right.
@EndMath
@EndifDoc
@EndDescription

@BeginParameter
@Item = c
FLOAT  Constant
@EndParameter

@EndOperator


@BeginOperator_setmisstoc

@Title     = Set missing value to constant
@Parameter = c

@BeginDescription
@IfMan
         / c      if i(t,x) EQ miss
o(t,x) = 
         \ i(t,x) if i(t,x) NE miss
@EndifMan
@IfDoc
@BeginMath
o(t,x) = \left\{
\begin{array}{cll}
 \mbox{\sl c} & \mbox{if} \;\; i(t,x) = \mbox{miss}     \\
 i(t,x)       & \mbox{if} \;\; i(t,x) \neq \mbox{miss}  \\
\end{array}   \right.
@EndMath
@EndifDoc
@EndDescription

@BeginParameter
@Item = c
FLOAT  Constant
@EndParameter

@EndOperator


@BeginOperator_setrtomiss

@Title     = Set range to missing value
@Parameter = rmin rmax

@BeginDescription
@IfMan
         / miss   if i(t,x) GE rmin AND i(t,x) LE rmax
o(t,x) = 
         \ i(t,x) if i(t,x) LT rmin AND i(t,x) GT rmax
@EndifMan
@IfDoc
@BeginMath
o(t,x) = \left\{
\begin{array}{cll}
 \mbox{miss}   & \mbox{if} \;\; i(t,x) \geq \mbox{\sl rmin} \wedge i(t,x) \leq \mbox{\sl rmax}  \\
 i(t,x) & \mbox{if} \;\; i(t,x) < \mbox{\sl rmin} \vee i(t,x) > \mbox{\sl rmax}  \\
\end{array}   \right.
@EndMath
@EndifDoc
@EndDescription

@BeginParameter rmax
@Item = rmin
FLOAT  Lower bound
@Item = rmax
FLOAT  Upper bound
@EndParameter

@EndOperator

@EndDoc
*/


void *Setmiss(void *argument)
{
  static char func[] = "Setmiss";
  int SETMISSVAL, SETCTOMISS, SETMISSTOC, SETRTOMISS;
  int operatorID;
  int streamID1, streamID2;
  int gridsize;
  int nrecs, recID;
  int nvars;
  int tsID;
  int varID, levelID;
  int vlistID1, vlistID2;
  int nmiss;
  int i;
  double missval, missval2 = 0;
  double rconst = 0, rmin = 0, rmax = 0;
  double *array;
  int taxisID1, taxisID2;

  cdoInitialize(argument);

  SETMISSVAL = cdoOperatorAdd("setmissval", 0, 0, "missing value");
  SETCTOMISS = cdoOperatorAdd("setctomiss", 0, 0, "constant");
  SETMISSTOC = cdoOperatorAdd("setmisstoc", 0, 0, "constant");
  SETRTOMISS = cdoOperatorAdd("setrtomiss", 0, 0, "range (min, max)");

  operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  if ( operatorID == SETMISSVAL )
    {
      operatorCheckArgc(1);
      missval2 = atof(operatorArgv()[0]);
    }
  else if ( operatorID == SETCTOMISS || operatorID == SETMISSTOC )
    {
      operatorCheckArgc(1);
      if ( operatorArgv()[0][0] == 'n' || operatorArgv()[0][0] == 'N' )
	{
#if ! defined  (HAVE_ISNAN)
	  cdoWarning("function >isnan< not available!");
#endif
	  rconst = 0.0/0.0;
	}
      else
	rconst = atof(operatorArgv()[0]);
    }
  else
    {
      operatorCheckArgc(2);
      rmin = atof(operatorArgv()[0]);
      rmax = atof(operatorArgv()[1]);
    }

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if ( operatorID == SETMISSVAL )
    {
      nvars = vlistNvars(vlistID2);
      for ( varID = 0; varID < nvars; varID++ )
	vlistDefVarMissval(vlistID2, varID, missval2);
    }

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);

  array = (double *) malloc(gridsize*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array, &nmiss);

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  missval = vlistInqVarMissval(vlistID1, varID);

	  if ( operatorID == SETMISSVAL )
	    {
	      nmiss = 0;
	      for ( i = 0; i < gridsize; i++ )
		if ( DBL_IS_EQUAL(array[i], missval)  || DBL_IS_EQUAL(array[i], (float)missval) ||
		     DBL_IS_EQUAL(array[i], missval2) || DBL_IS_EQUAL(array[i], (float)missval2) )
		  {
		    array[i] = missval2;
		    nmiss++;
		  }
	    }
	  else if ( operatorID == SETCTOMISS )
	    {
#if  defined  (HAVE_ISNAN)
	      if ( isnan(rconst) )
		{
		  for ( i = 0; i < gridsize; i++ )
		    if ( isnan(array[i]) )
		      {
			array[i] = missval;
			nmiss++;
		      }
		}
	      else
#endif
		{
		  for ( i = 0; i < gridsize; i++ )
		    if ( DBL_IS_EQUAL(array[i], rconst) || DBL_IS_EQUAL(array[i], (float)rconst) )
		      {
			array[i] = missval;
			nmiss++;
		      }
		}
	    }
	  else if ( operatorID == SETMISSTOC )
	    {
	      nmiss = 0;
	      for ( i = 0; i < gridsize; i++ )
		if ( DBL_IS_EQUAL(array[i], missval) || DBL_IS_EQUAL(array[i], (float)missval) )
		  {
		    array[i] = rconst;
		  }
	    }
	  else /* setrtomiss */
	    {
	      for ( i = 0; i < gridsize; i++ )
		if ( array[i] >= rmin && array[i] <= rmax )
		  {
		    array[i] = missval;
		    nmiss++;
		  }
	    }

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, array, nmiss);
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( array ) free(array);

  cdoFinish();

  return (0);
}

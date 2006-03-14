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

/*
@BeginDoc

@BeginModule

@Name      = Math
@Title     = Mathematical functions
@Section   = Mathematical functions
@Class     = Arithmetic
@Arguments = ifile ofile
@Operators = abs sqr sqrt exp log log10 sin cos tan asin acos atan

@EndModule


@BeginOperator_abs

@Title     = Absolute value

@BeginDescription
@IfDoc
@BeginMath
o(t,x) = abs{i(t,x)}
@EndMath
@EndifDoc
Calculates the absolute value of i(t,x).
@EndDescription

@EndOperator


@BeginOperator_sqr

@Title     = Square

@BeginDescription
@IfDoc
@BeginMath
o(t,x) = {i(t,x)}^2
@EndMath
@EndifDoc
Calculates the value of i(t,x) raised to the power of 2.
@EndDescription

@EndOperator


@BeginOperator_sqrt

@Title     = Square root

@BeginDescription
@IfDoc
@BeginMath
o(t,x) = \sqrt{i(t,x)}
@EndMath
@EndifDoc
Calculates the non-negative square root of i(t,x).
@EndDescription

@EndOperator


@BeginOperator_exp

@Title     = Exp

@BeginDescription
@IfDoc
@BeginMath
o(t,x) = e^{i(t,x)}
@EndMath
@EndifDoc
Calculates e (the base of natural logarithms) raised to the power of i(t,x).
@EndDescription

@EndOperator


@BeginOperator_log

@Title     = Logarithm

@BeginDescription
@IfDoc
@BeginMath
o(t,x) = log(i(t,x))
@EndMath
@EndifDoc
Calculates the natural logarithm of i(t,x).
@EndDescription

@EndOperator


@BeginOperator_log10

@Title     = Logarithm base 10

@BeginDescription
@IfDoc
@BeginMath
o(t,x) = log_{10}(i(t,x))
@EndMath
@EndifDoc
Calculates the base-10 logarithm of i(t,x).
@EndDescription

@EndOperator


@BeginOperator_sin

@Title     = Sine

@BeginDescription
@IfDoc
@BeginMath
o(t,x) = \sin(i(t,x))
@EndMath
@EndifDoc
Calculates the sine of i(t,x), where i(t,x) is given in radians.
@EndDescription

@EndOperator


@BeginOperator_cos

@Title     = Cosine

@BeginDescription
@IfDoc
@BeginMath
o(t,x) = \cos(i(t,x))
@EndMath
@EndifDoc
Calculates the cosine of i(t,x), where i(t,x) is given in radians.
@EndDescription

@EndOperator


@BeginOperator_tan

@Title     = Tangent

@BeginDescription
@IfDoc
@BeginMath
o(t,x) = \tan(i(t,x))
@EndMath
@EndifDoc
Calculates the tangent of i(t,x), where i(t,x) is given in radians.
@EndDescription

@EndOperator


@BeginOperator_asin

@Title     = Arcus sine

@BeginDescription
@IfDoc
@BeginMath
o(t,x) = asin(i(t,x))
@EndMath
@EndifDoc
Calculates the arcus sine of i(t,x); that is the value whose sine is i(t,x).
@EndDescription

@EndOperator


@BeginOperator_acos

@Title     = Arcus cosine

@BeginDescription
@IfDoc
@BeginMath
o(t,x) = acos(i(t,x))
@EndMath
@EndifDoc
Calculates the arcus cosine of i(t,x); that is the value whose cosine is i(t,x).
@EndDescription

@EndOperator


@BeginOperator_atan

@Title     = Arcus tangent

@BeginDescription
@IfDoc
@BeginMath
o(t,x) = atan(i(t,x))
@EndMath
@EndifDoc
Calculates the arcus tangent of i(t,x); that is the value whose tangent is i(t,x).
@EndDescription

@EndOperator

@EndDoc
*/


void *Math(void *argument)
{
  static char func[] = "Math";
  enum {ABS, SQR, SQRT, EXP, LOG, LOG10, SIN, COS, TAN, ASIN, ACOS, ATAN};
  int operatorID;
  int operfunc;
  int streamID1, streamID2;
  int gridsize;
  int nrecs, recID;
  int tsID;
  int varID, levelID;
  int vlistID1, vlistID2;
  int nmiss, nmiss2;
  int i;
  double missval1, missval2;
  double *array1, *array2;
  int taxisID1, taxisID2;

  cdoInitialize(argument);

  cdoOperatorAdd("abs",   ABS,   0, NULL);
  cdoOperatorAdd("sqr",   SQR,   0, NULL);
  cdoOperatorAdd("sqrt",  SQRT,  0, NULL);
  cdoOperatorAdd("exp",   EXP,   0, NULL);
  cdoOperatorAdd("log",   LOG,   0, NULL);
  cdoOperatorAdd("log10", LOG10, 0, NULL);
  cdoOperatorAdd("sin",   SIN,   0, NULL);
  cdoOperatorAdd("cos",   COS,   0, NULL);
  cdoOperatorAdd("tan",   TAN,   0, NULL);
  cdoOperatorAdd("asin",  ASIN,  0, NULL);
  cdoOperatorAdd("acos",  ACOS,  0, NULL);
  cdoOperatorAdd("atan",  ATAN,  0, NULL);
 
  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  nospec(vlistID1);

  gridsize = vlistGridsizeMax(vlistID1);

  array1 = (double *) malloc(gridsize*sizeof(double));
  array2 = (double *) malloc(gridsize*sizeof(double));

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array1, &nmiss);

	  missval1 = vlistInqVarMissval(vlistID1, varID);
	  missval2 = missval1;
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

	  switch ( operfunc )
	    {
	    case ABS:
	      for ( i = 0; i < gridsize; i++ )
		array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : fabs(array1[i]);
	      break;
	    case SQR:
	      for ( i = 0; i < gridsize; i++ )
		array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : MUL(array1[i], array1[i]);
	      break;
	    case SQRT:
	      for ( i = 0; i < gridsize; i++ )
		array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : ROOT(array1[i]);
	      break;
	    case EXP:
	      for ( i = 0; i < gridsize; i++ )
		array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : exp(array1[i]);
	      break;
	    case LOG:
	      for ( i = 0; i < gridsize; i++ )
		array2[i] = DBL_IS_EQUAL(array1[i], missval1) || array1[i] < 0 ? missval1 : log(array1[i]);
	      break;
	    case LOG10:
	      for ( i = 0; i < gridsize; i++ )
		array2[i] = DBL_IS_EQUAL(array1[i], missval1) || array1[i] < 0 ? missval1 : log10(array1[i]);
	      break;
	    case SIN:
	      for ( i = 0; i < gridsize; i++ )
		array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : sin(array1[i]);
	      break;
	    case COS:
	      for ( i = 0; i < gridsize; i++ )
		array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : cos(array1[i]);
	      break;
	    case TAN:
	      for ( i = 0; i < gridsize; i++ )
		array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : tan(array1[i]);
	      break;
	    case ASIN:
	      for ( i = 0; i < gridsize; i++ )
		array2[i] = DBL_IS_EQUAL(array1[i], missval1) || array1[i] < -1
		         || array1[i] > 1 ? missval1 : asin(array1[i]);
	      break;
	    case ACOS:
	      for ( i = 0; i < gridsize; i++ )
		array2[i] = DBL_IS_EQUAL(array1[i], missval1) || array1[i] < -1
		         || array1[i] > 1 ? missval1 : acos(array1[i]);
	      break;
	    case ATAN:
	      for ( i = 0; i < gridsize; i++ )
		array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : atan(array1[i]);
	      break;
	    default:
	      cdoAbort("operator not implemented!");
	    }

	  nmiss2 = 0;
	  for ( i = 0; i < gridsize; i++ )
	    if ( DBL_IS_EQUAL(array2[i], missval1) ) nmiss2++;

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, array2, nmiss2);
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( array2 ) free(array2);
  if ( array1 ) free(array1);

  cdoFinish();

  return (0);
}

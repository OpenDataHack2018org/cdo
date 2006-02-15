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

@Name      = Arithc
@Title     = Arithmetic with a constant
@Section   = Arithmetic
@Class     = Arithmetic
@Arguments = ifile ofile
@Operators = addc subc mulc divc

@EndModule


@BeginOperator_addc

@Title     = Add by constant
@Parameter = c

@BeginDesciption
@BeginMath
o(t,x) = i(t,x) + c
@EndMath
@EndDesciption

@BeginParameter
@Item = c
FLOAT  Constant
@EndParameter

@EndOperator


@BeginOperator_subc

@Title     = Subtract by constant
@Parameter = c

@BeginDesciption
@BeginMath
o(t,x) = i(t,x) - c
@EndMath
@EndDesciption

@BeginParameter
@Item = c
FLOAT  Constant
@EndParameter

@EndOperator


@BeginOperator_mulc

@Title     = Multiply by constant
@Parameter = c

@BeginDesciption
@BeginMath
o(t,x) = i(t,x) * c
@EndMath
@EndDesciption

@BeginParameter
@Item = c
FLOAT  Constant
@EndParameter

@EndOperator


@BeginOperator_divc

@Title     = Divide by constant
@Parameter = c

@BeginDesciption
@BeginMath
o(t,x) = i(t,x) / c
@EndMath
@EndDesciption

@BeginParameter
@Item = c
FLOAT  Constant
@EndParameter

@EndOperator

@EndDoc
*/


void *Arithc(void *argument)
{
  static char func[] = "Arithc";
  int operatorID;
  int operfunc;
  int streamID1, streamID2;
  int gridsize;
  int nrecs, recID;
  int tsID;
  int varID, levelID;
  int vlistID1, vlistID2;
  double rconst;
  FIELD field;
  int taxisID1, taxisID2;

  cdoInitialize(argument);

  cdoOperatorAdd("addc", func_add, 0, NULL);
  cdoOperatorAdd("subc", func_sub, 0, NULL);
  cdoOperatorAdd("mulc", func_mul, 0, NULL);
  cdoOperatorAdd("divc", func_div, 0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

  operatorInputArg("constant value");
  rconst = atof(operatorArgv()[0]);

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if ( operfunc == func_mul || operfunc == func_div )
    nospec(vlistID1);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);

  field.ptr    = (double *) malloc(gridsize*sizeof(double));
  field.weight = NULL;

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, field.ptr, &field.nmiss);

	  field.grid    = vlistInqVarGrid(vlistID1, varID);
	  field.missval = vlistInqVarMissval(vlistID1, varID);

	  farcfun(&field, rconst, operfunc);

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, field.ptr, field.nmiss);
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( field.ptr ) free(field.ptr);

  cdoFinish();

  return (0);
}

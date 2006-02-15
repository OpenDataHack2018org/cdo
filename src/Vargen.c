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

#include <string.h>
#include <math.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

/*
@BeginDoc

@BeginModule

@Name      = Vargen
@Title     = 
@Section   = Generation of variables
@Arguments = ofile
@Operators = const random

@EndModule


@BeginOperator_const

@Title     = Constant variable
@Parameter = const grid

@BeginDesciption
Generates a constant variable.
@EndDesciption

@BeginParameter const
@Item = const
FLOAT   Constant
@Item = grid
STRING  Grid description file or name
@EndParameter

@EndOperator


@BeginOperator_random

@Title     = Variable with random values
@Parameter = grid

@BeginDesciption
Generates a variable with rectangularly distrubuted random numbers in the interval [0,1].
@EndDesciption

@BeginParameter
@Item = grid
STRING  Grid description file or name
@EndParameter

@EndOperator

@EndDoc
*/


void *Vargen(void *argument)
{
  static char func[] = "Vargen";
  int RANDOM, CONST;
  int operatorID;
  int streamID;
  int nrecs;
  int tsID, recID, varID, levelID;
  int gridsize, i;
  int vlistID;
  int gridID, zaxisID, taxisID;
  const char *gridfile;
  double rconst = 0.0;
  double *array;

  cdoInitialize(argument);

  RANDOM = cdoOperatorAdd("random", 0, 0, "grid description file or name");
  CONST  = cdoOperatorAdd("const",  0, 0, "constant value, grid description file or name");

  operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  if ( operatorID == RANDOM )
    {
      operatorCheckArgc(1);
      gridfile = operatorArgv()[0];
    }
  else
    {
      operatorCheckArgc(2);
      rconst = atof(operatorArgv()[0]);
      gridfile = operatorArgv()[1];
    }

  gridID  = cdoDefineGrid(gridfile);

  zaxisID = zaxisNew(ZAXIS_SURFACE, 1);

  vlistID = vlistNew();
  varID   = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARIABLE);

  taxisID = taxisNew(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID, taxisID);

  streamID = streamOpenWrite(cdoStreamName(0), cdoFiletype());
  if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", cdoStreamName(0));

  streamDefVlist(streamID, vlistID);

  gridsize = gridInqSize(gridID);
  array = (double *) malloc(gridsize*sizeof(double));

  tsID = 0;
  streamDefTimestep(streamID, tsID);

  nrecs = 1;
  for ( recID = 0; recID < nrecs; recID++ )
    {
      levelID = 0;
      streamDefRecord(streamID, varID, levelID);

      if ( operatorID == RANDOM )
	{
	  for ( i = 0; i < gridsize; i++ )
	    array[i] = rand()/(RAND_MAX+1.0);
	}
      else
	{
	  for ( i = 0; i < gridsize; i++ )
	    array[i] = rconst;
	}

      streamWriteRecord(streamID, array, 0);
    }

  streamClose(streamID);

  if ( array ) free(array);

  cdoFinish();

  return (0);
}

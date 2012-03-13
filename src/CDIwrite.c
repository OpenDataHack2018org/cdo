/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2012 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *CDIwrite(void *argument)
{
  int nvars = 1, nlevs = 0, ntimesteps = 0;
  int operatorID;
  int streamID;
  int tsID, varID, levelID;
  int gridsize, i;
  int rval, rstart, rinc;
  int vlistID;
  int gridID = -1, zaxisID, taxisID;
  int vdate, vtime, julday;
  unsigned int seed = 1;
  const char *gridfile;
  double *levels = NULL;
  double ***vars = NULL;

  cdoInitialize(argument);

  operatorInputArg("grid, <nlevs, <ntimesteps, <nvars>>>");
  if ( operatorArgc() < 1 ) cdoAbort("Too few arguments!");
  if ( operatorArgc() > 4 ) cdoAbort("Too many arguments!");

  gridfile = operatorArgv()[0];
  gridID   = cdoDefineGrid(gridfile);
  if ( operatorArgc() >= 2 ) nlevs = atol(operatorArgv()[1]);
  if ( operatorArgc() >= 3 ) ntimesteps = atol(operatorArgv()[2]);
  if ( operatorArgc() >= 4 ) nvars = atol(operatorArgv()[3]);

  srand(seed);

  gridsize = gridInqSize(gridID);

  zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);

  if ( cdoVerbose )
    {
      cdoPrint("gridsize   : %d", gridInqSize);
      cdoPrint("nlevs      : %d", nlevs);
      cdoPrint("ntimesteps : %d", ntimesteps);
      cdoPrint("nvars      : %d", nvars);
    } 

  if ( nlevs <= 0 ) nlevs = 1;
  if ( ntimesteps <= 0 ) ntimesteps = 1;
  if ( nvars <= 0 ) nvars = 1;

  vars = (double ***) malloc(nvars*sizeof(double **));
  for ( varID = 0; varID < nvars; varID++ )
    {
      vars[varID] = (double **) malloc(nlevs*sizeof(double *));
      for ( levelID = 0; levelID < nlevs; levelID++ )
	{
	  vars[varID][levelID] = (double *) malloc(gridsize*sizeof(double));
	  for ( i = 0; i < gridsize; ++i )
	    vars[varID][levelID][i] = varID + rand()/(RAND_MAX+1.0);
	}
    }

  vlistID = vlistCreate();

  for ( i = 0; i < nvars; ++i )
    {
      varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARIABLE);
      vlistDefVarCode(vlistID, varID, varID+1);
      //    vlistDefVarName(vlistID, varID, );
    }

  taxisID = taxisCreate(TAXIS_RELATIVE);
  vlistDefTaxis(vlistID, taxisID);

  // vlistDefNtsteps(vlistID, 1);

  streamID = streamOpenWrite(cdoStreamName(0), cdoFiletype());

  streamDefVlist(streamID, vlistID);

  julday = date_to_julday(CALENDAR_PROLEPTIC, 10101);

  for ( tsID = 0; tsID < ntimesteps; tsID++ )
    {
      rval  = rstart + rinc*tsID;
      vdate = julday_to_date(CALENDAR_PROLEPTIC, julday + tsID);
      vtime = 0;
      taxisDefVdate(taxisID, vdate);
      taxisDefVtime(taxisID, vtime);
      streamDefTimestep(streamID, tsID);

      for ( varID = 0; varID < nvars; varID++ )
        {
          for ( levelID = 0; levelID < nlevs; levelID++ )
            {
              streamDefRecord(streamID, varID, levelID);
              streamWriteRecord(streamID, vars[varID][levelID], 0);
            }
        }
    }

  streamClose(streamID);

  vlistDestroy(vlistID);

  for ( varID = 0; varID < nvars; varID++ )
    {
      free(vars[varID]);
      for ( levelID = 0; levelID < nlevs; levelID++ ) free(vars[varID][levelID]);
    }
  free(vars);

  cdoFinish();

  return (0);
}

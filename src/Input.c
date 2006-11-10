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

      Input     input          ASCII input
      Input     inputsrv       SERVICE input
      Input     inputext       EXTRA input
*/


#include <stdio.h>
#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


static int input_iarray(int nval, int *array)
{
  int i, n;
  int ival = 0;

  for ( i = 0; i < nval; i++ )
    {
      n = scanf("%d", &array[i]);
      if ( n != 1 ) break;

      ival++;
    }

  return (ival);
}


static int input_darray(int nval, double *array)
{
  int i, n;
  int ival = 0;

  for ( i = 0; i < nval; i++ )
    {
      n = scanf("%lg", &array[i]);
      if ( n != 1 ) break;

      ival++;
    }

  return (ival);
}


void *Input(void *argument)
{
  static char func[] = "Input";
  int INPUT, INPUTSRV, INPUTEXT;
  int operatorID;
  int varID = 0;
  int nlevs = 1;
  int gridsize0 = 0, gridsize = 0;
  int gridID = -1, zaxisID, vdate = 0, vtime = 0;
  int nrecs;
  int levelID;
  int tsID, taxisID = 0;
  int streamID = -1;
  int vlistID = -1;
  int nmiss = 0;
  int code = 0, code0 = 0, level = 0, date = 0, time = 0, nlon = 0, nlat = 0;
  int output_filetype = FILETYPE_GRB;
  int rval;
  int ihead[8];
  double levels[1];
  double *array = NULL;

  cdoInitialize(argument);

  INPUT    = cdoOperatorAdd("input",    0, 0, NULL);
  INPUTSRV = cdoOperatorAdd("inputsrv", 0, 0, NULL);
  INPUTEXT = cdoOperatorAdd("inputext", 0, 0, NULL);

  operatorID = cdoOperatorID();

  if ( operatorID == INPUT )
    {
      operatorInputArg("grid description file or name");
      operatorCheckArgc(1);
      gridID = cdoDefineGrid(operatorArgv()[0]);
    }

  nrecs = 0;

  tsID = 0;
  while ( TRUE )
    {
      if ( operatorID == INPUT )
	{
	  output_filetype = cdoFiletype();

	  code     = -1;
	  level    = 0;
	  gridsize = gridInqSize(gridID);
	  date     = 0;
	  time     = 0;
	  
	  if ( nrecs == 0 )
	    array = (double *) malloc(gridsize*sizeof(double));
	  
	  cdoPrint("Enter all %d elements of record %d!\n", gridsize, nrecs+1);
	  
	  rval = input_darray(gridsize, array);

	  if ( nrecs > 0 && rval == 0 ) break;

	  if ( rval != gridsize )
	    cdoAbort("To few input elements (%d of %d)!", rval, gridsize);

	  if ( feof(stdin) ) break;
	}
      else if ( operatorID == INPUTEXT )
	{
	  output_filetype = FILETYPE_EXT;

	  cdoPrint("Enter header (code,level,date,time,nlon,nlat,dispo1,dispo2)"
		   " of record %d (or EOF(=^D))!\n", nrecs+1);

	  rval = input_iarray(4, ihead);
	  if ( feof(stdin) && nrecs == 0 )
	    cdoAbort("To few header elements (%d of %d)!", rval, 4);
	  if ( feof(stdin) ) break;
	  if ( rval != 4 ) cdoAbort("Invalid header input!");

	  date     = ihead[0];
	  code     = ihead[1];
	  level    = ihead[2];
	  gridsize = ihead[3];

	  time = 0;
	  
	  if ( nrecs == 0 )
	    {
	      code0 = code;
	      levels[0] = level;
	      gridsize0 = gridsize;

	      if ( gridsize < 0 )
		cdoAbort("Gridsize must not be negative!\n", gridsize);

	      array = (double *) malloc(gridsize*sizeof(double));

	      gridID = gridCreate(GRID_GENERIC, gridsize);
	    }
	  else
	    {
	      if ( gridsize != gridsize0 )
		cdoAbort("Gridsize must not change!");
	    }
	  
	  cdoPrint("Enter all %d elements of record %d!\n", gridsize, nrecs+1);
	  
	  rval = input_darray(gridsize, array);
	  if ( rval != gridsize ) cdoAbort("Invalid data input!");
	}
      else if ( operatorID == INPUTSRV )
	{
	  output_filetype = FILETYPE_SRV;
	  
	  cdoPrint("Enter header (code,level,date,time,nlon,nlat,dispo1,dispo2)"
		   " of record %d (or EOF(=^D))!\n", nrecs+1);

	  rval = input_iarray(8, ihead);
	  if ( feof(stdin) && nrecs == 0 )
	    cdoAbort("To few header elements (%d of %d)!", rval, 8);
	  if ( feof(stdin) ) break;
	  if ( rval != 8 ) cdoAbort("Invalid header input!");

	  code  = ihead[0];
	  level = ihead[1];
	  date  = ihead[2];
	  time  = ihead[3];
	  nlon  = ihead[4];
	  nlat  = ihead[5];

	  gridsize = nlon*nlat;

	  if ( nrecs == 0 )
	    {
	      code0 = code;
	      levels[0] = level;
	      gridsize0 = gridsize;
	  
	      if ( gridsize < 0 )
		cdoAbort("Gridsize must not be negative!\n", gridsize);

	      array = (double *) malloc(gridsize*sizeof(double));

	      gridID = gridCreate(GRID_GENERIC, gridsize);
	      gridDefXsize(gridID, nlon);
	      gridDefYsize(gridID, nlat);
	    }
	  else
	    {
	      if ( gridsize != gridsize0 )
		cdoAbort("Gridsize must not change!");
	    }
	  
	  cdoPrint("Enter all %d elements of record %d!\n", gridsize, nrecs+1);
	  
	  rval = input_darray(gridsize, array);
	  if ( rval != gridsize ) cdoAbort("Invalid data input!");
	}

      if ( nrecs == 0 )
	{
	  zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);

	  vlistID = vlistCreate();
	  varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARIABLE);
	  vlistDefVarCode(vlistID, varID, code);

	  taxisID = taxisCreate(TAXIS_ABSOLUTE);
	  vlistDefTaxis(vlistID, taxisID);

	  streamID = streamOpenWrite(cdoStreamName(0), output_filetype);
	  if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", cdoStreamName(0));

	  streamDefVlist(streamID, vlistID);
	}

      vdate = date;
      vtime = time;
      taxisDefVdate(taxisID, vdate);
      taxisDefVtime(taxisID, vtime);
      streamDefTimestep(streamID, tsID);

      for ( levelID = 0; levelID < nlevs; levelID++ )
	{
	  streamDefRecord(streamID, varID, levelID);
	  streamWriteRecord(streamID, array, nmiss);
	}

      nrecs++;
      tsID++;
    }

  if ( streamID >= 0 )
    {
      streamClose(streamID);
      vlistDestroy(vlistID);
    }

  if ( array ) free(array);

  cdoFinish();

  return (0);
}

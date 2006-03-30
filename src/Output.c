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

      Output     output          ASCII output
      Output     outputf         Formatted output
      Output     outputint       Integer output
      Output     outputsrv       SERVICE output
      Output     outputext       EXTRA output
*/


#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Output(void *argument)
{
  static char func[] = "Output";
  int OUTPUT, OUTPUTINT, OUTPUTSRV, OUTPUTEXT, OUTPUTF, OUTPUTTS, OUTPUTFLD;
  int operatorID;
  int i;
  int indf;
  int varID, recID;
  int gridsize = 0;
  int gridID, zaxisID, code, vdate, vtime;
  int gridtype, ngrids, gridID2;
  int nrecs;
  int levelID;
  int tsID, taxisID;
  int streamID = 0;
  int vlistID;
  int nmiss, nout;
  int nlon, nlat;
  int hour, minute;
  int year, month, day;
  int nelem = 0;
  int index;
  int ndiffgrids;
  const char *format = NULL;
  double level;
  double *xvals = NULL, *yvals = NULL;
  double *array = NULL;
  double xdate;
  double missval;

  cdoInitialize(argument);

  OUTPUT    = cdoOperatorAdd("output",    0, 0, NULL);
  OUTPUTINT = cdoOperatorAdd("outputint", 0, 0, NULL);
  OUTPUTSRV = cdoOperatorAdd("outputsrv", 0, 0, NULL);
  OUTPUTEXT = cdoOperatorAdd("outputext", 0, 0, NULL);
  OUTPUTF   = cdoOperatorAdd("outputf",   0, 0, NULL);
  OUTPUTTS  = cdoOperatorAdd("outputts",  0, 0, NULL);
  OUTPUTFLD = cdoOperatorAdd("outputfld", 0, 0, NULL);

  operatorID = cdoOperatorID();

  if ( operatorID == OUTPUTF )
    {
      operatorInputArg("format and number of elements");
      operatorCheckArgc(2);
      format = operatorArgv()[0];
      nelem  = atoi(operatorArgv()[1]);
    }

  for ( indf = 0; indf < cdoStreamCnt(); indf++ )
    {
      streamID = streamOpenRead(cdoStreamName(indf));
      if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", cdoStreamName(indf));

      vlistID = streamInqVlist(streamID);

      ngrids = vlistNgrids(vlistID);
      ndiffgrids = 0;
      for ( index = 1; index < ngrids; index++ )
	if ( vlistGrid(vlistID, 0) != vlistGrid(vlistID, index) )
	  ndiffgrids++;

      if ( ndiffgrids > 0 ) cdoAbort("Too many different grids!");

      gridID = vlistGrid(vlistID, 0);
      gridsize = gridInqSize(gridID);

      array = (double *) malloc(gridsize*sizeof(double));

      if ( operatorID == OUTPUTFLD )
	{
	  gridtype = gridInqType(gridID);

	  if ( gridtype == GRID_CURVILINEAR || gridtype == GRID_CELL ||
               gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN )
	    {
	      if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN )
		gridID2 = gridToCell(gridID);
	      else
		gridID2 = gridID;

	      xvals = (double *) malloc(gridsize*sizeof(double));
	      yvals = (double *) malloc(gridsize*sizeof(double));
	      gridInqXvals(gridID2, xvals);
	      gridInqYvals(gridID2, yvals);
	    }	  
	  else cdoAbort("Input grid unsupported!");
	}

      tsID = 0;
      taxisID = vlistInqTaxis(vlistID);
      while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	{
	  vdate = taxisInqVdate(taxisID);
	  vtime = taxisInqVtime(taxisID);

	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      streamInqRecord(streamID, &varID, &levelID);

	      code     = vlistInqVarCode(vlistID, varID);
	      gridID   = vlistInqVarGrid(vlistID, varID);
	      zaxisID  = vlistInqVarZaxis(vlistID, varID);
	      missval  = vlistInqVarMissval(vlistID, varID);
	      gridsize = gridInqSize(gridID);
	      nlon     = gridInqXsize(gridID);
	      nlat     = gridInqYsize(gridID);
	      level    = zaxisInqLevel(zaxisID, levelID);

	      streamReadRecord(streamID, array, &nmiss);

	      if ( operatorID == OUTPUTSRV )
		fprintf(stdout, "%4d %8g %8d %4d %8d %8d %d %d\n", code, level, vdate, vtime, nlon, nlat, 0, 0);

	      if ( operatorID == OUTPUTEXT )
		fprintf(stdout, "%8d %4d %8g %8d\n", vdate, code, level, gridsize);
		
	      if ( operatorID == OUTPUTINT )
		{
		  nout = 0;
		  for ( i = 0; i < gridsize; i++ )
		    {
		      if ( nout == 8 )
			{
			  nout = 0;
			  fprintf(stdout, "\n");
			}
		      fprintf(stdout, " %8d", (int) array[i]);
		      nout++;
		    }
		  fprintf(stdout, "\n");
		}
	      else if ( operatorID == OUTPUTF )
		{
		  nout = 0;
		  for ( i = 0; i < gridsize; i++ )
		    {
		      if ( nout == nelem )
			{
			  nout = 0;
			  fprintf(stdout, "\n");
			}
		      fprintf(stdout, format, array[i]);
		      nout++;
		    }
		  fprintf(stdout, "\n");
		}
	      else if ( operatorID == OUTPUTTS )
		{
		  if ( gridsize > 1 )
		    cdoAbort("operator works only with one gridpoint!");

		  decode_date(vdate, &year, &month, &day);
		  decode_time(vtime, &hour, &minute);

		  /*
		  xdate  = vdate - (vdate/100)*100 + (hour*60 + minute)/1440.;
		  */
		  fprintf(stdout, "%4.4d-%2.2d-%2.2d %2.2d:%2.2d %12.12g\n",
			  year, month, day, hour, minute, array[0]);
		}
	      else if ( operatorID == OUTPUTFLD )
		{
		  hour   = vtime / 100;
		  minute = vtime - hour*100;
		  xdate  = vdate - (vdate/100)*100 + (hour*60 + minute)/1440.;
		  for ( i = 0; i < gridsize; i++ )
		    if ( !DBL_IS_EQUAL(array[i], missval) )
		      fprintf(stdout, "%g\t%g\t%g\t%g\n", xdate, yvals[i], xvals[i], array[i]);
		}
	      else
		{
		  nout = 0;
		  for ( i = 0; i < gridsize; i++ )
		    {
		      if ( nout == 6 )
			{
			  nout = 0;
			  fprintf(stdout, "\n");
			}
		      fprintf(stdout, " %12.6g", array[i]);
		      nout++;
		    }
		  fprintf(stdout, "\n");
		}
	    }
	  tsID++;
	}
      streamClose(streamID);

      if ( array ) free(array);
      if ( xvals ) free(xvals);
      if ( yvals ) free(yvals);
    }

  cdoFinish();

  return (0);
}

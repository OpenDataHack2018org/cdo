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
#include <ctype.h>
#include <math.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "interpol.h"


int readnextpos(FILE *fp, double *julval, double *xpos, double *ypos)
{
  int year, month, day, hour, minute;
  int date, time;
  int stat;

  stat = fscanf(fp, "%d-%d-%d %d:%d %lf %lf", &year, &month, &day, &hour, &minute, xpos, ypos);

  if ( stat != EOF )
    {
      date = year*10000 + month*100 + day;
      time = hour*100 + minute;
      *julval = encode_julval(0, date, time);
    }

  return (stat);
}


void *Intgridtraj(void *argument)
{
  static char func[] = "Intgridtraj";
  int streamID1, streamID2;
  int nrecs, nvars, nlevel;
  int index, ngrids;
  int i, nrecords;
  int tsID, tsIDo, recID, varID, levelID;
  int gridsize;
  int vlistID1, vlistID2;
  int gridID1, gridID2;
  int taxisID1, taxisID2;
  int vdate, vtime;
  int offset;
  int nmiss;
  int *recVarID, *recLevelID;
  double julval1, julval2, julval;
  double fac1, fac2;
  double point;
  double *array, *single1, *single2;
  double **vardata1, **vardata2, *vardatap;
  double xpos, ypos;
  char *posfile;
  FILE *fp;

  cdoInitialize(argument);

  operatorInputArg("filename with grid trajectories");
  operatorCheckArgc(1);

  posfile = operatorArgv()[0];
  fp = fopen(posfile, "r");
  if ( fp == NULL ) cdoAbort("Open failed on %s!", posfile);

  readnextpos(fp, &julval, &xpos, &ypos);

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);

  nvars    = vlistNvars(vlistID1);
  nrecords = vlistNrecs(vlistID1);

  recVarID   = (int *) malloc(nrecords*sizeof(int));
  recLevelID = (int *) malloc(nrecords*sizeof(int));

  gridsize = vlistGridsizeMax(vlistID1);
  array = (double *) malloc(gridsize*sizeof(double));

  vardata1 = (double **) malloc(nvars*sizeof(double*));
  vardata2 = (double **) malloc(nvars*sizeof(double*));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      vardata1[varID] = (double *) malloc(gridsize*nlevel*sizeof(double));
      vardata2[varID] = (double *) malloc(gridsize*nlevel*sizeof(double));
    }

  gridID2 = gridNew(GRID_TRAJECTORY, 1);
  gridDefXsize(gridID2, 1);
  gridDefYsize(gridID2, 1);
  gridDefXvals(gridID2, &xpos);
  gridDefYvals(gridID2, &ypos);

  vlistID2 = vlistDuplicate(vlistID1);

  ngrids = vlistNgrids(vlistID1);
  for ( index = 0; index < ngrids; index++ )
    {
      gridID1 = vlistGrid(vlistID1, index);

      if ( gridInqType(gridID1) != GRID_LONLAT &&
	   gridInqType(gridID1) != GRID_GAUSSIAN )
	cdoAbort("Interpolation of %s data failed!", gridNamePtr(gridInqType(gridID1)) );

      vlistChangeGridIndex(vlistID2, index, gridID2);
    }

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  tsID = 0;
  nrecs = streamInqTimestep(streamID1, tsID++);
  julval1 = encode_julval(0, taxisInqVdate(taxisID1), taxisInqVtime(taxisID1));
  for ( recID = 0; recID < nrecs; recID++ )
    {
      streamInqRecord(streamID1, &varID, &levelID);
      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
      offset   = gridsize*levelID;
      single1  = vardata1[varID] + offset;
      streamReadRecord(streamID1, single1, &nmiss);
      if ( nmiss ) cdoAbort("missing values unsupported for this operator!");
    }

  tsIDo = 0;
  while ( julval1 <= julval )
    {
      nrecs = streamInqTimestep(streamID1, tsID++);
      if ( nrecs == 0 ) break;
      julval2 = encode_julval(0, taxisInqVdate(taxisID1), taxisInqVtime(taxisID1));

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  recVarID[recID]   = varID;
	  recLevelID[recID] = levelID;

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  offset   = gridsize*levelID;
	  single2  = vardata2[varID] + offset;
	  streamReadRecord(streamID1, single2, &nmiss);
	  if ( nmiss ) cdoAbort("missing values unsupported for this operator!");
	}

      while ( julval < julval2 )
	{
	  if ( julval >= julval1 && julval < julval2 )
	    {
	      decode_julval(0, julval, &vdate, &vtime);
	      taxisDefVdate(taxisID2, vdate);
	      taxisDefVtime(taxisID2, vtime);
	      streamDefTimestep(streamID2, tsIDo++);

	      fac1 = (julval2-julval) / (julval2-julval1);
	      fac2 = (julval-julval1) / (julval2-julval1);
	      /*
	      printf("      %f %f %f %f %f\n", julval, julval1, julval2, fac1, fac2);
	      */
	      for ( recID = 0; recID < nrecs; recID++ )
		{
		  varID    = recVarID[recID];
		  levelID  = recLevelID[recID];
		  gridID1  = vlistInqVarGrid(vlistID1, varID);
		  gridsize = gridInqSize(gridID1);
		  offset   = gridsize*levelID;
		  single1  = vardata1[varID] + offset;
		  single2  = vardata2[varID] + offset;

		  for ( i = 0; i < gridsize; i++ )
		    array[i] = single1[i]*fac1 + single2[i]*fac2;

		  intgrid(gridID1, array, gridID2, &point);

		  streamDefRecord(streamID2, varID, levelID);
		  streamWriteRecord(streamID2, &point, nmiss);
		}
	    }
	  if ( readnextpos(fp, &julval, &xpos, &ypos) == EOF ) break;
	  gridDefXvals(gridID2, &xpos);
	  gridDefYvals(gridID2, &ypos);
	}

      julval1 = julval2;
      for ( varID = 0; varID < nvars; varID++ )
	{
	  vardatap        = vardata1[varID];
	  vardata1[varID] = vardata2[varID];
	  vardata2[varID] = vardatap;
	}
    }

  fclose(fp);
  streamClose(streamID2);
  streamClose(streamID1);

  for ( varID = 0; varID < nvars; varID++ )
    {
      free(vardata1[varID]);
      free(vardata2[varID]);
    }
  free(vardata1);
  free(vardata2);
  if ( array )  free(array);

  cdoFinish();

  return (0);
}

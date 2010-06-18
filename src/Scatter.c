/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2010 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

static
void genGrids(int gridID1, int *gridIDs, int nxvals, int nyvals, int nsx, int nsy,
	      int **gridindex, int nsplit, int nvals)
{
  static char *func = "genGrids";
  int gridID2;
  int gridtype;
  int gridsize, nx, ny;
  int index, i, j, ix, iy;
  double *xvals = NULL, *yvals = NULL;

  gridtype = gridInqType(gridID1);
  if ( !(gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_GENERIC) )
    cdoAbort("Unsupported grid type: %s!", gridNamePtr(gridtype));

  gridsize = gridInqSize(gridID1);
  nx = gridInqXsize(gridID1);
  ny = gridInqYsize(gridID1);

  xvals = (double *) malloc(nx*sizeof(double));
  yvals = (double *) malloc(ny*sizeof(double));

  gridInqXvals(gridID1, xvals);
  gridInqYvals(gridID1, yvals);

  index = 0;
  for ( iy = 0; iy < nsy; ++iy )
    for ( ix = 0; ix < nsx; ++ix )
      {
	for ( j = 0; j < nyvals; ++j )
	  for ( i = 0; i < nxvals; ++i )
	    {
	      
	    }
      }

  for ( index = 0; index < nsplit; ++index )
    gridIDs[index] = gridDuplicate(gridID2);

  free(xvals);
  free(yvals);
}

typedef struct
{
  int gridID;
  int **gridindex;
} sgrid_t;


void *Scatter(void *argument)
{
  static char func[] = "Scatter";
  int nchars;
  int streamID1, streamID2;
  int *gridIDs = NULL, *vlistIDs = NULL, *streamIDs = NULL;
  int gridID1 = -1, varID;
  int nrecs, ngrids, nvals;
  int tsID, recID, levelID;
  int varID2, levelID2;
  int vlistID1, vlistID2;
  char filesuffix[32];
  char filename[4096];
  int index;
  int nsplit, nsx, nsy, ix, iy;
  int xinc, yinc;
  int gridsize;
  int gridtype = -1;
  int nmiss;
  int nx, ny, i;
  double *array = NULL;
  sgrid_t *grids;

  cdoInitialize(argument);

  operatorInputArg("xinc, yinc");
  operatorCheckArgc(2);
  xinc = atoi(operatorArgv()[0]);
  yinc = atoi(operatorArgv()[1]);

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);

  ngrids = vlistNgrids(vlistID1);

  for ( index = 0; index < ngrids; index++ )
    {
      gridID1 = vlistGrid(vlistID1, index);
      gridtype = gridInqType(gridID1);
      if ( gridtype == GRID_LONLAT   || gridtype == GRID_GAUSSIAN ||
	   (gridtype == GRID_GENERIC && gridInqXsize(gridID1) > 0 && gridInqYsize(gridID1) > 0) )
	   break;
    }

  if ( index == ngrids )
    cdoAbort("No Lon/Lat, Gaussian or Generic grid found (%s data unsupported)!",
	     gridNamePtr(gridtype));

  gridID1 = vlistGrid(vlistID1, 0);
  gridsize = gridInqSize(gridID1);
  nx = gridInqXsize(gridID1);
  ny = gridInqXsize(gridID1);
  for ( i = 1; i < ngrids; i++ )
    {
      gridID1 = vlistGrid(vlistID1, i);
      if ( gridsize != gridInqSize(gridID1) )
	cdoAbort("Gridsize must not change!");
    }

  array = (double *) malloc(gridsize*sizeof(double));
  
  if ( nx%xinc != 0 ) cdoAbort("xsize is not multiple of xinc!");
  nsx = nx/xinc;
  if ( ny%yinc != 0 ) cdoAbort("ysize is not multiple of yinc!");
  nsy = ny/yinc;

  nsplit = nsx*nsy;
  nvals  = gridsize/nsplit;

  gridIDs   = (int *) malloc(nsplit*sizeof(int));
  vlistIDs  = (int *) malloc(nsplit*sizeof(int));
  streamIDs = (int *) malloc(nsplit*sizeof(int));

  grids = (sgrid_t *) malloc(ngrids*sizeof(sgrid_t));
  for ( i = 0; i < ngrids; i++ )
    {  
      gridID1 = vlistGrid(vlistID1, i);
      grids[i].gridID = vlistGrid(vlistID1, i);
      grids[i].gridindex = (int **) malloc(nsplit*sizeof(int*));
      for ( index = 0; index < nsplit; index++ )
	grids[i].gridindex[index] = (int *) malloc(nvals*sizeof(int));
    }

  for ( i = 0; i < ngrids; i++ )
    {
      gridID1 = vlistGrid(vlistID1, i);
      genGrids(gridID1, gridIDs, xinc, yinc, nsx, nsy, grids[i].gridindex, nsplit, nvals);
      //   vlistChangeGridIndex(vlistIDs[index], index, gridIDs[index]);
    }


  strcpy(filename, cdoStreamName(1));
  nchars = strlen(filename);

  filesuffix[0] = 0;
  if ( cdoDisableFilesuffix == FALSE )
    {
      strcat(filesuffix, streamFilesuffix(cdoDefaultFileType));
      if ( cdoDefaultFileType == FILETYPE_GRB )
	if ( vlistIsSzipped(vlistID1) || cdoZtype == COMPRESS_SZIP )
	  strcat(filesuffix, ".sz");
    }


  for ( index = 0; index < nsplit; index++ )
    {
      iy = index/nsx;
      ix = index - iy*nsx;
      printf("index %d, ix %d, iy %d\n", index, ix, iy);

      vlistIDs[index] = vlistDuplicate(vlistID1);

      
      sprintf(filename+nchars, "%05d", index);
      if ( filesuffix[0] )
	sprintf(filename+nchars+5, "%s", filesuffix);
      streamIDs[index] = streamOpenWrite(filename, cdoFiletype());
      if ( streamIDs[index] < 0 ) cdiError(streamIDs[index], "Open failed on %s", filename);

      streamDefVlist(streamIDs[index], vlistIDs[index]);
    }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      for ( index = 0; index < nsplit; index++ )
	streamDefTimestep(streamIDs[index], tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array, &nmiss);

	  for ( index = 0; index < nsplit; index++ )
	    {
	      streamDefRecord(streamIDs[index], varID2, levelID2);
	      streamWriteRecord(streamIDs[index], array, nmiss);
	    }
	}

      tsID++;
    }

  streamClose(streamID1);

  for ( index = 0; index < nsplit; index++ )
    {
      streamClose(streamIDs[index]);
      vlistDestroy(vlistIDs[index]);
      //   gridDestroy(gridIDs[index]);
    }

  if ( array ) free(array);

  if ( gridIDs   ) free(gridIDs);
  if ( vlistIDs  ) free(vlistIDs);
  if ( streamIDs ) free(streamIDs);

  for ( i = 0; i < ngrids; i++ )
    {
      for ( index = 0; index < nsplit; index++ )
	free(grids[i].gridindex[index]);
      free(grids[i].gridindex);
    }

  free(grids);

  cdoFinish();

  return (0);
}

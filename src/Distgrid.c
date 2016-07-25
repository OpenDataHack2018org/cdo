/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2016 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
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

#define  MAX_BLOCKS  65536

static
void genGrids(int gridID1, int *gridIDs, int nxvals, int nyvals, int nxblocks, int nyblocks,
	      int **gridindex, int *ogridsize, int nsplit)
{
  int gridID2;
  int gridsize2;
  int i, j, ix, iy, offset;

  int gridtype = gridInqType(gridID1);
  int lregular = TRUE;
  if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_GENERIC )
    lregular = TRUE;
  else if ( gridtype == GRID_CURVILINEAR )
    lregular = FALSE;
  else
    cdoAbort("Unsupported grid type: %s!", gridNamePtr(gridtype));

  int nx = gridInqXsize(gridID1);
  int ny = gridInqYsize(gridID1);

  bool lxcoord = true;
  bool lycoord = true;
  if ( gridInqXvals(gridID1, NULL) == 0 ) lxcoord = false;
  if ( gridInqYvals(gridID1, NULL) == 0 ) lycoord = false;

  double *xvals = NULL, *yvals = NULL;
  double *xvals2 = NULL, *yvals2 = NULL;

  if ( lxcoord )
    {
      if ( lregular )
        {
          xvals = (double*) Malloc(nx*sizeof(double));
        }
      else
        {
          xvals = (double*) Malloc(nx*ny*sizeof(double));
          xvals2 = (double*) Malloc(nxvals*nyvals*sizeof(double));
        }
      
      gridInqXvals(gridID1, xvals);
    }

  if ( lycoord )
    {
      if ( lregular )
        {
          yvals = (double*) Malloc(ny*sizeof(double));
        }
      else
        {
          yvals = (double*) Malloc(nx*ny*sizeof(double));
          yvals2 = (double*) Malloc(nxvals*nyvals*sizeof(double));
        }

      gridInqYvals(gridID1, yvals);
    }
 
  int *xlsize = (int*) Malloc(nxblocks*sizeof(int));
  int *ylsize = (int*) Malloc(nyblocks*sizeof(int));

  for ( ix = 0; ix < nxblocks; ++ix ) xlsize[ix] = nxvals;
  if ( nx%nxblocks != 0 ) xlsize[nxblocks-1] = nx - (nxblocks-1)*nxvals;
  if ( cdoVerbose ) for ( ix = 0; ix < nxblocks; ++ix ) cdoPrint("xblock %d: %d", ix, xlsize[ix]);

  for ( iy = 0; iy < nyblocks; ++iy ) ylsize[iy] = nyvals;
  if ( ny%nyblocks != 0 ) ylsize[nyblocks-1] = ny - (nyblocks-1)*nyvals;
  if ( cdoVerbose ) for ( iy = 0; iy < nyblocks; ++iy ) cdoPrint("yblock %d: %d", iy, ylsize[iy]);

  int index = 0;
  for ( iy = 0; iy < nyblocks; ++iy )
    for ( ix = 0; ix < nxblocks; ++ix )
      {
	offset = iy*nyvals*nx + ix*nxvals;

	gridsize2 = xlsize[ix]*ylsize[iy];
	gridindex[index] = (int*) Malloc(gridsize2*sizeof(int));

	gridsize2 = 0;
        // printf("iy %d, ix %d offset %d\n", iy, ix,  offset);
	for ( j = 0; j < ylsize[iy]; ++j )
	  {
	    for ( i = 0; i < xlsize[ix]; ++i )
	      {
	       	// printf(">> %d %d %d\n", j, i, offset + j*nx + i);
                if ( !lregular )
                  {
                    if ( lxcoord ) xvals2[gridsize2] = xvals[offset + j*nx + i];
                    if ( lycoord ) yvals2[gridsize2] = yvals[offset + j*nx + i];
                  }
		gridindex[index][gridsize2++] = offset + j*nx + i;
	      }
	  }
	// printf("gridsize2 %d\n", gridsize2);

	gridID2 = gridCreate(gridtype, gridsize2);
	gridDefXsize(gridID2, xlsize[ix]);
	gridDefYsize(gridID2, ylsize[iy]);

        if ( lregular )
          {
            if ( lxcoord ) gridDefXvals(gridID2, xvals+ix*nxvals);
            if ( lycoord ) gridDefYvals(gridID2, yvals+iy*nyvals);
          }
        else
          {
            if ( lxcoord ) gridDefXvals(gridID2, xvals2);
            if ( lycoord ) gridDefYvals(gridID2, yvals2);
          }
        
	gridIDs[index] = gridID2;
	ogridsize[index] = gridsize2;

	index++;
	if ( index > nsplit )
	  cdoAbort("Internal problem, index exceeded bounds!");
      }

  if ( xvals2 ) Free(xvals2);
  if ( yvals2 ) Free(yvals2);
  Free(xvals);
  Free(yvals);
  Free(xlsize);
  Free(ylsize);
}

static
void window_cell(double *array1, double *array2, long gridsize2, int *cellidx)
{
  for ( long i = 0; i < gridsize2; ++i )
    array2[i] = array1[cellidx[i]];
}

typedef struct
{
  int gridID;
  int *gridIDs;
  int *gridsize;
  int **gridindex;
} sgrid_t;


void *Distgrid(void *argument)
{
  int gridID1;
  int varID, levelID;
  int nrecs;
  char filesuffix[32];
  char filename[8192];
  const char *refname;
  int index;
  int gridtype = -1;
  int nmiss;
  int i;
  double missval;

  cdoInitialize(argument);

  operatorInputArg("nxblocks, [nyblocks]");
  if ( operatorArgc() < 1 ) cdoAbort("Too few arguments!");
  if ( operatorArgc() > 2 ) cdoAbort("Too many arguments!");
  int nxblocks = parameter2int(operatorArgv()[0]);
  int nyblocks = 1;
  if ( operatorArgc() == 2 ) nyblocks = parameter2int(operatorArgv()[1]);

  if ( nxblocks <= 0 ) cdoAbort("nxblocks has to be greater than 0!");
  if ( nyblocks <= 0 ) cdoAbort("nyblocks has to be greater than 0!");

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);

  int ngrids = vlistNgrids(vlistID1);

  for ( index = 0; index < ngrids; index++ )
    {
      gridID1 = vlistGrid(vlistID1, index);
      gridtype = gridInqType(gridID1);
      if ( gridtype == GRID_LONLAT   || gridtype == GRID_GAUSSIAN || gridtype == GRID_CURVILINEAR ||
	  (gridtype == GRID_GENERIC && gridInqXsize(gridID1) > 0 && gridInqYsize(gridID1) > 0) )
	   break;
    }

  if ( index == ngrids )
    cdoAbort("No Lon/Lat, Gaussian, curvilinear or generic grid found (%s data unsupported)!", gridNamePtr(gridtype));

  gridID1 = vlistGrid(vlistID1, 0);
  int gridsize = gridInqSize(gridID1);
  int nx = gridInqXsize(gridID1);
  int ny = gridInqYsize(gridID1);
  for ( i = 1; i < ngrids; i++ )
    {
      gridID1 = vlistGrid(vlistID1, i);
      if ( gridsize != gridInqSize(gridID1) )
	cdoAbort("Gridsize must not change!");
    }

  if ( nxblocks > nx )
    {
      cdoPrint("nxblocks (%d) greater than nx (%d), set to %d!", nxblocks, nx, nx);
      nxblocks = nx;
    }
  if ( nyblocks > ny )
    {
      cdoPrint("nyblocks (%d) greater than ny (%d), set to %d!", nyblocks, ny, ny);
      nyblocks = ny;
    }

  int xinc = nx/nxblocks;
  int yinc = ny/nyblocks;

  if ( nx%xinc != 0 ) xinc++;
  if ( ny%yinc != 0 ) yinc++;

  int nsplit = nxblocks*nyblocks;
  if ( nsplit > MAX_BLOCKS ) cdoAbort("Too many blocks (max = %d)!", MAX_BLOCKS);

  double *array1 = (double*) Malloc(gridsize*sizeof(double));

  int *vlistIDs  = (int*) Malloc(nsplit*sizeof(int));
  int *streamIDs = (int*) Malloc(nsplit*sizeof(int));

  sgrid_t *grids = (sgrid_t*) Malloc(ngrids*sizeof(sgrid_t));
  for ( i = 0; i < ngrids; i++ )
    {  
      grids[i].gridID    = vlistGrid(vlistID1, i);
      grids[i].gridIDs   = (int*) Malloc(nsplit*sizeof(int));
      grids[i].gridsize  = (int*) Malloc(nsplit*sizeof(int));
      grids[i].gridindex = (int**) Malloc(nsplit*sizeof(int*));

      for ( index = 0; index < nsplit; index++ ) grids[i].gridindex[index] = NULL;
    }

  for ( index = 0; index < nsplit; index++ )
    vlistIDs[index] = vlistDuplicate(vlistID1);

  if ( cdoVerbose ) cdoPrint("ngrids=%d  nsplit=%d", ngrids, nsplit);

  for ( i = 0; i < ngrids; i++ )
    {
      gridID1 = vlistGrid(vlistID1, i);
      genGrids(gridID1, grids[i].gridIDs, xinc, yinc, nxblocks, nyblocks, grids[i].gridindex, grids[i].gridsize, nsplit);
      /*
      if ( cdoVerbose )
	for ( index = 0; index < nsplit; index++ )
	  cdoPrint("Block %d,  gridID %d,  gridsize %d", index+1, grids[i].gridIDs[index], gridInqSize(grids[i].gridIDs[index]));
      */
      for ( index = 0; index < nsplit; index++ )
	vlistChangeGridIndex(vlistIDs[index], i, grids[i].gridIDs[index]);
    }

  int gridsize2max = 0;
  for ( index = 0; index < nsplit; index++ )
    if ( grids[0].gridsize[index] > gridsize2max ) gridsize2max = grids[0].gridsize[index];

  double *array2 = (double*) Malloc(gridsize2max*sizeof(double));

  strcpy(filename, cdoStreamName(1)->args);
  int nchars = strlen(filename);

  refname = cdoStreamName(0)->argv[cdoStreamName(0)->argc-1];
  filesuffix[0] = 0;
  cdoGenFileSuffix(filesuffix, sizeof(filesuffix), streamInqFiletype(streamID1), vlistID1, refname);

  for ( index = 0; index < nsplit; index++ )
    {
      sprintf(filename+nchars, "%05d", index);
      if ( filesuffix[0] )
	sprintf(filename+nchars+5, "%s", filesuffix);

      argument_t *fileargument = file_argument_new(filename);
      streamIDs[index] = streamOpenWrite(fileargument, cdoFiletype());
      file_argument_free(fileargument);

      streamDefVlist(streamIDs[index], vlistIDs[index]);
    }

  if ( ngrids > 1 ) cdoPrint("Bausstelle: number of different grids > 1!");
  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      for ( index = 0; index < nsplit; index++ )
	streamDefTimestep(streamIDs[index], tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array1, &nmiss);

	  missval = vlistInqVarMissval(vlistID1, varID);

	  for ( index = 0; index < nsplit; index++ )
	    {
	      i = 0;
	      window_cell(array1, array2, grids[i].gridsize[index], grids[i].gridindex[index]);
	      streamDefRecord(streamIDs[index], varID, levelID);
	      if ( nmiss > 0 )
		{
		  nmiss = 0;
		  for ( int k = 0; k < grids[i].gridsize[index]; ++k )
		    if ( DBL_IS_EQUAL(array2[k], missval) ) nmiss++;
		}
	      streamWriteRecord(streamIDs[index], array2, nmiss);
	    }
	}

      tsID++;
    }

  streamClose(streamID1);

  for ( index = 0; index < nsplit; index++ )
    {
      streamClose(streamIDs[index]);
      vlistDestroy(vlistIDs[index]);
    }

  if ( array1 ) Free(array1);
  if ( array2 ) Free(array2);

  if ( vlistIDs  ) Free(vlistIDs);
  if ( streamIDs ) Free(streamIDs);

  for ( i = 0; i < ngrids; i++ )
    {
      for ( index = 0; index < nsplit; index++ )
	gridDestroy(grids[i].gridIDs[index]);
      Free(grids[i].gridIDs);
      Free(grids[i].gridsize);

      for ( index = 0; index < nsplit; index++ )
        Free(grids[i].gridindex[index]);
      Free(grids[i].gridindex);
    }
  Free(grids);

  cdoFinish();

  return 0;
}

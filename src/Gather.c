/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2014 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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
#include "util.h"


typedef struct
{
  int streamID;
  int vlistID;
  int gridID;
  double *array;
} ens_file_t;


typedef struct
{
  double x, y;
  int id;
}
xyinfo_t;

static
int cmpx(const void *s1, const void *s2)
{
  int cmp = 0;
  const xyinfo_t *xy1 = (const xyinfo_t *) s1;
  const xyinfo_t *xy2 = (const xyinfo_t *) s2;

  if      ( xy1->x < xy2->x ) cmp = -1;
  else if ( xy1->x > xy2->x ) cmp =  1;

  return (cmp);
}

static
int cmpxy_lt(const void *s1, const void *s2)
{
  int cmp = 0;
  const xyinfo_t *xy1 = (const xyinfo_t *) s1;
  const xyinfo_t *xy2 = (const xyinfo_t *) s2;

  if      ( xy1->y < xy2->y || (!(fabs(xy1->y - xy2->y) > 0) && xy1->x < xy2->x) ) cmp = -1;
  else if ( xy1->y > xy2->y || (!(fabs(xy1->y - xy2->y) > 0) && xy1->x > xy2->x) ) cmp =  1;

  return (cmp);
}

static
int cmpxy_gt(const void *s1, const void *s2)
{
  int cmp = 0;
  const xyinfo_t *xy1 = (const xyinfo_t *) s1;
  const xyinfo_t *xy2 = (const xyinfo_t *) s2;

  if      ( xy1->y > xy2->y || (!(fabs(xy1->y - xy2->y) > 0) && xy1->x < xy2->x) ) cmp = -1;
  else if ( xy1->y < xy2->y || (!(fabs(xy1->y - xy2->y) > 0) && xy1->x > xy2->x) ) cmp =  1;

  return (cmp);
}

static
int genGrid(int nfiles, ens_file_t *ef, int **gridindex, int igrid)
{
  int lsouthnorth = TRUE;
  int fileID;
  int gridID;
  int gridID2 = -1;
  int gridtype = -1;
  int *xsize, *ysize;
  int *xoff, *yoff;
  int xsize2, ysize2;
  int idx;
  int nx, ny, ix, iy, i, j, ij, offset;
  double **xvals, **yvals;
  double *xvals2, *yvals2;
  xyinfo_t *xyinfo;

  gridID   = vlistGrid(ef[0].vlistID, igrid);
  gridtype = gridInqType(gridID);
  if ( gridtype == GRID_GENERIC && gridInqXsize(gridID) == 0 && gridInqYsize(gridID) == 0 )
    return (gridID2);

  xsize = (int*) malloc(nfiles*sizeof(int));
  ysize = (int*) malloc(nfiles*sizeof(int));
  xyinfo = (xyinfo_t*) malloc(nfiles*sizeof(xyinfo_t));
  xvals = (double**) malloc(nfiles*sizeof(double*));
  yvals = (double**) malloc(nfiles*sizeof(double*));

  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      gridID   = vlistGrid(ef[fileID].vlistID, igrid);
      gridtype = gridInqType(gridID);
      if ( !(gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN ||
	    (gridtype == GRID_GENERIC && gridInqXsize(gridID) > 0 && gridInqYsize(gridID) > 0)) )
	cdoAbort("Unsupported grid type: %s!", gridNamePtr(gridtype));

      xsize[fileID] = gridInqXsize(gridID);
      ysize[fileID] = gridInqYsize(gridID);
      /*
      if ( xsize == 0 ) xsize = gridInqXsize(gridID);
      if ( ysize == 0 ) ysize = gridInqYsize(gridID);
      if ( xsize != gridInqXsize(gridID) ) cdoAbort("xsize differ!");
      if ( ysize != gridInqYsize(gridID) ) cdoAbort("ysize differ!");
      */
      xvals[fileID] = (double*) malloc(xsize[fileID]*sizeof(double));
      yvals[fileID] = (double*) malloc(ysize[fileID]*sizeof(double));
      gridInqXvals(gridID, xvals[fileID]);
      gridInqYvals(gridID, yvals[fileID]);

      // printf("fileID %d, gridID %d\n", fileID, gridID);

      xyinfo[fileID].x  = xvals[fileID][0];
      xyinfo[fileID].y  = yvals[fileID][0];
      xyinfo[fileID].id = fileID;

      if ( ysize[fileID] > 1 )
	{
	  if ( yvals[fileID][0] > yvals[fileID][ysize[fileID]-1] ) lsouthnorth = FALSE;
	}
    }
  /*
  if ( cdoVerbose )
    for ( fileID = 0; fileID < nfiles; fileID++ )
      printf("1 %d %g %g \n",  xyinfo[fileID].id, xyinfo[fileID].x, xyinfo[fileID].y);
  */
  qsort(xyinfo, nfiles, sizeof(xyinfo_t), cmpx);  	      
  /*
  if ( cdoVerbose )
    for ( fileID = 0; fileID < nfiles; fileID++ )
      printf("2 %d %g %g \n",  xyinfo[fileID].id, xyinfo[fileID].x, xyinfo[fileID].y);
  */
  if ( lsouthnorth )
    qsort(xyinfo, nfiles, sizeof(xyinfo_t), cmpxy_lt);  
  else
    qsort(xyinfo, nfiles, sizeof(xyinfo_t), cmpxy_gt);  	      
  /*
  if ( cdoVerbose )
    for ( fileID = 0; fileID < nfiles; fileID++ )
      printf("3 %d %g %g \n",  xyinfo[fileID].id, xyinfo[fileID].x, xyinfo[fileID].y);
  */
  nx = 1;
  for ( fileID = 1; fileID < nfiles; fileID++ )
    {
      if ( DBL_IS_EQUAL(xyinfo[0].y, xyinfo[fileID].y) ) nx++;
      else break;
    }
  ny = nfiles/nx;
  if ( cdoVerbose ) cdoPrint("nx %d  ny %d", nx, ny);

  xsize2 = 0;
  for ( i = 0; i < nx; ++i ) xsize2 += xsize[xyinfo[i].id];
  ysize2 = 0;
  for ( j = 0; j < ny; ++j ) ysize2 += ysize[xyinfo[j*nx].id];
  if ( cdoVerbose ) cdoPrint("xsize2 %d  ysize2 %d", xsize2, ysize2);

  xvals2 = (double*) malloc(xsize2*sizeof(double));
  yvals2 = (double*) malloc(ysize2*sizeof(double));

  xoff = (int*) malloc((nx+1)*sizeof(int));
  yoff = (int*) malloc((ny+1)*sizeof(int));

  xoff[0] = 0;
  for ( i = 0; i < nx; ++i )
    {
      idx = xyinfo[i].id;
      memcpy(xvals2+xoff[i], xvals[idx], xsize[idx]*sizeof(double));
      xoff[i+1] = xoff[i] + xsize[idx];
    }

  yoff[0] = 0;
  for ( j = 0; j < ny; ++j )
    {
      idx = xyinfo[j*nx].id;
      memcpy(yvals2+yoff[j], yvals[idx], ysize[idx]*sizeof(double));
      yoff[j+1] = yoff[j] + ysize[idx];
    }

  if ( gridindex != NULL )
    {
      for ( fileID = 0; fileID < nfiles; fileID++ )
	{
	  idx = xyinfo[fileID].id;
	  iy = fileID/nx;
	  ix = fileID - iy*nx;

          offset = yoff[iy]*xsize2 + xoff[ix];
	  /*
	  printf("fileID %d %d, iy %d, ix %d, offset %d\n",
		 fileID, xyinfo[fileID].id, iy, ix, offset);
	  */
	  ij = 0;
	  for ( j = 0; j < ysize[idx]; ++j )
	    for ( i = 0; i < xsize[idx]; ++i )
	      {
		gridindex[idx][ij++] = offset+j*xsize2+i;
	      }
	}
    }

  gridID2 = gridCreate(gridtype, xsize2*ysize2);
  gridDefXsize(gridID2, xsize2);
  gridDefYsize(gridID2, ysize2);
  gridDefXvals(gridID2, xvals2);
  gridDefYvals(gridID2, yvals2);

  free(xoff);
  free(yoff);
  free(xsize);
  free(ysize);
  free(xvals2);
  free(yvals2);

  char string[1024];
  string[0] = 0;
  gridID = vlistGrid(ef[0].vlistID, igrid);
  gridInqXname(gridID, string);
  gridDefXname(gridID2, string);
  gridInqYname(gridID, string);
  gridDefYname(gridID2, string);
  gridInqXlongname(gridID, string);
  gridDefXlongname(gridID2, string);
  gridInqYlongname(gridID, string);
  gridDefYlongname(gridID2, string);
  gridInqXunits(gridID, string);
  gridDefXunits(gridID2, string);
  gridInqYunits(gridID, string);
  gridDefYunits(gridID2, string);

  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      free(xvals[fileID]);
      free(yvals[fileID]);
    }
  free(xvals);
  free(yvals);
  free(xyinfo);

  return (gridID2);
}


void *Gather(void *argument)
{
  int i;
  int varID, recID;
  int nrecs, nrecs0;
  int levelID;
  int vlistID1, vlistID2;
  int nmiss;
  int taxisID1, taxisID2;
  double missval;
  int fileID;

  cdoInitialize(argument);
    
  int nfiles = cdoStreamCnt() - 1;

  const char *ofilename = cdoStreamName(nfiles)->args;

  if ( !cdoSilentMode && !cdoOverwriteMode )
    if ( fileExists(ofilename) )
      if ( !userFileOverwrite(ofilename) )
	cdoAbort("Outputfile %s already exists!", ofilename);

  ens_file_t *ef = (ens_file_t*) malloc(nfiles*sizeof(ens_file_t));

  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      ef[fileID].streamID = streamOpenRead(cdoStreamName(fileID));
      ef[fileID].vlistID  = streamInqVlist(ef[fileID].streamID);
    }

  int nvars = vlistNvars(ef[0].vlistID);
  int *vars  = (int*) malloc(nvars*sizeof(int));
  for ( varID = 0; varID < nvars; varID++ ) vars[varID] = FALSE;

  /* check that the contents is always the same */
  int cmpfunc;
  if ( nvars == 1 ) 
    cmpfunc = CMP_NAME | CMP_NLEVEL;
  else
    cmpfunc = CMP_NAME | CMP_NLEVEL;

  for ( fileID = 1; fileID < nfiles; fileID++ )
    vlistCompare(ef[0].vlistID, ef[fileID].vlistID, cmpfunc);

  vlistID1 = ef[0].vlistID;
  int gridsizemax = vlistGridsizeMax(vlistID1);
  int gridsize;
  for ( fileID = 1; fileID < nfiles; fileID++ )
    {
      gridsize = vlistGridsizeMax(ef[fileID].vlistID);
      if ( gridsize > gridsizemax ) gridsizemax = gridsize;
    }
  gridsize = gridsizemax;

  for ( fileID = 0; fileID < nfiles; fileID++ )
    ef[fileID].array = (double*) malloc(gridsizemax*sizeof(double));

  int ngrids = vlistNgrids(ef[0].vlistID);
  int *gridIDs = (int*) malloc(ngrids*sizeof(int));
  int **gridindex = (int **) malloc(nfiles*sizeof(int *));
  for ( fileID = 0; fileID < nfiles; fileID++ )
    gridindex[fileID] = (int*) malloc(gridsizemax*sizeof(int));

  int ginit = FALSE;
  for ( i = 0; i < ngrids; ++i )
    {
      if ( ginit == FALSE )
	{
	  gridIDs[i] = genGrid(nfiles, ef, gridindex, i);
	  if ( gridIDs[i] != -1 ) ginit = TRUE;
	}
      else
	gridIDs[i] = genGrid(nfiles, ef, NULL, i);
    }

  vlistID2 = vlistDuplicate(vlistID1);
  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int gridsize2 = 0;
  for ( i = 0; i < ngrids; ++i )
    {
      if ( gridIDs[i] != -1 ) 
	{
	  if ( gridsize2 == 0 ) gridsize2 = gridInqSize(gridIDs[i]);
	  if ( gridsize2 != gridInqSize(gridIDs[i]) ) cdoAbort("gridsize differ!");
	  vlistChangeGridIndex(vlistID2, i, gridIDs[i]);
	}
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      int gridID = vlistInqVarGrid(ef[0].vlistID, varID);

      for ( i = 0; i < ngrids; ++i )
	{
	  if ( gridIDs[i] != -1 ) 
	    {
	      if ( gridID == vlistGrid(ef[0].vlistID, i) )
	      vars[varID] = TRUE;
	      break;
	    }
	}
    }

  int streamID2 = streamOpenWrite(cdoStreamName(nfiles), cdoFiletype());
      
  streamDefVlist(streamID2, vlistID2);
	  
  double *array2 = (double*) malloc(gridsize2*sizeof(double));

  int tsID = 0;
  do
    {
      nrecs0 = streamInqTimestep(ef[0].streamID, tsID);
      for ( fileID = 1; fileID < nfiles; fileID++ )
	{
	  nrecs = streamInqTimestep(ef[fileID].streamID, tsID);
	  if ( nrecs != nrecs0 )
	    cdoAbort("Number of records at time step %d of %s and %s differ!", tsID+1, cdoStreamName(0)->args, cdoStreamName(fileID)->args);
	}

      taxisCopyTimestep(taxisID2, taxisID1);

      if ( nrecs0 > 0 ) streamDefTimestep(streamID2, tsID);
      
      for ( recID = 0; recID < nrecs0; recID++ )
	{
	  streamInqRecord(ef[0].streamID, &varID, &levelID);

	  missval = vlistInqVarMissval(vlistID1, varID);
	  for ( i = 0; i < gridsize2; i++ ) array2[i] = missval;

#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(fileID, nmiss, i)
#endif
	  for ( fileID = 0; fileID < nfiles; fileID++ )
	    {
	      int varIDx, levelIDx;
	      if ( fileID > 0 ) streamInqRecord(ef[fileID].streamID, &varIDx, &levelIDx);
	      streamReadRecord(ef[fileID].streamID, ef[fileID].array, &nmiss);

	      if ( vars[varID] )
		{
		  gridsize = gridInqSize(vlistInqVarGrid(ef[fileID].vlistID, varID));
		  for ( i = 0; i < gridsize; ++i )
		    array2[gridindex[fileID][i]] = ef[fileID].array[i];
		}
	    }

	  streamDefRecord(streamID2, varID, levelID);

	  if ( vars[varID] )
	    {
	      nmiss = 0;
	      for ( i = 0; i < gridsize2; i++ )
		if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss++;
	      
	      streamWriteRecord(streamID2, array2, nmiss);
	    }
	  else
	    streamWriteRecord(streamID2, ef[0].array, 0);
	}

      tsID++;
    }
  while ( nrecs0 > 0 );

  for ( fileID = 0; fileID < nfiles; fileID++ )
    streamClose(ef[fileID].streamID);

  streamClose(streamID2);

  for ( fileID = 0; fileID < nfiles; fileID++ )
    if ( ef[fileID].array ) free(ef[fileID].array);

  if ( ef ) free(ef);
  if ( array2 ) free(array2);

  free(gridIDs);
  if ( vars   ) free(vars);

  cdoFinish();

  return (0);
}

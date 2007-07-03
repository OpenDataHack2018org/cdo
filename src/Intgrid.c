/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2007 Uwe Schulzweida, schulzweida@dkrz.de
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

      Intgrid    interpolate     PINGO grid interpolation
      Intgrid    intgridbil      Bilinear grid interpolation
*/


#include <string.h>
#include <math.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "interpol.h"


int genBoxGrid(int gridID1, int xavg, int yavg)
{
  static char func[] = "genBoxGrid";
  int debug = 1;
  int i, j, i1;
  int gridID2, gridtype;
  int gridsize1, xsize1, ysize1;
  int gridsize2, xsize2, ysize2;
  double *xvals1, *yvals1, *xvals2, *yvals2;
  double *grid1_corner_lon = NULL, *grid1_corner_lat = NULL;
  double *grid2_corner_lon = NULL, *grid2_corner_lat = NULL;

  gridtype = gridInqType(gridID1);
  gridsize1 = gridInqSize(gridID1);
  xsize1 = gridInqXsize(gridID1);
  ysize1 = gridInqYsize(gridID1);

  if ( debug ) printf("grid1 %d %d %d\n", gridsize1, xsize1, ysize1);

  xsize2 = xsize1/xavg;
  ysize2 = ysize1/yavg;
  if ( xsize1%xavg ) xsize2++;
  if ( ysize1%yavg ) ysize2++;
  gridsize2 = xsize2*ysize2;

  if ( debug ) printf("grid2 %d %d %d\n", gridsize2, xsize2, ysize2);

  gridID2 = gridCreate(gridtype, gridsize2);
  gridDefXsize(gridID2, xsize2);
  gridDefYsize(gridID2, ysize2);

  if ( gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT )
    {
      xvals1 = (double *) malloc(xsize1*sizeof(double));
      yvals1 = (double *) malloc(ysize1*sizeof(double));
      xvals2 = (double *) malloc(xsize2*sizeof(double));
      yvals2 = (double *) malloc(ysize2*sizeof(double));
      gridInqXvals(gridID1, xvals1);
      gridInqYvals(gridID1, yvals1);

      if ( gridInqYbounds(gridID1, NULL) && gridInqXbounds(gridID1, NULL) )
	{
	  grid1_corner_lon = (double *) malloc(2*xsize1*sizeof(double));
	  grid1_corner_lat = (double *) malloc(2*ysize1*sizeof(double));
	  grid2_corner_lon = (double *) malloc(2*xsize2*sizeof(double));
	  grid2_corner_lat = (double *) malloc(2*ysize2*sizeof(double));
	  gridInqXbounds(gridID1, grid1_corner_lon);
	  gridInqYbounds(gridID1, grid1_corner_lat);
	}

      j = 0;
      for ( i = 0; i < xsize1; i += xavg )
	{
	  i1 = i+(xavg-1);
	  if ( i1 >= xsize1-1 ) i1 = xsize1-1; 
	  xvals2[j] = xvals1[i] + (xvals1[i1] - xvals1[i])/2;
	  if ( debug ) printf("x %d %d %d %g", i, i1, j, xvals2[j]);
	  if ( grid2_corner_lon )
	    {
	      grid2_corner_lon[2*j] = grid1_corner_lon[2*i];
	      grid2_corner_lon[2*j+1] = grid1_corner_lon[2*i1+1];
	      if ( debug ) printf(" %g %g", grid2_corner_lon[2*j], grid2_corner_lon[2*j+1]);
	    }
	  if ( debug ) printf("\n");
	  j++;
	}
      j = 0;
      for ( i = 0; i < ysize1; i += yavg )
	{
	  i1 = i+(yavg-1);
	  if ( i1 >= ysize1-1 ) i1 = ysize1-1; 
	  yvals2[j] = yvals1[i] + (yvals1[i1] - yvals1[i])/2;
	  if ( debug ) printf("y %d %d %d %g", i, i1, j, yvals2[j]);
	  if ( grid2_corner_lat )
	    {
	      grid2_corner_lat[2*j] = grid1_corner_lat[2*i];
	      grid2_corner_lat[2*j+1] = grid1_corner_lat[2*i1+1];
	      if ( debug ) printf(" %g %g", grid2_corner_lat[2*j], grid2_corner_lat[2*j+1]);
	    }
	  if ( debug ) printf("\n");
	  j++;
	}

      gridDefXvals(gridID2, xvals2);
      gridDefYvals(gridID2, yvals2);

      if ( grid2_corner_lon && grid2_corner_lat )
	{
	  gridDefNvertex(gridID2, 2);
	  gridDefXbounds(gridID2, grid2_corner_lon);
	  gridDefYbounds(gridID2, grid2_corner_lat);

	  free(grid2_corner_lon);
	  free(grid2_corner_lat);
	}
    }
  else
    {
      cdoAbort("Unsupported grid: %s", gridNamePtr(gridtype));
    }

  return gridID2;
}


void boxavg(FIELD *field1, FIELD *field2, int xavg, int yavg)
{
  static char func[] = "boxavg";
  int nlon1, nlat1;
  int nlon2, nlat2;
  int ilat, ilon;
  int gridID1, gridID2;
  int nmiss;
  double **xfield1;
  double *array1, *array2;
  double missval;
  /* static int index = 0; */

  gridID1  = field1->grid;
  gridID2 = field2->grid;
  array1   = field1->ptr;
  array2  = field2->ptr;
  missval   = field1->missval;

  nlon1 = gridInqXsize(gridID1);
  nlat1 = gridInqYsize(gridID1);

  nlon2 = gridInqXsize(gridID2);
  nlat2 = gridInqYsize(gridID2);

  xfield1 = (double **) malloc(nlat1*sizeof(double *));

  for ( ilat = 0; ilat < nlat1; ilat++ )
    xfield1[ilat] = array1 + ilat*nlon1;

    {
      int i, j, ii, jj, in;
      double **xfield2;

      xfield2 = (double **) malloc(nlat2 * sizeof(double *));

      for ( ilat = 0; ilat < nlat2; ilat++ )
	xfield2[ilat] = array2 + ilat*nlon2;

      for ( ilat = 0; ilat < nlat2; ilat++ )
	for ( ilon = 0; ilon < nlon2; ilon++ )
	  {
	    xfield2[ilat][ilon] = 0;

	    in = 0;
	    for ( j = 0; j < yavg; ++j )
	      {
		jj = ilat*yavg+j;
		if ( jj >= nlat1 ) break;
		for ( i = 0; i < xavg; ++i )
		  {
		    ii = ilon*xavg+i;
		    if ( ii >= nlon1 ) break;
		    in++;
		    xfield2[ilat][ilon] += xfield1[jj][ii];
		  }
	      }
	    xfield2[ilat][ilon] /= in;
	  }

      nmiss = 0;
      for ( i = 0; i < nlat2*nlon2; i++ )
	if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss++;

      field2->nmiss = nmiss;

      free(xfield2);
    }

  free(xfield1);
}



void *Intgrid(void *argument)
{
  static char func[] = "Intgrid";
  int INTGRID, INTPOINT, INTERPOLATE, BOXAVG;
  int operatorID;
  int streamID1, streamID2;
  int nrecs, ngrids;
  int index;
  int tsID, recID, varID, levelID;
  int gridsize;
  int vlistID1, vlistID2;
  int gridID1 = -1, gridID2 = -1;
  int nmiss;
  int xavg, yavg;
  double missval;
  double slon, slat;
  double *array1 = NULL, *array2 = NULL;
  FIELD field1, field2;
  int taxisID1, taxisID2;

  cdoInitialize(argument);

  INTGRID     = cdoOperatorAdd("intgridbil",  0, 0, NULL);
  INTPOINT    = cdoOperatorAdd("intpoint",    0, 0, NULL);
  INTERPOLATE = cdoOperatorAdd("interpolate", 0, 0, NULL);
  BOXAVG      = cdoOperatorAdd("boxavg",      0, 0, NULL);

  operatorID = cdoOperatorID();

  if ( operatorID == INTGRID || operatorID == INTERPOLATE )
    {
      operatorInputArg("grid description file or name");
      gridID2 = cdoDefineGrid(operatorArgv()[0]);
    }
  else if ( operatorID == INTPOINT )
    {
      operatorInputArg("longitude and latitude");
      operatorCheckArgc(2);
      slon = atof(operatorArgv()[0]);
      slat = atof(operatorArgv()[1]);
      gridID2 = gridCreate(GRID_LONLAT, 1);
      gridDefXsize(gridID2, 1);
      gridDefYsize(gridID2, 1);
      gridDefXvals(gridID2, &slon);
      gridDefYvals(gridID2, &slat);
    }
  else if ( operatorID == BOXAVG )
    {
      operatorInputArg("xavg, yavg");
      operatorCheckArgc(2);
      xavg = atoi(operatorArgv()[0]);
      yavg = atoi(operatorArgv()[1]);
    }

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  ngrids = vlistNgrids(vlistID1);
  for ( index = 0; index < ngrids; index++ )
    {
      gridID1 = vlistGrid(vlistID1, index);

      if ( operatorID == BOXAVG )
	{
	  if ( index == 0 )
	    {
	      if ( gridInqType(gridID1) != GRID_LONLAT && gridInqType(gridID1) != GRID_GAUSSIAN &&
		   gridInqType(gridID1) != GRID_CURVILINEAR   )
		cdoAbort("Interpolation of %s data unsupported!", gridNamePtr(gridInqType(gridID1)) );

	      gridID2 = genBoxGrid(gridID1, xavg, yavg);
	    }
	  else
	    cdoAbort("Too many different grids!");
	}
      else
	{
	  if ( gridInqType(gridID1) != GRID_LONLAT && gridInqType(gridID1) != GRID_GAUSSIAN )
	    cdoAbort("Interpolation of %s data unsupported!", gridNamePtr(gridInqType(gridID1)) );

	  if ( gridIsRotated(gridID1) )
	    cdoAbort("Rotated grids not supported!");
	}

      vlistChangeGridIndex(vlistID2, index, gridID2);
    }

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);
  array1   = (double *) malloc(gridsize*sizeof(double));

  gridsize = gridInqSize(gridID2);
  array2   = (double *) malloc(gridsize*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array1, &nmiss);

	  gridID1 = vlistInqVarGrid(vlistID1, varID);
	  missval = vlistInqVarMissval(vlistID1, varID);

	  field1.grid    = gridID1;
	  field1.nmiss   = nmiss;
	  field1.missval = missval;
	  field1.ptr     = array1;
	  field2.grid    = gridID2;
	  field2.ptr     = array2;
	  field2.nmiss   = 0;

	  if ( operatorID == INTGRID || operatorID == INTPOINT )
	    intgrid(&field1, &field2);
	  else if ( operatorID == INTERPOLATE )
	    interpolate(&field1, &field2);
	  else if ( operatorID == BOXAVG )
	    boxavg(&field1, &field2, xavg, yavg);

	  nmiss = field2.nmiss;

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, array2, nmiss);
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

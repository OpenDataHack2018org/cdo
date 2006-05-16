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

*/


#include <string.h>
#include <stdlib.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


static int gengrid(int gridID1, int lhalo, int rhalo)
{
  static char func[] = "gengrid";  
  int gridtype, gridID2;
  int nlon1, nlat1;
  int nlon2, nlat2;
  int i;
  int prec;
  int ilat, ilon;
  char xname[128], xlongname[128], xunits[128];
  char yname[128], ylongname[128], yunits[128];
  double *xvals1 = NULL, *yvals1 = NULL;
  double *xvals2 = NULL, *yvals2 = NULL;
  double *xbounds1 = NULL, *ybounds1 = NULL;
  double *xbounds2 = NULL, *ybounds2 = NULL;
  double *pxvals2 = NULL, *pyvals2 = NULL;
  double *pxbounds2 = NULL, *pybounds2 = NULL;

  nlon1 = gridInqXsize(gridID1);
  nlat1 = gridInqYsize(gridID1);

  nlon2 = nlon1 + lhalo + rhalo;
  nlat2 = nlat1;

  gridtype = gridInqType(gridID1);
  prec     = gridInqPrec(gridID1);

  gridID2 = gridCreate(gridtype, nlon2*nlat2);
  gridDefXsize(gridID2, nlon2);
  gridDefYsize(gridID2, nlat2);

  gridDefPrec(gridID2, prec);

  gridInqXname(gridID1, xname);
  gridInqXlongname(gridID1, xlongname);
  gridInqXunits(gridID1, xunits);
  gridInqYname(gridID1, yname);
  gridInqYlongname(gridID1, ylongname);
  gridInqYunits(gridID1, yunits);

  gridDefXname(gridID2, xname);
  gridDefXlongname(gridID2, xlongname);
  gridDefXunits(gridID2, xunits);
  gridDefYname(gridID2, yname);
  gridDefYlongname(gridID2, ylongname);
  gridDefYunits(gridID2, yunits);

  if ( gridInqXvals(gridID1, NULL) && gridInqYvals(gridID1, NULL) )
    {
      if ( gridtype == GRID_CURVILINEAR )
	{
	  xvals1 = (double *) malloc(nlon1*nlat1*sizeof(double));
	  yvals1 = (double *) malloc(nlon1*nlat1*sizeof(double));
	  xvals2 = (double *) malloc(nlon2*nlat2*sizeof(double));
	  yvals2 = (double *) malloc(nlon2*nlat2*sizeof(double));
	}
      else
	{
	  xvals1 = (double *) malloc(nlon1*sizeof(double));
	  yvals1 = (double *) malloc(nlat1*sizeof(double));
	  xvals2 = (double *) malloc(nlon2*sizeof(double));
	  yvals2 = (double *) malloc(nlat2*sizeof(double));
	}

      pxvals2 = xvals2;
      pyvals2 = yvals2;

      gridInqXvals(gridID1, xvals1);
      gridInqYvals(gridID1, yvals1);

      if ( gridtype == GRID_CURVILINEAR )
	{
	  for ( ilat = 0; ilat < nlat2; ilat++ )
	    {
	      for ( ilon = nlon1-lhalo; ilon < nlon1; ilon++ )
		{
		  *pxvals2++ = xvals1[ilat*nlon1 + ilon];
		  *pyvals2++ = yvals1[ilat*nlon1 + ilon];
		}

	      for ( ilon = 0; ilon < nlon1; ilon++ )
		{
		  *pxvals2++ = xvals1[ilat*nlon1 + ilon];
		  *pyvals2++ = yvals1[ilat*nlon1 + ilon];
		}

	      for ( ilon = 0; ilon < rhalo; ilon++ )
		{
		  *pxvals2++ = xvals1[ilat*nlon1 + ilon];
		  *pyvals2++ = yvals1[ilat*nlon1 + ilon];
		}
	    }
	}
      else
	{
	  for ( i = nlon1-lhalo; i < nlon1; i++ ) *pxvals2++ = xvals1[i] - 360;
	  for ( i = 0; i < nlon1; i++ ) *pxvals2++ = xvals1[i];
	  for ( i = 0; i < rhalo; i++ ) *pxvals2++ = xvals1[i] + 360;

	  for ( i = 0; i < nlat1; i++ ) yvals2[i] = yvals1[i];
	}
      /*
      for ( i = 0; i < nlat2; i++ ) printf("lat : %d %g\n", i+1, yvals2[i]);
      for ( i = 0; i < nlon2; i++ ) printf("lon : %d %g\n", i+1, xvals2[i]);
      */
      gridDefXvals(gridID2, xvals2);
      gridDefYvals(gridID2, yvals2);

      free(xvals1);
      free(yvals1);
      free(xvals2);
      free(yvals2);
    }

  if ( gridInqXbounds(gridID1, NULL) && gridInqYbounds(gridID1, NULL) )
    {
      if ( gridtype == GRID_CURVILINEAR )
	{
	  xbounds1 = (double *) malloc(4*nlon1*nlat1*sizeof(double));
	  ybounds1 = (double *) malloc(4*nlon1*nlat1*sizeof(double));
	  xbounds2 = (double *) malloc(4*nlon2*nlat2*sizeof(double));
	  ybounds2 = (double *) malloc(4*nlon2*nlat2*sizeof(double));
	}
      else
	{
	  xbounds1 = (double *) malloc(2*nlon1*sizeof(double));
	  ybounds1 = (double *) malloc(2*nlat1*sizeof(double));
	  xbounds2 = (double *) malloc(2*nlon2*sizeof(double));
	  ybounds2 = (double *) malloc(2*nlat2*sizeof(double));
	}

      pxbounds2 = xbounds2;
      pybounds2 = ybounds2;

      gridInqXbounds(gridID1, xbounds1);
      gridInqYbounds(gridID1, ybounds1);

      if ( gridtype == GRID_CURVILINEAR )
	{
	  gridDefNvertex(gridID2, 4);
	  for ( ilat = 0; ilat < nlat1; ilat++ )
	    {
	      for ( ilon = 4*(nlon1-lhalo); ilon < 4*nlon1; ilon++ )
		{
		  *pxbounds2++ = xbounds1[ilat*nlon1 + ilon];
		  *pybounds2++ = ybounds1[ilat*nlon1 + ilon];
		}

	      for ( ilon = 0; ilon < 4*nlon1; ilon++ )
		{
		  *pxbounds2++ = xbounds1[ilat*nlon1 + ilon];
		  *pybounds2++ = ybounds1[ilat*nlon1 + ilon];
		}

	      for ( ilon = 0; ilon < 4*rhalo; ilon++ )
		{
		  *pxbounds2++ = xbounds1[ilat*nlon1 + ilon];
		  *pybounds2++ = ybounds1[ilat*nlon1 + ilon];
		}
	    }
	}
      else
	{
	  gridDefNvertex(gridID2, 2);
	  for ( i = 2*(nlon1-lhalo); i < 2*nlon1; i++ ) *pxbounds2++ = xbounds1[i] - 360;
	  for ( i = 0; i < 2*nlon1; i++ ) *pxbounds2++ = xbounds1[i];
	  for ( i = 0; i < 2*rhalo; i++ ) *pxbounds2++ = xbounds1[i] + 360;

	  for ( i = 0; i < 2*nlat2; i++ ) ybounds2[i] = ybounds1[i];
	}

      gridDefXbounds(gridID2, xbounds2);
      gridDefYbounds(gridID2, ybounds2);

      free(xbounds1);
      free(ybounds1);
      free(xbounds2);
      free(ybounds2);
    }

  return (gridID2);
}


static int genindexgrid(int gridID1, int *lhalo, int *rhalo)
{
  int gridID2;
  int nlon1;

  operatorCheckArgc(2);

  *lhalo = atoi(operatorArgv()[0]);
  *rhalo = atoi(operatorArgv()[1]);

  nlon1 = gridInqXsize(gridID1);

  if ( *lhalo < 0 || *lhalo > nlon1 )
    {
      cdoWarning("left halo out of range. Set to 0.");
      *lhalo = 0;
    }

  if ( *rhalo < 0 || *rhalo > nlon1 )
    {
      cdoWarning("right halo out of range. Set to 0.");
      *rhalo = 0;
    }

  gridID2 = gengrid(gridID1, *lhalo, *rhalo);

  return (gridID2);
}


static void halo(double *array1, int gridID1, double *array2, int lhalo, int rhalo)
{
  int nlon1, nlat;
  int ilat, ilon;

  nlon1 = gridInqXsize(gridID1);
  nlat = gridInqYsize(gridID1);

  for ( ilat = 0; ilat < nlat; ilat++ )
    {
      for ( ilon = nlon1-lhalo; ilon < nlon1; ilon++ )
	*array2++ = array1[ilat*nlon1 + ilon];

      for ( ilon = 0; ilon < nlon1; ilon++ )
	*array2++ = array1[ilat*nlon1 + ilon];

      for ( ilon = 0; ilon < rhalo; ilon++ )
	*array2++ = array1[ilat*nlon1 + ilon];
    }
}


void *Sethalo(void *argument)
{
  static char func[] = "Sethalo";
  int streamID1, streamID2;
  int nrecs, nvars;
  int tsID, recID, varID, levelID;
  int gridsize, gridsize2;
  int vlistID1, vlistID2;
  int gridID1 = -1, gridID2;
  int index, ngrids, gridtype;
  int nmiss;
  int *vars;
  int i;
  int lhalo, rhalo;
  int ndiffgrids;
  double missval;
  double *array1 = NULL, *array2 = NULL;
  int taxisID1, taxisID2;

  cdoInitialize(argument);

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);

  ngrids = vlistNgrids(vlistID1);
  ndiffgrids = 0;
  for ( index = 1; index < ngrids; index++ )
    if ( vlistGrid(vlistID1, 0) != vlistGrid(vlistID1, index))
      ndiffgrids++;

  for ( index = 0; index < ngrids; index++ )
    {
      gridID1  = vlistGrid(vlistID1, index);
      gridtype = gridInqType(gridID1);
      if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN ) break;
      if ( gridtype == GRID_CURVILINEAR ) break;
      if ( gridtype == GRID_GENERIC &&
	   gridInqXsize(gridID1) > 0 && gridInqYsize(gridID1) > 0 ) break;
    }

  if ( gridInqType(gridID1) == GRID_GAUSSIAN_REDUCED )
    cdoAbort("Gaussian reduced grid found. Use option -R to convert it to a regular grid!");

  if ( index == ngrids ) cdoAbort("No regular grid found!");
  if ( ndiffgrids > 0 )  cdoAbort("Too many different grids!");

  operatorInputArg("left and right halo");
  gridID2 = genindexgrid(gridID1, &lhalo, &rhalo);

  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  ngrids = vlistNgrids(vlistID1);
  for ( index = 0; index < ngrids; index++ )
    {
      if ( gridID1 == vlistGrid(vlistID1, index) )
	{
	  vlistChangeGridIndex(vlistID2, index, gridID2);
	  break;
	}
    }

  nvars = vlistNvars(vlistID1);
  vars  = (int *) malloc(nvars*sizeof(int));
  for ( varID = 0; varID < nvars; varID++ )
    {
      if ( gridID1 == vlistInqVarGrid(vlistID1, varID) )
	vars[varID] = TRUE;
      else
	vars[varID] = FALSE;
    }

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  gridsize = gridInqSize(gridID1);
  array1 = (double *) malloc(gridsize*sizeof(double));

  gridsize2 = gridInqSize(gridID2);
  array2 = (double *) malloc(gridsize2*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  if ( vars[varID] )
	    {
	      streamReadRecord(streamID1, array1, &nmiss);

	      halo(array1, gridID1, array2, lhalo, rhalo);

	      if ( nmiss )
		{
		  nmiss = 0;
		  missval = vlistInqVarMissval(vlistID1, varID);
		  for ( i = 0; i < gridsize2; i++ )
		    if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss++;
		}

	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, array2, nmiss);
	    }
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( vars   ) free(vars);
  if ( array2 ) free(array2);
  if ( array1 ) free(array1);

  cdoFinish();

  return (0);
}

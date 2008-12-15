/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2008 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#if  defined  (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <stdarg.h> /* va_list */

#if defined (HAVE_LIBPROJ)
#  include "projects.h"
#endif

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "error.h"
#include "grid.h"


#define  deg2rad  (M_PI/180.)   /* conversion for deg to rad */
#define  rad2deg  (180./M_PI)   /* conversion for rad to deg */


int gridToZonal(int gridID1)
{
  static char func[] = "gridToZonal";
  int gridID2;
  int gridtype, gridsize;
  double  xval = 0;
  double *yvals;

  gridtype = gridInqType(gridID1);
  gridsize = gridInqYsize(gridID1);
  gridID2  = gridCreate(gridtype, gridsize);
	  
  if ( gridtype != GRID_LONLAT &&
       gridtype != GRID_GAUSSIAN &&
       (gridtype == GRID_GENERIC && gridsize <= 1) )
    {
      Error(func, "Gridtype %s unsupported!", gridNamePtr(gridtype));
    }
  else
    {
      gridDefXsize(gridID2, 1);
      gridDefYsize(gridID2, gridsize);

      gridDefXvals(gridID2, &xval);

      if ( gridInqYvals(gridID1, NULL) )
	{
	  yvals = (double *) malloc(gridsize*sizeof(double));

	  gridInqYvals(gridID1, yvals);
	  gridDefYvals(gridID2, yvals);

	  free(yvals);
	}
    }

  return (gridID2);
}


int gridToMeridional(int gridID1)
{
  static char func[] = "gridToMeridional";
  int gridID2;
  int gridtype, gridsize;
  double *xvals;
  double  yval = 0;

  gridtype = gridInqType(gridID1);
  gridsize = gridInqXsize(gridID1);
  gridID2  = gridCreate(gridtype, gridsize);
	  
  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
      {
	gridDefXsize(gridID2, gridsize);
	gridDefYsize(gridID2, 1);

	xvals = (double *) malloc(gridsize*sizeof(double));

	gridInqXvals(gridID1, xvals);
	gridDefXvals(gridID2, xvals);
	gridDefYvals(gridID2, &yval);

	free(xvals);

	break;
      }
    default:
      {
	Error(func, "Gridtype %s unsupported!", gridNamePtr(gridtype));
	break;
      }
    }

  return (gridID2);
}


void gridGenXbounds(int nx, double *xvals, double *xbounds)
{
  int i;

  for ( i = 0; i < nx-1; i++ )
    {
      xbounds[2*i+1]   = 0.5*(xvals[i] + xvals[i+1]);
      xbounds[2*(i+1)] = 0.5*(xvals[i] + xvals[i+1]);
    }

  xbounds[0]      = 2*xvals[0] - xbounds[1];
  xbounds[2*nx-1] = 2*xvals[nx-1] - xbounds[2*(nx-1)];
}


void gridGenYbounds(int ny, double *yvals, double *ybounds)
{
  int i;

  for ( i = 0; i < ny-1; i++ )
    {
      ybounds[2*i+1]   = 0.5*(yvals[i] + yvals[i+1]);
      ybounds[2*(i+1)] = 0.5*(yvals[i] + yvals[i+1]);
    }

  ybounds[0]      = 2*yvals[0] - ybounds[1];
  ybounds[2*ny-1] = 2*yvals[ny-1] - ybounds[2*(ny-1)];

  if ( yvals[0] > yvals[ny-1] )
    {
      if ( ybounds[0]      >  88 ) ybounds[0]      =  90;
      if ( ybounds[2*ny-1] < -88 ) ybounds[2*ny-1] = -90;
    }
  else
    {
      if ( ybounds[0]      < -88 ) ybounds[0]      = -90;
      if ( ybounds[2*ny-1] >  88 ) ybounds[2*ny-1] =  90;
    }
}


void gridGenYboundsM(int ny, double *yvals, double *ybounds)
{
  int i;

  for ( i = 0; i < ny-1; i++ )
    {
      ybounds[2*i+1]   = 0.5*(yvals[i] + yvals[i+1]);
      ybounds[2*(i+1)] = 0.5*(yvals[i] + yvals[i+1]);
    }

  ybounds[0]      = 2*yvals[0] - ybounds[1];
  ybounds[2*ny-1] = 2*yvals[ny-1] - ybounds[2*(ny-1)];
}


void gridGenRotBounds(int gridID, int nx, int ny,
		      double *xbounds, double *ybounds, double *xbounds2D, double *ybounds2D)
{
  int i, j, index;
  double minlon, maxlon;
  double minlat, maxlat;
  double xpole, ypole;

  xpole = gridInqXpole(gridID);
  ypole = gridInqYpole(gridID);

  for ( j = 0; j < ny; j++ )
    {
      if ( ybounds[0] > ybounds[1] )
	{
	  maxlat = ybounds[2*j];
	  minlat = ybounds[2*j+1];
	}
      else
	{
	  maxlat = ybounds[2*j+1];
	  minlat = ybounds[2*j];
	}

      for ( i = 0; i < nx; i++ )
	{
	  minlon = xbounds[2*i];
	  maxlon = xbounds[2*i+1];

	  index = j*4*nx + 4*i;
	  xbounds2D[index+0] = rls_to_rl(minlat, minlon, ypole, xpole);
	  xbounds2D[index+1] = rls_to_rl(minlat, maxlon, ypole, xpole);
	  xbounds2D[index+2] = rls_to_rl(maxlat, maxlon, ypole, xpole);
	  xbounds2D[index+3] = rls_to_rl(maxlat, minlon, ypole, xpole);

	  ybounds2D[index+0] = phs_to_ph(minlat, minlon, ypole);
	  ybounds2D[index+1] = phs_to_ph(minlat, maxlon, ypole);
	  ybounds2D[index+2] = phs_to_ph(maxlat, maxlon, ypole);
	  ybounds2D[index+3] = phs_to_ph(maxlat, minlon, ypole);
	}
    }
}


void gridGenXbounds2D(int nx, int ny, double *xbounds, double *xbounds2D)
{
  int i, j, index;
  double minlon, maxlon;

  for ( i = 0; i < nx; i++ )
    {
      minlon = xbounds[2*i];
      maxlon = xbounds[2*i+1];

      for ( j = 0; j < ny; j++ )
	{
	  index = j*4*nx + 4*i;
	  xbounds2D[index+0] = minlon;
	  xbounds2D[index+1] = maxlon;
	  xbounds2D[index+2] = maxlon;
	  xbounds2D[index+3] = minlon;
	}
    }
}


void gridGenYbounds2D(int nx, int ny, double *ybounds, double *ybounds2D)
{
  int i, j, index;
  double minlat, maxlat;

  for ( j = 0; j < ny; j++ )
    {
      if ( ybounds[0] > ybounds[1] )
	{
	  maxlat = ybounds[2*j];
	  minlat = ybounds[2*j+1];
	}
      else
	{
	  maxlat = ybounds[2*j+1];
	  minlat = ybounds[2*j];
	}

      for ( i = 0; i < nx; i++ )
	{
	  index = j*4*nx + 4*i;
	  ybounds2D[index+0] = minlat;
	  ybounds2D[index+1] = minlat;
	  ybounds2D[index+2] = maxlat;
	  ybounds2D[index+3] = maxlat;
	}
    }
}


char *gen_param(const char *fmt, ...)
{
  static char func[] = "get_param";
  va_list args;
  char str[256];
  char *rstr;
  int len;

  va_start(args, fmt);

  len = vsprintf(str, fmt, args);

  va_end(args);

  len++;
  rstr = (char *) malloc(len*sizeof(char));
  memcpy(rstr, str, len*sizeof(char));

  return (rstr);
}


int gridToCurvilinear(int gridID1)
{
  static char func[] = "gridToCurvilinear";
  int gridID2;
  int gridtype, gridsize;

  gridtype = gridInqType(gridID1);
  gridsize = gridInqSize(gridID1);
  gridID2  = gridCreate(GRID_CURVILINEAR, gridsize);
  gridDefPrec(gridID2, DATATYPE_FLT32);
	  
  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_LCC:
    case GRID_LAEA:
    case GRID_SINUSOIDAL:
      {
	int i, j, nx, ny;
	double *xvals = NULL, *yvals = NULL;
	double *xvals2D, *yvals2D;
	double *xbounds = NULL, *ybounds = NULL;
	double *xbounds2D, *ybounds2D;

	nx = gridInqXsize(gridID1);
	ny = gridInqYsize(gridID1);

	gridDefXsize(gridID2, nx);
	gridDefYsize(gridID2, ny);

	xvals2D = (double *) malloc(gridsize*sizeof(double));
	yvals2D = (double *) malloc(gridsize*sizeof(double));


	if ( gridtype == GRID_LCC )
	  {
	    double xi, xj;
	    double originLon, originLat, lonParY, lat1, lat2, xincm, yincm;
	    double zlat, zlon;
	    int projflag, scanflag;
	    int status;

	    gridInqLCC(gridID1, &originLon, &originLat, &lonParY, &lat1, &lat2, &xincm, &yincm,
		       &projflag, &scanflag);
	    /*
	    while ( originLon < 0 ) originLon += 360;
	    while ( lonParY   < 0 ) lonParY   += 360;
	    */
	    if ( ! DBL_IS_EQUAL(xincm, yincm) )
	      Warning(func, "X and Y increment must be equal on Lambert Conformal grid (Xinc = %g, Yinc = %g)\n", 
		      xincm, yincm);
	    if ( ! DBL_IS_EQUAL(lat1, lat2) )
	      Warning(func, "Lat1 and Lat2 must be equal on Lambert Conformal grid (Lat1 = %g, Lat2 = %g)\n", 
		      lat1, lat2);
		
	    for ( j = 0; j < ny; j++ )
	      for ( i = 0; i < nx; i++ )
		{
		  xi = i+1;
		  xj = j+1;
		  status = W3FB12(xi, xj, originLat, originLon, xincm, lonParY, lat1, &zlat, &zlon);
		  xvals2D[j*nx+i] = zlon;
		  yvals2D[j*nx+i] = zlat;
		}
	  }
	else
	  {
	    if ( ! (gridInqXvals(gridID1, NULL) && gridInqYvals(gridID1, NULL)) )
	      Error(func, "Grid has no values");

	    xvals = (double *) malloc(nx*sizeof(double));
	    yvals = (double *) malloc(ny*sizeof(double));

	    gridInqXvals(gridID1, xvals);
	    gridInqYvals(gridID1, yvals);

	    if ( gridIsRotated(gridID1) )
	      {
		double xpole, ypole;
		
		xpole = gridInqXpole(gridID1);
		ypole = gridInqYpole(gridID1);
		
		for ( j = 0; j < ny; j++ )
		  for ( i = 0; i < nx; i++ )
		    {
		      xvals2D[j*nx+i] = rls_to_rl(yvals[j], xvals[i], ypole, xpole);
		      yvals2D[j*nx+i] = phs_to_ph(yvals[j], xvals[i], ypole);
		    }	    
	      }
	    else
	      {
		if ( gridtype == GRID_SINUSOIDAL )
		  {
#if defined (HAVE_LIBPROJ)
		    int i;
		    PJ   *libProj;
		    char *params[20];
		    int nbpar=0;
		    projUV data, res;

		    nbpar = 0;
		    params[nbpar++] = (char*) "proj=sinu";
		    params[nbpar++] = (char*) "ellps=WGS84";

		    if ( cdoVerbose )
		      for ( i = 0; i < nbpar; ++i )
			cdoPrint("Proj.param[%d] = %s", i+1, params[i]);

		    libProj = pj_init(nbpar, params);
		    if ( !libProj )
		      cdoAbort("proj error: %s", pj_strerrno(pj_errno));

		    /* libProj->over = 1; */		// allow longitude > 180�

		    for ( j = 0; j < ny; j++ )
		      for ( i = 0; i < nx; i++ )
			{
			  data.u = xvals[i];
			  data.v = yvals[j];
			  res = pj_inv(data, libProj);
			  xvals2D[j*nx+i] = res.u*rad2deg;
			  yvals2D[j*nx+i] = res.v*rad2deg;
			}
#else
		    cdoAbort("proj4 support not compiled in!");
#endif
		  }
		else if ( gridtype == GRID_LAEA )
		  {
#if defined (HAVE_LIBPROJ)
		    int i;
		    PJ   *libProj;
		    char *params[20];
		    int nbpar=0;
		    projUV data, res;
		    double a, lon_0, lat_0;

		    gridInqLaea(gridID1, &a , &lon_0, &lat_0);

		    nbpar = 0;
		    params[nbpar++] = gen_param("proj=laea");
		    params[nbpar++] = gen_param("a=%g", a);
		    params[nbpar++] = gen_param("lon_0=%g", lon_0);
		    params[nbpar++] = gen_param("lat_0=%g", lat_0);

		    if ( cdoVerbose )
		      for ( i = 0; i < nbpar; ++i )
			cdoPrint("Proj.param[%d] = %s", i+1, params[i]);

		    libProj = pj_init(nbpar, &params[0]);
		    if ( !libProj )
		      cdoAbort("proj error: %s", pj_strerrno(pj_errno));

		    for ( i = 0; i < nbpar; ++i ) free(params[i]);

		    /* libProj->over = 1; */		// allow longitude > 180�

		    for ( j = 0; j < ny; j++ )
		      for ( i = 0; i < nx; i++ )
			{
			  data.u = xvals[i];
			  data.v = yvals[j];
			  res = pj_inv(data, libProj);
			  xvals2D[j*nx+i] = res.u*rad2deg;
			  yvals2D[j*nx+i] = res.v*rad2deg;
			}
#else
		    cdoAbort("proj4 support not compiled in!");
#endif
		  }
		else
		  {
		    for ( j = 0; j < ny; j++ )
		      for ( i = 0; i < nx; i++ )
			{
			  xvals2D[j*nx+i] = xvals[i];
			  yvals2D[j*nx+i] = yvals[j];
			}
		  }
	      }
	  }

	gridDefXvals(gridID2, xvals2D);
	gridDefYvals(gridID2, yvals2D);

	free(xvals2D);
	free(yvals2D);

	if ( gridtype == GRID_LCC )
	  {		
	    double xi, xj;
	    double originLon, originLat, lonParY, lat1, lat2, xincm, yincm;
	    double zlat, zlon;
	    int projflag, scanflag;
	    int status;
	    int index;

	    gridInqLCC(gridID1, &originLon, &originLat, &lonParY, &lat1, &lat2, &xincm, &yincm,
		       &projflag, &scanflag);

	    xbounds2D = (double *) malloc(4*gridsize*sizeof(double));
	    ybounds2D = (double *) malloc(4*gridsize*sizeof(double));

	    for ( j = 0; j < ny; j++ )
	      for ( i = 0; i < nx; i++ )
		{
		  index = j*4*nx + 4*i;

		  xi = i+1.5;
		  xj = j+1.5;
		  status = W3FB12(xi, xj, originLat, originLon, xincm, lonParY, lat1, &zlat,&zlon);
		  xbounds2D[index+0] = zlon;
		  ybounds2D[index+0] = zlat;

		  xi = i+0.5;
		  xj = j+1.5;
		  status = W3FB12(xi, xj, originLat, originLon, xincm, lonParY, lat1, &zlat,&zlon);
		  xbounds2D[index+1] = zlon;
		  ybounds2D[index+1] = zlat;

		  xi = i+0.5;
		  xj = j+0.5;
		  status = W3FB12(xi, xj, originLat, originLon, xincm, lonParY, lat1, &zlat,&zlon);
		  xbounds2D[index+2] = zlon;
		  ybounds2D[index+2] = zlat;

		  xi = i+1.5;
		  xj = j+0.5;
		  status = W3FB12(xi, xj, originLat, originLon, xincm, lonParY, lat1, &zlat,&zlon);
		  xbounds2D[index+3] = zlon;
		  ybounds2D[index+3] = zlat;
		}

	    gridDefXbounds(gridID2, xbounds2D);
	    gridDefYbounds(gridID2, ybounds2D);

	    free(xbounds2D);
	    free(ybounds2D);
	  }
	else
	  {
	    if ( gridInqXbounds(gridID1, NULL) )
	      {
		xbounds = (double *) malloc(2*nx*sizeof(double));
		gridInqXbounds(gridID1, xbounds);
	      }
	    else if ( nx > 1 )
	      {
		xbounds = (double *) malloc(2*nx*sizeof(double));
		gridGenXbounds(nx, xvals, xbounds);
	      }

	    if ( gridInqYbounds(gridID1, NULL) )
	      {
		ybounds = (double *) malloc(2*ny*sizeof(double));
		gridInqYbounds(gridID1, ybounds);
	      }
	    else if ( ny > 1 )
	      {
		ybounds = (double *) malloc(2*ny*sizeof(double));
		if ( gridtype == GRID_SINUSOIDAL || gridtype == GRID_LAEA )
		  gridGenYboundsM(ny, yvals, ybounds);
		else
		  gridGenYbounds(ny, yvals, ybounds);
	      }

	    free(xvals);
	    free(yvals);

	    if ( xbounds && ybounds )
	      {
		xbounds2D = (double *) malloc(4*gridsize*sizeof(double));
		ybounds2D = (double *) malloc(4*gridsize*sizeof(double));

		if ( gridIsRotated(gridID1) )
		  {
		    gridGenRotBounds(gridID1, nx, ny, xbounds, ybounds, xbounds2D, ybounds2D);
		  }
		else
		  {
		    if ( gridtype == GRID_SINUSOIDAL )
		      {
#if defined (HAVE_LIBPROJ)
			int index;
			int i;
			PJ   *libProj;
			char *params[20];
			int nbpar=0;
			projUV data, res;

			nbpar = 0;
			params[nbpar++] = (char*) "proj=sinu";
			params[nbpar++] = (char*) "ellps=WGS84";

			if ( cdoVerbose )
			  for ( i = 0; i < nbpar; ++i )
			    cdoPrint("Proj.param[%d] = %s", i+1, params[i]);

			libProj = pj_init(nbpar, params);
			if ( !libProj )
			  cdoAbort("proj error: %s", pj_strerrno(pj_errno));

			/* libProj->over = 1; */		// allow longitude > 180�
			/*
			for ( j = 0; j < 2*ny; j++ ) printf("ybounds %d %g\n", j, ybounds[j]);
			for ( i = 0; i < 2*nx; i++ ) printf("xbounds %d %g\n", i, xbounds[i]);
			*/
			for ( j = 0; j < ny; j++ )
			  for ( i = 0; i < nx; i++ )
			    {
			      index = j*4*nx + 4*i;

			      data.u = xbounds[2*i];
			      data.v = ybounds[2*j];
			      res = pj_inv(data, libProj);
			      xbounds2D[index+0] = res.u*rad2deg;
			      ybounds2D[index+0] = res.v*rad2deg;

			      data.u = xbounds[2*i];
			      data.v = ybounds[2*j+1];
			      res = pj_inv(data, libProj);
			      xbounds2D[index+1] = res.u*rad2deg;
			      ybounds2D[index+1] = res.v*rad2deg;

			      data.u = xbounds[2*i+1];
			      data.v = ybounds[2*j+1];
			      res = pj_inv(data, libProj);
			      xbounds2D[index+2] = res.u*rad2deg;
			      ybounds2D[index+2] = res.v*rad2deg;

			      data.u = xbounds[2*i+1];
			      data.v = ybounds[2*j];
			      res = pj_inv(data, libProj);
			      xbounds2D[index+3] = res.u*rad2deg;
			      ybounds2D[index+3] = res.v*rad2deg;
			    }
#else
			cdoAbort("proj4 support not compiled in!");
#endif
		      }
		    else if ( gridtype == GRID_LAEA )
		      {
#if defined (HAVE_LIBPROJ)
			int index;
			int i;
			PJ   *libProj;
			char *params[20];
			int nbpar=0;
			projUV data, res;
			double a, lon_0, lat_0;

			gridInqLaea(gridID1, &a , &lon_0, &lat_0);

			nbpar = 0;
			params[nbpar++] = gen_param("proj=laea");
			params[nbpar++] = gen_param("a=%g", a);
			params[nbpar++] = gen_param("lon_0=%g", lon_0);
			params[nbpar++] = gen_param("lat_0=%g", lat_0);

			if ( cdoVerbose )
			  for ( i = 0; i < nbpar; ++i )
			    cdoPrint("Proj.param[%d] = %s", i+1, params[i]);

			libProj = pj_init(nbpar, params);
			if ( !libProj )
			  cdoAbort("proj error: %s", pj_strerrno(pj_errno));

			for ( i = 0; i < nbpar; ++i ) free(params[i]);

			/* libProj->over = 1; */		// allow longitude > 180�
			/*
			for ( j = 0; j < 2*ny; j++ ) printf("ybounds %d %g\n", j, ybounds[j]);
			for ( i = 0; i < 2*nx; i++ ) printf("xbounds %d %g\n", i, xbounds[i]);
			*/
			for ( j = 0; j < ny; j++ )
			  for ( i = 0; i < nx; i++ )
			    {
			      index = j*4*nx + 4*i;

			      data.u = xbounds[2*i];
			      data.v = ybounds[2*j];
			      res = pj_inv(data, libProj);
			      xbounds2D[index+0] = res.u*rad2deg;
			      ybounds2D[index+0] = res.v*rad2deg;

			      data.u = xbounds[2*i];
			      data.v = ybounds[2*j+1];
			      res = pj_inv(data, libProj);
			      xbounds2D[index+1] = res.u*rad2deg;
			      ybounds2D[index+1] = res.v*rad2deg;

			      data.u = xbounds[2*i+1];
			      data.v = ybounds[2*j+1];
			      res = pj_inv(data, libProj);
			      xbounds2D[index+2] = res.u*rad2deg;
			      ybounds2D[index+2] = res.v*rad2deg;

			      data.u = xbounds[2*i+1];
			      data.v = ybounds[2*j];
			      res = pj_inv(data, libProj);
			      xbounds2D[index+3] = res.u*rad2deg;
			      ybounds2D[index+3] = res.v*rad2deg;
			    }
#else
			cdoAbort("proj4 support not compiled in!");
#endif
		      }
		    else
		      {
			gridGenXbounds2D(nx, ny, xbounds, xbounds2D);
			gridGenYbounds2D(nx, ny, ybounds, ybounds2D);
		      }
		  }
		
		gridDefXbounds(gridID2, xbounds2D);
		gridDefYbounds(gridID2, ybounds2D);
		
		free(xbounds);
		free(ybounds);
		free(xbounds2D);
		free(ybounds2D);
	      }
	  }

	break;
      }
    default:
      {
	Error(func, "Gridtype %s unsupported!", gridNamePtr(gridtype));
	break;
      }
    }

  return (gridID2);
}


int gridToCell(int gridID1)
{
  static char func[] = "gridToCell";
  int gridID2;
  int gridtype, gridsize;


  gridtype = gridInqType(gridID1);
  gridsize = gridInqSize(gridID1);
  gridID2  = gridCreate(GRID_CELL, gridsize);
  gridDefPrec(gridID2, DATATYPE_FLT32);
	  
  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
      {
	int i, j, nx, ny;
	double *xvals, *yvals;
	double *xvals2D, *yvals2D;
	double *xbounds = NULL, *ybounds = NULL;
	double *xbounds2D, *ybounds2D;

	gridDefXname(gridID2, "lon");
	gridDefYname(gridID2, "lat");
	gridDefXlongname(gridID2, "longitude");
	gridDefYlongname(gridID2, "latitude");
	gridDefXunits(gridID2, "degrees_east");
	gridDefYunits(gridID2, "degrees_north");

	gridDefNvertex(gridID2, 4);

	nx = gridInqXsize(gridID1);
	ny = gridInqYsize(gridID1);
	 
	gridDefXsize(gridID2, gridsize);
	gridDefYsize(gridID2, gridsize);

	xvals = (double *) malloc(nx*sizeof(double));
	yvals = (double *) malloc(ny*sizeof(double));

	xvals2D = (double *) malloc(gridsize*sizeof(double));
	yvals2D = (double *) malloc(gridsize*sizeof(double));

	gridInqXvals(gridID1, xvals);
	gridInqYvals(gridID1, yvals);

	for ( j = 0; j < ny; j++ )
	  for ( i = 0; i < nx; i++ )
	    {
	      xvals2D[j*nx+i] = xvals[i];
	      yvals2D[j*nx+i] = yvals[j];
	    }

	gridDefXvals(gridID2, xvals2D);
	gridDefYvals(gridID2, yvals2D);

	for ( j = 0; j < ny; j++ )
	  for ( i = 0; i < nx; i++ )
	    {
	      xvals2D[j*nx+i] = xvals[i];
	      yvals2D[j*nx+i] = yvals[j];
	    }

	gridDefXvals(gridID2, xvals2D);
	gridDefYvals(gridID2, yvals2D);

	free(xvals2D);
	free(yvals2D);

	if ( gridInqXbounds(gridID1, NULL) )
	  {
	    xbounds = (double *) malloc(2*nx*sizeof(double));
	    gridInqXbounds(gridID1, xbounds);
	  }
	else if ( nx > 1 )
	  {
	    xbounds = (double *) malloc(2*nx*sizeof(double));
	    gridGenXbounds(nx, xvals, xbounds);
	  }

	if ( gridInqYbounds(gridID1, NULL) )
	  {
	    ybounds = (double *) malloc(2*ny*sizeof(double));
	    gridInqYbounds(gridID1, ybounds);
	  }
	else if ( ny > 1 )
	  {
	    ybounds = (double *) malloc(2*ny*sizeof(double));
	    gridGenYbounds(ny, yvals, ybounds);
	  }

	free(xvals);
	free(yvals);

	if ( xbounds )
	  {
	    xbounds2D = (double *) malloc(4*gridsize*sizeof(double));
	    gridGenXbounds2D(nx, ny, xbounds, xbounds2D);
	    gridDefXbounds(gridID2, xbounds2D);

	    free(xbounds);
	    free(xbounds2D);
	  }

	if ( ybounds )
	  {
	    ybounds2D = (double *) malloc(4*gridsize*sizeof(double));
	    gridGenYbounds2D(nx, ny, ybounds, ybounds2D);
	    gridDefYbounds(gridID2, ybounds2D);

	    free(ybounds);
	    free(ybounds2D);
	  }

	break;
      }
    case GRID_GME:
      {
	int nd, ni, ni2, ni3, i, j;
	int nv = 6;
	int *imask;
	double *xvals, *yvals;
	double *xbounds, *ybounds;

	nd  = gridInqGMEnd(gridID1);
	ni  = gridInqGMEni(gridID1);
	ni2 = gridInqGMEni2(gridID1);
	ni3 = gridInqGMEni3(gridID1);

	imask   = (int *) malloc(gridsize*sizeof(int));
	xvals   = (double *) malloc(gridsize*sizeof(double));
	yvals   = (double *) malloc(gridsize*sizeof(double));
	xbounds = (double *) malloc(nv*gridsize*sizeof(double));
	ybounds = (double *) malloc(nv*gridsize*sizeof(double));

	gme_grid(gridsize, xvals, yvals, xbounds, ybounds, imask, ni, nd, ni2, ni3);
	
	for ( i = 0; i < gridsize; i++ )
	  {
	    xvals[i] *= RAD2DEG;
	    yvals[i] *= RAD2DEG;

	    for ( j = 0; j < nv; j++ )
	      {
		xbounds[i*nv + j] *= RAD2DEG;
		ybounds[i*nv + j] *= RAD2DEG;
	      }
	    /* printf("%d %g %g\n", i, xvals[i], yvals[i]); */
	  }
	
	gridDefXsize(gridID2, gridsize);
	gridDefYsize(gridID2, gridsize);

	gridDefXvals(gridID2, xvals);
	gridDefYvals(gridID2, yvals);

	gridDefMask(gridID2, imask);

	gridDefNvertex(gridID2, nv);

	gridDefXbounds(gridID2, xbounds);
	gridDefYbounds(gridID2, ybounds);

	gridDefXunits(gridID2, "degrees");
	gridDefYunits(gridID2, "degrees");

	free (imask);
	free (xvals);
	free (yvals);
	free (xbounds);
	free (ybounds);
	
	break;
      }
    default:
      {
	Error(func, "Gridtype %s unsupported!", gridNamePtr(gridtype));
	break;
      }
    }

  return (gridID2);
}


int gridGenArea(int gridID, double *area)
{
  static char func[] = "gridGenArea";
  int status = 0;
  int i, k;
  int gridtype;
  int nv, gridsize;
  int lgrid_gen_bounds = FALSE;
  double total_area;
  double *grid_center_lon = NULL;
  double *grid_center_lat = NULL;
  double *grid_corner_lon = NULL;
  double *grid_corner_lat = NULL;
  int *grid_mask = NULL;
  struct geo p1, p2, p3;
  struct cart c1, c2, c3;

  gridsize = gridInqSize(gridID);
  gridtype = gridInqType(gridID);

  if ( gridtype != GRID_LONLAT      &&
       gridtype != GRID_GAUSSIAN    &&
       gridtype != GRID_LCC         &&
       gridtype != GRID_GME         &&
       gridtype != GRID_CURVILINEAR &&
       gridtype != GRID_CELL )
    {
      cdoAbort("Internal error! Unsupported gridtype: %s", gridNamePtr(gridtype)); 
    }

  if ( gridtype != GRID_CELL && gridtype != GRID_CURVILINEAR )
    {
      if ( gridtype == GRID_GME )
	{
	  gridID = gridToCell(gridID);
	  grid_mask = (int *) malloc(gridsize*sizeof(int));
	  gridInqMask(gridID, grid_mask);
	}
      else
	{
	  gridID = gridToCurvilinear(gridID);
	  lgrid_gen_bounds = TRUE;
	}
    }

  gridtype = gridInqType(gridID);

  if ( gridtype == GRID_CELL )
    nv = gridInqNvertex(gridID);
  else
    nv = 4;
  
  grid_center_lon = (double *) malloc(gridsize*sizeof(double));
  grid_center_lat = (double *) malloc(gridsize*sizeof(double));

  gridInqXvals(gridID, grid_center_lon);
  gridInqYvals(gridID, grid_center_lat);

  grid_corner_lon = (double *) malloc(nv*gridsize*sizeof(double));
  grid_corner_lat = (double *) malloc(nv*gridsize*sizeof(double));

  if ( gridInqYbounds(gridID, NULL) && gridInqXbounds(gridID, NULL) )
    {
      gridInqXbounds(gridID, grid_corner_lon);
      gridInqYbounds(gridID, grid_corner_lat);
    }
  else
    {
      if ( lgrid_gen_bounds )
	{
	  int nlon = gridInqXsize(gridID);
	  int nlat = gridInqYsize(gridID);
	  genXbounds(nlon, nlat, grid_center_lon, grid_corner_lon);
	  genYbounds(nlon, nlat, grid_center_lat, grid_corner_lat);
	}
      else
	{
	  status = 1;
	  return (status);
	}
    }
  
  total_area = 0;
  for ( i = 0; i < gridsize; ++i )
    {
      area[i] = 0;
      
      p3.lon = grid_center_lon[i]*deg2rad; 
      p3.lat = grid_center_lat[i]*deg2rad;
      c3 = gc2cc(&p3);
      
      for ( k = 1; k < nv; ++k )
	{
	  p1.lon = grid_corner_lon[i*nv+k-1]*deg2rad; 
	  p1.lat = grid_corner_lat[i*nv+k-1]*deg2rad;
	  c1 = gc2cc(&p1);
	  p2.lon = grid_corner_lon[i*nv+k]*deg2rad; 
	  p2.lat = grid_corner_lat[i*nv+k]*deg2rad;
	  c2 = gc2cc(&p2);

	  if ( (fabs(p1.lon*rad2deg - p2.lon*rad2deg) > 179) ||
	       (fabs(p2.lon*rad2deg - p3.lon*rad2deg) > 179) ||
	       (fabs(p3.lon*rad2deg - p1.lon*rad2deg) > 179) ) return(2);

	  area[i] += areas(&c1, &c2, &c3);
	}

      p1.lon = grid_corner_lon[i*nv+0]*deg2rad; 
      p1.lat = grid_corner_lat[i*nv+0]*deg2rad;
      c1 = gc2cc(&p1);
      p2.lon = grid_corner_lon[i*nv+nv-1]*deg2rad; 
      p2.lat = grid_corner_lat[i*nv+nv-1]*deg2rad;
      c2 = gc2cc(&p2);

      if ( (fabs(p1.lon*rad2deg - p2.lon*rad2deg) > 179) ||
	   (fabs(p2.lon*rad2deg - p3.lon*rad2deg) > 179) ||
	   (fabs(p3.lon*rad2deg - p1.lon*rad2deg) > 179) ) return(2);

      area[i] += areas(&c1, &c2, &c3);

      total_area += area[i];
    }

  if ( cdoVerbose ) cdoPrint("Total area = %g\n", total_area);

  free(grid_center_lon);
  free(grid_center_lat);
  free(grid_corner_lon);
  free(grid_corner_lat);
  if ( grid_mask ) free(grid_mask);

  return (status);
}


int gridGenWeights(int gridID, double *grid_area, double *grid_wgts)
{
  static char func[] = "gridGenWeights";
  int i, nvals, gridsize, gridtype;
  int status = 0;
  int *grid_mask = NULL;
  double total_area;

  gridtype = gridInqType(gridID);
  gridsize = gridInqSize(gridID);
  
  if ( gridtype == GRID_GME )
    {
      gridID = gridToCell(gridID);	  
      grid_mask = (int *) malloc(gridsize*sizeof(int));
      gridInqMask(gridID, grid_mask);
    }

  total_area = 0;
  nvals = 0;
  for ( i = 0; i < gridsize; i++ )
    {
      if ( grid_mask )
	if ( grid_mask[i] == 0 ) continue;
      total_area += grid_area[i];
      nvals++;
    }

  if ( cdoVerbose ) cdoPrint("Total area = %g\n", total_area);

  for ( i = 0; i < gridsize; i++ )
    {
      if ( grid_mask )
	if ( grid_mask[i] == 0 )
	  {
	    grid_wgts[i] = 0;
	    continue;
	  }
      
      grid_wgts[i] = grid_area[i] / total_area;
    }
  
  if ( grid_mask ) free(grid_mask);

  return (status);
}


static
int gridWeightsOld(int gridID, double *weights)
{
  static char func[] = "gridWeightsOld";
  int status = FALSE;
  int i, j, len;

  len = gridInqSize(gridID);

  if ( gridHasArea(gridID) )
    {
      gridInqArea(gridID, weights);
    }
  else
    {
      int gridtype = gridInqType(gridID);

      if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN )
	{
	  int     nlat, nlon;
	  int     datapoint;
	  double *lats = NULL, *lons = NULL;
	  double sumw;
	  double phi1, phi2, theta1, theta2, sindphi;

	  nlon = gridInqXsize(gridID);
	  nlat = gridInqYsize(gridID);

	  lons = 1 + (double *) malloc((nlon+2)*sizeof(double));
	  lats = 1 + (double *) malloc((nlat+2)*sizeof(double));

	  gridInqXvals(gridID, lons);
	  gridInqYvals(gridID, lats);

	  /* Interpolate to find latitudes outside boundaries. */
	  lats[-1]   = 2*lats[0] - lats[1];
	  lats[nlat] = 2*lats[nlat-1] - lats[nlat-2];
	  lons[-1]   = 2*lons[0] - lons[1];
	  lons[nlon] = 2*lons[nlon-1] - lons[nlon-2];
  
	  /*  Calculate weights.  */
	  /*  phi 1 and 2 and theta 1 and 2 represent respectively the boundary */
	  /*  latitudes and longitudes of a particular grid square.             */
	  datapoint = 0;
	  sumw = 0;
	  for ( j = 0; j < nlat; j++ )
	    {
	      phi1 = (lats[j-1]+lats[j])/2*deg2rad;
	      phi2 = (lats[j+1]+lats[j])/2*deg2rad;
	      if ( phi1 < (-1*M_PI/2) ) phi1 = -1*M_PI/2;
	      if ( phi1 > (   M_PI/2) ) phi1 =    M_PI/2;
	      if ( phi2 > (   M_PI/2) ) phi2 =    M_PI/2;
	      if ( phi2 < (-1*M_PI/2) ) phi2 = -1*M_PI/2;
	      sindphi = sin(phi2)-sin(phi1);
	      for( i = 0; i < nlon; i++ )
		{
		  if ( lons[i] >= lons[0]+360 || fabs(lats[j]) > 90 )
		    weights[datapoint] = 0;
		  else
		    {
		      theta1 = (lons[i-1]+lons[i])/2*deg2rad;
		      theta2 = (lons[i+1]+lons[i])/2*deg2rad;
		      weights[datapoint] = fabs((theta2-theta1)*sindphi);
		      sumw += weights[datapoint];
		    }
		  datapoint++;
		}
	    }

	  /* Normalise weights.  */
	  if( !DBL_IS_EQUAL(sumw, 0) )
	    for( i = 0; i < datapoint; i++ ) weights[i] /= sumw;

	  if ( lons-1 ) free(lons-1);
	  if ( lats-1 ) free(lats-1);
	}
      else
	{
	  status = TRUE;

	  for ( i = 0; i < len; i++ ) weights[i] = 1./len;
	}
    }

  return (status);
}


int gridWeights(int gridID, double *grid_wgts)
{
  static char func[] = "gridWeights";
  int i, gridsize, gridtype;
  int a_status, w_status;
  double *grid_area;

  gridtype = gridInqType(gridID);
  gridsize = gridInqSize(gridID);
  
  grid_area = (double *) malloc(gridsize*sizeof(double));

  a_status = 0;

  if ( gridHasArea(gridID) )
    {
      if ( cdoVerbose ) cdoPrint("Using existing grid cell area!");
      gridInqArea(gridID, grid_area);
    }
  else
    {
      if ( gridtype != GRID_LONLAT      &&
	   gridtype != GRID_GAUSSIAN    &&
	   gridtype != GRID_LCC         &&
	   gridtype != GRID_GME         &&
	   gridtype != GRID_CURVILINEAR &&
	   gridtype != GRID_CELL )
	{
	  a_status = 1;
	}
      else
	{
	  a_status = gridGenArea(gridID, grid_area);
	}
    }

  if ( a_status == 0 )
    {
      w_status = gridGenWeights(gridID, grid_area, grid_wgts);
    }
  else
    {
      for ( i = 0; i < gridsize; ++i )
	grid_wgts[i] = 1./gridsize;

      w_status = 1;
    }

  free(grid_area);

  return (w_status);
}

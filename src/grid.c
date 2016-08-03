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

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#if defined(_OPENMP)
#include <omp.h>
#endif

#include <stdio.h>
#include <stdarg.h> /* va_list */

#if defined(HAVE_LIBPROJ)
#include "proj_api.h"
#endif

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "error.h"
#include "grid.h"


static
void scale_vec(double scalefactor, long nvals, double *restrict values)
{
#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(nvals, scalefactor, values)
#endif
  for ( long n = 0; n < nvals; ++n )
    {
      values[n] *= scalefactor;
    }
}


void grid_to_radian(const char *units, long nvals, double *restrict values, const char *description)
{
  if ( cmpstr(units, "degree") == 0 )
    {
      scale_vec(DEG2RAD, nvals, values);
    }
  else if ( cmpstr(units, "radian") == 0 )
    {
      /* No conversion necessary */
    }
  else
    {
      cdoWarning("Unknown units [%s] supplied for %s; proceeding assuming radians!", units, description);
    }
}


void grid_to_degree(const char *units, long nvals, double *restrict values, const char *description)
{
  if ( cmpstr(units, "radian") == 0 )
    {
      for ( long n = 0; n < nvals; ++n ) values[n] *= RAD2DEG;
    }
  else if ( cmpstr(units, "degree") == 0 )
    {
      /* No conversion necessary */
    }
  else
    {
      cdoWarning("Unknown units [%s] supplied for %s; proceeding assuming degress!", units, description);
    }
}


int gridToZonal(int gridID1)
{
  int gridtype = gridInqType(gridID1);
  int gridsize = gridInqYsize(gridID1);
  int gridID2  = gridCreate(gridtype, gridsize);
	  
  if ( gridtype == GRID_LONLAT   ||
       gridtype == GRID_GAUSSIAN ||
       gridtype == GRID_GENERIC )
    {
      gridDefXsize(gridID2, 1);
      gridDefYsize(gridID2, gridsize);

      double xval = 0;
      gridDefXvals(gridID2, &xval);

      if ( gridInqYvals(gridID1, NULL) )
	{
	  double *yvals = (double*) Malloc(gridsize*sizeof(double));
	  gridInqYvals(gridID1, yvals);
	  gridDefYvals(gridID2, yvals);
	  Free(yvals);
	}
    }
  else
    {
      Error("Gridtype %s unsupported!", gridNamePtr(gridtype));
    }

  return gridID2;
}


int gridToMeridional(int gridID1)
{
  int gridtype = gridInqType(gridID1);
  int gridsize = gridInqXsize(gridID1);
  int gridID2  = gridCreate(gridtype, gridsize);
	  
  if ( gridtype == GRID_LONLAT   ||
       gridtype == GRID_GAUSSIAN ||
       gridtype == GRID_GENERIC )
    {
      gridDefXsize(gridID2, gridsize);
      gridDefYsize(gridID2, 1);

      if ( gridInqXvals(gridID1, NULL) )
	{
	  double *xvals = (double*) Malloc(gridsize*sizeof(double));
	  gridInqXvals(gridID1, xvals);
	  gridDefXvals(gridID2, xvals);
	  Free(xvals);
	}

      double yval = 0;
      gridDefYvals(gridID2, &yval);
    }
  else
    {
      Error("Gridtype %s unsupported!", gridNamePtr(gridtype));
    }

  return gridID2;
}


void grid_gen_corners(int n, const double* restrict vals, double* restrict corners)
{
  if ( n == 1 )
    {
      corners[0] = vals[0];
      corners[1] = vals[0];
    }
  else
    {
      for ( int i = 0; i < n-1; ++i )
	corners[i+1] = 0.5*(vals[i] + vals[i+1]);

      corners[0] = 2*vals[0] - corners[1];
      corners[n] = 2*vals[n-1] - corners[n-1];
    }
}


void grid_gen_bounds(int n, const double *restrict vals, double *restrict bounds)
{
  for ( int i = 0; i < n-1; ++i )
    {
      bounds[2*i+1]   = 0.5*(vals[i] + vals[i+1]);
      bounds[2*(i+1)] = 0.5*(vals[i] + vals[i+1]);
    }

  bounds[0]     = 2*vals[0] - bounds[1];
  bounds[2*n-1] = 2*vals[n-1] - bounds[2*(n-1)];
}


void grid_check_lat_borders(int n, double *ybounds)
{
  if ( ybounds[0] > ybounds[n-1] )
    {
      if ( ybounds[0]   >  88 ) ybounds[0]   =  90;
      if ( ybounds[n-1] < -88 ) ybounds[n-1] = -90;
    }
  else
    {
      if ( ybounds[0]   < -88 ) ybounds[0]   = -90;
      if ( ybounds[n-1] >  88 ) ybounds[n-1] =  90;
    }
}


void grid_cell_center_to_bounds_X2D(const char* xunitstr, long xsize, long ysize, const double* restrict grid_center_lon, 
				    double* restrict grid_corner_lon, double dlon)
{
  if ( ! (dlon > 0) ) dlon = 360./xsize;
  /*
  if ( xsize == 1 || (grid_center_lon[xsize-1]-grid_center_lon[0]+dlon) < 359 )
    cdoAbort("Cannot calculate Xbounds for %d vals with dlon = %g", xsize, dlon);
  */
  for ( int i = 0; i < xsize; ++i )
    {
      double minlon = grid_center_lon[i] - 0.5*dlon;
      double maxlon = grid_center_lon[i] + 0.5*dlon;
      for ( int j = 0; j < ysize; ++j )
	{
	  int index = (j<<2)*xsize + (i<<2);
	  grid_corner_lon[index  ] = minlon;
	  grid_corner_lon[index+1] = maxlon;
	  grid_corner_lon[index+2] = maxlon;
	  grid_corner_lon[index+3] = minlon;
	}
    }
}

static
double genYmin(double y1, double y2)
{
  double dy = y2 - y1;
  double ymin = y1 - dy/2;

  if ( y1 < -85 && ymin < -87.5 ) ymin = -90;

  if ( cdoVerbose )
    cdoPrint("genYmin: y1 = %g  y2 = %g  dy = %g  ymin = %g", y1, y2, dy, ymin);

  return ymin;
}

static
double genYmax(double y1, double y2)
{
  double dy = y1 - y2;
  double ymax = y1 + dy/2;

  if ( y1 > 85 && ymax > 87.5 ) ymax = 90;

  if ( cdoVerbose )
    cdoPrint("genYmax: y1 = %g  y2 = %g  dy = %g  ymax = %g", y1, y2, dy, ymax);

  return ymax;
}



/*****************************************************************************/

void grid_cell_center_to_bounds_Y2D(const char* yunitstr, long xsize, long ysize, const double* restrict grid_center_lat, double* restrict grid_corner_lat)
{
  double minlat, maxlat;

  double firstlat = grid_center_lat[0];
  double lastlat  = grid_center_lat[xsize*ysize-1];

  // if ( ysize == 1 ) cdoAbort("Cannot calculate Ybounds for 1 value!");

  for ( int j = 0; j < ysize; ++j )
    {
      if ( ysize == 1 )
	{
	  minlat = grid_center_lat[0] - 360./ysize;
	  maxlat = grid_center_lat[0] + 360./ysize;
	}
      else
	{
	  int index = j*xsize;
	  if ( firstlat > lastlat )
	    {
	      if ( j == 0 )
		maxlat = genYmax(grid_center_lat[index], grid_center_lat[index+xsize]);
	      else
		maxlat = 0.5*(grid_center_lat[index]+grid_center_lat[index-xsize]);

	      if ( j == (ysize-1) )
		minlat = genYmin(grid_center_lat[index], grid_center_lat[index-xsize]);
	      else
		minlat = 0.5*(grid_center_lat[index]+grid_center_lat[index+xsize]);
	    }
	  else
	    {
	      if ( j == 0 )
		minlat = genYmin(grid_center_lat[index], grid_center_lat[index+xsize]);
	      else
		minlat = 0.5*(grid_center_lat[index]+grid_center_lat[index-xsize]);

	      if ( j == (ysize-1) )
		maxlat = genYmax(grid_center_lat[index], grid_center_lat[index-xsize]);
	      else
		maxlat = 0.5*(grid_center_lat[index]+grid_center_lat[index+xsize]);
	    }
	}

      for ( int i = 0; i < xsize; ++i )
	{
	  int index = (j<<2)*xsize + (i<<2);
	  grid_corner_lat[index  ] = minlat;
	  grid_corner_lat[index+1] = minlat;
	  grid_corner_lat[index+2] = maxlat;
	  grid_corner_lat[index+3] = maxlat;
	}
    }
}

static
void gridGenRotBounds(double xpole, double ypole, double angle, int nx, int ny,
		      double *xbounds, double *ybounds, double *xbounds2D, double *ybounds2D)
{
  double minlon, maxlon;
  double minlat, maxlat;

  for ( int j = 0; j < ny; j++ )
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

      for ( int i = 0; i < nx; i++ )
	{
	  minlon = xbounds[2*i];
	  maxlon = xbounds[2*i+1];

	  int index = j*4*nx + 4*i;
	  xbounds2D[index+0] = lamrot_to_lam(minlat, minlon, ypole, xpole, angle);
	  xbounds2D[index+1] = lamrot_to_lam(minlat, maxlon, ypole, xpole, angle);
	  xbounds2D[index+2] = lamrot_to_lam(maxlat, maxlon, ypole, xpole, angle);
	  xbounds2D[index+3] = lamrot_to_lam(maxlat, minlon, ypole, xpole, angle);

	  ybounds2D[index+0] = phirot_to_phi(minlat, minlon, ypole, angle);
	  ybounds2D[index+1] = phirot_to_phi(minlat, maxlon, ypole, angle);
	  ybounds2D[index+2] = phirot_to_phi(maxlat, maxlon, ypole, angle);
	  ybounds2D[index+3] = phirot_to_phi(maxlat, minlon, ypole, angle);
	}
    }
}


void grid_gen_xbounds2D(int nx, int ny, const double *restrict xbounds, double *restrict xbounds2D)
{
#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(nx, ny, xbounds, xbounds2D)
#endif
  for ( int i = 0; i < nx; ++i )
    {
      double minlon = xbounds[2*i  ];
      double maxlon = xbounds[2*i+1];

      for ( int j = 0; j < ny; ++j )
	{
	  int index = j*4*nx + 4*i;
	  xbounds2D[index  ] = minlon;
	  xbounds2D[index+1] = maxlon;
	  xbounds2D[index+2] = maxlon;
	  xbounds2D[index+3] = minlon;
	}
    }
}


void grid_gen_ybounds2D(int nx, int ny, const double *restrict ybounds, double *restrict ybounds2D)
{
#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(nx, ny, ybounds, ybounds2D)
#endif
  for ( int j = 0; j < ny; ++j )
    {
      double minlat, maxlat;
      if ( ybounds[0] > ybounds[1] )
	{
	  maxlat = ybounds[2*j  ];
	  minlat = ybounds[2*j+1];
	}
      else
	{
	  maxlat = ybounds[2*j+1];
	  minlat = ybounds[2*j  ];
	}

      for ( int i = 0; i < nx; ++i )
	{
	  int index = j*4*nx + 4*i;
	  ybounds2D[index  ] = minlat;
	  ybounds2D[index+1] = minlat;
	  ybounds2D[index+2] = maxlat;
	  ybounds2D[index+3] = maxlat;
	}
    }
}

static
char *gen_param(const char *fmt, ...)
{
  va_list args;
  char str[256];

  va_start(args, fmt);

  int len = vsprintf(str, fmt, args);

  va_end(args);

  len++;
  char *rstr = (char*) Malloc(len*sizeof(char));
  memcpy(rstr, str, len*sizeof(char));

  return rstr;
}

static
void lcc_to_geo(int gridID, int gridsize, double *xvals, double *yvals)
{
  double originLon, originLat, lonParY, lat1, lat2, xincm, yincm;
  int projflag, scanflag;
  proj_info_t proj;

  gridInqLCC(gridID, &originLon, &originLat, &lonParY, &lat1, &lat2, &xincm, &yincm, &projflag, &scanflag);
  /*
    while ( originLon < 0 ) originLon += 360;
    while ( lonParY   < 0 ) lonParY   += 360;
  */
  if ( IS_NOT_EQUAL(xincm, yincm) )
    Warning("X and Y increment must be equal on Lambert Conformal grid (Xinc = %g, Yinc = %g)\n", 
	    xincm, yincm);
  /*
  if ( IS_NOT_EQUAL(lat1, lat2) )
    Warning("Lat1 and Lat2 must be equal on Lambert Conformal grid (Lat1 = %g, Lat2 = %g)\n", 
	    lat1, lat2);
  */
  map_set(PROJ_LC, originLat, originLon, xincm, lonParY, lat1, lat2, &proj);

  double zlat, zlon;
  for ( int i = 0; i < gridsize; i++ )
    {
      double xi = xvals[i];
      double xj = yvals[i];
      // status = W3FB12(xi, xj, originLat, originLon, xincm, lonParY, lat1, &zlat, &zlon);
      ijll_lc(xi, xj, proj, &zlat, &zlon);
      xvals[i] = zlon;
      yvals[i] = zlat;
    }
}

static
void sinusoidal_to_geo(int gridsize, double *xvals, double *yvals)
{
#if defined(HAVE_LIBPROJ)
  char *params[20];
  projUV data, res;

  int nbpar = 0;
  params[nbpar++] = gen_param("proj=sinu");
  params[nbpar++] = gen_param("ellps=WGS84");

  if ( cdoVerbose )
    for ( int i = 0; i < nbpar; ++i )
      cdoPrint("Proj.param[%d] = %s", i+1, params[i]);

  projPJ proj = pj_init(nbpar, params);
  if ( !proj )
    cdoAbort("proj error: %s", pj_strerrno(pj_errno));

  for ( int i = 0; i < nbpar; ++i ) Free(params[i]);

  /* proj->over = 1; */		/* allow longitude > 180 */

  for ( int i = 0; i < gridsize; i++ )
    {
      data.u = xvals[i];
      data.v = yvals[i];
      res = pj_inv(data, proj);
      xvals[i] = res.u*RAD2DEG;
      yvals[i] = res.v*RAD2DEG;
      if ( xvals[i] < -9000. || xvals[i] > 9000. ) xvals[i] = -9999.;
      if ( yvals[i] < -9000. || yvals[i] > 9000. ) yvals[i] = -9999.;
    }

  pj_free(proj);
#else
  cdoAbort("proj4 support not compiled in!");
#endif
}

static
void grid_get_param_laea(int gridID, double *a, double *lon_0, double *lat_0)
{
  *a = 0; *lon_0 = 0; *lat_0 = 0;

  int gridtype = gridInqType(gridID);
  if ( gridtype == GRID_LAEA )
    gridInqLaea(gridID, a , lon_0, lat_0);
  else
    {
      const char *projection = "lambert_azimuthal_equal_area";
      char mapping[CDI_MAX_NAME]; mapping[0] = 0;
      cdiGridInqKeyStr(gridID, CDI_KEY_MAPPING, CDI_MAX_NAME, mapping);
      if ( mapping[0] && strcmp(mapping, projection) == 0 )
        {
          int atttype, attlen;
          char attname[CDI_MAX_NAME+1];

          int natts;
          vlistInqNatts(gridID, CDI_GLOBAL, &natts);

          for ( int iatt = 0; iatt < natts; ++iatt )
            {
              vlistInqAtt(gridID, CDI_GLOBAL, iatt, attname, &atttype, &attlen);

              if ( attlen != 1 ) continue;

              if ( atttype == DATATYPE_FLT32 || atttype == DATATYPE_FLT64 )
                {
                  double attflt;
                  vlistInqAttFlt(gridID, CDI_GLOBAL, attname, attlen, &attflt);
                  if      ( strcmp(attname, "earth_radius") == 0 )                    *a     = attflt;
                  else if ( strcmp(attname, "longitude_of_projection_origin")  == 0 ) *lon_0 = attflt;
                  else if ( strcmp(attname, "latitude_of_projection_origin") == 0 )   *lat_0 = attflt;
                }
            }
        }
      else
        Warning("%s mapping parameter missing!", projection);
    }
}

static
void laea_to_geo(int gridID, int gridsize, double *xvals, double *yvals)
{
#if defined(HAVE_LIBPROJ)
  char *params[20];
  projUV data, res;
  
  double a, lon_0, lat_0;
  grid_get_param_laea(gridID, &a , &lon_0, &lat_0);

  int nbpar = 0;
  params[nbpar++] = gen_param("proj=laea");
  if ( a > 0 ) params[nbpar++] = gen_param("a=%g", a);
  params[nbpar++] = gen_param("lon_0=%g", lon_0);
  params[nbpar++] = gen_param("lat_0=%g", lat_0);

  if ( cdoVerbose )
    for ( int i = 0; i < nbpar; ++i )
      cdoPrint("Proj.param[%d] = %s", i+1, params[i]);

  projPJ proj = pj_init(nbpar, &params[0]);
  if ( !proj ) cdoAbort("proj error: %s", pj_strerrno(pj_errno));

  for ( int i = 0; i < nbpar; ++i ) Free(params[i]);

  /* proj->over = 1; */		/* allow longitude > 180 */

  for ( int i = 0; i < gridsize; i++ )
    {
      data.u = xvals[i];
      data.v = yvals[i];
      res = pj_inv(data, proj);
      xvals[i] = res.u*RAD2DEG;
      yvals[i] = res.v*RAD2DEG;
      if ( xvals[i] < -9000. || xvals[i] > 9000. ) xvals[i] = -9999.;
      if ( yvals[i] < -9000. || yvals[i] > 9000. ) yvals[i] = -9999.;
    }

  pj_free(proj);
#else
  cdoAbort("proj4 support not compiled in!");
#endif
}

static
void lcc2_to_geo(int gridID, int gridsize, double *xvals, double *yvals)
{
#if defined(HAVE_LIBPROJ)
  char *params[20];
  projUV data, res;

  double a, lon_0, lat_0, lat_1, lat_2;
  gridInqLcc2(gridID, &a , &lon_0, &lat_0, &lat_1, &lat_2);

  int nbpar = 0;
  params[nbpar++] = gen_param("proj=lcc");
  if ( a > 0 ) params[nbpar++] = gen_param("a=%g", a);
  params[nbpar++] = gen_param("lon_0=%g", lon_0);
  params[nbpar++] = gen_param("lat_0=%g", lat_0);
  params[nbpar++] = gen_param("lat_1=%g", lat_1);
  params[nbpar++] = gen_param("lat_2=%g", lat_2);

  if ( cdoVerbose )
    for ( int i = 0; i < nbpar; ++i )
      cdoPrint("Proj.param[%d] = %s", i+1, params[i]);
  
  projPJ proj = pj_init(nbpar, &params[0]);
  if ( !proj )
    cdoAbort("proj error: %s", pj_strerrno(pj_errno));

  for ( int i = 0; i < nbpar; ++i ) Free(params[i]);

  /* proj->over = 1; */		/* allow longitude > 180 */
  
  for ( int i = 0; i < gridsize; i++ )
    {
      data.u = xvals[i];
      data.v = yvals[i];
      res = pj_inv(data, proj);
      xvals[i] = res.u*RAD2DEG;
      yvals[i] = res.v*RAD2DEG;
    }

  pj_free(proj);
#else
  cdoAbort("proj4 support not compiled in!");
#endif
}

/*
 * grib_get_reduced_row: code from GRIB_API 1.10.4
 *
 * Description:
 *   computes the number of points within the range lon_first->lon_last and the zero based indexes
 *   ilon_first,ilon_last of the first and last point given the number of points along a parallel (pl)
 *
 */
static
void grib_get_reduced_row(long pl,double lon_first,double lon_last,long* npoints,long* ilon_first, long* ilon_last )
{
  double dlon_first=0, dlon_last=0;

  double range = lon_last-lon_first;
  if ( range < 0 ) {range+=360;lon_first-=360;}

  /* computing integer number of points and coordinates without using floating point resolution*/
  *npoints=(range*pl)/360.0+1;
  *ilon_first=(lon_first*pl)/360.0;
  *ilon_last=(lon_last*pl)/360.0;

  long irange=*ilon_last-*ilon_first+1;

  if (irange != *npoints) {
    if (irange > *npoints) {
      /* checking if the first point is out of range*/
      dlon_first=((*ilon_first)*360.0)/pl;
      if (dlon_first < lon_first) {(*ilon_first)++;irange--;
      }

      /* checking if the last point is out of range*/
      dlon_last=((*ilon_last)*360.0)/pl;
      if (dlon_last > lon_last) {(*ilon_last)--;irange--;
      }
    } else {
      int ok=0;
      /* checking if the point before the first is in the range*/
      dlon_first=((*ilon_first-1)*360.0)/pl;
      if (dlon_first > lon_first) {(*ilon_first)--;irange++;ok=1;
      }

      /* checking if the point after the last is in the range*/
      dlon_last=((*ilon_last+1)*360.0)/pl;
      if (dlon_last < lon_last) {(*ilon_last)++;irange++;ok=1;
      }

      /* if neither of the two are triggered then npoints is too large */
      if (!ok) {(*npoints)--;
      }
    }

    //   assert(*npoints==irange);
  } else {
	  /* checking if the first point is out of range*/
	  dlon_first=((*ilon_first)*360.0)/pl;
	  if (dlon_first < lon_first) {
		  (*ilon_first)++;(*ilon_last)++;
	  }
  }

  if (*ilon_first<0) *ilon_first+=pl;
}

static
int qu2reg_subarea(int gridsize, int np, double xfirst, double xlast, 
		   double *array, int *rowlon, int ny, double missval, int *iret, int lmiss, int lperio, int lveggy)
{
  /* sub area (longitudes) */
  long ilon_firstx;
  long ilon_first, ilon_last;
  int i, j;
  long row_count;
  int rlon, nwork = 0;
  int np4 = np*4;
  int size = 0;
  int wlen;
  int ii;

  if ( np <= 0 ) cdoAbort("Number of values between pole and equator missing!");

  grib_get_reduced_row(np4, xfirst, xlast, &row_count, &ilon_firstx, &ilon_last);
  int nx = row_count;
  // printf("nx %d  %ld %ld lon1 %g lon2 %g\n", nx, ilon_firstx, ilon_last, (ilon_firstx*360.)/np4, (ilon_last*360.)/np4);

  for ( j = 0; j < ny; ++j ) nwork += rowlon[j];

  double **pwork = (double **) Malloc(ny*sizeof(double *));
  double *work = (double*) Malloc(ny*np4*sizeof(double));
  wlen = 0;
  pwork[0] = work;
  for ( j = 1; j < ny; ++j )
    {
      wlen += rowlon[j-1];
      pwork[j] = work + wlen;
    } 
  // printf(" ny, np4, nwork %d %d %d wlen %d\n", ny, np4, nwork, wlen);

  for ( j = 0; j < ny; ++j )
    {
      rlon = rowlon[j];
      for ( i = 0; i < rlon; ++i ) pwork[j][i] = missval;
    } 

  double *parray = array;
  for ( j = 0; j < ny; ++j )
    {
      rlon = rowlon[j];
      row_count = 0;
      grib_get_reduced_row(rlon, xfirst, xlast, &row_count, &ilon_first, &ilon_last);
      // printf("j %d xfirst %g xlast %g rowlon %d %ld %ld %ld %g %g\n", j, xfirst, xlast, rlon, row_count, ilon_first, ilon_last, (ilon_first*360.)/rlon, (ilon_last*360.)/rlon);

      for ( i = ilon_first; i < (ilon_first+row_count); ++i )
	{
	  ii = i;
	  if ( ii >= rlon ) ii -= rlon; 
	  pwork[j][ii] = *parray;
	  parray++;
	}
      size += row_count;
    }

  if ( gridsize != size ) cdoAbort("gridsize1 inconsistent!");

  (void) qu2reg3_double(work, rowlon, ny, np4, missval, iret, lmiss, lperio, lveggy);

  wlen = 0;
  pwork[0] = work;
  for ( j = 1; j < ny; ++j )
    {
      wlen += np4;
      pwork[j] = work + wlen;
    } 

  // printf("nx, ilon_firstx %d %ld\n", nx, ilon_firstx);
  parray = array;
  for ( j = 0; j < ny; ++j )
    {
      for ( i = ilon_firstx; i < (ilon_firstx+nx); ++i )
	{
	  ii = i;
	  if ( ii >= np4 ) ii -= np4; 
	  *parray = pwork[j][ii];
	  parray++;
	}
    }

  Free(work);
  Free(pwork);

  return nx;
}


void field2regular(int gridID1, int gridID2, double missval, double *array, int nmiss, int lnearest)
{
  int gridtype = gridInqType(gridID1);
  if ( gridtype != GRID_GAUSSIAN_REDUCED ) Error("Not a reduced Gaussian grid!");

  int lmiss = nmiss > 0;
  int lperio = 1;

  int ny = gridInqYsize(gridID1);
  int np = gridInqNP(gridID1);

  int *rowlon = (int*) Malloc(ny*sizeof(int));
  gridInqRowlon(gridID1, rowlon);

  double xfirstandlast[2] = {0, 0};
  gridInqXvals(gridID1, xfirstandlast);
  double xfirst = xfirstandlast[0];
  double xlast  = xfirstandlast[1];

  int iret;
  int nx = 0;
  if ( fabs(xfirst) > 0 || (np > 0 && fabs(xlast - (360.0-90.0/np)) > 90.0/np) )
    {
      nx = qu2reg_subarea(gridInqSize(gridID1), np, xfirst, xlast, array, rowlon, ny, missval, &iret, lmiss, lperio, lnearest);
    }
  else
    {
      nx = 2*ny;
      (void) qu2reg3_double(array, rowlon, ny, nx, missval, &iret, lmiss, lperio, lnearest);
    }

  if ( gridInqSize(gridID2) != nx*ny ) Error("Gridsize differ!");

  Free(rowlon);
}


int gridToRegular(int gridID1)
{
  int nx = 0;
  double *xvals = NULL;

  int gridtype = gridInqType(gridID1);
  if ( gridtype != GRID_GAUSSIAN_REDUCED ) Error("Not a reduced Gaussian grid!");

  int ny = gridInqYsize(gridID1);
  int np = gridInqNP(gridID1);

  double *yvals = (double*) Malloc(ny*sizeof(double));
  gridInqYvals(gridID1, yvals);

  double xfirstandlast[2] = {0, 0};
  gridInqXvals(gridID1, xfirstandlast);
  double xfirst = xfirstandlast[0];
  double xlast  = xfirstandlast[1];

  if ( fabs(xfirst) > 0 || (np > 0 && fabs(xlast - (360.0-90.0/np)) > 90.0/np) )
    {
      /* sub area (longitudes) */
      long ilon_first, ilon_last;
      long row_count;
      int np4 = np*4;
      int *rowlon = NULL;

      if ( np <= 0 ) cdoAbort("Number of values between pole and equator missing!");

      grib_get_reduced_row(np4, xfirst, xlast, &row_count, &ilon_first, &ilon_last);

      nx = row_count;
      xvals = (double*) Malloc(nx*sizeof(double));
      for ( int i = 0; i < nx; ++i )
	{
	  xvals[i] = ((ilon_first+i)*360.)/np4;
	  if ( xfirst > xlast ) xvals[i] -= 360.;
	}

      Free(rowlon);
    }
  else
    {
      nx = 2*ny;
      xvals = (double*) Malloc(nx*sizeof(double));
      for ( int i = 0; i < nx; ++i ) xvals[i] = i * 360./nx;
    }

  int gridsize = nx*ny;
  int gridID2  = gridCreate(GRID_GAUSSIAN, gridsize);
	  
  gridDefXsize(gridID2, nx);
  gridDefYsize(gridID2, ny);
  
  gridDefXvals(gridID2, xvals);
  gridDefYvals(gridID2, yvals);
  gridDefNP(gridID2, np);

  Free(xvals);
  Free(yvals);

  return gridID2;
}

static
void gridCopyMask(int gridID1, int gridID2, long gridsize)
{
  if ( gridInqMask(gridID1, NULL) )
    {
      int *mask = (int*) Malloc(gridsize*sizeof(int));
      gridInqMask(gridID1, mask);
      gridDefMask(gridID2, mask);
      Free(mask);
    }
}

static
bool check_range(long n, double *vals, double valid_min, double valid_max)
{
  bool status = false;

  for ( long i = 0; i < n; ++i )
    {
      if ( vals[i] < valid_min || vals[i] > valid_max )
	{
	  status = true;
	  break;
	}
    }

  return status;
}


int gridToCurvilinear(int gridID1, int lbounds)
{
  long index;
  int gridtype = gridInqType(gridID1);

  size_t gridsize = (size_t) gridInqSize(gridID1);
  int gridID2 = gridCreate(GRID_CURVILINEAR, (int) gridsize);
  gridDefPrec(gridID2, DATATYPE_FLT32);

  bool lproj_rll   = false;
  bool lproj_laea  = false;
  bool lproj_lcc   = false;
  bool lproj_lsinu = false;
  if ( gridtype == GRID_PROJECTION )
    {
      int projtype = gridInqProjType(gridID1);
      if      ( projtype == CDI_PROJ_RLL  ) lproj_rll   = true;
      else if ( projtype == CDI_PROJ_LAEA ) lproj_laea  = true;
      else if ( projtype == CDI_PROJ_LCC  ) lproj_lcc   = true;
      else if ( projtype == CDI_PROJ_SINU ) lproj_lsinu = true;

      if ( lproj_rll || lproj_laea || lproj_lcc || lproj_lsinu ) gridtype = GRID_LONLAT;
    }

  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_LCC:
    case GRID_LCC2:
    case GRID_LAEA:
    case GRID_SINUSOIDAL:
      {
        double xpole = 0, ypole = 0, angle = 0;
	double xscale = 1, yscale = 1;
	double *xvals = NULL, *yvals = NULL;
	double *xbounds = NULL, *ybounds = NULL;
	char xunits[CDI_MAX_NAME], yunits[CDI_MAX_NAME];

	int nx = gridInqXsize(gridID1);
	int ny = gridInqYsize(gridID1);

	gridInqXunits(gridID1, xunits);
	gridInqYunits(gridID1, yunits);

	if ( gridtype == GRID_LAEA || lproj_laea )
	  {
	    bool lvalid_xunits = false;
	    bool lvalid_yunits = false;
	    int len;
	    len = (int) strlen(xunits);
	    if ( len == 1 && memcmp(xunits, "m",  1) == 0 ) lvalid_xunits = true;
	    if ( len == 2 && memcmp(xunits, "km", 2) == 0 ) lvalid_xunits = true;
	    len = (int) strlen(yunits);
	    if ( len == 1 && memcmp(yunits, "m",  1) == 0 ) lvalid_yunits = true;
	    if ( len == 2 && memcmp(yunits, "km", 2) == 0 ) lvalid_yunits = true;

	    if ( lvalid_xunits == false )
	      cdoWarning("Possibly wrong result! Invalid x-coordinate units: \"%s\" (expected \"m\" or \"km\")", xunits);
	    if ( lvalid_yunits == false )
	      cdoWarning("Possibly wrong result! Invalid y-coordinate units: \"%s\" (expected \"m\" or \"km\")", yunits);
	  }

	if ( memcmp(xunits, "km", 2) == 0 ) xscale = 1000;
	if ( memcmp(yunits, "km", 2) == 0 ) yscale = 1000;

	gridDefXsize(gridID2, nx);
	gridDefYsize(gridID2, ny);

	double *xvals2D = (double*) Malloc(gridsize*sizeof(double));
	double *yvals2D = (double*) Malloc(gridsize*sizeof(double));

        if ( gridtype == GRID_LCC )
	  {
	    for ( int j = 0; j < ny; j++ )
	      for ( int i = 0; i < nx; i++ )
		{
		  xvals2D[j*nx+i] = i+1;
		  yvals2D[j*nx+i] = j+1;
		}

	    lcc_to_geo(gridID1, gridsize, xvals2D, yvals2D);
	  }
	else
	  {
            /*
	    if ( ! (gridInqXvals(gridID1, NULL) && gridInqYvals(gridID1, NULL)) )
	      Error("Grid has no coordinates");
            */
            if ( nx == 0 ) nx = 1;
            if ( ny == 0 ) ny = 1;

            xvals = (double*) Malloc(nx*sizeof(double));
            yvals = (double*) Malloc(ny*sizeof(double));
            
            if ( gridInqXvals(gridID1, NULL) )
              gridInqXvals(gridID1, xvals);
            else
              for ( int i = 0; i < nx; ++i ) xvals[i] = 0;
            
            if ( gridInqYvals(gridID1, NULL) )
              gridInqYvals(gridID1, yvals);
            else
              for ( int i = 0; i < ny; ++i ) yvals[i] = 0;
            
	    if ( gridIsRotated(gridID1) || lproj_rll )
	      {
                if ( lproj_rll )
                  {
                    gridInqProjParamRLL(gridID1, &xpole, &ypole, &angle);
                    gridDefProj(gridID2, gridID1);
                  }
                else
                  {
                    xpole = gridInqXpole(gridID1);
                    ypole = gridInqYpole(gridID1);
                    angle = gridInqAngle(gridID1);
                  }

		for ( int j = 0; j < ny; j++ )
		  for ( int i = 0; i < nx; i++ )
		    {
		      xvals2D[j*nx+i] = lamrot_to_lam(yvals[j], xvals[i], ypole, xpole, angle);
		      yvals2D[j*nx+i] = phirot_to_phi(yvals[j], xvals[i], ypole, angle);
		    }	    
	      }
	    else
	      {
		for ( int j = 0; j < ny; j++ )
		  for ( int i = 0; i < nx; i++ )
		    {
		      xvals2D[j*nx+i] = xscale*xvals[i];
		      yvals2D[j*nx+i] = yscale*yvals[j];
		    }

		if ( gridtype == GRID_SINUSOIDAL )
		  {
		    sinusoidal_to_geo(gridsize, xvals2D, yvals2D);
		    /* correct_sinxvals(nx, ny, xvals2D); */
		  }
		else if ( gridtype == GRID_LAEA || lproj_laea )
		  {
		    laea_to_geo(gridID1, gridsize, xvals2D, yvals2D);
		  }
		else if ( gridtype == GRID_LCC2 )
		  {
		    lcc2_to_geo(gridID1, gridsize, xvals2D, yvals2D);
		  }
	      }
	  }

	gridDefXvals(gridID2, xvals2D);
	gridDefYvals(gridID2, yvals2D);

	if ( xvals2D ) Free(xvals2D);
	if ( yvals2D ) Free(yvals2D);

	if ( !lbounds ) goto NO_BOUNDS;

	if ( gridtype == GRID_LCC )
	  {		
	    double *xbounds2D = (double*) Malloc(4*gridsize*sizeof(double));
	    double *ybounds2D = (double*) Malloc(4*gridsize*sizeof(double));

	    for ( int j = 0; j < ny; j++ )
	      for ( int i = 0; i < nx; i++ )
		{
		  index = j*4*nx + 4*i;

		  xbounds2D[index+0] = i+1.5;
		  ybounds2D[index+0] = j+1.5;

		  xbounds2D[index+1] = i+0.5;
		  ybounds2D[index+1] = j+1.5;

		  xbounds2D[index+2] = i+0.5;
		  ybounds2D[index+2] = j+0.5;

		  xbounds2D[index+3] = i+1.5;
		  ybounds2D[index+3] = j+0.5;
		}

	    lcc_to_geo(gridID1, 4*gridsize, xbounds2D, ybounds2D);

	    gridDefXbounds(gridID2, xbounds2D);
	    gridDefYbounds(gridID2, ybounds2D);

	    Free(xbounds2D);
	    Free(ybounds2D);
	  }
	else
	  {
	    if ( gridInqXbounds(gridID1, NULL) )
	      {
		xbounds = (double*) Malloc(2*nx*sizeof(double));
		gridInqXbounds(gridID1, xbounds);
		if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN )
		  if ( check_range(2*nx, xbounds, -720, 720) )
		    {
		      cdoWarning("longitude bounds out of range, skipped!");
		      Free(xbounds);
		      xbounds = NULL;
		    }
	      }
	    else if ( nx > 1 )
	      {
		xbounds = (double*) Malloc(2*nx*sizeof(double));
		grid_gen_bounds(nx, xvals, xbounds);
	      }

	    if ( gridInqYbounds(gridID1, NULL) )
	      {
		ybounds = (double*) Malloc(2*ny*sizeof(double));
		gridInqYbounds(gridID1, ybounds);
		if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN )
		  if ( check_range(2*ny, ybounds, -180, 180) )
		    {
		      cdoWarning("latitude bounds out of range, skipped!");
		      Free(ybounds);
		      ybounds = NULL;
		    }
	      }
	    else if ( ny > 1 )
	      {
		ybounds = (double*) Malloc(2*ny*sizeof(double));
		if ( gridtype == GRID_SINUSOIDAL || 
		     gridtype == GRID_LAEA || lproj_laea || 
		     gridtype == GRID_LCC2 )
		  grid_gen_bounds(ny, yvals, ybounds);
		else
		  {
		    grid_gen_bounds(ny, yvals, ybounds);
		    grid_check_lat_borders(2*ny, ybounds);
		  }
	      }

	    if ( xbounds && ybounds )
	      {
		double *xbounds2D = (double*) Malloc(4*gridsize*sizeof(double));
		double *ybounds2D = (double*) Malloc(4*gridsize*sizeof(double));

		if ( gridIsRotated(gridID1) || lproj_rll )
		  {
		    gridGenRotBounds(xpole, ypole, angle, nx, ny, xbounds, ybounds, xbounds2D, ybounds2D);
		  }
		else
		  {
		    if ( gridtype == GRID_SINUSOIDAL ||
			 gridtype == GRID_LAEA || lproj_laea || 
			 gridtype == GRID_LCC2 )
		      {
			for ( int j = 0; j < ny; j++ )
			  for ( int i = 0; i < nx; i++ )
			    {
			      index = j*4*nx + 4*i;

			      xbounds2D[index+0] = xscale*xbounds[2*i];
			      ybounds2D[index+0] = yscale*ybounds[2*j];

			      xbounds2D[index+1] = xscale*xbounds[2*i];
			      ybounds2D[index+1] = yscale*ybounds[2*j+1];

			      xbounds2D[index+2] = xscale*xbounds[2*i+1];
			      ybounds2D[index+2] = yscale*ybounds[2*j+1];

			      xbounds2D[index+3] = xscale*xbounds[2*i+1];
			      ybounds2D[index+3] = yscale*ybounds[2*j];
			    }
			
			if ( gridtype == GRID_SINUSOIDAL )
			  {
			    sinusoidal_to_geo(4*gridsize, xbounds2D, ybounds2D);
			    /*
			    xvals2D = (double*) Malloc(gridsize*sizeof(double));
			    for ( j = 0; j < 4; ++j )
			      {
				for ( i = 0; i < gridsize; ++i ) xvals2D[i] = xbounds2D[i*4+j];
				correct_sinxvals(nx, ny, xvals2D);
				for ( i = 0; i < gridsize; ++i ) xbounds2D[i*4+j] = xvals2D[i];
			      }
			    Free(xvals2D);
			    */
			  }
			else if ( gridtype == GRID_LAEA || lproj_laea )
			  laea_to_geo(gridID1, 4*gridsize, xbounds2D, ybounds2D);
			else if ( gridtype == GRID_LCC2 )
			  lcc2_to_geo(gridID1, 4*gridsize, xbounds2D, ybounds2D);
		      }
		    else
		      {
			grid_gen_xbounds2D(nx, ny, xbounds, xbounds2D);
			grid_gen_ybounds2D(nx, ny, ybounds, ybounds2D);
		      }
		  }
		
		gridDefXbounds(gridID2, xbounds2D);
		gridDefYbounds(gridID2, ybounds2D);
		
		if ( xbounds )  Free(xbounds);
		if ( ybounds )  Free(ybounds);
		if ( xbounds2D) Free(xbounds2D);
		if ( ybounds2D) Free(ybounds2D);
	      }
	  }

      NO_BOUNDS:

	if ( xvals ) Free(xvals);
	if ( yvals ) Free(yvals);

	gridCopyMask(gridID1, gridID2, gridsize);

	break;
      }
    default:
      {
	Error("Grid type >%s< unsupported!", gridNamePtr(gridtype));
	break;
      }
    }

  return gridID2;
}


int gridToUnstructuredSelecton(int gridID1, int selectionSize, int *selectionIndexList, int nocoords, int nobounds)
{
  /* transform input grid into a unstructured Version if necessary {{{ */
  int unstructuredGridID;
  if (GRID_UNSTRUCTURED == gridInqType(gridID1))
    {
      unstructuredGridID = gridID1;
    }
  else
    {
      unstructuredGridID = gridToUnstructured(gridID1,!nobounds);
    }
  int unstructuredGridSize = gridInqSize(unstructuredGridID);

  int unstructuredSelectionGridID = gridCreate(GRID_UNSTRUCTURED,selectionSize);

  if ( nocoords ) return (unstructuredSelectionGridID);
  /* }}} */

  /* copy meta data of coordinates {{{*/
  char xname[CDI_MAX_NAME], xlongname[CDI_MAX_NAME], xunits[CDI_MAX_NAME];
  char yname[CDI_MAX_NAME], ylongname[CDI_MAX_NAME], yunits[CDI_MAX_NAME];

  gridInqXname(unstructuredGridID, xname);
  gridInqXlongname(unstructuredGridID, xlongname);
  gridInqXunits(unstructuredGridID, xunits);
  gridInqYname(unstructuredGridID, yname);
  gridInqYlongname(unstructuredGridID, ylongname);
  gridInqYunits(unstructuredGridID, yunits);

  gridDefXname(unstructuredSelectionGridID, xname);
  gridDefXlongname(unstructuredSelectionGridID, xlongname);
  gridDefXunits(unstructuredSelectionGridID, xunits);
  gridDefYname(unstructuredSelectionGridID, yname);
  gridDefYlongname(unstructuredSelectionGridID, ylongname);
  gridDefYunits(unstructuredSelectionGridID, yunits);
  /* }}} */

  /* TODO: select bounds */

  /* copy relevant coordinate {{{ */
  double *xvalsUnstructured = (double*) Malloc(unstructuredGridSize*sizeof(double));
  double *yvalsUnstructured = (double*) Malloc(unstructuredGridSize*sizeof(double));
  gridInqXvals(unstructuredGridID, xvalsUnstructured);
  gridInqYvals(unstructuredGridID, yvalsUnstructured);


  gridDefXsize(unstructuredSelectionGridID, selectionSize);
  gridDefYsize(unstructuredSelectionGridID, selectionSize);
  double *xvals   = (double*) Malloc(selectionSize*sizeof(double));
  double *yvals   = (double*) Malloc(selectionSize*sizeof(double));

  for (int i = 0; i < selectionSize; i++)
  {
    xvals[i] = xvalsUnstructured[selectionIndexList[i]];
    yvals[i] = yvalsUnstructured[selectionIndexList[i]];
    /*
    cdoPrint("xval[%d](%d) = %g",i,selectionIndexList[i],xvals[i]);
    cdoPrint("yval[%d](%d) = %g",i,selectionIndexList[i],yvals[i]);
    */
  }
  gridDefXvals(unstructuredSelectionGridID,xvals);
  gridDefYvals(unstructuredSelectionGridID,yvals);
  /* }}} */

  /* copy bounds if requested {{{ */
  if ( ! nobounds )
  {
    int nvertex                 = gridInqNvertex(unstructuredGridID);
    double *xbounds             = (double*) Malloc(nvertex*selectionSize*sizeof(double));
    double *ybounds             = (double*) Malloc(nvertex*selectionSize*sizeof(double));
    double *xboundsUnstructured = (double*) Malloc(nvertex*unstructuredGridSize*sizeof(double));
    double *yboundsUnstructured = (double*) Malloc(nvertex*unstructuredGridSize*sizeof(double));
    gridInqXbounds(unstructuredGridID, xboundsUnstructured);
    gridInqYbounds(unstructuredGridID, yboundsUnstructured);
    for (int i = 0; i < selectionSize; i++)
    {
      for (int k = 0; k < nvertex; k++)
      {
        xbounds[i*nvertex+k] = xboundsUnstructured[selectionIndexList[i]*nvertex+k];
        ybounds[i*nvertex+k] = yboundsUnstructured[selectionIndexList[i]*nvertex+k];
      }
    }
    gridDefNvertex(unstructuredSelectionGridID,nvertex);
    gridDefXbounds(unstructuredSelectionGridID,xbounds);
    gridDefYbounds(unstructuredSelectionGridID,ybounds);

    Free(xboundsUnstructured);
    Free(yboundsUnstructured);
    Free(xbounds);
    Free(ybounds);
  }
  /* }}} */

  Free(xvalsUnstructured);
  Free(yvalsUnstructured);
  Free(xvals);
  Free(yvals);

  return (unstructuredSelectionGridID);
}

int gridToUnstructured(int gridID1, int lbounds)
{
  int gridtype = gridInqType(gridID1);
  int gridsize = gridInqSize(gridID1);
  int gridID2  = gridCreate(GRID_UNSTRUCTURED, gridsize);
  gridDefPrec(gridID2, DATATYPE_FLT32);
	  
  bool lproj_rll = false;
  if ( gridtype == GRID_PROJECTION && gridInqProjType(gridID1) == CDI_PROJ_RLL )
    {
      gridtype = GRID_LONLAT;
      lproj_rll = true;
    }

  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
      {
        double xpole = 0, ypole = 0, angle = 0;
	gridDefXname(gridID2, "lon");
	gridDefYname(gridID2, "lat");
	gridDefXlongname(gridID2, "longitude");
	gridDefYlongname(gridID2, "latitude");
	gridDefXunits(gridID2, "degrees_east");
	gridDefYunits(gridID2, "degrees_north");

	gridDefNvertex(gridID2, 4);

	int nx = gridInqXsize(gridID1);
	int ny = gridInqYsize(gridID1);
	 
	gridDefXsize(gridID2, gridsize);
	gridDefYsize(gridID2, gridsize);

	double *xvals = (double*) Malloc(nx*sizeof(double));
	double *yvals = (double*) Malloc(ny*sizeof(double));

	double *xvals2D = (double*) Malloc(gridsize*sizeof(double));
	double *yvals2D = (double*) Malloc(gridsize*sizeof(double));

	gridInqXvals(gridID1, xvals);
	gridInqYvals(gridID1, yvals);

	if ( gridIsRotated(gridID1) || lproj_rll )
	  {	    
            if ( lproj_rll )
              {
                gridInqProjParamRLL(gridID1, &xpole, &ypole, &angle);
              }
            else
              {
                xpole = gridInqXpole(gridID1);
                ypole = gridInqYpole(gridID1);
                angle = gridInqAngle(gridID1);
              }
		
	    for ( int j = 0; j < ny; j++ )
	      for ( int i = 0; i < nx; i++ )
		{
		  xvals2D[j*nx+i] = lamrot_to_lam(yvals[j], xvals[i], ypole, xpole, angle);
		  yvals2D[j*nx+i] = phirot_to_phi(yvals[j], xvals[i], ypole, angle);
		}
	  }
	else
	  {
	    for ( int j = 0; j < ny; j++ )
	      for ( int i = 0; i < nx; i++ )
		{
		  xvals2D[j*nx+i] = xvals[i];
		  yvals2D[j*nx+i] = yvals[j];
		}
	  }

	gridDefXvals(gridID2, xvals2D);
	gridDefYvals(gridID2, yvals2D);

        Free(xvals2D);
        Free(yvals2D);

	if ( lbounds )
	  {
	    double *xbounds = NULL, *ybounds = NULL;

	    if ( gridInqXbounds(gridID1, NULL) )
	      {
		xbounds = (double*) Malloc(2*nx*sizeof(double));
		gridInqXbounds(gridID1, xbounds);
	      }
	    else if ( nx > 1 )
	      {
		xbounds = (double*) Malloc(2*nx*sizeof(double));
		grid_gen_bounds(nx, xvals, xbounds);
	      }

	    if ( gridInqYbounds(gridID1, NULL) )
	      {
		ybounds = (double*) Malloc(2*ny*sizeof(double));
		gridInqYbounds(gridID1, ybounds);
	      }
	    else if ( ny > 1 )
	      {
		ybounds = (double*) Malloc(2*ny*sizeof(double));
		grid_gen_bounds(ny, yvals, ybounds);
		grid_check_lat_borders(2*ny, ybounds);
	      }

	    if ( xbounds && ybounds )
	      {
		double *xbounds2D = (double*) Malloc(4*gridsize*sizeof(double));
		double *ybounds2D = (double*) Malloc(4*gridsize*sizeof(double));

		if ( gridIsRotated(gridID1) || lproj_rll )
		  {
		    gridGenRotBounds(xpole, ypole, angle, nx, ny, xbounds, ybounds, xbounds2D, ybounds2D);
		  }
		else
		  {
		    grid_gen_xbounds2D(nx, ny, xbounds, xbounds2D);
		    grid_gen_ybounds2D(nx, ny, ybounds, ybounds2D);
		  }

		gridDefXbounds(gridID2, xbounds2D);
		gridDefYbounds(gridID2, ybounds2D);

                Free(xbounds);
                Free(ybounds);
                Free(xbounds2D);
                Free(ybounds2D);
	      }
	    }

        Free(xvals);
        Free(yvals);

	gridCopyMask(gridID1, gridID2, gridsize);

	break;
      }
    case GRID_CURVILINEAR:
      {
	gridID2 = gridDuplicate(gridID1);
	gridChangeType(gridID2, GRID_UNSTRUCTURED);
	gridDefXsize(gridID2, gridsize);
	gridDefYsize(gridID2, gridsize);

	break;
      }
    case GRID_GME:
      {
	int nv = 6;
	double *xbounds = NULL, *ybounds = NULL;

	int nd  = gridInqGMEnd(gridID1);
	int ni  = gridInqGMEni(gridID1);
	int ni2 = gridInqGMEni2(gridID1);
	int ni3 = gridInqGMEni3(gridID1);

	int *imask = (int*) Malloc(gridsize*sizeof(int));
	double *xvals = (double*) Malloc(gridsize*sizeof(double));
	double *yvals = (double*) Malloc(gridsize*sizeof(double));
	if ( lbounds )
	  {
	    xbounds = (double*) Malloc(nv*gridsize*sizeof(double));
	    ybounds = (double*) Malloc(nv*gridsize*sizeof(double));
	  }

	gme_grid(lbounds, gridsize, xvals, yvals, xbounds, ybounds, imask, ni, nd, ni2, ni3);
	
	for ( int i = 0; i < gridsize; i++ )
	  {
	    xvals[i] *= RAD2DEG;
	    yvals[i] *= RAD2DEG;

	    if ( lbounds )
	      for ( int j = 0; j < nv; j++ )
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

	gridDefMaskGME(gridID2, imask);

	gridDefNvertex(gridID2, nv);

	if ( lbounds )
	  {
	    gridDefXbounds(gridID2, xbounds);
	    gridDefYbounds(gridID2, ybounds);
	  }

	gridDefXunits(gridID2, "degrees_east");
	gridDefYunits(gridID2, "degrees_north");

        Free(imask);
        Free(xvals);
        Free(yvals);
	if ( xbounds ) Free(xbounds);
	if ( ybounds ) Free(ybounds);
	
	gridCopyMask(gridID1, gridID2, gridsize);

	break;
      }
    default:
      {
	Error("Grid type >%s< unsupported!", gridNamePtr(gridtype));
	break;
      }
    }

  return gridID2;
}


int gridCurvilinearToRegular(int gridID1)
{
  int gridID2 = -1;
  bool lx = true, ly = true;
	
  int gridtype = gridInqType(gridID1);
  int gridsize = gridInqSize(gridID1);

  if ( gridtype != GRID_CURVILINEAR ) return gridID2;

  int nx = gridInqXsize(gridID1);
  int ny = gridInqYsize(gridID1);
	
  double *xvals2D = (double*) Malloc(gridsize*sizeof(double));
  double *yvals2D = (double*) Malloc(gridsize*sizeof(double));

  gridInqXvals(gridID1, xvals2D);
  gridInqYvals(gridID1, yvals2D);

  double *xvals = (double*) Malloc(nx*sizeof(double));
  double *yvals = (double*) Malloc(ny*sizeof(double));

  for ( int i = 0; i < nx; i++ ) xvals[i] = xvals2D[i];
  for ( int j = 0; j < ny; j++ ) yvals[j] = yvals2D[j*nx];

  for ( int j = 1; j < ny; j++ )
    for ( int i = 0; i < nx; i++ )
      {
	if ( fabs(xvals[i] - xvals2D[j*nx+i]) > 1.e-6 )
	  {
	    lx = FALSE;
	    j = ny;
	    break;
	  }
      }
	
  for ( int i = 1; i < nx; i++ )
    for ( int j = 0; j < ny; j++ )
      {
	if ( fabs(yvals[j] - yvals2D[j*nx+i]) > 1.e-6 )
	  {
	    ly = FALSE;
	    i = nx;
	    break;
	  }
      }

  Free(xvals2D);
  Free(yvals2D);

  if ( lx && ly )
    {
      char xunits[CDI_MAX_NAME], yunits[CDI_MAX_NAME];
      
      gridID2  = gridCreate(GRID_LONLAT, gridsize);
      gridDefXsize(gridID2, nx);
      gridDefYsize(gridID2, ny);
      
      //  gridDefPrec(gridID2, DATATYPE_FLT32);

      gridInqXunits(gridID1, xunits);
      gridInqYunits(gridID1, yunits);

      grid_to_degree(xunits, nx, xvals, "grid1 center lon");
      grid_to_degree(yunits, ny, yvals, "grid1 center lat");

      gridDefXvals(gridID2, xvals);
      gridDefYvals(gridID2, yvals);
    }

  Free(xvals);
  Free(yvals);

  return gridID2;
}


int gridGenWeights(int gridID, double *grid_area, double *grid_wgts)
{
  int status = 0;
  int *grid_mask = NULL;

  int gridtype = gridInqType(gridID);
  int gridsize = gridInqSize(gridID);
  
  if ( gridtype == GRID_GME )
    {
      gridID = gridToUnstructured(gridID, 1);	  
      grid_mask = (int*) Malloc(gridsize*sizeof(int));
      gridInqMaskGME(gridID, grid_mask);
    }

  double total_area = 0;
  int nvals = 0;
  for ( int i = 0; i < gridsize; i++ )
    {
      if ( grid_mask )
	if ( grid_mask[i] == 0 ) continue;
      total_area += grid_area[i];
      nvals++;
    }

  if ( cdoVerbose ) cdoPrint("Total area = %g steradians", total_area);

  for ( int i = 0; i < gridsize; i++ )
    {
      if ( grid_mask )
	if ( grid_mask[i] == 0 )
	  {
	    grid_wgts[i] = 0;
	    continue;
	  }
      
      grid_wgts[i] = grid_area[i] / total_area;
    }
  
  if ( grid_mask ) Free(grid_mask);

  return status;
}


int gridWeightsOld(int gridID, double *weights)
{
  int status = FALSE;

  int len = gridInqSize(gridID);

  if ( gridHasArea(gridID) )
    {
      gridInqArea(gridID, weights);
    }
  else
    {
      int gridtype = gridInqType(gridID);
      if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN )
	{
	  int nlon = gridInqXsize(gridID);
	  int nlat = gridInqYsize(gridID);

	  double *lons = 1 + (double *) Malloc((nlon+2)*sizeof(double));
	  double *lats = 1 + (double *) Malloc((nlat+2)*sizeof(double));

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
	  int datapoint = 0;
	  double sumw = 0;
	  for ( int j = 0; j < nlat; j++ )
	    {
	      double phi1 = (lats[j-1]+lats[j])/2*DEG2RAD;
	      double phi2 = (lats[j+1]+lats[j])/2*DEG2RAD;
	      if ( phi1 < (-1*M_PI/2) ) phi1 = -1*M_PI/2;
	      if ( phi1 > (   M_PI/2) ) phi1 =    M_PI/2;
	      if ( phi2 > (   M_PI/2) ) phi2 =    M_PI/2;
	      if ( phi2 < (-1*M_PI/2) ) phi2 = -1*M_PI/2;
	      double sindphi = sin(phi2)-sin(phi1);
	      for ( int i = 0; i < nlon; i++ )
		{
		  if ( lons[i] >= lons[0]+360 || fabs(lats[j]) > 90 )
		    weights[datapoint] = 0;
		  else
		    {
		      double theta1 = (lons[i-1]+lons[i])/2*DEG2RAD;
		      double theta2 = (lons[i+1]+lons[i])/2*DEG2RAD;
		      weights[datapoint] = fabs((theta2-theta1)*sindphi);
		      sumw += weights[datapoint];
		    }
		  datapoint++;
		}
	    }

	  /* Normalise weights.  */
	  if( IS_NOT_EQUAL(sumw, 0) )
	    for( int i = 0; i < datapoint; i++ ) weights[i] /= sumw;

	  if ( lons-1 ) Free(lons-1);
	  if ( lats-1 ) Free(lats-1);
	}
      else
	{
	  status = TRUE;

	  for ( int i = 0; i < len; i++ ) weights[i] = 1./len;
	}
    }

  return status;
}


int gridWeights(int gridID, double *grid_wgts)
{
  int w_status = 1;
  int a_status = 0;

  int gridsize = gridInqSize(gridID);
  int gridtype = gridInqType(gridID);
  int projtype = (gridtype == GRID_PROJECTION) ? gridInqProjType(gridID) : -1; 
  
  double *grid_area = (double*) Malloc(gridsize*sizeof(double));

  if ( gridHasArea(gridID) )
    {
      if ( cdoVerbose ) cdoPrint("Using existing grid cell area!");
      gridInqArea(gridID, grid_area);
    }
  else
    {
      if ( gridtype == GRID_LONLAT      ||
	   gridtype == GRID_GAUSSIAN    ||
	   gridtype == GRID_LCC         ||
	   gridtype == GRID_LCC2        ||
	   gridtype == GRID_LAEA        ||
	   projtype == CDI_PROJ_LAEA    ||
	   gridtype == GRID_SINUSOIDAL  ||
	   gridtype == GRID_GME         ||
	   gridtype == GRID_CURVILINEAR ||
	   gridtype == GRID_UNSTRUCTURED )
	{
	  a_status = gridGenArea(gridID, grid_area);
	}
      else
	{
	  a_status = 1;
	}
    }

  if ( a_status == 0 )
    {
      w_status = gridGenWeights(gridID, grid_area, grid_wgts);
    }
  else
    {
      for ( int i = 0; i < gridsize; ++i ) grid_wgts[i] = 1./gridsize;
    }
  /*
  for ( i = 0; i < gridsize; ++i ) 
    printf("weights: %d %d %d %g %g\n", a_status, w_status, i, grid_area[i], grid_wgts[i]);
  */
  Free(grid_area);

  return w_status;
}


bool grid_is_distance_generic(int gridID)
{
  bool status = false;

  if ( gridInqType(gridID) == GRID_GENERIC )
    {
      char xunits[CDI_MAX_NAME];
      gridInqXunits(gridID, xunits);
      char yunits[CDI_MAX_NAME];
      gridInqXunits(gridID, yunits);

      if ( strcmp(xunits, "m") == 0 && strcmp(yunits, "m") == 0 &&
           gridInqXvals(gridID, NULL) && gridInqYvals(gridID, NULL) )
        status = true;
    }

  return status;
}

/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2017 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
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
#include "cdo_int.h"
#include "grid.h"

#define  PARAM_MISSVAL  -9.E33

int nfc_to_nlat(int nfc, int ntr)
{
  int nlat = nfc / (ntr+1);
  nlat /= 2;

  return nlat;
}


int nlat_to_ntr(int nlat)
{
  int ntr = (nlat*2 - 1) / 3;

  return ntr;
}


int nlat_to_ntr_linear(int nlat)
{
  int ntr = (nlat*2 - 1) / 2;

  return ntr;
}


int ntr_to_nlat(int ntr)
{
  int nlat = (int)lround((ntr*3.+1.)/2.);
  if ( (nlat % 2) > 0 )
    {
      nlat  = nlat + 1;
      /*
      int nlat2 = (int)lround(((ntr+1)*3.+1.)/2.);
      if ( nlat == nlat2 )
	cdoAbort("Computation of latitudes failed for truncation %d", ntr);
      */
    }

  return nlat;
}


int ntr_to_nlat_linear(int ntr)
{
  int nlat = (int)lround((ntr*2.+1.)/2.);
  if ( (nlat % 2) > 0 )
    {
      nlat  = nlat + 1;
      /*
      int nlat2 = (int)lround(((ntr+1)*2.+1.)/2.);
      if ( nlat == nlat2 )
	cdoAbort("Computation of latitudes failed for truncation %d", ntr);
      */
    }

  return nlat;
}


int nlat_to_nlon(int nlat)
{
  int nlon = 2 * nlat;

  /* check that FFT works with nlon */
  while ( 1 )
    {
      int n = nlon;
      if    ( n % 8 == 0 )  { n /= 8; }
      while ( n % 6 == 0 )  { n /= 6; }
      while ( n % 5 == 0 )  { n /= 5; }
      while ( n % 4 == 0 )  { n /= 4; }
      while ( n % 3 == 0 )  { n /= 3; }
      if    ( n % 2 == 0 )  { n /= 2; }

      if ( n <= 8 ) break;

      nlon = nlon + 2;
      /*
      if ( nlon > 9999 )
	{
	  nlon = 2 * nlat;
	  fprintf(stderr, "FFT does not work with len %d!\n", nlon);
	  break;
	}
      */
    }

  return nlon;
}


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


void grid_copy_attributes(int gridID1, int gridID2)
{
  char string[CDI_MAX_NAME];
  string[0] = 0;   cdiGridInqKeyStr(gridID1, CDI_KEY_XDIMNAME, CDI_MAX_NAME, string);
  if ( string[0] ) cdiGridDefKeyStr(gridID2, CDI_KEY_XDIMNAME, strlen(string)+1, string);
  string[0] = 0;   cdiGridInqKeyStr(gridID1, CDI_KEY_YDIMNAME, CDI_MAX_NAME, string);
  if ( string[0] ) cdiGridDefKeyStr(gridID2, CDI_KEY_YDIMNAME, strlen(string)+1, string);
  string[0] = 0;   cdiGridInqKeyStr(gridID1, CDI_KEY_VDIMNAME, CDI_MAX_NAME, string);
  if ( string[0] ) cdiGridDefKeyStr(gridID2, CDI_KEY_VDIMNAME, strlen(string)+1, string);
  string[0] = 0;   cdiGridInqKeyStr(gridID1, CDI_KEY_XNAME, CDI_MAX_NAME, string);
  if ( string[0] ) cdiGridDefKeyStr(gridID2, CDI_KEY_XNAME, strlen(string)+1, string);
  string[0] = 0;   cdiGridInqKeyStr(gridID1, CDI_KEY_YNAME, CDI_MAX_NAME, string);
  if ( string[0] ) cdiGridDefKeyStr(gridID2, CDI_KEY_YNAME, strlen(string)+1, string);
  string[0] = 0;   cdiGridInqKeyStr(gridID1, CDI_KEY_XLONGNAME, CDI_MAX_NAME, string);
  if ( string[0] ) cdiGridDefKeyStr(gridID2, CDI_KEY_XLONGNAME, strlen(string)+1, string);
  string[0] = 0;   cdiGridInqKeyStr(gridID1, CDI_KEY_YLONGNAME, CDI_MAX_NAME, string);
  if ( string[0] ) cdiGridDefKeyStr(gridID2, CDI_KEY_YLONGNAME, strlen(string)+1, string);
  string[0] = 0;   cdiGridInqKeyStr(gridID1, CDI_KEY_XUNITS, CDI_MAX_NAME, string);
  if ( string[0] ) cdiGridDefKeyStr(gridID2, CDI_KEY_XUNITS, strlen(string)+1, string);
  string[0] = 0;   cdiGridInqKeyStr(gridID1, CDI_KEY_YUNITS, CDI_MAX_NAME, string);
  if ( string[0] ) cdiGridDefKeyStr(gridID2, CDI_KEY_YUNITS, strlen(string)+1, string);
}


void grid_copy_mapping(int gridID1, int gridID2)
{
  char string[CDI_MAX_NAME];
  string[0] = 0;   cdiGridInqKeyStr(gridID1, CDI_KEY_MAPPING, CDI_MAX_NAME, string);
  if ( string[0] ) cdiGridDefKeyStr(gridID2, CDI_KEY_MAPPING, strlen(string)+1, string);
  string[0] = 0;   cdiGridInqKeyStr(gridID1, CDI_KEY_MAPNAME, CDI_MAX_NAME, string);
  if ( string[0] ) cdiGridDefKeyStr(gridID2, CDI_KEY_MAPNAME, strlen(string)+1, string);

  cdiCopyAtts(gridID1, CDI_GLOBAL, gridID2, CDI_GLOBAL);
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
      cdoAbort("Gridtype %s unsupported!", gridNamePtr(gridtype));
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
      cdoAbort("Gridtype %s unsupported!", gridNamePtr(gridtype));
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
  (void)xunitstr;

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
  (void)yunitstr;

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
void sinu_to_geo(int gridsize, double *xvals, double *yvals)
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


void grid_def_param_sinu(int gridID)
{
  const char *projection = "sinusoidal";
  cdiGridDefKeyStr(gridID, CDI_KEY_MAPNAME, (int)strlen(projection)+1, projection);
  const char *mapvarname = "Sinusoidal";
  cdiGridDefKeyStr(gridID, CDI_KEY_MAPPING, (int)strlen(mapvarname)+1, mapvarname);

  cdiDefAttTxt(gridID, CDI_GLOBAL, "grid_mapping_name", (int)strlen(projection), projection);
}


void grid_def_param_laea(int gridID, double a, double lon_0, double lat_0)
{
  const char *projection = "lambert_azimuthal_equal_area";
  cdiGridDefKeyStr(gridID, CDI_KEY_MAPNAME, (int)strlen(projection)+1, projection);
  const char *mapvarname = "Lambert_AEA";
  cdiGridDefKeyStr(gridID, CDI_KEY_MAPPING, (int)strlen(mapvarname)+1, mapvarname);

  cdiDefAttTxt(gridID, CDI_GLOBAL, "grid_mapping_name", (int)strlen(projection), projection);
  
  cdiDefAttFlt(gridID, CDI_GLOBAL, "earth_radius", CDI_DATATYPE_FLT64, 1, &a);
  cdiDefAttFlt(gridID, CDI_GLOBAL, "longitude_of_projection_origin", CDI_DATATYPE_FLT64, 1, &lon_0);
  cdiDefAttFlt(gridID, CDI_GLOBAL, "latitude_of_projection_origin", CDI_DATATYPE_FLT64, 1, &lat_0);
}

static
bool cdiInqAttConvertedToFloat(int gridID, int atttype, const char *attname, int attlen, double *attflt)
{
  bool status = true;

  if ( atttype == CDI_DATATYPE_INT32 )
    {
      int attint[attlen];
      cdiInqAttInt(gridID, CDI_GLOBAL, attname, attlen, attint);
      for ( int i = 0; i < attlen; ++i ) attflt[i] = (double)attint[i];
    }
  else if ( atttype == CDI_DATATYPE_FLT32 || atttype == CDI_DATATYPE_FLT64 )
    {
      cdiInqAttFlt(gridID, CDI_GLOBAL, attname, attlen, attflt);
    }
  else
    {
      status = false;
    }

  return status;
}

static
void grid_inq_param_laea(int gridID, double *a, double *lon_0, double *lat_0, double *x_0, double *y_0)
{
  *a = 0; *lon_0 = 0; *lat_0 = 0, *x_0 = 0, *y_0 = 0;

  int gridtype = gridInqType(gridID);
  if ( gridtype == GRID_PROJECTION )
    {
      const char *projection = "lambert_azimuthal_equal_area";
      char mapping[CDI_MAX_NAME]; mapping[0] = 0;
      cdiGridInqKeyStr(gridID, CDI_KEY_MAPNAME, CDI_MAX_NAME, mapping);
      if ( mapping[0] && strcmp(mapping, projection) == 0 )
        {
          int atttype, attlen;
          char attname[CDI_MAX_NAME+1];

          int natts;
          cdiInqNatts(gridID, CDI_GLOBAL, &natts);

          for ( int iatt = 0; iatt < natts; ++iatt )
            {
              cdiInqAtt(gridID, CDI_GLOBAL, iatt, attname, &atttype, &attlen);
              if ( attlen != 1 ) continue;

              double attflt;
              if ( cdiInqAttConvertedToFloat(gridID, atttype, attname, attlen, &attflt) )
                {
                  if      ( strcmp(attname, "earth_radius") == 0 )                    *a     = attflt;
                  else if ( strcmp(attname, "longitude_of_projection_origin") == 0 )  *lon_0 = attflt;
                  else if ( strcmp(attname, "latitude_of_projection_origin") == 0 )   *lat_0 = attflt;
                  else if ( strcmp(attname, "false_easting")  == 0 )                  *x_0   = attflt;
                  else if ( strcmp(attname, "false_northing") == 0 )                  *y_0   = attflt;
                }
            }
        }
      else
        cdoWarning("%s mapping parameter missing!", projection);
    }
}

#if defined(HAVE_LIBPROJ)
static
bool grid_inq_param_lcc(int gridID, double *a, double *rf, double *lon_0, double *lat_0, double *lat_1, double *lat_2, double *x_0, double *y_0)
{
  bool status = false;
  double xlon_0 = PARAM_MISSVAL, ylat_0 = PARAM_MISSVAL;
  *a = PARAM_MISSVAL; *rf = PARAM_MISSVAL; *lon_0 = PARAM_MISSVAL; *lat_0 = PARAM_MISSVAL;
  *lat_1 = PARAM_MISSVAL, *lat_2 = PARAM_MISSVAL;
  *x_0 = PARAM_MISSVAL, *y_0 = PARAM_MISSVAL;

  int gridtype = gridInqType(gridID);
  if ( gridtype == GRID_PROJECTION )
    {
      const char *projection = "lambert_conformal_conic";
      char mapping[CDI_MAX_NAME]; mapping[0] = 0;
      cdiGridInqKeyStr(gridID, CDI_KEY_MAPNAME, CDI_MAX_NAME, mapping);
      if ( mapping[0] && strcmp(mapping, projection) == 0 )
        {
          int atttype, attlen;
          char attname[CDI_MAX_NAME+1];

          int natts;
          cdiInqNatts(gridID, CDI_GLOBAL, &natts);

          for ( int iatt = 0; iatt < natts; ++iatt )
            {
              cdiInqAtt(gridID, CDI_GLOBAL, iatt, attname, &atttype, &attlen);
              if ( attlen > 2 ) continue;

              double attflt[2];
              if ( cdiInqAttConvertedToFloat(gridID, atttype, attname, attlen, attflt) )
                {
                  if      ( strcmp(attname, "earth_radius") == 0 )                       *a     = attflt[0];
                  else if ( strcmp(attname, "inverse_flattening") == 0 )                 *rf    = attflt[0];
                  else if ( strcmp(attname, "longitude_of_central_meridian") == 0 )      *lon_0 = attflt[0];
                  else if ( strcmp(attname, "latitude_of_projection_origin") == 0 )      *lat_0 = attflt[0];
                  else if ( strcmp(attname, "false_easting")  == 0 )                     *x_0   = attflt[0];
                  else if ( strcmp(attname, "false_northing") == 0 )                     *y_0   = attflt[0];
                  else if ( strcmp(attname, "longitudeOfFirstGridPointInDegrees") == 0 ) xlon_0 = attflt[0];
                  else if ( strcmp(attname, "latitudeOfFirstGridPointInDegrees")  == 0 ) ylat_0 = attflt[0];
                  else if ( strcmp(attname, "standard_parallel") == 0 )
                    {
                      *lat_1 = attflt[0];
                      *lat_2 = (attlen == 2) ? attflt[1] : attflt[0];
                    }
                }
            }

          status = true;
          if ( IS_EQUAL(*lon_0,PARAM_MISSVAL) ) { status = false; cdoWarning("%s mapping parameter %s missing!", projection, "longitude_of_central_meridian"); }
          if ( IS_EQUAL(*lat_0,PARAM_MISSVAL) ) { status = false; cdoWarning("%s mapping parameter %s missing!", projection, "latitude_of_central_meridian"); }
          if ( IS_EQUAL(*lat_1,PARAM_MISSVAL) ) { status = false; cdoWarning("%s mapping parameter %s missing!", projection, "standard_parallel"); }
          if ( status && IS_EQUAL(*x_0,PARAM_MISSVAL) && IS_EQUAL(*y_0,PARAM_MISSVAL) && IS_NOT_EQUAL(xlon_0,PARAM_MISSVAL) && IS_NOT_EQUAL(ylat_0,PARAM_MISSVAL) )
             { status = false; cdoWarning("%s mapping parameter %s missing!", projection, "false_easting and false_northing"); }
        }
      else
        cdoWarning("%s mapping parameter missing!", projection);
    }

  return status;
}
#endif


int proj_lonlat_to_lcc(double missval, double lon_0, double lat_0, double lat_1, double lat_2,
                       double a, double rf, int nvals, double *xvals, double *yvals)
{
  int status = 0;
#if defined(HAVE_LIBPROJ)
  char *params[20];

  int nbpar = 0;
  params[nbpar++] = gen_param("proj=lcc");
  if ( IS_NOT_EQUAL(a, missval) && a > 0 ) params[nbpar++] = gen_param("a=%g", a);
  if ( IS_NOT_EQUAL(rf, missval) && rf > 0 ) params[nbpar++] = gen_param("rf=%g", rf);
  params[nbpar++] = gen_param("lon_0=%g", lon_0);
  params[nbpar++] = gen_param("lat_0=%g", lat_0);
  params[nbpar++] = gen_param("lat_1=%g", lat_1);
  params[nbpar++] = gen_param("lat_2=%g", lat_2);
  params[nbpar++] = gen_param("units=m");
  //  params[nbpar++] = gen_param("no_defs");

  if ( cdoVerbose )
    for ( int i = 0; i < nbpar; ++i )
      cdoPrint("Proj.param[%d] = %s", i+1, params[i]);
  
  projPJ proj = pj_init(nbpar, &params[0]);
  if ( !proj ) status = -1;

  for ( int i = 0; i < nbpar; ++i ) Free(params[i]);

  /* proj->over = 1; */		/* allow longitude > 180 */

  if ( status == 0 )
    {
      projUV p;
      for ( int i = 0; i < nvals; i++ )
        {
          p.u = xvals[i]*DEG_TO_RAD;
          p.v = yvals[i]*DEG_TO_RAD;
          p = pj_fwd(p, proj);
          xvals[i] = p.u;
          yvals[i] = p.v;
        }

      pj_free(proj);
    }
#else
  status = -1;
#endif

  if ( status == -1 )
    for ( int i = 0; i < nvals; i++ )
      {
        xvals[i] = missval;
        yvals[i] = missval;
      }

  return status;
}

static
void lonlat_to_lcc(double missval, double lon_0, double lat_0, double lat_1, double lat_2,
                   double a, double rf, int nvals, double *xvals, double *yvals)
{
  int status = proj_lonlat_to_lcc(missval, lon_0, lat_0, lat_1, lat_2, a, rf, nvals, xvals, yvals);
#if defined(HAVE_LIBPROJ)
  if ( status == -1 ) cdoAbort("proj error: %s", pj_strerrno(pj_errno));
#else
  if ( status == -1 ) cdoAbort("proj4 support not compiled in!");
#endif
}


int cdo_lonlat_to_lcc(int gridID, int nvals, double *xvals, double *yvals)
{
  double lon_0, lat_0, lat_1, lat_2, a, rf, xval_0, yval_0, x_0, y_0;

  gridInqParamLCC(gridID, grid_missval, &lon_0, &lat_0, &lat_1, &lat_2, &a, &rf, &xval_0, &yval_0, &x_0, &y_0);

  lonlat_to_lcc(grid_missval, lon_0, lat_0, lat_1, lat_2, a, rf, nvals, xvals, yvals);

  return 0;
}

static
void laea_to_geo(int gridID, int gridsize, double *xvals, double *yvals)
{
#if defined(HAVE_LIBPROJ)
  char *params[20];
  
  double a, lon_0, lat_0, x_0, y_0;
  grid_inq_param_laea(gridID, &a, &lon_0, &lat_0, &x_0, &y_0);

  int nbpar = 0;
  params[nbpar++] = gen_param("proj=laea");
  if ( a > 0 ) params[nbpar++] = gen_param("a=%g", a);
  params[nbpar++] = gen_param("lon_0=%g", lon_0);
  params[nbpar++] = gen_param("lat_0=%g", lat_0);
  if ( IS_NOT_EQUAL(x_0,0) ) params[nbpar++] = gen_param("x_0=%g", x_0);
  if ( IS_NOT_EQUAL(y_0,0) ) params[nbpar++] = gen_param("y_0=%g", y_0);

  if ( cdoVerbose )
    for ( int i = 0; i < nbpar; ++i )
      cdoPrint("Proj.param[%d] = %s", i+1, params[i]);

  projPJ proj = pj_init(nbpar, &params[0]);
  if ( !proj ) cdoAbort("proj error: %s", pj_strerrno(pj_errno));

  for ( int i = 0; i < nbpar; ++i ) Free(params[i]);

  /* proj->over = 1; */		/* allow longitude > 180 */

  projUV p;
  for ( int i = 0; i < gridsize; i++ )
    {
      p.u = xvals[i];
      p.v = yvals[i];
      p = pj_inv(p, proj);
      xvals[i] = p.u*RAD_TO_DEG;
      yvals[i] = p.v*RAD_TO_DEG;
      if ( xvals[i] < -9000. || xvals[i] > 9000. ) xvals[i] = -9999.;
      if ( yvals[i] < -9000. || yvals[i] > 9000. ) yvals[i] = -9999.;
    }

  pj_free(proj);
#else
  cdoAbort("proj4 support not compiled in!");
#endif
}


int lcc_to_lonlat(int gridID, int gridsize, double *xvals, double *yvals)
{
  bool errstatus = false;
#if defined(HAVE_LIBPROJ)
  char *params[20];

  double a, rf, lon_0, lat_0, lat_1, lat_2, x_0, y_0;
  bool status = grid_inq_param_lcc(gridID, &a, &rf, &lon_0, &lat_0, &lat_1, &lat_2, &x_0, &y_0);
  if ( status == false ) cdoAbort("mapping parameter missing!");

  int nbpar = 0;
  params[nbpar++] = gen_param("proj=lcc");
  if ( IS_NOT_EQUAL(a,PARAM_MISSVAL)  && a  > 0 ) params[nbpar++] = gen_param("a=%g", a);
  if ( IS_NOT_EQUAL(rf,PARAM_MISSVAL) && rf > 0 ) params[nbpar++] = gen_param("rf=%g", rf);
  params[nbpar++] = gen_param("lon_0=%g", lon_0);
  params[nbpar++] = gen_param("lat_0=%g", lat_0);
  params[nbpar++] = gen_param("lat_1=%g", lat_1);
  params[nbpar++] = gen_param("lat_2=%g", lat_2);
  if ( IS_NOT_EQUAL(x_0,PARAM_MISSVAL) ) params[nbpar++] = gen_param("x_0=%g", x_0);
  if ( IS_NOT_EQUAL(y_0,PARAM_MISSVAL) ) params[nbpar++] = gen_param("y_0=%g", y_0);

  if ( cdoVerbose )
    for ( int i = 0; i < nbpar; ++i )
      cdoPrint("Proj.param[%d] = %s", i+1, params[i]);
  
  projPJ proj = pj_init(nbpar, &params[0]);
  if ( !proj ) cdoAbort("proj error: %s", pj_strerrno(pj_errno));

  for ( int i = 0; i < nbpar; ++i ) Free(params[i]);

  /* proj->over = 1; */		/* allow longitude > 180 */
  
  projUV p;
  for ( int i = 0; i < gridsize; i++ )
    {
      p.u = xvals[i];
      p.v = yvals[i];
      p = pj_inv(p, proj);
      xvals[i] = p.u*RAD_TO_DEG;
      yvals[i] = p.v*RAD_TO_DEG;
    }

  pj_free(proj);
#else
  errstatus = true;
  cdoAbort("proj4 support not compiled in!");
#endif
  return !errstatus;
}

static
void proj_to_geo(char *proj4param, int gridsize, double *xvals, double *yvals)
{
#if defined(HAVE_LIBPROJ)
  char *params[99];

  int nbpar;
  for ( nbpar = 0; nbpar < 99; ++nbpar )
    {
      while ( *proj4param == ' ' || *proj4param == '+' ) proj4param++;
      if ( *proj4param == 0 ) break;
      char *cstart = proj4param;
      while ( *proj4param != ' ' && *proj4param != 0 ) proj4param++;
      char *cend = proj4param;
      size_t len = cend - cstart;
      if ( len <= 0 ) break;
      bool lend = *cend == 0;
      if ( !lend ) *cend = 0;
      params[nbpar] = strdup(cstart);
      if ( !lend ) *cend = ' ';
    }

  if ( cdoVerbose )
    for ( int i = 0; i < nbpar; ++i )
      cdoPrint("Proj.param[%d] = %s", i+1, params[i]);

  projPJ proj = pj_init(nbpar, &params[0]);
  if ( !proj ) cdoAbort("proj error: %s", pj_strerrno(pj_errno));

  for ( int i = 0; i < nbpar; ++i ) Free(params[i]);

  projUV p;
  for ( int i = 0; i < gridsize; i++ )
    {
      p.u = xvals[i];
      p.v = yvals[i];
      p = pj_inv(p, proj);
      xvals[i] = p.u*RAD_TO_DEG;
      yvals[i] = p.v*RAD_TO_DEG;
      if ( xvals[i] < -9000. || xvals[i] > 9000. ) xvals[i] = -9999.;
      if ( yvals[i] < -9000. || yvals[i] > 9000. ) yvals[i] = -9999.;
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
  if ( gridtype != GRID_GAUSSIAN_REDUCED ) cdoAbort("Not a reduced Gaussian grid!");

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

  if ( gridInqSize(gridID2) != nx*ny ) cdoAbort("Gridsize differ!");

  Free(rowlon);
}


int gridToRegular(int gridID1)
{
  int nx = 0;
  double *xvals = NULL;

  int gridtype = gridInqType(gridID1);
  if ( gridtype != GRID_GAUSSIAN_REDUCED ) cdoAbort("Not a reduced Gaussian grid!");

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


char *grid_get_proj4param(int gridID)
{
  char *proj4param = NULL;

  int gridtype = gridInqType(gridID);
  if ( gridtype == GRID_PROJECTION )
    {
      int atttype, attlen, atttxtlen = 0;
      char *atttxt = NULL;
      char attname[CDI_MAX_NAME+1];

      int natts;
      cdiInqNatts(gridID, CDI_GLOBAL, &natts);

      for ( int iatt = 0; iatt < natts; ++iatt )
        {
          cdiInqAtt(gridID, CDI_GLOBAL, iatt, attname, &atttype, &attlen);

          if ( atttype == CDI_DATATYPE_TXT )
            {
              if ( attlen > atttxtlen )
                {
                  atttxt = (char*) Realloc(atttxt, (attlen+1));
                  atttxtlen = attlen;
                }
              cdiInqAttTxt(gridID, CDI_GLOBAL, attname, attlen, atttxt);
              atttxt[attlen] = 0;
              if ( strcmp(attname, "proj4_params") == 0 )
                {
                  proj4param = atttxt;
                  break;
                }
            }
        }
      if ( proj4param == NULL && atttxt ) Free(atttxt);
    }

  return proj4param;
}


int gridToCurvilinear(int gridID1, int lbounds)
{
  int index;
  int gridtype = gridInqType(gridID1);

  int nx = gridInqXsize(gridID1);
  int ny = gridInqYsize(gridID1);

  bool lxyvals = gridInqXvals(gridID1, NULL) && gridInqYvals(gridID1, NULL);
  if ( !lxyvals && gridtype != GRID_LCC ) cdoAbort("Grid coordinates missing!");

  int gridsize = gridInqSize(gridID1);
  int gridID2 = gridCreate(GRID_CURVILINEAR, gridsize);
  gridDefPrec(gridID2, CDI_DATATYPE_FLT32);

  char *proj4param = NULL;
  bool lproj4     = false;
  bool lproj_rll  = false;
  bool lproj_laea = false;
  bool lproj_lcc  = false;
  bool lproj_sinu = false;
  if ( gridtype == GRID_PROJECTION && gridsize == nx*ny )
    {
      proj4param = grid_get_proj4param(gridID1);
      if ( proj4param )
        {
          lproj4 = true;
          gridtype = GRID_LONLAT;
        }
      else
        {
          int projtype = gridInqProjType(gridID1);
          if      ( projtype == CDI_PROJ_RLL  ) lproj_rll  = true;
          else if ( projtype == CDI_PROJ_LAEA ) lproj_laea = true;
          else if ( projtype == CDI_PROJ_LCC  ) lproj_lcc  = true;
          else if ( projtype == CDI_PROJ_SINU ) lproj_sinu = true;

          if ( lproj_rll || lproj_laea || lproj_lcc || lproj_sinu ) gridtype = GRID_LONLAT;
        }
    }

  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_LCC:
      {
        double xpole = 0, ypole = 0, angle = 0;
	double xscale = 1, yscale = 1;
	double *xvals = NULL, *yvals = NULL;
	double *xbounds = NULL, *ybounds = NULL;
	char xunits[CDI_MAX_NAME], yunits[CDI_MAX_NAME];

	gridInqXunits(gridID1, xunits);
	gridInqYunits(gridID1, yunits);

	if ( lproj_laea || lproj_lcc || lproj_sinu )
	  {
	    int len;
	    len = (int) strlen(xunits);
            bool lvalid_xunits = (len == 1 && memcmp(xunits, "m",  1) == 0) ||
                                 (len == 2 && memcmp(xunits, "km", 2) == 0);
	    len = (int) strlen(yunits);
            bool lvalid_yunits = (len == 1 && memcmp(yunits, "m",  1) == 0) ||
                                 (len == 2 && memcmp(yunits, "km", 2) == 0);

	    if ( !lvalid_xunits )
	      cdoWarning("Possibly wrong result! Invalid x-coordinate units: \"%s\" (expected \"m\" or \"km\")", xunits);
	    if ( !lvalid_yunits )
	      cdoWarning("Possibly wrong result! Invalid y-coordinate units: \"%s\" (expected \"m\" or \"km\")", yunits);
	  }

	if ( memcmp(xunits, "km", 2) == 0 ) xscale = 1000;
	if ( memcmp(yunits, "km", 2) == 0 ) yscale = 1000;

	gridDefXsize(gridID2, nx);
	gridDefYsize(gridID2, ny);

	double *xvals2D = (double*) Malloc(gridsize*sizeof(double));
	double *yvals2D = (double*) Malloc(gridsize*sizeof(double));

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
            
        if ( lproj_rll )
          {
            gridInqParamRLL(gridID1, &xpole, &ypole, &angle);
            gridDefProj(gridID2, gridID1);

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

            if      ( lproj_sinu ) sinu_to_geo(gridsize, xvals2D, yvals2D);
            else if ( lproj_laea ) laea_to_geo(gridID1, gridsize, xvals2D, yvals2D);
            else if ( lproj_lcc  ) lcc_to_lonlat(gridID1, gridsize, xvals2D, yvals2D);
            else if ( lproj4     ) proj_to_geo(proj4param, gridsize, xvals2D, yvals2D);
          }

	gridDefXvals(gridID2, xvals2D);
	gridDefYvals(gridID2, yvals2D);

	if ( xvals2D ) Free(xvals2D);
	if ( yvals2D ) Free(yvals2D);

	if ( !lbounds ) goto NO_BOUNDS;

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
            if ( lproj_sinu || lproj_laea || lproj_lcc || lproj4 )
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

            if ( lproj_rll )
              {
                gridGenRotBounds(xpole, ypole, angle, nx, ny, xbounds, ybounds, xbounds2D, ybounds2D);
              }
            else
              {
                if ( lproj_sinu || lproj_laea || lproj_lcc || lproj4 )
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
			
                    if      ( lproj_sinu ) sinu_to_geo(4*gridsize, xbounds2D, ybounds2D);
                    else if ( lproj_laea ) laea_to_geo(gridID1, 4*gridsize, xbounds2D, ybounds2D);
                    else if ( lproj_lcc  ) lcc_to_lonlat(gridID1, 4*gridsize, xbounds2D, ybounds2D);
                    else if ( lproj4     ) proj_to_geo(proj4param, 4*gridsize, xbounds2D, ybounds2D);
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

      NO_BOUNDS:

	if ( xvals ) Free(xvals);
	if ( yvals ) Free(yvals);

	gridCopyMask(gridID1, gridID2, gridsize);

	break;
      }
    default:
      {
	cdoAbort("Grid type >%s< unsupported!", gridNamePtr(gridtype));
	break;
      }
    }

  if ( proj4param ) Free(proj4param);

  return gridID2;
}


int gridToUnstructuredSelecton(int gridID1, int selectionSize, int *selectionIndexList, int nocoords, int nobounds)
{
  /* transform input grid into a unstructured Version if necessary {{{ */
  int unstructuredGridID;
  if (GRID_UNSTRUCTURED == gridInqType(gridID1))
    unstructuredGridID = gridID1;
  else
    unstructuredGridID = gridToUnstructured(gridID1,!nobounds);

  int unstructuredGridSize = gridInqSize(unstructuredGridID);

  int unstructuredSelectionGridID = gridCreate(GRID_UNSTRUCTURED,selectionSize);

  if ( nocoords ) return unstructuredSelectionGridID;
  /* }}} */

  /* copy meta data of coordinates {{{*/
  grid_copy_attributes(unstructuredGridID, unstructuredSelectionGridID);
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

  return unstructuredSelectionGridID;
}


int gridToUnstructured(int gridID1, int lbounds)
{
  int gridtype = gridInqType(gridID1);
  int gridsize = gridInqSize(gridID1);
  int gridID2  = gridCreate(GRID_UNSTRUCTURED, gridsize);
  gridDefPrec(gridID2, CDI_DATATYPE_FLT32);
	  
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

	if ( lproj_rll )
	  {	    
            gridInqParamRLL(gridID1, &xpole, &ypole, &angle);
		
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

		if ( lproj_rll )
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

        int nd, ni, ni2, ni3;
        gridInqParamGME(gridID1, &nd, &ni, &ni2, &ni3);

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
	cdoAbort("Grid type >%s< unsupported!", gridNamePtr(gridtype));
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
      gridID2  = gridCreate(GRID_LONLAT, gridsize);
      gridDefXsize(gridID2, nx);
      gridDefYsize(gridID2, ny);
      
      //  gridDefPrec(gridID2, CDI_DATATYPE_FLT32);

      char xunits[CDI_MAX_NAME]; xunits[0] = 0;
      char yunits[CDI_MAX_NAME]; yunits[0] = 0;
      cdiGridInqKeyStr(gridID1, CDI_KEY_XUNITS, CDI_MAX_NAME, xunits);
      cdiGridInqKeyStr(gridID1, CDI_KEY_YUNITS, CDI_MAX_NAME, yunits);

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

  if ( cdoVerbose ) cdoPrint("Total area = %g", total_area);

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
	   projtype == CDI_PROJ_RLL     ||
	   projtype == CDI_PROJ_LAEA    ||
	   projtype == CDI_PROJ_LCC     ||
	   projtype == CDI_PROJ_SINU    ||
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

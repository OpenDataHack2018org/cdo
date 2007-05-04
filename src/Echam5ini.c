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

#if  defined  (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <string.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#if  defined  (HAVE_LIBNETCDF)
#  include "netcdf.h"
#endif

static int nvars_ml = 4;
static int nvars_sfc = 18;
static int filetype_ml = 1;
static int filetype_sfc = 2;
static const char strfiletype_ml[]  = "Initial file spectral";
static const char strfiletype_sfc[] = "Initial file surface";
static const char strfiletype_ini[] = "Initial file";
static const char strfiletype_res[] = "Restart history file";


typedef struct {
  int  gridtype;
  int  zaxistype;
  int  code;
  char *name;
  char *longname;
  char *units;
  int gridID;
  int zaxisID;
  int gridsize;
  int nlev;
  double *ptr;
} VAR;


void inivar(VAR *var, int gridtype, int zaxistype, int code, const char *name,
	       const char *longname, const char *units)
{
  static char func[] = "inivar_ml";
  
  var->gridtype  = gridtype;
  var->zaxistype = zaxistype;
  var->code      = code;
  var->name      = NULL;
  var->longname  = NULL;
  var->units     = NULL;
  if ( name )     var->name      = strdup(name);
  if ( longname ) var->longname  = strdup(longname);
  if ( units )    var->units     = strdup(units);
}

void inivars_ml(VAR **vars)
{
  static char func[] = "inivars_ml";

  *vars = (VAR *) malloc((nvars_ml+1)*sizeof(VAR));

  inivar(&(*vars)[0], GRID_GAUSSIAN, ZAXIS_HYBRID,  133, "Q",   "specific humidity", "kg/kg");
  inivar(&(*vars)[1], GRID_SPECTRAL, ZAXIS_HYBRID,  138, "SVO", "vorticity", "1/s");
  inivar(&(*vars)[2], GRID_SPECTRAL, ZAXIS_HYBRID,  155, "SD",  "divergence", "1/s");
  inivar(&(*vars)[3], GRID_SPECTRAL, ZAXIS_HYBRID,  130, "STP", "temperature", "K");
  /* Don't change the order (lsp must be the last one)! */
  inivar(&(*vars)[4], GRID_SPECTRAL, ZAXIS_SURFACE, 152, "LSP", "log surface pressure", "K");
}


#if  defined  (HAVE_LIBNETCDF)
static void nce(int istat)
{
  /*
    This routine provides a simple interface to netCDF error message routine.
  */

  if ( istat != NC_NOERR ) cdoAbort(nc_strerror(istat));
}
#endif


int read_e5ml(const char *filename, VAR **vars)
{
  static char func[] = "read_e5ml";
  int nvars = 0;
#if  defined  (HAVE_LIBNETCDF)
  int nc_dim_id, nc_var_id;
  size_t dimlen, nvals;
  size_t start[3];
  size_t count[3];
  int nlon, nlat, nlev, nlevp1, nvct, nsp, i, iv;
  int gridIDgp, gridIDsp, zaxisIDml, zaxisIDsfc;
  int gridtype, zaxistype;
  int nc_file_id;
  char filetype[256];
  size_t attlen;
  double *xvals, *yvals, *vct, *levs;

  /* open file and check file type */
  nce(nc_open(filename, NC_NOWRITE, &nc_file_id));

  nce(nc_get_att_text(nc_file_id, NC_GLOBAL, "file_type", filetype));
  nce(nc_inq_attlen(nc_file_id, NC_GLOBAL, "file_type", &attlen));
  filetype[attlen] = 0;

  if ( strcmp(filetype, strfiletype_ml) != 0 ) return (0);

  inivars_ml(vars);

  /* read dimensions */

  nce(nc_inq_dimid(nc_file_id, "lon", &nc_dim_id));
  nce(nc_inq_dimlen(nc_file_id, nc_dim_id, &dimlen));
  nlon = (int) dimlen;

  nce(nc_inq_dimid(nc_file_id, "lat", &nc_dim_id));
  nce(nc_inq_dimlen(nc_file_id, nc_dim_id, &dimlen));
  nlat = (int) dimlen;

  gridIDgp = gridCreate(GRID_GAUSSIAN, nlon*nlat);
  gridDefXsize(gridIDgp, nlon);
  gridDefYsize(gridIDgp, nlat);

  nce(nc_inq_dimid(nc_file_id, "nsp", &nc_dim_id));
  nce(nc_inq_dimlen(nc_file_id, nc_dim_id, &dimlen));
  nsp = (int) dimlen;

  gridIDsp = gridCreate(GRID_SPECTRAL, nsp*2);

  nce(nc_inq_dimid(nc_file_id, "nlev", &nc_dim_id));
  nce(nc_inq_dimlen(nc_file_id, nc_dim_id, &dimlen));
  nlev = (int) dimlen;
  nlevp1 = nlev + 1;
  nvct = nlevp1*2;

  zaxisIDsfc = zaxisCreate(ZAXIS_SURFACE, 1);
  zaxisIDml  = zaxisCreate(ZAXIS_HYBRID, nlev);

  levs = (double *) malloc(nlev*sizeof(double));
  for ( i = 0; i < nlev; i++ ) levs[i] = i+1;
  zaxisDefLevels(zaxisIDml, levs);
  free(levs);

  /* read variables */

  xvals = (double *) malloc(nlon*sizeof(double));
  yvals = (double *) malloc(nlat*sizeof(double));

  nce(nc_inq_varid(nc_file_id, "lon", &nc_var_id));
  nce(nc_get_var_double(nc_file_id, nc_var_id, xvals));

  nce(nc_inq_varid(nc_file_id, "lat", &nc_var_id));
  nce(nc_get_var_double(nc_file_id, nc_var_id, yvals));

  gridDefXvals(gridIDgp, xvals);
  gridDefYvals(gridIDgp, yvals);

  free(xvals);
  free(yvals);

  vct   = (double *) malloc(nvct*sizeof(double));

  nce(nc_inq_varid(nc_file_id, "vct_a", &nc_var_id));
  nce(nc_get_var_double(nc_file_id, nc_var_id, vct));

  nce(nc_inq_varid(nc_file_id, "vct_b", &nc_var_id));
  nce(nc_get_var_double(nc_file_id, nc_var_id, vct+nlevp1));

  zaxisDefVct(zaxisIDml, 2*nlevp1, vct);
  free(vct);

  for ( iv = 0; iv < nvars_ml; iv++ )
    {
      nvals = 0;

      gridtype  = (*vars)[iv].gridtype;
      zaxistype = (*vars)[iv].zaxistype;

      if ( gridtype == GRID_GAUSSIAN )
	{
	  (*vars)[iv].gridID = gridIDgp;
	  nvals += nlon*nlat;
	}
      else
	{
	  (*vars)[iv].gridID = gridIDsp;
	  nvals += nsp*2;
	}

      (*vars)[iv].zaxisID   = zaxisIDml;
      (*vars)[iv].gridsize  = nvals;
      (*vars)[iv].nlev      = nlev;

      (*vars)[iv].ptr = (double *) malloc(nlev*nvals*sizeof(double));
      
      for ( i = 0; i < nlev; i++ )
	{
	  if ( gridtype == GRID_GAUSSIAN )
	    {
	      start[0] = 0;     start[1] = i;  start[2] = 0;
	      count[0] = nlat;  count[1] = 1;  count[2] = nlon;     
	    }
	  else
	    {
	      start[0] = 0;     start[1] = 0;  start[2] = i;
	      count[0] = nsp;   count[1] = 2;  count[2] = 1;
	    }

	  nce(nc_inq_varid(nc_file_id, (*vars)[iv].name, &nc_var_id));
	  nce(nc_get_vara_double(nc_file_id, nc_var_id, start, count, (*vars)[iv].ptr+i*nvals));
	}
    }

  /* read lsp */

  (*vars)[nvars_ml].gridID    = gridIDsp;
  (*vars)[nvars_ml].zaxisID   = zaxisIDsfc;
  (*vars)[nvars_ml].gridsize  = nsp*2;
  (*vars)[nvars_ml].nlev      = 1;

  start[0] = 0;    start[1] = 0;  start[2] = nlev;
  count[0] = nsp;  count[1] = 2;  count[2] = 1;

  (*vars)[nvars_ml].ptr = (double *) malloc(nsp*2*sizeof(double));

  nce(nc_inq_varid(nc_file_id, "STP", &nc_var_id));
  nce(nc_get_vara_double(nc_file_id, nc_var_id, start, count, (*vars)[nvars_ml].ptr));

  /*close input file */
  nce(nc_close(nc_file_id));

  nvars = nvars_ml + 1;

#else
  cdoAbort("netCDF support not compiled in!");
#endif

  return (nvars);
}


int read_e5sfc(const char *filename, VAR **vars)
{
  static char func[] = "read_e5sfc";
  int nvars = 0;
#if  defined  (HAVE_LIBNETCDF)
  int nc_dim_id, nc_var_id;
  size_t dimlen, nvals;
  size_t start[3];
  size_t count[3];
  int nlon, nlat, nlev, nlevp1, nvct, nsp, i, iv;
  int gridIDgp, gridIDsp, zaxisIDml, zaxisIDsfc;
  int gridtype, zaxistype;
  int nc_file_id;
  char filetype[256];
  size_t attlen;
  double *xvals, *yvals, *vct, *levs;

  /* open file and check file type */
  nce(nc_open(filename, NC_NOWRITE, &nc_file_id));

  nce(nc_get_att_text(nc_file_id, NC_GLOBAL, "file_type", filetype));
  nce(nc_inq_attlen(nc_file_id, NC_GLOBAL, "file_type", &attlen));
  filetype[attlen] = 0;

  if ( strcmp(filetype, strfiletype_sfc) != 0 ) return (0);

  printf("%s\n", filetype);

  /*close input file */
  nce(nc_close(nc_file_id));

  nvars = 0;

#else
  cdoAbort("netCDF support not compiled in!");
#endif

  return (nvars);
}


int read_e5ini(const char *filename, VAR **vars)
{
  static char func[] = "read_e5ini";
  int nvars = 0;
#if  defined  (HAVE_LIBNETCDF)
  int nc_dim_id, nc_var_id;
  size_t dimlen, nvals;
  size_t start[3];
  size_t count[3];
  int nlon, nlat, nlev, nlevp1, nvct, nsp, i, iv;
  int gridIDgp, gridIDsp, zaxisIDml, zaxisIDsfc;
  int gridtype, zaxistype;
  int nc_file_id;
  char filetype[256];
  size_t attlen;
  double *xvals, *yvals, *vct, *levs;

  /* open file and check file type */
  nce(nc_open(filename, NC_NOWRITE, &nc_file_id));

  nce(nc_get_att_text(nc_file_id, NC_GLOBAL, "file_type", filetype));
  nce(nc_inq_attlen(nc_file_id, NC_GLOBAL, "file_type", &attlen));
  filetype[attlen] = 0;

  if ( strcmp(filetype, strfiletype_ini) != 0 ) return (0);

  printf("%s\n", filetype);

  /*close input file */
  nce(nc_close(nc_file_id));

  nvars = 0;

#else
  cdoAbort("netCDF support not compiled in!");
#endif

  return (nvars);
}


int read_e5res(const char *filename, VAR **vars)
{
  static char func[] = "read_e5res";
  int nvars = 0;
#if  defined  (HAVE_LIBNETCDF)
  int nc_dim_id, nc_var_id;
  int varid;
  size_t dimlen, nvals;
  size_t start[3];
  size_t count[3];
  int nlon, nlat, nlev, nvct, nsp, i, iv;
  int gridIDgp, gridIDsp, zaxisIDml, zaxisIDsfc;
  int gridtype, zaxistype;
  int nc_file_id;
  char filetype[256];
  size_t attlen;
  double *xvals, *yvals, *vct, *levs;
  int ncvarid;
  int ndims, ngatts, unlimdimid;
  int nvdims, nvatts;
  int dimidsp[9];
  int max_vars = 4096;
  nc_type xtype;
  char name[256];
  int lon_dimid, lat_dimid, nhgl_dimid, nlevp1_dimid, spc_dimid, nvclev_dimid;
  int complex_dimid, nmp1_dimid, belowsurface_dimid, lev_dimid, ilev_dimid;
  int surface_dimid, /*height2m_dimid, height10m_dimid,*/ n2_dimid;
  int lon, lat, nhgl, nlevp1, spc, nvclev;
  int complex, nmp1, belowsurface, lev, ilev;
  int surface, /*height2m, height10m,*/ n2;

  /* open file and check file type */
  nce(nc_open(filename, NC_NOWRITE, &nc_file_id));

  nce(nc_get_att_text(nc_file_id, NC_GLOBAL, "file_type", filetype));
  nce(nc_inq_attlen(nc_file_id, NC_GLOBAL, "file_type", &attlen));
  filetype[attlen] = 0;

  if ( strncmp(filetype, strfiletype_res, strlen(strfiletype_res)) != 0 ) return (0);

  printf("%s\n", filetype);

  nce(nc_inq_dimid(nc_file_id, "lon", &lon_dimid));
  nce(nc_inq_dimlen(nc_file_id, lon_dimid, &lon));

  nce(nc_inq_dimid(nc_file_id, "lat", &lat_dimid));
  nce(nc_inq_dimlen(nc_file_id, lat_dimid, &lat));

  nce(nc_inq_dimid(nc_file_id, "nhgl", &nhgl_dimid));
  nce(nc_inq_dimlen(nc_file_id, nhgl_dimid, &nhgl));

  nce(nc_inq_dimid(nc_file_id, "nlevp1", &nlevp1_dimid));
  nce(nc_inq_dimlen(nc_file_id, nlevp1_dimid, &nlevp1));

  nce(nc_inq_dimid(nc_file_id, "spc", &spc_dimid));
  nce(nc_inq_dimlen(nc_file_id, spc_dimid, &spc));

  nce(nc_inq_dimid(nc_file_id, "nvclev", &nvclev_dimid));
  nce(nc_inq_dimlen(nc_file_id, nvclev_dimid, &nvclev));

  nce(nc_inq_dimid(nc_file_id, "complex", &complex_dimid));
  nce(nc_inq_dimlen(nc_file_id, complex_dimid, &complex));

  nce(nc_inq_dimid(nc_file_id, "nmp1", &nmp1_dimid));
  nce(nc_inq_dimlen(nc_file_id, nmp1_dimid, &nmp1));

  nce(nc_inq_dimid(nc_file_id, "belowsurface", &belowsurface_dimid));
  nce(nc_inq_dimlen(nc_file_id, belowsurface_dimid, &belowsurface));

  nce(nc_inq_dimid(nc_file_id, "lev", &lev_dimid));
  nce(nc_inq_dimlen(nc_file_id, lev_dimid, &lev));

  nce(nc_inq_dimid(nc_file_id, "ilev", &ilev_dimid));
  nce(nc_inq_dimlen(nc_file_id, ilev_dimid, &ilev));

  nce(nc_inq_dimid(nc_file_id, "surface", &surface_dimid));
  nce(nc_inq_dimlen(nc_file_id, surface_dimid, &surface));
  /*
  nce(nc_inq_dimid(nc_file_id, "height2m", &height2m_dimid));
  nce(nc_inq_dimlen(nc_file_id, height2m_dimid, &height2m));

  nce(nc_inq_dimid(nc_file_id, "height10m", &height10m_dimid));
  nce(nc_inq_dimlen(nc_file_id, height10m_dimid, &height10m));
  */
  nce(nc_inq_dimid(nc_file_id, "n2", &n2_dimid));
  nce(nc_inq_dimlen(nc_file_id, n2_dimid, &n2));

  /* define gaussian grid */

  nlon = lon;
  nlat = lat;

  gridIDgp = gridCreate(GRID_GAUSSIAN, nlon*nlat);
  gridDefXsize(gridIDgp, nlon);
  gridDefYsize(gridIDgp, nlat);

  xvals = (double *) malloc(nlon*sizeof(double));
  yvals = (double *) malloc(nlat*sizeof(double));

  nce(nc_inq_varid(nc_file_id, "lon", &nc_var_id));
  nce(nc_get_var_double(nc_file_id, nc_var_id, xvals));

  nce(nc_inq_varid(nc_file_id, "lat", &nc_var_id));
  nce(nc_get_var_double(nc_file_id, nc_var_id, yvals));

  gridDefXvals(gridIDgp, xvals);
  gridDefYvals(gridIDgp, yvals);

  free(xvals);
  free(yvals);

  /* define surface level */

  zaxisIDsfc = zaxisCreate(ZAXIS_SURFACE, 1);

  /* define model level */

  nlev = lev;
  nvct = nvclev*2;

  zaxisIDml  = zaxisCreate(ZAXIS_HYBRID, nlev);

  levs = (double *) malloc(nlev*sizeof(double));
  for ( i = 0; i < nlev; i++ ) levs[i] = i+1;
  zaxisDefLevels(zaxisIDml, levs);
  free(levs);

  vct   = (double *) malloc(nvct*sizeof(double));

  nce(nc_inq_varid(nc_file_id, "vct_a", &nc_var_id));
  nce(nc_get_var_double(nc_file_id, nc_var_id, vct));

  nce(nc_inq_varid(nc_file_id, "vct_b", &nc_var_id));
  nce(nc_get_var_double(nc_file_id, nc_var_id, vct+nlevp1));

  zaxisDefVct(zaxisIDml, 2*nlevp1, vct);
  free(vct);


  nce(nc_inq(nc_file_id, &ndims, &nvars, &ngatts, &unlimdimid));

  *vars = (VAR *) malloc(max_vars*sizeof(VAR));

  varid = 0;
  for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {      
      nce(nc_inq_var(nc_file_id, ncvarid, name, &xtype, &nvdims, dimidsp, &nvatts));

      if ( nvdims == 4 )
	{
	  if ( dimidsp[0] == nhgl_dimid    && dimidsp[1] == nmp1_dimid &&
	       dimidsp[2] == complex_dimid && dimidsp[3] == lev_dimid )
	    {
	      /* varid++; */
	    }
	  else if ( dimidsp[0] == nhgl_dimid    && dimidsp[1] == nmp1_dimid &&
	            dimidsp[2] == complex_dimid && dimidsp[3] == ilev_dimid )
	    {
	      /* varid++; */
	    }
	  else
	    {
	      printf(" unsupported");
	    }
	}
      else if ( nvdims == 3 )
	{
	  if ( dimidsp[0] == lat_dimid && dimidsp[1] == lev_dimid &&
	       dimidsp[2] == lon_dimid)
	    {
	      nvals = nlon*nlat;
  
	      inivar(&(*vars)[varid], GRID_GAUSSIAN, ZAXIS_HYBRID,  0, name, NULL, NULL);

	      (*vars)[varid].gridID    = gridIDgp;
	      (*vars)[varid].zaxisID   = zaxisIDml;
	      (*vars)[varid].gridsize  = nvals;
	      (*vars)[varid].nlev      = nlev;

	      (*vars)[varid].ptr = (double *) malloc(nlev*nvals*sizeof(double));

	      for ( i = 0; i < nlev; i++ )
		{
		  start[0] = 0;     start[1] = i;  start[2] = 0;
		  count[0] = nlat;  count[1] = 1;  count[2] = nlon;     

		  nce(nc_inq_varid(nc_file_id, name, &nc_var_id));
		  nce(nc_get_vara_double(nc_file_id, nc_var_id, start, count, (*vars)[varid].ptr+i*nvals));
		}

	      varid++;
	    }
	  else if ( dimidsp[0] == lat_dimid && dimidsp[1] == ilev_dimid &&
		    dimidsp[2] == lon_dimid)
	    {
	      /* varid++; */
	    }
	  else if ( dimidsp[0] == lat_dimid && dimidsp[1] == belowsurface_dimid &&
		    dimidsp[2] == lon_dimid)
	    {
	      /* varid++; */
	    }
	  else if ( dimidsp[0] == lat_dimid && dimidsp[1] == n2_dimid &&
		    dimidsp[2] == lon_dimid)
	    {
	      /* varid++; */
	    }
	  else
	    {
	      printf(" unsupported");
	    }
	}
      else if ( nvdims == 2 )
	{
	  if ( dimidsp[0] == lat_dimid && dimidsp[1] == lon_dimid)
	    {
	      nvals = nlon*nlat;
  
	      inivar(&(*vars)[varid], GRID_GAUSSIAN, ZAXIS_SURFACE,  0, name, NULL, NULL);

	      (*vars)[varid].gridID    = gridIDgp;
	      (*vars)[varid].zaxisID   = zaxisIDsfc;
	      (*vars)[varid].gridsize  = nvals;
	      (*vars)[varid].nlev      = 1;

	      (*vars)[varid].ptr = (double *) malloc(nvals*sizeof(double));

	      for ( i = 0; i < nlev; i++ )
		{
		  nce(nc_inq_varid(nc_file_id, name, &nc_var_id));
		  nce(nc_get_var_double(nc_file_id, nc_var_id, (*vars)[varid].ptr));
		}

	      varid++;
	    }
	  else if ( dimidsp[0] == nhgl_dimid && dimidsp[1] == lev_dimid)
	    {
	      /* varid++; */
	    }
	  else
	    {
	      printf(" unsupported");
	    }
	}
      else if ( nvdims == 1 )
	{
	  if ( dimidsp[0] == lat_dimid && strcmp(name, "lat") == 0 )
	    {
	    }
	  else if ( dimidsp[0] == lon_dimid && strcmp(name, "lon") == 0 )
	    {
	    }
	  else if ( dimidsp[0] == nvclev_dimid && strcmp(name, "vct_a") == 0 )
	    {
	    }
	  else if ( dimidsp[0] == nvclev_dimid && strcmp(name, "vct_b") == 0 )
	    {
	    }
	  else
	    {
	      printf(" unsupported");
	    }
	}
      else
	{
	}
    }
  printf("nvars = %d\n", varid);

  /*close input file */
  nce(nc_close(nc_file_id));

  nvars = varid;

#else
  cdoAbort("netCDF support not compiled in!");
#endif

  return (nvars);
}


int write_e5res(const char *filename, VAR **vars)
{
  static char func[] = "write_e5res";
  int nvars = 0;
#if  defined  (HAVE_LIBNETCDF)
  int nc_dim_id, nc_var_id;
  int varid;
  size_t dimlen, nvals;
  size_t start[3];
  size_t count[3];
  int nlon, nlat, nlev, nvct, nsp, i, iv;
  int gridIDgp, gridIDsp, zaxisIDml, zaxisIDsfc;
  int gridtype, zaxistype;
  int nc_file_id;
  char filetype[256];
  size_t attlen;
  double *xvals, *yvals, *vct, *levs;
  int ncvarid;
  int ndims, ngatts, unlimdimid;
  int nvdims, nvatts;
  int dimidsp[9];
  int max_vars = 4096;
  nc_type xtype;
  char name[256];
  int lon_dimid, lat_dimid, nhgl_dimid, nlevp1_dimid, spc_dimid, nvclev_dimid;
  int complex_dimid, nmp1_dimid, belowsurface_dimid, lev_dimid, ilev_dimid;
  int surface_dimid, /*height2m_dimid, height10m_dimid,*/ n2_dimid;
  int lon, lat, nhgl, nlevp1, spc, nvclev;
  int complex, nmp1, belowsurface, lev, ilev;
  int surface, /*height2m, height10m,*/ n2;

  /* open file and check file type */
  nce(nc_open(filename, NC_NOWRITE, &nc_file_id));

  nce(nc_get_att_text(nc_file_id, NC_GLOBAL, "file_type", filetype));
  nce(nc_inq_attlen(nc_file_id, NC_GLOBAL, "file_type", &attlen));
  filetype[attlen] = 0;

  if ( strncmp(filetype, strfiletype_res, strlen(strfiletype_res)) != 0 ) return (0);

  printf("%s\n", filetype);

  nce(nc_inq_dimid(nc_file_id, "lon", &lon_dimid));
  nce(nc_inq_dimlen(nc_file_id, lon_dimid, &lon));

  nce(nc_inq_dimid(nc_file_id, "lat", &lat_dimid));
  nce(nc_inq_dimlen(nc_file_id, lat_dimid, &lat));

  nce(nc_inq_dimid(nc_file_id, "nhgl", &nhgl_dimid));
  nce(nc_inq_dimlen(nc_file_id, nhgl_dimid, &nhgl));

  nce(nc_inq_dimid(nc_file_id, "nlevp1", &nlevp1_dimid));
  nce(nc_inq_dimlen(nc_file_id, nlevp1_dimid, &nlevp1));

  nce(nc_inq_dimid(nc_file_id, "spc", &spc_dimid));
  nce(nc_inq_dimlen(nc_file_id, spc_dimid, &spc));

  nce(nc_inq_dimid(nc_file_id, "nvclev", &nvclev_dimid));
  nce(nc_inq_dimlen(nc_file_id, nvclev_dimid, &nvclev));

  nce(nc_inq_dimid(nc_file_id, "complex", &complex_dimid));
  nce(nc_inq_dimlen(nc_file_id, complex_dimid, &complex));

  nce(nc_inq_dimid(nc_file_id, "nmp1", &nmp1_dimid));
  nce(nc_inq_dimlen(nc_file_id, nmp1_dimid, &nmp1));

  nce(nc_inq_dimid(nc_file_id, "belowsurface", &belowsurface_dimid));
  nce(nc_inq_dimlen(nc_file_id, belowsurface_dimid, &belowsurface));

  nce(nc_inq_dimid(nc_file_id, "lev", &lev_dimid));
  nce(nc_inq_dimlen(nc_file_id, lev_dimid, &lev));

  nce(nc_inq_dimid(nc_file_id, "ilev", &ilev_dimid));
  nce(nc_inq_dimlen(nc_file_id, ilev_dimid, &ilev));

  nce(nc_inq_dimid(nc_file_id, "surface", &surface_dimid));
  nce(nc_inq_dimlen(nc_file_id, surface_dimid, &surface));
  /*
  nce(nc_inq_dimid(nc_file_id, "height2m", &height2m_dimid));
  nce(nc_inq_dimlen(nc_file_id, height2m_dimid, &height2m));

  nce(nc_inq_dimid(nc_file_id, "height10m", &height10m_dimid));
  nce(nc_inq_dimlen(nc_file_id, height10m_dimid, &height10m));
  */
  nce(nc_inq_dimid(nc_file_id, "n2", &n2_dimid));
  nce(nc_inq_dimlen(nc_file_id, n2_dimid, &n2));

  /* define gaussian grid */

  nlon = lon;
  nlat = lat;

  gridIDgp = gridCreate(GRID_GAUSSIAN, nlon*nlat);
  gridDefXsize(gridIDgp, nlon);
  gridDefYsize(gridIDgp, nlat);

  xvals = (double *) malloc(nlon*sizeof(double));
  yvals = (double *) malloc(nlat*sizeof(double));

  nce(nc_inq_varid(nc_file_id, "lon", &nc_var_id));
  nce(nc_get_var_double(nc_file_id, nc_var_id, xvals));

  nce(nc_inq_varid(nc_file_id, "lat", &nc_var_id));
  nce(nc_get_var_double(nc_file_id, nc_var_id, yvals));

  gridDefXvals(gridIDgp, xvals);
  gridDefYvals(gridIDgp, yvals);

  free(xvals);
  free(yvals);

  /* define surface level */

  zaxisIDsfc = zaxisCreate(ZAXIS_SURFACE, 1);

  /* define model level */

  nlev = lev;
  nvct = nvclev*2;

  zaxisIDml  = zaxisCreate(ZAXIS_HYBRID, nlev);

  levs = (double *) malloc(nlev*sizeof(double));
  for ( i = 0; i < nlev; i++ ) levs[i] = i+1;
  zaxisDefLevels(zaxisIDml, levs);
  free(levs);

  vct   = (double *) malloc(nvct*sizeof(double));

  nce(nc_inq_varid(nc_file_id, "vct_a", &nc_var_id));
  nce(nc_get_var_double(nc_file_id, nc_var_id, vct));

  nce(nc_inq_varid(nc_file_id, "vct_b", &nc_var_id));
  nce(nc_get_var_double(nc_file_id, nc_var_id, vct+nlevp1));

  zaxisDefVct(zaxisIDml, 2*nlevp1, vct);
  free(vct);


  nce(nc_inq(nc_file_id, &ndims, &nvars, &ngatts, &unlimdimid));

  *vars = (VAR *) malloc(max_vars*sizeof(VAR));

  varid = 0;
  for ( ncvarid = 0; ncvarid < nvars; ncvarid++ )
    {      
      nce(nc_inq_var(nc_file_id, ncvarid, name, &xtype, &nvdims, dimidsp, &nvatts));

      if ( nvdims == 4 )
	{
	  if ( dimidsp[0] == nhgl_dimid    && dimidsp[1] == nmp1_dimid &&
	       dimidsp[2] == complex_dimid && dimidsp[3] == lev_dimid )
	    {
	      /* varid++; */
	    }
	  else if ( dimidsp[0] == nhgl_dimid    && dimidsp[1] == nmp1_dimid &&
	            dimidsp[2] == complex_dimid && dimidsp[3] == ilev_dimid )
	    {
	      /* varid++; */
	    }
	  else
	    {
	      printf(" unsupported");
	    }
	}
      else if ( nvdims == 3 )
	{
	  if ( dimidsp[0] == lat_dimid && dimidsp[1] == lev_dimid &&
	       dimidsp[2] == lon_dimid)
	    {
	      nvals = nlon*nlat;
  
	      inivar(&(*vars)[varid], GRID_GAUSSIAN, ZAXIS_HYBRID,  0, name, NULL, NULL);

	      (*vars)[varid].gridID    = gridIDgp;
	      (*vars)[varid].zaxisID   = zaxisIDml;
	      (*vars)[varid].gridsize  = nvals;
	      (*vars)[varid].nlev      = nlev;

	      (*vars)[varid].ptr = (double *) malloc(nlev*nvals*sizeof(double));

	      for ( i = 0; i < nlev; i++ )
		{
		  start[0] = 0;     start[1] = i;  start[2] = 0;
		  count[0] = nlat;  count[1] = 1;  count[2] = nlon;     

		  nce(nc_inq_varid(nc_file_id, name, &nc_var_id));
		  nce(nc_get_vara_double(nc_file_id, nc_var_id, start, count, (*vars)[varid].ptr+i*nvals));
		}

	      varid++;
	    }
	  else if ( dimidsp[0] == lat_dimid && dimidsp[1] == ilev_dimid &&
		    dimidsp[2] == lon_dimid)
	    {
	      /* varid++; */
	    }
	  else if ( dimidsp[0] == lat_dimid && dimidsp[1] == belowsurface_dimid &&
		    dimidsp[2] == lon_dimid)
	    {
	      /* varid++; */
	    }
	  else if ( dimidsp[0] == lat_dimid && dimidsp[1] == n2_dimid &&
		    dimidsp[2] == lon_dimid)
	    {
	      /* varid++; */
	    }
	  else
	    {
	      printf(" unsupported");
	    }
	}
      else if ( nvdims == 2 )
	{
	  if ( dimidsp[0] == lat_dimid && dimidsp[1] == lon_dimid)
	    {
	      nvals = nlon*nlat;
  
	      inivar(&(*vars)[varid], GRID_GAUSSIAN, ZAXIS_SURFACE,  0, name, NULL, NULL);

	      (*vars)[varid].gridID    = gridIDgp;
	      (*vars)[varid].zaxisID   = zaxisIDsfc;
	      (*vars)[varid].gridsize  = nvals;
	      (*vars)[varid].nlev      = 1;

	      (*vars)[varid].ptr = (double *) malloc(nvals*sizeof(double));

	      for ( i = 0; i < nlev; i++ )
		{
		  nce(nc_inq_varid(nc_file_id, name, &nc_var_id));
		  nce(nc_get_var_double(nc_file_id, nc_var_id, (*vars)[varid].ptr));
		}

	      varid++;
	    }
	  else if ( dimidsp[0] == nhgl_dimid && dimidsp[1] == lev_dimid)
	    {
	      /* varid++; */
	    }
	  else
	    {
	      printf(" unsupported");
	    }
	}
      else if ( nvdims == 1 )
	{
	  if ( dimidsp[0] == lat_dimid && strcmp(name, "lat") == 0 )
	    {
	    }
	  else if ( dimidsp[0] == lon_dimid && strcmp(name, "lon") == 0 )
	    {
	    }
	  else if ( dimidsp[0] == nvclev_dimid && strcmp(name, "vct_a") == 0 )
	    {
	    }
	  else if ( dimidsp[0] == nvclev_dimid && strcmp(name, "vct_b") == 0 )
	    {
	    }
	  else
	    {
	      printf(" unsupported");
	    }
	}
      else
	{
	}
    }
  printf("nvars = %d\n", varid);

  /*close input file */
  nce(nc_close(nc_file_id));

  nvars = varid;

#else
  cdoAbort("netCDF support not compiled in!");
#endif

  return (nvars);
}


void *Echam5ini(void *argument)
{
  static char func[] = "Echam5ini";
  int operatorID;
  int operfunc;
  int READ_E5ML,  READ_E5SFC,  READ_E5INI,  READ_E5RES;
  int WRITE_E5ML, WRITE_E5SFC, WRITE_E5INI, WRITE_E5RES;
  int streamID1, streamID2 = CDI_UNDEFID;
  int nrecs = 0;
  int recID, varID, levelID;
  int vlistID1, vlistID2;
  int nvars = 0;
  int iv, nlev;
  int gridsize, nmiss;
  int taxisID, tsID;
  double *array = NULL;

  cdoInitialize(argument);

  READ_E5ML   = cdoOperatorAdd("read_e5ml",   func_read,  0, NULL);
  READ_E5SFC  = cdoOperatorAdd("read_e5sfc",  func_read,  0, NULL);
  READ_E5INI  = cdoOperatorAdd("read_e5ini",  func_read,  0, NULL);
  READ_E5RES  = cdoOperatorAdd("read_e5res",  func_read,  0, NULL);
  WRITE_E5ML  = cdoOperatorAdd("write_e5ml",  func_write, 0, NULL);
  WRITE_E5SFC = cdoOperatorAdd("write_e5sfc", func_write, 0, NULL);
  WRITE_E5INI = cdoOperatorAdd("write_e5ini", func_write, 0, NULL);
  WRITE_E5RES = cdoOperatorAdd("write_e5res", func_write, 0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

  if ( operfunc == func_read )
    {
      VAR *vars = NULL;

      if ( operatorID == READ_E5ML )
	nvars = read_e5ml(cdoStreamName(0), &vars);
      else if ( operatorID == READ_E5SFC )
	nvars = read_e5sfc(cdoStreamName(0), &vars);
      else if ( operatorID == READ_E5INI )
	nvars = read_e5ini(cdoStreamName(0), &vars);
      else if ( operatorID == READ_E5RES )
	nvars = read_e5res(cdoStreamName(0), &vars);
      else
	cdoAbort("Internal problem!");

      if ( nvars == 0 ) cdoAbort("Unsupported file type!");
      
      printf("nvars = %d\n", nvars);

      vlistID2 = vlistCreate();
      vlistDefNtsteps(vlistID2, 0);

      for ( iv = 0; iv < nvars; iv++ )
	{/*
	  fprintf(stderr, "%d %s %d %d %d %d\n", iv, vars[iv].name, vars[iv].gridID, vars[iv].zaxisID, gridInqSize(vars[iv].gridID), zaxisInqSize(vars[iv].zaxisID));
	 */
	  varID = vlistDefVar(vlistID2, vars[iv].gridID, vars[iv].zaxisID, TIME_CONSTANT);
	  if ( vars[iv].code > 0 ) vlistDefVarCode(vlistID2, varID, vars[iv].code);
	  if ( vars[iv].name )     vlistDefVarName(vlistID2, varID, vars[iv].name);
	  if ( vars[iv].longname ) vlistDefVarLongname(vlistID2, varID, vars[iv].longname);
          if ( vars[iv].units )    vlistDefVarUnits(vlistID2, varID, vars[iv].units);
	  vlistDefVarDatatype(vlistID2, varID, DATATYPE_FLT64);
	}

      taxisID = taxisCreate(TAXIS_ABSOLUTE);
      vlistDefTaxis(vlistID2, taxisID);

      if ( cdoDefaultFileType == CDI_UNDEFID )
	cdoDefaultFileType = FILETYPE_NC;

      streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
      if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

      streamDefVlist(streamID2, vlistID2);

      tsID = 0;
      streamDefTimestep(streamID2, tsID);

      for ( varID = 0; varID < nvars; varID++ )
	{
	  gridsize = vars[varID].gridsize;
	  nlev     = vars[varID].nlev;

	  for ( levelID = 0; levelID < nlev; levelID++ )
	    {
	      streamDefRecord(streamID2,  varID,  levelID);
	      streamWriteRecord(streamID2, vars[varID].ptr+levelID*gridsize, 0);
	    }
	}

      streamClose(streamID2);

      vlistDestroy(vlistID2);
    }
  else if ( operatorID == WRITE_E5ML )
    {
      VAR *vars = NULL;

      streamID1 = streamOpenRead(cdoStreamName(0));
      if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

      vlistID1 = streamInqVlist(streamID1);

      gridsize = vlistGridsizeMax(vlistID1);
      array = (double *) malloc(gridsize*sizeof(double));

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array, &nmiss);
	}
      /*
      if ( operatorID == WRITE_E5ML )
	write_e5ml(cdoStreamName(1), &vars, nvars);
      else if ( operatorID == WRITE_E5SFC )
	write_e5sfc(cdoStreamName(1), &vars, nvars);
      else if ( operatorID == WRITE_E5INI )
	write_e5ini(cdoStreamName(1), &vars, nvars);
      else if ( operatorID == WRITE_E5RES )
	write_e5res(cdoStreamName(1), &vars, nvars);
      else
	cdoAbort("Internal problem!");
      */
    }
  else
    {
      cdoAbort("Wrong operatorID!");
    }
  /*
  streamClose(streamID1);
  streamClose(streamID2);

  vlistDestroy(vlistID2);

  if ( array ) free(array);
  */
  cdoFinish();

  return (0);
}

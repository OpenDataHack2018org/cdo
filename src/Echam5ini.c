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
static const char strfiletype_ml[] = "Initial file spectral";
static const char strfiletype_sfc[] = "Initial file surface";


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


void inivar_ml(VAR *var, int gridtype, int zaxistype, int code, const char *name,
	       const char *longname, const char *units)
{
  static char func[] = "inivar_ml";
  
  var->gridtype  = gridtype;
  var->zaxistype = zaxistype;
  var->code      = code;
  var->name      = strdup(name);
  var->longname  = strdup(longname);
  var->units     = strdup(units);
}

void inivars_ml(VAR **vars)
{
  static char func[] = "inivars_ml";

  *vars = (VAR *) malloc((nvars_ml+1)*sizeof(VAR));

  inivar_ml(&(*vars)[0], GRID_GAUSSIAN, ZAXIS_HYBRID,  133, "Q",   "specific humidity", "kg/kg");
  inivar_ml(&(*vars)[1], GRID_SPECTRAL, ZAXIS_HYBRID,  138, "SVO", "vorticity", "1/s");
  inivar_ml(&(*vars)[2], GRID_SPECTRAL, ZAXIS_HYBRID,  155, "SD",  "divergence", "1/s");
  inivar_ml(&(*vars)[3], GRID_SPECTRAL, ZAXIS_HYBRID,  130, "STP", "temperature", "K");
  /* Don't change the order (lsp must be the last one)! */
  inivar_ml(&(*vars)[4], GRID_SPECTRAL, ZAXIS_SURFACE, 151, "LSP", "log surface pressure", "K");
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


#if  defined  (HAVE_LIBNETCDF)
static void read_echam5spec(int nc_file_id, VAR **vars)
{
  static char func[] = "read_echam5spec";
  int nc_dim_id, nc_var_id;
  size_t dimlen, nvals;
  size_t start[3];
  size_t count[3];
  int nlon, nlat, nlev, nlevp1, nvct, nsp, i, iv;
  int gridIDgp, gridIDsp, zaxisIDml, zaxisIDsfc;
  int gridtype, zaxistype;
  double *xvals, *yvals, *vct, *levs;

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
	      start[0] = 0;    start[1] = 0;  start[2] = i;
	      count[0] = nsp;  count[1] = 2;  count[2] = 1;
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

}
#endif


int read_echam5ini(const char *filename, VAR **vars)
{
  int filetype = 0;
#if  defined  (HAVE_LIBNETCDF)
  int nc_file_id;
  char inifiletype[256];
  size_t attlen;

  /* open file and check file type */
  nce(nc_open(filename, NC_NOWRITE, &nc_file_id));

  nce(nc_get_att_text(nc_file_id, NC_GLOBAL, "file_type", inifiletype));
  nce(nc_inq_attlen(nc_file_id, NC_GLOBAL, "file_type", &attlen));
  inifiletype[attlen] = 0;

  if ( strcmp(inifiletype, strfiletype_ml) == 0 )
    {
      filetype = filetype_ml;
      inivars_ml(vars);
    }
  else if ( strcmp(inifiletype, strfiletype_sfc) == 0 )
    {
      filetype = filetype_sfc;
    }
  else
    cdoAbort("Unsupported file type >%s<!", inifiletype);

  if ( filetype == filetype_ml )
    read_echam5spec(nc_file_id, vars);
      
  /*close input file */
  nce(nc_close(nc_file_id));

#else
  cdoAbort("netCDF library not available!");
#endif

  return (filetype);
}


void *Echam5ini(void *argument)
{
  static char func[] = "Echam5ini";
  int operatorID;
  int READ_E5INI, WRITE_E5INI;
  int streamID1, streamID2 = CDI_UNDEFID;
  int nrecs = 0;
  int recID, varID, levelID;
  int vlistID1, vlistID2;
  int nvars = 0;
  int iv, nlev;
  int filetype;
  int gridsize, nmiss;
  int taxisID, tsID;
  double *array = NULL;

  cdoInitialize(argument);

  READ_E5INI  = cdoOperatorAdd("read_e5ini",  0, 0, NULL);
  WRITE_E5INI = cdoOperatorAdd("write_e5ini", 0, 0, NULL);

  operatorID = cdoOperatorID();

  if ( operatorID == READ_E5INI )
    {
      VAR *vars = NULL;

      filetype = read_echam5ini(cdoStreamName(0), &vars);

      if ( vars == NULL ) cdoAbort("Unsupported file type!");
      
      if ( filetype == filetype_ml )
	nvars = nvars_ml+1;
      else if ( filetype == filetype_sfc )
	nvars = nvars_sfc;
      else
	cdoAbort("Unsupported file type!");

      vlistID2 = vlistCreate();

      for ( iv = 0; iv < nvars; iv++ )
	{	  
	  varID = vlistDefVar(vlistID2, vars[iv].gridID, vars[iv].zaxisID, TIME_CONSTANT);
	  vlistDefVarCode(vlistID2, varID, vars[iv].code);
	  vlistDefVarName(vlistID2, varID, vars[iv].name);
	  vlistDefVarLongname(vlistID2, varID, vars[iv].longname);
          vlistDefVarUnits(vlistID2, varID, vars[iv].units);
	}

      taxisID = taxisCreate(TAXIS_ABSOLUTE);
      vlistDefTaxis(vlistID2, taxisID);

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
  else if ( operatorID == WRITE_E5INI )
    {
      cdoAbort("Operator not implemented!");

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

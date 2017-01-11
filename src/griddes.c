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

#if defined(HAVE_LIBNETCDF)
#include "netcdf.h"
#endif

#include <ctype.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <limits.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "cdi_uuid.h"
#include "grid.h"
#include "griddes.h"
#include "error.h"


#define MAX_LINE_LEN 65536


int grid_read_pingo(FILE *gfp, const char *dname);


void gridInit(griddes_t *grid)
{
  grid->mask          = NULL;
  grid->xvals         = NULL;
  grid->yvals         = NULL;
  grid->xbounds       = NULL;
  grid->ybounds       = NULL;
  grid->area          = NULL;
  grid->type          = CDI_UNDEFID;
  grid->size          = 0;
  grid->xsize         = 0;
  grid->ysize         = 0;
  grid->np            = 0;
  grid->lcomplex      = 1;
  grid->prec          = 0;
  grid->ntr           = 0;
  grid->nvertex       = 0;
  grid->genBounds     = false;

  grid->originLon     = 0;
  grid->originLat     = 0;
  grid->lonParY       = 0;
  grid->lat1          = 0;
  grid->lat2          = 0;
  grid->projflag      = 0;
  grid->scanflag      = 64;
  grid->def_originLon = false;
  grid->def_originLat = false;
  grid->def_lonParY   = false;
  grid->def_lat1      = false;
  grid->def_lat2      = false;

  grid->a             = 0;
  grid->lon_0         = 0;
  grid->lat_0         = 0;
  grid->lat_1         = 0;
  grid->lat_2         = 0;
  grid->def_lon_0     = false;
  grid->def_lat_0     = false;
  grid->def_lat_1     = false;
  grid->def_lat_2     = false;

  grid->def_xfirst    = false;
  grid->def_yfirst    = false;
  grid->def_xlast     = false;
  grid->def_ylast     = false;
  grid->def_xinc      = false;
  grid->def_yinc      = false;
  grid->xfirst        = 0;
  grid->yfirst        = 0;
  grid->xlast         = 0;
  grid->ylast         = 0;
  grid->xinc          = 0;
  grid->yinc          = 0;
  grid->nd            = 0;
  grid->ni            = 0;
  grid->ni2           = 0;
  grid->ni3           = 0;
  grid->number        = 0;
  grid->position      = 0;
  grid->uuid[0]       = 0;
  grid->path[0]       = 0;
  grid->xname[0]      = 0;
  grid->xlongname[0]  = 0;
  grid->xunits[0]     = 0;
  grid->yname[0]      = 0;
  grid->ylongname[0]  = 0;
  grid->yunits[0]     = 0;
  grid->xdimname[0]   = 0;
  grid->ydimname[0]   = 0;
  grid->vdimname[0]   = 0;
}


int getoptname(char *optname, const char *optstring, int nopt)
{
  int nerr = 0;
  const char *pname = optstring;
  const char *pend  = optstring;

  for ( int i = 0; i < nopt; i++ )
    {
      pend = strchr(pname, ',');
      if ( pend == NULL )
	break;
      else
	pname = pend + 1;
    }

  if ( pend )
    {
      pend = strchr(pname, ',');
      size_t namelen;
      if ( pend == NULL )
	namelen = strlen(pname);
      else
	namelen = pend - pname;

      memcpy(optname, pname, namelen);
      optname[namelen] = '\0';
    }
  else
    nerr = 1;

  return nerr;
}


int gridDefine(griddes_t grid)
{
  int gridID = CDI_UNDEFID;
  int i;

  switch ( grid.type )
    {
    case GRID_GENERIC:
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_PROJECTION:
      {
	if ( grid.size != 1 )
	  {
	    if ( grid.xsize == 0 ) Error("xsize undefined!");
	    if ( grid.ysize == 0 ) Error("ysize undefined!");
	  }

	if ( grid.size == 0 ) grid.size = (long)grid.xsize*grid.ysize;

	if ( grid.size != (long)grid.xsize*grid.ysize )
	  Error("Inconsistent grid declaration: xsize*ysize!=gridsize (xsize=%d ysize=%d gridsize=%d)",
		grid.xsize, grid.ysize, grid.size);

	if ( grid.size < 0 || grid.size > INT_MAX ) Error("grid size (%ld) out of bounds (0 - %d)!", grid.size, INT_MAX);

	gridID = gridCreate(grid.type, grid.size);

	if ( grid.xsize > 0 ) gridDefXsize(gridID, grid.xsize);
	if ( grid.ysize > 0 ) gridDefYsize(gridID, grid.ysize);
	if ( grid.np    > 0 ) gridDefNP(gridID, grid.np);

	gridDefPrec(gridID, grid.prec);

	if ( (grid.def_xfirst || grid.def_xlast || grid.def_xinc) && grid.xvals == NULL )
	  {
	    grid.xvals = (double*) Malloc(grid.xsize*sizeof(double));
	    gridGenXvals(grid.xsize, grid.xfirst, grid.xlast, grid.xinc, grid.xvals);

	    if ( grid.genBounds && grid.xbounds == NULL && grid.xsize > 1 )
	      {
		grid.nvertex = 2;
		grid.xbounds = (double*) Malloc(grid.xsize*grid.nvertex*sizeof(double));
		for ( i = 0; i < (int) grid.xsize-1; i++ )
		  {
		    grid.xbounds[2*i+1]   = 0.5*(grid.xvals[i] + grid.xvals[i+1]);
		    grid.xbounds[2*(i+1)] = 0.5*(grid.xvals[i] + grid.xvals[i+1]);
		  }
		grid.xbounds[0] = 2*grid.xvals[0] - grid.xbounds[1];
		grid.xbounds[2*grid.xsize-1] = 2*grid.xvals[grid.xsize-1] - grid.xbounds[2*(grid.xsize-1)];
	      }
	  }

	if ( (grid.def_yfirst || grid.def_ylast || grid.def_yinc) && grid.yvals == NULL )
	  {
	    if ( ! grid.def_ylast ) grid.ylast = grid.yfirst;
	    grid.yvals = (double*) Malloc(grid.ysize*sizeof(double));
	    gridGenYvals(grid.type, grid.ysize, grid.yfirst, grid.ylast, grid.yinc, grid.yvals);

	    if ( grid.genBounds && grid.ybounds == NULL && grid.ysize > 1 )
	      {
		grid.nvertex = 2;
		grid.ybounds = (double*) Malloc(grid.ysize*grid.nvertex*sizeof(double));
		for ( i = 0; i < (int) grid.ysize-1; i++ )
		  {
		    grid.ybounds[2*i+1]   = 0.5*(grid.yvals[i] + grid.yvals[i+1]);
		    grid.ybounds[2*(i+1)] = 0.5*(grid.yvals[i] + grid.yvals[i+1]);
		  }

		if ( grid.yvals[0] > grid.yvals[grid.ysize-1] )
		  {
		    grid.ybounds[0] = 90;
		    grid.ybounds[grid.ysize*grid.nvertex-1] = -90;
		  }
		else
		  {
		    grid.ybounds[0] = -90;
		    grid.ybounds[grid.ysize*grid.nvertex-1] = 90;
		  }
	      }
	  }

	if ( grid.xvals )
	  {
	    gridDefXvals(gridID, grid.xvals);
	    Free(grid.xvals);
	  }

	if ( grid.yvals )
	  {
	    gridDefYvals(gridID, grid.yvals);
	    Free(grid.yvals);
	  }

	if ( grid.nvertex )
	  gridDefNvertex(gridID, grid.nvertex);

	if ( grid.xbounds )
	  {
	    gridDefXbounds(gridID, grid.xbounds);
	    Free(grid.xbounds);
	  }

	if ( grid.ybounds )
	  {
	    gridDefYbounds(gridID, grid.ybounds);
	    Free(grid.ybounds);
	  }

	if ( grid.mask )
	  {
	    gridDefMask(gridID, grid.mask);
	    Free(grid.mask);
	  }

	break;
      }
    case GRID_CURVILINEAR:
    case GRID_UNSTRUCTURED:
      {
	if ( grid.size == 0 )
	  {
	    if ( grid.type == GRID_CURVILINEAR )
	      grid.size = grid.xsize*grid.ysize;
	    else
	      grid.size = grid.xsize;
	  }

	gridID = gridCreate(grid.type, grid.size);

	gridDefPrec(gridID, grid.prec);

	if ( grid.type == GRID_CURVILINEAR )
	  {
	    if ( grid.xsize == 0 ) Error("xsize undefined!");
	    if ( grid.ysize == 0 ) Error("ysize undefined!");
	    gridDefXsize(gridID, grid.xsize);
	    gridDefYsize(gridID, grid.ysize);
	  }
	else
	  {
	    if ( grid.nvertex > 0 ) gridDefNvertex(gridID, grid.nvertex);
	    if ( grid.number > 0 )
	      {
		gridDefNumber(gridID, grid.number);
		if ( grid.position >= 0 ) gridDefPosition(gridID, grid.position);
	      }
	    if ( *grid.path ) gridDefReference(gridID, grid.path);
	  }

	if ( grid.xvals )
	  {
	    gridDefXvals(gridID, grid.xvals);
	    Free(grid.xvals);
	  }

	if ( grid.yvals )
	  {
	    gridDefYvals(gridID, grid.yvals);
	    Free(grid.yvals);
	  }

	if ( grid.area )
	  {
	    gridDefArea(gridID, grid.area);
	    Free(grid.area);
	  }

	if ( grid.xbounds )
	  {
	    gridDefXbounds(gridID, grid.xbounds);
	    Free(grid.xbounds);
	  }

	if ( grid.ybounds )
	  {
	    gridDefYbounds(gridID, grid.ybounds);
	    Free(grid.ybounds);
	  }

	if ( grid.mask )
	  {
	    gridDefMask(gridID, grid.mask);
	    Free(grid.mask);
	  }

	break;
      }
    case GRID_LCC:
      {
	if ( grid.xsize == 0 ) Error("xsize undefined!");
	if ( grid.ysize == 0 ) Error("ysize undefined!");

	if ( grid.size == 0 ) grid.size = grid.xsize*grid.ysize;

	gridID = gridCreate(grid.type, grid.size);

	gridDefPrec(gridID, grid.prec);

	gridDefXsize(gridID, grid.xsize);
	gridDefYsize(gridID, grid.ysize);

	if ( grid.def_originLon == false ) Error("originLon undefined!");
	if ( grid.def_originLat == false ) Error("originLat undefined!");
	if ( grid.def_lonParY   == false ) Error("lonParY undefined!");
	if ( grid.def_lat1      == false ) Error("lat1 undefined!");
	if ( grid.def_lat2      == false ) Error("lat2 undefined!");
	if ( grid.def_xinc      == false ) Error("xinc undefined!");
	if ( grid.def_yinc      == false ) Error("yinc undefined!");

	gridDefParamLCC(gridID, grid.originLon, grid.originLat, grid.lonParY,
		   grid.lat1, grid.lat2, grid.xinc, grid.yinc, grid.projflag, grid.scanflag);

	if ( grid.mask )
	  {
	    gridDefMask(gridID, grid.mask);
	    Free(grid.mask);
	  }

	break;
      }
    case GRID_SPECTRAL:
      {
	if ( grid.ntr == 0 )
	  Error("truncation undefined!");
	if ( grid.size == 0 )
	  grid.size = (grid.ntr+1) * (grid.ntr+2);

	gridID = gridCreate(grid.type, grid.size);

	gridDefPrec(gridID, grid.prec);

	gridDefTrunc(gridID, grid.ntr);

	gridDefComplexPacking(gridID, grid.lcomplex);

	break;
      }
    case GRID_GME:
      {
	if ( grid.nd   == 0 ) Error("nd undefined!");
	if ( grid.ni   == 0 ) Error("ni undefined!");
	if ( grid.size == 0 ) Error("size undefined!");

	gridID = gridCreate(grid.type, grid.size);

	gridDefPrec(gridID, grid.prec);

	gridDefParamGME(gridID, grid.nd, grid.ni, grid.ni2, grid.ni3);
	
	if ( grid.mask )
	  {
	    gridDefMask(gridID, grid.mask);
	    Free(grid.mask);
	  }

	break;
      }
    default:
      {
	if ( grid.type == -1 )
	  Error("Undefined grid type!");
	else
	  Error("Unsupported grid type: %s", gridNamePtr(grid.type));

	break;
      }
    }

  if ( grid.uuid[0] )      gridDefUUID(gridID, grid.uuid);

  if ( grid.xname[0]     ) cdiGridDefKeyStr(gridID, CDI_KEY_XNAME,     strlen(grid.xname)+1, grid.xname);
  if ( grid.xlongname[0] ) cdiGridDefKeyStr(gridID, CDI_KEY_XLONGNAME, strlen(grid.xlongname)+1, grid.xlongname);
  if ( grid.xunits[0]    ) cdiGridDefKeyStr(gridID, CDI_KEY_XUNITS,    strlen(grid.xunits)+1, grid.xunits);
  if ( grid.yname[0]     ) cdiGridDefKeyStr(gridID, CDI_KEY_YNAME,     strlen(grid.yname)+1, grid.yname);
  if ( grid.ylongname[0] ) cdiGridDefKeyStr(gridID, CDI_KEY_YLONGNAME, strlen(grid.ylongname)+1, grid.ylongname);
  if ( grid.yunits[0]    ) cdiGridDefKeyStr(gridID, CDI_KEY_YUNITS,    strlen(grid.yunits)+1, grid.yunits);
  if ( grid.xdimname[0]  ) cdiGridDefKeyStr(gridID, CDI_KEY_XDIMNAME,  strlen(grid.xdimname)+1, grid.xdimname);
  if ( grid.ydimname[0]  ) cdiGridDefKeyStr(gridID, CDI_KEY_YDIMNAME,  strlen(grid.ydimname)+1, grid.ydimname);
  if ( grid.vdimname[0]  ) cdiGridDefKeyStr(gridID, CDI_KEY_VDIMNAME,  strlen(grid.vdimname)+1, grid.vdimname);

  return gridID;
}

static
char *skipSeparator(char *pline)
{
  while ( isspace((int) *pline) ) pline++;
  if ( *pline == '=' || *pline == ':' ) pline++;
  while ( isspace((int) *pline) ) pline++;

  return pline;
}

static
char *read_att_name(char *pline, int len, char **attname, int *attlen)
{
  pline += len;
  if ( *pline == '_' )
    {
      pline++;
      *attlen = atol(pline);
      while ( isdigit(*pline) ) pline++;
    }
  else *attlen = 1;

  pline = skipSeparator(pline);
  *attname = pline;
  while ( *pline != ' ' && *pline != '=' && *pline != 0 ) pline++;
  char *endname = pline;
  pline = skipSeparator(pline);
  *endname = 0;

  return pline;
}

static
double read_value(const char *filename, const char *name, const char *pline)
{
  char *endptr;
  double val = strtod(pline, &endptr);
  if ( pline == endptr )
    cdoAbort("Couldn't read value for %s (grid description file: %s)!", name, filename);
 
  return val;
}


void cdo_read_field(const char *name, char *pline, int size, double *field, int *lineno, FILE *fp, const char *dname)
{
  char line[MAX_LINE_LEN];
  double fval;
  char *endptr;
  for ( int i = 0; i < size; i++ )
    {
      endptr = pline;
      fval = strtod(pline, &endptr);
      if ( pline == endptr )
        {
          (*lineno)++;
          if ( ! readline(fp, line, MAX_LINE_LEN) )
            cdoAbort("Incomplete command: >%s< (line: %d file: %s)", name, *lineno, dname);
          pline = line;
          fval = strtod(pline, &endptr);
        }
      field[i] = fval;
      pline = endptr;
    }
}


int gridFromFile(FILE *gfp, const char *dname)
{
  char line[MAX_LINE_LEN], *pline;
  int size;
  
  griddes_t grid;
  gridInit(&grid);

  int lineno = 0;
  while ( readline(gfp, line, MAX_LINE_LEN) )
    {
      lineno++;
      if ( line[0] == '#' ) continue;
      if ( line[0] == '\0' ) continue;
      size_t len = strlen(line);

      bool lerror = false;
      for ( size_t i = 0; i < len; ++i )
	if ( !(line[i] == 9 || (line[i] > 31 && line[i] < 127)) )
	  {
	    lerror = true;
	    line[i] = '#';
	  }
      if ( lerror ) cdoAbort("Grid description file >%s< contains illegal characters (line: %s)!", dname, line);

      pline = line;
      while ( isspace((int) *pline) ) pline++;
      if ( pline[0] == '\0' ) continue;
      if ( cmpstrlen(pline, "gridtype", len) == 0 )
	{
	  pline = skipSeparator(pline + len);
	  if ( cmpstrlen(pline, "lonlat", len)  == 0 ||
	       cmpstrlen(pline, "latlon", len)  == 0 )
	    {
	      grid.type = GRID_LONLAT;
	      grid.nvertex = 2;
	    }
	  else if ( cmpstrlen(pline, "gaussian", len)  == 0 )
	    {
	      grid.type = GRID_GAUSSIAN;
	      grid.nvertex = 2;
	    }
	  else if ( cmpstrlen(pline, "curvilinear", len)  == 0 )
	    {
	      grid.type = GRID_CURVILINEAR;
	      grid.nvertex = 4;
	    }
	  else if ( cmpstrlen(pline, "spectral", len)  == 0 )
	    grid.type = GRID_SPECTRAL;
	  else if ( cmpstrlen(pline, "unstructured", len)  == 0 )
	    grid.type = GRID_UNSTRUCTURED;
	  else if ( cmpstrlen(pline, "cell", len)  == 0 )
	    grid.type = GRID_UNSTRUCTURED;
	  else if ( cmpstrlen(pline, "gme", len)  == 0 )
	    grid.type = GRID_GME;
	  else if ( cmpstrlen(pline, "lcc", len)  == 0 )
	    grid.type = GRID_LCC;
	  else if ( cmpstrlen(pline, "lambert", len)  == 0 )
	    grid.type = GRID_LCC;
	  else if ( cmpstrlen(pline, "projection", len)  == 0 )
	    grid.type = GRID_PROJECTION;
	  else if ( cmpstrlen(pline, "generic", len)  == 0 )
	    grid.type = GRID_GENERIC;
	  else
	    cdoAbort("Invalid grid name : %s (grid description file: %s)", pline, dname);
	}
      else if ( cmpstrlen(pline, "gridprec", len)  == 0 )
	{
	  grid.prec = atol(skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "gridsize", len)  == 0 )
	{
	  grid.size = atol(skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "truncation", len)  == 0 )
	{
	  grid.ntr = atoi(skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "np", len)  == 0 )
	{
	  grid.np = atoi(skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "complexpacking", len)  == 0 )
	{
	  grid.lcomplex = atoi(skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "xname", len)  == 0 )
	{
	  strcpy(grid.xname, skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "xlongname", len)  == 0 )
	{
	  strcpy(grid.xlongname, skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "xunits", len)  == 0 )
	{
	  strcpy(grid.xunits, skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "xdimname", len)  == 0 )
	{
	  strcpy(grid.xdimname, skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "yname", len)  == 0 )
	{
	  strcpy(grid.yname, skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "ylongname", len)  == 0 )
	{
	  strcpy(grid.ylongname, skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "yunits", len)  == 0 )
	{
	  strcpy(grid.yunits, skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "ydimname", len)  == 0 )
	{
	  strcpy(grid.ydimname, skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "vdimname", len)  == 0 )
	{
	  strcpy(grid.vdimname, skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "nvertex", len)  == 0 )
	{
	  grid.nvertex = atol(skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "ni", len)  == 0 )
	{
	  grid.ni = atol(skipSeparator(pline + len));
          grid.nd = 10;
	}
      else if ( cmpstrlen(pline, "position", len)  == 0 )
	{
	  grid.position = atol(skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "number", len)  == 0 )
	{
	  grid.number = atol(skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "path", len)  == 0 )
	{
	  strcpy(grid.path, skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "uuid", len)  == 0 )
	{
	  char uuidOfHGridStr[256];
	  strcpy(uuidOfHGridStr, skipSeparator(pline + len));
	  cdiStr2UUID(uuidOfHGridStr, grid.uuid);
	}
      else if ( cmpstrlen(pline, "xsize", len)  == 0 )
	{
	  grid.xsize = atol(skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "nlon", len)  == 0 )
	{
	  grid.xsize = atol(skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "ysize", len)  == 0 )
	{
	  grid.ysize = atol(skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "nlat", len)  == 0 )
	{
	  grid.ysize = atol(skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "xfirst", len)  == 0 )
	{
	  grid.xfirst = read_value(dname, "xfirst", skipSeparator(pline + len));
	  grid.def_xfirst = true;
	}
      else if ( cmpstrlen(pline, "lonfirst", len)  == 0 )
	{
	  grid.xfirst = read_value(dname, "lonfirst", skipSeparator(pline + len));
	  grid.def_xfirst = true;
	}
      else if ( cmpstrlen(pline, "yfirst", len)  == 0 )
	{
	  grid.yfirst = read_value(dname, "yfirst", skipSeparator(pline + len));
	  grid.def_yfirst = true;
	}
      else if ( cmpstrlen(pline, "latfirst", len)  == 0 )
	{
	  grid.yfirst = read_value(dname, "latfirst", skipSeparator(pline + len));
	  grid.def_yfirst = true;
	}
      else if ( cmpstrlen(pline, "xlast", len)  == 0 )
	{
	  grid.xlast = read_value(dname, "xlast", skipSeparator(pline + len));
	  grid.def_xlast = true;
	}
      else if ( cmpstrlen(pline, "lonlast", len)  == 0 )
	{
	  grid.xlast = read_value(dname, "lonlast", skipSeparator(pline + len));
	  grid.def_xlast = true;
	}
      else if ( cmpstrlen(pline, "ylast", len)  == 0 )
	{
	  grid.ylast = read_value(dname, "ylast", skipSeparator(pline + len));
	  grid.def_ylast = true;
	}
      else if ( cmpstrlen(pline, "latlast", len)  == 0 )
	{
	  grid.ylast = read_value(dname, "latlast", skipSeparator(pline + len));
	  grid.def_ylast = true;
	}
      else if ( cmpstrlen(pline, "xinc", len)  == 0 )
	{
	  grid.xinc = read_value(dname, "xinc", skipSeparator(pline + len));
	  grid.def_xinc = true;
	}
      else if ( cmpstrlen(pline, "loninc", len)  == 0 )
	{
	  grid.xinc = read_value(dname, "loninc", skipSeparator(pline + len));
	  grid.def_xinc = true;
	}
      else if ( cmpstrlen(pline, "yinc", len)  == 0 )
	{
	  grid.yinc = read_value(dname, "yinc", skipSeparator(pline + len));
	  grid.def_yinc = true;
	}
      else if ( cmpstrlen(pline, "latinc", len)  == 0 )
	{
	  grid.yinc = read_value(dname, "latinc", skipSeparator(pline + len));
	  grid.def_yinc = true;
	}
      else if ( cmpstrlen(pline, "originLon", len)  == 0 )
	{
	  grid.originLon = read_value(dname, "originLon", skipSeparator(pline + len));
	  grid.def_originLon = true;
	}
      else if ( cmpstrlen(pline, "originLat", len)  == 0 )
	{
	  grid.originLat = read_value(dname, "originLat", skipSeparator(pline + len));
	  grid.def_originLat = true;
	}
      else if ( cmpstrlen(pline, "lonParY", len)  == 0 )
	{
	  grid.lonParY = read_value(dname, "lonParY", skipSeparator(pline + len));
	  grid.def_lonParY = true;
	}
      else if ( cmpstrlen(pline, "lat1", len)  == 0 )
	{
	  grid.lat1 = read_value(dname, "lat1", skipSeparator(pline + len));
	  grid.def_lat1 = true;
	}
      else if ( cmpstrlen(pline, "lat2", len)  == 0 )
	{
	  grid.lat2 = read_value(dname, "lat2", skipSeparator(pline + len));
	  grid.def_lat2 = true;
	}
      else if ( cmpstrlen(pline, "projection", len)  == 0 )
	{
	  pline = skipSeparator(pline + len);
	  if      ( cmpstrlen(pline, "north", len) == 0 )
	    {
	      grid.projflag = 0;
	      grid.scanflag = 64;
	    }
	  else if ( cmpstrlen(pline, "south", len) == 0 )
	    {
	      grid.projflag = 128;
	      grid.scanflag = 64;
	    }
	  else
	    cdoAbort("Invalid projection : %s (grid description file: %s)", pline, dname);
	}
      else if ( cmpstrlen(pline, "lon_0", len)  == 0 )
	{
	  grid.lon_0 = read_value(dname, "lon_0", skipSeparator(pline + len));
	  grid.def_lon_0 = true;
	}
      else if ( cmpstrlen(pline, "lat_0", len)  == 0 )
	{
	  grid.lat_0 = read_value(dname, "lat_0", skipSeparator(pline + len));
	  grid.def_lat_0 = true;
	}
      else if ( cmpstrlen(pline, "lat_1", len)  == 0 )
	{
	  grid.lat_1 = read_value(dname, "lat_1", skipSeparator(pline + len));
	  grid.def_lat_1 = true;
	}
      else if ( cmpstrlen(pline, "lat_2", len)  == 0 )
	{
	  grid.lat_2 = read_value(dname, "lat_2", skipSeparator(pline + len));
	  grid.def_lat_2 = true;
	}
      else if ( cmpstrlen(pline, "a", len)  == 0 )
	{
	  grid.a = read_value(dname, "a", skipSeparator(pline + len));
	}
      else if ( cmpstrlen(pline, "gridlatlon", len)  == 0 )
	{
	  double flat = 0, flon = 0;
	  if ( grid.size == 0 ) grid.size = grid.xsize * grid.ysize;
	  
	  grid.xvals = (double*) Malloc(grid.size*sizeof(double));
	  grid.yvals = (double*) Malloc(grid.size*sizeof(double));
	  for ( int i = 0; i < (int) grid.size; i++ )
	    {
              lineno++;
	      if ( ! readline(gfp, line, MAX_LINE_LEN) )
		cdoAbort("Incomplete command: >gridlatlon< (line: %d file: %s)", lineno, dname);

	      sscanf(line, "%lg %lg", &flat, &flon);
	      grid.yvals[i] = flat;
	      grid.xvals[i] = flon;
	    }
	}
      else if ( cmpstrlen(pline, "mask", len)  == 0 )
	{
	  size = grid.size;

	  if ( size > 0 )
	    {
	      long count = 0;
	      pline = skipSeparator(pline + len);
	      grid.mask = (int*) Malloc(size*sizeof(int));

	      for ( int i = 0; i < size; i++ )
		{
		  char *endptr = pline;
		  long lval = strtol(pline, &endptr, 10);
		  if ( pline == endptr )
		    {
                      lineno++;
		      if ( ! readline(gfp, line, MAX_LINE_LEN) )
			cdoAbort("Incomplete command: >mask< (line: %d file: %s)", lineno, dname);

		      pline = line;
		      lval = strtol(pline, &endptr, 10);
		    }
		  grid.mask[i] = (int)lval;
		  if ( grid.mask[i] == 1 ) count++;
		  pline = endptr;
		}

	      if ( count == size )
		{
		  Free(grid.mask);
		  grid.mask = NULL;
		}
	    }
	  else
	    cdoAbort("gridsize undefined (grid description file: %s)!", dname);
	}
      else if ( cmpstrlen(pline, "xvals", len)  == 0 )
	{
	  if ( grid.type == GRID_CURVILINEAR || grid.type == GRID_UNSTRUCTURED )
	    size = grid.size;
	  else
	    size = grid.xsize;

	  if ( size == 0 ) cdoAbort("xsize or gridsize undefined (grid description file: %s)!", dname);

          grid.xvals = (double*) Malloc(size*sizeof(double));
          pline = skipSeparator(pline + len);
          cdo_read_field("xvals", pline, size, grid.xvals, &lineno, gfp, dname);
	}
      else if ( cmpstrlen(pline, "yvals", len)  == 0 )
	{
	  if ( grid.type == GRID_CURVILINEAR || grid.type == GRID_UNSTRUCTURED )
	    size = grid.size;
	  else
	    size = grid.ysize;

	  if ( size == 0 ) cdoAbort("ysize or gridsize undefined (grid description file: %s)!", dname);

          grid.yvals = (double*) Malloc(size*sizeof(double));
          pline = skipSeparator(pline + len);
          cdo_read_field("yvals", pline, size, grid.yvals, &lineno, gfp, dname);
	}
      else if ( cmpstrlen(pline, "xbounds", len)  == 0 )
	{
	  if ( grid.nvertex == 0 )
	    {
	      if ( grid.type == GRID_LONLAT      ) grid.nvertex = 2;
	      if ( grid.type == GRID_CURVILINEAR ) grid.nvertex = 4;
	    }

	  if ( grid.type == GRID_CURVILINEAR || grid.type == GRID_UNSTRUCTURED )
	    size = grid.size;
	  else
	    size = grid.xsize;

          if ( size         == 0 ) cdoAbort("xsize or gridsize undefined (file: %s)!", dname);
          if ( grid.nvertex == 0 ) cdoAbort("nvertex undefined (file: %s)!", dname);

          grid.xbounds = (double*) Malloc(size*grid.nvertex*sizeof(double));
          pline = skipSeparator(pline + len);
          cdo_read_field("xbounds", pline, size*grid.nvertex, grid.xbounds, &lineno, gfp, dname);
	}
      else if ( cmpstrlen(pline, "ybounds", len)  == 0 )
	{
	  if ( grid.nvertex == 0 )
	    {
	      if ( grid.type == GRID_LONLAT      ) grid.nvertex = 2;
	      if ( grid.type == GRID_CURVILINEAR ) grid.nvertex = 4;
	    }

	  if ( grid.type == GRID_CURVILINEAR || grid.type == GRID_UNSTRUCTURED )
	    size = grid.size;
	  else
	    size = grid.ysize;

          if ( grid.ysize   == 0 ) cdoAbort("ysize or gridsize undefined (grid description file: %s)!", dname);
          if ( grid.nvertex == 0 ) cdoAbort("nvertex undefined!", dname);

          grid.ybounds = (double*) Malloc(size*grid.nvertex*sizeof(double));
          pline = skipSeparator(pline + len);
          cdo_read_field("ybounds", pline, size*grid.nvertex, grid.ybounds, &lineno, gfp, dname);
	}
      else if ( cmpstrlen(pline, "ATTR_TXT", len)  == 0 )
	{
          break;
	}
      else if ( cmpstrlen(pline, "ATTR_INT", len)  == 0 )
	{
          break;
	}
      else if ( cmpstrlen(pline, "ATTR_FLT", len)  == 0 )
	{
          break;
	}
      else
	{
	  if ( grid.type != CDI_UNDEFID )
	    cdoAbort("Invalid grid command : >%s< (line: %d file: %s)", pline, lineno, dname);
	}
    }
  /*
  printf("gridtype %d\n", grid.type);
  printf("gridsize %d\n", grid.size);
  printf("xsize %d\n", grid.xsize);
  printf("ysize %d\n", grid.ysize);
  */
  int gridID = (grid.type == CDI_UNDEFID ) ? -1 : gridDefine(grid);

  // define attributes
  int attlen;
  char *attname;
  do
    {
      if ( line[0] == '#' ) continue;
      if ( line[0] == '\0' ) continue;
      size_t len = strlen(line);

      bool lerror = false;
      for ( size_t i = 0; i < len; ++i )
	if ( !(line[i] == 9 || (line[i] > 31 && line[i] < 127)) )
	  {
	    lerror = true;
	    line[i] = '#';
	  }
      if ( lerror ) cdoAbort("Grid description file >%s< contains illegal characters (line: %s)!", dname, line);

      pline = line;
      while ( isspace((int) *pline) ) pline++;
      if ( pline[0] == '\0' ) continue;

      if ( cmpstrlen(pline, "ATTR_TXT", len)  == 0 )
	{
          pline = read_att_name(pline, len, &attname, &attlen);

          if ( *pline == '"' ) pline++;
          char *atttxt = pline;
          while ( *pline != 0 && *pline != '"' ) pline++;
          if ( *pline == '"' ) *pline = 0;

          if ( strcmp(attname, "grid_mapping_name") == 0 )
            cdiGridDefKeyStr(gridID, CDI_KEY_MAPPING, (int)strlen(atttxt)+1, atttxt);

          cdiDefAttTxt(gridID, CDI_GLOBAL, attname, (int)strlen(atttxt), atttxt);
	}
      else if ( cmpstrlen(pline, "ATTR_INT", len)  == 0 )
	{
          pline = read_att_name(pline, len, &attname, &attlen);

          int *attint = (int*) Malloc(attlen*sizeof(int));
          double *attflt = (double*) Malloc(attlen*sizeof(double));
          cdo_read_field("attint", pline, attlen, attflt, &lineno, gfp, dname);
          for ( int i = 0; i < attlen; ++i ) attint[i] = (int)lround(attflt[i]);
          cdiDefAttInt(gridID, CDI_GLOBAL, attname, CDI_DATATYPE_INT32, attlen, attint);
          free(attint);
          free(attflt);
	}
      else if ( cmpstrlen(pline, "ATTR_FLT", len)  == 0 )
	{
          pline = read_att_name(pline, len, &attname, &attlen);

          double *attflt = (double*) Malloc(attlen*sizeof(double));
          cdo_read_field("attflt", pline, attlen, attflt, &lineno, gfp, dname);
          cdiDefAttFlt(gridID, CDI_GLOBAL, attname, CDI_DATATYPE_FLT64, attlen, attflt);
          free(attflt);
	}
      else
	{
          cdoAbort("Invalid grid command : >%s< (line: %d file: %s)", pline, lineno, dname);
	}

      lineno++;
    }
  while ( readline(gfp, line, MAX_LINE_LEN) );

  return gridID;
}


int cdoDefineGrid(const char *gridfile)
{
  int gridID = -1;
  size_t len;
  bool isreg = false;
  bool lalloc = false;

  char *filename = expand_filename(gridfile);
  if ( filename )
    lalloc = true;
  else
    filename =  (char *) gridfile;

  int fileno = open(filename, O_RDONLY);
  if ( fileno >= 0 )
    {
      struct stat filestat;
      if ( fstat(fileno, &filestat) == 0 )
	isreg = S_ISREG(filestat.st_mode);
    }

  if ( fileno == -1 || !isreg )
    {
      if ( isreg ) close(fileno);

      gridID = grid_from_name(gridfile);

      if ( gridID == -1 ) cdoAbort("Open failed on %s!", gridfile);
    }
  else
    {
      char buffer[4];
      if ( read(fileno, buffer, 4) != 4 )
	SysError("Read grid from %s failed!", filename);

      close(fileno);

      if ( cmpstrlen(buffer, "CDF", len) == 0 )
	{
	  if ( cdoDebug ) cdoPrint("Grid from NetCDF file");
	  gridID = gridFromNCfile(filename);
	}

      if ( gridID == -1 )
	{
	  if ( cmpstrlen(buffer+1, "HDF", len) == 0 )
	    {
	      if ( cdoDebug ) cdoPrint("Grid from HDF5 file");
	      gridID = gridFromH5file(filename);
	    }
	}

      if ( gridID == -1 )
	{
	  if ( cmpstrlen(buffer+1, "HDF", len) == 0 )
	    {
	      if ( cdoDebug ) cdoPrint("Grid from NetCDF4 file");
	      gridID = gridFromNCfile(filename);
	    }
	}

      if ( gridID == -1 )
	{
	  if ( cdoDebug ) cdoPrint("Grid from CDI file");
	  openLock();
	  int streamID = streamOpenRead(filename);
	  openUnlock();
	  if ( streamID >= 0 )
	    {
	      int vlistID = streamInqVlist(streamID);
	      gridID  = vlistGrid(vlistID, 0);
	      streamClose(streamID);
	    }
	}

      if ( gridID == -1 )
	{
	  if ( cdoDebug ) cdoPrint("grid from ASCII file");
	  FILE *gfp = fopen(filename, "r");
	  //size_t buffersize = 20*1024*1024;
	  //char *buffer = (char*) Malloc(buffersize);
	  //setvbuf(gfp, buffer, _IOFBF, buffersize);
	  gridID = gridFromFile(gfp, filename);
	  fclose(gfp);
	  //free(buffer);
	}

      if ( gridID == -1 )
	{
	  if ( cdoDebug ) cdoPrint("grid from PINGO file");
	  FILE *gfp = fopen(filename, "r");
	  gridID = grid_read_pingo(gfp, filename);
	  fclose(gfp);
	}

      if ( gridID == -1 ) cdoAbort("Invalid grid description file %s!", filename);
    }

  if ( lalloc ) Free(filename);

  return gridID;
}


void defineGrid(const char *gridarg)
{
  char gridfile[4096];
  int nfile = 0;

  while ( getoptname(gridfile, gridarg, nfile++) == 0 )
    {      
      (void) cdoDefineGrid(gridfile);
    }
}

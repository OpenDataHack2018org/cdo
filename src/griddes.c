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
#  include "config.h"
#endif

#if defined(HAVE_LIBNETCDF)
#  include "netcdf.h"
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


#define UNDEFID -1

#define MAX_LINE_LEN 65536


void gridInit(griddes_t *grid)
{
  grid->mask          = NULL;
  grid->xvals         = NULL;
  grid->yvals         = NULL;
  grid->xbounds       = NULL;
  grid->ybounds       = NULL;
  grid->area          = NULL;
  grid->type          = UNDEFID;
  grid->size          = 0;
  grid->xsize         = 0;
  grid->ysize         = 0;
  grid->np            = 0;
  grid->lcomplex      = 1;
  grid->xpole         = 0;
  grid->ypole         = 0;
  grid->prec          = 0;
  grid->isRotated     = FALSE;
  grid->ntr           = 0;
  grid->nvertex       = 0;
  grid->genBounds     = FALSE;

  grid->originLon     = 0;
  grid->originLat     = 0;
  grid->lonParY       = 0;
  grid->lat1          = 0;
  grid->lat2          = 0;
  grid->projflag      = 0;
  grid->scanflag      = 64;
  grid->def_originLon = FALSE;
  grid->def_originLat = FALSE;
  grid->def_lonParY   = FALSE;
  grid->def_lat1      = FALSE;
  grid->def_lat2      = FALSE;

  grid->a             = 0;
  grid->lon_0         = 0;
  grid->lat_0         = 0;
  grid->lat_1         = 0;
  grid->lat_2         = 0;
  grid->def_lon_0     = FALSE;
  grid->def_lat_0     = FALSE;
  grid->def_lat_1     = FALSE;
  grid->def_lat_2     = FALSE;

  grid->def_xfirst    = FALSE;
  grid->def_yfirst    = FALSE;
  grid->def_xlast     = FALSE;
  grid->def_ylast     = FALSE;
  grid->def_xinc      = FALSE;
  grid->def_yinc      = FALSE;
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
  size_t namelen;
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
  int gridID = UNDEFID;
  int i;

  switch ( grid.type )
    {
    case GRID_GENERIC:
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_SINUSOIDAL:
    case GRID_LAEA:
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

	if ( grid.isRotated )
	  {
	    gridDefXpole(gridID, grid.xpole);
	    gridDefYpole(gridID, grid.ypole);
	    gridDefAngle(gridID, grid.angle);
	  }

	if ( grid.mask )
	  {
	    gridDefMask(gridID, grid.mask);
	    Free(grid.mask);
	  }

	if ( grid.type == GRID_LAEA )
	  {
	    if ( grid.a > 0 ) gridDefLaea(gridID, grid.a, grid.lon_0, grid.lat_0);
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

	if ( grid.def_originLon == FALSE ) Error("originLon undefined!");
	if ( grid.def_originLat == FALSE ) Error("originLat undefined!");
	if ( grid.def_lonParY   == FALSE ) Error("lonParY undefined!");
	if ( grid.def_lat1      == FALSE ) Error("lat1 undefined!");
	if ( grid.def_lat2      == FALSE ) Error("lat2 undefined!");
	if ( grid.def_xinc      == FALSE ) Error("xinc undefined!");
	if ( grid.def_yinc      == FALSE ) Error("yinc undefined!");

	gridDefLCC(gridID, grid.originLon, grid.originLat, grid.lonParY,
		   grid.lat1, grid.lat2, grid.xinc, grid.yinc, grid.projflag, grid.scanflag);

	if ( grid.mask )
	  {
	    gridDefMask(gridID, grid.mask);
	    Free(grid.mask);
	  }

	break;
      }
    case GRID_LCC2:
      {
	if ( grid.xsize == 0 ) Error("xsize undefined!");
	if ( grid.ysize == 0 ) Error("ysize undefined!");

	if ( grid.size == 0 ) grid.size = grid.xsize*grid.ysize;

	gridID = gridCreate(grid.type, grid.size);

	gridDefPrec(gridID, grid.prec);

	gridDefXsize(gridID, grid.xsize);
	gridDefYsize(gridID, grid.ysize);

	if ( grid.def_xfirst && grid.def_xinc && grid.xvals == NULL )
	  {
	    grid.xvals = (double*) Malloc(grid.xsize*sizeof(double));
	    for ( i = 0; i < grid.xsize; ++i )
	      grid.xvals[i] = grid.xfirst + i*grid.xinc;
	  }

	if ( grid.def_yfirst && grid.def_yinc && grid.yvals == NULL )
	  {
	    grid.yvals = (double*) Malloc(grid.ysize*sizeof(double));
	    for ( i = 0; i < grid.ysize; ++i )
	      grid.yvals[i] = grid.yfirst + i*grid.yinc;
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

	if ( grid.def_lon_0     == FALSE ) Error("lon_0 undefined!");
	if ( grid.def_lat_0     == FALSE ) Error("lat_0 undefined!");
	if ( grid.def_lat_1     == FALSE ) Error("lat_1 undefined!");
	if ( grid.def_lat_2     == FALSE ) grid.def_lat_2 = grid.def_lat_1;

	gridDefLcc2(gridID, grid.a, grid.lon_0, grid.lat_0, grid.lat_1, grid.lat_2);

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

	gridDefGMEnd(gridID, grid.nd);
	gridDefGMEni(gridID, grid.ni);
	gridDefGMEni2(gridID, grid.ni2);
	gridDefGMEni3(gridID, grid.ni3);
	
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

  if ( grid.xname[0]     ) gridDefXname(gridID, grid.xname);
  if ( grid.xlongname[0] ) gridDefXlongname(gridID, grid.xlongname);
  if ( grid.xunits[0]    ) gridDefXunits(gridID, grid.xunits);
  if ( grid.yname[0]     ) gridDefYname(gridID, grid.yname);
  if ( grid.ylongname[0] ) gridDefYlongname(gridID, grid.ylongname);
  if ( grid.yunits[0]    ) gridDefYunits(gridID, grid.yunits);
  if ( grid.xdimname[0]  ) cdiGridDefKeyStr(gridID, CDI_KEY_XDIMNAME, strlen(grid.xdimname)+1, grid.xdimname);
  if ( grid.ydimname[0]  ) cdiGridDefKeyStr(gridID, CDI_KEY_YDIMNAME, strlen(grid.ydimname)+1, grid.ydimname);
  if ( grid.vdimname[0]  ) cdiGridDefKeyStr(gridID, CDI_KEY_VDIMNAME, strlen(grid.vdimname)+1, grid.vdimname);

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


void fnmexp2(char *out, char *in1, const char *in2)
{
  const char *pos, *ch;
  char envv[20], *envr;
  int i,j;

  if ( *in1=='$' )
    {
      in1++;
      i = 0;
      while (*in1!='/' && *in1!='\0' && i<16)
	{
	  envv[i] = *in1;
	  i++; in1++;
	}
      envv[i] = '\0';
      envr = getenv(envv);
      if (envr)
	{
	  i = 0; j = 0;
	  while (*(envr+j))
	    {
	      *(out+i) = *(envr+j);
	      i++; j++;
	    }
	  while (*in1!='\0' && *in1!=' ' && *in1!='\n')
	    {
	      *(out+i) = *in1;
	      i++; in1++;
	    }
	  *(out+i) = '\0';
	}
      return;
    }
  ch = in2;
  pos=NULL;
  while (*ch!='\0' && *ch!=' ' && *ch!='\n')
    {
      if (*ch=='/') pos=ch;
      ch++;
    }
  if (pos) pos++;
  while (pos!=NULL && in2<pos)
    {
      *out = *in2;
      out++; in2++;
    }
  in1++;
  while (*in1!='\0' && *in1!=' ' && *in1!='\n')
    {
      *out = *in1;
      out++; in1++;
    }
  *out = '\0';
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

static
void read_field(const char *name, char *pline, int size, double *field, int *lineno, FILE *gfp, const char *dname)
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
          if ( ! readline(gfp, line, MAX_LINE_LEN) )
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
	  else if ( cmpstrlen(pline, "lcc2", len)  == 0 )
	    grid.type = GRID_LCC2;
	  else if ( cmpstrlen(pline, "lcc", len)  == 0 )
	    grid.type = GRID_LCC;
	  else if ( cmpstrlen(pline, "lambert", len)  == 0 )
	    grid.type = GRID_LCC;
	  else if ( cmpstrlen(pline, "sinusoidal", len)  == 0 )
	    grid.type = GRID_SINUSOIDAL;
	  else if ( cmpstrlen(pline, "laea", len)  == 0 )
	    grid.type = GRID_LAEA;
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
	  grid.def_xfirst = TRUE;
	}
      else if ( cmpstrlen(pline, "lonfirst", len)  == 0 )
	{
	  grid.xfirst = read_value(dname, "lonfirst", skipSeparator(pline + len));
	  grid.def_xfirst = TRUE;
	}
      else if ( cmpstrlen(pline, "yfirst", len)  == 0 )
	{
	  grid.yfirst = read_value(dname, "yfirst", skipSeparator(pline + len));
	  grid.def_yfirst = TRUE;
	}
      else if ( cmpstrlen(pline, "latfirst", len)  == 0 )
	{
	  grid.yfirst = read_value(dname, "latfirst", skipSeparator(pline + len));
	  grid.def_yfirst = TRUE;
	}
      else if ( cmpstrlen(pline, "xlast", len)  == 0 )
	{
	  grid.xlast = read_value(dname, "xlast", skipSeparator(pline + len));
	  grid.def_xlast = TRUE;
	}
      else if ( cmpstrlen(pline, "lonlast", len)  == 0 )
	{
	  grid.xlast = read_value(dname, "lonlast", skipSeparator(pline + len));
	  grid.def_xlast = TRUE;
	}
      else if ( cmpstrlen(pline, "ylast", len)  == 0 )
	{
	  grid.ylast = read_value(dname, "ylast", skipSeparator(pline + len));
	  grid.def_ylast = TRUE;
	}
      else if ( cmpstrlen(pline, "latlast", len)  == 0 )
	{
	  grid.ylast = read_value(dname, "latlast", skipSeparator(pline + len));
	  grid.def_ylast = TRUE;
	}
      else if ( cmpstrlen(pline, "xinc", len)  == 0 )
	{
	  grid.xinc = read_value(dname, "xinc", skipSeparator(pline + len));
	  grid.def_xinc = TRUE;
	}
      else if ( cmpstrlen(pline, "loninc", len)  == 0 )
	{
	  grid.xinc = read_value(dname, "loninc", skipSeparator(pline + len));
	  grid.def_xinc = TRUE;
	}
      else if ( cmpstrlen(pline, "yinc", len)  == 0 )
	{
	  grid.yinc = read_value(dname, "yinc", skipSeparator(pline + len));
	  grid.def_yinc = TRUE;
	}
      else if ( cmpstrlen(pline, "latinc", len)  == 0 )
	{
	  grid.yinc = read_value(dname, "latinc", skipSeparator(pline + len));
	  grid.def_yinc = TRUE;
	}
      else if ( cmpstrlen(pline, "originLon", len)  == 0 )
	{
	  grid.originLon = read_value(dname, "originLon", skipSeparator(pline + len));
	  grid.def_originLon = TRUE;
	}
      else if ( cmpstrlen(pline, "originLat", len)  == 0 )
	{
	  grid.originLat = read_value(dname, "originLat", skipSeparator(pline + len));
	  grid.def_originLat = TRUE;
	}
      else if ( cmpstrlen(pline, "lonParY", len)  == 0 )
	{
	  grid.lonParY = read_value(dname, "lonParY", skipSeparator(pline + len));
	  grid.def_lonParY = TRUE;
	}
      else if ( cmpstrlen(pline, "lat1", len)  == 0 )
	{
	  grid.lat1 = read_value(dname, "lat1", skipSeparator(pline + len));
	  grid.def_lat1 = TRUE;
	}
      else if ( cmpstrlen(pline, "lat2", len)  == 0 )
	{
	  grid.lat2 = read_value(dname, "lat2", skipSeparator(pline + len));
	  grid.def_lat2 = TRUE;
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
	  grid.def_lon_0 = TRUE;
	}
      else if ( cmpstrlen(pline, "lat_0", len)  == 0 )
	{
	  grid.lat_0 = read_value(dname, "lat_0", skipSeparator(pline + len));
	  grid.def_lat_0 = TRUE;
	}
      else if ( cmpstrlen(pline, "lat_1", len)  == 0 )
	{
	  grid.lat_1 = read_value(dname, "lat_1", skipSeparator(pline + len));
	  grid.def_lat_1 = TRUE;
	}
      else if ( cmpstrlen(pline, "lat_2", len)  == 0 )
	{
	  grid.lat_2 = read_value(dname, "lat_2", skipSeparator(pline + len));
	  grid.def_lat_2 = TRUE;
	}
      else if ( cmpstrlen(pline, "xnpole", len)  == 0 )
	{
	  grid.xpole = read_value(dname, "xnpole", skipSeparator(pline + len));
	  grid.isRotated = TRUE;
	}
      else if ( cmpstrlen(pline, "lonpole", len)  == 0 )
	{
	  grid.xpole = read_value(dname, "lonpole", skipSeparator(pline + len));
	  grid.isRotated = TRUE;
	}
      else if ( cmpstrlen(pline, "ynpole", len)  == 0 )
	{
	  grid.ypole = read_value(dname, "ynpole", skipSeparator(pline + len));
	  grid.isRotated = TRUE;
	}
      else if ( cmpstrlen(pline, "latpole", len)  == 0 )
	{
	  grid.ypole = read_value(dname, "latpole", skipSeparator(pline + len));
	  grid.isRotated = TRUE;
	}
      else if ( cmpstrlen(pline, "angle", len)  == 0 )
	{
	  grid.angle = read_value(dname, "angle", skipSeparator(pline + len));
	  grid.isRotated = TRUE;
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

	  if ( size > 0 )
	    {
	      grid.xvals = (double*) Malloc(size*sizeof(double));
	      pline = skipSeparator(pline + len);
              read_field("xvals", pline, size, grid.xvals, &lineno, gfp, dname);
	    }
	  else
	    cdoAbort("xsize or gridsize undefined (file: %s)!", dname);
	}
      else if ( cmpstrlen(pline, "yvals", len)  == 0 )
	{
	  if ( grid.type == GRID_CURVILINEAR || grid.type == GRID_UNSTRUCTURED )
	    size = grid.size;
	  else
	    size = grid.ysize;

	  if ( size > 0 )
	    {
	      grid.yvals = (double*) Malloc(size*sizeof(double));
	      pline = skipSeparator(pline + len);
              read_field("yvals", pline, size, grid.yvals, &lineno, gfp, dname);
	    }
	  else
	    cdoAbort("ysize or gridsize undefined (grid description file: %s)!", dname);
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

	  if ( size > 0 && grid.nvertex > 0 )
	    {	  
	      grid.xbounds = (double*) Malloc(size*grid.nvertex*sizeof(double));
	      pline = skipSeparator(pline + len);
              read_field("xbounds", pline, size*grid.nvertex, grid.xbounds, &lineno, gfp, dname);
	    }
	  else
	    {
	      if ( size         == 0 ) cdoAbort("xsize or gridsize undefined (file: %s)!", dname);
	      if ( grid.nvertex == 0 ) cdoAbort("nvertex undefined (file: %s)!", dname);
	    }
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

	  if ( size > 0 && grid.nvertex > 0 )
	    {	  
	      grid.ybounds = (double*) Malloc(size*grid.nvertex*sizeof(double));
	      pline = skipSeparator(pline + len);
              read_field("ybounds", pline, size*grid.nvertex, grid.ybounds, &lineno, gfp, dname);
	    }
	  else
	    {
	      if ( grid.ysize   == 0 ) cdoAbort("ysize or gridsize undefined (grid description file: %s)!", dname);
	      if ( grid.nvertex == 0 ) cdoAbort("nvertex undefined!", dname);
	    }
	}
      else
	{
	  if ( grid.type != UNDEFID )
	    cdoAbort("Invalid grid command : >%s< (line: %d file: %s)", pline, lineno, dname);
	}
    }
  /*
  printf("gridtype %d\n", grid.type);
  printf("gridsize %d\n", grid.size);
  printf("xsize %d\n", grid.xsize);
  printf("ysize %d\n", grid.ysize);
  */
  int gridID = (grid.type == UNDEFID ) ? -1 : gridDefine(grid);

  return gridID;
}


void skip_nondigit_lines(FILE *gfp)
{
  int c;

  if ( feof(gfp) ) return;

  while (1)
    {
      do
	c = fgetc(gfp);
      while ( (isspace(c) || c == ',') && c != EOF );

      if ( c == EOF || isdigit (c) || c == '.' || c == '+' || c == '-' ) break;
      else
	while ( c != '\n' && c != EOF )
	  c = fgetc(gfp);
    }

  ungetc(c, gfp);
}


int input_ival(FILE *gfp, int *ival)
{
  skip_nondigit_lines(gfp);

  if ( feof(gfp) ) return 0;

  *ival = 0;
  int read_items = fscanf(gfp, "%d", ival);

  return read_items;
}


int input_darray(FILE *gfp, int n_values, double *array)
{
  if ( n_values <= 0 ) return 0;

  int read_items = 0;
  for ( int i = 0; i < n_values; i++ )
    {
      skip_nondigit_lines(gfp);

      if ( feof(gfp) ) break;

      read_items += fscanf(gfp, "%lg", &array[i]);

      if ( feof(gfp) ) break;
    }

  return read_items;
}


int gridFromPingo(FILE *gfp, const char *dname)
{
  UNUSED(dname);
  int gridID = -1;
  int i;

  griddes_t grid;
  gridInit(&grid);

  int nlon, nlat;
  if ( ! input_ival(gfp, &nlon) ) return gridID;
  if ( ! input_ival(gfp, &nlat) ) return gridID;

  if ( nlon > 0 && nlon < 9999 && nlat > 0 && nlat < 9999 )
    {
      grid.xsize = nlon;
      grid.ysize = nlat;

      grid.xvals = (double*) Malloc(grid.xsize*sizeof(double));
      grid.yvals = (double*) Malloc(grid.ysize*sizeof(double));

      if ( ! input_ival(gfp, &nlon) ) return gridID;
      if ( nlon == 2 )
	{
	  if ( input_darray(gfp, 2, grid.xvals) != 2 ) return gridID;
	  grid.xvals[1] -= 360 * floor((grid.xvals[1] - grid.xvals[0]) / 360);

	  if ( grid.xsize > 1 )
	    if ( IS_EQUAL(grid.xvals[0], grid.xvals[1]) )
	      grid.xvals[1] += 360;

	  for ( i = 0; i < (int)grid.xsize; i++ )
	    grid.xvals[i] = grid.xvals[0] + i*(grid.xvals[1] - grid.xvals[0]);
	}
      else if ( nlon == (int)grid.xsize )
	{
	  if ( input_darray(gfp, nlon, grid.xvals) != nlon ) return gridID;
	  for ( i = 0; i < nlon - 1; i++ )
	    if ( grid.xvals[i+1] <= grid.xvals[i] ) break;

	  for ( i++; i < nlon; i++ )
	    {
	      grid.xvals[i] += 360;
	      if ( i < nlon - 1 && grid.xvals[i+1] + 360 <= grid.xvals[i] )
		{
		  Message("Longitudes are not in ascending order!");
		  return gridID;
		}
	    }
	}
      else
	return gridID;

      if ( ! input_ival(gfp, &nlat) ) return gridID;
      if ( nlat == 2 )
	{
	  if ( input_darray(gfp, 2, grid.yvals) != 2 ) return gridID;
	  for ( i = 0; i < (int)grid.ysize; i++ )
	    grid.yvals[i] = grid.yvals[0] + i*(grid.yvals[1] - grid.yvals[0]);
	}
      else if ( nlat == (int)grid.ysize )
	{
	  if ( input_darray(gfp, nlat, grid.yvals) != nlat ) return gridID;
	}
      else
	return gridID;

      if ( grid.yvals[0]      >  90.001  || 
	   grid.yvals[nlat-1] >  90.001  || 
	   grid.yvals[0]      < -90.001  || 
	   grid.yvals[nlat-1] < -90.001 )
	{
	  Message("Latitudes must be between 90 and -90!");
	  return gridID;
	}

      for ( i = 0; i < nlat - 1; i++ )
	if ( IS_EQUAL(grid.yvals[i+1], grid.yvals[i]) || (i < nlat - 2 &&
	    ((grid.yvals[i+1] > grid.yvals[i]) != (grid.yvals[i+2] > grid.yvals[i+1]))) )
	  {
	    Message("Latitudes must be in descending or ascending order!");
	    return gridID;
	  }
		    
      bool lgauss = false;
      if ( nlat > 2 ) /* check if gaussian */
	{
	  double *yvals, *yw;
	  yvals = (double*) Malloc(grid.ysize*sizeof(double));
	  yw    = (double*) Malloc(grid.ysize*sizeof(double));
	  gaussaw(yvals, yw, grid.ysize);
	  Free(yw);
	  for ( i = 0; i < (int) grid.ysize; i++ )
	    yvals[i] = asin(yvals[i])*RAD2DEG;

	  for ( i = 0; i < (int) grid.ysize; i++ )
	    if ( fabs(yvals[i] - grid.yvals[i]) > ((yvals[0] - yvals[1])/500) ) break;
		      
	  if ( i == (int) grid.ysize ) lgauss = true;

	  Free(yvals);
	}

      if ( lgauss )
	grid.type = GRID_GAUSSIAN;
      else
	grid.type = GRID_LONLAT;
    }
  
  if ( grid.type != UNDEFID ) gridID = gridDefine(grid);

  return gridID;
}


int nfc2nlat(int nfc, int ntr)
{
  int nlat = nfc / (ntr+1);
  nlat /= 2;

  return nlat;
}


int nlat2ntr(int nlat)
{
  int ntr = (nlat*2 - 1) / 3;

  return ntr;
}


int nlat2ntr_linear(int nlat)
{
  int ntr = (nlat*2 - 1) / 2;

  return ntr;
}


int ntr2nlat(int ntr)
{
  int nlat = (int)lround((ntr*3.+1.)/2.);
  if ( (nlat % 2) > 0 )
    {
      nlat  = nlat + 1;
      /*
      int nlat2 = (int)lround(((ntr+1)*3.+1.)/2.);
      if ( nlat == nlat2 )
	Error("Computation of latitudes failed for truncation %d", ntr);
      */
    }

  return nlat;
}


int ntr2nlat_linear(int ntr)
{
  int nlat = (int)lround((ntr*2.+1.)/2.);
  if ( (nlat % 2) > 0 )
    {
      nlat  = nlat + 1;
      /*
      int nlat2 = (int)lround(((ntr+1)*2.+1.)/2.);
      if ( nlat == nlat2 )
	Error("Computation of latitudes failed for truncation %d", ntr);
      */
    }

  return nlat;
}


int compNlon(int nlat)
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

      if ( nlon > 9999 )
	{
	  nlon = 2 * nlat;
	  fprintf(stderr, "FFT does not work with len %d!\n", nlon);
	  break;
	}
    }

  return nlon;
}

static
void gen_grid_lonlat(griddes_t *grid, const char *pline, double inc, double lon1, double lon2, double lat1, double lat2)
{
  int gridtype = GRID_LONLAT;
  bool lbounds = true;

  if ( *pline != 0 && (*pline == '+' || *pline == '-') && (isdigit((int) *(pline+1)) || ispunct((int) *(pline+1))) )
    {
      char *endptr = (char *) pline;
      double off = strtod(pline, &endptr);
      pline = endptr;
      
      lon1 -= off;
      lon2 += off;
      lat1 -= off;
      lat2 += off;
      if ( lat1 < -90 ) lat1 = -90;
      if ( lat2 >  90 ) lat2 =  90;
    }

  if ( *pline != 0 )
    {
      if ( *pline == '_' ) pline++;
      else return;

      if ( *pline == 0 ) return;

      if ( ! isdigit((int) *pline) && !ispunct((int) *pline) ) return;

      char *endptr = (char *) pline;
      inc = strtod(pline, &endptr);
      if ( *endptr != 0 )
        {
          pline = endptr;
          if ( *pline == '_' ) pline++;
          else return;
          
          if ( *pline == 0 ) return;
          if ( *pline == 'c' )
            {
              gridtype = GRID_CURVILINEAR;
              pline++;
              if ( *pline == '0' )
                {
                  lbounds = false;
                  pline++;
                }
            }
          else if ( *pline == 'u' )
            {
              gridtype = GRID_UNSTRUCTURED;
              pline++;
              if ( *pline == '0' )
                {
                  lbounds = false;
                  pline++;
                }
            }
          if ( *pline != 0 ) return;          
        }

      if ( inc < 1e-9 ) inc = 1;
    }

  grid->type = gridtype;

  if ( lon1 >= lon2 || lat1 >= lat2 )
    cdoAbort("Invalid grid box: lon1=%g lon2=%g lat1=%g lat2=%g", lon1, lon2, lat1, lat2);

  int nlon = (int) ((lon2 - lon1)/inc + 0.5);
  int nlat = (int) ((lat2 - lat1)/inc + 0.5);

  double *xvals = (double*) Malloc(nlon*sizeof(double));
  double *yvals = (double*) Malloc(nlat*sizeof(double));

  for ( int i = 0; i < nlon; ++i ) xvals[i] = lon1 + inc/2 + i*inc;
  for ( int i = 0; i < nlat; ++i ) yvals[i] = lat1 + inc/2 + i*inc;

  if ( gridtype == GRID_LONLAT )
    {
      grid->xsize = nlon;
      grid->ysize = nlat;
      grid->xvals = xvals;
      grid->yvals = yvals;
      xvals = NULL;
      yvals = NULL;
    }
  else
    {
      double gridsize = nlon*nlat;
      double *xvals2D = (double*) Malloc(gridsize*sizeof(double));
      double *yvals2D = (double*) Malloc(gridsize*sizeof(double));
      for ( int j = 0; j < nlat; j++ )
        for ( int i = 0; i < nlon; i++ )
          {
            xvals2D[j*nlon+i] = xvals[i];
            yvals2D[j*nlon+i] = yvals[j];
          }

      if ( gridtype == GRID_CURVILINEAR )
        {
          grid->xsize = nlon;
          grid->ysize = nlat;
        }
      else
        {
          grid->xsize = gridsize;
          grid->ysize = gridsize;
          if ( lbounds ) grid->nvertex = 4;
        }
      
      grid->xvals = xvals2D;
      grid->yvals = yvals2D;
      
      if ( lbounds && nlon > 1 && nlat > 1 )
        {
          double *xbounds = (double*) Malloc(2*nlon*sizeof(double));
          grid_gen_bounds(nlon, xvals, xbounds);
          
          double *ybounds = (double*) Malloc(2*nlat*sizeof(double));
          grid_gen_bounds(nlat, yvals, ybounds);
          grid_check_lat_borders(2*nlat, ybounds);

          double *xbounds2D = (double*) Malloc(4*gridsize*sizeof(double));
          double *ybounds2D = (double*) Malloc(4*gridsize*sizeof(double));

          grid_gen_xbounds2D(nlon, nlat, xbounds, xbounds2D);
          grid_gen_ybounds2D(nlon, nlat, ybounds, ybounds2D);

          Free(xbounds);
          Free(ybounds);
          grid->xbounds = xbounds2D;
          grid->ybounds = ybounds2D;
        }
   }

  if ( xvals ) Free(xvals);
  if ( yvals ) Free(yvals);
}


int gridFromName(const char *gridname)
{
  const char *pline;
  int gridID = UNDEFID;
  griddes_t grid;
  size_t len;
  char *endptr;

  gridInit(&grid);

  if ( gridname[0] == 't' && gridname[1] == 'l' ) /* tl<RES>grid or tl<RES>spec */
    {
      pline = &gridname[2];
      if ( isdigit((int) *pline) )
	{
	  grid.ntr = atoi(pline);
	  while ( isdigit((int) *pline) ) pline++;
	  if      ( cmpstrlen(pline, "grid", len) == 0 ) grid.type = GRID_GAUSSIAN;
	  else if ( cmpstrlen(pline, "zon",  len) == 0 ) grid.type = GRID_GAUSSIAN;
	  else if ( cmpstrlen(pline, "spec", len) == 0 ) grid.type = GRID_SPECTRAL;
	  else if ( cmpstrlen(pline, "",     len) == 0 ) grid.type = GRID_SPECTRAL;
      
	  if ( pline[len] != 0 ) return gridID;

	  if ( grid.type == GRID_GAUSSIAN )
	    {
	      grid.ysize = ntr2nlat_linear(grid.ntr);
	      grid.np    = grid.ysize/2;
	      if ( cmpstrlen(pline, "zon",  len) == 0 )
		grid.xsize = 1;
	      else
		grid.xsize = compNlon(grid.ysize);

	      grid.def_xfirst = TRUE;
	      grid.def_yfirst = TRUE;	      
	    }
	}
    }
  else if ( gridname[0] == 't' ) /* t<RES>grid or t<RES>spec */
    {
      pline = &gridname[1];
      if ( isdigit((int) *pline) )
	{
	  grid.ntr = atoi(pline);
	  while ( isdigit((int) *pline) ) pline++;
	  if      ( cmpstrlen(pline, "grid", len) == 0 ) grid.type = GRID_GAUSSIAN;
	  else if ( cmpstrlen(pline, "zon",  len) == 0 ) grid.type = GRID_GAUSSIAN;
	  else if ( cmpstrlen(pline, "spec", len) == 0 ) grid.type = GRID_SPECTRAL;
	  else if ( cmpstrlen(pline, "",     len) == 0 ) grid.type = GRID_SPECTRAL;
     
	  if ( pline[len] != 0 ) return gridID;

	  if ( grid.type == GRID_GAUSSIAN )
	    {
	      grid.ysize = ntr2nlat(grid.ntr);
	      grid.np    = grid.ysize/2;
	      if ( cmpstrlen(pline, "zon",  len) == 0 )
		grid.xsize = 1;
	      else
		grid.xsize = compNlon(grid.ysize);

	      grid.def_xfirst = TRUE;
	      grid.def_yfirst = TRUE;	      
	    }
	}
    }
  else if ( gridname[0] == 'r' ) /* r<LON>x<LAT> */
    {
      pline = &gridname[1];
      if ( isdigit((int) *pline) )
	{
	  grid.type = GRID_LONLAT;
	  grid.xsize = atoi(pline);
	  while ( isdigit((int) *pline) ) pline++;
	  pline++;
	  grid.ysize = atoi(pline);
	  while ( isdigit((int) *pline) ) pline++;

	  grid.def_xfirst = TRUE;
	  grid.def_yfirst = TRUE;
	}
    }
  else if ( gridname[0] == 'l' &&  gridname[1] == 'o' && gridname[2] == 'n' ) /* lon=<LON>_lat=<LAT> */
    {
      /* only one gridpoint */
      pline = &gridname[3];
      if ( *pline == '=' ) pline++;
      if ( isdigit((int) *pline) || ispunct((int) *pline) || *pline == '-' )
	{
	  grid.type = GRID_LONLAT;
	  grid.xsize = 1;
	  grid.ysize = 1;
	  grid.xvals = (double*) Malloc(sizeof(double));
	  grid.yvals = (double*) Malloc(sizeof(double));
	  grid.xvals[0] = atof(pline);
	  while ( isdigit((int) *pline) || ispunct((int) *pline) || *pline == '-' ) pline++;
	  if ( *pline == '_' ) pline++;
	  if ( ! (pline[0] == 'l' &&  pline[1] == 'a' && pline[2] == 't') ) return gridID;
	  pline += 3;
	  if ( *pline == '=' ) pline++;
	  if ( isdigit((int) *pline) || ispunct((int) *pline) || *pline == '-' )
	    grid.yvals[0] = atof(pline);
	  else
	    return gridID;
	}
    }
  else if ( gridname[0] == 'g' && gridname[1] == 'm' && gridname[2] == 'e' ) /* gme<NI> */
    {
      pline = &gridname[3];
      if ( isdigit((int) *pline) )
	{
	  long ni = strtol(pline, &endptr, 10);
	  if ( *endptr == 0 )
	    {
	      grid.type = GRID_GME;
	      grid.ni   = ni;
	      grid.nd   = 10;
	      factorni(grid.ni, &grid.ni2, &grid.ni3);
	      grid.size = (grid.ni+1)*(grid.ni+1)*10;
	    }
	}
    }
  else if ( gridname[0] == 'n' && gridname[1] == 'i' ) /* ni<NI> */
    {
      pline = &gridname[2];
      if ( isdigit((int) *pline) )
	{
	  long ni = strtol(pline, &endptr, 10);
	  if ( *endptr == 0 )
	    {
	      grid.type = GRID_GME;
	      grid.ni   = ni;
	      grid.nd   = 10;
	      factorni(grid.ni, &grid.ni2, &grid.ni3);
	      grid.size = (grid.ni+1)*(grid.ni+1)*10;
	    }
	}
    }
  else if ( gridname[0] == 'n' ) /* n<N> */
    {
      pline = &gridname[1];
      if ( isdigit((int) *pline) )
	{
	  long np = strtol(pline, &endptr, 10);
	  pline = endptr;

	  if ( cmpstrlen(pline, "zon",  len) == 0 )
	    {
	      grid.xsize = 1;
	      pline += 3;
	    }
	  else if ( *pline == 'b' )
	    {
	      grid.genBounds = TRUE;
	      pline++;
	    }

	  if ( *pline == 0 )
	    {
	      grid.type  = GRID_GAUSSIAN;
	      grid.np    = np;
	      grid.ysize = np*2;
	      if ( !grid.xsize ) grid.xsize = compNlon(grid.ysize);

	      grid.def_xfirst = TRUE;
	      grid.def_yfirst = TRUE;	      
	    }
	}
    }
  else if ( gridname[0] == 'g' && isdigit(gridname[1])) /* g<LON>x<LAT> or g<SIZE> */
    {
      pline = &gridname[1];
      if ( isdigit((int) *pline) )
	{
	  grid.type = GRID_GENERIC;
	  grid.xsize = atoi(pline);
	  while ( isdigit((int) *pline) ) pline++;
	  if ( *pline )
	    {
	      pline++;
	      grid.ysize = atoi(pline);
	      while ( isdigit((int) *pline) ) pline++;
	    }
	  else if ( grid.xsize == 1 )
	    {
	      grid.size  = 1;
	      grid.xsize = 0;
	    }
	}
    }
  else if ( strncmp(gridname, "germany", 7) == 0 ) /* germany_Xdeg */
    {
      double lon1 =   5.6, lon2 = 15.2;
      double lat1 =  47.1, lat2 = 55.1;
      double dll = 0.1;

      pline = &gridname[7];

      gen_grid_lonlat(&grid, pline, dll, lon1, lon2, lat1, lat2);
    }
  else if ( strncmp(gridname, "europe", 6) == 0 ) /* europe_Xdeg */
    {
      double lon1 = -30, lon2 = 60;
      double lat1 =  30, lat2 = 80;
      double dll = 1;

      pline = &gridname[6];

      gen_grid_lonlat(&grid, pline, dll, lon1, lon2, lat1, lat2);
    }
  else if ( strncmp(gridname, "africa", 6) == 0 ) /* africa_Xdeg */
    {
      double lon1 = -20, lon2 = 60;
      double lat1 = -40, lat2 = 40;
      double dll = 1;

      pline = &gridname[6];

      gen_grid_lonlat(&grid, pline, dll, lon1, lon2, lat1, lat2);
    }
  else if ( strncmp(gridname, "global", 6) == 0 ) /* global_Xdeg */
    {
      double lon1 = -180, lon2 = 180;
      double lat1 =  -90, lat2 =  90;
      double dll = 1;

      pline = &gridname[6];
  
      gen_grid_lonlat(&grid, pline, dll, lon1, lon2, lat1, lat2);
    }

  if ( grid.type != -1 ) gridID = gridDefine(grid);

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

      gridID = gridFromName(gridfile);

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
	  gridID = gridFromPingo(gfp, filename);
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

/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2012 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

#if  defined  (HAVE_LIBNETCDF)
#  include "netcdf.h"
#endif

#include <ctype.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "griddes.h"
#include "error.h"

#define  deg2rad  (M_PI/180.)   /* conversion for deg to rad */
#define  rad2deg  (180./M_PI)   /* conversion for rad to deg */

#define  cmpstr(s1, s2, len)  (memcmp(s1, s2, len = strlen(s2)))

/*
int  extInqPrec(int fileID);

int  extReadHeader(int fileID, int *header);
int  extReadDataDP(int fileID, double *data);

int  extOpen(const char *filename, const char *mode);
void extClose(int fileID);
*/

#define UNDEFID -1

#define MAX_LINE_LEN 65536


void gridInit(grid_t *grid)
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
  grid->path[0]       = 0;
  grid->xname[0]      = 0;
  grid->xlongname[0]  = 0;
  grid->xunits[0]     = 0;
  grid->yname[0]      = 0;
  grid->ylongname[0]  = 0;
  grid->yunits[0]     = 0;
}


int getoptname(char *optname, const char *optstring, int nopt)
{
  int i, nerr = 0;
  size_t namelen;
  const char *pname;
  const char *pend;

  pname = optstring;
  pend  = optstring;

  for ( i = 0; i < nopt; i++ )
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

  return (nerr);
}


int gridDefine(grid_t grid)
{
  int gridID = UNDEFID;
  int i;

  switch ( grid.type )
    {
    case GRID_GENERIC:
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_SINUSOIDAL:
      {
	if ( grid.size != 1 )
	  {
	    if ( grid.xsize == 0 ) Error("xsize undefined!");
	    if ( grid.ysize == 0 ) Error("ysize undefined!");
	  }

	if ( grid.size == 0 ) grid.size = grid.xsize*grid.ysize;

	gridID = gridCreate(grid.type, grid.size);

	if ( grid.xsize > 0 ) gridDefXsize(gridID, grid.xsize);
	if ( grid.ysize > 0 ) gridDefYsize(gridID, grid.ysize);
	if ( grid.np    > 0 ) gridDefNP(gridID, grid.np);

	gridDefPrec(gridID, grid.prec);

	if ( (grid.def_xfirst || grid.def_xlast || grid.def_xinc) && grid.xvals == NULL )
	  {
	    grid.xvals = (double *) malloc(grid.xsize*sizeof(double));
	    gridGenXvals(grid.xsize, grid.xfirst, grid.xlast, grid.xinc, grid.xvals);

	    if ( grid.genBounds && grid.xbounds == NULL && grid.xsize > 1 )
	      {
		grid.nvertex = 2;
		grid.xbounds = (double *) malloc(grid.xsize*grid.nvertex*sizeof(double));
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
	    grid.yvals = (double *) malloc(grid.ysize*sizeof(double));
	    gridGenYvals(grid.type, grid.ysize, grid.yfirst, grid.ylast, grid.yinc, grid.yvals);

	    if ( grid.genBounds && grid.ybounds == NULL && grid.ysize > 1 )
	      {
		grid.nvertex = 2;
		grid.ybounds = (double *) malloc(grid.ysize*grid.nvertex*sizeof(double));
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
	    free(grid.xvals);
	  }

	if ( grid.yvals )
	  {
	    gridDefYvals(gridID, grid.yvals);
	    free(grid.yvals);
	  }

	if ( grid.nvertex )
	  gridDefNvertex(gridID, grid.nvertex);

	if ( grid.xbounds )
	  {
	    gridDefXbounds(gridID, grid.xbounds);
	    free(grid.xbounds);
	  }

	if ( grid.ybounds )
	  {
	    gridDefYbounds(gridID, grid.ybounds);
	    free(grid.ybounds);
	  }

	if ( grid.isRotated )
	  {
	    gridDefXpole(gridID, grid.xpole);
	    gridDefYpole(gridID, grid.ypole);
	  }

	if ( grid.mask )
	  {
	    gridDefMask(gridID, grid.mask);
	    free(grid.mask);
	  }

	break;
      }
    case GRID_CURVILINEAR:
    case GRID_UNSTRUCTURED:
      {
	if ( grid.size == 0 ) grid.size = grid.xsize*grid.ysize;

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
	    gridDefNvertex(gridID, grid.nvertex);
	  }

	if ( grid.xvals )
	  {
	    gridDefXvals(gridID, grid.xvals);
	    free(grid.xvals);
	  }

	if ( grid.yvals )
	  {
	    gridDefYvals(gridID, grid.yvals);
	    free(grid.yvals);
	  }

	if ( grid.area )
	  {
	    gridDefArea(gridID, grid.area);
	    free(grid.area);
	  }

	if ( grid.xbounds )
	  {
	    gridDefXbounds(gridID, grid.xbounds);
	    free(grid.xbounds);
	  }

	if ( grid.ybounds )
	  {
	    gridDefYbounds(gridID, grid.ybounds);
	    free(grid.ybounds);
	  }

	if ( grid.mask )
	  {
	    gridDefMask(gridID, grid.mask);
	    free(grid.mask);
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
	    free(grid.mask);
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
	    grid.xvals = (double *) malloc(grid.xsize*sizeof(double));
	    for ( i = 0; i < grid.xsize; ++i )
	      grid.xvals[i] = grid.xfirst + i*grid.xinc;
	  }

	if ( grid.def_yfirst && grid.def_yinc && grid.yvals == NULL )
	  {
	    grid.yvals = (double *) malloc(grid.ysize*sizeof(double));
	    for ( i = 0; i < grid.ysize; ++i )
	      grid.yvals[i] = grid.yfirst + i*grid.yinc;
	  }

	if ( grid.xvals )
	  {
	    gridDefXvals(gridID, grid.xvals);
	    free(grid.xvals);
	  }

	if ( grid.yvals )
	  {
	    gridDefYvals(gridID, grid.yvals);
	    free(grid.yvals);
	  }	

	if ( grid.def_lon_0     == FALSE ) Error("lon_0 undefined!");
	if ( grid.def_lat_0     == FALSE ) Error("lat_0 undefined!");
	if ( grid.def_lat_1     == FALSE ) Error("lat_1 undefined!");
	if ( grid.def_lat_2     == FALSE ) grid.def_lat_2 = grid.def_lat_1;

	gridDefLcc2(gridID, grid.a, grid.lon_0, grid.lat_0, grid.lat_1, grid.lat_2);

	if ( grid.mask )
	  {
	    gridDefMask(gridID, grid.mask);
	    free(grid.mask);
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
	    free(grid.mask);
	  }

	break;
      }
    case GRID_REFERENCE:
      {
	gridID = gridCreate(grid.type, grid.size);

	gridDefNumber(gridID, grid.number);
	gridDefPosition(gridID, grid.position);
	if ( *grid.path ) gridDefReference(gridID, grid.path);

	break;
      }
    default:
      {
	if ( grid.type == -1 )
	  Error("Undefined grid type!");
	else
	  Error("Unsupported grid type: %s", gridNamePtr(grid.type));
      }
    }

  if ( grid.xname[0]     ) gridDefXname(gridID, grid.xname);
  if ( grid.xlongname[0] ) gridDefXlongname(gridID, grid.xlongname);
  if ( grid.xunits[0]    ) gridDefXunits(gridID, grid.xunits);
  if ( grid.yname[0]     ) gridDefYname(gridID, grid.yname);
  if ( grid.ylongname[0] ) gridDefYlongname(gridID, grid.ylongname);
  if ( grid.yunits[0]    ) gridDefYunits(gridID, grid.yunits);

  return (gridID);
}

static
char *skipSeparator(char *pline)
{
  while ( isspace((int) *pline) ) pline++;
  if ( *pline == '=' || *pline == ':' ) pline++;
  while ( isspace((int) *pline) ) pline++;

  return (pline);
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

/*
double *readfield(grid_t *grid, int record, char *format, char *filename)
{
  int fileID, rxysize, ierr, irec;
  double *vals;

  if ( grid->size == 0 )  Error("grid size = 0!");
  if ( format == NULL )   Error("format undefined!");
  if ( filename == NULL ) Error("file name undefined!");

  vals = (double *) malloc(grid->size*sizeof(double));

  if ( strcmp(format, "extra") == 0 )
    {
      int header[4];
      fileID = extOpen(filename, "r");
      if ( fileID == UNDEFID ) SysError(filename);

      for ( irec = 0; irec < record; irec++ )
	{
	  ierr = extReadHeader(fileID, header);
	  if ( ierr <= 0 ) Error("Record %d unexpected EOF in file %s", irec+1, filename);
	}
      grid->prec   = extInqPrec(fileID);
      rxysize = header[3];
      if ( rxysize != (int) grid->size ) Error("unexpected record size of %d!", rxysize);
      ierr = extReadDataDP(fileID, vals);
      extClose(fileID);
    }
  else
    Error("format %s unsupported!", format);

  return (vals);
}
*/
/*
double *readfield4(grid_t *grid, int record, char *format, char *filename)
{
  int fileID, rxysize, ierr, irec;
  double *vals;

  if ( grid->size == 0 )  Error("grid size = 0!");
  if ( format == NULL )   Error("format undefined!");
  if ( filename == NULL ) Error("file name undefined!");

  vals  = (double *) malloc(4*grid->size*sizeof(double));

  if ( strcmp(format, "extra") == 0 )
    {
      int header[4];
      fileID = extOpen(filename, "r");
      if ( fileID == UNDEFID ) SysError(filename);

      for ( irec = 0; irec < record; irec++ )
	{
	  ierr = extReadHeader(fileID, header);
	  if ( ierr <= 0 ) Error("Record %d unexpected EOF in file %s", irec+1, filename);
	}
      grid->prec   = extInqPrec(fileID);
      rxysize = header[3];
      if ( rxysize != (int) (4*grid->size) ) Error("unexpected record size of %d!", rxysize);
      ierr = extReadDataDP(fileID, vals);

      extClose(fileID);
    }
  else
    Error("format %s unsupported!", format);

  return (vals);
}
*/

double readflt(const char *name, const char *pline)
{
  double val;
  char *endptr;

  val = strtod(pline, &endptr);
  if ( pline == endptr )
    Warning("Couldn't read value for %s, set to zero!", name);

  return (val);
}


int gridFromFile(FILE *gfp, const char *dname)
{
  char line[MAX_LINE_LEN], *pline;
  int gridID = -1;
  int size;
  size_t len;
  grid_t grid;

  gridInit(&grid);

  while ( readline(gfp, line, MAX_LINE_LEN) )
    {
      if ( line[0] == '#' ) continue;
      if ( line[0] == '\0' ) continue;
      pline = line;
      while ( isspace((int) *pline) ) pline++;
      if ( pline[0] == '\0' ) continue;
      if ( cmpstr(pline, "gridtype", len) == 0 )
	{
	  pline = skipSeparator(pline + len);
	  if ( cmpstr(pline, "lonlat", len)  == 0 ||
	       cmpstr(pline, "latlon", len)  == 0 )
	    {
	      grid.type = GRID_LONLAT;
	      grid.nvertex = 2;
	    }
	  else if ( cmpstr(pline, "gaussian", len)  == 0 )
	    {
	      grid.type = GRID_GAUSSIAN;
	      grid.nvertex = 2;
	    }
	  else if ( cmpstr(pline, "curvilinear", len)  == 0 )
	    {
	      grid.type = GRID_CURVILINEAR;
	      grid.nvertex = 4;
	    }
	  else if ( cmpstr(pline, "spectral", len)  == 0 )
	    grid.type = GRID_SPECTRAL;
	  else if ( cmpstr(pline, "unstructured", len)  == 0 )
	    grid.type = GRID_UNSTRUCTURED;
	  else if ( cmpstr(pline, "cell", len)  == 0 )
	    grid.type = GRID_UNSTRUCTURED;
	  else if ( cmpstr(pline, "gme", len)  == 0 )
	    grid.type = GRID_GME;
	  else if ( cmpstr(pline, "reference", len)  == 0 )
	    grid.type = GRID_REFERENCE;
	  else if ( cmpstr(pline, "lcc2", len)  == 0 )
	    grid.type = GRID_LCC2;
	  else if ( cmpstr(pline, "lcc", len)  == 0 )
	    grid.type = GRID_LCC;
	  else if ( cmpstr(pline, "lambert", len)  == 0 )
	    grid.type = GRID_LCC;
	  else if ( cmpstr(pline, "sinusoidal", len)  == 0 )
	    grid.type = GRID_SINUSOIDAL;
	  else if ( cmpstr(pline, "laea", len)  == 0 )
	    grid.type = GRID_LAEA;
	  else if ( cmpstr(pline, "generic", len)  == 0 )
	    grid.type = GRID_GENERIC;
	  else
	    Warning("Invalid grid name : %s", pline);
	}
      else if ( cmpstr(pline, "gridprec", len)  == 0 )
	{
	  grid.prec = atol(skipSeparator(pline + len));
	}
      else if ( cmpstr(pline, "gridsize", len)  == 0 )
	{
	  grid.size = atol(skipSeparator(pline + len));
	}
      else if ( cmpstr(pline, "truncation", len)  == 0 )
	{
	  grid.ntr = atoi(skipSeparator(pline + len));
	}
      else if ( cmpstr(pline, "np", len)  == 0 )
	{
	  grid.np = atoi(skipSeparator(pline + len));
	}
      else if ( cmpstr(pline, "complexpacking", len)  == 0 )
	{
	  grid.lcomplex = atoi(skipSeparator(pline + len));
	}
      else if ( cmpstr(pline, "xname", len)  == 0 )
	{
	  strcpy(grid.xname, skipSeparator(pline + len));
	}
      else if ( cmpstr(pline, "xlongname", len)  == 0 )
	{
	  strcpy(grid.xlongname, skipSeparator(pline + len));
	}
      else if ( cmpstr(pline, "xunits", len)  == 0 )
	{
	  strcpy(grid.xunits, skipSeparator(pline + len));
	}
      else if ( cmpstr(pline, "yname", len)  == 0 )
	{
	  strcpy(grid.yname, skipSeparator(pline + len));
	}
      else if ( cmpstr(pline, "ylongname", len)  == 0 )
	{
	  strcpy(grid.ylongname, skipSeparator(pline + len));
	}
      else if ( cmpstr(pline, "yunits", len)  == 0 )
	{
	  strcpy(grid.yunits, skipSeparator(pline + len));
	}
      else if ( cmpstr(pline, "nvertex", len)  == 0 )
	{
	  grid.nvertex = atol(skipSeparator(pline + len));
	}
      else if ( cmpstr(pline, "ni", len)  == 0 )
	{
	  grid.ni = atol(skipSeparator(pline + len));
          grid.nd = 10;
	}
      else if ( cmpstr(pline, "position", len)  == 0 )
	{
	  grid.position = atol(skipSeparator(pline + len));
	}
      else if ( cmpstr(pline, "number", len)  == 0 )
	{
	  grid.number = atol(skipSeparator(pline + len));
	}
      else if ( cmpstr(pline, "path", len)  == 0 )
	{
	  strcpy(grid.path, skipSeparator(pline + len));
	}
      else if ( cmpstr(pline, "xsize", len)  == 0 )
	{
	  grid.xsize = atol(skipSeparator(pline + len));
	}
      else if ( cmpstr(pline, "nlon", len)  == 0 )
	{
	  grid.xsize = atol(skipSeparator(pline + len));
	}
      else if ( cmpstr(pline, "ysize", len)  == 0 )
	{
	  grid.ysize = atol(skipSeparator(pline + len));
	}
      else if ( cmpstr(pline, "nlat", len)  == 0 )
	{
	  grid.ysize = atol(skipSeparator(pline + len));
	}
      else if ( cmpstr(pline, "xfirst", len)  == 0 )
	{
	  grid.xfirst = readflt("xfirst", skipSeparator(pline + len));
	  grid.def_xfirst = TRUE;
	}
      else if ( cmpstr(pline, "lonfirst", len)  == 0 )
	{
	  grid.xfirst = readflt("lonfirst", skipSeparator(pline + len));
	  grid.def_xfirst = TRUE;
	}
      else if ( cmpstr(pline, "yfirst", len)  == 0 )
	{
	  grid.yfirst = readflt("yfirst", skipSeparator(pline + len));
	  grid.def_yfirst = TRUE;
	}
      else if ( cmpstr(pline, "latfirst", len)  == 0 )
	{
	  grid.yfirst = readflt("latfirst", skipSeparator(pline + len));
	  grid.def_yfirst = TRUE;
	}
      else if ( cmpstr(pline, "xlast", len)  == 0 )
	{
	  grid.xlast = readflt("xlast", skipSeparator(pline + len));
	  grid.def_xlast = TRUE;
	}
      else if ( cmpstr(pline, "lonlast", len)  == 0 )
	{
	  grid.xlast = readflt("lonlast", skipSeparator(pline + len));
	  grid.def_xlast = TRUE;
	}
      else if ( cmpstr(pline, "ylast", len)  == 0 )
	{
	  grid.ylast = readflt("ylast", skipSeparator(pline + len));
	  grid.def_ylast = TRUE;
	}
      else if ( cmpstr(pline, "latlast", len)  == 0 )
	{
	  grid.ylast = readflt("latlast", skipSeparator(pline + len));
	  grid.def_ylast = TRUE;
	}
      else if ( cmpstr(pline, "xinc", len)  == 0 )
	{
	  grid.xinc = readflt("xinc", skipSeparator(pline + len));
	  grid.def_xinc = TRUE;
	}
      else if ( cmpstr(pline, "loninc", len)  == 0 )
	{
	  grid.xinc = readflt("loninc", skipSeparator(pline + len));
	  grid.def_xinc = TRUE;
	}
      else if ( cmpstr(pline, "yinc", len)  == 0 )
	{
	  grid.yinc = readflt("yinc", skipSeparator(pline + len));
	  grid.def_yinc = TRUE;
	}
      else if ( cmpstr(pline, "latinc", len)  == 0 )
	{
	  grid.yinc = readflt("latinc", skipSeparator(pline + len));
	  grid.def_yinc = TRUE;
	}
      else if ( cmpstr(pline, "originLon", len)  == 0 )
	{
	  grid.originLon = readflt("originLon", skipSeparator(pline + len));
	  grid.def_originLon = TRUE;
	}
      else if ( cmpstr(pline, "originLat", len)  == 0 )
	{
	  grid.originLat = readflt("originLat", skipSeparator(pline + len));
	  grid.def_originLat = TRUE;
	}
      else if ( cmpstr(pline, "lonParY", len)  == 0 )
	{
	  grid.lonParY = readflt("lonParY", skipSeparator(pline + len));
	  grid.def_lonParY = TRUE;
	}
      else if ( cmpstr(pline, "lat1", len)  == 0 )
	{
	  grid.lat1 = readflt("lat1", skipSeparator(pline + len));
	  grid.def_lat1 = TRUE;
	}
      else if ( cmpstr(pline, "lat2", len)  == 0 )
	{
	  grid.lat2 = readflt("lat2", skipSeparator(pline + len));
	  grid.def_lat2 = TRUE;
	}
      else if ( cmpstr(pline, "projection", len)  == 0 )
	{
	  pline = skipSeparator(pline + len);
	  if      ( cmpstr(pline, "north", len) == 0 )
	    {
	      grid.projflag = 0;
	      grid.scanflag = 64;
	    }
	  else if ( cmpstr(pline, "south", len) == 0 )
	    {
	      grid.projflag = 128;
	      grid.scanflag = 64;
	    }
	  else
	    Warning("Invalid projection : %s", pline);
	}
      else if ( cmpstr(pline, "a", len)  == 0 )
	{
	  grid.a = readflt("a", skipSeparator(pline + len));
	}
      else if ( cmpstr(pline, "lon_0", len)  == 0 )
	{
	  grid.lon_0 = readflt("lon_0", skipSeparator(pline + len));
	  grid.def_lon_0 = TRUE;
	}
      else if ( cmpstr(pline, "lat_0", len)  == 0 )
	{
	  grid.lat_0 = readflt("lat_0", skipSeparator(pline + len));
	  grid.def_lat_0 = TRUE;
	}
      else if ( cmpstr(pline, "lat_1", len)  == 0 )
	{
	  grid.lat_1 = readflt("lat_1", skipSeparator(pline + len));
	  grid.def_lat_1 = TRUE;
	}
      else if ( cmpstr(pline, "lat_2", len)  == 0 )
	{
	  grid.lat_2 = readflt("lat_2", skipSeparator(pline + len));
	  grid.def_lat_2 = TRUE;
	}
      else if ( cmpstr(pline, "xnpole", len)  == 0 )
	{
	  grid.xpole = readflt("xnpole", skipSeparator(pline + len));
	  grid.isRotated = TRUE;
	}
      else if ( cmpstr(pline, "lonpole", len)  == 0 )
	{
	  grid.xpole = readflt("lonpole", skipSeparator(pline + len));
	  grid.isRotated = TRUE;
	}
      else if ( cmpstr(pline, "ynpole", len)  == 0 )
	{
	  grid.ypole = readflt("ynpole", skipSeparator(pline + len));
	  grid.isRotated = TRUE;
	}
      else if ( cmpstr(pline, "latpole", len)  == 0 )
	{
	  grid.ypole = readflt("latpole", skipSeparator(pline + len));
	  grid.isRotated = TRUE;
	}
      else if ( cmpstr(pline, "gridlatlon", len)  == 0 )
	{
	  int i;
	  double flat = 0, flon = 0;
	  if ( grid.size == 0 ) grid.size = grid.xsize * grid.ysize;
	  
	  grid.xvals = (double *) malloc(grid.size*sizeof(double));
	  grid.yvals = (double *) malloc(grid.size*sizeof(double));
	  for ( i = 0; i < (int) grid.size; i++ )
	    {
	      if ( ! readline(gfp, line, MAX_LINE_LEN) )
		{
		  Warning("Incomplete command: >gridlatlon<");
		  break;
		}
	      sscanf(line, "%lg %lg", &flat, &flon);
	      grid.yvals[i] = flat;
	      grid.xvals[i] = flon;
	    }
	}
      else if ( cmpstr(pline, "mask", len)  == 0 )
	{
	  int i = 0;
	  long lval;
	  char *endptr;

	  size = grid.size;

	  if ( size > 0 )
	    {
	      long count = 0;
	      pline = skipSeparator(pline + len);
	      grid.mask = (int *) malloc(size*sizeof(int));

	      for ( i = 0; i < size; i++ )
		{
		  endptr = pline;
		  lval = strtol(pline, &endptr, 10);
		  if ( pline == endptr )
		    {
		      if ( ! readline(gfp, line, MAX_LINE_LEN) )
			{
			  Warning("Incomplete command: >mask<");
			  break;
			}
		      pline = line;
		      lval = strtol(pline, &endptr, 10);
		    }
		  grid.mask[i] = (int)lval;
		  if ( grid.mask[i] == 1 ) count++;
		  pline = endptr;
		}

	      if ( count == size )
		{
		  free(grid.mask);
		  grid.mask = NULL;
		}
	    }
	  else
	    Warning("gridsize undefined!");
	}
      else if ( cmpstr(pline, "xvals", len)  == 0 )
	{
	  int i = 0;
	  double fval;
	  char *endptr;

	  if ( grid.type == GRID_CURVILINEAR || grid.type == GRID_UNSTRUCTURED )
	    size = grid.size;
	  else
	    size = grid.xsize;

	  if ( size > 0 )
	    {
	      pline = skipSeparator(pline + len);
	      grid.xvals = (double *) malloc(size*sizeof(double));

	      for ( i = 0; i < size; i++ )
		{
		  endptr = pline;
		  fval = strtod(pline, &endptr);
		  if ( pline == endptr )
		    {
		      if ( ! readline(gfp, line, MAX_LINE_LEN) )
			{
			  Warning("Incomplete command: >xvals<");
			  break;
			}
		      pline = line;
		      fval = strtod(pline, &endptr);
		    }
		  grid.xvals[i] = fval;
		  pline = endptr;
		}
	    }
	  else
	    Warning("xsize or gridsize undefined!");
	}
      else if ( cmpstr(pline, "yvals", len)  == 0 )
	{
	  int i = 0;
	  double fval;
	  char *endptr;

	  if ( grid.type == GRID_CURVILINEAR || grid.type == GRID_UNSTRUCTURED )
	    size = grid.size;
	  else
	    size = grid.ysize;

	  if ( size > 0 )
	    {
	      pline = skipSeparator(pline + len);
	      grid.yvals = (double *) malloc(size*sizeof(double));

	      for ( i = 0; i < size; i++ )
		{
		  endptr = pline;
		  fval = strtod(pline, &endptr);
		  if ( pline == endptr )
		    {
		      if ( ! readline(gfp, line, MAX_LINE_LEN) )
			{
			  Warning("Incomplete command: >yvals<");
			  break;
			}
		      pline = line;
		      fval = strtod(pline, &endptr);
		    }
		  grid.yvals[i] = fval;
		  pline = endptr;
		}
	    }
	  else
	    Warning("ysize or gridsize undefined!");
	}
      else if ( cmpstr(pline, "xbounds", len)  == 0 )
	{
	  int i = 0;
	  double fval;
	  char *endptr;

	  if ( grid.nvertex == 0 )
	    {
	      if ( grid.type == GRID_CURVILINEAR ) grid.nvertex = 4;
	    }

	  if ( grid.type == GRID_CURVILINEAR || grid.type == GRID_UNSTRUCTURED )
	    size = grid.size;
	  else
	    size = grid.xsize;

	  if ( size > 0 && grid.nvertex > 0 )
	    {	  
	      pline = skipSeparator(pline + len);
	      grid.xbounds = (double *) malloc(size*grid.nvertex*sizeof(double));

	      for ( i = 0; i < (int) (size*grid.nvertex); i++ )
		{
		  endptr = pline;
		  fval = strtod(pline, &endptr);
		  if ( pline == endptr )
		    {
		      if ( ! readline(gfp, line, MAX_LINE_LEN) )
			{
			  Warning("Incomplete command: >xbounds<");
			  break;
			}
		      pline = line;
		      fval = strtod(pline, &endptr);
		    }
		  grid.xbounds[i] = fval;
		  pline = endptr;
		}
	    }
	  else
	    {
	      if ( size         == 0 ) Warning("xsize or gridsize undefined!");
	      if ( grid.nvertex == 0 ) Warning("nvertex undefined!");
	    }
	}
      else if ( cmpstr(pline, "ybounds", len)  == 0 )
	{
	  int i = 0;
	  double fval;
	  char *endptr;

	  if ( grid.nvertex == 0 )
	    {
	      if ( grid.type == GRID_CURVILINEAR ) grid.nvertex = 4;
	    }

	  if ( grid.type == GRID_CURVILINEAR || grid.type == GRID_UNSTRUCTURED )
	    size = grid.size;
	  else
	    size = grid.ysize;

	  if ( size > 0 && grid.nvertex > 0 )
	    {	  
	      pline = skipSeparator(pline + len);
	      grid.ybounds = (double *) malloc(size*grid.nvertex*sizeof(double));

	      for ( i = 0; i < (int) (size*grid.nvertex); i++ )
		{
		  endptr = pline;
		  fval = strtod(pline, &endptr);
		  if ( pline == endptr )
		    {
		      if ( ! readline(gfp, line, MAX_LINE_LEN) )
			{
			  Warning("Incomplete command: >ybounds<");
			  break;
			}
		      pline = line;
		      fval = strtod(pline, &endptr);
		    }
		  grid.ybounds[i] = fval;
		  pline = endptr;
		}
	    }
	  else
	    {
	      if ( grid.ysize   == 0 ) Warning("ysize or gridsize undefined!");
	      if ( grid.nvertex == 0 ) Warning("nvertex undefined!");
	    }
	}
      else
	{
	  if ( grid.type != UNDEFID )
	    Warning("Invalid grid command : >%s<", pline);
	}
    }
  /*
  printf("gridtype %d\n", grid.type);
  printf("gridsize %d\n", grid.size);
  printf("xsize %d\n", grid.xsize);
  printf("ysize %d\n", grid.ysize);
  */
  if ( grid.type != UNDEFID ) gridID = gridDefine(grid);

  return (gridID);
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
  int read_items;

  skip_nondigit_lines(gfp);

  if ( feof(gfp) ) return(0);

  *ival = 0;
  read_items = fscanf(gfp, "%d", ival);

  return (read_items);
}


int input_darray(FILE *gfp, int n_values, double *array)
{
  int i;
  int read_items;

  if ( n_values <= 0 ) return (0);

  read_items = 0;
  for ( i = 0; i < n_values; i++ )
    {
      skip_nondigit_lines(gfp);

      if ( feof(gfp) ) break;

      read_items += fscanf(gfp, "%lg", &array[i]);

      if ( feof(gfp) ) break;
    }

  return (read_items);
}


int gridFromPingo(FILE *gfp, const char *dname)
{
  int gridID = -1;
  int i;
  int nlon, nlat;
  int lgauss = FALSE;
  grid_t grid;

  gridInit(&grid);

  if ( ! input_ival(gfp, &nlon) ) return (gridID);
  if ( ! input_ival(gfp, &nlat) ) return (gridID);

  if ( nlon > 0 && nlon < 9999 && nlat > 0 && nlat < 9999 )
    {
      grid.xsize = nlon;
      grid.ysize = nlat;

      grid.xvals = (double *) malloc(grid.xsize*sizeof(double));
      grid.yvals = (double *) malloc(grid.ysize*sizeof(double));

      if ( ! input_ival(gfp, &nlon) ) return (gridID);
      if ( nlon == 2 )
	{
	  if ( input_darray(gfp, 2, grid.xvals) != 2 ) return (gridID);
	  grid.xvals[1] -= 360 * floor((grid.xvals[1] - grid.xvals[0]) / 360);

	  if ( grid.xsize > 1 )
	    if ( IS_EQUAL(grid.xvals[0], grid.xvals[1]) )
	      grid.xvals[1] += 360;

	  for ( i = 0; i < (int)grid.xsize; i++ )
	    grid.xvals[i] = grid.xvals[0] + i*(grid.xvals[1] - grid.xvals[0]);
	}
      else if ( nlon == (int)grid.xsize )
	{
	  if ( input_darray(gfp, nlon, grid.xvals) != nlon ) return (gridID);
	  for ( i = 0; i < nlon - 1; i++ )
	    if ( grid.xvals[i+1] <= grid.xvals[i] ) break;

	  for ( i++; i < nlon; i++ )
	    {
	      grid.xvals[i] += 360;
	      if ( i < nlon - 1 && grid.xvals[i+1] + 360 <= grid.xvals[i] )
		{
		  Message("Longitudes are not in ascending order!");
		  return (gridID);
		}
	    }
	}
      else
	return (gridID);

      if ( ! input_ival(gfp, &nlat) ) return (gridID);
      if ( nlat == 2 )
	{
	  if ( input_darray(gfp, 2, grid.yvals) != 2 ) return (gridID);
	  for ( i = 0; i < (int)grid.ysize; i++ )
	    grid.yvals[i] = grid.yvals[0] + i*(grid.yvals[1] - grid.yvals[0]);
	}
      else if ( nlat == (int)grid.ysize )
	{
	  if ( input_darray(gfp, nlat, grid.yvals) != nlat ) return (gridID);
	}
      else
	return (gridID);

      if ( grid.yvals[0]      >  90.001  || 
	   grid.yvals[nlat-1] >  90.001  || 
	   grid.yvals[0]      < -90.001  || 
	   grid.yvals[nlat-1] < -90.001 )
	{
	  Message("Latitudes must be between 90 and -90!");
	  return (gridID);
	}

      for ( i = 0; i < nlat - 1; i++ )
	if ( IS_EQUAL(grid.yvals[i+1], grid.yvals[i]) || (i < nlat - 2 &&
	    ((grid.yvals[i+1] > grid.yvals[i]) != (grid.yvals[i+2] > grid.yvals[i+1]))) )
	  {
	    Message("Latitudes must be in descending or ascending order!");
	    return (gridID);
	  }
		    
      if ( nlat > 2 ) /* check if gaussian */
	{
	  double *yvals, *yw;
	  yvals = (double *) malloc(grid.ysize*sizeof(double));
	  yw    = (double *) malloc(grid.ysize*sizeof(double));
	  gaussaw(yvals, yw, grid.ysize);
	  free(yw);
	  for ( i = 0; i < (int) grid.ysize; i++ )
	    yvals[i] = asin(yvals[i])*rad2deg;

	  for ( i = 0; i < (int) grid.ysize; i++ )
	    if ( fabs(yvals[i] - grid.yvals[i]) > ((yvals[0] - yvals[1])/500) ) break;
		      
	  if ( i == (int) grid.ysize ) lgauss = TRUE;

	  free(yvals);
	}

      if ( lgauss )
	grid.type = GRID_GAUSSIAN;
      else
	grid.type = GRID_LONLAT;
    }
  
  if ( grid.type != UNDEFID ) gridID = gridDefine(grid);

  return (gridID);
}


int nfc2nlat(int nfc, int ntr)
{
  int nlat;

  nlat = nfc / (ntr+1);
  nlat /= 2;

  return (nlat);
}


int nlat2ntr(int nlat)
{
  int ntr;

  ntr = (nlat*2 - 1) / 3;

  return (ntr);
}


int nlat2ntr_linear(int nlat)
{
  int ntr;

  ntr = (nlat*2 - 1) / 2;

  return (ntr);
}


int ntr2nlat(int ntr)
{
  int nlat, nlat2;

  nlat = NINT((ntr*3.+1.)/2.);
  if ( (nlat % 2) > 0 )
    {
      nlat  = nlat + 1;
      nlat2 = NINT(((ntr+1)*3.+1.)/2.);
      /*
      if ( nlat == nlat2 )
	Error("Computation of latitudes failed for truncation %d", ntr);
      */
    }

  return (nlat);
}


int ntr2nlat_linear(int ntr)
{
  int nlat, nlat2;

  nlat = NINT((ntr*2.+1.)/2.);
  if ( (nlat % 2) > 0 )
    {
      nlat  = nlat + 1;
      nlat2 = NINT(((ntr+1)*2.+1.)/2.);
      /*
      if ( nlat == nlat2 )
	Error("Computation of latitudes failed for truncation %d", ntr);
      */
    }

  return (nlat);
}


int compNlon(int nlat)
{
  int nlon, n;

  nlon = 2 * nlat;

  /* check that FFT works with nlon */
  while ( 1 )
    {
      n = nlon;
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

  return (nlon);
}

static
void gen_grid_lonlat(grid_t *grid, const char *pline, double inc, double lon1, double lon2, double lat1, double lat2)
{
  int nlon, nlat, i;
  int gridtype = GRID_LONLAT;
  char *endptr;

  if ( *pline != 0 )
    {
      if ( *pline == '_' ) pline++;
      else return;

      if ( *pline == 0 ) return;

      if ( ! isdigit((int) *pline) && !ispunct((int) *pline) ) return;

      endptr = (char *) pline;
      inc = strtod(pline, &endptr);
      if ( *endptr != 0 ) return;

      if ( inc < 1e-9 ) inc = 1;
    }

  grid->type = gridtype;

  nlon = (int) ((lon2 - lon1)/inc + 0.5);
  nlat = (int) ((lat2 - lat1)/inc + 0.5);
  grid->xsize = nlon;
  grid->ysize = nlat;

  grid->xvals = (double *) malloc(grid->xsize*sizeof(double));
  grid->yvals = (double *) malloc(grid->ysize*sizeof(double));

  for ( i = 0; i < nlon; ++i ) grid->xvals[i] = lon1 + inc/2 + i*inc;
  for ( i = 0; i < nlat; ++i ) grid->yvals[i] = lat1 + inc/2 + i*inc;
}


int gridFromName(const char *gridname)
{
  const char *pline;
  int gridID = UNDEFID;
  grid_t grid;
  size_t len;

  gridInit(&grid);

  if ( gridname[0] == 't' && gridname[1] == 'l' ) /* tl<RES>grid or tl<RES>spec */
    {
      pline = &gridname[2];
      if ( isdigit((int) *pline) )
	{
	  grid.ntr = atoi(pline);
	  while ( isdigit((int) *pline) ) pline++;
	  if      ( cmpstr(pline, "grid", len) == 0 ) grid.type = GRID_GAUSSIAN;
	  else if ( cmpstr(pline, "zon",  len) == 0 ) grid.type = GRID_GAUSSIAN;
	  else if ( cmpstr(pline, "spec", len) == 0 ) grid.type = GRID_SPECTRAL;
	  else if ( cmpstr(pline, "",     len) == 0 ) grid.type = GRID_SPECTRAL;
      
	  if ( pline[len] != 0 ) return (gridID);

	  if ( grid.type == GRID_GAUSSIAN )
	    {
	      grid.ysize = ntr2nlat_linear(grid.ntr);
	      grid.np    = grid.ysize/2;
	      if ( cmpstr(pline, "zon",  len) == 0 )
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
	  if      ( cmpstr(pline, "grid", len) == 0 ) grid.type = GRID_GAUSSIAN;
	  else if ( cmpstr(pline, "zon",  len) == 0 ) grid.type = GRID_GAUSSIAN;
	  else if ( cmpstr(pline, "spec", len) == 0 ) grid.type = GRID_SPECTRAL;
	  else if ( cmpstr(pline, "",     len) == 0 ) grid.type = GRID_SPECTRAL;
     
	  if ( pline[len] != 0 ) return (gridID);

	  if ( grid.type == GRID_GAUSSIAN )
	    {
	      grid.ysize = ntr2nlat(grid.ntr);
	      grid.np    = grid.ysize/2;
	      if ( cmpstr(pline, "zon",  len) == 0 )
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
	  grid.xvals = (double *) malloc(sizeof(double));
	  grid.yvals = (double *) malloc(sizeof(double));
	  grid.xvals[0] = atof(pline);
	  while ( isdigit((int) *pline) || ispunct((int) *pline) || *pline == '-' ) pline++;
	  if ( *pline == '_' ) pline++;
	  if ( ! (pline[0] == 'l' &&  pline[1] == 'a' && pline[2] == 't') ) return(gridID);
	  pline += 3;
	  if ( *pline == '=' ) pline++;
	  if ( isdigit((int) *pline) || ispunct((int) *pline) || *pline == '-' )
	    grid.yvals[0] = atof(pline);
	  else
	    return (gridID);
	}
    }
  else if ( gridname[0] == 'g' && gridname[1] == 'm' && gridname[2] == 'e' ) /* gme<NI> */
    {
      pline = &gridname[3];
      if ( isdigit((int) *pline) )
	{
	  grid.type = GRID_GME;
	  grid.ni   = atoi(pline);
	  grid.nd   = 10;
	  factorni(grid.ni, &grid.ni2, &grid.ni3);
	  grid.size = (grid.ni+1)*(grid.ni+1)*10;
	}
    }
  else if ( gridname[0] == 'n' && gridname[1] == 'i' ) /* ni<NI> */
    {
      pline = &gridname[2];
      if ( isdigit((int) *pline) )
	{
	  grid.type = GRID_GME;
	  grid.ni   = atoi(pline);
	  grid.nd   = 10;
	  factorni(grid.ni, &grid.ni2, &grid.ni3);
	  grid.size = (grid.ni+1)*(grid.ni+1)*10;
	}
    }
  else if ( gridname[0] == 'n' ) /* n<N> */
    {
      pline = &gridname[1];
      if ( isdigit((int) *pline) )
	{
	  int np;
	  np = atoi(pline);
	  while ( isdigit((int) *pline) ) pline++;

	  grid.type = GRID_GAUSSIAN;
	  grid.np    = np;
	  grid.ysize = np*2;
	  grid.xsize = compNlon(grid.ysize);

	  if ( cmpstr(pline, "zon",  len) == 0 ) 
	    grid.xsize = 1;
	  else if ( *pline == 'b' )
	    grid.genBounds = TRUE;

	  grid.def_xfirst = TRUE;
	  grid.def_yfirst = TRUE;	      
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

  return (gridID);
}


int cdoDefineGrid(const char *gridfile)
{
  FILE *gfp;
  char buffer[4];
  int gridID = -1;
  size_t len;
  struct stat filestat;
  int fileno;
  int isreg = FALSE;

  fileno = open(gridfile, O_RDONLY);
  if ( fileno >= 0 )
    {
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
      if ( read(fileno, buffer, 4) != 4 )
	SysError("Read grid from %s failed!", gridfile);

      close(fileno);

      if ( cmpstr(buffer, "CDF", len) == 0 )
	{
	  if ( cdoDebug ) cdoPrint("Grid from netCDF file");
	  gridID = gridFromNCfile(gridfile);
	}

      if ( gridID == -1 )
	{
	  if ( cmpstr(buffer+1, "HDF", len) == 0 )
	    {
	      if ( cdoDebug ) cdoPrint("Grid from HDF5 file");
	      gridID = gridFromH5file(gridfile);
	    }
	}

      if ( gridID == -1 )
	{
	  int streamID;
	  if ( cdoDebug ) cdoPrint("Grid from CDI file");
	  openLock();
	  streamID = streamOpenRead(gridfile);
	  openUnlock();
	  if ( streamID >= 0 )
	    {
	      int vlistID;
	      vlistID = streamInqVlist(streamID);
	      gridID  = vlistGrid(vlistID, 0);
	      streamClose(streamID);
	    }
	}

      if ( gridID == -1 )
	{
	  if ( cdoDebug ) cdoPrint("grid from ASCII file");
	  gfp = fopen(gridfile, "r");
	  gridID = gridFromFile(gfp, gridfile);
	  fclose(gfp);
	}

      if ( gridID == -1 )
	{
	  if ( cdoDebug ) cdoPrint("grid from PINGO file");
	  gfp = fopen(gridfile, "r");
	  gridID = gridFromPingo(gfp, gridfile);
	  fclose(gfp);
	}

      if ( gridID == -1 ) cdoAbort("Invalid grid description file %s!", gridfile);
    }

  return (gridID);
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

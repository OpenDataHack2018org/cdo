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

#if  defined  (HAVE_LIBNETCDF)
#  include "netcdf.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "error.h"

#define  deg2rad  (M_PI/180.)   /* conversion for deg to rad */
#define  rad2deg  (180./M_PI)   /* conversion for rad to deg */

/*
int  extInqPrec(int fileID);

int  extReadHeader(int fileID, int *header);
int  extReadDataDP(int fileID, double *data);

int  extOpen(const char *filename, const char *mode);
void extClose(int fileID);
*/

#define UNDEFID -1

#define MAX_LINE_LEN 65536


void gridInit(GRID *grid)
{
  grid->xvals         = NULL;
  grid->yvals         = NULL;
  grid->xbounds       = NULL;
  grid->ybounds       = NULL;
  grid->area          = NULL;
  grid->type          = UNDEFID;
  grid->size          = 0;
  grid->xsize         = 0;
  grid->ysize         = 0;
  grid->xpole         = 0;
  grid->ypole         = 0;
  grid->prec          = 0;
  grid->isRotated     = FALSE;
  grid->ntr           = 0;
  grid->nvertex       = 0;

  grid->originLon     = 0;
  grid->originLat     = 0;
  grid->lonParY       = 0;
  grid->lat1          = 0;
  grid->lat2          = 0;
  grid->def_originLon = FALSE;
  grid->def_originLat = FALSE;
  grid->def_lonParY   = FALSE;
  grid->def_lat1      = FALSE;
  grid->def_lat2      = FALSE;

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
  grid->xname[0]      = 0;
  grid->xlongname[0]  = 0;
  grid->xunits[0]     = 0;
  grid->yname[0]      = 0;
  grid->ylongname[0]  = 0;
  grid->yunits[0]     = 0;
  grid->nd = 0;
  grid->ni = 0;
  grid->ni2 = 0;
  grid->ni3 = 0;
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

      strncpy(optname, pname, namelen);
      optname[namelen] = '\0';
    }
  else
    nerr = 1;

  return (nerr);
}


int gridDefine(GRID grid)
{
  static char func[] = "gridDefine";
  int gridID = UNDEFID;
  int i;

  switch ( grid.type )
    {
    case GRID_GENERIC:
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
      {
	if ( grid.xsize == 0 ) Error(func, "xsize undefined!");
	if ( grid.ysize == 0 ) Error(func, "ysize undefined!");

	if ( grid.size == 0 ) grid.size = grid.xsize*grid.ysize;

	gridID = gridCreate(grid.type, grid.size);

	gridDefPrec(gridID, grid.prec);

	gridDefXsize(gridID, grid.xsize);
	gridDefYsize(gridID, grid.ysize);

	if ( (grid.def_xfirst || grid.def_xlast || grid.def_xinc) && grid.xvals == NULL )
	  {
	    grid.xvals = (double *) malloc(grid.xsize*sizeof(double));
	    gridGenXvals(grid.xsize, grid.xfirst, grid.xlast, grid.xinc, grid.xvals);
	    /*
	    if ( grid.type == GRID_GAUSSIAN && grid.xbounds == NULL && grid.xsize > 1 )
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
	    */
	  }

	if ( (grid.def_yfirst || grid.def_ylast || grid.def_yinc) && grid.yvals == NULL )
	  {
	    if ( ! grid.def_ylast ) grid.ylast = grid.yfirst;
	    grid.yvals = (double *) malloc(grid.ysize*sizeof(double));
	    gridGenYvals(grid.type, grid.ysize, grid.yfirst, grid.ylast, grid.yinc, grid.yvals);
	    /*
	    if ( grid.type == GRID_GAUSSIAN && grid.ybounds == NULL && grid.ysize > 1 )
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
	    */
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
	    printf("gridDefXbounds\n");
	    gridDefXbounds(gridID, grid.xbounds);
	    free(grid.xbounds);
	  }

	if ( grid.ybounds )
	  {
	    printf("gridDefYbounds\n");
	    gridDefYbounds(gridID, grid.ybounds);
	    free(grid.ybounds);
	  }

	if ( grid.isRotated )
	  {
	    gridDefXpole(gridID, grid.xpole);
	    gridDefYpole(gridID, grid.ypole);
	  }
	break;
      }
    case GRID_CURVILINEAR:
    case GRID_CELL:
      {
	if ( grid.size == 0 ) grid.size = grid.xsize*grid.ysize;

	gridID = gridCreate(grid.type, grid.size);

	gridDefPrec(gridID, grid.prec);

	if ( grid.type == GRID_CURVILINEAR )
	  {
	    if ( grid.xsize == 0 ) Error(func, "xsize undefined!");
	    if ( grid.ysize == 0 ) Error(func, "ysize undefined!");
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

	break;
      }
    case GRID_LAMBERT:
      {
	if ( grid.xsize == 0 ) Error(func, "xsize undefined!");
	if ( grid.ysize == 0 ) Error(func, "ysize undefined!");

	if ( grid.size == 0 ) grid.size = grid.xsize*grid.ysize;

	gridID = gridCreate(grid.type, grid.size);

	gridDefPrec(gridID, grid.prec);

	gridDefXsize(gridID, grid.xsize);
	gridDefYsize(gridID, grid.ysize);

	if ( grid.def_originLon == FALSE ) Error(func, "originLon undefined!");
	if ( grid.def_originLat == FALSE ) Error(func, "originLat undefined!");
	if ( grid.def_lonParY   == FALSE ) Error(func, "lonParY undefined!");
	if ( grid.def_lat1      == FALSE ) Error(func, "lat1 undefined!");
	if ( grid.def_lat2      == FALSE ) Error(func, "lat2 undefined!");
	if ( grid.def_xinc      == FALSE ) Error(func, "xinc undefined!");
	if ( grid.def_yinc      == FALSE ) Error(func, "yinc undefined!");

	gridDefLambert(gridID, grid.originLon, grid.originLat, grid.lonParY,
		       grid.lat1, grid.lat2, grid.xinc, grid.yinc);

	break;
      }
    case GRID_SPECTRAL:
      {
	if ( grid.ntr == 0 )
	  Error(func, "truncation undefined!");
	if ( grid.size == 0 )
	  grid.size = (grid.ntr+1) * (grid.ntr+2);

	gridID = gridCreate(grid.type, grid.size);

	gridDefPrec(gridID, grid.prec);

	gridDefTrunc(gridID, grid.ntr);
	
	break;
      }
    case GRID_GME:
      {
	if ( grid.nd == 0 ) Error(func, "nd undefined!");
	if ( grid.ni == 0 ) Error(func, "ni undefined!");
	if ( grid.size == 0 ) Error(func, "size undefined!");

	gridID = gridCreate(grid.type, grid.size);

	gridDefPrec(gridID, grid.prec);

	gridDefGMEnd(gridID, grid.nd);
	gridDefGMEni(gridID, grid.ni);
	gridDefGMEni2(gridID, grid.ni2);
	gridDefGMEni3(gridID, grid.ni3);
	
	break;
      }
    default:
      {
	if ( grid.type == -1 )
	  Error(func, "gridtype undefined!");
	else
	  Error(func, "%s grid unsupported!", gridNamePtr(grid.type));
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

static char *skipSeparator(char *pline)
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
double *readfield(GRID *grid, int record, char *format, char *filename)
{
  static char func[] = "readfield";
  int fileID, rxysize, ierr, irec;
  double *vals;

  if ( grid->size == 0 )  Error(func, "grid size = 0!");
  if ( format == NULL )   Error(func, "format undefined!");
  if ( filename == NULL ) Error(func, "file name undefined!");

  vals = (double *) malloc(grid->size*sizeof(double));

  if ( strcmp(format, "extra") == 0 )
    {
      int header[4];
      fileID = extOpen(filename, "r");
      if ( fileID == UNDEFID ) SysError(func, filename);

      for ( irec = 0; irec < record; irec++ )
	{
	  ierr = extReadHeader(fileID, header);
	  if ( ierr <= 0 ) Error(func, "Record %d unexpected EOF in file %s", irec+1, filename);
	}
      grid->prec   = extInqPrec(fileID);
      rxysize = header[3];
      if ( rxysize != (int) grid->size ) Error(func, "unexpected record size of %d!", rxysize);
      ierr = extReadDataDP(fileID, vals);
      extClose(fileID);
    }
  else
    Error(func, "format %s unsupported!", format);

  return (vals);
}
*/
/*
double *readfield4(GRID *grid, int record, char *format, char *filename)
{
  static char func[] = "readfield4";
  int fileID, rxysize, ierr, irec;
  double *vals;

  if ( grid->size == 0 )  Error(func, "grid size = 0!");
  if ( format == NULL )   Error(func, "format undefined!");
  if ( filename == NULL ) Error(func, "file name undefined!");

  vals  = (double *) malloc(4*grid->size*sizeof(double));

  if ( strcmp(format, "extra") == 0 )
    {
      int header[4];
      fileID = extOpen(filename, "r");
      if ( fileID == UNDEFID ) SysError(func, filename);

      for ( irec = 0; irec < record; irec++ )
	{
	  ierr = extReadHeader(fileID, header);
	  if ( ierr <= 0 ) Error(func, "Record %d unexpected EOF in file %s", irec+1, filename);
	}
      grid->prec   = extInqPrec(fileID);
      rxysize = header[3];
      if ( rxysize != (int) (4*grid->size) ) Error(func, "unexpected record size of %d!", rxysize);
      ierr = extReadDataDP(fileID, vals);

      extClose(fileID);
    }
  else
    Error(func, "format %s unsupported!", format);

  return (vals);
}
*/

int gridFromFile(FILE *gfp, const char *dname)
{
  static char func[] = "gridFromFile";
  char line[MAX_LINE_LEN], *pline;
  /* char path[4096]; */
  int gridID = -1;
  int size;
  GRID grid;

  gridInit(&grid);

  while ( readline(gfp, line, MAX_LINE_LEN) )
    {
      if ( line[0] == '#' ) continue;
      if ( line[0] == '\0' ) continue;
      pline = line;
      while ( isspace((int) *pline) ) pline++;
      if ( pline[0] == '\0' ) continue;
      if ( strncmp(pline, "gridtype", 8) == 0 )
	{
	  pline = skipSeparator(pline + 8);
	  if ( strncmp(pline, "lonlat", 6)  == 0 ||
	       strncmp(pline, "latlon", 6)  == 0 )
	    {
	      grid.type = GRID_LONLAT;
	      grid.nvertex = 2;
	    }
	  else if ( strncmp(pline, "gaussian", 8)  == 0 )
	    {
	      grid.type = GRID_GAUSSIAN;
	      grid.nvertex = 2;
	    }
	  else if ( strncmp(pline, "curvilinear", 11)  == 0 )
	    {
	      grid.type = GRID_CURVILINEAR;
	      grid.nvertex = 4;
	    }
	  else if ( strncmp(pline, "cell", 4)  == 0 )
	    grid.type = GRID_CELL;
	  else if ( strncmp(pline, "gme", 3)  == 0 )
	    grid.type = GRID_GME;
	  else if ( strncmp(pline, "lambert", 7)  == 0 )
	    grid.type = GRID_LAMBERT;
	  else
	    Warning(func, "Invalid grid name : %s", pline);
	}
      else if ( strncmp(pline, "gridprec", 8)  == 0 )
	{
	  grid.prec = atol(skipSeparator(pline + 8));
	}
      else if ( strncmp(pline, "gridsize", 8)  == 0 )
	{
	  grid.size = atol(skipSeparator(pline + 8));
	}
      else if ( strncmp(pline, "xname", 5)  == 0 )
	{
	  strcpy(grid.xname, skipSeparator(pline + 5));
	}
      else if ( strncmp(pline, "xlongname", 9)  == 0 )
	{
	  strcpy(grid.xlongname, skipSeparator(pline + 9));
	}
      else if ( strncmp(pline, "xunits", 6)  == 0 )
	{
	  strcpy(grid.xunits, skipSeparator(pline + 6));
	}
      else if ( strncmp(pline, "yname", 5)  == 0 )
	{
	  strcpy(grid.yname, skipSeparator(pline + 5));
	}
      else if ( strncmp(pline, "ylongname", 9)  == 0 )
	{
	  strcpy(grid.ylongname, skipSeparator(pline + 9));
	}
      else if ( strncmp(pline, "yunits", 6)  == 0 )
	{
	  strcpy(grid.yunits, skipSeparator(pline + 6));
	}
      else if ( strncmp(pline, "nvertex", 7)  == 0 )
	{
	  grid.nvertex = atol(skipSeparator(pline + 7));
	}
      else if ( strncmp(pline, "ni", 2)  == 0 )
	{
	  grid.ni = atol(skipSeparator(pline + 2));
          grid.nd = 10;
	}
      else if ( strncmp(pline, "xsize", 5)  == 0 )
	{
	  grid.xsize = atol(skipSeparator(pline + 5));
	}
      else if ( strncmp(pline, "nlon", 4)  == 0 )
	{
	  grid.xsize = atol(skipSeparator(pline + 4));
	}
      else if ( strncmp(pline, "ysize", 5)  == 0 )
	{
	  grid.ysize = atol(skipSeparator(pline + 5));
	}
      else if ( strncmp(pline, "nlat", 4)  == 0 )
	{
	  grid.ysize = atol(skipSeparator(pline + 4));
	}
      else if ( strncmp(pline, "xfirst", 6)  == 0 )
	{
	  grid.xfirst = atof(skipSeparator(pline + 6));
	  grid.def_xfirst = TRUE;
	}
      else if ( strncmp(pline, "lonfirst", 8)  == 0 )
	{
	  grid.xfirst = atof(skipSeparator(pline + 8));
	  grid.def_xfirst = TRUE;
	}
      else if ( strncmp(pline, "yfirst", 6)  == 0 )
	{
	  grid.yfirst = atof(skipSeparator(pline + 6));
	  grid.def_yfirst = TRUE;
	}
      else if ( strncmp(pline, "latfirst", 8)  == 0 )
	{
	  grid.yfirst = atof(skipSeparator(pline + 8));
	  grid.def_yfirst = TRUE;
	}
      else if ( strncmp(pline, "xlast", 5)  == 0 )
	{
	  grid.xlast = atof(skipSeparator(pline + 5));
	  grid.def_xlast = TRUE;
	}
      else if ( strncmp(pline, "lonlast", 7)  == 0 )
	{
	  grid.xlast = atof(skipSeparator(pline + 7));
	  grid.def_xlast = TRUE;
	}
      else if ( strncmp(pline, "ylast", 5)  == 0 )
	{
	  grid.ylast = atof(skipSeparator(pline + 5));
	  grid.def_ylast = TRUE;
	}
      else if ( strncmp(pline, "latlast", 7)  == 0 )
	{
	  grid.ylast = atof(skipSeparator(pline + 7));
	  grid.def_ylast = TRUE;
	}
      else if ( strncmp(pline, "xinc", 4)  == 0 )
	{
	  grid.xinc = atof(skipSeparator(pline + 4));
	  grid.def_xinc = TRUE;
	}
      else if ( strncmp(pline, "loninc", 6)  == 0 )
	{
	  grid.xinc = atof(skipSeparator(pline + 6));
	  grid.def_xinc = TRUE;
	}
      else if ( strncmp(pline, "yinc", 4)  == 0 )
	{
	  grid.yinc = atof(skipSeparator(pline + 4));
	  grid.def_yinc = TRUE;
	}
      else if ( strncmp(pline, "latinc", 6)  == 0 )
	{
	  grid.yinc = atof(skipSeparator(pline + 6));
	  grid.def_yinc = TRUE;
	}
      else if ( strncmp(pline, "originLon", 9)  == 0 )
	{
	  grid.originLon = atof(skipSeparator(pline + 9));
	  grid.def_originLon = TRUE;
	}
      else if ( strncmp(pline, "originLat", 9)  == 0 )
	{
	  grid.originLat = atof(skipSeparator(pline + 9));
	  grid.def_originLat = TRUE;
	}
      else if ( strncmp(pline, "lonParY", 7)  == 0 )
	{
	  grid.lonParY = atof(skipSeparator(pline + 7));
	  grid.def_lonParY = TRUE;
	}
      else if ( strncmp(pline, "lat1", 4)  == 0 )
	{
	  grid.lat1 = atof(skipSeparator(pline + 4));
	  grid.def_lat1 = TRUE;
	}
      else if ( strncmp(pline, "lat2", 4)  == 0 )
	{
	  grid.lat2 = atof(skipSeparator(pline + 4));
	  grid.def_lat2 = TRUE;
	}
      else if ( strncmp(pline, "xnpole", 6)  == 0 )
	{
	  grid.xpole = atof(skipSeparator(pline + 6));
	  grid.isRotated = TRUE;
	}
      else if ( strncmp(pline, "lonpole", 7)  == 0 )
	{
	  grid.xpole = atof(skipSeparator(pline + 7));
	  grid.isRotated = TRUE;
	}
      else if ( strncmp(pline, "ynpole", 6)  == 0 )
	{
	  grid.ypole = atof(skipSeparator(pline + 6));
	  grid.isRotated = TRUE;
	}
      else if ( strncmp(pline, "latpole", 7)  == 0 )
	{
	  grid.ypole = atof(skipSeparator(pline + 7));
	  grid.isRotated = TRUE;
	}
      /*
      else if ( strncmp(pline, "xvalsf", 6)  == 0 )
	{
	  char *format = NULL;
	  char *file = NULL;
	  int record = 1;
	  pline = skipSeparator(pline + 6);
	  if ( strncmp(pline, "format", 6)  == 0 )
	    {
	      pline = skipSeparator(pline + 6);
	      format = pline;
	      while ( isalnum((int) *pline ) ) pline++;
	      *pline++ = 0;
	    }
	  pline = skipSeparator(pline);
	  if ( strncmp(pline, "record", 6)  == 0 )
	    {
	      pline = skipSeparator(pline + 6);
	      record = atoi(pline);
	      while ( isalnum((int) *pline ) ) pline++;
	      *pline++ = 0;
	    }
	  pline = skipSeparator(pline);
	  if ( strncmp(pline, "file", 4)  == 0 )
	    {
	      pline = skipSeparator(pline + 4);
	      file = pline;
	      while ( *pline != ' ' &&  *pline != 0 ) pline++;
	      *pline++ = 0;
	    }

	  if ( file == NULL ) Error(func, "file name undefined!");

	  if ( *file == '^' || *file == '$' )
	    fnmexp2(path, file, dname);
	  else
	    strcpy(path, file);

	  grid.xvals = readfield(&grid, record, format, path);
	}
      else if ( strncmp(pline, "yvalsf", 6)  == 0 )
	{
	  char *format = NULL;
	  char *file = NULL;
	  int record = 1;
	  pline = skipSeparator(pline + 6);
	  if ( strncmp(pline, "format", 6)  == 0 )
	    {
	      pline = skipSeparator(pline + 6);
	      format = pline;
	      while ( isalnum((int) *pline ) ) pline++;
	      *pline++ = 0;
	    }
	  pline = skipSeparator(pline);
	  if ( strncmp(pline, "record", 6)  == 0 )
	    {
	      pline = skipSeparator(pline + 6);
	      record = atoi(pline);
	      while ( isalnum((int) *pline ) ) pline++;
	      *pline++ = 0;
	    }
	  pline = skipSeparator(pline);
	  if ( strncmp(pline, "file", 4)  == 0 )
	    {
	      pline = skipSeparator(pline + 4);
	      file = pline;
	      while ( *pline != ' ' &&  *pline != 0 ) pline++;
	      *pline++ = 0;
	    }

	  if ( file == NULL ) Error(func, "file name undefined!");

	  if ( *file == '^' || *file == '$' )
	    fnmexp2(path, file, dname);
	  else
	    strcpy(path, file);

	  grid.yvals = readfield(&grid, record, format, path);
	}
      else if ( strncmp(pline, "gridareaf", 9)  == 0 )
	{
	  char *format = NULL;
	  char *file = NULL;
	  int record = 1;
	  pline = skipSeparator(pline + 9);
	  if ( strncmp(pline, "format", 6)  == 0 )
	    {
	      pline = skipSeparator(pline + 6);
	      format = pline;
	      while ( isalnum((int) *pline ) ) pline++;
	      *pline++ = 0;
	    }
	  pline = skipSeparator(pline);
	  if ( strncmp(pline, "record", 6)  == 0 )
	    {
	      pline = skipSeparator(pline + 6);
	      record = atoi(pline);
	      while ( isalnum((int) *pline ) ) pline++;
	      *pline++ = 0;
	    }
	  pline = skipSeparator(pline);
	  if ( strncmp(pline, "file", 4)  == 0 )
	    {
	      pline = skipSeparator(pline + 4);
	      file = pline;
	      while ( *pline != ' ' &&  *pline != 0 ) pline++;
	      *pline++ = 0;
	    }

	  if ( file == NULL ) Error(func, "file name undefined!");

	  if ( *file == '^' || *file == '$' )
	    fnmexp2(path, file, dname);
	  else
	    strcpy(path, file);

	  grid.area = readfield(&grid, record, format, path);
	}
      else if ( strncmp(pline, "xboundsf", 8)  == 0 )
	{
	  char *format = NULL;
	  char *file = NULL;
	  int record = 1;
	  pline = skipSeparator(pline + 8);
	  if ( strncmp(pline, "format", 6)  == 0 )
	    {
	      pline = skipSeparator(pline + 6);
	      format = pline;
	      while ( isalnum((int) *pline ) ) pline++;
	      *pline++ = 0;
	    }
	  pline = skipSeparator(pline);
	  if ( strncmp(pline, "record", 6)  == 0 )
	    {
	      pline = skipSeparator(pline + 6);
	      record = atoi(pline);
	      while ( isalnum((int) *pline ) ) pline++;
	      *pline++ = 0;
	    }
	  pline = skipSeparator(pline);
	  if ( strncmp(pline, "file", 4)  == 0 )
	    {
	      pline = skipSeparator(pline + 4);
	      file = pline;
	      while ( *pline != ' ' &&  *pline != 0 ) pline++;
	      *pline++ = 0;
	    }

	  if ( file == NULL ) Error(func, "file name undefined!");

	  if ( *file == '^' || *file == '$' )
	    fnmexp2(path, file, dname);
	  else
	    strcpy(path, file);

	  grid.xbounds = readfield4(&grid, record, format, path);
	}
      else if ( strncmp(pline, "yboundsf", 8)  == 0 )
	{
	  char *format = NULL;
	  char *file = NULL;
	  int record = 1;
	  pline = skipSeparator(pline + 8);
	  if ( strncmp(pline, "format", 6)  == 0 )
	    {
	      pline = skipSeparator(pline + 6);
	      format = pline;
	      while ( isalnum((int) *pline ) ) pline++;
	      *pline++ = 0;
	    }
	  pline = skipSeparator(pline);
	  if ( strncmp(pline, "record", 6)  == 0 )
	    {
	      pline = skipSeparator(pline + 6);
	      record = atoi(pline);
	      while ( isalnum((int) *pline ) ) pline++;
	      *pline++ = 0;
	    }
	  pline = skipSeparator(pline);
	  if ( strncmp(pline, "file", 4)  == 0 )
	    {
	      pline = skipSeparator(pline + 4);
	      file = pline;
	      while ( *pline != ' ' &&  *pline != 0 ) pline++;
	      *pline++ = 0;
	    }

	  if ( file == NULL ) Error(func, "file name undefined!");

	  if ( *file == '^' || *file == '$' )
	    fnmexp2(path, file, dname);
	  else
	    strcpy(path, file);

	  grid.ybounds = readfield4(&grid, record, format, path);
	}
      */
      else if ( strncmp(pline, "gridlatlon", 10)  == 0 )
	{
	  int i;
	  float flat, flon;
	  if ( grid.size == 0 ) grid.size = grid.xsize * grid.ysize;
	  
	  grid.xvals = (double *) malloc(grid.size*sizeof(double));
	  grid.yvals = (double *) malloc(grid.size*sizeof(double));
	  for ( i = 0; i < (int) grid.size; i++ )
	    {
	      if ( ! readline(gfp, line, MAX_LINE_LEN) )
		{
		  Warning(func, "Incomplete command: >gridlatlon<");
		  break;
		}
	      sscanf(line, "%g %g", &flat, &flon);
	      grid.yvals[i] = flat;
	      grid.xvals[i] = flon;
	    }
	}
      else if ( strncmp(pline, "xvals", 5)  == 0 )
	{
	  int i = 0;
	  double fval;
	  char *endptr;

	  if ( grid.type == GRID_CURVILINEAR || grid.type == GRID_CELL )
	    size = grid.size;
	  else
	    size = grid.xsize;

	  if ( size > 0 )
	    {
	      pline = skipSeparator(pline + 5);
	      grid.xvals = (double *) malloc(size*sizeof(double));

	      for ( i = 0; i < size; i++ )
		{
		  endptr = pline;
		  fval = strtod(pline, &endptr);
		  if ( pline == endptr )
		    {
		      if ( ! readline(gfp, line, MAX_LINE_LEN) )
			{
			  Warning(func, "Incomplete command: >xvals<");
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
	    Warning(func, "xsize or gridsize undefined!");
	}
      else if ( strncmp(pline, "yvals", 5)  == 0 )
	{
	  int i = 0;
	  double fval;
	  char *endptr;

	  if ( grid.type == GRID_CURVILINEAR || grid.type == GRID_CELL )
	    size = grid.size;
	  else
	    size = grid.ysize;

	  if ( size > 0 )
	    {
	      pline = skipSeparator(pline + 5);
	      grid.yvals = (double *) malloc(size*sizeof(double));

	      for ( i = 0; i < size; i++ )
		{
		  endptr = pline;
		  fval = strtod(pline, &endptr);
		  if ( pline == endptr )
		    {
		      if ( ! readline(gfp, line, MAX_LINE_LEN) )
			{
			  Warning(func, "Incomplete command: >yvals<");
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
	    Warning(func, "ysize or gridsize undefined!");
	}
      else if ( strncmp(pline, "xbounds", 7)  == 0 )
	{
	  int i = 0;
	  double fval;
	  char *endptr;

	  if ( grid.nvertex == 0 )
	    {
	      if ( grid.type == GRID_CURVILINEAR ) grid.nvertex = 4;
	    }

	  if ( grid.type == GRID_CURVILINEAR || grid.type == GRID_CELL )
	    size = grid.size;
	  else
	    size = grid.xsize;

	  if ( size > 0 && grid.nvertex > 0 )
	    {	  
	      pline = skipSeparator(pline + 7);
	      grid.xbounds = (double *) malloc(size*grid.nvertex*sizeof(double));

	      for ( i = 0; i < (int) (size*grid.nvertex); i++ )
		{
		  endptr = pline;
		  fval = strtod(pline, &endptr);
		  if ( pline == endptr )
		    {
		      if ( ! readline(gfp, line, MAX_LINE_LEN) )
			{
			  Warning(func, "Incomplete command: >xbounds<");
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
	      if ( size         == 0 ) Warning(func, "xsize or gridsize undefined!");
	      if ( grid.nvertex == 0 ) Warning(func, "nvertex undefined!");
	    }
	}
      else if ( strncmp(pline, "ybounds", 7)  == 0 )
	{
	  int i = 0;
	  double fval;
	  char *endptr;

	  if ( grid.nvertex == 0 )
	    {
	      if ( grid.type == GRID_CURVILINEAR ) grid.nvertex = 4;
	    }

	  if ( grid.type == GRID_CURVILINEAR || grid.type == GRID_CELL )
	    size = grid.size;
	  else
	    size = grid.ysize;

	  if ( size > 0 && grid.nvertex > 0 )
	    {	  
	      pline = skipSeparator(pline + 7);
	      grid.ybounds = (double *) malloc(size*grid.nvertex*sizeof(double));

	      for ( i = 0; i < (int) (size*grid.nvertex); i++ )
		{
		  endptr = pline;
		  fval = strtod(pline, &endptr);
		  if ( pline == endptr )
		    {
		      if ( ! readline(gfp, line, MAX_LINE_LEN) )
			{
			  Warning(func, "Incomplete command: >ybounds<");
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
	      if ( grid.ysize   == 0 ) Warning(func, "ysize or gridsize undefined!");
	      if ( grid.nvertex == 0 ) Warning(func, "nvertex undefined!");
	    }
	}
      else
	{
	  if ( grid.type != UNDEFID )
	    Warning(func, "Invalid grid command : >%s<", pline);
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
  static char func[] = "gridFromPingo";
  int gridID = -1;
  int i;
  int nlon, nlat;
  int lgauss = FALSE;
  GRID grid;

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
	    if ( DBL_IS_EQUAL(grid.xvals[0], grid.xvals[1]) )
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
		  Message(func, "Longitudes are not in ascending order!");
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
	  Message(func, "Latitudes must be between 90 and -90!");
	  return (gridID);
	}

      for ( i = 0; i < nlat - 1; i++ )
	if ( DBL_IS_EQUAL(grid.yvals[i+1], grid.yvals[i]) || (i < nlat - 2 &&
	    ((grid.yvals[i+1] > grid.yvals[i]) != (grid.yvals[i+2] > grid.yvals[i+1]))) )
	  {
	    Message(func, "Latitudes must be in descending or ascending order!");
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
  static char func[] = "ntr2nlat";
  int nlat, nlat2;

  nlat = NINT((ntr*3.+1.)/2.);
  if ( (nlat % 2) > 0 )
    {
      nlat  = nlat + 1;
      nlat2 = NINT(((ntr+1)*3.+1.)/2.);
      /*
      if ( nlat == nlat2 )
	Error(func, "Computation of latitudes failed for truncation %d", ntr);
      */
    }

  return (nlat);
}


int ntr2nlat_linear(int ntr)
{
  static char func[] = "ntr2nlat_linear";
  int nlat, nlat2;

  nlat = NINT((ntr*2.+1.)/2.);
  if ( (nlat % 2) > 0 )
    {
      nlat  = nlat + 1;
      nlat2 = NINT(((ntr+1)*2.+1.)/2.);
      /*
      if ( nlat == nlat2 )
	Error(func, "Computation of latitudes failed for truncation %d", ntr);
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

static void factorni(int kni, int *kni2, int *kni3)
{
  /**********************************************************************/
  /* factorni computes the factors of the integer input kni, assuming   */ 
  /* that kni decomposes into kni3 factors (kni3 either 0 or 1) of "3"  */
  /* and kni2 (kni2 > 0) factors of "2". The subroutine returns the     */
  /* number of factors of "3", kni3, number of factors of "2", kni2.    */
  /* Bails out in case of error.                                        */
  /**********************************************************************/
  /*  Author: D. Majewski, DWD, January 2000                            */
  /**********************************************************************/ 
  /*  Input                                                             */
  /*  kni     INT   number of intervals on a main triangle side         */
  /**********************************************************************/
  /*  Output                                                            */
  /*  kni2    INT      exponent of "2", kni2 > 0                        */ 
  /*  kni3    INT      exponent of "3", either 0 or 1                   */
  /**********************************************************************/

  int mx;

  mx    = kni;

  *kni2  = 0;
  *kni3  = 0;

  while (mx > 1)
    { 
      if (mx%2 == 0) 
	{
          *kni2 = *kni2 + 1;
          mx   = mx/2;
	}
      else if (mx%3 == 0) 
	{
          *kni3 = *kni3 + 1;
          mx   = mx/3;
	}
      else
	{
	  /* error return */
        }
    }

  /* kni3 must not be greater than */

  if (*kni3 > 1) 
    {
      /* error return */
    }
}


int gridFromName(const char *gridname)
{
  const char *pline;
  int gridID = UNDEFID;
  GRID grid;

  gridInit(&grid);

  if ( gridname[0] == 't' && gridname[1] == 'l' ) /* tl<RES>grid or tl<RES>spec */
    {
      pline = &gridname[2];
      if ( isdigit((int) *pline) )
	{
	  grid.ntr = atoi(pline);
	  while ( isdigit((int) *pline) ) pline++;
	  if      ( strncmp(pline, "grid", 4) == 0 ) grid.type = GRID_GAUSSIAN;
	  else if ( strncmp(pline, "spec", 4) == 0 ) grid.type = GRID_SPECTRAL;
      
	  grid.ysize = ntr2nlat_linear(grid.ntr);
	  grid.xsize = compNlon(grid.ysize);

	  if ( grid.type == GRID_GAUSSIAN )
	    {
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
	  if      ( strncmp(pline, "grid", 4) == 0 ) grid.type = GRID_GAUSSIAN;
	  else if ( strncmp(pline, "spec", 4) == 0 ) grid.type = GRID_SPECTRAL;
      
	  grid.ysize = ntr2nlat(grid.ntr);
	  grid.xsize = compNlon(grid.ysize);

	  if ( grid.type == GRID_GAUSSIAN )
	    {
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
  else if ( gridname[0] == 'g' ) /* g<LON>x<LAT> */
    {
      pline = &gridname[1];
      if ( isdigit((int) *pline) )
	{
	  grid.type = GRID_GENERIC;
	  grid.xsize = atoi(pline);
	  while ( isdigit((int) *pline) ) pline++;
	  pline++;
	  grid.ysize = atoi(pline);
	  while ( isdigit((int) *pline) ) pline++;
	}
    }

  if ( grid.type != -1 )
    gridID = gridDefine(grid);

  return (gridID);
}


int cdoDefineGrid(const char *gridfile)
{
  static char func[] = "cdoDefineGrid";
  FILE *gfp;
  char buffer[4];
  int gridID = -1;

  gfp = fopen(gridfile, "r");
  if ( gfp == NULL )
    {
      gridID = gridFromName(gridfile);

      if ( gridID == -1 ) cdoAbort("Open failed on %s!", gridfile);
    }
  else
    {
      if ( fread(buffer, 1, 4, gfp) != 4 )
	SysError(func, "Read grid from %s failed!", gridfile);

      fclose(gfp);

      if ( strncmp(buffer, "CDF", 3) == 0 )
	{
	  if ( cdoDebug ) cdoPrint("Grid from netCDF file");
	  gridID = gridFromNCfile(gridfile);
	}

      if ( gridID == -1 )
	{
	  if ( strncmp(buffer+1, "HDF", 3) == 0 )
	    {
	      if ( cdoDebug ) cdoPrint("Grid from HDF5 file");
	      gridID = gridFromH5file(gridfile);
	    }
	}

      if ( gridID == -1 )
	{
	  int streamID;
	  if ( cdoDebug ) cdoPrint("Grid from CDI file");
	  streamID = streamOpenRead(gridfile);
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


struct cart {
  double x[3];
};

struct geo {
  double lon;
  double lat;
};

struct cart gc2cc(struct geo *position);
double areas(struct cart *dv1, struct cart *dv2, struct cart *dv3);

int gridGenArea(int gridID, double *area)
{
  static char func[] = "gridArea";
  int status = FALSE;
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
       gridtype != GRID_LAMBERT     &&
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
	  status = TRUE;
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

	  area[i] += areas(&c1, &c2, &c3);
	}

      p1.lon = grid_corner_lon[i*nv+0]*deg2rad; 
      p1.lat = grid_corner_lat[i*nv+0]*deg2rad;
      c1 = gc2cc(&p1);
      p2.lon = grid_corner_lon[i*nv+nv-1]*deg2rad; 
      p2.lat = grid_corner_lat[i*nv+nv-1]*deg2rad;
      c2 = gc2cc(&p2);

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
  int status = FALSE;
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


int gridWeights(int gridID, double *grid_wgts)
{
  static char func[] = "gridWeights";
  int i, gridsize, gridtype;
  int a_status, w_status;
  double *grid_area;

  gridtype = gridInqType(gridID);
  gridsize = gridInqSize(gridID);
  
  grid_area = (double *) malloc(gridsize*sizeof(double));

  a_status = FALSE;

  if ( gridHasArea(gridID) )
    {
      if ( cdoVerbose ) cdoPrint("Using existing grid cell area!");
      gridInqArea(gridID, grid_area);
    }
  else
    {
      if ( gridtype != GRID_LONLAT      &&
	   gridtype != GRID_GAUSSIAN    &&
	   gridtype != GRID_LAMBERT     &&
	   gridtype != GRID_GME         &&
	   gridtype != GRID_CURVILINEAR &&
	   gridtype != GRID_CELL )
	{
	  a_status = TRUE;
	}
      else
	{
	  a_status = gridGenArea(gridID, grid_area);
	}
    }

  if ( a_status == FALSE )
    {
      w_status = gridGenWeights(gridID, grid_area, grid_wgts);
    }
  else
    {
      for ( i = 0; i < gridsize; ++i )
	grid_wgts[i] = 1./gridsize;

      w_status = TRUE;
    }

  free(grid_area);

  return (w_status);
}

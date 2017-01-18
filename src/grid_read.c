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

#include <cdi.h>
#include "cdi_uuid.h"
#include "cdo_int.h"
#include "griddes.h"


#define MAX_LINE_LEN 65536

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

#define TEST_NEWFORMAT

#ifdef TEST_NEWFORMAT

typedef struct {
  keyValues_t *kv;
  bool isValid;
} kvmap_t;


int grid_read(FILE *gfp, const char *dname)
{
  char uuidOfHGridStr[256];

  list_t *pmlist = namelist_to_pmlist(gfp, dname);
  if ( pmlist == NULL ) return -1;
  list_t *kvlist = *(list_t **)pmlist->head->data;
  if ( kvlist == NULL ) return -1;

  griddes_t grid;
  gridInit(&grid);

  size_t nkv = list_size(kvlist);
  if ( nkv == 0 ) return -1;
  kvmap_t *kvmap = (kvmap_t*) Malloc(nkv*sizeof(kvmap_t));
  for ( size_t i = 0; i < nkv; ++i ) kvmap[i].isValid = false;

  size_t ik = 0;
  for ( listNode_t *kvnode = kvlist->head; kvnode; kvnode = kvnode->next )
    {
      keyValues_t *kv = *(keyValues_t **)kvnode->data;
      if ( ik == 0 && !STR_IS_EQ(kv->key, "gridtype") )
        cdoAbort("First grid description parameter must be >gridtype< (found: %s)!", kv->key);

      if ( kv->nvalues == 0 )
        {
          cdoWarning("Grid description parameter %s has no values, skipped!", kv->key);
        }
      else
        {
          kvmap[ik].isValid = true;
          kvmap[ik].kv = kv;
        }
      ik++;
    }

  size_t igmap = 0;
  for ( size_t ik = 0; ik < nkv; ++ik )
    {
      if ( !kvmap[ik].isValid ) continue;

      keyValues_t *kv = kvmap[ik].kv;
      const char *key = kv->key;
      size_t nvalues = kv->nvalues;
      const char *value = (kv->nvalues > 0) ? kv->values[0] : NULL;
      //bool lv1 = (kv->nvalues == 1);

      // printf("%s = ", key); if ( value  ) printf("%s", kv->value); printf("\n");

      if ( STR_IS_EQ(key, "gridtype") )
        {
          const char *gridtype = parameter2word(value);

          if      ( STR_IS_EQ(gridtype, "lonlat") )       grid.type = GRID_LONLAT;
          else if ( STR_IS_EQ(gridtype, "latlon") )       grid.type = GRID_LONLAT;
          else if ( STR_IS_EQ(gridtype, "gaussian") )     grid.type = GRID_GAUSSIAN;
          else if ( STR_IS_EQ(gridtype, "curvilinear") )  grid.type = GRID_CURVILINEAR;
          else if ( STR_IS_EQ(gridtype, "unstructured") ) grid.type = GRID_UNSTRUCTURED;
          else if ( STR_IS_EQ(gridtype, "cell") )         grid.type = GRID_UNSTRUCTURED;
          else if ( STR_IS_EQ(gridtype, "spectral") )     grid.type = GRID_SPECTRAL;
          else if ( STR_IS_EQ(gridtype, "gme") )          grid.type = GRID_GME;
          else if ( STR_IS_EQ(gridtype, "lcc") )          grid.type = GRID_LCC;
          else if ( STR_IS_EQ(gridtype, "lambert") )      grid.type = GRID_LCC;
          else if ( STR_IS_EQ(gridtype, "projection") )   grid.type = GRID_PROJECTION;
          else if ( STR_IS_EQ(gridtype, "generic") )      grid.type = GRID_GENERIC;
	  else cdoAbort("Invalid gridtype : %s (grid description file: %s)", gridtype, dname);
            
          if ( grid.type == GRID_LONLAT || grid.type == GRID_GAUSSIAN ) grid.nvertex = 2;
          else if ( grid.type == GRID_CURVILINEAR ) grid.nvertex = 4;
        }
      else if ( STR_IS_EQ(key, "gridprec") )       grid.prec = parameter2int(value);
      else if ( STR_IS_EQ(key, "gridsize") )       grid.size = parameter2int(value);
      else if ( STR_IS_EQ(key, "xsize") )          grid.xsize = parameter2int(value);
      else if ( STR_IS_EQ(key, "nlon") )           grid.xsize = parameter2int(value);
      else if ( STR_IS_EQ(key, "ysize") )          grid.ysize = parameter2int(value);
      else if ( STR_IS_EQ(key, "nlat") )           grid.ysize = parameter2int(value);
      else if ( STR_IS_EQ(key, "truncation") )     grid.ntr = parameter2int(value);
      else if ( STR_IS_EQ(key, "np") )             grid.np = parameter2int(value);
      else if ( STR_IS_EQ(key, "complexpacking") ) grid.lcomplex = parameter2int(value);
      else if ( STR_IS_EQ(key, "nvertex") )        grid.nvertex = parameter2int(value);
      else if ( STR_IS_EQ(key, "ni") )           { grid.ni = parameter2int(value); grid.nd = 10; }
      else if ( STR_IS_EQ(key, "position") )       grid.position = parameter2int(value);
      else if ( STR_IS_EQ(key, "number") )         grid.number = parameter2int(value);
      else if ( STR_IS_EQ(key, "xname") )     strcpy(grid.xname, parameter2word(value));
      else if ( STR_IS_EQ(key, "yname") )     strcpy(grid.yname, parameter2word(value));
      else if ( STR_IS_EQ(key, "xdimname") )  strcpy(grid.xdimname, parameter2word(value));
      else if ( STR_IS_EQ(key, "ydimname") )  strcpy(grid.ydimname, parameter2word(value));
      else if ( STR_IS_EQ(key, "vdimname") )  strcpy(grid.vdimname, parameter2word(value));
      else if ( STR_IS_EQ(key, "xlongname") ) strcpy(grid.xlongname, value);
      else if ( STR_IS_EQ(key, "ylongname") ) strcpy(grid.ylongname, value);
      else if ( STR_IS_EQ(key, "xunits") )    strcpy(grid.xunits, value);
      else if ( STR_IS_EQ(key, "yunits") )    strcpy(grid.yunits, value);
      else if ( STR_IS_EQ(key, "path") )      strcpy(grid.path, value);
      else if ( STR_IS_EQ(key, "uuid") )      { strcpy(uuidOfHGridStr, value); cdiStr2UUID(uuidOfHGridStr, grid.uuid); }
      else if ( STR_IS_EQ(key, "xfirst") )    { grid.xfirst = parameter2double(value); grid.def_xfirst = true; }
      else if ( STR_IS_EQ(key, "yfirst") )    { grid.yfirst = parameter2double(value); grid.def_yfirst = true; }
      else if ( STR_IS_EQ(key, "xlast") )     { grid.xlast = parameter2double(value); grid.def_xlast = true; }
      else if ( STR_IS_EQ(key, "ylast") )     { grid.ylast = parameter2double(value); grid.def_ylast = true; }
      else if ( STR_IS_EQ(key, "xinc") )      { grid.xinc = parameter2double(value); grid.def_xinc = true; }
      else if ( STR_IS_EQ(key, "yinc") )      { grid.yinc = parameter2double(value); grid.def_yinc = true; }
      else if ( STR_IS_EQ(key, "originLon") ) { grid.originLon = parameter2double(value); grid.def_originLon = true; }
      else if ( STR_IS_EQ(key, "originLat") ) { grid.originLat = parameter2double(value); grid.def_originLat = true; }
      else if ( STR_IS_EQ(key, "lonParY") )   { grid.lonParY = parameter2double(value); grid.def_lonParY = true; }
      else if ( STR_IS_EQ(key, "lat1") )      { grid.lat1 = parameter2double(value); grid.def_lat1 = true; }
      else if ( STR_IS_EQ(key, "lat2") )      { grid.lat2 = parameter2double(value); grid.def_lat2 = true; }
      else if ( STR_IS_EQ(key, "lat_0") )     { grid.lat_0 = parameter2double(value); grid.def_lat_0 = true; }
      else if ( STR_IS_EQ(key, "lat_1") )     { grid.lat_1 = parameter2double(value); grid.def_lat_1 = true; }
      else if ( STR_IS_EQ(key, "lat_2") )     { grid.lat_2 = parameter2double(value); grid.def_lat_2 = true; }
      else if ( STR_IS_EQ(key, "a") )         grid.a = parameter2double(value);
      else if ( STR_IS_EQ(key, "projection") )
        {
          if      ( STR_IS_EQ(value, "north") ) { grid.projflag = 0;    grid.scanflag = 64; }
	  else if ( STR_IS_EQ(value, "south") ) { grid.projflag = 128;  grid.scanflag = 64; }
	  else cdoAbort("Invalid projection : %s (grid description file: %s)", value, dname);
        }
      else if ( STR_IS_EQ(key, "xvals") )
        {
          size_t size = (grid.type == GRID_CURVILINEAR || grid.type == GRID_UNSTRUCTURED) ? grid.size : grid.xsize;
          if ( size == 0 ) cdoAbort("xsize or gridsize undefined (grid description file: %s)!", dname);
          if ( size != nvalues )
            cdoAbort("Number of xvals=%zu and size of xvals=%zu differ (grid description file: %s)!", nvalues, size, dname);

          grid.xvals = (double*) Malloc(size*sizeof(double));
          for ( size_t i = 0; i < size; ++i ) grid.xvals[i] = parameter2double(kv->values[i]);
        }
      else if ( STR_IS_EQ(key, "yvals") )
        {
          size_t size = (grid.type == GRID_CURVILINEAR || grid.type == GRID_UNSTRUCTURED) ? grid.size : grid.ysize;
          if ( size == 0 ) cdoAbort("ysize or gridsize undefined (grid description file: %s)!", dname);
          if ( size != nvalues )
            cdoAbort("Number of yvals=%zu and size of yvals=%zu differ (grid description file: %s)!", nvalues, size, dname);

          grid.yvals = (double*) Malloc(size*sizeof(double));
          for ( size_t i = 0; i < size; ++i ) grid.yvals[i] = parameter2double(kv->values[i]);
        }
      else if ( STR_IS_EQ(key, "xbounds") )
        {
          size_t size = (grid.type == GRID_CURVILINEAR || grid.type == GRID_UNSTRUCTURED) ? grid.size : grid.xsize;
          if ( size == 0 ) cdoAbort("xsize or gridsize undefined (grid description file: %s)!", dname);
          if ( grid.nvertex == 0 ) cdoAbort("nvertex undefined (grid description file: %s)!", dname);
          if ( grid.nvertex*size != nvalues )
            cdoAbort("Number of xbounds=%zu and size of xbounds=%zu differ (grid description file: %s)!", nvalues, grid.nvertex*size, dname);

          grid.xbounds = (double*) Malloc(grid.nvertex*size*sizeof(double));
          for ( size_t i = 0; i < grid.nvertex*size; ++i ) grid.xbounds[i] = parameter2double(kv->values[i]);
        }
      else if ( STR_IS_EQ(key, "ybounds") )
        {
          size_t size = (grid.type == GRID_CURVILINEAR || grid.type == GRID_UNSTRUCTURED) ? grid.size : grid.ysize;
          if ( size == 0 ) cdoAbort("ysize or gridsize undefined (grid description file: %s)!", dname);
          if ( grid.nvertex == 0 ) cdoAbort("nvertex undefined (grid description file: %s)!", dname);
          if ( grid.nvertex*size != nvalues )
            cdoAbort("Number of ybounds=%zu and size of ybounds=%zu differ (grid description file: %s)!", nvalues, grid.nvertex*size, dname);

          grid.ybounds = (double*) Malloc(grid.nvertex*size*sizeof(double));
          for ( size_t i = 0; i < grid.nvertex*size; ++i ) grid.ybounds[i] = parameter2double(kv->values[i]);
        }
      else if ( STR_IS_EQ(key, "gridlatlon") )
        {
	  if ( grid.size == 0 ) grid.size = grid.xsize * grid.ysize;
          if ( grid.size == 0 ) cdoAbort("gridsize undefined (grid description file: %s)!", dname);
          if ( (size_t)grid.size*2 != nvalues )
            cdoAbort("Number of gridlonlat values=%zu and size of grid=%zu differ (grid description file: %s)!", nvalues, grid.size*2, dname);
	  grid.xvals = (double*) Malloc(grid.size*sizeof(double));
	  grid.yvals = (double*) Malloc(grid.size*sizeof(double));
          for ( size_t i = 0; i < (size_t)grid.size; ++i )
            {
	      grid.yvals[i] = parameter2double(kv->values[2*i]);
	      grid.xvals[i] = parameter2double(kv->values[2*i+1]);
            }
        }
      else if ( STR_IS_EQ(key, "mask") )
        {
          size_t size = grid.size;
          if ( grid.size == 0 ) cdoAbort("gridsize undefined (grid description file: %s)!", dname);
          if ( size != nvalues )
            cdoAbort("Number of mask values=%zu and size of grid=%zu differ (grid description file: %s)!", nvalues, size, dname);
          grid.mask = (int*) Malloc(size*sizeof(int));
          size_t count = 0;
          for ( size_t i = 0; i < size; ++i )
            {
              grid.mask[i] = parameter2int(kv->values[i]);
              if ( grid.mask[i] == 1 ) count++;
            }
          if ( count == size ) { Free(grid.mask); grid.mask = NULL; }
        }
      else if ( STR_IS_EQ(key, "grid_mapping") ) { igmap = ik; break; }
      else cdoAbort("Invalid parameter : >%s< (grid description file: %s)", key, dname);
    }

  int gridID = (grid.type == CDI_UNDEFID) ? -1 : gridDefine(grid);

  if ( igmap > 0 )
    {
      for ( size_t ik = igmap; ik < nkv; ++ik )
        {
          if ( !kvmap[ik].isValid ) continue;

          keyValues_t *kv = kvmap[ik].kv;
          const char *key = kv->key;
          size_t nvalues = kv->nvalues;
          const char *value = (kv->nvalues > 0) ? kv->values[0] : NULL;
          // bool lv1 = (kv->nvalues == 1);

          // printf("%s = ", key); if ( value  ) printf("%s", kv->value); printf("\n");

          if ( STR_IS_EQ(key, "grid_mapping") )
            {
              cdiGridDefKeyStr(gridID, CDI_KEY_MAPNAME, (int)strlen(value)+1, value);
              continue;
            }
          if ( STR_IS_EQ(key, "grid_mapping_name") )
            cdiGridDefKeyStr(gridID, CDI_KEY_MAPPING, (int)strlen(value)+1, value);

          int dtype = literals_find_datatype(nvalues, kv->values);

          if ( dtype == CDI_DATATYPE_INT8 || dtype == CDI_DATATYPE_INT16 || dtype == CDI_DATATYPE_INT32 )
            {
              int *ivals = (int*) Malloc(nvalues*sizeof(int));
              for ( size_t i = 0; i < nvalues; ++i ) ivals[i] = literal_to_int(kv->values[i]);
              cdiDefAttInt(gridID, CDI_GLOBAL, key, dtype, nvalues, ivals);
              Free(ivals);
            }
          else if ( dtype == CDI_DATATYPE_FLT32 || dtype == CDI_DATATYPE_FLT64 )
            {
              double *dvals = (double*) Malloc(nvalues*sizeof(double));
              for ( size_t i = 0; i < nvalues; ++i ) dvals[i] = literal_to_double(kv->values[i]);
              cdiDefAttFlt(gridID, CDI_GLOBAL, key, dtype, nvalues, dvals);
              Free(dvals);
            }
          else
            {
              int len = (value && *value) ? (int) strlen(value) : 0;
              cdiDefAttTxt(gridID, CDI_GLOBAL, key, len, value);
            }
        }
    }

  list_destroy(pmlist);

  Free(kvmap);

  return gridID;
}

#else

static
char *skipSeparator(char *pline)
{
  while ( isspace((int) *pline) ) pline++;
  if ( *pline == '=' || *pline == ':' ) pline++;
  while ( isspace((int) *pline) ) pline++;

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


int grid_read(FILE *gfp, const char *dname)
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
	  else if ( cmpstrlen(pline, "unstructured", len)  == 0 )
	    grid.type = GRID_UNSTRUCTURED;
	  else if ( cmpstrlen(pline, "cell", len)  == 0 )
	    grid.type = GRID_UNSTRUCTURED;
	  else if ( cmpstrlen(pline, "spectral", len)  == 0 )
	    grid.type = GRID_SPECTRAL;
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
#endif

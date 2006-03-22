/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2006 Uwe Schulzweida, schulzweida@dkrz.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "error.h"

#define UNDEFID -1


typedef struct {
  double *vals;
  double *lbounds;
  double *ubounds;
  double *vct;
  int     vctsize;
  int     type;
  int     size;
  char    name[128];
  char    longname[256];
  char    units[128];
}
ZAXIS;


void zaxisInit(ZAXIS *zaxis)
{
  zaxis->vals        = NULL;
  zaxis->lbounds     = NULL;
  zaxis->ubounds     = NULL;
  zaxis->vct         = NULL;
  zaxis->type        = UNDEFID;
  zaxis->vctsize     = 0;
  zaxis->size        = 0;
  zaxis->name[0]     = 0;
  zaxis->longname[0] = 0;
  zaxis->units[0]    = 0;
}


static int getoptname(char *optname, const char *optstring, int nopt)
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


int zaxisDefine(ZAXIS zaxis)
{
  static char func[] = "zaxisDefine";
  int zaxisID = UNDEFID;

  if ( zaxis.type == -1 ) Error(func, "zaxistype undefined!");

  if ( zaxis.size == 0 ) Error(func, "zaxis size undefined!");

  zaxisID = zaxisNew(zaxis.type, zaxis.size);

  if ( zaxis.vals )
    {
      zaxisDefLevels(zaxisID, zaxis.vals);
      free(zaxis.vals);
    }
  if ( zaxis.lbounds )
    {
      zaxisDefLbounds(zaxisID, zaxis.lbounds);
      free(zaxis.lbounds);
    }
  if ( zaxis.ubounds )
    {
      zaxisDefLbounds(zaxisID, zaxis.ubounds);
      free(zaxis.ubounds);
    }

  if ( zaxis.name[0] )     zaxisDefName(zaxisID, zaxis.name);
  if ( zaxis.longname[0] ) zaxisDefLongname(zaxisID, zaxis.longname);
  if ( zaxis.units[0] )    zaxisDefUnits(zaxisID, zaxis.units);

  if ( zaxis.type == ZAXIS_HYBRID )
    {
      if ( zaxis.vctsize && zaxis.vct )
	zaxisDefVct(zaxisID, zaxis.vctsize, zaxis.vct);
      else
	Warning(func, "vct undefined!");	    
    }

  return (zaxisID);
}


static char *skipSeparator(char *pline)
{
  while ( isspace((int) *pline) ) pline++;
  if ( *pline == '=' || *pline == ':' ) pline++;
  while ( isspace((int) *pline) ) pline++;

  return (pline);
}


int zaxisFromFile(FILE *gfp)
{
  static char func[] = "zaxisFromFile";
  char line[1024], *pline;
  int zaxisID;
  ZAXIS zaxis;

  zaxisInit(&zaxis);

  while ( readline(gfp, line, 1024) )
    {
      if ( line[0] == '#' ) continue;
      if ( line[0] == '\0' ) continue;
      pline = line;
      while ( isspace((int) *pline) ) pline++;
      if ( pline[0] == '\0' ) continue;
      if ( strncmp(pline, "zaxistype", 9) == 0 || 
	   strncmp(pline, "type", 4) == 0 )
	{
	  if ( *pline == 'z' )
	    pline = skipSeparator(pline + 9);
	  else
	    pline = skipSeparator(pline + 4);

	  if ( strncmp(pline, "pressure", 6) == 0 )
	    zaxis.type = ZAXIS_PRESSURE;
	  else if ( strncmp(pline, "hybrid", 6)  == 0 )
	    zaxis.type = ZAXIS_HYBRID;
	  else if ( strncmp(pline, "height", 6)  == 0 )
	    zaxis.type = ZAXIS_HEIGHT;
	  else if ( strncmp(pline, "depth below sea", 15)  == 0 )
	    zaxis.type = ZAXIS_DEPTH_BELOW_SEA;
	  else if ( strncmp(pline, "depth below land", 16)  == 0 )
	    zaxis.type = ZAXIS_DEPTH_BELOW_LAND;
	  else if ( strncmp(pline, "isentropic", 10)  == 0 )
	    zaxis.type = ZAXIS_ISENTROPIC;
	  else if ( strncmp(pline, "surface", 7)  == 0 )
	    zaxis.type = ZAXIS_SURFACE;
	  else
	    Warning(func, "Invalid zaxisname : %s", pline);
	}
      else if ( strncmp(pline, "size", 4)  == 0 )
	{
	  zaxis.size = atol(skipSeparator(pline + 4));
	}
      else if ( strncmp(pline, "vctsize", 7)  == 0 )
	{
	  zaxis.vctsize = atol(skipSeparator(pline + 7));
	}
      else if ( strncmp(pline, "name", 4)  == 0 )
	{
	  strcpy(zaxis.name, skipSeparator(pline + 4));
	}
      else if ( strncmp(pline, "longname", 8)  == 0 )
	{
	  strcpy(zaxis.longname, skipSeparator(pline + 8));
	}
      else if ( strncmp(pline, "units", 5)  == 0 )
	{
	  strcpy(zaxis.units, skipSeparator(pline + 5));
	}
      else if ( strncmp(pline, "levels", 6)  == 0 )
	{
	  int i;
	  float flev;

	  if ( zaxis.size > 0 )
	    {
	      pline = skipSeparator(pline + 6);
	  
	      zaxis.vals = (double *) malloc(zaxis.size*sizeof(double));
	      for ( i = 0; i < zaxis.size; i++ )
		{
		  pline = skipSeparator(pline);
		  if ( strlen(pline) == 0 )
		    {
		      if ( ! readline(gfp, line, 1024) )
			{
			  Warning(func, "Incomplete commmand: >levels<");
			  break;
			}
		      pline = line;
		      pline = skipSeparator(pline);
		    }
		  sscanf(pline, "%g", &flev);
		  zaxis.vals[i] = flev;
		  while ( isalnum((int) *pline) ||
			  isdigit((int) *pline) ||
			  ispunct((int) *pline) ) pline++;
		}
	    }
	  else
	    {
	      Warning(func, "size undefined!");
	    }
	}
      else if ( strncmp(pline, "vct", 3)  == 0 )
	{
	  int i;
	  float flev;

	  if ( zaxis.vctsize > 0 )
	    {
	      pline = skipSeparator(pline + 3);
	  
	      zaxis.vct = (double *) malloc(zaxis.vctsize*sizeof(double));
	      for ( i = 0; i < zaxis.vctsize; i++ )
		{
		  pline = skipSeparator(pline);
		  if ( strlen(pline) == 0 )
		    {
		      if ( ! readline(gfp, line, 1024) )
			{
			  Warning(func, "Incomplete commmand: >vct<");
			  break;
			}
		      pline = line;
		      pline = skipSeparator(pline);
		    }
		  sscanf(pline, "%g", &flev);
		  zaxis.vct[i] = flev;
		  while ( isalnum((int) *pline) ||
			  isdigit((int) *pline) ||
			  ispunct((int) *pline) ) pline++;
		}
	    }
	  else
	    {
	      Warning(func, "vctsize undefined!");
	    }
	}
      else if ( strncmp(pline, "lbounds", 7)  == 0 )
	{
	  int i;
	  float flev;

	  if ( zaxis.size > 0 )
	    {
	      pline = skipSeparator(pline + 7);
	  
	      zaxis.lbounds = (double *) malloc(zaxis.size*sizeof(double));
	      for ( i = 0; i < zaxis.size; i++ )
		{
		  pline = skipSeparator(pline);
		  if ( strlen(pline) == 0 )
		    {
		      if ( ! readline(gfp, line, 1024) )
			{
			  Warning(func, "Incomplete commmand: >lbounds<");
			  break;
			}
		      pline = line;
		      pline = skipSeparator(pline);
		    }
		  sscanf(pline, "%g", &flev);
		  zaxis.lbounds[i] = flev;
		  while ( isalnum((int) *pline) ||
			  isdigit((int) *pline) ||
			  ispunct((int) *pline) ) pline++;
		}
	    }
	  else
	    {
	      Warning(func, "size undefined!");
	    }
	}
      else if ( strncmp(pline, "ubounds", 7)  == 0 )
	{
	  int i;
	  float flev;

	  if ( zaxis.size > 0 )
	    {
	      pline = skipSeparator(pline + 7);
	  
	      zaxis.ubounds = (double *) malloc(zaxis.size*sizeof(double));
	      for ( i = 0; i < zaxis.size; i++ )
		{
		  pline = skipSeparator(pline);
		  if ( strlen(pline) == 0 )
		    {
		      if ( ! readline(gfp, line, 1024) )
			{
			  Warning(func, "Incomplete commmand: >ubounds<");
			  break;
			}
		      pline = line;
		      pline = skipSeparator(pline);
		    }
		  sscanf(pline, "%g", &flev);
		  zaxis.ubounds[i] = flev;
		  while ( isalnum((int) *pline) ||
			  isdigit((int) *pline) ||
			  ispunct((int) *pline) ) pline++;
		}
	    }
	  else
	    {
	      Warning(func, "size undefined!");
	    }
	}
      else
	Warning(func, "Invalid zaxis command : >%s<", pline);
    }

  zaxisID = zaxisDefine(zaxis);

  return (zaxisID);
}


int cdoDefineZaxis(const char *zaxisfile)
{
  static char func[] = "cdoDefineZaxis";
  FILE *zfp;
  int zaxisID = -1;

  if ( cdoDebug ) cdoPrint("zaxis from ASCII file");
  zfp = fopen(zaxisfile, "r");
  if ( zfp == 0 )
    {
      SysError(func, zaxisfile);
    }
  else
    {
      zaxisID = zaxisFromFile(zfp);
      fclose(zfp);
    }

  if ( zaxisID == -1 ) cdoAbort("Invalid zaxis description file %s!", zaxisfile);

  return (zaxisID);
}


void defineZaxis(const char *zaxisarg)
{
  char zaxisfile[1024];
  int nfile = 0;

  while ( getoptname(zaxisfile, zaxisarg, nfile++) == 0 )
    {      
      (void) cdoDefineZaxis(zaxisfile);
    }
}

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

/*
   This module contains the following operators:

      Exprf      expr            Evaluate expressions
      Exprf      exprf           Evaluate expressions from script file
      Exprf      aexpr           Append evaluated expressions
      Exprf      aexprf          Append evaluated expressions from script file
*/
/*
Operatoren: +, -, *, \, ^, ==, !=, >, <, >=, <=, <=>, &&, ||, ?:
Functions: sqrt, exp, log, log10, sin, cos, tan, asin, acos, atan
Functions: min, max, avg, std, var
Constansts: M_PI, M_E
*/

#include <sys/types.h> /* stat */
#include <sys/stat.h>  /* stat */
#include <unistd.h>    /* stat */

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"
#include "expr.h"

void grid_cell_area(int gridID, double *array);
int getSurfaceID(int vlistID);

static
char *exprs_from_arg(const char *arg)
{
  char *exprs = NULL;
  
  size_t slen = strlen(arg);
  exprs = (char*) Malloc(slen+2);
  strcpy(exprs, operatorArgv()[0]);
  if ( exprs[slen-1] != ';' )
    {
      exprs[slen]   = ';';
      exprs[slen+1] = 0;
    }

  return exprs;
}

static
char *exprs_from_file(const char *exprf)
{
  char *exprs = NULL;
  
  /* Open expr script file for reading */
  FILE *fp = fopen(exprf, "r");
  if( fp == NULL ) cdoAbort("Open failed on %s", exprf);

  struct stat filestat;
  if ( stat(exprf, &filestat) != 0 ) cdoAbort("Stat failed on %s", exprf);

  size_t fsize = (size_t) filestat.st_size;
  exprs = (char*) Malloc(fsize+1);

  int ichar, ipos = 0;
  while ( (ichar = fgetc(fp)) != EOF ) exprs[ipos++] = ichar;

  exprs[ipos] = 0;
  if ( ipos == 0 ) cdoAbort("%s is empty!", exprf);

  fclose(fp);

  return exprs;
}

#define MAX_PARAMS 4096


static
paramType *params_new(int vlistID1)
{
  paramType *params = (paramType*) Malloc(MAX_PARAMS*sizeof(paramType));
  memset(params, 0, MAX_PARAMS*sizeof(paramType));

  int nvars1 = vlistNvars(vlistID1);

  char name[CDI_MAX_NAME];
  char longname[CDI_MAX_NAME];
  char units[CDI_MAX_NAME];
  for ( int varID = 0; varID < nvars1; varID++ )
    {
      int gridID     = vlistInqVarGrid(vlistID1, varID);
      int zaxisID    = vlistInqVarZaxis(vlistID1, varID);
      int steptype   = vlistInqVarTsteptype(vlistID1, varID);
      int ngp        = gridInqSize(gridID);
      int nlev       = zaxisInqSize(zaxisID);
      double missval = vlistInqVarMissval(vlistID1, varID);

      vlistInqVarName(vlistID1, varID, name);
      vlistInqVarLongname(vlistID1, varID, longname);
      vlistInqVarUnits(vlistID1, varID, units);
      
      params[varID].select   = false;
      params[varID].remove   = false;
      params[varID].coord    = 0;
      params[varID].gridID   = gridID;
      params[varID].zaxisID  = zaxisID;
      params[varID].steptype = steptype;
      params[varID].ngp      = ngp;
      params[varID].nlev     = nlev;
      params[varID].missval  = missval;
      params[varID].nmiss    = 0;
      params[varID].data     = NULL;
      params[varID].name     = strdup(name);
      params[varID].longname = strdup(longname);
      params[varID].units    = strdup(units);
    }

  return params;
}

static
int params_add_ts(parse_param_t *parse_arg)
{
  int varID = -1;
  paramType *params = parse_arg->params;
  if ( params )
    {
      varID = parse_arg->nparams;
      if ( varID >= parse_arg->maxparams )
        cdoAbort("Too many parameter (limit=%d)", parse_arg->maxparams);

      params[varID].name     = strdup("_ts");
      params[varID].gridID   = parse_arg->pointID;
      params[varID].zaxisID  = parse_arg->surfaceID;
      params[varID].steptype = TIME_VARIABLE;
      params[varID].ngp      = 1;
      params[varID].nlev     = 1;
      
      parse_arg->nparams++;
    }

  return varID;
}

static
void params_delete(paramType *params)
{
  if ( params )
    {
      for ( int varID = 0; varID < MAX_PARAMS; varID++ )
        {
          if ( params[varID].data )     Free(params[varID].data);
          if ( params[varID].name )     Free(params[varID].name);
          if ( params[varID].longname ) Free(params[varID].longname);
          if ( params[varID].units )    Free(params[varID].units);
        }
      Free(params);
    }
}

void *Expr(void *argument)
{
  cdoInitialize(argument);

  parse_param_t parse_arg;
  void *scanner;
  int yy_scan_string(const char *str, void *scanner);

  yylex_init(&scanner);
  yyset_extra(&parse_arg, scanner);

#define REPLACES_VARIABLES(id) cdoOperatorF1(id)
#define READS_COMMAND_LINE(id) cdoOperatorF2(id)

  cdoOperatorAdd("expr",   1, 1, "expressions");
  cdoOperatorAdd("exprf",  1, 0, "expr script filename");
  cdoOperatorAdd("aexpr",  0, 1, "expressions");
  cdoOperatorAdd("aexprf", 0, 0, "expr script filename");

  int operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  char *exprs = NULL;
  if ( READS_COMMAND_LINE(operatorID) )
    exprs = exprs_from_arg(operatorArgv()[0]);
  else
    exprs = exprs_from_file(operatorArgv()[0]);

  if ( cdoVerbose ) cdoPrint(exprs);

  int streamID1 = streamOpenRead(cdoStreamName(0));
  int vlistID1 = streamInqVlist(streamID1);
  int nvars1 = vlistNvars(vlistID1);

  int pointID   = gridCreate(GRID_GENERIC, 1);
  int surfaceID = getSurfaceID(vlistID1);

  paramType *params = params_new(vlistID1);

  parse_arg.maxparams  = MAX_PARAMS;
  parse_arg.nparams    = nvars1;
  parse_arg.nvars1     = nvars1;
  parse_arg.init       = true;
  parse_arg.debug      = false;
  if ( cdoVerbose ) parse_arg.debug = true;
  parse_arg.params     = params;
  parse_arg.pointID    = pointID;
  parse_arg.surfaceID  = surfaceID;
  parse_arg.needed     = (bool*) Malloc(nvars1*sizeof(bool));

  /* Set all input variables to 'needed' if replacing is switched off */
  for ( int varID = 0; varID < nvars1; varID++ )
    parse_arg.needed[varID] = ! REPLACES_VARIABLES(operatorID);

  int vartsID = params_add_ts(&parse_arg);
                  
  yy_scan_string(exprs, scanner);
  yyparse(&parse_arg, scanner);

  parse_arg.init = false;

  if ( cdoVerbose )
    for ( int varID = 0; varID < nvars1; varID++ )
      if ( parse_arg.needed[varID] )
	cdoPrint("Needed var: %d %s", varID, params[varID].name);

  if ( cdoVerbose )
    for ( int varID = 0; varID < parse_arg.nparams; varID++ )
      cdoPrint("var: %d %s ngp=%zu nlev=%zu coord=%c",
               varID, params[varID].name, params[varID].ngp, params[varID].nlev, params[varID].coord);

  int *varIDmap = (int*) Malloc(parse_arg.nparams*sizeof(int));

  int vlistID2 = vlistCreate();
  if ( ! REPLACES_VARIABLES(operatorID) )
    {
      vlistClearFlag(vlistID1);
      int pidx = 0;
      for ( int varID = 0; varID < nvars1; varID++ )
        {
          params[varID].select = false;
          if ( params[varID].remove == false )
            {
              varIDmap[pidx++] = varID;
              int nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
              for ( int levID = 0; levID < nlevs; levID++ )
                vlistDefFlag(vlistID1, varID, levID, TRUE);
            }
        }
      vlistCopyFlag(vlistID2, vlistID1);
    }

  for ( int pidx = 0; pidx < parse_arg.nparams; pidx++ )
    {
      if ( pidx <  nvars1 && params[pidx].select == false ) continue;
      if ( pidx >= nvars1 && params[pidx].name[0] == '_' ) continue;
      if ( pidx >= nvars1 && params[pidx].remove == true ) continue;
      if ( pidx >= nvars1 && params[pidx].coord ) continue;

      int varID = vlistDefVar(vlistID2, params[pidx].gridID, params[pidx].zaxisID, params[pidx].steptype);
      vlistDefVarName(vlistID2, varID, params[pidx].name);
      vlistDefVarMissval(vlistID2, varID, params[pidx].missval);
      if ( params[pidx].units ) vlistDefVarUnits(vlistID2, varID, params[pidx].units);
      if ( params[pidx].longname ) vlistDefVarLongname(vlistID2, varID, params[pidx].longname);
      if ( memcmp(params[pidx].name, "var", 3) == 0 )
        {
          if ( strlen(params[pidx].name) > 3 && isdigit(params[pidx].name[3]) )
            {
              int code = atoi(params[pidx].name+3);
              vlistDefVarCode(vlistID2, varID, code);
            }
        }
      varIDmap[varID] = pidx;
    }

  int nvars2 = vlistNvars(vlistID2);
  if ( nvars2 == 0 ) cdoAbort("No output variable found!");

  for ( int varID = 0; varID < nvars1; varID++ )
    {
      if ( parse_arg.needed[varID] )
        {
          size_t ngp  = params[varID].ngp;
          size_t nlev = params[varID].nlev;
          params[varID].data = (double*) Malloc(ngp*nlev*sizeof(double));
        }
    }

  for ( int varID = parse_arg.nvars1; varID < parse_arg.nparams; varID++ )
    {
      size_t ngp  = params[varID].ngp;
      size_t nlev = params[varID].nlev;
      params[varID].data = (double*) Malloc(ngp*nlev*sizeof(double));
    }

  // cleanup needed!!!
  for ( int varID = parse_arg.nvars1; varID < parse_arg.nparams; varID++ )
    {
      int coord = params[varID].coord;
      if ( coord )
        {
          if ( coord == 'x' || coord == 'y' )
            {
              int gridID = params[varID].gridID;
              if ( gridInqType(gridID) == GRID_GENERIC )
                cdoAbort("%s: not a geographical coordinate!", params[varID].name);
              if ( gridInqType(gridID) == GRID_GME )
                gridID = gridToUnstructured(gridID, 0);
              if ( gridInqType(gridID) != GRID_UNSTRUCTURED && gridInqType(gridID) != GRID_CURVILINEAR )
                gridID = gridToCurvilinear(gridID, 0);

              if      ( coord == 'x' ) gridInqXvals(gridID, params[varID].data);
              else if ( coord == 'y' ) gridInqYvals(gridID, params[varID].data);
              if ( gridID != params[varID].gridID ) gridDestroy(gridID);
            }
          else if ( coord == 'a' )
            {
              int gridID = params[varID].gridID;
              grid_cell_area(gridID, params[varID].data);
            }
          else if ( coord == 'z' )
            {
              int zaxisID = params[varID].zaxisID;
              zaxisInqLevels(zaxisID, params[varID].data);
            }
          else
            cdoAbort("Computation of coordinate %c not implemented!", coord);
        }
    }
 
  if ( cdoVerbose ) vlistPrint(vlistID2);
    
  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  int nrecs;
  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      params[vartsID].data[0] = tsID+1;
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( int varID = 0; varID < nvars1; varID++ ) params[varID].nmiss = 0;

      for ( int recID = 0; recID < nrecs; recID++ )
	{
          int varID, levelID;
	  streamInqRecord(streamID1, &varID, &levelID);
	  if ( parse_arg.needed[varID] )
	    {
	      size_t offset = params[varID].ngp*levelID;
	      double *vardata = params[varID].data + offset;
              int nmiss;
	      streamReadRecord(streamID1, vardata, &nmiss);
	      params[varID].nmiss += nmiss;
	    }
	}

      for ( int varID = 0; varID < nvars2; varID++ )
	{
          int pidx = varIDmap[varID];
          if ( pidx < nvars1 ) continue;
          size_t ngp  = params[pidx].ngp;
          size_t nlev = params[pidx].nlev;

          params[pidx].nmiss = 0;
	  memset(params[pidx].data, 0, ngp*nlev*sizeof(double));
	}

      yy_scan_string(exprs, scanner);
      yyparse(&parse_arg, scanner);

      for ( int varID = 0; varID < nvars2; varID++ )
	{
          int pidx = varIDmap[varID];
	  double missval = vlistInqVarMissval(vlistID2, varID);

          size_t ngp = params[pidx].ngp;
          int nlev = (int) params[pidx].nlev;
	  for ( int levelID = 0; levelID < nlev; levelID++ )
	    {
              size_t offset = ngp*levelID;
	      double *vardata = params[pidx].data + offset;

	      int nmiss = 0;
	      for ( size_t i = 0; i < ngp; i++ )
		if ( DBL_IS_EQUAL(vardata[i], missval) ) nmiss++;

	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, vardata, nmiss);
	    }
	}

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  vlistDestroy(vlistID2);

  yylex_destroy(scanner);

  if ( exprs ) Free(exprs);

  params_delete(params);

  if ( parse_arg.needed ) Free(parse_arg.needed);
  if ( varIDmap ) Free(varIDmap);

  cdoFinish();

  return 0;
}

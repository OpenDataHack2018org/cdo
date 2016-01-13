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
#include "expr.h"


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

#define MAX_PARAM 4096

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

  int vlistID2 = -1;
  int nvars = 0;
  if ( REPLACES_VARIABLES(operatorID) )
    {
      vlistID2 = vlistCreate();
      nvars = 0;
    }
  else
    {
      vlistID2 = vlistDuplicate(vlistID1);
      nvars = nvars1;
    }

  int vlisttmp = vlistCreate();

  int nparam         = MAX_PARAM;
  paramType *param1  = (paramType*) Malloc(nvars1*sizeof(paramType));
  parse_arg.init     = true;
  parse_arg.debug    = false;
  if ( cdoVerbose ) parse_arg.debug = true;
  parse_arg.nvars1   = nvars1;
  parse_arg.param1   = param1;
  parse_arg.vlistID1 = vlistID1;
  parse_arg.vlistID2 = vlistID2;
  parse_arg.vlisttmp = vlisttmp;
  parse_arg.gridID2  = -1;
  parse_arg.zaxisID2 = -1;
  parse_arg.tsteptype2 = -1;
   
  for ( int varID = 0; varID < nvars1; varID++ )
    {
      int gridID     = vlistInqVarGrid(vlistID1, varID);
      int zaxisID    = vlistInqVarZaxis(vlistID1, varID);
      int ngp        = gridInqSize(gridID);
      int nlev       = zaxisInqSize(zaxisID);
      double missval = vlistInqVarMissval(vlistID1, varID);
      
      param1[varID].gridID  = gridID;
      param1[varID].zaxisID = zaxisID;
      param1[varID].ngp     = ngp;
      param1[varID].nlev    = nlev;
      param1[varID].missval = missval;
      param1[varID].nmiss   = 0;
      param1[varID].data    = NULL;
    }

  /* Set all input variables to 'needed' if replacing is switched off */
  for ( int varID = 0; varID < nvars1; varID++ )
    parse_arg.needed[varID] = ! REPLACES_VARIABLES(operatorID);

  yy_scan_string(exprs, scanner);
  yyparse(&parse_arg, scanner);

  parse_arg.init = false;

  int nvars2 = vlistNvars(vlistID2);
  if ( nvars2 == 0 ) cdoAbort("No output variable found!");

  if ( cdoVerbose ) vlistPrint(vlistID2);

  if ( cdoVerbose )
    for ( int varID = 0; varID < nvars1; varID++ )
      if ( parse_arg.needed[varID] )
	printf("Needed var: %d %s\n", varID, parse_arg.varname[varID]);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  parse_arg.vardata2 = (double**) Malloc(nvars2*sizeof(double*));

  for ( int varID = 0; varID < nvars1; varID++ )
    {
      if ( parse_arg.needed[varID] )
        {
          int ngp  = param1[varID].ngp;
          int nlev = param1[varID].nlev;
          param1[varID].data = (double*) Malloc(ngp*nlev*sizeof(double));
        }
    }

  for ( int varID = 0; varID < nvars; varID++ )
    {
      parse_arg.vardata2[varID] = param1[varID].data;
    }

  for ( int varID = nvars; varID < nvars2; varID++ )
    {
      int gridID  = vlistInqVarGrid(vlistID2, varID);
      int zaxisID = vlistInqVarZaxis(vlistID2, varID);

      int gridsize = gridInqSize(gridID);
      int nlevel   = zaxisInqSize(zaxisID);
      parse_arg.vardata2[varID] = (double*) Malloc(gridsize*nlevel*sizeof(double));
    }

  int nrecs;
  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( int varID = 0; varID < nvars1; varID++ ) param1[varID].nmiss = 0;

      for ( int recID = 0; recID < nrecs; recID++ )
	{
          int varID, levelID;
	  streamInqRecord(streamID1, &varID, &levelID);
	  if ( parse_arg.needed[varID] )
	    {
	      int offset = param1[varID].ngp*levelID;
              int nmiss;
	      double *vardata = param1[varID].data + offset;
	      streamReadRecord(streamID1, vardata, &nmiss);
	      param1[varID].nmiss += nmiss;
	    }
	}

      for ( int varID = nvars; varID < nvars2; varID++ )
	{
	  int gridID   = vlistInqVarGrid(vlistID2, varID);
	  int zaxisID  = vlistInqVarZaxis(vlistID2, varID);
	  int gridsize = gridInqSize(gridID);
	  int nlevel   = zaxisInqSize(zaxisID);

	  memset(parse_arg.vardata2[varID], 0, gridsize*nlevel*sizeof(double));
	}

      yy_scan_string(exprs, scanner);
      yyparse(&parse_arg, scanner);

      for ( int varID = 0; varID < nvars2; varID++ )
	{
	  int gridID  = vlistInqVarGrid(vlistID2, varID);
	  int zaxisID = vlistInqVarZaxis(vlistID2, varID);
	  double missval = vlistInqVarMissval(vlistID2, varID);

	  int gridsize = gridInqSize(gridID);
	  int nlevel   = zaxisInqSize(zaxisID);
	  for ( int levelID = 0; levelID < nlevel; levelID++ )
	    {
	      int offset = gridsize*levelID;
	      double *vardata = parse_arg.vardata2[varID] + offset;

	      int nmiss = 0;
	      for ( int i = 0; i < gridsize; i++ )
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

  if ( param1 )
    {
      for ( int varID = 0; varID < nvars1; varID++ )
        if ( param1[varID].data ) Free(param1[varID].data);
      
      Free(param1);
    }

  cdoFinish();

  return 0;
}

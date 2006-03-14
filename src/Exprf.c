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

#include <string.h>    /* memcpy */
#include <sys/types.h> /* stat */
#include <sys/stat.h>  /* stat */
#include <unistd.h>    /* stat */

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "expr.h"

/*
Operatoren: +, -, *, \
Functions: sqrt, exp, log, log10, sin, cos, tan, asin, acos, atan
Functions: min, max, avg, std, var
Constansts: M_PI, M_E
*/

/*
@BeginDoc

@BeginModule

@Name      = Expr
@Title     = Arithmetic processor
@Section   = Arithmetic processor
@Class     = Arithmetic
@Arguments = ifile ofile
@Operators = expr exprf

@EndModule


@BeginOperator_expr

@Title     = Evaluate expressions
@Parameter = instr

@BeginDescription
This Operator arithmetically processes every timestep of ifile.
The processing instructions are read from the parameter.
Each individual assignment statement must end with a semi-colon.
@EndDescription

@BeginParameter
@Item = instr
STRING  Processing instructions
@EndParameter

@BeginExample
cdo -expr,'var1=aprl+aprc;var2=ts-273.15;' ifile ofile
@EndExample

@EndOperator


@BeginOperator_exprf

@Title     = Evaluate expressions from script file
@Parameter = filename

@BeginDescription
This Operator arithmetically processes every timestep of ifile.
Contrary to expr the processing instructions are read from a file.
@EndDescription

@BeginParameter
@Item = filename
STRING  File with processing instructions
@EndParameter

@BeginExample
cdo -expr,myexpr ifile ofile

and the ASCII file myexpr contains e.g:
#############################
var1 = aprl + aprc;
var2 = ts - 273.15;
#############################
@EndExample

@EndOperator

@EndDoc
*/

void *Expr(void *argument)
{
  static char func[] = "Expr";
  int EXPR, EXPRF;
  int operatorID;
  char *exprs = NULL;
  const char *exprf = NULL;
  int streamID1, streamID2 = CDI_UNDEFID;
  int offset;
  int nrecs, nvars, nvars2;
  int gridID, zaxisID;
  int tsID, recID, varID, levelID;
  int vlistID1, vlistID2;
  int gridsize, nlevel;
  int nmiss;
  int taxisID1, taxisID2;
  int lwarn = TRUE;
  double missval;
  double *array = NULL;
  double *single1, *single2;
  prs_sct prs_arg;
  extern int yyparse(void *);
  extern int yy_scan_string(const char *);

  cdoInitialize(argument);

  EXPR  = cdoOperatorAdd("expr",  0, 0, "expressions");
  EXPRF = cdoOperatorAdd("exprf", 0, 0, "expr script filename");

  operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  if ( operatorID == EXPR )
    {
      exprs = operatorArgv()[0];
    }
  else
    {
      int ichar, ipos = 0;
      FILE *fp;
      size_t fsize;
      struct stat filestat;

      exprf = operatorArgv()[0];

      /* Open script file for reading */
      if( (fp = fopen(exprf, "r")) == NULL ) cdoAbort("Open failed on %s", exprf);

      if ( stat(exprf, &filestat) != 0 ) cdoAbort("Stat failed on %s", exprf);

      fsize = (size_t) filestat.st_size;
      exprs = (char *) malloc(fsize+1);

      while ( (ichar = fgetc(fp)) != EOF ) exprs[ipos++] = ichar;

      exprs[ipos] = 0;

      if ( ipos == 0 ) cdoAbort("%s is empty!", exprf);

      fclose(fp);
    }

  if ( cdoVerbose ) cdoPrint(exprs);


  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);

  nvars = vlistNvars(vlistID1);

  vlistID2 = vlistNew();

  prs_arg.init = 1;
  prs_arg.vlistID1 = vlistID1;
  prs_arg.vlistID2 = vlistID2;
  prs_arg.nvars1 = 0;
  prs_arg.debug  = 0;
  for ( varID = 0; varID < nvars; varID++ )
    prs_arg.var_needed[varID] = FALSE;

  yy_scan_string(exprs);
  yyparse((void *) &prs_arg);

  prs_arg.init = 0;

  nvars2 = vlistNvars(vlistID2);

  if ( nvars2 == 0 ) cdoAbort("No variable in output!");

  if ( cdoVerbose ) vlistPrint(vlistID2);

  for ( varID = 0; varID < nvars; varID++ )
    if ( prs_arg.var_needed[varID] && cdoVerbose )
      printf("var_needed: %d %s\n", varID, prs_arg.var[varID]);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  prs_arg.vardata1 = (double **) malloc(nvars*sizeof(double*));
  prs_arg.vardata2 = (double **) malloc(nvars2*sizeof(double*));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID  = vlistInqVarGrid(vlistID1, varID);
      zaxisID = vlistInqVarZaxis(vlistID1, varID);
      /*     prs_arg.missval = vlistInqVarMissval(vlistID1, varID); */

      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(zaxisID);
      if ( prs_arg.var_needed[varID] )
	prs_arg.vardata1[varID] = (double *) malloc(gridsize*nlevel*sizeof(double));
      else
	prs_arg.vardata1[varID] = NULL;
    }

  for ( varID = 0; varID < nvars2; varID++ )
    {
      gridID  = vlistInqVarGrid(vlistID2, varID);
      zaxisID = vlistInqVarZaxis(vlistID2, varID);

      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(zaxisID);
      prs_arg.vardata2[varID] = (double *) malloc(gridsize*nlevel*sizeof(double));
    }

  gridsize = vlistGridsizeMax(vlistID1);
  array = (double *) malloc(gridsize*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  if ( prs_arg.var_needed[varID] )
	    {
	      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      offset   = gridsize*levelID;
	      single1  = prs_arg.vardata1[varID] + offset;
	      streamReadRecord(streamID1, single1, &nmiss);
	      if ( nmiss && lwarn )
		{
		  cdoWarning("Missing values unsupported for this operator!");
		  lwarn = FALSE;
		}
	    }
	}

      for ( varID = 0; varID < nvars2; varID++ )
	{
	  gridID   = vlistInqVarGrid(vlistID2, varID);
	  zaxisID  = vlistInqVarZaxis(vlistID2, varID);
	  gridsize = gridInqSize(gridID);
	  nlevel   = zaxisInqSize(zaxisID);

	  memset(prs_arg.vardata2[varID], 0, gridsize*nlevel*sizeof(double));
	}

      yy_scan_string(exprs);
      yyparse((void *) &prs_arg);

      for ( varID = 0; varID < nvars2; varID++ )
	{
	  gridID   = vlistInqVarGrid(vlistID2, varID);
	  zaxisID  = vlistInqVarZaxis(vlistID2, varID);
	  missval  = vlistInqVarMissval(vlistID2, varID);

	  gridsize = gridInqSize(gridID);
	  nlevel   = zaxisInqSize(zaxisID);
	  /* nmiss    = prs_arg.nmiss; */
	  /* if ( nmiss ) fprintf(stdout, "out nmiss = %d\n", nmiss); */
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      offset   = gridsize*levelID;
	      single2  = prs_arg.vardata2[varID] + offset;
	      nmiss = 0;
	      if ( missval < -1.e30 || missval > 1.e30 )
		{
		  int i;
		  for ( i = 0; i < gridsize; i++ )
		    if ( single2[i] < -1.e30 || single2[i] > 1.e30 )
		      {
			single2[i] = missval;
			nmiss++;
		      }
		}
	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, single2, nmiss);
	    }
	}

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( array ) free(array);

  cdoFinish();

  return (0);
}

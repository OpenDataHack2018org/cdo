/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2009 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Select2     select         Select fields
*/


#include <string.h>
#include <math.h>   /* fabs */

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "error.h"
#include "util.h"
#include "functs.h"
#include "list.h"
#include "namelist.h"


#define  PML_INT         1
#define  PML_FLT         2
#define  PML_WORD        3
#define  PML_DATE        4
#define  PML_TIME        4

typedef struct {
  char *name;
  int maxpar;
  int numpar;
  char *par;
}
PARAMETER;

#define PML_DEF_INT(name, size, val)  int par##name[size]; int flag##name[size]; int npar##name = 0; int name = 0
#define PML_DEF_FLT(name, size, val)  double par##name[size]; int npar##name = 0; double name = 0
#define PML_ADD_INT(nml, name) pmlAdd(nml, #name, PML_INT, 0, par##name, sizeof(par##name)/sizeof(int))
#define PML_ADD_FLT(nml, name) pmlAdd(nml, #name, PML_FLT, 0, par##name, sizeof(par##name)/sizeof(double))
#define PML_NUM(nml, name)   npar##name = pmlNum(nml, #name)
#define PML_PAR(name)        npar##name, par##name, name

static PARAMETER Parameter[] =
{
  {"code", 1024, 0, NULL},
};

static int NumParameter = sizeof(Parameter) / sizeof(Parameter[0]);

#define MAX_PML_ENTRY  256

typedef struct
{
  char *name;
  size_t len;
  void *ptr;
  int type;
  int occ;
  int dis;
  size_t size;
} pml_entry_t;


typedef struct
{
  int size;
  int dis;
  char *name;
  /* PML_LINE line; */
  pml_entry_t *entry[MAX_PML_ENTRY];
} pml_t;


static void pml_init(pml_t *pml, const char *name)
{
  static char func[] = "pml_init";
  pml->size = 0;
  pml->dis  = 1;
  pml->name = strdup(name);
}


pml_t *pmlNew(const char *name)
{
  static char func[] = "pmlNew";
  pml_t *pml;

  pml = (pml_t *) malloc(sizeof(pml_t));

  pml_init(pml, name);

  return (pml);
}


void pmlPrint(pml_t *pml)
{
  pml_entry_t *entry;
  int i, j, nout;

  if ( pml == NULL ) return;

  fprintf(stdout, "Parameter list: %s\n", pml->name);
  fprintf(stdout, " Num  Name             Type  Size   Dis   Occ  Entries\n");

  for ( i = 0; i < pml->size; ++i )
    {
      entry = pml->entry[i];
      fprintf(stdout, "%4d  %-16s %4d  %4d  %4d  %4d ",
	      i+1, pml->entry[i]->name, pml->entry[i]->type, (int)pml->entry[i]->size,
	      pml->entry[i]->dis, pml->entry[i]->occ);
      nout = pml->entry[i]->occ;
      if ( nout > 8 ) nout = 8;

      if      ( entry->type == PML_WORD )
	for ( j = 0; j < nout; j++ )
	  fprintf(stdout, " %s", ((char **)entry->ptr)[j]);
      else if ( entry->type == PML_INT )
	for ( j = 0; j < nout; j++ )
	  fprintf(stdout, " %d", ((int *)entry->ptr)[j]);
      else if ( entry->type == PML_FLT )
	for ( j = 0; j < nout; j++ )
	  fprintf(stdout, " %g", ((double *)entry->ptr)[j]);
      
      fprintf(stdout, "\n");
    }
}


int pmlAdd(pml_t *pml, const char *name, int type, int dis, void *ptr, size_t size)
{
  static char func[] = "pmlAdd";
  pml_entry_t *pml_entry;
  int entry = 0;

  if ( pml->size >= MAX_PML_ENTRY )
    {
      fprintf(stderr, "Too many entries in parameter list %s! (Max = %d)\n", pml->name, MAX_PML_ENTRY);
      return (-1);
    }

  pml_entry = (pml_entry_t *) malloc(sizeof(pml_entry_t));

  pml_entry->name = strdup(name);
  pml_entry->len  = strlen(name);
  pml_entry->type = type;
  pml_entry->ptr  = ptr;
  pml_entry->size = size;
  pml_entry->dis  = dis;
  pml_entry->occ  = 0;

  entry = pml->size;
  pml->entry[pml->size++] = pml_entry;

  return (entry);
}


int pmlNum(pml_t *pml, const char *name)
{
  pml_entry_t *entry;
  int i, nocc = 0;

  if ( pml == NULL ) return (nocc);

  for ( i = 0; i < pml->size; i++ )
    {
      entry = pml->entry[i];
      if ( strcmp(name, entry->name) == 0 )
	{
	  nocc = entry->occ;
	  break;
	}
    }

  if ( i == pml->size )
    fprintf(stderr, "Parameter list entry %s not found in %s\n", name, pml->name);

  return (nocc);
}


int pml_add_entry(pml_entry_t *entry, char *arg)
{
  int status = 0;

  if ( entry->type == PML_INT )
    {
      if ( entry->occ < (int) entry->size )
	((int *) entry->ptr)[entry->occ++] = atoi(arg);
    }
  else if ( entry->type == PML_FLT )
    {
      if ( entry->occ < (int) entry->size )
	((double *) entry->ptr)[entry->occ++] = atof(arg);
    }
  else
    {
      fprintf(stderr, "unsupported type!\n");
      return (status);
    }
}


int pmlProcess(pml_entry_t *entry, int argc, char **argv)
{
  int i;
  int len;
  char *parg;
  char *epos;

  for ( i = 0; i < argc; ++i )
    {
      parg = argv[i];
      if ( i == 0 )
	{
	  epos = strchr(parg, '=');
	  if ( epos == NULL )
	    {
	      fprintf(stderr, "internal problem, keyword not found!\n");
	    }
	  parg += epos-parg+1;
	}

      printf("process: %s %d %s\n", entry->name, i, parg);
      pml_add_entry(entry, parg);
    }
}


int pmlRead(pml_t *pml, int argc, char **argv)
{
  static char func[] = "pmlRead";
  pml_entry_t *entry;
  pml_entry_t *pentry[MAX_PML_ENTRY];
  int params[MAX_PML_ENTRY];
  int num_par[MAX_PML_ENTRY];
  int len_par[MAX_PML_ENTRY];
  int nparams = 0;
  int i, istart;
  char *epos;
  size_t len;
  char *parbuf;
  int bufsize = 0;
  int status = 0;

  for ( i = 0; i < argc; ++i ) printf("pmlRead: %d %s\n", i, argv[i]);

  for ( i = 0; i < argc; ++i )
    {
      len = strlen(argv[i]);
      len_par[i] = (int)len;
      bufsize += len+1;
    }

  printf("bufsize %d\n", bufsize);
  parbuf = (char *) malloc(bufsize*sizeof(char));
  memset(parbuf, 0, bufsize*sizeof(char));

  istart = 0;
  while ( istart < argc )
    {

      epos = strchr(argv[istart], '=');
      if ( epos == NULL )
	{
	  fprintf(stderr, "Parameter >%s< has no keyword!\n", argv[istart]);
	  status = 1;
	  goto END_LABEL;
	}

      len = epos - argv[istart];
      printf ("len = %d\n", len);
      for ( i = 0; i < pml->size; ++i )
	{
	  entry = pml->entry[i];
	  if ( entry->len == len )
	    if ( strncmp(entry->name, argv[istart], len) == 0 ) break;
	}

      if ( i == pml->size )
	{
	  fprintf(stderr, "Parameter >%s< not available!\n", argv[istart]);
	  status = 2;
	  goto END_LABEL;
	}

      num_par[nparams] = 0;
      pentry[nparams] = entry;
      params[nparams] = istart;
      num_par[nparams] = 1;
      
      istart++;
      for ( i = istart; i < argc; ++i )
	{
	  epos = strchr(argv[i], '=');
	  if ( epos != NULL ) break;

	  num_par[nparams]++;
	}

      istart = i;

      nparams++;
    }

  for ( i = 0; i < nparams; ++i )
    {
      printf("param %d %s %d %s\n", i, pentry[i]->name, num_par[i],  argv[params[i]]);
      pmlProcess(pentry[i], num_par[i], &argv[params[i]]);
    }


 END_LABEL:

  free(parbuf);

  return (status);
}


int par_check_int(int npar, int *parlist, int *flaglist, int par)
{
  int i, found;

  if ( npar == 0 ) found = 1;
  else             found = 0;

  for ( i = 0; i < npar; i++ )
    if ( par == parlist[i] ) { found = 1; break; }

  return (found);
}


void *Select2(void *argument)
{
  const char func[] = "Select2";
  int SELECT, SELCODE, SELNAME, SELLEVEL, SELLEVIDX, SELGRID, SELGRIDNAME, SELZAXIS, SELZAXISNAME, SELLTYPE; 
  int SELTABNUM, DELCODE, DELNAME, SELSTDNAME;
  int operatorID;
  int streamID1, streamID2;
  int tsID, nrecs;
  int nvars, nlevs;
  int code, tabnum, gridID, zaxisID, levID;
  double level;
  int varID2, levelID2;
  int recID, varID, levelID;
  int *intarr = NULL, nsel = 0;
  int *selfound = NULL;
  double *fltarr = NULL;
  char varname[256];
  char stdname[256];
  char gridname[256];
  char zaxisname[256];
  char **argnames = NULL;
  int vlistID1 = -1, vlistID2 = -1;
  int isel;
  int i;
  int lcopy = FALSE;
  int gridsize;
  int nmiss;
  double *array = NULL;
  int taxisID1, taxisID2;
  int ltype;
  int sargc;
  char **sargv;
  LIST *ilist = listNew(INT_LIST);
  LIST *flist = listNew(FLT_LIST);
  pml_t *pml;
  PML_DEF_INT(xcode,   1024, 0);
  PML_DEF_FLT(xlevel,  1024, 0);

  cdoInitialize(argument);

  SELECT  = cdoOperatorAdd("select", 0, 0, "parameter list");

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  nsel     = operatorArgc();
  argnames = operatorArgv();

  if ( cdoVerbose )
    for ( i = 0; i < nsel; i++ )
      printf("name %d = %s\n", i+1, argnames[i]);

  sargc = nsel;
  sargv = (char **) malloc(sargc*sizeof(char *));

  for ( i = 0; i < nsel; i++ )
    {
      argnames = operatorArgv();
      /*
      if ( i == 0 )
	{
	  if ( strncmp(argnames[i], "xcode=", 6) == 0 )
	    {
	      sargv[i] = strdup(argnames[i]+5);
	    }
	  else
	    {
	      cdoAbort("Parameter >code< not found");
	    }
	}
      else
      */
	{
	  sargv[i] = strdup(argnames[i]);
	}
    }

  if ( cdoVerbose )
    for ( i = 0; i < sargc; i++ )
      printf("sargc %d = %s\n", i+1, sargv[i]);

  nsel = args2intlist(sargc, sargv, ilist);
  intarr = (int *) listArrayPtr(ilist);

  if ( cdoVerbose )
    for ( i = 0; i < nsel; i++ )
      printf("int %d = %d\n", i+1, intarr[i]);

  pml = pmlNew("SELECT");

  PML_ADD_INT(pml, xcode);
  PML_ADD_FLT(pml, xlevel);
  /*
  pmlAdd(pml, "i2",  PML_INT,    1, &i2,  sizeof(i2)/sizeof(int));
  pmlAdd(pml, "dm",  PML_FLT,    1, &dm,  sizeof(dm)/sizeof(double));
  pmlAdd(pml, "var", PML_WORD,   0, var,  sizeof(var)/sizeof(char *));
  */
  pmlRead(pml, nsel, argnames);

  pmlPrint(pml);

  printf("nparxcode: %d\n", PML_NUM(pml, xcode));
  /*
  pmlDelete(pml);
  */

  if ( operatorID == SELNAME || operatorID == DELNAME || operatorID == SELSTDNAME ||
       operatorID == SELGRIDNAME || operatorID == SELZAXISNAME )
    {
      nsel     = operatorArgc();
      argnames = operatorArgv();

      if ( cdoVerbose )
	for ( i = 0; i < nsel; i++ )
	  printf("name %d = %s\n", i+1, argnames[i]);
    }
  else if ( operatorID == SELLEVEL )
    {
      nsel = args2fltlist(operatorArgc(), operatorArgv(), flist);
      fltarr = (double *) listArrayPtr(flist);

      if ( cdoVerbose )
	for ( i = 0; i < nsel; i++ )
	  printf("flt %d = %g\n", i+1, fltarr[i]);
    }
  else
    {
      /*      nsel = args2intlist(operatorArgc(), operatorArgv(), ilist); */
      /*
      nsel = args2intlist(sargc, sargv, ilist);
      intarr = (int *) listArrayPtr(ilist);

      if ( cdoVerbose )
	for ( i = 0; i < nsel; i++ )
	  printf("int %d = %d\n", i+1, intarr[i]);
      */
    }

  operatorID = SELCODE;

  if ( nsel )
    {
      selfound = (int *) malloc(nsel*sizeof(int));
      for ( i = 0; i < nsel; i++ ) selfound[i] = FALSE;
    }

  /*
  if ( nsel == 0 )
    cdoAbort("missing code argument!");
  */
  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);

  vlistClearFlag(vlistID1);
  nvars = vlistNvars(vlistID1);
  for ( varID = 0; varID < nvars; varID++ )
    {
      vlistInqVarName(vlistID1, varID, varname);
      vlistInqVarStdname(vlistID1, varID, stdname);
      code    = vlistInqVarCode(vlistID1, varID);
      tabnum  = tableInqNum(vlistInqVarTable(vlistID1, varID));
      gridID  = vlistInqVarGrid(vlistID1, varID);
      zaxisID = vlistInqVarZaxis(vlistID1, varID);
      nlevs   = zaxisInqSize(zaxisID);
      gridName(gridInqType(gridID), gridname);
      zaxisName(zaxisInqType(zaxisID), zaxisname);

      for ( levID = 0; levID < nlevs; levID++ )
	{
	  level = zaxisInqLevel(zaxisID, levID);

	  if ( operatorID == DELCODE || operatorID == DELNAME )
	    vlistDefFlag(vlistID1, varID, levID, TRUE);

	  printf("xcode: %d %d\n", code, par_check_int(nparxcode, parxcode, flagxcode, code));

	  for ( isel = 0; isel < nsel; isel++ )
	    {
	      if ( operatorID == SELCODE )
		{
		  if ( intarr[isel] == code )
		    {
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == SELNAME )
		{
		  if ( strcmp(argnames[isel], varname) == 0 )
		    {
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == SELSTDNAME )
		{
		  if ( strcmp(argnames[isel], stdname) == 0 )
		    {
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == SELLEVEL )
		{
		  if ( fabs(fltarr[isel] - level) < 0.0001 )
		    {
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == SELLEVIDX )
		{
		  if ( intarr[isel] == (levID+1) )
		    {
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == SELGRID )
		{
		  if ( intarr[isel] == (gridID+1) )
		    {
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == SELGRIDNAME )
		{
		  if ( strncmp(argnames[isel], gridname, strlen(argnames[isel])) == 0 )
		    {
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == SELZAXIS )
		{
		  if ( intarr[isel] == (zaxisID+1) )
		    {
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == SELZAXISNAME )
		{
		  if ( strncmp(argnames[isel], zaxisname, strlen(argnames[isel])) == 0 )
		    {
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == SELTABNUM )
		{
		  if ( intarr[isel] == tabnum )
		    {
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == DELCODE )
		{
		  if ( intarr[isel] == code )
		    {
		      vlistDefFlag(vlistID1, varID, levID, FALSE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == DELNAME )
		{
		  if ( strcmp(argnames[isel], varname) == 0 )
		    {
		      vlistDefFlag(vlistID1, varID, levID, FALSE);
		      selfound[isel] = TRUE;
		    }
		}
	      else if ( operatorID == SELLTYPE )
		{
		  ltype = zaxis2ltype(zaxisID);

		  if ( intarr[isel] == ltype )
		    {
		      vlistDefFlag(vlistID1, varID, levID, TRUE);
		      selfound[isel] = TRUE;
		    }
		}
	    }
	}
    }

  for ( isel = 0; isel < nsel; isel++ )
    {
      if ( selfound[isel] == FALSE )
	{
	  if ( operatorID == SELCODE || operatorID == DELCODE )
	    {
	      cdoWarning("Code number %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELNAME || operatorID == DELNAME )
	    {
	      cdoWarning("Variable name %s not found!", argnames[isel]);
	    }
	  else if ( operatorID == SELSTDNAME )
	    {
	      cdoWarning("Variable with standard name %s not found!", argnames[isel]);
	    }
	  else if ( operatorID == SELLEVEL )
	    {
	      cdoWarning("Level %g not found!", fltarr[isel]);
	    }
	  else if ( operatorID == SELLEVIDX )
	    {
	      cdoWarning("Level index %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELGRID )
	    {
	      cdoWarning("Grid %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELGRIDNAME )
	    {
	      cdoWarning("Grid name %s not found!", argnames[isel]);
	    }
	  else if ( operatorID == SELZAXIS )
	    {
	      cdoWarning("Zaxis %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELZAXISNAME )
	    {
	      cdoWarning("Zaxis name %s not found!", argnames[isel]);
	    }
	  else if ( operatorID == SELTABNUM )
	    {
	      cdoWarning("Table number %d not found!", intarr[isel]);
	    }
	  else if ( operatorID == SELLTYPE )
	    {
	      cdoWarning("GRIB level type %d not found!", intarr[isel]);
	    }
	}
    }

  vlistID2 = vlistCreate();
  vlistCopyFlag(vlistID2, vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  nrecs = vlistNrecs(vlistID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(vlistID1);
      array = (double *) malloc(gridsize*sizeof(double));
    }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
     
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  if ( vlistInqFlag(vlistID1, varID, levelID) == TRUE )
	    {
	      varID2   = vlistFindVar(vlistID2, varID);
	      levelID2 = vlistFindLevel(vlistID2, varID, levelID);

	      streamDefRecord(streamID2, varID2, levelID2);
	      if ( lcopy )
		{
		  streamCopyRecord(streamID2, streamID1);
		}
	      else
		{
		  streamReadRecord(streamID1, array, &nmiss);
		  streamWriteRecord(streamID2, array, nmiss);
		}
     	    }
	}
       
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);
 
  vlistDestroy(vlistID2);

  if ( ! lcopy )
    if ( array ) free(array);

  if ( selfound ) free(selfound);

  listDelete(ilist);
  listDelete(flist);

  cdoFinish();

  return (NULL);
}

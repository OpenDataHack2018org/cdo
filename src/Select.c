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

      Select      select         Select fields
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "error.h"
#include "util.h"
//#include "list.h"

double datestr_to_double(const char *datestr, int opt);

#define  PML_INT         1
#define  PML_FLT         2
#define  PML_WORD        3
#define  PML_DATE        4
#define  PML_TIME        4

/*
typedef struct {
  char *name;
  int maxpar;
  int numpar;
  char *par;
}
PARAMETER;

static PARAMETER Parameter[] =
{
  {"code", 1024, 0, NULL},
};

static int NumParameter = sizeof(Parameter) / sizeof(Parameter[0]);
*/

#define PML_DEF(name, size, txt)      bool flag_##name[size]; int npar_##name = 0; int max_##name = size; const char str_##name[] = txt
#define PML_DEF_INT(name, size, txt)  int par_##name[size]; int name = 0; PML_DEF(name, size, txt)
#define PML_DEF_FLT(name, size, txt)  double par_##name[size]; double name = 0; PML_DEF(name, size, txt)
#define PML_DEF_WORD(name, size, txt) char *par_##name[size]; const char *name = 0; PML_DEF(name, size, txt)
#define PML_INIT_INT(name)            memset(flag_##name, 0, max_##name * sizeof(bool))
#define PML_INIT_FLT(name)            memset(flag_##name, 0, max_##name * sizeof(bool))
#define PML_INIT_WORD(name)           memset(flag_##name, 0, max_##name * sizeof(bool))
#define PML_ADD_INT(nml, name)        pmlAdd(nml, #name, PML_INT,  0, par_##name, sizeof(par_##name)/sizeof(int))
#define PML_ADD_FLT(nml, name)        pmlAdd(nml, #name, PML_FLT,  0, par_##name, sizeof(par_##name)/sizeof(double))
#define PML_ADD_WORD(nml, name)       pmlAdd(nml, #name, PML_WORD, 0, par_##name, sizeof(par_##name)/sizeof(char *))
#define PML_NUM(nml, name)            npar_##name = pmlNum(nml, #name)
#define PML_PAR(name)                 npar_##name, par_##name, name

#define PAR_CHECK_INT_FLAG(name)      par_check_int_flag(npar_##name, par_##name, flag_##name, str_##name)
#define PAR_CHECK_FLT_FLAG(name)      par_check_flt_flag(npar_##name, par_##name, flag_##name, str_##name)
#define PAR_CHECK_WORD_FLAG(name)     par_check_word_flag(npar_##name, par_##name, flag_##name, str_##name)
#define PAR_CHECK_INT(name)           par_check_int(npar_##name, par_##name, flag_##name, name)
#define PAR_CHECK_FLT(name)           par_check_flt(npar_##name, par_##name, flag_##name, name)
#define PAR_CHECK_WORD(name)          par_check_word(npar_##name, par_##name, flag_##name, name)
#define PAR_CHECK_DATE(name)          par_check_date(npar_##name, par_##name, flag_##name, name)
#define PAR_CHECK_SEASON(name, month) par_check_season(npar_##name, par_##name, flag_##name, month)

#define MAX_PLIST_ENTRY  256
#define MAX_PML_ENTRY    256

typedef struct
{
  char *text;
  size_t len;
  void *ptr;
  int type;
  int occ;
  int dis;
  size_t size;
} plist_entry_t;


typedef struct
{
  int size;
  plist_entry_t *entry[MAX_PLIST_ENTRY];
} plist_t;


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


static
void pml_init(pml_t *pml, const char *name)
{
  pml->size = 0;
  pml->dis  = 1;
  pml->name = strdup(name);
}


pml_t *pmlNew(const char *name)
{
  pml_t *pml = (pml_t*) Malloc(sizeof(pml_t));

  pml_init(pml, name);

  return pml;
}


void pmlDestroy(pml_t *pml)
{
  if ( pml == NULL ) return;

  for ( int i = 0; i < pml->size; ++i )
    {
      if ( pml->entry[i] ) Free(pml->entry[i]);
    }

  Free(pml);
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
  if ( pml->size >= MAX_PML_ENTRY )
    {
      fprintf(stderr, "Too many entries in parameter list %s! (Max = %d)\n", pml->name, MAX_PML_ENTRY);
      return -1;
    }

  pml_entry_t *pml_entry = (pml_entry_t*) Malloc(sizeof(pml_entry_t));

  pml_entry->name = strdup(name);
  pml_entry->len  = strlen(name);
  pml_entry->type = type;
  pml_entry->ptr  = ptr;
  pml_entry->size = size;
  pml_entry->dis  = dis;
  pml_entry->occ  = 0;

  int entry = pml->size;
  pml->entry[pml->size++] = pml_entry;

  return entry;
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

  return nocc;
}

void split_intstring(const char *intstr, int *first, int *last, int *inc);

int pml_add_entry(pml_entry_t *entry, char *arg)
{
  int status = 0;

  if ( entry->type == PML_INT )
    {
      int ival, first, last, inc;

      split_intstring(arg, &first, &last, &inc);

      if ( inc >= 0 )
	{
	  for ( ival = first; ival <= last; ival += inc )
	    if ( entry->occ < (int) entry->size )
	      ((int *) entry->ptr)[entry->occ++] = ival;
	}
      else
	{
	  for ( ival = first; ival >= last; ival += inc )
	    if ( entry->occ < (int) entry->size )
	      ((int *) entry->ptr)[entry->occ++] = ival;
	}
    }
  else if ( entry->type == PML_FLT )
    {
      if ( entry->occ < (int) entry->size )
	((double *) entry->ptr)[entry->occ++] = atof(arg);
    }
  else if ( entry->type == PML_WORD )
    {
      if ( entry->occ < (int) entry->size )
	((char **) entry->ptr)[entry->occ++] = strdupx(arg);
    }
  else
    {
      fprintf(stderr, "unsupported type!\n");
    }

  return status;
}


void pmlProcess(pml_entry_t *entry, int argc, char **argv)
{
  char *parg;
  char *epos;

  for ( int i = 0; i < argc; ++i )
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

      pml_add_entry(entry, parg);
    }
}


int pmlRead(pml_t *pml, int argc, char **argv)
{
  pml_entry_t *entry = NULL;
  pml_entry_t *pentry[MAX_PML_ENTRY];
  int params[MAX_PML_ENTRY];
  int num_par[MAX_PML_ENTRY];
  int nparams = 0;
  int i;
  char *epos;
  size_t len;
  int bufsize = 0;
  int status = 0;
  /*
  if ( cdoVerbose )
    for ( i = 0; i < argc; ++i ) printf("pmlRead: %d %s\n", i, argv[i]);
  */
  for ( i = 0; i < argc; ++i )
    {
      len = strlen(argv[i]);
      bufsize += len+1;
    }

  char *parbuf = (char*) Malloc(bufsize*sizeof(char));
  memset(parbuf, 0, bufsize*sizeof(char));

  int istart = 0;
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
      for ( i = 0; i < pml->size; ++i )
	{
	  entry = pml->entry[i];
	  if ( entry->len == len )
	    if ( memcmp(entry->name, argv[istart], len) == 0 ) break;
	}

      if ( i == pml->size )
	{
	  fprintf(stderr, "Parameter >%s< has an invalid keyword!\n", argv[istart]);
	  status = 2;
	  goto END_LABEL;
	}

      num_par[nparams] = 0;
      pentry[nparams]  = entry;
      params[nparams]  = istart;
      num_par[nparams] = 1;
      
      istart++;
      for ( i = istart; i < argc; ++i )
	{
	  if ( *argv[i] == 0 ) { i++; break;}
	  epos = strchr(argv[i], '=');
	  if ( epos != NULL ) break;

	  num_par[nparams]++;
	}

      istart = i;

      nparams++;
    }

  for ( i = 0; i < nparams; ++i )
    {
      pmlProcess(pentry[i], num_par[i], &argv[params[i]]);
    }


 END_LABEL:

  Free(parbuf);

  return status;
}


bool par_check_int(int npar, int *parlist, bool *flaglist, int par)
{
  bool found = false;
  for ( int i = 0; i < npar; i++ )
    if ( par == parlist[i] ) { found = true; flaglist[i] = true;/* break;*/}

  return found;
}


bool par_check_flt(int npar, double *parlist, bool *flaglist, double par)
{
  bool found = false;
  for ( int i = 0; i < npar; i++ )
    if ( fabs(par - parlist[i]) < 1.e-4 ) { found = true; flaglist[i] = true;/* break;*/}

  return found;
}


bool par_check_word(int npar, char **parlist, bool *flaglist, const char *par)
{
  bool found = false;
  for ( int i = 0; i < npar; i++ )
    if ( wildcardmatch(parlist[i], par) == 0 ) { found = true; flaglist[i] = true;/* break;*/}

  return found;
}


bool par_check_date(int npar, char **parlist, bool *flaglist, const char *par)
{
  bool found = false;
  char wcdate[512];

  if ( *par == ' ' ) ++par;

  for ( int i = 0; i < npar; i++ )
    {
      strcpy(wcdate, parlist[i]);
      strcat(wcdate, "*");
      if ( wildcardmatch(wcdate, par) == 0 ) { found = true; flaglist[i] = true;/* break;*/}
    }

  return found;
}

void season_to_months(const char *season, int *imonths);

bool par_check_season(int npar, char **parlist, bool *flaglist, int month)
{
  assert(month>=1&&month<=12);
  bool found = false;
  int imon[13]; /* 1-12 ! */

  for ( int i = 0; i < npar; i++ )
    {
      for ( int i = 0; i < 13; ++i ) imon[i] = 0;
      season_to_months(parlist[i], imon);
      if ( imon[month] ) { found = true; flaglist[i] = true;/* break;*/}
    }

  return found;
}


void par_check_int_flag(int npar, int *parlist, bool *flaglist, const char *txt)
{
  for ( int i = 0; i < npar; ++i )
    if ( flaglist[i] == false )
      cdoWarning("%s >%d< not found!", txt, parlist[i]);
}


void par_check_flt_flag(int npar, double *parlist, bool *flaglist, const char *txt)
{
  for ( int i = 0; i < npar; ++i )
    if ( flaglist[i] == false )
      cdoWarning("%s >%g< not found!", txt, parlist[i]);
}


void par_check_word_flag(int npar, char **parlist, bool *flaglist, const char *txt)
{
  for ( int i = 0; i < npar; ++i )
    if ( flaglist[i] == false )
      cdoWarning("%s >%s< not found!", txt, parlist[i]);
}


int vlist_get_psvarid(int vlistID, int zaxisID)
{
  int psvarid = -1;
  char name[CDI_MAX_NAME];
  char psname[CDI_MAX_NAME];
  psname[0] = 0;
  zaxisInqPsName(zaxisID, psname);

  if ( psname[0] )
    {
      int nvars = vlistNvars(vlistID);
      for ( int varID = 0; varID < nvars; ++varID )
        {
          vlistInqVarName(vlistID, varID, name);
          if ( strcmp(name, psname) == 0 )
            {
              psvarid = varID;
              break;
            }
        }
      if ( cdoVerbose && psvarid == -1 )
        cdoWarning("Surface pressure variable not found - %s", psname);
    }

  return psvarid;
}


void *Select(void *argument)
{
  int streamID2 = CDI_UNDEFID;
  int nrecs;
  int nvars, nvars2, nlevs;
  int zaxisID;
  int varID2, levelID2;
  int varID, levelID;
  int last_year = -999999999;
  char paramstr[32];
  char varname[CDI_MAX_NAME];
  char stdname[CDI_MAX_NAME];
  char gname[CDI_MAX_NAME];
  char zname[CDI_MAX_NAME];
  int vlistID0 = -1, vlistID2 = -1;
  int result = FALSE;
  int taxisID2 = CDI_UNDEFID;
  int ntsteps;
  bool ltimsel = false;
  bool *vars = NULL;
  double **vardata2 = NULL;
  double *array = NULL;
  double fstartdate = -99999999999.;
  double fenddate   = -99999999999.;

  PML_DEF_INT(timestep_of_year, 4096, "Timestep of year");
  PML_DEF_INT(timestep,         4096, "Timestep");
  PML_DEF_INT(year,             1024, "Year");
  PML_DEF_INT(month,              32, "Month");
  PML_DEF_INT(day,                32, "Day");
  PML_DEF_INT(hour,               24, "Hour");
  PML_DEF_INT(minute,             60, "Minute");
  PML_DEF_INT(code,             1024, "Code number");
  PML_DEF_INT(levidx,           1024, "Level index");
  PML_DEF_INT(ltype,             256, "Level type");
  PML_DEF_INT(zaxisnum,          256, "Zaxis number");
  PML_DEF_INT(gridnum,           256, "Grid number");
  PML_DEF_FLT(level,            1024, "Level");
  PML_DEF_WORD(name,            1024, "Variable name");
  PML_DEF_WORD(param,           1024, "Parameter");
  PML_DEF_WORD(zaxisname,        256, "Zaxis name");
  PML_DEF_WORD(gridname,         256, "Grid name");
  PML_DEF_WORD(steptype,          32, "Time step type");
  PML_DEF_WORD(startdate,          1, "Start date");
  PML_DEF_WORD(enddate,            1, "End date");
  PML_DEF_WORD(season,            12, "Season");
  PML_DEF_WORD(date,            1024, "Date");

  PML_INIT_INT(timestep_of_year);
  PML_INIT_INT(timestep);
  PML_INIT_INT(year);
  PML_INIT_INT(month);
  PML_INIT_INT(day);
  PML_INIT_INT(hour);
  PML_INIT_INT(minute);
  PML_INIT_INT(code);
  PML_INIT_INT(levidx);
  PML_INIT_INT(ltype);
  PML_INIT_INT(zaxisnum);
  PML_INIT_INT(gridnum);
  PML_INIT_FLT(level);
  PML_INIT_WORD(name);
  PML_INIT_WORD(param);
  PML_INIT_WORD(zaxisname);
  PML_INIT_WORD(gridname);
  PML_INIT_WORD(steptype);
  PML_INIT_WORD(startdate);
  PML_INIT_WORD(enddate);
  PML_INIT_WORD(season);
  PML_INIT_WORD(date);

  cdoInitialize(argument);

  int SELECT = cdoOperatorAdd("select", 0, 0, "parameter list");
  int DELETE = cdoOperatorAdd("delete", 0, 0, "parameter list");

  bool lcopy = false;
  if ( UNCHANGED_RECORD ) lcopy = true;

  int operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  int nsel = operatorArgc();
  char **argnames = operatorArgv();

  if ( cdoVerbose )
    for ( int i = 0; i < nsel; i++ )
      printf("name %d = %s\n", i+1, argnames[i]);

  pml_t *pml = pmlNew("SELECT");

  PML_ADD_INT(pml, timestep_of_year);
  PML_ADD_INT(pml, timestep);
  PML_ADD_INT(pml, year);
  PML_ADD_INT(pml, month);
  PML_ADD_INT(pml, day);
  PML_ADD_INT(pml, hour);
  PML_ADD_INT(pml, minute);
  PML_ADD_INT(pml, code);
  PML_ADD_INT(pml, levidx);
  PML_ADD_INT(pml, ltype);
  PML_ADD_INT(pml, zaxisnum);
  PML_ADD_INT(pml, gridnum);
  PML_ADD_FLT(pml, level);
  PML_ADD_WORD(pml, name);
  PML_ADD_WORD(pml, param);
  PML_ADD_WORD(pml, zaxisname);
  PML_ADD_WORD(pml, gridname);
  PML_ADD_WORD(pml, steptype);
  PML_ADD_WORD(pml, startdate);
  PML_ADD_WORD(pml, enddate);
  PML_ADD_WORD(pml, season);
  PML_ADD_WORD(pml, date);

  pmlRead(pml, nsel, argnames);

  if ( cdoVerbose ) pmlPrint(pml);

  PML_NUM(pml, timestep_of_year);
  PML_NUM(pml, timestep);
  PML_NUM(pml, year);
  PML_NUM(pml, month);
  PML_NUM(pml, day);
  PML_NUM(pml, hour);
  PML_NUM(pml, minute);
  PML_NUM(pml, code);
  PML_NUM(pml, levidx);
  PML_NUM(pml, ltype);
  PML_NUM(pml, zaxisnum);
  PML_NUM(pml, gridnum);
  PML_NUM(pml, level);
  PML_NUM(pml, name);
  PML_NUM(pml, param);
  PML_NUM(pml, zaxisname);
  PML_NUM(pml, gridname);
  PML_NUM(pml, steptype);
  PML_NUM(pml, startdate);
  PML_NUM(pml, enddate);
  PML_NUM(pml, season);
  PML_NUM(pml, date);

  int streamCnt = cdoStreamCnt();
  int nfiles = streamCnt - 1;

  if ( !cdoVerbose && nfiles > 1 ) progressInit();

  int tsID2 = 0;
  for ( int indf = 0; indf < nfiles; indf++ )
    {
      if ( !cdoVerbose && nfiles > 1 ) progressStatus(0, 1, (indf+1.)/nfiles);
      if ( cdoVerbose ) cdoPrint("Process file: %s", cdoStreamName(indf)->args);

      int streamID1 = streamOpenRead(cdoStreamName(indf));

      int vlistID1 = streamInqVlist(streamID1);
      int taxisID1 = vlistInqTaxis(vlistID1);

      bool lcopy_const = false;

      if ( indf == 0 )
	{
	  // vlistID0 = vlistDuplicate(vlistID1);

	  vlistClearFlag(vlistID1);
	  nvars = vlistNvars(vlistID1);
	  vars  = (bool*) Malloc(nvars*sizeof(bool));

	  if ( operatorID == DELETE )
	    {
	      result = FALSE;
	      for ( int varID = 0; varID < nvars; varID++ )
		{
		  zaxisID = vlistInqVarZaxis(vlistID1, varID);
		  nlevs   = zaxisInqSize(zaxisID);
		  for ( int levID = 0; levID < nlevs; levID++ )
		    vlistDefFlag(vlistID1, varID, levID, TRUE);
		}
	    }
	  else if ( operatorID == SELECT )
	    {
	      result = TRUE;
	    }

	  for ( int varID = 0; varID < nvars; varID++ )
	    {
	      int iparam = vlistInqVarParam(vlistID1, varID);
	      code = vlistInqVarCode(vlistID1, varID);
	      vlistInqVarName(vlistID1, varID, varname);
	      vlistInqVarStdname(vlistID1, varID, stdname);

	      cdiParamToString(iparam, paramstr, sizeof(paramstr));

	      name  = varname;
	      param = paramstr;

	      int gridID  = vlistInqVarGrid(vlistID1, varID);
	      zaxisID = vlistInqVarZaxis(vlistID1, varID);
	      nlevs   = zaxisInqSize(zaxisID);
	      ltype   = zaxis2ltype(zaxisID);

              zaxisnum = vlistZaxisIndex(vlistID1, zaxisID)+1;
              int zaxistype = zaxisInqType(zaxisID);
              zaxisName(zaxistype, zname);
              zaxisname = zname;

              gridnum = vlistGridIndex(vlistID1, gridID)+1;
              int gridtype = gridInqType(gridID);
              gridName(gridtype, gname);
              gridname = gname;

              int tsteptype = vlistInqVarTsteptype(vlistID1, varID);
              if      ( tsteptype == TSTEP_CONSTANT ) steptype = "constant";
              else if ( tsteptype == TSTEP_INSTANT  ) steptype = "instant";
              else if ( tsteptype == TSTEP_INSTANT2 ) steptype = "instant";
              else if ( tsteptype == TSTEP_INSTANT3 ) steptype = "instant";
              else if ( tsteptype == TSTEP_MIN      ) steptype = "min";
              else if ( tsteptype == TSTEP_MAX      ) steptype = "max";
              else if ( tsteptype == TSTEP_AVG      ) steptype = "avg";
              else if ( tsteptype == TSTEP_ACCUM    ) steptype = "accum";
              else if ( tsteptype == TSTEP_RANGE    ) steptype = "range";
              else if ( tsteptype == TSTEP_DIFF     ) steptype = "diff";
              else                                    steptype = "unknown";
              
	      vars[varID] = false;
              bool found_code  = npar_code      && PAR_CHECK_INT(code);
              bool found_name  = npar_name      && PAR_CHECK_WORD(name);
              bool found_param = npar_param     && PAR_CHECK_WORD(param);
              bool found_grid  = npar_gridnum   && PAR_CHECK_INT(gridnum);
              bool found_gname = npar_gridname  && PAR_CHECK_WORD(gridname);
              bool found_stype = npar_steptype  && PAR_CHECK_WORD(steptype);
              bool found_ltype = npar_ltype     && PAR_CHECK_INT(ltype);
              bool found_zaxis = npar_zaxisnum  && PAR_CHECK_INT(zaxisnum);
              bool found_zname = npar_zaxisname && PAR_CHECK_WORD(zaxisname);
              bool lvar  = found_code || found_name || found_param;
              bool lstep = npar_steptype ? found_stype : true;
              bool lgrid = (npar_gridnum || npar_gridname) ? (found_grid || found_gname) : true;
              bool lvert = (npar_ltype || npar_zaxisnum || npar_zaxisname) ? (found_ltype || found_zaxis || found_zname) : true;
	      
              if ( !vars[varID] && lgrid && lvar) vars[varID] = true;
              if ( !vars[varID] && lvert && lvar) vars[varID] = true;
              if ( !vars[varID] && lstep && lvar) vars[varID] = true;
              if ( !vars[varID] && !lvar )
                {
                  if      ( found_grid || found_gname ) vars[varID] = true;
                  else if ( found_stype ) vars[varID] = true;
                  else if ( found_ltype || found_zaxis || found_zname ) vars[varID] = true;
                  else if ( npar_levidx || npar_level )
                    {
                      for ( int levID = 0; levID < nlevs; levID++ )
                        {
                          levidx = levID + 1;
                          level = zaxisInqLevel(zaxisID, levID);
                          if ( !vars[varID] && npar_levidx && PAR_CHECK_INT(levidx) )  vars[varID] = true;
                          if ( !vars[varID] && npar_level  && PAR_CHECK_FLT(level)  )  vars[varID] = true;
                        }
                    }
                }
	    }

	  for ( int varID = 0; varID < nvars; varID++ )
	    {
	      if ( vars[varID] )
		{
		  int zaxisID = vlistInqVarZaxis(vlistID1, varID);
                  if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID )
                    {
                      int psvarid = vlist_get_psvarid(vlistID1, zaxisID);
                      if ( psvarid != -1 && !vars[psvarid] ) vars[psvarid] = true;
                    }
                }
            }

	  for ( int varID = 0; varID < nvars; varID++ )
	    {
	      if ( vars[varID] )
		{
		  zaxisID = vlistInqVarZaxis(vlistID1, varID);
		  nlevs   = zaxisInqSize(zaxisID);

		  for ( int levID = 0; levID < nlevs; levID++ )
		    {
		      levidx = levID + 1;
		      level = zaxisInqLevel(zaxisID, levID);
		      
		      if ( nlevs == 1 && IS_EQUAL(level, 0) )
			{
			  vlistDefFlag(vlistID1, varID, levID, result);
			}
		      else
			{
			  if ( npar_levidx )
			    {
			      if ( PAR_CHECK_INT(levidx) )
				vlistDefFlag(vlistID1, varID, levID, result);
			    }
			  else if ( npar_level )
			    {
			      if ( PAR_CHECK_FLT(level) )
				vlistDefFlag(vlistID1, varID, levID, result);
			    }
			  else
			    {
			      vlistDefFlag(vlistID1, varID, levID, result);
			    }
			}
		    }
		}
	    }

	  PAR_CHECK_INT_FLAG(code);
	  PAR_CHECK_INT_FLAG(levidx);
	  PAR_CHECK_INT_FLAG(ltype);
	  PAR_CHECK_INT_FLAG(zaxisnum);
	  PAR_CHECK_INT_FLAG(gridnum);
	  PAR_CHECK_FLT_FLAG(level);
	  PAR_CHECK_WORD_FLAG(name);
	  PAR_CHECK_WORD_FLAG(param);
	  PAR_CHECK_WORD_FLAG(zaxisname);
	  PAR_CHECK_WORD_FLAG(gridname);
	  PAR_CHECK_WORD_FLAG(steptype);

	  if ( npar_date || npar_startdate || npar_enddate || npar_season ) ltimsel = true;
	  if ( npar_timestep_of_year || npar_timestep || npar_year || npar_month || npar_day || npar_hour || npar_minute ) ltimsel = true;

	  int npar = 0;
	  for ( int varID = 0; varID < nvars; varID++ )
	    {
	      zaxisID = vlistInqVarZaxis(vlistID1, varID);
	      nlevs   = zaxisInqSize(zaxisID);
	      for ( int levID = 0; levID < nlevs; levID++ )
		if ( vlistInqFlag(vlistID1, varID, levID) == result )
                  {
                    npar++;
                    break;
                  }
	    }

	  if ( npar == 0 )
	    {
	      if ( ltimsel == true )
		{
                  lcopy_const = true;

		  for ( int varID = 0; varID < nvars; varID++ )
		    {
		      vars[varID] = true;
		      zaxisID = vlistInqVarZaxis(vlistID1, varID);
		      nlevs   = zaxisInqSize(zaxisID);
		      for ( int levID = 0; levID < nlevs; levID++ )
			vlistDefFlag(vlistID1, varID, levID, TRUE);
		    }
		}
	      else
		{
		  cdoAbort("No variable selected!");
		}
	    }

	  //if ( cdoVerbose ) vlistPrint(vlistID1);

	  vlistID0 = vlistDuplicate(vlistID1);
	  for ( int varID = 0; varID < nvars; varID++ )
	    {
	      zaxisID = vlistInqVarZaxis(vlistID1, varID);
	      nlevs   = zaxisInqSize(zaxisID);
	      for ( int levID = 0; levID < nlevs; levID++ )
		vlistDefFlag(vlistID0, varID, levID, vlistInqFlag(vlistID1, varID, levID));
	    }

	  //if ( cdoVerbose ) vlistPrint(vlistID0);

	  vlistID2 = vlistCreate();
	  vlistCopyFlag(vlistID2, vlistID0);

	  //if ( cdoVerbose ) vlistPrint(vlistID2);

	  taxisID2 = taxisDuplicate(taxisID1);
	  vlistDefTaxis(vlistID2, taxisID2);

	  ntsteps = vlistNtsteps(vlistID1);

	  nvars2 = vlistNvars(vlistID2);

	  if ( ntsteps == 1 && nfiles == 1 )
	    {
	      for ( varID = 0; varID < nvars2; ++varID )
		if ( vlistInqVarTsteptype(vlistID2, varID) != TSTEP_CONSTANT ) break;

	      if ( varID == nvars2 ) ntsteps = 0;
	    }

	  int ntsteps2 = ntsteps;
	  if ( operatorID == SELECT && npar_timestep == 1 ) ntsteps2 = 1;
	  
	  if ( ntsteps2 == 0 || ntsteps2 == 1 ) vlistDefNtsteps(vlistID2, ntsteps2);

	  if ( ntsteps2 == 0 && nfiles > 1 )
	    {	      
	      for ( int varID = 0; varID < nvars2; ++varID )
		vlistDefVarTsteptype(vlistID2, varID, TSTEP_INSTANT);
	    }

	  // support for negative timestep values
	  if ( npar_timestep > 0 && ntsteps > 0 && nfiles == 1 )
	    {
	      for ( int i = 0; i < npar_timestep; i++ )
		{
		  if ( par_timestep[i] < 0 )
		    {
		      if ( cdoVerbose )
			cdoPrint("timestep %d changed to %d", par_timestep[i], par_timestep[i] + ntsteps + 1 );
		      par_timestep[i] += ntsteps + 1;
		    }
		}
	    }

	  if ( ! lcopy )
	    {
	      int gridsize = vlistGridsizeMax(vlistID1);
	      if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
	      array = (double*) Malloc(gridsize*sizeof(double));
	    }

	  startdate = par_startdate[0];
	  enddate   = par_enddate[0];
	  if ( npar_startdate ) fstartdate = datestr_to_double(startdate, 0);
	  if ( npar_enddate   ) fenddate   = datestr_to_double(enddate, 1);
	}
      else
	{
	  vlistCompare(vlistID0, vlistID1, CMP_ALL);
	}


      if ( nvars2 == 0 )
	{
	  cdoWarning("No resulting variables available!");
	  goto END_LABEL;
	}

      if ( lcopy_const )
        {
          vardata2 = (double**) malloc(nvars2*sizeof(double));
          for ( int varID = 0; varID < nvars2; ++varID ) vardata2[varID] = NULL;
        }

      bool lstop = false;
      int tsID1 = 0;
      while ( (nrecs = streamInqTimestep(streamID1, tsID1)) )
	{
	  bool copytimestep = true;

	  if ( ltimsel == true )
	    {
	      copytimestep = false;
	      timestep = tsID1 + 1;

	      if ( operatorID == SELECT && npar_timestep > 0 && timestep > par_timestep[npar_timestep-1] )
		{
		  lstop = true;
		  break;
		}

	      int vdate = taxisInqVdate(taxisID1);
	      int vtime = taxisInqVtime(taxisID1);
              int second;
	      cdiDecodeDate(vdate, &year, &month, &day);
	      cdiDecodeTime(vtime, &hour, &minute, &second);
              UNUSED(season);

	      if ( year != last_year )
		{
		  timestep_of_year = 0;
		  last_year = year;
		}

	      timestep_of_year++;

	      if ( npar_timestep && PAR_CHECK_INT(timestep) ) copytimestep = true;
	      if ( npar_timestep_of_year && PAR_CHECK_INT(timestep_of_year) ) copytimestep = true;

	      if ( !copytimestep && npar_date == 0 && npar_timestep == 0 && npar_timestep_of_year == 0 )
		{
		  bool lseason = false, lyear = false, lmonth = false, lday = false, lhour = false, lminute = false;

		  if ( npar_season == 0 || (npar_season && PAR_CHECK_SEASON(season, month)) ) lseason   = true;
		  if ( npar_year   == 0 || (npar_year   && PAR_CHECK_INT(year))   ) lyear   = true;
		  if ( npar_month  == 0 || (npar_month  && PAR_CHECK_INT(month))  ) lmonth  = true;
		  if ( npar_day    == 0 || (npar_day    && PAR_CHECK_INT(day))    ) lday    = true;
		  if ( npar_hour   == 0 || (npar_hour   && PAR_CHECK_INT(hour))   ) lhour   = true;
		  if ( npar_minute == 0 || (npar_minute && PAR_CHECK_INT(minute)) ) lminute = true;

		  if ( lseason && lyear && lmonth && lday && lhour && lminute ) copytimestep = true;
		}

	      double fdate = ((double)vdate) + ((double)vtime)/1000000.;

	      if ( npar_enddate )
		{
		  if ( fdate > fenddate )
		    {
		      flag_enddate[0] = true;
		      copytimestep = false;
		      if ( operatorID == SELECT )
			{
			  lstop = true;
			  break;
			}
		    }
		  else
		    {
		      copytimestep = true;
		    }
		}

	      if ( npar_startdate )
		{
		  if ( fdate < fstartdate )
		    {
		      copytimestep = false;
		    }
		  else
		    {
		      flag_startdate[0] = true;
		      copytimestep = true;
		    }
		}

              
              if ( npar_date )
                {
                  char vdatetimestr[64];
                  datetime2str(vdate, vtime, vdatetimestr, sizeof(vdatetimestr));
                  date = vdatetimestr;
                  if ( PAR_CHECK_DATE(date) ) copytimestep = true;
                }

	      if ( operatorID == DELETE ) copytimestep = !copytimestep;

              if ( copytimestep && indf == 0 && tsID1 == 0 ) lcopy_const = false;
	    }

	  if ( copytimestep == true )
	    {
	      if ( streamID2 == CDI_UNDEFID )
		{
		  streamID2 = streamOpenWrite(cdoStreamName(nfiles), cdoFiletype());
		  streamDefVlist(streamID2, vlistID2);
		}

	      taxisCopyTimestep(taxisID2, taxisID1);
	      streamDefTimestep(streamID2, tsID2);

              if ( lcopy_const && tsID2 == 0 )
                {
                  for ( varID2 = 0; varID2 < nvars2; ++varID2 )
                    {
                      if ( vardata2[varID2] )
                        {
                          double missval = vlistInqVarMissval(vlistID2, varID2);
                          int gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID2));
                          int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID2));
                          for ( int levelID2 = 0; levelID2 < nlevel; ++levelID2 )
                            {
                              double *pdata = vardata2[varID2]+gridsize*levelID;
                              int nmiss = 0;
                              for ( int i = 0; i < gridsize; ++i )
                                if ( DBL_IS_EQUAL(pdata[i], missval) ) nmiss++;

                              streamDefRecord(streamID2, varID2, levelID2);
                              streamWriteRecord(streamID2, pdata, nmiss);
                            }
                        }
                    }
                }
              
	      for ( int recID = 0; recID < nrecs; recID++ )
		{
		  streamInqRecord(streamID1, &varID, &levelID);
		  if ( vlistInqFlag(vlistID0, varID, levelID) == TRUE )
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
                          int nmiss;
			  streamReadRecord(streamID1, array, &nmiss);
			  streamWriteRecord(streamID2, array, nmiss);
			}
		    }
		}

	      tsID2++;              
	    }
          else if ( lcopy_const && indf == 0 && tsID1 == 0 )
            {
	      for ( int recID = 0; recID < nrecs; recID++ )
		{
		  streamInqRecord(streamID1, &varID, &levelID);
		  if ( vlistInqFlag(vlistID0, varID, levelID) == TRUE )
		    {
		      varID2 = vlistFindVar(vlistID2, varID);
                      if ( vlistInqVarTsteptype(vlistID2, varID2) == TSTEP_CONSTANT )
                        {
                          levelID2 = vlistFindLevel(vlistID2, varID, levelID);
                          int gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID2));
                          if ( levelID == 0 )
                            {
                              int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID2));
                              vardata2[varID2] = (double*) malloc(gridsize*nlevel*sizeof(double));
                            }
                          int nmiss;
                          streamReadRecord(streamID1, vardata2[varID2]+gridsize*levelID, &nmiss);
                        }
		    }
		}
            }

	  tsID1++;
	}
      
      streamClose(streamID1);

      if ( lstop ) break;
    }

 END_LABEL:

  if ( !cdoVerbose && nfiles > 1 ) progressStatus(0, 1, 1);    

  PAR_CHECK_INT_FLAG(timestep_of_year);
  PAR_CHECK_INT_FLAG(timestep);
  PAR_CHECK_INT_FLAG(year);
  PAR_CHECK_INT_FLAG(month);
  PAR_CHECK_INT_FLAG(day);
  PAR_CHECK_INT_FLAG(hour);
  PAR_CHECK_INT_FLAG(minute);
  PAR_CHECK_WORD_FLAG(startdate);
  UNUSED(str_enddate);
  //  PAR_CHECK_WORD_FLAG(enddate);
  PAR_CHECK_WORD_FLAG(season);
  PAR_CHECK_WORD_FLAG(date);

  if ( streamID2 != CDI_UNDEFID ) streamClose(streamID2);

  vlistDestroy(vlistID0);
  vlistDestroy(vlistID2);

  pmlDestroy(pml);

  if ( array ) Free(array);
  if ( vars ) Free(vars);
  if ( vardata2 )
    {
      for ( int varID = 0; varID < nvars2; ++ varID )
        if ( vardata2[varID] ) free(vardata2[varID]);

      free(vardata2);
    }

  if ( tsID2 == 0 ) cdoAbort("No timesteps selected!");

  cdoFinish();

  return 0;
}

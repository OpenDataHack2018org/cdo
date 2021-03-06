/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2018 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
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

#include "cdo_int.h"
#include "pstream_int.h"
#include "grid.h"
#include "expr.h"
#include "datetime.h"

void grid_cell_area(int gridID, double *array);
int getSurfaceID(int vlistID);

static char *
exprs_from_arg(const char *arg)
{
  size_t slen = strlen(arg);
  char *exprs = (char *) Malloc(slen + 2);
  strcpy(exprs, operatorArgv()[0]);
  if (exprs[slen - 1] != ';')
    {
      exprs[slen] = ';';
      exprs[slen + 1] = 0;
    }

  return exprs;
}

static char *
exprs_from_file(const char *exprf)
{
  /* Open expr script file for reading */
  FILE *fp = fopen(exprf, "r");
  if (fp == NULL) cdoAbort("Open failed on %s", exprf);

  struct stat filestat;
  if (stat(exprf, &filestat) != 0) cdoAbort("Stat failed on %s", exprf);

  size_t fsize = (size_t) filestat.st_size;
  char *exprs = (char *) Malloc(fsize + 1);

  int ichar, ipos = 0;
  while ((ichar = fgetc(fp)) != EOF) exprs[ipos++] = ichar;

  exprs[ipos] = 0;
  if (ipos == 0) cdoAbort("%s is empty!", exprf);

  fclose(fp);

  return exprs;
}

#define MAX_PARAMS 4096

char *
str_replace(char *target, const char *needle, const char *replacement)
{
  char buffer[1024] = { 0 };
  char *insert_point = &buffer[0];
  const char *tmp = target;
  size_t needle_len = strlen(needle);
  size_t repl_len = strlen(replacement);

  while (1)
    {
      const char *p = strstr(tmp, needle);

      // walked past last occurrence of needle; copy remaining part
      if (p == NULL)
        {
          strcpy(insert_point, tmp);
          break;
        }

      // copy part before needle
      memcpy(insert_point, tmp, p - tmp);
      insert_point += p - tmp;

      // copy replacement string
      memcpy(insert_point, replacement, repl_len);
      insert_point += repl_len;

      // adjust pointers, move on
      tmp = p + needle_len;
    }

  // write altered string back to target
  strcpy(target, buffer);

  return target;
}

static char *
exprs_expand(char *exprs, int vlistID)
{
  size_t len = strlen(exprs);
  unsigned nequal = 0, nsemi = 0;
  for (size_t i = 0; i < len; ++i)
    {
      if (exprs[i] == '=') nequal++;
      if (exprs[i] == ';') nsemi++;
    }

  const char *needle = "_T";
  if (nequal == 0 && nsemi == 1 && strstr(exprs, needle))
    {
      char varname[CDI_MAX_NAME];
      int nvars = vlistNvars(vlistID);
      const size_t bufsize = 1024;
      char *sbuf = (char *) Malloc(bufsize);
      char *buf = (char *) Malloc(nvars * bufsize);
      buf[0] = 0;
      for (int varID = 0; varID < nvars; ++varID)
        {
          vlistInqVarName(vlistID, varID, varname);
          strcpy(sbuf, exprs);
          str_replace(sbuf, needle, varname);
          strcat(buf, varname);
          strcat(buf, "=");
          strcat(buf, sbuf);
        }
      Free(sbuf);
      Free(exprs);
      exprs = buf;
    }

  return exprs;
}

static paramType *
params_new(int vlistID)
{
  paramType *params = (paramType *) Malloc(MAX_PARAMS * sizeof(paramType));
  memset(params, 0, MAX_PARAMS * sizeof(paramType));

  int nvars1 = vlistNvars(vlistID);

  char name[CDI_MAX_NAME];
  char longname[CDI_MAX_NAME];
  char units[CDI_MAX_NAME];

  for (int varID = 0; varID < nvars1; varID++)
    {
      int gridID = vlistInqVarGrid(vlistID, varID);
      int zaxisID = vlistInqVarZaxis(vlistID, varID);
      int datatype = vlistInqVarDatatype(vlistID, varID);
      int steptype = vlistInqVarTimetype(vlistID, varID);
      size_t nlat = gridInqYsize(gridID);
      size_t ngp = gridInqSize(gridID);
      size_t nlev = zaxisInqSize(zaxisID);
      double missval = vlistInqVarMissval(vlistID, varID);

      vlistInqVarName(vlistID, varID, name);
      vlistInqVarLongname(vlistID, varID, longname);
      vlistInqVarUnits(vlistID, varID, units);

      params[varID].type = PARAM_VAR;
      params[varID].select = false;
      params[varID].remove = false;
      params[varID].lmiss = true;
      params[varID].coord = 0;
      params[varID].gridID = gridID;
      params[varID].zaxisID = zaxisID;
      params[varID].datatype = datatype;
      params[varID].steptype = steptype;
      params[varID].nlat = nlat;
      params[varID].ngp = ngp;
      params[varID].nlev = nlev;
      params[varID].missval = missval;
      params[varID].nmiss = 0;
      params[varID].data = NULL;
      params[varID].name = strdup(name);
      params[varID].longname = strdup(longname);
      params[varID].units = strdup(units);
    }

  return params;
}

static void
params_add_coord(parseParamType *parse_arg, int coord, int cdiID, size_t size, const char *units, const char *longname)
{
  int ncoords = parse_arg->ncoords;
  if (ncoords >= parse_arg->maxcoords) cdoAbort("Too many coordinates (limit=%d)", parse_arg->maxcoords);

  parse_arg->coords[ncoords].needed = false;
  parse_arg->coords[ncoords].coord = coord;
  parse_arg->coords[ncoords].cdiID = cdiID;
  parse_arg->coords[ncoords].size = size;
  parse_arg->coords[ncoords].units = NULL;
  parse_arg->coords[ncoords].longname = NULL;
  parse_arg->coords[ncoords].data = NULL;
  if (units) parse_arg->coords[ncoords].units = strdup(units);
  if (longname) parse_arg->coords[ncoords].longname = strdup(longname);

  parse_arg->ncoords++;
}

int
params_get_coordID(parseParamType *parse_arg, int coord, int cdiID)
{
  int ncoords = parse_arg->ncoords;
  for (int coordID = 0; coordID < ncoords; ++coordID)
    {
      if (parse_arg->coords[coordID].coord == coord && parse_arg->coords[coordID].cdiID == cdiID) return coordID;
    }

  cdoAbort("%s: coordinate %c not found!", __func__, coord);

  return -1;
}

static void
params_add_coordinates(int vlistID, parseParamType *parse_arg)
{
  char longname[CDI_MAX_NAME];
  char units[CDI_MAX_NAME];

  int ngrids = vlistNgrids(vlistID);
  for (int index = 0; index < ngrids; ++index)
    {
      int gridID = vlistGrid(vlistID, index);
      size_t size = gridInqSize(gridID);
      gridInqXunits(gridID, units);
      params_add_coord(parse_arg, 'x', gridID, size, units, "longitude");
      gridInqYunits(gridID, units);
      params_add_coord(parse_arg, 'y', gridID, size, units, "latitude");

      params_add_coord(parse_arg, 'a', gridID, size, "m^2", "grid cell area");
      params_add_coord(parse_arg, 'w', gridID, size, NULL, "grid cell area weights");
    }

  int nzaxis = vlistNzaxis(vlistID);
  for (int index = 0; index < nzaxis; ++index)
    {
      int zaxisID = vlistZaxis(vlistID, index);
      int size = zaxisInqSize(zaxisID);
      zaxisInqUnits(zaxisID, units);
      zaxisInqLongname(zaxisID, longname);
      params_add_coord(parse_arg, 'z', zaxisID, size, units, longname);
    }
}

static int
params_add_ts(parseParamType *parse_arg)
{
  int varID = -1;
  paramType *params = parse_arg->params;
  if (params)
    {
      varID = parse_arg->nparams;
      if (varID >= parse_arg->maxparams) cdoAbort("Too many parameter (limit=%d)", parse_arg->maxparams);

      params[varID].name = strdup("_ts");
      params[varID].gridID = parse_arg->pointID;
      params[varID].zaxisID = parse_arg->surfaceID;
      params[varID].steptype = TIME_VARYING;
      params[varID].ngp = CLEN;
      params[varID].nlev = 1;

      parse_arg->nparams++;
    }

  return varID;
}

static void
params_delete(paramType *params)
{
  if (params)
    {
      for (int varID = 0; varID < MAX_PARAMS; varID++)
        {
          if (params[varID].data) Free(params[varID].data);
          if (params[varID].name) Free(params[varID].name);
          if (params[varID].longname) Free(params[varID].longname);
          if (params[varID].units) Free(params[varID].units);
        }
      Free(params);
    }
}

static void
parseParamInit(parseParamType *parse_arg, int vlistID, int pointID, int zonalID, int surfaceID, paramType *params)
{
  int nvars = vlistNvars(vlistID);
  int ngrids = vlistNgrids(vlistID);
  int nzaxis = vlistNzaxis(vlistID);
  int maxcoords = ngrids * 4 + nzaxis;

  parse_arg->maxparams = MAX_PARAMS;
  parse_arg->nparams = nvars;
  parse_arg->nvars1 = nvars;
  parse_arg->init = true;
  parse_arg->debug = cdoVerbose != 0;
  parse_arg->params = params;
  parse_arg->pointID = pointID;
  parse_arg->zonalID = zonalID;
  parse_arg->surfaceID = surfaceID;
  parse_arg->needed = (bool *) Malloc(nvars * sizeof(bool));
  parse_arg->coords = (coordType *) Malloc(maxcoords * sizeof(coordType));
  parse_arg->maxcoords = maxcoords;
  parse_arg->ncoords = 0;
}

static
int genZonalID(int vlistID)
{
  int zonalID = -1;

  int ngrids = vlistNgrids(vlistID);
  for (int index = 0; index < ngrids; index++)
    {
      int gridID = vlistGrid(vlistID, index);
      int gridtype = gridInqType(gridID);
      if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_GENERIC)
        if (gridInqXsize(gridID) > 1 && gridInqYsize(gridID) >= 1)
          {
            zonalID = gridToZonal(gridID);
            break;
          }
    }
 
  return zonalID;
}

void *
Expr(void *process)
{
  cdoInitialize(process);

  void *scanner;
  int yy_scan_string(const char *str, void *scanner);
  yylex_init(&scanner);

  parseParamType parse_arg;
  yyset_extra(&parse_arg, scanner);

  // clang-format off
  cdoOperatorAdd("expr",   1, 1, "expressions");
  cdoOperatorAdd("exprf",  1, 0, "expr script filename");
  cdoOperatorAdd("aexpr",  0, 1, "expressions");
  cdoOperatorAdd("aexprf", 0, 0, "expr script filename");
  // clang-format on

  int operatorID = cdoOperatorID();
  bool replacesVariables = cdoOperatorF1(operatorID);
  bool readsCommandLine = cdoOperatorF2(operatorID);

  operatorInputArg(cdoOperatorEnter(operatorID));

  char *exprs = readsCommandLine ? exprs_from_arg(operatorArgv()[0]) : exprs_from_file(operatorArgv()[0]);

  int streamID1 = cdoStreamOpenRead(0);
  int vlistID1 = cdoStreamInqVlist(streamID1);

  exprs = exprs_expand(exprs, vlistID1);
  if (cdoVerbose) cdoPrint(exprs);

  int nvars1 = vlistNvars(vlistID1);

  int pointID = gridCreate(GRID_GENERIC, 1);
  int zonalID = genZonalID(vlistID1);
  int surfaceID = getSurfaceID(vlistID1);

  paramType *params = params_new(vlistID1);

  parseParamInit(&parse_arg, vlistID1, pointID, zonalID, surfaceID, params);

  // Set all input variables to 'needed' if replacing is switched off
  for (int varID = 0; varID < nvars1; varID++) parse_arg.needed[varID] = !replacesVariables;

  int vartsID = params_add_ts(&parse_arg);
  parse_arg.tsID = vartsID;
  params_add_coordinates(vlistID1, &parse_arg);

  CDO_parser_errorno = 0;
  yy_scan_string(exprs, scanner);
  yyparse(&parse_arg, scanner);
  if (CDO_parser_errorno != 0) cdoAbort("Syntax error!");

  parse_arg.init = false;

  if (cdoVerbose)
    for (int varID = 0; varID < nvars1; varID++)
      if (parse_arg.needed[varID]) cdoPrint("Needed var: %d %s", varID, params[varID].name);

  for (int varID = 0; varID < parse_arg.nparams; varID++)
    {
      if (params[varID].ngp == 0 && params[varID].nlev == 0)
        {
          params[varID].ngp = 1;
          params[varID].nlev = 1;
          params[varID].type = PARAM_CONST;
        }
    }

  if (cdoVerbose)
    for (int varID = 0; varID < parse_arg.nparams; varID++)
      cdoPrint("var: %d %s ngp=%zu nlev=%zu coord=%c", varID, params[varID].name, params[varID].ngp, params[varID].nlev,
               params[varID].coord == 0 ? ' ' : params[varID].coord);

  int *varIDmap = (int *) Malloc(parse_arg.nparams * sizeof(int));

  int vlistID2 = vlistCreate();
  if (!replacesVariables)
    {
      vlistClearFlag(vlistID1);
      int pidx = 0;
      for (int varID = 0; varID < nvars1; varID++)
        {
          params[varID].select = false;
          if (params[varID].remove == false)
            {
              varIDmap[pidx++] = varID;
              int nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
              for (int levID = 0; levID < nlevs; levID++) vlistDefFlag(vlistID1, varID, levID, TRUE);
            }
        }
      cdoVlistCopyFlag(vlistID2, vlistID1);
    }

  for (int pidx = 0; pidx < parse_arg.nparams; pidx++)
    {
      if (pidx < nvars1 && params[pidx].select == false) continue;
      if (pidx >= nvars1 && params[pidx].type == PARAM_CONST) continue;
      if (pidx >= nvars1 && params[pidx].name[0] == '_') continue;
      if (pidx >= nvars1 && params[pidx].remove == true) continue;
      if (pidx >= nvars1 && params[pidx].coord) continue;

      int varID = vlistDefVar(vlistID2, params[pidx].gridID, params[pidx].zaxisID, params[pidx].steptype);
      vlistDefVarName(vlistID2, varID, params[pidx].name);
      if (params[pidx].lmiss) vlistDefVarMissval(vlistID2, varID, params[pidx].missval);
      if (params[pidx].units) vlistDefVarUnits(vlistID2, varID, params[pidx].units);
      if (params[pidx].longname) vlistDefVarLongname(vlistID2, varID, params[pidx].longname);
      if (memcmp(params[pidx].name, "var", 3) == 0)
        {
          if (strlen(params[pidx].name) > 3 && isdigit(params[pidx].name[3]))
            {
              int code = atoi(params[pidx].name + 3);
              vlistDefVarCode(vlistID2, varID, code);
            }
        }
      varIDmap[varID] = pidx;
    }

  if (cdoVerbose)
    {
      for (int varID = 0; varID < nvars1; varID++)
        if (parse_arg.needed[varID]) printf("needed: %d %s\n", varID, parse_arg.params[varID].name);
      cdoPrint("vlistNvars(vlistID1)=%d, vlistNvars(vlistID2)=%d", vlistNvars(vlistID1), vlistNvars(vlistID2));
    }

  int nvars2 = vlistNvars(vlistID2);
  if (nvars2 == 0) cdoAbort("No output variable found!");

  for (int varID = 0; varID < nvars1; varID++)
    {
      if (parse_arg.needed[varID])
        {
          size_t ngp = params[varID].ngp;
          size_t nlev = params[varID].nlev;
          size_t nItems = MAX(4, ngp * nlev);
          params[varID].data = (double *) Malloc(nItems * sizeof(double));
        }
    }

  for (int varID = parse_arg.nvars1; varID < parse_arg.nparams; varID++)
    {
      size_t ngp = params[varID].ngp;
      size_t nlev = params[varID].nlev;
      size_t nItems = MAX(4, ngp * nlev);
      params[varID].data = (double *) Malloc(nItems * sizeof(double));
    }

  for (int i = 0; i < parse_arg.ncoords; i++)
    {
      if (parse_arg.coords[i].needed)
        {
          int coord = parse_arg.coords[i].coord;
          if (coord == 'x' || coord == 'y' || coord == 'a' || coord == 'w')
            {
              int gridID = parse_arg.coords[i].cdiID;
              size_t ngp = parse_arg.coords[i].size;
              double *data = (double *) Malloc(ngp * sizeof(double));
              parse_arg.coords[i].data = data;
              if (coord == 'x' || coord == 'y')
                {
                  if (gridInqType(gridID) == GRID_GENERIC) cdoAbort("Grid has no geographical coordinates!");
                  if (gridInqType(gridID) == GRID_GME) gridID = gridToUnstructured(gridID, 0);
                  if (gridInqType(gridID) != GRID_UNSTRUCTURED && gridInqType(gridID) != GRID_CURVILINEAR)
                    gridID = gridToCurvilinear(gridID, 0);

                  if (coord == 'x') gridInqXvals(gridID, data);
                  if (coord == 'y') gridInqYvals(gridID, data);

                  if (gridID != parse_arg.coords[i].cdiID) gridDestroy(gridID);
                }
              else if (coord == 'a')
                {
                  grid_cell_area(gridID, data);
                }
              else if (coord == 'w')
                {
                  data[0] = 1;
                  if (ngp > 1)
                    {
                      int wstatus = gridWeights(gridID, data);
                      if (wstatus) cdoWarning("Grid cell bounds not available, using constant grid cell area weights!");
                    }
                }
            }
          else if (coord == 'z')
            {
              int zaxisID = parse_arg.coords[i].cdiID;
              size_t nlev = parse_arg.coords[i].size;
              double *data = (double *) Malloc(nlev * sizeof(double));
              parse_arg.coords[i].data = data;
              cdoZaxisInqLevels(zaxisID, data);
            }
          else
            cdoAbort("Computation of coordinate %c not implemented!", coord);
        }
    }

  for (int varID = parse_arg.nvars1; varID < parse_arg.nparams; varID++)
    {
      int coord = params[varID].coord;
      if (coord)
        {
          char *varname = strdup(params[varID].name);
          varname[strlen(varname) - 2] = 0;
          if (coord == 'x' || coord == 'y' || coord == 'a' || coord == 'w')
            {
              int coordID = params_get_coordID(&parse_arg, coord, params[varID].gridID);
              int gridID = parse_arg.coords[coordID].cdiID;
              size_t ngp = parse_arg.coords[coordID].size;
              double *data = parse_arg.coords[coordID].data;
              assert(gridID == params[varID].gridID);
              assert(data != NULL);

              arrayCopy(ngp, data, params[varID].data);
            }
          else if (coord == 'z')
            {
              int coordID = params_get_coordID(&parse_arg, coord, params[varID].zaxisID);
              int zaxisID = parse_arg.coords[coordID].cdiID;
              size_t nlev = parse_arg.coords[coordID].size;
              double *data = parse_arg.coords[coordID].data;
              assert(zaxisID == params[varID].zaxisID);
              assert(data != NULL);

              arrayCopy(nlev, data, params[varID].data);
            }
          else
            cdoAbort("Computation of coordinate %c not implemented!", coord);

          free(varname);
        }
    }

  if (cdoVerbose) vlistPrint(vlistID2);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = cdoStreamOpenWrite(1, cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int64_t vdate0 = 0;
  int vtime0 = 0;
  int calendar = taxisInqCalendar(taxisID1);

  int nrecs;
  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      int64_t vdate = taxisInqVdate(taxisID1);
      int vtime = taxisInqVtime(taxisID1);

      double jdelta = 0;

      if (tsID)
        {
          juldate_t juldate0 = juldate_encode(calendar, vdate0, vtime0);
          juldate_t juldate = juldate_encode(calendar, vdate, vtime);
          jdelta = juldate_to_seconds(juldate_sub(juldate, juldate0));
        }

      vdate0 = vdate;
      vtime0 = vtime;

      params[vartsID].data[CTIMESTEP] = tsID + 1;
      params[vartsID].data[CDATE] = vdate;
      params[vartsID].data[CTIME] = vtime;
      params[vartsID].data[CDELTAT] = jdelta;

      int year, mon, day;
      int hour, minute, second;
      cdiDecodeDate(vdate, &year, &mon, &day);
      cdiDecodeTime(vtime, &hour, &minute, &second);

      params[vartsID].data[CDAY] = day;
      params[vartsID].data[CMONTH] = mon;
      params[vartsID].data[CYEAR] = year;
      params[vartsID].data[CSECOND] = second;
      params[vartsID].data[CMINUTE] = minute;
      params[vartsID].data[CHOUR] = hour;

      taxisCopyTimestep(taxisID2, taxisID1);

      pstreamDefTimestep(streamID2, tsID);

      for (int varID = 0; varID < nvars1; varID++)
        if (tsID == 0 || params[varID].steptype != TIME_CONSTANT) params[varID].nmiss = 0;

      for (int recID = 0; recID < nrecs; recID++)
        {
          int varID, levelID;
          pstreamInqRecord(streamID1, &varID, &levelID);
          if (parse_arg.needed[varID])
            {
              size_t offset = params[varID].ngp * levelID;
              double *vardata = params[varID].data + offset;
              size_t nmiss;
              pstreamReadRecord(streamID1, vardata, &nmiss);
              params[varID].nmiss += nmiss;
            }
        }

      for (int varID = 0; varID < nvars2; varID++)
        {
          int pidx = varIDmap[varID];
          if (pidx < nvars1) continue;
          size_t ngp = params[pidx].ngp;
          size_t nlev = params[pidx].nlev;

          params[pidx].nmiss = 0;
          arrayFill(ngp * nlev, params[pidx].data, 0.0);
        }

      yy_scan_string(exprs, scanner);
      yyparse(&parse_arg, scanner);

      for (int varID = 0; varID < nvars2; varID++)
        {
          int pidx = varIDmap[varID];

          if (tsID > 0 && params[pidx].steptype == TIME_CONSTANT) continue;

          double missval = vlistInqVarMissval(vlistID2, varID);

          size_t ngp = params[pidx].ngp;
          int nlev = (int) params[pidx].nlev;
          for (int levelID = 0; levelID < nlev; levelID++)
            {
              size_t offset = ngp * levelID;
              double *vardata = params[pidx].data + offset;
              size_t nmiss = arrayNumMV(ngp, vardata, missval);
              pstreamDefRecord(streamID2, varID, levelID);
              pstreamWriteRecord(streamID2, vardata, nmiss);
            }
        }

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  vlistDestroy(vlistID2);

  yylex_destroy(scanner);

  if (exprs) Free(exprs);

  params_delete(params);

  if (parse_arg.needed) Free(parse_arg.needed);
  if (parse_arg.coords)
    {
      for (int i = 0; i < parse_arg.ncoords; i++)
        {
          if (parse_arg.coords[i].data) Free(parse_arg.coords[i].data);
          if (parse_arg.coords[i].units) Free(parse_arg.coords[i].units);
          if (parse_arg.coords[i].longname) Free(parse_arg.coords[i].longname);
        }
      Free(parse_arg.coords);
    }
  if (varIDmap) Free(varIDmap);

  gridDestroy(pointID);
  if ( zonalID != -1 ) gridDestroy(zonalID);

  cdoFinish();

  return 0;
}

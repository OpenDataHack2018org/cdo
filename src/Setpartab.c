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

      Setpartab  setpartab       Set parameter table
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"
#include "pmlist.h"
#include "convert_units.h"

int stringToParam(const char *paramstr);
void paramToStringLong(int param, char *paramstr, int maxlen);

typedef enum {CODE_NUMBER, PARAMETER_ID, VARIABLE_NAME, STANDARD_NAME} pt_mode_t;

typedef struct
{
  bool convert;
  bool remove;
  // missing value
  bool changemissval;
  double missval_old;
  //
  bool lfactor;
  double factor;
  //
  bool checkvalid;
  double valid_min;
  double valid_max;
  //
  bool check_min_mean_abs;
  double ok_min_mean_abs;
  //
  bool check_max_mean_abs;
  double ok_max_mean_abs;
  // units
  bool changeunits;
  char units_old[CDI_MAX_NAME];
  char units[CDI_MAX_NAME];
  // varname
  char name[CDI_MAX_NAME];
  // converter
  void *ut_converter;
} var_t;


static
void defineVarAttText(int vlistID2, int varID, const char *attname, const char *atttext)
{
  int len = strlen(atttext);
  cdiDefAttTxt(vlistID2, varID, attname, len, atttext);
}

static
void defineVarUnits(var_t *vars, int vlistID2, int varID, const char *units)
{
  char units_old[CDI_MAX_NAME];

  vlistInqVarUnits(vlistID2, varID, units_old);
  size_t len1 = strlen(units_old);
  size_t len2 = strlen(units);

  if ( strcmp(units, units_old) != 0 )
    {
      if ( len1 > 0 && len2 > 0 )
	{
	  vars[varID].changeunits = true;
	  strcpy(vars[varID].units_old, units_old);
	  strcpy(vars[varID].units, units);
	}

      vlistDefVarUnits(vlistID2, varID, units);
      defineVarAttText(vlistID2, varID, "original_units", units_old);
    }
}

static
void read_partab(pt_mode_t ptmode, const char *filename, int nvars, int vlistID2, var_t *vars)
{
  const char *hentry[] = {"Header"};
  const char *ventry[] = {"variable_entry", "parameter"};
  int nventry = (int) sizeof(ventry)/sizeof(ventry[0]);
  int nhentry = (int) sizeof(hentry)/sizeof(hentry[0]);
  char valstr[CDI_MAX_NAME];
  char varname[CDI_MAX_NAME];
  int codenum;

  list_t *pml = cdo_parse_namelist(filename);
  if ( pml == NULL ) return;

  // search for global missing value
  bool lmissval = false;
  double missval;
  list_t *kvl = pml_get_kvl_ventry(pml, nhentry, hentry);
  if ( kvl )
    {
      keyValues_t *kv = kvlist_search(kvl, "missing_value");
      if ( kv && kv->nvalues > 0 )
        {
          lmissval = true;
          missval = parameter2double(kv->values[0]);
        }
    }

  for ( int varID = 0; varID < nvars; varID++ )
    {
      vlistInqVarName(vlistID2, varID, varname);

      strcpy(vars[varID].name, varname);
      if ( lmissval )
        {
          double missval_old = vlistInqVarMissval(vlistID2, varID);
          if ( ! DBL_IS_EQUAL(missval, missval_old) )
            {
              vars[varID].changemissval = true;
              vars[varID].missval_old = missval_old;
              vlistDefVarMissval(vlistID2, varID, missval);
            }
        }

      list_t *kvl = NULL;
      if ( ptmode == CODE_NUMBER )
        {
          codenum = vlistInqVarCode(vlistID2, varID);          
          snprintf(valstr, sizeof(valstr), "%d", codenum);
          kvl = pml_search_kvl_ventry(pml, "code", valstr, nventry, ventry);
          if ( kvl )
            {
              int tableID = vlistInqVarTable(vlistID2, varID);
              int tabnum  = tableInqNum(tableID);
              int levtype = zaxisInqLtype(vlistInqVarZaxis(vlistID2, varID));
              int table = tabnum;
              int ltype = levtype;
              keyValues_t *kv = kvlist_search(kvl, "table");
              if ( kv && kv->nvalues == 1 ) table = parameter2int(kv->values[0]);
              kv = kvlist_search(kvl, "ltype");
              if ( kv && kv->nvalues == 1 ) ltype = parameter2int(kv->values[0]);
              if ( !(tabnum == table && levtype == ltype) ) kvl = NULL;
            }
        }
      else if ( ptmode == PARAMETER_ID )
        {
          char paramstr[32];
          int param   = vlistInqVarParam(vlistID2, varID);
          paramToStringLong(param, paramstr, sizeof(paramstr));
          snprintf(valstr, sizeof(valstr), "%s", paramstr);
          kvl = pml_search_kvl_ventry(pml, "param", valstr, nventry, ventry);
          if ( kvl )
            {
              int levtype = zaxisInqLtype(vlistInqVarZaxis(vlistID2, varID));
              int ltype = levtype;
              keyValues_t *kv = kvlist_search(kvl, "ltype");
              if ( kv && kv->nvalues == 1 ) ltype = parameter2int(kv->values[0]);
              if ( !(levtype == ltype) ) kvl = NULL;
            }  
        }
      else if ( ptmode == VARIABLE_NAME )
        {
          kvl = pml_search_kvl_ventry(pml, "name", varname, nventry, ventry);
        }

      if ( kvl )
        {
          int pnum, ptab, pdum;
          cdiDecodeParam(vlistInqVarParam(vlistID2, varID), &pnum, &ptab, &pdum);
          bool lvalid_min = false, lvalid_max = false;

          for ( listNode_t *kvnode = kvl->head; kvnode; kvnode = kvnode->next )
            {
              keyValues_t *kv = *(keyValues_t **)kvnode->data;
              const char *key = kv->key;
              const char *value = (kv->nvalues == 1) ? kv->values[0] : NULL;
              if ( !value ) continue;
              
              // printf("key=%s  value=%s\n", key, value);

              if      ( STR_IS_EQ(key, "standard_name") ) vlistDefVarStdname(vlistID2, varID, value);
              else if ( STR_IS_EQ(key, "long_name")     ) vlistDefVarLongname(vlistID2, varID, value);
              else if ( STR_IS_EQ(key, "units")         ) defineVarUnits(vars, vlistID2, varID, value);
              else if ( STR_IS_EQ(key, "name")          ) /*vlistDefVarName(vlistID2, varID, parameter2word(value))*/;
              else if ( STR_IS_EQ(key, "out_name")      )
                {
                  vlistDefVarName(vlistID2, varID, parameter2word(value));
                  defineVarAttText(vlistID2, varID, "original_name", vars[varID].name);
                }
              else if ( STR_IS_EQ(key, "param")         ) vlistDefVarParam(vlistID2, varID, stringToParam(parameter2word(value)));
              else if ( STR_IS_EQ(key, "out_param")     ) vlistDefVarParam(vlistID2, varID, stringToParam(parameter2word(value)));
              else if ( STR_IS_EQ(key, "code")          ) vlistDefVarParam(vlistID2, varID, cdiEncodeParam(parameter2int(value), ptab, 255));
              else if ( STR_IS_EQ(key, "out_code")      ) vlistDefVarParam(vlistID2, varID, cdiEncodeParam(parameter2int(value), ptab, 255));
              else if ( STR_IS_EQ(key, "comment")       ) defineVarAttText(vlistID2, varID, "comment", value);
              else if ( STR_IS_EQ(key, "cell_methods")  ) defineVarAttText(vlistID2, varID, "cell_methods", value);
              else if ( STR_IS_EQ(key, "cell_measures") ) defineVarAttText(vlistID2, varID, "cell_measures", value);
              else if ( STR_IS_EQ(key, "delete")        ) vars[varID].remove = parameter2bool(value);
              else if ( STR_IS_EQ(key, "convert")       ) vars[varID].convert = parameter2bool(value);
              else if ( STR_IS_EQ(key, "factor")        )
                {
                  vars[varID].lfactor = true;
                  vars[varID].factor = parameter2double(value);
                  if ( cdoVerbose ) cdoPrint("%s - scale factor %g", varname, vars[varID].factor);
                }
              else if ( STR_IS_EQ(key, "missval") )
                {
                  double missval = parameter2double(value);
                  double missval_old = vlistInqVarMissval(vlistID2, varID);
                  if ( ! DBL_IS_EQUAL(missval, missval_old) )
                    {
                      if ( cdoVerbose ) cdoPrint("%s - change missval from %g to %g", varname, missval_old, missval);
                      vars[varID].changemissval = true;
                      vars[varID].missval_old = missval_old;
                      vlistDefVarMissval(vlistID2, varID, missval);
                    }
                }
              else if ( STR_IS_EQ(key, "valid_min") )
                {
                  lvalid_min = true;
                  vars[varID].valid_min = parameter2double(value);
                }
              else if ( STR_IS_EQ(key, "valid_max") )
                {
                  lvalid_max = true;
                  vars[varID].valid_max = parameter2double(value);
                }
              else if ( STR_IS_EQ(key, "ok_min_mean_abs") )
                {
                  vars[varID].check_min_mean_abs = true;
                  vars[varID].ok_min_mean_abs = parameter2double(value);
                }
              else if ( STR_IS_EQ(key, "ok_max_mean_abs") )
                {
                  vars[varID].check_max_mean_abs = true;
                  vars[varID].ok_max_mean_abs = parameter2double(value);
                }
              else if ( STR_IS_EQ(key, "datatype") || STR_IS_EQ(key, "type") )
                {
                  int datatype = str2datatype(parameter2word(value));
                  if ( datatype != -1 ) vlistDefVarDatatype(vlistID2, varID, datatype);
                }
              else
                {
                  if ( cdoVerbose ) cdoPrint("Attribute %s:%s not supported!", varname,  key);
                }
            }

          if ( lvalid_min && lvalid_max ) vars[varID].checkvalid = true;
        }
      else
        {
          cdoPrint("Variable %s not found!", varname);
        }
    }

  list_destroy(pml);
}

static
void check_data(int vlistID2, int varID2, int varID, var_t *vars, long gridsize, double missval, double *array)
{
  char varname[CDI_MAX_NAME];
  int nvals = 0;
  double amean = 0, aval;
  double amin  =  1.e300;
  double amax  = -1.e300;
  
  for ( long i = 0; i < gridsize; ++i )
    {
      aval = array[i];
      if ( !DBL_IS_EQUAL(aval, missval) )
	{
	  if ( aval < amin ) amin = aval;
	  if ( aval > amax ) amax = aval;
	  amean += aval;
	  nvals++;
	}
    }

  if ( nvals > 0 ) amean /= nvals;

  int n_lower_min = 0;
  int n_greater_max = 0;
  for ( long i = 0; i < gridsize; ++i )
    {
      aval = array[i];
      if ( !DBL_IS_EQUAL(aval, missval) )
	{
	  if ( aval < vars[varID].valid_min ) n_lower_min++;
	  if ( aval > vars[varID].valid_max ) n_greater_max++;
	}
    }

  vlistInqVarName(vlistID2, varID2, varname);

  if ( n_lower_min > 0 )
    cdoWarning("Invalid value(s) detected for variable '%s': %i values were lower than minimum valid value (%.4g).",
	       varname, n_lower_min, vars[varID].valid_min);
  if ( n_greater_max > 0 )
    cdoWarning("Invalid value(s) detected for variable '%s': %i values were greater than maximum valid value (%.4g).",
	       varname, n_greater_max, vars[varID].valid_max);

  amean = fabs(amean);

  if ( vars[varID].check_min_mean_abs )
    {
      if ( amean < .1*vars[varID].ok_min_mean_abs )
	cdoWarning("Invalid Absolute Mean for variable '%s' (%.5g) is lower by more than an order of magnitude than minimum allowed: %.4g",
		 varname, amean, vars[varID].ok_min_mean_abs);

      if ( amean < vars[varID].ok_min_mean_abs)
	cdoWarning("Invalid Absolute Mean for variable '%s' (%.5g) is lower than minimum allowed: %.4g",
		   varname, amean, vars[varID].ok_min_mean_abs);
    }

  if ( vars[varID].check_max_mean_abs )
    {
      if ( amean > 10.*vars[varID].ok_max_mean_abs )
	cdoWarning("Invalid Absolute Mean for variable '%s' (%.5g) is greater by more than an order of magnitude than maximum allowed: %.4g",
		 varname, amean, vars[varID].ok_max_mean_abs);
      
      if ( amean > vars[varID].ok_max_mean_abs )
	cdoWarning("Invalid Absolute Mean for variable '%s' (%.5g) is greater than maximum allowed: %.4g",
		   varname, amean, vars[varID].ok_max_mean_abs);
    }
}


void *Setpartab(void *argument)
{
  int nrecs;
  int varID, levelID;
  int nmiss;
  bool delvars = false;
  int tableID = -1;
  int tableformat = 0;
  double missval;

  cdoInitialize(argument);

  int SETCODETAB = cdoOperatorAdd("setcodetab",  0, 0, "parameter code table name");
  int SETPARTABC = cdoOperatorAdd("setpartabc",  0, 0, "parameter table name");
  int SETPARTABP = cdoOperatorAdd("setpartabp",  0, 0, "parameter table name");
  int SETPARTABN = cdoOperatorAdd("setpartabn",  0, 0, "parameter table name");

  int operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  if ( operatorArgc() < 1 ) cdoAbort("Too few arguments!");

  bool convert_data = false;
  if ( operatorArgc() == 2 )
    {
      if ( strcmp("convert", operatorArgv()[1]) == 0 ) convert_data = true;
      else cdoAbort("Unknown parameter: >%s<", operatorArgv()[1]); 
    }

  if ( operatorArgc() > 2 ) cdoAbort("Too many arguments!");

  pt_mode_t ptmode = CODE_NUMBER;
  if      ( operatorID == SETCODETAB ) ptmode = CODE_NUMBER;
  else if ( operatorID == SETPARTABC ) ptmode = CODE_NUMBER;
  else if ( operatorID == SETPARTABP ) ptmode = PARAMETER_ID;
  else if ( operatorID == SETPARTABN ) ptmode = VARIABLE_NAME;

  if ( ptmode == CODE_NUMBER )
    {
      char *partab = operatorArgv()[0];
      FILE *fp = NULL;
      if ( fileExists(partab) ) fp = fopen(partab, "r");
      if ( fp != NULL )
	{
	  fseek(fp, 0L, SEEK_END);
	  size_t fsize = (size_t) ftell(fp);
	  char *parbuf = (char *) Malloc(fsize+1);
	  fseek(fp, 0L, SEEK_SET);
	  fread(parbuf, fsize, 1, fp);
	  parbuf[fsize] = 0;
	  fseek(fp, 0L, SEEK_SET);

	  if ( atoi(parbuf) == 0 ) tableformat = 1;

	  fclose(fp);
	  Free(parbuf);
	}

      if ( tableformat == 0 ) tableID = defineTable(partab);
    }
  else if (  ptmode == PARAMETER_ID )
    {
      tableformat = 1;
    }
  else if (  ptmode == VARIABLE_NAME )
    {
      tableformat = 1;
    }

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);
  /* vlistPrint(vlistID2);*/

  int nvars = vlistNvars(vlistID2);
  var_t *vars = (var_t *) Malloc(nvars*sizeof(var_t));
  memset(vars, 0, nvars*sizeof(var_t));

  if ( convert_data )
    for ( varID = 0; varID < nvars; ++varID ) vars[varID].convert = true;

  if ( tableformat == 0 )
    {
      /*
        for ( int varID = 0; varID < nvars; varID++ )
	  vlistDefVarTable(vlistID2, varID, tableID);
      */
      char name[CDI_MAX_NAME], longname[CDI_MAX_NAME], units[CDI_MAX_NAME];
      for ( int varID = 0; varID < nvars; varID++ )
        {
          int param = vlistInqVarParam(vlistID2, varID);
          int pdis, pcat, pnum;
          cdiDecodeParam(param, &pnum, &pcat, &pdis);
          if ( pdis == 255 )
            {
              int code = pnum;
              if ( tableInqParName(tableID, code, name) == 0 )
                {
                  vlistDefVarName(vlistID2, varID, name);
                  longname[0] = 0;
                  tableInqParLongname(tableID, code, longname);
                  vlistDefVarLongname(vlistID2, varID, longname);
                  units[0] = 0;
                  tableInqParUnits(tableID, code, units);
                  vlistDefVarUnits(vlistID2, varID, units);
                }
            }
          vlistDefVarTable(vlistID2, varID, tableID);
        }
    }
  else
    {
      const char *filename = operatorArgv()[0];
      read_partab(ptmode, filename, nvars, vlistID2, vars);

      for ( int varID = 0; varID < nvars; ++varID )
	if ( vars[varID].remove )
          {
            delvars = true;
            break;
          }

      if ( delvars )
	{
	  vlistClearFlag(vlistID1);
	  vlistClearFlag(vlistID2);

	  for ( int varID = 0; varID < nvars; varID++ )
	    {
	      int zaxisID = vlistInqVarZaxis(vlistID2, varID);
	      int nlevs   = zaxisInqSize(zaxisID);
	      for ( int levID = 0; levID < nlevs; levID++ )
		{
		  vlistDefFlag(vlistID1, varID, levID, TRUE);
		  vlistDefFlag(vlistID2, varID, levID, TRUE);
		  if ( vars[varID].remove )
		    {
		      vlistDefFlag(vlistID1, varID, levID, FALSE);
		      vlistDefFlag(vlistID2, varID, levID, FALSE);
		    }
		}
	    }

	  int vlistIDx = vlistCreate();
	  vlistCopyFlag(vlistIDx, vlistID2);

	  vlistDestroy(vlistID2);

	  vlistID2 = vlistIDx;
          if ( vlistNvars(vlistID2) == 0 ) cdoAbort("No variable selected!");
	}

      for ( int varID = 0; varID < nvars; ++varID )
        {
          var_t *var = &vars[varID];
          if ( var->convert == false ) var->changeunits = false;
          if ( var->changeunits )
            cdoConvertUnits(&var->ut_converter, &var->changeunits, (char*)&var->units, (char*)&var->units_old, var->name);
        }
    }

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  /* vlistPrint(vlistID2);*/
  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  long gridsize = vlistGridsizeMax(vlistID1);
  if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
  double *array = (double *) Malloc(gridsize*sizeof(double));

  int tsID1 = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID1)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID1);
	       
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  int varID2 = varID;
	  int levelID2 = levelID;

	  if ( delvars )
	    {
	      if ( vars[varID].remove ) continue;

	      if ( vlistInqFlag(vlistID1, varID, levelID) == TRUE )
		{
		  varID2   = vlistFindVar(vlistID2, varID);
		  levelID2 = vlistFindLevel(vlistID2, varID, levelID);
		}
	    }

	  streamDefRecord(streamID2,  varID2,  levelID2);

	  streamReadRecord(streamID1, array, &nmiss);

	  missval = vlistInqVarMissval(vlistID2, varID2);
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID2));
	  if ( vlistInqVarNumber(vlistID2, varID2) != CDI_REAL ) gridsize *= 2;

	  if ( nmiss > 0 && vars[varID].changemissval )
	    {
	      for ( long i = 0; i < gridsize; ++i )
		{
		  if ( DBL_IS_EQUAL(array[i], vars[varID].missval_old) ) array[i] = missval;
		}
	    }

	  if ( vars[varID].lfactor )
	    {
	      for ( long i = 0; i < gridsize; ++i )
		{
		  if ( !DBL_IS_EQUAL(array[i], missval) ) array[i] *= vars[varID].factor;
		}
	    }

#if defined(HAVE_UDUNITS2)
	  if ( vars[varID].changeunits )
	    {
	      int nerr = 0;
	      for ( long i = 0; i < gridsize; ++i )
		{
		  if ( !DBL_IS_EQUAL(array[i], missval) )
		    {
		      array[i] = cv_convert_double((const cv_converter*)vars[varID].ut_converter, array[i]);
		      if ( ut_get_status() != UT_SUCCESS ) nerr++;
		    }
		}
	      if ( nerr )
		{
		  cdoWarning("Udunits: Error converting units from [%s] to [%s], parameter: %s",
			     vars[varID].units_old, vars[varID].units, vars[varID].name);
		  vars[varID].changeunits = false;
		}
	    }
#endif
	  
	  streamWriteRecord(streamID2, array, nmiss);

	  if ( vars[varID].checkvalid || vars[varID].check_min_mean_abs || vars[varID].check_max_mean_abs )
	    check_data(vlistID2, varID2, varID, vars, gridsize, missval, array);
	}
      tsID1++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

#if defined(HAVE_UDUNITS2)
  for ( int varID = 0; varID < nvars; varID++ )
    if ( vars[varID].ut_converter ) cdoConvertFree(vars[varID].ut_converter);

  cdoConvertDestroy();
#endif

  if ( array ) Free(array);
  if ( vars  ) Free(vars);

  cdoFinish();

  return 0;
}

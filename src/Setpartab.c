/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2012 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

#if  defined  (HAVE_CONFIG_H)
#  include "config.h"
#endif

#if defined (HAVE_LIBUDUNITS2)
#  include <udunits2/udunits2.h>
#endif

#include <errno.h>
#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"
#include "namelist.h"


#if defined (HAVE_LIBUDUNITS2)
ut_system *ut_read = NULL;

static
void *get_converter(char *src_unit_str, char *tgt_unit_str, int *rstatus)
{
  ut_unit *src_unit, *tgt_unit;
  cv_converter *ut_units_converter = NULL;
  int status;

  *rstatus = -1;

  if ( ut_read == NULL )
    {
      errno = 0;
      ut_read = ut_read_xml(NULL);
      status = ut_get_status();
      if ( status == UT_PARSE )
	{
	  if ( cdoVerbose ) cdoWarning("Udunits: Couldn't parse unit database!");
	}
      if ( status == UT_OPEN_ENV || status == UT_OPEN_DEFAULT || status == UT_OS )
	{
	  if ( cdoVerbose ) cdoWarning("Udunits: %s", strerror(errno));
	}
      errno = 0;
      if ( status != UT_SUCCESS )
	{
	  if ( cdoVerbose ) cdoWarning("Udunits: Error reading units system!");
	  return NULL;
	}
    }

  ut_trim(src_unit_str, UT_ASCII);
  src_unit = ut_parse(ut_read, src_unit_str, UT_ASCII);
  if ( ut_get_status() != UT_SUCCESS )
    {
      if ( cdoVerbose ) cdoWarning("Udunits: Error parsing units: [%s]", src_unit_str);
      return NULL;
    }

  ut_trim(tgt_unit_str, UT_ASCII);
  tgt_unit = ut_parse(ut_read, tgt_unit_str, UT_ASCII);
  if ( ut_get_status() != UT_SUCCESS )
    {
      if ( cdoVerbose ) cdoWarning("Udunits: Error parsing units: [%s]", tgt_unit_str);
      return NULL;
    }

  status = ut_compare(src_unit, tgt_unit);
  if ( status == 0 ) *rstatus = -2;

  if ( *rstatus == -1 )
    {
      status = ut_are_convertible(src_unit, tgt_unit);
      if ( status == 0 ) *rstatus = -3;
    }

  if ( *rstatus == -1 )
    {
      ut_units_converter = ut_get_converter(src_unit, tgt_unit);
      if ( ut_units_converter == NULL || ut_get_status() != UT_SUCCESS )
	{
	  if ( cdoVerbose ) cdoWarning("Udunits: Error getting converter from [%s] to [%s]", src_unit_str, tgt_unit_str);
	}
      else
	*rstatus = 0;
    }

  ut_free(src_unit);
  if ( ut_get_status() != UT_SUCCESS )
    {
      if ( cdoVerbose ) cdoWarning("Udunits: Error freeing units [%s]", src_unit_str);
      return NULL;
    }
     
  ut_free(tgt_unit);
  if ( ut_get_status() != UT_SUCCESS )
    {
      if ( cdoVerbose ) cdoWarning("Udunits: Error freeing units [%s]", tgt_unit_str);
      return NULL;
    }

  return ((void *) ut_units_converter);
}
#endif

typedef struct
{
  // missing value
  int changemissval;
  double missval_old;
  double missval;
  // units
  int changeunits;
  char units_old[CDI_MAX_NAME];
  char units[CDI_MAX_NAME];
  // varname
  char name[CDI_MAX_NAME];
  // converter
  void *ut_converter;
} var_t;

int lwarn_udunits = TRUE;

static
void defineVarUnits(var_t *vars, int vlistID2, int varID, char *units, char *name)
{
  char units_old[CDI_MAX_NAME];
  size_t len1, len2;
  vlistInqVarUnits(vlistID2, varID, units_old);
  len1 = strlen(units_old);
  len2 = strlen(units);

  if ( memcmp(units, units_old, len2) != 0 )
    {
      if ( len1 > 0 && len2 > 0 )
	{
	  int status;
	  vars[varID].changeunits = TRUE;
	  strcpy(vars[varID].units_old, units_old);
	  strcpy(vars[varID].units, units);
#if defined (HAVE_LIBUDUNITS2)
	  vars[varID].ut_converter = get_converter(units_old, units, &status);
	  if ( vars[varID].ut_converter == NULL )
	    {
	      if ( status == -2 )
		{
		  if ( cdoVerbose )
		    cdoPrint("%s - not converted from [%s] to [%s], units are equal!", name, units_old, units);
		}
	      else if ( status == -3 )
		{
		  cdoWarning("%s - converting units from [%s] to [%s] failed, not convertible!", name, units_old, units);
		}
	      else
		cdoWarning("%s - converting units from [%s] to [%s] failed!", name, units_old, units);
	      vars[varID].changeunits = FALSE;
	    }
	  else
	    {
	      if ( cdoVerbose )
		{
		  char buf[64];
		  cv_get_expression(vars[varID].ut_converter, buf, 64, name);
		  cdoPrint("%s - convert units from [%s] to [%s] (expression: %s).", name, units_old, units, buf);
		}
	    }
#else
	  if ( lwarn_udunits )
	    {
	      cdoWarning("Can not convert units, UDUNITS2 support not compiled in!");
	      vars[varID].changeunits = FALSE;
	      lwarn_udunits = FALSE;
	    }
#endif
	}
      vlistDefVarUnits(vlistID2, varID, units);
    }
}


void *Setpartab(void *argument)
{
  int SETPARTAB, SETPARTABV;
  int operatorID;
  int streamID1, streamID2 = CDI_UNDEFID;
  int nrecs, nvars, newval = -1, tabnum = 0;
  int tsID1, recID, varID, levelID;
  int vlistID1, vlistID2;
  int taxisID1, taxisID2;
  int nmiss;
  long gridsize;
  int index, zaxisID1, zaxisID2, nzaxis, nlevs;
  int tableID = -1;
  int tableformat = 0;
  int zaxistype;
  int newparam = 0;
  char *newname = NULL, *partab = NULL;
  double missval;
  double newlevel = 0;
  double *levels = NULL;
  double *array = NULL;
  var_t *vars = NULL;


  cdoInitialize(argument);

  SETPARTAB  = cdoOperatorAdd("setpartab",  0, 0, "parameter table");
  SETPARTABV = cdoOperatorAdd("setpartabv", 0, 0, "parameter table");

  operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));
  if ( operatorID == SETPARTAB )
    {
      FILE *fp;
      size_t fsize;
      char *parbuf = NULL;
      size_t nbytes;

      partab = operatorArgv()[0];
      fp = fopen(partab, "r");
      if ( fp != NULL )
	{
	  fseek(fp, 0L, SEEK_END);
	  fsize = (size_t) ftell(fp);
	  parbuf = (char *) malloc(fsize+1);
	  fseek(fp, 0L, SEEK_SET);
	  nbytes = fread(parbuf, fsize, 1, fp);
	  parbuf[fsize] = 0;
	  fseek(fp, 0L, SEEK_SET);

	  if ( atoi(parbuf) == 0 ) tableformat = 1;

	  fclose(fp);
	  free(parbuf);
	}

      if ( tableformat == 0 ) tableID = defineTable(partab);
    }
  else if ( operatorID == SETPARTABV )
    {
      tableformat = 1;
    }

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);
  /* vlistPrint(vlistID2);*/

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  nvars = vlistNvars(vlistID2);
  vars = (var_t *) malloc(nvars*sizeof(var_t));
  memset(vars, 0, nvars*sizeof(var_t));

  if ( operatorID == SETPARTAB || operatorID == SETPARTABV )
    {
      if ( tableformat == 0 )
	{
	  for ( varID = 0; varID < nvars; varID++ )
	    vlistDefVarTable(vlistID2, varID, tableID);
	}
      else
	{
	  FILE *fp;
	  namelist_t *nml;
	  int nml_code, nml_out_code, nml_table, nml_datatype, nml_name, nml_out_name, nml_stdname;
	  int nml_longname, nml_units, nml_ltype, nml_missval;
	  int locc, i;
	  int code, out_code, table, ltype;
	  int nml_index = 0;
	  int codenum, tabnum, levtype;
	  double missval;
	  char *datatypestr = NULL;
	  char *name = NULL, *out_name = NULL, *stdname = NULL, longname[CDI_MAX_NAME] = "", units[CDI_MAX_NAME] = "";
	  char varname[CDI_MAX_NAME];

	  partab = operatorArgv()[0];
	  fp = fopen(partab, "r");
	  if ( fp == NULL ) cdoAbort("Open failed on parameter table %s!", partab);

	  nml = namelistNew("parameter");
	  nml->dis = 0;

	  nml_code     = namelistAdd(nml, "code",            NML_INT,  0, &code, 1);
	  nml_out_code = namelistAdd(nml, "out_code",        NML_INT,  0, &out_code, 1);
	  nml_table    = namelistAdd(nml, "table",           NML_INT,  0, &table, 1);
	  nml_ltype    = namelistAdd(nml, "ltype",           NML_INT,  0, &ltype, 1);
	  nml_missval  = namelistAdd(nml, "missing_value",   NML_FLT,  0, &missval, 1);
	  nml_datatype = namelistAdd(nml, "datatype",        NML_WORD, 0, &datatypestr, 1);
	  nml_name     = namelistAdd(nml, "name",            NML_WORD, 0, &name, 1);
	  nml_out_name = namelistAdd(nml, "out_name",        NML_WORD, 0, &out_name, 1);
	  nml_stdname  = namelistAdd(nml, "standard_name",   NML_WORD, 0, &stdname, 1);
	  nml_longname = namelistAdd(nml, "long_name",       NML_TEXT, 0, longname, sizeof(longname));
	  nml_units    = namelistAdd(nml, "units",           NML_TEXT, 0, units, sizeof(units));
	      
	  while ( ! feof(fp) )
	    {
	      namelistReset(nml);

	      namelistRead(fp, nml);

	      locc = FALSE;
	      for ( i = 0; i < nml->size; i++ )
		{
		  if ( nml->entry[i]->occ ) { locc = TRUE; break; }
		}

	      if ( locc )
		{
		  /* namelistPrint(nml); */

		  nml_index++;

		  if ( operatorID == SETPARTAB )
		    {
		      if ( nml->entry[nml_code]->occ == 0 )
			{
			  cdoPrint("Parameter %d skipped, code number not found!", nml_index);
			  continue;
			}
		    }
		  else
		    {
		      if ( nml->entry[nml_name]->occ == 0 )
			{
			  cdoWarning("Parameter %d skipped, variable name not found!", nml_index);
			  continue;
			}
		    }

		  for ( varID = 0; varID < nvars; varID++ )
		    {
		      if ( operatorID == SETPARTAB )
			{
			  codenum = vlistInqVarCode(vlistID2, varID);
			  tableID = vlistInqVarTable(vlistID2, varID);
			  tabnum  = tableInqNum(tableID);
			  levtype = zaxisInqLtype(vlistInqVarZaxis(vlistID2, varID));
			  /*
			  printf("code = %d  tabnum = %d  ltype = %d\n", codenum, tabnum, levtype);
			  */
			  if ( nml->entry[nml_table]->occ == 0 ) table = tabnum;
			  if ( nml->entry[nml_ltype]->occ == 0 ) ltype = levtype;

			  if ( codenum == code && tabnum == table && levtype == ltype ) break;
			}
		      else
			{
			  vlistInqVarName(vlistID2, varID, varname);
			  if ( strcmp(varname, name) == 0 ) break;
			}
		    }

		  if ( varID < nvars )
		    {
		      if ( nml->entry[nml_code]->occ     ) vlistDefVarCode(vlistID2, varID, code);
		      if ( nml->entry[nml_out_code]->occ ) vlistDefVarCode(vlistID2, varID, out_code);
		      if ( nml->entry[nml_name]->occ     ) strcpy(vars[varID].name, name);
		      if ( nml->entry[nml_name]->occ     ) vlistDefVarName(vlistID2, varID, name);
		      if ( nml->entry[nml_out_name]->occ ) vlistDefVarName(vlistID2, varID, out_name);
		      if ( nml->entry[nml_stdname]->occ  ) vlistDefVarStdname(vlistID2, varID, stdname);
		      if ( nml->entry[nml_longname]->occ ) vlistDefVarLongname(vlistID2, varID, longname);
		      if ( nml->entry[nml_units]->occ    ) defineVarUnits(vars, vlistID2, varID, units, name);
		      if ( nml->entry[nml_datatype]->occ )
			{
			  int datatype = str2datatype(datatypestr);
			  if ( datatype != -1 ) vlistDefVarDatatype(vlistID2, varID, datatype);
			}
		      if ( nml->entry[nml_missval]->occ )
			{
			  double missval_old;
			  missval_old = vlistInqVarMissval(vlistID2, varID);
			  if ( ! DBL_IS_EQUAL(missval, missval_old) )
			    {
			      if ( cdoVerbose ) 
				cdoPrint("%s - change missval from %g to %g", name, missval_old, missval);
			      vars[varID].changemissval = TRUE;
			      vars[varID].missval_old = missval_old;
			      vars[varID].missval = missval;
			      vlistDefVarMissval(vlistID2, varID, missval);
			    }
			}
		    }
		  else
		    {
		      if ( cdoVerbose )
			{
			  if ( operatorID == SETPARTAB )
			    {
			      if ( nml->entry[nml_table]->occ == 0 )
				cdoPrint("Code %d not found!", code);
			      else
				cdoPrint("Code %d and table %d not found!", code, table);
			    }
			  else
			    cdoPrint("%s - not found!", name);
			}
		    }
		}
	      else
		break;
	    }
	  
	  namelistDelete(nml);

	  fclose(fp);
	}
    }

  /* vlistPrint(vlistID2);*/
  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);
  if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
  array = (double *) malloc(gridsize*sizeof(double));

  tsID1 = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID1)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID1);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamDefRecord(streamID2,  varID,  levelID);

	  streamReadRecord(streamID1, array, &nmiss);

	  missval = vlistInqVarMissval(vlistID2, varID);
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  if ( vlistInqVarNumber(vlistID2, varID) != CDI_REAL ) gridsize *= 2;

	  if ( nmiss > 0 && vars[varID].changemissval == TRUE )
	    {
	      for ( long i = 0; i < gridsize; ++i )
		{
		  if ( DBL_IS_EQUAL(array[i], vars[varID].missval_old) ) array[i] = vars[varID].missval;
		}
	    }

#if defined (HAVE_LIBUDUNITS2)
	  if ( vars[varID].changeunits == TRUE )
	    {
	      int nerr = 0;
	      for ( long i = 0; i < gridsize; ++i )
		{
		  if ( !DBL_IS_EQUAL(array[i], missval) )
		    {
		      array[i] = cv_convert_double(vars[varID].ut_converter, array[i]);
		      if ( ut_get_status() != UT_SUCCESS ) nerr++;
		    }
		}
	      if ( nerr )
		{
		  cdoWarning("Udunits: Error converting units from [%s] to [%s], parameter: %s",
			     vars[varID].units_old, vars[varID].units, vars[varID].name);
		  vars[varID].changeunits = FALSE;
		}
	    }
#endif
	  
	  streamWriteRecord(streamID2, array, nmiss);
	}
      tsID1++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

#if defined (HAVE_LIBUDUNITS2)
  for ( varID = 0; varID < nvars; varID++ )
    if ( vars[varID].ut_converter ) cv_free(vars[varID].ut_converter);
  if ( ut_read ) ut_free_system(ut_read);
#endif

  if ( array ) free(array);
  if ( vars  ) free(vars);

  cdoFinish();

  return (0);
}

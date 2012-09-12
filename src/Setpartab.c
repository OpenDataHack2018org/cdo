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

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "namelist.h"


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
  typedef struct
  {

    int changemissval;
    double missval_old;
    double missval;
    int changeunits;
    char units_old[CDI_MAX_NAME];
    char units[CDI_MAX_NAME];
  } var_t;
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
	  int nml_code, nml_new_code, nml_table, nml_datatype, nml_name, nml_new_name, nml_stdname;
	  int nml_longname, nml_units, nml_ltype, nml_missval;
	  int locc, i;
	  int code, new_code, table, ltype;
	  int nml_index = 0;
	  int codenum, tabnum, levtype;
	  double missval;
	  char *datatype = NULL;
	  char *name = NULL, *new_name = NULL, *stdname = NULL, longname[CDI_MAX_NAME] = "", units[CDI_MAX_NAME] = "";
	  char varname[CDI_MAX_NAME];

	  partab = operatorArgv()[0];
	  fp = fopen(partab, "r");
	  if ( fp == NULL ) cdoAbort("Internal problem! Parameter table %s not available", partab);

	  nml = namelistNew("parameter");
	  nml->dis = 0;

	  nml_code     = namelistAdd(nml, "code",          NML_INT,  0, &code, 1);
	  nml_new_code = namelistAdd(nml, "new_code",      NML_INT,  0, &new_code, 1);
	  nml_table    = namelistAdd(nml, "table",         NML_INT,  0, &table, 1);
	  nml_ltype    = namelistAdd(nml, "ltype",         NML_INT,  0, &ltype, 1);
	  nml_missval  = namelistAdd(nml, "missval",       NML_FLT,  0, &missval, 1);
	  nml_datatype = namelistAdd(nml, "datatype",      NML_WORD, 0, &datatype, 1);
	  nml_name     = namelistAdd(nml, "name",          NML_WORD, 0, &name, 1);
	  nml_new_name = namelistAdd(nml, "new_name",      NML_WORD, 0, &new_name, 1);
	  nml_stdname  = namelistAdd(nml, "standard_name", NML_WORD, 0, &stdname, 1);
	  nml_longname = namelistAdd(nml, "long_name",     NML_TEXT, 0, longname, sizeof(longname));
	  nml_units    = namelistAdd(nml, "units",         NML_TEXT, 0, units, sizeof(units));
	      
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
		      if ( nml->entry[nml_new_code]->occ ) vlistDefVarCode(vlistID2, varID, new_code);
		      if ( nml->entry[nml_name]->occ     ) vlistDefVarName(vlistID2, varID, name);
		      if ( nml->entry[nml_new_name]->occ ) vlistDefVarName(vlistID2, varID, new_name);
		      if ( nml->entry[nml_stdname]->occ  ) vlistDefVarStdname(vlistID2, varID, stdname);
		      if ( nml->entry[nml_longname]->occ ) vlistDefVarLongname(vlistID2, varID, longname);
		      if ( nml->entry[nml_units]->occ    )
			{
			  char units_old[CDI_MAX_NAME];
			  size_t len1, len2;
			  vlistInqVarUnits(vlistID2, varID, units_old);
			  len1 = strlen(units_old);
			  len2 = strlen(units);

			  if ( memcmp(units, units_old, len2) != 0 )
			    {
			      if ( cdoVerbose ) cdoPrint("%s - change units from [%s] to [%s]", name, units_old, units);
			      if ( len1 > 0 && len2 > 0 )
				{
				  vars[varID].changeunits = TRUE;
				  strcpy(vars[varID].units_old, units_old);
				  strcpy(vars[varID].units, units);
				}
			      vlistDefVarUnits(vlistID2, varID, units);
			    }
			}
		      if ( nml->entry[nml_datatype]->occ )
			{
			  int dtype = -1;
			  if ( strlen(datatype) == 3 )
			    {
			      if      ( memcmp(datatype, "F32", 3) == 0 ) dtype = DATATYPE_FLT32;
			      else if ( memcmp(datatype, "F64", 3) == 0 ) dtype = DATATYPE_FLT64;
			    }
			  if ( dtype != -1 ) vlistDefVarDatatype(vlistID2, varID, dtype);
			}
		      if ( nml->entry[nml_missval]->occ )
			{
			  double missval_old;
			  missval_old = vlistInqVarMissval(vlistID2, varID);
			  if ( ! DBL_IS_EQUAL(missval, missval_old) )
			    {
			      if ( cdoVerbose ) cdoPrint("%s - change missval from %g to %g", name, missval_old, missval);
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

	  if ( vars[varID].changeunits == TRUE )
	    {
	      for ( long i = 0; i < gridsize; ++i )
		{
		  //  if ( !DBL_IS_EQUAL(array[i], missval) ) array[i] = cv_convert_double(ut_cdo_converter, array[i])v;
		}	      
	    }
	  
	  streamWriteRecord(streamID2, array, nmiss);
	}
      tsID1++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( array ) free(array);
  if ( vars  ) free(vars);

  cdoFinish();

  return (0);
}

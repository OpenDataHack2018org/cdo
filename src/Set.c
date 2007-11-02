/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2007 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Set        setpartab       Set parameter table
      Set        setcode         Set code number
      Set        setname         Set variable name
      Set        setlevel        Set level
      Set        setltype        Set GRIB level type
*/


#include <string.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "namelist.h"


void *Set(void *argument)
{
  static char func[] = "Set";
  int SETPARTAB, SETPARTABV, SETCODE, SETNAME, SETLEVEL, SETLTYPE;
  int operatorID;
  int streamID1, streamID2 = CDI_UNDEFID;
  int nrecs, nvars, newval = -1;
  int tsID1, recID, varID, levelID;
  int vlistID1, vlistID2;
  int taxisID1, taxisID2;
  int nmiss;
  int gridsize;
  int index, zaxisID1, zaxisID2, nzaxis, nlevs;
  int tableID = -1;
  int tableformat = 0;
  int zaxistype;
  char *newname = NULL, *partab = NULL;
  double newlevel = 0;
  double *levels = NULL;
  double *array = NULL;

  cdoInitialize(argument);

  SETPARTAB  = cdoOperatorAdd("setpartab",  0, 0, "parameter table");
  SETPARTABV = cdoOperatorAdd("setpartabv", 0, 0, "parameter table");
  SETCODE    = cdoOperatorAdd("setcode",    0, 0, "code number");
  SETNAME    = cdoOperatorAdd("setname",    0, 0, "variable name");
  SETLEVEL   = cdoOperatorAdd("setlevel",   0, 0, "level");
  SETLTYPE   = cdoOperatorAdd("setltype",   0, 0, "GRIB level type");

  operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));
  if ( operatorID == SETCODE || operatorID == SETLTYPE )
    {
      newval = atoi(operatorArgv()[0]);
    }
  else if ( operatorID == SETNAME )
    {
      newname = operatorArgv()[0];
    }
  else if ( operatorID == SETPARTAB )
    {
      FILE *fp;
      size_t fsize;
      char *parbuf = NULL;

      partab = operatorArgv()[0];
      fp = fopen(partab, "r");
      if ( fp != NULL )
	{
	  fseek(fp, 0L, SEEK_END);
	  fsize = (size_t) ftell(fp);
	  parbuf = (char *) malloc(fsize+1);
	  fseek(fp, 0L, SEEK_SET);
	  fread(parbuf, fsize, 1, fp);
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
  else if ( operatorID == SETLEVEL )
    {
      newlevel = atof(operatorArgv()[0]);
    }

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);
  /* vlistPrint(vlistID2);*/

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if ( operatorID == SETCODE )
    {
      nvars = vlistNvars(vlistID2);
      for ( varID = 0; varID < nvars; varID++ )
	vlistDefVarCode(vlistID2, varID, newval);
    }
  else if ( operatorID == SETNAME )
    {
      vlistDefVarName(vlistID2, 0, newname);
    }
  else if ( operatorID == SETPARTAB || operatorID == SETPARTABV )
    {
      nvars = vlistNvars(vlistID2);

      if ( tableformat == 0 )
	{
	  for ( varID = 0; varID < nvars; varID++ )
	    vlistDefVarTable(vlistID2, varID, tableID);
	}
      else
	{
	  FILE *fp;
	  NAMELIST *nml;
	  int nml_code, nml_new_code, nml_table, nml_datatype, nml_name, nml_new_name, nml_stdname;
	  int nml_longname, nml_units;
	  int locc, i;
	  int code, new_code, table;
	  int nml_index = 0;
	  int tabnum;
	  char *datatype = NULL;
	  char *name = NULL, *new_name = NULL, *stdname = NULL, longname[256] = "", units[256] = "";
	  char varname[256];

	  partab = operatorArgv()[0];
	  fp = fopen(partab, "r");
	  if ( fp == NULL ) cdoAbort("Internal problem! Parameter table %s not available", partab);

	  nml = namelistNew("parameter");
	  nml->dis = 0;

	  nml_code     = namelistAdd(nml, "code",          NML_INT,  0, &code, 1);
	  nml_new_code = namelistAdd(nml, "new_code",      NML_INT,  0, &new_code, 1);
	  nml_table    = namelistAdd(nml, "table",         NML_INT,  0, &table, 1);
	  nml_datatype = namelistAdd(nml, "datatype",      NML_WORD, 0, &datatype, 1);
	  nml_name     = namelistAdd(nml, "name",          NML_WORD, 0, &name, 1);
	  nml_new_name = namelistAdd(nml, "new_name",      NML_WORD, 0, &new_name, 1);
	  nml_stdname  = namelistAdd(nml, "standard_name", NML_WORD, 0, &stdname, 1);
	  nml_longname = namelistAdd(nml, "long_name",     NML_TEXT, 0, longname, sizeof(longname));
	  nml_units    = namelistAdd(nml, "units",         NML_TEXT, 0, units, sizeof(units));
	      
	  while ( ! feof(fp) )
	    {
	      namelistClear(nml);

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

		      if ( nml->entry[nml_table]->occ == 0 )
			for ( varID = 0; varID < nvars; varID++ )
			  {
			    if ( vlistInqVarCode(vlistID2, varID) == code ) break;
			  }
		      else
			for ( varID = 0; varID < nvars; varID++ )
			  {
			    tableID = vlistInqVarTable(vlistID2, varID);
			    tabnum  = tableInqNum(tableID);
			    if ( vlistInqVarCode(vlistID2, varID) == code &&
				 tabnum == table ) break;
			  }
		    }
		  else
		    {
		      if ( nml->entry[nml_name]->occ == 0 )
			{
			  cdoWarning("Parameter %d skipped, variable name not found!", nml_index);
			  continue;
			}

		      for ( varID = 0; varID < nvars; varID++ )
			{
			  vlistInqVarName(vlistID2, varID, varname);
			  if ( strcmp(varname, name) == 0 ) break;
			}
		    }

		  if ( varID < nvars )
		    {
		      if ( nml->entry[nml_code]->occ )     vlistDefVarCode(vlistID2, varID, code);
		      if ( nml->entry[nml_new_code]->occ ) vlistDefVarCode(vlistID2, varID, new_code);
		      if ( nml->entry[nml_name]->occ )     vlistDefVarName(vlistID2, varID, name);
		      if ( nml->entry[nml_new_name]->occ ) vlistDefVarName(vlistID2, varID, new_name);
		      if ( nml->entry[nml_stdname]->occ )  vlistDefVarStdname(vlistID2, varID, stdname);
		      if ( nml->entry[nml_longname]->occ ) vlistDefVarLongname(vlistID2, varID, longname);
		      if ( nml->entry[nml_units]->occ )    vlistDefVarUnits(vlistID2, varID, units);
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
			    cdoPrint("Variable %s not found!", name);
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
  else if ( operatorID == SETLEVEL )
    {
      nzaxis = vlistNzaxis(vlistID2);
      for ( index = 0; index < nzaxis; index++ )
	{
	  zaxisID1 = vlistZaxis(vlistID2, index);
	  zaxisID2 = zaxisDuplicate(zaxisID1);
	  nlevs = zaxisInqSize(zaxisID2);
	  levels = (double *) malloc(nlevs*sizeof(double));
	  zaxisInqLevels(zaxisID2, levels);
	  levels[0] = newlevel;
	  zaxisDefLevels(zaxisID2, levels);
	  /* printf("lev %d %g\n", nlevs, levels[0]); */
	  vlistChangeZaxis(vlistID2, index, zaxisID2);
	  free(levels);
	}
    }
  else if ( operatorID == SETLTYPE )
    {
      nzaxis = vlistNzaxis(vlistID2);
      for ( index = 0; index < nzaxis; index++ )
	{
	  zaxisID1 = vlistZaxis(vlistID2, index);
	  zaxisID2 = zaxisDuplicate(zaxisID1);

	  zaxistype = ltype2ztype(newval);

	  zaxisChangeType(zaxisID2, zaxistype);
	  if ( zaxistype == ZAXIS_GENERIC ) zaxisDefLtype(zaxisID2, newval);
	  vlistChangeZaxis(vlistID2, index, zaxisID2);
	}
    }

  /* vlistPrint(vlistID2);*/
  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID2);
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
	  streamWriteRecord(streamID2, array, nmiss);
	}
      tsID1++;
    }

  streamClose(streamID1);
  streamClose(streamID2);

  if ( array ) free(array);

  cdoFinish();

  return (0);
}

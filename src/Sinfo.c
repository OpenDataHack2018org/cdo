/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2014 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Sinfo      sinfo           Short dataset information
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"

#include "printinfo.h"

const char * tunit2str(int tunits)
{
  if      ( tunits == TUNIT_YEAR )    return ("years");
  else if ( tunits == TUNIT_MONTH )   return ("months");
  else if ( tunits == TUNIT_DAY )     return ("days");
  else if ( tunits == TUNIT_12HOURS ) return ("12hours");
  else if ( tunits == TUNIT_6HOURS )  return ("6hours");
  else if ( tunits == TUNIT_3HOURS )  return ("3hours");
  else if ( tunits == TUNIT_HOUR )    return ("hours");
  else if ( tunits == TUNIT_MINUTE )  return ("minutes");
  else if ( tunits == TUNIT_SECOND )  return ("seconds");
  else                                return ("unknown");
}


const char * calendar2str(int calendar)
{
  if      ( calendar == CALENDAR_STANDARD )  return ("standard");
  else if ( calendar == CALENDAR_PROLEPTIC ) return ("proleptic_gregorian");
  else if ( calendar == CALENDAR_360DAYS )   return ("360_day");
  else if ( calendar == CALENDAR_365DAYS )   return ("365_day");
  else if ( calendar == CALENDAR_366DAYS )   return ("366_day");
  else                                       return ("unknown");
}


void *Sinfo(void *argument)
{
  enum {func_generic, func_param, func_name, func_code};
  int operatorID;
  int operfunc, lensemble;
  int indf;
  int varID;
  int gridsize = 0;
  int gridID, zaxisID, code, tabnum, param;
  int vdate, vtime;
  int nrecs, nvars, ntsteps;
  int levelsize;
  int tsID, ntimeout;
  int tsteptype, taxisID;
  char varname[CDI_MAX_NAME];
  char paramstr[32];
  char vdatestr[32], vtimestr[32];
  char *modelptr, *instptr;
  int streamID = 0;
  int vlistID;
  int datatype;
  char pstr[4];

  cdoInitialize(argument);

  cdoOperatorAdd("sinfo",   func_generic, 0, NULL);
  cdoOperatorAdd("sinfop",  func_param,   0, NULL);
  cdoOperatorAdd("sinfon",  func_name,    0, NULL);
  cdoOperatorAdd("sinfoc",  func_code,    0, NULL);
  cdoOperatorAdd("seinfo",  func_generic, 1, NULL);
  cdoOperatorAdd("seinfop", func_param,   1, NULL);
  cdoOperatorAdd("seinfon", func_name,    1, NULL);
  cdoOperatorAdd("seinfoc", func_code,    1, NULL);

  operatorID = cdoOperatorID();

  operfunc  = cdoOperatorF1(operatorID);
  lensemble = cdoOperatorF2(operatorID);

  for ( indf = 0; indf < cdoStreamCnt(); indf++ )
    {
      streamID = streamOpenRead(cdoStreamName(indf));

      vlistID = streamInqVlist(streamID);

      set_text_color(stdout, BRIGHT, BLACK);
      fprintf(stdout, "   File format");
      reset_text_color(stdout);
      fprintf(stdout, " : ");
      printFiletype(streamID, vlistID);

      set_text_color(stdout, BRIGHT, BLACK);
      if ( lensemble )
	fprintf(stdout, "%6d : Institut Source   Ttype    Einfo Levels Num  Gridsize Num Dtype : ",  -(indf+1));
      else
	fprintf(stdout, "%6d : Institut Source   Ttype    Levels Num  Gridsize Num Dtype : ",  -(indf+1));

      if      ( operfunc == func_name ) fprintf(stdout, "Parameter name");
      else if ( operfunc == func_code ) fprintf(stdout, "Table Code");
      else                              fprintf(stdout, "Parameter ID");

      if ( cdoVerbose ) fprintf(stdout, " : Extra" );              
      reset_text_color(stdout);
      fprintf(stdout, "\n" );              

      nvars = vlistNvars(vlistID);

      for ( varID = 0; varID < nvars; varID++ )
	{
	  param   = vlistInqVarParam(vlistID, varID);
	  code    = vlistInqVarCode(vlistID, varID);
	  tabnum  = tableInqNum(vlistInqVarTable(vlistID, varID));
	  gridID  = vlistInqVarGrid(vlistID, varID);
	  zaxisID = vlistInqVarZaxis(vlistID, varID);

	  set_text_color(stdout, BRIGHT, BLACK);
	  fprintf(stdout, "%6d", varID+1);
	  reset_text_color(stdout);
	  set_text_color(stdout, RESET, BLACK);
	  fprintf(stdout, " : ");
	  reset_text_color(stdout);
	      
	  set_text_color(stdout, RESET, BLUE);
	  /* institute info */
	  instptr = institutInqNamePtr(vlistInqVarInstitut(vlistID, varID));
	  if ( instptr )
	    fprintf(stdout, "%-8s ", instptr);
	  else
	    fprintf(stdout, "unknown  ");

	  /* source info */
	  modelptr = modelInqNamePtr(vlistInqVarModel(vlistID, varID));
	  if ( modelptr )
	    {
	      size_t len = strlen(modelptr);
	      if ( len > 10 )
		for ( size_t i = 3; i < len; ++i )
		  if ( modelptr[i] == ' ' )
		    {
		      modelptr[i] = 0;
		      break;
		    }
	      fprintf(stdout, "%-8s ", modelptr);
	    }
	  else
	    fprintf(stdout, "unknown  ");

	  /* tsteptype */
	  tsteptype = vlistInqVarTsteptype(vlistID, varID);
	  if      ( tsteptype == TSTEP_CONSTANT ) fprintf(stdout, "%-8s ", "constant");
	  else if ( tsteptype == TSTEP_INSTANT  ) fprintf(stdout, "%-8s ", "instant");
	  else if ( tsteptype == TSTEP_INSTANT2 ) fprintf(stdout, "%-8s ", "instant");
	  else if ( tsteptype == TSTEP_INSTANT3 ) fprintf(stdout, "%-8s ", "instant");
	  else if ( tsteptype == TSTEP_MIN      ) fprintf(stdout, "%-8s ", "min");
	  else if ( tsteptype == TSTEP_MAX      ) fprintf(stdout, "%-8s ", "max");
	  else if ( tsteptype == TSTEP_ACCUM    ) fprintf(stdout, "%-8s ", "accum");
	  else                                    fprintf(stdout, "%-8s ", "unknown");

	  /* ensemble information */
	  if ( lensemble )
	    {
	      int ensID, ensCount, forecast_type;
	      if ( vlistInqVarEnsemble(vlistID, varID, &ensID, &ensCount, &forecast_type) )
		fprintf(stdout, "%2d/%-2d ", ensID, ensCount);
	      else
		fprintf(stdout, "--/-- ");
	    }

	  /* layer info */
	  levelsize = zaxisInqSize(zaxisID);
	  fprintf(stdout, "%6d ", levelsize);
	  fprintf(stdout, "%3d ", vlistZaxisIndex(vlistID, zaxisID) + 1);

	  /* grid info */
	  gridsize = gridInqSize(gridID);
	  fprintf(stdout, "%9d ", gridsize);
	  fprintf(stdout, "%3d ", vlistGridIndex(vlistID, gridID) + 1);

	  /* datatype */
	  datatype = vlistInqVarDatatype(vlistID, varID);
	  datatype2str(datatype, pstr);

	  fprintf(stdout, " %-3s", pstr);

	  if ( vlistInqVarCompType(vlistID, varID) == COMPRESS_NONE )
	    fprintf(stdout, "  ");
	  else
	    fprintf(stdout, "z ");

	  reset_text_color(stdout);
	      
	  set_text_color(stdout, RESET, BLACK);
	  fprintf(stdout, ": ");
	  reset_text_color(stdout);

	  /* parameter info */
	  cdiParamToString(param, paramstr, sizeof(paramstr));

	  if ( operfunc == func_name ) vlistInqVarName(vlistID, varID, varname);

	  set_text_color(stdout, BRIGHT, GREEN);
	  if ( operfunc == func_name )
	    fprintf(stdout, "%-14s", varname);
	  else if ( operfunc == func_code )
	    fprintf(stdout, "%4d %4d   ", tabnum, code);
	  else
	    fprintf(stdout, "%-14s", paramstr);
	  reset_text_color(stdout);

	  if ( cdoVerbose )
	    {
	      char varextra[CDI_MAX_NAME];
	      vlistInqVarExtra(vlistID, varID, varextra);
	      fprintf(stdout, " : %s", varextra );              
	    }

	  fprintf(stdout, "\n");
	}

      set_text_color(stdout, BRIGHT, BLACK);
      fprintf(stdout, "   Grid coordinates");
      reset_text_color(stdout);
      fprintf(stdout, " :\n");

      printGridInfo(vlistID);

      set_text_color(stdout, BRIGHT, BLACK);
      fprintf(stdout, "   Vertical coordinates");
      reset_text_color(stdout);
      fprintf(stdout, " :\n");

      printZaxisInfo(vlistID);

      taxisID = vlistInqTaxis(vlistID);
      ntsteps = vlistNtsteps(vlistID);

      if ( ntsteps != 0 )
	{
	  set_text_color(stdout, BRIGHT, BLACK);
	  fprintf(stdout, "   Time coordinate");
	  reset_text_color(stdout);
	  if ( ntsteps == CDI_UNDEFID )
	    fprintf(stdout, " :  unlimited steps\n");
	  else
	    fprintf(stdout, " :  %d step%s\n", ntsteps, ntsteps == 1 ? "" : "s");

	  if ( taxisID != CDI_UNDEFID )
	    {
	      if ( taxisInqType(taxisID) == TAXIS_RELATIVE )
		{
		  int calendar, tunits;

		  vdate = taxisInqRdate(taxisID);
		  vtime = taxisInqRtime(taxisID);

		  date2str(vdate, vdatestr, sizeof(vdatestr));
		  time2str(vtime, vtimestr, sizeof(vtimestr));

		  fprintf(stdout, "     RefTime = %s %s", vdatestr, vtimestr);
		      
		  tunits = taxisInqTunit(taxisID);
		  if ( tunits != CDI_UNDEFID )  fprintf(stdout, "  Units = %s", tunit2str(tunits));
	      
		  calendar = taxisInqCalendar(taxisID);
		  if ( calendar != CDI_UNDEFID )  fprintf(stdout, "  Calendar = %s", calendar2str(calendar));

		  if ( taxisHasBounds(taxisID) )
		    fprintf(stdout, "  Bounds = true");

		  fprintf(stdout, "\n");
		}
	    }

	  fprintf(stdout, "  YYYY-MM-DD hh:mm:ss  YYYY-MM-DD hh:mm:ss  YYYY-MM-DD hh:mm:ss  YYYY-MM-DD hh:mm:ss\n");

	  set_text_color(stdout, RESET, MAGENTA);
	  ntimeout = 0;
	  tsID = 0;
	  while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	    {
	      if ( ntimeout == 4 )
		{
		  ntimeout = 0;
		  fprintf(stdout, "\n");
		}

	      vdate = taxisInqVdate(taxisID);
	      vtime = taxisInqVtime(taxisID);

	      date2str(vdate, vdatestr, sizeof(vdatestr));
	      time2str(vtime, vtimestr, sizeof(vtimestr));

	      fprintf(stdout, " %s %s", vdatestr, vtimestr);

	      ntimeout++;
	      tsID++;
	    }
	  reset_text_color(stdout);
	  fprintf(stdout, "\n");
	}

      streamClose(streamID);
    }

  cdoFinish();

  return (0);
}

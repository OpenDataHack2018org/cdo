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

      Sinfo      sinfo           Short dataset information
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

#include "cdi_uuid.h"
#include "printinfo.h"
#include "text.h"

const char *
tunit2str(int tunits)
{
  // clang-format off
  if      ( tunits == TUNIT_YEAR )       return ("years");
  else if ( tunits == TUNIT_MONTH )      return ("months");
  else if ( tunits == TUNIT_DAY )        return ("days");
  else if ( tunits == TUNIT_12HOURS )    return ("12hours");
  else if ( tunits == TUNIT_6HOURS )     return ("6hours");
  else if ( tunits == TUNIT_3HOURS )     return ("3hours");
  else if ( tunits == TUNIT_HOUR )       return ("hours");
  else if ( tunits == TUNIT_30MINUTES )  return ("30minutes");
  else if ( tunits == TUNIT_QUARTER )    return ("15minutes");
  else if ( tunits == TUNIT_MINUTE )     return ("minutes");
  else if ( tunits == TUNIT_SECOND )     return ("seconds");
  else                                   return ("unknown");
  // clang-format on
}

const char *
calendar2str(int calendar)
{
  // clang-format off
  if      ( calendar == CALENDAR_STANDARD )  return ("standard");
  else if ( calendar == CALENDAR_GREGORIAN ) return ("gregorian");
  else if ( calendar == CALENDAR_PROLEPTIC ) return ("proleptic_gregorian");
  else if ( calendar == CALENDAR_360DAYS )   return ("360_day");
  else if ( calendar == CALENDAR_365DAYS )   return ("365_day");
  else if ( calendar == CALENDAR_366DAYS )   return ("366_day");
  else                                       return ("unknown");
  // clang-format on
}

static void
limitStringLength(char *string, size_t maxlen)
{
  string[maxlen - 1] = 0;
  size_t len = strlen(string);

  if (len > 10)
    {
      for (size_t i = 3; i < len; ++i)
        if (string[i] == ' ' || string[i] == ',' || (i > 10 && string[i] == '.'))
          {
            string[i] = 0;
            break;
          }
    }
}

void *
Sinfo(void *process)
{
  enum
  {
    func_generic,
    func_param,
    func_name,
    func_code
  };
  char tmpname[CDI_MAX_NAME];
  char varname[CDI_MAX_NAME];
  char paramstr[32];
  char vdatestr[32], vtimestr[32];
  char pstr[4];

  cdoInitialize(process);

  // clang-format off
  cdoOperatorAdd("sinfo",   func_generic, 0, NULL);
  cdoOperatorAdd("sinfop",  func_param,   0, NULL);
  cdoOperatorAdd("sinfon",  func_name,    0, NULL);
  cdoOperatorAdd("sinfoc",  func_code,    0, NULL);
  cdoOperatorAdd("seinfo",  func_generic, 1, NULL);
  cdoOperatorAdd("seinfop", func_param,   1, NULL);
  cdoOperatorAdd("seinfon", func_name,    1, NULL);
  cdoOperatorAdd("seinfoc", func_code,    1, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();

  int operfunc = cdoOperatorF1(operatorID);
  int lensemble = cdoOperatorF2(operatorID);

  for (int indf = 0; indf < cdoStreamCnt(); indf++)
    {
      int streamID = cdoStreamOpenRead(cdoStreamName(indf));
      int vlistID = cdoStreamInqVlist(streamID);

      set_text_color(stdout, BRIGHT, BLACK);
      fprintf(stdout, "   File format");
      reset_text_color(stdout);
      fprintf(stdout, " : ");
      printFiletype(streamID, vlistID);

      int nvars = vlistNvars(vlistID);
      int nsubtypes = vlistNsubtypes(vlistID);

      set_text_color(stdout, BRIGHT, BLACK);
      if (lensemble)
        fprintf(stdout, "%6d : Institut Source   T Steptype Einfo Levels Num    Points Num Dtype : ", -(indf + 1));
      else if (nsubtypes > 1)
        fprintf(stdout, "%6d : Institut Source   T Steptype Subtypes Levels Num    Points Num Dtype : ", -(indf + 1));
      else
        fprintf(stdout, "%6d : Institut Source   T Steptype Levels Num    Points Num Dtype : ", -(indf + 1));

      if (operfunc == func_name)
        fprintf(stdout, "Parameter name");
      else if (operfunc == func_code)
        fprintf(stdout, "Table Code");
      else
        fprintf(stdout, "Parameter ID");

      if (cdoVerbose) fprintf(stdout, " : Extra");
      reset_text_color(stdout);
      fprintf(stdout, "\n");

      for (int varID = 0; varID < nvars; varID++)
        {
          int param = vlistInqVarParam(vlistID, varID);
          int code = vlistInqVarCode(vlistID, varID);
          int tabnum = tableInqNum(vlistInqVarTable(vlistID, varID));
          int gridID = vlistInqVarGrid(vlistID, varID);
          int zaxisID = vlistInqVarZaxis(vlistID, varID);

          set_text_color(stdout, BRIGHT, BLACK);
          fprintf(stdout, "%6d", varID + 1);
          reset_text_color(stdout);
          set_text_color(stdout, RESET, BLACK);
          fprintf(stdout, " : ");
          reset_text_color(stdout);

          set_text_color(stdout, RESET, BLUE);
          // institute info
          const char *instptr = institutInqNamePtr(vlistInqVarInstitut(vlistID, varID));
          strcpy(tmpname, "unknown");
          if (instptr) strncpy(tmpname, instptr, CDI_MAX_NAME);
          limitStringLength(tmpname, 32);
          fprintf(stdout, "%-8s ", tmpname);

          // source info
          const char *modelptr = modelInqNamePtr(vlistInqVarModel(vlistID, varID));
          strcpy(tmpname, "unknown");
          if (modelptr) strncpy(tmpname, modelptr, CDI_MAX_NAME);
          limitStringLength(tmpname, 32);
          fprintf(stdout, "%-8s ", tmpname);

          // timetype
          int timetype = vlistInqVarTimetype(vlistID, varID);
          fprintf(stdout, "%c ", timetype == TIME_CONSTANT ? 'c' : 'v');

          // tsteptype
          int tsteptype = vlistInqVarTsteptype(vlistID, varID);
          // clang-format off
	  if      ( tsteptype == TSTEP_INSTANT  ) fprintf(stdout, "%-8s ", "instant");
	  else if ( tsteptype == TSTEP_INSTANT2 ) fprintf(stdout, "%-8s ", "instant");
	  else if ( tsteptype == TSTEP_INSTANT3 ) fprintf(stdout, "%-8s ", "instant");
	  else if ( tsteptype == TSTEP_MIN      ) fprintf(stdout, "%-8s ", "min");
	  else if ( tsteptype == TSTEP_MAX      ) fprintf(stdout, "%-8s ", "max");
	  else if ( tsteptype == TSTEP_AVG      ) fprintf(stdout, "%-8s ", "avg");
	  else if ( tsteptype == TSTEP_ACCUM    ) fprintf(stdout, "%-8s ", "accum");
	  else if ( tsteptype == TSTEP_RANGE    ) fprintf(stdout, "%-8s ", "range");
	  else if ( tsteptype == TSTEP_DIFF     ) fprintf(stdout, "%-8s ", "diff");
	  else                                    fprintf(stdout, "%-8s ", "unknown");
          // clang-format on

          /* ensemble information */
          if (lensemble)
            {
              int perturbationNumber, numberOfForecastsInEnsemble;
              int r1 = cdiInqKeyInt(vlistID, varID, CDI_KEY_PERTURBATIONNUMBER, &perturbationNumber);
              int r2 = cdiInqKeyInt(vlistID, varID, CDI_KEY_NUMBEROFFORECASTSINENSEMBLE, &numberOfForecastsInEnsemble);
              if (r1 == 0 && r2 == 0)
                fprintf(stdout, "%2d/%-2d ", perturbationNumber, numberOfForecastsInEnsemble);
              else
                fprintf(stdout, "--/-- ");
            }

          if (nsubtypes > 1)
            {
              int subtypeID = vlistInqVarSubtype(vlistID, varID);
              int subtypesize = subtypeInqSize(subtypeID);
              fprintf(stdout, " %6d  ", subtypesize);
              fprintf(stdout, "%3d ", vlistSubtypeIndex(vlistID, subtypeID) + 1);
            }
          reset_text_color(stdout);

          /* layer info */
          int levelsize = zaxisInqSize(zaxisID);
          set_text_color(stdout, RESET, GREEN);
          fprintf(stdout, "%6d ", levelsize);
          reset_text_color(stdout);
          fprintf(stdout, "%3d ", vlistZaxisIndex(vlistID, zaxisID) + 1);

          /* grid info */
          size_t gridsize = gridInqSize(gridID);
          set_text_color(stdout, RESET, GREEN);
          fprintf(stdout, "%9zu ", gridsize);
          reset_text_color(stdout);
          fprintf(stdout, "%3d ", vlistGridIndex(vlistID, gridID) + 1);

          /* datatype */
          int datatype = vlistInqVarDatatype(vlistID, varID);
          datatype2str(datatype, pstr);

          set_text_color(stdout, RESET, BLUE);
          fprintf(stdout, " %-3s", pstr);

          int comptype = vlistInqVarCompType(vlistID, varID);
          if (comptype == CDI_COMPRESS_NONE)
            fprintf(stdout, "  ");
          else
            fprintf(stdout, "%c ", (int) comp_name(comptype)[0]);

          reset_text_color(stdout);

          set_text_color(stdout, RESET, BLACK);
          fprintf(stdout, ": ");
          reset_text_color(stdout);

          /* parameter info */
          cdiParamToString(param, paramstr, sizeof(paramstr));

          if (operfunc == func_name) vlistInqVarName(vlistID, varID, varname);

          set_text_color(stdout, BRIGHT, GREEN);
          if (operfunc == func_name)
            fprintf(stdout, "%-14s", varname);
          else if (operfunc == func_code)
            fprintf(stdout, "%4d %4d   ", tabnum, code);
          else
            fprintf(stdout, "%-14s", paramstr);
          reset_text_color(stdout);

          if (cdoVerbose)
            {
              char varextra[CDI_MAX_NAME];
              vlistInqVarExtra(vlistID, varID, varextra);
              fprintf(stdout, " : %s", varextra);
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

      if (nsubtypes > 1)
        {
          fprintf(stdout, "   Subtypes");
          fprintf(stdout, " :\n");

          printSubtypeInfo(vlistID);
        }

      int taxisID = vlistInqTaxis(vlistID);
      int ntsteps = vlistNtsteps(vlistID);

      if (ntsteps != 0)
        {
          set_text_color(stdout, BRIGHT, BLACK);
          fprintf(stdout, "   Time coordinate");
          reset_text_color(stdout);
          if (ntsteps == CDI_UNDEFID)
            fprintf(stdout, " :  unlimited steps\n");
          else
            fprintf(stdout, " :  %d step%s\n", ntsteps, ntsteps == 1 ? "" : "s");

          if (taxisID != CDI_UNDEFID)
            {
              if (taxisInqType(taxisID) != TAXIS_ABSOLUTE)
                {
                  int64_t vdate = taxisInqRdate(taxisID);
                  int vtime = taxisInqRtime(taxisID);

                  date2str(vdate, vdatestr, sizeof(vdatestr));
                  time2str(vtime, vtimestr, sizeof(vtimestr));

                  fprintf(stdout, "     RefTime = %s %s", vdatestr, vtimestr);

                  int tunits = taxisInqTunit(taxisID);
                  if (tunits != CDI_UNDEFID) fprintf(stdout, "  Units = %s", tunit2str(tunits));

                  int calendar = taxisInqCalendar(taxisID);
                  if (calendar != CDI_UNDEFID) fprintf(stdout, "  Calendar = %s", calendar2str(calendar));

                  if (taxisHasBounds(taxisID)) fprintf(stdout, "  Bounds = true");

                  fprintf(stdout, "\n");

                  if (taxisInqType(taxisID) == TAXIS_FORECAST)
                    {
                      vdate = taxisInqFdate(taxisID);
                      vtime = taxisInqFtime(taxisID);

                      date2str(vdate, vdatestr, sizeof(vdatestr));
                      time2str(vtime, vtimestr, sizeof(vtimestr));

                      fprintf(stdout, "     ForecastRefTime = %s %s", vdatestr, vtimestr);
                      fprintf(stdout, "\n");
                    }
                }
            }

          fprintf(stdout, "  YYYY-MM-DD hh:mm:ss  YYYY-MM-DD hh:mm:ss  YYYY-MM-DD hh:mm:ss  YYYY-MM-DD hh:mm:ss\n");

          set_text_color(stdout, RESET, MAGENTA);

          printTimesteps(streamID, taxisID, cdoVerbose);

          reset_text_color(stdout);
          fprintf(stdout, "\n");
        }

      pstreamClose(streamID);
    }

  cdoFinish();

  return 0;
}

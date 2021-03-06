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

      Splittime  splithour       Split hours
      Splittime  splitday        Split days
      Splittime  splitmon        Split months
      Splittime  splitseas       Split seasons
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

#include <time.h>

#define MAX_STREAMS 32

struct tm
datetime_to_tm(int64_t date, int time)
{
  int year, month, day, hour, minute, second;
  cdiDecodeDate(date, &year, &month, &day);
  cdiDecodeTime(time, &hour, &minute, &second);

  struct tm stime;
  memset(&stime, 0, sizeof(struct tm));

  stime.tm_sec = second;
  stime.tm_min = minute;
  stime.tm_hour = hour;
  stime.tm_mday = day;
  stime.tm_mon = month - 1;
  stime.tm_year = year - 1900;

  return stime;
}

void *
Splittime(void *process)
{
  int streamID2;
  int varID;
  int nrecs;
  int levelID;
  int streamIDs[MAX_STREAMS], tsIDs[MAX_STREAMS];
  int index = 0;
  size_t nmiss;
  char filesuffix[32];
  char filename[8192];

  cdoInitialize(process);

  if (processSelf().m_ID != 0) cdoAbort("This operator can't be combined with other operators!");

  bool lcopy = UNCHANGED_RECORD;

  // clang-format off
  int SPLITHOUR = cdoOperatorAdd("splithour", func_time, 10000, NULL);
  int SPLITDAY  = cdoOperatorAdd("splitday",  func_date,     1, NULL);
  int SPLITMON  = cdoOperatorAdd("splitmon",  func_date,   100, NULL);
  int SPLITSEAS = cdoOperatorAdd("splitseas", func_date,   100, NULL);
  // clang-format on

  UNUSED(SPLITDAY);
  UNUSED(SPLITHOUR);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);
  int operintval = cdoOperatorF2(operatorID);

  const char *format = NULL;
  if (operatorID == SPLITMON && operatorArgc() == 1) format = operatorArgv()[0];

  const char *seas_name[4];
  get_season_name(seas_name);

  for (int i = 0; i < MAX_STREAMS; i++) streamIDs[i] = -1;
  for (int i = 0; i < MAX_STREAMS; i++) tsIDs[i] = 0;

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  strcpy(filename, cdoGetObase());
  int nchars = strlen(filename);

  const char *refname = cdoGetObase();
  filesuffix[0] = 0;
  cdoGenFileSuffix(filesuffix, sizeof(filesuffix), pstreamInqFiletype(streamID1), vlistID1, refname);

  std::vector<double> array;
  //  if ( ! lcopy )
  {
    int gridsizemax = vlistGridsizeMax(vlistID1);
    if (vlistNumber(vlistID1) != CDI_REAL) gridsizemax *= 2;
    array.resize(gridsizemax);
  }

  int nvars = vlistNvars(vlistID1);
  int nconst = 0;
  for (varID = 0; varID < nvars; varID++)
    if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT) nconst++;

  std::vector<std::vector<Field>> vars;
  if (nconst)
    {
      vars.resize(nvars);

      for (varID = 0; varID < nvars; varID++)
        {
          if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT)
            {
              int gridID = vlistInqVarGrid(vlistID1, varID);
              int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
              size_t gridsize = gridInqSize(gridID);

              vars[varID].resize(nlevel);

              for (levelID = 0; levelID < nlevel; levelID++)
                {
                  field_init(&vars[varID][levelID]);
                  vars[varID][levelID].grid = gridID;
                  vars[varID][levelID].ptr = (double *) Malloc(gridsize * sizeof(double));
                }
            }
        }
    }

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      int64_t vdate = taxisInqVdate(taxisID1);
      int vtime = taxisInqVtime(taxisID1);

      if (operfunc == func_date)
        {
          index = (vdate / operintval) % 100;
          if (index < 0) index = -index;

          if (operatorID == SPLITSEAS) index = month_to_season(index);
        }
      else if (operfunc == func_time)
        {
          index = (vtime / operintval) % 100;
        }

      if (index < 0 || index >= MAX_STREAMS) cdoAbort("Index out of range!");

      streamID2 = streamIDs[index];
      if (streamID2 < 0)
        {
          if (operatorID == SPLITSEAS)
            {
              sprintf(filename + nchars, "%3s", seas_name[index]);
              if (filesuffix[0]) sprintf(filename + nchars + 3, "%s", filesuffix);
            }
          else
            {
              size_t slen;
              char oformat[32];
              strcpy(oformat, "%02d");

              if (operatorID == SPLITMON && format)
                {
                  char sbuf[32];
                  struct tm stime = datetime_to_tm(vdate, vtime);
                  slen = strftime(sbuf, 32, format, &stime);

                  if (slen) strcpy(oformat, sbuf);
                }

              slen = sprintf(filename + nchars, oformat, index);
              if (filesuffix[0]) sprintf(filename + nchars + slen, "%s", filesuffix);
            }

          if (cdoVerbose) cdoPrint("create file %s", filename);

          streamID2 = cdoStreamOpenWrite(filename, cdoFiletype());
          pstreamDefVlist(streamID2, vlistID2);
          streamIDs[index] = streamID2;
        }

      pstreamDefTimestep(streamID2, tsIDs[index]);

      if (tsID > 0 && tsIDs[index] == 0 && nconst)
        {
          for (varID = 0; varID < nvars; varID++)
            {
              if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT)
                {
                  int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
                  for (levelID = 0; levelID < nlevel; levelID++)
                    {
                      pstreamDefRecord(streamID2, varID, levelID);
                      nmiss = vars[varID][levelID].nmiss;
                      pstreamWriteRecord(streamID2, vars[varID][levelID].ptr, nmiss);
                    }
                }
            }
        }

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamDefRecord(streamID2, varID, levelID);

          if (lcopy && !(tsID == 0 && nconst))
            {
              pstreamCopyRecord(streamID2, streamID1);
            }
          else
            {
              pstreamReadRecord(streamID1, array.data(), &nmiss);
              pstreamWriteRecord(streamID2, array.data(), nmiss);

              if (tsID == 0 && nconst)
                {
                  if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT)
                    {
                      int gridID = vlistInqVarGrid(vlistID1, varID);
                      size_t gridsize = gridInqSize(gridID);
                      arrayCopy(gridsize, array.data(), vars[varID][levelID].ptr);
                      vars[varID][levelID].nmiss = nmiss;
                    }
                }
            }
        }

      tsIDs[index]++;
      tsID++;
    }

  pstreamClose(streamID1);

  if (nconst)
    {
      for (varID = 0; varID < nvars; varID++)
        {
          if (vlistInqVarTimetype(vlistID2, varID) == TIME_CONSTANT)
            {
              int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
              for (levelID = 0; levelID < nlevel; levelID++)
                if (vars[varID][levelID].ptr) Free(vars[varID][levelID].ptr);
            }
        }
    }

  for (index = 0; index < MAX_STREAMS; index++)
    {
      if (streamIDs[index] >= 0) pstreamClose(streamIDs[index]);
    }

  vlistDestroy(vlistID2);

  cdoFinish();

  return 0;
}

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

      Inttime    inttime         Time interpolation
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "datetime.h"

int get_tunits(const char *unit, int *incperiod, int *incunit, int *tunit);

void *
Inttime(void *process)
{
  int streamID2 = -1;
  int nlevel;
  int varID, levelID;
  int64_t vdate;
  int vtime;
  size_t gridsize;
  size_t offset;
  int incperiod = 0, incunit = 3600, tunit = TUNIT_HOUR;
  int year, month, day, hour, minute, second;
  double *single1, *single2;
  double *vardatap;

  cdoInitialize(process);

  operatorInputArg("date,time<,increment> (format YYYY-MM-DD,hh:mm:ss)");
  if (operatorArgc() < 2) cdoAbort("Too few arguments!");

  const char *datestr = operatorArgv()[0];
  const char *timestr = operatorArgv()[1];

  if (strchr(datestr, '-'))
    {
      year = 1;
      month = 1;
      day = 1;
      sscanf(datestr, "%d-%d-%d", &year, &month, &day);
      vdate = cdiEncodeDate(year, month, day);
    }
  else
    {
      vdate = parameter2long(datestr);
    }

  if (strchr(timestr, ':'))
    {
      hour = 0;
      minute = 0;
      second = 0;
      sscanf(timestr, "%d:%d:%d", &hour, &minute, &second);
      vtime = cdiEncodeTime(hour, minute, second);
    }
  else
    {
      vtime = parameter2int(timestr);
    }

  if (operatorArgc() == 3)
    {
      const char *timeunits = operatorArgv()[2];
      incperiod = (int) strtol(timeunits, NULL, 10);
      if (timeunits[0] == '-' || timeunits[0] == '+') timeunits++;
      while (isdigit((int) *timeunits)) timeunits++;

      get_tunits(timeunits, &incperiod, &incunit, &tunit);
    }
  /* increment in seconds */
  int64_t ijulinc = (int64_t) incperiod * incunit;

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  if (ijulinc == 0) vlistDefNtsteps(vlistID2, 1);

  int nvars = vlistNvars(vlistID1);

  int maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recinfo(maxrecs);

  size_t gridsizemax = vlistGridsizeMax(vlistID1);
  std::vector<double> array(gridsizemax);

  std::vector<std::vector<size_t>> nmiss1(nvars);
  std::vector<std::vector<size_t>> nmiss2(nvars);
  double **vardata1 = (double **) Malloc(nvars * sizeof(double *));
  double **vardata2 = (double **) Malloc(nvars * sizeof(double *));

  for (varID = 0; varID < nvars; varID++)
    {
      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      nmiss1[varID].resize(nlevel);
      nmiss2[varID].resize(nlevel);
      vardata1[varID] = (double *) Malloc(gridsize * nlevel * sizeof(double));
      vardata2[varID] = (double *) Malloc(gridsize * nlevel * sizeof(double));
    }

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  if (taxisHasBounds(taxisID2)) taxisDeleteBounds(taxisID2);
  vlistDefTaxis(vlistID2, taxisID2);

  int calendar = taxisInqCalendar(taxisID1);

  juldate_t juldate = juldate_encode(calendar, vdate, vtime);

  if (cdoVerbose)
    {
      cdoPrint("date %ld  time %d", vdate, vtime);
      cdoPrint("juldate  = %f", juldate_to_seconds(juldate));
      cdoPrint("ijulinc = %lld", ijulinc);
    }

  int tsID = 0;
  int nrecs = cdoStreamInqTimestep(streamID1, tsID++);
  juldate_t juldate1 = juldate_encode(calendar, taxisInqVdate(taxisID1), taxisInqVtime(taxisID1));
  for (int recID = 0; recID < nrecs; recID++)
    {
      pstreamInqRecord(streamID1, &varID, &levelID);
      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
      offset = gridsize * levelID;
      single1 = vardata1[varID] + offset;
      pstreamReadRecord(streamID1, single1, &nmiss1[varID][levelID]);
    }

  if (cdoVerbose)
    {
      cdoPrint("date %ld  time %d", taxisInqVdate(taxisID1), taxisInqVtime(taxisID1));
      cdoPrint("juldate1  = %f", juldate_to_seconds(juldate1));
    }

  if (juldate_to_seconds(juldate1) > juldate_to_seconds(juldate)) cdoWarning("start time %ld %d out of range!", vdate, vtime);

  int tsIDo = 0;
  while (juldate_to_seconds(juldate1) <= juldate_to_seconds(juldate))
    {
      nrecs = cdoStreamInqTimestep(streamID1, tsID++);
      if (nrecs == 0) break;

      juldate_t juldate2 = juldate_encode(calendar, taxisInqVdate(taxisID1), taxisInqVtime(taxisID1));
      if (cdoVerbose)
        {
          cdoPrint("date %ld  time %d", taxisInqVdate(taxisID1), taxisInqVtime(taxisID1));
          cdoPrint("juldate2  = %f", juldate_to_seconds(juldate2));
        }

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);

          recinfo[recID].varID = varID;
          recinfo[recID].levelID = levelID;
          recinfo[recID].lconst = vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT;

          gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
          offset = gridsize * levelID;
          single2 = vardata2[varID] + offset;
          pstreamReadRecord(streamID1, single2, &nmiss2[varID][levelID]);
        }

      while (juldate_to_seconds(juldate) <= juldate_to_seconds(juldate2))
        {
          if (juldate_to_seconds(juldate) >= juldate_to_seconds(juldate1)
              && juldate_to_seconds(juldate) <= juldate_to_seconds(juldate2))
            {
              juldate_decode(calendar, juldate, &vdate, &vtime);

              if (cdoVerbose)
                {
                  char vdatestr[32], vtimestr[32];
                  /*
                  cdoPrint("juldate1 %f", juldate_to_seconds(juldate1));
                  cdoPrint("juldate  %f", juldate_to_seconds(juldate));
                  cdoPrint("juldate2 %f", juldate_to_seconds(juldate2));
                  */
                  date2str(vdate, vdatestr, sizeof(vdatestr));
                  time2str(vtime, vtimestr, sizeof(vtimestr));
                  cdoPrint("%s %s  %f  %d", vdatestr, vtimestr, juldate_to_seconds(juldate), calendar);
                }

              if (streamID2 == -1)
                {
                  streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
                  pstreamDefVlist(streamID2, vlistID2);
                }

              taxisDefVdate(taxisID2, vdate);
              taxisDefVtime(taxisID2, vtime);
              pstreamDefTimestep(streamID2, tsIDo++);

              double fac1
                  = juldate_to_seconds(juldate_sub(juldate2, juldate)) / juldate_to_seconds(juldate_sub(juldate2, juldate1));
              double fac2
                  = juldate_to_seconds(juldate_sub(juldate, juldate1)) / juldate_to_seconds(juldate_sub(juldate2, juldate1));

              for (int recID = 0; recID < nrecs; recID++)
                {
                  varID = recinfo[recID].varID;
                  levelID = recinfo[recID].levelID;
                  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
                  offset = gridsize * levelID;
                  single1 = vardata1[varID] + offset;
                  single2 = vardata2[varID] + offset;

                  size_t nmiss3 = 0;

                  if (nmiss1[varID][levelID] > 0 || nmiss2[varID][levelID] > 0)
                    {
                      double missval1 = vlistInqVarMissval(vlistID1, varID);
                      double missval2 = vlistInqVarMissval(vlistID2, varID);

                      for (size_t i = 0; i < gridsize; i++)
                        {
                          if (!DBL_IS_EQUAL(single1[i], missval1) && !DBL_IS_EQUAL(single2[i], missval2))
                            array[i] = single1[i] * fac1 + single2[i] * fac2;
                          else if (DBL_IS_EQUAL(single1[i], missval1) && !DBL_IS_EQUAL(single2[i], missval2) && fac2 >= 0.5)
                            array[i] = single2[i];
                          else if (DBL_IS_EQUAL(single2[i], missval2) && !DBL_IS_EQUAL(single1[i], missval1) && fac1 >= 0.5)
                            array[i] = single1[i];
                          else
                            {
                              array[i] = missval1;
                              nmiss3++;
                            }
                        }
                    }
                  else
                    {
                      for (size_t i = 0; i < gridsize; i++) array[i] = single1[i] * fac1 + single2[i] * fac2;
                    }

                  pstreamDefRecord(streamID2, varID, levelID);
                  pstreamWriteRecord(streamID2, array.data(), nmiss3);
                }
            }

          if (ijulinc == 0) break;

          if (tunit == TUNIT_MONTH || tunit == TUNIT_YEAR)
            {
              juldate_decode(calendar, juldate, &vdate, &vtime);

              cdiDecodeDate(vdate, &year, &month, &day);

              month += (int) ijulinc;

              while (month > 12)
                {
                  month -= 12;
                  year++;
                }
              while (month < 1)
                {
                  month += 12;
                  year--;
                }

              vdate = cdiEncodeDate(year, month, day);

              juldate = juldate_encode(calendar, vdate, vtime);
            }
          else
            {
              juldate = juldate_add_seconds(ijulinc, juldate);
            }
        }

      juldate1 = juldate2;
      for (varID = 0; varID < nvars; varID++)
        {
          vardatap = vardata1[varID];
          vardata1[varID] = vardata2[varID];
          vardata2[varID] = vardatap;
        }
    }

  for (varID = 0; varID < nvars; varID++)
    {
      Free(vardata1[varID]);
      Free(vardata2[varID]);
    }

  Free(vardata1);
  Free(vardata2);

  if (streamID2 != -1) pstreamClose(streamID2);
  pstreamClose(streamID1);

  if (tsIDo == 0) cdoWarning("date/time out of time axis, no time step interpolated!");

  cdoFinish();

  return 0;
}

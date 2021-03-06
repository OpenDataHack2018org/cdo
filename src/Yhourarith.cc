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

      Yhourarith  yhouradd         Add multi-year hourly time series
      Yhourarith  yhoursub         Subtract multi-year hourly time series
      Yhourarith  yhourmul         Multiply multi-year hourly time series
      Yhourarith  yhourdiv         Divide multi-year hourly time series
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

#define MAX_HOUR 9301 /* 31*12*25 + 1 */

static int
hour_of_year(int64_t vdate, int vtime)
{
  int year, month, day, houroy;
  int hour, minute, second;

  cdiDecodeDate(vdate, &year, &month, &day);
  cdiDecodeTime(vtime, &hour, &minute, &second);

  if (month >= 1 && month <= 12 && day >= 1 && day <= 31 && hour >= 0 && hour < 24)
    houroy = ((month - 1) * 31 + day - 1) * 25 + hour + 1;
  else
    houroy = 0;

  if (houroy < 0 || houroy >= MAX_HOUR)
    {
      char vdatestr[32], vtimestr[32];
      date2str(vdate, vdatestr, sizeof(vdatestr));
      time2str(vtime, vtimestr, sizeof(vtimestr));
      cdoAbort("Hour of year %d out of range (%s %s)!", houroy, vdatestr, vtimestr);
    }

  return houroy;
}

void *
Yhourarith(void *process)
{
  int nrecs;
  int varID, levelID;
  size_t nmiss;

  cdoInitialize(process);

  cdoOperatorAdd("yhouradd", func_add, 0, NULL);
  cdoOperatorAdd("yhoursub", func_sub, 0, NULL);
  cdoOperatorAdd("yhourmul", func_mul, 0, NULL);
  cdoOperatorAdd("yhourdiv", func_div, 0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));
  int streamID2 = cdoStreamOpenRead(cdoStreamName(1));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = cdoStreamInqVlist(streamID2);
  int vlistID3 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, CMP_ALL);

  size_t gridsize = vlistGridsizeMax(vlistID1);

  Field field1, field2;
  field_init(&field1);
  field_init(&field2);
  field1.ptr = (double *) Malloc(gridsize * sizeof(double));
  field2.ptr = (double *) Malloc(gridsize * sizeof(double));

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = vlistInqTaxis(vlistID2);
  int taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID3);

  int streamID3 = cdoStreamOpenWrite(cdoStreamName(2), cdoFiletype());
  pstreamDefVlist(streamID3, vlistID3);

  int nvars = vlistNvars(vlistID2);

  std::vector<std::vector<std::vector<double>>> vardata2(MAX_HOUR);
  std::vector<std::vector<std::vector<size_t>>> varnmiss2(MAX_HOUR);

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID2, tsID)))
    {
      int64_t vdate = taxisInqVdate(taxisID2);
      int vtime = taxisInqVtime(taxisID2);

      int houroy = hour_of_year(vdate, vtime);
      if (vardata2[houroy].size() > 0) cdoAbort("Hour of year %d already allocatd!", houroy);

      vardata2[houroy].resize(nvars);
      varnmiss2[houroy].resize(nvars);

      for (varID = 0; varID < nvars; varID++)
        {
          size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
          size_t nlev = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
          vardata2[houroy][varID].resize(nlev * gridsize);
          varnmiss2[houroy][varID].resize(nlev);
        }

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID2, &varID, &levelID);
          size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
          size_t offset = gridsize * levelID;
          pstreamReadRecord(streamID2, &vardata2[houroy][varID][offset], &nmiss);
          varnmiss2[houroy][varID][levelID] = nmiss;
        }

      tsID++;
    }

  tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      int64_t vdate = taxisInqVdate(taxisID1);
      int vtime = taxisInqVtime(taxisID1);

      int houroy = hour_of_year(vdate, vtime);
      if (vardata2[houroy].size() == 0) cdoAbort("Hour of year %d not found!", houroy);

      taxisCopyTimestep(taxisID3, taxisID1);
      pstreamDefTimestep(streamID3, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamReadRecord(streamID1, field1.ptr, &nmiss);
          field1.nmiss = nmiss;
          field1.grid = vlistInqVarGrid(vlistID1, varID);
          field1.missval = vlistInqVarMissval(vlistID1, varID);

          size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
          size_t offset = gridsize * levelID;
          arrayCopy(gridsize, &vardata2[houroy][varID][offset], field2.ptr);
          field2.nmiss = varnmiss2[houroy][varID][levelID];
          field2.grid = vlistInqVarGrid(vlistID2, varID);
          field2.missval = vlistInqVarMissval(vlistID2, varID);

          farfun(&field1, field2, operfunc);

          pstreamDefRecord(streamID3, varID, levelID);
          pstreamWriteRecord(streamID3, field1.ptr, field1.nmiss);
        }

      tsID++;
    }

  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if (field1.ptr) Free(field1.ptr);
  if (field2.ptr) Free(field2.ptr);

  cdoFinish();

  return 0;
}

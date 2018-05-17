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

      Ymonarith  ymonadd         Add multi-year monthly time series
      Ymonarith  ymonsub         Subtract multi-year monthly time series
      Ymonarith  ymonmul         Multiply multi-year monthly time series
      Ymonarith  ymondiv         Divide multi-year monthly time series
      Ymonarith  yseasadd        Add multi-year seasonal time series
      Ymonarith  yseassub        Subtract multi-year seasonal time series
      Ymonarith  yseasmul        Multiply multi-year seasonal time series
      Ymonarith  yseasdiv        Divide multi-year seasonal time series
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

#define MAX_MON 12

void *
Ymonarith(void *process)
{
  enum
  {
    MONTHLY,
    SEASONAL
  };
  int nrecs;
  int varID, levelID;
  size_t nmiss;
  int year, mon, day;
  const char *seas_name[4];

  cdoInitialize(process);

  // clang-format off
  cdoOperatorAdd("ymonadd",  func_add, MONTHLY, NULL);
  cdoOperatorAdd("ymonsub",  func_sub, MONTHLY, NULL);
  cdoOperatorAdd("ymonmul",  func_mul, MONTHLY, NULL);
  cdoOperatorAdd("ymondiv",  func_div, MONTHLY, NULL);
  cdoOperatorAdd("yseasadd", func_add, SEASONAL, NULL);
  cdoOperatorAdd("yseassub", func_sub, SEASONAL, NULL);
  cdoOperatorAdd("yseasmul", func_mul, SEASONAL, NULL);
  cdoOperatorAdd("yseasdiv", func_div, SEASONAL, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);
  int opertype = cdoOperatorF2(operatorID);

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

  if (opertype == SEASONAL) get_season_name(seas_name);

  std::vector<std::vector<std::vector<double>>> vardata2(MAX_MON);
  std::vector<std::vector<std::vector<size_t>>> varnmiss2(MAX_MON);

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID2, tsID)))
    {
      int64_t vdate = taxisInqVdate(taxisID2);

      cdiDecodeDate(vdate, &year, &mon, &day);
      if (mon < 1 || mon > MAX_MON) cdoAbort("Month %d out of range!", mon);
      mon--;

      if (opertype == SEASONAL) mon = month_to_season(mon + 1);

      if (vardata2[mon].size())
        {
          if (opertype == SEASONAL)
            cdoAbort("Season %s already allocatd!", seas_name[mon]);
          else
            cdoAbort("Month %d already allocatd!", mon);
        }

      vardata2[mon].resize(nvars);
      varnmiss2[mon].resize(nvars);
      ;

      for (varID = 0; varID < nvars; varID++)
        {
          size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
          size_t nlev = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
          vardata2[mon][varID].resize(nlev * gridsize);
          varnmiss2[mon][varID].resize(nlev);
        }

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID2, &varID, &levelID);
          size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
          size_t offset = gridsize * levelID;
          pstreamReadRecord(streamID2, &vardata2[mon][varID][offset], &nmiss);
          varnmiss2[mon][varID][levelID] = nmiss;
        }

      tsID++;
    }

  tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      int64_t vdate = taxisInqVdate(taxisID1);

      cdiDecodeDate(vdate, &year, &mon, &day);
      if (mon < 1 || mon > MAX_MON) cdoAbort("Month %d out of range!", mon);
      mon--;

      if (opertype == SEASONAL) mon = month_to_season(mon + 1);

      if (vardata2[mon].size() == 0)
        {
          if (opertype == SEASONAL)
            cdoAbort("Season %s not found!", seas_name[mon]);
          else
            cdoAbort("Month %d not found!", mon);
        }

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

          arrayCopy(gridsize, &vardata2[mon][varID][offset], field2.ptr);
          field2.nmiss = varnmiss2[mon][varID][levelID];
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

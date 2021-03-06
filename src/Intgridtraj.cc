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

*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "interpol.h"
#include "datetime.h"

static int
readnextpos(FILE *fp, int calendar, juldate_t *juldate, double *xpos, double *ypos)
{
  int year = 0, month = 0, day = 0, hour = 0, minute = 0, second = 0;

  *xpos = 0;
  *ypos = 0;

  int stat = fscanf(fp, "%d-%d-%d %d:%d:%d %lf %lf", &year, &month, &day, &hour, &minute, &second, xpos, ypos);

  if (stat != EOF)
    {
      int date = cdiEncodeDate(year, month, day);
      int time = cdiEncodeTime(hour, minute, second);
      *juldate = juldate_encode(calendar, date, time);
    }

  return stat;
}

void *
Intgridtraj(void *process)
{
  int varID, levelID;
  int64_t vdate;
  int vtime;
  size_t nmiss;
  double point;
  double xpos, ypos;
  int calendar = CALENDAR_STANDARD;

  cdoInitialize(process);

  operatorInputArg("filename with grid trajectories");
  operatorCheckArgc(1);

  char *posfile = operatorArgv()[0];
  FILE *fp = fopen(posfile, "r");
  if (fp == NULL) cdoAbort("Open failed on %s!", posfile);

  juldate_t juldate;
  readnextpos(fp, calendar, &juldate, &xpos, &ypos);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);

  Field field1, field2;
  field_init(&field1);
  field_init(&field2);

  int nvars = vlistNvars(vlistID1);

  int maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recinfo(maxrecs);

  size_t gridsizemax = vlistGridsizeMax(vlistID1);
  std::vector<double> array(gridsizemax);

  double **vardata1 = (double **) Malloc(nvars * sizeof(double *));
  double **vardata2 = (double **) Malloc(nvars * sizeof(double *));

  for (varID = 0; varID < nvars; varID++)
    {
      size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
      size_t nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      vardata1[varID] = (double *) Malloc(gridsize * nlevel * sizeof(double));
      vardata2[varID] = (double *) Malloc(gridsize * nlevel * sizeof(double));
    }

  int gridID2 = gridCreate(GRID_TRAJECTORY, 1);
  gridDefXsize(gridID2, 1);
  gridDefYsize(gridID2, 1);
  gridDefXvals(gridID2, &xpos);
  gridDefYvals(gridID2, &ypos);

  int vlistID2 = vlistDuplicate(vlistID1);

  int ngrids = vlistNgrids(vlistID1);
  for (int index = 0; index < ngrids; index++)
    {
      int gridID1 = vlistGrid(vlistID1, index);

      if (gridInqType(gridID1) != GRID_LONLAT && gridInqType(gridID1) != GRID_GAUSSIAN)
        cdoAbort("Unsupported grid type: %s", gridNamePtr(gridInqType(gridID1)));

      vlistChangeGridIndex(vlistID2, index, gridID2);
    }

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = CDI_UNDEFID;

  int tsID = 0;
  int nrecs = cdoStreamInqTimestep(streamID1, tsID++);
  juldate_t juldate1 = juldate_encode(calendar, taxisInqVdate(taxisID1), taxisInqVtime(taxisID1));
  for (int recID = 0; recID < nrecs; recID++)
    {
      pstreamInqRecord(streamID1, &varID, &levelID);
      size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
      size_t offset = gridsize * levelID;
      double *single1 = vardata1[varID] + offset;
      pstreamReadRecord(streamID1, single1, &nmiss);
      if (nmiss) cdoAbort("Missing values unsupported for this operator!");
    }

  int tsIDo = 0;
  while (juldate_to_seconds(juldate1) <= juldate_to_seconds(juldate))
    {
      nrecs = cdoStreamInqTimestep(streamID1, tsID++);
      if (nrecs == 0) break;
      juldate_t juldate2 = juldate_encode(calendar, taxisInqVdate(taxisID1), taxisInqVtime(taxisID1));

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);

          recinfo[recID].varID = varID;
          recinfo[recID].levelID = levelID;
          recinfo[recID].lconst = vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT;

          size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
          size_t offset = gridsize * levelID;
          double *single2 = vardata2[varID] + offset;
          pstreamReadRecord(streamID1, single2, &nmiss);
          if (nmiss) cdoAbort("Missing values unsupported for this operator!");
        }

      while (juldate_to_seconds(juldate) < juldate_to_seconds(juldate2))
        {
          if (juldate_to_seconds(juldate) >= juldate_to_seconds(juldate1)
              && juldate_to_seconds(juldate) < juldate_to_seconds(juldate2))
            {
              if (streamID2 == CDI_UNDEFID)
                {
                  streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
                  pstreamDefVlist(streamID2, vlistID2);
                }

              juldate_decode(calendar, juldate, &vdate, &vtime);
              taxisDefVdate(taxisID2, vdate);
              taxisDefVtime(taxisID2, vtime);
              pstreamDefTimestep(streamID2, tsIDo++);

              double fac1
                  = juldate_to_seconds(juldate_sub(juldate2, juldate)) / juldate_to_seconds(juldate_sub(juldate2, juldate1));
              double fac2
                  = juldate_to_seconds(juldate_sub(juldate, juldate1)) / juldate_to_seconds(juldate_sub(juldate2, juldate1));
              /*
              printf("      %f %f %f %f %f\n", juldate_to_seconds(juldate),
                                               juldate_to_seconds(juldate1),
                                               juldate_to_seconds(juldate2),
              fac1, fac2);
              */
              for (int recID = 0; recID < nrecs; recID++)
                {
                  varID = recinfo[recID].varID;
                  levelID = recinfo[recID].levelID;
                  double missval = vlistInqVarMissval(vlistID1, varID);
                  int gridID1 = vlistInqVarGrid(vlistID1, varID);
                  size_t gridsize = gridInqSize(gridID1);
                  size_t offset = gridsize * levelID;
                  double *single1 = vardata1[varID] + offset;
                  double *single2 = vardata2[varID] + offset;

                  for (size_t i = 0; i < gridsize; i++) array[i] = single1[i] * fac1 + single2[i] * fac2;

                  field1.grid = gridID1;
                  field1.nmiss = nmiss;
                  field1.missval = missval;
                  field1.ptr = array.data();
                  field2.grid = gridID2;
                  field2.ptr = &point;
                  field2.nmiss = 0;

                  intgridbil(&field1, &field2);

                  pstreamDefRecord(streamID2, varID, levelID);
                  pstreamWriteRecord(streamID2, &point, nmiss);
                }
            }
          if (readnextpos(fp, calendar, &juldate, &xpos, &ypos) == EOF) break;
          gridDefXvals(gridID2, &xpos);
          gridDefYvals(gridID2, &ypos);
        }

      juldate1 = juldate2;
      for (varID = 0; varID < nvars; varID++)
        {
          double *vardatap = vardata1[varID];
          vardata1[varID] = vardata2[varID];
          vardata2[varID] = vardatap;
        }
    }

  if (tsIDo == 0)
    {
      juldate_decode(calendar, juldate, &vdate, &vtime);
      char vdatestr[32], vtimestr[32];
      date2str(vdate, vdatestr, sizeof(vdatestr));
      time2str(vtime, vtimestr, sizeof(vtimestr));
      cdoWarning("Date/time %s %s not found!", vdatestr, vtimestr);
    }

  fclose(fp);
  if (streamID2 != CDI_UNDEFID) pstreamClose(streamID2);
  pstreamClose(streamID1);

  for (varID = 0; varID < nvars; varID++)
    {
      Free(vardata1[varID]);
      Free(vardata2[varID]);
    }
  Free(vardata1);
  Free(vardata2);

  cdoFinish();

  return 0;
}

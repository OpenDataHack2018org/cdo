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

      Regres      regres           Regression
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

/* Same code as Trend ! */
void *
Regres(void *process)
{
  int nrecs;
  int varID, levelID;
  size_t nmiss;
  double temp1, temp2;
  enum
  {
    nwork = 5
  };
  Field **work[nwork];

  cdoInitialize(process);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  vlistDefNtsteps(vlistID2, 1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);
  /*
  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);
  */
  int streamID3 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID3, vlistID2);

  int maxrecs = vlistNrecs(vlistID1);
  std::vector<RecordInfo> recinfo(maxrecs);

  size_t gridsize = vlistGridsizeMax(vlistID1);

  Field field1, field2;
  field_init(&field1);
  field_init(&field2);
  field1.ptr = (double *) Malloc(gridsize * sizeof(double));
  field2.ptr = (double *) Malloc(gridsize * sizeof(double));

  for (int w = 0; w < nwork; w++) work[w] = field_calloc(vlistID1, FIELD_PTR);

  int64_t vdate = 0;
  int vtime = 0;
  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      tsID++;  // don't move this line !!!

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);

          if (tsID == 1)
            {
              recinfo[recID].varID = varID;
              recinfo[recID].levelID = levelID;
              recinfo[recID].lconst = vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT;
            }

          pstreamReadRecord(streamID1, field1.ptr, &nmiss);

          double missval = vlistInqVarMissval(vlistID1, varID);
          int gridID = vlistInqVarGrid(vlistID1, varID);
          size_t gridsize = gridInqSize(gridID);

          for (size_t i = 0; i < gridsize; i++)
            if (!DBL_IS_EQUAL(field1.ptr[i], missval))
              {
                work[0][varID][levelID].ptr[i] += tsID;
                work[1][varID][levelID].ptr[i] += tsID * tsID;
                work[2][varID][levelID].ptr[i] += tsID * field1.ptr[i];
                work[3][varID][levelID].ptr[i] += field1.ptr[i];
                work[4][varID][levelID].ptr[i]++;
              }
        }
    }

  taxisDefVdate(taxisID2, vdate);
  taxisDefVtime(taxisID2, vtime);
  /* pstreamDefTimestep(streamID2, 0); */
  pstreamDefTimestep(streamID3, 0);

  for (int recID = 0; recID < maxrecs; recID++)
    {
      int varID = recinfo[recID].varID;
      int levelID = recinfo[recID].levelID;

      double missval = vlistInqVarMissval(vlistID1, varID);
      int gridID = vlistInqVarGrid(vlistID1, varID);
      size_t gridsize = gridInqSize(gridID);

      double missval1 = missval;
      double missval2 = missval;

      for (size_t i = 0; i < gridsize; i++)
        {
          temp1 = SUBMN(work[2][varID][levelID].ptr[i], DIVMN(MULMN(work[0][varID][levelID].ptr[i], work[3][varID][levelID].ptr[i]),
                                                              work[4][varID][levelID].ptr[i]));
          temp2 = SUBMN(work[1][varID][levelID].ptr[i], DIVMN(MULMN(work[0][varID][levelID].ptr[i], work[0][varID][levelID].ptr[i]),
                                                              work[4][varID][levelID].ptr[i]));

          field2.ptr[i] = DIVMN(temp1, temp2);
          field1.ptr[i] = SUBMN(DIVMN(work[3][varID][levelID].ptr[i], work[4][varID][levelID].ptr[i]),
                                MULMN(DIVMN(work[0][varID][levelID].ptr[i], work[4][varID][levelID].ptr[i]), field2.ptr[i]));
        }

      /*
      nmiss = arrayNumMV(gridsize, field1.ptr, missval);
      pstreamDefRecord(streamID2, varID, levelID);
      pstreamWriteRecord(streamID2, field1.ptr, nmiss);
      */
      nmiss = arrayNumMV(gridsize, field2.ptr, missval);
      pstreamDefRecord(streamID3, varID, levelID);
      pstreamWriteRecord(streamID3, field2.ptr, nmiss);
    }

  for (int w = 0; w < nwork; w++) field_free(work[w], vlistID1);

  if (field1.ptr) Free(field1.ptr);
  if (field2.ptr) Free(field2.ptr);

  pstreamClose(streamID3);
  /* pstreamClose(streamID2); */
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

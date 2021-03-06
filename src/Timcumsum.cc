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

      Timcumsum    timcumsum         Cumulative sum over time
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

void *
Timcumsum(void *process)
{
  int nrecs;
  int varID, levelID;
  size_t nmiss;

  cdoInitialize(process);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  size_t gridsize = vlistGridsizeMax(vlistID1);
  if (vlistNumber(vlistID1) != CDI_REAL) gridsize *= 2;

  Field field;
  field_init(&field);
  field.ptr = (double *) Malloc(gridsize * sizeof(double));

  Field **vars1 = field_malloc(vlistID1, FIELD_PTR);

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);

          Field *pvars1 = &vars1[varID][levelID];

          gridsize = gridInqSize(pvars1->grid);

          if (tsID == 0)
            {
              pstreamReadRecord(streamID1, pvars1->ptr, &nmiss);
              // pvars1->nmiss = nmiss;
              if (nmiss)
                for (size_t i = 0; i < gridsize; ++i)
                  if (DBL_IS_EQUAL(pvars1->ptr[i], pvars1->missval)) pvars1->ptr[i] = 0;
            }
          else
            {
              pstreamReadRecord(streamID1, field.ptr, &nmiss);
              // field.nmiss   = nmiss;
              field.size = gridsize;
              field.grid = pvars1->grid;
              field.missval = pvars1->missval;

              if (nmiss)
                for (size_t i = 0; i < gridsize; ++i)
                  if (DBL_IS_EQUAL(field.ptr[i], pvars1->missval)) field.ptr[i] = 0;

              farfun(pvars1, field, func_sum);
            }

          pstreamDefRecord(streamID2, varID, levelID);
          pstreamWriteRecord(streamID2, pvars1->ptr, pvars1->nmiss);
        }

      tsID++;
    }

  field_free(vars1, vlistID1);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if (field.ptr) Free(field.ptr);

  cdoFinish();

  return 0;
}

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

     Deltat    deltat         Delta t
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

void *
Deltat(void *process)
{
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

  Field **vars = field_malloc(vlistID1, FIELD_PTR);

  size_t gridsizemax = vlistGridsizeMax(vlistID1);
  std::vector<double> array1(gridsizemax);
  std::vector<double> array2(gridsizemax);

  int tsID = 0;
  int nrecs = cdoStreamInqTimestep(streamID1, tsID);
  for (int recID = 0; recID < nrecs; ++recID)
    {
      pstreamInqRecord(streamID1, &varID, &levelID);
      pstreamReadRecord(streamID1, vars[varID][levelID].ptr, &nmiss);
      vars[varID][levelID].nmiss = nmiss;
    }

  tsID++;
  int tsID2 = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID2);

      for (int recID = 0; recID < nrecs; ++recID)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamReadRecord(streamID1, array1.data(), &nmiss);

          double missval = vars[varID][levelID].missval;
          double *array0 = vars[varID][levelID].ptr;
          size_t gridsize = vars[varID][levelID].size;
          if (nmiss || vars[varID][levelID].nmiss)
            {
              for (size_t i = 0; i < gridsize; ++i)
                {
                  if (DBL_IS_EQUAL(array0[i], missval) || DBL_IS_EQUAL(array1[i], missval))
                    array2[i] = missval;
                  else
                    array2[i] = array1[i] - array0[i];
                }

              nmiss = arrayNumMV(gridsize, array2.data(), missval);
            }
          else
            {
              for (size_t i = 0; i < gridsize; ++i) array2[i] = array1[i] - array0[i];
            }

          arrayCopy(gridsize, array1.data(), array0);

          pstreamDefRecord(streamID2, varID, levelID);
          pstreamWriteRecord(streamID2, array2.data(), nmiss);
        }

      tsID++;
      tsID2++;
    }

  field_free(vars, vlistID1);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

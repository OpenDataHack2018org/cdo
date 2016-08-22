/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2016 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
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
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Deltat(void *argument)
{
  int varID, levelID;
  int nmiss;

  cdoInitialize(argument);

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  field_t **vars = field_malloc(vlistID1, FIELD_PTR);
  
  int gridsizemax = vlistGridsizeMax(vlistID1);
  double *array1 = (double*) Malloc(gridsizemax*sizeof(double));
  double *array2 = (double*) Malloc(gridsizemax*sizeof(double));

  int tsID = 0;
  int nrecs = streamInqTimestep(streamID1, tsID);
  for ( int recID = 0; recID < nrecs; ++recID )
    {
      streamInqRecord(streamID1, &varID, &levelID);
      streamReadRecord(streamID1, vars[varID][levelID].ptr, &nmiss);
      vars[varID][levelID].nmiss = nmiss;
    }

  tsID++;
  int tsID2 = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      streamDefTimestep(streamID2, tsID2);

      for ( int recID = 0; recID < nrecs; ++recID )
        {
          streamInqRecord(streamID1, &varID, &levelID);
          streamReadRecord(streamID1, array1, &nmiss);

          double missval = vars[varID][levelID].missval;
          double *array0 = vars[varID][levelID].ptr;
          int gridsize = vars[varID][levelID].size;
          if ( nmiss || vars[varID][levelID].nmiss )
            {
              for ( int i = 0; i < gridsize; ++i )
                {
                  if ( DBL_IS_EQUAL(array0[i], missval) || DBL_IS_EQUAL(array1[i], missval) )
                    array2[i] = missval;
                  else
                    array2[i] = array1[i] - array0[i];
                }

              nmiss = 0;
              for ( int i = 0; i < gridsize; ++i )
                if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss++;
            }
          else
            {
              for ( int i = 0; i < gridsize; ++i )
                array2[i] = array1[i] - array0[i];
            }
          
          for ( int i = 0; i < gridsize; ++i ) array0[i] = array1[i];

          streamDefRecord(streamID2, varID, levelID);
          streamWriteRecord(streamID2, array2, nmiss);
        }
      
      tsID++;
      tsID2++;
    }

  Free(array1);
  Free(array2);
  field_free(vars, vlistID1);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return 0;
}

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

      Timcumsum    timcumsum         Cumulative sum over time
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Timcumsum(void *argument)
{
  int nrecs;
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

  int gridsize = vlistGridsizeMax(vlistID1);
  if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;

  field_t field;
  field_init(&field);
  field.ptr = (double*) Malloc(gridsize*sizeof(double));

  field_t **vars1 = field_malloc(vlistID1, FIELD_PTR);

  int tsID  = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
    
      for ( int recID = 0; recID < nrecs; recID++ )
        {
          streamInqRecord(streamID1, &varID, &levelID);

          field_t *pvar1 = &vars1[varID][levelID];
              
          gridsize = gridInqSize(pvar1->grid);

          if ( tsID == 0 )
            {
              streamReadRecord(streamID1, pvar1->ptr, &nmiss);
              // pvar1->nmiss = (size_t)nmiss;
              if ( nmiss )
                for ( int i = 0; i < gridsize; ++i )
                  if ( DBL_IS_EQUAL(pvar1->ptr[i], pvar1->missval) ) pvar1->ptr[i] = 0;
            }
          else
            {
              streamReadRecord(streamID1, field.ptr, &nmiss);
              // field.nmiss   = (size_t)nmiss;
              field.size    = gridsize;
              field.grid    = pvar1->grid;
              field.missval = pvar1->missval;

              if ( nmiss )
                for ( int i = 0; i < gridsize; ++i )
                  if ( DBL_IS_EQUAL(field.ptr[i], pvar1->missval) ) field.ptr[i] = 0;

              farfun(pvar1, field, func_sum);
            }
          
	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, pvar1->ptr, (int)pvar1->nmiss);
        }

      tsID++;
    }

  field_free(vars1, vlistID1);

  streamClose(streamID2);
  streamClose(streamID1);

  if ( field.ptr ) Free(field.ptr);

  cdoFinish();

  return 0;
}

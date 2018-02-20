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

     Subtrend   subtrend        Subtract trend
*/


#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"


void *Subtrend(void *process)
{
  int gridID, varID, levelID;
  size_t nmiss;

  cdoInitialize(process);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));
  int streamID2 = cdoStreamOpenRead(cdoStreamName(1));
  int streamID3 = cdoStreamOpenRead(cdoStreamName(2));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = cdoStreamInqVlist(streamID2);
  int vlistID3 = cdoStreamInqVlist(streamID3);
  int vlistID4 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, CMP_DIM);
  vlistCompare(vlistID1, vlistID3, CMP_DIM);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID4 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID4, taxisID4);

  int streamID4 = cdoStreamOpenWrite(cdoStreamName(3), cdoFiletype());
  pstreamDefVlist(streamID4, vlistID4);

  size_t gridsize = vlistGridsizeMax(vlistID1);

  field_type field1, field4;
  field_init(&field1);
  field_init(&field4);
  field1.ptr = (double*) Malloc(gridsize*sizeof(double));
  field4.ptr = (double*) Malloc(gridsize*sizeof(double));

  field_type **vars2 = field_malloc(vlistID1, FIELD_PTR);
  field_type **vars3 = field_malloc(vlistID1, FIELD_PTR);


  int tsID = 0;
  int nrecs = cdoStreamInqTimestep(streamID2, tsID);

  for ( int recID = 0; recID < nrecs; recID++ )
    {
      pstreamInqRecord(streamID2, &varID, &levelID);
      pstreamReadRecord(streamID2, vars2[varID][levelID].ptr, &nmiss);
    }

  tsID = 0;
  nrecs = cdoStreamInqTimestep(streamID3, tsID);

  for ( int recID = 0; recID < nrecs; recID++ )
    {
      pstreamInqRecord(streamID3, &varID, &levelID);
      pstreamReadRecord(streamID3, vars3[varID][levelID].ptr, &nmiss);
    }


  tsID = 0;
  while ( (nrecs = cdoStreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID4, taxisID1);
      pstreamDefTimestep(streamID4, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamReadRecord(streamID1, field1.ptr, &nmiss);

	  gridID   = vlistInqVarGrid(vlistID1, varID);
	  gridsize = gridInqSize(gridID);

	  double missval = vlistInqVarMissval(vlistID1, varID);
	  double missval1 = missval;
	  double missval2 = missval;
	  for ( size_t i = 0; i < gridsize; i++ )
	    field4.ptr[i] = SUBMN(field1.ptr[i], ADDMN(vars2[varID][levelID].ptr[i], MULMN(vars3[varID][levelID].ptr[i], tsID)));
    
	  nmiss = arrayNumMV(gridsize, field4.ptr, missval);
	  pstreamDefRecord(streamID4, varID, levelID);
	  pstreamWriteRecord(streamID4, field4.ptr, nmiss);
	}

      tsID++;
    }

  field_free(vars2, vlistID1);
  field_free(vars3, vlistID1);

  if ( field1.ptr ) Free(field1.ptr);
  if ( field4.ptr ) Free(field4.ptr);

  pstreamClose(streamID4);
  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2017 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
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

      Condc      ifthenc         If then constant
      Condc      ifnotthenc      If not then constant
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Condc(void *argument)
{
  int nrecs;
  int varID, levelID;
  int nmiss, nmiss2;
  int i;
  double missval;

  cdoInitialize(argument);

  // clang-format off
  int IFTHENC    = cdoOperatorAdd("ifthenc",    0, 0, NULL);
  int IFNOTTHENC = cdoOperatorAdd("ifnotthenc", 0, 0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();

  operatorInputArg("constant value");
  double rc = parameter2double(operatorArgv()[0]);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  nospec(vlistID1);

  int gridsize = vlistGridsizeMax(vlistID1);

  double *array1 = (double*) Malloc(gridsize*sizeof(double));
  double *array2 = (double*) Malloc(gridsize*sizeof(double));

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      pstreamDefTimestep(streamID2, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamReadRecord(streamID1, array1, &nmiss);

	  missval  = vlistInqVarMissval(vlistID1, varID);
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

	  if ( operatorID == IFTHENC )
	    {
	      for ( i = 0; i < gridsize; i++ )
		array2[i] = !DBL_IS_EQUAL(array1[i], missval) && !DBL_IS_EQUAL(array1[i], 0.) ? rc : missval;
	    }
	  else if ( operatorID == IFNOTTHENC )
	    {
	      for ( i = 0; i < gridsize; i++ )
		array2[i] = !DBL_IS_EQUAL(array1[i], missval) && DBL_IS_EQUAL(array1[i], 0.) ? rc : missval;
	    }
	  else
	    {
	      cdoAbort("Operator not implemented!");
	    }

	  nmiss2 = 0;
	  for ( i = 0; i < gridsize; i++ )
	    if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss2++;

	  pstreamDefRecord(streamID2, varID, levelID);
	  pstreamWriteRecord(streamID2, array2, nmiss2);
	}
      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if ( array2 ) Free(array2);
  if ( array1 ) Free(array1);

  cdoFinish();

  return 0;
}

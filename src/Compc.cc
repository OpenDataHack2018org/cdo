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

      Compc      eqc             Equal constant
      Compc      nec             Not equal constant
      Compc      lec             Less equal constant
      Compc      ltc             Less then constant
      Compc      gec             Greater equal constant
      Compc      gtc             Greater then constant
*/


#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"


void *Compc(void *process)
{
  int nrecs;
  int varID, levelID;
  size_t nmiss;

  cdoInitialize(process);

  int EQC = cdoOperatorAdd("eqc", 0, 0, NULL);
  int NEC = cdoOperatorAdd("nec", 0, 0, NULL);
  int LEC = cdoOperatorAdd("lec", 0, 0, NULL);
  int LTC = cdoOperatorAdd("ltc", 0, 0, NULL);
  int GEC = cdoOperatorAdd("gec", 0, 0, NULL);
  int GTC = cdoOperatorAdd("gtc", 0, 0, NULL);

  int operatorID = cdoOperatorID();

  operatorInputArg("constant value");
  double rc = parameter2double(operatorArgv()[0]);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  nospec(vlistID1);

  size_t gridsizemax = vlistGridsizeMax(vlistID1);

  double *array1 = (double*) Malloc(gridsizemax*sizeof(double));
  double *array2 = (double*) Malloc(gridsizemax*sizeof(double));

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
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

	  double missval  = vlistInqVarMissval(vlistID1, varID);
	  size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

	  bool rc_is_missval = DBL_IS_EQUAL(rc, missval);

          int datatype = vlistInqVarDatatype(vlistID1, varID);
          double rcv = (datatype == CDI_DATATYPE_FLT32) ? (float)rc : rc;

	  if ( operatorID == EQC )
	    {
	      for ( size_t i = 0; i < gridsize; i++ )
		array2[i] = DBL_IS_EQUAL(array1[i], missval) || rc_is_missval ? missval : IS_EQUAL(array1[i], rcv);
	    }
	  else if ( operatorID == NEC )
	    {
	      for ( size_t i = 0; i < gridsize; i++ )
		array2[i] = DBL_IS_EQUAL(array1[i], missval) || rc_is_missval ? missval : IS_NOT_EQUAL(array1[i], rcv);
	    }
	  else if ( operatorID == LEC )
	    {
	      for ( size_t i = 0; i < gridsize; i++ )
		array2[i] = DBL_IS_EQUAL(array1[i], missval) || rc_is_missval ? missval : array1[i] <= rcv;
	    }
	  else if ( operatorID == LTC )
	    {
	      for ( size_t i = 0; i < gridsize; i++ )
		array2[i] = DBL_IS_EQUAL(array1[i], missval) || rc_is_missval ? missval : array1[i] < rcv;
	    }
	  else if ( operatorID == GEC )
	    {
	      for ( size_t i = 0; i < gridsize; i++ )
		array2[i] = DBL_IS_EQUAL(array1[i], missval) || rc_is_missval ? missval : array1[i] >= rcv;
	    }
	  else if ( operatorID == GTC )
	    {
	      for ( size_t i = 0; i < gridsize; i++ )
		array2[i] = DBL_IS_EQUAL(array1[i], missval) || rc_is_missval ? missval : array1[i] > rcv;
	    }
	  else
	    {
	      cdoAbort("Operator not implemented!");
	    }

          size_t nmiss2 = 0;
	  for ( size_t i = 0; i < gridsize; i++ )
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

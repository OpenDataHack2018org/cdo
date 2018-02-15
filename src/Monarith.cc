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

      Monarith  monadd         Add monthly time series
      Monarith  monsub         Subtract monthly time series
      Monarith  monmul         Multiply monthly time series
      Monarith  mondiv         Divide monthly time series
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"


void *Monarith(void *process)
{
  int nrecs, nrecs2, nlev;
  int varID, levelID;
  int offset;
  size_t nmiss;
  int yearmon2 = -1;

  cdoInitialize(process);

  cdoOperatorAdd("monadd", func_add, 0, NULL);
  cdoOperatorAdd("monsub", func_sub, 0, NULL);
  cdoOperatorAdd("monmul", func_mul, 0, NULL);
  cdoOperatorAdd("mondiv", func_div, 0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));
  int streamID2 = cdoStreamOpenRead(cdoStreamName(1));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = cdoStreamInqVlist(streamID2);
  int vlistID3 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, CMP_ALL);
  
  size_t gridsize = vlistGridsizeMax(vlistID1);

  field_type field1, field2;
  field_init(&field1);
  field_init(&field2);
  field1.ptr = (double*) Malloc(gridsize*sizeof(double));
  field2.ptr = (double*) Malloc(gridsize*sizeof(double));

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = vlistInqTaxis(vlistID2);
  int taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID3);

  int streamID3 = cdoStreamOpenWrite(cdoStreamName(2), cdoFiletype());
  pstreamDefVlist(streamID3, vlistID3);

  int nvars  = vlistNvars(vlistID2);

  double **vardata2  = (double **) Malloc(nvars*sizeof(double *));
  size_t **varnmiss2 = (size_t **) Malloc(nvars*sizeof(size_t *));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
      nlev     = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
      vardata2[varID]  = (double*) Malloc(nlev*gridsize*sizeof(double));
      varnmiss2[varID] = (size_t*) Malloc(nlev*sizeof(size_t));
    }

  int tsID  = 0;
  int tsID2 = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      int vdate = taxisInqVdate(taxisID1);
      int yearmon1 = vdate / 100;

      if ( yearmon1 != yearmon2 )
	{
	  int year1, mon1;

	  year1 = yearmon1/100;
	  mon1  = yearmon1 - (yearmon1/100)*100;

	  if ( cdoVerbose ) cdoPrint("Process: Year = %4d  Month = %2d", year1, mon1);

	  nrecs2 = pstreamInqTimestep(streamID2, tsID2);
	  if ( nrecs2 == 0 )
	    cdoAbort("Missing year=%4d mon=%2d in %s!", year1, mon1, cdoGetStreamName(1).c_str());

	  vdate = taxisInqVdate(taxisID2);
	  yearmon2 = vdate / 100;

	  if ( yearmon1 != yearmon2 )
	    {
	      int year2 = yearmon2/100;
	      int mon2  = yearmon2 - (yearmon2/100)*100;
	      cdoAbort("Timestep %d in %s has wrong date! Current year=%4d mon=%2d, expected year=%4d mon=%2d",
		       tsID2+1, cdoGetStreamName(1).c_str(), year2, mon2, year1, mon1);
	    }

	  for ( int recID = 0; recID < nrecs2; recID++ )
	    {
	      pstreamInqRecord(streamID2, &varID, &levelID);

	      gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	      offset   = gridsize*levelID;

	      pstreamReadRecord(streamID2, vardata2[varID]+offset, &nmiss);
	      varnmiss2[varID][levelID] = nmiss;
	    }

	  tsID2++;
	}

      taxisCopyTimestep(taxisID3, taxisID1);
      pstreamDefTimestep(streamID3, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamReadRecord(streamID1, field1.ptr, &nmiss);
          field1.nmiss   = nmiss;
	  field1.grid    = vlistInqVarGrid(vlistID1, varID);
	  field1.missval = vlistInqVarMissval(vlistID1, varID);

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  offset   = gridsize*levelID;
	  arrayCopy(gridsize, vardata2[varID]+offset, field2.ptr);
	  field2.nmiss   = varnmiss2[varID][levelID];
	  field2.grid    = vlistInqVarGrid(vlistID2, varID);
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

  for ( varID = 0; varID < nvars; varID++ )
    {
      Free(vardata2[varID]);
      Free(varnmiss2[varID]);
    }

  Free(vardata2);
  Free(varnmiss2);

  if ( field1.ptr ) Free(field1.ptr);
  if ( field2.ptr ) Free(field2.ptr);

  cdoFinish();

  return 0;
}

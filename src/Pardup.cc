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

      Pardup     pardup          Duplicate parameters
      Pardup     parmul          Multiply parameters
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"


void *Pardup(void *process)
{
  int nrecs;
  int varID, varID2, levelID;
  long offset;
  int nmul = 0;
  size_t nmiss;
  int nlevel;
  double *single;

  cdoInitialize(process);

  int PARDUP = cdoOperatorAdd("pardup", 0, 0, NULL);
  int PARMUL = cdoOperatorAdd("parmul", 0, 0, NULL);

  int operatorID = cdoOperatorID();

  if ( operatorID == PARDUP )
    {
      nmul = 2;
    }
  else if ( operatorID == PARMUL )
    {
      operatorInputArg("number of multiply");
      nmul = parameter2int(operatorArgv()[0]);
    }
  else
    cdoAbort("operator not implemented!");

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int nvars    = vlistNvars(vlistID1);
  int nrecords = vlistNrecs(vlistID1);

  int *recVarID   = (int*) Malloc(nrecords*sizeof(int));
  int *recLevelID = (int*) Malloc(nrecords*sizeof(int));

  size_t gridsize = vlistGridsizeMax(vlistID1);
  double *array    = (double*) Malloc(gridsize*sizeof(double));
  double **vardata = (double **) Malloc(nvars*sizeof(double *));
  size_t **varnmiss = (size_t **) Malloc(nvars*sizeof(size_t *));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      vardata[varID]  = (double*) Malloc(gridsize*nlevel*sizeof(double));
      varnmiss[varID] = (size_t*) Malloc(nlevel*sizeof(size_t));
    }

  for ( int i = 1; i < nmul; i++ )
    {
      vlistCat(vlistID2, vlistID1);
      for ( varID = 0; varID < nvars; varID++ )
	vlistDefVarParam(vlistID2, varID+nvars*i, cdiEncodeParam(-(varID+nvars*i+1), 255, 255));
    }

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int tsID = 0;
  while ( (nrecs = cdoStreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);

	  recVarID[recID]   = varID;
	  recLevelID[recID] = levelID;

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  offset   = gridsize*levelID;
	  single   = vardata[varID] + offset;
  
	  pstreamReadRecord(streamID1, single, &nmiss);
	  varnmiss[varID][levelID] = nmiss;
	}

      for ( int i = 0; i < nmul; i++ )
	for ( int recID = 0; recID < nrecs; recID++ )
	  {
	    varID    = recVarID[recID];
	    varID2   = varID + i*nvars;
	    levelID  = recLevelID[recID];

	    gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	    offset   = gridsize*levelID;
	    single   = vardata[varID] + offset;
	    nmiss    = varnmiss[varID][levelID];

	    arrayCopy(gridsize, single, array);
	    pstreamDefRecord(streamID2,  varID2,  levelID);
	    pstreamWriteRecord(streamID2, array, nmiss);
	  }

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  for ( varID = 0; varID < nvars; varID++ ) Free(vardata[varID]);
  for ( varID = 0; varID < nvars; varID++ ) Free(varnmiss[varID]);
  Free(vardata);
  Free(varnmiss);
  Free(array);

  cdoFinish();

  return 0;
}

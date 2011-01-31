/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2011 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Isosurface(void *argument)
{
  int streamID1, streamID2;
  int vlistID1, vlistID2;
  int gridsize, nlevel;
  int recID, nrecs;
  int gridID;
  int i, offset;
  int tsID, varID, levelID;
  int nmiss, nvars;
  int zaxisID, nzaxis;
  double missval;
  int *vars = NULL;
  field_t *vars1 = NULL;
  field_t field;
  double *single;
  int taxisID1, taxisID2;

  cdoInitialize(argument);

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);

  vlistClearFlag(vlistID1);
  nvars = vlistNvars(vlistID1);
  for ( varID = 0; varID < nvars; varID++ )
    vlistDefFlag(vlistID1, varID, 0, TRUE);

  vlistID2 = vlistCreate();
  vlistCopyFlag(vlistID2, vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);
  nzaxis  = vlistNzaxis(vlistID1);
  for ( i = 0; i < nzaxis; i++ )
    if ( zaxisInqSize(vlistZaxis(vlistID1, i)) > 1 )
      vlistChangeZaxisIndex(vlistID2, i, zaxisID);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);

  field.ptr = (double *) malloc(gridsize*sizeof(double));

  vars  =     (int *) malloc(nvars*sizeof(int));
  vars1 = (field_t *) malloc(nvars*sizeof(field_t));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      zaxisID  = vlistInqVarZaxis(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(zaxisID);
      missval  = vlistInqVarMissval(vlistID1, varID);

      vars1[varID].grid    = gridID;
      vars1[varID].nmiss   = 0;
      vars1[varID].missval = missval;
      vars1[varID].ptr     = (double *) malloc(gridsize*nlevel*sizeof(double));
    }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( varID = 0; varID < nvars; varID++ )
	{
	  vars[varID] = FALSE;
	  vars1[varID].nmiss = 0;
	}

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  offset   = gridsize*levelID;
	  single   = vars1[varID].ptr + offset;
	  
	  streamReadRecord(streamID1, single, &nmiss);
	  vars1[varID].nmiss += nmiss;
	  vars[varID] = TRUE;
	}

      for ( varID = 0; varID < nvars; varID++ )
	{
	  streamDefRecord(streamID2, varID, 0);
	  streamWriteRecord(streamID2, vars1[varID].ptr, vars1[varID].nmiss);
	}

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vars[varID] )
	    {
	      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
		  offset   = gridsize*levelID;
		  single   = vars1[varID].ptr + offset;
		  streamDefRecord(streamID2, varID, levelID);
		  streamWriteRecord(streamID2, single, nmiss);
		}
	    }
	}

      tsID++;
    }

  for ( varID = 0; varID < nvars; varID++ ) free(vars1[varID].ptr);
  free(vars1);

  free(vars);

  free(field.ptr);

  streamClose(streamID2);
  streamClose(streamID1);

  vlistDestroy(vlistID2);

  cdoFinish();

  return (0);
}

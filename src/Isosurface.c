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
  int gridsize;
  int recID, nrecs;
  int gridID;
  int i;
  int tsID, varID, levelID;
  int nmiss, nvars;
  int zaxisID, nzaxis;
  double missval;
  field_t *vars1 = NULL, *vars2 = NULL, *samp1 = NULL;
  field_t field;
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

  vars1 = (field_t *) malloc(nvars*sizeof(field_t));
  samp1 = (field_t *) malloc(nvars*sizeof(field_t));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      missval  = vlistInqVarMissval(vlistID1, varID);

      vars1[varID].grid    = gridID;
      vars1[varID].nsamp   = 0;
      vars1[varID].nmiss   = 0;
      vars1[varID].missval = missval;
      vars1[varID].ptr     = (double *) malloc(gridsize*sizeof(double));
      samp1[varID].grid    = gridID;
      samp1[varID].nmiss   = 0;
      samp1[varID].missval = missval;
      samp1[varID].ptr     = NULL;
    }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
          vars1[varID].nsamp++;
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  streamReadRecord(streamID1, vars1[varID].ptr, &nmiss);
	  vars1[varID].nmiss = nmiss;
	}

      for ( varID = 0; varID < nvars; varID++ )
	{
	  streamDefRecord(streamID2, varID, 0);
	  streamWriteRecord(streamID2, vars1[varID].ptr, vars1[varID].nmiss);
	  vars1[varID].nsamp = 0;
	}

      tsID++;
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      free(vars1[varID].ptr);
      if ( samp1[varID].ptr ) free(samp1[varID].ptr);
    }

  free(vars1);
  free(samp1);

  if ( field.ptr ) free(field.ptr);

  streamClose(streamID2);
  streamClose(streamID1);

  vlistDestroy(vlistID2);

  cdoFinish();

  return (0);
}

/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult
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

      Ydaypctl   ydaypctl        Multi-year daily percentiles
*/


#include <stdio.h>
#include <math.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "field.h"
#include "percentiles.h"

#define  NDAY       373


void *Ydaypctl(void *argument)
{
  static const char *func = "Ydaypctl";
  int gridsize;
  int varID;
  int recID;
  int gridID;
  int vdate, vtime;
  int year, month, day, dayoy;
  int nrecs, nrecords;
  int levelID;
  int tsID;
  int otsID;
  long nsets[NDAY];
  int streamID1, streamID2, streamID3, streamID4;
  int vlistID1, vlistID2, vlistID3, vlistID4, taxisID1, taxisID2, taxisID3, taxisID4;
  int nmiss;
  int nvars, nlevels;
  int *recVarID, *recLevelID;
  int vdates1[NDAY], vtimes1[NDAY];
  int vdates2[NDAY], vtimes2[NDAY];
  double missval;
  field_t **vars1[NDAY];
  field_t field;
  double pn;
  HISTOGRAM_SET *hsets[NDAY];

  cdoInitialize(argument);
  cdoOperatorAdd("ydaypctl", func_pctl, 0, NULL);

  operatorInputArg("percentile number");
  pn = atof(operatorArgv()[0]);
      
  if ( !(pn > 0 && pn < 100) )
    cdoAbort("Illegal argument: percentile number %g is not in the range 0..100!", pn);

  for ( dayoy = 0; dayoy < NDAY; dayoy++ )
    {
      vars1[dayoy] = NULL;
      hsets[dayoy] = NULL;
      nsets[dayoy] = 0;
    }

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));
  streamID2 = streamOpenRead(cdoStreamName(1));
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));
  streamID3 = streamOpenRead(cdoStreamName(2));
  if ( streamID3 < 0 ) cdiError(streamID3, "Open failed on %s", cdoStreamName(2));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = streamInqVlist(streamID2);
  vlistID3 = streamInqVlist(streamID3);
  vlistID4 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, CMP_HRD);
  vlistCompare(vlistID1, vlistID3, CMP_HRD);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = vlistInqTaxis(vlistID2);
  taxisID3 = vlistInqTaxis(vlistID3);
  /* TODO - check that time axes 2 and 3 are equal */

  taxisID4 = taxisCreate(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID4, taxisID4);

  streamID4 = streamOpenWrite(cdoStreamName(3), cdoFiletype());
  if ( streamID4 < 0 ) cdiError(streamID4, "Open failed on %s", cdoStreamName(3));

  streamDefVlist(streamID4, vlistID4);

  nvars    = vlistNvars(vlistID1);
  nrecords = vlistNrecs(vlistID1);

  recVarID   = (int *) malloc(nrecords*sizeof(int));
  recLevelID = (int *) malloc(nrecords*sizeof(int));

  gridsize = vlistGridsizeMax(vlistID1);
  field.ptr = (double *) malloc(gridsize*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID2, tsID)) )
    {
      if ( nrecs != streamInqTimestep(streamID3, tsID) )
        cdoAbort("Number of records in time step %d of %s and %s are different!", tsID+1, cdoStreamName(1), cdoStreamName(2));
      
      vdate = taxisInqVdate(taxisID2);
      vtime = taxisInqVtime(taxisID2);
      
      if ( vdate != taxisInqVdate(taxisID3) || vtime != taxisInqVtime(taxisID3) )
        cdoAbort("Verification dates for time step %d of %s and %s are different!", tsID+1, cdoStreamName(1), cdoStreamName(2));
        
      if ( cdoVerbose ) cdoPrint("process timestep: %d %d %d", tsID+1, vdate, vtime);

      cdiDecodeDate(vdate, &year, &month, &day);

      if ( month >= 1 && month <= 12 )
	dayoy = (month-1)*31 + day;
      else
	dayoy = 0;

      if ( dayoy < 0 || dayoy >= NDAY )
	cdoAbort("Day %d out of range!", dayoy);

      vdates2[dayoy] = vdate;
      vtimes2[dayoy] = vtime;

      if ( vars1[dayoy] == NULL )
	{
	  vars1[dayoy] = (field_t **) malloc(nvars*sizeof(field_t *));
          hsets[dayoy] = hsetCreate(nvars);

	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      gridID   = vlistInqVarGrid(vlistID1, varID);
	      gridsize = gridInqSize(gridID);
	      nlevels  = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      missval  = vlistInqVarMissval(vlistID1, varID);

	      vars1[dayoy][varID] = (field_t *)  malloc(nlevels*sizeof(field_t));
              hsetCreateVarLevels(hsets[dayoy], varID, nlevels, gridID);
	      
	      for ( levelID = 0; levelID < nlevels; levelID++ )
		{
		  vars1[dayoy][varID][levelID].grid    = gridID;
		  vars1[dayoy][varID][levelID].nmiss   = 0;
		  vars1[dayoy][varID][levelID].missval = missval;
		  vars1[dayoy][varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
		}
	    }
	}
      
      for ( recID = 0; recID < nrecs; recID++ )
        {
          streamInqRecord(streamID2, &varID, &levelID);
	  streamReadRecord(streamID2, vars1[dayoy][varID][levelID].ptr, &nmiss);
          vars1[dayoy][varID][levelID].nmiss = nmiss;
        }
      for ( recID = 0; recID < nrecs; recID++ )
        {
          streamInqRecord(streamID3, &varID, &levelID);
	  streamReadRecord(streamID3, field.ptr, &nmiss);
          field.nmiss   = nmiss;
          field.grid    = vars1[dayoy][varID][levelID].grid;
	  field.missval = vars1[dayoy][varID][levelID].missval;
	  
	  hsetDefVarLevelBounds(hsets[dayoy], varID, levelID, &vars1[dayoy][varID][levelID], &field);
        }
      
      tsID++;
    }
  
  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      if ( cdoVerbose ) cdoPrint("process timestep: %d %d %d", tsID+1, vdate, vtime);

      cdiDecodeDate(vdate, &year, &month, &day);

      if ( month >= 1 && month <= 12 )
	dayoy = (month-1)*31 + day;
      else
	dayoy = 0;

      if ( dayoy < 0 || dayoy >= NDAY )
	cdoAbort("Day %d out of range!", dayoy);
	
      vdates1[dayoy] = vdate;
      vtimes1[dayoy] = vtime;
      
      if ( vars1[dayoy] == NULL )
        cdoAbort("No data for day %d in %s and %s", dayoy, cdoStreamName(1), cdoStreamName(2));
        
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  recVarID[recID]   = varID;
	  recLevelID[recID] = levelID;

	  streamReadRecord(streamID1, vars1[dayoy][varID][levelID].ptr, &nmiss);
	  vars1[dayoy][varID][levelID].nmiss = nmiss;
	      
	  hsetAddVarLevelValues(hsets[dayoy], varID, levelID, &vars1[dayoy][varID][levelID]);
	}

      nsets[dayoy]++;
      tsID++;
    }

  otsID = 0;
  for ( dayoy = 0; dayoy < NDAY; dayoy++ )
    if ( nsets[dayoy] )
      {
        if ( vdates1[dayoy] != vdates2[dayoy] )
          cdoAbort("Verification dates for day %d of %s, %s and %s are different!", dayoy, cdoStreamName(1), cdoStreamName(2), cdoStreamName(3));
        if ( vtimes1[dayoy] != vtimes2[dayoy] )
          cdoAbort("Verification times for day %d of %s, %s and %s are different!", dayoy, cdoStreamName(1), cdoStreamName(2), cdoStreamName(3));
        
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTime(vlistID1, varID) == TIME_CONSTANT ) continue;
	    nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      
	    for ( levelID = 0; levelID < nlevels; levelID++ )
	      hsetGetVarLevelPercentiles(&vars1[dayoy][varID][levelID], hsets[dayoy], varID, levelID, pn);
	  }

	taxisDefVdate(taxisID4, vdates1[dayoy]);
	taxisDefVtime(taxisID4, vtimes1[dayoy]);
	streamDefTimestep(streamID4, otsID++);

	for ( recID = 0; recID < nrecords; recID++ )
	  {
	    varID    = recVarID[recID];
	    levelID  = recLevelID[recID];

	    if ( otsID == 1 || vlistInqVarTime(vlistID1, varID) == TIME_VARIABLE )
	      {
		streamDefRecord(streamID4, varID, levelID);
		streamWriteRecord(streamID4, vars1[dayoy][varID][levelID].ptr, vars1[dayoy][varID][levelID].nmiss);
	      }
	  }
      }

  for ( dayoy = 0; dayoy < NDAY; dayoy++ )
    {
      if ( vars1[dayoy] != NULL )
	{
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      for ( levelID = 0; levelID < nlevels; levelID++ )
		free(vars1[dayoy][varID][levelID].ptr);
	      free(vars1[dayoy][varID]);
	    }
	  free(vars1[dayoy]); 
	  hsetDestroy(hsets[dayoy]);
	}
    }

  if ( field.ptr ) free(field.ptr);

  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);

  streamClose(streamID4);
  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (0);
}

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


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "percentiles_hist.h"

#define  NDAY       373

int getmonthday(int date);

void *Ydaypctl(void *argument)
{
  int varID;
  int gridID;
  int vdate, vtime;
  int year, month, day, dayoy;
  int nrecs;
  int levelID;
  int nmiss;
  int nlevels;
  int vdates1[NDAY], vtimes1[NDAY];
  int vdates2[NDAY];
  long nsets[NDAY];
  field_type **vars1[NDAY];
  HISTOGRAM_SET *hsets[NDAY];

  cdoInitialize(argument);
  cdoOperatorAdd("ydaypctl", func_pctl, 0, NULL);

  operatorInputArg("percentile number");
  double pn = parameter2double(operatorArgv()[0]);
  percentile_check_number(pn);

  for ( dayoy = 0; dayoy < NDAY; dayoy++ )
    {
      vars1[dayoy] = NULL;
      hsets[dayoy] = NULL;
      nsets[dayoy] = 0;
    }

  int streamID1 = pstreamOpenRead(cdoStreamName(0));
  int streamID2 = pstreamOpenRead(cdoStreamName(1));
  int streamID3 = pstreamOpenRead(cdoStreamName(2));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = pstreamInqVlist(streamID2);
  int vlistID3 = pstreamInqVlist(streamID3);
  int vlistID4 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, CMP_ALL);
  vlistCompare(vlistID1, vlistID3, CMP_ALL);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = vlistInqTaxis(vlistID2);
  int taxisID3 = vlistInqTaxis(vlistID3);
  /* TODO - check that time axes 2 and 3 are equal */

  int taxisID4 = taxisDuplicate(taxisID1);
  if ( taxisHasBounds(taxisID4) ) taxisDeleteBounds(taxisID4);
  vlistDefTaxis(vlistID4, taxisID4);

  int streamID4 = pstreamOpenWrite(cdoStreamName(3), cdoFiletype());
  pstreamDefVlist(streamID4, vlistID4);

  int nvars    = vlistNvars(vlistID1);
  int nrecords = vlistNrecs(vlistID1);

  int *recVarID   = (int*) Malloc(nrecords*sizeof(int));
  int *recLevelID = (int*) Malloc(nrecords*sizeof(int));

  int gridsize = vlistGridsizeMax(vlistID1);

  field_type field;
  field_init(&field);
  field.ptr = (double*) Malloc(gridsize*sizeof(double));

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID2, tsID)) )
    {
      if ( nrecs != pstreamInqTimestep(streamID3, tsID) )
        cdoAbort("Number of records at time step %d of %s and %s differ!", tsID+1, cdoStreamName(1)->args, cdoStreamName(2)->args);
      
      vdate = taxisInqVdate(taxisID2);
      vtime = taxisInqVtime(taxisID2);
      
      if ( vdate != taxisInqVdate(taxisID3) )
        cdoAbort("Verification dates at time step %d of %s and %s differ!", tsID+1, cdoStreamName(1)->args, cdoStreamName(2)->args);
        
      if ( cdoVerbose ) cdoPrint("process timestep: %d %d %d", tsID+1, vdate, vtime);

      cdiDecodeDate(vdate, &year, &month, &day);

      if ( month >= 1 && month <= 12 )
	dayoy = (month-1)*31 + day;
      else
	dayoy = 0;

      if ( dayoy < 0 || dayoy >= NDAY )
	cdoAbort("Day %d out of range!", dayoy);

      vdates2[dayoy] = vdate;

      if ( vars1[dayoy] == NULL )
	{
	  vars1[dayoy] = field_malloc(vlistID1, FIELD_PTR);
          hsets[dayoy] = hsetCreate(nvars);

	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      gridID   = vlistInqVarGrid(vlistID1, varID);
	      nlevels  = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));

              hsetCreateVarLevels(hsets[dayoy], varID, nlevels, gridID);
	    }
	}
      
      for ( int recID = 0; recID < nrecs; recID++ )
        {
          pstreamInqRecord(streamID2, &varID, &levelID);
	  pstreamReadRecord(streamID2, vars1[dayoy][varID][levelID].ptr, &nmiss);
          vars1[dayoy][varID][levelID].nmiss = nmiss;
        }
      for ( int recID = 0; recID < nrecs; recID++ )
        {
          pstreamInqRecord(streamID3, &varID, &levelID);
	  pstreamReadRecord(streamID3, field.ptr, &nmiss);
          field.nmiss   = nmiss;
          field.grid    = vars1[dayoy][varID][levelID].grid;
	  field.missval = vars1[dayoy][varID][levelID].missval;
	  
	  hsetDefVarLevelBounds(hsets[dayoy], varID, levelID, &vars1[dayoy][varID][levelID], &field);
        }
      
      tsID++;
    }
  
  tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
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
        cdoAbort("No data for day %d in %s and %s", dayoy, cdoStreamName(1)->args, cdoStreamName(2)->args);
        
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);

	  if ( tsID == 0 )
	    {
	      recVarID[recID]   = varID;
	      recLevelID[recID] = levelID;
	    }

	  pstreamReadRecord(streamID1, vars1[dayoy][varID][levelID].ptr, &nmiss);
	  vars1[dayoy][varID][levelID].nmiss = nmiss;
	      
	  hsetAddVarLevelValues(hsets[dayoy], varID, levelID, &vars1[dayoy][varID][levelID]);
	}

      nsets[dayoy]++;
      tsID++;
    }

  int otsID = 0;
  for ( dayoy = 0; dayoy < NDAY; dayoy++ )
    if ( nsets[dayoy] )
      {
        if ( getmonthday(vdates1[dayoy]) !=  getmonthday(vdates2[dayoy]) )
          cdoAbort("Verification dates for the day %d of %s and %s are different!",
                   dayoy, cdoStreamName(0)->args, cdoStreamName(1)->args);
        
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	    nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	      
	    for ( levelID = 0; levelID < nlevels; levelID++ )
	      hsetGetVarLevelPercentiles(&vars1[dayoy][varID][levelID], hsets[dayoy], varID, levelID, pn);
	  }

	taxisDefVdate(taxisID4, vdates1[dayoy]);
	taxisDefVtime(taxisID4, vtimes1[dayoy]);
	pstreamDefTimestep(streamID4, otsID);

	for ( int recID = 0; recID < nrecords; recID++ )
	  {
	    varID    = recVarID[recID];
	    levelID  = recLevelID[recID];

	    if ( otsID && vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

	    pstreamDefRecord(streamID4, varID, levelID);
	    pstreamWriteRecord(streamID4, vars1[dayoy][varID][levelID].ptr, vars1[dayoy][varID][levelID].nmiss);
	  }

	otsID++;
      }

  for ( dayoy = 0; dayoy < NDAY; dayoy++ )
    {
      if ( vars1[dayoy] != NULL )
	{
	  field_free(vars1[dayoy], vlistID1); 
	  hsetDestroy(hsets[dayoy]);
	}
    }

  if ( field.ptr ) Free(field.ptr);

  if ( recVarID   ) Free(recVarID);
  if ( recLevelID ) Free(recLevelID);

  pstreamClose(streamID4);
  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}
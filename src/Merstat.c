/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2006 Uwe Schulzweida, schulzweida@dkrz.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "functs.h"

/*
@BeginDoc

@BeginModule

@Name      = Merstat
@Title     = Meridional statistic
@Section   = Statistical description of the data
@Class     = Statistic
@Class     = Arithmetic
@Arguments = ifile ofile
@Operators = mermin mermax mersum mermean meravg merstd mervar

@EndModule


@BeginOperator_mermin

@Title     = Meridional minimum

@BeginDescription
For every longitude the minimum over all latitudes is computed.
@EndDescription

@EndOperator

@BeginOperator_mermax

@Title     = Meridional maximum

@BeginDescription
For every longitude the maximum over all latitudes is computed.
@EndDescription

@EndOperator

@BeginOperator_mersum

@Title     = Meridional sum

@BeginDescription
For every longitude the sum over all latitudes is computed.
@EndDescription

@EndOperator

@BeginOperator_mermean

@Title     = Meridional mean

@BeginDescription
For every longitude the mean over all latitudes is computed.
@EndDescription

@EndOperator

@BeginOperator_meravg

@Title     = Meridional average

@BeginDescription
For every longitude the average over all latitudes is computed.
@EndDescription

@EndOperator

@BeginOperator_mervar

@Title     = Meridional variance

@BeginDescription
For every longitude the variance over all latitudes is computed.
@EndDescription

@EndOperator

@BeginOperator_merstd

@Title     = Meridional standard deviation

@BeginDescription
For every longitude the standard deviation over all latitudes is computed.
@EndDescription

@EndOperator

@EndDoc
*/


void *Merstat(void *argument)
{
  static char func[] = "Merstat";
  int operatorID;
  int operfunc;
  int streamID1, streamID2;
  int vlistID1, vlistID2;
  int gridID1, gridID2, lastgrid = -1;
  int nlonmax;
  int index, ngrids;
  int recID, nrecs;
  int tsID, varID, levelID;
  int lim;
  FIELD field1, field2;
  int taxisID1, taxisID2;

  cdoInitialize(argument);

  cdoOperatorAdd("mermin",  func_min,  0, NULL);
  cdoOperatorAdd("mermax",  func_max,  0, NULL);
  cdoOperatorAdd("mersum",  func_sum,  0, NULL);
  cdoOperatorAdd("mermean", func_mean, 0, NULL);
  cdoOperatorAdd("meravg",  func_avg,  0, NULL);
  cdoOperatorAdd("mervar",  func_var,  0, NULL);
  cdoOperatorAdd("merstd",  func_std,  0, NULL);
 
  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  ngrids = vlistNgrids(vlistID1);
  index = 0;
  gridID1 = vlistGrid(vlistID1, index);
  if ( gridInqType(gridID1) != GRID_LONLAT &&
       gridInqType(gridID1) != GRID_GAUSSIAN )
    cdoAbort("Unsupported gridtype: %s", gridNamePtr(gridInqType(gridID1)));

  gridID2 = gridToMeridional(gridID1);
  vlistChangeGridIndex(vlistID2, index, gridID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  gridID1 = vlistInqVarGrid(vlistID1, 0);
  nlonmax = gridInqXsize(gridID1); /* max nlon ? */

  lim = vlistGridsizeMax(vlistID1);
  field1.ptr    = (double *) malloc(lim*sizeof(double));
  field1.weight = (double *) malloc(lim*sizeof(double));

  field2.ptr  = (double *) malloc(nlonmax*sizeof(double));
  field2.grid = gridID2;

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, field1.ptr, &field1.nmiss);

	  field1.grid    = vlistInqVarGrid(vlistID1, varID);
	  if ( field1.grid != lastgrid )
	    {
	      lastgrid = field1.grid;
	      gridWeights(field1.grid, field1.weight);
	    }
	  field1.missval = vlistInqVarMissval(vlistID1, varID);
	  field2.missval = vlistInqVarMissval(vlistID1, varID);

	  merfun(field1, &field2, operfunc);

	  streamDefRecord(streamID2, varID,  levelID);
	  streamWriteRecord(streamID2, field2.ptr, field2.nmiss);
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( field1.ptr )    free(field1.ptr);
  if ( field1.weight ) free(field1.weight);
  if ( field2.ptr )    free(field2.ptr);

  cdoFinish();

  return (0);
}

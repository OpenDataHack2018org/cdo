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

@Name      = Zonstat
@Title     = Zonal statistic
@Section   = Statistical description of the data
@Class     = Statistic
@Arguments = ifile ofile
@Operators = zonmin zonmax zonsum zonmean zonavg zonstd zonvar

@EndModule


@BeginOperator_zonmin

@Title     = Zonal minimum

@BeginDescription
For every latitude the minimum over all longitudes is computed.
@EndDescription

@EndOperator

@BeginOperator_zonmax

@Title     = Zonal maximum

@BeginDescription
For every latitude the maximum over all longitudes is computed.
@EndDescription

@EndOperator

@BeginOperator_zonsum

@Title     = Zonal sum

@BeginDescription
For every latitude the sum over all longitudes is computed.
@EndDescription

@EndOperator

@BeginOperator_zonmean

@Title     = Zonal mean

@BeginDescription
For every latitude the mean over all longitudes is computed.
@EndDescription

@EndOperator

@BeginOperator_zonavg

@Title     = Zonal average

@BeginDescription
For every latitude the average over all longitudes is computed.
@EndDescription

@EndOperator

@BeginOperator_zonvar

@Title     = Zonal variance

@BeginDescription
For every latitude the variance over all longitudes is computed.
@EndDescription

@EndOperator

@BeginOperator_zonstd

@Title     = Zonal standard deviation

@BeginDescription
For every latitude the standard deviation over all longitudes is computed.
@EndDescription

@EndOperator

@EndDoc
*/


void *Zonstat(void *argument)
{
  static char func[] = "Zonstat";
  int operatorID;
  int operfunc;
  int streamID1, streamID2;
  int vlistID1, vlistID2;
  int gridID1, gridID2;
  int nlatmax;
  int index, ngrids;
  int recID, nrecs;
  int tsID, varID, levelID;
  int lim;
  int taxisID1, taxisID2;
  FIELD field1, field2;

  cdoInitialize(argument);

  cdoOperatorAdd("zonmin",  func_min,  0, NULL);
  cdoOperatorAdd("zonmax",  func_max,  0, NULL);
  cdoOperatorAdd("zonsum",  func_sum,  0, NULL);
  cdoOperatorAdd("zonmean", func_mean, 0, NULL);
  cdoOperatorAdd("zonavg",  func_avg,  0, NULL);
  cdoOperatorAdd("zonvar",  func_var,  0, NULL);
  cdoOperatorAdd("zonstd",  func_std,  0, NULL);

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
  if ( ngrids > 1 ) cdoAbort("Too many different grids!");
  index = 0;
  gridID1 = vlistGrid(vlistID1, index);
  if ( gridInqType(gridID1) != GRID_LONLAT &&
       gridInqType(gridID1) != GRID_GAUSSIAN &&
       (gridInqType(gridID1) == GRID_GENERIC && gridInqYsize(gridID1) <= 1) )
    cdoAbort("Unsupported gridtype: %s", gridNamePtr(gridInqType(gridID1)));

  gridID2 = gridToZonal(gridID1);
  vlistChangeGridIndex(vlistID2, index, gridID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  gridID1 = vlistInqVarGrid(vlistID1, 0);
  nlatmax = gridInqYsize(gridID1); /* max nlat ? */

  lim = vlistGridsizeMax(vlistID1);
  field1.ptr  = (double *) malloc(lim*sizeof(double));
  field2.ptr  = (double *) malloc(nlatmax*sizeof(double));
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
	  field1.missval = vlistInqVarMissval(vlistID1, varID);
	  field2.missval = vlistInqVarMissval(vlistID1, varID);

	  zonfun(field1, &field2, operfunc);

	  streamDefRecord(streamID2, varID,  levelID);
	  streamWriteRecord(streamID2, field2.ptr, field2.nmiss);
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( field1.ptr ) free(field1.ptr);
  if ( field2.ptr ) free(field2.ptr);

  cdoFinish();

  return (0);
}

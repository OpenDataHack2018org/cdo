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

      Ydayarith  ydayadd         Add multi-year daily time series
      Ydayarith  ydaysub         Subtract multi-year daily time series
      Ydayarith  ydaymul         Multiply multi-year daily time series
      Ydayarith  ydaydiv         Divide multi-year daily time series
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  MAX_DOY   373

void *Ydayarith(void *argument)
{
  int nrecs;
  int varID, levelID;
  int nmiss;
  int year, month, day;
  int **varnmiss2[MAX_DOY];
  double **vardata2[MAX_DOY];

  cdoInitialize(argument);

  cdoOperatorAdd("ydayadd", func_add, 0, NULL);
  cdoOperatorAdd("ydaysub", func_sub, 0, NULL);
  cdoOperatorAdd("ydaymul", func_mul, 0, NULL);
  cdoOperatorAdd("ydaydiv", func_div, 0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  int streamID1 = streamOpenRead(cdoStreamName(0));
  int streamID2 = streamOpenRead(cdoStreamName(1));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = streamInqVlist(streamID2);
  int vlistID3 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, CMP_ALL);

  int gridsize = vlistGridsizeMax(vlistID1);

  field_t field1, field2;
  field_init(&field1);
  field_init(&field2);
  field1.ptr = (double*) Malloc(gridsize*sizeof(double));
  field2.ptr = (double*) Malloc(gridsize*sizeof(double));

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = vlistInqTaxis(vlistID2);
  int taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID3);

  int streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());

  streamDefVlist(streamID3, vlistID3);

  int nvars = vlistNvars(vlistID2);

  for ( int dayoy = 0; dayoy < MAX_DOY ; dayoy++ ) vardata2[dayoy] = NULL;

  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID2, tsID)) )
    {
      int vdate = taxisInqVdate(taxisID2);

      cdiDecodeDate(vdate, &year, &month, &day);

      int dayoy = 0;
      if ( month >= 1 && month <= 12 ) dayoy = (month-1)*31 + day;

      if ( dayoy < 0 || dayoy >= MAX_DOY )
	cdoAbort("Day of year %d out of range (date=%d)!", dayoy, vdate);

      if ( vardata2[dayoy] != NULL ) cdoAbort("Day of year %d already allocatd (date=%d)!", dayoy, vdate);

      vardata2[dayoy]  = (double **) Malloc(nvars*sizeof(double *));
      varnmiss2[dayoy] = (int **) Malloc(nvars*sizeof(int *));

      for ( varID = 0; varID < nvars; varID++ )
	{
          size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  size_t nlev     = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	  vardata2[dayoy][varID]  = (double*) Malloc(nlev*gridsize*sizeof(double));
	  varnmiss2[dayoy][varID] = (int*) Malloc(nlev*sizeof(int));
	}

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID2, &varID, &levelID);

          size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	  size_t offset   = gridsize*levelID;

	  streamReadRecord(streamID2, vardata2[dayoy][varID]+offset, &nmiss);
	  varnmiss2[dayoy][varID][levelID] = nmiss;
	}

      tsID++;
    }


  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      int vdate = taxisInqVdate(taxisID1);

      cdiDecodeDate(vdate, &year, &month, &day);
      
      int dayoy = 0;
      if ( month >= 1 && month <= 12 ) dayoy = (month-1)*31 + day;

      if ( dayoy < 0 || dayoy >= MAX_DOY )
	cdoAbort("Day of year %d out of range (date=%d)!", dayoy, vdate);

      taxisCopyTimestep(taxisID3, taxisID1);

      streamDefTimestep(streamID3, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, field1.ptr, &nmiss);
          field1.nmiss = (size_t) nmiss;

          size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
          size_t offset   = gridsize*levelID;
	  if ( vardata2[dayoy] == NULL ) cdoAbort("Day of year %d not found (date=%d)!", dayoy, vdate);
	  memcpy(field2.ptr, vardata2[dayoy][varID]+offset, gridsize*sizeof(double));
	  field2.nmiss = varnmiss2[dayoy][varID][levelID];

	  field1.grid    = vlistInqVarGrid(vlistID1, varID);
	  field1.missval = vlistInqVarMissval(vlistID1, varID);

	  field2.grid    = vlistInqVarGrid(vlistID2, varID);
	  field2.missval = vlistInqVarMissval(vlistID2, varID);

	  farfun(&field1, field2, operfunc);

          nmiss = (int) field1.nmiss;
	  streamDefRecord(streamID3, varID, levelID);
	  streamWriteRecord(streamID3, field1.ptr, nmiss);
	}
      tsID++;
    }

  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);

  for ( int dayoy = 0; dayoy < MAX_DOY; dayoy++ )
    if ( vardata2[dayoy] )
      {
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    Free(vardata2[dayoy][varID]);
	    Free(varnmiss2[dayoy][varID]);
	  }

        Free(vardata2[dayoy]);
        Free(varnmiss2[dayoy]);
      }

  if ( field1.ptr ) Free(field1.ptr);
  if ( field2.ptr ) Free(field2.ptr);

  cdoFinish();

  return 0;
}

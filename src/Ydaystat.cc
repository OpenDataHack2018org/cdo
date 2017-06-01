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

      Ydaystat   ydaymin         Multi-year daily minimum
      Ydaystat   ydaymax         Multi-year daily maximum
      Ydaystat   ydaysum         Multi-year daily sum
      Ydaystat   ydaymean        Multi-year daily mean
      Ydaystat   ydayavg         Multi-year daily average
      Ydaystat   ydayvar         Multi-year daily variance
      Ydaystat   ydayvar1        Multi-year daily variance [Normalize by (n-1)]
      Ydaystat   ydaystd         Multi-year daily standard deviation
      Ydaystat   ydaystd1        Multi-year daily standard deviation [Normalize by (n-1)]
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  MAX_DOY       373

typedef struct {
  short varID;
  short levelID;
} recinfo_t;


void *Ydaystat(void *argument)
{
  int varID, levelID;
  int year, month, day;
  int nrecs;
  int nsets[MAX_DOY];
  int nmiss;
  int vdates[MAX_DOY], vtimes[MAX_DOY];
  field_type **vars1[MAX_DOY], **vars2[MAX_DOY], **samp1[MAX_DOY];

  cdoInitialize(argument);

  cdoOperatorAdd("ydaymin",  func_min,  0, NULL);
  cdoOperatorAdd("ydaymax",  func_max,  0, NULL);
  cdoOperatorAdd("ydaysum",  func_sum,  0, NULL);
  cdoOperatorAdd("ydaymean", func_mean, 0, NULL);
  cdoOperatorAdd("ydayavg",  func_avg,  0, NULL);
  cdoOperatorAdd("ydayvar",  func_var,  0, NULL);
  cdoOperatorAdd("ydayvar1", func_var1, 0, NULL);
  cdoOperatorAdd("ydaystd",  func_std,  0, NULL);
  cdoOperatorAdd("ydaystd1", func_std1, 0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  bool lmean   = operfunc == func_mean || operfunc == func_avg;
  bool lstd    = operfunc == func_std || operfunc == func_std1;
  bool lvarstd = operfunc == func_std || operfunc == func_var || operfunc == func_std1 || operfunc == func_var1;
  int divisor  = operfunc == func_std1 || operfunc == func_var1;

  for ( int dayoy = 0; dayoy < MAX_DOY; dayoy++ )
    {
      vars1[dayoy] = NULL;
      vars2[dayoy] = NULL;
      samp1[dayoy] = NULL;
      nsets[dayoy] = 0;
    }

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  if ( taxisHasBounds(taxisID2) ) taxisDeleteBounds(taxisID2);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  int maxrecs = vlistNrecs(vlistID1);

  recinfo_t *recinfo = (recinfo_t *) Malloc(maxrecs*sizeof(recinfo_t));

  int gridsizemax = vlistGridsizeMax(vlistID1);

  field_type field;
  field_init(&field);
  field.ptr = (double*) Malloc(gridsizemax*sizeof(double));

  int tsID = 0;
  int otsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      int vdate = taxisInqVdate(taxisID1);
      int vtime = taxisInqVtime(taxisID1);

      if ( cdoVerbose ) cdoPrint("process timestep: %d %d %d", tsID+1, vdate, vtime);

      cdiDecodeDate(vdate, &year, &month, &day);

      int dayoy = 0;
      if ( month >= 1 && month <= 12 ) dayoy = (month-1)*31 + day;

      if ( dayoy < 0 || dayoy >= MAX_DOY )
	cdoAbort("Day of year %d out of range (date=%d)!", dayoy, vdate);

      vdates[dayoy] = vdate;
      vtimes[dayoy] = vtime;

      if ( vars1[dayoy] == NULL )
	{
	  vars1[dayoy] = field_malloc(vlistID1, FIELD_PTR);
	  samp1[dayoy] = field_malloc(vlistID1, FIELD_NONE);
	  if ( lvarstd )
	    vars2[dayoy] = field_malloc(vlistID1, FIELD_PTR);
	}

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  if ( tsID == 0 )
	    {
              recinfo[recID].varID   = varID;
              recinfo[recID].levelID = levelID;
	    }

          field_type *pvars1 = &vars1[dayoy][varID][levelID];
          field_type *pvars2 = vars2[dayoy] ? &vars2[dayoy][varID][levelID] : NULL;

	  int gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

	  if ( nsets[dayoy] == 0 )
	    {
	      streamReadRecord(streamID1, pvars1->ptr, &nmiss);
	      pvars1->nmiss = (size_t)nmiss;

	      if ( nmiss > 0 || samp1[dayoy][varID][levelID].ptr )
		{
		  if ( samp1[dayoy][varID][levelID].ptr == NULL )
		    samp1[dayoy][varID][levelID].ptr = (double*) Malloc(gridsize*sizeof(double));

		  for ( int i = 0; i < gridsize; i++ )
                    samp1[dayoy][varID][levelID].ptr[i] = !DBL_IS_EQUAL(pvars1->ptr[i], pvars1->missval);
		}
	    }
	  else
	    {
	      streamReadRecord(streamID1, field.ptr, &nmiss);
              field.nmiss   = (size_t)nmiss;
	      field.grid    = pvars1->grid;
	      field.missval = pvars1->missval;

	      if ( field.nmiss > 0 || samp1[dayoy][varID][levelID].ptr )
		{
		  if ( samp1[dayoy][varID][levelID].ptr == NULL )
		    {
		      samp1[dayoy][varID][levelID].ptr = (double*) Malloc(gridsize*sizeof(double));
		      for ( int i = 0; i < gridsize; i++ )
			samp1[dayoy][varID][levelID].ptr[i] = nsets[dayoy];
		    }
		  
		  for ( int i = 0; i < gridsize; i++ )
		    if ( !DBL_IS_EQUAL(field.ptr[i], pvars1->missval) )
		      samp1[dayoy][varID][levelID].ptr[i]++;
		}

	      if ( lvarstd )
		{
		  farsumq(pvars2, field);
		  farsum(pvars1, field);
		}
	      else
		{
		  farfun(pvars1, field, operfunc);
		}
	    }
	}

      if ( nsets[dayoy] == 0 && lvarstd )
        for ( int recID = 0; recID < maxrecs; recID++ )
          {
            int varID   = recinfo[recID].varID;
            int levelID = recinfo[recID].levelID;
            field_type *pvars1 = &vars1[dayoy][varID][levelID];
            field_type *pvars2 = &vars2[dayoy][varID][levelID];

	    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

            farmoq(pvars2, *pvars1);
	  }

      nsets[dayoy]++;
      tsID++;
    }

  // set the year to the minimum of years found on output timestep
  int outyear = 1e9;
  for ( int dayoy = 0; dayoy < MAX_DOY; dayoy++ )
    if ( nsets[dayoy] )
      {
        cdiDecodeDate(vdates[dayoy], &year, &month, &day);
        if ( year < outyear ) outyear = year;
      }
  for ( int dayoy = 0; dayoy < MAX_DOY; dayoy++ )
    if ( nsets[dayoy] )
      {
        cdiDecodeDate(vdates[dayoy], &year, &month, &day);
        if ( year > outyear ) vdates[dayoy] = cdiEncodeDate(outyear, month, day);
        //  printf("vdates[%d] = %d  nsets = %d\n", dayoy, vdates[dayoy], nsets[dayoy]);
      }

  for ( int dayoy = 0; dayoy < MAX_DOY; dayoy++ )
    if ( nsets[dayoy] )
      {
	if ( lmean )
          for ( int recID = 0; recID < maxrecs; recID++ )
            {
              int varID   = recinfo[recID].varID;
              int levelID = recinfo[recID].levelID;
              field_type *pvars1 = &vars1[dayoy][varID][levelID];

              if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

              if ( samp1[dayoy][varID][levelID].ptr == NULL )
                farcdiv(pvars1, (double)nsets[dayoy]);
              else
                fardiv(pvars1, samp1[dayoy][varID][levelID]);
	  }
	else if ( lvarstd )
          for ( int recID = 0; recID < maxrecs; recID++ )
            {
              int varID   = recinfo[recID].varID;
              int levelID = recinfo[recID].levelID;
              field_type *pvars1 = &vars1[dayoy][varID][levelID];
              field_type *pvars2 = &vars2[dayoy][varID][levelID];

              if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

              if ( samp1[dayoy][varID][levelID].ptr == NULL )
                {
                  if ( lstd ) farcstd(pvars1, *pvars2, nsets[dayoy], divisor);
                  else        farcvar(pvars1, *pvars2, nsets[dayoy], divisor);
                }
              else
                {
                  if ( lstd ) farstd(pvars1, *pvars2, samp1[dayoy][varID][levelID], divisor);
                  else        farvar(pvars1, *pvars2, samp1[dayoy][varID][levelID], divisor);
                }
            }

	taxisDefVdate(taxisID2, vdates[dayoy]);
	taxisDefVtime(taxisID2, vtimes[dayoy]);
	streamDefTimestep(streamID2, otsID);

        for ( int recID = 0; recID < maxrecs; recID++ )
          {
            int varID   = recinfo[recID].varID;
            int levelID = recinfo[recID].levelID;
            field_type *pvars1 = &vars1[dayoy][varID][levelID];

	    if ( otsID && vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

	    streamDefRecord(streamID2, varID, levelID);
	    streamWriteRecord(streamID2, pvars1->ptr, (int)pvars1->nmiss);
	  }

	otsID++;
      }

  for ( int dayoy = 0; dayoy < MAX_DOY; dayoy++ )
    {
      if ( vars1[dayoy] != NULL )
	{
	  field_free(vars1[dayoy], vlistID1);
	  field_free(samp1[dayoy], vlistID1);
	  if ( lvarstd ) field_free(vars2[dayoy], vlistID1);
	}
    }

  if ( field.ptr ) Free(field.ptr);

  Free(recinfo);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return 0;
}

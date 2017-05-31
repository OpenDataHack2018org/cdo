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

      Yseasstat  yseasmin        Multi-year seasonally minimum
      Yseasstat  yseasmax        Multi-year seasonally maximum
      Yseasstat  yseassum        Multi-year seasonally sum
      Yseasstat  yseasmean       Multi-year seasonally mean
      Yseasstat  yseasavg        Multi-year seasonally average
      Yseasstat  yseasvar        Multi-year seasonally variance
      Yseasstat  yseasvar1       Multi-year seasonally variance [Normalize by (n-1)]
      Yseasstat  yseasstd        Multi-year seasonally standard deviation
      Yseasstat  yseasstd1       Multi-year seasonally standard deviation [Normalize by (n-1)]
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"


#define  NSEAS       4

typedef struct {
  int vdate;
  int vtime;
}
date_time_t;

typedef struct {
  short varID;
  short levelID;
} recinfo_t;


void set_date(int vdate_new, int vtime_new, date_time_t *datetime)
{
  int year, month, day;
  cdiDecodeDate(vdate_new, &year, &month, &day);
  if ( month == 12 ) vdate_new = cdiEncodeDate(year-1, month, day);

  if ( vdate_new > datetime->vdate )
    {
      datetime->vdate = vdate_new;
      datetime->vtime = vtime_new;
    }
}


void *Yseasstat(void *argument)
{
  int varID;
  int year, month, day;
  int nrecs;
  int levelID;
  int nsets[NSEAS];
  int nmiss;
  date_time_t datetime[NSEAS];
  field_type **vars1[NSEAS], **vars2[NSEAS], **samp1[NSEAS];

  cdoInitialize(argument);

  cdoOperatorAdd("yseasmin",  func_min,  0, NULL);
  cdoOperatorAdd("yseasmax",  func_max,  0, NULL);
  cdoOperatorAdd("yseassum",  func_sum,  0, NULL);
  cdoOperatorAdd("yseasmean", func_mean, 0, NULL);
  cdoOperatorAdd("yseasavg",  func_avg,  0, NULL);
  cdoOperatorAdd("yseasvar",  func_var,  0, NULL);
  cdoOperatorAdd("yseasvar1", func_var1, 0, NULL);
  cdoOperatorAdd("yseasstd",  func_std,  0, NULL);
  cdoOperatorAdd("yseasstd1", func_std1, 0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  for ( int seas = 0; seas < NSEAS; seas++ )
    {
      vars1[seas]  = NULL;
      vars2[seas]  = NULL;
      samp1[seas]  = NULL;
      nsets[seas]  = 0;
      datetime[seas].vdate = 0;
      datetime[seas].vtime = 0;
    }

  bool lmean   = operfunc == func_mean || operfunc == func_avg;
  bool lstd    = operfunc == func_std || operfunc == func_std1;
  bool lvarstd = operfunc == func_std || operfunc == func_var || operfunc == func_std1 || operfunc == func_var1;
  int divisor  = operfunc == func_std1 || operfunc == func_var1;

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

  int gridsize = vlistGridsizeMax(vlistID1);

  field_type field;
  field_init(&field);
  field.ptr = (double*) Malloc(gridsize*sizeof(double));

  int tsID = 0;
  int otsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      int vdate = taxisInqVdate(taxisID1);
      int vtime = taxisInqVtime(taxisID1);
      cdiDecodeDate(vdate, &year, &month, &day);

      int seas = month_to_season(month);

      set_date(vdate, vtime, &datetime[seas]);

      if ( vars1[seas] == NULL )
	{
	  vars1[seas] = field_malloc(vlistID1, FIELD_PTR);
	  samp1[seas] = field_malloc(vlistID1, FIELD_NONE);
	  if ( lvarstd )
	    vars2[seas] = field_malloc(vlistID1, FIELD_PTR);
	}

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

	  if ( tsID == 0 )
	    {
              recinfo[recID].varID   = varID;
              recinfo[recID].levelID = levelID;
	    }

          field_type *pvars1 = &vars1[seas][varID][levelID];
          field_type *pvars2 = vars2[seas] ? &vars2[seas][varID][levelID] : NULL;

	  gridsize = pvars1->size;

	  if ( nsets[seas] == 0 )
	    {
	      streamReadRecord(streamID1, pvars1->ptr, &nmiss);
	      pvars1->nmiss = (size_t)nmiss;

	      if ( nmiss > 0 || samp1[seas][varID][levelID].ptr )
		{
		  if ( samp1[seas][varID][levelID].ptr == NULL )
		    samp1[seas][varID][levelID].ptr = (double*) Malloc(gridsize*sizeof(double));

		  for ( int i = 0; i < gridsize; i++ )
                    samp1[seas][varID][levelID].ptr[i] = !DBL_IS_EQUAL(pvars1->ptr[i], pvars1->missval);
		}
	    }
	  else
	    {
	      streamReadRecord(streamID1, field.ptr, &nmiss);
              field.nmiss   = (size_t)nmiss;
	      field.grid    = pvars1->grid;
	      field.missval = pvars1->missval;

	      if ( field.nmiss > 0 || samp1[seas][varID][levelID].ptr )
		{
		  if ( samp1[seas][varID][levelID].ptr == NULL )
		    {
		      samp1[seas][varID][levelID].ptr = (double*) Malloc(gridsize*sizeof(double));
		      for ( int i = 0; i < gridsize; i++ )
			samp1[seas][varID][levelID].ptr[i] = nsets[seas];
		    }
		  
		  for ( int i = 0; i < gridsize; i++ )
		    if ( !DBL_IS_EQUAL(field.ptr[i], pvars1->missval) )
		      samp1[seas][varID][levelID].ptr[i]++;
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

      if ( nsets[seas] == 0 && lvarstd )
        for ( int recID = 0; recID < maxrecs; recID++ )
          {
            int varID   = recinfo[recID].varID;
            int levelID = recinfo[recID].levelID;
            field_type *pvars1 = &vars1[seas][varID][levelID];
            field_type *pvars2 = &vars2[seas][varID][levelID];

	    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

            farmoq(pvars2, *pvars1);
	  }

      nsets[seas]++;
      tsID++;
    }

  for ( int seas = 0; seas < NSEAS; seas++ )
    if ( nsets[seas] )
      {
	if ( lmean )
          for ( int recID = 0; recID < maxrecs; recID++ )
            {
              int varID   = recinfo[recID].varID;
              int levelID = recinfo[recID].levelID;
              field_type *pvars1 = &vars1[seas][varID][levelID];

	      if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
              
              if ( samp1[seas][varID][levelID].ptr == NULL )
                farcdiv(pvars1, (double)nsets[seas]);
              else
                fardiv(pvars1, samp1[seas][varID][levelID]);
	    }
	else if ( lvarstd )
          for ( int recID = 0; recID < maxrecs; recID++ )
            {
              int varID   = recinfo[recID].varID;
              int levelID = recinfo[recID].levelID;
              field_type *pvars1 = &vars1[seas][varID][levelID];
              field_type *pvars2 = &vars2[seas][varID][levelID];

	      if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

              if ( samp1[seas][varID][levelID].ptr == NULL )
                {
                  if ( lstd ) farcstd(pvars1, *pvars2, nsets[seas], divisor);
                  else        farcvar(pvars1, *pvars2, nsets[seas], divisor);
                }
              else
                {
                  if ( lstd ) farstd(pvars1, *pvars2, samp1[seas][varID][levelID], divisor);
                  else        farvar(pvars1, *pvars2, samp1[seas][varID][levelID], divisor);
		}
	    }

	taxisDefVdate(taxisID2, datetime[seas].vdate);
	taxisDefVtime(taxisID2, datetime[seas].vtime);
	streamDefTimestep(streamID2, otsID);

	for ( int recID = 0; recID < maxrecs; recID++ )
	  {
            int varID   = recinfo[recID].varID;
            int levelID = recinfo[recID].levelID;
            field_type *pvars1 = &vars1[seas][varID][levelID];

	    if ( otsID && vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

	    streamDefRecord(streamID2, varID, levelID);
	    streamWriteRecord(streamID2, pvars1->ptr, (int)pvars1->nmiss);
	  }

	otsID++;
      }

  for ( int seas = 0; seas < NSEAS; seas++ )
    {
      if ( vars1[seas] != NULL )
	{
	  field_free(vars1[seas], vlistID1);
	  field_free(samp1[seas], vlistID1);
	  if ( lvarstd ) Free(vars2[seas]);
	}
    }

  if ( field.ptr ) Free(field.ptr);

  Free(recinfo);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return 0;
}

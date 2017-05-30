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

      Seasstat   seasrange       Seasonal range
      Seasstat   seasmin         Seasonal minimum
      Seasstat   seasmax         Seasonal maximum
      Seasstat   seassum         Seasonal sum
      Seasstat   seasmean        Seasonal mean
      Seasstat   seasavg         Seasonal average
      Seasstat   seasvar         Seasonal variance
      Seasstat   seasvar1        Seasonal variance [Normalize by (n-1)]
      Seasstat   seasstd         Seasonal standard deviation
      Seasstat   seasstd1        Seasonal standard deviation [Normalize by (n-1)]
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


typedef struct {
  short varID;
  short levelID;
} recinfo_t;


void *Seasstat(void *argument)
{
  int timestat_date = TIMESTAT_MEAN;
  int vdate0 = 0, vtime0 = 0;
  int vdate1 = 0, vtime1 = 0;
  int nrecs;
  int varID, levelID;
  int year, month, day, seas0 = 0;
  int nmiss;
  int oldmon = 0;
  int nseason = 0;
  const char *seas_name[4];

  cdoInitialize(argument);

  cdoOperatorAdd("seasrange", func_range, 0, NULL);
  cdoOperatorAdd("seasmin",   func_min,   0, NULL);
  cdoOperatorAdd("seasmax",   func_max,   0, NULL);
  cdoOperatorAdd("seassum",   func_sum,   0, NULL);
  cdoOperatorAdd("seasmean",  func_mean,  0, NULL);
  cdoOperatorAdd("seasavg",   func_avg,   0, NULL);
  cdoOperatorAdd("seasvar",   func_var,   0, NULL);
  cdoOperatorAdd("seasvar1",  func_var1,  0, NULL);
  cdoOperatorAdd("seasstd",   func_std,   0, NULL);
  cdoOperatorAdd("seasstd1",  func_std1,  0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  int season_start = get_season_start();
  get_season_name(seas_name);

  bool lrange  = operfunc == func_range;
  bool lmean   = operfunc == func_mean || operfunc == func_avg;
  bool lstd    = operfunc == func_std || operfunc == func_std1;
  bool lvarstd = operfunc == func_std || operfunc == func_var || operfunc == func_std1 || operfunc == func_var1;
  int  divisor = operfunc == func_std1 || operfunc == func_var1;

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  if ( taxisInqType(taxisID2) == TAXIS_FORECAST ) taxisDefType(taxisID2, TAXIS_RELATIVE);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  int maxrecs = vlistNrecs(vlistID1);

  recinfo_t *recinfo = (recinfo_t *) Malloc(maxrecs*sizeof(recinfo_t));

  dtlist_type *dtlist = dtlist_new();
  dtlist_set_stat(dtlist, timestat_date);
  dtlist_set_calendar(dtlist, taxisInqCalendar(taxisID1));

  int gridsize = vlistGridsizeMax(vlistID1);

  field_type field;
  field_init(&field);
  field.ptr = (double*) Malloc(gridsize*sizeof(double));

  field_type **samp1 = field_malloc(vlistID1, FIELD_NONE);
  field_type **vars1 = field_malloc(vlistID1, FIELD_PTR);
  field_type **vars2 = NULL;
  if ( lvarstd || lrange ) vars2 = field_malloc(vlistID1, FIELD_PTR);

  int tsID  = 0;
  int otsID = 0;
  while ( TRUE )
    {
      long nsets = 0;
      bool newseas = false;
      while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
	{
	  dtlist_taxisInqTimestep(dtlist, taxisID1, nsets);
	  int vdate = dtlist_get_vdate(dtlist, nsets);
	  int vtime = dtlist_get_vtime(dtlist, nsets);

	  cdiDecodeDate(vdate, &year, &month, &day);

	  int newmon = month;

	  if ( season_start == START_DEC && newmon == 12 ) newmon = 0;

          int seas = month_to_season(month);

	  if ( nsets == 0 )
	    {
	      nseason++;
	      vdate0 = vdate;
	      vtime0 = vtime;
	      seas0  = seas;
	      oldmon = newmon;
	    }

	  if ( newmon < oldmon ) newseas = true;

	  if ( (seas != seas0) || newseas ) break;

	  oldmon = newmon;

	  for ( int recID = 0; recID < nrecs; recID++ )
	    {
	      streamInqRecord(streamID1, &varID, &levelID);

	      if ( tsID == 0 )
		{
                  recinfo[recID].varID   = varID;
                  recinfo[recID].levelID = levelID;
		}

              field_type *pvar1 = &vars1[varID][levelID];

	      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

	      if ( nsets == 0 )
		{
		  streamReadRecord(streamID1, pvar1->ptr, &nmiss);
		  pvar1->nmiss = (size_t)nmiss;
                  if ( lrange )
                    {
                      field_type *pvar2 = &vars2[varID][levelID];
		      for ( int i = 0; i < gridsize; i++ )
                        pvar2->ptr[i] = pvar1->ptr[i];
                      pvar2->nmiss = (size_t)nmiss;
                    }

		  if ( nmiss > 0 || samp1[varID][levelID].ptr )
		    {
		      if ( samp1[varID][levelID].ptr == NULL )
			samp1[varID][levelID].ptr = (double*) Malloc(gridsize*sizeof(double));

		      for ( int i = 0; i < gridsize; i++ )
                        samp1[varID][levelID].ptr[i] = !DBL_IS_EQUAL(pvar1->ptr[i], pvar1->missval);
		    }
		}
	      else
		{
		  streamReadRecord(streamID1, field.ptr, &nmiss);
                  field.nmiss   = (size_t)nmiss;
		  field.grid    = pvar1->grid;
		  field.missval = pvar1->missval;

		  if ( field.nmiss > 0 || samp1[varID][levelID].ptr )
		    {
		      if ( samp1[varID][levelID].ptr == NULL )
			{
			  samp1[varID][levelID].ptr = (double*) Malloc(gridsize*sizeof(double));
			  for ( int i = 0; i < gridsize; i++ )
			    samp1[varID][levelID].ptr[i] = nsets;
			}

		      for ( int i = 0; i < gridsize; i++ )
			if ( !DBL_IS_EQUAL(field.ptr[i], pvar1->missval) )
			  samp1[varID][levelID].ptr[i]++;
		    }

		  if ( lvarstd )
		    {
                      field_type *pvar2 = &vars2[varID][levelID];
		      farsumq(pvar2, field);
		      farsum(pvar1, field);
		    }
                  else if ( lrange )
                    {
                      field_type *pvar2 = &vars2[varID][levelID];
                      farmin(pvar2, field);
                      farmax(pvar1, field);
                    }
		  else
		    {
		      farfun(pvar1, field, operfunc);
		    }
		}
	    }

	  if ( nsets == 0 && lvarstd )
            for ( int recID = 0; recID < maxrecs; recID++ )
              {
                int varID   = recinfo[recID].varID;
                int levelID = recinfo[recID].levelID;
                field_type *pvar1 = &vars1[varID][levelID];
                field_type *pvar2 = &vars2[varID][levelID];

		if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

                farmoq(pvar2, *pvar1);
	      }

	  vdate1 = vdate;
	  vtime1 = vtime;
	  nsets++;
	  tsID++;
	}

      if ( nrecs == 0 && nsets == 0 ) break;

      if ( lmean )
        for ( int recID = 0; recID < maxrecs; recID++ )
          {
            int varID   = recinfo[recID].varID;
            int levelID = recinfo[recID].levelID;
            field_type *pvar1 = &vars1[varID][levelID];

	    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

            if ( samp1[varID][levelID].ptr == NULL )
              farcdiv(pvar1, (double)nsets);
            else
              fardiv(pvar1, samp1[varID][levelID]);
	  }
      else if ( lvarstd )
        for ( int recID = 0; recID < maxrecs; recID++ )
          {
            int varID   = recinfo[recID].varID;
            int levelID = recinfo[recID].levelID;
            field_type *pvar1 = &vars1[varID][levelID];
            field_type *pvar2 = &vars2[varID][levelID];

            if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

            if ( samp1[varID][levelID].ptr == NULL )
              {
                if ( lstd ) farcstd(pvar1, *pvar2, nsets, divisor);
                else        farcvar(pvar1, *pvar2, nsets, divisor);
              }
            else
              {
                if ( lstd ) farstd(pvar1, *pvar2, samp1[varID][levelID], divisor);
                else        farvar(pvar1, *pvar2, samp1[varID][levelID], divisor);
	      }
	  }
      else if ( lrange )
        for ( int recID = 0; recID < maxrecs; recID++ )
          {
            int varID   = recinfo[recID].varID;
            int levelID = recinfo[recID].levelID;
            field_type *pvar1 = &vars1[varID][levelID];
            field_type *pvar2 = &vars2[varID][levelID];

            if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

            farsub(pvar1, *pvar2);
	  }

      if ( cdoVerbose )
	{
	  char vdatestr0[32], vtimestr0[32];
	  char vdatestr1[32], vtimestr1[32];
	  date2str(vdate0, vdatestr0, sizeof(vdatestr0));
	  time2str(vtime0, vtimestr0, sizeof(vtimestr0));
	  date2str(vdate1, vdatestr1, sizeof(vdatestr1));
	  time2str(vtime1, vtimestr1, sizeof(vtimestr1));
	  cdoPrint("season: %3d %3s  start: %s %s  end: %s %s ntimesteps: %d", 
		   nseason, seas_name[seas0], vdatestr0, vtimestr0, vdatestr1, vtimestr1, nsets);
	}

      dtlist_stat_taxisDefTimestep(dtlist, taxisID2, nsets);
      streamDefTimestep(streamID2, otsID);

      if ( nsets < 3 )
	{
	  char vdatestr[32];
	  date2str(vdate0, vdatestr, sizeof(vdatestr));
	  cdoWarning("Season %3d (%s) has only %d input time step%s!", 
		     otsID+1, vdatestr, nsets, nsets == 1 ? "" : "s");
	}

      for ( int recID = 0; recID < maxrecs; recID++ )
	{
          int varID   = recinfo[recID].varID;
          int levelID = recinfo[recID].levelID;
          field_type *pvar1 = &vars1[varID][levelID];

	  if ( otsID && vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, pvar1->ptr, (int)pvar1->nmiss);
	}

      if ( nrecs == 0 ) break;
      otsID++;
    }


  field_free(vars1, vlistID1);
  field_free(samp1, vlistID1);
  if ( lvarstd ) field_free(vars2, vlistID1);

  if ( field.ptr ) Free(field.ptr);

  Free(recinfo);

  dtlist_delete(dtlist);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return 0;
}

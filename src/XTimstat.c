/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2015 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
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

      Timstat    timmin          Time minimum
      Timstat    timmax          Time maximum
      Timstat    timsum          Time sum
      Timstat    timmean         Time mean
      Timstat    timavg          Time average
      Timstat    timvar          Time variance
      Timstat    timvar1         Time variance [Divisor is (n-1)]
      Timstat    timstd          Time standard deviation
      Timstat    timstd1         Time standard deviation [Divisor is (n-1)]
      Hourstat   hourmin         Hourly minimum
      Hourstat   hourmax         Hourly maximum
      Hourstat   hoursum         Hourly sum
      Hourstat   hourmean        Hourly mean
      Hourstat   houravg         Hourly average
      Hourstat   hourvar         Hourly variance
      Hourstat   hourvar1        Hourly variance [Divisor is (n-1)]
      Hourstat   hourstd         Hourly standard deviation
      Hourstat   hourstd1        Hourly standard deviation [Divisor is (n-1)]
      Daystat    daymin          Daily minimum
      Daystat    daymax          Daily maximum
      Daystat    daysum          Daily sum
      Daystat    daymean         Daily mean
      Daystat    dayavg          Daily average
      Daystat    dayvar          Daily variance
      Daystat    dayvar1         Daily variance [Divisor is (n-1)]
      Daystat    daystd          Daily standard deviation
      Daystat    daystd1         Daily standard deviation [Divisor is (n-1)]
      Monstat    monmin          Monthly minimum
      Monstat    monmax          Monthly maximum
      Monstat    monsum          Monthly sum
      Monstat    monmean         Monthly mean
      Monstat    monavg          Monthly average
      Monstat    monvar          Monthly variance
      Monstat    monvar1         Monthly variance [Divisor is (n-1)]
      Monstat    monstd          Monthly standard deviation
      Monstat    monstd1         Monthly standard deviation [Divisor is (n-1)]
      Yearstat   yearmin         Yearly minimum
      Yearstat   yearmax         Yearly maximum
      Yearstat   yearsum         Yearly sum
      Yearstat   yearmean        Yearly mean
      Yearstat   yearavg         Yearly average
      Yearstat   yearvar         Yearly variance
      Yearstat   yearvar1        Yearly variance [Divisor is (n-1)]
      Yearstat   yearstd         Yearly standard deviation
      Yearstat   yearstd1        Yearly standard deviation [Divisor is (n-1)]
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *XTimstat(void *argument)
{
  int timestat_date = TIMESTAT_MEAN;
  int gridsize;
  int vdate = 0, vtime = 0;
  int vdate0 = 0, vtime0 = 0;
  int nrecs;
  int varID, levelID, recID;
  long nsets;
  int i;
  int streamID3 = -1;
  int vlistID3, taxisID3 = -1;
  int nmiss;
  int nlevels;
  int lvfrac = FALSE;
  int nwpv; // number of words per value; real:1  complex:2
  char indate1[DATE_LEN+1], indate2[DATE_LEN+1];
  double vfrac = 1;
  double missval;

  cdoInitialize(argument);

  cdoOperatorAdd("xtimmin",    func_min,  DATE_LEN, NULL);
  cdoOperatorAdd("xtimmax",    func_max,  DATE_LEN, NULL);
  cdoOperatorAdd("xtimsum",    func_sum,  DATE_LEN, NULL);
  cdoOperatorAdd("xtimmean",   func_mean, DATE_LEN, NULL);
  cdoOperatorAdd("xtimavg",    func_avg,  DATE_LEN, NULL);
  cdoOperatorAdd("xtimvar",    func_var,  DATE_LEN, NULL);
  cdoOperatorAdd("xtimvar1",   func_var1, DATE_LEN, NULL);
  cdoOperatorAdd("xtimstd",    func_std,  DATE_LEN, NULL);
  cdoOperatorAdd("xtimstd1",   func_std1, DATE_LEN, NULL);
  cdoOperatorAdd("xmonmin",    func_min,   8, NULL);
  cdoOperatorAdd("xmonmax",    func_max,   8, NULL);
  cdoOperatorAdd("xmonsum",    func_sum,   8, NULL);
  cdoOperatorAdd("xmonmean",   func_mean,  8, NULL);
  cdoOperatorAdd("xmonavg",    func_avg,   8, NULL);
  cdoOperatorAdd("xmonvar",    func_var,   8, NULL);
  cdoOperatorAdd("xmonvar1",   func_var1,  8, NULL);
  cdoOperatorAdd("xmonstd",    func_std,   8, NULL);
  cdoOperatorAdd("xmonstd1",   func_std1,  8, NULL);

  int operatorID = cdoOperatorID();
  int operfunc   = cdoOperatorF1(operatorID);
  int comparelen = cdoOperatorF2(operatorID);

  int lmean   = operfunc == func_mean || operfunc == func_avg;
  int lstd    = operfunc == func_std || operfunc == func_std1;
  int lvarstd = operfunc == func_std || operfunc == func_var || operfunc == func_std1 || operfunc == func_var1;
  double divisor = operfunc == func_std1 || operfunc == func_var1;

  if ( operfunc == func_mean )
    {
      int oargc = operatorArgc();
      char **oargv = operatorArgv();

      if ( oargc == 1 )
	{
	  lvfrac = TRUE;
	  vfrac = atof(oargv[0]);
	  if ( cdoVerbose ) cdoPrint("Set vfrac to %g", vfrac);
	  if ( vfrac < 0 || vfrac > 1 ) cdoAbort("vfrac out of range!");
	}
      else if ( oargc > 1 )
	cdoAbort("Too many arguments!");
    }

  int cmplen = DATE_LEN - comparelen;

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  if ( cmplen == 0 ) vlistDefNtsteps(vlistID2, 1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  if ( taxisInqType(taxisID2) == TAXIS_FORECAST ) taxisDefType(taxisID2, TAXIS_RELATIVE);
  vlistDefTaxis(vlistID2, taxisID2);

  int nvars = vlistNvars(vlistID1);

  if ( cmplen == 0 && CDO_Reduce_Dim )
    for ( varID = 0; varID < nvars; ++varID )
      vlistDefVarTsteptype(vlistID2, varID, TSTEP_CONSTANT);

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  if ( cdoDiag )
    {
      char filename[8192];

      strcpy(filename, cdoOperatorName(operatorID));
      strcat(filename, "_");
      strcat(filename, cdoStreamName(1)->args);
      argument_t *fileargument = file_argument_new(filename);
      streamID3 = streamOpenWrite(fileargument, cdoFiletype());
      file_argument_free(fileargument);

      vlistID3 = vlistDuplicate(vlistID1);

      for ( varID = 0; varID < nvars; ++varID )
	{
	  vlistDefVarDatatype(vlistID3, varID, DATATYPE_INT32);
	  vlistDefVarMissval(vlistID3, varID, -1);
	  vlistDefVarUnits(vlistID3, varID, "");
	  vlistDefVarAddoffset(vlistID3, varID, 0);
	  vlistDefVarScalefactor(vlistID3, varID, 1);
	}

      taxisID3 = taxisDuplicate(taxisID1);
      vlistDefTaxis(vlistID3, taxisID3);

      streamDefVlist(streamID3, vlistID3);
    }

  dtlist_type *dtlist = dtlist_new();
  dtlist_set_stat(dtlist, timestat_date);
  dtlist_set_calendar(dtlist, taxisInqCalendar(taxisID1));

  gridsize = vlistGridsizeMax(vlistID1);
  if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;

  field_t **input_vars = field_malloc(vlistID1, FIELD_PTR);
  field_t **vars1 = field_malloc(vlistID1, FIELD_PTR);
  field_t **samp1 = field_malloc(vlistID1, FIELD_NONE);
  field_t **vars2 = NULL;
  if ( lvarstd ) vars2 = field_malloc(vlistID1, FIELD_PTR);

  int tsID  = 0;
  int otsID = 0;
  while ( TRUE )
    {
      nsets = 0;
      while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
	{
	  dtlist_taxisInqTimestep(dtlist, taxisID1, nsets);
	  vdate = dtlist_get_vdate(dtlist, nsets);
	  vtime = dtlist_get_vtime(dtlist, nsets);

	  if ( nsets == 0 ) SET_DATE(indate2, vdate, vtime);
	  SET_DATE(indate1, vdate, vtime);

	  if ( DATE_IS_NEQ(indate1, indate2, cmplen) ) break;

	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      streamInqRecord(streamID1, &varID, &levelID);
              streamReadRecord(streamID1, input_vars[varID][levelID].ptr, &nmiss);
              input_vars[varID][levelID].nmiss = nmiss;
            }

          for ( varID = 0; varID < nvars; varID++ )
            {
              nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
              for ( levelID = 0; levelID < nlevels; levelID++ )
                {
                  nwpv     = vars1[varID][levelID].nwpv;
                  gridsize = gridInqSize(vars1[varID][levelID].grid);

                  if ( nsets == 0 )
                    {
                      memcpy(vars1[varID][levelID].ptr, input_vars[varID][levelID].ptr, nwpv*gridsize*sizeof(double));
                      nmiss = input_vars[varID][levelID].nmiss;
                      vars1[varID][levelID].nmiss = nmiss;
                      if ( nmiss > 0 || samp1[varID][levelID].ptr )
                        {
                          if ( samp1[varID][levelID].ptr == NULL )
                            samp1[varID][levelID].ptr = (double*) Malloc(nwpv*gridsize*sizeof(double));

                          for ( i = 0; i < nwpv*gridsize; i++ )
                            if ( DBL_IS_EQUAL(vars1[varID][levelID].ptr[i], vars1[varID][levelID].missval) )
                              samp1[varID][levelID].ptr[i] = 0;
                            else
                              samp1[varID][levelID].ptr[i] = 1;
                        }
                    }
                  else
                    {
                      nmiss = input_vars[varID][levelID].nmiss;
                      if ( nmiss > 0 || samp1[varID][levelID].ptr )
                        {
                          if ( samp1[varID][levelID].ptr == NULL )
                            {
                              samp1[varID][levelID].ptr = (double*) Malloc(nwpv*gridsize*sizeof(double));
                              for ( i = 0; i < nwpv*gridsize; i++ )
                                samp1[varID][levelID].ptr[i] = nsets;
                            }
                          
                          for ( i = 0; i < nwpv*gridsize; i++ )
                            if ( !DBL_IS_EQUAL(input_vars[varID][levelID].ptr[i], vars1[varID][levelID].missval) )
                              samp1[varID][levelID].ptr[i]++;
                        }
                  
                      if ( lvarstd )
                        {
                          farsumq(&vars2[varID][levelID], input_vars[varID][levelID]);
                          farsum(&vars1[varID][levelID], input_vars[varID][levelID]);
                        }
                      else
                        {
                          farfun(&vars1[varID][levelID], input_vars[varID][levelID], operfunc);
                        }
                    }
		}
	    }

	  if ( nsets == 0 && lvarstd )
	    for ( varID = 0; varID < nvars; varID++ )
	      {
		if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
		nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
		for ( levelID = 0; levelID < nlevels; levelID++ )
		  farmoq(&vars2[varID][levelID], vars1[varID][levelID]);
	      }

	  vdate0 = vdate;
	  vtime0 = vtime;
	  nsets++;
	  tsID++;
	}

      if ( nrecs == 0 && nsets == 0 ) break;

      if ( lmean )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	    nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevels; levelID++ )
	      {
		if ( samp1[varID][levelID].ptr == NULL )
		  farcdiv(&vars1[varID][levelID], (double)nsets);
		else
                  {
                    //   farround(&samp1[varID][levelID]); not necessary
                    fardiv(&vars1[varID][levelID], samp1[varID][levelID]);
                  }
              }
	  }
      else if ( lvarstd )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	    nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevels; levelID++ )
	      {
		if ( samp1[varID][levelID].ptr == NULL )
		  {
		    if ( lstd )
		      farcstd(&vars1[varID][levelID], vars2[varID][levelID], nsets, divisor);
		    else
		      farcvar(&vars1[varID][levelID], vars2[varID][levelID], nsets, divisor);
		  }
		else
		  {
		    if ( lstd )
		      farstd(&vars1[varID][levelID], vars2[varID][levelID], samp1[varID][levelID], divisor);
		    else
		      farvar(&vars1[varID][levelID], vars2[varID][levelID], samp1[varID][levelID], divisor);
		  }
	      }
	  }

      if ( cdoVerbose )
	{
	  char vdatestr[32], vtimestr[32];
	  date2str(vdate0, vdatestr, sizeof(vdatestr));
	  time2str(vtime0, vtimestr, sizeof(vtimestr));
	  cdoPrint("%s %s  vfrac = %g, nsets = %d", vdatestr, vtimestr, vfrac, nsets);
	}

      if ( lvfrac && operfunc == func_mean )
	for ( varID = 0; varID < nvars; varID++ )
	  {
	    if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;
	    nwpv     = vars1[varID][levelID].nwpv;
	    gridsize = gridInqSize(vars1[varID][levelID].grid);
	    nlevels  = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	    for ( levelID = 0; levelID < nlevels; levelID++ )
	      {
		missval = vars1[varID][levelID].missval;
		if ( samp1[varID][levelID].ptr )
		  {
		    int irun = 0;
		    for ( i = 0; i < nwpv*gridsize; ++i )
		      {
			if ( (samp1[varID][levelID].ptr[i] / nsets) < vfrac )
			  {
			    vars1[varID][levelID].ptr[i] = missval;
			    irun++;
			  }
		      }

		    if ( irun )
		      {
			nmiss = 0;
			for ( i = 0; i < nwpv*gridsize; ++i )
			  if ( DBL_IS_EQUAL(vars1[varID][levelID].ptr[i], missval) ) nmiss++;
			vars1[varID][levelID].nmiss = nmiss;
		      }
		  }
	      }
	  }

      dtlist_stat_taxisDefTimestep(dtlist, taxisID2, nsets);
      streamDefTimestep(streamID2, otsID);

      if ( cdoDiag )
	{
	  dtlist_stat_taxisDefTimestep(dtlist, taxisID3, nsets);
	  streamDefTimestep(streamID3, otsID);
	}

      for ( varID = 0; varID < nvars; ++varID )
	{
	  if ( otsID && vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT ) continue;

          nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
          for ( levelID = 0; levelID < nlevels; ++ levelID )
            {
              streamDefRecord(streamID2, varID, levelID);
              streamWriteRecord(streamID2, vars1[varID][levelID].ptr,  vars1[varID][levelID].nmiss);
              if ( cdoDiag )
                {
                  if ( samp1[varID][levelID].ptr )
                    {
                      streamDefRecord(streamID3, varID, levelID);
                      streamWriteRecord(streamID3, samp1[varID][levelID].ptr, 0);
                    }
                }
            }
	}

      if ( nrecs == 0 ) break;
      otsID++;
    }


  field_free(input_vars, vlistID1);
  field_free(vars1, vlistID1);
  field_free(samp1, vlistID1);
  if ( lvarstd ) field_free(vars2, vlistID1);

  dtlist_delete(dtlist);

  if ( cdoDiag ) streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return 0;
}
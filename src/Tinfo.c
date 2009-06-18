/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2007-2009 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Tinfo      tinfo           Time information
*/


#include <stdio.h>
#include <string.h>
#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "dtypes.h"


#define MAX_GAPS  64
#define MAX_NTSM 128

enum {TU_SECONDS=0, TU_MINUTES, TU_HOURS, TU_DAYS, TU_MONTHS, TU_YEARS};
char *tunits[] = {"second", "minute", "hour", "day", "month", "year"};
int   iunits[] = {1, 60, 3600, 86400, 1, 12};

static
void printTunit(int unit)
{
  if ( unit == TUNIT_YEAR )
    fprintf(stdout, "  Units = years");
  else if ( unit == TUNIT_MONTH )
    fprintf(stdout, "  Units = months");
  else if ( unit == TUNIT_DAY )
    fprintf(stdout, "  Units = days");
  else if ( unit == TUNIT_HOUR )
    fprintf(stdout, "  Units = hours");
  else if ( unit == TUNIT_MINUTE )
    fprintf(stdout, "  Units = minutes");
  else if ( unit == TUNIT_SECOND )
    fprintf(stdout, "  Units = seconds");
  else
    fprintf(stdout, "  Units = unknown");
}


static
void printCalendar(int calendar)
{
  if      ( calendar == CALENDAR_STANDARD )
    fprintf(stdout, "  Calendar = STANDARD");
  else if ( calendar == CALENDAR_PROLEPTIC )
    fprintf(stdout, "  Calendar = PROLEPTIC");
  else if ( calendar == CALENDAR_360DAYS )
    fprintf(stdout, "  Calendar = 360DAYS");
  else if ( calendar == CALENDAR_365DAYS )
    fprintf(stdout, "  Calendar = 365DAYS");
  else if ( calendar == CALENDAR_366DAYS )
    fprintf(stdout, "  Calendar = 366DAYS");
  else
    fprintf(stdout, "  Calendar = unknown");
}


void getTimeInc(int lperiod, int deltam, int deltay, int *incperiod, int *incunit)
{
  *incperiod = 0;
  *incunit   = 0;

  if ( lperiod/60 > 0 && lperiod/60 < 60 )
    {
      *incperiod = lperiod/60;
      *incunit = TU_MINUTES;
    }
  else if ( lperiod/3600 > 0 && lperiod/3600 < 24 )
    {
      *incperiod = lperiod/3600;
      *incunit = TU_HOURS;
    }
  else if ( lperiod/(3600*24) > 0 && lperiod/(3600*24) < 32 )
    {
      *incperiod = lperiod/(3600*24);
      *incunit = TU_DAYS;
      if ( *incperiod > 27 && deltam == 1 )
	{
	  *incperiod = 1;
	  *incunit = TU_MONTHS;
	}
    }
  else if ( lperiod/(3600*24*30) > 0 && lperiod/(3600*24*30) < 12 )
    {
      *incperiod = deltam;
      *incunit = TU_MONTHS;
    }
  else if ( lperiod/(3600*24*30*12) > 0 )
    {
      *incperiod = deltay;
      *incunit = TU_YEARS;
    }
  else
    {
      *incperiod = lperiod;
      *incunit = TU_SECONDS;
    }
}

static
void printBounds(int taxisID, int calendar)
{
  int vdate0, vdate1;
  int vtime0, vtime1;
  int year0, month0, day0, hour0, minute0, second0;
  int year1, month1, day1, hour1, minute1, second1;
  INT64 lperiod;
  int incperiod = 0, incunit = 0;
  juldate_t juldate1, juldate0;
  double jdelta;
  int deltam, deltay;
  int i, len;

  taxisInqVdateBounds(taxisID, &vdate0, &vdate1);
  taxisInqVtimeBounds(taxisID, &vtime0, &vtime1);

  decode_date(vdate0, &year0, &month0, &day0);
  decode_time(vtime0, &hour0, &minute0, &second0);

  decode_date(vdate1, &year1, &month1, &day1);
  decode_time(vtime1, &hour1, &minute1, &second1);

  fprintf(stdout, "%5.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d",
	  year0, month0, day0, hour0, minute0, second0);
  fprintf(stdout, "%5.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d",
	  year1, month1, day1, hour1, minute1, second1);

  juldate0  = juldate_encode(calendar, vdate0, vtime0);
  juldate1  = juldate_encode(calendar, vdate1, vtime1);
  jdelta    = juldate_to_seconds(juldate_sub(juldate1, juldate0));
  lperiod   = (INT64)(jdelta+0.5);
  incperiod = (int) lperiod;

  deltay = year1-year0;
  deltam = deltay*12 + (month1-month0);

  getTimeInc(lperiod, deltam, deltay, &incperiod, &incunit);
  
  /* fprintf(stdout, "  %g  %g  %g  %d", jdelta, jdelta/3600, fmod(jdelta,3600), incperiod%3600);*/
  len = fprintf(stdout, " %3d %s%s", incperiod, tunits[incunit], incperiod>1?"s":"");
  for ( i = 0; i < 11-len; ++i ) fprintf(stdout, " ");
}


void *Tinfo(void *argument)
{
  int vdate0 = 0, vtime0 = 0;
  int vdate, vtime;
  int nrecs, ntsteps;
  int tsID, ntimeout;
  int taxisID;
  int streamID;
  int vlistID;
  int year0, month0, day0, hour0, minute0, second0;
  int year, month, day, hour, minute, second;
  int calendar, unit;
  int incperiod0 = 0, incunit0 = 0;
  int incperiod1 = 0, incunit1 = 0;
  INT64 lperiod;
  int incperiod = 0, incunit = 0;
  int its = 0, igap;
  int ngaps = 0;
  int ntsm[MAX_GAPS];
  int rangetsm[MAX_GAPS][2];
  int vdatem[MAX_GAPS][MAX_NTSM];
  int vtimem[MAX_GAPS][MAX_NTSM];
  juldate_t juldate, juldate0;
  double jdelta;
  int i, len;
	  

  cdoInitialize(argument);

  streamID = streamOpenRead(cdoStreamName(0));
  if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", cdoStreamName(0));

  vlistID = streamInqVlist(streamID);

  fprintf(stdout, "\n");

  taxisID = vlistInqTaxis(vlistID);
  ntsteps = vlistNtsteps(vlistID);

  if ( ntsteps != 0 )
    {
      if ( ntsteps == CDI_UNDEFID )
	fprintf(stdout, "   Time axis :  unlimited steps\n");
      else
	fprintf(stdout, "   Time axis :  %d step%s\n", ntsteps, ntsteps == 1 ? "" : "s");

      if ( taxisID != CDI_UNDEFID )
	{
	  if ( taxisInqType(taxisID) == TAXIS_RELATIVE )
	    {
	      vdate = taxisInqRdate(taxisID);
	      vtime = taxisInqRtime(taxisID);
	      
	      decode_date(vdate, &year, &month, &day);
	      decode_time(vtime, &hour, &minute, &second);

	      fprintf(stdout, "     RefTime = %4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d",
		      year, month, day, hour, minute, second);
		      
	      unit = taxisInqTunit(taxisID);
	      if ( unit != CDI_UNDEFID ) printTunit(unit);
	      
	      calendar = taxisInqCalendar(taxisID);
	      if ( calendar != CDI_UNDEFID ) printCalendar(calendar);

	      fprintf(stdout, "\n");
	    }
	}

      calendar = taxisInqCalendar(taxisID);

      if ( taxisHasBounds(taxisID) ) 
	fprintf(stdout, "\nTimestep YYYY-MM-DD hh:mm:ss   Inrement YYYY-MM-DD hh:mm:ss YYYY-MM-DD hh:mm:ss  Difference\n");
      else
	fprintf(stdout, "\nTimestep YYYY-MM-DD hh:mm:ss   Inrement\n");

      tsID = 0;
      while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	{  
	  vdate = taxisInqVdate(taxisID);
	  vtime = taxisInqVtime(taxisID);
	  
	  decode_date(vdate, &year, &month, &day);
	  decode_time(vtime, &hour, &minute, &second);

	  fprintf(stdout, "%6d  %5.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d",
		  tsID+1, year, month, day, hour, minute, second);
	  if ( tsID )
	    {
	      int deltam, deltay;

	      decode_date(vdate0, &year0, &month0, &day0);
	      decode_time(vtime0, &hour0, &minute0, &second0);

	      juldate0  = juldate_encode(calendar, vdate0, vtime0);
	      juldate   = juldate_encode(calendar, vdate, vtime);
	      jdelta    = juldate_to_seconds(juldate_sub(juldate, juldate0));
	      lperiod   = (INT64)(jdelta+0.5);
	      incperiod = (int) lperiod;

	      deltay = year-year0;
	      deltam = deltay*12 + (month-month0);

	      getTimeInc(lperiod, deltam, deltay, &incperiod, &incunit);

	      /* fprintf(stdout, "  %g  %g  %g  %d", jdelta, jdelta/3600, fmod(jdelta,3600), incperiod%3600);*/
	      len = fprintf(stdout, " %3d %s%s", incperiod, tunits[incunit], incperiod>1?"s":"");
	      for ( i = 0; i < 11-len; ++i ) fprintf(stdout, " ");
	    }
	  else
	    {
	      fprintf(stdout, "   --------");
	    }

	  if ( taxisHasBounds(taxisID) ) printBounds(taxisID, calendar);

	  if ( tsID > 1 )
	    {
	      if ( incperiod != incperiod0 || incunit != incunit0 )
		{
		  int ndate, ntime;
		  int ijulinc = incperiod0 * iunits[incunit0];
		  /*
		  printf("%d %s changed to %d %s\n", 
			 incperiod0, tunits[incunit0], incperiod, tunits[incunit]);
		  */
		  if ( ngaps < MAX_GAPS )
		    {
		      rangetsm[ngaps][0] = tsID;
		      rangetsm[ngaps][1] = tsID+1;

		      if ( incunit0 == TU_MONTHS || incunit0 == TU_YEARS )
			{
			  its = 0;
			  ndate = vdate0;
			  while ( TRUE )
			    {
			      decode_date(ndate, &year, &month, &day);
				  
			      month += ijulinc;
				  
			      while ( month > 12 ) { month -= 12; year++; }
			      while ( month <  1 ) { month += 12; year--; }
			      
			      if ( day0 == 31 )
				day = days_per_month(calendar, year, month);

			      ndate = encode_date(year, month, day);
			      ntime = vtime0;
			      if ( ndate >= vdate ) break;
			      /* printf("\n1 %d %d\n", ndate, ntime); */
			      if ( its < MAX_NTSM )
				{
				  vdatem[ngaps][its] = ndate;
				  vtimem[ngaps][its] = ntime;
				}
			      its++;
			    }
			}
		      else
			{
			  its = 0;
			  juldate0 = juldate_add_seconds(ijulinc, juldate0);
			  while ( juldate_to_seconds(juldate0) < juldate_to_seconds(juldate) )
			    {
			      juldate_decode(calendar, juldate0, &ndate, &ntime);
			      juldate0 = juldate_add_seconds(ijulinc, juldate0);
			      if ( its < MAX_NTSM )
				{
				  vdatem[ngaps][its] = ndate;
				  vtimem[ngaps][its] = ntime;
				}
			      its++;
			    }
			}			
		      ntsm[ngaps] = its;
		    }
		  ngaps++;
		  if ( cdoVerbose )
		    fprintf(stdout, "  <--- Gap %d, missing %d timestep%s",
			    ngaps, its, its>1?"s":"");
		}
	    }

	  if ( tsID )
	    {
	      if ( tsID == 1 )
		{
		  incperiod0 = incperiod;
		  incunit0   = incunit;
		}

	      incperiod1 = incperiod;
	      incunit1   = incunit;
	    }
	  fprintf(stdout, "\n");

	  vdate0 = vdate;
	  vtime0 = vtime;

	  tsID++;
	}
    }

  if ( cdoVerbose && ngaps )
    {
      fprintf(stdout, "\nFound potentially %d gap%s in the time series", ngaps, ngaps>1?"s":"");
      if ( ngaps >= MAX_GAPS )
	{
	  ngaps = MAX_GAPS;
	  fprintf(stdout, ", here are the first %d", ngaps);
	}
      fprintf(stdout, ":\n");
      for ( igap = 0; igap < ngaps; ++igap )
	{
	  fprintf(stdout, "  Gap %d between timestep %d and %d. Missing %d timestep%s",
		  igap+1, rangetsm[igap][0], rangetsm[igap][1], ntsm[igap], ntsm[igap]>1?"s":"");
	  if ( ntsm[igap] >= MAX_NTSM )
	    {
	      ntsm[igap] = MAX_NTSM;
	      fprintf(stdout, ", here are the first %d", ntsm[igap]);
	    }
	  fprintf(stdout, ":\n");
	  
	  ntimeout = 0;
	  for ( its = 0; its < ntsm[igap]; ++its )
	    {
	      if ( ntimeout == 4 )
		{
		  ntimeout = 0;
		  fprintf(stdout, "\n");
		}

	      vdate = vdatem[igap][its];
	      vtime = vtimem[igap][its];

	      decode_date(vdate, &year, &month, &day);
	      decode_time(vtime, &hour, &minute, &second);

	      fprintf(stdout, "   %4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d",
		      year, month, day, hour, minute, second);
	      ntimeout++;
	      tsID++;
	    }
	  fprintf(stdout, "\n");
	}
    }

  streamClose(streamID);

  cdoFinish();

  return (0);
}

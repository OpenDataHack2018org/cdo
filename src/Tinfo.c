/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2007-2008 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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
  if ( calendar == CALENDAR_STANDARD )
    fprintf(stdout, "  Calendar = STANDARD");
  else if ( calendar == CALENDAR_NONE )
    fprintf(stdout, "  Calendar = NONE");
  else if ( calendar == CALENDAR_360DAYS )
    fprintf(stdout, "  Calendar = 360DAYS");
  else if ( calendar == CALENDAR_365DAYS )
    fprintf(stdout, "  Calendar = 365DAYS");
  else if ( calendar == CALENDAR_366DAYS )
    fprintf(stdout, "  Calendar = 366DAYS");
  else
    fprintf(stdout, "  Calendar = unknown");
}


static
void getTimeInc(int lperiod, int deltam, int deltay, int *incperiod, int *incunit)
{
  *incperiod = 0;
  *incunit   = 0;

  if ( lperiod/60 > 0 && lperiod/60 < 60 )
    {
      *incperiod = lperiod/60;
      *incunit = 1;
    }
  else if ( lperiod/3600 > 0 && lperiod/3600 < 24 )
    {
      *incperiod = lperiod/3600;
      *incunit = 2;
    }
  else if ( lperiod/(3600*24) > 0 && lperiod/(3600*24) < 32 )
    {
      *incperiod = lperiod/(3600*24);
      *incunit = 3;
      if ( *incperiod > 27 && deltam == 1 )
	{
	  *incperiod = 1;
	  *incunit = 4;
	}
    }
  else if ( lperiod/(3600*24*30) > 0 && lperiod/(3600*24*30) < 12 )
    {
      *incperiod = deltam;
      *incunit = 4;
    }
  else if ( lperiod/(3600*24*30*12) > 0 )
    {
      *incperiod = deltay;
      *incunit = 5;
    }
}

static
void printBounds(int taxisID, int dpy)
{
  int vdate0, vdate1;
  int vtime0, vtime1;
  int year0, month0, day0, hour0, minute0;
  int year1, month1, day1, hour1, minute1;
  INT64 lperiod;
  int incperiod = 0, incunit = 0;
  double julval1 = 0, julval0 = 0, jdelta;
  int deltam, deltay;
  int i, len;

  taxisInqVdateBounds(taxisID, &vdate0, &vdate1);
  taxisInqVtimeBounds(taxisID, &vtime0, &vtime1);

  decode_date(vdate0, &year0, &month0, &day0);
  decode_time(vtime0, &hour0, &minute0);

  decode_date(vdate1, &year1, &month1, &day1);
  decode_time(vtime1, &hour1, &minute1);

  fprintf(stdout, "%6.4d-%2.2d-%2.2d %2.2d:%2.2d", year0, month0, day0, hour0, minute0);
  fprintf(stdout, "%6.4d-%2.2d-%2.2d %2.2d:%2.2d", year1, month1, day1, hour1, minute1);

  julval0 = encode_julval(dpy, vdate0, vtime0);
  julval1 = encode_julval(dpy, vdate1, vtime1);
  jdelta  = julval1 - julval0;
  lperiod = (INT64)(jdelta+0.5);
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
  int year0, month0, day0, hour0, minute0;
  int year, month, day, hour, minute;
  int calendar, unit, dpy;
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
  double julval = 0, julval0 = 0, jdelta;
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
	      decode_time(vtime, &hour, &minute);

	      fprintf(stdout, "     RefTime = %4.4d-%2.2d-%2.2d %2.2d:%2.2d",
		      year, month, day, hour, minute);
		      
	      unit = taxisInqTunit(taxisID);
	      if ( unit != CDI_UNDEFID ) printTunit(unit);
	      
	      calendar = taxisInqCalendar(taxisID);
	      if ( calendar != CDI_UNDEFID ) printCalendar(calendar);

	      fprintf(stdout, "\n");
	    }
	}

      calendar = taxisInqCalendar(taxisID);

      dpy = calendar_dpy(calendar);

      if ( taxisHasBounds(taxisID) ) 
	fprintf(stdout, "\nTimestep  YYYY-MM-DD hh:mm   Inrement  YYYY-MM-DD hh:mm  YYYY-MM-DD hh:mm  Difference\n");
      else
	fprintf(stdout, "\nTimestep  YYYY-MM-DD hh:mm   Inrement\n");

      tsID = 0;
      while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	{  
	  vdate = taxisInqVdate(taxisID);
	  vtime = taxisInqVtime(taxisID);
	  
	  decode_date(vdate, &year, &month, &day);
	  decode_time(vtime, &hour, &minute);

	  fprintf(stdout, "%6d  %6.4d-%2.2d-%2.2d %2.2d:%2.2d", tsID+1, year, month, day, hour, minute);
	  if ( tsID )
	    {
	      int deltam, deltay;

	      decode_date(vdate0, &year0, &month0, &day0);
	      decode_time(vtime0, &hour0, &minute0);

	      julval0 = encode_julval(dpy, vdate0, vtime0);
	      julval  = encode_julval(dpy, vdate, vtime);
	      jdelta  = julval - julval0;
	      lperiod = (INT64)(jdelta+0.5);
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

	  if ( taxisHasBounds(taxisID) ) printBounds(taxisID, dpy);

	  if ( tsID > 1 )
	    {
	      if ( incperiod != incperiod0 || incunit != incunit0 )
		{
		  int ndate, ntime;
		  int ijulinc = incperiod0 * iunits[incunit0];
		  if ( ngaps < MAX_GAPS )
		    {
		      rangetsm[ngaps][0] = tsID;
		      rangetsm[ngaps][1] = tsID+1;

		      if ( incunit0 == 4 || incunit0 == 5 )
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
				day = days_per_month(dpy, year, month);

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
			  julval0 += ijulinc;
			  while ( julval0 < julval )
			    {
			      decode_julval(dpy, julval0, &ndate, &ntime);
			      julval0 += ijulinc;
			      /* printf("\n2 %d %d %g %d\n", ndate, ntime, julval0, ijulinc); */
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
      fprintf(stdout, "\nFound potentially %d gaps in the time series", ngaps);
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
	      decode_time(vtime, &hour, &minute);

	      fprintf(stdout, "   %4.4d-%2.2d-%2.2d %2.2d:%2.2d", year, month, day, hour, minute);
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

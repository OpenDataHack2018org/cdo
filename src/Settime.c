/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2013 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Settime    setdate         Set date
      Settime    settime         Set time
      Settime    setday          Set day
      Settime    setmon          Set month
      Settime    setyear         Set year
      Settime    settunits       Set time units
      Settime    settaxis        Set time axis
      Settime    setreftime      Set reference time
      Settime    setcalendar     Set calendar
      Settime    shifttime       Shift timesteps
*/

#include <ctype.h>  /* isdigit */

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


int get_tunits(const char *unit, int *incperiod, int *incunit, int *tunit)
{
  size_t len;
	
  len = strlen(unit);
  
  if      ( memcmp(unit, "seconds", len) == 0 ) { *incunit =     1; *tunit = TUNIT_SECOND;  }
  else if ( memcmp(unit, "minutes", len) == 0 ) { *incunit =    60; *tunit = TUNIT_MINUTE;  }
  else if ( memcmp(unit, "hours", len)   == 0 ) { *incunit =  3600; *tunit = TUNIT_HOUR;    }
  else if ( memcmp(unit, "3hours", len)  == 0 ) { *incunit = 10800; *tunit = TUNIT_3HOURS;  }
  else if ( memcmp(unit, "6hours", len)  == 0 ) { *incunit = 21600; *tunit = TUNIT_6HOURS;  }
  else if ( memcmp(unit, "12hours", len) == 0 ) { *incunit = 43200; *tunit = TUNIT_12HOURS; }
  else if ( memcmp(unit, "days", len)    == 0 ) { *incunit = 86400; *tunit = TUNIT_DAY;     }
  else if ( memcmp(unit, "months", len)  == 0 ) { *incunit =     1; *tunit = TUNIT_MONTH;   }
  else if ( memcmp(unit, "years", len)   == 0 ) { *incunit =    12; *tunit = TUNIT_YEAR;    }
  else cdoAbort("time unit >%s< unsupported", unit);

  if ( *tunit == TUNIT_HOUR )
    {
      if      ( *incperiod ==  3 ) { *incperiod = 1; *incunit = 10800; *tunit = TUNIT_3HOURS;  }
      else if ( *incperiod ==  6 ) { *incperiod = 1; *incunit = 21600; *tunit = TUNIT_6HOURS;  }
      else if ( *incperiod == 12 ) { *incperiod = 1; *incunit = 43200; *tunit = TUNIT_12HOURS; }
    }

  return (0);
}

static
void shifttime(int calendar, int tunit, int ijulinc, int *pdate, int *ptime)
{
  int year, month, day;
  int vdate = *pdate;
  int vtime = *ptime;
  juldate_t juldate;

  if ( tunit == TUNIT_MONTH || tunit == TUNIT_YEAR )
    {
      cdiDecodeDate(vdate, &year, &month, &day);
	      
      month += ijulinc;

      while ( month > 12 ) { month -= 12; year++; }
      while ( month <  1 ) { month += 12; year--; }

      vdate = cdiEncodeDate(year, month, day);

      *pdate = vdate;
    }
  else
    {
      juldate = juldate_encode(calendar, vdate, vtime);
      juldate = juldate_add_seconds(ijulinc, juldate);
      juldate_decode(calendar, juldate, &vdate, &vtime);

      *pdate = vdate;
      *ptime = vtime;

      if ( cdoVerbose )
	cdoPrint("juldate, ijulinc, vdate, vtime: %g %d %d %d",
		 juldate_to_seconds(juldate), ijulinc, vdate, vtime);
    }
}


void *Settime(void *argument)
{
  int SETYEAR, SETMON, SETDAY, SETDATE, SETTIME, SETTUNITS;
  int SETTAXIS, SETREFTIME, SETCALENDAR, SHIFTTIME;
  int operatorID;
  int streamID1, streamID2 = CDI_UNDEFID;
  int nrecs, newval = 0, ntsteps, nvars;
  int tsID1, recID, varID, levelID;
  int vlistID1, vlistID2;
  int vdate, vtime;
  int vdateb[2], vtimeb[2];
  int sdate = 0, stime = 0;
  int taxisID1, taxisID2 = CDI_UNDEFID;
  int nmiss;
  int gridsize;
  int tunit = TUNIT_DAY;
  int ijulinc = 0, incperiod = 0, incunit = 0;
  int year = 1, month = 1, day = 1, hour = 0, minute = 0, second = 0;
  int day0;
  int taxis_has_bounds, copy_timestep = FALSE;
  int calendar;
  int newcalendar = CALENDAR_STANDARD;
  // int nargs;
  const char *datestr, *timestr;
  char *rstr;
  juldate_t juldate;
  double *array = NULL;

  cdoInitialize(argument);

  SETYEAR     = cdoOperatorAdd("setyear",     0,  1, "year");
  SETMON      = cdoOperatorAdd("setmon",      0,  1, "month");
  SETDAY      = cdoOperatorAdd("setday",      0,  1, "day");
  SETDATE     = cdoOperatorAdd("setdate",     0,  1, "date (format: YYYY-MM-DD)");
  SETTIME     = cdoOperatorAdd("settime",     0,  1, "time (format: hh:mm:ss)");
  SETTUNITS   = cdoOperatorAdd("settunits",   0,  1, "time units (seconds, minutes, hours, days, months, years)");
  SETTAXIS    = cdoOperatorAdd("settaxis",    0, -2, "date,time<,increment> (format YYYY-MM-DD,hh:mm:ss)");
  SETREFTIME  = cdoOperatorAdd("setreftime",  0, -2, "date,time<,units> (format YYYY-MM-DD,hh:mm:ss)");
  SETCALENDAR = cdoOperatorAdd("setcalendar", 0,  1, "calendar (standard, proleptic, 360days, 365days, 366days)");
  SHIFTTIME   = cdoOperatorAdd("shifttime",   0,  1, "shift value");

  operatorID = cdoOperatorID();
  // nargs = cdoOperatorF2(operatorID);

  operatorInputArg(cdoOperatorEnter(operatorID));

  //  if ( operatorArgc()

  if ( operatorID == SETTAXIS || operatorID == SETREFTIME )
    {
      if ( operatorArgc() < 2 ) cdoAbort("Too few arguments!");

      datestr = operatorArgv()[0];
      timestr = operatorArgv()[1];

      if ( strchr(datestr, '-') )
	{
	  sscanf(datestr, "%d-%d-%d", &year, &month, &day);
	  sdate = cdiEncodeDate(year, month, day);
	}
      else
	{
	  sdate = (int)strtol(datestr, &rstr, 10);
	  if ( *rstr != 0 ) cdoAbort("Parameter string contains invalid characters: %s", datestr);
	}

      if ( strchr(timestr, ':') )
	{
	  sscanf(timestr, "%d:%d:%d", &hour, &minute, &second);
	  stime = cdiEncodeTime(hour, minute, second);
	}
      else
	{
	  stime = (int)strtol(timestr, &rstr, 10);
	  if ( *rstr != 0 ) cdoAbort("Parameter string contains invalid characters: %s", timestr);
	}

      if ( operatorArgc() == 3 )
	{
	  const char *timeunits = operatorArgv()[2];
	  incperiod = (int)strtol(timeunits, NULL, 10);;
	  while ( isdigit((int) *timeunits) ) timeunits++;

	  get_tunits(timeunits, &incperiod, &incunit, &tunit);
	}
      /* increment in seconds */
      ijulinc = incperiod * incunit;
    }
  else if ( operatorID == SETDATE )
    {
      if ( operatorArgc() < 1 ) cdoAbort("Too few arguments!");
      datestr = operatorArgv()[0];
      if ( strchr(datestr, '-') )
	{
	  sscanf(datestr, "%d-%d-%d", &year, &month, &day);
	  newval = cdiEncodeDate(year, month, day);
	}
      else
	{
	  newval = (int)strtol(datestr, &rstr, 10);
	  if ( *rstr != 0 ) cdoAbort("Parameter string contains invalid characters: %s", datestr);
	}
    }
  else if ( operatorID == SETTIME )
    {
      if ( operatorArgc() < 1 ) cdoAbort("Too few arguments!");
      timestr = operatorArgv()[0];

      if ( strchr(timestr, ':') )
	{
	  sscanf(timestr, "%d:%d:%d", &hour, &minute, &second);
	  newval = cdiEncodeTime(hour, minute, second);
	}
      else
	{
	  newval = (int)strtol(timestr, &rstr, 10);
	  if ( *rstr != 0 ) cdoAbort("Parameter string contains invalid characters: %s", timestr);
	}
    }
  else if ( operatorID == SHIFTTIME )
    {
      const char *timeunits = operatorArgv()[0];
      incperiod = (int)strtol(timeunits, NULL, 10);;
      if ( timeunits[0] == '-' || timeunits[0] == '+' ) timeunits++;
      while ( isdigit((int) *timeunits) ) timeunits++;

      get_tunits(timeunits, &incperiod, &incunit, &tunit);

      /* increment in seconds */
      ijulinc = incperiod * incunit;
    }
  else if ( operatorID == SETTUNITS )
    {
      int idum;
      const char *timeunits = operatorArgv()[0];
      incperiod = 0;
      get_tunits(timeunits, &incperiod, &idum, &tunit);
    }
  else if ( operatorID == SETCALENDAR )
    {
      size_t len;
      char *cname = operatorArgv()[0];
      len = strlen(cname);      
      if      ( memcmp(cname, "standard" , len) == 0 ) { newcalendar = CALENDAR_STANDARD;}
      else if ( memcmp(cname, "proleptic", len) == 0 ) { newcalendar = CALENDAR_PROLEPTIC;}
      else if ( memcmp(cname, "360days",   len) == 0 ) { newcalendar = CALENDAR_360DAYS;}
      else if ( memcmp(cname, "365days",   len) == 0 ) { newcalendar = CALENDAR_365DAYS;}
      else if ( memcmp(cname, "366days",   len) == 0 ) { newcalendar = CALENDAR_366DAYS;}
      else cdoAbort("Calendar >%s< unsupported!", cname);
    }
  else
    {
      newval = (int)strtol(operatorArgv()[0], &rstr, 10);
      if ( *rstr != 0 ) cdoAbort("Parameter string contains invalid characters: %s", operatorArgv()[0]);
    }

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxis_has_bounds = taxisHasBounds(taxisID1);
  ntsteps  = vlistNtsteps(vlistID1);
  nvars    = vlistNvars(vlistID1);

  if ( ntsteps == 1 )
    {
      for ( varID = 0; varID < nvars; ++varID )
	if ( vlistInqVarTsteptype(vlistID1, varID) != TSTEP_CONSTANT ) break;

      if ( varID == nvars ) ntsteps = 0;
    }

  if ( ntsteps == 0 )
    {
      for ( varID = 0; varID < nvars; ++varID )
	vlistDefVarTsteptype(vlistID2, varID, TSTEP_INSTANT);
    }

  calendar = taxisInqCalendar(taxisID1);

  if ( cdoVerbose ) cdoPrint("calendar = %d", calendar);

  if ( operatorID == SETREFTIME )
    {
      copy_timestep = TRUE;

      if ( taxisInqType(taxisID1) == TAXIS_ABSOLUTE )
	{
	  cdoPrint("Changing absolute to relative time axis!");

	  taxisID2 = taxisCreate(TAXIS_RELATIVE);
	}
      else
	taxisID2 = taxisDuplicate(taxisID1);

      if ( operatorArgc() != 3 ) tunit = taxisInqTunit(taxisID1);
      taxisDefTunit(taxisID2, tunit);
    }
  else if ( operatorID == SETTUNITS )
    {
      copy_timestep = TRUE;

      if ( taxisInqType(taxisID1) == TAXIS_ABSOLUTE )
	{
	  cdoPrint("Changing absolute to relative time axis!");

	  taxisID2 = taxisCreate(TAXIS_RELATIVE);
	  taxisDefTunit(taxisID2, tunit);
	}
      else
	taxisID2 = taxisDuplicate(taxisID1);
    }
  else if ( operatorID == SETCALENDAR )
    {
      copy_timestep = TRUE;
      /*
      if ( ((char *)argument)[0] == '-' )
	cdoAbort("This operator does not work with pipes!");
      */
      if ( taxisInqType(taxisID1) == TAXIS_ABSOLUTE )
	{/*
	  if ( cdoFiletype() != FILETYPE_NC )
	    cdoAbort("This operator does not work on an absolute time axis!");
	 */
	  cdoPrint("Changing absolute to relative time axis!");
	  taxisID2 = taxisCreate(TAXIS_RELATIVE);
	}
      else
	taxisID2 = taxisDuplicate(taxisID1);
    }
  else
    taxisID2 = taxisDuplicate(taxisID1);

  if ( operatorID == SETTAXIS )
    {
      taxisDefTunit(taxisID2, tunit);
      taxisDefRdate(taxisID2, sdate);
      taxisDefRtime(taxisID2, stime);
      juldate = juldate_encode(calendar, sdate, stime);
    }
  else if ( operatorID == SETTUNITS )
    {
      taxisDefTunit(taxisID2, tunit);
    }
  else if ( operatorID == SETCALENDAR )
    {
      taxisDefCalendar(taxisID2, newcalendar);
    }

  if ( operatorID != SHIFTTIME )
    if ( taxis_has_bounds && copy_timestep == FALSE )
      {
	cdoWarning("Time bounds unsupported by this operator, removed!");
	taxisDeleteBounds(taxisID2);
	taxis_has_bounds = FALSE;
      }

  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);
  if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
  array = (double *) malloc(gridsize*sizeof(double));

  tsID1 = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID1)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      if ( operatorID == SETTAXIS )
	{
	  if ( tunit == TUNIT_MONTH || tunit == TUNIT_YEAR )
	    {
	      vtime = stime;
	      if ( tsID1 == 0 )
		{
		  vdate = sdate;
		  cdiDecodeDate(vdate, &year, &month, &day0);
		}
	      else
		{	      
		  month += ijulinc;

		  while ( month > 12 ) { month -= 12; year++; }
		  while ( month <  1 ) { month += 12; year--; }

		  if ( day0 == 31 )
		    day = days_per_month(calendar, year, month);
		  else
		    day = day0;

		  vdate = cdiEncodeDate(year, month, day);
		}
	    }
	  else
	    {
	      juldate_decode(calendar, juldate, &vdate, &vtime);
	      juldate = juldate_add_seconds(ijulinc, juldate);
	    }
	}
      else if ( operatorID == SHIFTTIME )
	{
	  shifttime(calendar, tunit, ijulinc, &vdate, &vtime);
	  if ( taxis_has_bounds )
	    {
	      taxisInqVdateBounds(taxisID1, &vdateb[0], &vdateb[1]);
	      taxisInqVtimeBounds(taxisID1, &vtimeb[0], &vtimeb[1]);	      
	      shifttime(calendar, tunit, ijulinc, &vdateb[0], &vtimeb[0]);
	      shifttime(calendar, tunit, ijulinc, &vdateb[1], &vtimeb[1]);
	    }
	}
      else if ( operatorID == SETREFTIME || operatorID == SETCALENDAR || operatorID == SETTUNITS )
	{
	  ;
	}
      else
	{
	  cdiDecodeDate(vdate, &year, &month, &day);

	  if ( operatorID == SETYEAR ) year  = newval;
	  if ( operatorID == SETMON  ) month = newval;
	  if ( operatorID == SETMON && (month < 0 || month > 16) ) cdoAbort("parameter month=%d out of range!", month);
	  if ( operatorID == SETDAY  ) day   = newval;
	  if ( operatorID == SETDAY && (day < 0 || day > 31) ) cdoAbort("parameter day=%d %d out of range!", day);
      
	  vdate = cdiEncodeDate(year, month, day);

	  if ( operatorID == SETDATE  ) vdate = newval;
	  if ( operatorID == SETTIME  ) vtime = newval;
	}

      if ( copy_timestep )
	{
	  taxisCopyTimestep(taxisID2, taxisID1);
	  if ( operatorID == SETREFTIME )
	    {
	      taxisDefRdate(taxisID2, sdate);
	      taxisDefRtime(taxisID2, stime);
	    }
	}
      else
	{
	  int numavg = taxisInqNumavg(taxisID1);
	  taxisDefNumavg(taxisID2, numavg);

	  taxisDefVdate(taxisID2, vdate);
	  taxisDefVtime(taxisID2, vtime);
	  if ( taxis_has_bounds )
	    {
	      taxisDefVdateBounds(taxisID2, vdateb[0], vdateb[1]);
	      taxisDefVtimeBounds(taxisID2, vtimeb[0], vtimeb[1]);
	    }
	}

      streamDefTimestep(streamID2, tsID1);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamDefRecord(streamID2,  varID,  levelID);
	  
	  streamReadRecord(streamID1, array, &nmiss);
	  streamWriteRecord(streamID2, array, nmiss);
	}
      
      tsID1++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( array ) free(array);

  cdoFinish();

  return (0);
}

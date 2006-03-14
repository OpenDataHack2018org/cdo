/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2005 Uwe Schulzweida, schulzweida@dkrz.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include <string.h>
#include <ctype.h>  /* isdigit */

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

/*
@BeginDoc

@BeginModule

@Name      = Settime
@Title     = Set time
@Section   = Manipulating the header/field
@Class     = Manipulation
@Arguments = ifile ofile
@Operators = setdate settime setday setmon setyear settunits settaxis setreftime setcalendar shifttime

@EndModule


@BeginOperator_setdate

@Title     = Set date
@Parameter = date

@BeginDescription
Sets the date in every timestep to the same given value.
@EndDescription

@BeginParameter date
@Item = date
STRING  New date (format YYYY-MM-DD)
@EndParameter

@EndOperator


@BeginOperator_settime

@Title     = Set time
@Parameter = time

@BeginDescription
Sets the time in every timestep to the same given value.
@EndDescription

@BeginParameter timr
@Item = time
STRING  New time (format HH:MM)
@EndParameter

@EndOperator


@BeginOperator_setday

@Title     = Set day
@Parameter = day

@BeginDescription
Sets the day in every timestep to the same given value.
@EndDescription

@BeginParameter day
@Item = day
INTEGER  Value of the new day
@EndParameter

@EndOperator


@BeginOperator_setmon

@Title     = Set month
@Parameter = month

@BeginDescription
Sets the month in every timestep to the same given value.
@EndDescription

@BeginParameter month
@Item = month
INTEGER  Value of the new month
@EndParameter

@EndOperator


@BeginOperator_setyear

@Title     = Set year
@Parameter = year

@BeginDescription
Sets the year in every timestep to the same given value.
@EndDescription

@BeginParameter year
@Item = year
INTEGER  Value of the new year
@EndParameter

@EndOperator


@BeginOperator_settunits

@Title     = Set time units
@Parameter = units

@BeginDescription
Sets the time units.
@EndDescription

@BeginParameter units
@Item = units
STRING  Base units of the time axis (minutes, hours, days, months, years).
@EndParameter

@EndOperator


@BeginOperator_settaxis

@Title     = Set time axis
@Parameter = date time [inc]

@BeginDescription
Sets the time axis.
@EndDescription

@BeginParameter date
@Item = date
STRING  Start date (format YYYY-MM-DD)
@Item = time
STRING  Start time (format HH:MM)
@Item = inc
STRING  Optional increment (e.g. 12hour) [default: 0hour]
@EndParameter

@EndOperator


@BeginOperator_setreftime

@Title     = Set reference time
@Parameter = date time

@BeginDescription
Sets the reference time of an relative time axis.
@EndDescription

@BeginParameter date
@Item = date
STRING  Reference date (format YYYY-MM-DD)
@Item = time
STRING  Reference time (format HH:MM)
@EndParameter

@EndOperator


@BeginOperator_setcalendar

@Title     = Set calendar
@Parameter = calendar

@BeginDescription
Sets the calendar.
@EndDescription

@BeginParameter calendar
@Item = calendar
STRING  Calendar (standard, 360days, 365days, 366days)
@EndParameter

@EndOperator


@BeginOperator_shifttime

@Title     = Shift timesteps
@Parameter = sval

@BeginDescription
Shifts all timesteps by the parameter sval.
@EndDescription

@BeginParameter
@Item = sval
STRING   Shift value (e.g. -3hour)
@EndParameter

@EndOperator

@EndDoc
*/


void *Settime(void *argument)
{
  static char func[] = "Settime";
  int SETYEAR, SETMON, SETDAY, SETDATE, SETTIME, SETTUNITS;
  int SETTAXIS, SETREFTIME, SETCALENDAR, SHIFTTIME;
  int operatorID;
  int streamID1, streamID2 = CDI_UNDEFID;
  int nrecs, newval = 0;
  int tsID1, recID, varID, levelID;
  int vlistID1, vlistID2;
  int vdate, vtime;
  int sdate = 0, stime = 0;
  int taxisID1, taxisID2 = CDI_UNDEFID;
  int nmiss;
  int gridsize;
  int tunit = TUNIT_DAY;
  int ijulinc = 0, incperiod = 0, incunit = 0;
  int year, month, day, hour, minute;
  int dpy, calendar;
  int newcalendar = CALENDAR_STANDARD;
  const char *datestr, *timestr;
  double julval = 0;
  double *array = NULL;

  cdoInitialize(argument);

  SETYEAR     = cdoOperatorAdd("setyear",     0, 0, "year");
  SETMON      = cdoOperatorAdd("setmon",      0, 0, "month");
  SETDAY      = cdoOperatorAdd("setday",      0, 0, "day");
  SETDATE     = cdoOperatorAdd("setdate",     0, 0, "date (format YYYY-MM-DD)");
  SETTIME     = cdoOperatorAdd("settime",     0, 0, "time (format HH:MM)");
  SETTUNITS   = cdoOperatorAdd("settunits",   0, 0, "time units (minutes, hours, days, months, years)");
  SETTAXIS    = cdoOperatorAdd("settaxis",    0, 0, "date,time<,increment> (format YYYY-MM-DD,HH:MM)");
  SETREFTIME  = cdoOperatorAdd("setreftime",  0, 0, "date,time (format YYYY-MM-DD,HH:MM)");
  SETCALENDAR = cdoOperatorAdd("setcalendar", 0, 0, "calendar (standard, 360days, 365days, 366days)");
  SHIFTTIME   = cdoOperatorAdd("shifttime",   0, 0, "shift value");

  operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  if ( operatorID == SETTAXIS || operatorID == SETREFTIME )
    {
      if ( operatorArgc() < 2 ) cdoAbort("Not enough arguments!");

      datestr = operatorArgv()[0];
      timestr = operatorArgv()[1];

      if ( strchr(datestr, '-') == NULL )
	{
	  sdate = atoi(datestr);
	}
      else
	{
	  sscanf(datestr, "%d-%d-%d", &year, &month, &day);
	  sdate = year*10000 + month*100 + day;
	}

      if ( strchr(timestr, ':') == NULL )
	{
	  stime = atoi(timestr);
	}
      else
	{
	  sscanf(timestr, "%d:%d", &hour, &minute);
	  stime = hour*100 + minute;
	}

      if ( operatorArgc() == 3 )
	{
	  size_t len;
	  char *unit = operatorArgv()[2];
	  incperiod = atoi(unit);
	  while ( isdigit((int) *unit) ) unit++;
	  len = strlen(unit);      
	  if      ( strncmp(unit, "minutes", len) == 0 ) { incunit =    60; tunit = TUNIT_MINUTE;}
	  else if ( strncmp(unit, "hours", len)   == 0 ) { incunit =  3600; tunit = TUNIT_HOUR;}
	  else if ( strncmp(unit, "days", len)    == 0 ) { incunit = 86400; tunit = TUNIT_DAY;}
	  else if ( strncmp(unit, "months", len)  == 0 ) { incunit =     1; tunit = TUNIT_MONTH;}
	  else if ( strncmp(unit, "years", len)   == 0 ) { incunit =    12; tunit = TUNIT_YEAR;}
	  else cdoAbort("time unit >%s< unsupported", unit);
	}
      /* increment in seconds */
      ijulinc = incperiod * incunit;
    }
  else if ( operatorID == SETDATE )
    {
      if ( operatorArgc() < 1 ) cdoAbort("Not enough arguments!");
      datestr = operatorArgv()[0];
      if ( strchr(datestr, '-') == NULL )
	{
	  newval = atoi(datestr);
	}
      else
	{
	  sscanf(datestr, "%d-%d-%d", &year, &month, &day);
	  newval = year*10000 + month*100 + day;
	}
    }
  else if ( operatorID == SETTIME )
    {
      if ( operatorArgc() < 1 ) cdoAbort("Not enough arguments!");
      timestr = operatorArgv()[0];

      if ( strchr(timestr, ':') == NULL )
	{
	  newval = atoi(timestr);
	}
      else
	{
	  sscanf(timestr, "%d:%d", &hour, &minute);
	  newval = hour*100 + minute;
	}
    }
  else if ( operatorID == SHIFTTIME )
    {
	  size_t len;
	  char *unit = operatorArgv()[0];
	  incperiod = atoi(unit);
	  if ( unit[0] == '-' ) unit++;
	  while ( isdigit((int) *unit) ) unit++;
	  len = strlen(unit);      
	  if      ( strncmp(unit, "minutes", len) == 0 ) { incunit =    60; tunit = TUNIT_MINUTE;}
	  else if ( strncmp(unit, "hours", len)   == 0 ) { incunit =  3600; tunit = TUNIT_HOUR;}
	  else if ( strncmp(unit, "days", len)    == 0 ) { incunit = 86400; tunit = TUNIT_DAY;}
	  else if ( strncmp(unit, "months", len)  == 0 ) { incunit =     1; tunit = TUNIT_MONTH;}
	  else if ( strncmp(unit, "years", len)   == 0 ) { incunit =    12; tunit = TUNIT_YEAR;}
	  else cdoAbort("time unit >%s< unsupported", unit);

      /* increment in seconds */
      ijulinc = incperiod * incunit;
    }
  else if ( operatorID == SETTUNITS )
    {
      size_t len;
      char *unit = operatorArgv()[0];
      len = strlen(unit);      
      if      ( strncmp(unit, "minutes", len) == 0 ) { tunit = TUNIT_MINUTE;}
      else if ( strncmp(unit, "hours", len)   == 0 ) { tunit = TUNIT_HOUR;}
      else if ( strncmp(unit, "days", len)    == 0 ) { tunit = TUNIT_DAY;}
      else if ( strncmp(unit, "months", len)  == 0 ) { tunit = TUNIT_MONTH;}
      else if ( strncmp(unit, "years", len)   == 0 ) { tunit = TUNIT_YEAR;}
      else cdoAbort("time unit >%s< unsupported", unit);
    }
  else if ( operatorID == SETCALENDAR )
    {
      size_t len;
      char *cname = operatorArgv()[0];
      len = strlen(cname);      
      if      ( strncmp(cname, "standard", len) == 0 ) { newcalendar = CALENDAR_STANDARD;}
      else if ( strncmp(cname, "360days", len)  == 0 ) { newcalendar = CALENDAR_360DAYS;}
      else if ( strncmp(cname, "365days", len)  == 0 ) { newcalendar = CALENDAR_365DAYS;}
      else if ( strncmp(cname, "366days", len)  == 0 ) { newcalendar = CALENDAR_366DAYS;}
      else cdoAbort("calendar >%s< unsupported", cname);
    }
  else
    {
      newval = atoi(operatorArgv()[0]);
    }

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);

  calendar = taxisInqCalendar(taxisID1);

  dpy = calendar_dpy(calendar);

  if ( cdoVerbose )
    cdoPrint("calendar = %d;  dpy = %d\n", calendar, dpy);

  if ( operatorID == SETREFTIME )
    {
      if ( taxisInqType(taxisID1) == TAXIS_ABSOLUTE )
	cdoPrint("Changing absolute to relative time axis!");

      taxisID2 = taxisNew(TAXIS_RELATIVE);
      if ( operatorArgc() != 3 ) tunit = taxisInqTunit(taxisID1);
      taxisDefTunit(taxisID2, tunit);
    }
  else if ( operatorID == SETCALENDAR )
    {
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
	  taxisID2 = taxisNew(TAXIS_RELATIVE);
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
      julval = encode_julval(dpy, sdate, stime);
    }
  else if ( operatorID == SETTUNITS )
    {
      taxisDefTunit(taxisID2, tunit);
    }
  else if ( operatorID == SETREFTIME )
    {
      taxisDefRdate(taxisID2, sdate);
      taxisDefRtime(taxisID2, stime);
    }
  else if ( operatorID == SETCALENDAR )
    {
      taxisDefCalendar(taxisID2, newcalendar);
    }

  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID2);
  array = (double *) malloc(gridsize*sizeof(double));

  tsID1 = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID1)) )
    {
      if ( operatorID == SETTAXIS )
	{
	  if ( incunit == 1 || incunit == 12 )
	    {
	      vtime = stime;
	      if ( tsID1 == 0 )
		{
		  vdate = sdate;
		}
	      else
		{
		  year  = vdate/10000;
		  month = (vdate - year*10000) / 100;
		  day   =  vdate - year*10000 - month*100;
	      
		  month += ijulinc;

		  while ( month > 12 ) { month -= 12; year++; }
		  while ( month <  1 ) { month += 12; year--; }

		  vdate = year*10000 + month*100 + day;
		}
	    }
	  else
	    {
	      decode_julval(dpy, julval, &vdate, &vtime);
	      julval += ijulinc;
	    }
	}
      else if ( operatorID == SHIFTTIME )
	{
	  vdate = taxisInqVdate(taxisID1);
	  vtime = taxisInqVtime(taxisID1);

	  if ( incunit == 1 || incunit == 12 )
	    {
	      decode_date(vdate, &year, &month, &day);
	      
	      month += ijulinc;

	      while ( month > 12 )
		{
		  month -= 12;
		  year++;
		}

	      while ( month < 1 )
		{
		  month += 12;
		  year--;
		}

	      vdate = year*10000 + month*100 + day;
	    }
	  else
	    {
	      julval = encode_julval(dpy, vdate, vtime);
	      julval += ijulinc;
	      decode_julval(dpy, julval, &vdate, &vtime);
	    }
	}
      else if ( operatorID == SETREFTIME )
	{
	  vdate = taxisInqVdate(taxisID1);
	  vtime = taxisInqVtime(taxisID1);
	}
      else if ( operatorID == SETCALENDAR )
	{
	  vdate = taxisInqVdate(taxisID1);
	  vtime = taxisInqVtime(taxisID1);
	}
      else
	{
	  vdate = taxisInqVdate(taxisID1);
	  vtime = taxisInqVtime(taxisID1);

	  decode_date(vdate, &year, &month, &day);

	  if ( operatorID == SETYEAR ) year  = newval;
	  if ( operatorID == SETMON  ) month = newval;
	  if ( operatorID == SETDAY  ) day   = newval;
      
	  vdate = year*10000 + month*100 + day;

	  if ( operatorID == SETDATE  ) vdate = newval;
	  if ( operatorID == SETTIME  ) vtime = newval;
	}

      taxisDefVdate(taxisID2, vdate);
      taxisDefVtime(taxisID2, vtime);

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

  streamClose(streamID1);
  streamClose(streamID2);

  if ( array ) free(array);

  cdoFinish();

  return (0);
}

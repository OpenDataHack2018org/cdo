/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2007 Uwe Schulzweida, schulzweida@dkrz.de
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




void *Tinfo(void *argument)
{
  int vdate, vtime;
  int nrecs, ntsteps;
  int tsID, ntimeout;
  int timeID, taxisID;
  int nbyte, nbyte0;
  int index;
  int streamID;
  int vlistID;
  int year, month, day, hour, minute;

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
	  int calendar, unit;
	  
	  if ( taxisInqType(taxisID) == TAXIS_RELATIVE )
	    {
	      vdate = taxisInqRdate(taxisID);
	      vtime = taxisInqRtime(taxisID);
	      
	      decode_date(vdate, &year, &month, &day);
	      decode_time(vtime, &hour, &minute);

	      fprintf(stdout, "     RefTime = %4.4d-%2.2d-%2.2d %2.2d:%2.2d",
		      year, month, day, hour, minute);
		      
	      unit = taxisInqTunit(taxisID);
	      if ( unit != CDI_UNDEFID )
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
	      
	      calendar = taxisInqCalendar(taxisID);
	      if ( calendar != CDI_UNDEFID )
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

	      fprintf(stdout, "\n");
	    }
	}

      fprintf(stdout, "Timestep   YYYY-MM-DD hh:mm\n");

      ntimeout = 0;
      tsID = 0;
      while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	{  
	  vdate = taxisInqVdate(taxisID);
	  vtime = taxisInqVtime(taxisID);
	  
	  decode_date(vdate, &year, &month, &day);
	  decode_time(vtime, &hour, &minute);

	  fprintf(stdout, "%5d   %4.4d-%2.2d-%2.2d %2.2d:%2.2d\n", tsID+1, year, month, day, hour, minute);
	  ntimeout++;
	  tsID++;
	}
    }

  streamClose(streamID);

  cdoFinish();

  return (0);
}

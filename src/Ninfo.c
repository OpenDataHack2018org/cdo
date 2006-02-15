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

#include <stdio.h>
#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

/*
@BeginDoc

@BeginModule

@Name      = Ninfo
@Title     = 
@Section   = Information
@Class     = Information
@Arguments = ifile
@Operators = nyear nmon ndate ntime ncode nvar nlevel

@EndModule


@BeginOperator_nyear

@Title     = Number of years

@BeginDesciption
Prints the number of different years.
@EndDesciption

@EndOperator


@BeginOperator_nmon

@Title     = Number of months

@BeginDesciption
Prints the number of different combinations of years and months.
@EndDesciption

@EndOperator


@BeginOperator_ndate

@Title     = Number of dates

@BeginDesciption
Prints the number of different dates.
@EndDesciption

@EndOperator


@BeginOperator_ntime

@Title     = Number of timesteps

@BeginDesciption
Prints the number of timesteps.
@EndDesciption

@EndOperator


@BeginOperator_ncode

@Title     = Number of codes

@BeginDesciption
Prints the number of different codes.
@EndDesciption

@EndOperator


@BeginOperator_nvar

@Title     = Number of variables

@BeginDesciption
Prints the number of different variables.
@EndDesciption

@EndOperator

@BeginOperator_nlevel

@Title     = Number of levels

@BeginDesciption
Prints the number of levels for each variable.
@EndDesciption

@EndOperator

@EndDoc
*/

void *Ninfo(void *argument)
{
  enum {NYEAR, NMON, NDATE, NTIME, NCODE, NVAR, NLEVEL};
  int operatorID;
  int operfunc;
  int varID, zaxisID;
  int vdate;
  int nrecs, nvars;
  int levelsize;
  int tsID, ndate, date0 = 0;
  int mon0 = 0, mon, nmon, year0 = 0, year, nyear;
  int taxisID;
  int streamID;
  int vlistID;

  cdoInitialize(argument);

  cdoOperatorAdd("nyear",  NYEAR,  0, NULL);
  cdoOperatorAdd("nmon",   NMON,   0, NULL);
  cdoOperatorAdd("ndate",  NDATE,  0, NULL);
  cdoOperatorAdd("ntime",  NTIME,  0, NULL);
  cdoOperatorAdd("ncode",  NCODE,  0, NULL);
  cdoOperatorAdd("nvar",   NVAR,   0, NULL);
  cdoOperatorAdd("nlevel", NLEVEL, 0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

  streamID = streamOpenRead(cdoStreamName(0));
  if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", cdoStreamName(0));

  vlistID = streamInqVlist(streamID);

  nvars   = vlistNvars(vlistID);
  taxisID = vlistInqTaxis(vlistID);

  switch ( operfunc )
    {
    case NYEAR:
      nyear = 0;
      tsID = 0;
      while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	{
	  vdate = taxisInqVdate(taxisID);

	  year = vdate/10000;
	 
	  if ( tsID == 0 || year0 != year )
	    {
	      year0 = year;
	      nyear++;
	    }

	  tsID++;
	}
      fprintf(stdout, "%d\n", nyear);
      break;
    case NMON:
      nmon = 0;
      tsID = 0;
      while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	{
	  vdate = taxisInqVdate(taxisID);

	  year = vdate/10000;
	  mon  = (vdate - year*10000) / 100;
	 
	  if ( tsID == 0 || mon0 != mon )
	    {
	      mon0 = mon;
	      nmon++;
	    }

	  tsID++;
	}
      fprintf(stdout, "%d\n", nmon);
      break;
    case NDATE:
      ndate = 0;
      tsID = 0;
      while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	{
	  vdate = taxisInqVdate(taxisID);
	 
	  if ( tsID == 0 || date0 != vdate )
	    {
	      date0 = vdate;
	      ndate++;
	    }

	  tsID++;
	}
      fprintf(stdout, "%d\n", ndate);
      break;
    case NTIME:
      tsID = 0;
      while ( (nrecs = streamInqTimestep(streamID, tsID)) ) tsID++;
      fprintf(stdout, "%d\n", tsID);
      break;
    case NCODE:
    case NVAR:
      fprintf(stdout, "%d\n", nvars);
      break;
    case NLEVEL:
      for ( varID = 0; varID < nvars; varID++ )
	{
	  zaxisID = vlistInqVarZaxis(vlistID, varID);
	  levelsize = zaxisInqSize(zaxisID);
	  fprintf(stdout, "%d\n", levelsize);
	}
      break;
    default:
      cdoAbort("operator not implemented!");
    }

  streamClose(streamID);

  cdoFinish();

  return (0);
}

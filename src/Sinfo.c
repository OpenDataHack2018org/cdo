/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2006 Uwe Schulzweida, schulzweida@dkrz.de
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
#include <string.h>
#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

/*
@BeginDoc

@BeginModule

@Name      = Sinfo
@Title     = Short information
@Section   = Information
@Class     = Information
@Arguments = ifile
@Operators = sinfo sinfov

@EndModule


@BeginOperator_sinfo

@Title     = Short file information

@BeginDesciption
Prints short information for each variable of a file.
For each variable the operator print in one line the:
@BeginItemize
@Item = variable number
@Item = institute and source
@Item = code and codetable
@Item = horizontal grid size and number
@Item = vertical grid size and number
@EndItemize
@EndDesciption

@EndOperator


@BeginOperator_sinfov

@Title     = Short file information

@BeginDesciption
The same as operator sinfo. Using the name instead of the code number
to identify the variable.
@EndDesciption

@EndOperator

@EndDoc
*/

#define  CDI_BIGENDIAN            0   /* Data type BIGENDIAN     */
#define  CDI_LITTLEENDIAN         1   /* Data type LITTLEENDIAN  */

static void printFiletype(int streamID)
{
  int filetype;

  filetype = streamInqFiletype(streamID);

  switch ( filetype )
    {
    case FILETYPE_GRB:
      printf("   File format: GRIB");
      break;
    case FILETYPE_NC:
      printf("   File format: netCDF");
      break;
    case FILETYPE_NC2:
      printf("   File format: netCDF2");
      break;
    case FILETYPE_SRV:
      printf("   File format: SERVICE");
      switch ( streamInqByteorder(streamID) )
	{
	case CDI_BIGENDIAN:
	  printf("  BIGENDIAN"); break;
	case CDI_LITTLEENDIAN:
	  printf("  LITTLEENDIAN"); break;
	default:
	  printf("  byteorder: %d undefined", streamInqByteorder(streamID)); break;
	}
      break;
    case FILETYPE_EXT:
      printf("   File format: EXTRA");
      switch ( streamInqByteorder(streamID) )
	{
	case CDI_BIGENDIAN:
	  printf("  BIGENDIAN"); break;
	case CDI_LITTLEENDIAN:
	  printf("  LITTLEENDIAN"); break;
	default:
	  printf("  byteorder: %d undefined", streamInqByteorder(streamID)); break;
	}
      break;
    case FILETYPE_IEG:
      printf("   File format: IEG");
      switch ( streamInqByteorder(streamID) )
	{
	case CDI_BIGENDIAN:
	  printf("  BIGENDIAN"); break;
	case CDI_LITTLEENDIAN:
	  printf("  LITTLEENDIAN"); break;
	default:
	  printf("  byteorder: %d undefined", streamInqByteorder(streamID)); break;
	}
      break;
    default:
      printf("  File format: unsupported filetype %d" , filetype);
    }

  printf("\n");
}


void *Sinfo(void *argument)
{
  int SINFO, SINFOV;
  int operatorID;
  int indf;
  int varID;
  int gridsize = 0;
  int gridID, zaxisID, code;
  int vdate, vtime;
  int nrecs, nvars, ngrids, nzaxis, ntsteps;
  int levelID, levelsize;
  int tsID, ntimeout;
  int xsize, ysize, trunc;
  int timeID, taxisID;
  int gridtype, nbyte, nbyte0;
  int index;
  char varname[128];
  char longname[128];
  char units[128];
  double level;
  char *modelptr, *instptr;
  int streamID = 0;
  int vlistID;
  int prec;
  int year, month, day, hour, minute;
  char pstr[4];

  cdoInitialize(argument);

  SINFO  = cdoOperatorAdd("sinfo",  0, 0, NULL);
  SINFOV = cdoOperatorAdd("sinfov", 0, 0, NULL);

  operatorID = cdoOperatorID();

  for ( indf = 0; indf < cdoStreamCnt(); indf++ )
    {
      streamID = streamOpenRead(cdoStreamName(indf));
      if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", cdoStreamName(indf));

      vlistID = streamInqVlist(streamID);

      printFiletype(streamID);

      if ( operatorID == SINFOV )
	fprintf(stdout,
		"%6d : Institut Source  Table   Time   Typ  Grid Size Num  Levels Num Varname\n",  -(indf+1));
      else
	fprintf(stdout,
		"%6d : Institut Source  Table Code   Time   Typ  Grid Size Num  Levels Num\n",  -(indf+1));

      nvars = vlistNvars(vlistID);

      for ( varID = 0; varID < nvars; varID++ )
	{
	  code    = vlistInqVarCode(vlistID, varID);
	  gridID  = vlistInqVarGrid(vlistID, varID);
	  zaxisID = vlistInqVarZaxis(vlistID, varID);

	  if ( operatorID == SINFOV ) vlistInqVarName(vlistID, varID, varname);

	  gridsize = gridInqSize(gridID);

	  fprintf(stdout, "%6d : ", varID + 1);

	  instptr = institutInqNamePtr(vlistInqVarInstitut(vlistID, varID));
	  if ( instptr )
	    fprintf(stdout, "%-9s", instptr);
	  else
	    fprintf(stdout, "unknown  ");

	  modelptr = modelInqNamePtr(vlistInqVarModel(vlistID, varID));
	  if ( modelptr )
	    fprintf(stdout, "%-9s", modelptr);
	  else
	    fprintf(stdout, "unknown  ");

	  fprintf(stdout, "%4d", tableInqNum(vlistInqVarTable(vlistID, varID)));

	  if ( operatorID != SINFOV )
	    fprintf(stdout, "%5d", code);

	  timeID = vlistInqVarTime(vlistID, varID);
	  if ( timeID == TIME_CONSTANT )
	    fprintf(stdout, " constant");
	  else
	    fprintf(stdout, " variable");

	  prec = vlistInqVarDatatype(vlistID, varID);
	  if      ( prec == DATATYPE_PACK  ) strcpy(pstr, "P0");
	  else if ( prec == DATATYPE_PACK1 ) strcpy(pstr, "P1");
	  else if ( prec == DATATYPE_PACK2 ) strcpy(pstr, "P2");
	  else if ( prec == DATATYPE_PACK3 ) strcpy(pstr, "P3");
	  else if ( prec == DATATYPE_REAL4 ) strcpy(pstr, "R4");
	  else if ( prec == DATATYPE_REAL8 ) strcpy(pstr, "R8");
	  else if ( prec == DATATYPE_INT1  ) strcpy(pstr, "I1");
	  else if ( prec == DATATYPE_INT2  ) strcpy(pstr, "I2");
	  else if ( prec == DATATYPE_INT4  ) strcpy(pstr, "I4");
	  else                               strcpy(pstr, "-1");

	  fprintf(stdout, " %-3s", pstr);

	  fprintf(stdout, " %9d", gridsize);

	  fprintf(stdout, " %3d ", gridID + 1);

	  levelsize = zaxisInqSize(zaxisID);
	  fprintf(stdout, " %6d", levelsize);
	  fprintf(stdout, " %3d", zaxisID + 1);
	    
	  if ( operatorID == SINFOV )
	    fprintf(stdout, "  %-10s", varname);

	  fprintf(stdout, "\n");
	}

      ngrids = vlistNgrids(vlistID);
      fprintf(stdout, "   Horizontal grids :\n");
      for ( index = 0; index < ngrids; index++ )
	{
	  gridID   = vlistGrid(vlistID, index);
	  gridtype = gridInqType(gridID);
	  trunc    = gridInqTrunc(gridID);
	  gridsize = gridInqSize(gridID);
	  xsize    = gridInqXsize(gridID);
	  ysize    = gridInqYsize(gridID);

	  /*	  nbyte0   = fprintf(stdout, "  %4d : %-23s : ",*/
	  nbyte0   = fprintf(stdout, "  %4d : %-12s > ",
			     gridID+1, gridNamePtr(gridtype));

	  if ( gridtype == GRID_LONLAT   ||
	       gridtype == GRID_GAUSSIAN ||
	       gridtype == GRID_GAUSSIAN_REDUCED )
	    {
	      double lonfirst = 0.0, lonlast = 0.0;
	      double latfirst = 0.0, latlast = 0.0;
	      double loninc = 0.0, latinc = 0.0;

	      latfirst = gridInqYval(gridID, 0);
	      latlast  = gridInqYval(gridID, ysize-1);
	      latinc   = gridInqYinc(gridID);
	      if ( gridtype == GRID_GAUSSIAN_REDUCED )
		{
		  fprintf(stdout, "size : dim = %d  nlat = %d\n", gridsize, ysize);
		  fprintf(stdout, "%*s", nbyte0, "");
		  fprintf(stdout, "longitude : reduced\n");
		}
	      else
		{
		  lonfirst = gridInqXval(gridID, 0);
		  lonlast  = gridInqXval(gridID, xsize-1);
		  loninc   = gridInqXinc(gridID);
		  fprintf(stdout, "size      : dim = %d  nlon = %d  nlat = %d\n", gridsize, xsize, ysize);
		  fprintf(stdout, "%*s", nbyte0, "");
		  fprintf(stdout, "longitude : first = %.9g  last = %.9g", lonfirst, lonlast);
		  if ( !DBL_IS_EQUAL(loninc, 0) )
		    fprintf(stdout, "  inc = %.9g", loninc);
		  fprintf(stdout, "\n");
		}
	      fprintf(stdout, "%*s", nbyte0, "");
	      fprintf(stdout, "latitude  : first = %.9g  last = %.9g", latfirst, latlast);
	      if ( !DBL_IS_EQUAL(latinc, 0) )
		fprintf(stdout, "  inc = %.9g", latinc);
	      fprintf(stdout, "\n");

	      if ( gridIsRotated(gridID) )
		{
		  double lonpole, latpole;
		  lonpole = gridInqXpole(gridID);
		  latpole = gridInqYpole(gridID);
		  fprintf(stdout, "%*s", nbyte0, "");
		  fprintf(stdout, "northpole : lon = %.9g  lat = %.9g\n", lonpole, latpole);
		}
		
	      if ( gridInqXbounds(gridID, NULL) || gridInqYbounds(gridID, NULL) )
		{
		  fprintf(stdout, "%*s", nbyte0, "");
		  fprintf(stdout, "available :");
		  if ( gridInqXbounds(gridID, NULL) ) fprintf(stdout, " xbounds");
		  if ( gridInqYbounds(gridID, NULL) ) fprintf(stdout, " ybounds");
		  fprintf(stdout, "\n");
		}
	    }
	  else if ( gridtype == GRID_SPECTRAL )
	    {
	      fprintf(stdout, "size      : dim = %d  truncation = %d  spc = %d\n",
		      gridsize, trunc, gridsize/2);
	    }
	  else if ( gridtype == GRID_GME )
	    {
	      int ni, nd;
	      ni = gridInqGMEni(gridID);
	      nd = gridInqGMEnd(gridID);
	      fprintf(stdout, "size      : dim = %d  nd = %d  ni = %d\n", gridsize, nd, ni);
	    }
	  else if ( gridtype == GRID_CURVILINEAR )
	    {
	      fprintf(stdout, "size      : dim = %d  nx = %d  ny = %d\n", gridsize, xsize, ysize);
	    }
	  else if ( gridtype == GRID_CELL )
	    {
	      fprintf(stdout, "size      : dim = %d  nvertex = %d\n", gridsize, gridInqNvertex(gridID));
	    }
	  else /* if ( gridtype == GRID_GENERIC ) */
	    {
	      if ( ysize == 0 )
		fprintf(stdout, "size      : dim = %d\n", gridsize);
	      else
		fprintf(stdout, "size      : dim = %d  nx = %d  ny = %d\n", gridsize, xsize, ysize);
	    }

	  if ( gridtype == GRID_CURVILINEAR || gridtype == GRID_CELL || gridtype == GRID_GENERIC )
	    {
	      if ( gridInqXvals(gridID, NULL) || gridInqYvals(gridID, NULL) || gridHasArea(gridID) ||
		   gridInqXbounds(gridID, NULL) || gridInqYbounds(gridID, NULL) )
		{
		  fprintf(stdout, "%*s", nbyte0, "");
		  fprintf(stdout, "available :");
		  if ( gridInqXvals(gridID, NULL) )   fprintf(stdout, " xvals");
		  if ( gridInqYvals(gridID, NULL) )   fprintf(stdout, " yvals");
		  if ( gridInqXbounds(gridID, NULL) ) fprintf(stdout, " xbounds");
		  if ( gridInqYbounds(gridID, NULL) ) fprintf(stdout, " ybounds");
		  if ( gridHasArea(gridID) )          fprintf(stdout, " area");
		  fprintf(stdout, "\n");
		}
	    }
	}

      nzaxis = vlistNzaxis(vlistID);
      fprintf(stdout, "   Vertical grids :\n");
      for ( index = 0; index < nzaxis; index++)
	{
	  zaxisID   = vlistZaxis(vlistID, index);
	  levelsize = zaxisInqSize(zaxisID);
	  /* zaxisInqLongname(zaxisID, longname); */
	  zaxisName(zaxisInqType(zaxisID), longname);
	  longname[16] = 0;
	  zaxisInqUnits(zaxisID, units);
	  units[12] = 0;
	  nbyte0    = fprintf(stdout, "  %4d : %-16s  %5s : ", zaxisID+1, longname, units);
	  nbyte = nbyte0;
	  for ( levelID = 0; levelID < levelsize; levelID++ )
	    {
	      if ( nbyte > 80 )
		{
		  fprintf(stdout, "\n");
		  fprintf(stdout, "%*s", nbyte0, "");
		  nbyte = nbyte0;
		}
	      level = zaxisInqLevel(zaxisID, levelID);
	      nbyte += fprintf(stdout, "%.9g ", level);
	    }
	  fprintf(stdout, "\n");
	  if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
	    {
	      double level1, level2;
	      nbyte = nbyte0;
	      nbyte0 = fprintf(stdout, "%32s : ", "bounds");
	      for ( levelID = 0; levelID < levelsize; levelID++ )
		{
		  if ( nbyte > 80 )
		    {
		      fprintf(stdout, "\n");
		      fprintf(stdout, "%*s", nbyte0, "");
		      nbyte = nbyte0;
		    }
		  level1 = zaxisInqLbound(zaxisID, levelID);
		  level2 = zaxisInqUbound(zaxisID, levelID);
		  nbyte += fprintf(stdout, "%.9g-%.9g ", level1, level2);
		}
	      fprintf(stdout, "\n");
	    }
	}

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

	  fprintf(stdout, "   YYYY-MM-DD HH:MM   YYYY-MM-DD HH:MM   YYYY-MM-DD HH:MM   YYYY-MM-DD HH:MM\n");

	  ntimeout = 0;
	  tsID = 0;
	  while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	    {
	      if ( ntimeout == 4 )
		{
		  ntimeout = 0;
		  fprintf(stdout, "\n");
		}

	      vdate = taxisInqVdate(taxisID);
	      vtime = taxisInqVtime(taxisID);

	      decode_date(vdate, &year, &month, &day);
	      decode_time(vtime, &hour, &minute);

	      fprintf(stdout, "   %4.4d-%2.2d-%2.2d %2.2d:%2.2d", year, month, day, hour, minute);
	      ntimeout++;
	      tsID++;
	    }
	  fprintf(stdout, "\n");
	}

      streamClose(streamID);
    }

  cdoFinish();

  return (0);
}

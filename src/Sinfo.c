/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2008 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Sinfo      sinfo           Short dataset information
*/


#include <stdio.h>
#include <string.h>
#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void printFiletype(int streamID, int vlistID)
{
  int filetype;

  filetype = streamInqFiletype(streamID);

  switch ( filetype )
    {
    case FILETYPE_GRB:
      printf("GRIB");
      break;
    case FILETYPE_NC:
      printf("netCDF");
      break;
    case FILETYPE_NC2:
      printf("netCDF2");
      break;
    case FILETYPE_NC4:
      printf("netCDF4");
      break;
    case FILETYPE_SRV:
      printf("SERVICE");
      break;
    case FILETYPE_EXT:
      printf("EXTRA");
      break;
    case FILETYPE_IEG:
      printf("IEG");
      break;
    default:
      printf("  File format: unsupported filetype %d" , filetype);
    }

  if ( filetype == FILETYPE_SRV || filetype == FILETYPE_EXT || filetype == FILETYPE_IEG )
    {
      switch ( streamInqByteorder(streamID) )
	{
	case CDI_BIGENDIAN:
	  printf("  BIGENDIAN"); break;
	case CDI_LITTLEENDIAN:
	  printf("  LITTLEENDIAN"); break;
	default:
	  printf("  byteorder: %d undefined", streamInqByteorder(streamID)); break;
	}
    }

  if ( filetype == FILETYPE_GRB || filetype == FILETYPE_NC4 )
    {
      int nvars, varID;
      int ztype;

      nvars = vlistNvars(vlistID);

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( ztype = vlistInqVarZtype(vlistID, varID) )
	    {
	      if ( ztype == COMPRESS_SZIP )
		printf(" SZIP");
	      else if ( ztype == COMPRESS_ZIP )
		printf(" ZIP");

	      break;
	    }
	}
    }

  printf("\n");
}


static void printGridInfo(int vlistID)
{
  static char func[] = "printGridInfo";
  int ngrids, index;
  int gridID, gridtype, trunc, gridsize, xsize, ysize;
  int nbyte0;

  ngrids = vlistNgrids(vlistID);
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
	      if ( gridIsCircular(gridID) )
		fprintf(stdout, "  circular");
	      fprintf(stdout, "\n");
	    }
	  fprintf(stdout, "%*s", nbyte0, "");
	  fprintf(stdout, "latitude  : first = %.9g  last = %.9g", latfirst, latlast);
	  if ( !DBL_IS_EQUAL(latinc, 0) && gridtype == GRID_LONLAT )
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
      else if ( gridtype == GRID_CURVILINEAR || gridtype == GRID_CELL )
	{
	  if ( gridtype == GRID_CURVILINEAR )
	    fprintf(stdout, "size      : dim = %d  nx = %d  ny = %d\n", gridsize, xsize, ysize);
	  else
	    fprintf(stdout, "size      : dim = %d  nvertex = %d\n", gridsize, gridInqNvertex(gridID));

	  if ( gridInqXvals(gridID, NULL) && gridInqYvals(gridID, NULL) )
	    {
	      int i;
	      char xunits[256];
	      char yunits[256];
	      double *xvals, *yvals;
	      double xfirst, xlast, yfirst, ylast;
	      xvals = (double *) malloc(gridsize*sizeof(double));
	      yvals = (double *) malloc(gridsize*sizeof(double));

	      gridInqXvals(gridID, xvals);
	      gridInqYvals(gridID, yvals);
	      gridInqXunits(gridID, xunits);
	      gridInqYunits(gridID, yunits);

	      xfirst = xvals[0];
	      xlast  = xvals[0];
	      yfirst = yvals[0];
	      ylast  = yvals[0];
	      for ( i = 1; i < gridsize; i++ )
		{
		  if ( xvals[i] < xfirst ) xfirst = xvals[i];
		  if ( xvals[i] > xlast )  xlast  = xvals[i];
		  if ( yvals[i] < yfirst ) yfirst = yvals[i];
		  if ( yvals[i] > ylast )  ylast  = yvals[i];
		}

	      fprintf(stdout, "%*s", nbyte0, "");
	      fprintf(stdout, "longitude : min = %.9g  max = %.9g  %s", xfirst, xlast, xunits);
	      if ( gridIsCircular(gridID) )
		fprintf(stdout, "  circular");
	      fprintf(stdout, "\n");
	      fprintf(stdout, "%*s", nbyte0, "");
	      fprintf(stdout, "latitude  : min = %.9g  max = %.9g  %s\n", yfirst, ylast, yunits);
	      
	      free(xvals);
	      free(yvals);
	    }
	}
      else if ( gridtype == GRID_LCC )
	{
	  double originLon, originLat, lonParY, lat1, lat2, xincm, yincm;
	  int projflag, scanflag;

	  gridInqLCC(gridID, &originLon, &originLat, &lonParY, &lat1, &lat2, &xincm, &yincm,
		     &projflag, &scanflag);

	  fprintf(stdout, "size      : dim = %d  nx = %d  ny = %d  ", gridsize, xsize, ysize);
	  if ( (projflag&128) == 0 )
	    fprintf(stdout, "North Pole\n");
	  else
	    fprintf(stdout, "South Pole\n");
	  fprintf(stdout, "%*s", nbyte0, "");	  
	  fprintf(stdout, "            originLon = %g  originLat = %g  lonParY = %g\n",
		  originLon, originLat, lonParY);
	  fprintf(stdout, "%*s", nbyte0, "");	  
	  fprintf(stdout, "            lat1 = %g  lat2 = %g  xinc = %gm  yinc = %gm\n", 
		  lat1, lat2, xincm, yincm);
	}
      else /* if ( gridtype == GRID_GENERIC ) */
	{
	  if ( ysize == 0 )
	    fprintf(stdout, "size      : dim = %d\n", gridsize);
	  else
	    {
	      fprintf(stdout, "size      : dim = %d  nx = %d  ny = %d\n", gridsize, xsize, ysize);
	      if ( gridIsCircular(gridID) )
		{
		  fprintf(stdout, "%*s", nbyte0, "");	  
		  fprintf(stdout, "longitude :  circular\n");
		}
	    }
	}

      if ( gridtype == GRID_CURVILINEAR || gridtype == GRID_CELL ||
	   gridtype == GRID_GENERIC || gridtype == GRID_LCC )
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
}


void *Sinfo(void *argument)
{
  int SINFO, SINFOV, SINFOP;
  int operatorID;
  int indf;
  int varID;
  int gridsize = 0;
  int gridID, zaxisID, code;
  int zaxistype, ltype;
  int vdate, vtime;
  int nrecs, nvars, nzaxis, ntsteps;
  int levelID, levelsize;
  int tsID, ntimeout;
  int timeID, taxisID;
  int nbyte, nbyte0;
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
  SINFOP = cdoOperatorAdd("sinfop", 0, 0, NULL);

  operatorID = cdoOperatorID();

  for ( indf = 0; indf < cdoStreamCnt(); indf++ )
    {
      streamID = streamOpenRead(cdoStreamName(indf));
      if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", cdoStreamName(indf));

      vlistID = streamInqVlist(streamID);

      printf("   File format: ");
      printFiletype(streamID, vlistID);

      if ( operatorID == SINFOV )
	fprintf(stdout,
		"%6d : Institut Source   Varname      Time   Typ  Grid Size Num  Levels Num\n",  -(indf+1));
      else if ( operatorID == SINFOP )
	fprintf(stdout,
		"%6d : Institut Source   Param     Time   Typ  Grid Size Num  Levels Num\n",  -(indf+1));
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

	  if ( operatorID == SINFOV )
	    fprintf(stdout, "%-10s", varname);
	  else if ( operatorID == SINFOP )
	    fprintf(stdout, "%03d.%03d ", tableInqNum(vlistInqVarTable(vlistID, varID)), code);
	  else
	    fprintf(stdout, "%4d %4d", tableInqNum(vlistInqVarTable(vlistID, varID)), code);

	  timeID = vlistInqVarTime(vlistID, varID);
	  if ( timeID == TIME_CONSTANT )
	    fprintf(stdout, " constant");
	  else
	    fprintf(stdout, " variable");

	  prec = vlistInqVarDatatype(vlistID, varID);
	  if      ( prec == DATATYPE_PACK   ) strcpy(pstr, "P0");
	  else if ( prec > 0 && prec <= 32  ) sprintf(pstr, "P%d", prec);
	  else if ( prec == DATATYPE_FLT32  ) strcpy(pstr, "F32");
	  else if ( prec == DATATYPE_FLT64  ) strcpy(pstr, "F64");
	  else if ( prec == DATATYPE_INT8   ) strcpy(pstr, "I8");
	  else if ( prec == DATATYPE_INT16  ) strcpy(pstr, "I16");
	  else if ( prec == DATATYPE_INT32  ) strcpy(pstr, "I32");
	  else if ( prec == DATATYPE_UINT8  ) strcpy(pstr, "U8");
	  else if ( prec == DATATYPE_UINT16 ) strcpy(pstr, "U16");
	  else if ( prec == DATATYPE_UINT32 ) strcpy(pstr, "U32");
	  else                                strcpy(pstr, "-1");

	  fprintf(stdout, " %-3s", pstr);

	  if ( vlistInqVarZtype(vlistID, varID) == COMPRESS_NONE )
	    fprintf(stdout, " ");
	  else
	    fprintf(stdout, "z");

	  fprintf(stdout, "%9d", gridsize);

	  fprintf(stdout, " %3d ", gridID + 1);

	  levelsize = zaxisInqSize(zaxisID);
	  fprintf(stdout, " %6d", levelsize);
	  fprintf(stdout, " %3d", zaxisID + 1);

	  fprintf(stdout, "\n");
	}

      fprintf(stdout, "   Horizontal grids :\n");
      printGridInfo(vlistID);

      nzaxis = vlistNzaxis(vlistID);
      fprintf(stdout, "   Vertical grids :\n");
      for ( index = 0; index < nzaxis; index++)
	{
	  zaxisID   = vlistZaxis(vlistID, index);
	  zaxistype = zaxisInqType(zaxisID);
	  ltype     = zaxisInqLtype(zaxisID);
	  levelsize = zaxisInqSize(zaxisID);
	  /* zaxisInqLongname(zaxisID, longname); */
	  zaxisName(zaxistype, longname);
	  longname[16] = 0;
	  zaxisInqUnits(zaxisID, units);
	  units[12] = 0;
	  if ( zaxistype == ZAXIS_GENERIC && ltype != 0 )
	    nbyte0    = fprintf(stdout, "  %4d : %-10s  (ltype=%3d) : ", zaxisID+1, longname, ltype);
	  else
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

	  fprintf(stdout, "   YYYY-MM-DD hh:mm   YYYY-MM-DD hh:mm   YYYY-MM-DD hh:mm   YYYY-MM-DD hh:mm\n");

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

	      fprintf(stdout, " %6.4d-%2.2d-%2.2d %2.2d:%2.2d", year, month, day, hour, minute);
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

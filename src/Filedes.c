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

/*
   This module contains the following operators:

      Filedes    vardes          Variable description
      Filedes    griddes         Grid description
      Filedes    vct             Vertical coordinate table
*/


#include <string.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Filedes(void *argument)
{
  int GRIDDES, ZAXISDES, VCT, VARDES, TAXISDES, FILEDES, VLIST, PARTAB;
  int operatorID;
  int streamID = 0;
  int zaxisID;
  int nvars, ngrids, nzaxis;
  int type, index;
  int vlistID;

  cdoInitialize(argument);

  GRIDDES  = cdoOperatorAdd("griddes",  0, 0, NULL);
  ZAXISDES = cdoOperatorAdd("zaxisdes", 0, 0, NULL);
  TAXISDES = cdoOperatorAdd("taxisdes", 0, 0, NULL);
  VCT      = cdoOperatorAdd("vct",      0, 0, NULL);
  VARDES   = cdoOperatorAdd("vardes",   0, 0, NULL);
  FILEDES  = cdoOperatorAdd("filedes",  0, 0, NULL);
  VLIST    = cdoOperatorAdd("vlist",    0, 0, NULL);
  PARTAB   = cdoOperatorAdd("partab",   0, 0, NULL);

  operatorID = cdoOperatorID();

  streamID = streamOpenRead(cdoStreamName(0));
  if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", cdoStreamName(0));

  vlistID = streamInqVlist(streamID);

  nvars  = vlistNvars(vlistID);
  ngrids = vlistNgrids(vlistID);
  nzaxis = vlistNzaxis(vlistID);

  if ( operatorID == GRIDDES )
    {
      for ( index = 0; index < ngrids; index++ )
	gridPrint(vlistGrid(vlistID, index));
    }
  else if ( operatorID == ZAXISDES )
    {
      for ( index = 0; index < nzaxis; index++ )
	zaxisPrint(vlistZaxis(vlistID, index));
    }
  else if ( operatorID == TAXISDES )
    {
      int vdate, vtime, ntsteps, nrecs;
      int year, month, day, hour, minute;
      int taxisID, tsID;

      taxisID = vlistInqTaxis(vlistID);
      ntsteps = vlistNtsteps(vlistID);

      fprintf(stdout, "   Time axis  : ");
      if ( taxisInqType(taxisID) == TAXIS_RELATIVE )
	fprintf(stdout, "relative");
      else if ( taxisInqType(taxisID) == TAXIS_ABSOLUTE )
	fprintf(stdout, "absolute");
      else
	fprintf(stdout, "unknown");
      fprintf(stdout, "\n");

      taxisID = vlistInqTaxis(vlistID);

      if ( ntsteps != 0 )
	{
	  if ( ntsteps == CDI_UNDEFID )
	    fprintf(stdout, "   Time steps :  unlimited\n");
	  else
	    fprintf(stdout, "   Time steps :  %d\n", ntsteps);

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

	  fprintf(stdout, "  time verification time        lower bound        upper bound\n");
	  fprintf(stdout, "  step  YYYY-MM-DD hh:mm   YYYY-MM-DD hh:mm   YYYY-MM-DD hh:mm\n");

	  tsID = 0;
	  while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	    {
	      vdate = taxisInqVdate(taxisID);
	      vtime = taxisInqVtime(taxisID);

	      decode_date(vdate, &year, &month, &day);
	      decode_time(vtime, &hour, &minute);

	      tsID++;
	      fprintf(stdout, " %5d %4.4d-%2.2d-%2.2d %2.2d:%2.2d", tsID, year, month, day, hour, minute);
	      fprintf(stdout, "\n");
	    }
	}
    }
  else if ( operatorID == VCT )
    {
      for ( index = 0; index < nzaxis; index++)
	{
	  zaxisID = vlistZaxis(vlistID, index);
	  type = zaxisInqType(zaxisID);
	  if ( type == ZAXIS_HYBRID )
	    {
	      int i, vctsize;
	      const double *vct;

	      vctsize = zaxisInqVctSize(zaxisID);
	      vct     = zaxisInqVctPtr(zaxisID);
		
	      for ( i = 0; i < vctsize/2; i++ )
		fprintf(stdout, "%5d %25.17f %25.17f\n", i, vct[i], vct[vctsize/2+i]);
	    }
	}
    }
  else if ( operatorID == VLIST )
    {
      vlistPrint(vlistID);
    }
  else if ( operatorID == VARDES )
    {
      int varID, code;
      char varname[128], varlongname[128], varunits[128];

      for ( varID = 0; varID < nvars; varID++ )
	{
	  varname[0]     = 0;
	  varlongname[0] = 0;
	  varunits[0]    = 0;
	  code     = vlistInqVarCode(vlistID, varID);
	  vlistInqVarName(vlistID, varID, varname);
	  vlistInqVarLongname(vlistID, varID, varlongname);
	  vlistInqVarUnits(vlistID, varID, varunits);
	  fprintf(stdout, "%4d  %-12s", code, varname);
	  if ( strlen(varlongname) )
	    {
	      fprintf(stdout, "  %s", varlongname);
	      if ( strlen(varunits) )
		fprintf(stdout, " [%s]", varunits);
	    }
	  fprintf(stdout, "\n");
	}   
    }
  else if ( operatorID == PARTAB )
    {
      int varID, code, tabnum, tableID, prec;
      char pstr[4];
      char varname[128], varlongname[128], varstdname[128], varunits[128];

      for ( varID = 0; varID < nvars; varID++ )
	{
	  fprintf(stdout, "&PARAMETER\n");

	  varname[0]     = 0;
	  varlongname[0] = 0;
	  varunits[0]    = 0;
	  code     = vlistInqVarCode(vlistID, varID);
	  tableID  = vlistInqVarTable(vlistID, varID);
	  tabnum   = tableInqNum(tableID);
	  vlistInqVarName(vlistID, varID, varname);
	  vlistInqVarStdname(vlistID, varID, varstdname);
	  vlistInqVarLongname(vlistID, varID, varlongname);
	  vlistInqVarUnits(vlistID, varID, varunits);

	  prec = vlistInqVarDatatype(vlistID, varID);
	  if      ( prec == DATATYPE_PACK   ) strcpy(pstr, "P0");
	  else if ( prec > 0 && prec <= 32  ) sprintf(pstr, "P%d", prec);
	  else if ( prec == DATATYPE_FLT32  ) strcpy(pstr, "F32");
	  else if ( prec == DATATYPE_FLT64  ) strcpy(pstr, "F64");
	  else if ( prec == DATATYPE_INT8   ) strcpy(pstr, "I8");
	  else if ( prec == DATATYPE_INT16  ) strcpy(pstr, "I16");
	  else if ( prec == DATATYPE_INT32  ) strcpy(pstr, "I32");
	  else                                strcpy(pstr, "-1");

	  if ( code   > 0 ) fprintf(stdout, "  CODE=%d\n", code);
	  if ( tabnum > 0 ) fprintf(stdout, "  TABLE=%d\n", tabnum);
	  fprintf(stdout, "  NAME=%s\n", varname);
	  if ( strlen(varstdname) )
	    fprintf(stdout, "  STANDARD_NAME=%s\n", varstdname);
	  if ( strlen(varlongname) )
	    fprintf(stdout, "  LONG_NAME=\"%s\"\n", varlongname);
	  if ( strlen(varunits) )
	    fprintf(stdout, "  UNITS=\"%s\"\n", varunits);

	  /* if ( pstr ) fprintf(stdout, "  DATATYPE=%s\n", pstr); */

	  fprintf(stdout, "/\n");
	}   
    }
  else if ( operatorID == FILEDES )
    {
      int filetype;

      printf("\n");
      filetype = streamInqFiletype(streamID);
      switch ( filetype )
	{
	case FILETYPE_GRB:
	  printf("  GRIB data\n");
	  break;
	case FILETYPE_NC:
	  printf("  netCDF data\n");
	  break;
	case FILETYPE_NC2:
	  printf("  netCDF2 data\n");
	  break;
	case FILETYPE_SRV:
	  printf("  SERVICE data\n");
	  switch ( streamInqByteorder(streamID) )
	    {
	    case CDI_BIGENDIAN:
	      printf("  byteorder is BIGENDIAN\n"); break;
	    case CDI_LITTLEENDIAN:
	      printf("  byteorder is LITTLEENDIAN\n"); break;
	    default:
	      printf("  byteorder %d undefined\n", streamInqByteorder(streamID)); break;
	    }
	   break;
	case FILETYPE_EXT:
	  printf("  EXTRA data\n");
	  switch ( streamInqByteorder(streamID) )
	    {
	    case CDI_BIGENDIAN:
	      printf("  byteorder is BIGENDIAN\n"); break;
	    case CDI_LITTLEENDIAN:
	      printf("  byteorder is LITTLEENDIAN\n"); break;
	    default:
	      printf("  byteorder %d undefined\n", streamInqByteorder(streamID)); break;
	    }
	   break;
	case FILETYPE_IEG:
	  printf("  IEG data\n");
	  switch ( streamInqByteorder(streamID) )
	    {
	    case CDI_BIGENDIAN:
	      printf("  byteorder is BIGENDIAN\n"); break;
	    case CDI_LITTLEENDIAN:
	      printf("  byteorder is LITTLEENDIAN\n"); break;
	    default:
	      printf("  byteorder %d undefined\n", streamInqByteorder(streamID)); break;
	    }
	   break;
	default:
	  printf("  unsupported filetype %d\n" , filetype);
	}

      printf("\n");
    }

  streamClose(streamID);

  cdoFinish();

  return (0);
}

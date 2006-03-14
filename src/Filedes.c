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

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

/*
@BeginDoc

@BeginModule

@Name      = Filedes
@Title     = File description
@Section   = Information
@Class     = Information
@Arguments = ifile
@Operators = vardes griddes vct

@EndModule


@BeginOperator_vardes

@Title     = Variable description

@BeginDescription
Prints a table with a description of all variables.
For each variable the operator print in one line the
code, name, description and units.
@EndDescription

@EndOperator


@BeginOperator_griddes

@Title     = Grid description

@BeginDescription
Prints the description of all grids in a file.
@EndDescription

@EndOperator

@BeginOperator_vct

@Title     = Vertical coordinate table

@BeginDescription
Prints the vertical coordinate table.
@EndDescription

@EndOperator

@EndDoc
*/

#define  CDI_BIGENDIAN            0   /* Data type BIGENDIAN     */
#define  CDI_LITTLEENDIAN         1   /* Data type LITTLEENDIAN  */

void *Filedes(void *argument)
{
  int GRIDDES, ZAXISDES, VCT, VARDES, TIMEDES, FILEDES, VLIST;
  int operatorID;
  int streamID = 0;
  int zaxisID;
  int nvars, ngrids, nzaxis;
  int type, index;
  int vlistID;

  cdoInitialize(argument);

  GRIDDES  = cdoOperatorAdd("griddes",  0, 0, NULL);
  ZAXISDES = cdoOperatorAdd("zaxisdes", 0, 0, NULL);
  VCT      = cdoOperatorAdd("vct",      0, 0, NULL);
  VARDES   = cdoOperatorAdd("vardes",   0, 0, NULL);
  TIMEDES  = cdoOperatorAdd("timedes",  0, 0, NULL);
  FILEDES  = cdoOperatorAdd("filedes",  0, 0, NULL);
  VLIST    = cdoOperatorAdd("vlist",    0, 0, NULL);

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
	default:
	  printf("  unsupported filetype %d\n" , filetype);
	}

      printf("\n");
    }

  streamClose(streamID);

  cdoFinish();

  return (0);
}

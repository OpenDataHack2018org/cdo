/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2012 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


static
const char *filetypestr(int filetype)
{
  switch ( filetype )
    {
    case FILETYPE_GRB:  return ("GRIB");            break;
    case FILETYPE_GRB2: return ("GRIB2");           break;
    case FILETYPE_NC:   return ("netCDF");          break;
    case FILETYPE_NC2:  return ("netCDF2");         break;
    case FILETYPE_NC4:  return ("netCDF4");         break;
    case FILETYPE_NC4C: return ("netCDF4 classic"); break;
    case FILETYPE_SRV:  return ("SERVICE");         break;
    case FILETYPE_EXT:  return ("EXTRA");           break;
    case FILETYPE_IEG:  return ("IEG");             break;
    default:            return ("");
    }
}

static
const char *datatypestr(int datatype)
{
  static char str[20];

  str[0] = 0;
  sprintf(str, "%d bit packed", datatype);

  if      ( datatype == DATATYPE_PACK   ) return ("P0");
  else if ( datatype > 0 && datatype <= 32 ) return (str);
  else if ( datatype == DATATYPE_CPX32  ) return ("C32");
  else if ( datatype == DATATYPE_CPX64  ) return ("C64");
  else if ( datatype == DATATYPE_FLT32  ) return ("32 bit floats");
  else if ( datatype == DATATYPE_FLT64  ) return ("64 bit floats");
  else if ( datatype == DATATYPE_INT8   ) return ("I8");
  else if ( datatype == DATATYPE_INT16  ) return ("I16");
  else if ( datatype == DATATYPE_INT32  ) return ("I32");
  else if ( datatype == DATATYPE_UINT8  ) return ("U8");
  else if ( datatype == DATATYPE_UINT16 ) return ("U16");
  else if ( datatype == DATATYPE_UINT32 ) return ("U32");
  else                                    return ("");
}

static
off_t filesize(const char *filename)
{
  FILE *fp;
  off_t pos = 0;

  fp = fopen(filename, "r");
  if ( fp == NULL )
    {
      fprintf(stderr, "Open failed on %s\n", filename);
    }
  else
    {
      fseek(fp, 0L, SEEK_END);
      pos = ftello(fp);
    }
  
  return pos;
}

void *CDIwrite(void *argument)
{
  int nvars = 10, nlevs = 0, ntimesteps = 30;
  char *defaultgrid = "global_.2";
  int operatorID;
  int streamID;
  int tsID, varID, levelID;
  int gridsize, i;
  int rval, rstart, rinc;
  int vlistID;
  int gridID = -1, zaxisID, taxisID;
  int vdate, vtime, julday;
  int filetype, datatype;
  off_t fsize;
  unsigned int seed = 1;
  const char *gridfile;
  double file_size, data_size = 0;
  double tw;
  double *levels = NULL;
  double ***vars = NULL;
  extern int timer_write;

  srand(seed);

  cdoInitialize(argument);

  if ( cdoVerbose ) cdoPrint("parameter: <grid, <nlevs, <ntimesteps, <nvars>>>>");
  // operatorInputArg("<grid, <nlevs, <ntimesteps, <nvars>>>>");

  if ( operatorArgc() > 4 ) cdoAbort("Too many arguments!");

  gridfile = defaultgrid;
  if ( operatorArgc() >= 1 ) gridfile = operatorArgv()[0];
  if ( operatorArgc() >= 2 ) nlevs = atol(operatorArgv()[1]);
  if ( operatorArgc() >= 3 ) ntimesteps = atol(operatorArgv()[2]);
  if ( operatorArgc() >= 4 ) nvars = atol(operatorArgv()[3]);

  gridID   = cdoDefineGrid(gridfile);
  gridsize = gridInqSize(gridID);

  zaxisID  = zaxisCreate(ZAXIS_SURFACE, 1);

  if ( cdoVerbose )
    {
      cdoPrint("gridsize   : %d", gridInqSize);
      cdoPrint("nlevs      : %d", nlevs);
      cdoPrint("ntimesteps : %d", ntimesteps);
      cdoPrint("nvars      : %d", nvars);
    } 

  if ( nlevs <= 0 ) nlevs = 1;
  if ( ntimesteps <= 0 ) ntimesteps = 1;
  if ( nvars <= 0 ) nvars = 1;

  vars = (double ***) malloc(nvars*sizeof(double **));
  for ( varID = 0; varID < nvars; varID++ )
    {
      vars[varID] = (double **) malloc(nlevs*sizeof(double *));
      for ( levelID = 0; levelID < nlevs; levelID++ )
	{
	  vars[varID][levelID] = (double *) malloc(gridsize*sizeof(double));
	  for ( i = 0; i < gridsize; ++i )
	    vars[varID][levelID][i] = varID + rand()/(RAND_MAX+1.0);
	}
    }

  vlistID = vlistCreate();

  for ( i = 0; i < nvars; ++i )
    {
      varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARIABLE);
      vlistDefVarCode(vlistID, varID, varID+1);
      //    vlistDefVarName(vlistID, varID, );
    }

  taxisID = taxisCreate(TAXIS_RELATIVE);
  vlistDefTaxis(vlistID, taxisID);

  // vlistDefNtsteps(vlistID, 1);

  streamID = streamOpenWrite(cdoStreamName(0), cdoFiletype());

  streamDefVlist(streamID, vlistID);

  filetype = streamInqFiletype(streamID);
  datatype = vlistInqVarDatatype(vlistID, 0);
	  
  julday = date_to_julday(CALENDAR_PROLEPTIC, 19870101);

  for ( tsID = 0; tsID < ntimesteps; tsID++ )
    {
      rval  = rstart + rinc*tsID;
      vdate = julday_to_date(CALENDAR_PROLEPTIC, julday + tsID);
      vtime = 0;
      taxisDefVdate(taxisID, vdate);
      taxisDefVtime(taxisID, vtime);
      streamDefTimestep(streamID, tsID);

      for ( varID = 0; varID < nvars; varID++ )
        {
          for ( levelID = 0; levelID < nlevs; levelID++ )
            {
              streamDefRecord(streamID, varID, levelID);
              streamWriteRecord(streamID, vars[varID][levelID], 0);
	      data_size += gridsize*8;
            }
        }
    }

  streamClose(streamID);

  tw = timer_val(timer_write);

  data_size /= 1024.*1024.*1024.;
  cdoPrint("Wrote %.1f GB of 64 bit floats to %s %s", data_size, datatypestr(datatype), filetypestr(filetype));

  fsize = filesize(cdoStreamName(0));
  file_size = fsize;
  file_size /= 1024.*1024.*1024.;
  cdoPrint("Wrote %.1f GB in %.1f seconds, total %.1f MB/s", file_size, tw, 1024*file_size/tw);

  vlistDestroy(vlistID);

  for ( varID = 0; varID < nvars; varID++ )
    {
      free(vars[varID]);
      for ( levelID = 0; levelID < nlevs; levelID++ ) free(vars[varID][levelID]);
    }
  free(vars);

  cdoFinish();

  return (0);
}

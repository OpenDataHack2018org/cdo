/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2018 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
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

#include "cdo_int.h"
#include "pstream_int.h"

#define NLON 1440
#define NLAT 720
#define MAX_VARS 6

static void
init_amsr_day(int vlistID, int gridID, int zaxisID, int nvars)
{
  /*
    Version-5 RSS AMSR-E or AMSR-J daily files

    filename  with path in form satname_yyyymmdd_v5.gz
    where satname  = name of satellite (amsre or amsr)
             yyyy	= year
               mm	= month
               dd	= day of month

    1:time	time of measurement in fractional hours GMT
    2:sst     	sea surface temperature in deg Celcius
    3:wind	10m surface wind in meters/second
    4:vapor	columnar water vapor in millimeters
    5:cloud	cloud liquid water in millimeters
    6:rain   	rain rate in millimeters/hour
  */
  const char *name[] = { "hours", "sst", "wind", "vapor", "cloud", "rain" };
  const char *units[] = { "h", "deg Celcius", "m/s", "mm", "mm", "mm/h" };
  double xscale[] = { 0.1, 0.15, 0.2, 0.3, 0.01, 0.1 };
  double xminval[] = { 0., -3., 0., 0., 0., 0. };
  int i, varID;

  for (i = 0; i < nvars; ++i)
    {
      varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARYING);
      vlistDefVarName(vlistID, varID, name[i]);
      vlistDefVarUnits(vlistID, varID, units[i]);
      vlistDefVarDatatype(vlistID, varID, CDI_DATATYPE_INT16);
      vlistDefVarMissval(vlistID, varID, 254);
      vlistDefVarScalefactor(vlistID, varID, xscale[i]);
      vlistDefVarAddoffset(vlistID, varID, xminval[i]);
    }
}

static void
init_amsr_averaged(int vlistID, int gridID, int zaxisID, int nvars)
{
  /*
    Version-5 AMSR-E or AMSR-J time-averaged files including:
          3-day		(average of 3 days ending on file date)
          weekly	(average of 7 days ending on Saturday of file date)
          monthly	(average of all days in month)


     filename
           format of file names are:
                3-day	satname_yyyymmddv5_d3d.gz
                weekly	satname_yyyymmddv5.gz
                monthly	satname_yyyymmv5.gz

        where	satname	=name of satellite (amsre or amsr)
                        yyyy	=year
                        mm	=month
                        dd	=day of month

    1:sst       sea surface temperature in deg Celcius
    2:wind      10m surface wind in meters/second
    3:vapor	columnar water vapor in millimeters
    4:cloud     cloud liquid water in millimeters
    5:rain	rain rate in millimeters/hour
  */
  const char *name[] = { "sst", "wind", "vapor", "cloud", "rain" };
  const char *units[] = { "deg Celcius", "m/s", "mm", "mm", "mm/h" };
  double xscale[] = { 0.15, 0.2, 0.3, 0.01, 0.1 };
  double xminval[] = { -3., 0., 0., 0., 0. };
  int i, varID;

  for (i = 0; i < nvars; ++i)
    {
      /* varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_CONSTANT); */
      varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARYING);
      vlistDefVarName(vlistID, varID, name[i]);
      vlistDefVarUnits(vlistID, varID, units[i]);
      vlistDefVarDatatype(vlistID, varID, CDI_DATATYPE_INT16);
      vlistDefVarMissval(vlistID, varID, 254);
      vlistDefVarScalefactor(vlistID, varID, xscale[i]);
      vlistDefVarAddoffset(vlistID, varID, xminval[i]);
    }
}

static void
read_amsr(FILE *fp, int vlistID, int nvars, double *data[], size_t *nmiss)
{
  int varID, i, gridsize;
  unsigned char *amsr_data = NULL;
  double xminval, xscale, missval;
  size_t size;

  for (varID = 0; varID < nvars; ++varID)
    {
      gridsize = gridInqSize(vlistInqVarGrid(vlistID, varID));
      amsr_data = (unsigned char *) Realloc(amsr_data, gridsize);
      size = fread(amsr_data, 1, gridsize, fp);
      if ((int) size != gridsize) cdoAbort("Read error!");

      missval = vlistInqVarMissval(vlistID, varID);
      xminval = vlistInqVarAddoffset(vlistID, varID);
      xscale = vlistInqVarScalefactor(vlistID, varID);

      nmiss[varID] = 0;
      for (i = 0; i < gridsize; ++i)
        {
          if (amsr_data[i] <= 250)
            {
              data[varID][i] = amsr_data[i] * xscale + xminval;
            }
          else
            {
              data[varID][i] = missval;
              nmiss[varID]++;
            }
        }
    }

  Free(amsr_data);
}

static void
write_data(int streamID, int nvars, double *data[], size_t *nmiss)
{
  for (int varID = 0; varID < nvars; ++varID)
    {
      pstreamDefRecord(streamID, varID, 0);
      pstreamWriteRecord(streamID, data[varID], nmiss[varID]);
    }
}

static int
getDate(const char *name)
{
  int date = 0;
  char *pname = (char *) strchr(name, '_');

  if (pname) date = atoi(pname + 1);

  return date;
}

void *
Importamsr(void *process)
{
  int tsID;
  int nvars;
  int vtime = 0;
  double xvals[NLON], yvals[NLAT];
  double *data[MAX_VARS];
  size_t nmiss[MAX_VARS];

  cdoInitialize(process);

  FILE *fp = fopen(cdoGetStreamName(0).c_str(), "r");
  if (fp == NULL)
    {
      perror(cdoGetStreamName(0).c_str());
      exit(EXIT_FAILURE);
    }

  fseek(fp, 0L, SEEK_END);
  size_t fsize = (size_t) ftell(fp);
  fseek(fp, 0L, SEEK_SET);

  int64_t vdate = getDate(cdoGetStreamName(0).c_str());
  if (vdate <= 999999) vdate = vdate * 100 + 1;

  int streamID = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());

  /*
    Longitude  is 0.25*xdim-0.125    degrees east
    Latitude   is 0.25*ydim-90.125
  */
  size_t gridsize = NLON * NLAT;
  int gridID = gridCreate(GRID_LONLAT, gridsize);
  gridDefXsize(gridID, NLON);
  gridDefYsize(gridID, NLAT);
  for (int i = 0; i < NLON; ++i) xvals[i] = 0.25 * (i + 1) - 0.125;
  for (int i = 0; i < NLAT; ++i) yvals[i] = 0.25 * (i + 1) - 90.125;
  gridDefXvals(gridID, xvals);
  gridDefYvals(gridID, yvals);

  int zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);

  int taxisID = taxisCreate(TAXIS_ABSOLUTE);

  int vlistID = vlistCreate();
  vlistDefTaxis(vlistID, taxisID);

  if (fsize == 12441600)
    {
      nvars = 6;
      for (int i = 0; i < nvars; ++i) data[i] = (double *) Malloc(gridsize * sizeof(double));

      init_amsr_day(vlistID, gridID, zaxisID, nvars);

      pstreamDefVlist(streamID, vlistID);

      vtime = 13000; /* 1:30:00 */
      for (tsID = 0; tsID < 2; ++tsID)
        {
          taxisDefVdate(taxisID, vdate);
          taxisDefVtime(taxisID, vtime);
          vtime += 120000; /* 13:30:00 */
          pstreamDefTimestep(streamID, tsID);

          read_amsr(fp, vlistID, nvars, data, nmiss);

          write_data(streamID, nvars, data, nmiss);
        }

      for (int i = 0; i < nvars; ++i) Free(data[i]);
    }
  else if (fsize == 5184000)
    {
      nvars = 5;
      for (int i = 0; i < nvars; ++i) data[i] = (double *) Malloc(gridsize * sizeof(double));

      init_amsr_averaged(vlistID, gridID, zaxisID, nvars);

      /* vlistDefNtsteps(vlistID, 0);*/
      pstreamDefVlist(streamID, vlistID);

      taxisDefVdate(taxisID, vdate);
      taxisDefVtime(taxisID, vtime);
      tsID = 0;
      pstreamDefTimestep(streamID, tsID);

      read_amsr(fp, vlistID, nvars, data, nmiss);

      write_data(streamID, nvars, data, nmiss);

      for (int i = 0; i < nvars; ++i) Free(data[i]);
    }
  else
    cdoAbort("Unexpected file size for AMSR data!");

  processDefVarNum(vlistNvars(vlistID));

  pstreamClose(streamID);

  fclose(fp);

  vlistDestroy(vlistID);
  gridDestroy(gridID);
  zaxisDestroy(zaxisID);
  taxisDestroy(taxisID);

  cdoFinish();

  return 0;
}

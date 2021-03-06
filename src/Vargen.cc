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

/*
   This module contains the following operators:

      Vargen     const           Create a constant field
      Vargen     random          Field with random values
      Vargen     stdatm          Field values for pressure and temperature for
                                 the standard atmosphere
*/

#ifdef HAVE_CONFIG_H
#include "config.h"  // ENABLE_DATA
#endif

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "listarray.h"
#include "grid.h"
#include "constants.h"
#include "stdnametable.h"

#if defined(ENABLE_DATA)
static const double etopo_scale = 3;
static const double etopo_offset = 11000;
static const unsigned short etopo[] = {
#include "etopo.h"
};

static const double temp_scale = 500;
static const double temp_offset = -220;
static const unsigned short temp[] = {
#include "temp.h"
};

static const double mask_scale = 1;
static const double mask_offset = 0;
static const unsigned short mask[] = {
#include "mask.h"
};
#endif

/*  some Constants for creating temperatur and pressure for the standard
 * atmosphere */
#define T_ZERO (213.0)
#define T_DELTA (75.0)
#define SCALEHEIGHT (10000.0) /* [m] */
#define P_ZERO (1013.25)      /* surface pressure [hPa] */
#define CC_R (287.05)         /* specific gas constant for air */
static double TMP4PRESSURE = (C_EARTH_GRAV * SCALEHEIGHT) / (CC_R * T_ZERO);

static double
std_atm_temperatur(double height)
{
  /*
    Compute the temperatur for the given height (in meters) according to the
    solution of the hydrostatic atmosphere
   */
  return (T_ZERO + T_DELTA * exp((-1) * (height / SCALEHEIGHT)));
}

static double
std_atm_pressure(double height)
{
  /*
    Compute the pressure for the given height (in meters) according to the
    solution of the hydrostatic atmosphere
   */
  return (P_ZERO * exp((-1) * TMP4PRESSURE * log((exp(height / SCALEHEIGHT) * T_ZERO + T_DELTA) / (T_ZERO + T_DELTA))));
}

static void
conv_generic_grid(int gridID, size_t gridsize, double *xvals2D, double *yvals2D)
{
  size_t xsize = gridInqXsize(gridID);
  size_t ysize = gridInqYsize(gridID);

  assert(gridsize == xsize * ysize);

  std::vector<double> xcoord(xsize);
  std::vector<double> ycoord(ysize);

  gridInqXvals(gridID, &xcoord[0]);
  gridInqYvals(gridID, &ycoord[0]);

  double xmin, xmax;
  arrayMinMax(xsize, &xcoord[0], &xmin, &xmax);

  double ymin, ymax;
  arrayMinMax(ysize, &ycoord[0], &ymin, &ymax);

  double xrange = xmax - xmin;
  double yrange = ymax - ymin;

  for (size_t j = 0; j < ysize; ++j)
    for (size_t i = 0; i < xsize; ++i)
      {
        xvals2D[j * xsize + i] = xcoord[i] * M_PI / xrange;
        yvals2D[j * xsize + i] = ycoord[j] * M_PI / yrange;
      }
}

static void
remap_nn_reg2d_reg2d(size_t nx, size_t ny, const double *restrict data, int gridID, double *restrict array)
{
  if (gridInqType(gridID) != GRID_LONLAT) cdoAbort("Internal error, wrong grid type!");

  size_t nxvals = gridInqXsize(gridID);
  size_t nyvals = gridInqYsize(gridID);
  std::vector<double> xvals(nxvals);
  std::vector<double> yvals(nyvals);

  gridInqXvals(gridID, &xvals[0]);
  gridInqYvals(gridID, &yvals[0]);

  /* Convert lat/lon units if required */
  char units[CDI_MAX_NAME];
  gridInqXunits(gridID, units);
  grid_to_degree(units, nxvals, &xvals[0], "grid center lon");
  gridInqYunits(gridID, units);
  grid_to_degree(units, nyvals, &yvals[0], "grid center lat");

  for (size_t j = 0; j < nyvals; j++)
    {
      double yval = yvals[j];
      for (size_t i = 0; i < nxvals; i++)
        {
          double xval = xvals[i];
          if (xval >= 180) xval -= 360;
          if (xval < -180) xval += 360;
          size_t ii = (xval + 180) * 2;
          size_t jj = (yval + 90) * 2;
          if (ii >= nx) ii = nx - 1;
          if (jj >= ny) jj = ny - 1;
          array[j * nxvals + i] = data[jj * nx + ii];
        }
    }
}

static void
remap_nn_reg2d_nonreg2d(size_t nx, size_t ny, const double *restrict data, int gridID, double *restrict array)
{
  int gridID2 = gridID;
  size_t gridsize = gridInqSize(gridID2);
  std::vector<double> xvals(gridsize);
  std::vector<double> yvals(gridsize);

  if (gridInqType(gridID2) == GRID_GME) gridID2 = gridToUnstructured(gridID2, 0);

  if (gridInqType(gridID2) != GRID_UNSTRUCTURED && gridInqType(gridID2) != GRID_CURVILINEAR)
    gridID2 = gridToCurvilinear(gridID2, 0);

  gridInqXvals(gridID2, &xvals[0]);
  gridInqYvals(gridID2, &yvals[0]);

  /* Convert lat/lon units if required */
  char units[CDI_MAX_NAME];
  gridInqXunits(gridID2, units);
  grid_to_degree(units, gridsize, &xvals[0], "grid center lon");
  gridInqYunits(gridID2, units);
  grid_to_degree(units, gridsize, &yvals[0], "grid center lat");

  for (size_t i = 0; i < gridsize; i++)
    {
      double xval = xvals[i];
      double yval = yvals[i];
      if (xval >= 180) xval -= 360;
      if (xval < -180) xval += 360;
      size_t ii = (xval + 180) * 2;
      size_t jj = (yval + 90) * 2;
      if (ii >= nx) ii = nx - 1;
      if (jj >= ny) jj = ny - 1;
      array[i] = data[jj * nx + ii];
    }

  if (gridID != gridID2) gridDestroy(gridID2);
}

static void
remap_nn_reg2d(size_t nx, size_t ny, const double *restrict data, int gridID, double *restrict array)
{
  if (gridInqType(gridID) == GRID_LONLAT)
    remap_nn_reg2d_reg2d(nx, ny, data, gridID, array);
  else
    remap_nn_reg2d_nonreg2d(nx, ny, data, gridID, array);
}

#define NLON 720
#define NLAT 360

void *
Vargen(void *process)
{
  int ntimesteps, nlevels = 1;
  int varID, varID2 = -1, levelID;
  int gridID = -1, gridIDdata = -1, zaxisID;
  double rstart = 0, rstop = 0, rinc = 0;
  double rconst = 0;
  double *levels = NULL;
  double lon[NLON], lat[NLAT];
  size_t nlon = NLON;
  size_t nlat = NLAT;
  ListArray<double> listArrayFlt;

  cdoInitialize(process);

  // clang-format off
  int RANDOM  = cdoOperatorAdd("random",  0, 0, "grid description file or name, <seed>");
  int SINCOS  = cdoOperatorAdd("sincos",  0, 0, "grid description file or name");
  int COSHILL = cdoOperatorAdd("coshill", 0, 0, "grid description file or name");
  int CONST   = cdoOperatorAdd("const",   0, 0, "constant value, grid description file or name");
  int FOR     = cdoOperatorAdd("for",     0, 0, "start, end, <increment>");
  int TOPO    = cdoOperatorAdd("topo",    0, 0, NULL);
  int TEMP    = cdoOperatorAdd("temp",    0, 0, NULL);
  int MASK    = cdoOperatorAdd("mask",    0, 0, NULL);
  int STDATM  = cdoOperatorAdd("stdatm",  0, 0, "height levels [m]");
  // clang-format on

  int operatorID = cdoOperatorID();

  if (operatorID == RANDOM)
    {
      unsigned int seed = 1;
      operatorInputArg(cdoOperatorEnter(operatorID));
      if (operatorArgc() < 1) cdoAbort("Too few arguments!");
      if (operatorArgc() > 2) cdoAbort("Too many arguments!");
      gridID = cdoDefineGrid(operatorArgv()[0]);
      if (operatorArgc() == 2)
        {
          int idum = parameter2int(operatorArgv()[1]);
          if (idum >= 0 && idum < 0x7FFFFFFF) seed = idum;
        }
      srand(seed);
    }
  else if (operatorID == SINCOS || operatorID == COSHILL)
    {
      operatorInputArg(cdoOperatorEnter(operatorID));
      operatorCheckArgc(1);
      gridID = cdoDefineGrid(operatorArgv()[0]);
    }
  else if (operatorID == CONST)
    {
      operatorInputArg(cdoOperatorEnter(operatorID));
      operatorCheckArgc(2);
      rconst = parameter2double(operatorArgv()[0]);
      gridID = cdoDefineGrid(operatorArgv()[1]);
    }
  else if (operatorID == TOPO || operatorID == TEMP || operatorID == MASK)
    {
      gridIDdata = gridCreate(GRID_LONLAT, nlon * nlat);
      gridDefXsize(gridIDdata, nlon);
      gridDefYsize(gridIDdata, nlat);

      for (size_t i = 0; i < nlon; i++) lon[i] = -179.75 + i * 0.5;
      for (size_t i = 0; i < nlat; i++) lat[i] = -89.75 + i * 0.5;

      gridDefXvals(gridIDdata, lon);
      gridDefYvals(gridIDdata, lat);

      gridID = gridIDdata;

      if (operatorArgc() == 1) gridID = cdoDefineGrid(operatorArgv()[0]);
    }
  else if (operatorID == FOR)
    {
      double lon = 0, lat = 0;
      operatorInputArg(cdoOperatorEnter(operatorID));
      if (operatorArgc() < 2) cdoAbort("Too few arguments!");
      if (operatorArgc() > 3) cdoAbort("Too many arguments!");

      rstart = parameter2double(operatorArgv()[0]);
      rstop = parameter2double(operatorArgv()[1]);
      if (operatorArgc() == 3)
        rinc = parameter2double(operatorArgv()[2]);
      else
        rinc = 1;

      if (DBL_IS_EQUAL(rinc, 0.0)) cdoAbort("Increment is zero!");

      gridID = gridCreate(GRID_LONLAT, 1);
      gridDefXsize(gridID, 1);
      gridDefYsize(gridID, 1);
      gridDefXvals(gridID, &lon);
      gridDefYvals(gridID, &lat);
    }
  else if (operatorID == STDATM)
    {
      double lon = 0, lat = 0;

      operatorInputArg(cdoOperatorEnter(operatorID));
      nlevels = listArrayFlt.argvToFlt(operatorArgc(), operatorArgv());
      levels = listArrayFlt.data();

      if (cdoVerbose)
        for (int i = 0; i < nlevels; ++i) printf("levels %d: %g\n", i, levels[i]);

      gridID = gridCreate(GRID_LONLAT, 1);
      gridDefXsize(gridID, 1);
      gridDefYsize(gridID, 1);
      gridDefXvals(gridID, &lon);
      gridDefYvals(gridID, &lat);
    }

  if (operatorID == STDATM)
    {
      zaxisID = zaxisCreate(ZAXIS_HEIGHT, nlevels);
      zaxisDefLevels(zaxisID, levels);
      zaxisDefName(zaxisID, "level");
      zaxisDefLongname(zaxisID, "Level");
      zaxisDefUnits(zaxisID, "m");
    }
  else
    {
      zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);
      nlevels = 1;
    }

  int vlistID = vlistCreate();

  int timetype = (operatorID == FOR) ? TIME_VARYING : TIME_CONSTANT;

  varID = vlistDefVar(vlistID, gridID, zaxisID, timetype);
  /*
     For the standard atmosphere two output variables are generated: pressure
     and temperatur. The first (varID) is pressure, second (varID2) is
     temperatur. Add an additional variable for the standard atmosphere.
   */
  if (operatorID == STDATM) varID2 = vlistDefVar(vlistID, gridID, zaxisID, TIME_CONSTANT);

  if (operatorID == MASK) vlistDefVarDatatype(vlistID, varID, CDI_DATATYPE_INT8);

  if (operatorID == STDATM)
    {
      vlistDefVarName(vlistID, varID, "P");
      vlistDefVarParam(vlistID, varID, cdiEncodeParam(1, 255, 255));
      vlistDefVarStdname(vlistID, varID, "air_pressure");
      vlistDefVarLongname(vlistID, varID, "pressure");
      vlistDefVarUnits(vlistID, varID, "hPa");

      vlistDefVarName(vlistID, varID2, "T");
      vlistDefVarParam(vlistID, varID2, cdiEncodeParam(130, 128, 255));
      vlistDefVarStdname(vlistID, varID2, var_stdname(air_temperature));
      vlistDefVarLongname(vlistID, varID2, "temperature");
      vlistDefVarUnits(vlistID, varID2, "K");
    }
  else
    {
      vlistDefVarName(vlistID, varID, cdoOperatorName(operatorID));
      if (operatorID == TOPO) vlistDefVarUnits(vlistID, varID, "m");
      if (operatorID == TEMP) vlistDefVarUnits(vlistID, varID, "K");
    }

  int taxisID = taxisCreate(TAXIS_RELATIVE);
  vlistDefTaxis(vlistID, taxisID);

  if (operatorID == RANDOM || operatorID == SINCOS || operatorID == COSHILL || operatorID == CONST || operatorID == TOPO
      || operatorID == TEMP || operatorID == MASK || operatorID == STDATM)
    vlistDefNtsteps(vlistID, 1);

  int streamID = cdoStreamOpenWrite(0, cdoFiletype());

  pstreamDefVlist(streamID, vlistID);

  size_t gridsize = gridInqSize(gridID);
  size_t datasize = gridsize;
  std::vector<double> array(gridsize);
  double *data = &array[0];
  if (gridID != gridIDdata && gridIDdata != -1)
    {
      datasize = gridInqSize(gridIDdata);
      data = (double *) Malloc(datasize * sizeof(double));
    }

  if (operatorID == FOR)
    ntimesteps = 1.001 + ((rstop - rstart) / rinc);
  else
    {
      vlistDefNtsteps(vlistID, 0);
      ntimesteps = 1;
    }

  int64_t julday = date_to_julday(CALENDAR_PROLEPTIC, 10101);

  int nvars = vlistNvars(vlistID);

  for (int tsID = 0; tsID < ntimesteps; tsID++)
    {
      double rval = rstart + rinc * tsID;
      int64_t vdate = julday_to_date(CALENDAR_PROLEPTIC, julday + tsID);
      int vtime = 0;
      taxisDefVdate(taxisID, vdate);
      taxisDefVtime(taxisID, vtime);
      pstreamDefTimestep(streamID, tsID);

      for (varID = 0; varID < nvars; varID++)
        {
          nlevels = zaxisInqSize(vlistInqVarZaxis(vlistID, varID));
          for (levelID = 0; levelID < nlevels; levelID++)
            {
              pstreamDefRecord(streamID, varID, levelID);

              if (operatorID == RANDOM)
                {
                  for (size_t i = 0; i < gridsize; i++) array[i] = ((double) rand()) / ((double) RAND_MAX);
                }
              else if (operatorID == SINCOS || operatorID == COSHILL)
                {
                  double *xvals = (double *) Malloc(gridsize * sizeof(double));
                  double *yvals = (double *) Malloc(gridsize * sizeof(double));

                  if (grid_is_distance_generic(gridID))
                    {
                      conv_generic_grid(gridID, gridsize, xvals, yvals);
                    }
                  else
                    {
                      if (gridInqType(gridID) == GRID_GME) gridID = gridToUnstructured(gridID, 0);

                      if (gridInqType(gridID) != GRID_UNSTRUCTURED && gridInqType(gridID) != GRID_CURVILINEAR)
                        gridID = gridToCurvilinear(gridID, 0);

                      gridInqXvals(gridID, xvals);
                      gridInqYvals(gridID, yvals);

                      /* Convert lat/lon units if required */
                      char units[CDI_MAX_NAME];
                      gridInqXunits(gridID, units);
                      grid_to_radian(units, gridsize, xvals, "grid center lon");
                      gridInqYunits(gridID, units);
                      grid_to_radian(units, gridsize, yvals, "grid center lat");
                    }

                  if (operatorID == SINCOS)
                    {
                      for (size_t i = 0; i < gridsize; i++) array[i] = cos(1.0 * xvals[i]) * sin(2.0 * yvals[i]);
                    }
                  else if (operatorID == COSHILL)
                    {
                      for (size_t i = 0; i < gridsize; i++) array[i] = 2 - cos(acos(cos(xvals[i]) * cos(yvals[i])) / 1.2);
                    }

                  Free(xvals);
                  Free(yvals);
                }
              else if (operatorID == CONST)
                {
                  for (size_t i = 0; i < gridsize; i++) array[i] = rconst;
                }
              else if (operatorID == TOPO)
                {
#if defined(ENABLE_DATA)
                  for (size_t i = 0; i < datasize; i++) data[i] = etopo[i] / etopo_scale - etopo_offset;
#else
                  cdoAbort("Operator support disabled!");
#endif
                }
              else if (operatorID == TEMP)
                {
#if defined(ENABLE_DATA)
                  for (size_t i = 0; i < datasize; i++) data[i] = temp[i] / temp_scale - temp_offset;
#else
                  cdoAbort("Operator support disabled!");
#endif
                }
              else if (operatorID == MASK)
                {
#if defined(ENABLE_DATA)
                  for (size_t i = 0; i < datasize; i++) data[i] = mask[i] / mask_scale - mask_offset;
#else
                  cdoAbort("Operator support disabled!");
#endif
                }
              else if (operatorID == FOR)
                {
                  array[0] = rval;
                }
              else if (operatorID == STDATM)
                {
                  array[0] = (varID == varID2) ? std_atm_temperatur(levels[levelID]) : std_atm_pressure(levels[levelID]);
                }

              if (gridID != gridIDdata && (operatorID == TOPO || operatorID == TEMP || operatorID == MASK))
                {
                  remap_nn_reg2d(nlon, nlat, data, gridID, &array[0]);
                }

              pstreamWriteRecord(streamID, &array[0], 0);
            }
        }
    }

  pstreamClose(streamID);

  vlistDestroy(vlistID);

  if (gridID != gridIDdata && gridIDdata != -1) Free(data);

  cdoFinish();

  return 0;
}

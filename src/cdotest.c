/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2006 Brockmann Consult
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "cdi.h"
#include "cdo_int.h"


#define TO_KELVIN(x) ((x) + 273.15)
#define MISSVAL -9.0E33F


static int equals(double expected, double actual, double eps)
{
  return fabs(expected - actual) < eps; 
}


static double humidityIndex(double t, double p, double h, double missval)
{
  static const double tmin = 26.0;
  static const double hmin = 0.40;
  
  if ( t < tmin || h < hmin )
    return missval;
    
  return t + (5.0 / 9.0) * ((0.01 * h * p * 6.112 * pow(10.0, (7.5 * t) / (237.7 + t))) - 10.0);
}


/* reads a NetCDF file containing data for a single grid point */ 
static void readNcFile(const char path[], double **vars, int nvars, int nts)
{
  int taxisID;
  int vlistID, varID, streamID, tsID;
  int nmiss, nrecs;
  
  streamID = streamOpenRead(path);
  if ( streamID < 0 )
    {
      fprintf(stderr, "%s\n", cdiStringError(streamID));
      exit(EXIT_FAILURE);
    }
  
  vlistID = streamInqVlist(streamID);
  
  assert(streamNtsteps(streamID) == nts);
  
  assert(vlistGridsizeMax(vlistID) == 1);
  assert(vlistNvars(vlistID) == nvars);

  taxisID = vlistInqTaxis(vlistID);
  
  for (tsID = 0; tsID < nts; ++tsID)
    {
      nrecs = streamInqTimestep(streamID, tsID);
      assert(nrecs == nvars);
      
      taxisInqVdate(taxisID);
      taxisInqVtime(taxisID);
    
      for (varID = 0; varID < nvars; ++varID)
        streamReadVar(streamID, varID, &vars[varID][tsID], &nmiss);
    }
  
  streamClose(streamID);
}


/* writes a NetCDF file containing data for a single grid point */
static void writeNcFile(const char path[], const double array[], int length)
{
  int gridID, zaxisID, taxisID;
  int vlistID, varID, streamID, tsID;
  int nmiss;
  
  double lons[] = {0.0};
  double lats[] = {0.0};
  double value;
  
  gridID = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID, 1);
  gridDefYsize(gridID, 1);
  gridDefXvals(gridID, lons);
  gridDefYvals(gridID, lats);
  
  zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);
  vlistID = vlistCreate();
  
  varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARIABLE);
  vlistDefVarName(vlistID, varID, "test_values");
  vlistDefVarMissval(vlistID, varID, MISSVAL);
  
  taxisID = taxisCreate(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID, taxisID);
  
  streamID = streamOpenWrite(path, FILETYPE_NC);
  if ( streamID < 0 ) 
    {
      fprintf(stderr, "%s\n", cdiStringError(streamID));
      exit(EXIT_FAILURE);
    }
  
  streamDefVlist(streamID, vlistID);
  
  for (tsID = 0; tsID < length; ++tsID)
    {
      taxisDefVdate(taxisID, 20060101 + tsID);
      taxisDefVtime(taxisID, 2359);
      streamDefTimestep(streamID, tsID);
      
      value = array[tsID];
      nmiss = DBL_IS_EQUAL(value, MISSVAL) ? 1 : 0;
      
      streamWriteVar(streamID, varID, &value, nmiss);
    }
    
  streamClose(streamID);
  
  vlistDestroy(vlistID);
  taxisDestroy(taxisID);
  zaxisDestroy(zaxisID);
  gridDestroy(gridID); 
}


/* creates the necessary storage for nvars variables with nts time steps on a single grid point */
static double **createVars(int nvars, int nts)
{
  static const char func[] = "createVars";
  
  double *array = (double *) malloc(nvars*nts*sizeof(double));
  double **vars = (double **) malloc(nvars*sizeof(double));
  
  int i;
  
  for (i = 0; i < nvars; ++i)
    vars[i] = &array[i * nts];
    
  return vars;  
}


/* destroys storage */
static void destroyVars(double **vars)
{
  static const char func[] = "destroyVars";
  
  if ( vars != NULL)
    {
      free(vars[0]);
      free(vars);
    }
}


/* gets the path of the CDO binary executable */
static char *getCdoPath()
{
  char *cdoPath = getenv("CDO_PATH");
  
  if ( cdoPath == NULL )
    {
      struct stat filestat;
      int status;
      status = stat("$HOME/bin/cdo", &filestat);
      if ( status == 0 )
	return "$HOME/bin/cdo";
      else
	{
	  fprintf(stderr, "cdo binary not found! Use CDO_PATH to set the path to cdo.\n");
	  exit(-1);
	}
    }
    
  return cdoPath;
}


/* submits a CDO command */
static int submitCdoCommand(const char *argument)
{
  static const char func[] = "submitCdoCommand";
  
  const char *cdoPath = getCdoPath();
  char *cdoCommand = (char *) malloc(strlen(cdoPath) + strlen(" ") + strlen(argument) + 1);
  
  int status;
  
  cdoCommand[0] = '\0';
  strcat(cdoCommand, cdoPath);
  strcat(cdoCommand, " ");
  strcat(cdoCommand, argument);
  
  status = system(cdoCommand);
  free(cdoCommand);
  
  return status;
}


static void testEcaFd()
{
  const double array[] = {MISSVAL, MISSVAL, TO_KELVIN(1.0), TO_KELVIN(1.0), 
    TO_KELVIN(-1.0), TO_KELVIN(-1.0)};
  
  int nvars = 1;
  int nts   = 1;
  
  double **vars = createVars(nvars, nts);
  
  writeNcFile("in.nc", array, 2);

  submitCdoCommand("eca_fd in.nc out.nc");
  readNcFile("out.nc", vars, nvars, nts);
  assert(DBL_IS_EQUAL(vars[0][0], MISSVAL));

  writeNcFile("in.nc", array, 3);

  submitCdoCommand("eca_fd in.nc out.nc");
  readNcFile("out.nc", vars, nvars, nts);
  assert(DBL_IS_EQUAL(vars[0][0], 0.0));

  writeNcFile("in.nc", array, 5);

  submitCdoCommand("eca_fd in.nc out.nc");
  readNcFile("out.nc", vars, nvars, nts);
  assert(DBL_IS_EQUAL(vars[0][0], 1.0));

  writeNcFile("in.nc", array, 6);

  submitCdoCommand("eca_fd in.nc out.nc");
  readNcFile("out.nc", vars, nvars, nts);
  assert(DBL_IS_EQUAL(vars[0][0], 2.0));
  
  destroyVars(vars);
}


static void testEcaSu()
{
  const double array[] = {MISSVAL, MISSVAL, TO_KELVIN(26.0), TO_KELVIN(24.0), 
    TO_KELVIN(26.0), TO_KELVIN(24.0)};
  
  int nvars = 1;
  int nts   = 1;
  
  double **vars = createVars(nvars, nts);
  
  writeNcFile("in.nc", array, 2);

  submitCdoCommand("eca_su in.nc out.nc");
  readNcFile("out.nc", vars, nvars, nts);
  assert(DBL_IS_EQUAL(vars[0][0], MISSVAL));

  writeNcFile("in.nc", array, 6);

  submitCdoCommand("eca_su in.nc out.nc");
  readNcFile("out.nc", vars, nvars, nts);
  assert(DBL_IS_EQUAL(vars[0][0], 2.0));

  submitCdoCommand("eca_su,20.0 in.nc out.nc");
  readNcFile("out.nc", vars, nvars, nts);
  assert(DBL_IS_EQUAL(vars[0][0], 4.0));

  submitCdoCommand("eca_su,30.0 in.nc out.nc");
  readNcFile("out.nc", vars, nvars, nts);
  assert(DBL_IS_EQUAL(vars[0][0], 0.0));
  
  destroyVars(vars);
}


static void testEcaFdns()
{
  const double array1[] = {MISSVAL, TO_KELVIN(1.0), TO_KELVIN(-1.0), 
    TO_KELVIN(-1.0), TO_KELVIN(-1.0)};
  const double array2[] = {0.0, 0.0, 1.0, 0.0, MISSVAL};
  
  int nvars = 1;
  int nts   = 1;
  
  double **vars = createVars(nvars, nts);
  
  writeNcFile("in1.nc", array1, 1);
  writeNcFile("in2.nc", array2, 1);

  submitCdoCommand("eca_fdns in1.nc in2.nc out.nc");
  readNcFile("out.nc", vars, nvars, nts);
  assert(DBL_IS_EQUAL(vars[0][0], MISSVAL));
  
  writeNcFile("in1.nc", array1, 2);
  writeNcFile("in2.nc", array2, 2);

  submitCdoCommand("eca_fdns in1.nc in2.nc out.nc");
  readNcFile("out.nc", vars, nvars, nts);
  assert(DBL_IS_EQUAL(vars[0][0], 0.0));

  writeNcFile("in1.nc", array1, 3);
  writeNcFile("in2.nc", array2, 3);

  submitCdoCommand("eca_fdns in1.nc in2.nc out.nc");
  readNcFile("out.nc", vars, nvars, nts);
  assert(DBL_IS_EQUAL(vars[0][0], 0.0));

  writeNcFile("in1.nc", array1, 4);
  writeNcFile("in2.nc", array2, 4);

  submitCdoCommand("eca_fdns in1.nc in2.nc out.nc");
  readNcFile("out.nc", vars, nvars, nts);
  assert(DBL_IS_EQUAL(vars[0][0], 1.0));

  writeNcFile("in1.nc", array1, 5);
  writeNcFile("in2.nc", array2, 5);

  submitCdoCommand("eca_fdns in1.nc in2.nc out.nc");
  readNcFile("out.nc", vars, nvars, nts);
  assert(DBL_IS_EQUAL(vars[0][0], 1.0));

  writeNcFile("in1.nc", array1 + 4, 1);
  writeNcFile("in2.nc", array2 + 4, 1);

  submitCdoCommand("eca_fdns in1.nc in2.nc out.nc");
  readNcFile("out.nc", vars, nvars, nts);
  assert(DBL_IS_EQUAL(vars[0][0], MISSVAL));

  destroyVars(vars);
}


static void testHi()
{
  const double array1[] = {MISSVAL, 70.0, 36.0, 46.0, 56.0};
  const double array2[] = {1.0, 1.0, 1.0, 1.0, 1.0};
  const double array3[] = {0.4, 0.4, 0.3, 0.5, 0.6};
  
  int nvars = 1;
  int nts   = 5;
  
  double **vars = createVars(nvars, nts);
  
  writeNcFile("in1.nc", array1, 5);
  writeNcFile("in2.nc", array2, 5);
  writeNcFile("in3.nc", array3, 5);

  submitCdoCommand("hi in1.nc in2.nc in3.nc out.nc");

  readNcFile("out.nc", vars, nvars, nts);
  assert(equals(vars[0][0], humidityIndex(array1[0], array2[0], array3[0], MISSVAL), 1.0e-5));
  assert(equals(vars[0][1], humidityIndex(array1[1], array2[1], array3[1], MISSVAL), 1.0e-5));
  assert(equals(vars[0][2], humidityIndex(array1[2], array2[2], array3[2], MISSVAL), 1.0e-5));
  assert(equals(vars[0][3], humidityIndex(array1[3], array2[3], array3[3], MISSVAL), 1.0e-5));
  assert(equals(vars[0][4], humidityIndex(array1[4], array2[4], array3[4], MISSVAL), 1.0e-5));

  destroyVars(vars);
}


int main(void)
{
  testEcaFd();
  testEcaSu();
  testEcaFdns();
  
  testHi();
  
  return 0;
}

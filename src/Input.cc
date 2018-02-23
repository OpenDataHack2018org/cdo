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

      Input     input          ASCII input
      Input     inputsrv       SERVICE input
      Input     inputext       EXTRA input
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

static size_t
input_iarray(size_t nval, int *array)
{
  size_t ival = 0;

  for (size_t i = 0; i < nval; i++)
    {
      size_t n = scanf("%d", &array[i]);
      if (n != 1) break;

      ival++;
    }

  return ival;
}

size_t input_darray(FILE *gfp, size_t nval, double *array);

void *
Input(void *process)
{
  int varID = 0;
  size_t gridsize0 = 0, gridsize = 0;
  int taxisID = 0;
  int streamID = -1;
  int vlistID = -1;
  int code = 0, level = 0, date = 0, time = 0, nlon = 0, nlat = 0;
  int output_filetype = CDI_FILETYPE_GRB;
  int ihead[8];
  double missval = 0;
  double levels[1];
  double *array = NULL;

  cdoInitialize(process);

  // clang-format off
  int INPUT    = cdoOperatorAdd("input",    0, 0, NULL);
  int INPUTSRV = cdoOperatorAdd("inputsrv", 0, 0, NULL);
  int INPUTEXT = cdoOperatorAdd("inputext", 0, 0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();

  int gridID = -1;
  int zaxisID = -1;
  if (operatorID == INPUT)
    {
      operatorInputArg("grid description file or name");
      if (operatorArgc() == 1)
        operatorCheckArgc(1);
      else
        operatorCheckArgc(2);

      gridID = cdoDefineGrid(operatorArgv()[0]);
      if (operatorArgc() == 2) zaxisID = cdoDefineZaxis(operatorArgv()[1]);
    }

  int nlevs = 1;
  if (zaxisID != -1) nlevs = zaxisInqSize(zaxisID);

  levels[0] = 0;
  int nrecs = 0;

  int tsID = 0;
  while (TRUE)
    {
      if (operatorID == INPUT)
        {
          output_filetype = cdoFiletype();

          code = -1;
          level = 0;
          gridsize = gridInqSize(gridID);
          date = 0;
          time = 0;

          if (nrecs == 0) array = (double *) Malloc(gridsize * nlevs * sizeof(double));

          cdoPrint("Enter all %zu elements of timestep %d!", gridsize * nlevs, nrecs + 1);

          size_t rval = input_darray(stdin, gridsize * nlevs, array);

          if (nrecs > 0 && rval == 0) break;

          if (rval != gridsize * nlevs) cdoAbort("Too few input elements (%zu of %zu)!", rval, gridsize * nlevs);

          if (feof(stdin)) break;
        }
      else if (operatorID == INPUTEXT)
        {
          output_filetype = cdoDefaultFileType;
          if (output_filetype == CDI_UNDEFID) output_filetype = CDI_FILETYPE_EXT;

          cdoPrint("Enter header (code,level,date,time,nlon,nlat,dispo1,dispo2)"
                   " of record %d (or EOF(=^D))!",
                   nrecs + 1);

          int rval = input_iarray(4, ihead);
          if (feof(stdin) && nrecs == 0) cdoAbort("Too few header elements (%d of %d)!", rval, 4);
          if (feof(stdin)) break;
          if (rval != 4) cdoAbort("Invalid header input!");

          date = ihead[0];
          code = ihead[1];
          level = ihead[2];
          gridsize = ihead[3];

          time = 0;

          if (nrecs == 0)
            {
              levels[0] = level;
              gridsize0 = gridsize;

              array = (double *) Malloc(gridsize * sizeof(double));

              gridID = gridCreate(GRID_GENERIC, gridsize);
            }
          else
            {
              if (gridsize != gridsize0) cdoAbort("Gridsize must not change!");
            }

          cdoPrint("Enter all %zu elements of record %d!", gridsize, nrecs + 1);

          rval = input_darray(stdin, gridsize, array);
          if (rval != (int) gridsize) cdoAbort("Invalid data input!");
        }
      else if (operatorID == INPUTSRV)
        {
          output_filetype = cdoDefaultFileType;
          if (output_filetype == CDI_UNDEFID) output_filetype = CDI_FILETYPE_SRV;

          cdoPrint("Enter header (code,level,date,time,nlon,nlat,dispo1,dispo2)"
                   " of record %d (or EOF(=^D))!",
                   nrecs + 1);

          int rval = input_iarray(8, ihead);
          if (feof(stdin) && nrecs == 0) cdoAbort("Too few header elements (%d of %d)!", rval, 8);
          if (feof(stdin)) break;
          if (rval != 8) cdoAbort("Invalid header input!");

          code = ihead[0];
          level = ihead[1];
          date = ihead[2];
          time = ihead[3];
          nlon = ihead[4];
          nlat = ihead[5];

          gridsize = nlon * nlat;

          if (nrecs == 0)
            {
              levels[0] = level;
              gridsize0 = gridsize;

              array = (double *) Malloc(gridsize * sizeof(double));

              gridID = gridCreate(GRID_GENERIC, gridsize);
              gridDefXsize(gridID, nlon);
              gridDefYsize(gridID, nlat);
            }
          else
            {
              if (gridsize != gridsize0) cdoAbort("Gridsize must not change!");
            }

          cdoPrint("Enter all %zu elements of record %d!", gridsize, nrecs + 1);

          rval = input_darray(stdin, gridsize, array);
          if (rval != (int) gridsize) cdoAbort("Invalid data input!");
        }

      if (nrecs == 0)
        {
          if (zaxisID == -1)
            {
              zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);
              zaxisDefLevels(zaxisID, levels);
            }

          vlistID = vlistCreate();
          varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARYING);
          vlistDefVarParam(vlistID, varID, cdiEncodeParam(code, 255, 255));

          missval = vlistInqVarMissval(vlistID, varID);

          taxisID = taxisCreate(TAXIS_RELATIVE);
          vlistDefTaxis(vlistID, taxisID);

          streamID = cdoStreamOpenWrite(cdoStreamName(0), output_filetype);

          pstreamDefVlist(streamID, vlistID);
        }

      int vdate = date;
      int vtime = time;
      taxisDefVdate(taxisID, vdate);
      taxisDefVtime(taxisID, vtime);
      pstreamDefTimestep(streamID, tsID);

      for (int levelID = 0; levelID < nlevs; levelID++)
        {
          int offset = gridsize * levelID;
          size_t nmiss = arrayNumMV(gridsize, array + offset, missval);
          pstreamDefRecord(streamID, varID, levelID);
          pstreamWriteRecord(streamID, array + offset, nmiss);
        }

      nrecs++;
      tsID++;
    }

  if (streamID >= 0)
    {
      pstreamClose(streamID);
      vlistDestroy(vlistID);
    }

  if (array) Free(array);

  cdoFinish();

  return 0;
}

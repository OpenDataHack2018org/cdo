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

double intlin(double x, double y1, double x1, double y2, double x2);

static void
isosurface(double isoval, long nlev1, double *lev1, Field *field3D, Field *field2D)
{
  size_t gridsize = gridInqSize(field3D->grid);
  size_t nmiss = field3D->nmiss;
  double missval = field3D->missval;
  double *data3D = field3D->ptr;
  double *data2D = field2D->ptr;

  for (size_t i = 0; i < gridsize; ++i)
    {
      data2D[i] = missval;

      for (long k = 0; k < (nlev1 - 1); ++k)
        {
          double val1 = data3D[k * gridsize + i];
          double val2 = data3D[(k + 1) * gridsize + i];

          if (nmiss > 0)
            {
              bool lmiss1 = DBL_IS_EQUAL(val1, missval);
              bool lmiss2 = DBL_IS_EQUAL(val2, missval);
              if (lmiss1 && lmiss2) continue;
              if (lmiss1 && IS_EQUAL(isoval, val2)) data2D[i] = lev1[k + 1];
              if (lmiss2 && IS_EQUAL(isoval, val1)) data2D[i] = lev1[k];
              if (lmiss1 || lmiss2) continue;
            }

          if ((isoval >= val1 && isoval <= val2) || (isoval >= val2 && isoval <= val1))
            {
              data2D[i] = IS_EQUAL(val1, val2) ? lev1[k] : intlin(isoval, lev1[k], val1, lev1[k + 1], val2);
              break;
            }
        }
    }

  field2D->missval = missval;
  field2D->nmiss = arrayNumMV(gridsize, data2D, missval);
}

void *
Isosurface(void *process)
{
  int nrecs;
  int i;
  int varID, levelID;
  int zaxisID, zaxisID1 = -1;
  size_t nmiss;

  cdoInitialize(process);

  operatorInputArg("isoval");

  operatorCheckArgc(1);

  double isoval = parameter2double(operatorArgv()[0]);

  if (cdoVerbose) cdoPrint("Isoval: %g", isoval);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int nzaxis = vlistNzaxis(vlistID1);
  int nlevel = 0;
  for (i = 0; i < nzaxis; i++)
    {
      zaxisID = vlistZaxis(vlistID1, i);
      nlevel = zaxisInqSize(zaxisID);
      if (zaxisInqType(zaxisID) != ZAXIS_HYBRID && zaxisInqType(zaxisID) != ZAXIS_HYBRID_HALF)
        if (nlevel > 1)
          {
            zaxisID1 = zaxisID;
            break;
          }
    }
  if (i == nzaxis) cdoAbort("No processable variable found!");

  int nlev1 = nlevel;
  std::vector<double> lev1(nlev1);
  cdoZaxisInqLevels(zaxisID1, &lev1[0]);

  int zaxisIDsfc = zaxisCreate(ZAXIS_SURFACE, 1);
  for (i = 0; i < nzaxis; i++)
    if (zaxisID1 == vlistZaxis(vlistID1, i)) vlistChangeZaxisIndex(vlistID2, i, zaxisIDsfc);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  size_t gridsizemax = vlistGridsizeMax(vlistID1);

  Field field;
  field_init(&field);
  field.ptr = (double *) Malloc(gridsizemax * sizeof(double));

  int nvars = vlistNvars(vlistID1);

  std::vector<bool> liso(nvars);
  std::vector<bool> vars(nvars);
  std::vector<Field> vars1(nvars);

  for (varID = 0; varID < nvars; varID++)
    {
      int gridID = vlistInqVarGrid(vlistID1, varID);
      int zaxisID = vlistInqVarZaxis(vlistID1, varID);
      size_t gridsize = gridInqSize(gridID);
      int nlevel = zaxisInqSize(zaxisID);
      double missval = vlistInqVarMissval(vlistID1, varID);

      liso[varID] = (zaxisID == zaxisID1);

      field_init(&vars1[varID]);
      vars1[varID].grid = gridID;
      vars1[varID].zaxis = zaxisID;
      vars1[varID].nmiss = 0;
      vars1[varID].missval = missval;
      vars1[varID].ptr = (double *) Malloc(gridsize * nlevel * sizeof(double));
    }

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for (varID = 0; varID < nvars; varID++)
        {
          vars[varID] = false;
          vars1[varID].nmiss = 0;
        }

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
          size_t offset = gridsize * levelID;
          double *single = vars1[varID].ptr + offset;

          pstreamReadRecord(streamID1, single, &nmiss);
          vars1[varID].nmiss += nmiss;
          vars[varID] = true;
        }

      for (varID = 0; varID < nvars; varID++)
        {
          if (vars[varID])
            {
              if (liso[varID])
                {
                  isosurface(isoval, nlev1, &lev1[0], &vars1[varID], &field);

                  pstreamDefRecord(streamID2, varID, 0);
                  pstreamWriteRecord(streamID2, field.ptr, field.nmiss);
                }
              else
                {
                  size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
                  int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
                  double missval = vlistInqVarMissval(vlistID2, varID);

                  for (levelID = 0; levelID < nlevel; levelID++)
                    {
                      size_t offset = gridsize * levelID;
                      double *single = vars1[varID].ptr + offset;
                      nmiss = arrayNumMV(gridsize, single, missval);
                      pstreamDefRecord(streamID2, varID, levelID);
                      pstreamWriteRecord(streamID2, single, nmiss);
                    }
                }
            }
        }

      tsID++;
    }

  for (varID = 0; varID < nvars; varID++) Free(vars1[varID].ptr);

  Free(field.ptr);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  vlistDestroy(vlistID2);

  cdoFinish();

  return 0;
}

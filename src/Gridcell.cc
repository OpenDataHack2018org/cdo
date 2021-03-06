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

      Gridcell   gridarea        Grid cell area in m^2
      Gridcell   gridweights     Grid cell weights
      Gridcell   gridmask        Grid mask
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "grid.h"
#include "constants.h"

static inline double
orthodrome(double px1, double py1, double px2, double py2)
{
  return acos(sin(py1) * sin(py2) + cos(py1) * cos(py2) * cos(px2 - px1));
}

void
grid_cell_area(int gridID, double *array)
{
  int gridtype = gridInqType(gridID);
  int projtype = (gridtype == GRID_PROJECTION) ? gridInqProjType(gridID) : -1;

  if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || projtype == CDI_PROJ_RLL || projtype == CDI_PROJ_LAEA
      || projtype == CDI_PROJ_SINU || projtype == CDI_PROJ_LCC || gridtype == GRID_GME || gridtype == GRID_CURVILINEAR
      || gridtype == GRID_UNSTRUCTURED)
    {
      if (gridHasArea(gridID))
        {
          if (cdoVerbose) cdoPrint("Using existing grid cell area!");
          gridInqArea(gridID, array);
        }
      else
        {
          int status = gridGenArea(gridID, array);
          if (status == 1)
            cdoAbort("%s: Grid corner missing!", __func__);
          else if (status == 2)
            cdoAbort("%s: Can't compute grid cell area for this grid!", __func__);

          size_t ngp = gridInqSize(gridID);
          for (size_t i = 0; i < ngp; ++i) array[i] *= PlanetRadius * PlanetRadius;
        }
    }
  else
    {
      if (gridtype == GRID_GAUSSIAN_REDUCED)
        cdoAbort("Unsupported grid type: %s, use CDO option -R to convert reduced to regular grid!", gridNamePtr(gridtype));
      else
        cdoAbort("%s: Unsupported grid type: %s", __func__, gridNamePtr(gridtype));
    }
}

void *
Gridcell(void *process)
{
  cdoInitialize(process);

  // clang-format off
  int GRIDAREA    = cdoOperatorAdd("gridarea",     1,  0, NULL);
  int GRIDWGTS    = cdoOperatorAdd("gridweights",  1,  0, NULL);
  int GRIDMASK    = cdoOperatorAdd("gridmask",     0,  0, NULL);
  int GRIDDX      = cdoOperatorAdd("griddx",       1,  0, NULL);
  int GRIDDY      = cdoOperatorAdd("griddy",       1,  0, NULL);
  int GRIDCELLIDX = cdoOperatorAdd("gridcellidx",  0,  0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();

  bool needRadius = cdoOperatorF1(operatorID) > 0;
  if (needRadius)
    {
      char *envstr = getenv("PLANET_RADIUS");
      if (envstr)
        {
          double fval = atof(envstr);
          if (fval > 0) PlanetRadius = fval;
        }
    }

  if (cdoVerbose) cdoPrint("PlanetRadius: %g", PlanetRadius);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);

  int ngrids = vlistNgrids(vlistID1);

  if (ngrids > 1) cdoWarning("Found more than 1 grid, using the first one!");

  int gridID = vlistGrid(vlistID1, 0);
  int zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);

  int vlistID2 = vlistCreate();
  int varID = vlistDefVar(vlistID2, gridID, zaxisID, TIME_CONSTANT);
  vlistDefNtsteps(vlistID2, 0);

  if (operatorID == GRIDAREA)
    {
      vlistDefVarName(vlistID2, varID, "cell_area");
      vlistDefVarStdname(vlistID2, varID, "area");
      vlistDefVarLongname(vlistID2, varID, "area of grid cell");
      vlistDefVarUnits(vlistID2, varID, "m2");
      vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_FLT64);
    }
  else if (operatorID == GRIDWGTS)
    {
      vlistDefVarName(vlistID2, varID, "cell_weights");
      vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_FLT64);
    }
  else if (operatorID == GRIDMASK)
    {
      vlistDefVarName(vlistID2, varID, "grid_mask");
      vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_UINT8);
    }
  else if (operatorID == GRIDDX)
    {
      vlistDefVarName(vlistID2, varID, "dx");
      vlistDefVarLongname(vlistID2, varID, "delta x");
      vlistDefVarUnits(vlistID2, varID, "m");
    }
  else if (operatorID == GRIDDY)
    {
      vlistDefVarName(vlistID2, varID, "dy");
      vlistDefVarLongname(vlistID2, varID, "delta y");
      vlistDefVarUnits(vlistID2, varID, "m");
    }
  else if (operatorID == GRIDCELLIDX)
    {
      vlistDefVarName(vlistID2, varID, "gridcellidx");
      vlistDefVarLongname(vlistID2, varID, "grid cell index");
    }

  int taxisID = taxisCreate(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID2, taxisID);

  size_t gridsize = gridInqSize(gridID);
  std::vector<double> array(gridsize);

  if (operatorID == GRIDAREA)
    {
      grid_cell_area(gridID, &array[0]);
    }
  else if (operatorID == GRIDWGTS)
    {
      int status = gridWeights(gridID, &array[0]);
      if (status != 0) cdoWarning("Using constant grid cell area weights!");
    }
  else if (operatorID == GRIDMASK)
    {
      std::vector<int> mask(gridsize);
      if (gridInqMask(gridID, NULL))
        {
          gridInqMask(gridID, &mask[0]);
        }
      else
        {
          for (size_t i = 0; i < gridsize; ++i) mask[i] = 1;
        }

      for (size_t i = 0; i < gridsize; ++i) array[i] = mask[i];
    }
  else if (operatorID == GRIDCELLIDX)
    {
      for (size_t i = 0; i < gridsize; ++i) array[i] = i + 1;
    }
  else if (operatorID == GRIDDX || operatorID == GRIDDY)
    {
      int gridtype = gridInqType(gridID);
      int projtype = (gridtype == GRID_PROJECTION) ? gridInqProjType(gridID) : -1;
      if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || projtype == CDI_PROJ_LCC || gridtype == GRID_CURVILINEAR)
        {
          double len1 = 0, len2 = 0;
          char units[CDI_MAX_NAME];

          if (gridtype != GRID_CURVILINEAR) gridID = gridToCurvilinear(gridID, 1);

          gridsize = gridInqSize(gridID);
          size_t xsize = gridInqXsize(gridID);
          size_t ysize = gridInqYsize(gridID);

          std::vector<double> xv(gridsize);
          std::vector<double> yv(gridsize);

          gridInqXvals(gridID, &xv[0]);
          gridInqYvals(gridID, &yv[0]);

          // Convert lat/lon units if required
          gridInqXunits(gridID, units);
          grid_to_radian(units, gridsize, &xv[0], "grid longitudes");
          grid_to_radian(units, gridsize, &yv[0], "grid latitudes");

          if (operatorID == GRIDDX)
            {
              for (size_t j = 0; j < ysize; ++j)
                for (size_t i = 0; i < xsize; ++i)
                  {
                    if (i == 0)
                      {
                        len2 = orthodrome(xv[j * xsize + i], yv[j * xsize + i], xv[j * xsize + i + 1], yv[j * xsize + i + 1]);
                        len1 = len2;
                      }
                    else if (i == (xsize - 1))
                      {
                        len1 = orthodrome(xv[j * xsize + i - 1], yv[j * xsize + i - 1], xv[j * xsize + i], yv[j * xsize + i]);
                        len2 = len1;
                      }
                    else
                      {
                        len1 = orthodrome(xv[j * xsize + i - 1], yv[j * xsize + i - 1], xv[j * xsize + i], yv[j * xsize + i]);
                        len2 = orthodrome(xv[j * xsize + i], yv[j * xsize + i], xv[j * xsize + i + 1], yv[j * xsize + i + 1]);
                      }

                    array[j * xsize + i] = 0.5 * (len1 + len2) * PlanetRadius;
                  }
            }
          else
            {
              for (size_t i = 0; i < xsize; ++i)
                for (size_t j = 0; j < ysize; ++j)
                  {
                    if (j == 0)
                      {
                        len2 = orthodrome(xv[j * xsize + i], yv[j * xsize + i], xv[(j + 1) * xsize + i], yv[(j + 1) * xsize + i]);
                        len1 = len2;
                      }
                    else if (j == (ysize - 1))
                      {
                        len1 = orthodrome(xv[(j - 1) * xsize + i], yv[(j - 1) * xsize + i], xv[j * xsize + i], yv[j * xsize + i]);
                        len2 = len1;
                      }
                    else
                      {
                        len1 = orthodrome(xv[(j - 1) * xsize + i], yv[(j - 1) * xsize + i], xv[j * xsize + i], yv[j * xsize + i]);
                        len2 = orthodrome(xv[j * xsize + i], yv[j * xsize + i], xv[(j + 1) * xsize + i], yv[(j + 1) * xsize + i]);
                      }

                    array[j * xsize + i] = 0.5 * (len1 + len2) * PlanetRadius;
                  }
            }
        }
      else
        {
          if (gridtype == GRID_GAUSSIAN_REDUCED)
            cdoAbort("Unsupported grid type: %s, use CDO option -R to convert reduced to regular grid!", gridNamePtr(gridtype));
          else
            cdoAbort("Unsupported grid type: %s", gridNamePtr(gridtype));
        }
    }

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);
  pstreamDefTimestep(streamID2, 0);
  pstreamDefRecord(streamID2, 0, 0);
  pstreamWriteRecord(streamID2, &array[0], 0);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

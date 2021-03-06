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

      Intgrid    interpolate     PINGO grid interpolation
      Intgrid    intgridbil      Bilinear grid interpolation
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "interpol.h"
#include "grid.h"
#include <vector>

int
genThinoutGrid(int gridID1, size_t xinc, size_t yinc)
{
  size_t nlon1 = gridInqXsize(gridID1);
  size_t nlat1 = gridInqYsize(gridID1);

  size_t nlon2 = nlon1 / xinc;
  size_t nlat2 = nlat1 / yinc;
  if (nlon1 % xinc) nlon2++;
  if (nlat1 % yinc) nlat2++;
  size_t gridsize2 = nlon2 * nlat2;

  int gridID2 = gridCreate(GRID_LONLAT, gridsize2);
  gridDefXsize(gridID2, nlon2);
  gridDefYsize(gridID2, nlat2);

  int gridtype = gridInqType(gridID1);
  if (gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT)
    {
      std::vector<double> xvals1(nlon1);
      std::vector<double> yvals1(nlat1);
      std::vector<double> xvals2(nlon2);
      std::vector<double> yvals2(nlat2);
      gridInqXvals(gridID1, &xvals1[0]);
      gridInqYvals(gridID1, &yvals1[0]);

      size_t olat = 0;
      for (size_t ilat = 0; ilat < nlat1; ilat += yinc) yvals2[olat++] = yvals1[ilat];

      size_t olon = 0;
      for (size_t ilon = 0; ilon < nlon1; ilon += xinc) xvals2[olon++] = xvals1[ilon];

      gridDefXvals(gridID2, &xvals2[0]);
      gridDefYvals(gridID2, &yvals2[0]);
    }
  else
    {
      cdoAbort("Unsupported grid: %s", gridNamePtr(gridtype));
    }

  return gridID2;
}

int
genBoxavgGrid(int gridID1, size_t xinc, size_t yinc)
{
  size_t i, j, i1;

  size_t nlon1 = gridInqXsize(gridID1);
  size_t nlat1 = gridInqYsize(gridID1);

  size_t nlon2 = nlon1 / xinc;
  size_t nlat2 = nlat1 / yinc;
  if (nlon1 % xinc) nlon2++;
  if (nlat1 % yinc) nlat2++;
  size_t gridsize2 = nlon2 * nlat2;

  int gridID2 = gridCreate(GRID_LONLAT, gridsize2);
  gridDefXsize(gridID2, nlon2);
  gridDefYsize(gridID2, nlat2);

  int gridtype = gridInqType(gridID1);
  if (gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT)
    {
      std::vector<double> xvals1(nlon1);
      std::vector<double> yvals1(nlat1);
      std::vector<double> xvals2(nlon2);
      std::vector<double> yvals2(nlat2);
      gridInqXvals(gridID1, &xvals1[0]);
      gridInqYvals(gridID1, &yvals1[0]);

      std::vector<double> grid1_corner_lon, grid1_corner_lat;
      std::vector<double> grid2_corner_lon, grid2_corner_lat;
      if (gridInqYbounds(gridID1, NULL) && gridInqXbounds(gridID1, NULL))
        {
          grid1_corner_lon.resize(2 * nlon1);
          grid1_corner_lat.resize(2 * nlat1);
          grid2_corner_lon.resize(2 * nlon2);
          grid2_corner_lat.resize(2 * nlat2);
          gridInqXbounds(gridID1, &grid1_corner_lon[0]);
          gridInqYbounds(gridID1, &grid1_corner_lat[0]);
        }

      j = 0;
      for (i = 0; i < nlon1; i += xinc)
        {
          i1 = i + (xinc - 1);
          if (i1 >= nlon1 - 1) i1 = nlon1 - 1;
          xvals2[j] = xvals1[i] + (xvals1[i1] - xvals1[i]) / 2;
          if (!grid2_corner_lon.empty())
            {
              grid2_corner_lon[2 * j] = grid1_corner_lon[2 * i];
              grid2_corner_lon[2 * j + 1] = grid1_corner_lon[2 * i1 + 1];
            }
          j++;
        }
      j = 0;
      for (i = 0; i < nlat1; i += yinc)
        {
          i1 = i + (yinc - 1);
          if (i1 >= nlat1 - 1) i1 = nlat1 - 1;
          yvals2[j] = yvals1[i] + (yvals1[i1] - yvals1[i]) / 2;
          if (!grid2_corner_lat.empty())
            {
              grid2_corner_lat[2 * j] = grid1_corner_lat[2 * i];
              grid2_corner_lat[2 * j + 1] = grid1_corner_lat[2 * i1 + 1];
            }
          j++;
        }

      gridDefXvals(gridID2, &xvals2[0]);
      gridDefYvals(gridID2, &yvals2[0]);

      if (!grid2_corner_lon.empty() && !grid2_corner_lat.empty())
        {
          gridDefNvertex(gridID2, 2);
          gridDefXbounds(gridID2, &grid2_corner_lon[0]);
          gridDefYbounds(gridID2, &grid2_corner_lat[0]);
        }
    }
  else
    {
      cdoAbort("Unsupported grid: %s", gridNamePtr(gridtype));
    }

  return gridID2;
}

static void
boxavg(Field *field1, Field *field2, size_t xinc, size_t yinc)
{
  int gridID1 = field1->grid;
  int gridID2 = field2->grid;
  double *array1 = field1->ptr;
  double *array2 = field2->ptr;
  double missval = field1->missval;

  size_t nlon1 = gridInqXsize(gridID1);
  size_t nlat1 = gridInqYsize(gridID1);

  size_t nlon2 = gridInqXsize(gridID2);
  size_t nlat2 = gridInqYsize(gridID2);

  std::vector<double *> xfield1(nlat1);
  for (size_t ilat = 0; ilat < nlat1; ilat++) xfield1[ilat] = array1 + ilat * nlon1;

  std::vector<double *> xfield2(nlat2);
  for (size_t ilat = 0; ilat < nlat2; ilat++) xfield2[ilat] = array2 + ilat * nlon2;

  for (size_t ilat = 0; ilat < nlat2; ilat++)
    for (size_t ilon = 0; ilon < nlon2; ilon++)
      {
        xfield2[ilat][ilon] = 0;

        size_t in = 0;
        for (size_t j = 0; j < yinc; ++j)
          {
            size_t jj = ilat * yinc + j;
            if (jj >= nlat1) break;
            for (size_t i = 0; i < xinc; ++i)
              {
                size_t ii = ilon * xinc + i;
                if (ii >= nlon1) break;
                in++;
                xfield2[ilat][ilon] += xfield1[jj][ii];
              }
          }
        xfield2[ilat][ilon] /= in;
      }

  field2->nmiss = arrayNumMV(nlat2 * nlon2, array2, missval);
}

static void
thinout(Field *field1, Field *field2, int xinc, int yinc)
{
  int gridID1 = field1->grid;
  int gridID2 = field2->grid;
  double *array1 = field1->ptr;
  double *array2 = field2->ptr;
  double missval = field1->missval;

  size_t nlon1 = gridInqXsize(gridID1);
  size_t nlat1 = gridInqYsize(gridID1);

  size_t nlon2 = gridInqXsize(gridID2);
  size_t nlat2 = gridInqYsize(gridID2);

  std::vector<double *> xfield1(nlat1);
  for (size_t ilat = 0; ilat < nlat1; ilat++) xfield1[ilat] = array1 + ilat * nlon1;

  std::vector<double *> xfield2(nlat2);
  for (size_t ilat = 0; ilat < nlat2; ilat++) xfield2[ilat] = array2 + ilat * nlon2;

  size_t olat = 0;
  for (size_t ilat = 0; ilat < nlat1; ilat += yinc)
    {
      size_t olon = 0;
      for (size_t ilon = 0; ilon < nlon1; ilon += xinc)
        {
          xfield2[olat][olon] = xfield1[ilat][ilon];
          olon++;
        }
      olat++;
    }

  field2->nmiss = arrayNumMV(nlat2 * nlon2, array2, missval);
}

void *
Intgrid(void *process)
{
  int nrecs;
  int varID, levelID;
  int gridID1 = -1, gridID2 = -1;
  int xinc = 0, yinc = 0;
  size_t nmiss;
  double missval;

  cdoInitialize(process);

  // clang-format off
  int INTGRIDBIL  = cdoOperatorAdd("intgridbil",  0, 0, NULL);
  int INTGRIDDIS  = cdoOperatorAdd("intgriddis",  0, 0, NULL);
  int INTGRIDNN   = cdoOperatorAdd("intgridnn",   0, 0, NULL);
  int INTERPOLATE = cdoOperatorAdd("interpolate", 0, 0, NULL);
  int BOXAVG      = cdoOperatorAdd("boxavg",      0, 0, NULL);
  int THINOUT     = cdoOperatorAdd("thinout",     0, 0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();

  if (operatorID == INTGRIDBIL || operatorID == INTERPOLATE || operatorID == INTGRIDDIS || operatorID == INTGRIDNN)
    {
      operatorInputArg("grid description file or name");
      gridID2 = cdoDefineGrid(operatorArgv()[0]);
    }
  else if (operatorID == THINOUT || operatorID == BOXAVG)
    {
      operatorInputArg("xinc, yinc");
      operatorCheckArgc(2);
      xinc = parameter2int(operatorArgv()[0]);
      yinc = parameter2int(operatorArgv()[1]);
    }

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int ngrids = vlistNgrids(vlistID1);
  for (int index = 0; index < ngrids; index++)
    {
      gridID1 = vlistGrid(vlistID1, index);
      int gridtype = gridInqType(gridID1);

      if (operatorID == BOXAVG || operatorID == THINOUT)
        {
          if (index == 0)
            {
              if (gridtype != GRID_LONLAT && gridtype != GRID_GAUSSIAN
                  /* && gridtype != GRID_CURVILINEAR */)
                cdoAbort("Interpolation of %s data unsupported!", gridNamePtr(gridtype));

              gridID2 = operatorID == BOXAVG ? genBoxavgGrid(gridID1, xinc, yinc) : genThinoutGrid(gridID1, xinc, yinc);
            }
          else
            cdoAbort("Too many different grids!");
        }
      else if (operatorID == INTGRIDBIL || operatorID == INTERPOLATE)
        {
          bool ldistgen = (grid_is_distance_generic(gridID1) && grid_is_distance_generic(gridID2));
          if (!ldistgen && gridtype != GRID_LONLAT && gridtype != GRID_GAUSSIAN)
            cdoAbort("Interpolation of %s data unsupported!", gridNamePtr(gridtype));
        }
      else if (operatorID == INTGRIDNN || operatorID == INTGRIDDIS)
        {
          int projtype = (gridtype == GRID_PROJECTION) ? gridInqProjType(gridID1) : -1;
          bool lproj4param = (gridtype == GRID_PROJECTION) && grid_has_proj4param(gridID1);
          if (gridtype != GRID_LONLAT && !lproj4param && projtype != CDI_PROJ_RLL && projtype != CDI_PROJ_LAEA
              && projtype != CDI_PROJ_SINU && projtype != CDI_PROJ_LCC && gridtype != GRID_GAUSSIAN && gridtype != GRID_GME
              && gridtype != GRID_CURVILINEAR && gridtype != GRID_UNSTRUCTURED)
            cdoAbort("Interpolation of %s data unsupported!", gridNamePtr(gridtype));
        }

      vlistChangeGridIndex(vlistID2, index, gridID2);
    }

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  size_t gridsize = vlistGridsizeMax(vlistID1);
  std::vector<double> array1(gridsize);

  gridsize = gridInqSize(gridID2);
  std::vector<double> array2(gridsize);

  Field field1, field2;
  field_init(&field1);
  field_init(&field2);

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamReadRecord(streamID1, array1.data(), &nmiss);

          gridID1 = vlistInqVarGrid(vlistID1, varID);
          missval = vlistInqVarMissval(vlistID1, varID);

          field1.grid = gridID1;
          field1.ptr = array1.data();
          field1.nmiss = nmiss;
          field1.missval = missval;
          field2.grid = gridID2;
          field2.ptr = array2.data();
          field2.missval = missval;
          field2.nmiss = 0;

          // clang-format off
	  if      ( operatorID == INTGRIDBIL )  intgridbil(&field1, &field2);
	  else if ( operatorID == INTGRIDNN )   intgriddis(&field1, &field2, 1);
	  else if ( operatorID == INTGRIDDIS )  intgriddis(&field1, &field2, 4);
	  else if ( operatorID == INTERPOLATE ) interpolate(&field1, &field2);
	  else if ( operatorID == BOXAVG )      boxavg(&field1, &field2, xinc, yinc);
	  else if ( operatorID == THINOUT )     thinout(&field1, &field2, xinc, yinc);
          // clang-format on

          pstreamDefRecord(streamID2, varID, levelID);
          pstreamWriteRecord(streamID2, array2.data(), field2.nmiss);
        }
      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

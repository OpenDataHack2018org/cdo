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

      Invert     invertlat       Invert latitude
      Invert     invertlon       Invert longitude
      Invert     invertlatdes    Invert latitude description
      Invert     invertlondes    Invert longitude description
      Invert     invertlatdata   Invert latitude data
      Invert     invertlondata   Invert longitude data
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

static void
invertLonDes(int vlistID)
{
  int ngrids = vlistNgrids(vlistID);
  for (int index = 0; index < ngrids; index++)
    {
      int gridID1 = vlistGrid(vlistID, index);
      int gridID2 = gridDuplicate(gridID1);

      int gridtype = gridInqType(gridID1);

      if (!(gridtype == GRID_GENERIC || gridtype == GRID_GAUSSIAN || gridtype == GRID_PROJECTION || gridtype == GRID_LONLAT
            || gridtype == GRID_CURVILINEAR))
        cdoAbort("Unsupported gridtype: %s!", gridNamePtr(gridtype));

      if (gridInqXvals(gridID1, NULL))
        {
          size_t nlon = gridInqXsize(gridID1);
          size_t nlat = gridInqYsize(gridID1);
          size_t size = (gridtype == GRID_CURVILINEAR) ? nlon * nlat : nlon;

          std::vector<double> xv1(size);
          std::vector<double> xv2(size);

          gridInqXvals(gridID1, xv1.data());

          if (gridtype == GRID_CURVILINEAR)
            {
              for (size_t ilat = 0; ilat < nlat; ilat++)
                for (size_t ilon = 0; ilon < nlon; ilon++) xv2[ilat * nlon + nlon - ilon - 1] = xv1[ilat * nlon + ilon];
            }
          else
            {
              for (size_t ilon = 0; ilon < nlon; ilon++) xv2[nlon - ilon - 1] = xv1[ilon];
            }

          gridDefXvals(gridID2, xv2.data());
        }

      if (gridInqXbounds(gridID1, NULL))
        {
          size_t nlon = gridInqXsize(gridID1);
          size_t nlat = gridInqYsize(gridID1);
          int nv = gridInqNvertex(gridID1);
          size_t size = (gridtype == GRID_CURVILINEAR) ? nv * nlon * nlat : nv * nlon;

          std::vector<double> xb1(size);
          std::vector<double> xb2(size);

          gridInqXbounds(gridID1, xb1.data());

          if (gridtype == GRID_CURVILINEAR)
            {
              for (size_t ilat = 0; ilat < nlat; ilat++)
                for (size_t ilon = 0; ilon < nlon; ilon++)
                  for (int iv = 0; iv < nv; iv++)
                    xb2[ilat * nlon * nv + (nlon - ilon - 1) * nv + iv] = xb1[ilat * nlon * nv + ilon * nv + iv];
            }
          else
            {
              for (size_t ilon = 0; ilon < nlon; ilon++)
                {
                  xb2[nlon * 2 - ilon * 2 - 1] = xb1[ilon * 2];
                  xb2[nlon * 2 - ilon * 2 - 2] = xb1[ilon * 2 + 1];
                }
            }

          gridDefXbounds(gridID2, xb2.data());
        }

      vlistChangeGrid(vlistID, gridID1, gridID2);
    }
}

static void
invertLatCoord(int gridID)
{
  int gridtype = gridInqType(gridID);

  if (gridInqYvals(gridID, NULL))
    {
      size_t nlon = gridInqXsize(gridID);
      size_t nlat = gridInqYsize(gridID);
      size_t size = (gridtype == GRID_CURVILINEAR) ? nlon * nlat : nlat;

      std::vector<double> yv1(size);
      std::vector<double> yv2(size);

      if (gridtype == GRID_CURVILINEAR)
        {
          gridInqXvals(gridID, yv1.data());

          for (size_t ilat = 0; ilat < nlat; ilat++)
            for (size_t ilon = 0; ilon < nlon; ilon++) yv2[(nlat - ilat - 1) * nlon + ilon] = yv1[ilat * nlon + ilon];

          gridDefXvals(gridID, yv2.data());

          gridInqYvals(gridID, yv1.data());

          for (size_t ilat = 0; ilat < nlat; ilat++)
            for (size_t ilon = 0; ilon < nlon; ilon++) yv2[(nlat - ilat - 1) * nlon + ilon] = yv1[ilat * nlon + ilon];

          gridDefYvals(gridID, yv2.data());
        }
      else
        {
          gridInqYvals(gridID, yv1.data());

          for (size_t ilat = 0; ilat < nlat; ilat++) yv2[nlat - ilat - 1] = yv1[ilat];

          gridDefYvals(gridID, yv2.data());
        }
    }

  if (gridInqYbounds(gridID, NULL))
    {
      size_t nlon = gridInqXsize(gridID);
      size_t nlat = gridInqYsize(gridID);
      int nv = gridInqNvertex(gridID);
      size_t size = (gridtype == GRID_CURVILINEAR) ? nv * nlon * nlat : nv * nlat;

      std::vector<double> yb1(size);
      std::vector<double> yb2(size);

      gridInqYbounds(gridID, yb1.data());

      if (gridtype == GRID_CURVILINEAR)
        {
          for (size_t ilat = 0; ilat < nlat; ilat++)
            for (size_t ilon = 0; ilon < nlon; ilon++)
              for (int iv = 0; iv < nv; iv++)
                yb2[(nlat - ilat - 1) * nlon * nv + ilon * nv + iv] = yb1[ilat * nlon * nv + ilon * nv + iv];
        }
      else
        {
          for (size_t ilat = 0; ilat < nlat; ilat++)
            {
              yb2[nlat * 2 - ilat * 2 - 1] = yb1[ilat * 2];
              yb2[nlat * 2 - ilat * 2 - 2] = yb1[ilat * 2 + 1];
            }
        }

      gridDefYbounds(gridID, yb2.data());
    }
}

static void
invertLatDes(int vlistID)
{
  int ngrids = vlistNgrids(vlistID);
  for (int index = 0; index < ngrids; index++)
    {
      int gridID1 = vlistGrid(vlistID, index);
      int gridID2 = gridDuplicate(gridID1);

      int gridtype = gridInqType(gridID1);

      if (!(gridtype == GRID_GENERIC || gridtype == GRID_GAUSSIAN || gridtype == GRID_PROJECTION || gridtype == GRID_LONLAT
            || gridtype == GRID_CURVILINEAR))
        cdoAbort("Unsupported gridtype: %s!", gridNamePtr(gridtype));

      invertLatCoord(gridID2);

      int projID = gridInqProj(gridID2);
      if (projID != CDI_UNDEFID) invertLatCoord(projID);

      vlistChangeGrid(vlistID, gridID1, gridID2);
    }
}

static void
invertLonData(double *array1, double *array2, int gridID1)
{
  size_t nlon = gridInqXsize(gridID1);
  size_t nlat = gridInqYsize(gridID1);

  if (nlat > 0)
    {
      double **field1 = (double **) Malloc(nlat * sizeof(double *));
      double **field2 = (double **) Malloc(nlat * sizeof(double *));

      for (size_t ilat = 0; ilat < nlat; ilat++)
        {
          field1[ilat] = array1 + ilat * nlon;
          field2[ilat] = array2 + ilat * nlon;
        }

      for (size_t ilat = 0; ilat < nlat; ilat++)
        for (size_t ilon = 0; ilon < nlon; ilon++) field2[ilat][nlon - ilon - 1] = field1[ilat][ilon];

      if (field1) Free(field1);
      if (field2) Free(field2);
    }
  else
    {
      array2[0] = array1[0];
    }
}

static void
invertLatData(double *array1, double *array2, int gridID1)
{
  size_t nlon = gridInqXsize(gridID1);
  size_t nlat = gridInqYsize(gridID1);

  if (nlat > 0)
    {
      double **field1 = (double **) Malloc(nlat * sizeof(double *));
      double **field2 = (double **) Malloc(nlat * sizeof(double *));

      for (size_t ilat = 0; ilat < nlat; ilat++)
        {
          field1[ilat] = array1 + ilat * nlon;
          field2[ilat] = array2 + ilat * nlon;
        }

      for (size_t ilat = 0; ilat < nlat; ilat++) arrayCopy(nlon, field1[ilat], field2[nlat - ilat - 1]);

      if (field1) Free(field1);
      if (field2) Free(field2);
    }
  else
    {
      array2[0] = array1[0];
    }
}

void *
Invert(void *process)
{
  int nrecs;
  int varID, levelID;
  int gridID1;
  size_t nmiss;

  cdoInitialize(process);

  // clang-format off
  cdoOperatorAdd("invertlat",     func_all, func_lat, NULL);
  cdoOperatorAdd("invertlon",     func_all, func_lon, NULL);
  cdoOperatorAdd("invertlatdes",  func_hrd, func_lat, NULL);
  cdoOperatorAdd("invertlondes",  func_hrd, func_lon, NULL);
  cdoOperatorAdd("invertlatdata", func_fld, func_lat, NULL);
  cdoOperatorAdd("invertlondata", func_fld, func_lon, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();
  int operfunc1 = cdoOperatorF1(operatorID);
  int operfunc2 = cdoOperatorF2(operatorID);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if (operfunc1 == func_all || operfunc1 == func_hrd)
    {
      if (operfunc2 == func_lat)
        invertLatDes(vlistID2);
      else
        invertLonDes(vlistID2);
    }

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());

  pstreamDefVlist(streamID2, vlistID2);

  size_t gridsize = vlistGridsizeMax(vlistID1);

  std::vector<double> array1(gridsize);
  std::vector<double> array2(gridsize);

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      pstreamDefTimestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamReadRecord(streamID1, array1.data(), &nmiss);

          pstreamDefRecord(streamID2, varID, levelID);

          if (operfunc1 == func_all || operfunc1 == func_fld)
            {
              gridID1 = vlistInqVarGrid(vlistID1, varID);

              if (operfunc2 == func_lat)
                invertLatData(array1.data(), array2.data(), gridID1);
              else
                invertLonData(array1.data(), array2.data(), gridID1);

              pstreamWriteRecord(streamID2, array2.data(), nmiss);
            }
          else
            {
              pstreamWriteRecord(streamID2, array1.data(), nmiss);
            }
        }
      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

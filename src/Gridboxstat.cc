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

      Gridboxstat    gridboxrange        Gridbox range
      Gridboxstat    gridboxmin          Gridbox minimum
      Gridboxstat    gridboxmax          Gridbox maximum
      Gridboxstat    gridboxsum          Gridbox sum
      Gridboxstat    gridboxmean         Gridbox mean
      Gridboxstat    gridboxavg          Gridbox average
      Gridboxstat    gridboxstd          Gridbox standard deviation
      Gridboxstat    gridboxstd1         Gridbox standard deviation [Normalize by (n-1)]
      Gridboxstat    gridboxvar          Gridbox variance
      Gridboxstat    gridboxvar1         Gridbox variance [Normalize by (n-1)]
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "grid.h"
#include "cdoOptions.h"

static void
genBoxGridReg2D(int gridID1, size_t xinc, size_t yinc, int gridID2)
{
  size_t i, j, i1;

  size_t nlon1 = gridInqXsize(gridID1);
  size_t nlat1 = gridInqYsize(gridID1);
  size_t nlon2 = gridInqXsize(gridID2);
  size_t nlat2 = gridInqYsize(gridID2);

  bool gridHasBounds = (gridInqXbounds(gridID1, NULL) && gridInqYbounds(gridID1, NULL));

  {
    std::vector<double> xvals1(nlon1);
    std::vector<double> yvals1(nlat1);
    std::vector<double> xvals2(nlon2);
    std::vector<double> yvals2(nlat2);
    gridInqXvals(gridID1, xvals1.data());
    gridInqYvals(gridID1, yvals1.data());

    j = 0;
    for (i = 0; i < nlon1; i += xinc)
      {
        i1 = i + (xinc - 1);
        if (i1 >= nlon1 - 1) i1 = nlon1 - 1;
        xvals2[j] = xvals1[i] + (xvals1[i1] - xvals1[i]) / 2.;
        j++;
      }

    j = 0;
    for (i = 0; i < nlat1; i += yinc)
      {
        i1 = i + (yinc - 1);
        if (i1 >= nlat1 - 1) i1 = nlat1 - 1;
        yvals2[j] = yvals1[i] + (yvals1[i1] - yvals1[i]) / 2;
        j++;
      }

    gridDefXvals(gridID2, xvals2.data());
    gridDefYvals(gridID2, yvals2.data());
  }

  if (gridHasBounds)
    {
      std::vector<double> grid1_corner_lon(2 * nlon1);
      std::vector<double> grid1_corner_lat(2 * nlat1);
      std::vector<double> grid2_corner_lon(2 * nlon2);
      std::vector<double> grid2_corner_lat(2 * nlat2);
      gridInqXbounds(gridID1, grid1_corner_lon.data());
      gridInqYbounds(gridID1, grid1_corner_lat.data());

      j = 0;
      for (i = 0; i < nlon1; i += xinc)
        {
          i1 = i + (xinc - 1);
          if (i1 >= nlon1 - 1) i1 = nlon1 - 1;
          grid2_corner_lon[2 * j] = grid1_corner_lon[2 * i];
          grid2_corner_lon[2 * j + 1] = grid1_corner_lon[2 * i1 + 1];
          j++;
        }

      j = 0;
      for (i = 0; i < nlat1; i += yinc)
        {
          i1 = i + (yinc - 1);
          if (i1 >= nlat1 - 1) i1 = nlat1 - 1;
          grid2_corner_lat[2 * j] = grid1_corner_lat[2 * i];
          grid2_corner_lat[2 * j + 1] = grid1_corner_lat[2 * i1 + 1];
          j++;
        }

      gridDefNvertex(gridID2, 2);
      gridDefXbounds(gridID2, grid2_corner_lon.data());
      gridDefYbounds(gridID2, grid2_corner_lat.data());
    }
}

static void
genBoxGridCurv2D(int gridID1, size_t xinc, size_t yinc, int gridID2)
{
  size_t nlon1 = gridInqXsize(gridID1);
  size_t nlat1 = gridInqYsize(gridID1);
  size_t nlon2 = gridInqXsize(gridID2);
  size_t nlat2 = gridInqYsize(gridID2);
  size_t gridsize1 = gridInqSize(gridID1);
  size_t gridsize2 = gridInqSize(gridID2);

  size_t x1, y1, x2, y2;
  size_t use_x1, use_y1;
  size_t g1_add, g2_add;
  int corner, add;
  int circular = gridIsCircular(gridID1);
  double on_up, on_lo, ox_up, ox_lo, an_le, an_ri, ax_le, ax_ri;
  double xvals2_0 = 0;
  double area_norm;

  bool gridHasBounds = (gridInqXbounds(gridID1, NULL) && gridInqYbounds(gridID1, NULL));

  std::vector<double> xvals1(gridsize1);
  std::vector<double> yvals1(gridsize1);
  std::vector<double> xvals2(gridsize2);
  std::vector<double> yvals2(gridsize2);
  gridInqXvals(gridID1, xvals1.data());
  gridInqYvals(gridID1, yvals1.data());

  // Convert lat/lon units if required
  {
    char units[CDI_MAX_NAME];
    gridInqXunits(gridID1, units);
    grid_to_degree(units, gridsize1, xvals1.data(), "grid center lon");
    gridInqYunits(gridID1, units);
    grid_to_degree(units, gridsize1, yvals1.data(), "grid center lat");
  }

  std::vector<double> grid1_corner_lon, grid1_corner_lat;
  std::vector<double> grid2_corner_lon, grid2_corner_lat;
  if (gridHasBounds)
    {
      grid1_corner_lon.resize(4 * gridsize1);
      grid1_corner_lat.resize(4 * gridsize1);
      grid2_corner_lon.resize(4 * gridsize2);
      grid2_corner_lat.resize(4 * gridsize2);
      gridInqXbounds(gridID1, grid1_corner_lon.data());
      gridInqYbounds(gridID1, grid1_corner_lat.data());

      // Convert lat/lon units if required
      {
        char units[CDI_MAX_NAME];
        gridInqXunits(gridID1, units);
        grid_to_degree(units, 4 * gridsize1, grid1_corner_lon.data(), "grid corner lon");
        gridInqYunits(gridID1, units);
        grid_to_degree(units, 4 * gridsize1, grid1_corner_lat.data(), "grid corner lat");
      }
    }

  // Process grid2 bounds
  area_norm = xinc * yinc;
  for (y2 = 0; y2 < nlat2; y2++)
    {
      for (x2 = 0; x2 < nlon2; x2++)
        {
          g2_add = (y2 * nlon2 + x2);
          on_up = on_lo = 360.;
          ox_up = ox_lo = -360.;
          an_ri = an_le = 90.;
          ax_ri = ax_le = -90.;

          for (y1 = y2 * yinc; y1 < yinc * (y2 + 1); y1++)
            {
              use_y1 = (y1 >= nlat1) ? nlat1 - 1 : y1;
              for (x1 = x2 * xinc; x1 < xinc * (x2 + 1); x1++)
                {
                  use_x1 = x1;
                  if (x1 >= nlon1)
                    {
                      if (circular && use_y1 == y1)
                        use_y1 -= 1;
                      else
                        use_x1 = nlon1 - 1;
                    }

                  g1_add = (use_y1 * nlon1) + use_x1;

                  if (y1 == y2 * yinc && x1 == x2 * xinc)
                    {
                      xvals2_0 = xvals1[g1_add];
                      xvals2[g2_add] = xvals1[g1_add] / area_norm;
                      yvals2[g2_add] = yvals1[g1_add] / area_norm;
                    }
                  else if (fabs(xvals1[g1_add] - xvals2_0) > 270.)
                    {
                      if ((xvals1[g1_add] - xvals2_0) > 270.)
                        xvals2[g2_add] += (xvals1[g1_add] - 360.) / area_norm;
                      else if ((xvals1[g1_add] - xvals2_0) < -270.)
                        xvals2[g2_add] += (xvals1[g1_add] + 360.) / area_norm;
                      yvals2[g2_add] += yvals1[g1_add] / area_norm;
                    }
                  else
                    {
                      xvals2[g2_add] += xvals1[g1_add] / area_norm;
                      yvals2[g2_add] += yvals1[g1_add] / area_norm;
                    }

                  if (gridHasBounds)
                    {
                      int c_flag[4], corner2, g1_add2, g1_add3;
                      double lon, lat, lon2, lat2, lon3, lat3;
                      /* upper left cell */
                      if (y1 == y2 * yinc && x1 == x2 * xinc)
                        {
                          c_flag[0] = c_flag[1] = c_flag[2] = c_flag[3] = 0;
                          for (corner = 0; corner < 4; corner++)
                            {
                              add = 4 * g1_add + corner;
                              lon = grid1_corner_lon[add];
                              lat = grid1_corner_lat[add];
                              g1_add2 = g1_add + 1;
                              if (g1_add + nlon1 > gridsize1)
                                {
                                  cdoWarning("Can't find cell below upper left");
                                  continue;
                                }
                              g1_add3 = g1_add + nlon1;
                              for (corner2 = 0; corner2 < 4; corner2++)
                                {
                                  lon2 = grid1_corner_lon[4 * g1_add2 + corner2];
                                  lat2 = grid1_corner_lat[4 * g1_add2 + corner2];
                                  lon3 = grid1_corner_lon[4 * g1_add3 + corner2];
                                  lat3 = grid1_corner_lat[4 * g1_add3 + corner2];
                                  if ((IS_EQUAL(lon2, lon) && IS_EQUAL(lat2, lat)) || (IS_EQUAL(lon3, lon) && IS_EQUAL(lat3, lat)))
                                    c_flag[corner] = 1;
                                }
                            }
                          for (corner = 0; corner < 4; corner++)
                            if (!c_flag[corner]) break;
                          on_up = grid1_corner_lon[4 * g1_add + corner];
                          ax_le = grid1_corner_lat[4 * g1_add + corner];
                          if (c_flag[0] + c_flag[1] + c_flag[2] + c_flag[3] < 3) cdoWarning("found two matching corners!");
                        }

                      /* upper right cell */
                      if ((y1 == y2 * yinc) && (x1 == (x2 + 1) * xinc - 1))
                        {
                          c_flag[0] = c_flag[1] = c_flag[2] = c_flag[3] = 0;
                          for (corner = 0; corner < 4; corner++)
                            {
                              add = 4 * g1_add + corner;
                              lon = grid1_corner_lon[add];
                              lat = grid1_corner_lat[add];
                              g1_add2 = g1_add - 1;
                              if (g1_add + nlon1 > gridsize1)
                                {
                                  cdoWarning("Can't find cell below upper left");
                                  continue;
                                }
                              g1_add3 = g1_add + nlon1;
                              for (corner2 = 0; corner2 < 4; corner2++)
                                {
                                  lon2 = grid1_corner_lon[4 * g1_add2 + corner2];
                                  lat2 = grid1_corner_lat[4 * g1_add2 + corner2];
                                  lon3 = grid1_corner_lon[4 * g1_add3 + corner2];
                                  lat3 = grid1_corner_lat[4 * g1_add3 + corner2];
                                  if ((IS_EQUAL(lon2, lon) && IS_EQUAL(lat2, lat)) || (IS_EQUAL(lon3, lon) && IS_EQUAL(lat3, lat)))
                                    c_flag[corner] = 1;
                                }
                            }
                          for (corner = 0; corner < 4; corner++)
                            if (!c_flag[corner]) break;
                          ox_up = grid1_corner_lon[4 * g1_add + corner];
                          ax_ri = grid1_corner_lat[4 * g1_add + corner];
                          if (c_flag[0] + c_flag[1] + c_flag[2] + c_flag[3] < 3) cdoWarning("found two matching corners!");
                        }

                      /* lower right cell */
                      if ((y1 == (y2 + 1) * yinc - 1) && (x1 == (x2 + 1) * xinc - 1))
                        {
                          c_flag[0] = c_flag[1] = c_flag[2] = c_flag[3] = 0;
                          for (corner = 0; corner < 4; corner++)
                            {
                              add = 4 * g1_add + corner;
                              lon = grid1_corner_lon[add];
                              lat = grid1_corner_lat[add];
                              g1_add2 = g1_add - 1;
                              if (g1_add < nlon1)
                                {
                                  cdoWarning("Can't find cell above lower right left");
                                  continue;
                                }
                              g1_add3 = g1_add - nlon1;
                              for (corner2 = 0; corner2 < 4; corner2++)
                                {
                                  lon2 = grid1_corner_lon[4 * g1_add2 + corner2];
                                  lat2 = grid1_corner_lat[4 * g1_add2 + corner2];
                                  lon3 = grid1_corner_lon[4 * g1_add3 + corner2];
                                  lat3 = grid1_corner_lat[4 * g1_add3 + corner2];
                                  if ((IS_EQUAL(lon2, lon) && IS_EQUAL(lat2, lat)) || (IS_EQUAL(lon3, lon) && IS_EQUAL(lat3, lat)))
                                    c_flag[corner] = 1;
                                }
                            }
                          for (corner = 0; corner < 4; corner++)
                            if (!c_flag[corner]) break;
                          ox_lo = grid1_corner_lon[4 * g1_add + corner];
                          an_ri = grid1_corner_lat[4 * g1_add + corner];
                          if (c_flag[0] + c_flag[1] + c_flag[2] + c_flag[3] < 3) cdoWarning("found two matching corners!");
                        }

                      /* lower left cell */
                      if ((y1 == (y2 + 1) * yinc - 1) && (x1 == x2 * xinc))
                        {
                          c_flag[0] = c_flag[1] = c_flag[2] = c_flag[3] = 0;
                          for (corner = 0; corner < 4; corner++)
                            {
                              add = 4 * g1_add + corner;
                              lon = grid1_corner_lon[add];
                              lat = grid1_corner_lat[add];
                              g1_add2 = g1_add + 1;
                              if (g1_add < nlon1)
                                {
                                  cdoWarning("Can't find cell above lower right left");
                                  continue;
                                }
                              g1_add3 = g1_add - nlon1;
                              for (corner2 = 0; corner2 < 4; corner2++)
                                {
                                  lon2 = grid1_corner_lon[4 * g1_add2 + corner2];
                                  lat2 = grid1_corner_lat[4 * g1_add2 + corner2];
                                  lon3 = grid1_corner_lon[4 * g1_add3 + corner2];
                                  lat3 = grid1_corner_lat[4 * g1_add3 + corner2];
                                  if ((IS_EQUAL(lon2, lon) && IS_EQUAL(lat2, lat)) || (IS_EQUAL(lon3, lon) && IS_EQUAL(lat3, lat)))
                                    c_flag[corner] = 1;
                                }
                            }
                          for (corner = 0; corner < 4; corner++)
                            if (!c_flag[corner]) break;
                          on_lo = grid1_corner_lon[4 * g1_add + corner];
                          an_le = grid1_corner_lat[4 * g1_add + corner];
                          if (c_flag[0] + c_flag[1] + c_flag[2] + c_flag[3] < 3) cdoWarning("found two matching corners!");
                        }
                    } /* if ( gridHasBounds) */
                }     /* for ( x1 = x2*xinc; x1 < xinc*(x2+1) ; x1++ ) */
            }         /* for ( y1 = y2*yinc; y1 < yinc*(y2+1); y1++ ) */

          if (gridHasBounds)
            {
              /* upper left corner */
              grid2_corner_lon[4 * g2_add + 3] = on_up;
              grid2_corner_lat[4 * g2_add + 3] = ax_le;
              /* upper right corner */
              grid2_corner_lon[4 * g2_add + 2] = ox_up;
              grid2_corner_lat[4 * g2_add + 2] = ax_ri;
              /* lower right corner */
              grid2_corner_lon[4 * g2_add + 1] = ox_lo;
              grid2_corner_lat[4 * g2_add + 1] = an_ri;
              /* lower left corner */
              grid2_corner_lon[4 * g2_add + 0] = on_lo;
              grid2_corner_lat[4 * g2_add + 0] = an_le;
            }

          //  while ( xvals2[g2_add] >  180. ) xvals2[g2_add] -= 360.;
          //  while ( xvals2[g2_add] < -180. ) xvals2[g2_add] += 360.;
        } /* for ( x2 = 0; x2 < nlon2; x2++ ) */
    }     /* for ( y2 = 0; y2 < nlat2; y2++ ) */

  gridDefXvals(gridID2, xvals2.data());
  gridDefYvals(gridID2, yvals2.data());

  if (gridHasBounds)
    {
      gridDefNvertex(gridID2, 4);
      gridDefXbounds(gridID2, grid2_corner_lon.data());
      gridDefYbounds(gridID2, grid2_corner_lat.data());
    }
}

static int
genBoxGrid(int gridID1, size_t xinc, size_t yinc)
{
  if (xinc < 1 || yinc < 1) cdoAbort("xinc and yinc must not be smaller than 1!");

  int gridID2 = -1;
  int gridtype = gridInqType(gridID1);
  if (gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT || gridtype == GRID_CURVILINEAR || gridtype == GRID_GENERIC)
    {
      size_t nlon1 = gridInqXsize(gridID1);
      size_t nlat1 = gridInqYsize(gridID1);
      if (xinc > nlon1 || yinc > nlat1) cdoAbort("xinc and/or yinc exceeds gridsize!");

      size_t nlon2 = nlon1 / xinc;
      size_t nlat2 = nlat1 / yinc;
      if (nlon1 % xinc) nlon2++;
      if (nlat1 % yinc) nlat2++;
      size_t gridsize2 = nlon2 * nlat2;

      gridID2 = gridCreate(gridtype, gridsize2);
      gridDefXsize(gridID2, nlon2);
      gridDefYsize(gridID2, nlat2);
    }
  else
    {
      cdoAbort("Unsupported grid: %s", gridNamePtr(gridtype));
    }

  if (gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT)
    {
      genBoxGridReg2D(gridID1, xinc, yinc, gridID2);
    }
  else if (gridtype == GRID_GENERIC)
    {
    }
  else if (gridtype == GRID_CURVILINEAR)
    {
      genBoxGridCurv2D(gridID1, xinc, yinc, gridID2);
    }
  else
    {
      cdoAbort("Unsupported grid: %s", gridNamePtr(gridtype));
    }

  return gridID2;
}

static void
gridboxstat(Field *field1, Field *field2, size_t xinc, size_t yinc, int statfunc)
{
  bool useWeight = (field1->weight != NULL);
  /*
  double findex = 0;

  progressInit();
  */

  size_t gridsize = xinc * yinc;
  Field *field = (Field *) Malloc(Threading::ompNumThreads * sizeof(Field));
  for (int i = 0; i < Threading::ompNumThreads; i++)
    {
      field[i].size = gridsize;
      field[i].ptr = (double *) Malloc(gridsize * sizeof(double));
      field[i].weight = (useWeight) ? (double *) Malloc(gridsize * sizeof(double)) : NULL;
      field[i].missval = field1->missval;
      field[i].nmiss = 0;
    }

  int gridID1 = field1->grid;
  int gridID2 = field2->grid;
  double *array1 = field1->ptr;
  double *array2 = field2->ptr;
  double missval = field1->missval;

  size_t nlon1 = gridInqXsize(gridID1);
  size_t nlat1 = gridInqYsize(gridID1);

  size_t nlon2 = gridInqXsize(gridID2);
  size_t nlat2 = gridInqYsize(gridID2);

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t ig = 0; ig < nlat2 * nlon2; ++ig)
    {
      int ompthID = cdo_omp_get_thread_num();

      /*
      int lprogress = 1;
#ifdef  _OPENMP
      if ( ompthID != 0 ) lprogress = 0;
#endif
      findex++;
      if ( lprogress ) progressStatus(0, 1, findex/nlat2*nlon2);
      */
      size_t ilat = ig / nlon2;
      size_t ilon = ig - ilat * nlon2;

      size_t isize = 0;
      field[ompthID].nmiss = 0;
      for (size_t j = 0; j < yinc; ++j)
        {
          size_t jj = ilat * yinc + j;
          if (jj >= nlat1) break;
          for (size_t i = 0; i < xinc; ++i)
            {
              size_t ii = ilon * xinc + i;
              size_t index = jj * nlon1 + ii;
              if (ii >= nlon1) break;
              field[ompthID].ptr[isize] = array1[index];
              if (useWeight) field[ompthID].weight[isize] = field1->weight[index];
              if (DBL_IS_EQUAL(field[ompthID].ptr[isize], field[ompthID].missval)) field[ompthID].nmiss++;
              isize++;
            }
        }

      field[ompthID].size = isize;
      field2->ptr[ig] = fldfun(field[ompthID], statfunc);
    }

  field2->nmiss = arrayNumMV(nlat2 * nlon2, array2, missval);

  for (int i = 0; i < Threading::ompNumThreads; i++)
    {
      if (field[i].ptr) Free(field[i].ptr);
      if (field[i].weight) Free(field[i].weight);
    }

  if (field) Free(field);
}

void *
Gridboxstat(void *process)
{
  int lastgrid = -1;
  int nrecs;
  int varID, levelID;
  bool wstatus = false;
  char varname[CDI_MAX_NAME];

  cdoInitialize(process);

  operatorInputArg("xinc, yinc");
  operatorCheckArgc(2);
  int xinc = parameter2int(operatorArgv()[0]);
  int yinc = parameter2int(operatorArgv()[1]);

  // clang-format off
  cdoOperatorAdd("gridboxrange", func_range, 0, NULL);
  cdoOperatorAdd("gridboxmin",   func_min,   0, NULL);
  cdoOperatorAdd("gridboxmax",   func_max,   0, NULL);
  cdoOperatorAdd("gridboxsum",   func_sum,   0, NULL);
  cdoOperatorAdd("gridboxmean",  func_meanw, 1, NULL);
  cdoOperatorAdd("gridboxavg",   func_avgw,  1, NULL);
  cdoOperatorAdd("gridboxvar",   func_varw,  1, NULL);
  cdoOperatorAdd("gridboxvar1",  func_var1w, 1, NULL);
  cdoOperatorAdd("gridboxstd",   func_stdw,  1, NULL);
  cdoOperatorAdd("gridboxstd1",  func_std1w, 1, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);
  bool needWeights = cdoOperatorF2(operatorID) != 0;

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int ngrids = vlistNgrids(vlistID1);
  if (ngrids > 1) cdoAbort("Too many different grids!");

  int gridID1 = vlistGrid(vlistID1, 0);
  if (gridInqType(gridID1) == GRID_GAUSSIAN_REDUCED)
    cdoAbort("Gaussian reduced grid found. Use option -R to convert it to a regular grid!");

  int gridID2 = genBoxGrid(gridID1, xinc, yinc);
  for (int index = 0; index < ngrids; index++) vlistChangeGridIndex(vlistID2, index, gridID2);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  Field field1, field2;
  field_init(&field1);
  field_init(&field2);

  size_t gridsize1 = gridInqSize(gridID1);
  field1.ptr = (double *) Malloc(gridsize1 * sizeof(double));
  field1.weight = needWeights ? (double *) Malloc(gridsize1 * sizeof(double)) : NULL;

  size_t gridsize2 = gridInqSize(gridID2);
  field2.ptr = (double *) Malloc(gridsize2 * sizeof(double));
  field2.weight = NULL;

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          size_t nmiss;
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamReadRecord(streamID1, field1.ptr, &nmiss);
          field1.nmiss = nmiss;

          field1.grid = vlistInqVarGrid(vlistID1, varID);
          field1.size = gridInqSize(field1.grid);
          field1.missval = vlistInqVarMissval(vlistID1, varID);

          field2.grid = gridID2;
          field2.size = gridsize2;
          field2.missval = field1.missval;

          if (needWeights && field1.grid != lastgrid)
            {
              lastgrid = field1.grid;
              wstatus = gridWeights(field1.grid, field1.weight);
            }
          if (wstatus != 0 && tsID == 0 && levelID == 0)
            {
              vlistInqVarName(vlistID1, varID, varname);
              cdoWarning("Using constant grid cell area weights for variable %s!", varname);
            }

          gridboxstat(&field1, &field2, xinc, yinc, operfunc);

          pstreamDefRecord(streamID2, varID, levelID);
          pstreamWriteRecord(streamID2, field2.ptr, field2.nmiss);
        }
      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if (field1.ptr) Free(field1.ptr);
  if (field1.weight) Free(field1.weight);

  if (field2.ptr) Free(field2.ptr);
  if (field2.weight) Free(field2.weight);

  cdoFinish();

  return 0;
}

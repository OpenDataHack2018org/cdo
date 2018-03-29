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
#include "grid.h"
#include "util.h" /* progressStatus */
#include "cdoOptions.h"

/**
* Find the interval i-1 .. i in which an element x fits and return i, the
* bigger one of the interval borders or x itself if it is an interval border.
*
* If no interval can be found return the length of the array.

* @param *array ascending or descending sorted list
* @param nelem  length of the sorted list
* @param x      the element to find a position for
*/
static long
find_element(double x, long nelem, const double *restrict array)
{
  long ii;
  long mid = 0;
  long first = 1;
  long last = nelem;

  if (array[0] < array[nelem - 1])  // ascending order
    {
      // return the length of the array if x is out of bounds
      if (x < array[0] || x > array[nelem - 1]) return nelem;

      // search for the interval in which x fits
      // implementation: binary search algorithm
      for (ii = 1; ii < nelem; ++ii)
        {
          // binary search: divide search room in the middle
          // mid = first + ((last - first) >> 1);
          // faster!
          mid = (first + last) >> 1;

          // return the bigger interval border of the interval in which x fits
          // if ( x >= array[mid-1] && x <= array[mid] ) break;
          // faster!
          if (!(x < array[mid - 1] || x > array[mid])) break;

          // binary search: ignore half of the search room
          if (x > array[mid])
            first = mid;
          else
            last = mid;
        }
    }
  else
    {
      // return the length of the array if x is out of bounds
      if (x < array[nelem - 1] || x > array[0]) return nelem;

      // search for the interval in which x fits
      // implementation: binary search algorithm
      for (ii = 1; ii < nelem; ++ii)
        {
          // binary search: divide search room in the middle
          // mid = first + ((last - first) >> 1);
          // faster!
          mid = (first + last) >> 1;

          // return the bigger interval border of the interval in which x fits
          // if ( x >= array[mid] && x <= array[mid-1] ) break;
          // faster!
          if (!(x < array[mid] || x > array[mid - 1])) break;

          // binary search: ignore half of the search room
          if (x < array[mid])
            first = mid;
          else
            last = mid;
        }
    }

  if (mid > 1 && IS_EQUAL(x, array[mid - 1])) mid--;

  return mid;
}
/*
static
long find_element(double x, long nelem, const double *array)
{
  long ii;

  if ( array[0] < array[nelem-1] )
    {
      for ( ii = 1; ii < nelem; ii++ )
        if ( x >= array[ii-1] && x <= array[ii] ) break;
    }
  else
    {
      for ( ii = 1; ii < nelem; ii++ )
        if ( x >= array[ii] && x <= array[ii-1] ) break;
    }

  return ii;
}
*/

int
rect_grid_search(size_t *ii, size_t *jj, double x, double y, size_t nxm, size_t nym, const double *restrict xm,
                 const double *restrict ym)
{
  constexpr double rtol = 1.e-12;
  int lfound = 0;

  *jj = find_element(y, nym, ym);
  if (*jj >= nym && fabs(ym[0] - y) < rtol) *jj = 1;  // fix rounding errors

  if (*jj < nym)
    {
      *ii = find_element(x, nxm, xm);
      if (*ii >= nxm && fabs(xm[0] - x) < rtol) *ii = 1;  // fix rounding errors

      if (*ii < nxm) lfound = 1;
    }

  return lfound;
}

int
rect_grid_search2(long *imin, long *imax, double xmin, double xmax, long nxm, const double *restrict xm)
{
  int lfound = 0;
  *imin = nxm;
  *imax = -1;

  bool lascend = (xm[0] < xm[nxm - 1]);

  long i1 = find_element(xmin, nxm, xm);
  long i2 = find_element(xmax, nxm, xm);

  if (i1 > 0 && i1 < nxm)
    {
      lfound = 1;

      if (lascend)
        {
          if (i1 > 1 && xmin <= xm[i1 - 1]) i1--;
          *imin = i1 - 1;
          *imax = i1 - 1;
        }
      else
        {
          if (i1 < nxm - 1 && xmin <= xm[i1]) i1++;
          *imin = i1 - 1;
          *imax = i1 - 1;
        }
    }

  if (i2 > 0 && i2 < nxm)
    {
      lfound = 1;

      if (lascend)
        {
          if (i2 < nxm - 1 && xmax >= xm[i2]) i2++;
          *imax = i2 - 1;
          if (*imin == nxm) *imin = *imax;
        }
      else
        {
          if (i2 > 1 && xmax >= xm[i2 - 1]) i2--;
          *imin = i2 - 1;
          if (*imax == -1) *imax = *imin;
        }
    }

  return lfound;
}

double
intlinarr2p(long nxm, long nym, double **fieldm, const double *xm, const double *ym, double x, double y)
{
  long ii, jj;
  double value = 0;

  for (jj = 1; jj < nym; jj++)
    if (y >= MIN(ym[jj - 1], ym[jj]) && y <= MAX(ym[jj - 1], ym[jj])) break;

  for (ii = 1; ii < nxm; ii++)
    if (x >= xm[ii - 1] && x <= xm[ii]) break;

  if (jj < nym && ii < nxm)
    {
      value = fieldm[jj - 1][ii - 1] * (x - xm[ii]) * (y - ym[jj]) / ((xm[ii - 1] - xm[ii]) * (ym[jj - 1] - ym[jj]))
              + fieldm[jj - 1][ii] * (x - xm[ii - 1]) * (y - ym[jj]) / ((xm[ii] - xm[ii - 1]) * (ym[jj - 1] - ym[jj]))
              + fieldm[jj][ii - 1] * (x - xm[ii]) * (y - ym[jj - 1]) / ((xm[ii - 1] - xm[ii]) * (ym[jj] - ym[jj - 1]))
              + fieldm[jj][ii] * (x - xm[ii - 1]) * (y - ym[jj - 1]) / ((xm[ii] - xm[ii - 1]) * (ym[jj] - ym[jj - 1]));
    }

  return value;
}

static void
intlinarr2(double missval, int lon_is_circular, size_t nxm, size_t nym, const double *restrict fieldm,
           const double *restrict xm, const double *restrict ym, size_t gridsize2, double *field, const double *restrict x,
           const double *restrict y)
{
  size_t nlon1 = nxm;
  double findex = 0;

  if (lon_is_circular) nlon1--;
  size_t gridsize1 = nlon1 * nym;

  bool *grid1_mask = (bool *) Malloc(gridsize1 * sizeof(bool));
  for (size_t jj = 0; jj < nym; ++jj)
    for (size_t ii = 0; ii < nlon1; ++ii) grid1_mask[jj * nlon1 + ii] = !DBL_IS_EQUAL(fieldm[jj * nlon1 + ii], missval);

  progressInit();

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(findex, \
    field, fieldm, x, y, xm, ym, nxm, nym, gridsize2, missval, nlon1, lon_is_circular, grid1_mask)
#endif
  for (size_t i = 0; i < gridsize2; ++i)
    {
      size_t src_add[4]; /*  address for the four source points    */
      int lprogress = 1;
      if (cdo_omp_get_thread_num() != 0) lprogress = 0;

      field[i] = missval;

#ifdef _OPENMP
#pragma omp atomic
#endif
      findex++;
      if (lprogress) progressStatus(0, 1, findex / gridsize2);

      size_t ii, jj;
      int lfound = rect_grid_search(&ii, &jj, x[i], y[i], nxm, nym, xm, ym);

      if (lfound)
        {
          size_t iix = ii;
          if (lon_is_circular && iix == (nxm - 1)) iix = 0;
          src_add[0] = (jj - 1) * nlon1 + (ii - 1);
          src_add[1] = (jj - 1) * nlon1 + (iix);
          src_add[2] = (jj) *nlon1 + (ii - 1);
          src_add[3] = (jj) *nlon1 + (iix);

          /* Check to see if points are missing values */
          for (int n = 0; n < 4; ++n)
            if (!grid1_mask[src_add[n]]) lfound = 0;
        }

      if (lfound)
        {
          double wgts[4];
          wgts[0] = (x[i] - xm[ii]) * (y[i] - ym[jj]) / ((xm[ii - 1] - xm[ii]) * (ym[jj - 1] - ym[jj]));
          wgts[1] = (x[i] - xm[ii - 1]) * (y[i] - ym[jj]) / ((xm[ii] - xm[ii - 1]) * (ym[jj - 1] - ym[jj]));
          wgts[3] = (x[i] - xm[ii - 1]) * (y[i] - ym[jj - 1]) / ((xm[ii] - xm[ii - 1]) * (ym[jj] - ym[jj - 1]));
          wgts[2] = (x[i] - xm[ii]) * (y[i] - ym[jj - 1]) / ((xm[ii - 1] - xm[ii]) * (ym[jj] - ym[jj - 1]));

          // printf("%2ld %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",
          // dst_add, plon, plat, wgts[0], wgts[1], wgts[2], wgts[3], iw, jw);

          field[i] = 0;
          for (int n = 0; n < 4; ++n) field[i] += fieldm[src_add[n]] * wgts[n];
        }
    }

  progressStatus(0, 1, 1);

  if (grid1_mask) Free(grid1_mask);
}

double
intlin(double x, double y1, double x1, double y2, double x2)
{
  /*
    xlin - lineare interpolation

    Uwe Schulzweida  04/05/1995
  */
  double value = (y2 * (x - x1) + y1 * (x2 - x)) / (x2 - x1);

  return value;
}

void
intlinarr(long nxm, double *ym, double *xm, int nx, double *y, double *x)
{
  /*
    xlinarr - lineare interpolation over 1D array

    Uwe Schulzweida  04/05/1995
  */
  for (long jj = 1; jj < nxm; jj++)
    for (long j = 0; j < nx; j++)
      if (x[j] >= xm[jj - 1] && x[j] <= xm[jj]) y[j] = intlin(x[j], ym[jj - 1], xm[jj - 1], ym[jj], xm[jj]);
}

void
intgridbil(field_type *field1, field_type *field2)
{
  char xunits[CDI_MAX_NAME];
  int gridID1 = field1->grid;
  int gridID2 = field2->grid;
  double *array1 = field1->ptr;
  double *array2 = field2->ptr;
  double missval = field1->missval;

  size_t nlon1 = gridInqXsize(gridID1);
  size_t nlat1 = gridInqYsize(gridID1);

  int lon_is_circular = 0;

  bool lgeorefgrid = true;
  if (grid_is_distance_generic(gridID1) && grid_is_distance_generic(gridID2)) lgeorefgrid = false;

  double **array1_2D = (double **) Malloc(nlat1 * sizeof(double *));
  for (size_t ilat = 0; ilat < nlat1; ilat++) array1_2D[ilat] = array1 + ilat * nlon1;

  if (lgeorefgrid)
    {
      if (!(gridInqXvals(gridID1, NULL) && gridInqYvals(gridID1, NULL))) cdoAbort("Source grid has no coordinate values!");

      lon_is_circular = gridIsCircular(gridID1);

      if (lon_is_circular) nlon1 += 1;
    }

  double *lon1 = (double *) Malloc(nlon1 * sizeof(double));
  double *lat1 = (double *) Malloc(nlat1 * sizeof(double));
  gridInqXvals(gridID1, lon1);
  gridInqYvals(gridID1, lat1);

  if (lgeorefgrid)
    {
      if (lon_is_circular) lon1[nlon1 - 1] = 0;

      gridInqXunits(gridID1, xunits);

      grid_to_radian(xunits, nlon1, lon1, "grid1 center lon");
      grid_to_radian(xunits, nlat1, lat1, "grid1 center lat");

      if (lon_is_circular) lon1[nlon1 - 1] = lon1[0] + 2 * M_PI;
    }

  size_t xsize2 = gridInqXsize(gridID2);
  size_t ysize2 = gridInqYsize(gridID2);

  if (xsize2 == 1 && ysize2 == 1)
    {
      double lon2, lat2;
      gridInqXvals(gridID2, &lon2);
      gridInqYvals(gridID2, &lat2);

      gridInqXunits(gridID2, xunits);

      grid_to_radian(xunits, xsize2, &lon2, "grid2 center lon");
      grid_to_radian(xunits, ysize2, &lat2, "grid2 center lat");

      if (lon2 < lon1[0]) lon2 += 2 * M_PI;

      if (lon2 > lon1[nlon1 - 1])
        {
          double **field = array1_2D;
          array1_2D = (double **) Malloc(nlat1 * sizeof(double *));
          lon1 = (double *) Realloc(lon1, (nlon1 + 1) * sizeof(double));
          double *array = (double *) Malloc(nlat1 * (nlon1 + 1) * sizeof(double));

          for (size_t ilat = 0; ilat < nlat1; ilat++)
            {
              array1_2D[ilat] = array + ilat * (nlon1 + 1);
              memcpy(array1_2D[ilat], field[ilat], nlon1 * sizeof(double));
              array1_2D[ilat][nlon1] = array1_2D[ilat][0];
              lon1[nlon1] = lon1[0] + 2 * M_PI;
            }
          nlon1++;
          Free(field);
          Free(array);
        }

      if (lon2 < lon1[0] || lon2 > lon1[nlon1 - 1])
        cdoAbort("Longitude %f out of bounds (%f to %f)!", lon2, lon1[0], lon1[nlon1 - 1]);

      if (lat2 < MIN(lat1[0], lat1[nlat1 - 1]) || lat2 > MAX(lat1[0], lat1[nlat1 - 1]))
        cdoAbort("Latitude %f out of bounds (%f to %f)!", lat2, lat1[0], lat1[nlat1 - 1]);

      *array2 = intlinarr2p(nlon1, nlat1, array1_2D, lon1, lat1, lon2, lat2);
      /*
      printf("%5d %f %f %f\n", index++, lon2, lat2, *array2);
      */
    }
  else
    {
      if (lgeorefgrid)
        {
          if (gridInqType(gridID2) == GRID_GME) gridID2 = gridToUnstructured(gridID2, 0);

          if (gridInqType(gridID2) != GRID_UNSTRUCTURED && gridInqType(gridID2) != GRID_CURVILINEAR)
            gridID2 = gridToCurvilinear(gridID2, 0);

          if (!(gridInqXvals(gridID2, NULL) && gridInqYvals(gridID2, NULL))) cdoAbort("Target grid has no coordinate values!");
        }

      size_t gridsize2 = gridInqSize(gridID2);

      double *xvals2 = (double *) Malloc(gridsize2 * sizeof(double));
      double *yvals2 = (double *) Malloc(gridsize2 * sizeof(double));

      if (lgeorefgrid)
        {
          gridInqXvals(gridID2, xvals2);
          gridInqYvals(gridID2, yvals2);

          gridInqXunits(gridID2, xunits);

          grid_to_radian(xunits, gridsize2, xvals2, "grid2 center lon");
          grid_to_radian(xunits, gridsize2, yvals2, "grid2 center lat");

          for (size_t i = 0; i < gridsize2; ++i)
            {
              if (xvals2[i] < lon1[0]) xvals2[i] += 2 * M_PI;
              if (xvals2[i] > lon1[nlon1 - 1]) xvals2[i] -= 2 * M_PI;
            }
        }
      else
        {
          double *xcoord = (double *) Malloc(xsize2 * sizeof(double));
          double *ycoord = (double *) Malloc(ysize2 * sizeof(double));

          gridInqXvals(gridID2, xcoord);
          gridInqYvals(gridID2, ycoord);

          for (size_t j = 0; j < ysize2; ++j)
            for (size_t i = 0; i < xsize2; ++i)
              {
                xvals2[j * xsize2 + i] = xcoord[i];
                yvals2[j * xsize2 + i] = ycoord[j];
              }

          Free(xcoord);
          Free(ycoord);
        }

      intlinarr2(missval, lon_is_circular, nlon1, nlat1, array1, lon1, lat1, gridsize2, array2, xvals2, yvals2);

      field2->nmiss = arrayNumMV(gridsize2, array2, missval);

      Free(xvals2);
      Free(yvals2);
    }

  Free(lon1);
  Free(lat1);
  Free(array1_2D);
}

/* source code from pingo */
void
interpolate(field_type *field1, field_type *field2)
{
  int i;
  double *lono_array, *lato_array, *lono, *lato;
  double *lon_array, *lat_array, *lon, *lat;
  // size_t gridsize_i
  size_t gridsize_o;
  int gridIDi;
  double *arrayIn;
  int gridIDo;
  double *arrayOut;
  double missval;
  long ilon, ilat, nxlon, nxlat, olon, olat;
  long l11, l12, l21, l22, l1, l2;
  double volon1, volon2, volat1, volat2;
  double *volon11, *volon12;
  double *volon21, *volon22;
  double vilon1, vilon2, vilat1, vilat2;
  double vlon1, vlon2, vlat1, vlat2;
  long *ilon11, *ilon12;
  long *ilon21, *ilon22;
  long *ilat1, *ilat2;
  long ilon1, ilon2;
  double sum, wsum;
  int k, n;
  double *xin_array, *xlon, *xlat;
  double **in0, **xin, **xout;
  int wrap_around, xlat_is_ascending;
  double a11, a12, a21, a22, b11, b12, b21, b22, t;
  double faclon1, faclon2, faclat1, faclat2;
  int nlon, nlat, out_nlon, out_nlat;

  gridIDi = field1->grid;
  gridIDo = field2->grid;
  arrayIn = field1->ptr;
  arrayOut = field2->ptr;
  missval = field1->missval;

  // gridsize_i = gridInqSize(gridIDi);
  gridsize_o = gridInqSize(gridIDo);

  nlon = gridInqXsize(gridIDi);
  nlat = gridInqYsize(gridIDi);
  out_nlon = gridInqXsize(gridIDo);
  out_nlat = gridInqYsize(gridIDo);

  lon_array = (double *) Malloc((nlon + 2) * sizeof(double));
  lat_array = (double *) Malloc((nlat + 2) * sizeof(double));
  lon = lon_array + 1;
  lat = lat_array + 1;

  if (!(gridInqXvals(gridIDi, NULL) && gridInqYvals(gridIDi, NULL))) cdoAbort("Source grid has no values");

  if (!(gridInqXvals(gridIDo, NULL) && gridInqYvals(gridIDo, NULL))) cdoAbort("Target grid has no values");

  gridInqXvals(gridIDi, lon);
  gridInqYvals(gridIDi, lat);

  /* Convert lat/lon units if required */
  {
    char units[CDI_MAX_NAME];
    gridInqXunits(gridIDi, units);
    grid_to_degree(units, nlon, lon, "grid1 center lon");
    gridInqYunits(gridIDi, units);
    grid_to_degree(units, nlat, lat, "grid1 center lat");
  }

  if (nlon > 1)
    {
      lon[-1] = lon[nlon - 1] - 360 > 2 * lon[0] - lon[1] ? lon[nlon - 1] - 360 : 2 * lon[0] - lon[1];
      lon[nlon] = lon[0] + 360 < 2 * lon[nlon - 1] - lon[nlon - 2] ? lon[0] + 360 : 2 * lon[nlon - 1] - lon[nlon - 2];
    }
  else
    {
      lon[-1] = lon[0] - 360;
      lon[1] = lon[0] + 360;
    }

  if (nlat > 1)
    {
      lat[-1] = 2 * lat[0] - lat[1];
      lat[nlat] = 2 * lat[nlat - 1] - lat[nlat - 2];
    }
  else
    {
      lat[-1] = lat[0] - 10;
      lat[1] = lat[nlat - 1] + 10;
    }

  if (lat[-1] < -90) lat[-1] = -99;
  if (lat[-1] > 90) lat[-1] = 99;
  if (lat[nlat] < -90) lat[nlat] = -99;
  if (lat[nlat] > 90) lat[nlat] = 99;

  lono_array = (double *) Malloc((out_nlon < 2 ? 4 : out_nlon + 2) * sizeof(double));
  lono = lono_array + 1;
  lato_array = (double *) Malloc((out_nlat < 2 ? 4 : out_nlat + 2) * sizeof(double));
  lato = lato_array + 1;

  gridInqXvals(gridIDo, lono);
  gridInqYvals(gridIDo, lato);

  /* Convert lat/lon units if required */
  {
    char units[CDI_MAX_NAME];
    gridInqXunits(gridIDo, units);
    grid_to_degree(units, out_nlon, lono, "grid2 center lon");
    gridInqYunits(gridIDo, units);
    grid_to_degree(units, out_nlat, lato, "grid2 center lat");
  }

  for (i = 0; i < out_nlon - 1; i++)
    if (lono[i + 1] <= lono[i]) break;

  for (i++; i < out_nlon; i++)
    {
      lono[i] += 360;
      if (i < out_nlon - 1 && lono[i + 1] + 360 <= lono[i]) cdoAbort("Longitudes of output grid are not in ascending order!");
    }

  if (lono[out_nlon - 1] - lono[0] >= 360) cdoAbort("The area covered by the longitudes of output grid must not overlap!");

  if (lato[0] > 90.001 || lato[out_nlat - 1] > 90.001 || lato[0] < -90.001 || lato[out_nlat - 1] < -90.001)
    {
      cdoAbort("Latitudes of output grid must be between 90 and -90!");
    }

  for (i = 0; i < out_nlat - 1; i++)
    if (IS_EQUAL(lato[i + 1], lato[i]) || (i < out_nlat - 2 && ((lato[i + 1] > lato[i]) != (lato[i + 2] > lato[i + 1]))))
      {
        cdoAbort("Latitudes of output grid must be in descending or ascending "
                 "order!");
      }

  if (out_nlon > 1)
    {
      lono[-1] = lono[out_nlon - 1] - 360 > 2 * lono[0] - lono[1] ? lono[out_nlon - 1] - 360 : 2 * lono[0] - lono[1];
      lono[out_nlon] = lono[0] + 360 < 2 * lono[out_nlon - 1] - lono[out_nlon - 2]
                           ? lono[0] + 360
                           : 2 * lono[out_nlon - 1] - lono[out_nlon - 2];
    }
  else
    {
      lono[-1] = lono[0] - 360;
      lono[1] = lono[0] + 360;
    }

  if (out_nlat > 1)
    {
      lato[-1] = 2 * lato[0] - lato[1];
      lato[out_nlat] = 2 * lato[out_nlat - 1] - lato[out_nlat - 2];
    }
  else
    {
      lato[-1] = lato[0] - 10;
      lato[1] = lato[out_nlat - 1] + 10;
    }

  if (lato[-1] < -90) lato[-1] = -99;
  if (lato[-1] > 90) lato[-1] = 99;
  if (lato[out_nlat] < -90) lato[out_nlat] = -99;
  if (lato[out_nlat] > 90) lato[out_nlat] = 99;

  nxlon = 2 * nlon + 1;
  nxlat = 2 * nlat + 1;
  xin_array = (double *) Malloc(nxlon * nxlat * sizeof(double));
  xin = (double **) Malloc(nxlat * sizeof(double *));

  for (ilat = 0; ilat < nxlat; ilat++) xin[ilat] = xin_array + ilat * nxlon;

  xlon = (double *) Malloc(nxlon * sizeof(double));
  for (ilon = 0; ilon < nlon; ilon++)
    {
      xlon[2 * ilon + 1] = lon[ilon];
      xlon[2 * ilon] = (lon[ilon - 1] + lon[ilon]) / 2;
    }
  xlon[2 * nlon] = (lon[nlon - 1] + lon[nlon]) / 2;

  xlat = (double *) Malloc((2 * nlat + 1) * sizeof(double));
  for (ilat = 0; ilat < nlat; ilat++)
    {
      xlat[2 * ilat + 1] = lat[ilat];
      xlat[2 * ilat] = (lat[ilat - 1] + lat[ilat]) / 2;
    }
  xlat[2 * nlat] = (lat[nlat - 1] + lat[nlat]) / 2;

  in0 = (double **) Malloc(nlat * sizeof(double *));
  for (ilat = 0; ilat < nlat; ilat++) in0[ilat] = arrayIn + ilat * nlon;

  ilon11 = (long *) Malloc(out_nlon * sizeof(long));
  ilon12 = (long *) Malloc(out_nlon * sizeof(long));
  ilon21 = (long *) Malloc(out_nlon * sizeof(long));
  ilon22 = (long *) Malloc(out_nlon * sizeof(long));
  volon11 = (double *) Malloc(out_nlon * sizeof(double));
  volon12 = (double *) Malloc(out_nlon * sizeof(double));
  volon21 = (double *) Malloc(out_nlon * sizeof(double));
  volon22 = (double *) Malloc(out_nlon * sizeof(double));

  for (olon = 0; olon < out_nlon; olon++)
    {
      volon1 = (lono[olon - 1] + lono[olon]) / 2;
      volon2 = (lono[olon] + lono[olon + 1]) / 2;
      if (IS_EQUAL(volon1, volon2)) volon2 += 360;
      volon2 -= 360 * floor((volon1 - xlon[0]) / 360);
      volon1 -= 360 * floor((volon1 - xlon[0]) / 360);
      volon21[olon] = volon1;
      volon22[olon] = volon2;
      for (l21 = 0; l21 < nxlon && xlon[l21] < volon1; l21++)
        ;
      for (l22 = l21; l22 < nxlon && xlon[l22] < volon2; l22++)
        ;
      volon1 -= 360;
      volon2 -= 360;
      volon11[olon] = volon1;
      volon12[olon] = volon2;
      for (l11 = 0; xlon[l11] < volon1; l11++)
        ;
      for (l12 = l11; l12 < nxlon && xlon[l12] < volon2; l12++)
        ;
      ilon11[olon] = l11;
      ilon12[olon] = l12;
      ilon21[olon] = l21;
      ilon22[olon] = l22;
    }

  ilat1 = (long *) Malloc(out_nlat * sizeof(long));
  ilat2 = (long *) Malloc(out_nlat * sizeof(long));

  xlat_is_ascending = xlat[0] <= xlat[nxlat - 1];
  for (olat = 0; olat < out_nlat; olat++)
    {
      volat1 = (lato[olat - 1] + lato[olat]) / 2;
      volat2 = (lato[olat] + lato[olat + 1]) / 2;
      if (!xlat_is_ascending)
        {
          if (volat1 > volat2)
            {
              for (l1 = 0; l1 < nxlat && xlat[l1] > volat1; l1++)
                ;
              for (l2 = l1; l2 < nxlat && xlat[l2] > volat2; l2++)
                ;
            }
          else
            {
              for (l1 = 0; l1 < nxlat && xlat[l1] > volat2; l1++)
                ;
              for (l2 = l1; l2 < nxlat && xlat[l2] > volat1; l2++)
                ;
            }
        }
      else
        {
          if (volat1 < volat2)
            {
              for (l1 = 0; l1 < nxlat && xlat[l1] < volat1; l1++)
                ;
              for (l2 = l1; l2 < nxlat && xlat[l2] < volat2; l2++)
                ;
            }
          else
            {
              for (l1 = 0; l1 < nxlat && xlat[l1] < volat2; l1++)
                ;
              for (l2 = l1; l2 < nxlat && xlat[l2] < volat1; l2++)
                ;
            }
        }

      ilat1[olat] = l1;
      ilat2[olat] = l2;
    }

  xout = (double **) Malloc(out_nlat * sizeof(double *));
  for (olat = 0; olat < out_nlat; olat++) xout[olat] = arrayOut + olat * out_nlon;

  wrap_around = nlon > 1 && (lon[nlon - 1] >= lon[-1] + 360 - 0.001 || lon[nlon] >= lon[0] + 360 - 0.001);

  for (ilat = 0; ilat < nlat; ilat++)
    for (ilon = 0; ilon < nlon; ilon++) xin[2 * ilat + 1][2 * ilon + 1] = in0[ilat][ilon];

  for (ilat = 0; ilat < nxlat; ilat += 2)
    for (ilon = 1; ilon < nxlon; ilon += 2)
      {
        sum = 0;
        n = 0;
        if (ilat > 0 && !DBL_IS_EQUAL(xin[ilat - 1][ilon], missval))
          {
            sum += xin[ilat - 1][ilon];
            n++;
          }
        if (ilat < nxlat - 1 && !DBL_IS_EQUAL(xin[ilat + 1][ilon], missval))
          {
            sum += xin[ilat + 1][ilon];
            n++;
          }
        xin[ilat][ilon] = n ? sum / n : missval;
      }

  for (ilat = 1; ilat < nxlat; ilat += 2)
    for (ilon = 0; ilon < nxlon; ilon += 2)
      {
        sum = 0;
        n = 0;
        if (ilon > 0 && !DBL_IS_EQUAL(xin[ilat][ilon - 1], missval))
          {
            sum += xin[ilat][ilon - 1];
            n++;
          }
        if (ilon == 0 && wrap_around && !DBL_IS_EQUAL(xin[ilat][2 * nlon - 1], missval))
          {
            sum += xin[ilat][2 * nlon - 1];
            n++;
          }
        if (ilon < nxlon - 1 && !DBL_IS_EQUAL(xin[ilat][ilon + 1], missval))
          {
            sum += xin[ilat][ilon + 1];
            n++;
          }
        if (ilon == nxlon - 1 && wrap_around && !DBL_IS_EQUAL(xin[ilat][1], missval))
          {
            sum += xin[ilat][1];
            n++;
          }
        xin[ilat][ilon] = n ? sum / n : missval;
      }

  for (ilat = 0; ilat < nxlat; ilat += 2)
    for (ilon = 0; ilon < nxlon; ilon += 2)
      {
        sum = 0;
        n = 0;
        if (ilon > 0 && !DBL_IS_EQUAL(xin[ilat][ilon - 1], missval))
          {
            sum += xin[ilat][ilon - 1];
            n++;
          }
        if (ilon == 0 && wrap_around && !DBL_IS_EQUAL(xin[ilat][2 * nlon - 1], missval))
          {
            sum += xin[ilat][2 * nlon - 1];
            n++;
          }
        if (ilon < nxlon - 1 && !DBL_IS_EQUAL(xin[ilat][ilon + 1], missval))
          {
            sum += xin[ilat][ilon + 1];
            n++;
          }
        if (ilon == nxlon - 1 && wrap_around && !DBL_IS_EQUAL(xin[ilat][1], missval))
          {
            sum += xin[ilat][1];
            n++;
          }
        if (ilat > 0 && !DBL_IS_EQUAL(xin[ilat - 1][ilon], missval))
          {
            sum += xin[ilat - 1][ilon];
            n++;
          }
        if (ilat < nxlat - 1 && !DBL_IS_EQUAL(xin[ilat + 1][ilon], missval))
          {
            sum += xin[ilat + 1][ilon];
            n++;
          }
        xin[ilat][ilon] = n ? sum / n : missval;
      }

  for (olat = 0; olat < out_nlat; olat++)
    {
      if (lato[-1] < lato[out_nlat])
        {
          volat1 = (lato[olat - 1] + lato[olat]) / 2;
          volat2 = (lato[olat] + lato[olat + 1]) / 2;
        }
      else
        {
          volat2 = (lato[olat - 1] + lato[olat]) / 2;
          volat1 = (lato[olat] + lato[olat + 1]) / 2;
        }

      for (olon = 0; olon < out_nlon; olon++)
        {
          sum = 0;
          wsum = 0;
          for (k = 0; k < 2; k++)
            {
              if (k == 0)
                {
                  ilon1 = ilon11[olon];
                  ilon2 = ilon12[olon];
                  volon1 = volon11[olon];
                  volon2 = volon12[olon];
                }
              else
                {
                  ilon1 = ilon21[olon];
                  ilon2 = ilon22[olon];
                  volon1 = volon21[olon];
                  volon2 = volon22[olon];
                }

              for (ilon = ilon1; ilon <= ilon2; ilon++)
                {
                  if (ilon == 0 || ilon == nxlon) continue;
                  vilon1 = xlon[ilon - 1];
                  vilon2 = xlon[ilon];
                  for (ilat = ilat1[olat]; ilat <= ilat2[olat]; ilat++)
                    {
                      if (ilat == 0 || ilat == nxlat) continue;
                      if (xlat_is_ascending)
                        {
                          vilat1 = xlat[ilat - 1];
                          vilat2 = xlat[ilat];
                          a11 = xin[ilat - 1][ilon - 1];
                          a12 = xin[ilat - 1][ilon];
                          a21 = xin[ilat][ilon - 1];
                          a22 = xin[ilat][ilon];
                        }
                      else
                        {
                          vilat1 = xlat[ilat];
                          vilat2 = xlat[ilat - 1];
                          a11 = xin[ilat][ilon - 1];
                          a12 = xin[ilat][ilon];
                          a21 = xin[ilat - 1][ilon - 1];
                          a22 = xin[ilat - 1][ilon];
                        }
                      if (DBL_IS_EQUAL(a11, missval) || DBL_IS_EQUAL(a12, missval) || DBL_IS_EQUAL(a21, missval)
                          || DBL_IS_EQUAL(a22, missval))
                        {
                          continue;
                        }
                      if (volon1 <= vilon1 && vilon2 <= volon2 && volat1 <= vilat1 && vilat2 <= volat2)
                        {
                          vlon1 = vilon1 * M_PI / 180;
                          vlon2 = vilon2 * M_PI / 180;
                          vlat1 = vilat1 * M_PI / 180;
                          vlat2 = vilat2 * M_PI / 180;
                          b11 = a11;
                          b12 = a12;
                          b21 = a21;
                          b22 = a22;
                        }
                      else
                        {
                          vlon1 = (volon1 <= vilon1 ? vilon1 : volon1);
                          vlon2 = (vilon2 <= volon2 ? vilon2 : volon2);
                          vlat1 = (volat1 <= vilat1 ? vilat1 : volat1);
                          vlat2 = (vilat2 <= volat2 ? vilat2 : volat2);
                          if (vlon1 >= vlon2 - (volon2 - volon1) * 1e-5 || vlat1 >= vlat2 - (volat2 - volat1) * 1e-5)
                            {
                              continue;
                            }
                          faclon1 = (vlon1 - vilon1) / (vilon2 - vilon1);
                          faclon2 = (vlon2 - vilon1) / (vilon2 - vilon1);
                          faclat1 = (vlat1 - vilat1) / (vilat2 - vilat1);
                          faclat2 = (vlat2 - vilat1) / (vilat2 - vilat1);
                          vlon1 *= M_PI / 180;
                          vlon2 *= M_PI / 180;
                          vlat1 *= M_PI / 180;
                          vlat2 *= M_PI / 180;
                          b11 = a11 + (a12 - a11) * faclon1 + (a21 - a11) * faclat1
                                + (a22 - a12 - a21 + a11) * faclon1 * faclat1;
                          b12 = a11 + (a12 - a11) * faclon2 + (a21 - a11) * faclat1
                                + (a22 - a12 - a21 + a11) * faclon2 * faclat1;
                          b21 = a11 + (a12 - a11) * faclon1 + (a21 - a11) * faclat2
                                + (a22 - a12 - a21 + a11) * faclon1 * faclat2;
                          b22 = a11 + (a12 - a11) * faclon2 + (a21 - a11) * faclat2
                                + (a22 - a12 - a21 + a11) * faclon2 * faclat2;
                        }
                      wsum += (vlon2 - vlon1) * (sin(vlat2) - sin(vlat1));
                      t = 2 * sin((vlat2 + vlat1) / 2) * sin((vlat2 - vlat1) / 2) / (vlat2 - vlat1);
                      sum += (vlon2 - vlon1) / 2 * ((b11 + b12) * (t - sin(vlat1)) + (b21 + b22) * (sin(vlat2) - t));
                    }
                }
            }
          xout[olat][olon] = IS_NOT_EQUAL(wsum, 0) ? sum / wsum : missval;
        }
    }

  field2->nmiss = arrayNumMV(gridsize_o, arrayOut, missval);

  Free(lon_array);
  Free(lat_array);
  Free(lono_array);
  Free(lato_array);
  Free(xin);
  Free(xin_array);
  Free(xlon);
  Free(xlat);
  Free(in0);
  Free(ilon11);
  Free(ilon12);
  Free(ilon21);
  Free(ilon22);
  Free(volon11);
  Free(volon12);
  Free(volon21);
  Free(volon22);
  Free(ilat1);
  Free(ilat2);
  Free(xout);
}

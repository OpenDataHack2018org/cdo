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
#include "cdo_int.h"
#include "remap.h"

int
gridSearchSquareReg2dNN(size_t nx, size_t ny, size_t *restrict nbr_add, double *restrict nbr_dist, double plat, double plon,
                        const double *restrict src_center_lat, const double *restrict src_center_lon)
{
  int search_result = 0;

  const double coslat_dst = cos(plat);
  const double sinlat_dst = sin(plat);
  const double coslon_dst = cos(plon);
  const double sinlon_dst = sin(plon);

  double dist_min = DBL_MAX;
  for (unsigned n = 0; n < 4; ++n) nbr_dist[n] = DBL_MAX;

  size_t jjf = 0, jjl = ny - 1;
  if (plon >= src_center_lon[0] && plon <= src_center_lon[nx - 1])
    {
      if (src_center_lat[0] < src_center_lat[ny - 1])
        {
          if (plat <= src_center_lat[0])
            {
              jjf = 0;
              jjl = (ny == 1) ? 0 : 1;
            }
          else
            {
              jjf = (ny == 1) ? 0 : ny - 2;
              jjl = ny - 1;
            }
        }
      else
        {
          if (plat >= src_center_lat[0])
            {
              jjf = 0;
              jjl = (ny == 1) ? 0 : 1;
            }
          else
            {
              jjf = (ny == 1) ? 0 : ny - 2;
              jjl = ny - 1;
            }
        }
    }

  std::vector<double> sincoslon(nx);

  for (size_t ii = 0; ii < nx; ++ii) sincoslon[ii] = coslon_dst * cos(src_center_lon[ii]) + sinlon_dst * sin(src_center_lon[ii]);

  for (size_t jj = jjf; jj <= jjl; ++jj)
    {
      const double coslat = coslat_dst * cos(src_center_lat[jj]);
      const double sinlat = sinlat_dst * sin(src_center_lat[jj]);

      size_t jjskip = jj > 1 && jj < (ny - 2);

      for (size_t ii = 0; ii < nx; ++ii)
        {
          if (jjskip && ii > 1 && ii < (nx - 2)) continue;

          double distance = acos(coslat * sincoslon[ii] + sinlat);
          if (distance < dist_min)
            {
              size_t srch_add = jj * nx + ii;
              for (unsigned n = 0; n < 4; ++n)
                {
                  if (distance < nbr_dist[n])
                    {
                      for (unsigned i = 3; i > n; --i)
                        {
                          nbr_add[i] = nbr_add[i - 1];
                          nbr_dist[i] = nbr_dist[i - 1];
                        }
                      search_result = -1;
                      nbr_add[n] = srch_add;
                      nbr_dist[n] = distance;
                      dist_min = nbr_dist[3];
                      break;
                    }
                }
            }
        }
    }

  for (unsigned n = 0; n < 4; ++n) nbr_dist[n] = 1.0 / (nbr_dist[n] + TINY);
  double distance = 0.0;
  for (unsigned n = 0; n < 4; ++n) distance += nbr_dist[n];
  for (unsigned n = 0; n < 4; ++n) nbr_dist[n] /= distance;

  return search_result;
}

int
gridSearchSquareReg2d(RemapGrid *src_grid, size_t *restrict src_add, double *restrict src_lats, double *restrict src_lons,
                      double plat, double plon)
{
  /*
    Input variables:

      double plat                    ! latitude  of the search point
      double plon                    ! longitude of the search point

    Output variables:

      int    src_add[4]              ! address of each corner point enclosing P
      double src_lats[4]             ! latitudes  of the four corner points
      double src_lons[4]             ! longitudes of the four corner points
  */
  int search_result = 0;
  const double *restrict src_center_lat = src_grid->reg2d_center_lat;
  const double *restrict src_center_lon = src_grid->reg2d_center_lon;

  for (unsigned n = 0; n < 4; ++n) src_add[n] = 0;

  const size_t nx = src_grid->dims[0];
  const size_t ny = src_grid->dims[1];

  const size_t nxm = src_grid->is_cyclic ? nx + 1 : nx;

  if (plon < src_center_lon[0]) plon += PI2;
  if (plon > src_center_lon[nxm - 1]) plon -= PI2;

  size_t ii, jj;
  int lfound = rect_grid_search(&ii, &jj, plon, plat, nxm, ny, src_center_lon, src_center_lat);
  if (lfound)
    {
      size_t iix = ii;
      if (src_grid->is_cyclic && iix == (nxm - 1)) iix = 0;
      src_add[0] = (jj - 1) * nx + (ii - 1);
      src_add[1] = (jj - 1) * nx + (iix);
      src_add[2] = (jj) *nx + (iix);
      src_add[3] = (jj) *nx + (ii - 1);

      src_lons[0] = src_center_lon[ii - 1];
      src_lons[1] = src_center_lon[iix];
      /* For consistency, we must make sure all lons are in same 2pi interval */
      if (src_lons[0] > PI2) src_lons[0] -= PI2;
      if (src_lons[0] < 0) src_lons[0] += PI2;
      if (src_lons[1] > PI2) src_lons[1] -= PI2;
      if (src_lons[1] < 0) src_lons[1] += PI2;
      src_lons[2] = src_lons[1];
      src_lons[3] = src_lons[0];

      src_lats[0] = src_center_lat[jj - 1];
      src_lats[1] = src_lats[0];
      src_lats[2] = src_center_lat[jj];
      src_lats[3] = src_lats[2];

      search_result = 1;

      return search_result;
    }

  /*
    If no cell found, point is likely either in a box that straddles either pole
    or is outside the grid. Fall back to a distance-weighted average of the four
    closest points. Go ahead and compute weights here, but store in src_lats and
    return -add to prevent the parent routine from computing bilinear weights.
  */
  if (!src_grid->lextrapolate) return search_result;

  search_result = gridSearchSquareReg2dNN(nx, ny, src_add, src_lats, plat, plon, src_center_lat, src_center_lon);

  return search_result;
}

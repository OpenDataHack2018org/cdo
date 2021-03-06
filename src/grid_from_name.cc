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
#include "griddes.h"
#include "util_string.h"

size_t genIcosphereCoords(int subdivisions, bool lbounds, std::vector<double> &xvals, std::vector<double> &yvals,
                          std::vector<double> &xbounds, std::vector<double> &ybounds);

static void
gen_grid_icosphere(GridDesciption &grid, const char *pline)
{
  int gridtype = GRID_UNSTRUCTURED;
  bool lbounds = true;
  long b = 0;

  if (*pline != 0)
    {
      if (*pline == 'r')
        pline++;
      else
        return;

      if (*pline == 0) return;
      if (!isdigit((int) *pline)) return;

      char *endptr = (char *) pline;
      long r = strtol(pline, &endptr, 10);
      if (*endptr == 0 || r != 2) return;
      pline = endptr;

      if (*pline == 'b')
        pline++;
      else
        return;

      if (*pline == 0) return;
      if (!isdigit((int) *pline)) return;

      endptr = (char *) pline;
      b = strtol(pline, &endptr, 10);

      if (*endptr != 0)
        {
          pline = endptr;
          if (*pline != '_') return;
          pline++;
          if (*pline == 0) return;
          if (*pline == '0')
            {
              lbounds = false;
              pline++;
            }
          if (*pline != 0) return;
        }
    }

  grid.type = gridtype;
  if (lbounds) grid.nvertex = 3;

  size_t ncells = genIcosphereCoords(b + 1, lbounds, grid.xvals, grid.yvals, grid.xbounds, grid.ybounds);
  grid.xsize = ncells;
  grid.ysize = ncells;
  strcpy(grid.xname, "clon");
  strcpy(grid.yname, "clat");
  strcpy(grid.xunits, "radian");
  strcpy(grid.yunits, "radian");
}

static void
gen_grid_lonlat(GridDesciption &grid, const char *pline, double inc, double lon1, double lon2, double lat1, double lat2)
{
  int gridtype = GRID_LONLAT;
  bool lbounds = true;

  if (*pline != 0 && (*pline == '+' || *pline == '-') && (isdigit((int) *(pline + 1)) || ispunct((int) *(pline + 1))))
    {
      char *endptr = (char *) pline;
      double off = strtod(pline, &endptr);
      pline = endptr;

      lon1 -= off;
      lon2 += off;
      lat1 -= off;
      lat2 += off;
      if (lat1 < -90) lat1 = -90;
      if (lat2 > 90) lat2 = 90;
    }

  if (*pline != 0)
    {
      if (*pline == '_')
        pline++;
      else
        return;

      if (*pline == 0) return;

      if (!isdigit((int) *pline) && !ispunct((int) *pline)) return;

      char *endptr = (char *) pline;
      inc = strtod(pline, &endptr);
      if (*endptr != 0)
        {
          pline = endptr;
          if (*pline == '_')
            pline++;
          else
            return;

          if (*pline == 0) return;
          if (*pline == 'c')
            {
              gridtype = GRID_CURVILINEAR;
              pline++;
              if (*pline == '0')
                {
                  lbounds = false;
                  pline++;
                }
            }
          else if (*pline == 'u')
            {
              gridtype = GRID_UNSTRUCTURED;
              pline++;
              if (*pline == '0')
                {
                  lbounds = false;
                  pline++;
                }
            }
          if (*pline != 0) return;
        }

      if (inc < 1e-9) inc = 1;
    }

  grid.type = gridtype;

  if (lon1 >= lon2 || lat1 >= lat2) cdoAbort("Invalid grid box: lon1=%g lon2=%g lat1=%g lat2=%g", lon1, lon2, lat1, lat2);

  size_t nlon = (size_t)((lon2 - lon1) / inc + 0.5);
  size_t nlat = (size_t)((lat2 - lat1) / inc + 0.5);

  grid.xvals.resize(nlon);
  grid.yvals.resize(nlat);

  for (size_t i = 0; i < nlon; ++i) grid.xvals[i] = lon1 + inc / 2 + i * inc;
  for (size_t i = 0; i < nlat; ++i) grid.yvals[i] = lat1 + inc / 2 + i * inc;

  if (gridtype == GRID_LONLAT)
    {
      grid.xsize = nlon;
      grid.ysize = nlat;
    }
  else
    {
      std::vector<double> yvals(nlat);
      for (size_t j = 0; j < nlat; ++j) yvals[j] = grid.yvals[j];
      size_t gridsize = nlon * nlat;
      grid.xvals.resize(gridsize);
      grid.yvals.resize(gridsize);
      for (size_t j = 0; j < nlat; ++j)
        for (size_t i = 0; i < nlon; ++i)
          {
            grid.xvals[j * nlon + i] = grid.xvals[i];
            grid.yvals[j * nlon + i] = yvals[j];
          }

      if (gridtype == GRID_CURVILINEAR)
        {
          grid.xsize = nlon;
          grid.ysize = nlat;
        }
      else
        {
          grid.xsize = gridsize;
          grid.ysize = gridsize;
          if (lbounds) grid.nvertex = 4;
        }

      if (lbounds && nlon > 1 && nlat > 1)
        {
          std::vector<double> xbounds(2 * nlon);
          std::vector<double> ybounds(2 * nlat);

          grid_gen_bounds(nlon, grid.xvals.data(), xbounds.data());
          grid_gen_bounds(nlat, yvals.data(), ybounds.data());
          grid_check_lat_borders(2 * nlat, ybounds.data());

          grid.xbounds.resize(4 * gridsize);
          grid.ybounds.resize(4 * gridsize);
          grid_gen_xbounds2D(nlon, nlat, xbounds.data(), grid.xbounds.data());
          grid_gen_ybounds2D(nlon, nlat, ybounds.data(), grid.ybounds.data());
        }
    }
}

int
grid_from_name(const char *gridnameptr)
{
  const char *pline;
  int gridID = CDI_UNDEFID;
  GridDesciption grid;
  size_t len;
  char *endptr;

  char *gridname = strdup(gridnameptr);
  strtolower(gridname);

  if (gridname[0] == 't' && gridname[1] == 'l') /* tl<RES>grid or tl<RES>spec */
    {
      pline = &gridname[2];
      if (isdigit((int) *pline))
        {
          grid.ntr = atoi(pline);
          while (isdigit((int) *pline)) pline++;
          if (cmpstrlen(pline, "grid", len) == 0)
            grid.type = GRID_GAUSSIAN;
          else if (cmpstrlen(pline, "zon", len) == 0)
            grid.type = GRID_GAUSSIAN;
          else if (cmpstrlen(pline, "spec", len) == 0)
            grid.type = GRID_SPECTRAL;
          else if (cmpstrlen(pline, "", len) == 0)
            grid.type = GRID_SPECTRAL;

          if (pline[len] != 0) return gridID;

          if (grid.type == GRID_GAUSSIAN)
            {
              grid.ysize = ntr_to_nlat_linear(grid.ntr);
              grid.np = grid.ysize / 2;
              if (cmpstrlen(pline, "zon", len) == 0)
                grid.xsize = 1;
              else
                grid.xsize = nlat_to_nlon(grid.ysize);

              grid.def_xfirst = true;
              grid.def_yfirst = true;
            }
        }
    }
  else if (gridname[0] == 't') /* t<RES>grid or t<RES>spec */
    {
      pline = &gridname[1];
      if (isdigit((int) *pline))
        {
          grid.ntr = atoi(pline);
          while (isdigit((int) *pline)) pline++;
          if (cmpstrlen(pline, "grid", len) == 0)
            grid.type = GRID_GAUSSIAN;
          else if (cmpstrlen(pline, "zon", len) == 0)
            grid.type = GRID_GAUSSIAN;
          else if (cmpstrlen(pline, "spec", len) == 0)
            grid.type = GRID_SPECTRAL;
          else if (cmpstrlen(pline, "", len) == 0)
            grid.type = GRID_SPECTRAL;

          if (pline[len] != 0) return gridID;

          if (grid.type == GRID_GAUSSIAN)
            {
              grid.ysize = ntr_to_nlat(grid.ntr);
              grid.np = grid.ysize / 2;
              if (cmpstrlen(pline, "zon", len) == 0)
                grid.xsize = 1;
              else
                grid.xsize = nlat_to_nlon(grid.ysize);

              grid.def_xfirst = true;
              grid.def_yfirst = true;
            }
        }
    }
  else if (gridname[0] == 'r') /* r<LON>x<LAT> */
    {
      pline = &gridname[1];
      if (isdigit((int) *pline))
        {
          grid.type = GRID_LONLAT;
          grid.xsize = atoi(pline);
          while (isdigit((int) *pline)) pline++;
          pline++;
          grid.ysize = atoi(pline);
          while (isdigit((int) *pline)) pline++;

          grid.def_xfirst = true;
          grid.def_yfirst = true;
        }
    }
  else if (gridname[0] == 'l' && gridname[1] == 'o' && gridname[2] == 'n') /* lon=<LON>_lat=<LAT> */
    {
      /* only one gridpoint */
      pline = &gridname[3];
      if (*pline == '=') pline++;
      if (isdigit((int) *pline) || ispunct((int) *pline) || *pline == '-')
        {
          grid.type = GRID_LONLAT;
          grid.xsize = 1;
          grid.ysize = 1;
          grid.xvals.resize(1);
          grid.yvals.resize(1);
          grid.xvals[0] = atof(pline);
          while (isdigit((int) *pline) || ispunct((int) *pline) || *pline == '-') pline++;
          if (*pline == '_') pline++;
          if (!(pline[0] == 'l' && pline[1] == 'a' && pline[2] == 't')) return gridID;
          pline += 3;
          if (*pline == '=') pline++;
          if (isdigit((int) *pline) || ispunct((int) *pline) || *pline == '-')
            grid.yvals[0] = atof(pline);
          else
            return gridID;
        }
    }
  else if (gridname[0] == 'g' && gridname[1] == 'm' && gridname[2] == 'e') /* gme<NI> */
    {
      pline = &gridname[3];
      if (isdigit((int) *pline))
        {
          long ni = strtol(pline, &endptr, 10);
          if (*endptr == 0)
            {
              grid.type = GRID_GME;
              grid.ni = ni;
              grid.nd = 10;
              factorni(grid.ni, &grid.ni2, &grid.ni3);
              grid.size = (grid.ni + 1) * (grid.ni + 1) * 10;
            }
        }
    }
  else if (gridname[0] == 'n' && gridname[1] == 'i') /* ni<NI> */
    {
      pline = &gridname[2];
      if (isdigit((int) *pline))
        {
          long ni = strtol(pline, &endptr, 10);
          if (*endptr == 0)
            {
              grid.type = GRID_GME;
              grid.ni = ni;
              grid.nd = 10;
              factorni(grid.ni, &grid.ni2, &grid.ni3);
              grid.size = (grid.ni + 1) * (grid.ni + 1) * 10;
            }
        }
    }
  else if (gridname[0] == 'n') /* n<N> */
    {
      pline = &gridname[1];
      if (isdigit((int) *pline))
        {
          long np = strtol(pline, &endptr, 10);
          pline = endptr;

          if (cmpstrlen(pline, "zon", len) == 0)
            {
              grid.xsize = 1;
              pline += 3;
            }
          else if (*pline == 'b')
            {
              grid.genBounds = true;
              pline++;
            }

          if (*pline == 0)
            {
              grid.type = GRID_GAUSSIAN;
              grid.np = np;
              grid.ysize = np * 2;
              if (!grid.xsize) grid.xsize = nlat_to_nlon(grid.ysize);

              grid.def_xfirst = true;
              grid.def_yfirst = true;
            }
        }
    }
  else if (gridname[0] == 'g' && isdigit(gridname[1])) /* g<LON>x<LAT> or g<SIZE> */
    {
      pline = &gridname[1];
      if (isdigit((int) *pline))
        {
          grid.type = GRID_GENERIC;
          grid.xsize = atoi(pline);
          while (isdigit((int) *pline)) pline++;
          if (*pline)
            {
              pline++;
              grid.ysize = atoi(pline);
              while (isdigit((int) *pline)) pline++;
            }
          else if (grid.xsize == 1)
            {
              grid.size = 1;
              grid.xsize = 0;
            }
        }
    }
  else if (cmpstrlen(gridname, "germany", len) == 0)  // germany_Xdeg
    {
      double lon1 = 5.6, lon2 = 15.2;
      double lat1 = 47.1, lat2 = 55.1;
      double dll = 0.1;

      pline = &gridname[len];
      gen_grid_lonlat(grid, pline, dll, lon1, lon2, lat1, lat2);
    }
  else if (cmpstrlen(gridname, "europe", len) == 0)  // europe_Xdeg
    {
      double lon1 = -30, lon2 = 60;
      double lat1 = 30, lat2 = 80;
      double dll = 1;

      pline = &gridname[len];
      gen_grid_lonlat(grid, pline, dll, lon1, lon2, lat1, lat2);
    }
  else if (cmpstrlen(gridname, "africa", len) == 0)  // africa_Xdeg
    {
      double lon1 = -20, lon2 = 60;
      double lat1 = -40, lat2 = 40;
      double dll = 1;

      pline = &gridname[len];
      gen_grid_lonlat(grid, pline, dll, lon1, lon2, lat1, lat2);
    }
  else if (cmpstrlen(gridname, "global", len) == 0)  // global_Xdeg
    {
      double lon1 = -180, lon2 = 180;
      double lat1 = -90, lat2 = 90;
      double dll = 1;

      pline = &gridname[len];
      gen_grid_lonlat(grid, pline, dll, lon1, lon2, lat1, lat2);
    }
  else if (cmpstrlen(gridname, "ico", len) == 0)  // icoR02BXX
    {
      pline = &gridname[len];
      gen_grid_icosphere(grid, pline);
    }

  if (grid.type != -1) gridID = gridDefine(grid);

  free(gridname);

  return gridID;
}

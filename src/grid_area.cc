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

static void
lonlat_to_xyz(double lon, double lat, double *xyz)
{
  double coslat = cos(lat);
  xyz[0] = coslat * cos(lon);
  xyz[1] = coslat * sin(lon);
  xyz[2] = sin(lat);
}

/*
extern "C" {
#include "lib/yac/grid.h"
#include "lib/yac/grid_cell.h"
#include "lib/yac/area.h"
}

static
double yac_huiliers_area(int num_corners, double *cell_corner_lon, double
*cell_corner_lat)
{
  if ( num_corners < 3 ) return 0;

  double coordinates_xyz[num_corners*3];
  enum edge_type edge_types[num_corners];
  struct grid_cell cell =
    {.coordinates_x   = cell_corner_lon,
     .coordinates_y   = cell_corner_lat,
     .coordinates_xyz = coordinates_xyz,
     .edge_type       = edge_types,
     .num_corners     = num_corners};

  for ( int i = 0; i < num_corners; ++i ) edge_types[i] = GREAT_CIRCLE;
  for ( int i = 0; i < num_corners; ++i )
    lonlat_to_xyz(cell_corner_lon[i], cell_corner_lat[i], coordinates_xyz+i*3);

  double area = huiliers_area(cell);
  area /= (EarthRadius*EarthRadius);

  return area;
}
*/

static void
cross_product(const double *restrict a, const double *restrict b, double *restrict c)
{
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

/*
static
double scalar_product(const double *restrict a, const double *restrict b)
{
  return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}
*/

static double
norm(const double *restrict a)
{
  return (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}
/*
static
double mod_cell_area(int num_corners, double *cell_corner_lon, double
*cell_corner_lat)
{
  if ( num_corners < 3 ) return 0;

  // generalised version based on the ICON code, mo_base_geometry.f90
  // provided by Luis Kornblueh, MPI-M.

  int M = num_corners; // number of vertices
  int m; // loop index over number of vertices M
  int i; // loop index over the three dimensions

  double area = 0.0;
  double s[M];
  double ca[M];
  double a[M];

  double p[M][3];
  double u[M][3];

  // Convert into cartesian coordinates
  for ( m = 0; m < M; m++ )
    lonlat_to_xyz(cell_corner_lon[m], cell_corner_lat[m], p[m]);

  // First, compute cross products Uij = Vi x Vj.
  for ( m = 0; m < M; m++ )
    cross_product(p[m], p[(m+1)%M], u[m]);

  //  Normalize Uij to unit vectors.
  for ( m = 0; m < M; m++ )
    {
      s[m] = norm(u[m]);
      area += s[m];
    }

  // Test for a degenerated cells associated with collinear vertices.

  if ( fabs(area) > 0.0 )
    {
      for ( m = 0; m < M; m++ )
        s[m] = sqrt(s[m]);

      for ( m = 0; m < M; m++ )
        for ( i = 0; i < 3; i++ )
          u[m][i] = u[m][i]/s[m];

      //  Compute interior angles Ai as the dihedral angles between planes
      //  by using the definition of the scalar product
      //
      //	    ab = |a| |b| cos (phi)
      //
      //  As a and b are already normalised this reduces to
      //
      //            ab = cos (phi)

      //  There is no explanation so far for the - in the loop below.
      //  But otherwise we don't get the correct results for triangles
      //  and cells. Must have something to do with the theorem.

      for ( m = 0; m < M; m++ )
        {
          ca[m] = - scalar_product(u[m], u[(m+1)%M]);
          if ( ca[m] < -1.0 ) ca[m] = -1.0;
          if ( ca[m] >  1.0 ) ca[m] =  1.0;
          a[m] = acos(ca[m]);
        }

      //  Compute areas = a1 + a2 + a3 - (M-2) * pi. here for a unit sphere:
      area = - (double) (M-2) * M_PI;

      for ( m = 0; m < M; m++ )
        area += a[m];

      // area *= EarthRadius * EarthRadius;
      if ( area < 0.0 ) area = 0.0;
    }

  return area;
}
*/

/** area of a spherical triangle based on L'Huilier's Theorem
 *
 * source code is taken from code by Robert Oehmke of Earth System Modeling
 * Framework (www.earthsystemmodeling.org)
 *
 * the license statement for this routine is as follows:
 * Earth System Modeling Framework
 * Copyright 2002-2013, University Corporation for Atmospheric Research,
 * Massachusetts Institute of Technology, Geophysical Fluid Dynamics
 * Laboratory, University of Michigan, National Centers for Environmental
 * Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
 * NASA Goddard Space Flight Center.
 * Licensed under the University of Illinois-NCSA License.
 */
static double
mod_tri_area(const double *restrict u, const double *restrict v, const double *restrict w)
{
  double tmp_vec[3];

  cross_product(u, v, tmp_vec);
  double sina = sqrt(norm(tmp_vec));
  double a = asin(sina);

  cross_product(u, w, tmp_vec);
  double sinb = sqrt(norm(tmp_vec));
  double b = asin(sinb);

  cross_product(w, v, tmp_vec);
  double sinc = sqrt(norm(tmp_vec));
  double c = asin(sinc);

  double s = 0.5 * (a + b + c);

  double t = tan(s * 0.5) * tan((s - a) * 0.5) * tan((s - b) * 0.5) * tan((s - c) * 0.5);

  double area = fabs(4.0 * atan(sqrt(fabs(t))));

  return area;
}

/*
 * source code is taken from code by Robert Oehmke of Earth System Modeling
 * Framework (www.earthsystemmodeling.org) and adjusted to CDO data structures
 *
 * the license statement for this routine is as follows:
 * Earth System Modeling Framework
 * Copyright 2002-2013, University Corporation for Atmospheric Research,
 * Massachusetts Institute of Technology, Geophysical Fluid Dynamics
 * Laboratory, University of Michigan, National Centers for Environmental
 * Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
 * NASA Goddard Space Flight Center.
 * Licensed under the University of Illinois-NCSA License.
 */
static double
mod_huiliers_area(int num_corners, double *cell_corner_lon, double *cell_corner_lat)
{
  if (num_corners < 3) return 0;

  // sum areas around cell
  double sum = 0.0;
  double pnt1[3], pnt2[3], pnt3[3];

  lonlat_to_xyz(cell_corner_lon[0], cell_corner_lat[0], pnt1);
  lonlat_to_xyz(cell_corner_lon[1], cell_corner_lat[1], pnt2);

  for (int i = 2; i < num_corners; i++)
    {
      // points that make up a side of cell
      lonlat_to_xyz(cell_corner_lon[i], cell_corner_lat[i], pnt3);

      // compute angle for pnt2
      sum += mod_tri_area(pnt1, pnt2, pnt3);

      if (i < (num_corners - 1))
        {
          pnt2[0] = pnt3[0];
          pnt2[1] = pnt3[1];
          pnt2[2] = pnt3[2];
        }
    }

  return sum;
}

static double
mod_huiliers_area2(int num_corners, double *cell_corner_lon, double *cell_corner_lat, double cell_center_lon,
                   double cell_center_lat)
{
  if (num_corners < 3) return 0;

  // sum areas around cell
  double sum = 0.0;
  double pnt1[3], pnt2[3], pnt3[3];

  lonlat_to_xyz(cell_center_lon, cell_center_lat, pnt1);
  lonlat_to_xyz(cell_corner_lon[0], cell_corner_lat[0], pnt2);

  for (int i = 1; i < num_corners; i++)
    {
      if (IS_EQUAL(cell_corner_lon[i], cell_corner_lon[i - 1]) && IS_EQUAL(cell_corner_lat[i], cell_corner_lat[i - 1])) continue;

      // points that make up a side of cell
      lonlat_to_xyz(cell_corner_lon[i], cell_corner_lat[i], pnt3);

      // compute angle for pnt2
      sum += mod_tri_area(pnt1, pnt2, pnt3);

      pnt2[0] = pnt3[0];
      pnt2[1] = pnt3[1];
      pnt2[2] = pnt3[2];
    }

  if (!(IS_EQUAL(cell_corner_lon[0], cell_corner_lon[num_corners - 1])
        && IS_EQUAL(cell_corner_lat[0], cell_corner_lat[num_corners - 1])))
    {
      lonlat_to_xyz(cell_corner_lon[0], cell_corner_lat[0], pnt3);
      sum += mod_tri_area(pnt1, pnt2, pnt3);
    }

  return sum;
}

static void
getLonLatCorner(size_t nx, size_t idx, const double *grid_corner_lon, const double *grid_corner_lat, double *lons, double *lats)
{
  size_t j = idx / nx;
  size_t i = idx - j * nx;

  lons[0] = grid_corner_lon[2*i];
  lons[1] = grid_corner_lon[2*i+1];
  lons[2] = grid_corner_lon[2*i+1];
  lons[3] = grid_corner_lon[2*i];

  if ( grid_corner_lat[2*j+1] > grid_corner_lat[2*j] )
    {
      lats[0] = grid_corner_lat[2*j];
      lats[1] = grid_corner_lat[2*j];
      lats[2] = grid_corner_lat[2*j+1];
      lats[3] = grid_corner_lat[2*j+1];
    }
  else
    {
      lats[0] = grid_corner_lat[2*j+1];
      lats[1] = grid_corner_lat[2*j+1];
      lats[2] = grid_corner_lat[2*j];
      lats[3] = grid_corner_lat[2*j];
    }
}

static int
gridGenAreaReg2D(int gridID, double *area)
{
  int status = 0;

  size_t gridsize = gridInqSize(gridID);
  size_t nlon = gridInqXsize(gridID);
  size_t nlat = gridInqYsize(gridID);

  if (gridInqYvals(gridID, NULL) == 0 || gridInqXvals(gridID, NULL) == 0)
    {
      cdoWarning("Computation of grid cell area weights failed, grid cell center coordinates missing!");
      return 1;
    }

  char xunitstr[CDI_MAX_NAME];
  char yunitstr[CDI_MAX_NAME];
  gridInqXunits(gridID, xunitstr);
  gridInqYunits(gridID, yunitstr);

  std::vector<double> grid_center_lon(nlon);
  std::vector<double> grid_center_lat(nlat);

  gridInqXvals(gridID, grid_center_lon.data());
  gridInqYvals(gridID, grid_center_lat.data());

  std::vector<double> grid_corner_lon(2 * nlon);
  std::vector<double> grid_corner_lat(2 * nlat);

  if (gridInqYbounds(gridID, NULL) && gridInqXbounds(gridID, NULL))
    {
      gridInqXbounds(gridID, grid_corner_lon.data());
      gridInqYbounds(gridID, grid_corner_lat.data());
    }
  else
    {
      grid_gen_bounds(nlon, grid_center_lon.data(), grid_corner_lon.data());
      grid_gen_bounds(nlat, grid_center_lat.data(), grid_corner_lat.data());
      grid_check_lat_borders(2 * nlat, grid_corner_lat.data());
    }

  grid_to_radian(xunitstr, nlon, grid_center_lon.data(), "grid1 center longitudes");
  grid_to_radian(xunitstr, nlon * 2, grid_corner_lon.data(), "grid1 corner longitudes");

  grid_to_radian(yunitstr, nlat, grid_center_lat.data(), "grid1 center latitudes");
  grid_to_radian(yunitstr, nlat * 2, grid_corner_lat.data(), "grid1 corner latitudes");

  int nv = 4;
  double findex = 0;

  progressInit();

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(findex, gridsize, nlon, area, nv, grid_corner_lon, grid_corner_lat)
#endif
  for (size_t i = 0; i < gridsize; ++i)
    {
      int lprogress = 1;
      if (cdo_omp_get_thread_num() != 0) lprogress = 0;

#ifdef _OPENMP
#pragma omp atomic
#endif
      findex++;
      if (lprogress) progressStatus(0, 1, findex / gridsize);

      double lons[4], lats[4];
      getLonLatCorner(nlon, i, grid_corner_lon.data(), grid_corner_lat.data(), lons, lats);
      // area[i] = mod_cell_area(nv, grid_corner_lon+i*nv, grid_corner_lat+i*nv);
      area[i] = mod_huiliers_area(nv, lons, lats);
    }

  progressStatus(0, 1, 1);

  return status;
}

static int
gridGenAreaUnstruct(int gridID, double *area)
{
  int status = 0;
  bool lgrid_gen_bounds = false;
  bool lgriddestroy = false;

  size_t gridsize = gridInqSize(gridID);
  int gridtype = gridInqType(gridID);

  if (gridtype != GRID_UNSTRUCTURED && gridtype != GRID_CURVILINEAR)
    {
      if (gridtype == GRID_GME)
        {
          lgriddestroy = true;
          gridID = gridToUnstructured(gridID, 1);
          /*
          grid_mask = (int*) Malloc(gridsize*sizeof(int));
          gridInqMaskGME(gridID, grid_mask);
          if ( grid_mask ) Free(grid_mask);
          */
        }
      else
        {
          lgriddestroy = true;
          gridID = gridToCurvilinear(gridID, 1);
          lgrid_gen_bounds = true;
        }
    }

  if (gridtype == GRID_UNSTRUCTURED)
    {
      if (gridInqYvals(gridID, NULL) == 0 || gridInqXvals(gridID, NULL) == 0)
        {
          if (gridInqNumber(gridID) > 0)
            {
              lgriddestroy = true;
              gridID = referenceToGrid(gridID);
              if (gridID == -1) return 1;
            }
        }
    }

  gridtype = gridInqType(gridID);

  size_t nv = (gridtype == GRID_UNSTRUCTURED) ? gridInqNvertex(gridID) : 4;

  if (gridInqYvals(gridID, NULL) == 0 || gridInqXvals(gridID, NULL) == 0)
    {
      cdoWarning("Computation of grid cell area weights failed, grid cell center coordinates missing!");
      return 1;
    }

  if (nv == 0)
    {
      cdoWarning("Computation of grid cell area weights failed, grid cell corner coordinates missing!");
      return 1;
    }

  char xunitstr[CDI_MAX_NAME];
  char yunitstr[CDI_MAX_NAME];
  gridInqXunits(gridID, xunitstr);
  gridInqYunits(gridID, yunitstr);

  std::vector<double> grid_center_lon(gridsize);
  std::vector<double> grid_center_lat(gridsize);

  gridInqXvals(gridID, grid_center_lon.data());
  gridInqYvals(gridID, grid_center_lat.data());

  std::vector<double> grid_corner_lon(nv * gridsize);
  std::vector<double> grid_corner_lat(nv * gridsize);

  if (gridInqYbounds(gridID, NULL) && gridInqXbounds(gridID, NULL))
    {
      gridInqXbounds(gridID, grid_corner_lon.data());
      gridInqYbounds(gridID, grid_corner_lat.data());
    }
  else
    {
      if (lgrid_gen_bounds)
        {
          size_t nlon = gridInqXsize(gridID);
          size_t nlat = gridInqYsize(gridID);
          double dlon = 0;
          if (nlon == 1) dlon = 1;

          grid_cell_center_to_bounds_X2D(xunitstr, nlon, nlat, grid_center_lon.data(), grid_corner_lon.data(), dlon);
          grid_cell_center_to_bounds_Y2D(yunitstr, nlon, nlat, grid_center_lat.data(), grid_corner_lat.data());
        }
      else
        {
          status = 1;
          return status;
        }
    }

  grid_to_radian(xunitstr, gridsize, grid_center_lon.data(), "grid1 center longitudes");
  grid_to_radian(xunitstr, gridsize * nv, grid_corner_lon.data(), "grid1 corner longitudes");

  grid_to_radian(yunitstr, gridsize, grid_center_lat.data(), "grid1 center latitudes");
  grid_to_radian(yunitstr, gridsize * nv, grid_corner_lat.data(), "grid1 corner latitudes");

  if (lgriddestroy) gridDestroy(gridID);

  double findex = 0;

  progressInit();

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(findex, gridsize, area, nv, grid_corner_lon, grid_corner_lat, grid_center_lon, \
                                              grid_center_lat)
#endif
  for (size_t i = 0; i < gridsize; ++i)
    {
      int lprogress = 1;
      if (cdo_omp_get_thread_num() != 0) lprogress = 0;

#ifdef _OPENMP
#pragma omp atomic
#endif
      findex++;
      if (lprogress) progressStatus(0, 1, findex / gridsize);

      // area[i] = mod_cell_area(nv, grid_corner_lon+i*nv, grid_corner_lat+i*nv);
      if (nv <= 4)
        area[i] = mod_huiliers_area(nv, &grid_corner_lon[i * nv], &grid_corner_lat[i * nv]);
      else
        area[i]
            = mod_huiliers_area2(nv, &grid_corner_lon[i * nv], &grid_corner_lat[i * nv], grid_center_lon[i], grid_center_lat[i]);
    }

  progressStatus(0, 1, 1);

  return status;
}

int
gridGenArea(int gridID, double *area)
{
  int status = 0;

  size_t gridsize = gridInqSize(gridID);
  int gridtype = gridInqType(gridID);
  int projtype = (gridtype == GRID_PROJECTION) ? gridInqProjType(gridID) : -1;

  if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN)
    {
      status = gridGenAreaReg2D(gridID, area);
    }
  else if (projtype == CDI_PROJ_RLL || projtype == CDI_PROJ_LAEA || projtype == CDI_PROJ_SINU || projtype == CDI_PROJ_LCC
           || gridtype == GRID_GME || gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED)
    {
      status = gridGenAreaUnstruct(gridID, area);
    }
  else
    {
      cdoAbort("Internal error! Unsupported gridtype: %s", gridNamePtr(gridtype));
    }

  if (cdoVerbose)
    {
      double total_area = 0;
      for (size_t i = 0; i < gridsize; ++i) total_area += area[i];
      cdoPrint("Total area = %g steradians", total_area);
    }

  if (gridsize < 20)
    {
      double total_area = 0;
      for (size_t i = 0; i < gridsize; ++i) total_area += area[i];
      int nzero = 0;
      for (size_t i = 0; i < gridsize; ++i)
        if (IS_EQUAL(area[i], 0.)) nzero++;
      if (IS_EQUAL(total_area, 0.)) status = 2;
    }

  return status;
}

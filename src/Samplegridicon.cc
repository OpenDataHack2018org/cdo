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
#include "grid_search.h"

extern "C"
{
#include "clipping/geometry.h"
}

constexpr int MAX_CHILDS = 9;

typedef struct
{
  long ncells;
  long *neighbor;  // neighbor cell index
  long *parent;    // parent cell index
  long *child;     // child cell index
  const char *filename;
} cellindex_type;

static void
copy_data_to_index(long ncells, const double *restrict data, long *restrict cellindex)
{
  for (long i = 0; i < ncells; ++i)
    cellindex[i] = lround(data[i]);
}

static void
free_cellindex(cellindex_type *cellindex)
{
  if (cellindex->neighbor) Free(cellindex->neighbor);
  if (cellindex->parent) Free(cellindex->parent);
  if (cellindex->child) Free(cellindex->child);
  Free(cellindex);
}

static cellindex_type *
read_cellindex(const char *filename)
{
  openLock();
  int streamID = streamOpenRead(filename);
  openUnlock();

  if (streamID < 0) cdiOpenError(streamID, "Open failed on >%s<", filename);

  int vlistID = streamInqVlist(streamID);
  int ngrids = vlistNgrids(vlistID);
  int gridID = -1;
  for (int index = 0; index < ngrids; ++index)
    {
      gridID = vlistGrid(vlistID, index);
      if (gridInqType(gridID) == GRID_UNSTRUCTURED && gridInqNvertex(gridID) == 3) break;
    }

  if (gridID == -1) cdoAbort("No ICON grid found in %s!", filename);

  int nid = CDI_UNDEFID;
  int pid = CDI_UNDEFID;
  // int cid = CDI_UNDEFID;
  int nvars = vlistNvars(vlistID);
  char varname[CDI_MAX_NAME];
  for (int varID = 0; varID < nvars; ++varID)
    {
      vlistInqVarName(vlistID, varID, varname);
      if (strcmp(varname, "neighbor_cell_index") == 0)
        nid = varID;
      else if (strcmp(varname, "parent_cell_index") == 0)
        pid = varID;
      // else if ( strcmp(varname, "child_cell_index") == 0 ) cid = varID;
    }

  if (nid == CDI_UNDEFID) cdoAbort("neighbor_cell_index not found in %s!", filename);
  if (pid == CDI_UNDEFID) cdoAbort("parent_cell_index not found in %s!", filename);
  // if ( cid == CDI_UNDEFID ) cdoAbort("child_cell_index not found in %s!",
  // filename);

  long ncells = gridInqSize(gridID);

  cellindex_type *cellindex = (cellindex_type *) Malloc(sizeof(cellindex_type));
  cellindex->ncells = ncells;

  cellindex->neighbor = NULL;
  // cellindex->neighbor = (long*) Malloc(3*ncells*sizeof(long));
  cellindex->parent = (long *) Malloc(ncells * sizeof(long));
  cellindex->child = NULL;
  // cellindex->child    = (cid != CDI_UNDEFID) ? (int*)
  // Malloc(MAX_CHILDS*ncells*sizeof(int)) : NULL;
  double *data = (double *) Malloc(ncells * sizeof(double));

  int nrecs = streamInqTimestep(streamID, 0);
  for (int recID = 0; recID < nrecs; recID++)
    {
      int varID, levelID;
      size_t nmiss;
      streamInqRecord(streamID, &varID, &levelID);
      if (varID == pid /* || varID == nid || varID == cid */)
        {
          streamReadRecord(streamID, data, &nmiss);
          if (varID == pid) copy_data_to_index(ncells, data, cellindex->parent);
          // else if ( varID == nid ) copy_data_to_index(ncells, data,
          // cellindex->neighbor+levelID*ncells); else if ( varID == cid )
          // copy_data_to_index(ncells, data, cellindex->child+levelID*ncells);
        }
    }

  // Fortran to C index
  for (long i = 0; i < ncells; ++i)
    cellindex->parent[i] -= 1;
  // for ( long i = 0; i < 3*ncells; ++i ) cellindex->neighbor[i] -= 1;

  streamClose(streamID);

  Free(data);

  return cellindex;
}

static int
read_grid(const char *filename)
{
  openLock();
  int streamID = streamOpenRead(filename);
  openUnlock();

  if (streamID < 0) cdiOpenError(streamID, "Open failed on >%s<", filename);

  int vlistID = streamInqVlist(streamID);
  int ngrids = vlistNgrids(vlistID);
  int gridID = -1;
  for (int index = 0; index < ngrids; ++index)
    {
      gridID = vlistGrid(vlistID, index);
      if (gridInqType(gridID) == GRID_UNSTRUCTURED && gridInqNvertex(gridID) == 3) break;
    }

  if (gridID == -1) cdoAbort("No ICON grid found in %s!", filename);

  int gridID2 = gridDuplicate(gridID);

  streamClose(streamID);

  return gridID2;
}

/**
* Return the first index of element x fits.
*
* If no interval can be found return -1.

* @param *array ascending sorted list
* @param n      length of the sorted list
* @param search the element to find a position for
*/
static long
find_index(int search, long n, const long *restrict array)
{
  long first = 0;
  long last = n - 1;
  long middle = (first + last) / 2;

  while (first <= last)
    {
      if (array[middle] < search)
        first = middle + 1;
      else if (array[middle] == search)
        {
          for (long i = middle; i >= 0; i--)
            {
              if (array[i] == search)
                middle = i;
              else
                break;
            }
          return middle;
        }
      else
        last = middle - 1;

      middle = (first + last) / 2;
    }

  return -1;
}

typedef struct
{
  int p, i;
} sinfo_t;

static int
cmpsinfo(const void *s1, const void *s2)
{
  int cmp = 0;
  const sinfo_t *x = (const sinfo_t *) s1;
  const sinfo_t *y = (const sinfo_t *) s2;

  if (x->p < y->p)
    cmp = -1;
  else if (x->p > y->p)
    cmp = 1;

  return cmp;
}

static void
compute_child_from_parent(cellindex_type *cellindex1, cellindex_type *cellindex2)
{
  long ncells1 = cellindex1->ncells;
  long *parent1 = cellindex1->parent;

  long *idx1 = (long *) Malloc(ncells1 * sizeof(long));
  for (long i = 0; i < ncells1; ++i)
    idx1[i] = i;
  for (long i = 1; i < ncells1; ++i)
    if (parent1[i] < parent1[i - 1])
      {
        if (cdoVerbose) cdoPrint("Sort parent index of %s!", cellindex1->filename);
        sinfo_t *sinfo = (sinfo_t *) Malloc(ncells1 * sizeof(sinfo_t));
        for (long j = 0; j < ncells1; ++j)
          {
            sinfo[j].p = parent1[j];
            sinfo[j].i = idx1[j];
          }
        qsort(sinfo, ncells1, sizeof(sinfo_t), cmpsinfo);
        for (long j = 0; j < ncells1; ++j)
          {
            parent1[j] = sinfo[j].p;
            idx1[j] = sinfo[j].i;
          }
        Free(sinfo);
        break;
      }

  long ncells2 = cellindex2->ncells;
  long *child2 = (long *) Malloc(MAX_CHILDS * ncells2 * sizeof(long));
  cellindex2->child = child2;
  for (long i = 0; i < ncells2; ++i)
    {
      for (long k = 0; k < MAX_CHILDS; ++k)
        child2[i * MAX_CHILDS + k] = -1;
      long j = find_index(i, ncells1, parent1);
      if (j < 0) continue;
      for (long k = 0; k < MAX_CHILDS; ++k)
        {
          if (i != parent1[j + k]) break;
          //  child2[i*MAX_CHILDS+k] = j+k;
          child2[i * MAX_CHILDS + k] = idx1[j + k];
        }
      // if ( i%10000 == 0 ) printf("%d %d %d %d %d %d\n", i, j, parent1[j],
      // parent1[j+1], parent1[j+2], parent1[j+3]);
    }
  Free(idx1);
}

static void
read_coordinates(const char *filename, long n, double *lon, double *lat, int nv, double *lon_bnds, double *lat_bnds)
{
  openLock();
  int streamID = streamOpenRead(filename);
  openUnlock();

  if (streamID < 0) cdiOpenError(streamID, "Open failed on >%s<", filename);

  int vlistID = streamInqVlist(streamID);
  int ngrids = vlistNgrids(vlistID);
  int gridID = -1;
  for (int index = 0; index < ngrids; ++index)
    {
      gridID = vlistGrid(vlistID, index);
      if (gridInqType(gridID) == GRID_UNSTRUCTURED && (long) gridInqSize(gridID) == n && gridInqNvertex(gridID) == 3)
        break;
    }

  if (gridID == -1) cdoAbort("No ICON grid with %ld cells found in %s!", n, filename);

  gridInqXvals(gridID, lon);
  gridInqYvals(gridID, lat);

  char units[CDI_MAX_NAME];
  /* Convert lat/lon units if required */
  gridInqXunits(gridID, units);
  grid_to_radian(units, n, lon, "grid center lon");
  gridInqYunits(gridID, units);
  grid_to_radian(units, n, lat, "grid center lat");

  if (nv == 3 && lon_bnds && lat_bnds)
    {
      gridInqXbounds(gridID, lon_bnds);
      gridInqYbounds(gridID, lat_bnds);

      gridInqXunits(gridID, units);
      grid_to_radian(units, n * 3, lon_bnds, "grid corner lon");
      gridInqYunits(gridID, units);
      grid_to_radian(units, n * 3, lat_bnds, "grid corner lat");
    }

  streamClose(streamID);
}

int grid_search_nbr(struct gridsearch *gs, nbrWeightsType &nbrWeights, double plon, double plat);

int find_coordinate_to_ignore(double *cell_corners_xyz);
double calculate_the_polygon_area(double cell_corners[], int number_corners);
bool are_polygon_vertices_arranged_in_clockwise_order(double cell_area);
int winding_numbers_algorithm(double cell_corners[], int number_corners, double point[]);

#define MAX_SEARCH 128  // the triangles are distorted!

static void
compute_child_from_bounds(cellindex_type *cellindex2, long ncells2, double *grid_center_lon2, double *grid_center_lat2,
                          double *grid_corner_lon2, double *grid_corner_lat2, long ncells1, double *grid_center_lon1,
                          double *grid_center_lat1)
{
  bool xIsCyclic = false;
  size_t dims[2];
  dims[0] = ncells1;
  dims[1] = 0;
  struct gridsearch *gs = gridsearch_create(xIsCyclic, dims, ncells1, grid_center_lon1, grid_center_lat1);
  nbrWeightsType nbrWeights(MAX_SEARCH);
  size_t *nbr_addr = &nbrWeights.m_addr[0];

  int ncorner = 3;
  double corner_coordinates[3];
  double center_point_xyz[3];
  double cell_corners_xyz[12];
  double cell_corners_plane_projection[8];
  double center_point_plane_projection[2];

  long *child2 = (long *) Malloc(MAX_CHILDS * ncells2 * sizeof(long));
  cellindex2->child = child2;
  for (long cell_no2 = 0; cell_no2 < ncells2; ++cell_no2)
    {
      for (int k = 0; k < MAX_CHILDS; ++k)
        child2[cell_no2 * MAX_CHILDS + k] = -1;

      for (int corner_no = 0; corner_no < ncorner; corner_no++)
        {
          /* Conversion of corner spherical coordinates to Cartesian
           * coordinates. */

          LLtoXYZ(grid_corner_lon2[cell_no2 * ncorner + corner_no], grid_corner_lat2[cell_no2 * ncorner + corner_no],
                  corner_coordinates);

          /* The components of the result vector are appended to the list of
           * cell corner coordinates. */

          int off = corner_no * 3;
          cell_corners_xyz[off + 0] = corner_coordinates[0];
          cell_corners_xyz[off + 1] = corner_coordinates[1];
          cell_corners_xyz[off + 2] = corner_coordinates[2];
        }
      cell_corners_xyz[ncorner * 3 + 0] = cell_corners_xyz[0];
      cell_corners_xyz[ncorner * 3 + 1] = cell_corners_xyz[1];
      cell_corners_xyz[ncorner * 3 + 2] = cell_corners_xyz[2];

      int coordinate_to_ignore = find_coordinate_to_ignore(cell_corners_xyz);

      bool invert_result = false;
      if (cell_corners_xyz[coordinate_to_ignore - 1] < 0) invert_result = true;

      switch (coordinate_to_ignore)
        {
        case 1:
          for (int corner_no = 0; corner_no <= ncorner; corner_no++)
            {
              cell_corners_plane_projection[corner_no * 2 + 0] = cell_corners_xyz[corner_no * 3 + 1];
              cell_corners_plane_projection[corner_no * 2 + 1] = cell_corners_xyz[corner_no * 3 + 2];
            }
          break;
        case 2:
          for (int corner_no = 0; corner_no <= ncorner; corner_no++)
            {
              cell_corners_plane_projection[corner_no * 2 + 0] = cell_corners_xyz[corner_no * 3 + 2];
              cell_corners_plane_projection[corner_no * 2 + 1] = cell_corners_xyz[corner_no * 3 + 0];
            }
          break;
        case 3:
          for (int corner_no = 0; corner_no <= ncorner; corner_no++)
            {
              cell_corners_plane_projection[corner_no * 2 + 0] = cell_corners_xyz[corner_no * 3 + 0];
              cell_corners_plane_projection[corner_no * 2 + 1] = cell_corners_xyz[corner_no * 3 + 1];
            }
          break;
        }

      double polygon_area = calculate_the_polygon_area(cell_corners_plane_projection, ncorner + 1);
      bool is_clockwise = are_polygon_vertices_arranged_in_clockwise_order(polygon_area);

      if (invert_result) is_clockwise = !is_clockwise;
      if (is_clockwise) continue;

      grid_search_nbr(gs, nbrWeights, grid_center_lon2[cell_no2], grid_center_lat2[cell_no2]);
      int k = 0;

      for (int i = 0; i < MAX_SEARCH; ++i)
        {
          size_t cell_no1 = nbr_addr[i];
          if (cell_no1 < SIZE_MAX)
            {
              LLtoXYZ(grid_center_lon1[cell_no1], grid_center_lat1[cell_no1], center_point_xyz);

              switch (coordinate_to_ignore)
                {
                case 1:
                  center_point_plane_projection[0] = center_point_xyz[1];
                  center_point_plane_projection[1] = center_point_xyz[2];
                  break;
                case 2:
                  center_point_plane_projection[0] = center_point_xyz[2];
                  center_point_plane_projection[1] = center_point_xyz[0];
                  break;
                case 3:
                  center_point_plane_projection[0] = center_point_xyz[0];
                  center_point_plane_projection[1] = center_point_xyz[1];
                  break;
                }

              int winding_number = winding_numbers_algorithm(cell_corners_plane_projection, ncorner + 1,
                                                             center_point_plane_projection);
              // printf("%d %g %g %g %g %d\n", cell_no2,
              // RAD2DEG*grid_center_lon2[cell_no2],
              // RAD2DEG*grid_center_lat2[cell_no2],
              // RAD2DEG*grid_center_lon1[cell_no1],
              // RAD2DEG*grid_center_lat1[cell_no1], winding_number);
              if (winding_number != 0)
                {
                  if (k >= MAX_CHILDS)
                    cdoAbort("Internal problem, limit of MAX_CHILDS reached "
                             "(limit=9).");
                  child2[cell_no2 * MAX_CHILDS + k++] = (long) cell_no1;
                }
            }
        }
      // printf("%d k = %d\n", cell_no2, k);
    }
}

static void
compute_child_from_coordinates(cellindex_type *cellindex1, cellindex_type *cellindex2)
{
  long ncells1 = cellindex1->ncells;
  long ncells2 = cellindex2->ncells;

  double *lon1 = (double *) Malloc(ncells1 * sizeof(double));
  double *lat1 = (double *) Malloc(ncells1 * sizeof(double));
  double *lon2 = (double *) Malloc(ncells2 * sizeof(double));
  double *lat2 = (double *) Malloc(ncells2 * sizeof(double));
  double *lon2_bnds = (double *) Malloc(3 * ncells2 * sizeof(double));
  double *lat2_bnds = (double *) Malloc(3 * ncells2 * sizeof(double));

  read_coordinates(cellindex1->filename, ncells1, lon1, lat1, 0, NULL, NULL);
  read_coordinates(cellindex2->filename, ncells2, lon2, lat2, 3, lon2_bnds, lat2_bnds);

  compute_child_from_bounds(cellindex2, ncells2, lon2, lat2, lon2_bnds, lat2_bnds, ncells1, lon1, lat1);

  Free(lon1);
  Free(lat1);
  Free(lon2);
  Free(lat2);
  Free(lon2_bnds);
  Free(lat2_bnds);
}

static void
compute_child(cellindex_type *cellindex1, cellindex_type *cellindex2)
{
  bool lparent = true;
  long ncells1 = cellindex1->ncells;
  long *parent1 = cellindex1->parent;
  {
    long i;
    for (i = 0; i < ncells1; ++i)
      if (parent1[i] >= 0) break;
    if (i == ncells1) lparent = false;
  }

  if (lparent)
    compute_child_from_parent(cellindex1, cellindex2);
  else
    {
      compute_child_from_coordinates(cellindex1, cellindex2);
      // cdoAbort("Missing parent index of %s!", cellindex1->filename);
    }
}

static void
compute_sum(long i, long *n, double *sum, double *sumq, long kci, cellindex_type **cellindex, double *array)
{
  // printf("compute: i, kci %d %d\n", i, kci);
  long ncells2 = cellindex[kci]->ncells;
  if (i < 0 || i > ncells2) cdoAbort("Child grid cell index %ld out of bounds %ld!", i, ncells2);

  for (int k = 0; k < MAX_CHILDS; ++k)
    {
      long index = cellindex[kci]->child[i * MAX_CHILDS + k];
      if (index == -1) break;
      if (kci == 1)
        {
          *sum += array[index];
          *sumq += array[index] * array[index];
          *n += 1;
        }
      else
        compute_sum(index, n, sum, sumq, kci - 1, cellindex, array);
    }
}

static void
samplegrid(double missval, long nci, cellindex_type **cellindex, double *array1, double *array2, double *array3)
{
  static bool lstat = true;
  long kci = nci - 1;
  long ncells2 = cellindex[kci]->ncells;
  long nx = 0;
  double x = 0;
#ifdef _OPENMP
  //#pragma omp parallel for default(none) shared(missval, ncells2, kci,
  // cellindex, array1, array2, array3)
#endif
  for (long i = 0; i < ncells2; ++i)
    {
      long n = 0;
      double sum = 0, sumq = 0;
      compute_sum(i, &n, &sum, &sumq, kci, cellindex, array1);
      array2[i] = n ? sum / n : missval;  // mean
      double var1 = (n * n > n) ? (sumq * n - sum * sum) / (n * n - n) : missval;
      if (var1 < 0 && var1 > -1.e-5) var1 = 0;
      array3[i] = varToStd(var1, missval);  // std1
      if (lstat && n)
        {
          nx++;
          x += n;
        }
    }
  if (cdoVerbose && lstat)
    {
      lstat = false;
      cdoPrint("Mean number of childs %g", nx ? x / nx : 0);
    }
}

#include "pstream_int.h"

void *
Samplegridicon(void *process)
{
  int nrecs;
  int varID, levelID;

  cdoInitialize(process);

  cdoOperatorAdd("samplegridicon", 0, 0, "sample grids");

  int nsamplegrids = operatorArgc();
  if (nsamplegrids < 2) cdoAbort("Parameter missing!");

  std::vector<cellindex_type *> cellindex(nsamplegrids);

  for (int i = 0; i < nsamplegrids; ++i)
    {
      cellindex[i] = read_cellindex(operatorArgv()[i]);
      cellindex[i]->filename = operatorArgv()[i];
      if (cdoVerbose) cdoPrint("Found %ld grid cells in %s", cellindex[i]->ncells, cellindex[i]->filename);
    }

  for (int i = 0; i < nsamplegrids - 1; ++i)
    compute_child(cellindex[i], cellindex[i + 1]);

  int gridID2 = read_grid(operatorArgv()[nsamplegrids - 1]);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);

  long gridsize = vlistGridsizeMax(vlistID1);
  if (cdoVerbose) cdoPrint("Source gridsize = %zu", gridsize);
  if (gridsize != cellindex[0]->ncells)
    cdoAbort("Gridsize (%ld) of input stream and first grid (%ld) differ!", gridsize, cellindex[0]->ncells);
  if (vlistNumber(vlistID1) != CDI_REAL) gridsize *= 2;
  double *array1 = (double *) Malloc(gridsize * sizeof(double));

  int vlistID2 = vlistDuplicate(vlistID1);
  int vlistID3 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  int taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);
  vlistDefTaxis(vlistID3, taxisID3);

  int ngrids = vlistNgrids(vlistID1);
  for (int index = 0; index < ngrids; ++index)
    {
      int gridID = vlistGrid(vlistID1, index);
      int gridtype = gridInqType(gridID);
      if (!(gridtype == GRID_UNSTRUCTURED && gridInqNvertex(gridID) == 3))
        cdoAbort("Unsupported gridtype: %s with %d corners", gridNamePtr(gridtype), gridInqNvertex(gridID));

      vlistChangeGridIndex(vlistID2, index, gridID2);
      vlistChangeGridIndex(vlistID3, index, gridID2);
    }

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int streamID3 = cdoStreamOpenWrite(cdoStreamName(2), cdoFiletype());
  pstreamDefVlist(streamID3, vlistID3);

  long gridsize2 = gridInqSize(gridID2);
  if (cdoVerbose) cdoPrint("Target gridsize = %ld", gridsize2);
  if (vlistNumber(vlistID2) != CDI_REAL) gridsize2 *= 2;
  double *array2 = (double *) Malloc(gridsize2 * sizeof(double));
  double *array3 = (double *) Malloc(gridsize2 * sizeof(double));

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      taxisCopyTimestep(taxisID3, taxisID1);
      pstreamDefTimestep(streamID2, tsID);
      pstreamDefTimestep(streamID3, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          size_t nmiss;
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamReadRecord(streamID1, array1, &nmiss);

          double missval = vlistInqVarMissval(vlistID1, varID);

          samplegrid(missval, nsamplegrids, &cellindex[0], array1, array2, array3);

          nmiss = arrayNumMV(gridsize2, array2, missval);
          pstreamDefRecord(streamID2, varID, levelID);
          pstreamWriteRecord(streamID2, array2, nmiss);

          nmiss = arrayNumMV(gridsize2, array3, missval);
          pstreamDefRecord(streamID3, varID, levelID);
          pstreamWriteRecord(streamID3, array3, nmiss);
        }

      tsID++;
    }

  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  vlistDestroy(vlistID2);
  gridDestroy(gridID2);

  if (array3) Free(array3);
  if (array2) Free(array2);
  if (array1) Free(array1);

  for (int i = 0; i < nsamplegrids; ++i)
    free_cellindex(cellindex[i]);

  cdoFinish();

  return 0;
}

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
#include "pstream_int.h"

#define MAX_BLOCKS 65536

static void
genGrids(int gridID1, int *gridIDs, size_t nxvals, size_t nyvals, size_t nxblocks, size_t nyblocks, size_t **gridindex,
         size_t *ogridsize, size_t nsplit)
{
  double *xpvals = NULL, *ypvals = NULL;
  int gridtype = gridInqType(gridID1);
  bool lregular = true;
  if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_GENERIC)
    lregular = true;
  else if (gridtype == GRID_CURVILINEAR)
    lregular = false;
  else
    cdoAbort("Unsupported grid type: %s!", gridNamePtr(gridtype));

  size_t nx = gridInqXsize(gridID1);
  size_t ny = gridInqYsize(gridID1);

  bool lxcoord = gridInqXvals(gridID1, NULL) > 0;
  bool lycoord = gridInqYvals(gridID1, NULL) > 0;

  std::vector<double> xvals, yvals;
  std::vector<double> xvals2, yvals2;

  if (lxcoord)
    {
      if (lregular)
        {
          xvals.resize(nx);
        }
      else
        {
          xvals.resize(nx * ny);
          xvals2.resize(nxvals * nyvals);
        }

      gridInqXvals(gridID1, &xvals[0]);
    }

  if (lycoord)
    {
      if (lregular)
        {
          yvals.resize(ny);
        }
      else
        {
          yvals.resize(nx * ny);
          yvals2.resize(nxvals * nyvals);
        }

      gridInqYvals(gridID1, &yvals[0]);
    }

  std::vector<size_t> xlsize(nxblocks);
  std::vector<size_t> ylsize(nyblocks);

  for (size_t ix = 0; ix < nxblocks; ++ix) xlsize[ix] = nxvals;
  if (nx % nxblocks != 0) xlsize[nxblocks - 1] = nx - (nxblocks - 1) * nxvals;
  if (cdoVerbose)
    for (size_t ix = 0; ix < nxblocks; ++ix) cdoPrint("xblock %zu: %zu", ix, xlsize[ix]);

  for (size_t iy = 0; iy < nyblocks; ++iy) ylsize[iy] = nyvals;
  if (ny % nyblocks != 0) ylsize[nyblocks - 1] = ny - (nyblocks - 1) * nyvals;
  if (cdoVerbose)
    for (size_t iy = 0; iy < nyblocks; ++iy) cdoPrint("yblock %zu: %zu", iy, ylsize[iy]);

  size_t index = 0;
  for (size_t iy = 0; iy < nyblocks; ++iy)
    for (size_t ix = 0; ix < nxblocks; ++ix)
      {
        size_t offset = iy * nyvals * nx + ix * nxvals;
        size_t gridsize2 = xlsize[ix] * ylsize[iy];
        gridindex[index] = (size_t *) Malloc(gridsize2 * sizeof(size_t));

        gridsize2 = 0;
        // printf("iy %d, ix %d offset %d\n", iy, ix,  offset);
        for (size_t j = 0; j < ylsize[iy]; ++j)
          for (size_t i = 0; i < xlsize[ix]; ++i)
            {
              // printf(">> %d %d %d\n", j, i, offset + j*nx + i);
              if (!lregular)
                {
                  if (lxcoord) xvals2[gridsize2] = xvals[offset + j * nx + i];
                  if (lycoord) yvals2[gridsize2] = yvals[offset + j * nx + i];
                }
              gridindex[index][gridsize2++] = offset + j * nx + i;
            }
        // printf("gridsize2 %d\n", gridsize2);

        int gridID2 = gridCreate(gridtype, gridsize2);
        gridDefXsize(gridID2, xlsize[ix]);
        gridDefYsize(gridID2, ylsize[iy]);

        gridDefNP(gridID2, gridInqNP(gridID1));
        gridDefDatatype(gridID2, gridInqDatatype(gridID1));

        grid_copy_attributes(gridID1, gridID2);

        if (gridtype == GRID_PROJECTION) grid_copy_mapping(gridID1, gridID2);

        if (lregular)
          {
            if (lxcoord) gridDefXvals(gridID2, &xvals[ix * nxvals]);
            if (lycoord) gridDefYvals(gridID2, &yvals[iy * nyvals]);
          }
        else
          {
            if (lxcoord) gridDefXvals(gridID2, &xvals2[0]);
            if (lycoord) gridDefYvals(gridID2, &yvals2[0]);
          }

        int projID1 = gridInqProj(gridID1);
        if (projID1 != CDI_UNDEFID && gridInqType(projID1) == GRID_PROJECTION)
          {
            int projID2 = gridCreate(GRID_PROJECTION, gridsize2);
            gridDefXsize(projID2, xlsize[ix]);
            gridDefYsize(projID2, ylsize[iy]);

            grid_copy_attributes(projID1, projID2);
            grid_copy_mapping(projID1, projID2);

            bool lxpcoord = gridInqXvals(projID1, NULL) > 0;
            if (lxpcoord)
              {
                if (!xpvals)
                  {
                    xpvals = (double *) Malloc(nx * sizeof(double));
                    gridInqXvals(projID1, xpvals);
                  }
                gridDefXvals(projID2, xpvals + ix * nxvals);
              }
            bool lypcoord = gridInqYvals(projID1, NULL) > 0;
            if (lypcoord)
              {
                if (!ypvals)
                  {
                    ypvals = (double *) Malloc(ny * sizeof(double));
                    gridInqYvals(projID1, ypvals);
                  }
                gridDefYvals(projID2, ypvals + iy * nyvals);
              }

            gridDefProj(gridID2, projID2);
          }

        gridIDs[index] = gridID2;
        ogridsize[index] = gridsize2;

        index++;
        if (index > nsplit) cdoAbort("Internal problem, index exceeded bounds!");
      }

  if (xpvals) Free(xpvals);
  if (ypvals) Free(ypvals);
}

static void
window_cell(double *array1, double *array2, size_t gridsize2, size_t *cellidx)
{
  for (size_t i = 0; i < gridsize2; ++i) array2[i] = array1[cellidx[i]];
}

struct sgridType
{
  int gridID;
  int *gridIDs;
  size_t *gridsize;
  size_t **gridindex;
};

void *
Distgrid(void *process)
{
  int gridID1;
  int varID, levelID;
  int nrecs;
  char filesuffix[32];
  char filename[8192];
  int index;
  int gridtype = -1;
  size_t nmiss;

  cdoInitialize(process);

  operatorInputArg("nxblocks, [nyblocks]");
  if (operatorArgc() < 1) cdoAbort("Too few arguments!");
  if (operatorArgc() > 2) cdoAbort("Too many arguments!");
  size_t nxblocks = parameter2int(operatorArgv()[0]);
  size_t nyblocks = 1;
  if (operatorArgc() == 2) nyblocks = parameter2int(operatorArgv()[1]);

  if (nxblocks == 0) cdoAbort("nxblocks has to be greater than 0!");
  if (nyblocks == 0) cdoAbort("nyblocks has to be greater than 0!");

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));
  int vlistID1 = cdoStreamInqVlist(streamID1);

  int ngrids = vlistNgrids(vlistID1);

  for (index = 0; index < ngrids; index++)
    {
      gridID1 = vlistGrid(vlistID1, index);
      gridtype = gridInqType(gridID1);
      if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_CURVILINEAR
          || (gridtype == GRID_GENERIC && gridInqXsize(gridID1) > 0 && gridInqYsize(gridID1) > 0))
        break;
    }

  if (index == ngrids)
    cdoAbort("No Lon/Lat, Gaussian, curvilinear or generic grid found (%s data unsupported)!", gridNamePtr(gridtype));

  gridID1 = vlistGrid(vlistID1, 0);
  size_t gridsize = gridInqSize(gridID1);
  size_t nx = gridInqXsize(gridID1);
  size_t ny = gridInqYsize(gridID1);
  for (int i = 1; i < ngrids; i++)
    {
      gridID1 = vlistGrid(vlistID1, i);
      if (gridsize != gridInqSize(gridID1)) cdoAbort("Gridsize must not change!");
    }

  if (nxblocks > nx)
    {
      cdoPrint("nxblocks (%zu) greater than nx (%zu), set to %zu!", nxblocks, nx, nx);
      nxblocks = nx;
    }
  if (nyblocks > ny)
    {
      cdoPrint("nyblocks (%zu) greater than ny (%zu), set to %zu!", nyblocks, ny, ny);
      nyblocks = ny;
    }

  size_t xinc = nx / nxblocks;
  size_t yinc = ny / nyblocks;

  if (nx % xinc && nx % (xinc + 1)) xinc++;
  if (ny % yinc && ny % (yinc + 1)) yinc++;

  size_t nsplit = nxblocks * nyblocks;
  if (nsplit > MAX_BLOCKS) cdoAbort("Too many blocks (max = %d)!", MAX_BLOCKS);

  std::vector<double> array1(gridsize);
  std::vector<int> vlistIDs(nsplit);
  std::vector<int> streamIDs(nsplit);

  sgridType *grids = (sgridType *) Malloc(ngrids * sizeof(sgridType));
  for (int i = 0; i < ngrids; i++)
    {
      grids[i].gridID = vlistGrid(vlistID1, i);
      grids[i].gridIDs = (int *) Malloc(nsplit * sizeof(int));
      grids[i].gridsize = (size_t *) Malloc(nsplit * sizeof(size_t));
      grids[i].gridindex = (size_t **) Malloc(nsplit * sizeof(size_t *));

      for (size_t index = 0; index < nsplit; index++) grids[i].gridindex[index] = NULL;
    }

  for (size_t index = 0; index < nsplit; index++) vlistIDs[index] = vlistDuplicate(vlistID1);

  if (cdoVerbose) cdoPrint("ngrids=%d  nsplit=%zu", ngrids, nsplit);

  for (int i = 0; i < ngrids; i++)
    {
      gridID1 = vlistGrid(vlistID1, i);
      genGrids(gridID1, &grids[i].gridIDs[0], xinc, yinc, nxblocks, nyblocks, grids[i].gridindex, grids[i].gridsize, nsplit);
      /*
      if ( cdoVerbose )
        for ( size_t index = 0; index < nsplit; index++ )
          cdoPrint("Block %d,  gridID %d,  gridsize %zu", index+1,
      grids[i].gridIDs[index], gridInqSize(grids[i].gridIDs[index]));
      */
      for (size_t index = 0; index < nsplit; index++) vlistChangeGridIndex(vlistIDs[index], i, grids[i].gridIDs[index]);
    }

  size_t gridsize2max = 0;
  for (size_t index = 0; index < nsplit; index++)
    if (grids[0].gridsize[index] > gridsize2max) gridsize2max = grids[0].gridsize[index];

  std::vector<double> array2(gridsize2max);

  strcpy(filename, cdoGetObase());
  int nchars = strlen(filename);

  const char *refname = cdoGetObase();
  filesuffix[0] = 0;
  cdoGenFileSuffix(filesuffix, sizeof(filesuffix), pstreamInqFiletype(streamID1), vlistID1, refname);

  for (size_t index = 0; index < nsplit; index++)
    {
      sprintf(filename + nchars, "%05ld", (long) index);
      if (filesuffix[0]) sprintf(filename + nchars + 5, "%s", filesuffix);

      streamIDs[index] = cdoStreamOpenWrite(filename, cdoFiletype());
      pstreamDefVlist(streamIDs[index], vlistIDs[index]);
    }

  if (ngrids > 1) cdoPrint("Baustelle: number of different grids > 1!");
  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      for (size_t index = 0; index < nsplit; index++) pstreamDefTimestep(streamIDs[index], tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamReadRecord(streamID1, &array1[0], &nmiss);

          double missval = vlistInqVarMissval(vlistID1, varID);

          for (size_t index = 0; index < nsplit; index++)
            {
              int i = 0;
              window_cell(&array1[0], &array2[0], grids[i].gridsize[index], grids[i].gridindex[index]);
              pstreamDefRecord(streamIDs[index], varID, levelID);
              if (nmiss > 0) nmiss = arrayNumMV(grids[i].gridsize[index], &array2[0], missval);
              pstreamWriteRecord(streamIDs[index], &array2[0], nmiss);
            }
        }

      tsID++;
    }

  pstreamClose(streamID1);

  for (size_t index = 0; index < nsplit; index++)
    {
      pstreamClose(streamIDs[index]);
      vlistDestroy(vlistIDs[index]);
    }

  for (int i = 0; i < ngrids; i++)
    {
      for (size_t index = 0; index < nsplit; index++) gridDestroy(grids[i].gridIDs[index]);
      Free(grids[i].gridIDs);
      Free(grids[i].gridsize);

      for (size_t index = 0; index < nsplit; index++) Free(grids[i].gridindex[index]);
      Free(grids[i].gridindex);
    }
  Free(grids);

  cdoFinish();

  return 0;
}

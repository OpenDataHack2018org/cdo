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
#include "pstream_int.h"
#include "grid.h"
#include "util_files.h"

struct ensfileType
{
  int streamID;
  int vlistID;
  int gridID;
  size_t nmiss;
  size_t gridsize;
  std::vector<long> gridindex;
  std::vector<double> array;
};

struct xyinfoType
{
  double x, y;
  int id;
};

static int
cmpx(const void *a, const void *b)
{
  const double x = ((const xyinfoType *) a)->x;
  const double y = ((const xyinfoType *) b)->x;
  return ((x > y) - (x < y)) * 2 + (x > y) - (x < y);
}

static int
cmpxy_lt(const void *a, const void *b)
{
  const double x1 = ((const xyinfoType *) a)->x;
  const double x2 = ((const xyinfoType *) b)->x;
  const double y1 = ((const xyinfoType *) a)->y;
  const double y2 = ((const xyinfoType *) b)->y;

  int cmp = 0;
  if (y1 < y2 || (!(fabs(y1 - y2) > 0) && x1 < x2))
    cmp = -1;
  else if (y1 > y2 || (!(fabs(y1 - y2) > 0) && x1 > x2))
    cmp = 1;

  return cmp;
}

static int
cmpxy_gt(const void *a, const void *b)
{
  const double x1 = ((const xyinfoType *) a)->x;
  const double x2 = ((const xyinfoType *) b)->x;
  const double y1 = ((const xyinfoType *) a)->y;
  const double y2 = ((const xyinfoType *) b)->y;

  int cmp = 0;
  if (y1 > y2 || (!(fabs(y1 - y2) > 0) && x1 < x2))
    cmp = -1;
  else if (y1 < y2 || (!(fabs(y1 - y2) > 0) && x1 > x2))
    cmp = 1;

  return cmp;
}

static int
genGrid(int ngrids, int nfiles, std::vector<ensfileType> &ef, bool ginit, int igrid, long nxblocks)
{
  bool lsouthnorth = true;
  bool lregular = false;
  bool lcurvilinear = false;
  int gridID2 = -1;

  long nx = -1;
  if (nxblocks != -1) nx = nxblocks;

  int gridID = vlistGrid(ef[0].vlistID, igrid);
  int gridtype = gridInqType(gridID);
  if (ngrids > 1 && gridtype == GRID_GENERIC && gridInqXsize(gridID) == 0 && gridInqYsize(gridID) == 0) return gridID2;

  std::vector<long> xsize(nfiles);
  std::vector<long> ysize(nfiles);
  xyinfoType *xyinfo = (xyinfoType *) Malloc(nfiles * sizeof(xyinfoType));
  std::vector<std::vector<double>> xvals(nfiles);
  std::vector<std::vector<double>> yvals(nfiles);

  for (int fileID = 0; fileID < nfiles; fileID++)
    {
      gridID = vlistGrid(ef[fileID].vlistID, igrid);
      int gridtype = gridInqType(gridID);
      if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN)
        lregular = true;
      else if (gridtype == GRID_CURVILINEAR)
        lcurvilinear = true;
      else if (gridtype == GRID_GENERIC /*&& gridInqXsize(gridID) > 0 && gridInqYsize(gridID) > 0*/)
        ;
      else
        cdoAbort("Unsupported grid type: %s!", gridNamePtr(gridtype));

      xsize[fileID] = gridInqXsize(gridID);
      ysize[fileID] = gridInqYsize(gridID);
      if (xsize[fileID] == 0) xsize[fileID] = 1;
      if (ysize[fileID] == 0) ysize[fileID] = 1;

      if (lregular)
        {
          xvals[fileID].resize(xsize[fileID]);
          yvals[fileID].resize(ysize[fileID]);
        }
      else if (lcurvilinear)
        {
          xvals[fileID].resize(xsize[fileID] * ysize[fileID]);
          yvals[fileID].resize(xsize[fileID] * ysize[fileID]);
        }

      if (lregular || lcurvilinear)
        {
          gridInqXvals(gridID, &xvals[fileID][0]);
          gridInqYvals(gridID, &yvals[fileID][0]);
        }
      // printf("fileID %d, gridID %d\n", fileID, gridID);

      if (lregular)
        {
          xyinfo[fileID].x = xvals[fileID][0];
          xyinfo[fileID].y = yvals[fileID][0];
          xyinfo[fileID].id = fileID;

          if (ysize[fileID] > 1)
            {
              if (yvals[fileID][0] > yvals[fileID][ysize[fileID] - 1]) lsouthnorth = false;
            }
        }
      else
        {
          xyinfo[fileID].x = 0;
          xyinfo[fileID].y = 0;
          xyinfo[fileID].id = fileID;
        }
    }

  if (cdoVerbose && lregular)
    for (int fileID = 0; fileID < nfiles; fileID++) printf("1 %d %g %g \n", xyinfo[fileID].id, xyinfo[fileID].x, xyinfo[fileID].y);

  if (lregular)
    {
      qsort(xyinfo, nfiles, sizeof(xyinfoType), cmpx);

      if (cdoVerbose)
        for (int fileID = 0; fileID < nfiles; fileID++)
          printf("2 %d %g %g \n", xyinfo[fileID].id, xyinfo[fileID].x, xyinfo[fileID].y);

      if (lsouthnorth)
        qsort(xyinfo, nfiles, sizeof(xyinfoType), cmpxy_lt);
      else
        qsort(xyinfo, nfiles, sizeof(xyinfoType), cmpxy_gt);

      if (cdoVerbose)
        for (int fileID = 0; fileID < nfiles; fileID++)
          printf("3 %d %g %g \n", xyinfo[fileID].id, xyinfo[fileID].x, xyinfo[fileID].y);

      if (nx <= 0)
        {
          nx = 1;
          for (int fileID = 1; fileID < nfiles; fileID++)
            {
              if (DBL_IS_EQUAL(xyinfo[0].y, xyinfo[fileID].y))
                nx++;
              else
                break;
            }
        }
    }
  else
    {
      if (nx <= 0) nx = nfiles;
    }

  long ny = nfiles / nx;
  if (nx * ny != nfiles) cdoAbort("Number of input files (%ld) and number of blocks (%ldx%ld) differ!", nfiles, nx, ny);

  long xsize2 = 0;
  for (long i = 0; i < nx; ++i) xsize2 += xsize[xyinfo[i].id];
  long ysize2 = 0;
  for (long j = 0; j < ny; ++j) ysize2 += ysize[xyinfo[j * nx].id];
  if (cdoVerbose) cdoPrint("xsize2 %ld  ysize2 %ld", xsize2, ysize2);

  std::vector<double> xvals2, yvals2;
  if (lregular)
    {
      xvals2.resize(xsize2);
      yvals2.resize(ysize2);
    }
  else if (lcurvilinear)
    {
      xvals2.resize(xsize2 * ysize2);
      yvals2.resize(xsize2 * ysize2);
    }

  std::vector<long> xoff(nx + 1);
  std::vector<long> yoff(ny + 1);

  xoff[0] = 0;
  for (long i = 0; i < nx; ++i)
    {
      long idx = xyinfo[i].id;
      if (lregular) arrayCopy(xsize[idx], &xvals[idx][0], &xvals2[xoff[i]]);
      xoff[i + 1] = xoff[i] + xsize[idx];
    }

  yoff[0] = 0;
  for (long j = 0; j < ny; ++j)
    {
      long idx = xyinfo[j * nx].id;
      if (lregular) arrayCopy(ysize[idx], &yvals[idx][0], &yvals2[yoff[j]]);
      yoff[j + 1] = yoff[j] + ysize[idx];
    }

  if (ginit == false)
    {
      for (int fileID = 0; fileID < nfiles; fileID++)
        {
          long idx = xyinfo[fileID].id;
          long iy = fileID / nx;
          long ix = fileID - iy * nx;

          long offset = yoff[iy] * xsize2 + xoff[ix];
          /*
          printf("fileID %d %d, iy %d, ix %d, offset %d\n",
                 fileID, xyinfo[fileID].id, iy, ix, offset);
          */
          long ij = 0;
          for (long j = 0; j < ysize[idx]; ++j)
            for (long i = 0; i < xsize[idx]; ++i)
              {
                if (lcurvilinear)
                  {
                    xvals2[offset + j * xsize2 + i] = xvals[idx][ij];
                    yvals2[offset + j * xsize2 + i] = yvals[idx][ij];
                  }
                ef[idx].gridindex[ij++] = offset + j * xsize2 + i;
              }
        }
    }

  gridID2 = gridCreate(gridtype, xsize2 * ysize2);
  gridDefXsize(gridID2, xsize2);
  gridDefYsize(gridID2, ysize2);
  if (lregular || lcurvilinear)
    {
      gridDefXvals(gridID2, &xvals2[0]);
      gridDefYvals(gridID2, &yvals2[0]);
    }

  Free(xyinfo);

  gridID = vlistGrid(ef[0].vlistID, igrid);

  grid_copy_attributes(gridID, gridID2);

  if (gridtype == GRID_PROJECTION) grid_copy_mapping(gridID, gridID2);

  return gridID2;
}

void *
Collgrid(void *process)
{
  int nxblocks = -1;
  int varID, levelID;
  int nrecs0;

  cdoInitialize(process);

  int nfiles = cdoStreamCnt() - 1;
  const char *ofilename = cdoGetStreamName(nfiles).c_str();

  if (!cdoOverwriteMode && fileExists(ofilename) && !userFileOverwrite(ofilename))
    cdoAbort("Outputfile %s already exists!", ofilename);

  std::vector<ensfileType> ef(nfiles);

  for (int fileID = 0; fileID < nfiles; fileID++)
    {
      ef[fileID].streamID = cdoStreamOpenRead(cdoStreamName(fileID));
      ef[fileID].vlistID = cdoStreamInqVlist(ef[fileID].streamID);
    }

  int vlistID1 = ef[0].vlistID;
  vlistClearFlag(vlistID1);

  /* check that the contents is always the same */
  for (int fileID = 1; fileID < nfiles; fileID++) vlistCompare(vlistID1, ef[fileID].vlistID, CMP_NAME | CMP_NLEVEL);

  int nvars = vlistNvars(vlistID1);
  std::vector<bool> vars(nvars);
  for (varID = 0; varID < nvars; varID++) vars[varID] = false;
  std::vector<bool> vars1(nvars);
  for (varID = 0; varID < nvars; varID++) vars1[varID] = false;

  int nsel = operatorArgc();
  int noff = 0;

  if (nsel > 0)
    {
      int len = (int) strlen(operatorArgv()[0]);
      while (--len >= 0 && isdigit(operatorArgv()[0][len]))
        ;

      if (len == -1)
        {
          nsel--;
          noff++;
          nxblocks = parameter2int(operatorArgv()[0]);
        }
    }

  if (nsel == 0)
    {
      for (varID = 0; varID < nvars; varID++) vars1[varID] = true;
    }
  else
    {
      char **argnames = operatorArgv() + noff;

      if (cdoVerbose)
        for (int i = 0; i < nsel; i++) fprintf(stderr, "name %d = %s\n", i + 1, argnames[i]);

      std::vector<bool> selfound(nsel);
      for (int i = 0; i < nsel; i++) selfound[i] = false;

      char varname[CDI_MAX_NAME];
      for (varID = 0; varID < nvars; varID++)
        {
          vlistInqVarName(vlistID1, varID, varname);

          for (int isel = 0; isel < nsel; isel++)
            {
              if (strcmp(argnames[isel], varname) == 0)
                {
                  selfound[isel] = true;
                  vars1[varID] = true;
                }
            }
        }

      for (int isel = 0; isel < nsel; isel++)
        if (selfound[isel] == false) cdoAbort("Variable name %s not found!", argnames[isel]);
    }

  for (varID = 0; varID < nvars; varID++)
    {
      if (vars1[varID])
        {
          int zaxisID = vlistInqVarZaxis(vlistID1, varID);
          int nlevs = zaxisInqSize(zaxisID);
          for (int levID = 0; levID < nlevs; levID++) vlistDefFlag(vlistID1, varID, levID, TRUE);
        }
    }

  for (int fileID = 0; fileID < nfiles; fileID++)
    {
      size_t gridsize = vlistGridsizeMax(ef[fileID].vlistID);
      ef[fileID].gridsize = gridsize;
      ef[fileID].gridindex.resize(gridsize);
      ef[fileID].array.resize(gridsize);
    }

  int vlistID2 = vlistCreate();
  cdoVlistCopyFlag(vlistID2, vlistID1);
  /*
  if ( cdoVerbose )
    {
      vlistPrint(vlistID1);
      vlistPrint(vlistID2);
    }
  */
  // int vlistID2 = vlistDuplicate(vlistID1);
  int nvars2 = vlistNvars(vlistID2);

  int ngrids1 = vlistNgrids(vlistID1);
  int ngrids2 = vlistNgrids(vlistID2);

  std::vector<int> gridIDs(ngrids2);

  bool ginit = false;
  for (int i2 = 0; i2 < ngrids2; ++i2)
    {
      int i1;
      for (i1 = 0; i1 < ngrids1; ++i1)
        if (vlistGrid(vlistID1, i1) == vlistGrid(vlistID2, i2)) break;

      //   printf("i1 %d i2 %d\n", i1, i2);

      if (!ginit)
        {
          gridIDs[i2] = genGrid(ngrids2, nfiles, ef, ginit, i1, nxblocks);
          if (gridIDs[i2] != -1) ginit = true;
        }
      else
        gridIDs[i2] = genGrid(ngrids2, nfiles, ef, ginit, i1, nxblocks);
    }

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  size_t gridsize2 = 0;
  for (int i = 0; i < ngrids2; ++i)
    {
      if (gridIDs[i] != -1)
        {
          if (gridsize2 == 0) gridsize2 = gridInqSize(gridIDs[i]);
          if (gridsize2 != gridInqSize(gridIDs[i])) cdoAbort("gridsize differ!");
          vlistChangeGridIndex(vlistID2, i, gridIDs[i]);
        }
    }

  for (varID = 0; varID < nvars2; varID++)
    {
      int gridID = vlistInqVarGrid(vlistID2, varID);

      for (int i = 0; i < ngrids2; ++i)
        {
          if (gridIDs[i] != -1)
            {
              if (gridID == vlistGrid(vlistID2, i)) vars[varID] = true;
              break;
            }
        }
    }

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(nfiles), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  std::vector<double> array2;
  if (gridsize2 > 0) array2.resize(gridsize2);

  int tsID = 0;
  do
    {
      nrecs0 = cdoStreamInqTimestep(ef[0].streamID, tsID);
      for (int fileID = 1; fileID < nfiles; fileID++)
        {
          int nrecs = cdoStreamInqTimestep(ef[fileID].streamID, tsID);
          if (nrecs != nrecs0)
            cdoAbort("Number of records at time step %d of %s and %s differ!", tsID + 1, cdoGetStreamName(0).c_str(),
                     cdoGetStreamName(fileID).c_str());
        }

      taxisCopyTimestep(taxisID2, taxisID1);

      if (nrecs0 > 0) pstreamDefTimestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs0; recID++)
        {
          pstreamInqRecord(ef[0].streamID, &varID, &levelID);
          if (cdoVerbose && tsID == 0) printf(" tsID, recID, varID, levelID %d %d %d %d\n", tsID, recID, varID, levelID);

          for (int fileID = 1; fileID < nfiles; fileID++)
            {
              int varIDx, levelIDx;
              pstreamInqRecord(ef[fileID].streamID, &varIDx, &levelIDx);
            }

          if (vlistInqFlag(vlistID1, varID, levelID) == TRUE)
            {
              int varID2 = vlistFindVar(vlistID2, varID);
              int levelID2 = vlistFindLevel(vlistID2, varID, levelID);
              if (cdoVerbose && tsID == 0) printf("varID %d %d levelID %d %d\n", varID, varID2, levelID, levelID2);

              double missval = vlistInqVarMissval(vlistID2, varID2);
              for (size_t i = 0; i < gridsize2; i++) array2[i] = missval;

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
              for (int fileID = 0; fileID < nfiles; fileID++)
                {
                  pstreamReadRecord(ef[fileID].streamID, &ef[fileID].array[0], &ef[fileID].nmiss);

                  if (vars[varID2])
                    {
                      size_t gridsize = ef[fileID].gridsize;
                      for (size_t i = 0; i < gridsize; ++i) array2[ef[fileID].gridindex[i]] = ef[fileID].array[i];
                    }
                }

              pstreamDefRecord(streamID2, varID2, levelID2);

              if (vars[varID2])
                {
                  size_t nmiss = arrayNumMV(gridsize2, &array2[0], missval);
                  pstreamWriteRecord(streamID2, &array2[0], nmiss);
                }
              else
                pstreamWriteRecord(streamID2, &ef[0].array[0], ef[0].nmiss);
            }
        }

      tsID++;
    }
  while (nrecs0 > 0);

  for (int fileID = 0; fileID < nfiles; fileID++) pstreamClose(ef[fileID].streamID);

  pstreamClose(streamID2);

  cdoFinish();

  return 0;
}

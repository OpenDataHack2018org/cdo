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

#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <limits.h>

#include <cdi.h>
#include "cdo_int.h"
#include "grid.h"
#include "griddes.h"
#include "error.h"
#include "cdoDebugOutput.h"
#include "util_wildcards.h"

int grid_read(FILE *gfp, const char *dname);
int grid_read_pingo(FILE *gfp, const char *dname);

void
GridDesciption::init()
{
  this->type = CDI_UNDEFID;
  this->datatype = CDI_UNDEFID;
  this->size = 0;
  this->xsize = 0;
  this->ysize = 0;
  this->np = 0;
  this->lcomplex = 1;
  this->ntr = 0;
  this->nvertex = 0;
  this->genBounds = false;
  this->def_xfirst = false;
  this->def_yfirst = false;
  this->def_xlast = false;
  this->def_ylast = false;
  this->def_xinc = false;
  this->def_yinc = false;
  this->xfirst = 0;
  this->yfirst = 0;
  this->xlast = 0;
  this->ylast = 0;
  this->xinc = 0;
  this->yinc = 0;
  this->nd = 0;
  this->ni = 0;
  this->ni2 = 0;
  this->ni3 = 0;
  this->number = 0;
  this->position = 0;
  this->uuid[0] = 0;
  this->path[0] = 0;
  this->xname[0] = 0;
  this->xlongname[0] = 0;
  this->xunits[0] = 0;
  this->yname[0] = 0;
  this->ylongname[0] = 0;
  this->yunits[0] = 0;
  this->xdimname[0] = 0;
  this->ydimname[0] = 0;
  this->vdimname[0] = 0;
  this->uvRelativeToGrid = false;
  this->scanningMode = 64;
  /*
    scanningMode  = 128 * iScansNegatively + 64 * jScansPositively + 32 * jPointsAreConsecutive;
              64  = 128 * 0                + 64 *        1         + 32 * 0 
              00  = 128 * 0                + 64 *        0         + 32 * 0
              96  = 128 * 0                + 64 *        1         + 32 * 1
    Default  implicit scanning mode is 64: i and j scan positively, i points are consecutive (row-major)
  */
}

int
getoptname(char *optname, const char *optstring, int nopt)
{
  int nerr = 0;
  const char *pname = optstring;
  const char *pend = optstring;

  for (int i = 0; i < nopt; i++)
    {
      pend = strchr(pname, ',');
      if (pend == NULL)
        break;
      else
        pname = pend + 1;
    }

  if (pend)
    {
      pend = strchr(pname, ',');
      size_t namelen = (pend == NULL) ? strlen(pname) : pend - pname;
      memcpy(optname, pname, namelen);
      optname[namelen] = '\0';
    }
  else
    nerr = 1;

  return nerr;
}

int
gridDefine(GridDesciption &grid)
{
  int gridID = CDI_UNDEFID;

  switch (grid.type)
    {
    case GRID_GENERIC:
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_GAUSSIAN_REDUCED:
    case GRID_PROJECTION:
      {
        if (grid.size != 1)
          {
            if (grid.xsize == 0 && grid.type != GRID_GAUSSIAN_REDUCED) Error("xsize undefined!");
            if (grid.ysize == 0) Error("ysize undefined!");
          }

        if (grid.size == 0) grid.size = grid.xsize * grid.ysize;

        if (grid.xsize && grid.size != grid.xsize * grid.ysize)
          Error("Inconsistent grid declaration: xsize*ysize!=gridsize (xsize=%zu ysize=%zu gridsize=%zu)", grid.xsize, grid.ysize,
                grid.size);

        // if ( grid.size < 0 || grid.size > INT_MAX ) Error("grid size (%ld)
        // out of bounds (0 - %d)!", grid.size, INT_MAX);

        gridID = gridCreate(grid.type, grid.size);

        if (grid.xsize > 0) gridDefXsize(gridID, grid.xsize);
        if (grid.ysize > 0) gridDefYsize(gridID, grid.ysize);
        if (grid.np > 0) gridDefNP(gridID, grid.np);

        if (grid.uvRelativeToGrid) gridDefUvRelativeToGrid(gridID, 1);
        if (grid.nvertex) gridDefNvertex(gridID, grid.nvertex);

        if ((grid.def_xfirst || grid.def_xlast || grid.def_xinc) && grid.xvals.size() == 0)
          {
            grid.xvals.resize(grid.xsize);
            gridGenXvals(grid.xsize, grid.xfirst, grid.xlast, grid.xinc, grid.xvals.data());

            if (grid.genBounds && grid.xbounds.size() == 0 && grid.xsize > 1)
              {
                grid.nvertex = 2;
                grid.xbounds.resize(grid.xsize * grid.nvertex);
                for (size_t i = 0; i < grid.xsize - 1; ++i)
                  {
                    grid.xbounds[2 * i + 1] = 0.5 * (grid.xvals[i] + grid.xvals[i + 1]);
                    grid.xbounds[2 * (i + 1)] = 0.5 * (grid.xvals[i] + grid.xvals[i + 1]);
                  }
                grid.xbounds[0] = 2 * grid.xvals[0] - grid.xbounds[1];
                grid.xbounds[2 * grid.xsize - 1] = 2 * grid.xvals[grid.xsize - 1] - grid.xbounds[2 * (grid.xsize - 1)];
              }
          }

        if ((grid.def_yfirst || grid.def_ylast || grid.def_yinc) && grid.yvals.size() == 0)
          {
            if (!grid.def_ylast) grid.ylast = grid.yfirst;
            grid.yvals.resize(grid.ysize);
            gridGenYvals(grid.type, grid.ysize, grid.yfirst, grid.ylast, grid.yinc, grid.yvals.data());

            if (grid.genBounds && grid.ybounds.size() == 0 && grid.ysize > 1)
              {
                grid.nvertex = 2;
                grid.ybounds.resize(grid.ysize * grid.nvertex);
                for (size_t i = 0; i < grid.ysize - 1; ++i)
                  {
                    grid.ybounds[2 * i + 1] = 0.5 * (grid.yvals[i] + grid.yvals[i + 1]);
                    grid.ybounds[2 * (i + 1)] = 0.5 * (grid.yvals[i] + grid.yvals[i + 1]);
                  }

                if (grid.yvals[0] > grid.yvals[grid.ysize - 1])
                  {
                    grid.ybounds[0] = 90;
                    grid.ybounds[grid.ysize * grid.nvertex - 1] = -90;
                  }
                else
                  {
                    grid.ybounds[0] = -90;
                    grid.ybounds[grid.ysize * grid.nvertex - 1] = 90;
                  }
              }
          }

        if (grid.xvals.size()) gridDefXvals(gridID, grid.xvals.data());
        if (grid.yvals.size()) gridDefYvals(gridID, grid.yvals.data());
        if (grid.xbounds.size()) gridDefXbounds(gridID, grid.xbounds.data());
        if (grid.ybounds.size()) gridDefYbounds(gridID, grid.ybounds.data());
        if (grid.mask.size()) gridDefMask(gridID, grid.mask.data());
        if (grid.rowlon.size()) gridDefRowlon(gridID, grid.ysize, grid.rowlon.data());

        break;
      }
    case GRID_CURVILINEAR:
    case GRID_UNSTRUCTURED:
      {
        if (grid.size == 0) grid.size = (grid.type == GRID_CURVILINEAR) ? grid.xsize * grid.ysize : grid.xsize;

        gridID = gridCreate(grid.type, grid.size);

        if (grid.type == GRID_CURVILINEAR)
          {
            if (grid.xsize == 0) Error("xsize undefined!");
            if (grid.ysize == 0) Error("ysize undefined!");
            gridDefXsize(gridID, grid.xsize);
            gridDefYsize(gridID, grid.ysize);
          }
        else
          {
            if (grid.nvertex > 0) gridDefNvertex(gridID, grid.nvertex);
            if (grid.number > 0)
              {
                gridDefNumber(gridID, grid.number);
                if (grid.position >= 0) gridDefPosition(gridID, grid.position);
              }
            if (*grid.path) gridDefReference(gridID, grid.path);
          }

        if (grid.xvals.size()) gridDefXvals(gridID, grid.xvals.data());
        if (grid.yvals.size()) gridDefYvals(gridID, grid.yvals.data());
        if (grid.area.size()) gridDefArea(gridID, grid.area.data());
        if (grid.xbounds.size()) gridDefXbounds(gridID, grid.xbounds.data());
        if (grid.ybounds.size()) gridDefYbounds(gridID, grid.ybounds.data());
        if (grid.mask.size()) gridDefMask(gridID, grid.mask.data());

        break;
      }
    case GRID_SPECTRAL:
      {
        if (grid.ntr == 0) Error("truncation undefined!");
        if (grid.size == 0) grid.size = (grid.ntr + 1) * (grid.ntr + 2);

        gridID = gridCreate(grid.type, grid.size);

        gridDefTrunc(gridID, grid.ntr);
        gridDefComplexPacking(gridID, grid.lcomplex);

        break;
      }
    case GRID_GME:
      {
        if (grid.nd == 0) Error("nd undefined!");
        if (grid.ni == 0) Error("ni undefined!");
        if (grid.size == 0) Error("size undefined!");

        gridID = gridCreate(grid.type, grid.size);

        gridDefParamGME(gridID, grid.nd, grid.ni, grid.ni2, grid.ni3);

        if (grid.mask.size()) gridDefMask(gridID, grid.mask.data());

        break;
      }
    default:
      {
        if (grid.type == -1)
          Error("Undefined grid type!");
        else
          Error("Unsupported grid type: %s", gridNamePtr(grid.type));

        break;
      }
    }

  if (grid.datatype != CDI_UNDEFID) gridDefDatatype(gridID, grid.datatype);

  if (grid.uuid[0]) gridDefUUID(gridID, grid.uuid);

  if (grid.xname[0]) cdiGridDefKeyStr(gridID, CDI_KEY_XNAME, strlen(grid.xname) + 1, grid.xname);
  if (grid.xlongname[0]) cdiGridDefKeyStr(gridID, CDI_KEY_XLONGNAME, strlen(grid.xlongname) + 1, grid.xlongname);
  if (grid.xunits[0]) cdiGridDefKeyStr(gridID, CDI_KEY_XUNITS, strlen(grid.xunits) + 1, grid.xunits);
  if (grid.yname[0]) cdiGridDefKeyStr(gridID, CDI_KEY_YNAME, strlen(grid.yname) + 1, grid.yname);
  if (grid.ylongname[0]) cdiGridDefKeyStr(gridID, CDI_KEY_YLONGNAME, strlen(grid.ylongname) + 1, grid.ylongname);
  if (grid.yunits[0]) cdiGridDefKeyStr(gridID, CDI_KEY_YUNITS, strlen(grid.yunits) + 1, grid.yunits);
  if (grid.xdimname[0]) cdiGridDefKeyStr(gridID, CDI_KEY_XDIMNAME, strlen(grid.xdimname) + 1, grid.xdimname);
  if (grid.ydimname[0]) cdiGridDefKeyStr(gridID, CDI_KEY_YDIMNAME, strlen(grid.ydimname) + 1, grid.ydimname);
  if (grid.vdimname[0]) cdiGridDefKeyStr(gridID, CDI_KEY_VDIMNAME, strlen(grid.vdimname) + 1, grid.vdimname);

  return gridID;
}

int
cdoDefineGrid(const char *pgridfile)
{
  int gridID = -1;
  bool isreg = false;
  bool lalloc = false;
  char *gridfile = strdup(pgridfile);
  size_t len = strlen(gridfile);
  int gridno = 1;
  if ( len > 2 && gridfile[len-2] == ':' && isdigit(gridfile[len-1]) )
    {
      gridno = gridfile[len-1] - '0';
      gridfile[len-2] = 0;
    }

  char *filename = expand_filename(gridfile);
  if (filename)
    lalloc = true;
  else
    filename = (char *) gridfile;

  int fileno = open(filename, O_RDONLY);
  if (fileno >= 0)
    {
      struct stat filestat;
      if (fstat(fileno, &filestat) == 0) isreg = S_ISREG(filestat.st_mode);
    }

  if (fileno == -1 || !isreg)
    {
      if (isreg) close(fileno);

      gridID = grid_from_name(gridfile);

      if (gridID == -1) cdoAbort("Open failed on %s!", gridfile);
    }
  else
    {
      char buffer[4];
      if (read(fileno, buffer, 4) != 4) SysError("Read grid from %s failed!", filename);

      close(fileno);

      if (cmpstrlen(buffer, "CDF", len) == 0)
        {
          if (CdoDebug::cdoDebug) cdoPrint("Grid from NetCDF file");
          gridID = gridFromNCfile(filename);
        }

      if (gridID == -1)
        {
          if (cmpstrlen(buffer + 1, "HDF", len) == 0)
            {
              if (CdoDebug::cdoDebug) cdoPrint("Grid from HDF5 file");
              gridID = gridFromH5file(filename);
            }
        }

      if (gridID == -1)
        {
          if (cmpstrlen(buffer + 1, "HDF", len) == 0)
            {
              if (CdoDebug::cdoDebug) cdoPrint("Grid from NetCDF4 file");
              gridID = gridFromNCfile(filename);
            }
        }

      if (gridID == -1)
        {
          if (CdoDebug::cdoDebug) cdoPrint("Grid from CDI file");
          openLock();
          int streamID = streamOpenRead(filename);
          openUnlock();
          if (streamID >= 0)
            {
              int vlistID = streamInqVlist(streamID);
              int ngrids = vlistNgrids(vlistID);
              if ( gridno < 1 || gridno > ngrids ) cdoAbort("Grid number %d not available in %s!", gridno, filename);
              gridID = vlistGrid(vlistID, gridno-1);
              streamClose(streamID);
            }
        }

      if (gridID == -1)
        {
          if (CdoDebug::cdoDebug) cdoPrint("grid from ASCII file");
          FILE *gfp = fopen(filename, "r");
          // size_t buffersize = 20*1024*1024;
          // char *buffer = (char*) Malloc(buffersize);
          // setvbuf(gfp, buffer, _IOFBF, buffersize);
          gridID = grid_read(gfp, filename);
          fclose(gfp);
          // free(buffer);
        }

      if (gridID == -1)
        {
          if (CdoDebug::cdoDebug) cdoPrint("grid from PINGO file");
          FILE *gfp = fopen(filename, "r");
          gridID = grid_read_pingo(gfp, filename);
          fclose(gfp);
        }

      if (gridID == -1) cdoAbort("Invalid grid description file %s!", filename);
    }

  if (lalloc) Free(filename);
  free(gridfile);

  return gridID;
}

void
cdo_set_grids(const char *gridarg)
{
  char gridfile[4096];
  int nfile = 0;

  while (getoptname(gridfile, gridarg, nfile++) == 0)
    {
      (void) cdoDefineGrid(gridfile);
    }
}

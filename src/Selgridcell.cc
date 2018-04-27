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

      Selindex     selindex    Select grid indices
*/

#include <cdi.h>

#include "cdo_int.h"
#include "grid.h"
#include "listarray.h"
#include "pstream_int.h"
#include "util_files.h"

int gengridcell(int gridID1, size_t gridsize2, long *cellidx);

static int
genindexgrid(int gridID1, size_t gridsize2, long *cellidx)
{
  int gridID0 = gridID1;
  int gridtype1 = gridInqType(gridID1);

  if (gridtype1 == GRID_LONLAT || gridtype1 == GRID_GAUSSIAN || gridtype1 == GRID_PROJECTION)
    {
      gridID1 = gridToCurvilinear(gridID1, 0);
      gridtype1 = GRID_CURVILINEAR;
    }

  int gridID2 = -1;
  if (gridtype1 == GRID_UNSTRUCTURED || gridtype1 == GRID_CURVILINEAR) gridID2 = gengridcell(gridID1, gridsize2, cellidx);

  if (gridID0 != gridID1) gridDestroy(gridID1);

  return gridID2;
}

static void
sel_index(double *array1, double *array2, long nind, long *indarr)
{
  for (long i = 0; i < nind; ++i)
    {
      array2[i] = array1[indarr[i]];
    }
}

void *
Selgridcell(void *process)
{
  int nrecs;
  int varID;
  int gridID1 = -1, gridID2;
  int index, gridtype = -1;
  typedef struct
  {
    int gridID1, gridID2;
  } sindex_t;
  ListArray<int> listArrayInt;

  cdoInitialize(process);

  cdoOperatorAdd("selgridcell", 0, 0, "grid cell indices (1-N)");
  int DELGRIDCELL = cdoOperatorAdd("delgridcell", 0, 0, "grid cell indices (1-N)");

  operatorInputArg(cdoOperatorEnter(0));

  int operatorID = cdoOperatorID();

  int nind;
  int *indarr = NULL;
  if (operatorArgc() == 1 && fileExists(operatorArgv()[0]))
    {
      bool *cdo_read_mask(const char *maskfile, size_t *n);
      size_t n = 0;
      bool *mask = cdo_read_mask(operatorArgv()[0], &n);
      nind = 0;
      for (size_t i = 0; i < n; ++i)
        if (mask[i]) nind++;
      if (nind == 0)
        cdoAbort("Mask is empty!");
      else
        {
          indarr = (int *) Malloc(nind * sizeof(int));
          nind = 0;
          for (size_t i = 0; i < n; ++i)
            if (mask[i]) indarr[nind++] = i;
        }
      if (mask) Free(mask);
    }
  else
    {
      nind = listArrayInt.argvToInt(operatorArgc(), operatorArgv());
      indarr = listArrayInt.data();

      if (cdoVerbose)
        for (int i = 0; i < nind; i++) cdoPrint("int %d = %d", i + 1, indarr[i]);

      for (int i = 0; i < nind; i++) indarr[i] -= 1;
    }

  int indmin = indarr[0];
  int indmax = indarr[0];
  for (int i = 1; i < nind; i++)
    {
      if (indmax < indarr[i]) indmax = indarr[i];
      if (indmin > indarr[i]) indmin = indarr[i];
    }

  if (indmin < 0) cdoAbort("Index < 1 not allowed!");

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int nvars = vlistNvars(vlistID1);
  std::vector<bool> vars(nvars);
  for (varID = 0; varID < nvars; varID++) vars[varID] = false;

  int ngrids = vlistNgrids(vlistID1);
  std::vector<sindex_t> sindex(ngrids);

  long ncells = nind;
  std::vector<long> cellidx;
  if (operatorID == DELGRIDCELL)
    {
      size_t gridsize = vlistGridsizeMax(vlistID1);
      ncells = gridsize - nind;
      cellidx.resize(gridsize);
      for (size_t i = 0; i < gridsize; ++i) cellidx[i] = 1;
      for (long i = 0; i < nind; ++i) cellidx[indarr[i]] = 0;
      long j = 0;
      for (size_t i = 0; i < gridsize; ++i)
        if (cellidx[i] == 1) cellidx[j++] = i;
      if (j != ncells) cdoAbort("Internal error; number of cells differ");
    }
  else
    {
      cellidx.resize(nind);
      for (int i = 0; i < nind; ++i) cellidx[i] = indarr[i];
    }

  if (ncells == 0) cdoAbort("Mask is empty!");

  for (index = 0; index < ngrids; index++)
    {
      gridID1 = vlistGrid(vlistID1, index);
      gridtype = gridInqType(gridID1);

      size_t gridsize = gridInqSize(gridID1);
      if (gridsize == 1) continue;
      if (indmax >= (int) gridsize)
        {
          cdoWarning("Max grid index is greater than grid size, skipped grid %d!", index + 1);
          continue;
        }

      gridID2 = genindexgrid(gridID1, ncells, cellidx.data());

      if (gridID2 == -1)
        {
          cdoWarning("Unsupported grid type >%s<, skipped grid %d!", gridNamePtr(gridtype), index + 1);
          continue;
        }

      sindex[index].gridID1 = gridID1;
      sindex[index].gridID2 = gridID2;

      vlistChangeGridIndex(vlistID2, index, gridID2);

      for (varID = 0; varID < nvars; varID++)
        if (gridID1 == vlistInqVarGrid(vlistID1, varID)) vars[varID] = true;
    }

  for (varID = 0; varID < nvars; varID++)
    if (vars[varID]) break;

  if (varID >= nvars) cdoAbort("No variables selected!");

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  size_t gridsize = vlistGridsizeMax(vlistID1);
  if (vlistNumber(vlistID1) != CDI_REAL) gridsize *= 2;
  std::vector<double> array1(gridsize);

  size_t gridsize2 = vlistGridsizeMax(vlistID2);
  if (vlistNumber(vlistID2) != CDI_REAL) gridsize2 *= 2;
  std::vector<double> array2(gridsize2);

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          size_t nmiss;
          int varID, levelID;
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamReadRecord(streamID1, array1.data(), &nmiss);

          pstreamDefRecord(streamID2, varID, levelID);

          if (vars[varID])
            {
              gridID1 = vlistInqVarGrid(vlistID1, varID);

              for (index = 0; index < ngrids; index++)
                if (gridID1 == sindex[index].gridID1) break;

              if (index == ngrids) cdoAbort("Internal problem, grid not found!");

              gridsize2 = gridInqSize(sindex[index].gridID2);

              sel_index(array1.data(), array2.data(), ncells, cellidx.data());

              if (nmiss)
                {
                  double missval = vlistInqVarMissval(vlistID2, varID);
                  nmiss = arrayNumMV(gridsize2, array2.data(), missval);
                }

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

  vlistDestroy(vlistID2);

  cdoFinish();

  return 0;
}

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

      Setbox     setclonlatbox   Set lon/lat box to constant
      Setbox     setcindexbox    Set index box to constant
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "grid.h"

void genlonlatbox(int argc_offset, int gridID1, long *lat1, long *lat2, long *lon11, long *lon12, long *lon21, long *lon22);

void genindexbox(int argc_offset, int gridID1, long *lat1, long *lat2, long *lon11, long *lon12, long *lon21, long *lon22);

static void
setcbox(double constant, double *array, int gridID, long lat1, long lat2, long lon11, long lon12, long lon21, long lon22)
{
  long nlon = gridInqXsize(gridID);
  long nlat = gridInqYsize(gridID);

  for (long ilat = 0; ilat < nlat; ilat++)
    for (long ilon = 0; ilon < nlon; ilon++)
      if ((lat1 <= ilat && ilat <= lat2 && ((lon11 <= ilon && ilon <= lon12) || (lon21 <= ilon && ilon <= lon22))))
        {
          array[nlon * ilat + ilon] = constant;
        }
}

void *
Setbox(void *process)
{
  int nrecs;
  int varID, levelID;
  int gridID = -1;
  int index, gridtype;
  size_t nmiss;
  long lat1, lat2, lon11, lon12, lon21, lon22;

  cdoInitialize(process);

  // clang-format off
  int SETCLONLATBOX = cdoOperatorAdd("setclonlatbox", 0, 0, "constant, western and eastern longitude and southern and northern latitude");
  int SETCINDEXBOX = cdoOperatorAdd("setcindexbox", 0, 0, "constant, index of first and last longitude and index of first and last latitude");
  // clang-format on

  int operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  double constant = parameter2double(operatorArgv()[0]);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);

  int ngrids = vlistNgrids(vlistID1);
  int ndiffgrids = 0;
  for (index = 1; index < ngrids; index++)
    if (vlistGrid(vlistID1, 0) != vlistGrid(vlistID1, index)) ndiffgrids++;

  for (index = 0; index < ngrids; index++)
    {
      gridID = vlistGrid(vlistID1, index);
      gridtype = gridInqType(gridID);
      if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN) break;
      if (gridtype == GRID_CURVILINEAR) break;
      if (operatorID == SETCINDEXBOX && gridtype == GRID_GENERIC && gridInqXsize(gridID) > 0 && gridInqYsize(gridID) > 0) break;
    }

  if (gridInqType(gridID) == GRID_GAUSSIAN_REDUCED)
    cdoAbort("Gaussian reduced grid found. Use option -R to convert it to a regular grid!");

  if (index == ngrids) cdoAbort("No regular grid found!");
  if (ndiffgrids > 0) cdoAbort("Too many different grids!");

  operatorInputArg(cdoOperatorEnter(operatorID));

  if (operatorID == SETCLONLATBOX)
    genlonlatbox(1, gridID, &lat1, &lat2, &lon11, &lon12, &lon21, &lon22);
  else
    genindexbox(1, gridID, &lat1, &lat2, &lon11, &lon12, &lon21, &lon22);

  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int nvars = vlistNvars(vlistID1);
  std::vector<bool> vars(nvars);
  for (varID = 0; varID < nvars; varID++) vars[varID] = (gridID == vlistInqVarGrid(vlistID1, varID));

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());

  pstreamDefVlist(streamID2, vlistID2);

  size_t gridsize = gridInqSize(gridID);
  std::vector<double> array(gridsize);

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);

          if (vars[varID])
            {
              pstreamReadRecord(streamID1, array.data(), &nmiss);

              setcbox(constant, array.data(), gridID, lat1, lat2, lon11, lon12, lon21, lon22);

              double missval = vlistInqVarMissval(vlistID1, varID);
              nmiss = arrayNumMV(gridsize, array.data(), missval);
              pstreamDefRecord(streamID2, varID, levelID);
              pstreamWriteRecord(streamID2, array.data(), nmiss);
            }
        }

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

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

      Setgrid    setgrid         Set grid
      Setgrid    setgridtype     Set grid type
      Setgrid    setgridarea     Set grid area
      Setgrid    setgridmask     Set grid mask
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "grid.h"

void *
Setgrid(void *process)
{
  int nrecs;
  int varID, levelID;
  int gridID2 = -1;
  int gridtype = -1;
  size_t nmiss;
  size_t areasize = 0;
  size_t masksize = 0;
  bool lregular = false;
  bool lregularnn = false;
  bool ldereference = false;
  int number = 0, position = 0;
  int grid2_nvgp;
  int lbounds = TRUE;
  char *gridname = NULL;
  char *griduri = NULL;
  std::vector<int> grid2_vgpm;
  std::vector<double> gridmask, areaweight;

  cdoInitialize(process);

  // clang-format off
  int SETGRID       = cdoOperatorAdd("setgrid",       0, 0, "grid description file or name");
  int SETGRIDTYPE   = cdoOperatorAdd("setgridtype",   0, 0, "grid type");
  int SETGRIDAREA   = cdoOperatorAdd("setgridarea",   0, 0, "filename with area weights");
  int SETGRIDMASK   = cdoOperatorAdd("setgridmask",   0, 0, "filename with grid mask");
  int UNSETGRIDMASK = cdoOperatorAdd("unsetgridmask", 0, 0, NULL);
  int SETGRIDNUMBER = cdoOperatorAdd("setgridnumber", 0, 0, "grid number and optionally grid position");
  int SETGRIDURI    = cdoOperatorAdd("setgriduri",    0, 0, "reference URI of the horizontal grid");
  int USEGRIDNUMBER = cdoOperatorAdd("usegridnumber", 0, 0, "use existing grid identified by grid number");
  // clang-format on

  int operatorID = cdoOperatorID();

  if (operatorID != UNSETGRIDMASK) operatorInputArg(cdoOperatorEnter(operatorID));

  if (operatorID == SETGRID)
    {
      operatorCheckArgc(1);
      gridID2 = cdoDefineGrid(operatorArgv()[0]);
    }
  else if (operatorID == SETGRIDTYPE)
    {
      operatorCheckArgc(1);
      gridname = operatorArgv()[0];

      // clang-format off
      if      ( strcmp(gridname, "curvilinear0") == 0 )  {gridtype = GRID_CURVILINEAR; lbounds = 0;}
      else if ( strcmp(gridname, "curvilinear") == 0 )   {gridtype = GRID_CURVILINEAR; lbounds = 1;}
      else if ( strcmp(gridname, "cell") == 0 )           gridtype = GRID_UNSTRUCTURED;
      else if ( strcmp(gridname, "unstructured0") == 0 ) {gridtype = GRID_UNSTRUCTURED; lbounds = 0;}
      else if ( strcmp(gridname, "unstructured") == 0 )  {gridtype = GRID_UNSTRUCTURED; lbounds = 1;}
      else if ( strcmp(gridname, "generic") == 0 )        gridtype = GRID_GENERIC;
      else if ( strcmp(gridname, "dereference") == 0 )    ldereference = true;
      else if ( strcmp(gridname, "lonlat") == 0 )         gridtype = GRID_LONLAT;
      else if ( strcmp(gridname, "gaussian") == 0 )       gridtype = GRID_GAUSSIAN;
      else if ( strcmp(gridname, "regularnn") == 0 )     {gridtype = GRID_GAUSSIAN; lregularnn = true;}
      else if ( strcmp(gridname, "regular") == 0 )       {gridtype = GRID_GAUSSIAN; lregular = true;}
      else cdoAbort("Unsupported grid name: %s", gridname);
      // clang-format on
    }
  else if (operatorID == SETGRIDAREA)
    {
      operatorCheckArgc(1);
      char *areafile = operatorArgv()[0];

      int streamID = streamOpenRead(areafile);
      int vlistID = streamInqVlist(streamID);

      nrecs = streamInqTimestep(streamID, 0);
      streamInqRecord(streamID, &varID, &levelID);

      int gridID = vlistInqVarGrid(vlistID, varID);
      areasize = gridInqSize(gridID);
      areaweight.resize(areasize);

      streamReadRecord(streamID, areaweight.data(), &nmiss);
      streamClose(streamID);

      if (cdoVerbose)
        {
          double arrmean, arrmin, arrmax;
          arrayMinMaxMean(areasize, areaweight.data(), &arrmin, &arrmax, &arrmean);
          cdoPrint("areaweights: %zu %#12.5g%#12.5g%#12.5g", areasize, arrmin, arrmean, arrmax);
        }
    }
  else if (operatorID == SETGRIDMASK)
    {
      operatorCheckArgc(1);
      char *maskfile = operatorArgv()[0];
      int streamID = streamOpenRead(maskfile);

      int vlistID = streamInqVlist(streamID);

      nrecs = streamInqTimestep(streamID, 0);
      streamInqRecord(streamID, &varID, &levelID);

      double missval = vlistInqVarMissval(vlistID, varID);
      int gridID = vlistInqVarGrid(vlistID, varID);
      masksize = gridInqSize(gridID);
      gridmask.resize(masksize);

      streamReadRecord(streamID, gridmask.data(), &nmiss);
      streamClose(streamID);

      for (size_t i = 0; i < masksize; i++)
        if (DBL_IS_EQUAL(gridmask[i], missval)) gridmask[i] = 0;
    }
  else if (operatorID == USEGRIDNUMBER)
    {
      operatorCheckArgc(1);
      number = parameter2int(operatorArgv()[0]);
    }
  else if (operatorID == SETGRIDNUMBER)
    {
      if (operatorArgc() >= 1 && operatorArgc() <= 2)
        {
          number = parameter2int(operatorArgv()[0]);
          if (operatorArgc() == 2) position = parameter2int(operatorArgv()[1]);
        }
      else
        {
          operatorCheckArgc(1);
        }
    }
  else if (operatorID == SETGRIDURI)
    {
      operatorCheckArgc(1);
      griduri = operatorArgv()[0];
    }

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if (operatorID == SETGRID)
    {
      int found = 0;
      int ngrids = vlistNgrids(vlistID1);
      for (int index = 0; index < ngrids; index++)
        {
          int gridID1 = vlistGrid(vlistID1, index);

          if (gridInqSize(gridID1) == gridInqSize(gridID2))
            {
              vlistChangeGridIndex(vlistID2, index, gridID2);
              found++;
            }
        }
      if (!found) cdoWarning("No grid with %zu points found!", gridInqSize(gridID2));
    }
  else if (operatorID == SETGRIDNUMBER || operatorID == SETGRIDURI || operatorID == USEGRIDNUMBER)
    {
      if (operatorID == SETGRIDNUMBER)
        {
          int gridID1 = vlistGrid(vlistID1, 0);
          gridID2 = gridCreate(GRID_UNSTRUCTURED, gridInqSize(gridID1));
          gridDefNumber(gridID2, number);
          gridDefPosition(gridID2, position);
        }
      else if (operatorID == USEGRIDNUMBER)
        {
          if (number < 1 || number > vlistNgrids(vlistID1))
            cdoAbort("Invalid grid number: %d (max = %d)!", number, vlistNgrids(vlistID1));

          gridID2 = vlistGrid(vlistID1, number - 1);
        }
      else
        {
          int gridID1 = vlistGrid(vlistID1, 0);
          gridID2 = gridDuplicate(gridID1);
          gridDefReference(gridID2, griduri);
        }

      int found = 0;
      int ngrids = vlistNgrids(vlistID1);
      for (int index = 0; index < ngrids; index++)
        {
          int gridID1 = vlistGrid(vlistID1, index);

          if (gridInqSize(gridID1) == gridInqSize(gridID2))
            {
              vlistChangeGridIndex(vlistID2, index, gridID2);
              found++;
            }
        }
      if (!found) cdoWarning("No horizontal grid with %zu cells found!", gridInqSize(gridID2));
    }
  else if (operatorID == SETGRIDTYPE)
    {
      bool lrgrid = false;
      int ngrids = vlistNgrids(vlistID1);
      for (int index = 0; index < ngrids; index++)
        {
          int gridID1 = vlistGrid(vlistID1, index);
          int gridtype1 = gridInqType(gridID1);
          gridID2 = -1;

          if (gridtype1 == GRID_GENERIC && gridInqSize(gridID1) == 1) continue;

          if (lregular || lregularnn)
            {
              if (gridtype1 == GRID_GAUSSIAN_REDUCED) gridID2 = gridToRegular(gridID1);
            }
          else if (ldereference)
            {
              gridID2 = referenceToGrid(gridID1);
              if (gridID2 == -1) cdoAbort("Reference to horizontal grid not found!");
            }
          else
            {
              if (gridtype == GRID_CURVILINEAR)
                {
                  gridID2 = (gridtype1 == GRID_CURVILINEAR) ? gridID1 : gridToCurvilinear(gridID1, lbounds);
                }
              else if (gridtype == GRID_UNSTRUCTURED)
                {
                  bool ligme = false;
                  if (gridtype1 == GRID_GME) ligme = true;
                  gridID2 = gridToUnstructured(gridID1, 1);

                  if (ligme)
                    {
                      grid2_nvgp = gridInqSize(gridID2);
                      grid2_vgpm.resize(grid2_nvgp);
                      gridInqMaskGME(gridID2, grid2_vgpm.data());
                      gridCompress(gridID2);
                    }
                }
              else if (gridtype == GRID_LONLAT && gridtype1 == GRID_CURVILINEAR)
                {
                  gridID2 = gridCurvilinearToRegular(gridID1);
                  if (gridID2 == -1) cdoWarning("Conversion of curvilinear grid to regular grid failed!");
                }
              else if (gridtype == GRID_LONLAT && gridtype1 == GRID_UNSTRUCTURED)
                {
                  gridID2 = -1;
                  if (gridID2 == -1) cdoWarning("Conversion of unstructured grid to regular grid failed!");
                }
              else if (gridtype == GRID_LONLAT && gridtype1 == GRID_GENERIC)
                {
                  gridID2 = -1;
                  if (gridID2 == -1) cdoWarning("Conversion of generic grid to regular grid failed!");
                }
              else if (gridtype == GRID_LONLAT && gridtype1 == GRID_LONLAT)
                {
                  gridID2 = gridID1;
                }
              else
                cdoAbort("Unsupported grid name: %s", gridname);
            }

          if (gridID2 == -1)
            {
              if (!(lregular || lregularnn)) cdoAbort("Unsupported grid type!");
            }

          if (gridID2 != -1)
            {
              if (lregular || lregularnn) lrgrid = true;
              vlistChangeGridIndex(vlistID2, index, gridID2);
            }
        }

      if ((lregular || lregularnn) && !lrgrid) cdoWarning("No reduced Gaussian grid found!");
    }
  else if (operatorID == SETGRIDAREA)
    {
      int ngrids = vlistNgrids(vlistID1);
      for (int index = 0; index < ngrids; index++)
        {
          int gridID1 = vlistGrid(vlistID1, index);
          size_t gridsize = gridInqSize(gridID1);
          if (gridsize == areasize)
            {
              gridID2 = gridDuplicate(gridID1);
              gridDefArea(gridID2, areaweight.data());
              vlistChangeGridIndex(vlistID2, index, gridID2);
            }
        }
    }
  else if (operatorID == SETGRIDMASK)
    {
      int ngrids = vlistNgrids(vlistID1);
      for (int index = 0; index < ngrids; index++)
        {
          int gridID1 = vlistGrid(vlistID1, index);
          size_t gridsize = gridInqSize(gridID1);
          if (gridsize == masksize)
            {
              std::vector<int> mask(masksize);
              for (size_t i = 0; i < masksize; i++)
                {
                  if (gridmask[i] < 0 || gridmask[i] > 255)
                    mask[i] = 0;
                  else
                    mask[i] = (int) lround(gridmask[i]);
                }
              gridID2 = gridDuplicate(gridID1);
              gridDefMask(gridID2, mask.data());
              vlistChangeGridIndex(vlistID2, index, gridID2);
            }
        }
    }
  else if (operatorID == UNSETGRIDMASK)
    {
      int ngrids = vlistNgrids(vlistID1);
      for (int index = 0; index < ngrids; index++)
        {
          int gridID1 = vlistGrid(vlistID1, index);
          gridID2 = gridDuplicate(gridID1);
          gridDefMask(gridID2, NULL);
          vlistChangeGridIndex(vlistID2, index, gridID2);
        }
    }

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());

  pstreamDefVlist(streamID2, vlistID2);
  // vlistPrint(vlistID2);

  size_t gridsize = (lregular || lregularnn) ? vlistGridsizeMax(vlistID2) : vlistGridsizeMax(vlistID1);

  if (vlistNumber(vlistID1) != CDI_REAL) gridsize *= 2;
  std::vector<double> array(gridsize);

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamDefRecord(streamID2, varID, levelID);

          pstreamReadRecord(streamID1, array.data(), &nmiss);

          int gridID1 = vlistInqVarGrid(vlistID1, varID);
          if (lregular || lregularnn)
            {
              gridID2 = vlistInqVarGrid(vlistID2, varID);
              if (gridInqType(gridID1) == GRID_GAUSSIAN_REDUCED)
                {
                  double missval = vlistInqVarMissval(vlistID1, varID);
                  int lnearst = lregularnn ? 1 : 0;
                  field2regular(gridID1, gridID2, missval, array.data(), nmiss, lnearst);
                }
            }
          else if (gridInqType(gridID1) == GRID_GME)
            {
              size_t gridsize = gridInqSize(gridID1);
              size_t j = 0;
              for (size_t i = 0; i < gridsize; i++)
                if (grid2_vgpm[i]) array[j++] = array[i];
            }

          pstreamWriteRecord(streamID2, array.data(), nmiss);
        }

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

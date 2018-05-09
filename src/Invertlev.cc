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

      Invertlev     invertlev       Invert level
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

static void
invertLevDes(int vlistID)
{
  int nzaxis = vlistNzaxis(vlistID);
  for (int index = 0; index < nzaxis; index++)
    {
      int zaxisID1 = vlistZaxis(vlistID, index);
      int zaxisID2 = zaxisDuplicate(zaxisID1);
      int zaxistype = zaxisInqType(zaxisID1);

      int nlev = zaxisInqSize(zaxisID1);
      if (nlev <= 1) continue;

      if (zaxisInqLevels(zaxisID1, NULL))
        {
          std::vector<double> yv1(nlev);
          std::vector<double> yv2(nlev);
          zaxisInqLevels(zaxisID1, yv1.data());
          for (int ilev = 0; ilev < nlev; ++ilev) yv2[nlev - ilev - 1] = yv1[ilev];
          zaxisDefLevels(zaxisID2, yv2.data());
        }

      if (zaxisInqLbounds(zaxisID1, NULL) && zaxisInqUbounds(zaxisID1, NULL))
        {
          std::vector<double> yb1(nlev);
          std::vector<double> yb2(nlev);
          zaxisInqLbounds(zaxisID1, yb1.data());
          for (int ilev = 0; ilev < nlev; ++ilev) yb2[nlev - ilev - 1] = yb1[ilev];
          zaxisDefLbounds(zaxisID2, yb2.data());

          zaxisInqUbounds(zaxisID1, yb1.data());
          for (int ilev = 0; ilev < nlev; ++ilev) yb2[nlev - ilev - 1] = yb1[ilev];
          zaxisDefUbounds(zaxisID2, yb2.data());
        }

      if (zaxistype == ZAXIS_HYBRID || zaxistype == ZAXIS_HYBRID_HALF)
        {
          int vctsize = zaxisInqVctSize(zaxisID1);
          if (vctsize && vctsize % 2 == 0)
            {
              std::vector<double> vct1(vctsize);
              std::vector<double> vct2(vctsize);
              zaxisInqVct(zaxisID1, vct1.data());
              for (int i = 0; i < vctsize / 2; ++i)
                {
                  vct2[vctsize / 2 - 1 - i] = vct1[i];
                  vct2[vctsize - 1 - i] = vct1[vctsize / 2 + i];
                }
              zaxisDefVct(zaxisID2, vctsize, vct2.data());
            }
        }

      vlistChangeZaxis(vlistID, zaxisID1, zaxisID2);
    }
}

void *
Invertlev(void *process)
{
  int nrecs;
  int varID, levelID;
  size_t nmiss;
  int nlev, nlevel;
  int gridID, zaxisID;
  bool linvert = false;

  cdoInitialize(process);

  bool lcopy = UNCHANGED_RECORD;

  cdoOperatorAdd("invertlev", func_all, 0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if (operfunc == func_all || operfunc == func_hrd) invertLevDes(vlistID2);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());

  pstreamDefVlist(streamID2, vlistID2);

  size_t gridsizemax = vlistGridsizeMax(vlistID1);
  std::vector<double> array(gridsizemax);

  int nvars = vlistNvars(vlistID1);

  std::vector<std::vector<double>> vardata(nvars);
  std::vector<std::vector<size_t>> varnmiss(nvars);

  for (varID = 0; varID < nvars; varID++)
    {
      gridID = vlistInqVarGrid(vlistID1, varID);
      zaxisID = vlistInqVarZaxis(vlistID1, varID);
      size_t gridsize = gridInqSize(gridID);
      nlev = zaxisInqSize(zaxisID);

      if (nlev > 1)
        {
          linvert = true;
          vardata[varID].resize(gridsize * nlev);
          varnmiss[varID].resize(nlev);
        }
    }

  if (linvert == false) cdoWarning("No variables with invertable levels found!");

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      pstreamDefTimestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);

          if (vardata[varID].size())
            {
              gridID = vlistInqVarGrid(vlistID1, varID);
              size_t gridsize = gridInqSize(gridID);
              size_t offset = gridsize * levelID;

              pstreamReadRecord(streamID1, &vardata[varID][offset], &nmiss);
              varnmiss[varID][levelID] = nmiss;
            }
          else
            {
              pstreamDefRecord(streamID2, varID, levelID);
              if (lcopy)
                {
                  pstreamCopyRecord(streamID2, streamID1);
                }
              else
                {
                  pstreamReadRecord(streamID1, array.data(), &nmiss);
                  pstreamWriteRecord(streamID2, array.data(), nmiss);
                }
            }
        }

      for (varID = 0; varID < nvars; varID++)
        {
          if (vardata[varID].size())
            {
              gridID = vlistInqVarGrid(vlistID1, varID);
              zaxisID = vlistInqVarZaxis(vlistID1, varID);
              size_t gridsize = gridInqSize(gridID);
              nlevel = zaxisInqSize(zaxisID);
              for (levelID = 0; levelID < nlevel; levelID++)
                {
                  pstreamDefRecord(streamID2, varID, levelID);

                  size_t offset = gridsize * (nlevel - levelID - 1);
                  nmiss = varnmiss[varID][nlevel - levelID - 1];

                  pstreamWriteRecord(streamID2, &vardata[varID][offset], nmiss);
                }
            }
        }

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

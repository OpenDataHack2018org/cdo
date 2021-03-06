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

void *
Seloperator(void *process)
{
  int nrecs;
  int varID, levelID;
  int nlevs, code, zaxisID, selfound = FALSE;
  int levID, ltype = 0;
  int varID2, levelID2;
  bool sellevel, selcode, selltype;
  size_t nmiss;
  double slevel = 0, level;

  cdoInitialize(process);

  bool lcopy = UNCHANGED_RECORD;

  operatorInputArg("code, ltype, level");

  int scode = parameter2int(operatorArgv()[0]);
  int sltype = parameter2int(operatorArgv()[1]);

  if (operatorArgc() == 3) slevel = parameter2double(operatorArgv()[2]);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);

  int nvars = vlistNvars(vlistID1);
  for (varID = 0; varID < nvars; varID++)
    {
      code = vlistInqVarCode(vlistID1, varID);
      zaxisID = vlistInqVarZaxis(vlistID1, varID);
      nlevs = zaxisInqSize(zaxisID);

      ltype = zaxis2ltype(zaxisID);

      for (levID = 0; levID < nlevs; levID++)
        {
          level = cdoZaxisInqLevel(zaxisID, levID);

          if (operatorArgc() == 3)
            sellevel = IS_EQUAL(level, slevel);
          else
            sellevel = true;

          if (scode == -1 || scode == code)
            selcode = true;
          else
            selcode = false;

          if (sltype == -1 || sltype == ltype)
            selltype = true;
          else
            selltype = false;

          if (selcode && selltype && sellevel)
            {
              vlistDefFlag(vlistID1, varID, levID, true);
              selfound = true;
            }
        }
    }

  if (selfound == FALSE) cdoWarning("Code %d, ltype %d, level %g not found!", scode, sltype, slevel);

  int vlistID2 = vlistCreate();
  cdoVlistCopyFlag(vlistID2, vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  std::vector<double> array;
  if (!lcopy)
    {
      size_t gridsizemax = vlistGridsizeMax(vlistID1);
      array.resize(gridsizemax);
    }

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          if (vlistInqFlag(vlistID1, varID, levelID) == TRUE)
            {
              varID2 = vlistFindVar(vlistID2, varID);
              levelID2 = vlistFindLevel(vlistID2, varID, levelID);

              pstreamDefRecord(streamID2, varID2, levelID2);

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

      tsID++;
    }

  pstreamClose(streamID1);
  pstreamClose(streamID2);

  vlistDestroy(vlistID2);

  cdoFinish();

  return 0;
}

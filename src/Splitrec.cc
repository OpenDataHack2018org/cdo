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

      Split      splitrec        Split records
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

void *
Splitrec(void *process)
{
  int varID;
  int levelID;
  char filesuffix[32];
  char filename[8192];
  size_t nmiss;
  double *array = NULL;

  cdoInitialize(process);

  if (processSelf().m_ID != 0)
    cdoAbort("This operator can't be combined with other operators!");

  bool lcopy = UNCHANGED_RECORD;

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));
  int vlistID1 = cdoStreamInqVlist(streamID1);

  int nrecs = vlistNrecs(vlistID1);

  strcpy(filename, cdoGetObase());
  int nchars = strlen(filename);

  const char *refname = cdoGetObase();
  filesuffix[0] = 0;
  cdoGenFileSuffix(filesuffix, sizeof(filesuffix),
                   pstreamInqFiletype(streamID1), vlistID1, refname);

  if (!lcopy)
    {
      size_t gridsizemax = vlistGridsizeMax(vlistID1);
      if (vlistNumber(vlistID1) != CDI_REAL) gridsizemax *= 2;
      array = (double *) Malloc(gridsizemax * sizeof(double));
    }

  int index = 0;
  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);

          vlistClearFlag(vlistID1);
          vlistDefFlag(vlistID1, varID, levelID, TRUE);

          int vlistID2 = vlistCreate();
          cdoVlistCopyFlag(vlistID2, vlistID1);

          index++;
          sprintf(filename + nchars, "%06d", index);
          if (filesuffix[0]) sprintf(filename + nchars + 6, "%s", filesuffix);

          if (cdoVerbose) cdoPrint("create file %s", filename);

          int streamID2 = cdoStreamOpenWrite(filename, cdoFiletype());

          pstreamDefVlist(streamID2, vlistID2);

          int varID2 = vlistFindVar(vlistID2, varID);
          int levelID2 = vlistFindLevel(vlistID2, varID, levelID);

          pstreamDefTimestep(streamID2, 0);
          pstreamDefRecord(streamID2, varID2, levelID2);
          if (lcopy)
            {
              pstreamCopyRecord(streamID2, streamID1);
            }
          else
            {
              pstreamReadRecord(streamID1, array, &nmiss);
              pstreamWriteRecord(streamID2, array, nmiss);
            }

          pstreamClose(streamID2);
          vlistDestroy(vlistID2);
        }

      tsID++;
    }

  pstreamClose(streamID1);

  if (!lcopy)
    if (array) Free(array);

  cdoFinish();

  return 0;
}

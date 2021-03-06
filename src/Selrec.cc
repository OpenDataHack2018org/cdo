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

      Selrec     selrec          Select records
*/

#include <cdi.h>

#include "cdo_int.h" /* processSelf */
#include "pstream_int.h"
#include "listarray.h"

void *
Selrec(void *process)
{
  int nrecs;
  int varID, levelID;
  ListArray<int> listArrayInt;

  cdoInitialize(process);

  if (processSelf().m_ID != 0) cdoAbort("This operator can't be combined with other operators!");

  operatorInputArg("records");

  int nsel = listArrayInt.argvToInt(operatorArgc(), operatorArgv());

  int *intarr = listArrayInt.data();

  if (cdoVerbose)
    {
      for (int i = 0; i < nsel; i++) cdoPrint("intarr entry: %d %d", i, intarr[i]);
    }

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int filetype = pstreamInqFiletype(streamID1);

  if (filetype == CDI_FILETYPE_NC || filetype == CDI_FILETYPE_NC2 || filetype == CDI_FILETYPE_NC4 || filetype == CDI_FILETYPE_NC4C
      || filetype == CDI_FILETYPE_NC5)
    cdoAbort("This operator does not work on NetCDF data!");

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int recordID = 0;
  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          recordID++;
          pstreamInqRecord(streamID1, &varID, &levelID);

          for (int i = 0; i < nsel; i++)
            {
              if (recordID == intarr[i])
                {
                  pstreamDefRecord(streamID2, varID, levelID);
                  pstreamCopyRecord(streamID2, streamID1);

                  break;
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

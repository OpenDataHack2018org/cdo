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

      Merge      mergetime       Merge datasets sorted by date and time
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "util_files.h"

void *
Mergetime(void *process)
{
  int tsID2 = 0;
  int taxisID2 = CDI_UNDEFID;
  int64_t last_vdate = -1;
  int last_vtime = -1;
  bool skip_same_time = false;
  struct sfile_t
  {
    int streamID;
    int vlistID;
    int taxisID;
    int tsID;
    int64_t vdate;
    int vtime;
    int nrecs;
  };

  cdoInitialize(process);

  char *envstr = getenv("SKIP_SAME_TIME");
  if (envstr)
    {
      int ival = atoi(envstr);
      if (ival == 1)
        {
          skip_same_time = true;
          if (cdoVerbose) cdoPrint("Set SKIP_SAME_TIME to %d", ival);
        }
    }

  bool lcopy = UNCHANGED_RECORD;

  int nfiles = cdoStreamCnt() - 1;
  std::vector<sfile_t> sf(nfiles);

  for (int fileID = 0; fileID < nfiles; fileID++)
    {
      if (cdoVerbose) cdoPrint("process: %s", cdoGetStreamName(fileID).c_str());

      sf[fileID].streamID = cdoStreamOpenRead(cdoStreamName(fileID));
      sf[fileID].vlistID = cdoStreamInqVlist(sf[fileID].streamID);
      sf[fileID].taxisID = vlistInqTaxis(sf[fileID].vlistID);
    }

  // check that the contents is always the same
  for (int fileID = 1; fileID < nfiles; fileID++) vlistCompare(sf[0].vlistID, sf[fileID].vlistID, CMP_ALL);

  // read the first time step
  for (int fileID = 0; fileID < nfiles; fileID++)
    {
      sf[fileID].tsID = 0;
      sf[fileID].nrecs = cdoStreamInqTimestep(sf[fileID].streamID, sf[fileID].tsID);
      if (sf[fileID].nrecs == 0)
        {
          pstreamClose(sf[fileID].streamID);
          sf[fileID].streamID = -1;
        }
      else
        {
          sf[fileID].vdate = taxisInqVdate(sf[fileID].taxisID);
          sf[fileID].vtime = taxisInqVtime(sf[fileID].taxisID);
        }
    }

  const char *ofilename = cdoGetStreamName(nfiles).c_str();

  if (!cdoOverwriteMode && fileExists(ofilename) && !userFileOverwrite(ofilename))
    cdoAbort("Outputfile %s already exists!", ofilename);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(nfiles), cdoFiletype());

  std::vector<double> array;
  if (!lcopy)
    {
      size_t gridsizemax = vlistGridsizeMax(sf[0].vlistID);
      array.resize(gridsizemax);
    }

  while (true)
    {
      bool process_timestep = true;

      int next_fileID = -1;
      int64_t vdate = 0;
      int vtime = 0;
      for (int fileID = 0; fileID < nfiles; fileID++)
        {
          if (sf[fileID].streamID != -1)
            if (next_fileID == -1 || sf[fileID].vdate < vdate || (sf[fileID].vdate == vdate && sf[fileID].vtime < vtime))
              {
                next_fileID = fileID;
                vdate = sf[fileID].vdate;
                vtime = sf[fileID].vtime;
              }
        }

      int fileID = next_fileID;

      if (cdoVerbose) cdoPrint("nextstep = %d  vdate = %ld  vtime = %d", fileID, vdate, vtime);

      if (fileID == -1) break;

      if (skip_same_time)
        if (vdate == last_vdate && vtime == last_vtime)
          {
            char vdatestr[32], vtimestr[32];
            date2str(vdate, vdatestr, sizeof(vdatestr));
            time2str(vtime, vtimestr, sizeof(vtimestr));
            cdoPrint("Timestep %4d in stream %d (%s %s) already exists, skipped!", sf[fileID].tsID + 1, sf[fileID].streamID,
                     vdatestr, vtimestr);
            process_timestep = false;
          }

      if (process_timestep)
        {
          if (tsID2 == 0)
            {
              int vlistID1 = sf[fileID].vlistID;
              int vlistID2 = vlistDuplicate(vlistID1);
              int taxisID1 = vlistInqTaxis(vlistID1);
              taxisID2 = taxisDuplicate(taxisID1);
              vlistDefTaxis(vlistID2, taxisID2);

              pstreamDefVlist(streamID2, vlistID2);
            }

          last_vdate = vdate;
          last_vtime = vtime;

          taxisCopyTimestep(taxisID2, sf[fileID].taxisID);

          pstreamDefTimestep(streamID2, tsID2);

          for (int recID = 0; recID < sf[fileID].nrecs; recID++)
            {
              int varID, levelID;
              pstreamInqRecord(sf[fileID].streamID, &varID, &levelID);

              if (tsID2 > 0 && sf[fileID].tsID == 0)
                if (vlistInqVarTimetype(sf[fileID].vlistID, varID) == TIME_CONSTANT) continue;

              pstreamDefRecord(streamID2, varID, levelID);

              if (lcopy)
                {
                  pstreamCopyRecord(streamID2, sf[fileID].streamID);
                }
              else
                {
                  size_t nmiss;
                  pstreamReadRecord(sf[fileID].streamID, array.data(), &nmiss);
                  pstreamWriteRecord(streamID2, array.data(), nmiss);
                }
            }

          tsID2++;
        }

      sf[fileID].nrecs = cdoStreamInqTimestep(sf[fileID].streamID, ++sf[fileID].tsID);
      if (sf[fileID].nrecs == 0)
        {
          pstreamClose(sf[fileID].streamID);
          sf[fileID].streamID = -1;
        }
      else
        {
          sf[fileID].vdate = taxisInqVdate(sf[fileID].taxisID);
          sf[fileID].vtime = taxisInqVtime(sf[fileID].taxisID);
        }
    }

  pstreamClose(streamID2);

  cdoFinish();

  return 0;
}

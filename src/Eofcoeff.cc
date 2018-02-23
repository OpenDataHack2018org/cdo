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

 Eofcoeff             eofcoeff             process eof coefficients
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "grid.h"

// NO MISSING VALUE SUPPORT ADDED SO FAR

void *
Eofcoeff(void *process)
{
  char eof_name[16], oname[1024], filesuffix[32];
  double missval1 = -999, missval2;
  field_type in;
  field_type out;
  int i, varID, levelID;
  int nrecs;
  size_t nmiss;

  cdoInitialize(process);

  if (processSelf().m_ID != 0) cdoAbort("This operator can't be combined with other operators!");

  cdoOperatorAdd("eofcoeff", 0, 0, NULL);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));
  int streamID2 = cdoStreamOpenRead(cdoStreamName(1));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = cdoStreamInqVlist(streamID2);
  int vlistID3 = vlistDuplicate(vlistID2);

  // taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = vlistInqTaxis(vlistID2);
  int taxisID3 = taxisDuplicate(taxisID2);

  int gridID1 = vlistInqVarGrid(vlistID1, 0);
  int gridID2 = vlistInqVarGrid(vlistID2, 0);

  size_t gridsize = vlistGridsizeMax(vlistID1);
  if (gridsize != vlistGridsizeMax(vlistID2))
    cdoAbort("Gridsize of input files does not match! %zu and %zu", gridsize, vlistGridsizeMax(vlistID2));

  if (vlistNgrids(vlistID2) > 1 || vlistNgrids(vlistID1) > 1) cdoAbort("Too many different grids in input!");

  int nvars = vlistNvars(vlistID1) == vlistNvars(vlistID2) ? vlistNvars(vlistID1) : -1;
  int nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, 0));

  if (gridID1 != gridID2) cdoCompareGrids(gridID1, gridID2);

  strcpy(oname, cdoGetObase());
  int nchars = strlen(oname);

  const char *refname = cdoGetObase();
  filesuffix[0] = 0;
  cdoGenFileSuffix(filesuffix, sizeof(filesuffix), pstreamInqFiletype(streamID1), vlistID1, refname);

  field_type ***eof = (field_type ***) Malloc(nvars * sizeof(field_type **));
  for (varID = 0; varID < nvars; varID++)
    eof[varID] = (field_type **) Malloc(nlevs * sizeof(field_type *));

  int eofID = 0;
  while (1)
    {
      nrecs = cdoStreamInqTimestep(streamID1, eofID);
      if (nrecs == 0) break;

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          missval1 = vlistInqVarMissval(vlistID1, varID);
          if (eofID == 0)
            eof[varID][levelID] = (field_type *) Malloc(1 * sizeof(field_type));
          else
            eof[varID][levelID] = (field_type *) Realloc(eof[varID][levelID], (eofID + 1) * sizeof(field_type));
          eof[varID][levelID][eofID].grid = gridID1;
          eof[varID][levelID][eofID].nmiss = 0;
          eof[varID][levelID][eofID].missval = missval1;
          eof[varID][levelID][eofID].ptr = (double *) Malloc(gridsize * sizeof(double));
          arrayFill(gridsize, &eof[varID][levelID][eofID].ptr[0], missval1);

          if (varID >= nvars) cdoAbort("Internal error - too high varID");
          if (levelID >= nlevs) cdoAbort("Internal error - too high levelID");

          pstreamReadRecord(streamID1, eof[varID][levelID][eofID].ptr, &nmiss);
          eof[varID][levelID][eofID].nmiss = nmiss;
        }
      eofID++;
    }

  int neof = eofID;

  if (cdoVerbose) cdoPrint("%s contains %i eof's", cdoGetStreamName(0).c_str(), neof);
  // Create 1x1 Grid for output
  int gridID3 = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID3, 1);
  gridDefYsize(gridID3, 1);
  double xvals = 0;
  double yvals = 0;
  gridDefXvals(gridID3, &xvals);
  gridDefYvals(gridID3, &yvals);

  // Create var-list and time-axis for output

  int ngrids = vlistNgrids(vlistID3);
  for (i = 0; i < ngrids; i++)
    vlistChangeGridIndex(vlistID3, i, gridID3);

  vlistDefTaxis(vlistID3, taxisID3);
  for (varID = 0; varID < nvars; varID++)
    vlistDefVarTimetype(vlistID3, varID, TIME_VARYING);

  // open streams for eofcoeff output
  int *streamIDs = (int *) Malloc(neof * sizeof(int));
  for (eofID = 0; eofID < neof; eofID++)
    {
      oname[nchars] = '\0';

      sprintf(eof_name, "%5.5i", eofID);
      strcat(oname, eof_name);
      if (filesuffix[0]) strcat(oname, filesuffix);

      streamIDs[eofID] = cdoStreamOpenWrite(oname, cdoFiletype());

      if (cdoVerbose) cdoPrint("opened %s ('w')  as stream%i for %i. eof", oname, streamIDs[eofID], eofID + 1);

      pstreamDefVlist(streamIDs[eofID], vlistID3);
    }

  // ALLOCATE temporary fields for data read and write
  in.ptr = (double *) Malloc(gridsize * sizeof(double));
  in.grid = gridID1;
  out.missval = missval1;
  out.nmiss = 0;
  out.ptr = (double *) Malloc(1 * sizeof(double));

  int tsID = 0;
  while (1)
    {
      nrecs = cdoStreamInqTimestep(streamID2, tsID);
      if (nrecs == 0) break;

      taxisCopyTimestep(taxisID3, taxisID2);
      /*for ( eofID=0; eofID<neof; eofID++)
        {
          fprintf(stderr, "defining ts %i\n", tsID);
          pstreamDefTimestep(streamIDs[eofID],tsID);
        }
      */
      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID2, &varID, &levelID);
          missval2 = vlistInqVarMissval(vlistID2, varID);
          pstreamReadRecord(streamID2, in.ptr, &nmiss);
          in.nmiss = nmiss;

          for (eofID = 0; eofID < neof; eofID++)
            {
              if (recID == 0) pstreamDefTimestep(streamIDs[eofID], tsID);
              // if ( recID == 0 ) fprintf(stderr, "ts%i rec%i eof%i\n", tsID,
              // recID, eofID);
              out.ptr[0] = 0;
              out.grid = gridID3;
              out.missval = missval2;
              for (size_t i = 0; i < gridsize; i++)
                {
                  if (!DBL_IS_EQUAL(in.ptr[i], missval2) && !DBL_IS_EQUAL(eof[varID][levelID][eofID].ptr[i], missval1))
                    {
                      // double tmp =
                      // w[i]*in.ptr[i]*eof[varID][levelID][eofID].ptr[i];
                      double tmp = in.ptr[i] * eof[varID][levelID][eofID].ptr[i];
                      out.ptr[0] += tmp;
                    }
                }
              if (!DBL_IS_EQUAL(out.ptr[0], 0.))
                nmiss = 0;
              else
                {
                  nmiss = 1;
                  out.ptr[0] = missval2;
                }

              pstreamDefRecord(streamIDs[eofID], varID, levelID);
              // fprintf(stderr, "%d %d %d %d %d %g\n", streamIDs[eofID],tsID,
              // recID, varID, levelID,*out.ptr);
              pstreamWriteRecord(streamIDs[eofID], out.ptr, nmiss);
            }
          if (varID >= nvars) cdoAbort("Internal error - too high varID");
          if (levelID >= nlevs) cdoAbort("Internal error - too high levelID");
        }

      tsID++;
    }

  for (eofID = 0; eofID < neof; eofID++)
    pstreamClose(streamIDs[eofID]);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

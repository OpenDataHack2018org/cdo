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

      Splitsel   splitsel        Split time selection
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

void *
Splitsel(void *process)
{
  int nrecs = 0;
  int varID, levelID;
  int tsID;
  size_t nmiss;
  int i2 = 0;
  int nchars;
  char filesuffix[32];
  char filename[8192];

  cdoInitialize(process);

  if (processSelf().m_ID != 0) cdoAbort("This operator can't be combined with other operators!");

  bool lcopy = UNCHANGED_RECORD;

  cdoOperatorAdd("splitsel", 0, 0, NULL);

  /*  operatorInputArg("nsets <noffset <nskip>>"); */

  int nargc = operatorArgc();
  if (nargc < 1) cdoAbort("Too few arguments! Need %d found %d.", 1, nargc);

  /*   ndates = parameter2int(operatorArgv()[0]); */
  /*   if ( nargc > 1 ) noffset = parameter2int(operatorArgv()[1]); */
  /*   if ( nargc > 2 ) nskip   = parameter2int(operatorArgv()[2]); */
  /*   printf("%s %s %s\n", operatorArgv()[0],operatorArgv()[1],operatorArgv()[2]); */
  double noffset = 0.0, nskip = 0.0;
  double ndates = parameter2double(operatorArgv()[0]);
  if (nargc > 1) noffset = parameter2double(operatorArgv()[1]);
  if (nargc > 2) nskip = parameter2double(operatorArgv()[2]);

  if (cdoVerbose) cdoPrint("nsets = %f, noffset = %f, nskip = %f", ndates, noffset, nskip);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  /*   taxisID2 = taxisCreate(TAXIS_ABSOLUTE); */
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  strcpy(filename, cdoGetObase());
  nchars = strlen(filename);

  const char *refname = cdoGetObase();
  filesuffix[0] = 0;
  cdoGenFileSuffix(filesuffix, sizeof(filesuffix), pstreamInqFiletype(streamID1), vlistID1, refname);

  std::vector<double> array;
  //  if ( ! lcopy )
  {
    size_t gridsizemax = vlistGridsizeMax(vlistID1);
    if (vlistNumber(vlistID1) != CDI_REAL) gridsizemax *= 2;
    array.resize(gridsizemax);
  }

  int nvars = vlistNvars(vlistID1);
  int nconst = 0;
  for (varID = 0; varID < nvars; varID++)
    if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT) nconst++;

  std::vector<std::vector<Field>> vars;
  if (nconst)
    {
      vars.resize(nvars);

      for (varID = 0; varID < nvars; varID++)
        {
          if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT)
            {
              int gridID = vlistInqVarGrid(vlistID1, varID);
              int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
              size_t gridsize = gridInqSize(gridID);

              vars[varID].resize(nlevel);

              for (levelID = 0; levelID < nlevel; levelID++)
                {
                  field_init(&vars[varID][levelID]);
                  vars[varID][levelID].grid = gridID;
                  vars[varID][levelID].ptr = (double *) Malloc(gridsize * sizeof(double));
                }
            }
        }
    }

  int index = 0;
  int nsets = 0;

  /* offset */
  for (tsID = 0; tsID < noffset; tsID++)
    {
      nrecs = cdoStreamInqTimestep(streamID1, tsID);
      if (nrecs == 0)
        {
          cdoWarning("noffset is larger than number of timesteps!");
          goto LABEL_END;
        }

      if (tsID == 0 && nconst)
        for (int recID = 0; recID < nrecs; recID++)
          {
            pstreamInqRecord(streamID1, &varID, &levelID);
            if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT)
              {
                pstreamReadRecord(streamID1, vars[varID][levelID].ptr, &nmiss);
                vars[varID][levelID].nmiss = nmiss;
              }
          }
    }

  while (TRUE)
    {
      sprintf(filename + nchars, "%06d", index);
      sprintf(filename + nchars + 6, "%s", filesuffix);

      if (cdoVerbose) cdoPrint("create file %s", filename);

      int streamID2 = -1;
      int tsID2 = 0;

      for (; nsets < (int) (ndates * (index + 1)); nsets++)
        {
          nrecs = cdoStreamInqTimestep(streamID1, tsID);
          if (nrecs == 0) break;

          taxisCopyTimestep(taxisID2, taxisID1);

          if (streamID2 < 0)
            {
              streamID2 = cdoStreamOpenWrite(filename, cdoFiletype());
              pstreamDefVlist(streamID2, vlistID2);
            }

          pstreamDefTimestep(streamID2, tsID2);

          if (tsID > 0 && tsID2 == 0 && nconst)
            {
              for (varID = 0; varID < nvars; varID++)
                {
                  if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT)
                    {
                      int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
                      for (levelID = 0; levelID < nlevel; levelID++)
                        {
                          pstreamDefRecord(streamID2, varID, levelID);
                          nmiss = vars[varID][levelID].nmiss;
                          pstreamWriteRecord(streamID2, vars[varID][levelID].ptr, nmiss);
                        }
                    }
                }
            }

          for (int recID = 0; recID < nrecs; recID++)
            {

              pstreamInqRecord(streamID1, &varID, &levelID);
              pstreamDefRecord(streamID2, varID, levelID);
              if (lcopy && !(tsID == 0 && nconst))
                {
                  pstreamCopyRecord(streamID2, streamID1);
                }
              else
                {
                  pstreamReadRecord(streamID1, array.data(), &nmiss);
                  pstreamWriteRecord(streamID2, array.data(), nmiss);

                  if (tsID == 0 && nconst)
                    {
                      if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT)
                        {
                          int gridID = vlistInqVarGrid(vlistID1, varID);
                          size_t gridsize = gridInqSize(gridID);
                          arrayCopy(gridsize, array.data(), vars[varID][levelID].ptr);
                          vars[varID][levelID].nmiss = nmiss;
                        }
                    }
                }
            }

          tsID++;
          tsID2++;
        }

      pstreamClose(streamID2);
      if (nrecs == 0) break;

      nrecs = cdoStreamInqTimestep(streamID1, tsID);
      if (nrecs == 0) break;

      for (; i2 < (int) (nskip * (index + 1)); i2++)
        {
          nrecs = cdoStreamInqTimestep(streamID1, tsID);
          if (nrecs == 0) break;
          tsID++;
        }

      nrecs = cdoStreamInqTimestep(streamID1, tsID);
      if (nrecs == 0) break;

      index++;
    }

LABEL_END:

  pstreamClose(streamID1);

  if (nconst)
    {
      for (varID = 0; varID < nvars; varID++)
        {
          if (vlistInqVarTimetype(vlistID2, varID) == TIME_CONSTANT)
            {
              int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
              for (levelID = 0; levelID < nlevel; levelID++)
                if (vars[varID][levelID].ptr) Free(vars[varID][levelID].ptr);
            }
        }
    }

  vlistDestroy(vlistID2);

  cdoFinish();

  return 0;
}

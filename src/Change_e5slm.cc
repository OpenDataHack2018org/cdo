/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2018 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied 1warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

/*
   This module contains the following operators:

      Change_e5slm      change_e5slm          Change ECHAM5 sea land mask
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

void *
Change_e5slm(void *process)
{
  char name[CDI_MAX_NAME];
  int nrecs;
  int varID, levelID;
  size_t nmiss;

  cdoInitialize(process);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int taxisID1 = vlistInqTaxis(vlistID1);

  int vlistID2 = vlistDuplicate(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  /* get filename of SLM */
  operatorInputArg("filename of the sea land mask");
  operatorCheckArgc(1);
  const char *fn_slm = operatorArgv()[0];

  /* read SLM */
  int streamIDslm = streamOpenRead(fn_slm);

  int vlistIDslm = streamInqVlist(streamIDslm);

  size_t gridsize = gridInqSize(vlistInqVarGrid(vlistIDslm, 0));

  double *array = (double *) Malloc(gridsize * sizeof(double));
  double *cland = (double *) Malloc(gridsize * sizeof(double));
  bool *lsea = (bool *) Malloc(gridsize * sizeof(bool));

  streamInqTimestep(streamIDslm, 0);

  streamInqRecord(streamIDslm, &varID, &levelID);
  streamReadRecord(streamIDslm, cland, &nmiss);

  if (nmiss > 0) cdoAbort("SLM with missing values are unsupported!");

  double minval, maxval;
  arrayMinMaxMask(gridsize, cland, NULL, &minval, &maxval);
  if (minval < 0 || maxval > 1) cdoWarning("Values of SLM out of bounds! (minval=%g, maxval=%g)", minval, maxval);

  streamClose(streamIDslm);

  for (size_t i = 0; i < gridsize; ++i)
    lsea[i] = !(cland[i] > 0);

  int nvars = vlistNvars(vlistID1);
  short *codes = (short *) Malloc(nvars * sizeof(short));

  for (varID = 0; varID < nvars; ++varID)
    {
      if (gridsize != gridInqSize(vlistInqVarGrid(vlistID1, varID))) cdoAbort("gridsize differ!");

      int code = vlistInqVarCode(vlistID1, varID);
      vlistInqVarName(vlistID1, varID, name);

      if (code < 0)
        {
          if (strcmp(name, "SLM") == 0)
            code = 172;
          else if (strcmp(name, "ALAKE") == 0)
            code = 99;
          else if (strcmp(name, "WS") == 0)
            code = 140;
          else if (strcmp(name, "AZ0") == 0)
            code = 173;
          else if (strcmp(name, "ALB") == 0)
            code = 174;
          else if (strcmp(name, "VGRAT") == 0)
            code = 198;
          else if (strcmp(name, "FOREST") == 0)
            code = 212;
          else if (strcmp(name, "FAO") == 0)
            code = 226;
          else if (strcmp(name, "WSMX") == 0)
            code = 229;
          else if (strcmp(name, "GLAC") == 0)
            code = 232;
          else if (strcmp(name, "VLTCLIM") == 0)
            code = 71;
          else if (strcmp(name, "VGRATCLIM") == 0)
            code = 70;
        }

      codes[varID] = code;
    }

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      pstreamDefTimestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamReadRecord(streamID1, array, &nmiss);

          int code = codes[varID];
          if (code == 172)
            {
              cdoPrint("SLM changed!");
              for (size_t i = 0; i < gridsize; ++i)
                array[i] = cland[i];
            }
          else if (code == 99)
            {
              cdoPrint("ALAKE set all values to zero!");
              for (size_t i = 0; i < gridsize; ++i)
                array[i] = 0;
            }
          else if (code == 232)
            {
              cdoPrint("GLAC set sea points to %g!", array[0]);
              for (size_t i = 0; i < gridsize; ++i)
                if (cland[i] < 0.5) array[i] = array[0];
            }
          else if (code == 70 || code == 71 || code == 140 || code == 173 || code == 174 || code == 198 || code == 200
                   || code == 212 || code == 226 || code == 229)
            {
              cdoPrint("Code %d set sea points to %g!", code, array[0]);
              for (size_t i = 0; i < gridsize; ++i)
                if (lsea[i]) array[i] = array[0];
            }

          pstreamDefRecord(streamID2, varID, levelID);
          pstreamWriteRecord(streamID2, array, nmiss);
        }

      tsID++;
    }

  pstreamClose(streamID1);
  pstreamClose(streamID2);

  Free(array);
  Free(cland);
  Free(lsea);
  Free(codes);

  cdoFinish();

  return 0;
}

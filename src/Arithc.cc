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

      Arithc     addc            Add by constant
      Arithc     subc            Subtract by constant
      Arithc     mulc            Multiply by constant
      Arithc     divc            Divide by constant
      Arithc     mod             Modulo operator
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

void fill_vars(int vlistID, std::vector<bool> &vars)
{
  int varID;
  int nvars = vlistNvars(vlistID);
  vars.resize(nvars);

  if (cdoNumVarnames)
    {
      bool lfound = false;
      char varname[CDI_MAX_NAME];

      for (varID = 0; varID < nvars; ++varID)
        {
          vars[varID] = false;

          vlistInqVarName(vlistID, varID, varname);

          for (int i = 0; i < cdoNumVarnames; ++i)
            if (strcmp(varname, cdoVarnames[i]) == 0)
              {
                vars[varID] = true;
                lfound = true;
                break;
              }
        }

      if (!lfound) cdoAbort("Variable %s%s not found!", cdoVarnames[0], cdoNumVarnames > 1 ? ",..." : "");
    }
  else
    {
      for (varID = 0; varID < nvars; ++varID) vars[varID] = true;
    }
}

void *
Arithc(void *process)
{
  int nrecs;
  int varID, levelID;

  cdoInitialize(process);

  // clang-format off
  cdoOperatorAdd("addc", func_add, 1, "constant value");
  cdoOperatorAdd("subc", func_sub, 1, "constant value");
  cdoOperatorAdd("mulc", func_mul, 1, "constant value");
  cdoOperatorAdd("divc", func_div, 1, "constant value");
  cdoOperatorAdd("mod",  func_mod, 0, "divisor");
  // clang-format on

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);
  bool opercplx = cdoOperatorF2(operatorID);

  operatorInputArg(cdoOperatorEnter(operatorID));
  if (operatorArgc() < 1) cdoAbort("Too few arguments!");
  if (operatorArgc() > 2) cdoAbort("Too many arguments!");
  double rconst = parameter2double(operatorArgv()[0]);
  double rconstcplx[2];
  rconstcplx[0] = rconst;
  rconstcplx[1] = 0;
  if (operatorArgc() == 2) rconstcplx[1] = parameter2double(operatorArgv()[1]);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  std::vector<bool> vars;
  fill_vars(vlistID1, vars);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  size_t nwpv = (vlistNumber(vlistID1) == CDI_COMP) ? 2 : 1;
  if (nwpv == 2 && opercplx == false) cdoAbort("Fields with complex numbers are not supported by this operator!");
  size_t gridsizemax = nwpv * vlistGridsizeMax(vlistID1);

  Field field;
  field_init(&field);
  field.ptr = (double *) Malloc(gridsizemax * sizeof(double));
  field.weight = NULL;

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamReadRecord(streamID1, field.ptr, &field.nmiss);

          if (vars[varID])
            {
              int gridID = vlistInqVarGrid(vlistID1, varID);
              size_t gridsize = nwpv * gridInqSize(gridID);
              field.grid = gridID;
              field.missval = vlistInqVarMissval(vlistID1, varID);

              if (nwpv == 2)
                farcfuncplx(&field, rconstcplx, operfunc);
              else
                farcfun(&field, rconst, operfunc);

              // recalculate number of missing values
              field.nmiss = arrayNumMV(gridsize, field.ptr, field.missval);
            }

          pstreamDefRecord(streamID2, varID, levelID);
          pstreamWriteRecord(streamID2, field.ptr, field.nmiss);
        }

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if (field.ptr) Free(field.ptr);

  vlistDestroy(vlistID2);

  cdoFinish();

  return 0;
}

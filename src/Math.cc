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

      Math       abs             Absolute value
      Math       sqr             Square
      Math       sqrt            Square root
      Math       exp             Exponential
      Math       ln              Natural logarithm
      Math       log10           Base 10 logarithm
      Math       sin             Sine
      Math       cos             Cosine
      Math       tan             Tangent
      Math       asin            Arc sine
      Math       acos            Arc cosine
      Math       atan            Arc tangent
      Math       pow             Power
      Math       reci            Reciprocal
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

void *
Math(void *process)
{
  enum
  {
    ABS,
    FINT,
    FNINT,
    SQR,
    SQRT,
    EXP,
    LN,
    LOG10,
    SIN,
    COS,
    TAN,
    ASIN,
    ACOS,
    ATAN,
    POW,
    RECI,
    NOT,
    CONJ,
    RE,
    IM,
    ARG
  };
  int nrecs;
  int varID, levelID;
  size_t nmiss;
  size_t i;

  cdoInitialize(process);

  // clang-format off
  cdoOperatorAdd("abs",   ABS,   0, NULL);
  cdoOperatorAdd("int",   FINT,  0, NULL);
  cdoOperatorAdd("nint",  FNINT, 0, NULL);
  cdoOperatorAdd("sqr",   SQR,   0, NULL);
  cdoOperatorAdd("sqrt",  SQRT,  0, NULL);
  cdoOperatorAdd("exp",   EXP,   0, NULL);
  cdoOperatorAdd("ln",    LN,    0, NULL);
  cdoOperatorAdd("log10", LOG10, 0, NULL);
  cdoOperatorAdd("sin",   SIN,   0, NULL);
  cdoOperatorAdd("cos",   COS,   0, NULL);
  cdoOperatorAdd("tan",   TAN,   0, NULL);
  cdoOperatorAdd("asin",  ASIN,  0, NULL);
  cdoOperatorAdd("acos",  ACOS,  0, NULL);
  cdoOperatorAdd("atan",  ATAN,  0, NULL);
  cdoOperatorAdd("pow",   POW,   0, NULL);
  cdoOperatorAdd("reci",  RECI,  0, NULL);
  cdoOperatorAdd("not",   NOT,   0, NULL);
  cdoOperatorAdd("conj",  CONJ,  0, NULL);
  cdoOperatorAdd("re",    RE,    0, NULL);
  cdoOperatorAdd("im",    IM,    0, NULL);
  cdoOperatorAdd("arg",   ARG,   0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  if (operfunc == FNINT) cdo_check_round();

  double rc = 0;
  if (operfunc == POW)
    {
      operatorInputArg("value");
      rc = parameter2double(operatorArgv()[0]);
    }

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  if (operfunc == RE || operfunc == IM || operfunc == ABS || operfunc == ARG)
    {
      int nvars = vlistNvars(vlistID2);
      for (int varID = 0; varID < nvars; ++varID)
        {
          if (vlistInqVarDatatype(vlistID2, varID) == CDI_DATATYPE_CPX32) vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_FLT32);
          if (vlistInqVarDatatype(vlistID2, varID) == CDI_DATATYPE_CPX64) vlistDefVarDatatype(vlistID2, varID, CDI_DATATYPE_FLT64);
        }
    }

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  size_t gridsizemax = vlistGridsizeMax(vlistID1);
  if (vlistNumber(vlistID1) != CDI_REAL) gridsizemax *= 2;

  std::vector<double> array1(gridsizemax);
  std::vector<double> array2(gridsizemax);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamReadRecord(streamID1, &array1[0], &nmiss);

          double missval1 = vlistInqVarMissval(vlistID1, varID);
          double missval2 = missval1;
          size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
          int number = vlistInqVarNumber(vlistID1, varID);

          if (number == CDI_REAL)
            {
              switch (operfunc)
                {
                case ABS:
                  for (i = 0; i < gridsize; i++) array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : fabs(array1[i]);
                  break;
                case FINT:
                  for (i = 0; i < gridsize; i++) array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : (int) (array1[i]);
                  break;
                case FNINT:
                  for (i = 0; i < gridsize; i++) array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : round(array1[i]);
                  break;
                case SQR:
                  for (i = 0; i < gridsize; i++) array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : array1[i] * array1[i];
                  break;
                case SQRT:
                  for (i = 0; i < gridsize; i++) array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : SQRTMN(array1[i]);
                  break;
                case EXP:
                  for (i = 0; i < gridsize; i++) array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : exp(array1[i]);
                  break;
                case LN:
                  for (i = 0; i < gridsize; i++)
                    array2[i] = DBL_IS_EQUAL(array1[i], missval1) || array1[i] < 0 ? missval1 : log(array1[i]);
                  break;
                case LOG10:
                  for (i = 0; i < gridsize; i++)
                    array2[i] = DBL_IS_EQUAL(array1[i], missval1) || array1[i] < 0 ? missval1 : log10(array1[i]);
                  break;
                case SIN:
                  for (i = 0; i < gridsize; i++) array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : sin(array1[i]);
                  break;
                case COS:
                  for (i = 0; i < gridsize; i++) array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : cos(array1[i]);
                  break;
                case TAN:
                  for (i = 0; i < gridsize; i++) array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : tan(array1[i]);
                  break;
                case ASIN:
                  for (i = 0; i < gridsize; i++)
                    array2[i] = DBL_IS_EQUAL(array1[i], missval1) || array1[i] < -1 || array1[i] > 1 ? missval1 : asin(array1[i]);
                  break;
                case ACOS:
                  for (i = 0; i < gridsize; i++)
                    array2[i] = DBL_IS_EQUAL(array1[i], missval1) || array1[i] < -1 || array1[i] > 1 ? missval1 : acos(array1[i]);
                  break;
                case ATAN:
                  for (i = 0; i < gridsize; i++) array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : atan(array1[i]);
                  break;
                case POW:
                  for (i = 0; i < gridsize; i++) array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : pow(array1[i], rc);
                  break;
                case RECI:
                  for (i = 0; i < gridsize; i++)
                    array2[i] = DBL_IS_EQUAL(array1[i], missval1) || DBL_IS_EQUAL(array1[i], 0.) ? missval1 : 1 / array1[i];
                  break;
                case NOT:
                  for (i = 0; i < gridsize; i++) array2[i] = DBL_IS_EQUAL(array1[i], missval1) ? missval1 : IS_EQUAL(array1[i], 0);
                  break;
                case RE:
                case ARG:
                  for (i = 0; i < gridsize; i++) array2[i] = array1[i];
                  break;
                default: cdoAbort("Operator not implemented for real data!"); break;
                }

              nmiss = arrayNumMV(gridsize, &array2[0], missval1);
            }
          else
            {
              switch (operfunc)
                {
                case SQR:
                  for (i = 0; i < gridsize; i++)
                    {
                      array2[i * 2] = array1[i * 2] * array1[i * 2] + array1[i * 2 + 1] * array1[i * 2 + 1];
                      array2[i * 2 + 1] = 0;
                    }
                  break;
                case SQRT:
                  for (i = 0; i < gridsize; i++)
                    {
                      double abs = SQRTMN(ADDMN(MULMN(array1[2 * i], array1[2 * i]), MULMN(array1[2 * i + 1], array1[2 * i + 1])));
                      array2[i * 2] = MULMN(1 / sqrt(2.), SQRTMN(ADDMN(array1[i * 2], abs)));
                      array2[i * 2 + 1] = MULMN(1 / sqrt(2.), DIVMN(array1[2 * i + 1], SQRTMN(ADDMN(array1[2 * i], abs))));
                      ;
                    }
                  break;
                case CONJ:
                  for (i = 0; i < gridsize; i++)
                    {
                      array2[i * 2] = array1[i * 2];
                      array2[i * 2 + 1] = -array1[i * 2 + 1];
                    }
                  break;
                case RE:
                  for (i = 0; i < gridsize; i++) array2[i] = array1[i * 2];
                  break;
                case IM:
                  for (i = 0; i < gridsize; i++) array2[i] = array1[i * 2 + 1];
                  break;
                case ABS:
                  for (i = 0; i < gridsize; i++)
                    array2[i] = SQRTMN(ADDMN(MULMN(array1[2 * i], array1[2 * i]), MULMN(array1[2 * i + 1], array1[2 * i + 1])));
                  break;
                case ARG:
                  for (i = 0; i < gridsize; i++)
                    array2[i] = (DBL_IS_EQUAL(array1[2 * i], missval1) || DBL_IS_EQUAL(array1[2 * i + 1], missval1))
                                    ? missval1
                                    : atan2(array1[2 * i + 1], array1[2 * i]);
                  break;
                default: cdoAbort("Fields with complex numbers are not supported by this operator!"); break;
                }

              nmiss = 0;
            }

          pstreamDefRecord(streamID2, varID, levelID);
          pstreamWriteRecord(streamID2, &array2[0], nmiss);
        }

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  vlistDestroy(vlistID2);

  cdoFinish();

  return 0;
}

/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include <assert.h>
#include <string.h>

#include <cdi.h>

#include "cdo_int.h"
#include "grid.h"
#include "pstream_int.h"
#include "ecacore.h"
#include "ecautil.h"

#define FIRST_VAR_ID 0

#define IS_NOT_SET(x) (x == NULL)
#define IS_SET(x) (x != NULL)

void
eca1(const ECA_REQUEST_1 *request)
{
  const int operatorID = cdoOperatorID();

  size_t nmiss;
  int cmplen;
  char indate1[DATE_LEN + 1], indate2[DATE_LEN + 1];
  int64_t ivdate = 0, ovdate = 0;
  int ivtime = 0, ovtime = 0;
  int nrecs, nrecords;
  int gridID, zaxisID, varID, levelID;
  int itsID;
  int otsID;
  long nsets;
  int istreamID, ostreamID;
  int ivlistID, ovlistID, itaxisID, otaxisID;
  int nlevels;
  int *recVarID, *recLevelID;
  double missval;
  Field *var12 = NULL, *samp1 = NULL, *samp2 = NULL, *var13 = NULL, *var21 = NULL, *var23 = NULL, *var;
  Field field1, field2;

  cmplen = DATE_LEN - cdoOperatorF2(operatorID);

  istreamID = cdoStreamOpenRead(cdoStreamName(0));

  ivlistID = cdoStreamInqVlist(istreamID);
  ovlistID = vlistCreate();

  gridID = vlistInqVarGrid(ivlistID, FIRST_VAR_ID);
  zaxisID = vlistInqVarZaxis(ivlistID, FIRST_VAR_ID);
  missval = vlistInqVarMissval(ivlistID, FIRST_VAR_ID);

  varID = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARYING);

  vlistDefVarMissval(ovlistID, varID, missval);

  if (IS_SET(request->var1.name)) vlistDefVarName(ovlistID, varID, request->var1.name);
  if (IS_SET(request->var1.longname)) vlistDefVarLongname(ovlistID, varID, request->var1.longname);
  if (IS_SET(request->var1.units)) vlistDefVarUnits(ovlistID, varID, request->var1.units);

  if (IS_SET(request->var2.h2) || IS_SET(request->var2.h3))
    {
      varID = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARYING);

      vlistDefVarMissval(ovlistID, varID, missval);

      if (IS_SET(request->var2.name)) vlistDefVarName(ovlistID, varID, request->var2.name);
      if (IS_SET(request->var2.longname)) vlistDefVarLongname(ovlistID, varID, request->var2.longname);
      if (IS_SET(request->var2.units)) vlistDefVarUnits(ovlistID, varID, request->var2.units);
    }

  if (cdoOperatorF2(operatorID) == 16) vlistDefNtsteps(ovlistID, 1);

  itaxisID = vlistInqTaxis(ivlistID);
  otaxisID = taxisCreate(TAXIS_RELATIVE);
  taxisDefTunit(otaxisID, TUNIT_MINUTE);
  //  taxisDefCalendar(otaxisID, CALENDAR_PROLEPTIC);
  taxisDefCalendar(otaxisID, taxisInqCalendar(itaxisID));
  taxisDefRdate(otaxisID, 19550101);
  taxisDefRtime(otaxisID, 0);
  vlistDefTaxis(ovlistID, otaxisID);

  ostreamID = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(ostreamID, ovlistID);

  nrecords = vlistNrecs(ivlistID);
  recVarID = (int *) Malloc(nrecords * sizeof(int));
  recLevelID = (int *) Malloc(nrecords * sizeof(int));

  size_t gridsize = gridInqSize(gridID);

  field_init(&field1);
  field_init(&field2);

  field1.ptr = (double *) Malloc(gridsize * sizeof(double));
  if (IS_SET(request->var2.h2) || IS_SET(request->var2.h3))
    field2.ptr = (double *) Malloc(gridsize * sizeof(double));
  else
    field2.ptr = NULL;

  nlevels = zaxisInqSize(zaxisID);

  var12 = (Field *) Malloc(nlevels * sizeof(Field));
  samp1 = (Field *) Malloc(nlevels * sizeof(Field));
  samp2 = (Field *) Malloc(nlevels * sizeof(Field));
  if (IS_SET(request->var1.f3)) var13 = (Field *) Malloc(nlevels * sizeof(Field));

  if (IS_SET(request->var2.h2)) var21 = (Field *) Malloc(nlevels * sizeof(Field));
  if (IS_SET(request->var2.h3)) var23 = (Field *) Malloc(nlevels * sizeof(Field));

  for (levelID = 0; levelID < nlevels; levelID++)
    {
      field_init(&var12[levelID]);
      var12[levelID].grid = gridID;
      var12[levelID].nmiss = 0;
      var12[levelID].missval = missval;
      var12[levelID].ptr = (double *) Malloc(gridsize * sizeof(double));

      field_init(&samp1[levelID]);
      samp1[levelID].grid = gridID;
      samp1[levelID].nmiss = 0;
      samp1[levelID].missval = missval;
      samp1[levelID].ptr = (double *) Malloc(gridsize * sizeof(double));

      field_init(&samp2[levelID]);
      samp2[levelID].grid = gridID;
      samp2[levelID].nmiss = 0;
      samp2[levelID].missval = missval;
      samp2[levelID].ptr = NULL;

      if (IS_SET(request->var1.f3))
        {
          field_init(&var13[levelID]);
          var13[levelID].grid = gridID;
          var13[levelID].nmiss = 0;
          var13[levelID].missval = missval;
          var13[levelID].ptr = (double *) Malloc(gridsize * sizeof(double));
        }
      if (IS_SET(request->var2.h2))
        {
          field_init(&var21[levelID]);
          var21[levelID].grid = gridID;
          var21[levelID].nmiss = 0;
          var21[levelID].missval = missval;
          var21[levelID].ptr = (double *) Malloc(gridsize * sizeof(double));
        }
      if (IS_SET(request->var2.h3))
        {
          field_init(&var23[levelID]);
          var23[levelID].grid = gridID;
          var23[levelID].nmiss = 0;
          var23[levelID].missval = missval;
          var23[levelID].ptr = (double *) Malloc(gridsize * sizeof(double));
        }
    }

  itsID = 0;
  otsID = 0;
  while (TRUE)
    {
      nsets = 0;
      while ((nrecs = cdoStreamInqTimestep(istreamID, itsID)) > 0)
        {
          ivdate = taxisInqVdate(itaxisID);
          ivtime = taxisInqVtime(itaxisID);

          if (nsets == 0) SET_DATE(indate2, ivdate, ivtime);
          SET_DATE(indate1, ivdate, ivtime);

          if (DATE_IS_NEQ(indate1, indate2, cmplen)) break;

          for (int recID = 0; recID < nrecs; recID++)
            {
              pstreamInqRecord(istreamID, &varID, &levelID);

              if (itsID == 0)
                {
                  recVarID[recID] = varID;
                  recLevelID[recID] = levelID;
                }
              if (varID != FIRST_VAR_ID) continue;

              if (nsets == 0)
                {
                  for (size_t i = 0; i < gridsize; i++)
                    {
                      var12[levelID].ptr[i] = missval;
                      samp1[levelID].ptr[i] = missval;
                      if (IS_SET(samp2[levelID].ptr)) samp2[levelID].ptr[i] = 0.0;
                      if (IS_SET(request->var1.f3)) var13[levelID].ptr[i] = missval;
                      if (IS_SET(request->var2.h2)) var21[levelID].ptr[i] = missval;
                      if (IS_SET(request->var2.h3)) var23[levelID].ptr[i] = missval;
                    }
                  var12[levelID].nmiss = gridsize;
                  samp1[levelID].nmiss = gridsize;
                  if (IS_SET(request->var1.f3)) var13[levelID].nmiss = gridsize;
                  if (IS_SET(request->var2.h2)) var21[levelID].nmiss = gridsize;
                  if (IS_SET(request->var2.h3)) var23[levelID].nmiss = gridsize;
                }

              pstreamReadRecord(istreamID, field1.ptr, &nmiss);
              field1.nmiss = nmiss;
              field1.grid = var12[levelID].grid;
              field1.missval = var12[levelID].missval;

              farnum(&samp1[levelID], field1);

              if (IS_SET(request->var2.h2))
                {
                  arrayCopy(gridsize, field1.ptr, field2.ptr);
                  field2.nmiss = field1.nmiss;
                  field2.grid = field1.grid;
                  field2.missval = field1.missval;
                }

              if (IS_SET(request->var1.f1)) request->var1.f1(&field1, request->var1.f1arg);

              if (field1.nmiss > 0 || IS_SET(samp2[levelID].ptr))
                {
                  if (IS_NOT_SET(samp2[levelID].ptr))
                    {
                      samp2[levelID].ptr = (double *) Malloc(gridsize * sizeof(double));
                      for (size_t i = 0; i < gridsize; i++) samp2[levelID].ptr[i] = nsets;
                    }
                  for (size_t i = 0; i < gridsize; i++)
                    {
                      if (DBL_IS_EQUAL(field1.ptr[i], field1.missval)) continue;
                      samp2[levelID].ptr[i]++;
                    }
                }

              if (IS_NOT_EQUAL(request->var1.mulc, 0.0)) farcmul(&field1, request->var1.mulc);
              if (IS_NOT_EQUAL(request->var1.addc, 0.0)) farcadd(&field1, request->var1.addc);

              request->var1.f2(&var12[levelID], field1);

              if (IS_SET(request->var2.h2) || IS_SET(request->var2.h3))
                {
                  /* if h2 is null, use the output of f2 as input for h1 */
                  if (IS_NOT_SET(request->var2.h2))
                    {
                      arrayCopy(gridsize, var12[levelID].ptr, field2.ptr);
                      field2.nmiss = var12[levelID].nmiss;
                      field2.grid = var12[levelID].grid;
                      field2.missval = var12[levelID].missval;
                    }

                  if (IS_SET(request->var2.h1)) request->var2.h1(&field2, request->var2.h1arg);

                  if (IS_NOT_SET(request->var2.h2))
                    request->var2.h3(&var23[levelID], field2);
                  else
                    {
                      request->var2.h2(&var21[levelID], field2);
                      if (IS_SET(request->var2.h3)) request->var2.h3(&var23[levelID], var21[levelID]);
                    }
                }

              if (IS_SET(request->var1.f3)) request->var1.f3(&var13[levelID], var12[levelID]);
            }

          ovdate = ivdate;
          ovtime = ivtime;
          nsets++;
          itsID++;
        }

      if (nrecs == 0 && nsets == 0) break;

      if (request->var1.epilog == MEAN || request->var1.epilog == PERCENT_OF_TIME)
        for (levelID = 0; levelID < nlevels; levelID++)
          {
            if (IS_SET(request->var1.f3))
              var = &var13[levelID];
            else
              var = &var12[levelID];

            if (IS_NOT_SET(samp2[levelID].ptr))
              farcdiv(var, nsets);
            else
              fardiv(var, samp2[levelID]);

            if (request->var1.epilog == PERCENT_OF_TIME) farcmul(var, 100.0);
          }

      taxisDefVdate(otaxisID, ovdate);
      taxisDefVtime(otaxisID, ovtime);
      pstreamDefTimestep(ostreamID, otsID);

      if (otsID && vlistInqVarTimetype(ivlistID, FIRST_VAR_ID) == TIME_CONSTANT) continue;

      varID = 0;
      for (levelID = 0; levelID < nlevels; levelID++)
        {
          if (IS_SET(request->var1.f3))
            var = &var13[levelID];
          else
            var = &var12[levelID];

          farsel(var, samp1[levelID]);

          pstreamDefRecord(ostreamID, varID, levelID);
          pstreamWriteRecord(ostreamID, var->ptr, var->nmiss);
        }

      if (IS_SET(request->var2.h2) || IS_SET(request->var2.h3))
        {
          varID = 1;
          for (levelID = 0; levelID < nlevels; levelID++)
            {
              if (IS_SET(request->var2.h3))
                var = &var23[levelID];
              else
                var = &var21[levelID];

              farsel(var, samp1[levelID]);

              pstreamDefRecord(ostreamID, varID, levelID);
              pstreamWriteRecord(ostreamID, var->ptr, var->nmiss);
            }
        }

      if (nrecs == 0) break;
      otsID++;
    }

  for (levelID = 0; levelID < nlevels; levelID++)
    {
      Free(var12[levelID].ptr);
      Free(samp1[levelID].ptr);
      if (IS_SET(samp2[levelID].ptr)) Free(samp2[levelID].ptr);
    }
  Free(var12);
  Free(samp1);
  Free(samp2);

  if (IS_SET(var13))
    {
      for (levelID = 0; levelID < nlevels; levelID++) Free(var13[levelID].ptr);
      Free(var13);
    }
  if (IS_SET(var21))
    {
      for (levelID = 0; levelID < nlevels; levelID++) Free(var21[levelID].ptr);
      Free(var21);
    }
  if (IS_SET(var23))
    {
      for (levelID = 0; levelID < nlevels; levelID++) Free(var23[levelID].ptr);
      Free(var23);
    }

  if (IS_SET(field1.ptr)) Free(field1.ptr);
  if (IS_SET(field2.ptr)) Free(field2.ptr);

  if (IS_SET(recVarID)) Free(recVarID);
  if (IS_SET(recLevelID)) Free(recLevelID);

  pstreamClose(ostreamID);
  pstreamClose(istreamID);
}

void
eca2(const ECA_REQUEST_2 *request)
{
  const int operatorID = cdoOperatorID();

  size_t nmiss;
  int cmplen;
  char indate1[DATE_LEN + 1], indate2[DATE_LEN + 1];
  int64_t ivdate = 0, ovdate = 0;
  int ivtime = 0, ovtime = 0;
  int nrecs, nrecords;
  int gridID, zaxisID, varID, levelID;
  int itsID;
  int otsID;
  long nsets;
  int istreamID1, istreamID2, ostreamID;
  int ivlistID1, ivlistID2, ovlistID, itaxisID1, otaxisID;
  int nlevels;
  int *recVarID, *recLevelID;
  double missval1, missval2;
  Field *var14 = NULL, *samp1 = NULL, *samp2 = NULL, *samp3 = NULL, *total = NULL, *var15 = NULL, *var22 = NULL, *var;
  Field field1, field2;

  cmplen = DATE_LEN - cdoOperatorF2(operatorID);

  istreamID1 = cdoStreamOpenRead(cdoStreamName(0));
  istreamID2 = cdoStreamOpenRead(cdoStreamName(1));

  ivlistID1 = cdoStreamInqVlist(istreamID1);
  ivlistID2 = cdoStreamInqVlist(istreamID2);
  ovlistID = vlistCreate();

  vlistCompare(ivlistID1, ivlistID2, CMP_ALL);

  gridID = vlistInqVarGrid(ivlistID1, FIRST_VAR_ID);
  zaxisID = vlistInqVarZaxis(ivlistID1, FIRST_VAR_ID);
  missval1 = vlistInqVarMissval(ivlistID1, FIRST_VAR_ID);
  missval2 = vlistInqVarMissval(ivlistID2, FIRST_VAR_ID);

  varID = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARYING);

  vlistDefVarMissval(ovlistID, varID, missval1);

  if (IS_SET(request->var1.name)) vlistDefVarName(ovlistID, varID, request->var1.name);
  if (IS_SET(request->var1.longname)) vlistDefVarLongname(ovlistID, varID, request->var1.longname);
  if (IS_SET(request->var1.units)) vlistDefVarUnits(ovlistID, varID, request->var1.units);

  if (IS_SET(request->var2.h2))
    {
      varID = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARYING);

      vlistDefVarMissval(ovlistID, varID, missval1);

      if (IS_SET(request->var2.name)) vlistDefVarName(ovlistID, varID, request->var2.name);
      if (IS_SET(request->var2.longname)) vlistDefVarLongname(ovlistID, varID, request->var2.longname);
      if (IS_SET(request->var2.units)) vlistDefVarUnits(ovlistID, varID, request->var2.units);
    }

  if (cdoOperatorF2(operatorID) == 16) vlistDefNtsteps(ovlistID, 1);

  itaxisID1 = vlistInqTaxis(ivlistID1);
  otaxisID = taxisCreate(TAXIS_RELATIVE);
  taxisDefTunit(otaxisID, TUNIT_MINUTE);
  //  taxisDefCalendar(otaxisID, CALENDAR_PROLEPTIC);
  taxisDefCalendar(otaxisID, taxisInqCalendar(itaxisID1));
  taxisDefRdate(otaxisID, 19550101);
  taxisDefRtime(otaxisID, 0);
  vlistDefTaxis(ovlistID, otaxisID);

  ostreamID = cdoStreamOpenWrite(cdoStreamName(2), cdoFiletype());
  pstreamDefVlist(ostreamID, ovlistID);

  nrecords = vlistNrecs(ivlistID1);
  recVarID = (int *) Malloc(nrecords * sizeof(int));
  recLevelID = (int *) Malloc(nrecords * sizeof(int));

  size_t gridsize = gridInqSize(gridID);

  field_init(&field1);
  field_init(&field2);
  field1.ptr = (double *) Malloc(gridsize * sizeof(double));
  field2.ptr = (double *) Malloc(gridsize * sizeof(double));

  nlevels = zaxisInqSize(zaxisID);

  var14 = (Field *) Malloc(nlevels * sizeof(Field));
  samp1 = (Field *) Malloc(nlevels * sizeof(Field));
  samp2 = (Field *) Malloc(nlevels * sizeof(Field));
  samp3 = (Field *) Malloc(nlevels * sizeof(Field));

  if (request->var1.epilog == PERCENT_OF_TOTAL_AMOUNT) total = (Field *) Malloc(nlevels * sizeof(Field));
  if (IS_SET(request->var1.f5)) var15 = (Field *) Malloc(nlevels * sizeof(Field));
  if (IS_SET(request->var2.h2)) var22 = (Field *) Malloc(nlevels * sizeof(Field));

  for (levelID = 0; levelID < nlevels; levelID++)
    {
      field_init(&var14[levelID]);
      var14[levelID].grid = gridID;
      var14[levelID].nmiss = 0;
      var14[levelID].missval = missval1;
      var14[levelID].ptr = (double *) Malloc(gridsize * sizeof(double));

      field_init(&samp1[levelID]);
      samp1[levelID].grid = gridID;
      samp1[levelID].nmiss = 0;
      samp1[levelID].missval = missval1;
      samp1[levelID].ptr = (double *) Malloc(gridsize * sizeof(double));

      field_init(&samp2[levelID]);
      samp2[levelID].grid = gridID;
      samp2[levelID].nmiss = 0;
      samp2[levelID].missval = missval1;
      samp2[levelID].ptr = (double *) Malloc(gridsize * sizeof(double));

      field_init(&samp3[levelID]);
      samp3[levelID].grid = gridID;
      samp3[levelID].nmiss = 0;
      samp3[levelID].missval = missval1;
      samp3[levelID].ptr = NULL;

      if (request->var1.epilog == PERCENT_OF_TOTAL_AMOUNT)
        {
          field_init(&total[levelID]);
          total[levelID].grid = gridID;
          total[levelID].nmiss = 0;
          total[levelID].missval = missval1;
          total[levelID].ptr = (double *) Malloc(gridsize * sizeof(double));
        }
      if (IS_SET(request->var1.f5))
        {
          field_init(&var15[levelID]);
          var15[levelID].grid = gridID;
          var15[levelID].nmiss = 0;
          var15[levelID].missval = missval1;
          var15[levelID].ptr = (double *) Malloc(gridsize * sizeof(double));
        }
      if (IS_SET(request->var2.h2))
        {
          field_init(&var22[levelID]);
          var22[levelID].grid = gridID;
          var22[levelID].nmiss = 0;
          var22[levelID].missval = missval1;
          var22[levelID].ptr = (double *) Malloc(gridsize * sizeof(double));
        }
    }

  itsID = 0;
  otsID = 0;
  while (TRUE)
    {
      nsets = 0;
      while ((nrecs = cdoStreamInqTimestep(istreamID1, itsID)) > 0)
        {
          if (!cdoStreamInqTimestep(istreamID2, itsID)) cdoAbort("Input streams have different number of time steps!");

          ivdate = taxisInqVdate(itaxisID1);
          ivtime = taxisInqVtime(itaxisID1);

          if (nsets == 0) SET_DATE(indate2, ivdate, ivtime);
          SET_DATE(indate1, ivdate, ivtime);

          if (DATE_IS_NEQ(indate1, indate2, cmplen)) break;

          for (int recID = 0; recID < nrecs; recID++)
            {
              pstreamInqRecord(istreamID1, &varID, &levelID);
              pstreamInqRecord(istreamID2, &varID, &levelID);

              if (itsID == 0)
                {
                  recVarID[recID] = varID;
                  recLevelID[recID] = levelID;
                }
              if (varID != FIRST_VAR_ID) continue;

              if (nsets == 0)
                {
                  for (size_t i = 0; i < gridsize; i++)
                    {
                      var14[levelID].ptr[i] = missval1;
                      samp1[levelID].ptr[i] = missval1;
                      samp2[levelID].ptr[i] = missval1;
                      if (IS_SET(samp3[levelID].ptr)) samp3[levelID].ptr[i] = 0.0;
                      if (request->var1.epilog == PERCENT_OF_TOTAL_AMOUNT) total[levelID].ptr[i] = 0.0;
                      if (IS_SET(request->var1.f5)) var15[levelID].ptr[i] = missval1;
                      if (IS_SET(request->var2.h2)) var22[levelID].ptr[i] = missval1;
                    }
                  var14[levelID].nmiss = gridsize;
                  samp1[levelID].nmiss = gridsize;
                  samp2[levelID].nmiss = gridsize;
                  if (request->var1.epilog == PERCENT_OF_TOTAL_AMOUNT) total[levelID].nmiss = gridsize;
                  if (IS_SET(request->var1.f5)) var15[levelID].nmiss = gridsize;
                  if (IS_SET(request->var2.h2)) var22[levelID].nmiss = gridsize;
                }

              pstreamReadRecord(istreamID1, field1.ptr, &nmiss);
              field1.nmiss = nmiss;
              field1.grid = gridID;
              field1.missval = missval1;

              pstreamReadRecord(istreamID2, field2.ptr, &nmiss);
              field2.nmiss = nmiss;
              field2.grid = gridID;
              field2.missval = missval2;

              farnum(&samp1[levelID], field1);
              farnum(&samp2[levelID], field2);

              if (request->var1.epilog == PERCENT_OF_TOTAL_AMOUNT) farsum(&total[levelID], field1);

              if (IS_SET(request->var1.f1)) request->var1.f1(&field1, request->var1.f1arg);

              if (IS_SET(request->var1.f2)) request->var1.f2(&field2, request->var1.f2arg);

              if (field1.nmiss > 0 || IS_SET(samp3[levelID].ptr))
                {
                  if (IS_NOT_SET(samp3[levelID].ptr))
                    {
                      samp3[levelID].ptr = (double *) Malloc(gridsize * sizeof(double));
                      for (size_t i = 0; i < gridsize; i++) samp3[levelID].ptr[i] = nsets;
                    }
                  for (size_t i = 0; i < gridsize; i++)
                    {
                      if (DBL_IS_EQUAL(field1.ptr[i], field1.missval)) continue;
                      samp3[levelID].ptr[i]++;
                    }
                }

              request->var1.f3(&field1, field2);
              request->var1.f4(&var14[levelID], field1);

              if (IS_SET(request->var2.h2))
                {
                  arrayCopy(gridsize, var14[levelID].ptr, field2.ptr);
                  field2.nmiss = var14[levelID].nmiss;
                  field2.grid = var14[levelID].grid;
                  field2.missval = var14[levelID].missval;

                  if (IS_SET(request->var2.h1)) request->var2.h1(&field2, request->var2.h1arg);

                  request->var2.h2(&var22[levelID], field2);
                }

              if (IS_SET(request->var1.f5)) request->var1.f5(&var15[levelID], var14[levelID], request->var1.f5arg);
            }

          ovdate = ivdate;
          ovtime = ivtime;
          nsets++;
          itsID++;
        }

      if (nrecs == 0 && nsets == 0) break;

      if (request->var1.epilog == MEAN || request->var1.epilog == PERCENT_OF_TIME)
        for (levelID = 0; levelID < nlevels; levelID++)
          {
            if (IS_SET(request->var1.f5))
              var = &var15[levelID];
            else
              var = &var14[levelID];

            if (IS_NOT_SET(samp3[levelID].ptr))
              farcdiv(var, nsets);
            else
              fardiv(var, samp3[levelID]);

            if (request->var1.epilog == PERCENT_OF_TIME) farcmul(var, 100.0);
          }
      else if (request->var1.epilog == PERCENT_OF_TOTAL_AMOUNT)
        for (levelID = 0; levelID < nlevels; levelID++)
          {
            if (IS_SET(request->var1.f5))
              var = &var15[levelID];
            else
              var = &var14[levelID];

            fardiv(var, total[levelID]);
            farcmul(var, 100.0);
          }

      taxisDefVdate(otaxisID, ovdate);
      taxisDefVtime(otaxisID, ovtime);
      pstreamDefTimestep(ostreamID, otsID);

      if (otsID && vlistInqVarTimetype(ivlistID1, FIRST_VAR_ID) == TIME_CONSTANT) continue;

      varID = 0;
      for (levelID = 0; levelID < nlevels; levelID++)
        {
          if (IS_SET(request->var1.f5))
            var = &var15[levelID];
          else
            var = &var14[levelID];

          farsel(var, samp1[levelID]);
          farsel(var, samp2[levelID]);

          pstreamDefRecord(ostreamID, varID, levelID);
          pstreamWriteRecord(ostreamID, var->ptr, var->nmiss);
        }
      if (IS_SET(request->var2.h2))
        {
          varID = 1;
          for (levelID = 0; levelID < nlevels; levelID++)
            {
              var = &var22[levelID];

              farsel(var, samp1[levelID]);
              farsel(var, samp2[levelID]);

              pstreamDefRecord(ostreamID, varID, levelID);
              pstreamWriteRecord(ostreamID, var->ptr, var->nmiss);
            }
        }

      if (nrecs == 0) break;
      otsID++;
    }

  for (levelID = 0; levelID < nlevels; levelID++)
    {
      Free(var14[levelID].ptr);
      Free(samp1[levelID].ptr);
      Free(samp2[levelID].ptr);
      if (IS_SET(samp3[levelID].ptr)) Free(samp3[levelID].ptr);
    }
  Free(var14);
  Free(samp1);
  Free(samp2);
  Free(samp3);

  if (IS_SET(total))
    {
      for (levelID = 0; levelID < nlevels; levelID++) Free(total[levelID].ptr);
      Free(total);
    }
  if (IS_SET(var15))
    {
      for (levelID = 0; levelID < nlevels; levelID++) Free(var15[levelID].ptr);
      Free(var15);
    }
  if (IS_SET(var22))
    {
      for (levelID = 0; levelID < nlevels; levelID++) Free(var22[levelID].ptr);
      Free(var22);
    }

  if (IS_SET(field1.ptr)) Free(field1.ptr);
  if (IS_SET(field2.ptr)) Free(field2.ptr);

  if (IS_SET(recVarID)) Free(recVarID);
  if (IS_SET(recLevelID)) Free(recLevelID);

  pstreamClose(ostreamID);
  pstreamClose(istreamID2);
  pstreamClose(istreamID1);
}

void
eca3(const ECA_REQUEST_3 *request)
{
  const int operatorID = cdoOperatorID();

  size_t nmiss;
  int cmplen;
  char indate1[DATE_LEN + 1], indate2[DATE_LEN + 1];
  int64_t ivdate1 = 0, ivdate2 = 0;
  int ivtime1 = 0, ivtime2 = 0;
  int64_t ovdate = 0;
  int ovtime = 0;
  int nrecs, nrecords;
  int gridID, zaxisID, varID, levelID;
  int itsID;
  int otsID;
  long nsets;
  int istreamID1, istreamID2, ostreamID;
  int ivlistID1, ivlistID2, ovlistID, itaxisID1, itaxisID2, otaxisID;
  int nlevels;
  int *recVarID, *recLevelID;
  double missval;
  Field *var1 = NULL, *var2 = NULL;
  Field field1, field2;

  cmplen = DATE_LEN - cdoOperatorF2(operatorID);

  istreamID1 = cdoStreamOpenRead(cdoStreamName(0));
  istreamID2 = cdoStreamOpenRead(cdoStreamName(1));

  ivlistID1 = cdoStreamInqVlist(istreamID1);
  ivlistID2 = cdoStreamInqVlist(istreamID2);
  ovlistID = vlistCreate();

  vlistCompare(ivlistID1, ivlistID2, CMP_ALL);

  gridID = vlistInqVarGrid(ivlistID1, FIRST_VAR_ID);
  zaxisID = vlistInqVarZaxis(ivlistID1, FIRST_VAR_ID);
  missval = vlistInqVarMissval(ivlistID1, FIRST_VAR_ID);

  varID = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARYING);

  vlistDefVarMissval(ovlistID, varID, missval);

  if (IS_SET(request->name)) vlistDefVarName(ovlistID, varID, request->name);
  if (IS_SET(request->longname)) vlistDefVarLongname(ovlistID, varID, request->longname);
  if (IS_SET(request->units)) vlistDefVarUnits(ovlistID, varID, request->units);

  if (cdoOperatorF2(operatorID) == 16) vlistDefNtsteps(ovlistID, 1);

  itaxisID1 = vlistInqTaxis(ivlistID1);
  itaxisID2 = vlistInqTaxis(ivlistID2);
  otaxisID = taxisCreate(TAXIS_RELATIVE);
  taxisDefTunit(otaxisID, TUNIT_MINUTE);
  //  taxisDefCalendar(otaxisID, CALENDAR_PROLEPTIC);
  taxisDefCalendar(otaxisID, taxisInqCalendar(itaxisID1));
  taxisDefRdate(otaxisID, 19550101);
  taxisDefRtime(otaxisID, 0);
  vlistDefTaxis(ovlistID, otaxisID);

  ostreamID = cdoStreamOpenWrite(cdoStreamName(2), cdoFiletype());
  pstreamDefVlist(ostreamID, ovlistID);

  nrecords = vlistNrecs(ivlistID1);
  recVarID = (int *) Malloc(nrecords * sizeof(int));
  recLevelID = (int *) Malloc(nrecords * sizeof(int));

  size_t gridsize = gridInqSize(gridID);

  field_init(&field1);
  field_init(&field2);
  field1.ptr = (double *) Malloc(gridsize * sizeof(double));
  field2.ptr = (double *) Malloc(gridsize * sizeof(double));

  nlevels = zaxisInqSize(zaxisID);

  var1 = (Field *) Malloc(nlevels * sizeof(Field));
  var2 = (Field *) Malloc(nlevels * sizeof(Field));

  for (levelID = 0; levelID < nlevels; levelID++)
    {
      field_init(&var1[levelID]);
      var1[levelID].grid = gridID;
      var1[levelID].nmiss = 0;
      var1[levelID].missval = missval;
      var1[levelID].ptr = (double *) Malloc(gridsize * sizeof(double));

      field_init(&var2[levelID]);
      var2[levelID].grid = gridID;
      var2[levelID].nmiss = 0;
      var2[levelID].missval = missval;
      var2[levelID].ptr = (double *) Malloc(gridsize * sizeof(double));
    }

  itsID = 0;
  otsID = 0;
  while (TRUE)
    {
      nsets = 0;
      while ((nrecs = cdoStreamInqTimestep(istreamID1, itsID)) > 0)
        {
          if (!cdoStreamInqTimestep(istreamID2, itsID)) cdoAbort("Input streams have different number of time steps!");

          ivdate1 = taxisInqVdate(itaxisID1);
          ivdate2 = taxisInqVdate(itaxisID2);
          if (ivdate1 != ivdate2)
            cdoAbort("Input streams have different verification dates at time step %d!", itsID + 1);

          ivtime1 = taxisInqVtime(itaxisID1);
          ivtime2 = taxisInqVtime(itaxisID2);
          if (ivtime1 != ivtime2)
            cdoAbort("Input streams have different verification times at time step %d!", itsID + 1);

          if (nsets == 0) SET_DATE(indate2, ivdate1, ivtime1);
          SET_DATE(indate1, ivdate1, ivtime1);

          if (DATE_IS_NEQ(indate1, indate2, cmplen)) break;

          for (int recID = 0; recID < nrecs; recID++)
            {
              pstreamInqRecord(istreamID1, &varID, &levelID);
              pstreamInqRecord(istreamID2, &varID, &levelID);

              if (itsID == 0)
                {
                  recVarID[recID] = varID;
                  recLevelID[recID] = levelID;
                }
              if (varID != FIRST_VAR_ID) continue;

              if (nsets == 0)
                {
                  for (size_t i = 0; i < gridsize; i++)
                    {
                      var1[levelID].ptr[i] = missval;
                      var2[levelID].ptr[i] = missval;
                    }
                  var1[levelID].nmiss = gridsize;
                  var2[levelID].nmiss = gridsize;
                }

              pstreamReadRecord(istreamID1, field1.ptr, &nmiss);
              field1.nmiss = nmiss;
              field1.grid = var1[levelID].grid;
              field1.missval = var1[levelID].missval;

              pstreamReadRecord(istreamID2, field2.ptr, &nmiss);
              field2.nmiss = nmiss;
              field2.grid = var1[levelID].grid;
              field2.missval = var1[levelID].missval;

              request->f1(&var1[levelID], field1);
              request->f2(&var2[levelID], field2);
            }

          ovdate = ivdate1;
          ovtime = ivtime1;
          nsets++;
          itsID++;
        }

      if (nrecs == 0 && nsets == 0) break;

      for (levelID = 0; levelID < nlevels; levelID++) request->f3(&var1[levelID], var2[levelID]);

      taxisDefVdate(otaxisID, ovdate);
      taxisDefVtime(otaxisID, ovtime);
      pstreamDefTimestep(ostreamID, otsID);

      if (otsID && vlistInqVarTimetype(ivlistID1, FIRST_VAR_ID) == TIME_CONSTANT) continue;

      varID = 0;
      for (levelID = 0; levelID < nlevels; levelID++)
        {
          pstreamDefRecord(ostreamID, varID, levelID);
          pstreamWriteRecord(ostreamID, var1[levelID].ptr, var1[levelID].nmiss);
        }

      if (nrecs == 0) break;
      otsID++;
    }

  for (levelID = 0; levelID < nlevels; levelID++)
    {
      Free(var1[levelID].ptr);
      Free(var2[levelID].ptr);
    }
  Free(var1);
  Free(var2);

  if (IS_SET(field1.ptr)) Free(field1.ptr);
  if (IS_SET(field2.ptr)) Free(field2.ptr);

  if (IS_SET(recVarID)) Free(recVarID);
  if (IS_SET(recLevelID)) Free(recLevelID);

  pstreamClose(ostreamID);
  pstreamClose(istreamID2);
  pstreamClose(istreamID1);
}

void
eca4(const ECA_REQUEST_4 *request)
{
  const int operatorID = cdoOperatorID();

  size_t nmiss;
  int yearcnt = 0;
  int nrecs;
  int varID, levelID;
  int resetAtJan = FALSE, resetAtJul = FALSE;
  int isFirstYear = TRUE;
  int64_t ivdate = 0, ovdate = 0;
  int ivtime = 0, ovtime = 0;
  char indate1[DATE_LEN + 1], indate2[DATE_LEN + 1];
  Field *startDateWithHist[2], *endDateWithHist[2];

  int cmplen = DATE_LEN - cdoOperatorF2(operatorID);

  int istreamID1 = cdoStreamOpenRead(cdoStreamName(0));
  int istreamID2 = cdoStreamOpenRead(cdoStreamName(1));

  int ivlistID1 = cdoStreamInqVlist(istreamID1);
  int ivlistID2 = cdoStreamInqVlist(istreamID2);
  int ovlistID = vlistCreate();

  int gridID = vlistInqVarGrid(ivlistID1, FIRST_VAR_ID);
  if (gridInqSize(gridID) != gridInqSize(vlistInqVarGrid(ivlistID2, FIRST_VAR_ID)))
    cdoAbort("Grid sizes of the input fields do not match!");

  int zaxisID = vlistInqVarZaxis(ivlistID1, FIRST_VAR_ID);
  double missval = vlistInqVarMissval(ivlistID1, FIRST_VAR_ID);

  int ovarID1 = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARYING);

  vlistDefVarMissval(ovlistID, ovarID1, missval);

  if (IS_SET(request->name)) vlistDefVarName(ovlistID, ovarID1, request->name);
  if (IS_SET(request->longname)) vlistDefVarLongname(ovlistID, ovarID1, request->longname);
  if (IS_SET(request->units)) vlistDefVarUnits(ovlistID, ovarID1, request->units);

  int ovarID2 = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARYING);

  vlistDefVarMissval(ovlistID, ovarID2, missval);

  if (IS_SET(request->name2)) vlistDefVarName(ovlistID, ovarID2, request->name2);
  if (IS_SET(request->longname2)) vlistDefVarLongname(ovlistID, ovarID2, request->longname2);
  if (IS_SET(request->units2)) vlistDefVarUnits(ovlistID, ovarID2, request->units2);

  if (cdoOperatorF2(operatorID) == 16) vlistDefNtsteps(ovlistID, 1);

  int itaxisID = vlistInqTaxis(ivlistID1);
  int otaxisID = taxisCreate(TAXIS_RELATIVE);
  taxisDefTunit(otaxisID, TUNIT_MINUTE);
  //  taxisDefCalendar(otaxisID, CALENDAR_PROLEPTIC);
  taxisDefCalendar(otaxisID, taxisInqCalendar(itaxisID));
  taxisDefRdate(otaxisID, 19550101);
  taxisDefRtime(otaxisID, 0);
  vlistDefTaxis(ovlistID, otaxisID);

  int ostreamID = cdoStreamOpenWrite(cdoStreamName(2), cdoFiletype());
  pstreamDefVlist(ostreamID, ovlistID);

  int nrecords = vlistNrecs(ivlistID1);
  int *recVarID = (int *) Malloc(nrecords * sizeof(int));
  int *recLevelID = (int *) Malloc(nrecords * sizeof(int));

  bool lyvals = true;
  int gridtype = gridInqType(gridID);
  if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_PROJECTION)
    {
      gridID = gridToCurvilinear(gridID, 1);
    }
  else if (gridtype == GRID_GME)
    {
      gridID = gridToUnstructured(gridID, 1);
    }
  else
    {
      lyvals = false;
    }
  size_t gridsize = gridInqSize(gridID);
  /* for later check on northern\southern hemisphere */
  double *yvals = (double *) Malloc(gridsize * sizeof(double));
  if (lyvals)
    {
      gridInqYvals(gridID, yvals);
    }
  else
    {
      for (size_t i = 0; i < gridsize; ++i) yvals[i] = 20;  // Northern hemisphere
    }

  /* Two fields are needed because of the definition of gsl for northern and
   * southern hemisphere                                                      */
  Field fieldGt, fieldLt;
  field_init(&fieldGt);
  field_init(&fieldLt);
  fieldGt.ptr = (double *) Malloc(gridsize * sizeof(double));
  fieldLt.ptr = (double *) Malloc(gridsize * sizeof(double));

  /* field for the land-water-distribution */
  Field mask;
  field_init(&mask);
  mask.ptr = (double *) Malloc(gridsize * sizeof(double));

  int nlevels = zaxisInqSize(zaxisID);

  Field *startCount = (Field *) Malloc(nlevels * sizeof(Field));
  Field *endCount = (Field *) Malloc(nlevels * sizeof(Field));
  Field *gslDuration = (Field *) Malloc(nlevels * sizeof(Field));
  Field *gslFirstDay = (Field *) Malloc(nlevels * sizeof(Field));

  /* because of the different definitions for northern and southern hemisphere,
   * the values of the last year have to be present THE LAST YEAR HAS THE INDEX
   * 1 */
  for (int h = 0; h < 2; h++)
    {
      startDateWithHist[h] = (Field *) Malloc(nlevels * sizeof(Field));
      endDateWithHist[h] = (Field *) Malloc(nlevels * sizeof(Field));
    }

  for (levelID = 0; levelID < nlevels; levelID++)
    {
      field_init(&startCount[levelID]);
      startCount[levelID].grid = gridID;
      startCount[levelID].size = gridsize;
      startCount[levelID].nmiss = 0;
      startCount[levelID].missval = missval;
      startCount[levelID].ptr = (double *) Calloc(gridsize, sizeof(double));

      field_init(&endCount[levelID]);
      endCount[levelID].grid = gridID;
      endCount[levelID].size = gridsize;
      endCount[levelID].nmiss = 0;
      endCount[levelID].missval = missval;
      endCount[levelID].ptr = (double *) Calloc(gridsize, sizeof(double));

      field_init(&gslDuration[levelID]);
      gslDuration[levelID].grid = gridID;
      gslDuration[levelID].size = gridsize;
      gslDuration[levelID].nmiss = 0;
      gslDuration[levelID].missval = missval;
      gslDuration[levelID].ptr = (double *) Malloc(gridsize * sizeof(double));

      field_init(&gslFirstDay[levelID]);
      gslFirstDay[levelID].grid = gridID;
      gslFirstDay[levelID].size = gridsize;
      gslFirstDay[levelID].nmiss = 0;
      gslFirstDay[levelID].missval = missval;
      gslFirstDay[levelID].ptr = (double *) Malloc(gridsize * sizeof(double));

      for (int h = 0; h < 2; h++)
        {
          field_init(&startDateWithHist[h][levelID]);
          startDateWithHist[h][levelID].grid = gridID;
          startDateWithHist[h][levelID].size = gridsize;
          startDateWithHist[h][levelID].nmiss = 0;
          startDateWithHist[h][levelID].missval = missval;
          startDateWithHist[h][levelID].ptr = (double *) Malloc(gridsize * sizeof(double));

          field_init(&endDateWithHist[h][levelID]);
          endDateWithHist[h][levelID].grid = gridID;
          endDateWithHist[h][levelID].size = gridsize;
          endDateWithHist[h][levelID].nmiss = 0;
          endDateWithHist[h][levelID].missval = missval;
          endDateWithHist[h][levelID].ptr = (double *) Malloc(gridsize * sizeof(double));
        }
    }

  int itsID = 0;
  int otsID = 0;

  if (cdoStreamInqTimestep(istreamID2, itsID))
    {
      pstreamInqRecord(istreamID2, &varID, &levelID);
      pstreamReadRecord(istreamID2, mask.ptr, &nmiss);
      mask.nmiss = nmiss;
      mask.grid = gridID;
      mask.missval = vlistInqVarMissval(ivlistID2, 0);

      request->s3(&mask, request->s3arg);
    }
  else
    cdoAbort("Could not read land-water mask!");

  while (TRUE)
    {
      long nsets = 0;
      while ((nrecs = cdoStreamInqTimestep(istreamID1, itsID)) > 0)
        {
          ivdate = taxisInqVdate(itaxisID);
          ivtime = taxisInqVtime(itaxisID);

          int month = (ivdate % 10000) / 100;
          if (month < 1 || month > 12) cdoAbort("month %d out of range!", month);

          if (nsets == 0) SET_DATE(indate2, ivdate, ivtime);
          SET_DATE(indate1, ivdate, ivtime);

          if (DATE_IS_NEQ(indate1, indate2, cmplen))
            {
              resetAtJan = FALSE;
              resetAtJul = FALSE;
              break;
            }

          for (int recID = 0; recID < nrecs; recID++)
            {
              pstreamInqRecord(istreamID1, &varID, &levelID);

              if (itsID == 0)
                {
                  recVarID[recID] = varID;
                  recLevelID[recID] = levelID;
                }
              if (varID != FIRST_VAR_ID) continue;

              if (nsets == 0)
                {
                  for (size_t i = 0; i < gridsize; i++)
                    {
                      gslDuration[levelID].ptr[i] = missval;
                      gslFirstDay[levelID].ptr[i] = missval;
                      /* reinitialize the current year */
                      startDateWithHist[0][levelID].ptr[i] = missval;
                      endDateWithHist[0][levelID].ptr[i] = missval;
                    }
                  gslDuration[levelID].nmiss = 0;
                  gslFirstDay[levelID].nmiss = 0;
                  /* reinitialize the current year */
                  startDateWithHist[0][levelID].nmiss = gridsize;
                  endDateWithHist[0][levelID].nmiss = gridsize;
                }
              /* init the history ONCE */
              if (0 == itsID)
                {
                  for (size_t i = 0; i < gridsize; i++)
                    {
                      startDateWithHist[1][levelID].ptr[i] = missval;
                      endDateWithHist[1][levelID].ptr[i] = missval;
                    }
                  startDateWithHist[1][levelID].nmiss = gridsize;
                  endDateWithHist[1][levelID].nmiss = gridsize;
                }

              pstreamReadRecord(istreamID1, fieldGt.ptr, &nmiss);
              fieldGt.nmiss = nmiss;
              arrayCopy(gridsize, fieldGt.ptr, fieldLt.ptr);
              fieldLt.nmiss = fieldGt.nmiss;
              fieldGt.grid = startCount[levelID].grid;
              fieldGt.missval = startCount[levelID].missval;
              fieldLt.grid = startCount[levelID].grid;
              fieldLt.missval = startCount[levelID].missval;

              /* Reinitialization of (start|end)Count variables has to be done
               * different for norther and southern hemisphere */
              if (1 == month && !resetAtJan)
                {
                  /* reset northern startCount */
                  for (size_t i = 0; i < gridsize; i++)
                    {
                      if (yvals[i] >= 0.0)
                        if (!DBL_IS_EQUAL(startCount[levelID].ptr[i], missval))
                          {
                            startCount[levelID].ptr[i] = missval;
                            startCount[levelID].nmiss++;
                          }
                    }
                  /* reset southern endCount */
                  for (size_t i = 0; i < gridsize; i++)
                    {
                      if (yvals[i] < 0.0)
                        if (!DBL_IS_EQUAL(endCount[levelID].ptr[i], missval))
                          {
                            endCount[levelID].ptr[i] = missval;
                            endCount[levelID].nmiss++;
                          }
                    }

                  resetAtJan = TRUE;
                }
              if (7 == month && !resetAtJul)
                {
#ifdef _OPENMP
#pragma omp sections
#endif
                  {
#ifdef _OPENMP
#pragma omp section
#endif
                    {
                      /* reset northern endCount  */
                      for (size_t i = 0; i < gridsize; i++)
                        {
                          if (yvals[i] >= 0.0)
                            {
                              if (!DBL_IS_EQUAL(endCount[levelID].ptr[i], missval))
                                {
                                  endCount[levelID].ptr[i] = missval;
                                  endCount[levelID].nmiss++;
                                }
                            }
                        }
                    }
#ifdef _OPENMP
#pragma omp section
#endif
                    {
                      /* reset southern startCount */
                      for (size_t i = 0; i < gridsize; i++)
                        {
                          if (yvals[i] < 0.0)
                            {
                              if (!DBL_IS_EQUAL(startCount[levelID].ptr[i], missval))
                                {
                                  startCount[levelID].ptr[i] = missval;
                                  startCount[levelID].nmiss++;
                                }
                            }
                        }
                    }
                  }
                  resetAtJul = TRUE;
                }

/* count the day with temperature larger/smaller than the
 * given limit */
#ifdef _OPENMP
#pragma omp sections
#endif
              {
#ifdef _OPENMP
#pragma omp section
#endif
                {
                  farsel(&fieldGt, mask);
                  request->s1(&fieldGt, request->s1arg);
                  farnum2(&startCount[levelID], fieldGt);
                }
#ifdef _OPENMP
#pragma omp section
#endif
                {
                  farsel(&fieldLt, mask);
                  request->s2(&fieldLt, request->s1arg);
                  farnum2(&endCount[levelID], fieldLt);
                }
              }

              if (month < 7)
                {
                  for (size_t i = 0; i < gridsize; i++) /* dictinct between northern and southern sphere */
                    /* start with south */
                    if (yvals[i] < 0)
                      {
                        /* south: periods can also start in the first half of
                         * the year, but this date has already gone into the
                         * history */
                        if (DBL_IS_EQUAL(startDateWithHist[1][levelID].ptr[i], missval)
                            && IS_EQUAL(startCount[levelID].ptr[i], request->consecutiveDays))
                          {
                            startDateWithHist[1][levelID].ptr[i] = ivdate;
                            /* reset the endCount, because we are only
                             * interessted in the end of the eriod, if a start
                             * was found */
                            endCount[levelID].ptr[i] = missval;
                            endDateWithHist[0][levelID].ptr[i] = missval;
                          }
                        if (DBL_IS_EQUAL(endDateWithHist[0][levelID].ptr[i], missval)
                            && IS_EQUAL(endCount[levelID].ptr[i], request->consecutiveDays))
                          {
                            endDateWithHist[0][levelID].ptr[i] = ivdate;
                          }
                      }
                    else
                      {
                        if (DBL_IS_EQUAL(startDateWithHist[0][levelID].ptr[i], missval)
                            && IS_EQUAL(startCount[levelID].ptr[i], request->consecutiveDays))
                          {
                            startDateWithHist[0][levelID].ptr[i] = ivdate;
                          }
                      }
                }
              else
                {
                  for (size_t i = 0; i < gridsize; i++)
                    {
                      if (yvals[i] < 0)
                        {
                          if (DBL_IS_EQUAL(startDateWithHist[0][levelID].ptr[i], missval)
                              && IS_EQUAL(startCount[levelID].ptr[i], request->consecutiveDays))
                            {
                              startDateWithHist[0][levelID].ptr[i] = ivdate;
                            }
                        }
                      else
                        {
                          /* north: periods can also start in the second half of
                           * the year */
                          if (DBL_IS_EQUAL(startDateWithHist[0][levelID].ptr[i], missval)
                              && IS_EQUAL(startCount[levelID].ptr[i], request->consecutiveDays))
                            {
                              startDateWithHist[0][levelID].ptr[i] = ivdate;
                              /* reset the endCount, because we are only
                               * interessted in the end of the eriod, if a start
                               * was found */
                              endCount[levelID].ptr[i] = missval;
                              endDateWithHist[0][levelID].ptr[i] = missval;
                            }
                          if (DBL_IS_EQUAL(endDateWithHist[0][levelID].ptr[i], missval)
                              && IS_EQUAL(endCount[levelID].ptr[i], request->consecutiveDays))
                            {
                              endDateWithHist[0][levelID].ptr[i] = ivdate;
                            }
                        }
                    }
                }
              /* update nmiss for saving data in GRIB */
              fldunm(&startCount[levelID]);
              fldunm(&endCount[levelID]);
              fldunm(&startDateWithHist[1][levelID]);
              fldunm(&startDateWithHist[0][levelID]);
              fldunm(&endDateWithHist[1][levelID]);
              fldunm(&endDateWithHist[0][levelID]);
            }

          ovdate = ivdate;
          ovtime = ivtime;
          nsets++;
          itsID++;
        }

      if (nrecs == 0 && nsets == 0) break;

      adjustEndDate(nlevels, gridsize, yvals, missval, ovdate, startDateWithHist, endDateWithHist);

      /*  compute and write GSL for the previous year
       *  AND
       *  write the current start/end dates into the history
       *
       *  this is the default action if more than a year is available */
      if (yearcnt != 0)
        {
          computeGsl(nlevels, gridsize, yvals, missval, startDateWithHist, endDateWithHist, gslDuration, gslFirstDay, FALSE);

          /* values of the privous year */
          {
            writeGslStream(ostreamID, otaxisID, otsID, ovarID1, ovarID2, ivlistID1, FIRST_VAR_ID, gslDuration, gslFirstDay,
                           cdiEncodeDate(ovdate / 10000 - 1, 12, 31), ovtime, nlevels);
            otsID++;
          }
        }

      /*  if there is a previous year */
      if (ovdate != ivdate)
        {
          /*  if the first year of data was processed, the history has to
           *  be checked befor it get's updated. This is necessary, if a
           *  growing period on the southern hemisphere was found. Otherwise,
           *  it would get overwritten. */
          if (isFirstYear)
            {
              /*  Check for non missing values, i.e. is there any data for the
               *  previous year? */
              if (fldhvs(startDateWithHist[1], nlevels))
                {
                  computeGsl(nlevels, gridsize, yvals, missval, startDateWithHist, endDateWithHist, gslDuration, gslFirstDay,
                             FALSE);
                  {
                    writeGslStream(ostreamID, otaxisID, otsID, ovarID1, ovarID2, ivlistID1, FIRST_VAR_ID, gslDuration, gslFirstDay,
                                   cdiEncodeDate(ovdate / 10000 - 1, 12, 31), ovtime, nlevels);
                    otsID++;
                  }
                }
              isFirstYear = FALSE;
            }
#ifdef _OPENMP
#pragma omp sections
#endif
          {
            updateHist(startDateWithHist, nlevels, gridsize, yvals, FALSE);
#ifdef _OPENMP
#pragma omp section
#endif
            updateHist(endDateWithHist, nlevels, gridsize, yvals, TRUE);
          }
        }
      else /* process the current year, this only happens, if the last timestep
              is reached OR if data for only one year is present */
        {
          computeGsl(nlevels, gridsize, yvals, missval, startDateWithHist, endDateWithHist, gslDuration, gslFirstDay, TRUE);
          {
            writeGslStream(ostreamID, otaxisID, otsID, ovarID1, ovarID2, ivlistID1, FIRST_VAR_ID, gslDuration, gslFirstDay, ovdate,
                           ovtime, nlevels);
            otsID++;
          }
        }
      yearcnt++;

      if (nrecs == 0) break;
    }

  for (levelID = 0; levelID < nlevels; levelID++)
    {
      Free(startCount[levelID].ptr);
      Free(endCount[levelID].ptr);
      Free(gslDuration[levelID].ptr);
      Free(gslFirstDay[levelID].ptr);
      for (int h = 0; h < 2; h++)
        {
          Free(startDateWithHist[h][levelID].ptr);
          Free(endDateWithHist[h][levelID].ptr);
        }
    }
  for (int h = 0; h < 2; h++)
    {
      Free(startDateWithHist[h]);
      Free(endDateWithHist[h]);
    }
  Free(startCount);
  Free(endCount);
  Free(gslDuration);
  Free(gslFirstDay);

  if (IS_SET(fieldGt.ptr)) Free(fieldGt.ptr);
  if (IS_SET(fieldLt.ptr)) Free(fieldLt.ptr);

  if (IS_SET(recVarID)) Free(recVarID);
  if (IS_SET(recLevelID)) Free(recLevelID);

  pstreamClose(ostreamID);
  pstreamClose(istreamID1);
}

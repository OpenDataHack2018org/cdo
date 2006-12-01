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

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "dtypes.h"
#include "ecacore.h"
#include "ecautil.h"

#define FIRST_VAR_ID 0

#define IS_NOT_SET(x) (x == NULL)
#define IS_SET(x)     (x != NULL)


void eca1(const ECA_REQUEST_1 *request)
{
  static const char func[] = "eca1";
  const int operatorID = cdoOperatorID();
  
  INT64 intvdat;
  INT64 indate1 = 0, indate2;
  int gridsize;
  int ivdate = 0, ivtime = 0;
  int ovdate = 0, ovtime = 0;
  int nrecs, nrecords;
  int gridID, zaxisID, varID, levelID, recID;
  int itsID;
  int otsID;
  long nsets;
  int i;
  int istreamID, ostreamID;
  int ivlistID, ovlistID, itaxisID, otaxisID;
  int nlevels;
  int *recVarID, *recLevelID;
  double missval;
  FIELD *var12 = NULL, *samp1 = NULL, *var13 = NULL, *var21 = NULL, *var23 = NULL, *var;
  FIELD field1, field2;
  
  intvdat = (INT64) pow(10.0, (double) cdoOperatorIntval(operatorID));

  istreamID = streamOpenRead(cdoStreamName(0));
  if ( istreamID < 0 ) cdiError(istreamID, "Open failed on %s", cdoStreamName(0));

  ivlistID = streamInqVlist(istreamID);
  ovlistID = vlistCreate();
  
  gridID  = vlistInqVarGrid(ivlistID, FIRST_VAR_ID);
  zaxisID = vlistInqVarZaxis(ivlistID, FIRST_VAR_ID);
  varID   = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARIABLE);
  
  if ( IS_SET(request->var1.name) )
    vlistDefVarName(ovlistID, varID, request->var1.name);
  if ( IS_SET(request->var1.longname) ) 
    vlistDefVarLongname(ovlistID, varID, request->var1.longname);
  if ( IS_SET(request->var1.units) ) 
    vlistDefVarUnits(ovlistID, varID, request->var1.units);

  if ( IS_SET(request->var2.h2) || IS_SET(request->var2.h3) )
    {
      varID = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARIABLE);
  
      if ( IS_SET(request->var2.name) ) 
        vlistDefVarName(ovlistID, varID, request->var2.name);
      if ( IS_SET(request->var2.longname) ) 
        vlistDefVarLongname(ovlistID, varID, request->var2.longname);
    }
    
  if ( cdoOperatorIntval(operatorID) == 17 ) vlistDefNtsteps(ovlistID, 1);

  itaxisID = vlistInqTaxis(ivlistID);
  otaxisID = taxisCreate(TAXIS_ABSOLUTE);
  vlistDefTaxis(ovlistID, otaxisID);

  ostreamID = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( ostreamID < 0 ) cdiError(ostreamID, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(ostreamID, ovlistID);

  nrecords   = vlistNrecs(ivlistID);
  recVarID   = (int *) malloc(nrecords*sizeof(int));
  recLevelID = (int *) malloc(nrecords*sizeof(int));

  gridsize = gridInqSize(gridID);

  field1.ptr = (double *) malloc(gridsize*sizeof(double));
  if ( IS_SET(request->var2.h2) || IS_SET(request->var2.h3) ) 
    field2.ptr = (double *) malloc(gridsize*sizeof(double));
  else
    field2.ptr = NULL;

  nlevels = zaxisInqSize(zaxisID);
  missval = vlistInqVarMissval(ivlistID, FIRST_VAR_ID);

  var12 = (FIELD *) malloc(nlevels*sizeof(FIELD));
  samp1 = (FIELD *) malloc(nlevels*sizeof(FIELD));
  if ( IS_SET(request->var1.f3) ) 
    var13 = (FIELD *) malloc(nlevels*sizeof(FIELD));
    
  if ( IS_SET(request->var2.h2) ) 
    var21 = (FIELD *) malloc(nlevels*sizeof(FIELD));
  if ( IS_SET(request->var2.h3) ) 
    var23 = (FIELD *) malloc(nlevels*sizeof(FIELD));
      
  for ( levelID = 0; levelID < nlevels; levelID++ )
    {
      var12[levelID].grid    = gridID;
      var12[levelID].nmiss   = 0;
      var12[levelID].missval = missval;
      var12[levelID].ptr     = (double *) malloc(gridsize*sizeof(double));

      samp1[levelID].grid    = gridID;
      samp1[levelID].nmiss   = 0;
      samp1[levelID].missval = missval;
      samp1[levelID].ptr     = NULL;

      if ( IS_SET(request->var1.f3) )
        {
          var13[levelID].grid    = gridID;
          var13[levelID].nmiss   = 0;
          var13[levelID].missval = missval;
          var13[levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
        }
      if ( IS_SET(request->var2.h2) )
        {
          var21[levelID].grid    = gridID;
          var21[levelID].nmiss   = 0;
          var21[levelID].missval = missval;
          var21[levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
        }
      if ( IS_SET(request->var2.h3) )
        {
          var23[levelID].grid    = gridID;
          var23[levelID].nmiss   = 0;
          var23[levelID].missval = missval;
          var23[levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
        }
    }

  indate2 = 0;
  itsID   = 0;
  otsID   = 0;
  while ( TRUE )
    {
      nsets = 0;
      while ( (nrecs = streamInqTimestep(istreamID, itsID)) > 0 )
        {
          ivdate = taxisInqVdate(itaxisID);
          ivtime = taxisInqVtime(itaxisID);

          if ( nsets == 0 ) indate1 = (INT64)ivdate*10000 + ivtime;
          indate2 = (INT64)ivdate*10000 + ivtime;

          if ( indate2/intvdat != indate1/intvdat ) break;

          for ( recID = 0; recID < nrecs; recID++ )
            {
              streamInqRecord(istreamID, &varID, &levelID);

              if ( itsID == 0 )
                {
                  recVarID[recID]   = varID;
                  recLevelID[recID] = levelID;
                }
              if ( varID != FIRST_VAR_ID ) continue;

              if ( nsets == 0 )
                {
                  for ( i = 0; i < gridsize; i++ )
                    {
                      var12[levelID].ptr[i] = missval;
                      if ( IS_SET(samp1[levelID].ptr) )
                        samp1[levelID].ptr[i] = 0.0;
                      if ( IS_SET(request->var1.f3) ) 
                        var13[levelID].ptr[i] = missval;
                      if ( IS_SET(request->var2.h2) )
                        var21[levelID].ptr[i] = missval;
                      if ( IS_SET(request->var2.h3) )
                        var23[levelID].ptr[i] = missval;
                    }
                  var12[levelID].nmiss = gridsize;
                  if ( IS_SET(request->var1.f3) )
                    var13[levelID].nmiss = gridsize;
                  if ( IS_SET(request->var2.h2) )
                    var21[levelID].nmiss = gridsize; 
                  if ( IS_SET(request->var2.h3) )
                    var23[levelID].nmiss = gridsize; 
                }

              streamReadRecord(istreamID, field1.ptr, &field1.nmiss);
              field1.grid    = var12[levelID].grid;
              field1.missval = var12[levelID].missval;

              if ( IS_SET(request->var2.h2) )
                {
                  memcpy(field2.ptr, field1.ptr, gridsize*sizeof(double));
                  field2.nmiss   = field1.nmiss;
                  field2.grid    = field1.grid;
                  field2.missval = field1.missval;
                }
		
              if ( IS_SET(request->var1.f1) ) 
                request->var1.f1(&field1, request->var1.f1arg);
              
              if ( field1.nmiss > 0 || IS_SET(samp1[levelID].ptr) )
                {
                  if ( IS_NOT_SET(samp1[levelID].ptr) )
                    {
                      samp1[levelID].ptr = (double *) malloc(gridsize*sizeof(double));
                      for ( i = 0; i < gridsize; i++ )
                        samp1[levelID].ptr[i] = nsets;
                    }
                  for ( i = 0; i < gridsize; i++ )
                    {
                      if ( DBL_IS_EQUAL(field1.ptr[i], field1.missval) )
                        continue;
                      samp1[levelID].ptr[i]++;
                    }
                }
              
              if ( !DBL_IS_EQUAL(request->var1.mulc, 0.0) )
                farcmul(&field1, request->var1.mulc);
              if ( !DBL_IS_EQUAL(request->var1.addc, 0.0) )
                farcadd(&field1, request->var1.addc);
                
              request->var1.f2(&var12[levelID], field1);
              
              if ( IS_SET(request->var2.h2) || IS_SET(request->var2.h3) )
                { 
                  /* if h2 is null, use the output of f2 as input for h1 */
                  if ( IS_NOT_SET(request->var2.h2) )
                    {
                      memcpy(field2.ptr, var12[levelID].ptr, gridsize*sizeof(double));
                      field2.nmiss   = var12[levelID].nmiss;
                      field2.grid    = var12[levelID].grid;
                      field2.missval = var12[levelID].missval;
                    }
                  
                  if ( IS_SET(request->var2.h1) ) 
                    request->var2.h1(&field2, request->var2.h1arg);
                        
                  if ( IS_NOT_SET(request->var2.h2) )    
                    request->var2.h3(&var23[levelID], field2);
                  else
                    {
                      request->var2.h2(&var21[levelID], field2);
                      if ( IS_SET(request->var2.h3) )
                        request->var2.h3(&var23[levelID], var21[levelID]);
                    }
                }

              if ( IS_SET(request->var1.f3) )
                request->var1.f3(&var13[levelID], var12[levelID]);
            }

          ovdate = ivdate;
          ovtime = ivtime;
          nsets++;
          itsID++;
        }

      if ( nrecs == 0 && nsets == 0 ) break;
      
      if ( request->var1.epilog == MEAN || request->var1.epilog == PERCENT_OF_TIME )
        for ( levelID = 0; levelID < nlevels; levelID++ )
          {
            if ( IS_SET(request->var1.f3) )
              var = &var13[levelID];
            else 
              var = &var12[levelID];
              
            if ( IS_NOT_SET(samp1[levelID].ptr) )
              farcdiv(var, nsets);
            else
              fardiv(var, samp1[levelID]);
              
            if ( request->var1.epilog == PERCENT_OF_TIME )
              farcmul(var, 100.0);
          }

      taxisDefVdate(otaxisID, ovdate);
      taxisDefVtime(otaxisID, ovtime);
      streamDefTimestep(ostreamID, otsID++);

      if ( otsID == 1 || vlistInqVarTime(ivlistID, FIRST_VAR_ID) == TIME_VARIABLE )
        {
          varID = 0;
          for ( levelID = 0; levelID < nlevels; levelID++ )
            {
              if ( IS_SET(request->var1.f3) )
                var = &var13[levelID];
              else 
                var = &var12[levelID];
               
              streamDefRecord(ostreamID, varID, levelID);
              streamWriteRecord(ostreamID, var->ptr, var->nmiss);
            }
          if ( IS_SET(request->var2.h2) || IS_SET(request->var2.h3) )
            {
              varID = 1;
              for ( levelID = 0; levelID < nlevels; levelID++ )
                {
                  if ( IS_SET(request->var2.h3) )
                    var = &var23[levelID];
                  else 
                    var = &var21[levelID];
                  
                  streamDefRecord(ostreamID, varID, levelID);
                  streamWriteRecord(ostreamID, var->ptr, var->nmiss);
               }
            }
        }

      if ( nrecs == 0 ) break;
    }

  for ( levelID = 0; levelID < nlevels; levelID++ )
    {
      free(var12[levelID].ptr);
      if ( IS_SET(samp1[levelID].ptr) ) free(samp1[levelID].ptr);
    }  
  free(var12);
  free(samp1);
  
  if ( IS_SET(var13) ) 
    {
      for ( levelID = 0; levelID < nlevels; levelID++ )
        free(var13[levelID].ptr);
      free(var13);
    }
  if ( IS_SET(var21) ) 
    {
      for ( levelID = 0; levelID < nlevels; levelID++ )
        free(var21[levelID].ptr);
      free(var21);
    }
  if ( IS_SET(var23) ) 
    {
      for ( levelID = 0; levelID < nlevels; levelID++ )
        free(var23[levelID].ptr);
      free(var23);
    }
    
  if ( IS_SET(field1.ptr) ) free(field1.ptr);
  if ( IS_SET(field2.ptr) ) free(field2.ptr);

  if ( IS_SET(recVarID) )   free(recVarID);
  if ( IS_SET(recLevelID) ) free(recLevelID);

  streamClose(ostreamID);
  streamClose(istreamID);
}


void eca2(const ECA_REQUEST_2 *request)
{
  static const char func[] = "eca2";
  const int operatorID = cdoOperatorID();
  
  INT64 intvdat;
  INT64 indate1 = 0, indate2;
  int gridsize;
  int ivdate1 = 0, ivtime1 = 0;
  int ivdate2 = 0, ivtime2 = 0;
  int ovdate = 0, ovtime = 0;
  int nrecs, nrecords;
  int gridID, zaxisID, varID, levelID, recID;
  int itsID;
  int otsID;
  long nsets;
  int i;
  int istreamID1, istreamID2, ostreamID;
  int ivlistID1, ivlistID2, ovlistID, itaxisID1, itaxisID2, otaxisID;
  int nlevels;
  int *recVarID, *recLevelID;
  double missval;
  FIELD *var14 = NULL, *samp1 = NULL, *total = NULL, *var15 = NULL, *var22 = NULL, *var;
  FIELD field1, field2;
  
  intvdat = (INT64) pow(10.0, (double) cdoOperatorIntval(operatorID));

  istreamID1 = streamOpenRead(cdoStreamName(0));
  if ( istreamID1 < 0 ) cdiError(istreamID1, "Open failed on %s", cdoStreamName(0));
  istreamID2 = streamOpenRead(cdoStreamName(1));
  if ( istreamID2 < 0 ) cdiError(istreamID2, "Open failed on %s", cdoStreamName(1));

  ivlistID1 = streamInqVlist(istreamID1);
  ivlistID2 = streamInqVlist(istreamID2);
  ovlistID  = vlistCreate();
  
  vlistCompare(ivlistID1, ivlistID2, func_hrd);
  vlistCompare(ivlistID1, ivlistID2, func_sft);
  
  gridID  = vlistInqVarGrid(ivlistID1, FIRST_VAR_ID);
  zaxisID = vlistInqVarZaxis(ivlistID1, FIRST_VAR_ID);
  varID   = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARIABLE);
  
  if ( IS_SET(request->var1.name) ) 
    vlistDefVarName(ovlistID, varID, request->var1.name);
  if ( IS_SET(request->var1.longname) ) 
    vlistDefVarLongname(ovlistID, varID, request->var1.longname);
  if ( IS_SET(request->var1.units) ) 
    vlistDefVarUnits(ovlistID, varID, request->var1.units);

  if ( IS_SET(request->var2.h2) )
    {
      varID = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARIABLE);
  
      if ( IS_SET(request->var2.name) ) 
        vlistDefVarName(ovlistID, varID, request->var2.name);
      if ( IS_SET(request->var2.longname) ) 
        vlistDefVarLongname(ovlistID, varID, request->var2.longname);
    }

  if ( cdoOperatorIntval(operatorID) == 17 ) vlistDefNtsteps(ovlistID, 1);

  itaxisID1 = vlistInqTaxis(ivlistID1);
  itaxisID2 = vlistInqTaxis(ivlistID2);
  otaxisID  = taxisCreate(TAXIS_ABSOLUTE);
  vlistDefTaxis(ovlistID, otaxisID);

  ostreamID = streamOpenWrite(cdoStreamName(2), cdoFiletype());
  if ( ostreamID < 0 ) cdiError(ostreamID, "Open failed on %s", cdoStreamName(2));

  streamDefVlist(ostreamID, ovlistID);

  nrecords   = vlistNrecs(ivlistID1);
  recVarID   = (int *) malloc(nrecords*sizeof(int));
  recLevelID = (int *) malloc(nrecords*sizeof(int));

  gridsize = gridInqSize(gridID);

  field1.ptr = (double *) malloc(gridsize*sizeof(double));
  field2.ptr = (double *) malloc(gridsize*sizeof(double));

  nlevels = zaxisInqSize(zaxisID);
  missval = vlistInqVarMissval(ivlistID1, FIRST_VAR_ID);

  var14 = (FIELD *) malloc(nlevels*sizeof(FIELD));
  samp1 = (FIELD *) malloc(nlevels*sizeof(FIELD));
  
  if ( request->var1.epilog == PERCENT_OF_TOTAL_AMOUNT )
    total  = (FIELD *) malloc(nlevels*sizeof(FIELD));
  if ( IS_SET(request->var1.f5) )
    var15 = (FIELD *) malloc(nlevels*sizeof(FIELD));
  if ( IS_SET(request->var2.h2) )
    var22 = (FIELD *) malloc(nlevels*sizeof(FIELD));
      
  for ( levelID = 0; levelID < nlevels; levelID++ )
    {
      var14[levelID].grid    = gridID;
      var14[levelID].nmiss   = 0;
      var14[levelID].missval = missval;
      var14[levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
      
      samp1[levelID].grid    = gridID;
      samp1[levelID].nmiss   = 0;
      samp1[levelID].missval = missval;
      samp1[levelID].ptr     = NULL;
      
      if ( request->var1.epilog == PERCENT_OF_TOTAL_AMOUNT )
        {
          total[levelID].grid    = gridID;
          total[levelID].nmiss   = 0;
          total[levelID].missval = missval;
          total[levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
        }
      if ( IS_SET(request->var1.f5) )
        {
          var15[levelID].grid    = gridID;
          var15[levelID].nmiss   = 0;
          var15[levelID].missval = missval;
          var15[levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
        }
      if ( IS_SET(request->var2.h2) )
        {
          var22[levelID].grid    = gridID;
          var22[levelID].nmiss   = 0;
          var22[levelID].missval = missval;
          var22[levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
        }
    }

  indate2 = 0;
  itsID   = 0;
  otsID   = 0;
  while ( TRUE )
    {
      nsets = 0;
      while ( (nrecs = streamInqTimestep(istreamID1, itsID)) > 0 )
        {
	  if ( !streamInqTimestep(istreamID2, itsID) )
	    cdoAbort("Input streams have different number of time steps!");
          /*
          ivdate1 = taxisInqVdate(itaxisID1);
          ivdate2 = taxisInqVdate(itaxisID2);
          if ( ivdate1 != ivdate2 )
            cdoAbort("Input streams have different verification dates for time step %d!", itsID+1);
            
          ivtime1 = taxisInqVtime(itaxisID1);
          ivtime2 = taxisInqVtime(itaxisID2);
          if ( ivtime1 != ivtime2 )
            cdoAbort("Input streams have different verification times for time step %d!", itsID+1);
	  */
          if ( nsets == 0 ) indate1 = (INT64)ivdate1*10000 + ivtime1;
          indate2 = (INT64)ivdate1*10000 + ivtime1;

          if ( indate2/intvdat != indate1/intvdat ) break;

          for ( recID = 0; recID < nrecs; recID++ )
            {
              streamInqRecord(istreamID1, &varID, &levelID);
              streamInqRecord(istreamID2, &varID, &levelID);

              if ( itsID == 0 )
                {
                  recVarID[recID]   = varID;
                  recLevelID[recID] = levelID;
                }
              if ( varID != FIRST_VAR_ID ) continue;

              if ( nsets == 0 )
                {
                  for ( i = 0; i < gridsize; i++ )
                    {
                      var14[levelID].ptr[i] = missval;
                      if ( IS_SET(samp1[levelID].ptr) )
                        samp1[levelID].ptr[i] = 0.0;
                      if ( request->var1.epilog == PERCENT_OF_TOTAL_AMOUNT )
                        total[levelID].ptr[i] = 0.0;
                      if ( IS_SET(request->var1.f5) ) 
                        var15[levelID].ptr[i] = missval;
                      if ( IS_SET(request->var2.h2) ) 
                        var22[levelID].ptr[i] = missval;
                    }
                  var14[levelID].nmiss = gridsize;
                  if ( request->var1.epilog == PERCENT_OF_TOTAL_AMOUNT )
                    total[levelID].nmiss = gridsize;
                  if ( IS_SET(request->var1.f5) )
                    var15[levelID].nmiss = gridsize;
                  if ( IS_SET(request->var2.h2) )
                    var22[levelID].nmiss = gridsize;
                }

              streamReadRecord(istreamID1, field1.ptr, &field1.nmiss);
              field1.grid    = var14[levelID].grid;
              field1.missval = var14[levelID].missval;

              streamReadRecord(istreamID2, field2.ptr, &field2.nmiss);
              field2.grid    = var14[levelID].grid;
              field2.missval = var14[levelID].missval;

              if ( request->var1.epilog == PERCENT_OF_TOTAL_AMOUNT )
                farsum(&total[levelID], field1);
                             
              if ( IS_SET(request->var1.f1) ) 
                request->var1.f1(&field1, request->var1.f1arg);
              
              if ( IS_SET(request->var1.f2) )
                request->var1.f2(&field2, request->var1.f2arg);

              if ( field1.nmiss > 0 || IS_SET(samp1[levelID].ptr) )
                {
                  if ( IS_NOT_SET(samp1[levelID].ptr) )
                    {
                      samp1[levelID].ptr = (double *) malloc(gridsize*sizeof(double));
                      for ( i = 0; i < gridsize; i++ )
                        samp1[levelID].ptr[i] = nsets;
                    }
                  for ( i = 0; i < gridsize; i++ )
                    {
                      if ( DBL_IS_EQUAL(field1.ptr[i], field1.missval) )
                        continue;
                      samp1[levelID].ptr[i]++;
                    }
                }
              
              request->var1.f3(&field1, field2);
              request->var1.f4(&var14[levelID], field1);

              if ( IS_SET(request->var2.h2) )
                {
                  memcpy(field2.ptr, var14[levelID].ptr, gridsize*sizeof(double));
                  field2.nmiss   = var14[levelID].nmiss;
                  field2.grid    = var14[levelID].grid;
                  field2.missval = var14[levelID].missval;

                  if ( IS_SET(request->var2.h1) )
                    request->var2.h1(&field2, request->var2.h1arg);
                    
                  request->var2.h2(&var22[levelID], field2);
                }

              if ( IS_SET(request->var1.f5) )
                request->var1.f5(&var15[levelID], var14[levelID], request->var1.f5arg);              
            }

          ovdate = ivdate1;
          ovtime = ivtime1;
          nsets++;
          itsID++;
        }

      if ( nrecs == 0 && nsets == 0 ) break;
      
      if ( request->var1.epilog == MEAN || request->var1.epilog == PERCENT_OF_TIME )
        for ( levelID = 0; levelID < nlevels; levelID++ )
          {
	    if ( IS_SET(request->var1.f5) )
	      var = &var15[levelID];
	    else
	      var = &var14[levelID];
	  
            if ( IS_NOT_SET(samp1[levelID].ptr) )
              farcdiv(var, nsets);
            else
              fardiv(var, samp1[levelID]);

	    if ( request->var1.epilog == PERCENT_OF_TIME )
	      farcmul(var, 100.0);
          }
      else if ( request->var1.epilog == PERCENT_OF_TOTAL_AMOUNT )
        for ( levelID = 0; levelID < nlevels; levelID++ )
          {
	    if ( IS_SET(request->var1.f5) )
	      var = &var15[levelID];
	    else
	      var = &var14[levelID];
	  
            fardiv(var, total[levelID]);
	    farcmul(var, 100.0);
          }

      taxisDefVdate(otaxisID, ovdate);
      taxisDefVtime(otaxisID, ovtime);
      streamDefTimestep(ostreamID, otsID++);

      if ( otsID == 1 || vlistInqVarTime(ivlistID1, FIRST_VAR_ID) == TIME_VARIABLE )
        {
          varID = 0;
          for ( levelID = 0; levelID < nlevels; levelID++ )
            {
              if ( IS_SET(request->var1.f5) )
                var = &var15[levelID];
              else 
                var = &var14[levelID];
               
              streamDefRecord(ostreamID, varID, levelID);
              streamWriteRecord(ostreamID, var->ptr, var->nmiss);
            }
          if ( IS_SET(request->var2.h2) )
            {
              varID = 1;
              for ( levelID = 0; levelID < nlevels; levelID++ )
                {
                  var = &var22[levelID];
                  
                  streamDefRecord(ostreamID, varID, levelID);
                  streamWriteRecord(ostreamID, var->ptr, var->nmiss);
               }
            }
        }

      if ( nrecs == 0 ) break;
    }

  for ( levelID = 0; levelID < nlevels; levelID++ )
    {
      free(var14[levelID].ptr);
      if ( IS_SET(samp1[levelID].ptr) ) free(samp1[levelID].ptr);
    }
  free(var14);
  free(samp1);
  
  if ( IS_SET(total) ) 
    {
      for ( levelID = 0; levelID < nlevels; levelID++ )
        free(total[levelID].ptr);
      free(total);
    }
  if ( IS_SET(var15) ) 
    {
      for ( levelID = 0; levelID < nlevels; levelID++ )
        free(var15[levelID].ptr);
      free(var15);
    }
  if ( IS_SET(var22) ) 
    {
      for ( levelID = 0; levelID < nlevels; levelID++ )
        free(var22[levelID].ptr);
      free(var22);
    }
  
  if ( IS_SET(field1.ptr) ) free(field1.ptr);
  if ( IS_SET(field2.ptr) ) free(field2.ptr);

  if ( IS_SET(recVarID) )   free(recVarID);
  if ( IS_SET(recLevelID) ) free(recLevelID);

  streamClose(ostreamID);
  streamClose(istreamID2);
  streamClose(istreamID1);
}


void eca3(const ECA_REQUEST_3 *request)
{
  static const char func[] = "eca3";
  const int operatorID = cdoOperatorID();
  
  INT64 intvdat;
  INT64 indate1 = 0, indate2;
  int gridsize;
  int ivdate1 = 0, ivtime1 = 0;
  int ivdate2 = 0, ivtime2 = 0;
  int ovdate = 0, ovtime = 0;
  int nrecs, nrecords;
  int gridID, zaxisID, varID, levelID, recID;
  int itsID;
  int otsID;
  long nsets;
  int i;
  int istreamID1, istreamID2, ostreamID;
  int ivlistID1, ivlistID2, ovlistID, itaxisID1, itaxisID2, otaxisID;
  int nlevels;
  int *recVarID, *recLevelID;
  double missval;
  FIELD *var1 = NULL, *var2 = NULL;
  FIELD field1, field2;
  
  intvdat = (INT64) pow(10.0, (double) cdoOperatorIntval(operatorID));

  istreamID1 = streamOpenRead(cdoStreamName(0));
  if ( istreamID1 < 0 ) cdiError(istreamID1, "Open failed on %s", cdoStreamName(0));
  istreamID2 = streamOpenRead(cdoStreamName(1));
  if ( istreamID2 < 0 ) cdiError(istreamID2, "Open failed on %s", cdoStreamName(1));

  ivlistID1 = streamInqVlist(istreamID1);
  ivlistID2 = streamInqVlist(istreamID2);
  ovlistID  = vlistCreate();
  
  vlistCompare(ivlistID1, ivlistID2, func_hrd);
  vlistCompare(ivlistID1, ivlistID2, func_sft);
  
  gridID  = vlistInqVarGrid(ivlistID1, FIRST_VAR_ID);
  zaxisID = vlistInqVarZaxis(ivlistID1, FIRST_VAR_ID);
  varID   = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARIABLE);
  
  if ( IS_SET(request->name) ) 
    vlistDefVarName(ovlistID, varID, request->name);
  if ( IS_SET(request->longname) ) 
    vlistDefVarLongname(ovlistID, varID, request->longname);
  if ( IS_SET(request->units) ) 
    vlistDefVarUnits(ovlistID, varID, request->units);

  if ( cdoOperatorIntval(operatorID) == 17 ) vlistDefNtsteps(ovlistID, 1);

  itaxisID1 = vlistInqTaxis(ivlistID1);
  itaxisID2 = vlistInqTaxis(ivlistID2);
  otaxisID  = taxisCreate(TAXIS_ABSOLUTE);
  vlistDefTaxis(ovlistID, otaxisID);

  ostreamID = streamOpenWrite(cdoStreamName(2), cdoFiletype());
  if ( ostreamID < 0 ) cdiError(ostreamID, "Open failed on %s", cdoStreamName(2));

  streamDefVlist(ostreamID, ovlistID);

  nrecords   = vlistNrecs(ivlistID1);
  recVarID   = (int *) malloc(nrecords*sizeof(int));
  recLevelID = (int *) malloc(nrecords*sizeof(int));

  gridsize = gridInqSize(gridID);

  field1.ptr = (double *) malloc(gridsize*sizeof(double));
  field2.ptr = (double *) malloc(gridsize*sizeof(double));

  nlevels = zaxisInqSize(zaxisID);
  missval = vlistInqVarMissval(ivlistID1, FIRST_VAR_ID);

  var1 = (FIELD *) malloc(nlevels*sizeof(FIELD));
  var2 = (FIELD *) malloc(nlevels*sizeof(FIELD));
        
  for ( levelID = 0; levelID < nlevels; levelID++ )
    {
      var1[levelID].grid    = gridID;
      var1[levelID].nmiss   = 0;
      var1[levelID].missval = missval;
      var1[levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
            
      var2[levelID].grid    = gridID;
      var2[levelID].nmiss   = 0;
      var2[levelID].missval = missval;
      var2[levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
    }

  indate2 = 0;
  itsID   = 0;
  otsID   = 0;
  while ( TRUE )
    {
      nsets = 0;
      while ( (nrecs = streamInqTimestep(istreamID1, itsID)) > 0 )
        {
	  if ( !streamInqTimestep(istreamID2, itsID) )
	    cdoAbort("Input streams have different number of time steps!");
            
          ivdate1 = taxisInqVdate(itaxisID1);
          ivdate2 = taxisInqVdate(itaxisID2);
          if ( ivdate1 != ivdate2 )
            cdoAbort("Input streams have different verification dates for time step %d!", itsID+1);
            
          ivtime1 = taxisInqVtime(itaxisID1);
          ivtime2 = taxisInqVtime(itaxisID2);
          if ( ivtime1 != ivtime2 )
            cdoAbort("Input streams have different verification times for time step %d!", itsID+1);

          if ( nsets == 0 ) indate1 = (INT64)ivdate1*10000 + ivtime1;
          indate2 = (INT64)ivdate1*10000 + ivtime1;

          if ( indate2/intvdat != indate1/intvdat ) break;

          for ( recID = 0; recID < nrecs; recID++ )
            {
              streamInqRecord(istreamID1, &varID, &levelID);
              streamInqRecord(istreamID2, &varID, &levelID);

              if ( itsID == 0 )
                {
                  recVarID[recID]   = varID;
                  recLevelID[recID] = levelID;
                }
              if ( varID != FIRST_VAR_ID ) continue;

              if ( nsets == 0 )
                {
                  for ( i = 0; i < gridsize; i++ )
                    {
                      var1[levelID].ptr[i] = missval;
                      var2[levelID].ptr[i] = missval;
                    }
                  var1[levelID].nmiss = gridsize;
                  var2[levelID].nmiss = gridsize;
                }

              streamReadRecord(istreamID1, field1.ptr, &field1.nmiss);
              field1.grid    = var1[levelID].grid;
              field1.missval = var1[levelID].missval;

              streamReadRecord(istreamID2, field2.ptr, &field2.nmiss);
              field2.grid    = var1[levelID].grid;
              field2.missval = var1[levelID].missval;

              request->f1(&var1[levelID], field1);
              request->f2(&var2[levelID], field2);
            }

          ovdate = ivdate1;
          ovtime = ivtime1;
          nsets++;
          itsID++;
        }

      if ( nrecs == 0 && nsets == 0 ) break;
      
      for ( levelID = 0; levelID < nlevels; levelID++ )
        request->f3(&var1[levelID], var2[levelID]);

      taxisDefVdate(otaxisID, ovdate);
      taxisDefVtime(otaxisID, ovtime);
      streamDefTimestep(ostreamID, otsID++);

      if ( otsID == 1 || vlistInqVarTime(ivlistID1, FIRST_VAR_ID) == TIME_VARIABLE )
        {
          varID = 0;
          for ( levelID = 0; levelID < nlevels; levelID++ )
            {
              streamDefRecord(ostreamID, varID, levelID);
              streamWriteRecord(ostreamID, var1[levelID].ptr, var1[levelID].nmiss);
            }
        }

      if ( nrecs == 0 ) break;
    }

  for ( levelID = 0; levelID < nlevels; levelID++ )
    {
      free(var1[levelID].ptr);
      free(var2[levelID].ptr);
    }
  free(var1);
  free(var2);
    
  if ( IS_SET(field1.ptr) ) free(field1.ptr);
  if ( IS_SET(field2.ptr) ) free(field2.ptr);

  if ( IS_SET(recVarID) )   free(recVarID);
  if ( IS_SET(recLevelID) ) free(recLevelID);

  streamClose(ostreamID);
  streamClose(istreamID2);
  streamClose(istreamID1);
}


void eca4(const ECA_REQUEST_4 *request)
{
  static const char func[] = "eca4";
  const int operatorID = cdoOperatorID();
  
  INT64 intvdat;
  INT64 indate1 = 0, indate2;
  int gridsize;
  int ivdate = 0, ivtime = 0;
  int ovdate = 0, ovtime = 0;
  int month;
  int nrecs, nrecords;
  int gridID, zaxisID, varID, levelID, recID;
  int itsID;
  int otsID;
  long nsets;
  int i;
  int istreamID, ostreamID;
  int ivlistID, ovlistID, itaxisID, otaxisID;
  int nlevels;
  int *recVarID, *recLevelID;
  double missval;
  FIELD *var1 = NULL, *var2 = NULL, *sdat = NULL, *edat = NULL, *var3 = NULL, *done = NULL;
  FIELD field;
  
  intvdat = (INT64) pow(10.0, (double) cdoOperatorIntval(operatorID));

  istreamID = streamOpenRead(cdoStreamName(0));
  if ( istreamID < 0 ) cdiError(istreamID, "Open failed on %s", cdoStreamName(0));

  ivlistID = streamInqVlist(istreamID);
  ovlistID = vlistCreate();
  
  gridID  = vlistInqVarGrid(ivlistID, FIRST_VAR_ID);
  zaxisID = vlistInqVarZaxis(ivlistID, FIRST_VAR_ID);
  varID   = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARIABLE);
  
  if ( IS_SET(request->name) )
    vlistDefVarName(ovlistID, varID, request->name);
  if ( IS_SET(request->longname) ) 
    vlistDefVarLongname(ovlistID, varID, request->longname);
  if ( IS_SET(request->units) ) 
    vlistDefVarUnits(ovlistID, varID, request->units);

  varID = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARIABLE);
  
  if ( IS_SET(request->name2) ) 
    vlistDefVarName(ovlistID, varID, request->name2);
  if ( IS_SET(request->longname2) ) 
    vlistDefVarLongname(ovlistID, varID, request->longname2);
  if ( IS_SET(request->units2) ) 
    vlistDefVarUnits(ovlistID, varID, request->units2);

  varID = vlistDefVar(ovlistID, gridID, zaxisID, TIME_VARIABLE);
  
  if ( IS_SET(request->name3) ) 
    vlistDefVarName(ovlistID, varID, request->name3);
  if ( IS_SET(request->longname3) ) 
    vlistDefVarLongname(ovlistID, varID, request->longname3);
  if ( IS_SET(request->units3) ) 
    vlistDefVarUnits(ovlistID, varID, request->units3);
    
  if ( cdoOperatorIntval(operatorID) == 17 ) vlistDefNtsteps(ovlistID, 1);

  itaxisID = vlistInqTaxis(ivlistID);
  otaxisID = taxisCreate(TAXIS_ABSOLUTE);
  vlistDefTaxis(ovlistID, otaxisID);

  ostreamID = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( ostreamID < 0 ) cdiError(ostreamID, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(ostreamID, ovlistID);

  nrecords   = vlistNrecs(ivlistID);
  recVarID   = (int *) malloc(nrecords*sizeof(int));
  recLevelID = (int *) malloc(nrecords*sizeof(int));

  gridsize = gridInqSize(gridID);

  field.ptr = (double *) malloc(gridsize*sizeof(double));

  nlevels = zaxisInqSize(zaxisID);
  missval = vlistInqVarMissval(ivlistID, FIRST_VAR_ID);

  var1 = (FIELD *) malloc(nlevels*sizeof(FIELD));
  var2 = (FIELD *) malloc(nlevels*sizeof(FIELD));
  var3 = (FIELD *) malloc(nlevels*sizeof(FIELD));
  sdat = (FIELD *) malloc(nlevels*sizeof(FIELD));
  edat = (FIELD *) malloc(nlevels*sizeof(FIELD));
  done = (FIELD *) malloc(nlevels*sizeof(FIELD));
    
  for ( levelID = 0; levelID < nlevels; levelID++ )
    {
      var1[levelID].grid    = gridID;
      var1[levelID].nmiss   = 0;
      var1[levelID].missval = missval;
      var1[levelID].ptr     = (double *) malloc(gridsize*sizeof(double));

      var2[levelID].grid    = gridID;
      var2[levelID].nmiss   = 0;
      var2[levelID].missval = missval;
      var2[levelID].ptr     = (double *) malloc(gridsize*sizeof(double));

      var3[levelID].grid    = gridID;
      var3[levelID].nmiss   = 0;
      var3[levelID].missval = missval;
      var3[levelID].ptr     = (double *) malloc(gridsize*sizeof(double));

      sdat[levelID].grid    = gridID;
      sdat[levelID].nmiss   = 0;
      sdat[levelID].missval = missval;
      sdat[levelID].ptr     = (double *) malloc(gridsize*sizeof(double));

      edat[levelID].grid    = gridID;
      edat[levelID].nmiss   = 0;
      edat[levelID].missval = missval;
      edat[levelID].ptr     = (double *) malloc(gridsize*sizeof(double));

      done[levelID].grid    = gridID;
      done[levelID].nmiss   = 0;
      done[levelID].missval = missval;
      done[levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
    }

  indate2 = 0;
  itsID   = 0;
  otsID   = 0;
  while ( TRUE )
    {
      nsets = 0;
      while ( (nrecs = streamInqTimestep(istreamID, itsID)) > 0 )
        {
          ivdate = taxisInqVdate(itaxisID);
          ivtime = taxisInqVtime(itaxisID);

          if ( nsets == 0 ) indate1 = (INT64)ivdate*10000 + ivtime;
          indate2 = (INT64)ivdate*10000 + ivtime;

          if ( indate2/intvdat != indate1/intvdat ) break;

          for ( recID = 0; recID < nrecs; recID++ )
            {
              streamInqRecord(istreamID, &varID, &levelID);

              if ( itsID == 0 )
                {
                  recVarID[recID]   = varID;
                  recLevelID[recID] = levelID;
                }
              if ( varID != FIRST_VAR_ID ) continue;

              if ( nsets == 0 )
                {
                  for ( i = 0; i < gridsize; i++ )
                    {
                      var1[levelID].ptr[i] = missval;
                      var2[levelID].ptr[i] = missval;
                      var3[levelID].ptr[i] = missval;
                      sdat[levelID].ptr[i] = missval;
                      edat[levelID].ptr[i] = missval;
                      done[levelID].ptr[i] = 0.0;
                    }
                  var1[levelID].nmiss = gridsize;
                  var2[levelID].nmiss = gridsize;
                  var3[levelID].nmiss = gridsize;
                  sdat[levelID].nmiss = gridsize;
                  edat[levelID].nmiss = gridsize;
                  done[levelID].nmiss = 0;
                }

              streamReadRecord(istreamID, field.ptr, &field.nmiss);
              field.grid    = var1[levelID].grid;
              field.missval = var1[levelID].missval;
              
              month = (ivdate % 10000) / 100;
              if ( month < 1 || month > 12 )
	        cdoAbort("month %d out of range!", month);
              
              if ( month < 7 ) 
                {
                  request->s1(&field, request->s1arg);
                  farnum2(&var1[levelID], field);
                  
                  for ( i = 0; i < gridsize; i++ )
                    {
                      if ( done[levelID].ptr[i] != 1.0 )
                        {
                          if ( var1[levelID].ptr[i] == request->consecutiveDays )
                            {
                              var3[levelID].ptr[i] = 0.0;
                              done[levelID].ptr[i] = 1.0;
                            } 
                        }  
                      else 
                        {
                          if ( DBL_IS_EQUAL(sdat[levelID].ptr[i], missval) )
                            sdat[levelID].ptr[i] = ivdate;
                            
                          var3[levelID].ptr[i] += 1.0;
                        }  
                    }
                }
              else
                {
                  request->s2(&field, request->s2arg);
                  farnum2(&var2[levelID], field);

                  for ( i = 0; i < gridsize; i++ )
                    {
                      if ( done[levelID].ptr[i] != 1.0 )
                        continue;
                        
                      if ( DBL_IS_EQUAL(sdat[levelID].ptr[i], missval) )
                        sdat[levelID].ptr[i] = ivdate;
                         
                      if ( var2[levelID].ptr[i] == 0.0 )
                        edat[levelID].ptr[i] = ivdate;
                        
                      if ( var2[levelID].ptr[i] < request->consecutiveDays )
                        var3[levelID].ptr[i] += 1.0;
                      else
                        {
                          assert(var2[levelID].ptr[i] == request->consecutiveDays);
                          var3[levelID].ptr[i] = var3[levelID].ptr[i] - (request->consecutiveDays - 1);
                          done[levelID].ptr[i] = 0.0;
                        } 
                    }
                }
                
                var3[levelID].nmiss = 0;
                sdat[levelID].nmiss = 0;
                edat[levelID].nmiss = 0;
                for ( i = 0; i < gridsize; i++ )
                  {
                    if ( DBL_IS_EQUAL(var3[levelID].ptr[i], missval) )
                      var3[levelID].nmiss++;
                    if ( DBL_IS_EQUAL(sdat[levelID].ptr[i], missval) )
                      sdat[levelID].nmiss++;
                    if ( DBL_IS_EQUAL(edat[levelID].ptr[i], missval) )
                      edat[levelID].nmiss++;
                  }
            }

          ovdate = ivdate;
          ovtime = ivtime;
          nsets++;
          itsID++;
        }

      if ( nrecs == 0 && nsets == 0 ) break;
      
      taxisDefVdate(otaxisID, ovdate);
      taxisDefVtime(otaxisID, ovtime);
      streamDefTimestep(ostreamID, otsID++);

      if ( otsID == 1 || vlistInqVarTime(ivlistID, FIRST_VAR_ID) == TIME_VARIABLE )
        {
          varID = 0;
          for ( levelID = 0; levelID < nlevels; levelID++ )
            {
              streamDefRecord(ostreamID, varID, levelID);
              streamWriteRecord(ostreamID, var3[levelID].ptr, var3[levelID].nmiss);
            }
          varID = 1;
          for ( levelID = 0; levelID < nlevels; levelID++ )
            {
              streamDefRecord(ostreamID, varID, levelID);
              streamWriteRecord(ostreamID, sdat[levelID].ptr, sdat[levelID].nmiss);
            }
          varID = 2;
          for ( levelID = 0; levelID < nlevels; levelID++ )
            {
              streamDefRecord(ostreamID, varID, levelID);
              streamWriteRecord(ostreamID, edat[levelID].ptr, edat[levelID].nmiss);
            }
        }

      if ( nrecs == 0 ) break;
    }

  for ( levelID = 0; levelID < nlevels; levelID++ )
    {
      free(var1[levelID].ptr);
      free(var2[levelID].ptr);
      free(var3[levelID].ptr);
      free(sdat[levelID].ptr);
      free(edat[levelID].ptr);
      free(done[levelID].ptr);
    }  
  free(var1);
  free(var2);
  free(var3);
  free(sdat);
  free(edat);
  free(done);
  
  if ( IS_SET(field.ptr) ) free(field.ptr);

  if ( IS_SET(recVarID) )   free(recVarID);
  if ( IS_SET(recLevelID) ) free(recLevelID);

  streamClose(ostreamID);
  streamClose(istreamID);
}

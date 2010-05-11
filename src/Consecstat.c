/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2010 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Consectstep  consecsum  For each timestep, the current number of 
                              onsecutive timsteps is counted
      Consectstep  consects   For each period of consecutive timesteps, only
                              count its lenght + last contributing timesteps

   =============================================================================
   Created:  04/08/2010 11:58:01 AM
    Author:  Ralf Mueller (ram), ralf.mueller@zmaw.de
   Company:  Max-Planck-Institute for Meteorology
   =============================================================================
 */

#if defined (_OPENMP)
#  include <omp.h>
#endif

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "field.h"
#include "pstream.h"

enum {CONSECSUM, CONSECTS};
#define SWITCHWARN "Hit default case!This should never happen (%s).\n"

static void selEndOfPeriod(field_t *periods, field_t history, field_t current, int isLastTimestep)
{
  static char func[] = "selEndOfPeriod";
  long   i, len;
  double pmissval = periods->missval;
  double  *parray = periods->ptr;
  double hmissval = history.missval;
  double  *harray = history.ptr;
  double cmissval = current.missval;
  double  *carray = current.ptr;

  len = gridInqSize(periods->grid);
  if ( len != gridInqSize(current.grid) || (gridInqSize(current.grid) != gridInqSize(history.grid)) )
    cdoAbort("Fields have different gridsize (%s)", func);

  if (!isLastTimestep)
  {
    if ( current.nmiss > 0 || history.nmiss > 0 )
    {
#if defined (_OPENMP)
#pragma omp parallel for default(shared) schedule(static)
#endif
      for ( i = 0; i < len; i++ )
      {
        if ( !DBL_IS_EQUAL(harray[i], hmissval) )
        {
          if ( !DBL_IS_EQUAL(carray[i], cmissval) )
          {
            parray[i] = ( DBL_IS_EQUAL(carray[i], 0.0)  && IS_NOT_EQUAL(harray[i], 0.0) ) ? harray[i] : pmissval;
          }
          else /* DBL_IS_EQUAL(carray[i], cmissval) */
          {
            parray[i] = ( IS_NOT_EQUAL(harray[i], 0.0) ) ? harray[i] : pmissval;
          }
        }
        else /* DBL_IS_EQUAL(harray[i], hmissval) */
        {
          parray[i] = pmissval;
        }
      }
    }
    else
    {
#if defined (_OPENMP)
#pragma omp parallel for default(shared) schedule(static)
#endif
      for ( i = 0; i < len; i++ )
        parray[i] = ( DBL_IS_EQUAL(carray[i], 0.0)  && IS_NOT_EQUAL(harray[i], 0.0) ) ? harray[i] : pmissval;
    }
  }
  else
  {
    if ( current.nmiss > 0 )
    {
#if defined (_OPENMP)
#pragma omp parallel for default(shared) schedule(static)
#endif
      for ( i = 0; i < len; i++ )
        if ( !DBL_IS_EQUAL(carray[i], cmissval) )
        {
          parray[i] = ( DBL_IS_EQUAL(carray[i], 0.0) ) ? pmissval : carray[i];
        }
        else /* DBL_IS_EQUAL(carray[i], cmissval) */
        {
          parray[i] = pmissval;
        }
    }
    else
    {
#if defined (_OPENMP)
#pragma omp parallel for default(shared) schedule(static)
#endif
      for ( i = 0; i < len; i++ )
        parray[i] = ( DBL_IS_EQUAL(carray[i], 0.0) ) ? pmissval : carray[i];
    }
  }

  periods->nmiss = 0;
  for ( i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(parray[i], pmissval) ) periods->nmiss++;
}


void *Consecstat (void *argument)
{
  static char func[] = "Consecstat";
  int operatorID;

  int i;
  int istreamID, itaxisID, ivlistID, itsID;
  int ostreamID, otaxisID, ovlistID, otsID;
  int vdate, vtime;
  int histvdate = 0, histvtime = 0;
  int recID, nrecs;
  int varID, nvars;
  int levelID, nlevels; 
  int gridID, gridsize;
  double missval;
  double refval = 0.0;

  field_t **vars = NULL, **hist = NULL, **periods = NULL;
  field_t field;

  cdoInitialize(argument);
  cdoOperatorAdd("consecsum",CONSECSUM, 0, "refval");
  cdoOperatorAdd("consects" ,CONSECTS , 0, NULL);
  operatorID = cdoOperatorID();

  if ( operatorID == CONSECSUM )
    if ( operatorArgc() > 0 ) refval = atof(operatorArgv()[0]);

  istreamID = streamOpenRead(cdoStreamName(0));
  if ( istreamID < 0 ) cdiError(istreamID, "Open failed on %s", cdoStreamName(0));

  ivlistID = streamInqVlist(istreamID);
  itaxisID = vlistInqTaxis(ivlistID);
  ovlistID = vlistDuplicate(ivlistID);
  otaxisID = taxisDuplicate(itaxisID);
  vlistDefTaxis(ovlistID, otaxisID);

  field.ptr = (double *) malloc(vlistGridsizeMax(ovlistID)*sizeof(double));
  nvars     = vlistNvars(ivlistID);
  vars      = field_palloc(ivlistID, FIELD_PTR, FIELD_PINIT_STATIC);
  hist      = field_palloc(ivlistID, FIELD_PTR, FIELD_PINIT_STATIC);
  if ( operatorID == CONSECTS )
    periods = field_palloc(ivlistID, FIELD_PTR, FIELD_PINIT_STATIC);

  for ( varID = 0; varID < nvars; varID++ )
    vlistDefVarUnits(ovlistID, varID, "steps");

  ostreamID = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( ostreamID < 0 ) cdiError(ostreamID, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(ostreamID, ovlistID);

  itsID = 0;
  otsID = 0;
  while ( (nrecs = streamInqTimestep(istreamID, itsID)) )
  {
    vdate = taxisInqVdate(itaxisID);
    vtime = taxisInqVtime(itaxisID);
    switch (operatorID)
    {
      case CONSECSUM:
        taxisDefVdate(otaxisID, vdate);
        taxisDefVtime(otaxisID, vtime);
        streamDefTimestep(ostreamID, otsID);
        break;
      case CONSECTS:
        if (itsID != 0 )
        {
          taxisDefVdate(otaxisID, histvdate);
          taxisDefVtime(otaxisID, histvtime);
          streamDefTimestep(ostreamID, otsID-1);
        }
        break;
      default:
        printf (SWITCHWARN,func);
    }

    for ( recID = 0; recID < nrecs; recID++ )
    {
      streamInqRecord(istreamID, &varID, &levelID);
      streamReadRecord(istreamID, field.ptr, &field.nmiss);
      field.grid    = vlistInqVarGrid(ovlistID, varID);
      field.missval = vlistInqVarMissval(ovlistID, varID);
/* Localize 4 OpenMP */
#if defined (_OPENMP)
#pragma omp parallel for default(shared) schedule(static)
      for ( i = 0; i < gridInqSize(field.grid); i++ )
        field.ptr[i] = field.ptr[i];
#endif
      farsumtr(&vars[varID][levelID], field, refval);
      switch (operatorID)
      {
        case CONSECSUM:
          streamDefRecord(ostreamID, varID, levelID);
          streamWriteRecord(ostreamID, vars[varID][levelID].ptr, vars[varID][levelID].nmiss);
          break;
        case CONSECTS:
          if ( itsID != 0 )
          {
            selEndOfPeriod(&periods[varID][levelID], hist[varID][levelID], vars[varID][levelID], FALSE);
            streamDefRecord(ostreamID, varID, levelID);
            streamWriteRecord(ostreamID, periods[varID][levelID].ptr, periods[varID][levelID].nmiss);
          }
          memcpy(hist[varID][levelID].ptr,
                 vars[varID][levelID].ptr,
                 gridInqSize(vars[varID][levelID].grid)*sizeof(double));
          break;
        default:
          printf (SWITCHWARN,func);
      }
    }
    histvdate = vdate;
    histvtime = vtime;
    itsID++;
    otsID++;
  }

  if ( operatorID == CONSECTS ) /* Save the last timestep */
  {
    taxisDefVdate(otaxisID, vdate);
    taxisDefVtime(otaxisID, vtime);
    streamDefTimestep(ostreamID, otsID-1);
    for ( varID = 0; varID < nvars; varID++ )
    {
      nlevels = zaxisInqSize(vlistInqVarZaxis(ovlistID, varID));
      for ( levelID = 0; levelID < nlevels; levelID++ )
      {
        selEndOfPeriod(&periods[varID][levelID], hist[varID][levelID], vars[varID][levelID], TRUE);
        streamDefRecord(ostreamID, varID, levelID);
        streamWriteRecord(ostreamID, periods[varID][levelID].ptr, periods[varID][levelID].nmiss);
      }
    }
  }
  
  field_free(vars, ivlistID);
  field_free(hist, ivlistID);
  if ( operatorID == CONSECTS )
    field_free(periods, ivlistID);

  streamClose(istreamID);
  streamClose(ostreamID);

  cdoFinish();

  return (0);
}

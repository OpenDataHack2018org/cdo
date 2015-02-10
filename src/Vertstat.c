/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2015 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
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

      Vertstat   vertmin         Vertical minimum
      Vertstat   vertmax         Vertical maximum
      Vertstat   vertsum         Vertical sum
      Vertstat   vertmean        Vertical mean
      Vertstat   vertavg         Vertical average
      Vertstat   vertvar         Vertical variance
      Vertstat   vertstd         Vertical standard deviation
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define IS_SURFACE_LEVEL(zaxisID)  (zaxisInqType(zaxisID) == ZAXIS_SURFACE && zaxisInqSize(zaxisID) == 1)

static
int getSurfaceID(int vlistID)
{
  int surfID = -1;
  int zaxisID;
  int nzaxis = vlistNzaxis(vlistID);

  for ( int index = 0; index < nzaxis; ++index )
    {
      zaxisID = vlistZaxis(vlistID, index);
      if ( IS_SURFACE_LEVEL(zaxisID) )
	{
	  surfID = vlistZaxis(vlistID, index);
	  break;
	}
    }

  if ( surfID == -1 ) surfID = zaxisCreate(ZAXIS_SURFACE, 1);

  return surfID;
}

static
void setSurfaceID(int vlistID, int surfID)
{
  int zaxisID;
  int nzaxis = vlistNzaxis(vlistID);

  for ( int index = 0; index < nzaxis; ++index )
    {
      zaxisID = vlistZaxis(vlistID, index);
      if ( zaxisID != surfID || !IS_SURFACE_LEVEL(zaxisID) )
	vlistChangeZaxisIndex(vlistID, index, surfID);
    }
}

static
void getLayerThickness(int zaxisID, int nlev, double *weights)
{
  int i;
  double *levels  = (double *) malloc(nlev*sizeof(double));
  double *lbounds = (double *) malloc(nlev*sizeof(double));
  double *ubounds = (double *) malloc(nlev*sizeof(double));

  zaxisInqLevels(zaxisID, levels);
  if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
    {
      zaxisInqLbounds(zaxisID, lbounds);
      zaxisInqUbounds(zaxisID, ubounds);  
    }
  else
    {
      lbounds[0] = levels[0];
      ubounds[nlev-1] = levels[nlev-1];
      for ( i = 0; i < nlev-1; ++i )
	{
	  double bound = 0.5*(levels[i] + levels[i+1]);
	  lbounds[i+1] = bound;
	  ubounds[i]   = bound;
	}		  
    }

  for ( i = 0; i < nlev; ++i ) weights[i] = fabs(ubounds[i]-lbounds[i]);

  double wsum = 0;
  for ( i = 0; i < nlev; ++i ) wsum += weights[i];
  /*
    for ( i = 0; i < nlev; ++i ) weights[i] /= wsum;
  */

  if ( cdoVerbose )
    {
      printf("zaxisID=%d  nlev=%d\n", zaxisID, nlev);
      printf("levels=");
      for ( i = 0; i < nlev; ++i ) printf("%g ", levels[i]);
      printf("\n");
      printf("bounds=");
      for ( i = 0; i < nlev; ++i ) printf("%g/%g ", lbounds[i], ubounds[i]);
      printf("\n");
      printf("weights=");
      for ( i = 0; i < nlev; ++i ) printf("%g ", weights[i]);
      printf("\n");
      printf("weight sum: %g\n", wsum);
    }

  free(levels);
  free(lbounds);
  free(ubounds);
}

void *Vertstat(void *argument)
{
  int recID, nrecs;
  int gridID;
  int i;
  int varID, levelID;
  int nmiss;
  double missval;
  field_t *vars1 = NULL, *vars2 = NULL, *samp1 = NULL;
  field_t field;
  int needWeights = FALSE;
  typedef struct {
    int zaxisID;
    int numlevel;
    double *weights;
  }
  vert_t;

  cdoInitialize(argument);

                cdoOperatorAdd("vertmin",  func_min,  0, NULL);
                cdoOperatorAdd("vertmax",  func_max,  0, NULL);
                cdoOperatorAdd("vertsum",  func_sum,  0, NULL);
  int VERTINT = cdoOperatorAdd("vertint",  func_sum,  0, NULL);
                cdoOperatorAdd("vertmean", func_mean, 0, NULL);
                cdoOperatorAdd("vertavg",  func_avg,  0, NULL);
                cdoOperatorAdd("vertvar",  func_var,  0, NULL);
                cdoOperatorAdd("vertstd",  func_std,  0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc   = cdoOperatorF1(operatorID);

  if ( operatorID == VERTINT ) needWeights = TRUE;

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);

  vlistClearFlag(vlistID1);
  int nvars = vlistNvars(vlistID1);
  for ( varID = 0; varID < nvars; varID++ )
    vlistDefFlag(vlistID1, varID, 0, TRUE);

  int vlistID2 = vlistCreate();
  vlistCopyFlag(vlistID2, vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int surfID = getSurfaceID(vlistID1);
  setSurfaceID(vlistID2, surfID);

  int nzaxis = vlistNzaxis(vlistID1);
  int nlev, zaxisID;
  vert_t vert[nzaxis];
  if ( needWeights )
    {
      for ( int index = 0; index < nzaxis; ++index )
	{
	  zaxisID = vlistZaxis(vlistID1, index);
	  nlev = zaxisInqSize(zaxisID);
	  vert[index].numlevel = 0;
	  vert[index].zaxisID  = zaxisID;
	  if ( nlev > 1 )
	    {
	      vert[index].numlevel = nlev;
	      vert[index].weights = (double *) malloc(nlev*sizeof(double));
	      getLayerThickness(zaxisID, nlev, vert[index].weights); 
	    }
	}
    }

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  int gridsize = vlistGridsizeMax(vlistID1);

  field_init(&field);
  field.ptr = (double*) malloc(gridsize*sizeof(double));

  vars1 = (field_t*) malloc(nvars*sizeof(field_t));
  samp1 = (field_t*) malloc(nvars*sizeof(field_t));
  if ( operfunc == func_std || operfunc == func_var )
    vars2 = (field_t*) malloc(nvars*sizeof(field_t));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      zaxisID  = vlistInqVarZaxis(vlistID1, varID);
      missval  = vlistInqVarMissval(vlistID1, varID);

      field_init(&vars1[varID]);
      field_init(&samp1[varID]);
      vars1[varID].grid    = gridID;
      vars1[varID].zaxis   = zaxisID;
      vars1[varID].nsamp   = 0;
      vars1[varID].nmiss   = 0;
      vars1[varID].missval = missval;
      vars1[varID].ptr     = (double*) malloc(gridsize*sizeof(double));
      samp1[varID].grid    = gridID;
      samp1[varID].nmiss   = 0;
      samp1[varID].missval = missval;
      samp1[varID].ptr     = NULL;
      if ( operfunc == func_std || operfunc == func_var )
	{
	  field_init(&vars2[varID]);
	  vars2[varID].grid    = gridID;
	  vars2[varID].nmiss   = 0;
	  vars2[varID].missval = missval;
	  vars2[varID].ptr     = (double*) malloc(gridsize*sizeof(double));
	}
    }

  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);

          vars1[varID].nsamp++;
	  gridsize = gridInqSize(vars1[varID].grid);
	  zaxisID = vars1[varID].zaxis;

	  double vweight = 1.e30;
	  if ( operatorID == VERTINT && !IS_SURFACE_LEVEL(zaxisID) )
	    for ( int index = 0; index < nzaxis; ++index )
	      if ( vert[index].zaxisID == zaxisID )
		{
		  vweight = vert[index].weights[levelID];
		  break;
		}

	  if ( levelID == 0 )
	    {
	      streamReadRecord(streamID1, vars1[varID].ptr, &nmiss);
	      vars1[varID].nmiss = nmiss;

	      if ( vweight < 1.e30 ) farcmul(&vars1[varID], vweight);

	      if ( operfunc == func_std || operfunc == func_var )
		farmoq(&vars2[varID], vars1[varID]);

	      if ( nmiss > 0 || samp1[varID].ptr )
		{
		  if ( samp1[varID].ptr == NULL )
		    samp1[varID].ptr = (double*) malloc(gridsize*sizeof(double));

		  for ( i = 0; i < gridsize; i++ )
		    if ( DBL_IS_EQUAL(vars1[varID].ptr[i], vars1[varID].missval) )
		      samp1[varID].ptr[i] = 0;
		    else
		      samp1[varID].ptr[i] = 1;
		}
	    }
	  else
	    {
	      streamReadRecord(streamID1, field.ptr, &field.nmiss);
	      field.grid    = vars1[varID].grid;
	      field.missval = vars1[varID].missval;

	      if ( vweight < 1.e30 ) farcmul(&field, vweight);

	      if ( field.nmiss > 0 || samp1[varID].ptr )
		{
		  if ( samp1[varID].ptr == NULL )
		    {
		      samp1[varID].ptr = (double*) malloc(gridsize*sizeof(double));
		      for ( i = 0; i < gridsize; i++ )
			samp1[varID].ptr[i] = vars1[varID].nsamp;
		    }

		  for ( i = 0; i < gridsize; i++ )
		    if ( !DBL_IS_EQUAL(field.ptr[i], vars1[varID].missval) )
		      samp1[varID].ptr[i]++;
		}

	      if ( operfunc == func_std || operfunc == func_var )
		{
		  farsumq(&vars2[varID], field);
		  farsum(&vars1[varID], field);
		}
	      else
		{
		  farfun(&vars1[varID], field, operfunc);
		}
	    }
	}

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vars1[varID].nsamp )
	    {
	      if ( operfunc == func_mean || operfunc == func_avg )
		{
		  if ( samp1[varID].ptr == NULL )
		    farcmul(&vars1[varID], 1.0/vars1[varID].nsamp);
		  else
		    fardiv(&vars1[varID], samp1[varID]);
		}
	      else if ( operfunc == func_std || operfunc == func_var )
		{
		  if ( samp1[varID].ptr == NULL )
		    {
		      if ( operfunc == func_std )
			farcstd(&vars1[varID], vars2[varID], 1.0/vars1[varID].nsamp);
		      else
			farcvar(&vars1[varID], vars2[varID], 1.0/vars1[varID].nsamp);
		    }
		  else
		    {
		      farinv(&samp1[varID]);
		      if ( operfunc == func_std )
			farstd(&vars1[varID], vars2[varID], samp1[varID]);
		      else
			farvar(&vars1[varID], vars2[varID], samp1[varID]);
		    }
		}

	      streamDefRecord(streamID2, varID, 0);
	      streamWriteRecord(streamID2, vars1[varID].ptr, vars1[varID].nmiss);
	      vars1[varID].nsamp = 0;
	    }
	}

      tsID++;
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      free(vars1[varID].ptr);
      if ( samp1[varID].ptr ) free(samp1[varID].ptr);
      if ( operfunc == func_std || operfunc == func_var ) free(vars2[varID].ptr);
    }

  free(vars1);
  free(samp1);
  if ( operfunc == func_std || operfunc == func_var ) free(vars2);

  if ( field.ptr ) free(field.ptr);

  if ( needWeights )
    for ( int index = 0; index < nzaxis; ++index )
      if ( vert[index].numlevel > 1 )  free(vert[index].weights);

  streamClose(streamID2);
  streamClose(streamID1);

  vlistDestroy(vlistID2);

  cdoFinish();

  return (0);
}

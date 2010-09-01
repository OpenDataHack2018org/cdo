/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2006 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Mastrfu    mastrfu         Mass stream function
*/


#include <math.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


static void mastrfu(int gridID, int zaxisID, double *array1, double *array2)
{
  static const char *func = "mastrfu";
  int nlev;
  int nlat;
  int ilev, ilat, n;
  double *plevel, *phi, *cosphi, *dummy;
  double fact =  4*atan(1.0) * 6371000 / 9.81;
  double **field1, **field2;

  nlat = gridInqSize(gridID);
  nlev = zaxisInqSize(zaxisID);
  phi    = (double *) malloc(nlat*sizeof(double));
  dummy  = (double *) malloc(nlat*sizeof(double));
  cosphi = (double *) malloc(nlat*sizeof(double));
  plevel = (double *) malloc(nlev*sizeof(double));
  field1 = (double **) malloc(nlev*sizeof(double*));
  field2 = (double **) malloc(nlev*sizeof(double*));

  zaxisInqLevels(zaxisID, plevel);

  gaussaw(phi, dummy, nlat);

  for ( ilat = 0; ilat < nlat; ilat++ )
    cosphi[ilat] = sqrt(1.0 - phi[ilat]*phi[ilat]);

  for ( ilev = 0; ilev < nlev; ilev++ )
    {
      field1[ilev] = array1 + ilev*nlat;
      field2[ilev] = array2 + ilev*nlat;
    }

  for ( ilev = 0; ilev < nlev; ilev++ )
    for ( ilat = 0; ilat < nlat; ilat++ )
      field2[ilev][ilat] = 0.0;

  for ( ilev = nlev-1; ilev >= 0; ilev-- )
    for ( n = ilev; n < nlev-1; n++ )
      for ( ilat = 0; ilat < nlat; ilat++ )
	{
	  field2[ilev][ilat] = field2[ilev][ilat] +
	                       fact*(field1[n][ilat]+field1[n+1][ilat])*
	                       cosphi[ilat]*(plevel[n]-plevel[n+1]);
	}

  free(field2);
  free(field1);
  free(plevel);
  free(cosphi);
  free(dummy);
  free(phi);
}


void *Mastrfu(void *argument)
{
  static const char *func = "Mastrfu";
  int streamID1, streamID2;
  int nrecs;
  int tsID, recID, varID, levelID;
  int gridsize;
  int nvars, code, gridID, zaxisID, nlev;
  int vlistID1, vlistID2;
  int offset;
  int nmiss;
  double *array1, *array2;
  int taxisID1, taxisID2;

  cdoInitialize(argument);

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);

  nvars = vlistNvars(vlistID1);
  if ( nvars != 1 )
    cdoAbort("This operator works only with one variable!");

  code = vlistInqVarCode(vlistID1, 0);
  if ( code != 132 )
    cdoWarning("Unexpected code %d!", code);

  zaxisID = vlistInqVarZaxis(vlistID1, 0);
  if ( zaxisInqType(zaxisID) != ZAXIS_PRESSURE &&
       zaxisInqType(zaxisID) != ZAXIS_GENERIC )
    {
      char longname[128];
      zaxisInqLongname(zaxisID, longname);
      cdoWarning("Unexpected vertical grid %s!", longname);
    }

  gridID = vlistInqVarGrid(vlistID1, 0);
  if ( gridInqXsize(gridID) > 1 )
    cdoAbort("Grid must be a zonal mean!");

  gridsize = gridInqSize(gridID);
  nlev = zaxisInqSize(zaxisID);

  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  vlistDefVarCode(vlistID2, 0, 272);
  vlistDefVarName(vlistID2, 0, "mastrfu");
  vlistDefVarLongname(vlistID2, 0, "mass stream function");
  vlistDefVarUnits(vlistID2, 0, "kg/s");

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  array1 = (double *) malloc(gridsize*nlev*sizeof(double));
  array2 = (double *) malloc(gridsize*nlev*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  offset  = gridsize*levelID;
	  streamReadRecord(streamID1, array1+offset, &nmiss);
	  if ( nmiss ) cdoAbort("missing values unsupported for this operator!");
	}

      mastrfu(gridID, zaxisID, array1, array2);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  varID = 0;
	  levelID = recID;
	  streamDefRecord(streamID2, varID,  levelID);
	  offset  = gridsize*levelID;
	  streamWriteRecord(streamID2, array2+offset, nmiss);     
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( array1 ) free(array1);
  if ( array2 ) free(array2);

  cdoFinish();

  return (0);
}

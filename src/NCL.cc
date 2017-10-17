/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2017 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include <limits.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "libncl.h"


void uv2dv_cfd_W(double *u, double *v, double *lon, double *lat, size_t nlon, size_t nlat, size_t nlev, int iopt, double *div)
{
  int ierror;
  int bound_opt = iopt;
  int has_missing_u;
  double missing_u = 0, missing_du = 0, missing_ru = 0;

  size_t nlatnlon = nlat * nlon;

  // Test dimension sizes.
  if ( (nlon > INT_MAX) || (nlat > INT_MAX) )
    cdoAbort("nlat and/or nlon is greater than INT_MAX!");

  int inlon = (int) nlon;
  int inlat = (int) nlat;

  size_t gridsize_uv = nlatnlon;

  // Check for missing values.
  // coerce_missing(type_u,has_missing_u,&missing_u,&missing_du,&missing_ru);

  for ( size_t k = 0; k < nlev; ++k )
    {
      printf("level: %zu\n", k+1);
      double *tmp_u = u + k*gridsize_uv;
      double *tmp_v = v + k*gridsize_uv;
      double *tmp_div = div + k*gridsize_uv;
      // Init output array.
      memset(tmp_div, 0, gridsize_uv*sizeof(double));
      // Call the Fortran routine.
#ifdef HAVE_CF_INTERFACE
      DDVFIDF(tmp_u, tmp_v, lat, lon, inlon, inlat, missing_du, bound_opt, tmp_div, ierror);
#else
      cdoAbort("Fortran support not compiled in!");
#endif
    }
}


void *NCL(void *argument)
{
  int iopt = 1;
  int nrecs;
  int varID, levelID;
  size_t nmiss, nmissu, nmissv;

  cdoInitialize(argument);

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int varIDu = 1;
  int varIDv = 2;

  int gridIDu = vlistInqVarGrid(vlistID1, varIDu);
  int gridIDv = vlistInqVarGrid(vlistID1, varIDv);
  int gridtype = gridInqType(gridIDu);
  size_t gridsizeuv = gridInqSize(gridIDu);

  if ( !((gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN) && gridtype == gridInqType(gridIDv)) )
    cdoAbort("u and v must be on a regular lonlat or Gaussian grid!");
  
  if ( !(gridsizeuv == gridInqSize(gridIDv)) )
    cdoAbort("u and v must have the same grid size!");

  size_t nlon = gridInqXsize(gridIDu);
  size_t nlat = gridInqYsize(gridIDu);

  double *lon = (double*) Malloc(nlon*sizeof(double));
  double *lat = (double*) Malloc(nlat*sizeof(double));

  gridInqXvals(gridIDu, lon);
  gridInqYvals(gridIDu, lat);

  size_t nlev = zaxisInqSize(vlistInqVarZaxis(vlistID1, varIDu));

  if ( nlev != (size_t)zaxisInqSize(vlistInqVarZaxis(vlistID1, varIDv)) )
    cdoAbort("u and v must have the same number of level!");

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());

  pstreamDefVlist(streamID2, vlistID2);

  size_t gridsizemax = vlistGridsizeMax(vlistID1);
  double *array = (double*) Malloc(gridsizemax*sizeof(double));
  // level missing
  double *arrayu = (double*) Malloc(nlev*gridsizeuv*sizeof(double));
  double *arrayv = (double*) Malloc(nlev*gridsizeuv*sizeof(double));
  double *arrayd = (double*) Malloc(nlev*gridsizeuv*sizeof(double));

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);
	       
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);         
          pstreamReadRecord(streamID1, array, &nmiss);

          if ( varID == varIDu || varID == varIDv )
            {
              if ( varID == varIDu ) { memcpy(arrayu+levelID*gridsizeuv, array, gridsizeuv*sizeof(double)); nmissu = nmiss; }
              if ( varID == varIDv ) { memcpy(arrayv+levelID*gridsizeuv, array, gridsizeuv*sizeof(double)); nmissv = nmiss; }
            }
          else
            {
              pstreamDefRecord(streamID2,  varID,  levelID);
	      pstreamWriteRecord(streamID2, array, nmiss);
	    }
	}

      uv2dv_cfd_W(arrayu, arrayv, lon, lat, nlon, nlat, nlev, iopt, arrayd);

      for ( levelID = 0; levelID < nlev; ++levelID )
        {
          nmiss = 0;
          pstreamDefRecord(streamID2,  varIDu,  levelID);
          pstreamWriteRecord(streamID2, arrayd+levelID*gridsizeuv, nmiss);
        }

      for ( levelID = 0; levelID < nlev; ++levelID )
        {
          pstreamDefRecord(streamID2,  varIDv,  levelID);
          pstreamWriteRecord(streamID2, arrayv+levelID*gridsizeuv, nmiss);
        }

      tsID++;
    }

  pstreamClose(streamID1);
  pstreamClose(streamID2);

  vlistDestroy(vlistID2);

  if ( array ) Free(array);
  if ( arrayu ) Free(arrayu);
  if ( arrayv ) Free(arrayv);
  if ( arrayd ) Free(arrayd);

  cdoFinish();

  return 0;
}

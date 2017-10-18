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


void uv2dv_cfd_W(double missval, double *u, double *v, double *lon, double *lat, size_t nlon, size_t nlat, size_t nlev, int boundOpt, double *div)
{
  int ierror;
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
      double *tmp_u = u + k*gridsize_uv;
      double *tmp_v = v + k*gridsize_uv;
      double *tmp_div = div + k*gridsize_uv;
      // Init output array.
      memset(tmp_div, 0, gridsize_uv*sizeof(double));
      // Call the Fortran routine.
#ifdef HAVE_CF_INTERFACE
      DDVFIDF(tmp_u, tmp_v, lat, lon, inlon, inlat, missval, boundOpt, tmp_div, ierror);
#else
      cdoAbort("Fortran support not compiled in!");
#endif
    }
}

int find_name(int vlistID, char *name)
{
  char varname[CDI_MAX_NAME];

  int nvars = vlistNvars(vlistID);
  for ( int varID = 0; varID < nvars; ++varID )
    {
      vlistInqVarName(vlistID, varID, varname);
      if ( STR_IS_EQ(name, varname) ) return varID;
    }

  return CDI_UNDEFID;
}

static
void set_parameter(char *name_u, char *name_v, int *boundOpt)
{
  name_u[0] = 0;
  name_v[0] = 0;

  int pargc = operatorArgc();
  if ( pargc )
    { 
      char **pargv = operatorArgv();

      list_t *kvlist = list_new(sizeof(keyValues_t *), free_keyval, "PARAMETER");
      if ( kvlist_parse_cmdline(kvlist, pargc, pargv) != 0 ) cdoAbort("Parse error!");
      if ( cdoVerbose ) kvlist_print(kvlist);

      for ( listNode_t *kvnode = kvlist->head; kvnode; kvnode = kvnode->next )
        {
          keyValues_t *kv = *(keyValues_t **)kvnode->data;
          const char *key = kv->key;
          if ( kv->nvalues > 1 ) cdoAbort("Too many values for parameter key >%s<!", key);
          if ( kv->nvalues < 1 ) cdoAbort("Missing value for parameter key >%s<!", key);
          const char *value = kv->values[0];
          
          if      ( STR_IS_EQ(key, "u") ) strcpy(name_u, value);
          else if ( STR_IS_EQ(key, "v") ) strcpy(name_v, value);
          else if ( STR_IS_EQ(key, "boundOpt") ) *boundOpt = parameter2int(value);
          else cdoAbort("Invalid parameter key >%s<!", key);
        }          
          
      list_destroy(kvlist);
    }
      
  if ( cdoVerbose )
    cdoPrint("u = %s, v = %s", name_u, name_v);
}


void *NCL(void *argument)
{
  int boundOpt = -1;
  int nrecs;
  int varID, levelID;
  size_t nmiss;
  char name_u[CDI_MAX_NAME], name_v[CDI_MAX_NAME];

  cdoInitialize(argument);

  set_parameter(name_u, name_v, &boundOpt);

  if ( !name_u[0] ) strcpy(name_u, "u");
  if ( !name_v[0] ) strcpy(name_v, "v");

  int streamID1 = pstreamOpenRead(cdoStreamName(0));

  int vlistID1 = pstreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int varIDu = find_name(vlistID1, name_u);
  int varIDv = find_name(vlistID1, name_v);

  if ( varIDu == CDI_UNDEFID ) cdoAbort("%s not found!", name_u);
  if ( varIDv == CDI_UNDEFID ) cdoAbort("%s not found!", name_v);

  int gridIDu = vlistInqVarGrid(vlistID1, varIDu);
  int gridIDv = vlistInqVarGrid(vlistID1, varIDv);
  int gridtype = gridInqType(gridIDu);
  size_t gridsizeuv = gridInqSize(gridIDu);

  if ( !((gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN) && gridtype == gridInqType(gridIDv)) )
    cdoAbort("u and v must be on a regular lonlat or Gaussian grid!");
  
  if ( !(gridsizeuv == gridInqSize(gridIDv)) )
    cdoAbort("u and v must have the same grid size!");

  if ( boundOpt == -1 ) boundOpt = gridIsCircular(gridIDu) ? 1 : 0;
  if ( boundOpt < 0 || boundOpt > 3 ) cdoAbort("Parameter boundOpt=%d out of bounds (0-3)!", boundOpt);

  size_t nlon = gridInqXsize(gridIDu);
  size_t nlat = gridInqYsize(gridIDu);

  int zaxisIDu = vlistInqVarZaxis(vlistID1, varIDu);
  int nlev = zaxisInqSize(zaxisIDu);

  if ( nlev != zaxisInqSize(vlistInqVarZaxis(vlistID1, varIDv)) )
    cdoAbort("u and v must have the same number of level!");


  double missvalu = vlistInqVarMissval(vlistID1, varIDu);
  double missvalv = vlistInqVarMissval(vlistID1, varIDv);

  
  int timetype = vlistInqVarTimetype(vlistID1, varIDu);
  int varIDdiv = vlistDefVar(vlistID2, gridIDu, zaxisIDu, timetype);
  vlistDefVarName(vlistID2, varIDdiv, "d");
  vlistDefVarLongname(vlistID2, varIDdiv, "divergence");
  vlistDefVarUnits(vlistID2, varIDdiv, "1/s");
  vlistDefVarMissval(vlistID2, varIDdiv, missvalu);
  
  double *lon = (double*) Malloc(nlon*sizeof(double));
  double *lat = (double*) Malloc(nlat*sizeof(double));

  gridInqXvals(gridIDu, lon);
  gridInqYvals(gridIDu, lat);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = pstreamOpenWrite(cdoStreamName(1), cdoFiletype());

  pstreamDefVlist(streamID2, vlistID2);

  size_t gridsizemax = vlistGridsizeMax(vlistID1);
  double *array = (double*) Malloc(gridsizemax*sizeof(double));
  double *arrayu = (double*) Malloc(nlev*gridsizeuv*sizeof(double));
  double *arrayv = (double*) Malloc(nlev*gridsizeuv*sizeof(double));
  double *arrayd = (double*) Malloc(nlev*gridsizeuv*sizeof(double));

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      size_t nmissu = 0, nmissv = 0;

      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);
	       
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);         
          pstreamReadRecord(streamID1, array, &nmiss);

          if ( varID == varIDu || varID == varIDv )
            {
              if ( varID == varIDu ) { memcpy(arrayu+levelID*gridsizeuv, array, gridsizeuv*sizeof(double)); nmissu += nmiss; }
              if ( varID == varIDv ) { memcpy(arrayv+levelID*gridsizeuv, array, gridsizeuv*sizeof(double)); nmissv += nmiss; }
            }

          pstreamDefRecord(streamID2,  varID,  levelID);
          pstreamWriteRecord(streamID2, array, nmiss);
	}

      if ( nmissu != nmissv )
        {
          cdoAbort("u and v have different number of missing values!");
          if ( nmissu > 0 && !DBL_IS_EQUAL(missvalu, missvalv) )
            {
              for ( levelID = 0; levelID < nlev; ++levelID )
                {
                  double *parray = arrayv+levelID*gridsizeuv;
                  for ( size_t i = 0; i < gridsizeuv; ++i )
                    if ( DBL_IS_EQUAL(parray[i], missvalv) ) parray[i] = missvalu;
                  
                }
            }
        }

      uv2dv_cfd_W(missvalu, arrayu, arrayv, lon, lat, nlon, nlat, nlev, boundOpt, arrayd);

      for ( levelID = 0; levelID < nlev; ++levelID )
        {
          double *parray = arrayd+levelID*gridsizeuv;
          nmiss = 0;
          for ( size_t i = 0; i < gridsizeuv; ++i )
            if ( DBL_IS_EQUAL(parray[i], missvalu) ) nmiss++;

          pstreamDefRecord(streamID2, varIDdiv, levelID);
          pstreamWriteRecord(streamID2, parray, nmiss);
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

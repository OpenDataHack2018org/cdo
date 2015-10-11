#include <stdio.h>
#include <string.h>

#include <cdi.h>
#include <cdo.h>
#include "dmemory.h"
#include "field.h"
#include "util.h"


void field_init(field_t *field)
{
  memset(field, 0, sizeof(field_t));
}


field_t **field_allocate(int vlistID, int ptype, int init)
{
  int nvars = vlistNvars(vlistID);

  field_t **field = (field_t **) Malloc(nvars*sizeof(field_t *));

  for ( int varID = 0; varID < nvars; ++varID )
    {
      int nwpv     = vlistInqNWPV(vlistID, varID); // number of words per value; real:1  complex:2
      int gridID   = vlistInqVarGrid(vlistID, varID);
      int gridsize = gridInqSize(gridID);
      int zaxisID  = vlistInqVarZaxis(vlistID, varID);
      int nlevel   = zaxisInqSize(zaxisID);
      double missval  = vlistInqVarMissval(vlistID, varID);

      field[varID] = (field_t*) Malloc(nlevel*sizeof(field_t));

      for ( int levelID = 0; levelID < nlevel; ++levelID )
	{
	  field_init(&field[varID][levelID]);

	  field[varID][levelID].nwpv    = nwpv;
	  field[varID][levelID].grid    = gridID;
	  field[varID][levelID].nsamp   = 0;
	  field[varID][levelID].nmiss   = 0;
	  field[varID][levelID].nmiss2  = 0;
	  field[varID][levelID].missval = missval;
	  field[varID][levelID].ptr     = NULL;
	  field[varID][levelID].ptr2    = NULL;
	  field[varID][levelID].weight  = NULL;

	  if ( ptype & FIELD_PTR )
	    {
	      field[varID][levelID].ptr = (double*) Malloc(nwpv*gridsize*sizeof(double));
	      if ( init ) memset(field[varID][levelID].ptr, 0, nwpv*gridsize*sizeof(double));
	    }

	  if ( ptype & FIELD_PTR2 )
	    {
	      field[varID][levelID].ptr2 = (double*) Malloc(nwpv*gridsize*sizeof(double));
	      if ( init ) memset(field[varID][levelID].ptr2, 0, nwpv*gridsize*sizeof(double));
	    }

	  if ( ptype & FIELD_WGT )
	    {
	      field[varID][levelID].weight = (double*) Malloc(nwpv*gridsize*sizeof(double));
	      if ( init ) memset(field[varID][levelID].weight, 0, nwpv*gridsize*sizeof(double));
	    }    
	}
    }

  return field;
}


field_t **field_malloc(int vlistID, int ptype)
{
  return field_allocate(vlistID, ptype, 0);
}


field_t **field_calloc(int vlistID, int ptype)
{
  return field_allocate(vlistID, ptype, 1);
}


void field_free(field_t **field, int vlistID)
{
  int nvars = vlistNvars(vlistID);
  for ( int varID = 0; varID < nvars; ++varID )
    {
      int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID, varID));
      for ( int levelID = 0; levelID < nlevel; ++levelID )
	{
	  if ( field[varID][levelID].ptr )    Free(field[varID][levelID].ptr);
	  if ( field[varID][levelID].ptr2 )   Free(field[varID][levelID].ptr2);
       	  if ( field[varID][levelID].weight ) Free(field[varID][levelID].weight);
	}

      Free(field[varID]);
    }

  Free(field);
}

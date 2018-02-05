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
#include <cdi.h>
#include "cdo_int.h"


bool *cdo_read_timestepmask(const char *maskfile, int *n)
{
  *n = 0;

  int streamID = streamOpenRead(maskfile);
  if ( streamID == CDI_UNDEFID ) cdoAbort("Open failed on %s!", maskfile);

  int vlistID = streamInqVlist(streamID);

  int nvars = vlistNvars(vlistID);
  if ( nvars > 1 ) cdoAbort("timestepmask %s contains more than one variable!", maskfile);

  size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID, 0));
  if ( gridsize > 1 ) cdoAbort("timestepmask %s has more than one gridpoint!", maskfile);

  int nlev = zaxisInqSize(vlistInqVarZaxis(vlistID, 0));
  if ( nlev > 1 ) cdoAbort("timestepmask %s has more than one level!", maskfile);

  int nts = vlistNtsteps(vlistID);
  if ( nts == -1 )
    {
      nts = 0;
      while ( streamInqTimestep(streamID, nts) ) nts++;

      if ( cdoVerbose ) cdoPrint("%s: counted %i timeSteps in %s", __func__, nts, maskfile);

      streamClose(streamID);
      streamID = streamOpenRead(maskfile);

    }
  else
    if ( cdoVerbose ) cdoPrint("%s: found %i timeSteps in %s", __func__, nts, maskfile);

  *n = nts;
  bool *imask = (bool*) Malloc(nts*sizeof(bool));

  int nrecs;
  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID, tsID)) )
    {
      if ( nrecs != 1 ) cdoAbort("Internal error; unexprected number of records!");

      int varID, levelID;
      size_t nmiss;
      double value;
      streamInqRecord(streamID, &varID, &levelID);
      streamReadRecord(streamID, &value, &nmiss);
      
      imask[tsID] = !(nmiss || IS_EQUAL(value, 0));
      
      tsID++;  
    }
  
  streamClose(streamID);

  return imask;
}


bool *cdo_read_mask(const char *maskfile, size_t *n)
{
  *n = 0;

  int streamID = streamOpenRead(maskfile);
  if ( streamID == CDI_UNDEFID ) cdoAbort("Open failed on %s!", maskfile);

  int vlistID = streamInqVlist(streamID);

  int nvars = vlistNvars(vlistID);
  if ( nvars > 1 ) cdoAbort("Mask %s contains more than one variable!", maskfile);

  size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID, 0));

  int nlev = zaxisInqSize(vlistInqVarZaxis(vlistID, 0));
  if ( nlev > 1 ) cdoAbort("Mask %s has more than one level!", maskfile);

  *n = gridsize;
  bool *imask = (bool*) Malloc(gridsize*sizeof(bool));
  double *dmask = (double*) Malloc(gridsize*sizeof(double));

  int nrecs = streamInqTimestep(streamID, 0);
  if ( nrecs != 1 ) cdoAbort("Internal error; unexprected number of records!");

  int varID, levelID;
  size_t nmiss;
  streamInqRecord(streamID, &varID, &levelID);
  streamReadRecord(streamID, dmask, &nmiss);

  for ( size_t i = 0; i < gridsize; ++i )
    imask[i] = IS_NOT_EQUAL(dmask[i], 0);
      
      Free(dmask);
  
  streamClose(streamID);

  return imask;
}

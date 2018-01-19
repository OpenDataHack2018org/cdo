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
#include "pstream_int.h"
#include "cdoDebugOutput.h"

    void
pstreamClose(int pstreamID)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);
  if (pstreamptr == NULL)
    ERROR("Internal problem, stream ", pstreamID ," not open!");

  pstreamptr->close();
}
int
pstreamInqVlist(int pstreamID)
{
  if(CdoDebug::PSTREAM) MESSAGE("Inquiring Vlist from pstream ", pstreamID);
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);
  int vlistID =  pstreamptr->inqVlist();
  processDefVarNum(vlistNvars(vlistID));
  return vlistID;
}

void
pstreamDefVlist(int pstreamID, int vlistID)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);
  pstreamptr->defVlist(vlistID);
}

int
pstreamInqTimestep(int pstreamID, int tsID)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);
  int nrecs = -1;
#ifdef  HAVE_LIBPTHREAD
  if(pstreamptr->isPipe())
  {
      if (CdoDebug::PSTREAM){
          MESSAGE(pstreamptr->pipe->name.c_str(), " pstreamID ",  pstreamptr->self);}
      nrecs = pstreamptr->pipe->pipeInqTimestep(tsID);
  }
  else
#endif
  {
      if (pstreamptr->mfiles){tsID -= pstreamptr->tsID0;}
      nrecs = pstreamptr->inqTimestep(tsID);
      if (nrecs && tsID != pstreamptr->tsID){
          processDefTimesteps(pstreamID);
          pstreamptr->tsID = tsID;
        }
  }
  return nrecs;
}

void
pstreamDefTimestep(int pstreamID, int tsID)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);
  pstreamptr->defTimestep(tsID);
}


int
pstreamInqFiletype(int pstreamID)
{
  return pstream_to_pointer(pstreamID)->inqFileType();
}

void
pstreamInqGRIBinfo(int pstreamID, int *intnum, float *fltnum, off_t *bignum)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

  streamInqGRIBinfo(pstreamptr->m_fileID, intnum, fltnum, bignum);
}

int
pstreamFileID(int pstreamID)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

  return pstreamptr->m_fileID;
}

int
pstreamInqByteorder(int pstreamID)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);
  return pstreamptr->inqByteorder();
}



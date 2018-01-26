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
#include "process_int.h"
#include "exception.h"

    void
pstreamClose(int pstreamID)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);
  if (pstreamptr == NULL)
    ERROR("Internal problem, stream ", pstreamID ," not open!");

  processSelf().addNvals(pstreamptr->getNvals());
  pstreamptr->close();
}
int
pstreamInqVlist(int pstreamID)
{
  if(CdoDebug::PSTREAM) MESSAGE("Inquiring Vlist from pstream ", pstreamID);
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);

  int vlistID =  pstreamptr->inqVlist();

  if (vlistNumber(vlistID) == CDI_COMP && cdoStreamNumber() == CDI_REAL)
    cdoAbort("Complex fields are not supported by this operator!");

   if (vlistNumber(vlistID) == CDI_REAL && cdoStreamNumber() == CDI_COMP)
    cdoAbort("This operator needs complex fields!");

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

void
pstreamReadRecord(int pstreamID, double *data, size_t *nmiss)
{
  if (data == NULL)
    cdoAbort("Data pointer not allocated (pstreamReadRecord)!");

  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);
  pstreamptr->readRecord(data, nmiss);
}

void
pstreamReadRecordF(int pstreamID, float *data, size_t *nmiss)
{
  if (data == NULL)
    cdoAbort("Data pointer not allocated (pstreamReadRecord)!");

  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);
  pstreamptr->readRecordF(data, nmiss);
}

void
pstreamCopyRecord(int pstreamIDdest, int pstreamIDsrc)
{
  if (CdoDebug::PSTREAM)
    MESSAGE("pstreamIDdest = ",pstreamIDdest,"  pstreamIDsrc = ", pstreamIDsrc);

  pstream_t *pstreamptr_dest = pstream_to_pointer(pstreamIDdest);

  pstreamptr_dest->copyRecord(pstream_to_pointer(pstreamIDsrc));
}

int
pstreamInqRecord(int pstreamID, int *varID, int *levelID)
{
  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);
  return pstreamptr->inqRecord(varID, levelID);
}

void
pstreamWriteRecordF(int pstreamID, float *data, size_t nmiss)
{
  if (data == NULL)
    cdoAbort("Data pointer not allocated (%s)!", __func__);

  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);
  pstreamptr->writeRecordF(data,nmiss);
}

void
pstreamWriteRecord(int pstreamID, double *data, size_t nmiss)
{
  if (data == NULL)
    cdoAbort("Data pointer not allocated (%s)!", __func__);

  pstream_t *pstreamptr = pstream_to_pointer(pstreamID);
  pstreamptr->writeRecord(data,nmiss);
}

void
pstreamDefRecord(int pstreamID, int varID, int levelID)
{
  pstream_t *pstreamptr;

  pstreamptr = pstream_to_pointer(pstreamID);
  pstreamptr->defRecord(varID, levelID);
}


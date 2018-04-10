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

#ifdef HAVE_CONFIG_H
#include "config.h" /* _FILE_OFFSET_BITS influence off_t */
#endif

#include <cdi.h>
#include "pstream_int.h"
#include "cdoDebugOutput.h"
#include "process_int.h"
#include "exception.h"

void
pstreamClose(int pstreamID)
{
  PstreamType *pstreamptr = pstreamToPointer(pstreamID);
  if (pstreamptr == NULL) ERROR("Internal problem, stream ", pstreamID, " not open!");

  Cdo_Debug(CdoDebug::PSTREAM, "Adding ", pstreamptr->getNvals(), " to pstream ", pstreamptr->self, " ", pstreamptr->m_name);
  if (pstreamptr->rthreadID == pthread_self()) processSelf().addNvals(pstreamptr->getNvals());
  pstreamptr->close();
}

int
pstreamInqVlist(int pstreamID)
{
  if (CdoDebug::PSTREAM) MESSAGE("Inquiring Vlist from pstream ", pstreamID);
  PstreamType *pstreamptr = pstreamToPointer(pstreamID);
  int vlistID = pstreamptr->inqVlist();
  return vlistID;
}

void
pstreamDefVlist(int pstreamID, int vlistID)
{
  PstreamType *pstreamptr = pstreamToPointer(pstreamID);
  pstreamptr->defVlist(vlistID);
}

int
pstreamInqTimestep(PstreamType *p_pstreamptr, int tsID)
{
  int nrecs = -1;
#ifdef HAVE_LIBPTHREAD
  if (p_pstreamptr->isPipe())
    {
      if (CdoDebug::PSTREAM)
        {
          MESSAGE(p_pstreamptr->pipe->name.c_str(), " pstreamID ", p_pstreamptr->self);
        }
      nrecs = p_pstreamptr->pipe->pipeInqTimestep(tsID);
    }
  else
#endif
    {
      if (p_pstreamptr->mfiles)
        {
          tsID -= p_pstreamptr->tsID0;
        }
      nrecs = p_pstreamptr->inqTimestep(tsID);
      if (nrecs && tsID != p_pstreamptr->tsID)
        {
          p_pstreamptr->tsID = tsID;
        }
    }
  return nrecs;
}

void
pstreamDefTimestep(int pstreamID, int tsID)
{
  PstreamType *pstreamptr = pstreamToPointer(pstreamID);
  pstreamptr->defTimestep(tsID);
}

int
pstreamInqFiletype(int pstreamID)
{
  return pstreamToPointer(pstreamID)->inqFileType();
}

void
pstreamInqGRIBinfo(int pstreamID, int *intnum, float *fltnum, off_t *bignum)
{
  PstreamType *pstreamptr = pstreamToPointer(pstreamID);

  streamInqGRIBinfo(pstreamptr->m_fileID, intnum, fltnum, bignum);
}

int
pstreamFileID(int pstreamID)
{
  PstreamType *pstreamptr = pstreamToPointer(pstreamID);

  return pstreamptr->m_fileID;
}

int
pstreamInqByteorder(int pstreamID)
{
  PstreamType *pstreamptr = pstreamToPointer(pstreamID);
  return pstreamptr->inqByteorder();
}

void
pstreamReadRecord(int pstreamID, double *data, size_t *nmiss)
{
  if (data == NULL) cdoAbort("Data pointer not allocated (pstreamReadRecord)!");

  PstreamType *pstreamptr = pstreamToPointer(pstreamID);
  pstreamptr->readRecord(data, nmiss);
}

void
pstreamReadRecordF(int pstreamID, float *data, size_t *nmiss)
{
  if (data == NULL) cdoAbort("Data pointer not allocated (pstreamReadRecord)!");

  PstreamType *pstreamptr = pstreamToPointer(pstreamID);
  pstreamptr->readRecordF(data, nmiss);
}

void
pstreamCopyRecord(int pstreamIDdest, int pstreamIDsrc)
{
  if (CdoDebug::PSTREAM) MESSAGE("pstreamIDdest = ", pstreamIDdest, "  pstreamIDsrc = ", pstreamIDsrc);

  PstreamType *pstreamptr_dest = pstreamToPointer(pstreamIDdest);

  pstreamptr_dest->copyRecord(pstreamToPointer(pstreamIDsrc));
}

void
pstreamInqRecord(int pstreamID, int *varID, int *levelID)
{
  PstreamType *pstreamptr = pstreamToPointer(pstreamID);
  pstreamptr->inqRecord(varID, levelID);
}

void
pstreamWriteRecordF(int pstreamID, float *data, size_t nmiss)
{
  if (data == NULL) cdoAbort("Data pointer not allocated (%s)!", __func__);

  PstreamType *pstreamptr = pstreamToPointer(pstreamID);
  pstreamptr->writeRecordF(data, nmiss);
}

void
pstreamWriteRecord(int pstreamID, double *data, size_t nmiss)
{
  if (data == NULL) cdoAbort("Data pointer not allocated (%s)!", __func__);

  PstreamType *pstreamptr = pstreamToPointer(pstreamID);
  pstreamptr->writeRecord(data, nmiss);
}

void
pstreamDefRecord(int pstreamID, int varID, int levelID)
{
  PstreamType *pstreamptr;

  pstreamptr = pstreamToPointer(pstreamID);
  pstreamptr->defRecord(varID, levelID);
}

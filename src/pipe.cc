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

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600 /* struct timespec */
#endif

#include <stdio.h>
#include <string.h>
#include <time.h>  // time()
#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "error.h"
#include "dmemory.h"
#include "pipe.h"
#include "pstream_int.h"

#if defined(HAVE_LIBPTHREAD)

static int PipeDebug = 0;

static void
pipe_init(pipe_t *pipe)
{
  pthread_mutexattr_t m_attr;
  pthread_condattr_t c_attr;

  pthread_mutexattr_init(&m_attr);
  pthread_condattr_init(&c_attr);
  /*
#if defined(_POSIX_THREAD_PROCESS_SHARED)
  if ( PipeDebug )
    {
      Message("setpshared mutexattr to PTHREAD_PROCESS_SHARED");
      Message("setpshared condattr to PTHREAD_PROCESS_SHARED");
    }

  pthread_mutexattr_setpshared(&m_attr, PTHREAD_PROCESS_SHARED);
  pthread_condattr_setpshared(&c_attr, PTHREAD_PROCESS_SHARED);

  if ( PipeDebug )
    {
      int pshared;
      pthread_mutexattr_getpshared(&m_attr, &pshared);
      if ( pshared == PTHREAD_PROCESS_SHARED )
        Message("getpshared mutexattr is PTHREAD_PROCESS_SHARED");
      else if ( pshared == PTHREAD_PROCESS_PRIVATE )
        Message("getpshared mutexattr is PTHREAD_PROCESS_PRIVATE");

      pthread_condattr_getpshared(&c_attr, &pshared);
      if ( pshared == PTHREAD_PROCESS_SHARED )
        Message("getpshared condattr is PTHREAD_PROCESS_SHARED");
      else if ( pshared == PTHREAD_PROCESS_PRIVATE )
        Message("getpshared condattr is PTHREAD_PROCESS_PRIVATE");
    }
#else
  if ( PipeDebug )
    Message("_POSIX_THREAD_PROCESS_SHARED undefined");
#endif
  */
  pipe->EOP = false;

  pipe->recIDr = -1;
  pipe->recIDw = -1;
  pipe->tsIDr = -1;
  pipe->tsIDw = -1;

  pipe->nvals = 0;
  pipe->nmiss = 0;
  pipe->data = NULL;
  pipe->hasdata = 0;
  pipe->usedata = true;
  pipe->pstreamptr_in = 0;

  pipe->mutex = (pthread_mutex_t *) Malloc(sizeof(pthread_mutex_t));
  pthread_mutex_init(pipe->mutex, &m_attr);

  pipe->tsDef = (pthread_cond_t *) Malloc(sizeof(pthread_cond_t));
  pthread_cond_init(pipe->tsDef, &c_attr);
  pipe->tsInq = (pthread_cond_t *) Malloc(sizeof(pthread_cond_t));
  pthread_cond_init(pipe->tsInq, &c_attr);

  pipe->recDef = (pthread_cond_t *) Malloc(sizeof(pthread_cond_t));
  pthread_cond_init(pipe->recDef, &c_attr);
  pipe->recInq = (pthread_cond_t *) Malloc(sizeof(pthread_cond_t));
  pthread_cond_init(pipe->recInq, &c_attr);

  pipe->vlistDef = (pthread_cond_t *) Malloc(sizeof(pthread_cond_t));
  pthread_cond_init(pipe->vlistDef, &c_attr);
  pipe->isclosed = (pthread_cond_t *) Malloc(sizeof(pthread_cond_t));
  pthread_cond_init(pipe->isclosed, &c_attr);

  pipe->writeCond = (pthread_cond_t *) Malloc(sizeof(pthread_cond_t));
  pthread_cond_init(pipe->writeCond, &c_attr);

  pipe->readCond = (pthread_cond_t *) Malloc(sizeof(pthread_cond_t));
  pthread_cond_init(pipe->readCond, &c_attr);

  pthread_mutexattr_destroy(&m_attr);
  pthread_condattr_destroy(&c_attr);
}

pipe_t *
pipeNew()
{
  pipe_t *pipe = (pipe_t *) Malloc(sizeof(pipe_t));

  pipe_init(pipe);

  return pipe;
}

void
pipeDelete(pipe_t *pipe)
{
  if (pipe)
    {
      if (pipe->mutex)
        Free(pipe->mutex);
      if (pipe->tsDef)
        Free(pipe->tsDef);
      if (pipe->tsInq)
        Free(pipe->tsInq);
      if (pipe->recDef)
        Free(pipe->recDef);
      if (pipe->recInq)
        Free(pipe->recInq);
      if (pipe->vlistDef)
        Free(pipe->vlistDef);
      if (pipe->isclosed)
        Free(pipe->isclosed);
      if (pipe->writeCond)
        Free(pipe->writeCond);
      if (pipe->readCond)
        Free(pipe->readCond);
      Free(pipe);
    }
}

void
pipeDefVlist(pstream_t *pstreamptr, int vlistID)
{
  char *pname = pstreamptr->name;
  pipe_t *pipe = pstreamptr->pipe;

  if (PipeDebug)
    Message("%s pstreamID %d", pname, pstreamptr->self);

  // LOCK
  pthread_mutex_lock(pipe->mutex);
  pstreamptr->vlistID = vlistID;
  pthread_mutex_unlock(pipe->mutex);
  // UNLOCK

  pthread_cond_signal(pipe->vlistDef);
}

#define TIMEOUT 1  // wait 1 seconds
#define MIN_WAIT_CYCLES 10
#define MAX_WAIT_CYCLES 3600
int processNumsActive(void);

int
pipeInqVlist(pstream_t *pstreamptr)
{
  char *pname = pstreamptr->name;
  pipe_t *pipe = pstreamptr->pipe;
  int vlistID = -1;
  struct timespec time_to_wait;
  int retcode = 0;
  int nwaitcycles = 0;

  time_to_wait.tv_sec = 0;
  time_to_wait.tv_nsec = 0;

  if (PipeDebug)
    Message("%s pstreamID %d", pname, pstreamptr->self);

  // LOCK
  pthread_mutex_lock(pipe->mutex);
  time_to_wait.tv_sec = time(NULL);
  while (pstreamptr->vlistID == -1 && retcode == 0)
    {
      time_to_wait.tv_sec += TIMEOUT;
      // fprintf(stderr, "tvsec %g\n", (double) time_to_wait.tv_sec);
      if (PipeDebug)
        Message("%s wait of vlistDef", pname);
      // pthread_cond_wait(pipe->vlistDef, pipe->mutex);
      retcode = pthread_cond_timedwait(pipe->vlistDef, pipe->mutex, &time_to_wait);
      // fprintf(stderr, "self %d retcode %d %d %d\n", pstreamptr->self, retcode, processNumsActive(),
      // pstreamptr->vlistID);
      if (retcode != 0 && nwaitcycles++ < MAX_WAIT_CYCLES)
        {
          if (processNumsActive() > 1 || (processNumsActive() == 1 && nwaitcycles < MIN_WAIT_CYCLES))
            retcode = 0;
        }
    }

  if (retcode == 0)
    vlistID = pstreamptr->vlistID;
  else if (PipeDebug)
    Message("%s timeout!", pname);

  pthread_mutex_unlock(pipe->mutex);
  // UNLOCK

  return vlistID;
}

int
pipeInqTimestep(pstream_t *pstreamptr, int tsID)
{
  char *pname = pstreamptr->name;
  pipe_t *pipe = pstreamptr->pipe;
  int nrecs;

  if (PipeDebug)
    Message("%s pstreamID %d", pname, pstreamptr->self);

  // LOCK
  pthread_mutex_lock(pipe->mutex);
  pipe->usedata = false;
  pipe->recIDr = -1;
  if (tsID != pipe->tsIDr + 1)
    {
      if (!(tsID == pipe->tsIDr && pipe->tsIDr == pipe->tsIDw && pipe->recIDr == -1))
        Error("%s unexpected tsID %d %d %d", pname, tsID, pipe->tsIDr + 1, pipe->tsIDw);
    }

  pipe->tsIDr = tsID;
  while (pipe->tsIDw != tsID)
    {
      if (pipe->EOP)
        {
          if (PipeDebug)
            Message("%s EOP", pname);
          break;
        }
      if (pipe->hasdata)
        {
          if (PipeDebug)
            Message("%s has data", pname);
          pipe->hasdata = 0;
          pipe->data = NULL;
          pthread_cond_signal(pipe->readCond);
        }
      else if (PipeDebug)
        Message("%s has no data", pname);

      pthread_cond_signal(pipe->recInq); /* o.k. ??? */

      if (PipeDebug)
        Message("%s wait of tsDef", pname);
      pthread_cond_wait(pipe->tsDef, pipe->mutex);
    }

  if (pipe->EOP)
    nrecs = 0;
  else
    nrecs = pipe->nrecs;

  pthread_mutex_unlock(pipe->mutex);
  // UNLOCK

  pthread_cond_signal(pipe->tsInq);

  return nrecs;
}

void
pipeDefTimestep(pstream_t *pstreamptr, int tsID)
{
  char *pname = pstreamptr->name;
  pipe_t *pipe = pstreamptr->pipe;
  int nrecs;

  if (PipeDebug)
    Message("%s pstreamID %d", pname, pstreamptr->self);

  // LOCK
  pthread_mutex_lock(pipe->mutex);
  pipe->recIDw = -1;
  pipe->tsIDw++;
  if (tsID != pipe->tsIDw)
    Error("unexpected tsID %d(%d) for %s", tsID, pipe->tsIDw, pname);

  if (tsID == 0)
    nrecs = vlistNrecs(pstreamptr->vlistID);
  else
    {
      int vlistID, varID;
      vlistID = pstreamptr->vlistID;
      nrecs = 0;
      for (varID = 0; varID < vlistNvars(vlistID); varID++)
        if (vlistInqVarTsteptype(vlistID, varID) != TSTEP_CONSTANT)
          nrecs += zaxisInqSize(vlistInqVarZaxis(vlistID, varID));
      // Message("nrecs = %d nvars = %d", nrecs, vlistNvars(vlistID));
    }

  pipe->nrecs = nrecs;
  if (PipeDebug)
    Message("%s nrecs %d tsID %d %d %d", pname, nrecs, tsID, pipe->tsIDw, pipe->tsIDr);
  if (nrecs == 0)
    pipe->EOP = true;
  pthread_mutex_unlock(pipe->mutex);
  // UNLOCK

  pthread_cond_signal(pipe->tsDef);
  // sleep(1);

  // LOCK
  pthread_mutex_lock(pipe->mutex);
  while (pipe->tsIDr < tsID)
    {
      if (pipe->EOP)
        {
          if (PipeDebug)
            Message("EOP");
          break;
        }
      if (PipeDebug)
        Message("%s wait of tsInq (tsID %d %d)", pname, tsID, pipe->tsIDr);
      pthread_cond_wait(pipe->tsInq, pipe->mutex);
    }
  pthread_mutex_unlock(pipe->mutex);
  // UNLOCK
}

int
pipeInqRecord(pstream_t *pstreamptr, int *varID, int *levelID)
{
  char *pname = pstreamptr->name;
  pipe_t *pipe = pstreamptr->pipe;
  bool condSignal = false;

  if (PipeDebug)
    Message("%s pstreamID %d", pname, pstreamptr->self);

  // LOCK
  pthread_mutex_lock(pipe->mutex);
  if (PipeDebug)
    Message("%s has no data %d %d", pname, pipe->recIDr, pipe->recIDw);
  if (pipe->hasdata || pipe->usedata)
    {
      pipe->hasdata = 0;
      pipe->data = NULL;
      pipe->usedata = false;
      condSignal = true;
    }
  pthread_mutex_unlock(pipe->mutex);
  // UNLOCK

  if (condSignal)
    pthread_cond_signal(pipe->readCond);

  // LOCK
  pthread_mutex_lock(pipe->mutex);
  pipe->usedata = true;
  pipe->recIDr++;

  if (PipeDebug)
    Message("%s recID %d %d", pname, pipe->recIDr, pipe->recIDw);

  while (pipe->recIDw != pipe->recIDr)
    {
      if (pipe->EOP)
        {
          if (PipeDebug)
            Message("EOP");
          break;
        }
      if (PipeDebug)
        Message("%s wait of recDef", pname);
      pthread_cond_wait(pipe->recDef, pipe->mutex);
    }

  if (pipe->EOP)
    {
      *varID = -1;
      *levelID = -1;
    }
  else
    {
      *varID = pipe->varID;
      *levelID = pipe->levelID;
    }

  pthread_mutex_unlock(pipe->mutex);
  // UNLOCK

  pthread_cond_signal(pipe->recInq);

  return 0;
}

void
pipeDefRecord(pstream_t *pstreamptr, int varID, int levelID)
{
  char *pname = pstreamptr->name;
  pipe_t *pipe = pstreamptr->pipe;
  bool condSignal = false;

  if (PipeDebug)
    Message("%s pstreamID %d", pname, pstreamptr->self);

  // LOCK
  pthread_mutex_lock(pipe->mutex);
  if (PipeDebug)
    Message("%s has data %d %d", pname, pipe->recIDr, pipe->recIDw);
  if (pipe->hasdata)
    {
      pipe->hasdata = 0;
      pipe->data = NULL;
      condSignal = true;
    }
  pthread_mutex_unlock(pipe->mutex);
  // UNLOCK

  if (condSignal)
    pthread_cond_signal(pipe->readCond);

  // LOCK
  pthread_mutex_lock(pipe->mutex);
  pipe->usedata = true;
  pipe->recIDw++;
  pipe->varID = varID;
  pipe->levelID = levelID;
  if (PipeDebug)
    Message("%s recID %d %d", pname, pipe->recIDr, pipe->recIDw);
  pthread_mutex_unlock(pipe->mutex);
  // UNLOCK

  pthread_cond_signal(pipe->recDef);

  // LOCK
  pthread_mutex_lock(pipe->mutex);
  while (pipe->recIDr < pipe->recIDw)
    {
      if (pipe->tsIDw != pipe->tsIDr)
        break;
      if (pipe->EOP)
        break;
      if (PipeDebug)
        Message("%s wait of recInq %d", pname, pipe->recIDr);
      pthread_cond_wait(pipe->recInq, pipe->mutex);
    }
  pthread_mutex_unlock(pipe->mutex);
  // UNLOCK
}

void
pipeCopyRecord(pstream_t *pstreamptr_out, pstream_t *pstreamptr_in)
{
  char *ipname = pstreamptr_in->name;
  char *opname = pstreamptr_out->name;
  pipe_t *pipe = pstreamptr_out->pipe;

  if (PipeDebug)
    Message("%s pstreamIDin %d", ipname, pstreamptr_in->self);
  if (PipeDebug)
    Message("%s pstreamIDout %d", opname, pstreamptr_out->self);

  // LOCK
  pthread_mutex_lock(pipe->mutex);
  pipe->hasdata = 2; /* pipe */
  pipe->pstreamptr_in = pstreamptr_in;
  pthread_mutex_unlock(pipe->mutex);
  // UNLOCK

  pthread_cond_signal(pipe->writeCond);

  // LOCK
  pthread_mutex_lock(pipe->mutex);
  while (pipe->hasdata)
    {
      if (!pipe->usedata)
        break;

      if (pipe->recIDw != pipe->recIDr)
        break;

      if (pipe->EOP)
        {
          if (PipeDebug)
            Message("EOP");
          break;
        }
      if (PipeDebug)
        Message("%s wait of readCond", opname);
      pthread_cond_wait(pipe->readCond, pipe->mutex);
    }
  pthread_mutex_unlock(pipe->mutex);
  // UNLOCK
}

/***
 * copys data from a pipe to data
 *
 * @param data destination for the record data
 * @param pipe pipe that has the wanted data
 */
void
pipeReadPipeRecord(pipe_t *pipe, double *data, char *pname, int vlistID, int *nmiss)
{
  int datasize;

  if (!pipe->data)
    Error("No data pointer for %s", pname);

  datasize = gridInqSize(vlistInqVarGrid(vlistID, pipe->varID));
  pipe->nvals += datasize;
  if (vlistNumber(vlistID) != CDI_REAL)
    datasize *= 2;
  memcpy(data, pipe->data, datasize * sizeof(double));
  *nmiss = pipe->nmiss;
}

void
pipeGetReadTarget(pstream_t *pstreamptr, pstream_t *pstreamptr_in, char *pname)
{

  pstreamptr_in = pstreamptr->pipe->pstreamptr_in;
  pstreamptr = pstreamptr_in;
  while (pstreamptr_in->ispipe)
    {
      if (PipeDebug)
        fprintf(stderr, "%s: istream %d is pipe\n", __func__, pstreamptr_in->self);
      pstreamptr = pstreamptr_in;
      pstreamptr_in = pstreamptr_in->pipe->pstreamptr_in;
      if (pstreamptr_in == 0)
        break;
    }

  if (pstreamptr->pipe->hasdata != 1)
    {
      Error("Internal problem! istream undefined");
    }
  else if (!pstreamptr->pipe->data)
    {
      Error("No data pointer for %s", pname);
    }
}

void
pipeReadRecord(pstream_t *pstreamptr, double *data, int *nmiss)
{
  char *pname = pstreamptr->name;
  pipe_t *pipe = pstreamptr->pipe;

  *nmiss = 0;
  if (PipeDebug)
    Message("%s pstreamID %d", pname, pstreamptr->self);

  // LOCK
  pthread_mutex_lock(pipe->mutex);
  while (pipe->hasdata == 0)
    {
      if (PipeDebug)
        Message("%s wait of writeCond", pname);
      pthread_cond_wait(pipe->writeCond, pipe->mutex);
    }

  if (pipe->hasdata == 2)
    {
      pstream_t *pstreamptr_in = 0;
      pipeGetReadTarget(pstreamptr, pstreamptr_in, pname);
      if (pstreamptr_in == 0)
        {
          if (PipeDebug)
            fprintf(stderr, "pstreamID = %d\n", pstreamptr->self);

          pipeReadPipeRecord(pstreamptr->pipe, data, pname, pstreamptr->vlistID, nmiss);
        }
      else
        {
          if (PipeDebug)
            fprintf(stderr, "%s: istream %d is file\n", __func__, pstreamptr_in->self);

          streamReadRecord(pstreamptr_in->fileID, data, nmiss);
        }
    }
  else if (pipe->hasdata == 1)
    {
      pipeReadPipeRecord(pipe, data, pname, pstreamptr->vlistID, nmiss);
    }
  else
    {
      Error("data type %d not implemented", pipe->hasdata);
    }

  if (PipeDebug)
    Message("%s read record %d", pname, pipe->recIDr);

  pipe->hasdata = 0;
  pipe->data = NULL;
  pthread_mutex_unlock(pipe->mutex);
  // UNLOCK

  pthread_cond_signal(pipe->readCond);
}

void
pipeWriteRecord(pstream_t *pstreamptr, double *data, int nmiss)
{
  char *pname = pstreamptr->name;
  pipe_t *pipe = pstreamptr->pipe;

  if (PipeDebug)
    Message("%s pstreamID %d", pname, pstreamptr->self);

  /*
  if ( ! pipe->usedata ) return;
  */
  // LOCK
  pthread_mutex_lock(pipe->mutex);
  pipe->hasdata = 1; /* data pointer */
  pipe->data = data;
  pipe->nmiss = nmiss;
  pthread_mutex_unlock(pipe->mutex);
  // UNLOCK

  pthread_cond_signal(pipe->writeCond);

  if (PipeDebug)
    Message("%s write record %d", pname, pipe->recIDw);

  // LOCK
  pthread_mutex_lock(pipe->mutex);
  while (pipe->hasdata)
    {
      if (!pipe->usedata)
        break;
      /*
      printf("ts ids %d %d\n", pipe->tsIDw, pipe->tsIDr);
      printf("rec ids %d %d\n", pipe->recIDw, pipe->recIDr);
      */
      if (pipe->recIDw != pipe->recIDr)
        break;

      if (pipe->EOP)
        {
          if (PipeDebug)
            Message("EOP");
          break;
        }
      if (PipeDebug)
        Message("%s wait of readCond", pname);
      pthread_cond_wait(pipe->readCond, pipe->mutex);
    }
  pthread_mutex_unlock(pipe->mutex);
  // UNLOCK
}

void
pipeDebug(int debug)
{
  PipeDebug = debug;
}

#endif

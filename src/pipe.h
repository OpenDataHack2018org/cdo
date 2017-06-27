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

#ifndef PIPE_H
#define PIPE_H

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <stdbool.h>
#include <sys/types.h>

#if defined(HAVE_LIBPTHREAD)

#include <pthread.h>
#include "pthread_debug.h"
#include "pstream.h"

#endif

#if defined(HAVE_LIBPTHREAD)

struct pipe_s {
  bool    EOP;
  bool    usedata;
  short   hasdata;
  int     nrecs;
  int     varID, levelID;
  int     recIDr, recIDw, tsIDr, tsIDw;
  int     nmiss;
  double *data;
  pstream_t *pstreamptr_in;
  /* unsigned long */ off_t nvals;
  pthread_mutex_t *mutex;
  pthread_cond_t *tsDef, *tsInq, *vlistDef, *isclosed;
  pthread_cond_t *recDef, *recInq;
  pthread_cond_t *writeCond, *readCond;
};

typedef struct pipe_s pipe_t;

pipe_t *pipeNew(void);
void    pipeDelete(pipe_t *pipe);

void  pipeDebug(int debug);

void  pipeDefVlist(pstream_t *pstreamptr, int vlistID);
int   pipeInqVlist(pstream_t *pstreamptr);

void  pipeDefTimestep(pstream_t *pstreamptr, int tsID);
int   pipeInqTimestep(pstream_t *pstreamptr, int tsID);

void  pipeDefRecord(pstream_t *pstreamptr, int  varID, int  levelID);
int   pipeInqRecord(pstream_t *pstreamptr, int *varID, int *levelID);

void  pipeReadRecord(pstream_t *pstreamptr, double *data, int *nmiss);
void  pipeWriteRecord(pstream_t *pstreamptr, double *data, int nmiss);
void  pipeCopyRecord(pstream_t *pstreamptr_dest, pstream_t *pstreamptr_src);

#endif

#endif  /* PIPE_H */

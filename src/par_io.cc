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
#include "config.h"
#endif

#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
#endif

#include "cdo_int.h"
#include "par_io.h"
#include "pstream_int.h"

void *
readRecord(void *arg)
{
  int streamID;
  int *varID, *levelID;
  size_t *nmiss;
  double *array;
  read_arg_t *read_arg = (read_arg_t *) arg;

  streamID = read_arg->streamID;
  varID = read_arg->varID;
  levelID = read_arg->levelID;
  nmiss = read_arg->nmiss;
  array = read_arg->array;

  /* fprintf(stderr, "pstreamInqRecord: streamID = %d\n", streamID); */
  pstreamInqRecord(streamID, varID, levelID);
  pstreamReadRecord(streamID, array, nmiss);
  /* fprintf(stderr, "pstreamReadRecord: varID %d levelID %d\n", *varID,
   * *levelID); */

  return (NULL);
}

void
parReadRecord(int streamID, int *varID, int *levelID, double *array, size_t *nmiss, par_io_t *parIO)
{
  int lpario = FALSE;
  int recID = 0, nrecs = 0;
#ifdef HAVE_LIBPTHREAD
  pthread_t thrID;
  /* pthread_attr_t attr; */
  int rval;
#endif

#ifdef HAVE_LIBPTHREAD
  if (parIO)
    {
      lpario = TRUE;
      recID = parIO->recID;
      nrecs = parIO->nrecs;
      thrID = parIO->thrID;
    }
#endif

  if (recID == 0 || lpario == FALSE)
    {
      read_arg_t read_arg;
      read_arg.streamID = streamID;
      read_arg.varID = varID;
      read_arg.levelID = levelID;
      read_arg.nmiss = nmiss;
      read_arg.array = array;

      readRecord(&read_arg);
    }
#ifdef HAVE_LIBPTHREAD
  else
    {
      /* fprintf(stderr, "parIO1: %ld streamID %d %d %d\n", (long)thrID,
       * streamID, recID, nrecs); */
      rval = pthread_join(thrID, NULL);
      if (rval != 0) cdoAbort("pthread_join failed!");

      *varID = parIO->varID;
      *levelID = parIO->levelID;
      *nmiss = parIO->nmiss;
      /* fprintf(stderr, "parIO2: %ld streamID %d %d %d\n", (long)thrID,
       * streamID, *varID, *levelID); */
      arrayCopy(parIO->array_size, parIO->array, array);
    }

  if (lpario && nrecs > 1)
    {
      read_arg_t *read_arg = &(parIO->read_arg);
      if ((recID + 1) < nrecs)
        {
          if (recID == 0)
            {
              pthread_attr_init(&parIO->attr);
              pthread_attr_setdetachstate(&parIO->attr, PTHREAD_CREATE_JOINABLE);
            }

          read_arg->streamID = streamID;
          read_arg->varID = &parIO->varID;
          read_arg->levelID = &parIO->levelID;
          read_arg->nmiss = &parIO->nmiss;
          read_arg->array = parIO->array;

          /* fprintf(stderr, "pthread_create: streamID %d %d\n",
           * read_arg->streamID,streamID); */
          rval = pthread_create(&thrID, &parIO->attr, readRecord, read_arg);
          if (rval != 0) cdoAbort("pthread_create failed!");

          /* fprintf(stderr, "thrID = %ld\n", (long) thrID); */
          parIO->thrID = thrID;
        }
      else
        pthread_attr_destroy(&parIO->attr);
    }
#endif
}

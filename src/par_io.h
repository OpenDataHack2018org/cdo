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
#ifndef PAR_IO_H
#define PAR_IO_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
#endif

struct read_arg_t
{
  int streamID;
  int *varID, *levelID;
  size_t *nmiss;
  double *array;
};

struct par_io_t
{
  int varID, levelID;
  size_t nmiss;
  double *array;
  int array_size;
  int recID, nrecs;
  read_arg_t read_arg;
#ifdef HAVE_LIBPTHREAD
  pthread_t thrID;
  pthread_attr_t attr;
#endif
};

void parReadRecord(int streamID, int *varID, int *levelID, double *array, size_t *nmiss, par_io_t *parIO);

#endif /* PAR_IO_H */

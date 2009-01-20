#ifndef _PAR_IO_H
#define _PAR_IO_H

#if  defined  (HAVE_CONFIG_H)
#  include "config.h"
#endif

#if  defined  (HAVE_LIBPTHREAD)
#  include <pthread.h>
#endif

typedef struct {
  int varID, levelID, nmiss;
  double *array;
  int array_size;
  int recID, nrecs;
#if  defined  (HAVE_LIBPTHREAD)
  pthread_t thrID;
  pthread_attr_t attr;
#endif
}
par_io_t;


void parReadRecord(int streamID, int *varID, int *levelID, double *array, int *nmiss, par_io_t *parIO);

#endif  /* _PAR_IO_H */

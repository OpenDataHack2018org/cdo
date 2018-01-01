#ifndef PAR_IO_H
#define PAR_IO_H

#ifdef  HAVE_CONFIG_H
#  include "config.h"
#endif

#ifdef  HAVE_LIBPTHREAD
#  include <pthread.h>
#endif


typedef struct {
  int streamID;
  int *varID, *levelID;
  size_t *nmiss;
  double *array;
}
read_arg_t;


typedef struct {
  int varID, levelID;
  size_t nmiss;
  double *array;
  int array_size;
  int recID, nrecs;
  read_arg_t read_arg;
#ifdef  HAVE_LIBPTHREAD
  pthread_t thrID;
  pthread_attr_t attr;
#endif
}
par_io_t;


void parReadRecord(int streamID, int *varID, int *levelID, double *array, size_t *nmiss, par_io_t *parIO);

#endif  /* PAR_IO_H */

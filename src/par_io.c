#if  defined  (HAVE_CONFIG_H)
#  include "config.h"
#endif

#if  defined  (HAVE_LIBPTHREAD)
#  include <pthread.h>
#endif

#include <string.h> /* memcpy */

#include "cdo.h"
#include "par_io.h"
#include "pstream.h"

typedef struct {
  int streamID;
  int *varID, *levelID, *nmiss;
  double *array;
}
READ_ARG;


void *readRecord(void *arg)
{
  int streamID;
  int *varID, *levelID, *nmiss;
  double *array;
  READ_ARG *read_arg = (READ_ARG *) arg;

  streamID = read_arg->streamID;
  varID    = read_arg->varID;
  levelID  = read_arg->levelID;
  nmiss    = read_arg->nmiss;
  array    = read_arg->array;

  streamInqRecord(streamID, varID, levelID);
  streamReadRecord(streamID, array, nmiss);
  fprintf(stderr, "1 varID %d levelID %d\n", *varID, *levelID);

  return (NULL);
}


void parReadRecord(int streamID, int *varID, int *levelID, double *array, int *nmiss, par_io_t *parIO)
{
  int lpario = FALSE;
  READ_ARG read_arg;
  void *statusp;
  int recID = 0, nrecs = 0;
#if  defined  (HAVE_LIBPTHREAD)
  pthread_t thrID;
  /* pthread_attr_t attr; */
  int rval;
#endif

  read_arg.streamID = streamID;
  read_arg.varID    = varID;
  read_arg.levelID  = levelID;
  read_arg.nmiss    = nmiss;
  read_arg.array    = array;

#if  defined  (HAVE_LIBPTHREAD)
  if ( parIO )
    {
      lpario = TRUE;
      recID = parIO->recID;
      nrecs = parIO->nrecs;
      thrID = parIO->thrID;
    }
#endif

  if ( recID == 0 || lpario == FALSE )
    {
      statusp = readRecord(&read_arg);
    }
#if  defined  (HAVE_LIBPTHREAD)
  else
    {
      fprintf(stderr, "parIO: %ld streamID %d %d %d\n", (long)thrID, streamID, recID, nrecs);
      rval = pthread_join(thrID, &statusp);
      if ( rval != 0 ) cdoAbort("pthread_join failed!");
      /*
      if ( *(int *)statusp < 0 )
	cdoAbort("pthread_join failed! (status = %d)", *(int *)statusp);
      */

      *varID    = parIO->varID;
      *levelID  = parIO->levelID;
      *nmiss    = parIO->nmiss;
      memcpy(array, parIO->array, parIO->array_size*sizeof(double));
    }

  if ( lpario && nrecs > 1 )
    {
      if ( (recID+1) < nrecs )
	{
	  if ( recID == 0 )
	    {
	      pthread_attr_init(&parIO->attr);
	      pthread_attr_setdetachstate(&parIO->attr, PTHREAD_CREATE_JOINABLE);
	    }

	  read_arg.streamID = streamID;
	  read_arg.varID    = &parIO->varID;
	  read_arg.levelID  = &parIO->levelID;
	  read_arg.nmiss    = &parIO->nmiss;
	  read_arg.array    = parIO->array;

	  rval = pthread_create(&thrID, &parIO->attr, readRecord, &read_arg);
	  if ( rval != 0 ) cdoAbort("pthread_create failed!");

	  fprintf(stderr, "thrID = %ld\n", (long) thrID);
	  parIO->thrID = thrID;
	}
      else
	pthread_attr_destroy(&parIO->attr);
    }
#endif
}

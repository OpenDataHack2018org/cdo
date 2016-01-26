#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

#if defined(HAVE_LIBCMOR)
#include "cmor.h"
#endif

void *CMOR(void *argument)
{
  cdoInitialize(argument);

#if defined(HAVE_LIBCMOR)
  int nparams = operatorArgc();
  char **params = operatorArgv();
  char tab[CMOR_MAX_STRING];
  char *option, *value;
  char *var, *chunk, *expinfo, *modinfo, *userinfo;
  if ( cdoVerbose )
    for ( int i = 0; i < nparams; ++i )
      printf("param %d: %s\n", i, params[i]);

  if ( nparams < 1 ) cdoAbort("Too few arguments!");
  strncpy(tab, params[0], CMOR_MAX_STRING);
  for ( int i = 1; i < nparams; ++i )
    {
      option = strtok(params[i], "=");
      value = strtok(NULL, "=");
      if ( value == NULL )
        cdoAbort("Missing value!");
      if ( strncasecmp(option, "var", 3) == 0 )
        var = value;
      else if ( strncasecmp(option, "chunk", 5) == 0 )
        chunk = value;
      else if ( strncasecmp(option, "expinfo", 7) == 0 )
        expinfo = value;
      else if ( strncasecmp(option, "modinfo", 7) == 0 )
        modinfo = value;
      else if ( strncasecmp(option, "userinfo", 8) == 0 )
        userinfo = value;
      else
        cdoAbort("Unknown option!");
    }

  /*
  if ( cdoVerbose )
    cdoPrint("Using CMOR version %d.%d.%d.", CMOR_VERSION_MAJOR, CMOR_VERSION_MINOR, CMOR_VERSION_PATCH);
  */
  int streamID = streamOpenRead(cdoStreamName(0));

  int vlistID = streamInqVlist(streamID);
  int taxisID = vlistInqTaxis(vlistID);

  int gridsize = vlistGridsizeMax(vlistID);
  double *array = (double*) Malloc(gridsize*sizeof(double));

  int recID, varID, levelID, nrecs, nmiss;
  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID, tsID)) )
    {
      for ( recID = 0; recID < nrecs; recID++ )
        {
          streamInqRecord(streamID, &varID, &levelID);
          streamReadRecord(streamID, array, &nmiss);
        }

      tsID++;
    }

  streamClose(streamID);

  if ( array ) Free(array);
#else
  cdoWarning("CMOR support not compiled in!");
#endif

  cdoFinish();

  return 0;
}
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */

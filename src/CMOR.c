#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "par_io.h"
#include "pstream.h"

#if defined(HAVE_LIBCMOR)
//#include "cmor.h"
#endif

void *CMOR(void *argument)
{
  cdoInitialize(argument);

#if defined(HAVE_LIBCMOR)
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

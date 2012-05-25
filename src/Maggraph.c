#if  defined  (HAVE_CONFIG_H)
#  include "config.h" /* HAVE_LIBMAGICS */
#endif

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "pstream.h"

#if  defined  (HAVE_LIBMAGICS)
#include "magics_api.h"
#endif


#if  defined  (HAVE_LIBXML)

#include<libxml/parser.h>
#include<libxml/tree.h>
#include "template_parser.h"
#include "magics_template_parser.h"
#include "results_template_parser.h"

//xmlDoc *param_doc = NULL;
//xmlNode *root_node = NULL, *magics_node = NULL, *results_node = NULL;

#endif

static
void maggraph(const char *plotfile, const char *varname, long nfiles, long nts, int *vdate, int *vtime, double **datatab)
{
  long tsID, fileID;

  for ( tsID = 0; tsID < nts; ++tsID )
    {
      printf("%6d %6d", vdate[tsID], vtime[tsID]);
      for ( fileID = 0; fileID < nfiles; ++fileID )
	printf(" %6g", datatab[fileID][tsID]);
      printf("\n");
    }

#if  defined  (HAVE_LIBMAGICS)
#endif

}

#if  defined  (HAVE_LIBMAGICS)

static
void init_MAGICS( )

{
	mag_open();
}

static
void quit_MAGICS( )

{

  mag_close ();
  fprintf( stdout,"Exiting From MAGICS\n" );

}

#endif

#define NINC_ALLOC 1024

void *Maggraph(void *argument)
{
  int operatorID;
  int varID, levelID, recID;
  int gridID;
  int nrecs;
  int tsID;
  int streamID;
  int vlistID, vlistID0 = -1;
  int nmiss;
  int zaxisID, taxisID;
  int *vdate = NULL, *vtime = NULL;
  int fileID, nfiles;
  int nts = 0, nts_alloc = 0;
  char varname[CDI_MAX_NAME];
  double missval;
  double **datatab = NULL;
  double val;
  const char *ofilename;
  char  *Filename = "combined.xml";

  cdoInitialize(argument);

  // operatorID = cdoOperatorID();

  nfiles = cdoStreamCnt() - 1;
  ofilename = cdoStreamName(nfiles);

  datatab = (double **) malloc(nfiles*sizeof(double *));

  for ( fileID = 0; fileID < nfiles; fileID++ )
    datatab[fileID] = NULL;

  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      streamID = streamOpenRead(cdoStreamName(fileID));

      vlistID = streamInqVlist(streamID);
      taxisID = vlistInqTaxis(vlistID);

      if ( fileID == 0 )
	{
	  vlistInqVarName(vlistID, 0, varname);
	  gridID = vlistInqVarGrid(vlistID, 0);

	  if ( gridInqSize(gridID) != 1 ) cdoAbort("Variable has more than one grid point!");

	  vlistID0 = vlistDuplicate(vlistID);
	}
      else
	{
	  vlistCompare(vlistID0, vlistID, CMP_ALL);
	}

      tsID = 0;
      while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	{
	  if ( nrecs != 1 ) cdoAbort("Input streams have more than one record!\n");
	  if ( fileID == 0 )
	    {
	      nts++;

	      if ( nts > nts_alloc )
		{
		  nts_alloc += NINC_ALLOC;
		  datatab[fileID] = (double *) realloc(datatab[fileID], nts_alloc*sizeof(double));
		  vdate = (int *) realloc(vdate, nts_alloc*sizeof(int));
		  vtime = (int *) realloc(vtime, nts_alloc*sizeof(int));
		}
	      vdate[tsID] = taxisInqVdate(taxisID);
	      vtime[tsID] = taxisInqVtime(taxisID);
	    }
	  else
	    {
	      if ( (tsID+1) > nts ) cdoAbort("Too many timesteps in stream %s", cdoStreamName(fileID));

	      if ( tsID == 0 )
		{
		  datatab[fileID] = (double *) malloc(nts*sizeof(double));
		}
	    }
	  
	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      streamInqRecord(streamID, &varID, &levelID);
	      streamReadRecord(streamID, &val, &nmiss);	
	      datatab[fileID][tsID] = val;
	    }

	  tsID++;
	}

      streamClose(streamID);
    }
  
#if  defined  (HAVE_LIBXML)
  /* HARDCODED THE FILE NAME .. TO BE SENT AS COMMAND LINE ARGUMENT FOR THE MAGICS OPERATOR */
  init_XMLtemplate_parser( Filename );
  updatemagics_and_results_nodes( );
#endif


#if  defined  (HAVE_LIBMAGICS)
  init_MAGICS( );
#endif

  cdoPrint(" Creating PLOT for %s", varname);
  maggraph(ofilename, varname, nfiles, nts, vdate, vtime, datatab);

#if  defined  (HAVE_LIBXML)
  quit_XMLtemplate_parser( );
#endif

#if  defined  (HAVE_LIBMAGICS)
  quit_MAGICS( );
#endif

  if ( vlistID0 != -1 ) vlistDestroy(vlistID0);

  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      if ( datatab[fileID] ) free(datatab[fileID]);
    }

  free(datatab);

  if ( vdate ) free(vdate);
  if ( vtime ) free(vtime);

  cdoFinish();

  return (0);
}

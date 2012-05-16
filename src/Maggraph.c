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
void maggraph( const char *plotfile, const char *varname, long nlon, long nlat, double *grid_center_lon, double *grid_center_lat, double *array )
{

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


void *Maggraph(void *argument)
{
  int operatorID;
  int varID, recID;
  int gridsize;
  int gridID;
  int nrecs;
  int levelID;
  int tsID;
  int streamID;
  int vlistID;
  int nmiss;
  int nlon, nlat;
  int nlev;
  int zaxisID, taxisID;
  int vdate, vtime;
  char varname[CDI_MAX_NAME];
  double missval;
  double *array = NULL;
  double *grid_center_lat = NULL, *grid_center_lon = NULL;
  char units[CDI_MAX_NAME];
  char vdatestr[32], vtimestr[32];

  char  *Filename = "combined.xml";

  cdoInitialize(argument);

  // operatorID = cdoOperatorID();

  streamID = streamOpenRead(cdoStreamName(0));

  vlistID = streamInqVlist(streamID);
  taxisID = vlistInqTaxis(vlistID);

  varID = 0;
  gridID  = vlistInqVarGrid(vlistID, varID);
  zaxisID = vlistInqVarZaxis(vlistID, varID);
  missval = vlistInqVarMissval(vlistID, varID);

  if ( gridInqType(gridID) == GRID_GME          ) cdoAbort("GME grid unspported!");
  if ( gridInqType(gridID) == GRID_UNSTRUCTURED ) cdoAbort("Unstructured grid unspported!");

  if ( gridInqType(gridID) != GRID_CURVILINEAR )
    gridID = gridToCurvilinear(gridID, 1);

  gridsize = gridInqSize(gridID);
  nlon     = gridInqXsize(gridID);
  nlat     = gridInqYsize(gridID);
  nlev     = zaxisInqSize(zaxisID);

  array           = (double *) malloc(gridsize*sizeof(double));
  grid_center_lat = (double *) malloc(gridsize*sizeof(double));
  grid_center_lon = (double *) malloc(gridsize*sizeof(double));

  gridInqYvals(gridID, grid_center_lat);
  gridInqXvals(gridID, grid_center_lon);

  /* Convert lat/lon units if required */
  gridInqXunits(gridID, units);
  gridToDegree(units, "grid center lon", gridsize, grid_center_lon);
  gridInqYunits(gridID, units);
  gridToDegree(units, "grid center lat", gridsize, grid_center_lat);
					
  tsID = 0;

#if  defined  (HAVE_LIBXML)
  /* HARDCODED THE FILE NAME .. TO BE SENT AS COMMAND LINE ARGUMENT FOR THE MAGICS OPERATOR */
  init_XMLtemplate_parser( Filename );
  updatemagics_and_results_nodes( );
#endif


#if  defined  (HAVE_LIBMAGICS)
  init_MAGICS( );
#endif

  while ( (nrecs = streamInqTimestep(streamID, tsID)) )
    {
      vdate = taxisInqVdate(taxisID);
      vtime = taxisInqVtime(taxisID);
	      
      date2str(vdate, vdatestr, sizeof(vdatestr));
      time2str(vtime, vtimestr, sizeof(vtimestr));

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID, &varID, &levelID);
	  streamReadRecord(streamID, array, &nmiss);

	  vlistInqVarName(vlistID, varID, varname);

          fprintf( stderr," Creating PLOT for %s\n",varname );

	  maggraph(cdoStreamName(1), varname, nlon, nlat, grid_center_lon, grid_center_lat, array);

	  //	  break;
	}

      break;

      tsID++;
    }

  streamClose(streamID);

  if ( array  ) free(array);
  if ( grid_center_lon ) free(grid_center_lon);
  if ( grid_center_lat ) free(grid_center_lat);

#if  defined  (HAVE_LIBXML)
  quit_XMLtemplate_parser( );
#endif

#if  defined  (HAVE_LIBMAGICS)
  quit_MAGICS( );
#endif

  cdoFinish();


  return (0);
}


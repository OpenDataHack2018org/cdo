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

static
void magplot(const char *plotfile, long nlon, long nlat, double *grid_center_lon, double *grid_center_lat, double *array)
{
  long i;

#if  defined  (HAVE_LIBMAGICS)

  // open magics
  mag_open ();

  // set the output device 
  mag_setc ("output_format",    "pdf");
  mag_setc ("output_name",      plotfile);

  // Set the input data arrays to magics++
   
  mag_set2r("input_field", array, nlon, nlat);

  // mag_set2r("input_field_latitudes", grid_center_lat, nlon, nlat);
  // mag_set2r("input_field_longitudes", grid_center_lon, nlon, nlat);
    
  mag_setr("input_field_initial_latitude", -89.75);
  mag_setr("input_field_latitude_step", 0.5);

  mag_setr("input_field_initial_longitude", -179.75);
  mag_setr("input_field_longitude_step", 0.5);


  /* Area specification (SOUTH, WEST, NORTH, EAST ) */
  mag_setr ("SUBPAGE_LOWER_LEFT_LATITUDE",   -90.0);
  mag_setr ("SUBPAGE_LOWER_LEFT_LONGITUDE", -180.0);
  mag_setr ("SUBPAGE_UPPER_RIGHT_LATITUDE",   90.0);
  mag_setr ("SUBPAGE_UPPER_RIGHT_LONGITUDE", 180.0);


  /* set up the coastline attributes */
  mag_setc ("map_coastline_colour", "khaki");
  mag_setc ("map_grid_colour",      "grey");     

  /* define the contouring parameters */
  mag_setc ("contour",                  "on");
  mag_setc ("contour_line_colour",      "sky");
  mag_setc ("CONTOUR_HIGHLIGHT_COLOUR", "GREEN");
  mag_setc ("contour_label",            "on");
  mag_cont ();

  /* plot the title text and the coastlines */
  mag_text  ();
  mag_coast ();

  mag_close ();

#else
  cdoAbort("MAGICS support not compiled in!");
#endif

}


void *Magplot(void *argument)
{
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

  cdoInitialize(argument);

  streamID = streamOpenRead(cdoStreamName(0));

  vlistID = streamInqVlist(streamID);
  taxisID = vlistInqTaxis(vlistID);

  varID = 0;
  vlistInqVarName(vlistID, varID, varname);
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

	  magplot(cdoStreamName(1), nlon, nlat, grid_center_lon, grid_center_lat, array);

	  break;
	}

      break;

      tsID++;
    }

  streamClose(streamID);

  if ( array  ) free(array);
  if ( grid_center_lon ) free(grid_center_lon);
  if ( grid_center_lat ) free(grid_center_lat);

  cdoFinish();

  return (0);
}

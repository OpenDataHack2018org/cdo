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

extern xmlNode  *magics_node;

#endif

#define DBG 0

char *line_colours[] = {

     "RGB(1.,0.,0.)",
     "RGB(0.,1.,0.)",
     "RGB(0.,0.,1.)",
     "RGB(0.,1.,1.)",
     "RGB(0.,0.,0.)",
     "RGB(1.,1.,0.)",
};

int num_colours = sizeof( line_colours )/sizeof( char* );

static
void maggraph(const char *plotfile, const char *varname, long nfiles, long nts, int *vdate, int *vtime, double **datatab)
{
  long tsID, fileID, i;
  char   *lines = "My Graph";
  double *date_time;
  char *date_time_str[nts];
  char vdatestr[32], vtimestr[32], legend_text_data[256];
  double min_val = 1.0e+200, max_val = -1.0e+200;

  date_time = (double *) malloc(nts*sizeof(double));

  if( DBG )
  {
  	printf(" %6d %6d\n", nfiles, nts );
	printf("\n");
  }

  for ( tsID = 0; tsID < nts; ++tsID )
    {
      date_time[tsID] = tsID+1;
      date2str(vdate[tsID], vdatestr, sizeof(vdatestr));
      time2str(vtime[tsID], vtimestr, sizeof(vtimestr));
      date_time_str[tsID] = (char *)malloc(256);
      sprintf(date_time_str[tsID], "%s %s", vdatestr, vtimestr);

      if( DBG )
      {
      	printf("%d: %s\n", tsID, date_time_str[tsID]);
      	printf("%6d %6d", vdate[tsID], vtime[tsID]);
      }
      for ( fileID = 0; fileID < nfiles; ++fileID )
      {
        if( DBG )
	  printf("%d\n", fileID );
	if( datatab[fileID][tsID] < min_val )
            min_val = datatab[ fileID ][ tsID ];	
	if( datatab[fileID][tsID] > max_val )
            max_val = datatab[ fileID ][ tsID ];	

        if( DBG )
        {
	   printf(" %6g", datatab[fileID][tsID]);
           printf("\n");
        }
      }
    }

    if( DBG )
    {
      printf(" %6g %6g\n", min_val, max_val );
      printf(" %s %s\n", date_time_str[0], date_time_str[ nts-1 ] );
      printf("\n");
    }

  /* 
	1. Loop over the Files
	2. Loop over the number of time steps 
	3. Set the attributes for the magics data and plot
  */  
   
#if  defined  (HAVE_LIBMAGICS)

  magics_template_parser( magics_node );
  mag_setc("output_name", plotfile);
  mag_setc("subpage_map_projection", "cartesian"); 
  mag_setr("subpage_y_length", 14.);
  mag_setr("subpage_y_position", 1.5);


  /* Horizontal Axis attributes */
  mag_setc("axis_orientation","horizontal");
  mag_setc("axis_grid", "on");
  mag_setc("axis_grid_colour", "grey");
  mag_seti("axis_grid_thickness", 1);
  mag_setc("axis_grid_line_style", "dot");
  mag_setc("axis_type", "date");
  //mag_setc("axis_date_type", "automatic");
  mag_setc("axis_date_type", "months");
  mag_setc("axis_date_min_value", date_time_str[0]);
  mag_setc("axis_date_max_value", date_time_str[nts-1]);
  mag_axis();

  /* Vertical Axis attributes */
  mag_setc("axis_orientation", "vertical");
  mag_setc("axis_grid", "on");
  mag_setc("axis_type", "regular");
  mag_setc("axis_grid_colour", "grey");
  mag_seti("axis_grid_thickness", 1);
  mag_setc("axis_grid_line_style", "dot");
  //mag_setc("graph_axis_control", "automatic");
  mag_setr("axis_min_value", min_val);
  mag_setr("axis_max_value", max_val);
  mag_axis();

  /* To automatically set the min, max for the axes, based on the input data */
  //mag_setc("graph_axis_control", "automatic");
  //mag_graph ();
  

  /* Legend */
  mag_setc("legend", "on");
  mag_setc("legend_text_colour", "black");

  for ( i = 0; i < nfiles; ++i )
    {
          sprintf(legend_text_data, "data_%d", i+1);
  	  mag_setc("graph_line_colour", line_colours[ i%num_colours ]);
	  mag_seti("graph_line_thickness", 8 );
	  mag_setc("graph_symbol", "on");
	  mag_setc("legend_user_text", legend_text_data);
	  mag_seti("graph_symbol_marker_index", 1);
	  mag_setr("graph_symbol_height", 0.5);
	  mag_set1c("graph_curve_date_x_values", date_time_str, nts);
	  mag_set1r("graph_curve_y_values", datatab[i], nts);
          mag_graph ();
    }

  mag_set1c("text_lines", &lines, 1);
  mag_setc("text_html", "true");
  mag_setc("text_colour", "black");
  mag_setr("text_font_size", 0.6);
  mag_setc("text_mode", "positional");
  mag_setr("text_box_x_position", 1.5);
  mag_setr("text_box_y_position", 16.5);
  mag_setr("text_box_x_length", 20.);
  mag_setr("text_box_y_length", 2.5);
  mag_setc("text_border", "off");
  mag_setc("text_justification", "left");
  mag_text();

  free(date_time);

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

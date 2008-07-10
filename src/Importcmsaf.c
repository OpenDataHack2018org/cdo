#if  defined  (HAVE_CONFIG_H)
#  include "config.h"
#endif

#if  defined  (HAVE_LIBHDF5)
#  include "hdf5.h"
#endif

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  NLON    1440
#define  NLAT     720
#define  MAX_VARS   6


void init_amsr_day(int vlistID, int gridID, int zaxisID, int nvars)
{
  /*
    Version-5 RSS AMSR-E or AMSR-J daily files 

    filename  with path in form satname_yyyymmdd_v5.gz
    where satname  = name of satellite (amsre or amsr)
             yyyy	= year
	       mm	= month
	       dd	= day of month

    1:time	time of measurement in fractional hours GMT
    2:sst     	sea surface temperature in deg Celcius
    3:wind	10m surface wind in meters/second
    4:vapor	columnar water vapor in millimeters
    5:cloud	cloud liquid water in millimeters
    6:rain   	rain rate in millimeters/hour
  */
  char *name[]  = {"hours", "sst", "wind", "vapor", "cloud", "rain"};
  char *units[] = {"h", "deg Celcius", "m/s", "mm", "mm", "mm/h"};
  double xscale[]  = {0.1, 0.15, 0.2, 0.3, 0.01, 0.1};
  double xminval[] = {0., -3., 0., 0., 0., 0.};
  int i, varID;

  for ( i = 0; i < nvars; ++i )
    {
      varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARIABLE);
      vlistDefVarName(vlistID, varID, name[i]);
      vlistDefVarUnits(vlistID, varID, units[i]);
      vlistDefVarDatatype(vlistID, varID, DATATYPE_INT16);
      vlistDefVarMissval(vlistID, varID, 254);
      vlistDefVarScalefactor(vlistID, varID, xscale[i]);
      vlistDefVarAddoffset(vlistID, varID, xminval[i]);
    }
}


void init_amsr_averaged(int vlistID, int gridID, int zaxisID, int nvars)
{
  /*
    Version-5 AMSR-E or AMSR-J time-averaged files including:
	  3-day		(average of 3 days ending on file date)
	  weekly	(average of 7 days ending on Saturday of file date)
	  monthly	(average of all days in month)


     filename
	   format of file names are:
	     	3-day	satname_yyyymmddv5_d3d.gz
	       	weekly	satname_yyyymmddv5.gz
	       	monthly	satname_yyyymmv5.gz

       	where	satname	=name of satellite (amsre or amsr)
		      	yyyy	=year
		       	mm	=month
		       	dd	=day of month

    1:sst       sea surface temperature in deg Celcius
    2:wind      10m surface wind in meters/second
    3:vapor	columnar water vapor in millimeters
    4:cloud     cloud liquid water in millimeters
    5:rain	rain rate in millimeters/hour
  */
  char *name[]  = {"sst", "wind", "vapor", "cloud", "rain"};
  char *units[] = {"deg Celcius", "m/s", "mm", "mm", "mm/h"};
  double xscale[]  = {0.15, 0.2, 0.3, 0.01, 0.1};
  double xminval[] = {-3., 0., 0., 0., 0.};
  int i, varID;

  for ( i = 0; i < nvars; ++i )
    {
      /* varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_CONSTANT); */
      varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARIABLE);
      vlistDefVarName(vlistID, varID, name[i]);
      vlistDefVarUnits(vlistID, varID, units[i]);
      vlistDefVarDatatype(vlistID, varID, DATATYPE_INT16);
      vlistDefVarMissval(vlistID, varID, 254);
      vlistDefVarScalefactor(vlistID, varID, xscale[i]);
      vlistDefVarAddoffset(vlistID, varID, xminval[i]);
    }
}


void read_amsr(FILE *fp, int vlistID, int nvars, double *data[], int *nmiss)
{
  static char func[] = "read_amsr_averaged";
  int varID, i, gridsize;
  unsigned char *amsr_data = NULL;
  double xminval, xscale, missval;
  size_t size;

  for ( varID = 0; varID < nvars; ++varID )
    {
      gridsize = gridInqSize(vlistInqVarGrid(vlistID, varID));
      amsr_data = (unsigned char *) realloc(amsr_data, gridsize);
      size = fread(amsr_data, 1, gridsize, fp);
      if ( (int)size != gridsize ) cdoAbort("Read error!");

      missval = vlistInqVarMissval(vlistID, varID);
      xminval = vlistInqVarAddoffset(vlistID, varID);
      xscale  = vlistInqVarScalefactor(vlistID, varID);
      
      nmiss[varID] = 0;
      for ( i = 0; i < gridsize; ++i )
	{
	  if ( amsr_data[i] <= 250 )
	    {
	      data[varID][i] = amsr_data[i]*xscale + xminval;
	    }
	  else
	    {
	      data[varID][i] = missval;
	      nmiss[varID]++;
	    }
	}
    } 

  free(amsr_data);
}


void write_data(int streamID, int nvars, double *data[], int *nmiss)
{
  int varID;

  for ( varID = 0; varID < nvars; ++varID )
    {
      streamDefRecord(streamID, varID, 0);
      streamWriteRecord(streamID, data[varID], nmiss[varID]);
    }
}


int getDate(const char *name)
{
  int date = 0;
  size_t len;
  char *pname;

  len = strlen(name);

  pname = strchr(name, '_');

  if ( pname ) date = atoi(pname+1);

  return(date);
}


#if  defined  (HAVE_LIBHDF5)
static
herr_t file_info(hid_t loc_id, const char *name, void *opdata)
{
    H5G_stat_t statbuf;

    /*
     * Get type of the object and display its name and type.
     * The name of the object is passed to this function by 
     * the Library. Some magic :-)
     */
    H5Gget_objinfo(loc_id, name, FALSE, &statbuf);
    switch (statbuf.type) {
    case H5G_GROUP: 
         printf(" Object with name %s is a group \n", name);
         break;
    case H5G_DATASET: 
         printf(" Object with name %s is a dataset \n", name);
         break;
    case H5G_TYPE: 
         printf(" Object with name %s is a named datatype \n", name);
         break;
    default:
         printf(" Unable to identify an object ");
    }
    return 0;
}
#endif


void *Importcmsaf(void *argument)
{
#if  defined  (HAVE_LIBHDF5)
  static char func[] = "Importcmsaf";
  int streamID;
  int tsID;
  int gridID, zaxisID, taxisID, vlistID;
  int gridsize;
  int i;
  int nvars;
  int vdate = 0, vtime = 0;
  double xvals[NLON], yvals[NLAT];
  double *data[MAX_VARS];
  int nmiss[MAX_VARS];
  FILE *fp;
  size_t fsize;
  hid_t	  file_id;	/* HDF5 File ID	        	*/
  herr_t  status;	/* Generic return value		*/

  cdoInitialize(argument);

  /* Open an existing file. */
  file_id = H5Fopen(cdoStreamName(0), H5F_ACC_RDONLY, H5P_DEFAULT);

  H5Giterate(file_id, "/", NULL, file_info, NULL);

  streamID = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", cdoStreamName(1));

  /*
    Longitude  is 0.25*xdim-0.125    degrees east
    Latitude   is 0.25*ydim-90.125
  */
  /*
  gridsize = NLON*NLAT;
  gridID = gridCreate(GRID_LONLAT, gridsize);
  gridDefXsize(gridID, NLON);
  gridDefYsize(gridID, NLAT);
  for ( i = 0; i < NLON; ++i ) xvals[i] = 0.25*(i+1) - 0.125;
  for ( i = 0; i < NLAT; ++i ) yvals[i] = 0.25*(i+1) - 90.125;
  gridDefXvals(gridID, xvals);
  gridDefYvals(gridID, yvals);

  zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);

  taxisID = taxisCreate(TAXIS_ABSOLUTE);

  vlistID = vlistCreate();
  vlistDefTaxis(vlistID, taxisID);

  if ( fsize == 12441600 )
    {
      nvars = 6;
      for ( i = 0; i < nvars; ++i ) data[i] = (double *) malloc(gridsize*sizeof(double));

      init_amsr_day(vlistID, gridID, zaxisID, nvars);

      streamDefVlist(streamID, vlistID);

      vtime = 130;
      for ( tsID = 0; tsID < 2; ++tsID )
	{
	  taxisDefVdate(taxisID, vdate);
	  taxisDefVtime(taxisID, vtime);
	  vtime += 1200;
	  streamDefTimestep(streamID, tsID);
	  processDefTimesteps(streamID);

	  read_amsr(fp, vlistID, nvars, data, nmiss);

	  write_data(streamID, nvars, data, nmiss);
	}

      for ( i = 0; i < nvars; ++i ) free(data[i]);
    }
  else if ( fsize == 5184000 )
    {
      nvars = 5;
      for ( i = 0; i < nvars; ++i ) data[i] = (double *) malloc(gridsize*sizeof(double));

      init_amsr_averaged(vlistID, gridID, zaxisID, nvars);

      streamDefVlist(streamID, vlistID);
      
      taxisDefVdate(taxisID, vdate);
      taxisDefVtime(taxisID, vtime);
      tsID = 0;
      streamDefTimestep(streamID, tsID);
      processDefTimesteps(streamID);
      
      read_amsr(fp, vlistID, nvars, data, nmiss);

      write_data(streamID, nvars, data, nmiss);

      for ( i = 0; i < nvars; ++i ) free(data[i]);
    }
  else
    cdoAbort("Unexpected file size for AMSR data!");

  processDefVarNum(vlistNvars(vlistID), streamID);
  */
  streamClose(streamID);

  /* Close file */
  status = H5Fclose(file_id);

  /*
  vlistDestroy(vlistID);
  gridDestroy(gridID);
  zaxisDestroy(zaxisID);
  taxisDestroy(taxisID);
  */
  cdoFinish();
#else
  cdoWarning("HDF5 support not compiled in!");
#endif

  return (0);
}

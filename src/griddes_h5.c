#if  defined  (HAVE_CONFIG_H)
#  include "config.h"
#endif

#define H5_USE_16_API

#if  defined  (HAVE_LIBHDF5)
#  include "hdf5.h"
#endif

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "griddes.h"
#include "error.h"


#if  defined  (HAVE_LIBHDF5)
static
void nce(int istat)
{
  /*
    This routine provides a simple interface to netCDF error message routine.
  */
  /*
  if ( istat != NC_NOERR ) cdoAbort(nc_strerror(istat));
  */
}
#endif


#if  defined  (HAVE_LIBHDF5)
static herr_t
obj_info(hid_t loc_id, const char *name, void *objname)
{
  H5G_obj_t obj_type;
  H5G_stat_t statbuf;
  herr_t lexist = 0;

  H5Gget_objinfo(loc_id, name, FALSE, &statbuf);

  if ( strcmp(name, (char *) objname) == 0 )
    {
      lexist = 1;

      obj_type = statbuf.type;

      switch (obj_type) {
      case H5G_GROUP:
	if ( cdoVerbose ) cdoPrint(" Object with name %s is a group", name);
	break;
      case H5G_DATASET: 
	if ( cdoVerbose ) cdoPrint(" Object with name %s is a dataset", name);
	break;
      case H5G_TYPE: 
	if ( cdoVerbose ) cdoPrint(" Object with name %s is a named datatype", name);
	break;
      default:
	/*cdoAbort(" Unable to identify an object %s", name);*/
	break;
      }
  }

  return lexist;
}
#endif


#if  defined  (HAVE_LIBHDF5)
static
int h5find_object(hid_t file_id, char *name)
{
  int lexist = 0;

  lexist = (int) H5Giterate(file_id, "/", NULL, obj_info, (void *) name);

  return lexist;
}
#endif


int gridFromH5file(const char *gridfile)
{
  static char func[] = "gridFromH5file";
  int gridID = -1;
#if  defined  (HAVE_LIBHDF5)
  hid_t	  file_id;	/* HDF5 File ID	        	*/
  hid_t	  lon_id = -1;	/* Dataset ID	        	*/
  hid_t	  lat_id = -1;	/* Dataset ID	        	*/
  hid_t   dataspace;   
  hsize_t dims_out[9];  /* dataset dimensions           */
  herr_t  status;	/* Generic return value		*/
  int     rank;
  GRID    grid;


  gridInit(&grid);

  /* Open an existing file. */
  file_id = H5Fopen(gridfile, H5F_ACC_RDONLY, H5P_DEFAULT);

  if ( h5find_object(file_id, "lon") > 0 && 
       h5find_object(file_id, "lat") > 0 )
    {
      lon_id = H5Dopen(file_id, "/lon");
      lat_id = H5Dopen(file_id, "/lat");
    }
  else if ( h5find_object(file_id, "Longitude") > 0 && 
	    h5find_object(file_id, "Latitude") > 0 )
    {
      lon_id = H5Dopen(file_id, "/Longitude");
      lat_id = H5Dopen(file_id, "/Latitude");
    }
  
  if ( lon_id >= 0 && lat_id >= 0 )
    {
      dataspace = H5Dget_space(lon_id);    /* dataspace handle */
      rank      = H5Sget_simple_extent_ndims(dataspace);
      status    = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);

      if ( rank != 2 )
	{
	  cdoWarning("Unexpected rank = %d!", rank);
	  goto RETURN;
	}
      /*
      printf("\nRank: %d\nDimensions: %lu x %lu \n", rank,
	     (unsigned long)(dims_out[1]), (unsigned long)(dims_out[0]));
      */

      grid.xsize = (int)dims_out[1];
      grid.ysize = (int)dims_out[0];
      grid.size  = grid.xsize*grid.ysize;

      grid.xvals = (double *) malloc(grid.size*sizeof(double));
      grid.yvals = (double *) malloc(grid.size*sizeof(double));

      status = H5Dread(lon_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, grid.xvals);
      status = H5Dread(lat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, grid.yvals);

      status = H5Sclose(dataspace);

      /* Close the dataset. */
      status = H5Dclose(lon_id);
      status = H5Dclose(lat_id);

      grid.type = GRID_CURVILINEAR;
      grid.prec = DATATYPE_FLT32;

      gridID = gridDefine(grid);
    }
  else  if ( h5find_object(file_id, "where") > 0 )
    {
      double xscale = 1, yscale = 1;
      double xoffset = 0, yoffset = 0;
      hid_t grp_id;
      hid_t att_id;
      int i;

      grp_id = H5Gopen(file_id, "/where/lon/what");
      if ( grp_id >= 0 )
	{
	  att_id = H5Aopen_name(grp_id, "gain");
	  if ( att_id >= 0 )
	    {
	      status = H5Aread(att_id, H5T_NATIVE_DOUBLE, &xscale);
	      H5Aclose(att_id);
	    }

	  att_id = H5Aopen_name(grp_id, "offset");
	  if ( att_id >= 0 )
	    {
	      status = H5Aread(att_id, H5T_NATIVE_DOUBLE, &xoffset);
	      H5Aclose(att_id);
	    }
	  
	  H5Gclose(grp_id);
	}

      grp_id = H5Gopen(file_id, "/where/lat/what");
      if ( grp_id >= 0 )
	{
	  att_id = H5Aopen_name(grp_id, "gain");
	  if ( att_id >= 0 )
	    {
	      status = H5Aread(att_id, H5T_NATIVE_DOUBLE, &yscale);
	      H5Aclose(att_id);
	    }

	  att_id = H5Aopen_name(grp_id, "offset");
	  if ( att_id >= 0 )
	    {
	      status = H5Aread(att_id, H5T_NATIVE_DOUBLE, &yoffset);
	      H5Aclose(att_id);
	    }
	  
	  H5Gclose(grp_id);
	}

      /* Open an existing dataset. */
      lon_id = H5Dopen(file_id, "/where/lon/data");
      if ( lon_id >= 0 )
	lat_id = H5Dopen(file_id, "/where/lat/data");

      if ( lon_id >= 0 && lat_id >= 0 )
	{
	  dataspace = H5Dget_space(lon_id);    /* dataspace handle */
	  rank      = H5Sget_simple_extent_ndims(dataspace);
	  status    = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);

	  if ( rank != 2 )
	    {
	      cdoWarning("Unexpected rank = %d!", rank);
	      goto RETURN;
	    }
	  /*
	  printf("\nRank: %d\nDimensions: %lu x %lu \n", rank,
		 (unsigned long)(dims_out[1]), (unsigned long)(dims_out[0]));
	  */

	  grid.xsize = (int)dims_out[1];
	  grid.ysize = (int)dims_out[0];
	  grid.size  = grid.xsize*grid.ysize;

	  grid.xvals = (double *) malloc(grid.size*sizeof(double));
	  grid.yvals = (double *) malloc(grid.size*sizeof(double));

	  status = H5Dread(lon_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, grid.xvals);
	  status = H5Dread(lat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, grid.yvals);

	  status = H5Sclose(dataspace);
	  
	  /* Close the dataset. */
	  status = H5Dclose(lon_id);
	  status = H5Dclose(lat_id);

	  for ( i = 0; i < grid.size; ++i ) grid.xvals[i] = grid.xvals[i]*xscale + xoffset;
	  for ( i = 0; i < grid.size; ++i ) grid.yvals[i] = grid.yvals[i]*yscale + yoffset;

	  grid.type = GRID_CURVILINEAR;
	  grid.prec = DATATYPE_FLT32;

	  gridID = gridDefine(grid);
	}
    }

  /* Close file */
  status = H5Fclose(file_id);
#else
  cdoWarning("HDF5 support not compiled in!");
#endif

 RETURN:
  return (gridID);
}

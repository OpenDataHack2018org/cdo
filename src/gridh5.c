#if  defined  (HAVE_CONFIG_H)
#  include "config.h"
#endif

#if  defined  (HAVE_LIBHDF5)
#  include "hdf5.h"
#endif

#include <stdio.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "error.h"


#if  defined  (HAVE_LIBHDF5)
static void nce(int istat)
{
  /*
    This routine provides a simple interface to netCDF error message routine.
  */
  /*
  if ( istat != NC_NOERR ) cdoAbort(nc_strerror(istat));
  */
}
#endif


/*
 * Operator function.
 */
#if  defined  (HAVE_LIBHDF5)
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


int gridFromH5file(const char *gridfile)
{
  static char func[] = "gridFromH5file";
  int gridID = -1;
#if  defined  (HAVE_LIBHDF5)
  hid_t	  file_id;	/* HDF5 File ID	        	*/
  hid_t	  lon_id;	/* Dataset ID	        	*/
  hid_t	  lat_id;	/* Dataset ID	        	*/
  hid_t   dataspace;   
  hsize_t dims_out[9];  /* dataset dimensions           */
  herr_t  status;	/* Generic return value		*/
  int     rank;
  GRID    grid;


  gridInit(&grid);

  /* Open an existing file. */
  file_id = H5Fopen(gridfile, H5F_ACC_RDWR, H5P_DEFAULT);

  /* H5Giterate(file_id, "/", NULL, file_info, NULL); */

  /* Open an existing dataset. */
  lon_id = H5Dopen(file_id, "/lon", 0);
  if ( lon_id >= 0 )
    lat_id = H5Dopen(file_id, "/lat", 0);
  
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
	     (unsigned long)(dims_out[0]), (unsigned long)(dims_out[1]));
      */

      grid.xsize = (int)dims_out[0];
      grid.ysize = (int)dims_out[1];
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
  else
    {
      double xscale = 1, yscale = 1;
      double xoffset = 0, yoffset = 0;
      hid_t grp_id;
      hid_t att_id;
      int i;

      grp_id = H5Gopen(file_id, "/where/lon/what", H5P_DEFAULT);
      if ( grp_id >= 0 )
	{
	  att_id = H5Aopen(grp_id, "gain", H5P_DEFAULT);
	  if ( att_id >= 0 )
	    {
	      status = H5Aread(att_id, H5T_NATIVE_DOUBLE, &xscale);
	      H5Aclose(att_id);
	    }

	  att_id = H5Aopen(grp_id, "offset", H5P_DEFAULT);
	  if ( att_id >= 0 )
	    {
	      status = H5Aread(att_id, H5T_NATIVE_DOUBLE, &xoffset);
	      H5Aclose(att_id);
	    }
	  
	  H5Gclose(grp_id);
	}

      grp_id = H5Gopen(file_id, "/where/lat/what", H5P_DEFAULT);
      if ( grp_id >= 0 )
	{
	  att_id = H5Aopen(grp_id, "gain", H5P_DEFAULT);
	  if ( att_id >= 0 )
	    {
	      status = H5Aread(att_id, H5T_NATIVE_DOUBLE, &yscale);
	      H5Aclose(att_id);
	    }

	  att_id = H5Aopen(grp_id, "offset", H5P_DEFAULT);
	  if ( att_id >= 0 )
	    {
	      status = H5Aread(att_id, H5T_NATIVE_DOUBLE, &yoffset);
	      H5Aclose(att_id);
	    }
	  
	  H5Gclose(grp_id);
	}

      /* Open an existing dataset. */
      lon_id = H5Dopen(file_id, "/where/lon/data", H5P_DEFAULT);
      if ( lon_id >= 0 )
	lat_id = H5Dopen(file_id, "/where/lat/data", H5P_DEFAULT);

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
		 (unsigned long)(dims_out[0]), (unsigned long)(dims_out[1]));
	  */

	  grid.xsize = (int)dims_out[0];
	  grid.ysize = (int)dims_out[1];
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

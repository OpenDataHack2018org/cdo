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
  hid_t   dapl_id;
  hid_t	  file_id;	/* HDF5 File ID	        	*/
  hid_t	  dataset_id;	/* Dataset ID	        	*/
  hid_t   dataspace;   
  hsize_t dims_out[9];  /* dataset dimensions           */
  herr_t  status;	/* Generic return value		*/
  int     rank;
  int *dset_data;
  int dset_len = 9999*9999;

  /* Open an existing file. */
  file_id = H5Fopen(gridfile, H5F_ACC_RDWR, H5P_DEFAULT);

  H5Giterate(file_id, "/", NULL, file_info, NULL);

  /* Open an existing dataset. */
  dataset_id = H5Dopen(file_id, "/lon", 0);

  if ( dataset_id < 0 )
    {
    }

  if ( dataset_id >= 0 )
    {
      dataspace = H5Dget_space(dataset_id);    /* dataspace handle */
      rank      = H5Sget_simple_extent_ndims(dataspace);
      status    = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
      printf("\nRank: %d\nDimensions: %lu x %lu \n", rank,
	     (unsigned long)(dims_out[0]), (unsigned long)(dims_out[1]));

      dset_data = (int *) malloc(dset_len*sizeof(int));
      status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);
      free(dset_data);

      status = H5Sclose(dataspace);

      /* Close the dataset. */
      status = H5Dclose(dataset_id);
    }

  /* Close file */
  status = H5Fclose(file_id);
#else
  cdoWarning("HDF5 support not compiled in!");
#endif

  return (gridID);
}

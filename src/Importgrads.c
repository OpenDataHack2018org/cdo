#if  defined  (HAVE_CONFIG_H)
#  include "config.h"
#endif

#define H5_USE_16_API

#if  defined  (HAVE_LIBHDF5)
#  include "hdf5.h"
#endif

#include <ctype.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

#include "gradsdeslib.h"


void *Importgrads(void *argument)
{
#if  defined  (HAVE_LIBHDF5)
  static char func[] = "Importgrads";
  int streamID;
  int gridID = -1, zaxisID, taxisID, vlistID;
  int i, offset;
  int nmiss;
  int ivar;
  int varID, levelID, tsID;
  int nx, ny, nz, nt, gridsize;
  double *array;
  double missval, minval, maxval;
  hid_t	  file_id;	/* HDF5 File ID	        	*/
  int  status;
  dsets_t dsets;
  int vdate, vtime;
  int *vtimes = NULL;

  cdoInitialize(argument);

  if ( cdoDefaultFileType == CDI_UNDEFID )
    cdoDefaultFileType = FILETYPE_NC;

  dsets_init(&dsets);

  status = read_gradsdes((char *)cdoStreamName(0), &dsets);
  if ( status ) return(0);


  return (0);

  if ( dsets.nsets == 0 ) cdoAbort("No dataset found!");

  gridsize = dsets.obj[0].gridsize;
  nx = dsets.obj[0].nx;
  ny = dsets.obj[0].ny;
  nz = dsets.obj[0].nz;
  nt = dsets.obj[0].nt;

  for ( ivar = 0; ivar < dsets.nsets; ++ivar )
    if ( dsets.obj[ivar].nt > 1 )
      {
	nt = dsets.obj[ivar].nt;
	break;
      }

  if ( nt > 1 )
    {
      vtimes = (int *) malloc(nt*sizeof(int));
      
      for ( i = 0; i < nt; ++i ) vtimes[i] = i*100 + 45;

      if ( dsets.obj[ivar].time )
	{
	  long itime;
	  char *pline = dsets.obj[ivar].time;

	  for ( i = 0; i < nt; ++i )
	    {
	      itime = (int)strtol(pline, &pline, 10);
	      if ( itime < 0 || itime > 2400 )
		{
		  cdoWarning("Wrong time string!");
		  break;
		}
	      vtimes[i] = (int) itime;
	    }
	}
    }

  if ( cdoVerbose )
    for ( ivar = 0; ivar < dsets.nsets; ++ivar )
      cdoPrint(" Var %d %-20s %dx%d nlev = %d nts = %d", 
	       ivar, dsets.obj[ivar].name, nx, ny, nz, dsets.obj[ivar].nt);

  for ( ivar = 1; ivar < dsets.nsets; ++ivar )
    {
      if ( nx != dsets.obj[0].nx || ny != dsets.obj[0].ny )
	cdoAbort("Gridsize must not change!");
      if ( nz != dsets.obj[0].nz )
	cdoAbort("Number of levels must not change!");
    }
  /*
  if ( dsets.lgeoloc )
    {
      gridID = read_geolocation(file_id, nx, ny, dsets.lprojtype);
    }
  else if ( dsets.lregion )
    {
      gridID = read_region(file_id, nx, ny);
    }
  */
  if ( gridID == -1 )
    {
      gridID = gridCreate(GRID_GENERIC, gridsize);
      gridDefXsize(gridID, nx);
      gridDefYsize(gridID, ny);
    }

  if ( nz == 1 )
    zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);
  else
    {
      double *levels;
      levels = (double *) malloc(nz*sizeof(double));
      for ( i = 0; i < nz; ++i ) levels[i] = i+1;
      zaxisID = zaxisCreate(ZAXIS_GENERIC, nz);
      zaxisDefLevels(zaxisID, levels);
      free(levels);
    }

  vlistID = vlistCreate();

  if ( nt > 1 )
    taxisID = taxisCreate(TAXIS_RELATIVE);
  else
    taxisID = taxisCreate(TAXIS_ABSOLUTE);

  taxisDefCalendar(taxisID, CALENDAR_STANDARD);

  vlistDefTaxis(vlistID, taxisID);

  for ( ivar = 0; ivar < dsets.nsets; ++ivar )
    {
      if ( dsets.obj[ivar].nt > 1 )
	varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARIABLE);
      else
	varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_CONSTANT);

      vlistDefVarName(vlistID, varID,  dsets.obj[ivar].name);
      if ( dsets.obj[ivar].description )
	vlistDefVarLongname(vlistID, varID,  dsets.obj[ivar].description);
      if ( dsets.obj[ivar].units )
	vlistDefVarUnits(vlistID, varID,  dsets.obj[ivar].units);
      if ( dsets.obj[ivar].title )
	vlistDefAttTxt(vlistID, varID, "title", (int)strlen(dsets.obj[ivar].title), 
		       dsets.obj[ivar].title);
	
      /*
      vlistDefVarUnits(vlistID, varID, units[i]);
      */
      vlistDefVarDatatype(vlistID, varID, dsets.obj[ivar].dtype);
      if ( dsets.obj[ivar].lmissval )
	vlistDefVarMissval(vlistID, varID, dsets.obj[ivar].missval);
      if ( dsets.obj[ivar].lscale )
	vlistDefVarScalefactor(vlistID, varID, dsets.obj[ivar].scale);
      if ( dsets.obj[ivar].loffset )
	vlistDefVarAddoffset(vlistID, varID, dsets.obj[ivar].offset);
    }
  /*
  get_global_att(file_id, "/", vlistID);
  if ( dsets.lmetadata ) get_global_att(file_id, "Metadata", vlistID);

  vdate = get_vdate(vlistID);
  if ( vdate == 0 ) vdate = 10101;
  */
  streamID = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID, vlistID);

  for ( tsID = 0; tsID < nt; ++tsID )
    {
      taxisDefVdate(taxisID, vdate);
      vtime = 0;
      if ( vtimes ) vtime = vtimes[tsID];
      taxisDefVtime(taxisID, vtime);
      streamDefTimestep(streamID, tsID);

      for ( ivar = 0; ivar < dsets.nsets; ++ivar )
	{
	  varID   = ivar;

	  if ( tsID > 0 && dsets.obj[ivar].nt == 1 ) continue;
	  
	  gridsize = dsets.obj[ivar].gridsize;
	  missval  = dsets.obj[ivar].missval;

	  for ( levelID = 0; levelID < nz; ++levelID )
	    {
	      offset = gridsize*levelID;
	      if ( nz == 1 ) offset = gridsize*tsID;
	      array  = dsets.obj[ivar].array+offset;
	      
	      nmiss  = 0;
	      minval =  1e35;
	      maxval = -1e35;

	      for ( i = 0; i < gridsize; i++ )
		{
		  if ( !DBL_IS_EQUAL(array[i], missval) )
		    {
		      if ( array[i] < minval ) minval = array[i];
		      if ( array[i] > maxval ) maxval = array[i];
		    }
		}

	      for ( i = 0; i < gridsize; i++ )
		if ( DBL_IS_EQUAL(array[i], missval) ) nmiss++;

	      if ( cdoVerbose )
		cdoPrint(" Write var %d,  level %d, nmiss %d, missval %g, minval %g, maxval %g",
			 varID, levelID, nmiss, missval, minval, maxval);
	      /*
		if ( ! (missval < minval || missval > maxval) )
		cdoWarning(" Missval is inside of valid values! Name: %s  Range: %g - %g  Missval: %g\n",
		dsets.obj[ivar].name, minval, maxval, missval);
	      */
	      streamDefRecord(streamID,  varID,  levelID);
	      streamWriteRecord(streamID, array, nmiss);
	    }
	}
    }

  /* Close file */
  status = H5Fclose(file_id);

  processDefVarNum(vlistNvars(vlistID), streamID);

  streamClose(streamID);

  vlistDestroy(vlistID);
  gridDestroy(gridID);
  zaxisDestroy(zaxisID);
  taxisDestroy(taxisID);

  for ( ivar = 0; ivar < dsets.nsets; ++ivar )
    free(dsets.obj[ivar].array);

  if ( vtimes ) free(vtimes);

  cdoFinish();
#else
  cdoWarning("HDF5 support not compiled in!");
#endif

  return (0);
}

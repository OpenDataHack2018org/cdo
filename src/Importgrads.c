#if  defined  (HAVE_CONFIG_H)
#  include "config.h"
#endif


#include <ctype.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

#include "gradsdeslib.h"


void *Importgrads(void *argument)
{
  static char func[] = "Importgrads";
  int streamID;
  int gridID = -1, zaxisID, taxisID, vlistID;
  int i, offset;
  int nmiss;
  int ivar;
  int varID, levelID, tsID;
  int nx, ny, nz, nt, gridsize;
  double missval, minval, maxval;
  int  status;
  dsets_t pfi;
  int vdate, vtime;
  int *vtimes = NULL;
  int tcur, told,fnum;
  int tmin=0,tmax=0;
  char *ch=NULL;
  off_t flen;
  int nvars;
  int recID;
  int e, flag;
  size_t rc, recsize;
  int recoffset;
  char *rec = NULL;
  struct gavar *pvar;
  struct dt dtim, dtimi;
  double fmin, fmax;
  float *farray;
  double *array;

  cdoInitialize(argument);

  dsets_init(&pfi);

  status = read_gradsdes((char *)cdoStreamName(0), &pfi);
  if ( cdoVerbose ) fprintf(stderr, "status %d\n", status);
  if ( status ) cdoAbort("read_gradsdes failed");
  /*
  printf("filename: %s\n", pfi.name);
  pfi.infile = fopen(pfi.name, "rb");
  if ( pfi.infile==NULL )  cdoAbort("Open failed on %s!", pfi.name);
  fclose (pfi.infile);
  */

  nvars = pfi.vnum;
  pvar = pfi.pvar1;
  printf("gridsize %d %d %d %d\n", pfi.gsiz, pfi.xyhdr, pfi.dnum[0], pfi.dnum[1]);
  printf("dims %d %d %d %d\n", pfi.dnum[0], pfi.dnum[1], pfi.dnum[2], pfi.dnum[3]);

  gridsize = pfi.dnum[0]*pfi.dnum[1];
  recoffset = pfi.xyhdr*4;
  if ( pfi.seqflg ) recoffset += 4;

  recsize = pfi.gsiz*4;
  rec = (char *) malloc(recsize);

  array = (double *) malloc(gridsize*sizeof(double));

  for ( ivar = 0; ivar < nvars; ++ivar )
    {
      fprintf(stderr, "%s %s %d %d %d %d %s %d \n", 
	      pvar->abbrv, pvar->longnm, pvar->offset, pvar->recoff, pvar->levels, pvar->nvardims, pvar->varnm, pvar->var_t);
      pvar++;
    }

  if (pfi.tmplat)
    for ( i = 0; i <  pfi.dnum[3]; ++i )
      printf("%d %d\n", i, pfi.fnums[i]);

  pfi.infile = NULL;
  tcur = 0;
  e = 1;
  while (1)
    {    /* loop over all times for this ensemble */
      if (pfi.tmplat)
	{
	  /* make sure no file is open */
	  if (pfi.infile!=NULL) {
	    fclose(pfi.infile);
	    pfi.infile=NULL;
	  }
	  /* advance to first valid time step for this ensemble */
	  if (tcur==0) {
	    told = 0;
	    tcur = 1;
	    while (pfi.fnums[tcur-1] == -1) tcur++;  
	  }
	  else {  /* tcur!=0 */
	    told = pfi.fnums[tcur-1];
	    /* increment time step until fnums changes */
	    while (told==pfi.fnums[tcur-1] && tcur<=pfi.dnum[3]) {
	      tcur++;
	      if ( tcur > pfi.dnum[3] ) break;
	    }
	  }

	  /* make sure we haven't advanced past end of time axis */
	  if (tcur>pfi.dnum[3]) break;

	  /* check if we're past all valid time steps for this ensemble */
	  if ((told != -1) && (pfi.fnums[tcur-1] == -1)) break;

	  /* Find the range of t indexes that have the same fnums value.
	     These are the times that are contained in this particular file */
	  tmin = tcur;
	  tmax = tcur-1;
	  fnum = pfi.fnums[tcur-1];
	  if (fnum != -1) {
	    while (fnum == pfi.fnums[tmax])
	      {
		tmax++; 
		if (tmax == pfi.dnum[3]) break;
	      }
	    gr2t(pfi.grvals[3], (gadouble)tcur, &dtim); 
	    gr2t(pfi.grvals[3], (gadouble)1, &dtimi);
	    ch = gafndt(pfi.name, &dtim, &dtimi, pfi.abvals[3], pfi.pchsub1, NULL,tcur,e,&flag);
	    if (ch==NULL) {
	      printf(" grib1map error: couldn't determine data file name for e=%d t=%d\n",e,tcur);
	      /* return(1);*/
	      exit (-1);
	    }
	  }
	}
      else { 
	/* Data set is not templated */
	ch = pfi.name;
	tmin = 1;
	tmax = pfi.dnum[3];
      }
       
      fprintf(stderr, "tmin %d tmax %d\n", tmin, tmax);
      /* Open this file and position to start of first record */
      if ( cdoVerbose) cdoPrint(" opening file: %s",ch);
      pfi.infile = fopen(ch,"rb");
      if (pfi.infile==NULL) {
	if (pfi.tmplat) {
	  if ( cdoVerbose ) cdoPrint(" Could not open file: %s",ch);
	  continue;
	} else {
	  cdoAbort(" Could not open file: %s",ch);
	}
      }
      if (pfi.tmplat) gree(ch,"312");
       
      /* Get file size */
      fseeko(pfi.infile,0L,2);
      flen = ftello(pfi.infile);

      printf("flen %d tsiz %d\n", flen, pfi.tsiz);
       
      fseeko (pfi.infile,0,0);

      for ( tsID = tmin-1; tsID < tmax; ++tsID )
	{
	  gr2t(pfi.grvals[3], (gadouble)(tsID+1), &dtim); 
	  printf("tsID = %d %d-%d-%d %d:%d\n", tsID, dtim.yr, dtim.mo, dtim.dy, dtim.hr, dtim.mn);
	  for ( recID = 0; recID < pfi.trecs; ++recID )
	    {
	      rc = fread (rec, 1, recsize, pfi.infile);
	      if ( rc < recsize )
		{
		  fprintf (stderr, "I/O error reading record!\n");
		  exit(-1);
		}
	      if ( pfi.bswap ) gabswp(rec+recoffset, gridsize);
	      farray = (float *) (rec+recoffset);
	      fmin =  1.e33;
	      fmax = -1.e33;
	      nmiss = 0;
	      for ( i = 0; i < gridsize; ++i )
		{
		  array[i] = (double) farray[i];
		  if ( array[i] > pfi.ulow && array[i] < pfi.uhi )
		    {
		      array[i] = pfi.undef;
		      nmiss++;
		    }
		  else
		    {
		      if ( array[i] < fmin ) fmin = array[i];
		      if ( array[i] > fmax ) fmax = array[i];
		    }
		}
	      printf("%d %d %g %g %d %d\n", tsID, recID, fmin, fmax, recoffset, nmiss);
	    }
	}

      /* break out if not templating */
      if (!pfi.tmplat) break;
      
    } /* end of while (1) loop */

  return (0);

  if ( pfi.nsets == 0 ) cdoAbort("No dataset found!");

  gridsize = pfi.obj[0].gridsize;
  nx = pfi.obj[0].nx;
  ny = pfi.obj[0].ny;
  nz = pfi.obj[0].nz;
  nt = pfi.obj[0].nt;

  for ( ivar = 0; ivar < pfi.nsets; ++ivar )
    if ( pfi.obj[ivar].nt > 1 )
      {
	nt = pfi.obj[ivar].nt;
	break;
      }

  if ( nt > 1 )
    {
      vtimes = (int *) malloc(nt*sizeof(int));
      
      for ( i = 0; i < nt; ++i ) vtimes[i] = i*100 + 45;

      if ( pfi.obj[ivar].time )
	{
	  long itime;
	  char *pline = pfi.obj[ivar].time;

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
    for ( ivar = 0; ivar < pfi.nsets; ++ivar )
      cdoPrint(" Var %d %-20s %dx%d nlev = %d nts = %d", 
	       ivar, pfi.obj[ivar].name, nx, ny, nz, pfi.obj[ivar].nt);

  for ( ivar = 1; ivar < pfi.nsets; ++ivar )
    {
      if ( nx != pfi.obj[0].nx || ny != pfi.obj[0].ny )
	cdoAbort("Gridsize must not change!");
      if ( nz != pfi.obj[0].nz )
	cdoAbort("Number of levels must not change!");
    }
  /*
  if ( pfi.lgeoloc )
    {
      gridID = read_geolocation(file_id, nx, ny, pfi.lprojtype);
    }
  else if ( pfi.lregion )
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

  for ( ivar = 0; ivar < pfi.nsets; ++ivar )
    {
      if ( pfi.obj[ivar].nt > 1 )
	varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARIABLE);
      else
	varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_CONSTANT);

      vlistDefVarName(vlistID, varID,  pfi.obj[ivar].name);
      if ( pfi.obj[ivar].description )
	vlistDefVarLongname(vlistID, varID,  pfi.obj[ivar].description);
      if ( pfi.obj[ivar].units )
	vlistDefVarUnits(vlistID, varID,  pfi.obj[ivar].units);
      if ( pfi.obj[ivar].title )
	vlistDefAttTxt(vlistID, varID, "title", (int)strlen(pfi.obj[ivar].title), 
		       pfi.obj[ivar].title);
	
      /*
      vlistDefVarUnits(vlistID, varID, units[i]);
      */
      vlistDefVarDatatype(vlistID, varID, pfi.obj[ivar].dtype);
      if ( pfi.obj[ivar].lmissval )
	vlistDefVarMissval(vlistID, varID, pfi.obj[ivar].missval);
      if ( pfi.obj[ivar].lscale )
	vlistDefVarScalefactor(vlistID, varID, pfi.obj[ivar].scale);
      if ( pfi.obj[ivar].loffset )
	vlistDefVarAddoffset(vlistID, varID, pfi.obj[ivar].offset);
    }
  /*
  get_global_att(file_id, "/", vlistID);
  if ( pfi.lmetadata ) get_global_att(file_id, "Metadata", vlistID);

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

      for ( ivar = 0; ivar < pfi.nsets; ++ivar )
	{
	  varID   = ivar;

	  if ( tsID > 0 && pfi.obj[ivar].nt == 1 ) continue;
	  
	  gridsize = pfi.obj[ivar].gridsize;
	  missval  = pfi.obj[ivar].missval;

	  for ( levelID = 0; levelID < nz; ++levelID )
	    {
	      offset = gridsize*levelID;
	      if ( nz == 1 ) offset = gridsize*tsID;
	      array  = pfi.obj[ivar].array+offset;
	      
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
		pfi.obj[ivar].name, minval, maxval, missval);
	      */
	      streamDefRecord(streamID,  varID,  levelID);
	      streamWriteRecord(streamID, array, nmiss);
	    }
	}
    }

  /* Close file */
  /* status = H5Fclose(file_id);*/

  processDefVarNum(vlistNvars(vlistID), streamID);

  streamClose(streamID);

  vlistDestroy(vlistID);
  gridDestroy(gridID);
  zaxisDestroy(zaxisID);
  taxisDestroy(taxisID);

  for ( ivar = 0; ivar < pfi.nsets; ++ivar )
    free(pfi.obj[ivar].array);

  if ( vtimes ) free(vtimes);

  cdoFinish();

  return (0);
}

// This file is used in CDI and CDO !!!

#if defined (HAVE_CONFIG_H)
#  include "../src/config.h"
#endif

#include <stdio.h>

#define DATE_FORMAT "%5.4d-%2.2d-%2.2d"
#define TIME_FORMAT "%2.2d:%2.2d:%2.2d"

void datetime2str(int date, int time, char *datetimestr, int maxlen)
{
  int year, month, day;
  cdiDecodeDate(date, &year, &month, &day);
  int hour, minute, second;
  cdiDecodeTime(time, &hour, &minute, &second);

  int len = sprintf(datetimestr, DATE_FORMAT "T" TIME_FORMAT, year, month, day, hour, minute, second);
  if ( len > ( maxlen-1) )
    fprintf(stderr, "Internal problem (%s): sizeof input string is too small!\n", __func__);
}


void date2str(int date, char *datestr, int maxlen)
{
  int year, month, day;
  cdiDecodeDate(date, &year, &month, &day);

  int len = sprintf(datestr, DATE_FORMAT, year, month, day);
  if ( len > ( maxlen-1) )
    fprintf(stderr, "Internal problem (%s): sizeof input string is too small!\n", __func__);
}


void time2str(int time, char *timestr, int maxlen)
{
  int hour, minute, second;
  cdiDecodeTime(time, &hour, &minute, &second);

  int len = sprintf(timestr, TIME_FORMAT, hour, minute, second);
  if ( len > ( maxlen-1) )
    fprintf(stderr, "Internal problem (%s): sizeof input string is too small!\n", __func__);
}


void printFiletype(int streamID, int vlistID)
{
  int filetype = streamInqFiletype(streamID);

  switch ( filetype )
    {
    case FILETYPE_GRB:
      printf("GRIB");
      break;
    case FILETYPE_GRB2:
      printf("GRIB2");
      break;
    case FILETYPE_NC:
      printf("NetCDF");
      break;
    case FILETYPE_NC2:
      printf("NetCDF2");
      break;
    case FILETYPE_NC4:
      printf("NetCDF4");
      break;
    case FILETYPE_NC4C:
      printf("NetCDF4 classic");
      break;
    case FILETYPE_SRV:
      printf("SERVICE");
      break;
    case FILETYPE_EXT:
      printf("EXTRA");
      break;
    case FILETYPE_IEG:
      printf("IEG");
      break;
    default:
      printf("  File format: unsupported filetype %d" , filetype);
      break;
    }

  if ( filetype == FILETYPE_SRV || filetype == FILETYPE_EXT || filetype == FILETYPE_IEG )
    {
      switch ( streamInqByteorder(streamID) )
	{
	case CDI_BIGENDIAN:
	  printf("  BIGENDIAN"); break;
	case CDI_LITTLEENDIAN:
	  printf("  LITTLEENDIAN"); break;
	default:
	  printf("  byteorder: %d undefined", streamInqByteorder(streamID)); break;
	}
    }

  if ( filetype == FILETYPE_GRB || filetype == FILETYPE_NC4 || filetype == FILETYPE_NC4C )
    {
      int nvars, varID;
      int comptype;

      nvars = vlistNvars(vlistID);

      for ( varID = 0; varID < nvars; varID++ )
	{
	  comptype = vlistInqVarCompType(vlistID, varID);
	  if ( comptype )
	    {
	      if ( comptype == COMPRESS_SZIP )
		printf(" SZIP");
	      else if ( comptype == COMPRESS_ZIP )
		printf(" ZIP");

	      break;
	    }
	}
    }

  if ( filetype == FILETYPE_GRB2 )
    {
      int comptype;
      int nvars = vlistNvars(vlistID);
      for ( int varID = 0; varID < nvars; varID++ )
	{
	  comptype = vlistInqVarCompType(vlistID, varID);
	  if ( comptype )
	    {
	      if ( comptype == COMPRESS_JPEG )
		printf(" JPEG");

	      break;
	    }
	}
    }

  printf("\n");
}

static
void printGridInfo(int vlistID)
{
  char xname[CDI_MAX_NAME], yname[CDI_MAX_NAME], xunits[CDI_MAX_NAME], yunits[CDI_MAX_NAME];
  unsigned char uuidOfHGrid[CDI_UUID_SIZE];

  int ngrids = vlistNgrids(vlistID);
  for ( int index = 0; index < ngrids; index++ )
    {
      int gridID   = vlistGrid(vlistID, index);
      int gridtype = gridInqType(gridID);
      int trunc    = gridInqTrunc(gridID);
      int gridsize = gridInqSize(gridID);
      int xsize    = gridInqXsize(gridID);
      int ysize    = gridInqYsize(gridID);
      int xysize   = xsize*ysize;
      int prec     = gridInqPrec(gridID);

      int dig = (prec == DATATYPE_FLT64) ? 15 : 7;

      gridInqXname(gridID, xname);
      gridInqYname(gridID, yname);
      gridInqXunits(gridID, xunits);
      gridInqYunits(gridID, yunits);

      fprintf(stdout, "  %4d : %-24s", index+1, gridNamePtr(gridtype));

      if ( gridtype == GRID_LONLAT   ||
	   gridtype == GRID_LCC2 ||
	   gridtype == GRID_LAEA ||
	   gridtype == GRID_SINUSOIDAL ||
	   gridtype == GRID_GENERIC ||
	   gridtype == GRID_GAUSSIAN ||
	   gridtype == GRID_GAUSSIAN_REDUCED )
	{
	  double yfirst = gridInqYval(gridID, 0);
	  double ylast  = gridInqYval(gridID, ysize-1);
	  double yinc   = gridInqYinc(gridID);

          fprintf(stdout, " : points=%d", gridsize);
	  if ( gridtype == GRID_GAUSSIAN_REDUCED )
	    fprintf(stdout, "  nlat=%d", ysize);
	  else if ( xysize )
	    fprintf(stdout, " (%dx%d)", xsize, ysize);

	  if ( gridtype == GRID_GAUSSIAN || gridtype == GRID_GAUSSIAN_REDUCED )
	    fprintf(stdout, "  np=%d", gridInqNP(gridID));

	  fprintf(stdout, "\n");

          bool lxcoord = true, lycoord = true;
          if ( gridInqXvals(gridID, NULL) == 0 ) lxcoord = false;
          if ( gridInqYvals(gridID, NULL) == 0 ) lycoord = false;

	  if ( xsize > 0 && lxcoord )
	    {
              double xfirst = gridInqXval(gridID, 0);
              double xlast  = gridInqXval(gridID, xsize-1);
              double xinc   = gridInqXinc(gridID);
              fprintf(stdout, "%33s : %.*g", xname, dig, xfirst);
              if ( xsize > 1 )
                {
                  fprintf(stdout, " to %.*g", dig, xlast);
                  if ( IS_NOT_EQUAL(xinc, 0) )
                    fprintf(stdout, " by %.*g", dig, xinc);
                }
              fprintf(stdout, " %s", xunits);
              if ( gridIsCircular(gridID) ) fprintf(stdout, "  circular");
              fprintf(stdout, "\n");
	    }

	  if ( ysize > 0 && lycoord )
	    {
	      fprintf(stdout, "%33s : %.*g", yname, dig, yfirst);
	      if ( ysize > 1 )
                {
                  fprintf(stdout, " to %.*g", dig, ylast);
                  if ( IS_NOT_EQUAL(yinc, 0) && gridtype != GRID_GAUSSIAN && gridtype != GRID_GAUSSIAN_REDUCED )
                    fprintf(stdout, " by %.*g", dig, yinc);
                }
              fprintf(stdout, " %s", yunits);
	      fprintf(stdout, "\n");
	    }

	  if ( gridIsRotated(gridID) )
	    {
	      double lonpole = gridInqXpole(gridID);
	      double latpole = gridInqYpole(gridID);
	      double angle   = gridInqAngle(gridID);
	      fprintf(stdout, "%33s : lon=%.*g  lat=%.*g", "northpole", dig, lonpole, dig, latpole);
	      if ( IS_NOT_EQUAL(angle, 0) ) fprintf(stdout, "  angle=%.*g", dig, angle);
	      fprintf(stdout, "\n");
	    }

	  if ( gridInqXbounds(gridID, NULL) || gridInqYbounds(gridID, NULL) )
	    {
	      fprintf(stdout, "%33s :", "available");
	      if ( gridInqXbounds(gridID, NULL) && gridInqYbounds(gridID, NULL) ) fprintf(stdout, " cellbounds");
	      if ( gridHasArea(gridID) )          fprintf(stdout, " area");
	      if ( gridInqMask(gridID, NULL) )    fprintf(stdout, " mask");
	      fprintf(stdout, "\n");
	    }

	  if ( gridtype == GRID_LAEA )
	    {
	      double a, lon_0, lat_0;
	      gridInqLaea(gridID, &a, &lon_0, &lat_0);
	      fprintf(stdout, "%33s : a=%g  lon_0=%g  lat_0=%g\n", "projpar", a, lon_0, lat_0);
	    }

	  if ( gridtype == GRID_LCC2 )
	    {
	      double a, lon_0, lat_0, lat_1, lat_2;
	      gridInqLcc2(gridID, &a, &lon_0, &lat_0, &lat_1, &lat_2);
	      fprintf(stdout, "%33s : a=%7.0f  lon_0=%g  lat_0=%g  lat_1=%g  lat_2=%g\n",
                      "projpar", a, lon_0, lat_0, lat_1, lat_2);
	    }
	}
      else if ( gridtype == GRID_SPECTRAL )
	{
	  fprintf(stdout, " : points=%d  nsp=%d  truncation=%d", gridsize, gridsize/2, trunc);
          if ( gridInqComplexPacking(gridID) ) fprintf(stdout, "  complexPacking");
          fprintf(stdout, "\n");
	}
      else if ( gridtype == GRID_FOURIER )
	{
	  fprintf(stdout, " : points=%d  nfc=%d  truncation=%d\n", gridsize, gridsize/2, trunc);
	}
      else if ( gridtype == GRID_GME )
	{
	  int ni = gridInqGMEni(gridID);
	  int nd = gridInqGMEnd(gridID);
	  fprintf(stdout, " : points=%d  nd=%d  ni=%d\n", gridsize, nd, ni);
	}
      else if ( gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED )
	{
	  if ( gridtype == GRID_CURVILINEAR )
	    fprintf(stdout, " : points=%d (%dx%d)", gridsize, xsize, ysize);
	  else
	    fprintf(stdout, " : points=%d", gridsize);

          if ( gridtype == GRID_UNSTRUCTURED && gridInqNvertex(gridID) > 0 )
	    fprintf(stdout, "  nvertex=%d", gridInqNvertex(gridID));

          fprintf(stdout, "\n");

          if ( gridtype == GRID_UNSTRUCTURED )
            {
              int number   = gridInqNumber(gridID);
              int position = gridInqPosition(gridID);

              if ( number > 0 )
                {
                  fprintf(stdout, "%33s : number=%d  position=%d\n", "grid", number, position);
                }

              if ( gridInqReference(gridID, NULL) )
                {
                  char reference_link[8192];
                  gridInqReference(gridID, reference_link);
                  fprintf(stdout, "%33s : %s\n", "uri", reference_link);
                }
            }

	  if ( gridInqXvals(gridID, NULL) && gridInqYvals(gridID, NULL) )
	    {
	      double *xvals = (double*) malloc((size_t)gridsize*sizeof(double));
	      double *yvals = (double*) malloc((size_t)gridsize*sizeof(double));

	      gridInqXvals(gridID, xvals);
	      gridInqYvals(gridID, yvals);

	      double xfirst = xvals[0];
	      double xlast  = xvals[0];
	      double yfirst = yvals[0];
	      double ylast  = yvals[0];
	      for ( int i = 1; i < gridsize; i++ )
		{
		  if ( xvals[i] < xfirst ) xfirst = xvals[i];
		  if ( xvals[i] > xlast  ) xlast  = xvals[i];
		  if ( yvals[i] < yfirst ) yfirst = yvals[i];
		  if ( yvals[i] > ylast  ) ylast  = yvals[i];
		}

	      fprintf(stdout, "%33s : %.*g to %.*g %s", xname, dig, xfirst, dig, xlast, xunits);
	      if ( gridIsCircular(gridID) ) fprintf(stdout, "  circular");
	      fprintf(stdout, "\n");
	      fprintf(stdout, "%33s : %.*g to %.*g %s\n", yname, dig, yfirst, dig, ylast, yunits);

	      free(xvals);
	      free(yvals);
	    }
	}
      else if ( gridtype == GRID_LCC )
	{
	  double originLon, originLat, lonParY, lat1, lat2, xincm, yincm;
	  int projflag, scanflag;

	  gridInqLCC(gridID, &originLon, &originLat, &lonParY, &lat1, &lat2, &xincm, &yincm,
		     &projflag, &scanflag);

	  fprintf(stdout, " : points=%d (%dx%d)  ", gridsize, xsize, ysize);
	  if ( (projflag&128) == 0 )
	    fprintf(stdout, "North Pole\n");
	  else
	    fprintf(stdout, "South Pole\n");

	  fprintf(stdout, "%33s : originLon=%g  originLat=%g  lonParY=%g\n", " ", originLon, originLat, lonParY);
	  fprintf(stdout, "%33s : lat1=%g  lat2=%g  xinc=%g m  yinc=%g m\n", " ", lat1, lat2, xincm, yincm);
	}
      else /* if ( gridtype == GRID_GENERIC ) */
	{
	  if ( ysize == 0 )
	    fprintf(stdout, " : points=%d\n", gridsize);
	  else
            fprintf(stdout, " : points=%d (%dx%d)\n", gridsize, xsize, ysize);
	}

      if ( gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED || gridtype == GRID_LCC )
	{
	  if ( gridHasArea(gridID) ||
	       gridInqXbounds(gridID, NULL) || gridInqYbounds(gridID, NULL) )
	    {
	      fprintf(stdout, "%33s :", "available");
	      if ( gridInqXbounds(gridID, NULL) && gridInqYbounds(gridID, NULL) ) fprintf(stdout, " cellbounds");
	      if ( gridHasArea(gridID) )          fprintf(stdout, " area");
	      if ( gridInqMask(gridID, NULL) )    fprintf(stdout, " mask");
	      fprintf(stdout, "\n");
	    }
	}

      gridInqUUID(gridID, uuidOfHGrid);
      if ( !cdiUUIDIsNull(uuidOfHGrid) )
        {
          char uuidOfHGridStr[37];
          cdiUUID2Str(uuidOfHGrid, uuidOfHGridStr);
          if ( uuidOfHGridStr[0] != 0  && strlen(uuidOfHGridStr) == 36 )
            {
	      fprintf(stdout, "%33s : %s\n", "uuid", uuidOfHGridStr);
            }
        }
    }
}

static
void printZaxisInfo(int vlistID)
{
  char zaxisname[CDI_MAX_NAME], zname[CDI_MAX_NAME], zunits[CDI_MAX_NAME];

  int nzaxis = vlistNzaxis(vlistID);
  for ( int index = 0; index < nzaxis; index++ )
    {
      double zinc = 0;
      int zaxisID   = vlistZaxis(vlistID, index);
      int zaxistype = zaxisInqType(zaxisID);
      int ltype     = zaxisInqLtype(zaxisID);
      int levelsize = zaxisInqSize(zaxisID);
      int prec      = zaxisInqPrec(zaxisID);

      int dig = (prec == DATATYPE_FLT64) ? 15 : 7;

      zaxisName(zaxistype, zaxisname);
      zaxisInqName(zaxisID, zname);
      zaxisInqUnits(zaxisID, zunits);
      zunits[12] = 0;

      if ( zaxistype == ZAXIS_GENERIC && ltype != 0 )
        fprintf(stdout, "  %4d : %-12s (ltype=%3d) :", vlistZaxisIndex(vlistID, zaxisID)+1, zaxisname, ltype);
      else
        fprintf(stdout, "  %4d : %-24s :", vlistZaxisIndex(vlistID, zaxisID)+1, zaxisname);

      fprintf(stdout, " levels=%d", levelsize);
      fprintf(stdout, "\n");

      double *levels = (double*) malloc((size_t)levelsize*sizeof(double));
      zaxisInqLevels(zaxisID, levels);

      if ( !(zaxistype == ZAXIS_SURFACE && levelsize == 1 && !(fabs(levels[0]) > 0)) )
        {
          double zfirst = levels[0];
          double zlast  = levels[levelsize-1];
          if ( levelsize > 2 )
            {
              int levelID;
              zinc = (levels[levelsize-1] - levels[0]) / (levelsize-1);
              for ( levelID = 2; levelID < levelsize; ++levelID )
                if ( fabs(fabs(levels[levelID] - levels[levelID-1]) - zinc) > 0.001*zinc ) break;

              if ( levelID < levelsize ) zinc = 0;
            }

          fprintf(stdout, "%33s : %.*g", zname, dig, zfirst);
          if ( levelsize > 1 )
            {
              fprintf(stdout, " to %.*g", dig, zlast);
              if ( IS_NOT_EQUAL(zinc, 0) )
                fprintf(stdout, " by %.*g", dig, zinc);
            }
          fprintf(stdout, " %s", zunits);
          fprintf(stdout, "\n");
        }

      free(levels);

      if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
        {
          double level1, level2;
          fprintf(stdout, "%33s : ", "bounds");

          level1 = zaxisInqLbound(zaxisID, 0);
          level2 = zaxisInqUbound(zaxisID, 0);
          fprintf(stdout, "%.*g-%.*g", dig, level1, dig, level2);
          if ( levelsize > 1 )
            {
              level1 = zaxisInqLbound(zaxisID, levelsize-1);
              level2 = zaxisInqUbound(zaxisID, levelsize-1);
              fprintf(stdout, " to %.*g-%.*g", dig, level1, dig, level2);
              if ( IS_NOT_EQUAL(zinc, 0) )
                fprintf(stdout, " by %.*g", dig, zinc);
            }
          fprintf(stdout, " %s", zunits);
          fprintf(stdout, "\n");
        }

      if ( zaxistype == ZAXIS_HYBRID )
        {
          char psname[CDI_MAX_NAME];
          psname[0] = 0;
          zaxisInqPsName(zaxisID, psname);
          int vctsize = zaxisInqVctSize(zaxisID);
          if ( vctsize || psname[0] )
            {
	      fprintf(stdout, "%33s :", "available");
              if ( vctsize   ) fprintf(stdout, " vct");
              if ( psname[0] ) fprintf(stdout, "  ps: %s", psname);
              fprintf(stdout, "\n");
            }
        }

      if ( zaxistype == ZAXIS_REFERENCE )
        {
          int number   = zaxisInqNumber(zaxisID);

          if ( number > 0 )
            {
              fprintf(stdout, "%33s : ", "zaxis");
              fprintf(stdout, "number = %d\n", number);
            }

          unsigned char uuidOfVGrid[CDI_UUID_SIZE];
          zaxisInqUUID(zaxisID, uuidOfVGrid);
          if ( !cdiUUIDIsNull(uuidOfVGrid) )
            {
              char uuidOfVGridStr[37];
              cdiUUID2Str(uuidOfVGrid, uuidOfVGridStr);
              if ( uuidOfVGridStr[0] != 0  && strlen(uuidOfVGridStr) == 36 )
                {
                  fprintf(stdout, "%33s : ", "uuid");
                  fprintf(stdout, "%s\n", uuidOfVGridStr);
                }
            }
        }
    }
}

static
void printSubtypeInfo(int vlistID)
{
  int nsubtypes = vlistNsubtypes(vlistID);
  for ( int index = 0; index < nsubtypes; index++)
    {
      int subtypeID = vlistSubtype(vlistID, index);
      int subtypesize = subtypeInqSize(subtypeID);
      // subtypePrint(subtypeID);
      fprintf(stdout, "  %4d : %-24s :", vlistSubtypeIndex(vlistID, subtypeID)+1, "tiles");
      fprintf(stdout, " ntiles=%d", subtypesize);
      fprintf(stdout, "\n");
    }
}

static
int printDateTime(int ntimeout, int vdate, int vtime)
{
  char vdatestr[32], vtimestr[32];

  if ( ntimeout == 4 )
    {
      ntimeout = 0;
      fprintf(stdout, "\n");
    }

  date2str(vdate, vdatestr, sizeof(vdatestr));
  time2str(vtime, vtimestr, sizeof(vtimestr));

  fprintf(stdout, " %s %s", vdatestr, vtimestr);

  return (++ntimeout);
}

#define NUM_TIMESTEP 60
#define MAX_DOTS     80

static
int printDot(int ndotout, int *nfact, int *ncout)
{
  //printf("ncout %d %d %d\n",*ncout, (*ncout)%(*nfact), *nfact);
  if ( (*ncout)%(*nfact) == 0 )
    {
      if ( ndotout == MAX_DOTS )
	{
	  *ncout = 0;
	  ndotout = 0;
	  fprintf(stdout, "\n   ");
	  (*nfact) *= 10;
	}

      fprintf(stdout, ".");
      fflush(stdout);
      ndotout++;
    }

  (*ncout)++;

  return ndotout;
}

static
void printTimesteps(int streamID, int taxisID, int verbose)
{
  int nrecs;
  struct datetime {
    int vdate;
    int vtime;
    struct datetime *next;
  };
  struct datetime vdatetime[NUM_TIMESTEP];
  struct datetime *next_vdatetime = vdatetime;

  for ( int i = 0; i < NUM_TIMESTEP-1; ++i ) vdatetime[i].next = &vdatetime[i+1];
  vdatetime[NUM_TIMESTEP-1].next = &vdatetime[0];

  int ntimeout = 0;
  int ndotout = 0;
  int nvdatetime = 0;
  int ncout = 0;
  int nfact = 1;
  int tsID = 0;

#ifdef CDO
  dtlist_type *dtlist = dtlist_new();
#endif
  while ( (nrecs = streamInqTimestep(streamID, tsID)) )
    {
#ifdef CDO
      dtlist_taxisInqTimestep(dtlist, taxisID, 0);
      int vdate = dtlist_get_vdate(dtlist, 0);
      int vtime = dtlist_get_vtime(dtlist, 0);
#else
      int vdate = taxisInqVdate(taxisID);
      int vtime = taxisInqVtime(taxisID);
#endif

      if ( verbose || tsID < NUM_TIMESTEP )
	{
	  ntimeout = printDateTime(ntimeout, vdate, vtime);
	}
      else
	{
	  if ( tsID == 2*NUM_TIMESTEP ) fprintf(stdout, "\n   ");
	  if ( tsID >= 2*NUM_TIMESTEP ) ndotout = printDot(ndotout, &nfact, &ncout);

	  if ( nvdatetime < NUM_TIMESTEP )
	    {
	      vdatetime[nvdatetime].vdate = vdate;
	      vdatetime[nvdatetime].vtime = vtime;
	      nvdatetime++;
	    }
	  else
	    {
	      next_vdatetime->vdate = vdate;
	      next_vdatetime->vtime = vtime;
	      next_vdatetime = next_vdatetime->next;
	    }
	}

      tsID++;
    }

#ifdef CDO
  dtlist_delete(dtlist);
#endif
  if ( nvdatetime )
    {
      fprintf(stdout, "\n");

      ntimeout = 0;
      int toff = 0;
      if ( tsID > 2*NUM_TIMESTEP )
        {
          toff = tsID%4;
          if ( toff > 0 ) toff = 4 - toff;
          for ( int i = 0; i < toff; ++i ) next_vdatetime = next_vdatetime->next;
        }
      for ( int i = toff; i < nvdatetime; ++i )
	{
	  int vdate = next_vdatetime->vdate;
	  int vtime = next_vdatetime->vtime;
	  ntimeout = printDateTime(ntimeout, vdate, vtime);
	  next_vdatetime = next_vdatetime->next;
	}
    }
}
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */

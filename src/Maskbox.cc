/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2018 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

/*
   This module contains the following operators:

      Maskbox    masklonlatbox   Mask lon/lat box
      Maskbox    maskindexbox    Mask index box
      Maskbox    maskregion      Mask regions
*/


#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "grid.h"

#define MAX_LINE 256
#define MAX_VALS 1048576

static
int ReadCoords(double *xvals, double *yvals, const char *polyfile, FILE *fp)
{
  int z = 0, number = 0, jumpedlines = 0;
  char line[MAX_LINE];

  while ( readline(fp, line, MAX_LINE) )
    {
      int i = 0;
      if ( line[0] == '#' ) 
         { 
           jumpedlines++;
           continue;
         }  
      if ( line[0] == '\0' )
         { 
           jumpedlines++;
           continue;
         }  
      if ( line[0] == '&' ) break;
	 
      char *linep = &line[0];
     
      double xcoord = strtod(linep, &linep);
      
      if ( ! (fabs(xcoord) > 0) ) jumpedlines++;
      
      while( ( ( isdigit( (int) *linep ) == FALSE ) && ( *linep!='-' )) && ( i < 64) )
	{	  
	  if ( *linep == 0 ) 
	    {
	      cdoAbort(" y value missing in file %s at  line %d", polyfile, (number+jumpedlines+1));
	      break;
	    }
          if ( ( isspace( (int) *linep) == FALSE && ( *linep!='-' ) ) && ( linep != NULL ) ) 
	    cdoWarning("unknown character in file %s at line %d", polyfile, (number+jumpedlines+1) );

          linep++;
	  i++;
	}    
      if ( ( i >= 63 ) && ( number != 0 ) ) cdoAbort( "Wrong value format in file %s at line %d", polyfile, (number+jumpedlines+1) );
    
      double ycoord = strtod(linep, NULL);
     
      xvals[number] = xcoord;
      yvals[number] = ycoord;

      number++;
    }

 
  if ( ( number != 0 )&& ( ! (IS_EQUAL(xvals[0], xvals[number-1]) && IS_EQUAL(yvals[0], yvals[number-1])) ) )
    {
      xvals[number] = xvals[0];
      yvals[number] = yvals[0];
      number++;
    }
  

  if ( cdoVerbose ) 
    {
      for ( z = 0; z < number; z++ )  fprintf(stderr, "%d %g %g\n",  (z+1),  xvals[z], yvals[z]);
    }

  return number;
}


void genlonlatbox(int argc_offset, int gridID1, long *lat1, long *lat2, long *lon11, long *lon12, long *lon21, long *lon22);

void genindexbox(int argc_offset, int gridID1, long *lat1, long *lat2, long *lon11, long *lon12, long *lon21, long *lon22);


static
void maskbox(std::vector<bool> &mask, int gridID, size_t lat1, size_t lat2, size_t lon11, size_t lon12, size_t lon21, size_t lon22)
{
  size_t nlon = gridInqXsize(gridID);
  size_t nlat = gridInqYsize(gridID);

  for ( size_t ilat = 0; ilat < nlat; ilat++ )
    for ( size_t ilon = 0; ilon < nlon; ilon++ )
      if (  (lat1 <= ilat && ilat <= lat2 && 
	      ((lon11 <= ilon && ilon <= lon12) || (lon21 <= ilon && ilon <= lon22))) )
	mask[nlon*ilat + ilon] = false;
}

void getlonlatparams(int argc_offset, double *xlon1, double *xlon2, double *xlat1, double *xlat2);

static
void maskbox_cell(std::vector<bool> &mask, int gridID)
{
  double xlon1 = 0, xlon2 = 0, xlat1 = 0, xlat2 = 0;
  getlonlatparams(0, &xlon1, &xlon2, &xlat1, &xlat2);

  size_t gridsize = gridInqSize(gridID);

  std::vector<double> xvals(gridsize);
  std::vector<double> yvals(gridsize);

  gridInqXvals(gridID, &xvals[0]);
  gridInqYvals(gridID, &yvals[0]);

  char xunits[CDI_MAX_NAME];
  char yunits[CDI_MAX_NAME];
  gridInqXunits(gridID, xunits);
  gridInqYunits(gridID, yunits);

  double xfact = 1, yfact = 1;
  if ( strncmp(xunits, "radian", 6) == 0 ) xfact = RAD2DEG;
  if ( strncmp(yunits, "radian", 6) == 0 ) yfact = RAD2DEG;

  if ( xlon1 > xlon2 ) 
    cdoAbort("The second longitude have to be greater than the first one!");

  if ( xlat1 > xlat2 )
    {
      double xtemp = xlat1;
      xlat1 = xlat2;
      xlat2 = xtemp;
    }
  
  for ( size_t i = 0; i < gridsize; ++i )
    {
      mask[i] = true;

      double xval = xfact*xvals[i];
      double yval = yfact*yvals[i];
      if ( yval >= xlat1 && yval <= xlat2 )
        {
          if ( ((xval     >= xlon1 && xval     <= xlon2) ||
                (xval-360 >= xlon1 && xval-360 <= xlon2) ||
                (xval+360 >= xlon1 && xval+360 <= xlon2)) )
            {
              mask[i] = false;
            }
        }
    }
}

static
void maskregion(std::vector<bool> &mask, int gridID, double *xcoords, double *ycoords, size_t nofcoords)
{
  size_t nlon = gridInqXsize(gridID);
  size_t nlat = gridInqYsize(gridID);

  std::vector<double> xvals(nlon);
  std::vector<double> yvals(nlat);

  gridInqXvals(gridID, &xvals[0]);
  gridInqYvals(gridID, &yvals[0]);  

  /* Convert lat/lon units if required */
  {
    char units[CDI_MAX_NAME];
    gridInqXunits(gridID, units);
    grid_to_degree(units, nlon, &xvals[0], "grid center lon");
    gridInqYunits(gridID, units);
    grid_to_degree(units, nlat, &yvals[0], "grid center lat");
  }

  double xmin, xmax;
  arrayMinMax(nlon, &xvals[0], &xmin, &xmax);

  double ymin, ymax;
  arrayMinMax(nlat, &yvals[0], &ymin, &ymax);

  for ( size_t ilat = 0; ilat < nlat; ilat++ )
    {
      double yval = yvals[ilat];
      for ( size_t ilon = 0; ilon < nlon; ilon++ )
	{
          size_t i, j;
          int c = 0;
	  double xval = xvals[ilon];
	  if (!( ( ( xval > xmin ) || ( xval < xmax ) ) || ( (yval > ymin) || (yval < ymax) ) ) ) c = !c;
	  	  
          if ( c == 0 )
	    {
	      for ( i = 0, j = nofcoords-1; i < nofcoords; j = i++ )
	    
	      if ((((ycoords[i]<=yval) && (yval<ycoords[j])) ||
		   ((ycoords[j]<=yval) && (yval<ycoords[i]))) &&
		  ((xval) < (xcoords[j] - (xcoords[i])) * (yval - ycoords[i]) / (ycoords[j] - ycoords[i]) +(xcoords[i])))
		c = !c;
	    }

	  if ( c == 0 )
	    {
	      for ( i = 0, j = nofcoords-1; i < nofcoords; j = i++ )
		{
		  if ( xvals[ilon] > 180 )
		    {
		      if ((((ycoords[i]<=yval) && (yval<ycoords[j])) ||
			   ((ycoords[j]<=yval) && (yval<ycoords[i]))) &&
			  ((xval-360) < (xcoords[j] - (xcoords[i])) * (yval - ycoords[i]) / (ycoords[j] - ycoords[i]) +(xcoords[i])))
			c = !c;
		    }
		}
	    }
 
	  if ( c == 0 )
	    {
	      for ( i = 0, j = nofcoords-1; i< nofcoords; j = i++ )
		{		
		  if ( xval < 0 )
		    {
		      if ((((ycoords[i]<=yval) && (yval<ycoords[j])) ||
			   ((ycoords[j]<=yval) && (yval<ycoords[i]))) &&
			  ((xval+360) < (xcoords[j] - (xcoords[i])) * (yval - ycoords[i]) / (ycoords[j] - ycoords[i]) +(xcoords[i])))
			c = !c;
		    }
		}
	    }
	     
	  if ( c != 0 ) mask[nlon*ilat+ilon] =  false;
	}
    }
}


void *Maskbox(void *process)
{
  int nrecs;
  int varID, levelID;
  int gridID = -1;
  int index, gridtype = -1;
  long lat1, lat2, lon11, lon12, lon21, lon22;

  cdoInitialize(process);

  // clang-format off
  int MASKLONLATBOX = cdoOperatorAdd("masklonlatbox", 0, 0, "western and eastern longitude and southern and northern latitude");
  int MASKINDEXBOX  = cdoOperatorAdd("maskindexbox",  0, 0, "index of first and last longitude and index of first and last latitude");
  int MASKREGION    = cdoOperatorAdd("maskregion",    0, 0, "limiting coordinates of the region");
  // clang-format on

  int operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));
  int vlistID1 = cdoStreamInqVlist(streamID1);

  int ngrids = vlistNgrids(vlistID1);
  int ndiffgrids = 0;
  for ( index = 1; index < ngrids; index++ )
    if ( vlistGrid(vlistID1, 0) != vlistGrid(vlistID1, index))
      ndiffgrids++;

  for ( index = 0; index < ngrids; index++ )
    {
      gridID   = vlistGrid(vlistID1, index);
      gridtype = gridInqType(gridID);

      if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN ) break;
      if ( operatorID != MASKREGION && gridtype == GRID_CURVILINEAR ) break;
      if ( operatorID != MASKREGION && gridtype == GRID_UNSTRUCTURED ) break;
      if ( operatorID == MASKINDEXBOX && gridtype == GRID_GENERIC &&
	  gridInqXsize(gridID) > 0 && gridInqYsize(gridID) > 0 ) break;
    }

  if ( gridtype == GRID_GAUSSIAN_REDUCED )
    cdoAbort("Gaussian reduced grid found. Use option -R to convert it to a regular grid!");

  if ( index == ngrids ) cdoAbort("No regular lon/lat grid found!");
  if ( ndiffgrids > 0 )  cdoAbort("Too many different grids!");

  operatorInputArg(cdoOperatorEnter(operatorID));

  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int nvars = vlistNvars(vlistID1);
  std::vector<bool> vars(nvars);
  for ( varID = 0; varID < nvars; varID++ )
    {
      vars[varID] = (gridID == vlistInqVarGrid(vlistID1, varID));
    }

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  size_t gridsize = gridInqSize(gridID);
  std::vector<double> array(gridsize);
  std::vector<bool> mask(gridsize);
  for ( size_t i = 0; i < gridsize; ++i ) mask[i] = true;
 
  if ( operatorID == MASKLONLATBOX )
    {
      if ( gridtype == GRID_CURVILINEAR ||
           gridtype == GRID_UNSTRUCTURED )
        {
          maskbox_cell(mask, gridID);
        }
      else
        {
          genlonlatbox(0, gridID, &lat1, &lat2, &lon11, &lon12, &lon21, &lon22);
          maskbox(mask, gridID, lat1, lat2, lon11, lon12, lon21, lon22);
        }
    }
  else if ( operatorID == MASKINDEXBOX )
    {
      genindexbox(0, gridID, &lat1, &lat2, &lon11, &lon12, &lon21, &lon22);
      maskbox(mask, gridID, lat1, lat2, lon11, lon12, lon21, lon22);
    }
  else if ( operatorID == MASKREGION )
    {
      double *xcoords = (double*) Malloc(MAX_VALS*sizeof(double));
      double *ycoords = (double*) Malloc(MAX_VALS*sizeof(double));
      int nfiles = operatorArgc();
     
      for ( int i = 0; i < nfiles; i++ )
	{
	  const char *polyfile = operatorArgv()[i];
	  FILE *fp = fopen(polyfile, "r");
	 
	  if ( fp == 0 ) cdoAbort("Open failed on %s", polyfile);   
	  while ( TRUE )
	    {
	      int number = ReadCoords(xcoords, ycoords, polyfile, fp );
	      if ( number == 0 ) break;
	      if ( number < 3 ) cdoAbort( "Too less values in file %s", polyfile );
	      maskregion(mask, gridID, xcoords, ycoords, number);
	    }
	  fclose(fp); 
	}
      if ( xcoords ) Free( xcoords );
      if ( ycoords ) Free( ycoords );
    }

  int tsID = 0;
  while ( (nrecs = cdoStreamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);
	       
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);

	  if ( vars[varID] )
	    {
              size_t nmiss;
	      pstreamReadRecord(streamID1, &array[0], &nmiss);

              double missval = vlistInqVarMissval(vlistID1, varID);             
	      for ( size_t i = 0; i < gridsize; i++ )
                if ( mask[i] ) array[i] = missval;
		
	      nmiss = arrayNumMV(gridsize, &array[0], missval);
	      pstreamDefRecord(streamID2, varID, levelID);
	      pstreamWriteRecord(streamID2, &array[0], nmiss);
	    }
	}
      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

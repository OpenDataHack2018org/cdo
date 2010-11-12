/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2010 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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
   Ensstat3       ensrkhist_space  Ensemble ranked histogram averaged over time
   Ensstat3       ensrkhist_time   Ensemble ranked histogram averaged over space
   Ensstat3       ensroccurve      Ensamble Receiver Operating Characteristics
*/

#if defined (_OPENMP)
#  include <omp.h>
#endif

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"

enum TDATA_TYPE {TIME, SPACE};

#define time_data TIME
#define space_data SPACE

void *Ensstat3(void *argument)
{
  int operatorID;
  int operfunc, datafunc;
  int i,j,binID;
  int nvars;
  int cmpflag;
  int varID, recID;
  int gridsize = 0;
  int gridID, gridID2;
  int have_miss;
  int nrecs, nrecs0;
  int levelID;
  int tsID;
  int streamID = 0, streamID2;
  int vlistID, vlistID1, vlistID2;
  int nmiss;
  int taxisID1, taxisID2;
  int zaxisID,zaxisID2;
  int ompthID;
  int *varID2;
  int time_mode;
  double missval;
  double *levs;
  double val;
  int **array2 = NULL;
  field_t *field;
  int fileID, nfiles;
  const char *ofilename;

  typedef struct
  {
    int streamID;
    int vlistID;
    double *array;
  } ens_file_t;
  ens_file_t *ef = NULL;

  cdoInitialize(argument);

  cdoOperatorAdd("ensroc",          func_roc,  0,          NULL);
  cdoOperatorAdd("ensrkhist_space", func_rank, space_data, NULL);
  cdoOperatorAdd("ensrkhist_time",  func_rank, time_data,  NULL);
  
  operatorID = cdoOperatorID();
  operfunc = cdoOperatorF1(operatorID);
  datafunc = cdoOperatorF2(operatorID);
  
  nfiles = cdoStreamCnt() - 1;

  if ( cdoVerbose )
    cdoPrint("Ensemble over %d files.", nfiles);

  ofilename = cdoStreamName(nfiles);

  if ( !cdoSilentMode )
    if ( fileExist(ofilename) )
      if ( !userFileOverwrite(ofilename) )
	cdoAbort("Outputfile %s already exist!", ofilename);

  ef = (ens_file_t *) malloc(nfiles*sizeof(ens_file_t));


  /* *************************************************** */
  /* should each thread be allocating memory locally???? */
  /* ("first touch strategy")                            */
  /* --> #pragma omp parallel for ...                    */
  /* *************************************************** */
#if defined (_OPENMP)
  field = (field_t *) malloc(omp_get_max_threads()*sizeof(field_t));
  for ( i = 0; i < omp_get_max_threads(); i++ )
#else
  field = (field_t *) malloc(1*sizeof(field_t));
  for ( i = 0; i < 1; i++ )
#endif
    {
      field[i].size   = nfiles;
      field[i].ptr    = (double *) malloc(nfiles*sizeof(double));
      field[i].weight = (double *) malloc(nfiles*sizeof(double));
      for ( fileID = 0; fileID < nfiles; fileID++ )
	field[i].weight[fileID] = 1;
    }

  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      streamID = streamOpenRead(cdoStreamName(fileID));
      fprintf(stderr,"opened %s\n",cdoStreamName(fileID));
      if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", cdoStreamName(fileID));

      vlistID = streamInqVlist(streamID);

      ef[fileID].streamID = streamID;
      ef[fileID].vlistID = vlistID;
    }

  /* check that the contents is always the same */
  nvars = vlistNvars(ef[0].vlistID);
  if ( nvars == 1 ) 
    cmpflag = CMP_NAME | CMP_GRIDSIZE | CMP_NLEVEL | CMP_GRID;
  else // What is this supposed to do different? - is there missing a bracket?
    cmpflag = CMP_NAME | CMP_GRIDSIZE | CMP_NLEVEL | CMP_GRID;

  for ( fileID = 1; fileID < nfiles; fileID++ )
    vlistCompare(ef[0].vlistID, ef[fileID].vlistID, cmpflag);

  vlistID1 = ef[0].vlistID;
  vlistID2 = vlistCreate();
  varID2 = (int *) malloc( nvars*sizeof(int));

  levs = (double*) calloc ( nfiles, sizeof(double) );
  zaxisID2 = zaxisCreate(ZAXIS_GENERIC, nfiles);
  for ( i=0; i<nfiles; i++ )
    levs[i] = i;
  zaxisDefLevels(zaxisID2,levs);
  zaxisDefName(zaxisID2, "histogram_binID");
  time_mode = datafunc == TIME? TIME_VARIABLE : TIME_CONSTANT;
  for ( varID=0; varID<nvars; varID++) {

    /* **************************************************************** */
    /* nfiles includes the observation, so there are nfiles-1 ensembles */
    /* and exactly nfiles bins, in which the observation could fall     */
    /* **************************************************************** */

    if ( datafunc == TIME ) 
      {
	val = 0;
	gridID2 = gridCreate(GRID_LONLAT, 1);
	gridDefXsize(gridID2, 1);
	gridDefYsize(gridID2, 1);
	gridDefXvals(gridID2, &val);
	gridDefYvals(gridID2, &val);
      }
    else // datafunc == SPACE
      gridID2 = vlistInqVarGrid(vlistID1,varID);

    varID2[varID] = vlistDefVar(vlistID2, gridID2, zaxisID2, time_mode);
  }

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);
  

  for ( varID=0; varID< nvars; varID++ ){
    if ( zaxisInqSize(vlistInqVarZaxis(vlistID1, varID)) > 1 ) {
      cdoWarning("More than one level not supported when processing ranked histograms.");
      cdoWarning("Try to use `cdo splitlevel` to split the dataset into levels and apply");
      cdoWarning("the operator seperately to each level.");
      cdoAbort("Exit due to unsupported file structure");
    }
  } 

  streamID2 = streamOpenWrite(ofilename, cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", ofilename);

  streamDefVlist(streamID2, vlistID2);
	  
  gridsize = vlistGridsizeMax(vlistID1);

  for ( fileID = 0; fileID < nfiles; fileID++ )
    ef[fileID].array = (double *) malloc(gridsize*sizeof(double));

  if ( datafunc == SPACE ) 
    { /*need to memorize data for entire grid before writing          */
      array2 = (int **) malloc((nfiles+1)*sizeof(int*));
      for ( binID=0; binID<nfiles; binID++ ) 
	array2[binID] = (int*) calloc ( gridsize, sizeof(int) );
    }
  else /* data_func == TIME                                           */
    {  /* can process data separately for each timestep and only need */
       /* to cumulate values over the grid                            */
      array2    = (int**) malloc ( (nfiles+1)*sizeof(int *));
      for ( binID=0; binID<nfiles; binID++ )
	array2[binID] = (int *) calloc ( 1, sizeof(int) );
    }
  
  
  
  tsID = 0;
  do
    {
      nrecs0 = streamInqTimestep(ef[0].streamID, tsID);
      for ( fileID = 1; fileID < nfiles; fileID++ )
	{
	  streamID = ef[fileID].streamID;
	  nrecs = streamInqTimestep(streamID, tsID);
	  if ( nrecs != nrecs0 )
	    cdoAbort("Number of records changed from %d to %d", nrecs0, nrecs);
	}

      if ( datafunc == TIME || tsID == 0 ) {
	taxisCopyTimestep(taxisID2, taxisID1);
	if ( nrecs0 > 0 ) streamDefTimestep(streamID2, tsID);
      }

      fprintf(stderr,"TIMESTEP %i varID %i rec %i\n",tsID,varID,recID);
      
      for ( recID = 0; recID < nrecs0; recID++ )
	{
#if defined (_OPENMP)
#pragma omp parallel for default(shared) private(fileID, streamID, nmiss) \
  lastprivate(binID, varID, levelID)
#endif
	  for ( fileID = 0; fileID < nfiles; fileID++ )
	    {
	      streamID = ef[fileID].streamID;
	      streamInqRecord(streamID, &varID, &levelID);
	      streamReadRecord(streamID, ef[fileID].array, &nmiss);
	    }

	  gridID   = vlistInqVarGrid(vlistID1, varID);
	  gridsize = gridInqSize(gridID);
	  missval  = vlistInqVarMissval(vlistID1, varID);

	  nmiss = 0;
	  if ( datafunc == TIME ) 
	    for ( binID=0;binID<nfiles;binID++ )
	      array2[binID][0] = 0;
#if defined (_OPENMP)
#pragma omp parallel for default(shared) private(i, ompthID, fileID)
#endif
	  for ( i = 0; i < gridsize; i++ )
	    {
#if defined (_OPENMP)
	      ompthID = omp_get_thread_num();
#else
	      ompthID = 0;
#endif
	      field[ompthID].missval = missval;
	      field[ompthID].nmiss = 0;
	      have_miss = 0;
	      for ( fileID = 0; fileID < nfiles; fileID++ )
		{
		  field[ompthID].ptr[fileID] = ef[fileID].array[i];
		  if ( DBL_IS_EQUAL(field[ompthID].ptr[fileID], missval) ) 
		    {
		      have_miss = 1;
		      break;
		    }
		}
	      
	      // need to ignore all data for a gridpoint if a single ensemble
	      // has a missing value at that gridpoint.

	      for ( j=0; j<nfiles; j++ )
		fprintf(stderr,"%5.2g ",field[ompthID].ptr[j]);
	      binID = (int) fldfun(field[ompthID], operfunc);
	      fprintf(stderr,"-->%i\n",binID);

	      if ( datafunc == SPACE && ! have_miss) 
		array2[binID][i]++;
	      else if ( ! have_miss ) 
		array2[binID][0]++;
	    }

	  if ( datafunc == TIME ) {
	    for ( binID=0; binID<nfiles; binID++ ) {
	      val = (double)array2[binID][0];
	      fprintf(stderr,"%i ",(int)val);
	      streamDefRecord(streamID2, varID2[varID], binID);
	      streamWriteRecord(streamID2,&val, nmiss);
	    }
	  fprintf(stderr,"\n");
	  }

	}
      tsID++;
    }
  while ( nrecs0 > 0 );
  
  for ( binID=0; binID<nfiles; binID++ ) {
    double *tmpdoub = (double *)malloc (gridsize*sizeof(double));
    for(i=0; i<gridsize; i++ )
      tmpdoub[i] = (double) array2[binID][i];
    streamDefRecord(streamID2,varID2[varID],binID);
    streamWriteRecord(streamID2,tmpdoub,nmiss);
  }

  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      streamID = ef[fileID].streamID;
      streamClose(streamID);
    }

  streamClose(streamID2);

  for ( fileID = 0; fileID < nfiles; fileID++ )
    if ( ef[fileID].array ) free(ef[fileID].array);

  if ( ef ) free(ef);
  if ( array2 ) {
    for (binID=0; binID<nfiles; binID++ )
      free(array2[binID]);
    free(array2);
  }

  for ( i = 0; i < ompNumThreads; i++ )
    {
      if ( field[i].ptr    ) free(field[i].ptr);
      if ( field[i].weight ) free(field[i].weight);
    }
  
  if ( field ) free(field);
  
  cdoFinish();
  
  return (0);
}

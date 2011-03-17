/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2011 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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
   Ensstat5       enscrps          Ensemble cumulative ranked probability score
   Ensstat5       ensbrs           Ensemble Brier score
*/

#if defined (_OPENMP)
#  include <omp.h>
#endif

#include <cdi.h>
#include "cdo.h"
#include "statistic.h"
#include "cdo_int.h"
#include "pstream.h"
#include "util.h"
#include "merge_sort2.h"

enum OPERTYPE {CRPS, BRS};

void *Ensstat5(void *argument)
{
  int operatorID;
  int operfunc, datafunc;
  int i,k;
  int nvars,nrecs = 0, nrecs0, nmiss, nens, nfiles,nlevs,valcount;
  int cmpflag;
  int cum;
  int levelID, varID, recID, tsID, binID, ensID;
  int gridsize = 0;
  int gridID, gridID2;
  int have_miss = 0;
  int streamID = 0, streamID2;
  int vlistID, vlistID1, vlistID2;
  int taxisID1, taxisID2;
  int zaxisID1,zaxisID2;
  int ompthID;
  int *varID2;
  int xsize,ysize;
  double missval;
  double *alpha, *beta, *alpha_weights, *beta_weights;
  double *brs_g, *brs_o, *brs_g_weights, *brs_o_weights;
  double xval=0; double yval=0;
  double xa, *x;
  double *val;
  double *weights, sum_weights;
  double crps_reli, crps_pot,sprd, crps;
  double heavyside0, heavysideN;
  double brs_reli, brs_resol, brs_uncty, brs_thresh;
  double g,o,p;
  
  int fileID;
  const char *ofilename;

  typedef struct
  {
    int streamID;
    int vlistID;
    double *array;
  } ens_file_t;
  ens_file_t *ef = NULL;

  cdoInitialize(argument);
  
  cdoOperatorAdd("enscrps",   CRPS,  0,   NULL);
  cdoOperatorAdd("ensbrs",    BRS,   0,   NULL);
  
  operatorID = cdoOperatorID();
  operfunc = cdoOperatorF1(operatorID);
  datafunc = cdoOperatorF2(operatorID);

  nfiles = cdoStreamCnt() - 1;
  nens = nfiles-1;

  if ( operfunc == CRPS ) 
    cdoPrint("ensstat - crps");
  else if ( operfunc == BRS )  {
    operatorInputArg("Threshold for Brier score?");
    operatorCheckArgc(1);
    brs_thresh = atof(operatorArgv()[0]);
  }

  val = (double *) calloc ( nfiles,sizeof(double) );
  
  if ( operfunc == CRPS ) {
    alpha=(double *) calloc ( nens+1,sizeof(double) );
    beta =(double *) calloc ( nens+1,sizeof(double) );
    alpha_weights=(double *) calloc ( nens+1,sizeof(double) );
    beta_weights =(double *) calloc ( nens+1,sizeof(double) );
  }
  else if ( operfunc == BRS ) {
    brs_g = (double *) calloc ( nens+1,sizeof(double) );
    brs_o = (double *) calloc ( nens+1,sizeof(double) );
  }
  if ( cdoVerbose )
    cdoPrint("Ensemble over %d files (Ensstat5).", nfiles-1);

  ofilename = cdoStreamName(nfiles);

  if ( !cdoSilentMode && !cdoOverwriteMode )
    if ( fileExist(ofilename) )
      if ( !userFileOverwrite(ofilename) )
	cdoAbort("Outputfile %s already exist!", ofilename);

  ef = (ens_file_t *) malloc(nfiles*sizeof(ens_file_t));
  
  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      streamID = streamOpenRead(cdoStreamName(fileID));

      vlistID = streamInqVlist(streamID);
      
      ef[fileID].streamID = streamID;
      ef[fileID].vlistID = vlistID;
    }

  if ( cdoVerbose ) 
    cdoPrint("Opened %i Input Files for Ensemble Operator",nfiles);
  
  /* check for identical contents of all ensemble members */
  nvars = vlistNvars(ef[0].vlistID);
  fprintf(stderr,"nvars %i\n",nvars);
  if ( nvars == 1 ) 
    cmpflag = CMP_NAME | CMP_GRIDSIZE | CMP_NLEVEL | CMP_GRID;
  else 
    cdoAbort("Only single-variable files supported");

  for ( fileID = 1; fileID < nfiles; fileID++ )
    vlistCompare(ef[0].vlistID, ef[fileID].vlistID, cmpflag);

  vlistID1 = ef[0].vlistID;

  zaxisID1 = vlistInqVarZaxis(vlistID1,0);
  nlevs    = zaxisInqSize(zaxisID1);
  zaxisID2 = zaxisDuplicate(zaxisID1);

  vlistID2 = vlistCreate();
  varID2 = (int *) malloc( 2*sizeof(int));
  gridID2 = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID2, 1);
  gridDefYsize(gridID2, 1);
  gridDefXvals(gridID2, &xval);
  gridDefYvals(gridID2, &yval);

  if ( cdoVerbose )
    cdoPrint("Setup Grid for vlistID2");

  varID2[0] = vlistDefVar(vlistID2, gridID2, zaxisID2, TIME_VARIABLE);
  vlistDefVarName(vlistID2, varID2[0], "Reli");
  vlistDefVarLongname(vlistID2, varID2[0], "CRPS Reliability");
  vlistDefVarUnits(vlistID2, varID2[0], "1");
  
  varID2[1] = vlistDefVar(vlistID2, gridID2, zaxisID2, TIME_VARIABLE);
  vlistDefVarName(vlistID2, varID2[0], "CRPS_pot");
  vlistDefVarLongname(vlistID2, varID2[0], "Potential CRPS score");
  vlistDefVarUnits(vlistID2, varID2[0], "1");
  
  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);
  
  streamID2 = streamOpenWrite(ofilename, cdoFiletype());
  streamDefVlist(streamID2, vlistID2);
  
  gridsize = vlistGridsizeMax(vlistID1);
  weights=(double*) malloc (gridsize*sizeof(double));

  gridID = vlistInqVarGrid(vlistID1, 0);
  xsize = gridInqXsize(gridID);
  ysize = gridInqYsize(gridID);

  /*  if ( xsize > 1 && ysize > 1 )  {
    gridWeights(gridID, weights);
    sum_weights=0;
    for ( i=0; i<gridsize; i++ )  
      sum_weights += weights[i];
  }
  else*/ {
    for ( i=0; i< gridsize; i++ )
      weights[i] = 1./gridsize;
    sum_weights=1.;
  }

  if ( cdoVerbose ) 
    cdoPrint(" sum_weights %10.6f\n",sum_weights);

  for ( fileID = 0; fileID < nfiles; fileID++ )
    ef[fileID].array = (double *) malloc(gridsize*sizeof(double));
  
  tsID = 0;
  do
    {
      nrecs0 = streamInqTimestep(ef[0].streamID, tsID);
      for ( fileID = 1; fileID < nfiles; fileID++ )
	{
	  streamID = ef[fileID].streamID;
	  nrecs = streamInqTimestep(streamID, tsID);
	  if ( nrecs != nrecs0 )
	    cdoAbort("Number of records changed from %d to %d at time Step", nrecs0, nrecs, tsID);
	}

      taxisCopyTimestep(taxisID2, taxisID1);
      if ( nrecs0 > 0 ) streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs0; recID++ )
	{
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
	  valcount = 0;
	  heavyside0 = 0;
	  heavysideN = 0;

	  for ( i = 0; i < gridsize; i++ )
	    {
	      have_miss = 0;
	      for ( fileID = 0; fileID < nfiles; fileID++ )
		{
		  val[fileID] = ef[fileID].array[i];
		  if ( DBL_IS_EQUAL(val[fileID], missval) ) 
		    { have_miss = 1; break; }
		}

	      xa=val[0];                                     /* 1st file contains reference  */
	      x = &val[1];                                   /* Ensembles start at 2nd   file*/
	      sort_iter_single(nens,x,1);                    /* Sort The Ensemble Array      */

	      // only process if no missing value in ensemble
	      if ( ! have_miss && operfunc == CRPS )  
		{

		  if ( xa < x[0] ) {                             /* Consider outliers            */
		    beta[0] += (x[0]-xa)*weights[i];
		    heavyside0 += 1.;
		  }
		  if ( xa > x[nens-1] )  {
		    alpha[nens] += (xa-x[nens-1])*weights[i];
		    alpha_weights[nens] += weights[i];
		    heavysideN += 1.;
		  }

 		  /* Loop start at zero ==> 1st ensemble (c-indexing) */
		  for ( k=0; k<nens-1; k++ ) {                   /* Cumulate alpha and beta      */
		    if ( xa > x[k+1] )                           /* left of heavyside            */
		      alpha[k+1]+= (x[k+1]-x[k]) * weights[i];   /*                              */
		    else if ( xa < x[k] )                        /* right of heavyside           */
		      beta[k+1] += (x[k+1]-x[k]) * weights[i];   /*                              */
		    else if ( x[k+1] >= xa && xa >= x[k] ) {     /* hitting jump pf heavyside    */
		      alpha[k+1]+= (xa    -x[k]) * weights[i];   /* (occurs exactly once!)       */
		      beta[k+1] += (x[k+1]-xa)   * weights[i];   /* **************************** */
		    }
		  }
		}
	      else if ( operfunc == BRS ) 
		{
		  int occ = xa>brs_thresh? 1 : 0;
		  
		  if ( x[0] > brs_thresh ) 
		    brs_g[0] += weights[i];
		  else if ( x[nens-1] < brs_thresh ) 
		    brs_g[nens] += weights[i];		  
		  else
		    for ( k=0; k<nens-1; k++ ) {
		      if ( x[k+1] >= brs_thresh && brs_thresh >= x[k] ) {
			brs_g[k+1] += weights[i];
			break;
		      }
		    }

		  if ( x[0] > xa )
		    brs_o[0] += weights[i];
		  else if ( x[nens-1] < xa ) 
		    brs_o[nens] += weights[i];
		  else
		    for ( k=0; k<nens-1; k++ ) {
		      if ( x[k+1] >= xa && xa >= x[k] ) {
			brs_o[k+1] += weights[i];
			break;
		      }
		    }

		}
	    }        // for ( i=0; i<gridsize; i++ )
	  
	  if ( operfunc == CRPS ) {

	    // First Bin
	    p=0.; g=0.;
	    o = heavyside0/gridsize;
	    if ( o > 0. ) {
	      g = beta[0]/o;
	    }
	    crps_reli= g * (o - p) * (o - p);
	    crps_pot = g*o*(1.-o);
	    crps     = g*( (1.-o)*p*p  +  o*(1.-p)*(1.-p) );	  
	    
	    // Middle Bins
	    for ( k=1; k<nens; k++ ) {
	      p = (double)k/(double)nens;
	      
	      if ( ! DBL_IS_EQUAL(sum_weights,1.) ) {
		alpha[k] /= sum_weights; 
		beta[k] /= sum_weights;
	      }
	      
	      g = alpha[k]+beta[k];
	      o = beta[k] / ( alpha[k] + beta[k] ); 
	      
	      crps_reli    += g * (o - p) * (o - p);
	      crps_pot+= g*o*(1.-o);
	      crps    += g*( (1.-o)*p*p  +  o*(1.-p)*(1.-p) );	  
	    }
	    
	    // Last Bin
	    p=1.; g=0.;
	    o = 1. - heavysideN/gridsize; 
	    if ( o != 1. ) {
	      g = alpha[nens] / (1-o);
	      
	      crps_reli    += g * (o-p) * (o-p);
	      crps_pot+= g*o*(1-o);
	      crps    += g*( (1-o)*p*p  +  o*(1-p)*(1-p) );	  
	    }
	  } else if ( operfunc == BRS ) {
	    double gsum=0;
	    double obar=0; 
	    double osum=0;
	    double o,g,p;

	    brs_reli=0;
	    brs_resol=0;
	    brs_uncty=0;

	    for ( k=0; k<=nens; k++ ) {
	      obar += brs_g[k]*brs_o[k];
	      gsum += brs_g[k];
	      osum += brs_o[k];
	    }

	    if ( abs(osum)-1 > 1e-06 || abs(gsum)-1 > 1e-06 )  {
	      cdoAbort("Internal error - normalization constraint of problem not fulfilled");
	      cdoAbort("This is likely due to missing values");
	    }
	    o=0; p=0; g=0;
	    brs_uncty = obar * (1-obar);

	    for ( k=0; k<=nens; k++ ) {

	      g = brs_g[k];
	      o = brs_o[k];
	      p = k/(float)nens;

	      brs_reli += g * ( o-p ) * ( o-p );
	      brs_resol+= g * (o-obar) * (o-obar);
	      //fprintf(stderr,"%12.6g %12.6g %12.6g %12.6g %12.6g\n",obar,o,g,p,osum);
	      //fprintf(stderr,"%3i %12.6g %12.6g %12.6g %12.6g\n",k,brs_g[k], brs_reli,brs_resol, brs_uncty);
	    }

	    if ( cdoVerbose ) {
	      cdoPrint("Brier score for var %i level %i calculated",varID, levelID);
	      cdoPrint("obar %12.6g osum %12.6g gsum %12.6g",obar,osum,gsum);
	      cdoPrint("brs  %12.6g reli %12.6g resol %12.6g u %12.6g",
		       brs_reli-brs_resol+brs_uncty,brs_reli,brs_resol,brs_uncty);
	    }
	  }

	  
	  if ( cdoVerbose && operfunc == CRPS ) 
	    cdoPrint("CRPS:%12.6g reli:%12.6g crps_pot:%12.6g crps:%12.6g\n", 
		     crps,crps_reli,crps_pot, crps_reli+crps_pot);
	  if ( cdoVerbose && operfunc == BRS ) 
	    cdoPrint("BRS:");
	  
	  switch ( operfunc ) {
	  case ( CRPS ):
	    streamDefRecord(streamID2,varID2[0],levelID);
	    streamWriteRecord(streamID2,&crps_reli,have_miss);
	    streamDefRecord(streamID2,varID2[1],levelID);
	    streamWriteRecord(streamID2,&crps_pot,have_miss);
	    memset(alpha, 0, (nens+1) * sizeof(double) );
	    memset(beta,  0, (nens+1) * sizeof(double) ); 
	    heavyside0=0;
	    heavysideN=0;	  
	    break;
	  case ( BRS ):
	    cdoPrint("Write Brier score here\n");

	    memset(brs_o, 0, (nens+1)*sizeof(double) );
	    memset(brs_g, 0, (nens+1)*sizeof(double) );
	  }

	  
	}   // for ( recID = 0; recID < nrecs; recID++ ) 
      tsID++;
    }  while ( nrecs );
  
  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      streamID = ef[fileID].streamID;
      streamClose(streamID);
    }
  
  if ( operfunc != func_roc ) 
    streamClose(streamID2);
  
  for ( fileID = 0; fileID < nfiles; fileID++ )
    if ( ef[fileID].array ) free(ef[fileID].array);
  
  if ( ef ) free(ef);
  
  cdoFinish();
  
  return (0);
}

/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2015 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
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

     Timeof        eof             EOF in spatial or time space
     Timeof        eofspatial      EOF in spatial space
     Timeof        eoftime         EOF in time space
*/
/*
 * TODO: 
 * Role of the weights for eofs. Should not be mixed up with division with
 * number of contributing values during summation.
 */

//#define OLD_IMPLEMENTATION
#define WEIGHTS 1

#include <limits.h>  // LONG_MAX
#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"
#include "statistic.h"

enum T_EIGEN_MODE {JACOBI, DANIELSON_LANCZOS};

// NO MISSING VALUE SUPPORT ADDED SO FAR

void *EOFs(void * argument)
{
  char *envstr;

  enum {EOF_, EOF_TIME, EOF_SPATIAL};

  long i, ii, j, i1, i2, j1, j2;
  int nlevs=0 ;
  int nmiss;
  int tsID;
  int varID, recID, levelID;
  int nts=0;
  int n=0;
  int grid_space=0, time_space=0;
  int missval_warning=0;
  int timer_cov = 0, timer_eig = 0;

  int calendar = CALENDAR_STANDARD;
  juldate_t juldate;

  double sum;
  double missval=0;
  double xvals, yvals;
  double *df1p, *df2p;

  enum T_EIGEN_MODE eigen_mode = JACOBI;

  if ( cdoTimer )
    {
      timer_cov  = timer_new("Timeof cov");
      timer_eig  = timer_new("Timeof eig");
    }
  
  cdoInitialize(argument);

  cdoOperatorAdd("eof",        EOF_,        0, NULL);
  cdoOperatorAdd("eoftime",    EOF_TIME,    0, NULL);
  cdoOperatorAdd("eofspatial", EOF_SPATIAL, 0, NULL);

  int operatorID  = cdoOperatorID();
  int operfunc    = cdoOperatorF1(operatorID);

  operatorInputArg("Number of eigen functions to write out");
  int n_eig       = parameter2int(operatorArgv()[0]);

  envstr = getenv("CDO_SVD_MODE");
  
  if ( envstr && !strncmp(envstr,"danielson_lanczos",17) ) 
    eigen_mode = DANIELSON_LANCZOS;
  else if ( envstr && ! strncmp(envstr,"jacobi",6 ) )
    eigen_mode = JACOBI;
  else if ( envstr ) {
    cdoWarning("Unknown environmental setting %s for CDO_SVD_MODE. Available options are",envstr);
    cdoWarning("  - 'jacobi' for a one-sided parallelized jacobi algorithm");
    cdoWarning("  - 'danielson_lanzcos' for the D/L algorithm");
    envstr=NULL;
  }

  if ( cdoVerbose ) 
    cdoPrint("Using CDO_SVD_MODE '%s' from %s",
	     eigen_mode==JACOBI?"jacobi":"danielson_lanczos",
	     envstr?"Environment":" default");
  

  int streamID1  = streamOpenRead(cdoStreamName(0));
  int vlistID1   = streamInqVlist(streamID1);
  int taxisID1   = vlistInqTaxis(vlistID1);
  int gridID1    = vlistInqVarGrid(vlistID1, 0);
  long gridsize  = vlistGridsizeMax(vlistID1);
  int nvars      = vlistNvars(vlistID1);
  int nrecs      = vlistNrecs(vlistID1);

  double *weight = (double*) malloc(gridsize*sizeof(double));
  if ( WEIGHTS )
    gridWeights(gridID1, weight);
  else
    for( i = 0; i < gridsize; i++) weight[i] = 1;

  /*  eigenvalues */

  tsID = 0;

  /* COUNT NUMBER OF TIMESTEPS if EOF_ or EOF_TIME */
  if ( operfunc == EOF_ || operfunc == EOF_TIME )
    {
      if ( cdoVerbose )
	cdoPrint("Counting timesteps in ifile");
      
      while ( TRUE )
        {
          nrecs = streamInqTimestep(streamID1, tsID);
          if ( nrecs == 0 )  break;
          tsID++;
        }

      if ( cdoVerbose ) cdoPrint("Counted %i timeSteps", tsID);

      nts         = tsID;

      //TODO close on streamID1 ??  streamClose(streamID1);
      streamClose(streamID1);
      streamID1   = streamOpenRead(cdoStreamName(0));
      vlistID1    = streamInqVlist(streamID1);
      taxisID1    = vlistInqTaxis(vlistID1);
      if ( nts < gridsize || operfunc == EOF_TIME )
	{
	  time_space = 1;
	  grid_space = 0;
	}
      else
        {
          time_space = 0;
          grid_space = 1;
        }
    }
  else if ( operfunc == EOF_SPATIAL )
    {
      time_space = 0;
      grid_space = 1;
    }

  /* reset the requested number of eigen-function to the maximum if neccessary */
  if ( time_space )
    {
      if ( n_eig > nts )
        {
          cdoWarning("Solving in time-space:");
          cdoWarning("Number of eigen-functions to write out is bigger than number of time-steps.");
          cdoWarning("Setting n_eig to %i.", nts);
          cdoWarning("If You want to force a solution in grid-space use operator eofspatial");
          n_eig = nts;
        }
      n = nts;
    }
  else if ( grid_space )
    {
      if ( ((double)gridsize)*gridsize > (double)LONG_MAX )
	cdoAbort("Grid space to large!");

      if ( n_eig > gridsize )
        {
          cdoWarning("Solving in spatial space");
          cdoWarning("Number of eigen-functions to write out is bigger than grid size");
          cdoWarning("Setting n_eig to %i", gridsize);
          cdoWarning("If You want to force a solution in time-space use operator eoftime");
          n_eig = gridsize;
        }
      n = gridsize;
    }

  if ( cdoVerbose ) 
    cdoPrint("Calculating %i eigenvectors and %i eigenvalues in %s",
	     n_eig,n,grid_space==1?"grid_space" : "time_space");

  /* allocation of temporary fields and output structures */
  double *in              = (double *) malloc(gridsize*sizeof(double));
  int ***datacounts       = (int ***) malloc(nvars*sizeof(int **));
  double ****datafields   = (double ****) malloc(nvars*sizeof(double ***));
  double ****eigenvectors = (double ****) malloc(nvars*sizeof(double ***));
  double ****eigenvalues  = (double ****) malloc(nvars*sizeof(double ***));

  for ( varID = 0; varID < nvars; ++varID )
    {
      gridID1             = vlistInqVarGrid(vlistID1, varID);
      gridsize            = vlistGridsizeMax(vlistID1);
      nlevs               = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      missval             = vlistInqVarMissval(vlistID1, varID);

      datacounts[varID]   = (int **) malloc(nlevs*sizeof(int* ));
      datafields[varID]   = (double ***) malloc(nlevs*sizeof(double **));
      eigenvectors[varID] = (double ***) malloc(nlevs*sizeof(double **));
      eigenvalues[varID]  = (double ***) malloc(nlevs*sizeof(double **));

      for ( levelID = 0; levelID < nlevs; ++levelID )
        {
          if ( grid_space )
            {
              datafields[varID][levelID]    = (double **) malloc(1*sizeof(double *));
              datafields[varID][levelID][0] = (double*) malloc(gridsize*gridsize*sizeof(double));

              datacounts[varID][levelID]    = (int*) malloc(gridsize*gridsize*sizeof(int));
	      for ( i = 0; i < gridsize*gridsize; i++ )
		{
		  datacounts[varID][levelID][i] = 0;
		  datafields[varID][levelID][0][i] = 0;            
		}
	    }
          else if ( time_space )
            {
              datafields[varID][levelID] = (double **) malloc(nts*sizeof(double *));
              for ( tsID = 0; tsID < nts; tsID++ )
                {
                  datafields[varID][levelID][tsID] = (double *) malloc(gridsize*sizeof(double));
                  for ( i = 0; i < gridsize; ++i )
                    datafields[varID][levelID][tsID][i] = 0;
                }
              datacounts[varID][levelID] = (int *) malloc(gridsize*sizeof(int));	      
	      for ( i = 0; i < gridsize; i++ )
		datacounts[varID][levelID][i] = 0;
            }

          eigenvectors[varID][levelID] = (double **) malloc(n_eig*sizeof(double *));
          eigenvalues[varID][levelID]  = (double **) malloc(n*sizeof(double *));

          for ( i = 0; i < n; i++ )
            {
              if ( i < n_eig )
                {
                  eigenvectors[varID][levelID][i] = (double*) malloc(gridsize*sizeof(double));
                  for ( ii = 0; ii < gridsize; ++ii )
                    eigenvectors[varID][levelID][i][ii] = missval;
                }

              eigenvalues[varID][levelID][i] = (double*) malloc(1*sizeof(double));
              eigenvalues[varID][levelID][i][0]  = missval;
            }
        }
    }

  if ( cdoVerbose )
    cdoPrint("Allocated eigenvalue/eigenvector structures with nts=%i gridsize=%i", nts, gridsize);

  tsID = 0;

  /* read the data and create covariance matrices for each var & level */
  while ( TRUE )
    {
      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) break;

      for ( recID = 0; recID < nrecs; recID++ )
        {
          int i2;

          streamInqRecord(streamID1, &varID, &levelID);
          gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

          missval = vlistInqVarMissval(vlistID1, varID);
          streamReadRecord(streamID1, in, &nmiss);

          if ( grid_space )
            {
	      // This could be done in parallel to save lots of time
#if defined(_OPENMP)
#pragma omp parallel for private(i1,i2) default(shared)
#endif
              for ( i1 = 0; i1 < gridsize; i1++ )
                {
                  for ( i2 = i1; i2 < gridsize; i2++ )
                    {
                      if ( nmiss == 0 ||
			   (( ! DBL_IS_EQUAL(in[i1], missval) ) &&
			    ( ! DBL_IS_EQUAL(in[i2], missval) )) )
                        {
                          datafields[varID][levelID][0][i1*gridsize+i2] += in[i1]*in[i2];
                          datacounts[varID][levelID][i1*gridsize+i2]++;
                        }
                      else if ( missval_warning == 0 )
                        {
			  cdoWarning("Missing value support not checked for this operator!\n");
			  missval_warning = 1; 
			}
                    }
                }
            }
          else if ( time_space )
	    {
	      for ( i = 0; i < gridsize; ++i )
		{
		  if ( ! DBL_IS_EQUAL(in[i], missval ) )
		    {
		      datafields[varID][levelID][tsID][i] = in[i];
		      datacounts[varID][levelID][i]++;
		    }
		  else
		    {
		      if ( missval_warning == 0 )
			{
			  cdoWarning("Missing Value Support not Checked for this Operator!");
			  cdoWarning("Does not work with changing locations of missing values in time.");
			  missval_warning = 1;
			}
		      datafields[varID][levelID][tsID][i] = 0;
		    }
		}
	    }
        }
      tsID++;
    }

  if ( tsID == 1 ) cdoAbort("File consists of only one timestep!");

  if ( grid_space )
    for ( i1 = 0; i1 < gridsize; ++i1 )
      for ( i2 = 0; i2 < i1; ++i2 )
        {
          datafields[varID][levelID][0][i1*gridsize+i2] = datafields[varID][levelID][0][i2*gridsize+i1];
          datacounts[varID][levelID][i1*gridsize+i2]    = datacounts[varID][levelID][i2*gridsize+i1];
        }

  /*
  pack = (int*) malloc(gridsize*sizeof(int)); //TODO
  */

  for ( varID = 0; varID < nvars; varID++ )
    {
      char vname[256];
      vlistInqVarName(vlistID1, varID, vname);
      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
      nlevs    = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));

      if ( cdoVerbose )
	cdoPrint("Calculating cov matrices for %i levels of var%i (%s)", nlevs, varID, vname);

      for ( levelID = 0; levelID < nlevs; levelID++ )
        {
	  if ( cdoTimer ) timer_start(timer_cov);

	  if ( cdoVerbose ) cdoPrint("processing level %i",levelID);

          int i2;

	  double **datafieldv = datafields[varID][levelID];
	  int *datacountv = datacounts[varID][levelID];

          double **cov  = NULL; // TODO covariance matrix / eigenvectors after solving
          double *eigv  = NULL; // TODO eigenvalues

	  int *pack = (int *) malloc(gridsize*sizeof(int));
          int npack = 0;    // TODO already set to 0
          double sum_w = 0;

          if ( grid_space )
            {
              for ( i1 = 0; i1 < gridsize; i1++ )
                {
		  if ( datacountv[i1*gridsize + i1] > 1 ) 
		    pack[npack++] = i1;
                }

              for ( i1 = 0; i1 < npack; i1++ )
                sum_w += weight[pack[i1]];

	      n = npack;
	      if ( npack )
		{
		  cov = (double **) malloc(npack*sizeof(double *));
		  for ( i1 = 0; i1 < npack; i1++ )
		    cov[i1] = (double *) malloc(npack*sizeof(double));
		  eigv = (double *) malloc(npack*sizeof(double));
		}

              for ( i1 = 0; i1 < npack; i1++ )
		for ( i2 = i1; i2 < npack; i2++ )
		  if ( datacountv[pack[i1]*gridsize+pack[i2]] )
		    cov[i2][i1] = cov[i1][i2] =
		      datafieldv[0][pack[i1]*gridsize+pack[i2]]*   // covariance
		      sqrt(weight[pack[i1]]) * sqrt(weight[pack[i2]]) / sum_w /       // weights
		      (datacountv[pack[i1]*gridsize+pack[i2]]);   // number of data contributing
            }
          else if ( time_space )
            {
              for ( i = 0; i < gridsize; i++ )
                {
		  if ( datacountv[i] )
		    {
		      pack[npack] = i;
		      npack++;
		      sum_w += weight[i];
		    }
		}

	      if ( cdoVerbose )
		cdoPrint("allocating cov with %i x %i elements | npack=%i", nts, nts, npack);

              cov = (double **) malloc(nts*sizeof(double *));
              for ( j1 = 0; j1 < nts; j1++ )
                cov[j1] = (double *) malloc(nts*sizeof(double));
	      eigv = (double *) malloc(nts*sizeof(double));

#if defined(_OPENMP)
#pragma omp parallel for private(j1,j2,i,sum, df1p, df2p) default(shared) schedule(dynamic)
#endif
              for ( j1 = 0; j1 < nts; j1++ )
		for ( j2 = j1; j2 < nts; j2++ )
		  {
		    sum = 0;
		    df1p = datafieldv[j1];
		    df2p = datafieldv[j2];
		    for ( i = 0; i < npack; i++ )
		      {
			sum += weight[pack[i]]*df1p[pack[i]]*df2p[pack[i]];
		      }
		    cov[j2][j1] = cov[j1][j2] = sum / sum_w / nts;
		  }

	      if ( cdoVerbose )
		cdoPrint("finished calculation of cov-matrix for var %s",&vname[0]);
            }

	  if ( cdoTimer ) timer_stop(timer_cov);

          /* SOLVE THE EIGEN PROBLEM */
	  if ( cdoTimer ) timer_start(timer_eig);

	  if ( eigen_mode == JACOBI ) 
	    // TODO: use return status (>0 okay, -1 did not converge at all) 
	    parallel_eigen_solution_of_symmetric_matrix(cov, eigv, n, __func__);
	  else 
	    eigen_solution_of_symmetric_matrix(cov, eigv, n, __func__);

	  if ( cdoTimer ) timer_stop(timer_eig);
	  /* NOW: cov contains the eigenvectors, eigv the eigenvalues */
	  
          for ( i = 0; i < n; i++ ) 
            eigenvalues[varID][levelID][i][0] = eigv[i]*sum_w;

          for ( i = 0; i < n_eig; i++ )
            {
	      double *eigenvec = eigenvectors[varID][levelID][i];

              if ( grid_space )
		{
		  /*
		  double sumw= 0;
		  for ( j = 0; j < npack; j++ ) sumw += sqrt(weight[pack[j]]);
		  printf("sumw %g\n", sumw);
		  sumw= 0;
		  for ( j = 0; j < npack; j++ ) sumw += datafieldv[0][pack[j]];
		  printf("sumw %g\n", sumw);
		  */
		  for ( j = 0; j < npack; j++ )
		    eigenvec[pack[j]] = 
#ifdef OLD_IMPLEMENTATION
		      cov[i][j] / sqrt(weight[pack[j]]);
#else
		      cov[i][j] /*/ sqrt(weight[pack[j]])*/;
		  // cov[i][j] / (sqrt(weight[pack[j]])*29.15);
#endif
		}
              else if ( time_space )
                {
#if defined(_OPENMP)
#pragma omp parallel for private(i2,j,sum) shared(datafieldv, eigenvec)
#endif
                  for ( i2 = 0; i2 < npack; i2++ )
                    {
                      sum = 0;
                      for ( j = 0; j < nts; j++ )
                        sum += datafieldv[j][pack[i2]] * cov[i][j];

                      eigenvec[pack[i2]] = sum;
                    }

                  // NORMALIZING
                  sum = 0;

#if defined(_OPENMP)
#pragma omp parallel for private(i2) default(none) reduction(+:sum)	\
  shared(eigenvec,weight,pack,npack)
#endif
                  for ( i2 = 0; i2 < npack; i2++ )
		    {
		      /* 
		      ** do not need to account for weights as eigenvectors are non-weighted                                   
		      */ 
#ifdef OLD_IMPLEMENTATION
		      sum += weight[pack[i2]] *
#else
		      sum += /*weight[pack[i2]] **/
#endif
			eigenvec[pack[i2]] * eigenvec[pack[i2]];
		    }

                  if ( sum > 0 )
                    {
                      sum = sqrt(sum);
#if defined(_OPENMP)
#pragma omp parallel for private(i2) default(none)  shared(npack,pack,sum,eigenvec)
#endif
                      for( i2 = 0; i2 < npack; i2++ ) eigenvec[pack[i2]] /= sum;
                    }
                  else
		    {
#if defined(_OPENMP)
#pragma omp parallel for private(i2) default(none)  shared(npack,pack,eigenvec,missval)
#endif
		      for( i2 = 0; i2 < npack; i2++ ) eigenvec[pack[i2]] = missval;
		    }
                } // else if ( time_space )
            } // for ( i = 0; i < n_eig; i++ )

	  if ( eigv ) free(eigv);
	  for ( i=0; i<n; i++ ) if ( cov[i] ) free(cov[i]);
	  if ( cov ) free(cov);
	  if ( pack ) free(pack);

        } // for ( levelID = 0; levelID < nlevs; levelID++ )
    } // for ( varID = 0; varID < nvars; varID++ )


  /* write files with eigenvalues (ID3) and eigenvectors (ID2) */

  /* eigenvalues */
  int streamID2   = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  int vlistID2    = vlistDuplicate(vlistID1);
  int taxisID2    = taxisDuplicate(taxisID1);
  taxisDefRdate(taxisID2, 0);
  taxisDefRtime(taxisID2, 0);
  vlistDefTaxis(vlistID2, taxisID2);
  int gridID2     = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID2, 1);
  gridDefYsize(gridID2, 1);
  xvals    = 0;
  yvals    = 0;
  gridDefXvals(gridID2, &xvals);
  gridDefYvals(gridID2, &yvals);
  int ngrids   = vlistNgrids(vlistID2);
  for ( i = 0; i < ngrids; i++ )
    vlistChangeGridIndex(vlistID2, i, gridID2);

  /*  eigenvectors */
  int streamID3   = streamOpenWrite(cdoStreamName(2), cdoFiletype());

  int vlistID3    = vlistDuplicate(vlistID1);
  int taxisID3    = taxisDuplicate(taxisID1);
  taxisDefRdate(taxisID3, 0);
  taxisDefRtime(taxisID3, 0);
  vlistDefTaxis(vlistID3, taxisID3);

  streamDefVlist(streamID2, vlistID2);
  streamDefVlist(streamID3, vlistID3);

  int vdate = 10101;
  int vtime = 0;
  juldate = juldate_encode(calendar, vdate, vtime);
  for ( tsID = 0; tsID < n; tsID++ )
    {
      juldate = juldate_add_seconds(60, juldate);
      juldate_decode(calendar, juldate, &vdate, &vtime);

      taxisDefVdate(taxisID2, vdate);
      taxisDefVtime(taxisID2, vtime);
      streamDefTimestep(streamID2, tsID);

      if ( tsID < n_eig )
        {
          taxisDefVdate(taxisID3, vdate);
          taxisDefVtime(taxisID3, vtime);
          streamDefTimestep(streamID3, tsID);
        }

      for ( varID = 0; varID < nvars; varID++ )
        {
          nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
          for ( levelID = 0; levelID < nlevs; levelID++ )
            {
              if ( tsID < n_eig )
                {
                  nmiss = 0;
                  for ( i = 0; i < gridsize; i++ )
                    if ( DBL_IS_EQUAL(eigenvectors[varID][levelID][tsID][i], missval) ) nmiss++;

                  streamDefRecord(streamID3, varID, levelID);
                  streamWriteRecord(streamID3, eigenvectors[varID][levelID][tsID], nmiss);
                }

              if ( DBL_IS_EQUAL(eigenvalues[varID][levelID][tsID][0], missval) ) nmiss = 1;
              else nmiss = 0;
              streamDefRecord(streamID2, varID, levelID);
              streamWriteRecord(streamID2, eigenvalues[varID][levelID][tsID], nmiss);
            }
        }
    }
  
  for ( varID = 0; varID < nvars; varID++)
    {
      nlevs    = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      gridsize =  gridInqSize(vlistInqVarGrid(vlistID1, varID));
      
      for( levelID = 0; levelID < nlevs; levelID++ )
        {
	  int n_use = time_space == 1? nts : gridsize;
          for(i = 0; i < n_use; i++)
            {
              if ( i < n_eig ) 
                if (eigenvectors[varID][levelID][i])
		  free(eigenvectors[varID][levelID][i]);
	      if (eigenvalues[varID][levelID][i])
		free(eigenvalues[varID][levelID][i]);
	    }
	  
	  if ( grid_space ) 
	    free(datafields[varID][levelID][0]);
	  else if ( time_space )
	    for (tsID=0; tsID<nts; tsID++ )
	      free(datafields[varID][levelID][tsID]);

          free(datacounts[varID][levelID]);
	  free(datafields[varID][levelID]);
          free(eigenvectors[varID][levelID]);
          free(eigenvalues[varID][levelID]);
        }
      
      free(datafields[varID]);
      free(datacounts[varID]);
      free(eigenvectors[varID]);
      free(eigenvalues[varID]);
    }

  free(datafields);
  free(datacounts);
  free(eigenvectors);
  free(eigenvalues);
  free(in);
  free(weight);


  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);

  vlistDestroy(vlistID2);
  vlistDestroy(vlistID3);

  gridDestroy(gridID1);
  gridDestroy(gridID2);

  //  taxisDestroy(taxisID2);
  //  taxisDestroy(taxisID3);

  cdoFinish();

  return (0);
}

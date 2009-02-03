/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2007 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

*/


#include <string.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "statistic.h"

void *Timeof(void *argument)
{
  static char func[] = "Timeof";
  int operatorID;
  int streamID1, streamID2, streamID3;
  int gridsize;
  int vdate = 0, vtime = 0;
  int nrecs, nvars, nlevs ;
  int i, j, i1, i2;
  int nmiss;
  int tsID;
  int varID, recID, levelID, gridID;
  int vlistID1, vlistID2 = -1, vlistID3 = -1;
  int taxisID1, taxisID2, taxisID3;
  int gridID1, gridID2, gridID3;
  int ngrids;
  int reached_eof;
  int npack;
  int *pack;
  int ***iwork;
  double *w;
  double sum_w;
  double missval;
  double *xvals, *yvals;
  FIELD **fwork;
  FIELD ***o, ***o2; 
  FIELD in;     

  cdoInitialize(argument);
  cdoOperatorAdd("timeof", 0, 0, NULL);  
  operatorID = cdoOperatorID();
  printf("initialized operator %i\n", operatorID);
  
  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));
  vlistID1 = streamInqVlist(streamID1); 
  taxisID1 = vlistInqTaxis(vlistID1);
  printf("opened stream1 (%i) for reading\n", streamID1);
  
  gridsize = vlistGridsizeMax(vlistID1);
  nvars = vlistNvars(vlistID1);
  nrecs = vlistNrecs(vlistID1); 
  taxisID1 = vlistInqTaxis(vlistID1);
  w = (double*)malloc(gridsize*sizeof(double));
  gridWeights(gridID1, &w[0]);
  
  printf("got stream1 properties\n");
  
  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));
  vlistID2 = vlistDuplicate(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1); 
  vlistDefTaxis(vlistID2, taxisID2);  
  
  printf("set stream to properties\n");
  
  streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());
  if ( streamID3 < 0 ) cdiError(streamID3, "Open failed on %s", cdoStreamName(2));
  vlistID3 = vlistDuplicate(vlistID1); 
  taxisID3 = taxisDuplicate(taxisID1);  
  vlistDefTaxis(vlistID3, taxisID3);  
  
  gridID3 = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID3, 1);
  gridDefYsize(gridID3, 1);
  xvals=(double*)malloc(1*sizeof(double));
  yvals=(double*)malloc(1*sizeof(double));
  xvals[0]=0;
  yvals[0]=0;
  gridDefXvals(gridID3, xvals);
  gridDefYvals(gridID3, yvals);
  
  ngrids = vlistNgrids(vlistID3);
  
  for ( i = 0; i < ngrids; i++ )
    vlistChangeGridIndex(vlistID3, i, gridID3);
  
  reached_eof=0;
  
  in.ptr = (double *) malloc(gridsize*sizeof(double));

  fwork = (FIELD **) malloc(nvars*sizeof(FIELD*));
  o     = (FIELD ***)malloc(nvars*sizeof(FIELD**));
  o2    = (FIELD ***)malloc(nvars*sizeof(FIELD**));
  iwork = (int ***) malloc(nvars*sizeof(int **));
  

  for ( varID = 0; varID < nvars; ++varID )
    {
      gridID1  = vlistInqVarGrid(vlistID1, varID);      
      gridsize = vlistGridsizeMax(vlistID1);
      nlevs    = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      missval  = vlistInqVarMissval(vlistID1, varID);      
      
      fwork[varID] = (FIELD *)  malloc(nlevs*sizeof(FIELD));
      o[varID]     = (FIELD **) malloc(nlevs*sizeof(FIELD*));     
      o2[varID]    = (FIELD **) malloc(nlevs*sizeof(FIELD*));
      iwork[varID] = (int ** )  malloc(nlevs*sizeof(int* ));

      for ( levelID = 0; levelID < nlevs; ++levelID )
        {  
          fwork[varID][levelID].grid    = gridID1;
          fwork[varID][levelID].nmiss   = 0;
          fwork[varID][levelID].missval = missval;
          fwork[varID][levelID].ptr     = (double *) malloc(gridsize*gridsize*sizeof(double));
          iwork[varID][levelID] = (int *) malloc(gridsize*gridsize*sizeof(int));
          memset(fwork[varID][levelID].ptr, missval, gridsize*gridsize*sizeof(double));
          memset(iwork[varID][levelID], 0, gridsize*gridsize*sizeof(int)); 
          
          o[varID][levelID] = (FIELD *) malloc(gridsize*sizeof(FIELD));
          o2[varID][levelID]= (FIELD *) malloc(gridsize*sizeof(FIELD));
          for ( i = 0; i < gridsize; i++ )
            {
              o[varID][levelID][i].grid   = gridID2;
              o[varID][levelID][i].nmiss  = 0;
              o[varID][levelID][i].missval= missval;
              o[varID][levelID][i].ptr    = (double *)malloc(gridsize*sizeof(double));
              memset(o[varID][levelID][i].ptr, missval, gridsize*sizeof(double));
              
              o2[varID][levelID][i].grid    = gridID3;
              o2[varID][levelID][i].nmiss   = 0;
              o2[varID][levelID][i].missval = missval;
              o2[varID][levelID][i].ptr     = (double *)malloc(1*sizeof(double));
              o2[varID][levelID][i].ptr[0]  = missval;              
            }
        }      
    }  
  
  tsID=0; 
  /* Read the data and reate covariance matrices for each var & level */
  while ( TRUE )
    {      
      if ( reached_eof ) continue;
      
      nrecs = streamInqTimestep(streamID1, tsID);
      //printf("nrecs %i\n",nrecs);
      if ( nrecs == 0 )
        {
          reached_eof = 1;
          break;
        }
          
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);    
      
      for ( recID = 0; recID < nrecs; recID++ )
        {
          int i2;
          streamInqRecord(streamID1, &varID, &levelID); 
          gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));          
          missval = in.missval = vlistInqVarMissval(vlistID1, varID);           
          
          streamReadRecord(streamID1, in.ptr, &in.nmiss);   
          for ( i = 0; i < gridsize; ++i )
            {
              double tmp;
              for ( i2 = 0; i2 <= i; i2++ )
                {                
                  if ( ( ! DBL_IS_EQUAL(in.ptr[i], missval) ) && 
                      ( ! DBL_IS_EQUAL(in.ptr[i2], missval) ) )                            
                    {
                      tmp = in.ptr[i]*in.ptr[i2];                    
                      fwork[varID][levelID].ptr[i*gridsize+i2] += tmp;                      
                      iwork[varID][levelID][i*gridsize+i2]++;                                        
                    }
                }              
            }	   
          // creater lower triangular of covariance matrix;
          for (i=0;i<gridsize;i++)
            for(i2=0;i2<i;i2++)
              {
                fwork[varID][levelID].ptr[i2*gridsize+i]=fwork[varID][levelID].ptr[i*gridsize+i2];
                iwork[varID][levelID][i2*gridsize+i]    =iwork[varID][levelID][i*gridsize+i2];
              }
              
        }                                  
      tsID++;
    }
  printf("read data\n");
  
  pack = (int *)malloc(gridsize*sizeof(int));   
  
  for ( varID = 0; varID < nvars; varID++ )
    {
      for ( levelID = 0; levelID < nlevs; levelID++ )
        {
          double **cov;
          double *eigv;
          npack = 0;
          int i2;
          sum_w = 0;
          memset(pack, 0, gridsize*sizeof(int));
          for ( i = 0; i < gridsize; i++ )
            {        
              for ( i2 = 0; i2 < gridsize; i2++ )              
                  if (!iwork[varID][levelID][i2*gridsize+i])
                    break;
              
              if ( i2 == gridsize )
                {
                  pack[npack] = i;                  
                  npack++;        
                  sum_w += w[i];
                }            
            }          
          cov = (double **)malloc(npack*sizeof(double *));
          eigv = (double *)malloc(npack*sizeof(double));
          
          for (i1=0;i1<npack;i1++) 
            { 
              cov[i1]=(double*)malloc(npack*sizeof(double));
              for (i2 = 0; i2 < npack; i2++ )                               
                if ( iwork[varID][levelID][i1*gridsize+i2] )
                  {                    
                    cov[i1][i2] = fwork[varID][levelID].ptr[pack[i2]*gridsize+pack[i1]];                                   
                    // weight and normalize the covariances 
                    cov[i1][i2] *= sqrt (w[pack[i1]]) * sqrt (w[pack[i2]]) / sum_w / (iwork[varID][levelID][i1*gridsize+i2] - 1);
                  }
            }   
          printf("sumw %7.5f, weight[0] %7.5f weights[npack-1] %7.5f\n", sum_w, w[0], w[npack-1]);
          
           for(i1=0;i1<npack;i1++)
            {
              printf("\n");
              for(i2=0;i2<npack;i2++)
                printf("%7.5f ", cov[i1][i2]);
              printf("\n");
            }
          
          
          eigen_solution_of_symmetric_matrix(&cov[0], &eigv[0], npack, func);
          // cov contains the eigenvectors, eigv the eigenvalues
          for (i=0;i<npack;i++)
            {             
              for (j=0;j<npack;j++)
                o[varID][levelID][i].ptr[pack[j]] = cov[j][i]/sqrt(w[pack[j]]/sum_w);                                
              o2[varID][levelID][i].ptr[0] = eigv[i];
            }
        }
    }
    
  streamDefVlist(streamID3, vlistID3);
  streamDefVlist(streamID2, vlistID2);
  
  for ( tsID=0; tsID<npack; tsID++ )
    {
      taxisDefVdate(taxisID3, 0.0); 
      taxisDefVtime(taxisID3, 0.0);
      streamDefTimestep(streamID3, tsID); 
     
      taxisDefVdate(taxisID2, 0.0);
      taxisDefVtime(taxisID2, 0.0);  
      streamDefTimestep(streamID2, tsID); 
       
      for ( varID = 0; varID < nvars; varID++ )
        {
          nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
          for ( levelID=0; levelID < nlevs;levelID++ )
            {           
              nmiss = 0;              
              for ( i = 0; i < gridsize; i++ )
                if ( DBL_IS_EQUAL(o[varID][levelID][tsID].ptr[i], missval) ) nmiss++;              
              streamDefRecord(streamID2, varID, levelID);              
              streamWriteRecord(streamID2, o[varID][levelID][tsID].ptr, nmiss);
              
              if ( DBL_IS_EQUAL(o2[varID][levelID][tsID].ptr[i], missval) ) nmiss = 1;
              else nmiss = 0;
              streamDefRecord(streamID3, varID, levelID);
              streamWriteRecord(streamID3, o2[varID][levelID][tsID].ptr,nmiss);
              
            }     
        }
    }
          
  for ( varID=0;varID<nvars;varID++)
    {
      for(levelID=0;levelID<nlevs;levelID++)
        {
          for(i=0;i<gridsize;i++)
            {
              free(o[varID][levelID][i].ptr);
              free(o2[varID][levelID][i].ptr);
            }
          free(o[varID][levelID]);
          free(o2[varID][levelID]);
          free(iwork[varID][levelID]);
        }
      free(o[varID]);
      free(o2[varID]);
      free(fwork[varID]);
      free(iwork[varID]);

    }
  
  free(o);
  free(o2);
  free(fwork);
  free(iwork);
  free(in.ptr);
  
  
  
  
  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);
  
  cdoFinish();   
  
  return (0);  
}

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
 
 Eofcoeff             eofcoeff             process eof coefficients
*/
#define WEIGHTS 1

#include <string.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

// NO MISSING VALUE SUPPORT ADDED SO FAR

void *Eofcoeff(void * argument)
{
  static char func[] = "Eofcoeff";
  int operatorID;
  int operfunc;  
  int gridsize;
  int vdate = 0, vtime = 0;
  int i, ii, j, i1, i2, j1, j2;
  int nts, nrecs, nvars, nlevs, neof, n, nchars, nmiss, ngrids; 
  int varID, recID, levelID, tsID, eofID;
  int vlistID1, vlistID2, taxisID1, taxisID2, gridID1, gridID2, gridID3, streamID1, streamID2;
  int  *vlistIDs, *taxisIDs, *streamIDs, *varIDs;
  int reached_eof;
  char eof_name[5], oname[1024], filesuffix[32];
  double *w;
  double sum_w;
  double missval1, missval2;
  double *xvals, *yvals;
  FIELD ***eof;  
  FIELD in;  
  FIELD out;
  
  cdoInitialize(argument);
  cdoOperatorAdd("eofcoeff",  0,       0, NULL);
  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);
     
  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));
  vlistID1 = streamInqVlist(streamID1); 
  taxisID1 = vlistInqTaxis(vlistID1);  
  gridsize = vlistGridsizeMax(vlistID1);
  gridID1 = vlistInqVarGrid(vlistID1, 0);
  nvars = vlistNvars(vlistID1);
  nrecs = vlistNrecs(vlistID1); 
  nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, 0));
  taxisID1 = vlistInqTaxis(vlistID1);
  w = (double*)malloc(gridsize*sizeof(double));
  gridWeights(gridID1, &w[0]);
  
  streamID2 = streamOpenRead(cdoStreamName(1));
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));
  vlistID2 = streamInqVlist(vlistID1);
  taxisID2 = vlistInqTaxis(vlistID1); 
  if ( vlistGridsizeMax(vlistID2) != gridsize ||
      vlistInqVarGrid(vlistID2, 0) || gridID1 )
    cdoAbort("EOFs (%s) and data (%s) defined on different grids", cdoStreamName(0), cdoStreamName(1));    
 
  
  strcpy(oname, cdoStreamName(2));
  nchars = strlen(oname);
  
  filesuffix[0] = 0;
  if ( cdoDisableFilesuffix == FALSE )
    {
      strcat(filesuffix, streamFilesuffix(cdoDefaultFileType));
      if ( cdoDefaultFileType == FILETYPE_GRB )
        if ( vlistIsSzipped(vlistID1) || cdoZtype == COMPRESS_SZIP )
          strcat(filesuffix, ".sz");
    }
  reached_eof=0;
  
  eof = (FIELD ***) malloc (nvars * sizeof(FIELD**) );
  for ( varID=0; varID<nvars; varID++)
    eof[varID] = (FIELD **) malloc(nlevs*sizeof(FIELD*));
  
  while ( 1 )       
   {     
     nrecs = streamInqTimestep(streamID1, eofID);
     if ( nrecs == 0)
       {
         reached_eof = 1;
         break;
       }
     printf("reading %i. eof\n", eofID+1); 
     for ( recID = 0; recID < nrecs; recID++ )
       {         
         streamInqRecord(streamID1, &varID, &levelID);
         missval1 = vlistInqVarMissval(vlistID1, varID);
         if ( eofID == 0 )
           eof[varID][levelID] = (FIELD*) malloc (1*sizeof(FIELD));
         else
           eof[varID][levelID] = (FIELD*) realloc (eof[varID][levelID], (eofID+1)*sizeof(FIELD));
         eof[varID][levelID][eofID].grid   = gridID1;
         eof[varID][levelID][eofID].nmiss  = 0;
         eof[varID][levelID][eofID].missval= missval1;
         eof[varID][levelID][eofID].ptr    = (double*)malloc(gridsize*sizeof(double));
         memset(&eof[varID][levelID][eofID].ptr[0], missval1, gridsize*sizeof(double));
         if ( varID >= nvars )
           cdoAbort("Internal error - too high varID");
         if ( levelID >= nlevs )
           cdoAbort("Internal error - too high levelID");
         
         streamReadRecord(streamID1, eof[varID][levelID][eofID].ptr, 
                          &eof[varID][levelID][eofID].nmiss);
       }
     eofID++;
   }
  neof = eofID;
  printf("read EOFs\n");
  
  gridID3 = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID3, 1);
  gridDefYsize(gridID3, 1);
  xvals=(double*)malloc(1*sizeof(double));
  yvals=(double*)malloc(1*sizeof(double));
  xvals[0]=0;
  yvals[0]=0;
  gridDefXvals(gridID3, xvals);
  gridDefYvals(gridID3, yvals);
  
  streamIDs = (int *) malloc (neof*sizeof(int));
  taxisIDs  = (int *) malloc (neof*sizeof(int));
  vlistIDs  = (int *) malloc (neof*sizeof(int));  
  eofID = 0;
  for ( eofID = 0; eofID < neof; eofID++)
    {
      oname[nchars] = '\0';                       
      for(varID=0; varID<nvars; varID++) 
      sprintf(eof_name, "%5.5i", eofID);
      strcat(oname, eof_name);
      if ( filesuffix[0] )
        strcat(oname, filesuffix);
      
      streamIDs[eofID] = streamOpenWrite(oname, cdoFiletype());
      if ( streamIDs[eofID] < 0) cdiError(streamIDs[eofID], "Open failed on %s", oname);
      else fprintf(stderr, "opened %s for writing\n", oname);
      
      streamDefVlist(streamIDs[eofID], vlistID2);
      vlistIDs[eofID] = vlistDuplicate(vlistID2); 
      taxisIDs[eofID] = taxisDuplicate(taxisID2);  
      vlistDefTaxis(vlistIDs[eofID], taxisIDs[eofID]);  
                   
      ngrids = vlistNgrids(vlistIDs[eofID]);
      for ( i = 0; i < ngrids; i++ )
        vlistChangeGridIndex(vlistIDs[eofID], i, gridID3);
      printf("eofID%i gridID%i", eofID, gridID3);
 
      
    }
  
  // ALLOCATE temporary field for data reading
  in.ptr = (double*) malloc(gridsize*sizeof(double));
  in.grid = gridID1;
  
  out.missval = missval1;
  out.nmiss = 0;
  out.ptr = (double *) malloc (1*sizeof(double));

  printf("allocated field for reading\n");
  reached_eof=0;
  tsID=0;
  while ( 1 )
    {
      printf("processing timestep %i\n", tsID+1);
      nrecs = streamInqTimestep(streamID2, tsID);
      if ( nrecs == 0 )
        {
          reached_eof=1;
          break;
        }
      
      for ( recID =0; recID< nrecs; recID ++ )
        {
          streamInqRecord(streamID2, &varID, &levelID);
          missval2 = vlistInqVarMissval(vlistID2, varID);
          streamReadRecord(streamID2, in.ptr, &in.nmiss);  
          
          for ( eofID = 0; eofID < neof; eofID++ )
            {
              out.ptr[0]  = 0;
              out.grid    = gridID3;
              out.missval = missval2;            
              for(i=0;i<gridsize;i++)
                {
                  //printf("%i\n",i);
                  if ( in.ptr[i] != missval2 && 
                      eof[varID][levelID][eofID].ptr[i] != missval1 )
                    {
                      double tmp = w[i]*in.ptr[i]*eof[varID][levelID][eofID].ptr[i];
                      out.ptr[0] += tmp;
                      //printf("%7.5f %7.5f %7.5f\n", in.ptr[i], eof[varID][levelID][eofID].ptr[i], tmp);
                    }
                }
              printf("grid %i\n", out.grid);
              printf("stream%i ts%i var%i level%i eof%i %17.15f\n", streamIDs[eofID], tsID, varID, levelID, eofID, out.ptr[0]); 
              if ( out.ptr[0] ) nmiss=0;
              else { nmiss=1; out.ptr[0]=missval2; }
              
              //printf("writing nmiss=%i streamID %i \n",nmiss, streamIDs[eofID]);
              vdate = taxisInqVdate(taxisID2);
              vtime = taxisInqVtime(taxisID2);
              taxisDefVdate(taxisIDs[eofID], vdate);
              taxisDefVtime(taxisIDs[eofID], vtime);
              
              streamDefTimestep(streamIDs[eofID],tsID);
              streamDefRecord(streamIDs[eofID], varID, levelID);
              
              // TEST FOR GRID
              int tmp = vlistInqVarGrid(vlistIDs[eofID], varID);
              printf("grid %i\n", tmp);
              streamWriteRecord(streamIDs[eofID],out.ptr,nmiss);
              //sprintf("finished\n");
            }
          if ( varID >= nvars )
            cdoAbort("Internal error - too high varID");
          if ( levelID >= nlevs )
            cdoAbort("Internal error - too high levelID");                              
        }
      tsID++;
    }
}


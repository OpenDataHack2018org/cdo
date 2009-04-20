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
  char eof_name[5], oname[1024], filesuffix[32];
  double *w;
  double missval1, missval2;
  double *xvals, *yvals;  
  FIELD ***eof;  
  FIELD in;  
  FIELD out;
  int operatorID;
  int operfunc;  
  int gridsize;
  int i, varID, recID, levelID, tsID, eofID;    
  int gridID1, gridID2, gridID3;
  int nrecs, nvars, nlevs, neof, nchars, nmiss, ngrids; 
  int reached_eof;
  int streamID1, streamID2, *streamIDs;
  int taxisID1, taxisID2, taxisID3;
  int vlistID1, vlistID2, vlistID3;
  int vdate = 0, vtime = 0;

 
  
  cdoInitialize(argument);
  cdoOperatorAdd("eofcoeff",  0,       0, NULL);
  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);
     
  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));
  streamID2 = streamOpenRead(cdoStreamName(1));
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));
  
  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = streamInqVlist(vlistID1);
  
  taxisID1 = vlistInqTaxis(vlistID1);  
  taxisID2 = vlistInqTaxis(vlistID1); 

  gridID1 = vlistInqVarGrid(vlistID1, 0);
  gridID2 = vlistInqVarGrid(vlistID2, 0);
  
  gridsize ?= vlistGridsizeMax(vlistID1)==vlistGridsizeMax(vlistID2)? vlistInqGridsizeMax(vlistID1) : -1;  
  nvars    ?= vlistNvars(vlistID1)==vlistNvars(vlistID2)? vlistNvars(vlistID1): -1;
  nrecs = vlistNrecs(vlistID1); 
  nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, 0));
  taxisID1 = vlistInqTaxis(vlistID1);
  w = (double*)malloc(gridsize*sizeof(double));
  gridWeights(gridID2, &w[0]);
  
  
  
   if (vlistGridsizeMax(vlistID2)   != gridsize ||
      vlistInqVarGrid(vlistID2, 0) != gridID1 )
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
 
  
  eof = (FIELD ***) malloc (nvars * sizeof(FIELD**) );
  for ( varID=0; varID<nvars; varID++)
    eof[varID] = (FIELD **) malloc(nlevs*sizeof(FIELD*));
  reached_eof=0;
  while ( 1 )       
   {     
     nrecs = streamInqTimestep(streamID1, eofID);
     if ( nrecs == 0)
       {
         reached_eof = 1;
         break;
       }
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
  cdoPrint("%s contains %i eof's", cdoStreamName(1), neof);
  // Create 1x1 Grid for output
  gridID3 = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID3, 1);
  gridDefYsize(gridID3, 1);
  xvals=(double*)malloc(1*sizeof(double));
  yvals=(double*)malloc(1*sizeof(double));
  xvals[0]=0;
  yvals[0]=0;
  gridDefXvals(gridID3, xvals);
  gridDefYvals(gridID3, yvals);
  
  // Create var-list and time-axis for output
  vlistID3 = vlistDuplicate(vlistID2); 
  taxisID3 = taxisDuplicate(taxisID2);  
  
  // open streams for eofcoeff output
  streamIDs = (int *) malloc (neof*sizeof(int)); 
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
      if ( streamIDs[eofID] < 0) 
        cdiError(streamIDs[eofID], "Open failed on %s", oname);
      else if (cdoVerbose) 
        cdoPrint("opened %s ('w')  as stream%i for %i. eof", oname, streamIDs[eofID], eofID+1);
      
      streamDefVlist(streamIDs[eofID], vlistID3);
      vlistDefTaxis(vlistID3, taxisID3);  
                   
      ngrids = vlistNgrids(vlistID3);
      if (cdoVerbose)
        cdoPrint("ngrids for %i.eof: %i", eofID+1, ngrids);
      for ( i = 0; i < ngrids; i++ )
        vlistChangeGridIndex(vlistID3, i, gridID3);            
    }
  
  // ALLOCATE temporary fields for data read and write
  in.ptr = (double*) malloc(gridsize*sizeof(double));
  in.grid = gridID1;  
  out.missval = missval1;
  out.nmiss = 0;
  out.ptr = (double *) malloc (1*sizeof(double));
 
  reached_eof=0;
  tsID=0;
  while ( 1 )
    {      
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
                  if ( in.ptr[i] != missval2 && 
                      eof[varID][levelID][eofID].ptr[i] != missval1 )
                    {
                      double tmp = w[i]*in.ptr[i]*eof[varID][levelID][eofID].ptr[i];
                      out.ptr[0] += tmp;                   
                    }
                }
              //printf("grid %i\n", out.grid);
              //printf("stream%i ts%i var%i level%i eof%i %17.15f\n", streamIDs[eofID], tsID, varID, levelID, eofID, out.ptr[0]); 
              if ( out.ptr[0] ) nmiss=0;
              else { nmiss=1; out.ptr[0]=missval2; }
              
              //printf("writing nmiss=%i streamID %i \n",nmiss, streamIDs[eofID]);
              vdate = taxisInqVdate(taxisID2);
              vtime = taxisInqVtime(taxisID2);
              taxisDefVdate(taxisID3, vdate);
              taxisDefVtime(taxisID3, vtime);
              
              streamDefTimestep(streamIDs[eofID],tsID);
              streamDefRecord(streamIDs[eofID], varID, levelID);
              
              // TEST FOR GRID
              int tmp = vlistInqVarGrid(vlistID3, varID);
              //printf("grid %i\n", tmp);
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


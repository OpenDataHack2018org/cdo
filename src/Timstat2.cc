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

        Timstat2        timcor      correlates two data files on the same grid
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

// correlation in time
static
void correlationInit(size_t gridsize, const double *array1, const double *array2,
                     double missval1, double missval2, size_t *nofvals, 
                     double *work0, double *work1, double *work2, double *work3, double *work4)
{
  for ( size_t i = 0; i < gridsize; ++i )
    {
      if ( ( ! DBL_IS_EQUAL(array1[i], missval1) ) && 
           ( ! DBL_IS_EQUAL(array2[i], missval2) ) )
        {
          work0[i] += array1[i];
          work1[i] += array2[i];
          work2[i] += array1[i]*array1[i];
          work3[i] += array2[i]*array2[i];
          work4[i] += array1[i]*array2[i];
          nofvals[i]++;
        }
    }	 
}

static
size_t correlation(size_t gridsize, double missval1, double missval2, size_t *nofvals, 
                   double *work0, double *work1, double *work2, double *work3, double *work4)
{
  size_t nmiss = 0;
  double cor;

  for ( size_t i = 0; i < gridsize; ++i )
    {	  
      size_t nvals = nofvals[i];
      if ( nvals > 0 )
	{
	  double temp0 = MULMN(work0[i], work1[i]);
	  double temp1 = SUBMN(work4[i], DIVMN(temp0, nvals));
	  double temp2 = MULMN(work0[i], work0[i]);
	  double temp3 = MULMN(work1[i], work1[i]);
	  double temp4 = SUBMN(work2[i], DIVMN(temp2, nvals));
	  double temp5 = SUBMN(work3[i], DIVMN(temp3, nvals));
	  double temp6 = MULMN(temp4, temp5);

	  cor = DIVMN(temp1, SQRTMN(temp6));
          if      ( cor < -1 )  cor = -1;
          else if ( cor >  1 )  cor =  1;

	  if ( DBL_IS_EQUAL(cor, missval1) ) nmiss++;
	}
      else
	{
	  nmiss++;
	  cor = missval1;
	}

      work0[i] = cor;
    }

  return nmiss;
}

// covariance in time
static
void covarianceInit(size_t gridsize, const double *array1, const double *array2,
                    double missval1, double missval2, size_t *nofvals, 
                    double *work0, double *work1, double *work2)
{
  for ( size_t i = 0; i < gridsize; ++i )
    {
      if ( ( ! DBL_IS_EQUAL(array1[i], missval1) ) && 
           ( ! DBL_IS_EQUAL(array2[i], missval2) ) )
        {
          work0[i] += array1[i];
          work1[i] += array2[i];
          work2[i] += array1[i]*array2[i];
          nofvals[i]++;
        }
    }	 
}

static
size_t covariance(size_t gridsize, double missval1, double missval2, size_t *nofvals, 
                  double *work0, double *work1, double *work2)
{
  size_t nmiss = 0;
  double covar;

  for ( size_t i = 0; i < gridsize; ++i )
    {	  
      size_t nvals = nofvals[i];
      if ( nvals > 0 )
	{
          double dnvals = nvals;
	  double temp = DIVMN( MULMN(work0[i], work1[i]), dnvals*dnvals);
	  covar = SUBMN( DIVMN(work2[i], dnvals), temp);
	  if ( DBL_IS_EQUAL(covar, missval1) ) nmiss++;
	}
      else
	{
	  nmiss++;
	  covar = missval1;
	}

      work0[i] = covar;
    }

  return nmiss;
}


void *Timstat2(void *process)
{
  int vdate = 0, vtime = 0;
  int varID, levelID;
  size_t nmiss;

  cdoInitialize(process);

  // clang-format off
  cdoOperatorAdd("timcor",   func_cor,   0, NULL);
  cdoOperatorAdd("timcovar", func_covar, 0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();
  int operfunc   = cdoOperatorF1(operatorID);

  int nwork = 0;
  if      ( operfunc == func_cor   ) nwork = 5;
  else if ( operfunc == func_covar ) nwork = 3;

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));
  int streamID2 = cdoStreamOpenRead(cdoStreamName(1));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = cdoStreamInqVlist(streamID2);
  int vlistID3 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, CMP_ALL);
 
  int nvars  = vlistNvars(vlistID1);
  int nrecs  = vlistNrecs(vlistID1);
  int nrecs3 = nrecs;
  std::vector<int> recVarID(nrecs);
  std::vector<int> recLevelID(nrecs);

  int taxisID1 = vlistInqTaxis(vlistID1);
  //int taxisID2 = vlistInqTaxis(vlistID2);
  int taxisID3 = taxisDuplicate(taxisID1);
 
  vlistDefTaxis(vlistID3, taxisID3);
  int streamID3 = cdoStreamOpenWrite(cdoStreamName(2), cdoFiletype());
  pstreamDefVlist(streamID3, vlistID3);
 
  size_t gridsizemax = vlistGridsizeMax(vlistID1);
  std::vector<double> array1(gridsizemax);
  std::vector<double> array2(gridsizemax);
  				 
  std::vector<std::vector<std::vector<std::vector<double>>>> work(nvars);
  std::vector<std::vector<std::vector<size_t>>> nofvals(nvars);

  for ( varID = 0; varID < nvars; varID++ )
    {
      int gridID = vlistInqVarGrid(vlistID1, 0);  
      size_t gridsize = gridInqSize(gridID);
      int nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));

      work[varID].resize(nlevs);
      nofvals[varID].resize(nlevs);  

      for ( levelID = 0; levelID < nlevs; levelID++ )
	{
          nofvals[varID][levelID].resize(gridsize, 0);
	  work[varID][levelID].resize(nwork);
	  for ( int i = 0; i < nwork; i++ )
            work[varID][levelID][i].resize(gridsize);
	}
    }
 
  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      int nrecs2 = pstreamInqTimestep(streamID2, tsID);
      if ( nrecs != nrecs2 )
        cdoWarning("Input streams have different number of records!");

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  pstreamInqRecord(streamID2, &varID, &levelID);

	  if ( tsID == 0 )
	    {
	      recVarID[recID]   = varID;
	      recLevelID[recID] = levelID;	     	     
	    }	 

	  size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

	  double missval1 = vlistInqVarMissval(vlistID1, varID);
	  double missval2 = vlistInqVarMissval(vlistID2, varID);

	  pstreamReadRecord(streamID1, &array1[0], &nmiss);
	  pstreamReadRecord(streamID2, &array2[0], &nmiss);

	  if ( operfunc == func_cor )
	    {
              correlationInit(gridsize, &array1[0], &array2[0], missval1, missval2, &nofvals[varID][levelID][0],
                              &work[varID][levelID][0][0], &work[varID][levelID][1][0],
                              &work[varID][levelID][2][0], &work[varID][levelID][3][0], 
                              &work[varID][levelID][4][0]);
	    }
	  else if ( operfunc == func_covar )
	    {
              covarianceInit(gridsize, &array1[0], &array2[0], missval1, missval2, &nofvals[varID][levelID][0],
                             &work[varID][levelID][0][0], &work[varID][levelID][1][0],
                             &work[varID][levelID][2][0]);
	    }
	}

      tsID++;
    }

  tsID = 0;
  taxisDefVdate(taxisID3, vdate);
  taxisDefVtime(taxisID3, vtime);
  pstreamDefTimestep(streamID3, tsID);

  for ( int recID = 0; recID < nrecs3; recID++ )
    {
      varID   = recVarID[recID];
      levelID = recLevelID[recID];
   
      size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

      double missval1 = vlistInqVarMissval(vlistID1, varID);
      double missval2 = vlistInqVarMissval(vlistID2, varID);

      if ( operfunc == func_cor )
	{
	  nmiss = correlation(gridsize, missval1, missval2, &nofvals[varID][levelID][0],
                              &work[varID][levelID][0][0], &work[varID][levelID][1][0],
                              &work[varID][levelID][2][0], &work[varID][levelID][3][0], 
                              &work[varID][levelID][4][0]);
	}
      else if ( operfunc == func_covar )
	{
	  nmiss = covariance(gridsize, missval1, missval2, &nofvals[varID][levelID][0],
                             &work[varID][levelID][0][0], &work[varID][levelID][1][0],
                             &work[varID][levelID][2][0]);
	}

      pstreamDefRecord(streamID3, varID, levelID);
      pstreamWriteRecord(streamID3, &work[varID][levelID][0][0], nmiss);
    }

  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);
    
  cdoFinish();   
 
  return 0;
}

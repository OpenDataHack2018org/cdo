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

      Intlevel   intlevel3d      Linear level interpolation on a 3d vertical coordinates variable
      Intlevel   intlevelx3d     Linear level interpolation on a 3d vertical coordinates variable with extrapolation
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "listarray.h"
#include "after_vertint.h"


bool levelDirUp(int nlev, double *lev);
bool levelDirDown(int nlev, double *lev);
double vert_interp_lev_kernel(double w1, double w2, double var1L1, double var1L2, double missval);
void vert_gen_weights(int expol, int nlev1, double *lev1, int nlev2, double *lev2,
		      int *lev_idx1, int *lev_idx2, double *lev_wgt1, double *lev_wgt2);

/*
 * 3d vertical interpolation routine (see vert_interp_lev() in src/Intlevel.cc)
 */
static
void vert_interp_lev3d(size_t gridsize, double missval, double *vardata1, double *vardata2,
		       int nlev2, int *lev_idx1, int *lev_idx2, double *lev_wgt1, double *lev_wgt2)
{
  for ( int ilev = 0; ilev < nlev2; ilev++ )
    {
      size_t offset = ilev*gridsize;
      double *var2 = vardata2 + offset;

      for ( size_t i = 0; i < gridsize; i++ )
	{
          int idx1 = lev_idx1[offset+i];
          int idx2 = lev_idx2[offset+i];
          double wgt1 = lev_wgt1[offset+i];
          double wgt2 = lev_wgt2[offset+i];

          /* upper/lower values from input field */
          double var1L1 = *(vardata1+idx1);
          double var1L2 = *(vardata1+idx2);

          /* if (cdoVerbose) printf("i:%d level %d: idx1=%d idx2=%d (offset+i:%d) wgt1=%g wgt2=%g var1L1:%g var1L2:%g ",
           *                         i,       ilev, idx1,   idx2,    offset+i,    wgt1,   wgt2,   var1L1,   var1L2);
           */
          var2[i] = vert_interp_lev_kernel(wgt1, wgt2, var1L1, var1L2, missval);
	}
    }
}

/*
 * Create weights for the 3d vertical coordinate
 *
 * The resulting index sets lev_idx1 and lev_idx2 contain absolute numbers,i.e.
 * wrt. the given gridsize. They can directly be used to read values from 3d
 * data fields.
 *
 * 3d version of vert_gen_weights() (src/Intlevel.cc)
 */
static
void vert_gen_weights3d(bool expol, int nlev1, size_t gridsize, double *xlev1, int nlev2, double *xlev2,
			int *xlev_idx1, int *xlev_idx2, double *xlev_wgt1, double *xlev_wgt2)
{
  std::vector<double> lev1(nlev1);
  std::vector<double> lev2(nlev2);
  std::vector<int> lev_idx1(nlev2);
  std::vector<int> lev_idx2(nlev2);
  std::vector<double> lev_wgt1(nlev2);
  std::vector<double> lev_wgt2(nlev2);

  for ( size_t i = 0; i < gridsize; i++ )
    {
      for ( int k = 0; k < nlev1; ++k ) lev1[k] = xlev1[k*gridsize+i];
      for ( int k = 0; k < nlev2; ++k ) lev2[k] = xlev2[k*gridsize+i];

      vert_gen_weights(expol, nlev1, &lev1[0], nlev2, &lev2[0], &lev_idx1[0], &lev_idx2[0], &lev_wgt1[0], &lev_wgt2[0]);

      for ( int k = 0; k < nlev2; ++k ) xlev_idx1[k*gridsize+i] = lev_idx1[k]*gridsize+i;
      for ( int k = 0; k < nlev2; ++k ) xlev_idx2[k*gridsize+i] = lev_idx2[k]*gridsize+i;
      for ( int k = 0; k < nlev2; ++k ) xlev_wgt1[k*gridsize+i] = lev_wgt1[k];
      for ( int k = 0; k < nlev2; ++k ) xlev_wgt2[k*gridsize+i] = lev_wgt2[k];
    }
}


void *Intlevel3d(void *argument)
{
  size_t gridsize, gridSize, gridsizei, gridsizeo;
  int nrecs;
  int i, offset;
  int varID, levelID;
  int nvars,nvct;
  size_t nmiss;
  int zaxisID1 = -1, zaxisID3;
  int gridID3 = -1, gridID, zaxisID;
  int nlevi, nlevo, nlevel = 0, maxlev;
  double missval;
  double *lev2 = NULL;
  double *single1, *single2;
  int taxisID1, taxisID3;
  double *zlevels_in, *zlevels_out;
  size_t zlevels_in_miss, zlevels_out_miss;
  char varname[10]; 

  cdoInitialize(argument);

  // clang-format off
  int INTLEVEL3D  = cdoOperatorAdd("intlevel3d",  0, 0, NULL);
  int INTLEVELX3D = cdoOperatorAdd("intlevelx3d", 0, 0, NULL);
  // clang-format on

  UNUSED(INTLEVEL3D);

  int operatorID = cdoOperatorID();

  bool expol = (operatorID == INTLEVELX3D);

  operatorInputArg("icoordinate");

  int streamID1 = pstreamOpenRead(cdoStreamName(0));                 /*  input data */
  int streamID2 = pstreamOpenRead(cdoStreamName(1));                 /*  3d target vertical coordinate */
  int streamID3 = pstreamOpenWrite(cdoStreamName(2),cdoFiletype());  /*  output stream */

  /*  Read filename from Parameter */
  operatorInputArg("filename for vertical source coordinates variable");
  operatorCheckArgc(1);
  argument_t *fileargument = file_argument_new(operatorArgv()[0]);
  int streamID0 = pstreamOpenRead(fileargument);                     /*  3d vertical input coordinate */
  file_argument_free(fileargument);

  int vlistID0 = pstreamInqVlist(streamID0);
  int vlistID1 = pstreamInqVlist(streamID1); taxisID1 = vlistInqTaxis(vlistID1);
  int vlistID2 = pstreamInqVlist(streamID2);
  int vlistID3 = vlistDuplicate(vlistID1);  taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID3, taxisID1);

  /*
   * Read 3d input coordinate (streamID0)
   * * two additional levels are added (top + bottom) for later extrapolation checking
   */
  {
    int nvars = vlistNvars(vlistID0);
    if ( nvars != 1 ) cdoAbort("Only one single variable is allowed!");

    gridID     = vlistInqVarGrid(vlistID0, 0);
    zaxisID    = vlistInqVarZaxis(vlistID0, 0);
    gridsize   = gridInqSize(gridID);
    nlevel     = zaxisInqSize(zaxisID);

    zlevels_in = (double*) Malloc(gridsize*(nlevel+2)*sizeof(double));
    nlevi      = nlevel;   /* number of input levels for later use */
    gridsizei  = gridsize; /* horizontal gridsize of input z coordinate */
    nrecs      = pstreamInqTimestep(streamID0, 0);
    if (cdoVerbose) cdoPrint("%d records input 3d vertical height",nrecs);

    for ( int recID = 0; recID < nrecs; recID++ )
      {
        pstreamInqRecord(streamID0, &varID, &levelID);
        gridsize = gridInqSize(vlistInqVarGrid(vlistID0, varID));
        offset   = gridsize + gridsize*levelID;
        single1  = zlevels_in + offset;
        pstreamReadRecord(streamID0, single1, &zlevels_in_miss);
      }
  }

  /*
   * Read 3d output coordinate (streamID2)
   */
  {
    int nvars = vlistNvars(vlistID2);
    if (nvars != 1) cdoAbort("Only one single variable is allowed!");
    gridID      = vlistInqVarGrid(vlistID2, varID);
    gridID3     = gridID;
    zaxisID     = vlistInqVarZaxis(vlistID2, varID);
    gridsize    = gridInqSize(gridID);
    nlevel      = zaxisInqSize(zaxisID);

    zlevels_out = (double*) Malloc(gridsize*nlevel*sizeof(double));
    nlevo       = nlevel;  /* number of output levels for later use */
    gridsizeo   = gridsize;/* horizontal gridsize of output z coordinate */
    nrecs       = pstreamInqTimestep(streamID2, 0);
    if (cdoVerbose) cdoPrint("%d records target 3d vertical height and gridsize %d",nrecs,gridsize);

    for ( int recID = 0; recID < nrecs; recID++ )
      {
	pstreamInqRecord(streamID2, &varID, &levelID);
	gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	offset   = gridsize*levelID;
	single1  = zlevels_out + offset;
	pstreamReadRecord(streamID2, single1, &zlevels_out_miss);
      }
  }

  /* Missing values are not allowed for coordinate variables */
  if ( 0 != zlevels_in_miss  )
    cdoAbort("Input vertical coordinate variables are not allowed to contain missing values.");
  else
    {
      if ( cdoVerbose ) cdoPrint("Input vertical coordinate has no missing values.");
    }

  if ( 0 != zlevels_out_miss )
    cdoAbort("Output vertical coordinate variables are not allowd to contain missing values.");
  else
    {
      if ( cdoVerbose ) cdoPrint("Output vertical coordinate has no missing values.");
    }


  /*
   * gridsize of input and output vertical coordinate must be equal
   * (later use of gridsizeo ONLY)
   */
  if ( gridsizei != gridsizeo )
    cdoAbort("Input and output vertical coordinate must have the same gridsize!");

  gridSize = gridsizeo;

   /* input and output vertical coordinates must have exactly the same horizontal grid */
   if ( gridsizei != gridsizeo  )
     {
       /* i =0; printf ( "lonIn:%g latIn:%g lonOut:%g latOut:%g\n",lonIn[i],latIn[i],lonOut[i],latOut[i] ); */
       cdoAbort("Input and output vertical coordinates do NOT exactly have the same horizontal grid.");
     }

  /*
   * Check for the correct vertical axis in the input: Variables with the same
   * number of levels as the input vertical levels from operators parameter
   * (streamID0). Variables with a different z-axis should be copied into output.
   */
  int nzaxis = vlistNzaxis(vlistID1);
  for ( i = 0; i < nzaxis; ++i )
    {
      zaxisID = vlistZaxis(vlistID1, i);
      nlevel  = zaxisInqSize(zaxisID);
      if ( nlevel == nlevi )
        {
          zaxisID1 = zaxisID;
          break;
        }
    }
  if ( i == nzaxis ) cdoAbort("No processable variable found (vertical coordinate differ)!");

  int ngrids = vlistNgrids(vlistID1);
  for ( i = 0; i < ngrids; ++i )
    {
      gridID = vlistGrid(vlistID1, i);
      gridsize = gridInqSize(gridID);
      if ( gridsize == gridSize ) break;
    }
  if ( i == nzaxis ) cdoAbort("No processable variable found (grid coordinate differ)!");
  
  /*
   * Check monotony of vertical levels
   */
  std::vector<double> lev1(nlevi);
  for ( int i = 0; i < nlevi; ++i ) lev1[i] = zlevels_in[(i+1)*gridSize]; 
  bool lup = levelDirUp(nlevi, &lev1[0]);
  bool ldown = levelDirDown(nlevi, &lev1[0]);

  /* Add artificial values for intication of extrapolation areas (lowermost + upmost levels) */
  if ( lup )
    {
      for ( size_t i = 0; i < gridSize ;i++ )
      {
        zlevels_in[i]                      = -1.e33;
        zlevels_in[(nlevi+1)*gridSize + i] =  1.e33;
      }
    }
  else if ( ldown )
    {
      for ( size_t i = 0; i < gridSize ;i++ )
      {
        zlevels_in[i]                      =  1.e33;
        zlevels_in[(nlevi+1)*gridSize + i] = -1.e33;
      }
    }
  else
    cdoWarning("Non monotonic zaxis!");

  /*
   * Create weights for later interpolation - assumption: input vertical correct is constant in time
   */
  int *lev_idx1 = (int*) Malloc(nlevo*gridSize*sizeof(int));
  int *lev_idx2 = (int*) Malloc(nlevo*gridSize*sizeof(int));
  double *lev_wgt1 = (double*) Malloc(nlevo*gridSize*sizeof(double));
  double *lev_wgt2 = (double*) Malloc(nlevo*gridSize*sizeof(double));

  vert_gen_weights3d(expol, nlevi+2, gridSize, zlevels_in, nlevo, zlevels_out, lev_idx1, lev_idx2, lev_wgt1, lev_wgt2);

  /*
   * Copy z-axis information to output z-axis
   */
  zaxisID3 = zaxisCreate(zaxisInqType(zaxisID1), nlevo);
  lev2 = (double*) Malloc(nlevo*sizeof(double));
  /* fill values with its indices */
  for (i=0;i<nlevo;i++)
    lev2[i] = (double) i+1;
  zaxisDefLevels(zaxisID3, lev2);
  zaxisDefName(zaxisID3, "lev");
  /*  copy VCT from input vlistID1 to output vlistID3 if there is one */
  nvct = zaxisInqVctSize(zaxisID1);
  if ( nvct > 0 ) zaxisDefVct(zaxisID3,zaxisInqVctSize(zaxisID1), zaxisInqVctPtr(zaxisID1));

  for ( i = 0; i < nzaxis; i++ )
    if ( zaxisID1 == vlistZaxis(vlistID1, i) )
      vlistChangeZaxisIndex(vlistID3, i, zaxisID3);
  /* add the vertical output field to the output stream */
  int oz3dvarID = vlistDefVar(vlistID3, gridID3, zaxisID3, TIME_VARYING);
  {
    char str[256];
    str[0] = 0;
    vlistInqVarName(vlistID2,0,str);
    vlistDefVarName(vlistID3,oz3dvarID,str);
    str[0] = 0;
    vlistInqVarLongname(vlistID2,0,str);
    if ( str[0] ) vlistDefVarLongname(vlistID3,oz3dvarID, str);
    str[0] = 0;
    vlistInqVarUnits(vlistID2,0, str);
    if ( str[0] ) vlistDefVarUnits(vlistID3,oz3dvarID, str);
  }

  pstreamDefVlist(streamID3, vlistID3);

  maxlev = nlevi > nlevo ? nlevi : nlevo;
  nvars = vlistNvars(vlistID1);
  bool *vars = (bool*) Malloc(nvars*sizeof(bool));
  bool *varinterp = (bool*) Malloc(nvars*sizeof(bool));   /* marker for variables to be interpolated       */
  size_t **varnmiss = (size_t**) Malloc(nvars*sizeof(size_t*)); /* can for missing values of arbitrary variables */
  double **vardata1 = (double**) Malloc(nvars*sizeof(double*)); /* input                                         */
  double **vardata2 = (double**) Malloc(nvars*sizeof(double*)); /* output                                        */

  /* by default no variable should be interpolated */
  for ( i = 0; i < nvars; i++ ) varinterp[varID] = false;

  for ( varID = 0; varID < nvars; varID++ )
    {
      int gridID = vlistInqVarGrid(vlistID1, varID);
      int zaxisID = vlistInqVarZaxis(vlistID1, varID);
      size_t gridsize = gridInqSize(gridID);
      int nlevel = zaxisInqSize(zaxisID);

      vlistInqVarName(vlistID1, varID, varname);

      vardata1[varID] = (double*) Malloc(gridsize*nlevel*sizeof(double));

      /*  variabls for interpolation:
       *  * have the required vertical axis, i.e. the correct number of levels (nlevi)
       *  * have the same number of horizontal grid points (i.e. same gridSize) like the two vertical coordinates
       *  * are NOT the output vertical coordinates itself
       */
      if ( zaxisID == zaxisID1 && varID != oz3dvarID && gridsize == gridSize )
        {
          if ( gridsizeo != gridsize )
            {
              varinterp[varID] = false;
              vardata2[varID]  = vardata1[varID];
              varnmiss[varID]  = (size_t*) Malloc(nlevel*sizeof(size_t));
              if ( cdoVerbose ) cdoPrint("Ignore variable %s (levels=%d gridsize=%zu)!", varname, nlevel, gridsize);
            }
          else
            {
              varinterp[varID] = true;
              vardata2[varID]  = (double*) Malloc(gridsize*nlevo*sizeof(double));
              varnmiss[varID]  = (size_t*) Malloc(maxlev*sizeof(size_t));
              memset(varnmiss[varID], 0, maxlev*sizeof(size_t));
            }
        }
      else
        {
          varinterp[varID] = false;
          vardata2[varID]  = vardata1[varID];
          varnmiss[varID]  = (size_t*) Malloc(nlevel*sizeof(size_t));
          if ( cdoVerbose ) cdoPrint("Ignore variable %s (levels=%d gridsize=%zu)!", varname, nlevel, gridsize);
        }
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      if ( varinterp[varID] ) break;
    }
  if ( varID == nvars ) cdoAbort("No processable variable found!");

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      for ( varID = 0; varID < nvars; ++varID ) vars[varID] = false;

      taxisCopyTimestep(taxisID3, taxisID1);
      pstreamDefTimestep(streamID3, tsID);

      /*
       * Read the whole 3d data field
       */
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
          vlistInqVarName(vlistID1, varID, varname); 
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  offset   = gridsize*levelID;
	  single1  = vardata1[varID] + offset;
          pstreamReadRecord(streamID1, single1, &varnmiss[varID][levelID]);
	  vars[varID] = true;
	}

      /* Perform the interpolation on all valid data variables */
      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vars[varID] && varinterp[varID] )
	    {
	      gridID   = vlistInqVarGrid(vlistID1, varID);
	      missval  = vlistInqVarMissval(vlistID1, varID);
	      gridsize = gridInqSize(gridID);

	      vert_interp_lev3d(gridsize, missval, vardata1[varID], vardata2[varID],
				nlevo, lev_idx1, lev_idx2, lev_wgt1, lev_wgt2);

	      for ( levelID = 0; levelID < nlevo; levelID++ )
		{
		  gridsize = gridInqSize(vlistInqVarGrid(vlistID3, varID));
		  offset   = gridsize*levelID;
		  single2  = vardata2[varID] + offset;
		  nmiss    = 0;
		  for ( size_t i = 0; i < gridsize; ++i )
		    if ( DBL_IS_EQUAL(single2[i], missval) ) nmiss++;
		  varnmiss[varID][levelID] = nmiss;
		}
	    }
          else
          { 
            vlistInqVarName(vlistID1, varID, varname); 
            if ( cdoVerbose && tsID <= 1 ) cdoPrint("Perform no interpolation on variable %s", varname);
          }
	}

      /* write the output */
      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vars[varID] )
	    {
	      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID3, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  gridsize = gridInqSize(vlistInqVarGrid(vlistID3, varID));
		  offset   = gridsize*levelID;

		  single2  = vardata2[varID] + offset;
		  pstreamDefRecord(streamID3, varID, levelID);
		  pstreamWriteRecord(streamID3, single2, varnmiss[varID][levelID]);
		}
	    }
	}

      /* copy output z coordinate to output stream */
      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID3, oz3dvarID));
      for ( levelID = 0; levelID < nlevel; levelID++ )
        {
          gridsize = gridInqSize(vlistInqVarGrid(vlistID3, oz3dvarID));
          offset   = gridsize*levelID;
          single2  = zlevels_out + offset;
          pstreamDefRecord(streamID3, oz3dvarID, levelID);
          pstreamWriteRecord(streamID3, single2, 0);

        }
      tsID++;
    }

  nvars = vlistNvars(vlistID1);
  for ( varID = 0; varID < nvars; varID++ )
    {
      Free(varnmiss[varID]);
      Free(vardata1[varID]);
      if ( varinterp[varID] ) Free(vardata2[varID]);
    } 


  Free(varinterp);
  Free(varnmiss);
  Free(vardata2);
  Free(vardata1);
  Free(vars);

  Free(lev_idx1);
  Free(lev_idx2);
  Free(lev_wgt1);
  Free(lev_wgt2);

  pstreamClose(streamID0);
  pstreamClose(streamID1);
  pstreamClose(streamID2);
  pstreamClose(streamID3);

  vlistDestroy(vlistID3);

  cdoFinish();

  return 0;
}

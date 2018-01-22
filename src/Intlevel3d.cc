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


/*
 * 3d vertical interpolation routine (see vert_interp_lev() in src/Intlevel.cc)
 */
static
void vert_interp_lev3d(int gridsize, double missval, double *vardata1, double *vardata2,
		       int nlev2, int *lev_idx1, int *lev_idx2, double *lev_wgt1, double *lev_wgt2)
{
  int idx1, idx2;
  int offset;
  double wgt1, wgt2;
  double w1, w2;
  double var1L1, var1L2, *var2;

  for ( int ilev = 0; ilev < nlev2; ilev++ )
    {
      offset = ilev*gridsize;
      var2 = vardata2 + offset;

      for ( int i = 0; i < gridsize; i++ )
	{
          idx1 = lev_idx1[offset+i];
          idx2 = lev_idx2[offset+i];
          wgt1 = lev_wgt1[offset+i];
          wgt2 = lev_wgt2[offset+i];

          /* upper/lower values from input field */
          var1L1 = *(vardata1+idx1);
          var1L2 = *(vardata1+idx2);

          /* if (cdoVerbose) printf("i:%d level %d: idx1=%d idx2=%d (offset+i:%d) wgt1=%g wgt2=%g var1L1:%g var1L2:%g ",
           *                         i,       ilev, idx1,   idx2,    offset+i,    wgt1,   wgt2,   var1L1,   var1L2);
           */
	  w1 = wgt1;
	  w2 = wgt2;
	  if ( DBL_IS_EQUAL(var1L1, missval) ) w1 = 0;
	  if ( DBL_IS_EQUAL(var1L2, missval) ) w2 = 0;

	  if ( IS_EQUAL(w1, 0) && IS_EQUAL(w2, 0) )
	    {
	      var2[i] = missval;
	    }
	  else if ( IS_EQUAL(w1, 0) )
	    {
              var2[i] = (w2 >= 0.5) ? var1L2 : missval;	      
	    }
	  else if ( IS_EQUAL(w2, 0) )
	    {
              var2[i] = (w1 >= 0.5) ? var1L1 : missval;	      
	    }
	  else
	    {
	      var2[i] = var1L1*w1 + var1L2*w2;
	    }
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
void vert_gen_weights3d(bool expol, int nlev1, int gridsize, double *lev1, int nlev2, double *lev2,
			int *lev_idx1, int *lev_idx2, double *lev_wgt1, double *lev_wgt2)
{
  int i,i1, i2;
  int idx1 = 0, idx2 = 0;
  double val1, val2 = 0;

  for ( i = 0; i < gridsize; i++ )
    {
      for ( i2 = 0; i2 < nlev2; i2++ )
        {
          /* Because 2 levels were added to the source vertical coordinate (one on
           * top, one at the bottom), its loop starts at 1 */
          for ( i1 = 1; i1 < nlev1; i1++ )
            {
              if ( lev1[(i1-1)*gridsize+i] < lev1[i1*gridsize+i] )
                {
                  idx1 = (i1-1)*gridsize+i;
                  idx2 = i1*gridsize+i;
                }
              else
                {
                  idx1 = i1*gridsize+i;
                  idx2 = (i1-1)*gridsize+i;
                }
              val1 = lev1[idx1];
              val2 = lev1[idx2];

              if ( lev2[i2*gridsize+i] > val1 && lev2[i2*gridsize+i] <= val2 ) break;
            }

          if ( i1 == nlev1 ) 
            {
              if ( expol )
                cdoAbort("Level %g at index %d not found! Use extrapolation", lev2[i2*gridsize],i2);
              else
                cdoAbort("Level %g at index %d not found!");
            }

          if ( i1-1 == 0 ) /* destination levels ios not covert by the first two input z levels */
            {
              lev_idx1[i2*gridsize+i] = gridsize+i;
              lev_idx2[i2*gridsize+i] = gridsize+i;
              lev_wgt1[i2*gridsize+i] = 0;
              lev_wgt2[i2*gridsize+i] = (expol || IS_EQUAL(lev2[i2*gridsize+i], val2));
            }
          else if ( i1 == nlev1-1 ) /* destination level is beyond the last value of the input z field */
            {
              lev_idx1[i2*gridsize+i] = (nlev1-2)*gridsize+i;
              lev_idx2[i2*gridsize+i] = (nlev1-2)*gridsize+i;
              lev_wgt1[i2*gridsize+i] = (expol || IS_EQUAL(lev2[i2*gridsize+i], val2));
              lev_wgt2[i2*gridsize+i] = 0;
            }
          else /* target z values has two bounday values in input z field */
            {
              lev_idx1[i2*gridsize+i] = idx1;
              lev_idx2[i2*gridsize+i] = idx2;
              lev_wgt1[i2*gridsize+i] = (lev1[idx2]        - lev2[i2*gridsize+i]) / (lev1[idx2] - lev1[idx1]);
              lev_wgt2[i2*gridsize+i] = (lev2[i2*gridsize+i] - lev1[idx1])        / (lev1[idx2] - lev1[idx1]);

            }
  /*         if (cdoVerbose)
   *         {
   *           printf("i:%d i2:%d\ti2*gridsize+i:%d\tlev2[i2*gridsize+i]:%g\tidx1:%d\tidx2:%d\tlev1[idx1]:%g\tlev1[idx2]:%g\t",
   *                   i, i2, i2*gridsize+i,         lev2[i2*gridsize+i],    idx1,    idx2,    lev1[idx1],    lev1[idx2]);
   *           printf("\tlev_wgt1:%g\tlev_wgt2:%g\n", lev_wgt1[i2*gridsize+i], lev_wgt2[i2*gridsize+i]);
   *         }
   */
          /* backshift of the indices because of the two additional levels in input vertical coordinate */
          lev_idx1[i2*gridsize+i] -= gridsize;
          lev_idx2[i2*gridsize+i] -= gridsize;

        }
    }
}


/*
 * Create weights for the 1d vertical coordinate from a 3d vertical coordinate
 *
 * The resulting index sets lev_idx1 and lev_idx2 contain absolute numbers,i.e.
 * wrt. the given gridsize. They can directly be used to read values from 3d
 * data fields.
 *
 * 3d1d version of vert_gen_weights() (src/Intlevel.cc)
 */
static
void vert_gen_weights3d1d(bool expol, int nlev1, int gridsize, double *lev1, int nlev2, double *lev2,
			  int *lev_idx1, int *lev_idx2, double *lev_wgt1, double *lev_wgt2)
{
  int i,i1, i2;
  int idx1 = 0, idx2 = 0;
  double val1, val2 = 0;

  for ( i = 0; i < gridsize; i++ )
    {
      for ( i2 = 0; i2 < nlev2; i2++ )
        {
          /* Because 2 levels were added to the source vertical coordinate (one on
           * top, one at the bottom), its loop starts at 1 */
          for ( i1 = 1; i1 < nlev1; i1++ )
            {
              if ( lev1[(i1-1)*gridsize+i] < lev1[i1*gridsize+i] )
                {
                  idx1 = (i1-1)*gridsize+i;
                  idx2 = i1*gridsize+i;
                }
              else
                {
                  idx1 = i1*gridsize+i;
                  idx2 = (i1-1)*gridsize+i;
                }
              val1 = lev1[idx1];
              val2 = lev1[idx2];

              if ( lev2[i2] > val1 && lev2[i2] <= val2 ) break;
            }

          if ( i1 == nlev1 ) 
            {
              if ( expol )
                cdoAbort("Level %g at index %d not found! Use extrapolation", lev2[i2],i2);
              else
                cdoAbort("Level %g at index %d not found!");
            }

          if ( i1-1 == 0 ) /* destination levels ios not covert by the first two input z levels */
            {
              lev_idx1[i2*gridsize+i] = gridsize+i;
              lev_idx2[i2*gridsize+i] = gridsize+i;
              lev_wgt1[i2*gridsize+i] = 0;
              lev_wgt2[i2*gridsize+i] = (expol || IS_EQUAL(lev2[i2], val2));
            }
          else if ( i1 == nlev1-1 ) /* destination level is beyond the last value of the input z field */
            {
              lev_idx1[i2*gridsize+i] = (nlev1-2)*gridsize+i;
              lev_idx2[i2*gridsize+i] = (nlev1-2)*gridsize+i;
              lev_wgt1[i2*gridsize+i] = (expol || IS_EQUAL(lev2[i2], val2));
              lev_wgt2[i2*gridsize+i] = 0;
            }
          else /* target z values has two bounday values in input z field */
            {
              lev_idx1[i2*gridsize+i] = idx1;
              lev_idx2[i2*gridsize+i] = idx2;
              lev_wgt1[i2*gridsize+i] = (lev1[idx2]  - lev2[i2]) / (lev1[idx2] - lev1[idx1]);
              lev_wgt2[i2*gridsize+i] = (lev2[i2] - lev1[idx1])  / (lev1[idx2] - lev1[idx1]);

            }
  /*         if (cdoVerbose)
   *         {
   *           printf("i:%d i2:%d\ti2*gridsize+i:%d\tlev2[i2]:%g\tidx1:%d\tidx2:%d\tlev1[idx1]:%g\tlev1[idx2]:%g\t",
   *                   i, i2, i2*gridsize+i,         lev2[i2],    idx1,    idx2,    lev1[idx1],    lev1[idx2]);
   *           printf("\tlev_wgt1:%g\tlev_wgt2:%g\n", lev_wgt1[i2*gridsize+i], lev_wgt2[i2*gridsize+i]);
   *         }
   */
          /* backshift of the indices because of the two additional levels in input vertical coordinate */
          lev_idx1[i2*gridsize+i] -= gridsize;
          lev_idx2[i2*gridsize+i] -= gridsize;
        }
    }
}


void *Intlevel3d(void *argument)
{
  int gridsize, gridSize, gridsizei, gridsizeo;
  int nrecs;
  int i, offset;
  int tsID, varID, levelID;
  int nvars,nvct;
  int nzaxis;
  size_t nmiss;

  int nlonIn, nlatIn, nlonOut, nlatOut;
  //double *lonIn, *latIn, *lonOut, *latOut;

  int zaxisID1 = -1, zaxisID3;
  int gridID3 = -1, gridID, zaxisID;
  int nlevi, nlevo, nlevel = 0, maxlev;
  double missval;
  double *lev1 = NULL, *lev2 = NULL;
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
    nvars = vlistNvars(vlistID0);
    if ( nvars != 1 ) cdoAbort("Only one single variable is allowed!");

    gridID     = vlistInqVarGrid(vlistID0, 0);
    zaxisID    = vlistInqVarZaxis(vlistID0, 0);
    gridsize   = gridInqSize(gridID);
    nlevel     = zaxisInqSize(zaxisID);

    nlonIn  = gridInqXsize(gridID);
    nlatIn  = gridInqYsize(gridID);
    /*
    lonIn   = (double*) Malloc(nlonIn*sizeof(double));
    latIn   = (double*) Malloc(nlatIn*sizeof(double));
    gridInqXvals(gridID, lonIn);
    gridInqYvals(gridID, latIn);
    */
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
    nvars = vlistNvars(vlistID2);
    if (nvars != 1) cdoAbort("Only one single variable is allowed!");
    gridID      = vlistInqVarGrid(vlistID2, varID);
    gridID3     = gridID;
    zaxisID     = vlistInqVarZaxis(vlistID2, varID);
    gridsize    = gridInqSize(gridID);
    nlevel      = zaxisInqSize(zaxisID);

    nlonOut = gridInqXsize(gridID);
    nlatOut = gridInqYsize(gridID);
    /*
    lonOut  = (double*) Malloc(nlonOut*sizeof(double));
    latOut  = (double*) Malloc(nlatOut*sizeof(double));
    gridInqXvals(gridID, lonOut);
    gridInqYvals(gridID, latOut);
    */
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
   if ( nlonIn != nlonOut || 
        nlatIn != nlatOut /*||
        memcmp(lonIn,lonOut,nlonIn*sizeof(double)) ||
        memcmp(latIn,latOut,nlatIn*sizeof(double))*/ )
     {
       /* i =0; printf ( "lonIn:%g latIn:%g lonOut:%g latOut:%g\n",lonIn[i],latIn[i],lonOut[i],latOut[i] ); */
       cdoAbort("Input and output vertical coordinates do NOT exactly have the same horizontal grid.");
     }

  /*
   * Check for the correct vertical axis in the input: Variables with the same
   * number of levels as the input vertical levels from operators parameter
   * (streamID0). Variables with a different z-axis should be copied into output.
   */
  nzaxis = vlistNzaxis(vlistID1);
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
  bool lup = false;
  bool ldown = false;
  lev1 = zlevels_in + gridSize;
  for ( i = 0; i < nlevi-1; i++ )
    {
      if ( i == 1 )
	{
	  if ( lev1[(i+1)*gridSize] > lev1[i*gridSize] )
	    lup = true;
	  else if ( lev1[(i+1)*gridSize] < lev1[i*gridSize] )
	    ldown = true;	
	}
      else
	{
	  if ( lup )
	    {
	      if ( !(lev1[(i+1)*gridSize] > lev1[i*gridSize]) ) lup = false;
	    }
	  else if ( ldown )
	    {
	      if ( !(lev1[(i+1)*gridSize] < lev1[i*gridSize]) ) ldown = false;
	    }
	}
    }

  /* Add artificial values for intication of extrapolation areas (lowermost + upmost levels) */
  if ( lup )
    {
      for ( i = 0; i < gridSize ;i++)
      {
        zlevels_in[i]                      = -1.e33;
        zlevels_in[(nlevi+1)*gridSize + i] =  1.e33;
      }
    }
  else if ( ldown )
    {
      for ( i = 0; i < gridSize ;i++)
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
    lev2[i] = (double) i;
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

  maxlev    = nlevi > nlevo ? nlevi : nlevo;
  nvars     = vlistNvars(vlistID1);
  bool *vars = (bool*) Malloc(nvars*sizeof(bool));
  bool *varinterp = (bool*) Malloc(nvars*sizeof(bool));   /* marker for variables to be interpolated       */
  size_t **varnmiss = (size_t**) Malloc(nvars*sizeof(size_t*)); /* can for missing values of arbitrary variables */
  double **vardata1 = (double**) Malloc(nvars*sizeof(double*)); /* input                                         */
  double **vardata2 = (double**) Malloc(nvars*sizeof(double*)); /* output                                        */

  /* by default no variable should be interpolated */
  for ( i = 0; i < nvars; i++ ) varinterp[varID] = false;

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      zaxisID  = vlistInqVarZaxis(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(zaxisID);

      vlistInqVarName(vlistID1, varID, varname);

      vardata1[varID] = (double*) Malloc(gridsize*nlevel*sizeof(double));

      /*  variabls for interpolation:
       *  * have the required vertical axis, i.e. the correct number of levels (nlevi)
       *  * have the same number of horizontal grid points (i.e. same gridSize) like the two vertical coordinates
       *  * are NOT the output vertical coordinates itself
       */
      if ( zaxisID == zaxisID1 && varID != oz3dvarID && gridsize == gridSize )
        {
          nlonIn  = gridInqXsize(gridID);
          nlatIn  = gridInqYsize(gridID);
	  /*
          lonIn   = (double*) Malloc(nlonIn*sizeof(double));
          latIn   = (double*) Malloc(nlatIn*sizeof(double));
          gridInqXvals(gridID, lonIn);
          gridInqYvals(gridID, latIn);
	  */
          if ( nlonIn != nlonOut || 
               nlatIn != nlatOut /*||
               memcmp(lonIn,lonOut,nlonIn*sizeof(double)) ||
               memcmp(latIn,latOut,nlatIn*sizeof(double))*/ )
            {
              varinterp[varID] = false;
              vardata2[varID]  = vardata1[varID];
              varnmiss[varID]  = (size_t*) Malloc(nlevel*sizeof(size_t));
              if ( cdoVerbose ) cdoPrint("Ignore variable %s with %d levels",varname,nlevel);
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
          if ( cdoVerbose ) cdoPrint("Ignore variable %s with %d levels",varname,nlevel);
        }
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      if ( varinterp[varID] ) break;
    }
  if ( varID == nvars ) cdoAbort("No processable variable found!");

  tsID = 0;
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
		  for ( i = 0; i < gridsize; ++i )
		    if ( DBL_IS_EQUAL(single2[i], missval) ) nmiss++;
		  varnmiss[varID][levelID] = nmiss;
		}
	    }
          else
          { 
            vlistInqVarName(vlistID1, varID, varname); 
            if ( cdoVerbose ) cdoPrint("Perform no interpolation on variable %s",varname);
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

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

      Intlevel   intlevel        Linear level interpolation
      Intlevel   intlevel3d      Linear level interpolation on a 3d vertical coordinates variable
*/


#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "listarray.h"


static
double vert_interp_lev_kernel(double w1, double w2, double var1L1, double var1L2, double missval)
{
  double var2 = missval;

  if ( DBL_IS_EQUAL(var1L1, missval) ) w1 = 0;
  if ( DBL_IS_EQUAL(var1L2, missval) ) w2 = 0;

  if ( IS_EQUAL(w1, 0) && IS_EQUAL(w2, 0) )
    {
      var2 = missval;
    }
  else if ( IS_EQUAL(w1, 0) )
    {
      var2 = (w2 >= 0.5) ? var1L2 : missval;	      
    }
  else if ( IS_EQUAL(w2, 0) )
    {
      var2 = (w1 >= 0.5) ? var1L1 : missval;	      
    }
  else
    {
      var2 = var1L1*w1 + var1L2*w2;
    }

  return var2;
}

static
void vert_interp_lev(size_t gridsize, double missval, double *vardata1, double *vardata2,
		     int nlev2, int *lev_idx1, int *lev_idx2, double *lev_wgt1, double *lev_wgt2)
{
  for ( int ilev = 0; ilev < nlev2; ++ilev )
    {
      int idx1 = lev_idx1[ilev];
      int idx2 = lev_idx2[ilev];
      double wgt1 = lev_wgt1[ilev];
      double wgt2 = lev_wgt2[ilev];
      double *var2 = vardata2+gridsize*ilev;

      // if ( cdoVerbose ) cdoPrint("level %d: idx1=%d idx2=%d wgt1=%g wgt2=%g", ilev, idx1, idx2, wgt1, wgt2);

      double *var1L1 = vardata1+gridsize*idx1;
      double *var1L2 = vardata1+gridsize*idx2;

#ifdef  _OPENMP
#pragma omp parallel for default(none) shared(gridsize, var2, var1L1, var1L2, wgt1, wgt2, missval)
#endif
      for ( size_t i = 0; i < gridsize; ++i )
	{
          var2[i] = vert_interp_lev_kernel(wgt1, wgt2, var1L1[i], var1L2[i], missval);
	}
    }
}

/*
 * 3d vertical interpolation routine (see vert_interp_lev() in src/Intlevel.cc)
 */
void vert_interp_lev3d(size_t gridsize, double missval, double *vardata1, double *vardata2,
		       int nlev2, int *lev_idx1, int *lev_idx2, double *lev_wgt1, double *lev_wgt2)
{
  for ( int ilev = 0; ilev < nlev2; ilev++ )
    {
      size_t offset = ilev*gridsize;
      double *var2 = vardata2 + offset;

#ifdef  _OPENMP
#pragma omp parallel for default(none) shared(offset, gridsize, vardata1, var2, lev_idx1, lev_idx2, lev_wgt1, lev_wgt2, missval)
#endif
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


void vert_gen_weights(int expol, int nlev1, double *lev1, int nlev2, double *lev2,
		      int *lev_idx1, int *lev_idx2, double *lev_wgt1, double *lev_wgt2)
{
  int i1;
  int idx1 = 0, idx2 = 0;
  double val1, val2 = 0;

  for ( int i2 = 0; i2 < nlev2; ++i2 )
    {
      // Because 2 levels were added to the source vertical coordinate (one on top, one at the bottom), its loop starts at 1
      for ( i1 = 1; i1 < nlev1; ++i1 )
	{
	  if ( lev1[i1-1] < lev1[i1] )
	    {
	      idx1 = i1-1;
	      idx2 = i1;
	    }
	  else
	    {
	      idx1 = i1;
	      idx2 = i1-1;
	    }
	  val1 = lev1[idx1];
	  val2 = lev1[idx2];

	  if ( lev2[i2] > val1 && lev2[i2] <= val2 ) break;
	}

      if ( i1 == nlev1 ) cdoAbort("Level %g not found!", lev2[i2]);

      if ( i1-1 == 0 ) // destination levels ios not covert by the first two input z levels
        {
          lev_idx1[i2] = 1;
          lev_idx2[i2] = 1;
          lev_wgt1[i2] = 0;
          lev_wgt2[i2] = (expol || IS_EQUAL(lev2[i2], val2));
        }
      else if ( i1 == nlev1-1 ) // destination level is beyond the last value of the input z field
        {
          lev_idx1[i2] = nlev1-2;
          lev_idx2[i2] = nlev1-2;
          lev_wgt1[i2] = (expol || IS_EQUAL(lev2[i2], val2));
          lev_wgt2[i2] = 0;
        }
      else // target z values has two bounday values in input z field
        {
          lev_idx1[i2] = idx1;
          lev_idx2[i2] = idx2;
          lev_wgt1[i2] = (lev1[idx2] - lev2[i2]) / (lev1[idx2] - lev1[idx1]);
          lev_wgt2[i2] = (lev2[i2] - lev1[idx1]) / (lev1[idx2] - lev1[idx1]);
        }
      // backshift of the indices because of the two additional levels in input vertical coordinate
      lev_idx1[i2]--;
      lev_idx2[i2]--;
      /*
        printf("%d %g %d %d %g %g %d %d %g %g\n",
        i2, lev2[i2], idx1, idx2, lev1[idx1], lev1[idx2], 
        lev_idx1[i2], lev_idx2[i2], lev_wgt1[i2], lev_wgt2[i2]);
      */
    }
}

static
void vert_gen_weights3d1d(bool expol, int nlev1, size_t gridsize, double *xlev1, int nlev2, double *lev2,
                          int *xlev_idx1, int *xlev_idx2, double *xlev_wgt1, double *xlev_wgt2)
{
  std::vector<double> lev1(nlev1);
  std::vector<int> lev_idx1(nlev2);
  std::vector<int> lev_idx2(nlev2);
  std::vector<double> lev_wgt1(nlev2);
  std::vector<double> lev_wgt2(nlev2);

  for ( size_t i = 0; i < gridsize; i++ )
    {
      for ( int k = 0; k < nlev1; ++k ) lev1[k] = xlev1[k*gridsize+i];

      vert_gen_weights(expol, nlev1, &lev1[0], nlev2, &lev2[0], &lev_idx1[0], &lev_idx2[0], &lev_wgt1[0], &lev_wgt2[0]);

      for ( int k = 0; k < nlev2; ++k ) xlev_idx1[k*gridsize+i] = lev_idx1[k]*gridsize+i;
      for ( int k = 0; k < nlev2; ++k ) xlev_idx2[k*gridsize+i] = lev_idx2[k]*gridsize+i;
      for ( int k = 0; k < nlev2; ++k ) xlev_wgt1[k*gridsize+i] = lev_wgt1[k];
      for ( int k = 0; k < nlev2; ++k ) xlev_wgt2[k*gridsize+i] = lev_wgt2[k];
    }
}


bool levelDirUp(int nlev, double *lev)
{
  bool lup = (nlev > 1 && lev[1] > lev[0]);
  for ( int i = 1; i < nlev-1; ++i )
    if ( lup && !(lev[i+1] > lev[i]) ) return false;
        
  return lup;
}


bool levelDirDown(int nlev, double *lev)
{
  bool ldown = (nlev > 1 && lev[1] < lev[0]);
  for ( int i = 1; i < nlev-1; ++i )
    if ( ldown && !(lev[i+1] < lev[i]) ) return false;

  return ldown;
}


void vert_init_level_0_and_N(int nlev, size_t gridsize, double *zlevels)
{
  /*
   * Check monotony of vertical levels
   */
  std::vector<double> zlev(nlev);
  for ( int i = 0; i < nlev; ++i ) zlev[i] = zlevels[(i+1)*gridsize]; 
  bool lup = levelDirUp(nlev, &zlev[0]);
  bool ldown = levelDirDown(nlev, &zlev[0]);

  /* Add artificial values for indication of extrapolation areas (lowermost + upmost levels) */
  if ( lup )
    {
      for ( size_t i = 0; i < gridsize; i++ )
        {
          zlevels[i]                     = -1.e33;
          zlevels[(nlev+1)*gridsize + i] =  1.e33;
        }
    }
  else if ( ldown )
    {
      for ( size_t i = 0; i < gridsize; i++ )
        {
          zlevels[i]                     =  1.e33;
          zlevels[(nlev+1)*gridsize + i] = -1.e33;
        }
    }
  else
    cdoWarning("Non monotonic zaxis!");

  if ( cdoVerbose ) for ( int i = 0; i < nlev+2; ++i ) cdoPrint("lev1 %d: %g", i, zlevels[i*gridsize]);
}


void *Intlevel(void *process)
{
  int nrecs;
  int varID, levelID;
  int zaxisID1 = -1;
  int nlevel = 0;

  cdoInitialize(process);

  // clang-format off
  int INTLEVEL   = cdoOperatorAdd("intlevel",  0, 0, NULL);
  int INTLEVELX  = cdoOperatorAdd("intlevelx", 0, 0, NULL);
  // clang-format on

  UNUSED(INTLEVEL);

  int operatorID = cdoOperatorID();

  bool expol = (operatorID == INTLEVELX);

  operatorInputArg("<zvar> levels");

  int argc = operatorArgc();
  char **argv = operatorArgv();
  const char *zvarname = NULL;
  if ( argc > 1 && isalpha(*argv[0]) )
    {
      zvarname = argv[0];
      argc--;
      argv++;
      if ( cdoVerbose ) cdoPrint("zvarname = %s", zvarname);
    }
  lista_t *flista = lista_new(FLT_LISTA);
  int nlev2 = args2flt_lista(argc, argv, flista);
  double *lev2 = (double *) lista_dataptr(flista);

  if ( cdoVerbose ) for ( int i = 0; i < nlev2; ++i ) cdoPrint("lev2 %d: %g", i, lev2[i]);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int nvars = vlistNvars(vlistID1);

  // Find z-variable
  int zaxisID2 = CDI_UNDEFID;
  int nlev1 = 0;
  double *lev1 = NULL;
  int zvarID = CDI_UNDEFID;
  bool zvarIsVarying = false;
  size_t zvarGridsize = 0;
  size_t wisize = 0;
  if ( zvarname )
    {
      for ( varID = 0; varID < nvars; varID++ )
        {
          char varname[CDI_MAX_NAME]; 
          vlistInqVarName(vlistID1, varID, varname);
          if ( STR_IS_EQ(zvarname, varname) )
            {
              zvarID = varID;
              break;
            }
        }

      if ( zvarID == CDI_UNDEFID ) cdoAbort("Variable %s not found!", zvarname);
      zvarIsVarying = vlistInqVarTimetype(vlistID1, zvarID) == TIME_VARYING;
      zvarGridsize = gridInqSize(vlistInqVarGrid(vlistID1, zvarID));
      zaxisID1 = vlistInqVarZaxis(vlistID1, zvarID);
      nlev1 = zaxisInqSize(zaxisID1);

      lev1 = (double*) Malloc(zvarGridsize*(nlev1+2)*sizeof(double));

      zaxisID2 = zaxisCreate(ZAXIS_GENERIC, nlev2);

      char str[CDI_MAX_NAME];
      strcpy(str, "zlev");
      zaxisDefName(zaxisID2, str);
      str[0] = 0;
      vlistInqVarLongname(vlistID1, zvarID, str);
      if ( str[0] ) zaxisDefLongname(zaxisID2, str);
      str[0] = 0;
      vlistInqVarUnits(vlistID1, zvarID, str);
      if ( str[0] ) zaxisDefUnits(zaxisID2, str);

      wisize = zvarGridsize*nlev2;
    }
  else
    {
      int i;
      int nzaxis = vlistNzaxis(vlistID1);
      for ( i = 0; i < nzaxis; i++ )
        {
          int zaxisID = vlistZaxis(vlistID1, i);
          nlevel = zaxisInqSize(zaxisID);
          if ( zaxisInqType(zaxisID) != ZAXIS_HYBRID && zaxisInqType(zaxisID) != ZAXIS_HYBRID_HALF )
            if ( nlevel > 1 )
              {
                zaxisID1 = zaxisID;
                break;
              }
        }
      if ( i == nzaxis ) cdoAbort("No processable variable found!");

      zaxisID2 = zaxisCreate(zaxisInqType(zaxisID1), nlev2);

      char str[CDI_MAX_NAME];
      zaxisInqName(zaxisID1, str);
      zaxisDefName(zaxisID2, str);
      str[0] = 0;
      zaxisInqLongname(zaxisID1, str);
      if ( str[0] ) zaxisDefLongname(zaxisID2, str);
      str[0] = 0;
      zaxisInqUnits(zaxisID1, str);
      if ( str[0] ) zaxisDefUnits(zaxisID2, str);

      zaxisDefDatatype(zaxisID2, zaxisInqDatatype(zaxisID1));

      nlev1 = nlevel;
      lev1 = (double*) Malloc((nlev1+2)*sizeof(double));
      cdoZaxisInqLevels(zaxisID1, lev1+1);

      bool lup = levelDirUp(nlev1, lev1+1);
      bool ldown = levelDirDown(nlev1, lev1+1);

      if ( lup )
        {
          lev1[0]       = -1.e33; 
          lev1[nlev1+1] =  1.e33;
        }
      else if ( ldown )
        {
          lev1[0]       =  1.e33; 
          lev1[nlev1+1] = -1.e33;
        }
      else
        cdoWarning("Non monotonic zaxis!");

      if ( cdoVerbose ) for ( int i = 0; i < nlev1+2; ++i ) cdoPrint("lev1 %d: %g", i, lev1[i]);

      wisize = nlev2;
    }

  zaxisDefLevels(zaxisID2, lev2);
  int nzaxis = vlistNzaxis(vlistID1);
  for ( int i = 0; i < nzaxis; i++ )
    if ( zaxisID1 == vlistZaxis(vlistID1, i) )
      vlistChangeZaxisIndex(vlistID2, i, zaxisID2);

  int *lev_idx1 = (int*) Malloc(wisize*sizeof(int));
  int *lev_idx2 = (int*) Malloc(wisize*sizeof(int));
  double *lev_wgt1 = (double*) Malloc(wisize*sizeof(double));
  double *lev_wgt2 = (double*) Malloc(wisize*sizeof(double));

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());

  pstreamDefVlist(streamID2, vlistID2);

  std::vector<bool> vars(nvars);
  std::vector<bool> varinterp(nvars);
  std::vector<size_t *> varnmiss(nvars);
  std::vector<double *> vardata1(nvars);
  std::vector<double *> vardata2(nvars);

  int maxlev = nlev1 > nlev2 ? nlev1 : nlev2;

  for ( varID = 0; varID < nvars; varID++ )
    {
      int gridID = vlistInqVarGrid(vlistID1, varID);
      int zaxisID = vlistInqVarZaxis(vlistID1, varID);
      size_t gridsize = gridInqSize(gridID);
      int nlevel = zaxisInqSize(zaxisID);

      vardata1[varID] = (double*) Malloc(gridsize*nlevel*sizeof(double));

      if ( zaxisID == zaxisID1 )
	{
	  varinterp[varID] = true;
	  vardata2[varID]  = (double*) Malloc(gridsize*nlev2*sizeof(double));
	  varnmiss[varID]  = (size_t*) Calloc(maxlev, sizeof(size_t));
	}
      else
	{
	  varinterp[varID] = false;
	  vardata2[varID]  = vardata1[varID];
	  varnmiss[varID]  = (size_t*) Malloc(nlevel*sizeof(size_t));
	}
    }

  int tsID = 0;
  while ( (nrecs = pstreamInqTimestep(streamID1, tsID)) )
    {
      for ( varID = 0; varID < nvars; ++varID ) vars[varID] = false;

      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  pstreamInqRecord(streamID1, &varID, &levelID);
	  size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  size_t offset = gridsize*levelID;
	  double *single1 = vardata1[varID] + offset;

	  pstreamReadRecord(streamID1, single1, &varnmiss[varID][levelID]);
	  vars[varID] = true;
	}

      if ( tsID == 0 || zvarIsVarying )
        {
          if ( zvarname )
            {
              for ( levelID = 0; levelID < nlev1; ++levelID )
                {
                  size_t offset1 = zvarGridsize*levelID;
                  size_t offset2 = zvarGridsize*(levelID+1);
                  arrayCopy(zvarGridsize, vardata1[zvarID]+offset1, lev1+offset2);
                }
                  
              vert_init_level_0_and_N(nlev1, zvarGridsize, lev1);
              vert_gen_weights3d1d(expol, nlev1+2, zvarGridsize, lev1, nlev2, lev2, lev_idx1, lev_idx2, lev_wgt1, lev_wgt2);
            }
          else
            vert_gen_weights(expol, nlev1+2, lev1, nlev2, lev2, lev_idx1, lev_idx2, lev_wgt1, lev_wgt2);
        }

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vars[varID] && varinterp[varID] )
	    {
	      int gridID = vlistInqVarGrid(vlistID1, varID);
	      double missval = vlistInqVarMissval(vlistID1, varID);
	      size_t gridsize = gridInqSize(gridID);

              if ( zvarname )
                vert_interp_lev3d(gridsize, missval, vardata1[varID], vardata2[varID],
                                  nlev2, lev_idx1, lev_idx2, lev_wgt1, lev_wgt2);
              else
                vert_interp_lev(gridsize, missval, vardata1[varID], vardata2[varID],
                                nlev2, lev_idx1, lev_idx2, lev_wgt1, lev_wgt2);

	      for ( levelID = 0; levelID < nlev2; levelID++ )
		{
		  size_t offset = gridsize*levelID;
		  double *single2 = vardata2[varID] + offset;
		  varnmiss[varID][levelID] = arrayNumMV(gridsize, single2, missval);
		}
	    }
	}

      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( vars[varID] )
	    {
              size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));
	      int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID));
	      for ( levelID = 0; levelID < nlevel; levelID++ )
		{
		  size_t offset = gridsize*levelID;
		  double *single2 = vardata2[varID] + offset;
		  pstreamDefRecord(streamID2, varID, levelID);
		  pstreamWriteRecord(streamID2, single2, varnmiss[varID][levelID]);
		}
	    }
	}

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  for ( varID = 0; varID < nvars; varID++ )
    {
      Free(varnmiss[varID]);
      Free(vardata1[varID]);
      if ( varinterp[varID] ) Free(vardata2[varID]);
    }

  Free(lev_idx1);
  Free(lev_idx2);
  Free(lev_wgt1);
  Free(lev_wgt2);

  Free(lev1);

  lista_destroy(flista);

  cdoFinish();

  return 0;
}

/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2016 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
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

      Select      select         Select fields
*/
//#define TEST_KVL

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "error.h"
#include "util.h"
#include "pml.h"
#include "pmlist.h"
#ifdef TEST_KVL
#include "sellist.h"
#endif

double datestr_to_double(const char *datestr, int opt);

static
void write_const_vars(int streamID2, int vlistID2, int nvars, double **vardata2)
{
  for ( int varID2c = 0; varID2c < nvars; ++varID2c )
    {
      if ( vardata2[varID2c] )
        {
          double missval = vlistInqVarMissval(vlistID2, varID2c);
          int gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID2c));
          int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID2c));
          for ( int levelID2c = 0; levelID2c < nlevel; ++levelID2c )
            {
              double *pdata = vardata2[varID2c]+gridsize*levelID2c;
              int nmiss = 0;
              for ( int i = 0; i < gridsize; ++i )
                if ( DBL_IS_EQUAL(pdata[i], missval) ) nmiss++;

              // if ( levelID2c == 0 ) printf("Write varID %d\n", varID2c);
              streamDefRecord(streamID2, varID2c, levelID2c);
              streamWriteRecord(streamID2, pdata, nmiss);
            }
          Free(vardata2[varID2c]);
          vardata2[varID2c] = NULL;
        }
    }
}  


void *Select(void *argument)
{
  bool lconstvars = true;
  int streamID2 = CDI_UNDEFID;
  int nrecs;
  int nvars, nvars2;
  int varID, levelID;
  int last_year = -999999999;
  char paramstr[32];
  char varname[CDI_MAX_NAME];
  char stdname[CDI_MAX_NAME];
  char gname[CDI_MAX_NAME];
  char zname[CDI_MAX_NAME];
  int vlistID0 = -1, vlistID2 = -1;
  int taxisID2 = CDI_UNDEFID;
  int ntsteps2 = 0;
  bool ltimsel = false;
  bool *vars = NULL;
  double **vardata2 = NULL;
  double *array = NULL;
  double fstartdate = -99999999999.;
  double fenddate   = -99999999999.;

  cdoInitialize(argument);

  int SELECT = cdoOperatorAdd("select", 0, 0, "parameter list");
  int DELETE = cdoOperatorAdd("delete", 0, 0, "parameter list");

  bool lcopy = UNCHANGED_RECORD;

  int operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  int nsel = operatorArgc();
  char **argnames = operatorArgv();

  if ( nsel == 0 ) cdoAbort("Parameter missing!");

  if ( cdoVerbose )
    for ( int i = 0; i < nsel; ++i )
      cdoPrint("name %d = %s", i+1, argnames[i]);

#ifdef TEST_KVL
  //src/cdo infon -select,name=var129,var152,level=0 ../../test/data/pl_data
  list_t *kvlist = kvlist_new("SELECT");
  if ( kvlist_parse_cmdline(kvlist, nsel, argnames) != 0 ) cdoAbort("Parse error!");
  if ( cdoVerbose ) kvlist_print(kvlist);

  sellist_t *sellist = sellist_create(kvlist);

  SELLIST_ADD_INT(sellist, timestep_of_year, "Timestep of year");
  SELLIST_ADD_INT(sellist, timestep,         "Timestep");
  SELLIST_ADD_INT(sellist, year,             "Year");
  SELLIST_ADD_INT(sellist, month,            "Month");
  SELLIST_ADD_INT(sellist, day,              "Day");
  SELLIST_ADD_INT(sellist, hour,             "Hour");
  SELLIST_ADD_INT(sellist, minute,           "Minute");
  SELLIST_ADD_INT(sellist, code,             "Code number");
  SELLIST_ADD_INT(sellist, levidx,           "Level index");
  SELLIST_ADD_INT(sellist, ltype,            "Level type");
  SELLIST_ADD_INT(sellist, zaxisnum,         "Zaxis number");
  SELLIST_ADD_INT(sellist, gridnum,          "Grid number");
  SELLIST_ADD_FLT(sellist, level,            "Level");
  SELLIST_ADD_WORD(sellist, name,            "Variable name");
  SELLIST_ADD_WORD(sellist, param,           "Parameter");
  SELLIST_ADD_WORD(sellist, zaxisname,       "Zaxis name");
  SELLIST_ADD_WORD(sellist, gridname,        "Grid name");
  SELLIST_ADD_WORD(sellist, steptype,        "Time step type");
  SELLIST_ADD_WORD(sellist, startdate,       "Start date");
  SELLIST_ADD_WORD(sellist, enddate,         "End date");
  SELLIST_ADD_WORD(sellist, season,          "Season");
  SELLIST_ADD_WORD(sellist, date,            "Date");
#else
  pml_t *pml = pml_create("SELECT");

  PML_ADD_INT(pml, timestep_of_year, 4096, "Timestep of year");
  PML_ADD_INT(pml, timestep,         4096, "Timestep");
  PML_ADD_INT(pml, year,             1024, "Year");
  PML_ADD_INT(pml, month,              32, "Month");
  PML_ADD_INT(pml, day,                32, "Day");
  PML_ADD_INT(pml, hour,               24, "Hour");
  PML_ADD_INT(pml, minute,             60, "Minute");
  PML_ADD_INT(pml, code,             1024, "Code number");
  PML_ADD_INT(pml, levidx,           1024, "Level index");
  PML_ADD_INT(pml, ltype,             256, "Level type");
  PML_ADD_INT(pml, zaxisnum,          256, "Zaxis number");
  PML_ADD_INT(pml, gridnum,           256, "Grid number");
  PML_ADD_FLT(pml, level,            1024, "Level");
  PML_ADD_WORD(pml, name,            1024, "Variable name");
  PML_ADD_WORD(pml, param,           1024, "Parameter");
  PML_ADD_WORD(pml, zaxisname,        256, "Zaxis name");
  PML_ADD_WORD(pml, gridname,         256, "Grid name");
  PML_ADD_WORD(pml, steptype,          32, "Time step type");
  PML_ADD_WORD(pml, startdate,          1, "Start date");
  PML_ADD_WORD(pml, enddate,            1, "End date");
  PML_ADD_WORD(pml, season,            12, "Season");
  PML_ADD_WORD(pml, date,            1024, "Date");

  int status = pml_read(pml, nsel, argnames);
  if ( cdoVerbose ) pml_print(pml);
  if ( status != 0 ) cdoAbort("Parameter read error!");
#endif
  int streamCnt = cdoStreamCnt();
  int nfiles = streamCnt - 1;

  if ( !cdoVerbose && nfiles > 1 ) progressInit();

  timestep = 0;
  int tsID2 = 0;
  for ( int indf = 0; indf < nfiles; ++indf )
    {
      if ( !cdoVerbose && nfiles > 1 ) progressStatus(0, 1, (indf+1.)/nfiles);
      if ( cdoVerbose ) cdoPrint("Process file: %s", cdoStreamName(indf)->args);

      int streamID1 = streamOpenRead(cdoStreamName(indf));

      int vlistID1 = streamInqVlist(streamID1);
      int taxisID1 = vlistInqTaxis(vlistID1);

      bool lcopy_const = false;

      if ( indf == 0 )
	{
          bool xresult = true;

	  // vlistID0 = vlistDuplicate(vlistID1);

	  vlistClearFlag(vlistID1);
	  nvars = vlistNvars(vlistID1);
	  vars  = (bool*) Malloc(nvars*sizeof(bool));

	  if ( operatorID == DELETE )
	    {
	      xresult = false;
	      for ( varID = 0; varID < nvars; ++varID )
		{
		  int zaxisID = vlistInqVarZaxis(vlistID1, varID);
		  int nlevs   = zaxisInqSize(zaxisID);
		  for ( int levID = 0; levID < nlevs; ++levID )
		    vlistDefFlag(vlistID1, varID, levID, TRUE);
		}
	    }

#ifdef TEST_KVL
          bool lvarsel = SELLIST_NOCC(code) || SELLIST_NOCC(ltype) || SELLIST_NOCC(zaxisnum) ||
            SELLIST_NOCC(gridnum) || SELLIST_NOCC(name) || SELLIST_NOCC(param) ||
            SELLIST_NOCC(zaxisname) || SELLIST_NOCC(gridname) || SELLIST_NOCC(steptype);

          bool llevsel = SELLIST_NOCC(level) || SELLIST_NOCC(levidx);

	  ltimsel = SELLIST_NOCC(date) || SELLIST_NOCC(startdate) || SELLIST_NOCC(enddate) || SELLIST_NOCC(season) ||
            SELLIST_NOCC(timestep_of_year) || SELLIST_NOCC(timestep) || SELLIST_NOCC(year) || SELLIST_NOCC(month) ||
            SELLIST_NOCC(day) || SELLIST_NOCC(hour) || SELLIST_NOCC(minute);
#else
          bool lvarsel = PML_NOCC(pml, code) || PML_NOCC(pml, ltype) || PML_NOCC(pml, zaxisnum) ||
            PML_NOCC(pml, gridnum) || PML_NOCC(pml, name) || PML_NOCC(pml, param) ||
            PML_NOCC(pml, zaxisname) || PML_NOCC(pml, gridname) || PML_NOCC(pml, steptype);

          bool llevsel = PML_NOCC(pml, level) || PML_NOCC(pml, levidx);

	  ltimsel = PML_NOCC(pml, date) || PML_NOCC(pml, startdate) || PML_NOCC(pml, enddate) || PML_NOCC(pml, season) ||
            PML_NOCC(pml, timestep_of_year) || PML_NOCC(pml, timestep) || PML_NOCC(pml, year) || PML_NOCC(pml, month) ||
            PML_NOCC(pml, day) || PML_NOCC(pml, hour) || PML_NOCC(pml, minute);
#endif
          
	  for ( varID = 0; varID < nvars; ++varID )
	    {
	      int iparam = vlistInqVarParam(vlistID1, varID);
	      code = vlistInqVarCode(vlistID1, varID);
	      vlistInqVarName(vlistID1, varID, varname);
	      vlistInqVarStdname(vlistID1, varID, stdname);

	      cdiParamToString(iparam, paramstr, sizeof(paramstr));

	      name  = varname;
	      param = paramstr;

	      int gridID  = vlistInqVarGrid(vlistID1, varID);
	      int zaxisID = vlistInqVarZaxis(vlistID1, varID);
	      int nlevs   = zaxisInqSize(zaxisID);
	      ltype = zaxis2ltype(zaxisID);

              zaxisnum = vlistZaxisIndex(vlistID1, zaxisID)+1;
              int zaxistype = zaxisInqType(zaxisID);
              zaxisName(zaxistype, zname);
              zaxisname = zname;

              gridnum = vlistGridIndex(vlistID1, gridID)+1;
              int gridtype = gridInqType(gridID);
              gridName(gridtype, gname);
              gridname = gname;

              int tsteptype = vlistInqVarTsteptype(vlistID1, varID);
              if      ( tsteptype == TSTEP_CONSTANT ) steptype = "constant";
              else if ( tsteptype == TSTEP_INSTANT  ) steptype = "instant";
              else if ( tsteptype == TSTEP_INSTANT2 ) steptype = "instant";
              else if ( tsteptype == TSTEP_INSTANT3 ) steptype = "instant";
              else if ( tsteptype == TSTEP_MIN      ) steptype = "min";
              else if ( tsteptype == TSTEP_MAX      ) steptype = "max";
              else if ( tsteptype == TSTEP_AVG      ) steptype = "avg";
              else if ( tsteptype == TSTEP_ACCUM    ) steptype = "accum";
              else if ( tsteptype == TSTEP_RANGE    ) steptype = "range";
              else if ( tsteptype == TSTEP_DIFF     ) steptype = "diff";
              else                                    steptype = "unknown";
              
	      vars[varID] = false;
#ifdef TEST_KVL
              bool found_code  = SELLIST_CHECK(code);
              bool found_name  = SELLIST_CHECK(name);
              bool found_param = SELLIST_CHECK(param);
              bool found_grid  = SELLIST_CHECK(gridnum);
              bool found_gname = SELLIST_CHECK(gridname);
              bool found_stype = SELLIST_CHECK(steptype);
              bool found_ltype = SELLIST_CHECK(ltype);
              bool found_zaxis = SELLIST_CHECK(zaxisnum);
              bool found_zname = SELLIST_CHECK(zaxisname);

              bool lvar  = found_code || found_name || found_param;
              bool lstep = SELLIST_NOCC(steptype) ? found_stype : true;
              bool lgrid = (SELLIST_NOCC(gridnum) || SELLIST_NOCC(gridname)) ? (found_grid || found_gname) : true;
              bool lvert = (SELLIST_NOCC(ltype) || SELLIST_NOCC(zaxisnum) || SELLIST_NOCC(zaxisname)) ? (found_ltype || found_zaxis || found_zname) : true;
#else
              bool found_code  = PML_NOCC(pml, code)      && PML_CHECK_INT(pml, code);
              bool found_name  = PML_NOCC(pml, name)      && PML_CHECK_WORD(pml, name);
              bool found_param = PML_NOCC(pml, param)     && PML_CHECK_WORD(pml, param);
              bool found_grid  = PML_NOCC(pml, gridnum)   && PML_CHECK_INT(pml, gridnum);
              bool found_gname = PML_NOCC(pml, gridname)  && PML_CHECK_WORD(pml, gridname);
              bool found_stype = PML_NOCC(pml, steptype)  && PML_CHECK_WORD(pml, steptype);
              bool found_ltype = PML_NOCC(pml, ltype)     && PML_CHECK_INT(pml, ltype);
              bool found_zaxis = PML_NOCC(pml, zaxisnum)  && PML_CHECK_INT(pml, zaxisnum);
              bool found_zname = PML_NOCC(pml, zaxisname) && PML_CHECK_WORD(pml, zaxisname);

              bool lvar  = found_code || found_name || found_param;
              bool lstep = PML_NOCC(pml, steptype) ? found_stype : true;
              bool lgrid = (PML_NOCC(pml, gridnum) || PML_NOCC(pml, gridname)) ? (found_grid || found_gname) : true;
              bool lvert = (PML_NOCC(pml, ltype) || PML_NOCC(pml, zaxisnum) || PML_NOCC(pml, zaxisname)) ? (found_ltype || found_zaxis || found_zname) : true;
#endif
	     
              if ( !vars[varID] && lgrid && lvar ) vars[varID] = true;
              if ( !vars[varID] && lvert && lvar ) vars[varID] = true;
              if ( !vars[varID] && lstep && lvar ) vars[varID] = true;

              if ( !vars[varID] && !lvar )
                {
                  if      ( found_grid || found_gname ) vars[varID] = true;
                  else if ( found_stype ) vars[varID] = true;
                  else if ( found_ltype || found_zaxis || found_zname ) vars[varID] = true;
#ifdef TEST_KVL
                  else if ( !lvarsel && (SELLIST_NOCC(levidx) || SELLIST_NOCC(level)) )
#else
                  else if ( !lvarsel && (PML_NOCC(pml, levidx) || PML_NOCC(pml, level)) )
#endif
                    {
                      for ( int levID = 0; levID < nlevs; ++levID )
                        {
                          levidx = levID + 1;
                          level = cdoZaxisInqLevel(zaxisID, levID);
#ifdef TEST_KVL
                          if ( !vars[varID] && SELLIST_CHECK(levidx) ) vars[varID] = true;
                          if ( !vars[varID] && SELLIST_CHECK(level)  ) vars[varID] = true;
#else
                          if ( !vars[varID] && PML_NOCC(pml, levidx) && PML_CHECK_INT(pml, levidx) ) vars[varID] = true;
                          if ( !vars[varID] && PML_NOCC(pml, level)  && PML_CHECK_FLT(pml, level)  ) vars[varID] = true;
#endif
                        }
                    }
                }
	    }

	  for ( varID = 0; varID < nvars; ++varID )
	    {
	      if ( vars[varID] )
		{
		  int zaxisID = vlistInqVarZaxis(vlistID1, varID);
                  if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID )
                    {
                      int psvarid = vlist_get_psvarid(vlistID1, zaxisID);
                      if ( psvarid != -1 && !vars[psvarid] ) vars[psvarid] = true;
                    }
                }
            }

	  for ( varID = 0; varID < nvars; ++varID )
	    {
	      if ( vars[varID] )
		{
		  int zaxisID = vlistInqVarZaxis(vlistID1, varID);
		  int nlevs   = zaxisInqSize(zaxisID);
		  for ( int levID = 0; levID < nlevs; ++levID )
		    {
		      levidx = levID + 1;
                      level = cdoZaxisInqLevel(zaxisID, levID);
		      
		      if ( nlevs == 1 && IS_EQUAL(level, 0) )
			{
#ifdef TEST_KVL
                          SELLIST_CHECK(level);
#else
                          PML_CHECK_FLT(pml, level);
#endif
			  vlistDefFlag(vlistID1, varID, levID, xresult);
			}
		      else
			{
#ifdef TEST_KVL
			  if ( SELLIST_NOCC(levidx) )
			    {
			      if ( SELLIST_CHECK(levidx) )
				vlistDefFlag(vlistID1, varID, levID, xresult);
			    }
			  else if ( SELLIST_NOCC(level) )
			    {
			      if ( SELLIST_CHECK(level) )
				vlistDefFlag(vlistID1, varID, levID, xresult);
			    }
#else
			  if ( PML_NOCC(pml, levidx) )
			    {
			      if ( PML_CHECK_INT(pml, levidx) )
				vlistDefFlag(vlistID1, varID, levID, xresult);
			    }
			  else if ( PML_NOCC(pml, level) )
			    {
			      if ( PML_CHECK_FLT(pml, level) )
				vlistDefFlag(vlistID1, varID, levID, xresult);
			    }
#endif
			  else
			    {
			      vlistDefFlag(vlistID1, varID, levID, xresult);
			    }
			}
		    }
		}
	    }

#ifdef TEST_KVL
	  SELLIST_CHECK_FLAG(code);
	  SELLIST_CHECK_FLAG(levidx);
	  SELLIST_CHECK_FLAG(ltype);
	  SELLIST_CHECK_FLAG(zaxisnum);
	  SELLIST_CHECK_FLAG(gridnum);
	  SELLIST_CHECK_FLAG(level);
	  SELLIST_CHECK_FLAG(name);
	  SELLIST_CHECK_FLAG(param);
	  SELLIST_CHECK_FLAG(zaxisname);
	  SELLIST_CHECK_FLAG(gridname);
	  SELLIST_CHECK_FLAG(steptype);
#else
	  PML_CHECK_INT_FLAG(pml, code);
	  PML_CHECK_INT_FLAG(pml, levidx);
	  PML_CHECK_INT_FLAG(pml, ltype);
	  PML_CHECK_INT_FLAG(pml, zaxisnum);
	  PML_CHECK_INT_FLAG(pml, gridnum);
	  PML_CHECK_FLT_FLAG(pml, level);
	  PML_CHECK_WORD_FLAG(pml, name);
	  PML_CHECK_WORD_FLAG(pml, param);
	  PML_CHECK_WORD_FLAG(pml, zaxisname);
	  PML_CHECK_WORD_FLAG(pml, gridname);
	  PML_CHECK_WORD_FLAG(pml, steptype);
#endif

	  int npar = 0;
	  for ( varID = 0; varID < nvars; ++varID )
	    {
	      int zaxisID = vlistInqVarZaxis(vlistID1, varID);
	      int nlevs   = zaxisInqSize(zaxisID);
	      for ( int levID = 0; levID < nlevs; ++levID )
		if ( vlistInqFlag(vlistID1, varID, levID) == TRUE )
                  {
                    npar++;
                    break;
                  }
	    }

	  if ( npar == 0 )
	    {
	      if ( (! lvarsel) && (! llevsel) && ltimsel )
		{
                  lcopy_const = true;

		  for ( varID = 0; varID < nvars; ++varID )
		    {
		      vars[varID] = true;
		      int zaxisID = vlistInqVarZaxis(vlistID1, varID);
		      int nlevs   = zaxisInqSize(zaxisID);
		      for ( int levID = 0; levID < nlevs; ++levID )
			vlistDefFlag(vlistID1, varID, levID, TRUE);
		    }
		}
	      else
		{
		  cdoAbort("No variable selected!");
		}
	    }
          else
            {
              if ( lvarsel && ltimsel )
                {
                  for ( varID = 0; varID < nvars; ++varID )
                    {
                      if ( vars[varID] == true && vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT )
                        {
                          lcopy_const = true;
                          break;
                        }
                    }
                }
            }

	  //if ( cdoVerbose ) vlistPrint(vlistID1);

	  vlistID0 = vlistDuplicate(vlistID1);
	  for ( varID = 0; varID < nvars; ++varID )
	    {
	      int zaxisID = vlistInqVarZaxis(vlistID1, varID);
	      int nlevs   = zaxisInqSize(zaxisID);
	      for ( int levID = 0; levID < nlevs; ++levID )
		vlistDefFlag(vlistID0, varID, levID, vlistInqFlag(vlistID1, varID, levID));
	    }

	  //if ( cdoVerbose ) vlistPrint(vlistID0);

	  vlistID2 = vlistCreate();
	  vlistCopyFlag(vlistID2, vlistID0);

	  //if ( cdoVerbose ) vlistPrint(vlistID2);

	  taxisID2 = taxisDuplicate(taxisID1);
	  vlistDefTaxis(vlistID2, taxisID2);

	  int ntsteps = vlistNtsteps(vlistID1);

	  nvars2 = vlistNvars(vlistID2);

	  if ( ntsteps == 1 && nfiles == 1 )
	    {
	      for ( varID = 0; varID < nvars2; ++varID )
		if ( vlistInqVarTsteptype(vlistID2, varID) != TSTEP_CONSTANT ) break;

	      if ( varID == nvars2 ) ntsteps = 0;
	    }

	  ntsteps2 = ntsteps;
#ifdef TEST_KVL
	  if ( operatorID == SELECT && SELLIST_NOCC(timestep) == 1 ) ntsteps2 = 1;
#else
	  if ( operatorID == SELECT && PML_NOCC(pml, timestep) == 1 ) ntsteps2 = 1;
#endif
	  
	  if ( ntsteps2 == 0 && nfiles > 1 )
	    {
              lconstvars = false;
	      for ( varID = 0; varID < nvars2; ++varID )
		vlistDefVarTsteptype(vlistID2, varID, TSTEP_INSTANT);
	    }

	  // support for negative timestep values
#ifdef TEST_KVL
	  if ( SELLIST_NOCC(timestep) > 0 && ntsteps > 0 && nfiles == 1 )
#else
	  if ( PML_NOCC(pml, timestep) > 0 && ntsteps > 0 && nfiles == 1 )
#endif
	    {
#ifdef TEST_KVL
	      for ( int i = 0; i < SELLIST_NOCC(timestep); ++i )
#else
	      for ( int i = 0; i < PML_NOCC(pml, timestep); ++i )
#endif
		{
#ifdef TEST_KVL
                  int ptimestep;
                  SELLIST_GET_PAR(timestep, i, &ptimestep);
		  if ( ptimestep < 0 )
		    {
		      if ( cdoVerbose )
			cdoPrint("timestep %d changed to %d", ptimestep, ptimestep + ntsteps + 1);
		      ptimestep += ntsteps + 1;
                      SELLIST_DEF_PAR(timestep, i, &ptimestep);
		    }
#else
		  if ( par_timestep[i] < 0 )
		    {
		      if ( cdoVerbose )
			cdoPrint("timestep %d changed to %d", par_timestep[i], par_timestep[i] + ntsteps + 1);
		      par_timestep[i] += ntsteps + 1;
		    }
#endif
		}
	    }

	  if ( ! lcopy )
	    {
	      int gridsize = vlistGridsizeMax(vlistID1);
	      if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
	      array = (double*) Malloc(gridsize*sizeof(double));
	    }

#ifdef TEST_KVL
	  SELLIST_GET_PAR(startdate, 0, &startdate);
	  SELLIST_GET_PAR(enddate, 0, &enddate);
	  if ( SELLIST_NOCC(startdate) ) fstartdate = datestr_to_double(startdate, 0);
	  if ( SELLIST_NOCC(enddate)   ) fenddate   = datestr_to_double(enddate, 1);
#else
	  startdate = par_startdate[0];
	  enddate   = par_enddate[0];
	  if ( PML_NOCC(pml, startdate) ) fstartdate = datestr_to_double(startdate, 0);
	  if ( PML_NOCC(pml, enddate)   ) fenddate   = datestr_to_double(enddate, 1);
#endif
	}
      else
	{
	  vlistCompare(vlistID0, vlistID1, CMP_ALL);
	}


      if ( nvars2 == 0 )
	{
	  cdoWarning("No resulting variables available!");
	  goto END_LABEL;
	}

      if ( lcopy_const )
        {
          vardata2 = (double**) Malloc(nvars2*sizeof(double));
          for ( varID = 0; varID < nvars2; ++varID ) vardata2[varID] = NULL;
        }

      bool lstop = false;
      int tsID1 = 0;
      while ( (nrecs = streamInqTimestep(streamID1, tsID1)) )
	{
          timestep++;
	  bool copytimestep = true;

	  if ( ltimsel == true )
	    {
	      copytimestep = false;

#ifdef TEST_KVL
	      if ( operatorID == SELECT && SELLIST_NOCC(timestep) > 0 )
		{
                  int ptimestep;
                  SELLIST_GET_PAR(timestep, SELLIST_NOCC(timestep)-1, &ptimestep);
                  if ( timestep > ptimestep )
                    {
                      lstop = true;
                      break;
                    }
                }
#else
	      if ( operatorID == SELECT && PML_NOCC(pml, timestep) > 0 && timestep > par_timestep[PML_NOCC(pml, timestep)-1] )
		{
		  lstop = true;
		  break;
		}
#endif

	      int vdate = taxisInqVdate(taxisID1);
	      int vtime = taxisInqVtime(taxisID1);
              int second;
	      cdiDecodeDate(vdate, &year, &month, &day);
	      cdiDecodeTime(vtime, &hour, &minute, &second);
              UNUSED(season);

	      if ( year != last_year )
		{
		  timestep_of_year = 0;
		  last_year = year;
		}

	      timestep_of_year++;

#ifdef TEST_KVL
	      if ( SELLIST_CHECK(timestep) ) copytimestep = true;
	      if ( SELLIST_CHECK(timestep_of_year) ) copytimestep = true;

	      if ( !copytimestep && SELLIST_NOCC(date) == 0 && SELLIST_NOCC(timestep) == 0 && SELLIST_NOCC(timestep_of_year) == 0 )
		{
		  bool lseason = false, lyear = false, lmonth = false, lday = false, lhour = false, lminute = false;

		  if ( SELLIST_NOCC(season) == 0 || SELLIST_CHECK_SEASON(season, month) ) lseason   = true;
		  if ( SELLIST_NOCC(year)   == 0 || SELLIST_CHECK(year)   ) lyear   = true;
		  if ( SELLIST_NOCC(month)  == 0 || SELLIST_CHECK(month)  ) lmonth  = true;
		  if ( SELLIST_NOCC(day)    == 0 || SELLIST_CHECK(day)    ) lday    = true;
		  if ( SELLIST_NOCC(hour)   == 0 || SELLIST_CHECK(hour)   ) lhour   = true;
		  if ( SELLIST_NOCC(minute) == 0 || SELLIST_CHECK(minute) ) lminute = true;

		  if ( lseason && lyear && lmonth && lday && lhour && lminute ) copytimestep = true;
		}
#else
	      if ( PML_NOCC(pml, timestep) && PML_CHECK_INT(pml, timestep) ) copytimestep = true;
	      if ( PML_NOCC(pml, timestep_of_year) && PML_CHECK_INT(pml, timestep_of_year) ) copytimestep = true;

	      if ( !copytimestep && PML_NOCC(pml, date) == 0 && PML_NOCC(pml, timestep) == 0 && PML_NOCC(pml, timestep_of_year) == 0 )
		{
		  bool lseason = false, lyear = false, lmonth = false, lday = false, lhour = false, lminute = false;

		  if ( PML_NOCC(pml, season) == 0 || (PML_NOCC(pml, season) && PML_CHECK_SEASON(pml, season, month)) ) lseason   = true;
		  if ( PML_NOCC(pml, year)   == 0 || (PML_NOCC(pml, year)   && PML_CHECK_INT(pml, year))   ) lyear   = true;
		  if ( PML_NOCC(pml, month)  == 0 || (PML_NOCC(pml, month)  && PML_CHECK_INT(pml, month))  ) lmonth  = true;
		  if ( PML_NOCC(pml, day)    == 0 || (PML_NOCC(pml, day)    && PML_CHECK_INT(pml, day))    ) lday    = true;
		  if ( PML_NOCC(pml, hour)   == 0 || (PML_NOCC(pml, hour)   && PML_CHECK_INT(pml, hour))   ) lhour   = true;
		  if ( PML_NOCC(pml, minute) == 0 || (PML_NOCC(pml, minute) && PML_CHECK_INT(pml, minute)) ) lminute = true;

		  if ( lseason && lyear && lmonth && lday && lhour && lminute ) copytimestep = true;
		}
#endif

	      double fdate = ((double)vdate) + ((double)vtime)/1000000.;

#ifdef TEST_KVL
	      if ( SELLIST_NOCC(enddate) )
#else
	      if ( PML_NOCC(pml, enddate) )
#endif
		{
                  copytimestep = (fdate <= fenddate);
		  if ( fdate > fenddate )
		    {
#ifdef TEST_KVL
                      SELLIST_DEF_FLAG(enddate, 0, true);
#else
		      flag_enddate[0] = true;
#endif
		      if ( operatorID == SELECT )
			{
			  lstop = true;
			  break;
			}
		    }
		}

#ifdef TEST_KVL
	      if ( SELLIST_NOCC(startdate) )
#else
	      if ( PML_NOCC(pml, startdate) )
#endif
		{
                  copytimestep = (fdate >= fstartdate);
#ifdef TEST_KVL
		  if ( fdate >= fstartdate ) SELLIST_DEF_FLAG(startdate, 0, true);
#else
		  if ( fdate >= fstartdate ) flag_startdate[0] = true;
#endif
		}

              
#ifdef TEST_KVL
              if ( SELLIST_NOCC(date) )
#else
              if ( PML_NOCC(pml, date) )
#endif
                {
                  char vdatetimestr[64];
                  datetime2str(vdate, vtime, vdatetimestr, sizeof(vdatetimestr));
                  date = vdatetimestr;
#ifdef TEST_KVL
                  if ( SELLIST_CHECK_DATE(date) ) copytimestep = true;
#else
                  if ( PML_CHECK_DATE(pml, date) ) copytimestep = true;
#endif
                }

	      if ( operatorID == DELETE ) copytimestep = !copytimestep;

              if ( copytimestep && indf == 0 && tsID1 == 0 ) lcopy_const = false;
	    }

	  if ( copytimestep == true )
	    {
	      if ( streamID2 == CDI_UNDEFID )
		{
                  bool lasttimestep = (nfiles == 1) && (ntsteps2 > 1) && (ntsteps2 == (tsID1+1));
                  if ( lasttimestep && tsID2 == 0 ) ntsteps2 = 1;
                  if ( ntsteps2 == 0 || ntsteps2 == 1 ) vlistDefNtsteps(vlistID2, ntsteps2);
		  streamID2 = streamOpenWrite(cdoStreamName(nfiles), cdoFiletype());
		  streamDefVlist(streamID2, vlistID2);
		}

	      taxisCopyTimestep(taxisID2, taxisID1);
	      streamDefTimestep(streamID2, tsID2);
              
	      for ( int recID = 0; recID < nrecs; ++recID )
		{
		  streamInqRecord(streamID1, &varID, &levelID);
		  if ( vlistInqFlag(vlistID0, varID, levelID) == TRUE )
		    {
                      if ( lconstvars && tsID2 > 0 && tsID1 == 0 )
                        if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT )
                          continue;

		      int varID2   = vlistFindVar(vlistID2, varID);
		      int levelID2 = vlistFindLevel(vlistID2, varID, levelID);
		      
                      if ( lcopy_const && tsID2 == 0 )
                        write_const_vars(streamID2, vlistID2, varID2, vardata2);

                      streamDefRecord(streamID2, varID2, levelID2);
                      // if ( levelID2 == 0 ) printf("Write varID %d\n", varID2);
		      if ( lcopy )
			{
			  streamCopyRecord(streamID2, streamID1);
			}
		      else
			{
                          int nmiss;
			  streamReadRecord(streamID1, array, &nmiss);
			  streamWriteRecord(streamID2, array, nmiss);
			}
		    }
		}

              if ( lcopy_const && tsID2 == 0 )
                write_const_vars(streamID2, vlistID2, nvars2, vardata2);

	      tsID2++;              
	    }
          else if ( lcopy_const && indf == 0 && tsID1 == 0 )
            {
	      for ( int recID = 0; recID < nrecs; ++recID )
		{
		  streamInqRecord(streamID1, &varID, &levelID);
		  if ( vlistInqFlag(vlistID0, varID, levelID) == TRUE )
		    {
		      int varID2 = vlistFindVar(vlistID2, varID);
                      if ( vlistInqVarTsteptype(vlistID2, varID2) == TSTEP_CONSTANT )
                        {
                          int levelID2 = vlistFindLevel(vlistID2, varID, levelID);
                          int gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID2));
                          if ( levelID == 0 )
                            {
                              int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID2));
                              vardata2[varID2] = (double*) Malloc(gridsize*nlevel*sizeof(double));
                            }
                          int nmiss;
                          streamReadRecord(streamID1, vardata2[varID2]+gridsize*levelID2, &nmiss);
                        }
		    }
		}
            }

	  tsID1++;
	}
      
      streamClose(streamID1);

      if ( lstop ) break;
    }

 END_LABEL:

  if ( !cdoVerbose && nfiles > 1 ) progressStatus(0, 1, 1);    

#ifdef TEST_KVL
  SELLIST_CHECK_FLAG(timestep_of_year);
  SELLIST_CHECK_FLAG(timestep);
  SELLIST_CHECK_FLAG(year);
  SELLIST_CHECK_FLAG(month);
  SELLIST_CHECK_FLAG(day);
  SELLIST_CHECK_FLAG(hour);
  SELLIST_CHECK_FLAG(minute);
  SELLIST_CHECK_FLAG(startdate);
  //  SELLIST_CHECK_FLAG(enddate);
  SELLIST_CHECK_FLAG(season);
  SELLIST_CHECK_FLAG(date);
#else
  PML_CHECK_INT_FLAG(pml, timestep_of_year);
  PML_CHECK_INT_FLAG(pml, timestep);
  PML_CHECK_INT_FLAG(pml, year);
  PML_CHECK_INT_FLAG(pml, month);
  PML_CHECK_INT_FLAG(pml, day);
  PML_CHECK_INT_FLAG(pml, hour);
  PML_CHECK_INT_FLAG(pml, minute);
  PML_CHECK_WORD_FLAG(pml, startdate);
  //  PML_CHECK_WORD_FLAG(pml, enddate);
  PML_CHECK_WORD_FLAG(pml, season);
  PML_CHECK_WORD_FLAG(pml, date);
#endif

  if ( streamID2 != CDI_UNDEFID ) streamClose(streamID2);

  vlistDestroy(vlistID0);
  vlistDestroy(vlistID2);


#ifdef TEST_KVL
  sellist_destroy(sellist);
  kvlist_destroy(kvlist);
#else
  pml_destroy(pml);
#endif

  if ( array ) Free(array);
  if ( vars ) Free(vars);
  if ( vardata2 )
    {
      for ( varID = 0; varID < nvars2; ++varID )
        if ( vardata2[varID] ) Free(vardata2[varID]);

      Free(vardata2);
    }

  if ( tsID2 == 0 ) cdoAbort("No timesteps selected!");

  cdoFinish();

  return 0;
}

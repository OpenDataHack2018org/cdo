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

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "error.h"
#include "util.h"
#include "pmlist.h"

double datestr_to_double(const char *datestr, int opt);

int vlist_get_psvarid(int vlistID, int zaxisID)
{
  int psvarid = -1;
  char name[CDI_MAX_NAME];
  char psname[CDI_MAX_NAME];
  psname[0] = 0;
  zaxisInqPsName(zaxisID, psname);

  if ( psname[0] )
    {
      int nvars = vlistNvars(vlistID);
      for ( int varID = 0; varID < nvars; ++varID )
        {
          vlistInqVarName(vlistID, varID, name);
          if ( strcmp(name, psname) == 0 )
            {
              psvarid = varID;
              break;
            }
        }
      if ( cdoVerbose && psvarid == -1 )
        cdoWarning("Surface pressure variable not found - %s", psname);
    }

  return psvarid;
}


void *Select(void *argument)
{
  bool lconstvars = true;
  int streamID2 = CDI_UNDEFID;
  int nrecs;
  int nvars, nvars2, nlevs;
  int zaxisID;
  int varID2, levelID2;
  int varID, levelID;
  int last_year = -999999999;
  char paramstr[32];
  char varname[CDI_MAX_NAME];
  char stdname[CDI_MAX_NAME];
  char gname[CDI_MAX_NAME];
  char zname[CDI_MAX_NAME];
  int vlistID0 = -1, vlistID2 = -1;
  int result = FALSE;
  int taxisID2 = CDI_UNDEFID;
  int ntsteps;
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

  if ( cdoVerbose )
    for ( int i = 0; i < nsel; i++ )
      printf("name %d = %s\n", i+1, argnames[i]);

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

  int streamCnt = cdoStreamCnt();
  int nfiles = streamCnt - 1;

  if ( !cdoVerbose && nfiles > 1 ) progressInit();

  timestep = 0;
  int tsID2 = 0;
  for ( int indf = 0; indf < nfiles; indf++ )
    {
      if ( !cdoVerbose && nfiles > 1 ) progressStatus(0, 1, (indf+1.)/nfiles);
      if ( cdoVerbose ) cdoPrint("Process file: %s", cdoStreamName(indf)->args);

      int streamID1 = streamOpenRead(cdoStreamName(indf));

      int vlistID1 = streamInqVlist(streamID1);
      int taxisID1 = vlistInqTaxis(vlistID1);

      bool lcopy_const = false;

      if ( indf == 0 )
	{
	  // vlistID0 = vlistDuplicate(vlistID1);

	  vlistClearFlag(vlistID1);
	  nvars = vlistNvars(vlistID1);
	  vars  = (bool*) Malloc(nvars*sizeof(bool));

	  if ( operatorID == DELETE )
	    {
	      result = FALSE;
	      for ( varID = 0; varID < nvars; varID++ )
		{
		  zaxisID = vlistInqVarZaxis(vlistID1, varID);
		  nlevs   = zaxisInqSize(zaxisID);
		  for ( int levID = 0; levID < nlevs; levID++ )
		    vlistDefFlag(vlistID1, varID, levID, TRUE);
		}
	    }
	  else if ( operatorID == SELECT )
	    {
	      result = TRUE;
	    }

	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      int iparam = vlistInqVarParam(vlistID1, varID);
	      code = vlistInqVarCode(vlistID1, varID);
	      vlistInqVarName(vlistID1, varID, varname);
	      vlistInqVarStdname(vlistID1, varID, stdname);

	      cdiParamToString(iparam, paramstr, sizeof(paramstr));

	      name  = varname;
	      param = paramstr;

	      int gridID  = vlistInqVarGrid(vlistID1, varID);
	      zaxisID = vlistInqVarZaxis(vlistID1, varID);
	      nlevs   = zaxisInqSize(zaxisID);
	      ltype   = zaxis2ltype(zaxisID);

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
	     
              if ( !vars[varID] && lgrid && lvar) vars[varID] = true;
              if ( !vars[varID] && lvert && lvar) vars[varID] = true;
              if ( !vars[varID] && lstep && lvar) vars[varID] = true;
              if ( !vars[varID] && !lvar )
                {
                  if      ( found_grid || found_gname ) vars[varID] = true;
                  else if ( found_stype ) vars[varID] = true;
                  else if ( found_ltype || found_zaxis || found_zname ) vars[varID] = true;
                  else if ( PML_NOCC(pml, levidx) || PML_NOCC(pml, level) )
                    {
                      for ( int levID = 0; levID < nlevs; levID++ )
                        {
                          levidx = levID + 1;
                          level = zaxisInqLevel(zaxisID, levID);
                          if ( !vars[varID] && PML_NOCC(pml, levidx) && PML_CHECK_INT(pml, levidx) )  vars[varID] = true;
                          if ( !vars[varID] && PML_NOCC(pml, level)  && PML_CHECK_FLT(pml, level)  )  vars[varID] = true;
                        }
                    }
                }
	    }

	  for ( varID = 0; varID < nvars; varID++ )
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

	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      if ( vars[varID] )
		{
		  zaxisID = vlistInqVarZaxis(vlistID1, varID);
		  nlevs   = zaxisInqSize(zaxisID);

		  for ( int levID = 0; levID < nlevs; levID++ )
		    {
		      levidx = levID + 1;
		      level = zaxisInqLevel(zaxisID, levID);
		      
		      if ( nlevs == 1 && IS_EQUAL(level, 0) )
			{
			  vlistDefFlag(vlistID1, varID, levID, result);
			}
		      else
			{
			  if ( PML_NOCC(pml, levidx) )
			    {
			      if ( PML_CHECK_INT(pml, levidx) )
				vlistDefFlag(vlistID1, varID, levID, result);
			    }
			  else if ( PML_NOCC(pml, level) )
			    {
			      if ( PML_CHECK_FLT(pml, level) )
				vlistDefFlag(vlistID1, varID, levID, result);
			    }
			  else
			    {
			      vlistDefFlag(vlistID1, varID, levID, result);
			    }
			}
		    }
		}
	    }

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

	  if ( PML_NOCC(pml, date) || PML_NOCC(pml, startdate) || PML_NOCC(pml, enddate) || PML_NOCC(pml, season) ) ltimsel = true;
	  if ( PML_NOCC(pml, timestep_of_year) || PML_NOCC(pml, timestep) || PML_NOCC(pml, year) || PML_NOCC(pml, month) || PML_NOCC(pml, day) || PML_NOCC(pml, hour) || PML_NOCC(pml, minute) ) ltimsel = true;

	  int npar = 0;
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      zaxisID = vlistInqVarZaxis(vlistID1, varID);
	      nlevs   = zaxisInqSize(zaxisID);
	      for ( int levID = 0; levID < nlevs; levID++ )
		if ( vlistInqFlag(vlistID1, varID, levID) == result )
                  {
                    npar++;
                    break;
                  }
	    }

	  if ( npar == 0 )
	    {
	      if ( ltimsel == true )
		{
                  lcopy_const = true;

		  for ( varID = 0; varID < nvars; varID++ )
		    {
		      vars[varID] = true;
		      zaxisID = vlistInqVarZaxis(vlistID1, varID);
		      nlevs   = zaxisInqSize(zaxisID);
		      for ( int levID = 0; levID < nlevs; levID++ )
			vlistDefFlag(vlistID1, varID, levID, TRUE);
		    }
		}
	      else
		{
		  cdoAbort("No variable selected!");
		}
	    }

	  //if ( cdoVerbose ) vlistPrint(vlistID1);

	  vlistID0 = vlistDuplicate(vlistID1);
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      zaxisID = vlistInqVarZaxis(vlistID1, varID);
	      nlevs   = zaxisInqSize(zaxisID);
	      for ( int levID = 0; levID < nlevs; levID++ )
		vlistDefFlag(vlistID0, varID, levID, vlistInqFlag(vlistID1, varID, levID));
	    }

	  //if ( cdoVerbose ) vlistPrint(vlistID0);

	  vlistID2 = vlistCreate();
	  vlistCopyFlag(vlistID2, vlistID0);

	  //if ( cdoVerbose ) vlistPrint(vlistID2);

	  taxisID2 = taxisDuplicate(taxisID1);
	  vlistDefTaxis(vlistID2, taxisID2);

	  ntsteps = vlistNtsteps(vlistID1);

	  nvars2 = vlistNvars(vlistID2);

	  if ( ntsteps == 1 && nfiles == 1 )
	    {
	      for ( varID = 0; varID < nvars2; ++varID )
		if ( vlistInqVarTsteptype(vlistID2, varID) != TSTEP_CONSTANT ) break;

	      if ( varID == nvars2 ) ntsteps = 0;
	    }

	  int ntsteps2 = ntsteps;
	  if ( operatorID == SELECT && PML_NOCC(pml, timestep) == 1 ) ntsteps2 = 1;
	  
	  if ( ntsteps2 == 0 || ntsteps2 == 1 ) vlistDefNtsteps(vlistID2, ntsteps2);

	  if ( ntsteps2 == 0 && nfiles > 1 )
	    {
              lconstvars = false;
	      for ( varID = 0; varID < nvars2; ++varID )
		vlistDefVarTsteptype(vlistID2, varID, TSTEP_INSTANT);
	    }

	  // support for negative timestep values
	  if ( PML_NOCC(pml, timestep) > 0 && ntsteps > 0 && nfiles == 1 )
	    {
	      for ( int i = 0; i < PML_NOCC(pml, timestep); i++ )
		{
		  if ( par_timestep[i] < 0 )
		    {
		      if ( cdoVerbose )
			cdoPrint("timestep %d changed to %d", par_timestep[i], par_timestep[i] + ntsteps + 1 );
		      par_timestep[i] += ntsteps + 1;
		    }
		}
	    }

	  if ( ! lcopy )
	    {
	      int gridsize = vlistGridsizeMax(vlistID1);
	      if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
	      array = (double*) Malloc(gridsize*sizeof(double));
	    }

	  startdate = par_startdate[0];
	  enddate   = par_enddate[0];
	  if ( PML_NOCC(pml, startdate) ) fstartdate = datestr_to_double(startdate, 0);
	  if ( PML_NOCC(pml, enddate)   ) fenddate   = datestr_to_double(enddate, 1);
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
          vardata2 = (double**) malloc(nvars2*sizeof(double));
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

	      if ( operatorID == SELECT && PML_NOCC(pml, timestep) > 0 && timestep > par_timestep[PML_NOCC(pml, timestep)-1] )
		{
		  lstop = true;
		  break;
		}

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

	      double fdate = ((double)vdate) + ((double)vtime)/1000000.;

	      if ( PML_NOCC(pml, enddate) )
		{
		  if ( fdate > fenddate )
		    {
		      flag_enddate[0] = true;
		      copytimestep = false;
		      if ( operatorID == SELECT )
			{
			  lstop = true;
			  break;
			}
		    }
		  else
		    {
		      copytimestep = true;
		    }
		}

	      if ( PML_NOCC(pml, startdate) )
		{
		  if ( fdate < fstartdate )
		    {
		      copytimestep = false;
		    }
		  else
		    {
		      flag_startdate[0] = true;
		      copytimestep = true;
		    }
		}

              
              if ( PML_NOCC(pml, date) )
                {
                  char vdatetimestr[64];
                  datetime2str(vdate, vtime, vdatetimestr, sizeof(vdatetimestr));
                  date = vdatetimestr;
                  if ( PML_CHECK_DATE(pml, date) ) copytimestep = true;
                }

	      if ( operatorID == DELETE ) copytimestep = !copytimestep;

              if ( copytimestep && indf == 0 && tsID1 == 0 ) lcopy_const = false;
	    }

	  if ( copytimestep == true )
	    {
	      if ( streamID2 == CDI_UNDEFID )
		{
		  streamID2 = streamOpenWrite(cdoStreamName(nfiles), cdoFiletype());
		  streamDefVlist(streamID2, vlistID2);
		}

	      taxisCopyTimestep(taxisID2, taxisID1);
	      streamDefTimestep(streamID2, tsID2);

              if ( lcopy_const && tsID2 == 0 )
                {
                  for ( varID2 = 0; varID2 < nvars2; ++varID2 )
                    {
                      if ( vardata2[varID2] )
                        {
                          double missval = vlistInqVarMissval(vlistID2, varID2);
                          int gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID2));
                          int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID2));
                          for ( int levelID2 = 0; levelID2 < nlevel; ++levelID2 )
                            {
                              double *pdata = vardata2[varID2]+gridsize*levelID;
                              int nmiss = 0;
                              for ( int i = 0; i < gridsize; ++i )
                                if ( DBL_IS_EQUAL(pdata[i], missval) ) nmiss++;

                              streamDefRecord(streamID2, varID2, levelID2);
                              streamWriteRecord(streamID2, pdata, nmiss);
                            }
                        }
                    }
                }
              
	      for ( int recID = 0; recID < nrecs; recID++ )
		{
		  streamInqRecord(streamID1, &varID, &levelID);
		  if ( vlistInqFlag(vlistID0, varID, levelID) == TRUE )
		    {
                      if ( lconstvars && tsID2 > 0 && tsID1 == 0 )
                        if ( vlistInqVarTsteptype(vlistID1, varID) == TSTEP_CONSTANT )
                          continue;

		      varID2   = vlistFindVar(vlistID2, varID);
		      levelID2 = vlistFindLevel(vlistID2, varID, levelID);
		      
		      streamDefRecord(streamID2, varID2, levelID2);
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

	      tsID2++;              
	    }
          else if ( lcopy_const && indf == 0 && tsID1 == 0 )
            {
	      for ( int recID = 0; recID < nrecs; recID++ )
		{
		  streamInqRecord(streamID1, &varID, &levelID);
		  if ( vlistInqFlag(vlistID0, varID, levelID) == TRUE )
		    {
		      varID2 = vlistFindVar(vlistID2, varID);
                      if ( vlistInqVarTsteptype(vlistID2, varID2) == TSTEP_CONSTANT )
                        {
                          levelID2 = vlistFindLevel(vlistID2, varID, levelID);
                          int gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID2));
                          if ( levelID == 0 )
                            {
                              int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID2, varID2));
                              vardata2[varID2] = (double*) malloc(gridsize*nlevel*sizeof(double));
                            }
                          int nmiss;
                          streamReadRecord(streamID1, vardata2[varID2]+gridsize*levelID, &nmiss);
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

  if ( streamID2 != CDI_UNDEFID ) streamClose(streamID2);

  vlistDestroy(vlistID0);
  vlistDestroy(vlistID2);

  pml_destroy(pml);

  if ( array ) Free(array);
  if ( vars ) Free(vars);
  if ( vardata2 )
    {
      for ( varID = 0; varID < nvars2; ++ varID )
        if ( vardata2[varID] ) free(vardata2[varID]);

      free(vardata2);
    }

  if ( tsID2 == 0 ) cdoAbort("No timesteps selected!");

  cdoFinish();

  return 0;
}

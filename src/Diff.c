/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2008 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Diff       diff            Compare two datasets
*/


#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Diff(void *argument)
{
  static char func[] = "Diff";
  int DIFF, DIFFV, SDIFF;
  int operatorID;
  int i;
  int indg;
  int varID, recID;
  int gridsize;
  int gridID, zaxisID, code, vdate, vtime;
  int nrecs, nrecs2;
  int levelID;
  int tsID;
  int dsgn, zero;
  int streamID1, streamID2;
  int vlistID1, vlistID2;
  int taxisID;
  int nmiss1, nmiss2;
  int ndrec = 0, nd2rec = 0, ngrec = 0;
  char varname[128];
  int year, month, day, hour, minute;
  double *array1, *array2;
  double absm, relm;
  double missval1, missval2;

  cdoInitialize(argument);

  DIFF  = cdoOperatorAdd("diff",  0, 0, NULL);
  DIFFV = cdoOperatorAdd("diffv", 0, 0, NULL);
  SDIFF = cdoOperatorAdd("sdiff", 0, 0, NULL);

  operatorID = cdoOperatorID();

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  streamID2 = streamOpenRead(cdoStreamName(1));
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = streamInqVlist(streamID2);

  vlistCompare(vlistID1, vlistID2, func_sft);

  gridsize = vlistGridsizeMax(vlistID1);

  array1 = (double *) malloc(gridsize*sizeof(double));
  array2 = (double *) malloc(gridsize*sizeof(double));

  if ( operatorID == DIFFV )
    fprintf(stdout, "               Date  Time Varname     Level    Size    Miss :"
	    " S Z      Absdiff     Reldiff\n");
  else if ( operatorID == DIFF )
    fprintf(stdout, "               Date  Time Code  Level    Size    Miss :"
	    " S Z      Absdiff     Reldiff\n");

  indg = 0;
  tsID = 0;
  taxisID = vlistInqTaxis(vlistID1);
  while ( TRUE )
    {
      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs > 0 )
	{
	  vdate = taxisInqVdate(taxisID);
	  vtime = taxisInqVtime(taxisID);

	  decode_date(vdate, &year, &month, &day);
	  decode_time(vtime, &hour, &minute);
	}

      nrecs2 = streamInqTimestep(streamID2, tsID);

      if ( nrecs == 0 || nrecs2 == 0 ) break;

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamInqRecord(streamID2, &varID, &levelID);

	  indg += 1;

	  code     = vlistInqVarCode(vlistID1, varID);
	  gridID   = vlistInqVarGrid(vlistID1, varID);
	  zaxisID  = vlistInqVarZaxis(vlistID1, varID);
	  gridsize = gridInqSize(gridID);
	  missval1 = vlistInqVarMissval(vlistID1, varID);
	  missval2 = vlistInqVarMissval(vlistID2, varID);

	  if ( operatorID == DIFFV || operatorID == DIFF )
	    {
	      if ( operatorID == DIFFV ) vlistInqVarName(vlistID1, varID, varname);

	      if ( operatorID == DIFFV )
		fprintf(stdout, "%6d : %4.4d-%2.2d-%2.2d %2.2d:%2.2d %-8s ",
			indg, year, month, day, hour, minute, varname);
	      else if ( operatorID == DIFF )
		fprintf(stdout, "%6d : %4.4d-%2.2d-%2.2d %2.2d:%2.2d %3d",
			indg, year, month, day, hour, minute, code);

	      fprintf(stdout, " %7g ", zaxisInqLevel(zaxisID, levelID));
	    }

	  streamReadRecord(streamID1, array1, &nmiss1);
	  streamReadRecord(streamID2, array2, &nmiss2);

          absm = 0.0;
	  relm = 0.0;
	  dsgn = FALSE;
          zero = FALSE;

	  for ( i = 0; i < gridsize; i++ )
	    {
	      if ( !DBL_IS_EQUAL(array1[i], missval1) && !DBL_IS_EQUAL(array2[i], missval2) )
		{
		  absm = MAX(absm, fabs(array1[i]-array2[i]));
		  if ( array1[i]*array2[i] < 0 )
		    dsgn = TRUE;
		  else if ( DBL_IS_EQUAL(array1[i]*array2[i], 0) )
		    zero = TRUE;
		  else
		    relm = MAX(relm, fabs(array1[i]-array2[i]) /
			   MAX(fabs(array1[i]), fabs(array2[i])));
		}
	      else if ( (DBL_IS_EQUAL(array1[i], missval1) && !DBL_IS_EQUAL(array2[i], missval2)) || 
			(!DBL_IS_EQUAL(array1[i], missval1) && DBL_IS_EQUAL(array2[i], missval2)) )
		{
		  relm = 1.0;
		}
	    }

	  if ( operatorID == DIFFV || operatorID == DIFF )
	    {
	      fprintf(stdout, "%7d %7d :", gridsize, MAX(nmiss1, nmiss2));

	      fprintf(stdout, " %c %c ", dsgn ? 'T' : 'F', zero ? 'T' : 'F');
	      fprintf(stdout, "%#12.5g%#12.5g\n", absm, relm);
	    }

	  ngrec++;
	  if ( absm > 0     || relm > 0     ) ndrec++;
	  if ( absm > 1.e-3 || relm > 1.e-3 ) nd2rec++;
	}
      tsID++;
    }

  fprintf(stdout, "  %d of %d records differ\n", ndrec, ngrec);
  if ( ndrec != nd2rec )
    fprintf(stdout, "  %d of %d records differ more than 0.001\n", nd2rec, ngrec);
  /*  fprintf(stdout, "  %d of %d records differ more then one thousandth\n", nprec, ngrec); */
  if ( nrecs == 0 && nrecs2 > 0 )
    cdoWarning("stream2 has more time steps than stream1!");
  if ( nrecs > 0 && nrecs2 == 0 )
    cdoWarning("stream1 has more time steps than stream2!");

  streamClose(streamID1);
  streamClose(streamID2);

  if ( array1 ) free(array1);
  if ( array2 ) free(array2);

  cdoFinish();

  return (0);
}

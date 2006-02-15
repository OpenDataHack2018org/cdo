/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2006 Uwe Schulzweida, schulzweida@dkrz.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "functs.h"

/*
@BeginDoc

@BeginModule

@Name      = Fldstat
@Title     = Field statistic
@Section   = Statistical description of the data
@Class     = Statistic
@Arguments = ifile ofile
@Operators = fldmin fldmax fldsum fldmean fldavg fldstd fldvar

@BeginDesciption
In this program there is the different notion of "mean" and "average"
to distinguish two different kinds of treatment of missing values:
While computing the mean, only the not missing values are considered
to belong to the sample with the side effect of a probably reduced sample
size. Computing the average is just adding the sample members and divide
the result by the sample size. For example, the mean of 1, 2, miss and 3
is (1+2+3)/3 = 2, whereas the average is (1+2+miss+3)/4 = miss/4 = miss.
If there are no missing values in the sample, the average and the mean are
identical.
In this chapter the abbreviations as in the following table are used:

\vspace{3mm}

\fbox{\parbox{15cm}{
\begin{eqnarray*}
\begin{array}{l}
\makebox[3cm][l]{{\bf mean} resp. {\bf avg}} \\
\end{array}
 &  &
 n^{-1} \sum_{i=1}^{n} x_i 
\\
\begin{array}{l}
\makebox[3cm][l]{{\bf mean} resp. {\bf avg}} \\
\mbox{weighted by} \\
\{w_i, i=1, ..., n\}  \\
\end{array}
 & &
  \left ( \sum_{j=1}^{n} w_j \right )^{-1} \sum\limits_{i=1}^{n} w_i \, x_i \\
\\
\begin{array}{l}
%\makebox[3cm][l]{Standard deviation} \\
\makebox[3cm][l]{Variance} \\
\makebox[3cm][l]{{\bf var}} \\
\end{array}
 &  &
 n^{-1} \sum_{i=1}^{n} (x_i - \overline{x})^2
\\
\begin{array}{l}
\makebox[3cm][l]{{\bf var} weighted by} \\
\{w_i, i=1, ..., n\}  \\
\end{array}
 & &
  \left ( \sum_{j=1}^{n} w_j \right )^{-1} \sum\limits_{i=1}^{n} w_i \, 
  \left ( x_i - \left ( \sum_{j=1}^{n} w_j \right )^{-1} \sum\limits_{j=1}^{n} w_j \, x_j \right)^2 \\
\end{eqnarray*}
}}

@EndDesciption

@EndModule


@BeginOperator_fldmin

@Title     = Field minimum

@BeginDesciption
@IfMan
o(t,1) = min{i(t',x'), t'=t}
@EndifMan
@IfDoc
@BeginMath
o(t,1) = \mbox{\bf min}\{i(t',x'), t' = t\}
@EndMath
@EndifDoc
@EndDesciption

@EndOperator

@BeginOperator_fldmax

@Title     = Field maximum

@BeginDesciption
@IfMan
o(t,1) = max{i(t',x'), t'=t}
@EndifMan
@IfDoc
@BeginMath
o(t,1) = \mbox{\bf max}\{i(t',x'), t' = t\}
@EndMath
@EndifDoc
@EndDesciption

@EndOperator

@BeginOperator_fldsum

@Title     = Field sum

@BeginDesciption
@IfMan
o(t,1) = sum{i(t,x)}
@EndifMan
@IfDoc
@BeginMath
o(t,1) = \sum\limits_x i(t,x)
@EndMath
@EndifDoc
@EndDesciption

@EndOperator

@BeginOperator_fldmean

@Title     = Field mean

@BeginDesciption
@IfMan
o(t,1) = mean{i(t',x'), t'=t}
@EndifMan
@IfDoc
@BeginMath
o(t,1) = \mbox{\bf mean}\{i(t',x'), t' = t\}
@EndMath
@EndifDoc
weighted by area weights obtained by the input field.
@EndDesciption

@EndOperator

@BeginOperator_fldavg

@Title     = Field average

@BeginDesciption
@IfMan
o(t,1) = avg{i(t',x'), t'=t}
@EndifMan
@IfDoc
@BeginMath
o(t,1) = \mbox{\bf avg}\{i(t',x'), t' = t\}
@EndMath
@EndifDoc
weighted by area weights obtained by the input field.
@EndDesciption

@EndOperator

@BeginOperator_fldvar

@Title     = Field variance

@BeginDesciption
@IfMan
o(t,1) = var{i(t',x'), t'=t}
@EndifMan
@IfDoc
@BeginMath
o(t,1) = \mbox{\bf var}\{i(t',x'), t' = t\}
@EndMath
@EndifDoc
weighted by area weights obtained by the input field.
@EndDesciption

@EndOperator

@BeginOperator_fldstd

@Title     = Field standard deviation

@BeginDesciption
@IfMan
o(t,1) = sqrt{var{i(t',x'), t'=t}}
@EndifMan
@IfDoc
@BeginMath
o(t,1) = \sqrt{\mbox{\bf var}\{i(t',x'), t' = t\}}
@EndMath
@EndifDoc
weighted by area weights obtained by the input field.
@EndDesciption

@EndOperator

@EndDoc
*/


void *Fldstat(void *argument)
{
  static char func[] = "Fldstat";
  int operatorID;
  int operfunc;
  int streamID1, streamID2;
  int vlistID1, vlistID2;
  int gridID2, lastgrid = -1;
  int index, ngrids;
  int recID, nrecs;
  int tsID, varID, levelID;
  int lim;
  int needWeights = FALSE;
  int nmiss;
  double slon, slat;
  double sglval;
  FIELD field;
  int taxisID1, taxisID2;

  cdoInitialize(argument);

  cdoOperatorAdd("fldmin",  func_min,  0, NULL);
  cdoOperatorAdd("fldmax",  func_max,  0, NULL);
  cdoOperatorAdd("fldsum",  func_sum,  0, NULL);
  cdoOperatorAdd("fldmean", func_mean, 0, NULL);
  cdoOperatorAdd("fldavg",  func_avg,  0, NULL);
  cdoOperatorAdd("fldvar",  func_var,  0, NULL);
  cdoOperatorAdd("fldstd",  func_std,  0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

  if ( operfunc == func_mean || operfunc == func_avg ||
       operfunc == func_var  || operfunc == func_std )
    needWeights = TRUE;

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  slon = 0;
  slat = 0;
  gridID2 = gridNew(GRID_LONLAT, 1);
  gridDefXsize(gridID2, 1);
  gridDefYsize(gridID2, 1);
  gridDefXvals(gridID2, &slon);
  gridDefYvals(gridID2, &slat);

  ngrids = vlistNgrids(vlistID1);
  index = 0;

  vlistChangeGridIndex(vlistID2, index, gridID2);
  if ( ngrids > 1 ) cdoAbort("Too many different grids!");

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  lim = vlistGridsizeMax(vlistID1);
  field.ptr    = (double *) malloc(lim*sizeof(double));
  field.weight = NULL;
  if ( needWeights )
    field.weight = (double *) malloc(lim*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, field.ptr, &field.nmiss);

	  field.grid = vlistInqVarGrid(vlistID1, varID);
	  field.size = gridInqSize(field.grid);
	  if ( field.grid != lastgrid )
	    {
	      lastgrid = field.grid;
	      if ( needWeights ) gridWeights(field.grid, field.weight);
	    }
	  field.missval = vlistInqVarMissval(vlistID1, varID);

	  sglval = fldfun(field, operfunc);

	  if ( DBL_IS_EQUAL(sglval, field.missval) )
	    nmiss = 1;
	  else
	    nmiss = 0;

	  streamDefRecord(streamID2, varID,  levelID);
	  streamWriteRecord(streamID2, &sglval, nmiss);
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( field.ptr )    free(field.ptr);
  if ( field.weight ) free(field.weight);

  cdoFinish();

  return (0);
}

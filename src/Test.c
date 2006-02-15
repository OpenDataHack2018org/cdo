/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2005 Uwe Schulzweida, schulzweida@dkrz.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include <stdio.h>

#include "cdi.h"
#include "cdo.h"
#include "pstream.h"

/*
@BeginDoc

@BeginModule

@Name      = 
@Title     = 
@Section   = 
@Arguments = ifile ofile
@Operators = 

@EndModule


@BeginOperator_

@Title     = 
@Parameter = 

@BeginDesciption
@EndDesciption

@BeginParameter
@Item = 
@EndParameter

@EndOperator

@EndDoc
*/

void *Test(void *argument)
{
  /*
  int streamID1, streamID2;
  */

  cdoInitialize(argument);

  /*
  streamID1 = streamOpenRead(cdoStreamName(0));
  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamClose(streamID2);
  streamClose(streamID1);
  */
  cdoFinish();

  return (0);
}

void *Test2(void *argument)
{
  /*
  int streamID1, streamID2, streamID3;
  */

  cdoInitialize(argument);

  /*
  streamID1 = streamOpenRead(cdoStreamName(0));
  streamID2 = streamOpenRead(cdoStreamName(1));
  streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());

  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);
  */
  cdoFinish();

  return (0);
}

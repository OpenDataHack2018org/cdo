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

#ifndef _PROCESS_H
#define _PROCESS_H

int  processSelf(void);
int  processCreate(void);
void processDelete(void);
int  processInqTimesteps(void);
void processDefTimesteps(int streamID);
int  processInqVarNum(void);
int  processInqStreamNum(void);
int  processInqStreamID(int streamindex);
void processAddStream(int streamID);
void processDefVarNum(int nvars, int streamID);
void processDefArgument(const char *argument);

void processStartTime(double *utime, double *stime);
void processEndTime(double *utime, double *stime);
void processAccuTime(double utime, double stime);

int  processInqChildNum(void);

const char *processOperatorArg(void);
const char *processInqOpername(void);
const char *processInqPrompt(void);

#endif  /* _PROCESS_H */

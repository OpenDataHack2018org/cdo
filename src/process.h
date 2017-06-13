/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2017 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
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

#include <sys/types.h> /* off_t */
#include "util.h"

constexpr int MAX_PROCESS  =   128;
constexpr int MAX_STREAM   =    64;
constexpr int MAX_OPERATOR =   128;
constexpr int MAX_OARGC    =  4096;
constexpr int MAX_FILES    = 65536;


typedef struct {
  int         f1;
  int         f2;
  const char *name;
  const char *enter;
}
oper_t;

typedef struct {
#if defined(HAVE_LIBPTHREAD)
  pthread_t   threadID;
  int         l_threadID;
#endif
  short       nchild;
  short       nInStream;
  short       nOutStream;
  short       inputStreams[MAX_STREAM];
  short       outputStreams[MAX_STREAM];
  double      s_utime;
  double      s_stime;
  double      a_utime;
  double      a_stime;
  double      cputime;

  off_t       nvals;
  short       nvars;
  int         ntimesteps;
  short       streamCnt;
  argument_t *streamNames;
  char       *xoperator;
  const char *operatorName;
  char       *operatorArg;
  int         oargc;
  char       *oargv[MAX_OARGC];
  char        prompt[64];
  short       noper;
  oper_t      oper[MAX_OPERATOR];
}
process_t;



int  processSelf(void);
int  processCreate(void);
void processDelete(void);
int  processInqTimesteps(void);
void processDefTimesteps(int streamID);
int  processInqVarNum(void);
int  processInqInputStreamNum(void);
int  processInqOutputStreamNum(void);
int  processInqInputStreamID(int streamindex);
int  processInqOutputStreamID(int streamindex);
void processAddInputStream(int streamID);
void processAddOutputStream(int streamID);
void processDelStream(int streamID);
void processDefVarNum(int nvars, int streamID);
void processDefArgument(void *vargument);

void processStartTime(double *utime, double *stime);
void processEndTime(double *utime, double *stime);
void processAccuTime(double utime, double stime);

void processDefCputime(int processID, double cputime);
double processInqCputime(int processID);

void processAddNvals(off_t nvals);
off_t processInqNvals(int processID);
int processNums(void);

int  processInqChildNum(void);

const char *processOperatorArg(void);
const char *processInqOpername(void);
const char *processInqOpername2(int processID);
const char *processInqPrompt(void);

void print_process(int p_process_id);
#endif  /* _PROCESS_H */

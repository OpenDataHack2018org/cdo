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

#ifndef PROCESS_INT_H
#define PROCESS_INT_H

#include "util.h"
#include "process.h"

extern std::map<int, process_t> Process;

void cdoInitialize(void *process);
void cdoFinish();
void printEndTimes(cdoTimes p_times, char *p_memstring);

int cdoOperatorAdd(const char *name, int f1, int f2, const char *enter);
int cdoOperatorID(void);
int cdoOperatorF1(int operID);
int cdoOperatorF2(int operID);

const char *cdoOperatorName(int operID);
const char *cdoOperatorEnter(int operID);
std::string cdoGetInStreamName(int p_inStream);
std::string cdoGetOutStreamName(int p_outStream);
std::string cdoGetStreamName(int p_streamIndex);
char* cdoGetObase();

int cdoStreamNumber();
int cdoStreamCnt(void);
int cdoStreamName(int cnt);

int cdoStreamOpenRead(int inStreamIDX);
int cdoStreamOpenWrite(int outStreamIDX, int filetype);
int cdoStreamOpenWrite(std::string p_filename, int filetype);
int cdoStreamOpenAppend(int outStreamIDX);

bool cdoOutFileExists(int outStreamIDX);
bool cdoInFileExists(int inStreamIDX);
int checkStreamCnt();


void operatorCheckArgc(int numargs);
void operatorInputArg(const char *enter);
char ** operatorArgv(void);
int operatorArgc(void);

void processDefVarNum(int nvars);
void processDefTimesteps(int streamID);

int processInqVarNum();
int processInqTimesteps(void);
const char *processInqPrompt(void);

int processNums(void);
int processNumsActive();

process_t &processSelf(void);
process_t* getProcess(int p_processID);

void clearProcesses();

void processAccuTime(double utime, double stime);
void processStartTime(double *utime, double *stime);

void createProcesses(int argc, const char **argv);
process_t *processCreate(const char *command);

int cdoStreamInqVlist(int pstreamID);

#endif

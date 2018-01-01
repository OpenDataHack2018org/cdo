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

#ifndef _PROCESS_H
#define _PROCESS_H

#include <vector>
#include "util.h"
#include "pstream.h"
#include "modules.h"

#include <vector>
#include <iostream>

constexpr int MAX_PROCESS = 128;
constexpr int MAX_STREAM = 64;
constexpr int MAX_OPERATOR = 128;
constexpr int MAX_OARGC = 4096;
constexpr int MAX_FILES = 65536;

enum class ProcessCheckResult{UNLIMITED_STREAM_COUNTS, INPUT_STREAM_MISSING, OUTPUT_STREAM_MISSING, TOO_MANY_STREAMS, TOO_FEW_STREAMS, FILENAME_HAS_OPERATOR_MARKER, OUTFILE_IS_INFILE, SUCCESS};

typedef struct
{
  int f1;
  int f2;
  const char *name;
  const char *enter;
} oper_t;

class process_t
{
public:
  int m_ID;
#ifdef  HAVE_LIBPTHREAD
  pthread_t threadID;
  int l_threadID;
#endif
  short nchild;
  std::vector<process_t *> childProcesses;
  std::vector<process_t *> parentProcesses;
  std::vector<pstream_t *> inputStreams;
  std::vector<pstream_t *> outputStreams;
  short m_cntIn;
  short m_cntOut;

  double s_utime;
  double s_stime;
  double a_utime;
  double a_stime;
  double cputime;

  size_t nvals;
  short nvars;

  int ntimesteps;
  int m_streamCnt;
  std::vector<argument_t> streamArguments;
  const char *m_operatorCommand;
  const char *operatorName;
  char *operatorArg;
  char prompt[64];
  short noper;

  modules_t m_module;
  std::vector<char *> oargv;
  int oargc;
  oper_t oper[MAX_OPERATOR];

  int getInStreamCnt();
  int getOutStreamCnt();
  void initProcess();
  void print_process();
  void defArgument();
  void setOperatorArgv(const char *operatorArguments);
  void setStreams(int argc, std::vector<char *> &argv);
  void addChild(process_t *child_process);
  void addParent(process_t *parent_process);
  bool hasAllInputs();
  process_t(int p_ID, const char *operatorCommand);

private:
  void defPrompt();
  process_t();
  void OpenRead(int p_input_idx);
  void OpenWrite(int p_input_idx);
  void OpenAppend(int p_input_idx);
  void setStreamNames(int argc, std::vector<char *> &argv);
  int expand_wildcards(int streamCnt);
  int checkStreamCnt();
};

extern std::map<int, process_t> Process;

pstream_t *processInqInputStream(int streamindex);
pstream_t *processInqOutputStream(int streamindex);
process_t &processSelf(void);
process_t *processCreate(void);
process_t *processCreate(const char *command);
void processDelete(void);
int processInqTimesteps(void);
void processDefTimesteps(int streamID);
int processInqVarNum(void);
int processInqInputStreamNum(void);
int processInqOutputStreamNum(void);
void processAddInputStream(pstream_t *p_pstream_ptr);
void processAddOutputStream(pstream_t *p_pstream_ptr);
void processDelStream(int streamID);
void processDefVarNum(int nvars);
void processDefArgument(void *vargument);

void processStartTime(double *utime, double *stime);
void processEndTime(double *utime, double *stime);
void processAccuTime(double utime, double stime);

void processDefCputime(int processID, double cputime);
double processInqCputime(int processID);

void processAddNvals(size_t nvals);
size_t processInqNvals(int processID);
int processNums(void);

int processInqChildNum(void);

const char *processOperatorArg(void);
const char *processInqOpername(void);
const char *processInqOpername2(int processID);
const char *processInqPrompt(void);

const argument_t *cdoStreamName(int cnt);
void createProcesses(int argc, const char **argv);
void clearProcesses();
int processNumsActive();

  int checkStreamCnt();
#endif /* _PROCESS_H */

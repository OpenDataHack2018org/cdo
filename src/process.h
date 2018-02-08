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

#include "pstream.h"
#include "modules.h"

#include <vector>
#include <iostream>
#include <string>

constexpr int MAX_PROCESS = 128;
constexpr int MAX_STREAM = 64;
constexpr int MAX_OPERATOR = 128;
constexpr int MAX_OARGC = 4096;
constexpr int MAX_FILES = 65536;

struct cdoTimes
{
  double s_utime, s_stime;
  double e_utime, e_stime;
  double c_cputime = 0, c_usertime = 0, c_systime = 0;
  double p_cputime = 0, p_usertime = 0, p_systime = 0;
};

struct oper_t
{
  int f1;
  int f2;
  const char *name;
  const char *enter;
};

class ProcessType
{
public:
  int m_ID;
  int m_posInParent;
#ifdef  HAVE_LIBPTHREAD
  pthread_t threadID;
  int l_threadID;
#endif
  short nchild;
  std::vector<ProcessType *> childProcesses;
  std::vector<ProcessType *> parentProcesses;
  std::vector<PstreamType *> inputStreams;
  std::vector<PstreamType *> outputStreams;
  int nChildActive = 0;
  short m_cntIn;
  short m_cntOut;

  double s_utime;
  double s_stime;
  double a_utime;
  double a_stime;
  double cputime;

  size_t m_nvals;
  short nvars;

  int ntimesteps;
  int m_streamCnt;
  const char *m_operatorCommand;
  const char *operatorName;
  char *operatorArg;
  char prompt[64];
  short m_noper;
  bool m_isActive;  /*TEMP*/ //not used right now, maybe later (12.Jan.2018)

  module_t m_module;
  std::vector<char *> m_oargv;
  int m_oargc;
  oper_t oper[MAX_OPERATOR];

  ProcessType(int p_ID, const char* p_operatorNamme, const char *operatorCommand);

  pthread_t run();

  int getInStreamCnt();
  int getOutStreamCnt();
  void initProcess();
  void print_process();
  void setOperatorArgv(const char *operatorArguments);
  void addChild(ProcessType *child_process);
  void addParent(ProcessType *parent_process);
  bool hasAllInputs();
  void addFileInStream(std::string file);
  void addFileOutStream(std::string file);
  void addPipeInStream();
  void addPipeOutStream();
  void addNvals(size_t p_nvals);
  void query_user_exit(const char *argument);
  void inqUserInputForOpArg(const char *enter);
  int operatorAdd(const char *name, int f1, int f2, const char *enter);
  int getOperatorID();
  void setInactive();
  const char *inqPrompt();
  cdoTimes getTimes(int p_processNums);
  void printProcessedValues();
  void printBenchmarks(cdoTimes p_times, char *p_memstring);
  int checkStreamCnt();

private:
  ProcessType();
  void defPrompt();
};

#endif /* _PROCESS_H */


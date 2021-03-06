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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_PTHREAD_H
#include <pthread.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#include <stdio.h>
#include <string.h>

#include "cdo_int.h"
#include "error.h"
#include "modules.h"
#include "util.h"
#include "pstream_int.h"
#include "dmemory.h"
#include "pthread.h"
#include "cdoDebugOutput.h"
#include "util_string.h"
#include "util_files.h"
#include "text.h"
#include "cdoOptions.h"
#include <iostream>
#include <string>

ProcessType::ProcessType(int p_ID, const char *p_operatorName, const char *operatorCommand)
    : m_ID(p_ID), operatorName(p_operatorName)
{
  initProcess();
  setOperatorArgv(operatorCommand);
  m_operatorCommand = operatorCommand;

  defPrompt();  // has to be called after get operatorName

  m_module = getModule(p_operatorName);
}

void
ProcessType::setOperatorArgv(const char *operatorArguments)
{
  if (operatorArguments)
    {
      char *operatorArg = strdup(operatorArguments);
      // fprintf(stderr, "processDefArgument: %d %s\n", oargc, operatorArg);

      while ((operatorArg = strchr(operatorArg, ',')) != NULL)
        {
          *operatorArg = '\0';
          operatorArg++;
          if (strlen(operatorArg))
            {
              m_oargv.push_back(operatorArg);
            }
        }
    }
  m_oargc = m_oargv.size();
}

void
ProcessType::initProcess()
{
#ifdef HAVE_LIBPTHREAD
  threadID = pthread_self();
#endif
  m_isActive = true;
  nchild = 0;

  cdoProcessTime(&s_utime, &s_stime);
  a_utime = 0;
  a_stime = 0;
  cputime = 0;
  m_nvals = 0;
  nvars = 0;
  ntimesteps = 0;

  m_streamCnt = 0;

  m_oargc = 0;
  m_operatorCommand = "UNINITALIZED";
  operatorArg = "UNINITALIZED";

  m_noper = 0;
}

int
ProcessType::getInStreamCnt()
{
  return inputStreams.size();
}
int
ProcessType::getOutStreamCnt()
{
  return outputStreams.size();
}

void
ProcessType::addNvals(size_t p_nvals)
{
  m_nvals += p_nvals;
}

void
ProcessType::defPrompt()
{
  if (m_ID == 0)
    sprintf(prompt, "%s %s", Cdo::progname, operatorName);
  else
    sprintf(prompt, "%s(%d) %s", Cdo::progname, m_ID + 1, operatorName);
}

const char *
ProcessType::inqPrompt()
{
  const char *newPrompt = "cdo";
  if (prompt[0]) newPrompt = prompt;

  return newPrompt;
}

void
ProcessType::handleProcessErr(ProcessStatus p_proErr)
{
  switch (p_proErr)
    {
    case ProcessStatus::UnlimitedIOCounts:
      {
        CdoError::Abort(Cdo::progname, "I/O stream counts unlimited no allowed!");
        break;
      }
    case ProcessStatus::MissInput:
      {
        CdoError::Abort(Cdo::progname, "Input streams missing!");
        break;
      }
    case ProcessStatus::MissOutput:
      {
        CdoError::Abort(Cdo::progname, "Output streams missing!");
        break;
      }
    case ProcessStatus::TooManyStreams:
      {
        CdoError::Abort(Cdo::progname, "Too many streams specified! Operator ", m_operatorCommand, " needs ", m_module.streamInCnt,
                        " input and ", m_module.streamOutCnt, " output streams.");
        break;
      }
    case ProcessStatus::TooFewStreams:
      {
        CdoError::Abort(Cdo::progname, "Too few streams specified! Operator ", m_operatorCommand, " needs ", m_module.streamInCnt,
                        " input and ", m_module.streamOutCnt, " output streams.");
        break;
      }
    }
}

void
ProcessType::validate()
{

  ProcessStatus processStatus = checkStreamCnt();
  if (processStatus != ProcessStatus::Ok)
    {
      handleProcessErr(processStatus);
    }
  int errorIdx = checkInFileStreams();
  if (errorIdx != -1)
    {
      CdoError::Abort(Cdo::progname, "Input file ", inputStreams[errorIdx]->m_mfnames[0].c_str(), " of process ", operatorName,
                      " does not exists");
    }
}
ProcessStatus
ProcessType::checkStreamCnt()
{
  int wantedStreamInCnt, wantedStreamOutCnt;
  int streamInCnt0;
  int streamCnt = 0;
  int obase = FALSE;

  wantedStreamInCnt = m_module.streamInCnt;
  wantedStreamOutCnt = m_module.streamOutCnt;

  streamInCnt0 = wantedStreamInCnt;

  if (wantedStreamOutCnt == -1)
    {
      wantedStreamOutCnt = 1;
      obase = TRUE;
    }

  if (wantedStreamInCnt == -1 && wantedStreamOutCnt == -1) return ProcessStatus::UnlimitedIOCounts;

  // printf(" wantedStreamInCnt,wantedStreamOutCnt %d %d\n",
  // wantedStreamInCnt,wantedStreamOutCnt);
  if (wantedStreamInCnt == -1)
    {
      wantedStreamInCnt = m_streamCnt - wantedStreamOutCnt;
      if (wantedStreamInCnt < 1) return ProcessStatus::MissInput;
    }

  if (wantedStreamOutCnt == -1)
    {
      wantedStreamOutCnt = m_streamCnt - wantedStreamInCnt;
      if (wantedStreamOutCnt < 1) return ProcessStatus::MissOutput;
    }
  // printf(" wantedStreamInCnt,wantedStreamOutCnt %d %d\n",
  // wantedStreamInCnt,wantedStreamOutCnt);

  streamCnt = wantedStreamInCnt + wantedStreamOutCnt;
  // printf(" streamCnt %d %d\n", m_streamCnt, streamCnt);

  if (m_streamCnt > streamCnt) return ProcessStatus::TooManyStreams;

  if (m_streamCnt < streamCnt && !obase) return ProcessStatus::TooFewStreams;

  if (wantedStreamInCnt == 1 && streamInCnt0 == -1) return ProcessStatus::Ok;

  return ProcessStatus::Ok;
}

int
ProcessType::checkInFileStreams()
{
  for (unsigned int i = 0; i < inputStreams.size(); i++)
    {
      if (!inputStreams[i]->isPipe() && !fileExists(inputStreams[i]->m_mfnames[0].c_str()))
        {
          return i;
        }
    }
  return -1;
}

bool
ProcessType::hasAllInputs()
{
  if (m_module.streamInCnt == -1)
    {
      return false;
    }
  return m_module.streamInCnt == static_cast<short>(inputStreams.size());
}

void
ProcessType::setInactive()
{
  m_isActive = false;
}

void
ProcessType::inqUserInputForOpArg(const char *enter)
{
  int oargc = m_oargc;

  if (oargc == 0)
    {
      char line[1024];
      char *pline = line;
      size_t pos, len, linelen;
      int lreadline = 1;

      if (enter)
        {
          set_text_color(stderr, BRIGHT, MAGENTA);
          fprintf(stderr, "%-16s : ", prompt);
          reset_text_color(stderr);
          // set_text_color(stderr, BLINK, BLACK);
          fprintf(stderr, "Enter %s > ", enter);
          // reset_text_color(stderr);
        }

      while (lreadline)
        {
          readline(stdin, pline, 1024);

          lreadline = 0;
          while (1)
            {

              pos = 0;
              while (pline[pos] == ' ' || pline[pos] == ',') pos++;
              pline += pos;
              linelen = strlen(pline);
              if (linelen > 0)
                {
                  if (pline[0] == '\\')
                    {
                      lreadline = 1;
                      break;
                    }
                  len = 0;
                  while (pline[len] != ' ' && pline[len] != ',' && pline[len] != '\\' && len < linelen) len++;

                  m_oargv.push_back((char *) Malloc(len + 1));
                  memcpy(m_oargv[oargc], pline, len);
                  m_oargv[oargc][len] = '\0';
                  oargc++;

                  pline += len;
                }
              else
                break;
            }
        }

      m_oargc = oargc;
    }
}

int
ProcessType::operatorAdd(const char *name, int f1, int f2, const char *enter)
{
  int operID = m_noper;

  if (operID >= MAX_OPERATOR) cdoAbort("Maximum number of %d operators reached!", MAX_OPERATOR);

  oper[m_noper].f1 = f1;
  oper[m_noper].f2 = f2;
  oper[m_noper].name = name;
  oper[m_noper].enter = enter;

  m_noper++;

  return operID;
}

int
ProcessType::getOperatorID()
{
  int operID = -1;

  if (m_noper > 0)
    {
      for (operID = 0; operID < m_noper; operID++)
        {
          if (oper[operID].name)
            {
              if (strcmp(operatorName, oper[operID].name) == 0)
                {
                  break;
                }
            }
        }
      if (operID == m_noper)
        {
          CdoError::Abort(Cdo::progname, prompt, "Operator not callable by this name! Name is:", operatorName);
        }
    }
  else
    {
      cdoAbort("Operator not initialized!");
    }

  return operID;
}

void
ProcessType::addFileInStream(std::string file)
{
  inputStreams.push_back(create_pstream(file));
  m_streamCnt++;
}

void
ProcessType::addFileOutStream(std::string file)
{
  if (file[0] == '-')
    {
      CdoError::Abort(Cdo::progname, "Missing output file. Found an operator instead of filename: ", file);
    }
  outputStreams.push_back(create_pstream(file));
  m_streamCnt++;
}

void
ProcessType::addChild(ProcessType *childProcess)
{
  childProcesses.push_back(childProcess);
  nchild = childProcesses.size();
  addPipeInStream();
}
void
ProcessType::addPipeInStream()
{
#ifdef HAVE_LIBPTHREAD
  inputStreams.push_back(create_pstream(m_ID, inputStreams.size()));
  m_streamCnt++;
#else
  cdoAbort("Cannot use pipes, pthread support not compiled in!");
#endif
}

void
ProcessType::addParent(ProcessType *parentProcess)
{
  parentProcesses.push_back(parentProcess);
  m_posInParent = parentProcess->inputStreams.size() - 1;
  addPipeOutStream();
}
void
ProcessType::addPipeOutStream()
{
  outputStreams.push_back(parentProcesses[0]->inputStreams[m_posInParent]);
  m_streamCnt++;
}

pthread_t
ProcessType::run()
{
  Cdo_Debug(CdoDebug::PROCESS, "starting new thread for process ", m_ID);
  pthread_attr_t attr;
  int status = pthread_attr_init(&attr);
  if (status) SysError("pthread_attr_init failed for '%s'", operatorName);
  status = pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  if (status) SysError("pthread_attr_setdetachstate failed for '%s'", operatorName);
  /*
    param.sched_priority = 0;
    status = pthread_attr_setschedparam(&attr, &param);
    if ( status ) SysError("pthread_attr_setschedparam failed for '%s'",
    newarg+1);
  */
  /* status = pthread_attr_setinheritsched(&attr, PTHREAD_EXPLICIT_SCHED); */
  /* if ( status ) SysError("pthread_attr_setinheritsched failed for '%s'",
   * newarg+1); */

  pthread_attr_getscope(&attr, &pthreadScope);

  /* status = pthread_attr_setscope(&attr, PTHREAD_SCOPE_PROCESS); */
  /* if ( status ) SysError("pthread_attr_setscope failed for '%s'", newarg+1);
   */
  /* If system scheduling scope is specified, then the thread is scheduled
   * against all threads in the system */
  /* pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM); */

  size_t stacksize = 0;
  status = pthread_attr_getstacksize(&attr, &stacksize);
  if (stacksize < 2097152)
    {
      stacksize = 2097152;
      pthread_attr_setstacksize(&attr, stacksize);
    }

  pthread_t thrID;
  int rval = pthread_create(&thrID, &attr, m_module.func, this);
  if (rval != 0)
    {
      errno = rval;
      SysError("pthread_create failed for '%s'", operatorName);
    }
  m_isActive = true;
  return thrID;
}

void
ProcessType::query_user_exit(const char *argument)
{
/* modified code from NCO */
#define USR_RPL_MAX_LNG 10 /* Maximum length for user reply */
#define USR_RPL_MAX_NBR 10 /* Maximum number of chances for user to reply */
  char usr_rpl[USR_RPL_MAX_LNG];
  int usr_rpl_int;
  short nbr_itr = 0;
  size_t usr_rpl_lng = 0;

  /* Initialize user reply string */
  usr_rpl[0] = 'z';
  usr_rpl[1] = '\0';

  while (!(usr_rpl_lng == 1 && (*usr_rpl == 'o' || *usr_rpl == 'O' || *usr_rpl == 'e' || *usr_rpl == 'E')))
    {
      if (nbr_itr++ > USR_RPL_MAX_NBR)
        {
          (void) fprintf(stdout, "\n%s: ERROR %d failed attempts to obtain valid interactive input.\n", prompt, nbr_itr - 1);
          exit(EXIT_FAILURE);
        }

      if (nbr_itr > 1) (void) fprintf(stdout, "%s: ERROR Invalid response.\n", prompt);
      (void) fprintf(stdout, "%s: %s exists ---`e'xit, or `o'verwrite (delete existing file) (e/o)? ", prompt, argument);
      (void) fflush(stdout);
      if (fgets(usr_rpl, USR_RPL_MAX_LNG, stdin) == NULL) continue;

      /* Ensure last character in input string is \n and replace that with \0 */
      usr_rpl_lng = strlen(usr_rpl);
      if (usr_rpl_lng >= 1)
        if (usr_rpl[usr_rpl_lng - 1] == '\n')
          {
            usr_rpl[usr_rpl_lng - 1] = '\0';
            usr_rpl_lng--;
          }
    }

  /* Ensure one case statement for each exit condition in preceding while loop */
  usr_rpl_int = (int) usr_rpl[0];
  switch (usr_rpl_int)
    {
    case 'E':
    case 'e': exit(EXIT_SUCCESS); break;
    case 'O':
    case 'o': break;
    default: exit(EXIT_FAILURE); break;
    } /* end switch */
}

/** function for operators with obase usage, will be called while operator execution */

cdoTimes
ProcessType::getTimes(int p_processNums)
{
  cdoTimes times;
  times.s_utime = s_utime;
  times.s_stime = s_stime;

  // both variables are set in cdoProcessTime
  cdoProcessTime(&times.e_utime, &times.e_stime);

  times.c_usertime = times.e_utime - times.s_utime;
  times.c_systime = times.e_stime - times.s_stime;
  times.c_cputime = times.c_usertime + times.c_systime;

#ifdef HAVE_LIBPTHREAD
  if (pthreadScope == PTHREAD_SCOPE_PROCESS)
    {
      times.c_usertime /= p_processNums;
      times.c_systime /= p_processNums;
      times.c_cputime /= p_processNums;
    }
#endif

  cputime = times.c_cputime;

  return times;
}

void
ProcessType::printBenchmarks(cdoTimes p_times, char *p_memstring)
{
#ifdef HAVE_SYS_TIMES_H
  if (m_ID == 0)
    {
      if (!Options::silentMode) fprintf(stderr, " [%.2fs%s]\n", p_times.c_cputime, p_memstring);
    }
  else
    {
      fprintf(stderr, "\n");
    }
#else
  fprintf(stderr, "\n");
#endif
}

void
ProcessType::printProcessedValues()
{
  set_text_color(stderr, RESET, GREEN);
  fprintf(stderr, "%s: ", processInqPrompt());
  reset_text_color(stderr);

  int64_t nvals = m_nvals;

  if (nvals > 0)
    {
      if (sizeof(int64_t) > sizeof(size_t))
#ifdef _WIN32
        fprintf(stderr, "Processed %I64d value%s from %d variable%s",
#else
        fprintf(stderr, "Processed %jd value%s from %d variable%s",
#endif
                (intmax_t) nvals, ADD_PLURAL(nvals), nvars, ADD_PLURAL(nvars));
      else
        fprintf(stderr, "Processed %zu value%s from %d variable%s", (size_t) nvals, ADD_PLURAL(nvals), nvars, ADD_PLURAL(nvars));
    }
  else if (nvars > 0)
    {
      fprintf(stderr, "Processed %d variable%s", nvars, ADD_PLURAL(nvars));
    }

  if (ntimesteps > 0) fprintf(stderr, " over %d timestep%s", ntimesteps, ADD_PLURAL(ntimesteps));

  //  fprintf(stderr, ".");
}
bool
ProcessType::hasOutStream(PstreamType *p_streamPtr)
{
  for (PstreamType *streamPtr : outputStreams)
    {
      if (streamPtr == p_streamPtr)
        {
          return true;
        }
    }
  return false;
}
bool
ProcessType::hasInStream(PstreamType *p_streamPtr)
{
  for (PstreamType *streamPtr : inputStreams)
    {
      if (streamPtr == p_streamPtr)
        {
          return true;
        }
    }
  return false;
}

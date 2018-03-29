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

#if defined(HAVE_PTHREAD_H)
#include <pthread.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
#include <sys/stat.h> /* stat */

#include <map>
#include <stack>
#include <vector>

// Debug and message includes
#include "cdoDebugOutput.h"
#include "exception.h"
#include <cdi.h>

#include "cdoOptions.h"
#include "exception.h"
#include "process_int.h"
#include "pstream_int.h"
#include "text.h"
#include "util_files.h"
#include "util_operatorStrings.h"

std::map<int, ProcessType> Process;
std::map<int, char *> obase;
std::vector<pthread_t> threadIDs;
/*TEMP*/  // Possibly not the best solution (19.Jan.2018)

static int NumProcess = 0;
static int NumProcessActive = 0;

#ifdef HAVE_LIBPTHREAD
pthread_mutex_t processMutex = PTHREAD_MUTEX_INITIALIZER;
#endif

ProcessType &
processSelf(void)
{
#ifdef HAVE_LIBPTHREAD
  pthread_t thID = pthread_self();

  pthread_mutex_lock(&processMutex);

  for (auto &id_process_pair : Process)
    {
      if (id_process_pair.second.l_threadID)
        {
          if (pthread_equal(id_process_pair.second.threadID, thID))
            {
              pthread_mutex_unlock(&processMutex);
              return id_process_pair.second;
            }
        }
    }

  pthread_mutex_unlock(&processMutex);
#endif

  return Process.find(0)->second;
}

void
processDefVarNum(int nvars)
{
  ProcessType &process = processSelf();
  process.nvars += nvars;
}

int
processInqVarNum(void)
{
  return processSelf().nvars;
}

int
cdoStreamInqTimestep(int pstreamID, int tsID)
{
  int nrecs = 0;
  Cdo_Debug(CdoDebug::PROCESS, "Inquiring Timestep ", tsID, " for stream ", pstreamID);
  ProcessType &process = processSelf();
  PstreamType *pstreamptr = pstreamToPointer(pstreamID);
  nrecs = pstreamInqTimestep(pstreamptr, tsID);
  if (nrecs && tsID != pstreamptr->tsID)
    {
      process.timesteps.insert(tsID);
      process.ntimesteps = process.timesteps.size();
    }
  return nrecs;
}

int
operatorArgc(void)
{
  return processSelf().m_oargc;
}

char **
operatorArgv(void)
{
  if (CdoDebug::PROCESS)
    {
      std::string oargv_str = "";
      for (auto entry : processSelf().m_oargv)
        {
          oargv_str += std::string(entry) + " ";
        }
      Cdo_Debug(CdoDebug::PROCESS, "Getting ", processSelf().m_oargv.size(), " operator arguments: ", oargv_str);
    }

  return &processSelf().m_oargv[0];
}

void
operatorCheckArgc(int numargs)
{
  int argc = processSelf().m_oargc;

  if (argc < numargs)
    cdoAbort("Too few arguments! Need %d found %d.", numargs, argc);
  else if (argc > numargs)
    cdoAbort("Too many arguments! Need %d found %d.", numargs, argc);
}

void
operatorInputArg(const char *enter)
{
  ProcessType &process = processSelf();
  process.inqUserInputForOpArg(enter);
}

int
cdoOperatorAdd(const char *name, int f1, int f2, const char *enter)
{
  ProcessType &process = processSelf();
  return process.operatorAdd(name, f1, f2, enter);
}

int
cdoOperatorID(void)
{
  ProcessType &process = processSelf();
  return process.getOperatorID();
}

int
cdoOperatorF1(int operID)
{
  return processSelf().oper[operID].f1;
}

int
cdoOperatorF2(int operID)
{
  return processSelf().oper[operID].f2;
}

const char *
cdoOperatorName(int operID)
{
  return processSelf().oper[operID].name;
}

const char *
cdoOperatorEnter(int operID)
{
  return processSelf().oper[operID].enter;
}

int
cdoStreamNumber()
{
  return operatorStreamNumber(processSelf().operatorName);
}
int
cdoStreamCnt(void)
{
  int cnt = processSelf().m_streamCnt;
  return cnt;
}

int
cdoStreamName(int cnt)
{
  return cnt;
}

ProcessType *
processCreate(const char *command)
{
#ifdef HAVE_LIBPTHREAD
  pthread_mutex_lock(&processMutex);
#endif

  Cdo_Debug(CdoDebug::PROCESS, "Creating new process for command: ", command);
  int processID = NumProcess++;

  const char *operatorName = get_original(getOperatorName(command));
  auto success = Process.insert(std::make_pair(processID, ProcessType(processID, operatorName, command)));
  if (success.second == false)
    {
      ERROR("Process ", processID, " could not be created");
    }

  NumProcessActive++;
  Cdo_Debug(CdoDebug::PROCESS, "NumProcessActive: ", NumProcessActive);

#ifdef HAVE_LIBPTHREAD
  pthread_mutex_unlock(&processMutex);
#endif

  if (processID >= MAX_PROCESS) ERROR("Limit of ", MAX_PROCESS, " processes reached!");

  return &success.first->second;
}

int
processNums(void)
{
#ifdef HAVE_LIBPTHREAD
  pthread_mutex_lock(&processMutex);
#endif

  int pnums = Process.size();

#ifdef HAVE_LIBPTHREAD
  pthread_mutex_unlock(&processMutex);
#endif

  return pnums;
}

int
processNumsActive(void)
{
#ifdef HAVE_LIBPTHREAD
  pthread_mutex_lock(&processMutex);
#endif

  int pnums = NumProcessActive;

#ifdef HAVE_LIBPTHREAD
  pthread_mutex_unlock(&processMutex);
#endif

  return pnums;
}

/**
 * Checks if \p p_argvEntry is a name of a existing file.
 * If no such file exists \p p_argvEntry is added to obase.
 * Else Cdo exits with error message
 * @return returns a pointer to the newly generated process
 */
void
handleObase(const char *p_argvEntry)
{
  if (!fileExists(p_argvEntry))
    {
      obase.insert({ 0, strdup(p_argvEntry) });
    }
  else
    {
      CdoError::Abort(Cdo::progname, "Obase missing. Found existing file: ", p_argvEntry, "instead");
    }
}

/**
 * Handles process Creation for cdo command p_argvEntry
 * Adds \p p_parentProcess as parent to the newly created process.
 * Adds newly created process to \p p_parentProcess.
 * Also checks if \p p_parentProcess accepts another processes streams
 * as input and exits with error message if not.
 */
ProcessType *
createNewProcess(ProcessType *p_parentProces, const char *argvEntry)
{
  ProcessType *newProcess = processCreate(argvEntry);
  if (newProcess->m_module.streamOutCnt == 0)
    {
      CdoError::Abort(Cdo::progname, "operator -", p_parentProces->operatorName, " can not take -", newProcess->operatorName,
                      "  with 0 outputs as input");
      exit(EXIT_FAILURE);
    }

  p_parentProces->addChild(newProcess);
  newProcess->addParent(p_parentProces);
  return newProcess;
}

void
createProcesses(int argc, const char **argv)
{
  Cdo_Debug(CdoDebug::PROCESS, "== Process Creation Start ==");
  Cdo_Debug(CdoDebug::PROCESS, "operators:  ", CdoDebug::argvToString(argc, argv));
  ProcessType *root_process = processCreate(argv[0]);

  ProcessType *currentProcess;
  ProcessType *lastAdded;

  int idx = 1;
  std::stack<ProcessType *> call_stack;

  call_stack.push(root_process);
  currentProcess = call_stack.top();
  int cntOutFiles = (int) currentProcess->m_module.streamOutCnt;

  int temp_argc = argc - cntOutFiles;

  for (int i = 0; i < cntOutFiles; i++)
    {
      Cdo_Debug(CdoDebug::PROCESS, "Creating new pstream for output file: ", argv[temp_argc + i]);
      root_process->addFileOutStream(argv[temp_argc + i]);
    }

  if (cntOutFiles == -1)
    {
      handleObase(argv[argc - 1]);
      cntOutFiles = 1;
    }

  lastAdded = root_process;
  if (idx < argc - cntOutFiles)
    {
      const char *argvEntry;
      do
        {
          argvEntry = argv[idx];
          Cdo_Debug(CdoDebug::PROCESS, "iteration ", idx, ", current argv: ", argv[idx],
                    ",  currentProcess: ", currentProcess->operatorName);
          if (argvEntry[0] == '-')
            {
              Cdo_Debug(CdoDebug::PROCESS, "Found new Operator: ", argvEntry);
              currentProcess = createNewProcess(currentProcess, argvEntry);
              lastAdded = currentProcess;
              call_stack.push(currentProcess);
            }
          else if (currentProcess->m_module.streamInCnt != 0)
            {
              Cdo_Debug(CdoDebug::PROCESS, "adding in file to ", currentProcess->operatorName);
              currentProcess->addFileInStream(argvEntry);
            }
          while (call_stack.top()->hasAllInputs() && call_stack.top() != root_process)
            {
              Cdo_Debug(CdoDebug::PROCESS, "Removing ", call_stack.top()->operatorName, " from stack");
              call_stack.pop();
            }
          currentProcess = call_stack.top();
          idx++;
        }
      while (!currentProcess->hasAllInputs() && idx < argc - cntOutFiles);
    }

  if (idx != argc - cntOutFiles)
    {
      CdoError::Abort(Cdo::progname, " Too many input streams for operator '", lastAdded->operatorName, "'!");
    }

  while (!call_stack.empty())
    {
      call_stack.pop();
    }
  Cdo_Debug(CdoDebug::PROCESS, "== Process Creation End ==");

  setProcessNum(Process.size());
}

void
cdoValidateProcesses()
{
  for (auto &process : Process)
    {
      process.second.checkStreamCnt();
    }
}

void
cdoStreamClose(int p_pstreamIDX)
{
  ProcessType &process = processSelf();
  if (p_pstreamIDX >= process.getInStreamCnt())
    {
      process.outputStreams[p_pstreamIDX - process.getInStreamCnt()]->close();
    }
  else if (p_pstreamIDX < process.getInStreamCnt() + process.getOutStreamCnt())
    {
      PstreamType *pstream = process.inputStreams[p_pstreamIDX];
      pstream->close();
      process.addNvals(pstream->getNvals());
    }
}

void
clearProcesses()
{
  Process.clear();
  NumProcess = 0;
  NumProcessActive = 0;
}

int
cdoStreamOpenRead(int inStreamIDX)
{
  Cdo_Debug(CdoDebug::PROCESS, "Getting in stream ", inStreamIDX, " of process ", processSelf().m_ID);
  ProcessType &process = processSelf();
  if (process.getInStreamCnt() < inStreamIDX || inStreamIDX < 0)
    {
      ERROR("instream ", inStreamIDX, " of process ", process.m_ID, " not found");
    }
  PstreamType *inStream = process.inputStreams[inStreamIDX];

  if (inStream->isPipe())
    {
      Cdo_Debug(CdoDebug::PROCESS, "Trying to open pipe: ", inStream->pipe->name);
      inStream->pstreamOpenReadPipe();
      process.nChildActive++;
    }
  else
    {
      Cdo_Debug(CdoDebug::PROCESS, "Trying to open file: ", inStream->m_mfnames[0]);
      inStream->pstreamOpenReadFile(inStream->m_mfnames[0].c_str());
    }

  return inStream->self;  // return ID
}

int
cdoStreamOpenWrite(int p_outStreamIDX, int filetype)
{
  Cdo_Debug(CdoDebug::PROCESS, "Getting out stream ", p_outStreamIDX, " of process ", processSelf().m_ID);

  ProcessType &process = processSelf();
  int outStreamIDX = p_outStreamIDX - process.inputStreams.size();
  if (outStreamIDX > process.getOutStreamCnt() || outStreamIDX < 0)
    {
      ERROR("outstream ", outStreamIDX, " of ", process.m_ID, " not found.", "Was called with streamIDX = ", p_outStreamIDX);
    }
  PstreamType *outStream = process.outputStreams[outStreamIDX];

  if (outStream->ispipe)
    {
      outStream->pstreamOpenWritePipe(outStream->pipe->name.c_str(), filetype);
    }
  else
    {
      Cdo_Debug(CdoDebug::PROCESS, "Trying to open: ", outStream->m_mfnames[0]);

      if (Options::cdoInteractive)
        {
          struct stat stbuf;

          int rstatus = stat(outStream->m_name.c_str(), &stbuf);
          /* If permanent file already exists, query user whether to overwrite or exit */
          if (rstatus != -1) process.query_user_exit(outStream->m_name.c_str());
        }

      outStream->pstreamOpenWriteFile(filetype);
    }
  return outStream->self;
}

int
cdoStreamOpenWrite(std::string p_filename, int filetype)
{
  int pstreamID = -1;
  ProcessType &process = processSelf();
  process.addFileOutStream(p_filename);
  pstreamID = process.outputStreams.back()->pstreamOpenWriteFile(filetype);

  if (pstreamID == -1)
    {
      ERROR("Could not create pstream for file: ", p_filename);
    }

  return pstreamID;
}

bool
cdoInFileExists(int inStreamIDX)
{
  PstreamType *inStream = processSelf().inputStreams[inStreamIDX];
  return fileExists(inStream->m_mfnames[0].c_str());
}
bool
cdoOutFileExists(int outStreamIDX)
{
  PstreamType *outStream = processSelf().outputStreams[outStreamIDX];
  return fileExists(outStream->m_mfnames[0].c_str());
}
int
cdoStreamOpenAppend(int p_outFileIndex)
{
  ProcessType &process = processSelf();
  int streamIndex = p_outFileIndex - process.inputStreams.size();
  PstreamType *outStream = process.outputStreams[streamIndex];
  int pstreamID = -1;
  if (outStream->ispipe)
    {
      Cdo_Debug(CdoDebug::PROCESS, "pipe ", outStream->pipe->name.c_str());
      cdoAbort("This operator doesn't work with pipes!");
    }
  else
    {
      outStream->openAppend(outStream->m_mfnames[0].c_str());
      pstreamID = outStream->self;
    }
  if (pstreamID == -1)
    {
      ERROR("could not append to ", outStream->m_mfnames[0]);
    }
  return pstreamID;
}

ProcessType *
getProcess(int p_processID)
{
  auto process = Process.find(p_processID);
  if (process == Process.end())
    {
      ERROR("Process with ID: ", p_processID, "not found");
    }
  return &(process->second);
}

std::string
cdoGetOutStreamName(int p_outStream)
{
  ProcessType &process = processSelf();
  return process.outputStreams[p_outStream]->m_name;
}

std::string
cdoGetInStreamName(int p_inStream)
{
  ProcessType &process = processSelf();
  return process.inputStreams[p_inStream]->m_name;
}
std::string
cdoGetStreamName(int p_streamIndex)
{
  std::string streamName;
  ProcessType &process = processSelf();
  Cdo_Debug(CdoDebug::PROCESS, "stridx ", p_streamIndex);
  if (p_streamIndex >= static_cast<int>(process.inputStreams.size()))
    {
      Cdo_Debug(CdoDebug::PROCESS, "Getting output stream name", p_streamIndex);
      streamName = cdoGetOutStreamName(p_streamIndex - process.inputStreams.size());
    }
  else
    {
      Cdo_Debug(CdoDebug::PROCESS, "Getting input stream name", p_streamIndex);
      streamName = cdoGetInStreamName(p_streamIndex);
    }
  Cdo_Debug(CdoDebug::PROCESS, "StreamName is:", streamName);
  return streamName;
}

bool
cdoStreamIsPipe(int p_streamIndex)
{
  ProcessType &process = processSelf();
  if (p_streamIndex >= static_cast<int>(process.inputStreams.size()))
    {
      return process.outputStreams[p_streamIndex - process.inputStreams.size()]->isPipe();
    }
  else
    {
      return process.inputStreams[p_streamIndex]->isPipe();
    }
}
char *
cdoGetObase()
{
  ProcessType &process = processSelf();
  if (obase.find(process.m_ID) == obase.end())
    {
      ERROR("No obase found, please check the module if this operator is defined for obase usage");
    }

  return obase[process.m_ID];
}

void
cdoInitialize(void *p_process)
{
#if defined(_OPENMP)
  omp_set_num_threads(Threading::ompNumThreads);  // Has to be called for every module (pthread)!
#endif
  ProcessType *process = (ProcessType *) p_process;

  // std::cout << arg->processID << std::endl;
  Cdo_Debug(CdoDebug::PROCESS, "Initializing process: ", process->m_operatorCommand);
  process->threadID = pthread_self();
// std::cout << "SomeMarker" << Process.size() << std::endl;

#if defined(HAVE_LIBPTHREAD)
  if (CdoDebug::PSTREAM) Cdo_Debug(CdoDebug::PROCESS, "process ", processSelf().m_ID, " thread ", pthread_self());
#endif
}

void
cdoFinish(void)
{
  ProcessType &process = processSelf();

  Cdo_Debug(CdoDebug::PROCESS, "Finishing process: ", process.m_ID);

#ifdef HAVE_LIBPTHREAD
  Cdo_Debug(CdoDebug::PROCESS, "process ", process.m_ID, " thread ", pthread_self());
#endif
  if (!Options::silentMode) process.printProcessedValues();

  cdoTimes times = process.getTimes(processNums());

  processAccuTime(times.c_usertime, times.c_systime);

  process.setInactive();
  NumProcessActive--;

  char memstring[32] = { "" };

  if (process.m_ID == 0)
    {
      printEndTimes(memstring);
    }
  process.printBenchmarks(times, memstring);
}

const char *
processInqPrompt(void)
{

  ProcessType &process = processSelf();
  return process.inqPrompt();
}

extern "C" {
size_t getPeakRSS();
}

void
printEndTimes(char *p_memstring)
{
  size_t memmax = getPeakRSS();
  if (memmax)
    {
      size_t muindex = 0;
      const char *mu[] = { "B", "KB", "MB", "GB", "TB", "PB" };
      const size_t nmu = sizeof(mu) / sizeof(char *);
      while (memmax > 9999 && muindex < nmu - 1)
        {
          memmax /= 1024;
          muindex++;
        }
      snprintf(p_memstring, sizeof(p_memstring), " %zu%s", memmax, mu[muindex]);
    }
}

void
processAccuTime(double utime, double stime)
{
  ProcessType *process = getProcess(0);
  process->a_utime += utime;
  process->a_stime += stime;
}

void
processStartTime(double *utime, double *stime)
{
  // used in: Command.cc, CDItest.cc, process.cc
  ProcessType &process = processSelf();

  *utime = process.s_utime;
  *stime = process.s_stime;
}

int
cdoStreamInqVlist(int pstreamID)
{
  int vlistID = pstreamInqVlist(pstreamID);
  if (vlistNumber(vlistID) == CDI_COMP && cdoStreamNumber() == CDI_REAL)
    cdoAbort("Fields with complex numbers are not supported by this operator!");

  if (vlistNumber(vlistID) == CDI_REAL && cdoStreamNumber() == CDI_COMP)
    cdoAbort("This operator needs fields with complex numbers!");
  processDefVarNum(vlistNvars(vlistID));
  return vlistID;
}

void
runProcesses()
{
  for (auto &idProcessPair : Process)
    {
      if (idProcessPair.first != 0)
        {
          /*TEMP*/
          if (Options::silentMode == 0)
            {
              set_text_color(stderr, RESET, GREEN);
              fprintf(stderr, "%s: ", idProcessPair.second.prompt);
              reset_text_color(stderr);
              std::cerr << "Process started" << std::endl;
            }
          threadIDs.push_back(idProcessPair.second.run());
        }
    }
  threadIDs.push_back(pthread_self());
  getProcess(0)->m_module.func(getProcess(0));
}

void
killProcesses()
{
  for (auto threadID : threadIDs)
    {
      if (threadID != pthread_self())
        {
          pthread_cancel(threadID);
        }
    }
}

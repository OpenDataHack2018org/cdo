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

#ifdef  HAVE_CONFIG_H
#include "config.h"
#endif

#if defined(HAVE_PTHREAD_H)
#include <pthread.h>
#endif

#include <stdio.h>
#include <string.h>

#if defined(HAVE_GLOB_H)
#include <glob.h>
#endif
#if defined(HAVE_WORDEXP_H)
#include <wordexp.h>
#endif

#include "cdo.h"
#include "cdo_int.h"
#include "error.h"
#include "modules.h"
#include "util.h"
#include "pstream.h"
#include "dmemory.h"
#include "pthread.h"
#include "cdoDebugOutput.h"

#include <algorithm>
#include <map>
#include <stack>

#ifdef  HAVE_LIBPTHREAD
pthread_mutex_t processMutex = PTHREAD_MUTEX_INITIALIZER;
#endif

static process_t *root_process;
std::map<int, process_t> Process;

static int NumProcess = 0;
static int NumProcessActive = 0;

process_t::process_t(int p_ID, const char *operatorCommand) : m_ID(p_ID)
{
  initProcess();
  operatorName = getOperatorName(operatorCommand);
  setOperatorArgv(operatorCommand);
  m_operatorCommand = operatorCommand;

  defPrompt();  // has to be called after get operatorName

  m_module = getModule(operatorName);
}

void
process_t::setOperatorArgv(const char *operatorArguments)
{
  if (operatorArguments)
    {
      char *operatorArg = strdup(operatorArguments);
      // fprintf(stderr, "processDefArgument: %d %s\n", oargc, operatorArg);

      while ((operatorArg = strchr(operatorArg, ',')) != NULL)
        {
          *operatorArg = '\0';
          *operatorArg++;
          if (strlen(operatorArg))
            {
              oargv.push_back(operatorArg);
            }
        }
    }
  oargc = oargv.size();
}

void
process_t::initProcess()
{
#ifdef  HAVE_LIBPTHREAD
  threadID = pthread_self();
  l_threadID = 1;
#endif
  nchild = 0;

  cdoProcessTime(&s_utime, &s_stime);
  a_utime = 0;
  a_stime = 0;
  cputime = 0;
  nvals = 0;
  nvars = 0;
  ntimesteps = 0;

  m_streamCnt = 0;

  oargc = 0;
  m_operatorCommand = "UNINITALIZED";
  operatorArg = "UNINITALIZED";

  noper = 0;
}

int
process_t::getInStreamCnt()
{
  return inputStreams.size();
}
int
process_t::getOutStreamCnt()
{
  return outputStreams.size();
}

process_t *
processCreate(const char *command)
{
#ifdef  HAVE_LIBPTHREAD
  pthread_mutex_lock(&processMutex);
#endif

  if(CdoDebug::PROCESS){
    MESSAGE("Creating new process for command: ", command);
  }
  int processID = NumProcess++;
  auto success = Process.insert(std::make_pair(processID, process_t(processID, command)));
  if(success.second == false)
  {
    ERROR("Process ", processID," could not be created");
  }

  NumProcessActive++;
  if(CdoDebug::PROCESS){
    MESSAGE("NumProcessActive: ", NumProcessActive);
  }

#ifdef  HAVE_LIBPTHREAD
  pthread_mutex_unlock(&processMutex);
#endif

  if (processID >= MAX_PROCESS)
    Error("Limit of %d processes reached!", MAX_PROCESS);

  return &success.first->second;
}

process_t &
processSelf(void)
{
#ifdef  HAVE_LIBPTHREAD
  pthread_t thID = pthread_self();

  pthread_mutex_lock(&processMutex);

  for (auto &id_process_pair : Process)
    if (id_process_pair.second.l_threadID)
      {
        if (pthread_equal(id_process_pair.second.threadID, thID))
          {
            pthread_mutex_unlock(&processMutex);
            return id_process_pair.second;
          }
      }

  pthread_mutex_unlock(&processMutex);

#endif
  return Process.find(0)->second;
}

int
processNums(void)
{
#ifdef  HAVE_LIBPTHREAD
  pthread_mutex_lock(&processMutex);
#endif

  int pnums = Process.size();

#ifdef  HAVE_LIBPTHREAD
  pthread_mutex_unlock(&processMutex);
#endif

  return pnums;
}

int
processNumsActive(void)
{
#ifdef  HAVE_LIBPTHREAD
  pthread_mutex_lock(&processMutex);
#endif

  int pnums = NumProcessActive;

#ifdef  HAVE_LIBPTHREAD
  pthread_mutex_unlock(&processMutex);
#endif

  return pnums;
}

void
processAddNvals(size_t nvals)
{
  processSelf().nvals += nvals;
}

size_t
processInqNvals(int processID)
{
  return Process.find(processID)->second.nvals;
}

void
processAddOutputStream(pstream_t *p_pstream_ptr)
{
  process_t &process = processSelf();

  int sindex = process.getOutStreamCnt();

  if (sindex >= MAX_STREAM)
    Error("Limit of %d output streams per process reached (processID = %d)!", MAX_STREAM, process.m_ID);

  process.outputStreams.push_back(p_pstream_ptr);
}

void
processAddInputStream(pstream_t *p_pstream_ptr)
{
  process_t &process = processSelf();

  if (p_pstream_ptr->isPipe())
    {
      process.nchild++;
    }
  else
  {
    }
  int sindex = process.getInStreamCnt();

  if (sindex >= MAX_STREAM)
    Error("Limit of %d input streams per process reached (processID = %d)!", MAX_STREAM, process.m_ID);

  process.inputStreams.push_back(p_pstream_ptr);
}

void
processDelStream(int streamID)
{
  UNUSED(streamID);
}

void
processDefCputime(int processID, double cputime)
{
  Process.find(processID)->second.cputime = cputime;
}

double
processInqCputime(int processID)
{
  return Process.find(processID)->second.cputime;
}

void
processStartTime(double *utime, double *stime)
{
  process_t &process = processSelf();

  *utime = process.s_utime;
  *stime = process.s_stime;
}

void
processEndTime(double *utime, double *stime)
{
  process_t &process = Process.find(0)->second;
  *utime = process.a_utime;
  *stime = process.a_stime;
}

void
processAccuTime(double utime, double stime)
{
  process_t &process = Process.find(0)->second;
  process.a_utime += utime;
  process.a_stime += stime;
}

int
processInqOutputStreamNum(void)
{
  return processSelf().getOutStreamCnt();
}

int
processInqInputStreamNum(void)
{
  return processSelf().getInStreamCnt();
}

int
processInqChildNum(void)
{
  return processSelf().nchild;
}

pstream_t *
processInqOutputStream(int streamindex)
{
  return (processSelf().outputStreams[streamindex]);
}
pstream_t *
processInqInputStream(int streamindex)
{
  return (processSelf().inputStreams[streamindex]);
}

const char *
processInqOpername2(process_t &process)
{
  return process.operatorName;
}

const char *
processInqOpername(void)
{
  return processSelf().operatorName;
}

void
process_t::defPrompt()
{
  if (m_ID == 0)
    sprintf(prompt, "%s %s", CDO_progname, operatorName);
  else
    sprintf(prompt, "%s(%d) %s", CDO_progname, m_ID + 1, operatorName);
}

const char *
processInqPrompt(void)
{
  process_t &process = processSelf();

  const char *prompt = "cdo";
  if (process.prompt[0])
    prompt = process.prompt;

  return prompt;
}

#if defined(HAVE_GLOB_H)
static int
get_glob_flags(void)
{
  int glob_flags = 0;

#if defined(GLOB_NOCHECK)
  glob_flags |= GLOB_NOCHECK;
#endif
#if defined(GLOB_TILDE)
  glob_flags |= GLOB_TILDE;
#endif

  return glob_flags;
}
#endif

#if defined(HAVE_WORDEXP_H)
/* Convert a shell pattern into a list of filenames. */
static argument_t *
glob_pattern(const char *restrict string)
{
  size_t cnt, length = 0;
  int flags = WRDE_UNDEF;
  char **p;

  wordexp_t glob_results;
  glob_results.we_wordc = 0;
  argument_t *argument = NULL;

  // glob the input argument or do even more shell magic
  int status = wordexp(string, &glob_results, flags);

  // How much space do we need?
  for (p = glob_results.we_wordv, cnt = glob_results.we_wordc; cnt; p++, cnt--)
    {
      length += strlen(*p) + 1;
    }

  // Allocate the space and generate the list.
  argument = argument_new(glob_results.we_wordc, length);

  // put all generated filenames into the argument_t data structure
  for (cnt = 0; cnt < glob_results.we_wordc; cnt++)
    {
      argument->argv[cnt] = strdupx(glob_results.we_wordv[cnt]);
      strcat(argument->args, glob_results.we_wordv[cnt]);
      if (cnt < glob_results.we_wordc - 1)
        strcat(argument->args, " ");
    }

  if (status == 0)
    wordfree(&glob_results);

  return argument;
}
#endif

int
cdoStreamCnt(void)
{
  int cnt = processSelf().m_streamCnt;
  return cnt;
}

const argument_t *
cdoStreamName(int cnt)
{
  process_t &process = processSelf();

  if (cnt > process.m_streamCnt || cnt < 0)
    Error("count %d out of range!", cnt);

  return &(process.streamArguments[cnt]);
}

const char *
processOperator(void)
{
  return processSelf().m_operatorCommand;
}

static int skipInputStreams(int argc, std::vector<char *> &argv, int globArgc, int nstreams);

static int
getGlobArgc(int argc, std::vector<char *> &argv, int globArgc)
{
  char *opername = &argv[globArgc][1];
  char *comma_position = strchr(opername, ',');
  if (comma_position)
    *comma_position = 0;

  int streamInCnt = operatorStreamInCnt(opername);
  int streamOutCnt = operatorStreamOutCnt(opername);

  if (streamInCnt == -1)
    streamInCnt = 1;

  if (streamOutCnt > 1)
    cdoAbort("More than one output stream not allowed in CDO pipes (Operator %s)!", opername);

  globArgc++;

  if (streamInCnt > 0)
    globArgc = skipInputStreams(argc, argv, globArgc, streamInCnt);
  if (comma_position)
    *comma_position = ',';

  return globArgc;
}

static int
skipInputStreams(int argc, std::vector<char *> &argv, int globArgc, int nstreams)
{
  while (nstreams > 0)
    {
      if (globArgc >= argc)
        {
          cdoAbort("Too few arguments. Check command line!");
          break;
        }
      if (argv[globArgc][0] == '-')
        {
          globArgc = getGlobArgc(argc, argv, globArgc);
        }
      else
        globArgc++;

      nstreams--;
    }

  return globArgc;
}

static int
getStreamCnt(int argc, std::vector<char *> &argv)
{
  int streamCnt = 0;
  int globArgc = 1;

  while (globArgc < argc)
    {
      if (argv[globArgc][0] == '-')
        {
          globArgc = getGlobArgc(argc, argv, globArgc);
        }
      else
        globArgc++;

      streamCnt++;
    }

  return streamCnt;
}

/*
static void setStreamNames(int argc, std::vector<char *> *argv)
{
    //check for output of first
    int current_argv_entry;
    std::vector<pstream_t*> &current_instreams = root_process->inputStreams;
    std::vector<pstream_t*> &current_outstreams = root_process->outputStreams;

}
*/
void
process_t::setStreams(int argc, std::vector<char *> &argv)
{
  int streamCnt = getStreamCnt(argc, argv);

  nvals = 0;
  nvars = 0;
  ntimesteps = 0;

  m_streamCnt = 0; /* filled in setStreamNames */
  if (streamCnt)
    streamArguments = std::vector<argument_t>(streamCnt);
  for (int i = 0; i < streamCnt; i++)
    {
      streamArguments[i].argc = 0;
      streamArguments[i].args = NULL;
    }

  setStreamNames(argc, argv);

  int status = checkStreamCnt();

  if (status == 0 && streamCnt != streamCnt)
    Error("Internal problem with stream count %d %d", streamCnt, streamCnt);
  /*
  for ( i = 0; i < streamCnt; i++ )
    fprintf(stderr, "setStreams: stream %d %s\n", i+1, process.streamArguments[i].args);
  */
}

void
process_t::setStreamNames(int argc, std::vector<char *> &argv)
{
  int i, ac;
  int globArgc = 1;
  int globArgcStart;
  char *streamname;
  int len;

  while (globArgc < argc)
    {
      if (argv[globArgc][0] == '-')
        {
          globArgcStart = globArgc;

          globArgc = getGlobArgc(argc, argv, globArgc);
          len = 0;
          for (i = globArgcStart; i < globArgc; i++)
            {
              len += strlen(argv[i]) + 1;
            }
          streamname = (char *) Calloc(1, len);
          for (i = globArgcStart; i < globArgc; i++)
            {
              strcat(streamname, argv[i]);
              if (i < globArgc - 1)
                {
                  strcat(streamname, " ");
                }
            }
          for (i = 1; i < len - 1; i++)
            {
              if (streamname[i] == '\0')
                {
                  streamname[i] = ' ';
                }
            }

          streamArguments[m_streamCnt].args = streamname;
          ac = globArgc - globArgcStart;
          // printf("setStreamNames:  ac %d  streamname1: %s\n", ac, streamname);
          streamArguments[m_streamCnt].argv.resize(ac);
          for (i = 0; i < ac; ++i)
            streamArguments[m_streamCnt].argv[i] = argv[i + globArgcStart];
          streamArguments[m_streamCnt].argc = ac;
          m_streamCnt++;
          // printf("setStreamNames:  streamname1: %s\n", streamname);
        }
      else
        {
          len = strlen(argv[globArgc]) + 1;
          streamname = (char *) Malloc(len);
          strcpy(streamname, argv[globArgc]);
          streamArguments[m_streamCnt].args = streamname;
          ac = 1;
          streamArguments[m_streamCnt].argv.resize(ac);
          streamArguments[m_streamCnt].argv[0] = argv[globArgc];
          streamArguments[m_streamCnt].argc = ac;
          streamArguments[m_streamCnt].args = streamname;
          m_streamCnt++;
          // printf("setStreamNames:  streamname2: %s\n", streamname);
          globArgc++;
        }
    }
}

static int
find_wildcard(const char *string, size_t len)
{
  int status = 0;

  if (len > 0)
    {
      if (string[0] == '~')
        status = 1;

      if (status == 0)
        {
          for (size_t i = 0; i < len; ++i)
            if (string[i] == '?' || string[i] == '*' || string[i] == '[')
              {
                status = 1;
                break;
              }
        }
    }

  return status;
}

char *
expand_filename(const char *string)
{
  char *filename = NULL;

  if (find_wildcard(string, strlen(string)))
    {
#if defined(HAVE_GLOB_H)
      int glob_flags = get_glob_flags();
      glob_t glob_results;

      glob(string, glob_flags, 0, &glob_results);

      if (glob_results.gl_pathc == 1)
        filename = strdupx(glob_results.gl_pathv[0]);

      globfree(&glob_results);
#endif
    }

  return filename;
}

 int
process_t::expand_wildcards(int streamCnt)
{
  const char *streamname0 = streamArguments[0].args;

  if (streamname0[0] == '-')
    return 1;

#if defined(HAVE_WORDEXP_H)
  argument_t *glob_arg = glob_pattern(streamname0);

  // skip if the input argument starts with an operator (starts with -)
  // otherwise adapt streams if there are several files (>1)
  // in case of one filename skip, no adaption needed
  if (glob_arg->argc > 1 && glob_arg->argv[0][0] != '-')
    {
      if (cdoVerbose)
        cdoPrint("Replaced >%s< by", streamArguments[0].args);

      streamCnt = streamCnt - 1 + glob_arg->argc;

      Free(streamArguments[0].args);

      streamArguments.resize(streamCnt);

      // move output streams to the end
      for (int i = 1; i < m_streamCnt; ++i)
        streamArguments[i + glob_arg->argc - 1] = streamArguments[i];

      for (int i = 0; i < glob_arg->argc; ++i)
        {
          argument_t &current_argument = streamArguments[i];
          current_argument.argc = 1;
          current_argument.argv.resize(current_argument.argc);
          current_argument.argv[0] = strdupx(glob_arg->argv[i]);
          current_argument.args = strdupx(glob_arg->argv[i]);
          if (cdoVerbose)
            cdoPrint("         >%s<", glob_arg->argv[i]);
        }

      m_streamCnt = streamCnt;
    }

  Free(glob_arg);
#endif

  return 1;
}

int process_t::checkStreamCnt(void)
{
  int wantedStreamInCnt,wantedStreamOutCnt;
  int streamInCnt0;
  int streamCnt = 0;
  int i, j;
  int obase = FALSE;
  int status = 0;

  wantedStreamInCnt= operatorStreamInCnt(operatorName);
  wantedStreamOutCnt = operatorStreamOutCnt(operatorName);

  streamInCnt0 = wantedStreamInCnt;

  if (wantedStreamOutCnt == -1)
    {
      wantedStreamOutCnt = 1;
      obase = TRUE;
    }

  if (wantedStreamInCnt== -1 && wantedStreamOutCnt == -1)
    cdoAbort("I/O stream counts unlimited no allowed!");

  // printf(" wantedStreamInCnt,wantedStreamOutCnt %d %d\n", wantedStreamInCnt,wantedStreamOutCnt);
  if (wantedStreamInCnt== -1)
    {
      wantedStreamInCnt = m_streamCnt - wantedStreamOutCnt;
      if (wantedStreamInCnt< 1)
        cdoAbort("Input streams missing!");
    }

  if (wantedStreamOutCnt == -1)
    {
      wantedStreamOutCnt = m_streamCnt - wantedStreamInCnt;
      if (wantedStreamOutCnt < 1)
        cdoAbort("Output streams missing!");
    }
  // printf(" wantedStreamInCnt,wantedStreamOutCnt %d %d\n", wantedStreamInCnt,wantedStreamOutCnt);

  streamCnt = wantedStreamInCnt+ wantedStreamOutCnt;
  // printf(" streamCnt %d %d\n", m_streamCnt, streamCnt);

  if (m_streamCnt > streamCnt)
    cdoAbort("Too many streams!"
             " Operator needs %d input and %d output streams.",
             wantedStreamInCnt,
             wantedStreamOutCnt);

  if (m_streamCnt < streamCnt)
    cdoAbort("Too few streams specified!"
             " Operator needs %d input and %d output streams.",
             wantedStreamInCnt,
             wantedStreamOutCnt);

  for (i = wantedStreamInCnt;i < streamCnt; i++)
    {
      if (streamArguments[i].args[0] == '-')
        {
          cdoAbort("Output file name %s must not begin with \"-\"!", streamArguments[i].args);
        }
      else if (!obase)
        {
          for (j = 0; j < wantedStreamInCnt;j++) /* does not work with files in pipes */
            if (strcmp(streamArguments[i].args, streamArguments[j].args) == 0)
              cdoAbort("Output file name %s is equal to input file name"
                       " on position %d!\n",
                       streamArguments[i].args,
                       j + 1);
        }
    }

  if (wantedStreamInCnt == 1 && streamInCnt0 == -1)
    status = expand_wildcards( streamCnt);

  return status;
}

bool
process_t::hasAllInputs()
{
  if(m_module.streamInCnt == -1)
  {
      return false;
  }
  return  m_module.streamInCnt == (inputStreams.size());
}

#include <fstream>
void print_creation_results(std::ofstream &p_outfile)
{
 p_outfile << std::endl << "RESULTS:" << std::endl;
  for (auto &process : Process)
    {
      p_outfile << "process: " << process.second.operatorName << " has children: " << std::endl;
      for (auto child : process.second.childProcesses)
        {
          p_outfile << child->m_ID << ", ";
        }
      for (auto outstream : process.second.inputStreams)
        {
          p_outfile << "S: " << outstream->self << " ";
        }
    }
  p_outfile << std::endl;

}

void
createProcesses(int argc, const char **argv)
{
  if(CdoDebug::PROCESS){
  std::string input_string = "";
     MESSAGE("== Process Creation Start ==");
    MESSAGE("operators:  ",CdoDebug::argvToString(argc, argv));
  }
  root_process = processCreate(argv[0]);

  process_t *current_process;
  process_t *parent_process;

  int idx = 1;
  std::stack<process_t *> call_stack;

  call_stack.push(root_process);
  current_process = call_stack.top();
  int cntOutFiles = std::max(0, (int)current_process->m_module.streamOutCnt);
  for(int i = 0; i < cntOutFiles; i++)
  {
    if(CdoDebug::PROCESS)
    {
        MESSAGE("Creating new pstream for output file: ", argv[argc - (i + 1)]);
    }
    root_process->outputStreams.push_back(create_pstream(std::string(argv[argc - (i + 1)])));
  }
  do
    {
      if(CdoDebug::PROCESS){
      MESSAGE(
              "iteration " , idx , ", current argv: " , argv[idx] ,
              ",  current_process: " , current_process->operatorName
              );
      }
    if (argv[idx][0] == '-')
        {
          if(CdoDebug::PROCESS){
            MESSAGE("Found new Operator: ", argv[idx]);
          }
          parent_process = current_process;
          current_process = processCreate(argv[idx]);
          parent_process->addChild(current_process);
          current_process->addParent(parent_process);
          call_stack.push(current_process);
        }
      else
        {
          if(CdoDebug::PROCESS){
            MESSAGE("adding file to ", current_process->operatorName);
          }
          MESSAGE(argv[idx]);
          current_process->inputStreams.push_back(create_pstream(argv[idx]));
        }

    while (current_process->hasAllInputs() && current_process != root_process)
      {
        if(CdoDebug::PROCESS)
        {
            MESSAGE("Removing ", current_process->operatorName, " from stack");
        }
        call_stack.pop();
        current_process = call_stack.top();
      }
      idx++;
    }
  while ((current_process != root_process || !root_process->hasAllInputs()) && idx < argc - cntOutFiles);

  if(CdoDebug::PROCESS)
  {
    MESSAGE("== Process Creation End ==");
  }
}

void
processDefArgument(void *vargument)
{
  /*
process_t &process = processSelf();
char *operatorArg;
char *commapos;
std::vector< char*> &oargv = process.oargv;
*/
  process_t &process = processSelf();
  std::vector<char*> &oargv = process.oargv;
  /*

  process.m_operatorCommand = argv[0];
  process.operatorName = getOperatorName(process.m_operatorCommand);
  process.operatorArg = getOperatorArg(process.m_operatorCommand);
  operatorArg = process.operatorArg;

  if (operatorArg)
    {
      orgv.push_back(operatorArg);
      // fprintf(stderr, "processDefArgument: %d %s\n", oargc, operatorArg);

      char *commapos = operatorArg;
      while ((commapos = strchr(commapos, ',')) != NULL)
        {
          *commapos = '\0';
          commapos++;
          if (strlen(commapos))
            {
              oargv.push_back(commapos);
            }
        }
      process.oargc = oargv.size();
    }

  processDefPrompt(process.operatorName);

*/
}

void
processDefVarNum(int nvars)
{
  process_t &process = processSelf();
  /*  if ( streamID == process.streams[0] ) */
  process.nvars += nvars;
}

int
processInqVarNum(void)
{
  return processSelf().nvars;
}

void
processDefTimesteps(int streamID)
{
  process_t &process = processSelf();

  UNUSED(streamID);
  /*
  int i;
  printf("streamID %d %d %d %d\n", streamID, Process[processID].streams[0], Process[processID].streams[1], processID);

  for ( i = 0; i < Process[processID].nstream; i++)
    printf("streamID %d %d %d %d << \n", processID, Process[processID].nstream, i, Process[processID].streams[i]);
  */
  /*  if ( streamID == Process[processID].streams[0] )*/
  process.ntimesteps++;
}

int
processInqTimesteps(void)
{
  return processSelf().ntimesteps;
}

void
processDelete(void)
{
  process_t &process = processSelf();

// fprintf(stderr, "delete processID %d\n", processID);
#ifdef  HAVE_LIBPTHREAD
  pthread_mutex_lock(&processMutex);

  process.l_threadID = 0;
#endif
  NumProcessActive--;

#ifdef  HAVE_LIBPTHREAD
  pthread_mutex_unlock(&processMutex);
#endif
}

int
operatorArgc(void)
{
  return processSelf().oargc;
}

char **
operatorArgv(void)
{
  return &processSelf().oargv[0];
}

void
operatorCheckArgc(int numargs)
{
  int argc = processSelf().oargc;

  if (argc < numargs)
    cdoAbort("Too few arguments! Need %d found %d.", numargs, argc);
  else if (argc > numargs)
    cdoAbort("Too many arguments! Need %d found %d.", numargs, argc);
}

void
operatorInputArg(const char *enter)
{
  process_t &process = processSelf();

  int oargc = process.oargc;

  if (oargc == 0)
    {
      char line[1024];
      char *pline = line;
      size_t pos, len, linelen;
      int lreadline = 1;

      if (enter)
        {
          set_text_color(stderr, BRIGHT, MAGENTA);
          fprintf(stderr, "%-16s : ", processInqPrompt());
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
              while (pline[pos] == ' ' || pline[pos] == ',')
                pos++;
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
                  while (pline[len] != ' ' && pline[len] != ',' && pline[len] != '\\' && len < linelen)
                    len++;

                  process.oargv.push_back((char*)Malloc(len + 1));
                  memcpy(process.oargv[oargc], pline, len);
                  process.oargv[oargc][len] = '\0';
                  oargc++;

                  pline += len;
                }
              else
                break;
            }
        }

      process.oargc = oargc;
    }
}

int
cdoOperatorAdd(const char *name, int f1, int f2, const char *enter)
{
  process_t &process = processSelf();
  int operID = process.noper;

  if (operID >= MAX_OPERATOR)
    cdoAbort("Maximum number of %d operators reached!", MAX_OPERATOR);

  process.oper[operID].f1 = f1;
  process.oper[operID].f2 = f2;
  process.oper[operID].name = name;
  process.oper[operID].enter = enter;

  process.noper++;

  return operID;
}

int
cdoOperatorID(void)
{
  process_t &process = processSelf();
  int operID = -1;

  if (process.noper > 0)
    {
      for (operID = 0; operID < process.noper; operID++)
        {
          if (process.oper[operID].name)
            if (strcmp(process.operatorName, process.oper[operID].name) == 0)
              break;
        }
      if (operID == process.noper)
        cdoAbort("Operator not callable by this name!");
    }
  else
    {
      cdoAbort("Operator not initialized!");
    }

  return operID;
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

void
process_t::print_process()
{
#ifdef  HAVE_LIBPTHREAD
  std::cout << " processID       : " << m_ID << std::endl;
  std::cout << " threadID        : " << threadID << std::endl;
  std::cout << " l_threadID      : " << l_threadID << std::endl;
#endif
  std::cout << " nchild          : " << nchild << std::endl;
  int nInStream = getInStreamCnt();
  int nOutStream = getOutStreamCnt();
  std::cout << " nInStream       : " << nInStream << std::endl;
  std::cout << " nOutStream      : " << nOutStream << std::endl;
  for (int i = 0; i < nInStream; i++)
    {
      std::cout << "    " << childProcesses[i]->m_ID << std::endl;
    }
  for (int i = 0; i < nOutStream; i++)
    {
      std::cout << "    " << parentProcesses[i]->m_ID << std::endl;
    }
  if ( s_utime > 0 )
    {
      std::cout << " s_utime         : " << s_utime << std::endl;
    }
  else
    {
      std::cout << " s_utime         : " << "UNINITALIZED" << std::endl;
    }
  if ( s_stime > 0 )
    {
      std::cout << " s_stime         : " << s_stime << std::endl;
    }
  else
    {
      std::cout << " s_stime         : "
                << "UNINITALIZED" << std::endl;
    }

  std::cout << " a_utime         : " << a_utime << std::endl;
  std::cout << " a_stime         : " << a_stime << std::endl;
  std::cout << " cputime         : " << cputime << std::endl;

  if (nvals)
    {
      std::cout << " nvals           : " << nvals << std::endl;
    }
  else
    {

      std::cout << " nvals           : "
                << "UNINITALIZED" << std::endl;
    }
  if (nvars)
    {
      std::cout << " nvars           : " << nvars << std::endl;
    }
  else
    {

      std::cout << " nvars           : "
                << "UNINITALIZED" << std::endl;
    }
  if (ntimesteps)
    {
      std::cout << " ntimesteps      : " << ntimesteps << std::endl;
    }
  else
    {

      std::cout << " ntimesteps      : "
                << "UNINITALIZED" << std::endl;
    }
  std::cout << " ntimesteps      : " << ntimesteps << std::endl;
  std::cout << " streamCnt       : " << m_streamCnt << std::endl;
  // std::cout << " streamArguments     : " << streamArguments                  <<  std::endl;
  std::cout << " m_operatorCommand       : " << m_operatorCommand << std::endl;
  std::cout << " operatorName    : " << operatorName << std::endl;
  std::cout << " operatorArg     : " << operatorArg << std::endl;
  std::cout << " oargc           : " << oargc << std::endl;
  std::cout << " noper           : " << noper << std::endl;
}

static void
processClosePipes(void)
{
  int nstream = processInqInputStreamNum();
  for (int sindex = 0; sindex < nstream; sindex++)
    {
      pstream_t *pstreamptr = processInqInputStream(sindex);

      if (CdoDebug::PROCESS)
        MESSAGE("process ",processSelf().m_ID,"  instream ",sindex,"  close streamID ", pstreamptr->self);

      if (pstreamptr)
        pstreamptr->close();
    }

  nstream = processInqOutputStreamNum();
  for (int sindex = 0; sindex < nstream; sindex++)
    {
      pstream_t *pstreamptr = processInqOutputStream(sindex);

      if (CdoDebug::PROCESS)
        MESSAGE("process ",processSelf().m_ID,"  outstream ",sindex,"  close streamID ", pstreamptr->self);


      if (pstreamptr)
        pstreamptr->close();
    }
}

extern "C" {
size_t getPeakRSS( );
}

void cdoFinish(void)
{
  int processID = processSelf().m_ID;
  int nvars, ntimesteps;
  char memstring[32] = { "" };
  double s_utime, s_stime;
  double e_utime, e_stime;
  double c_cputime = 0, c_usertime = 0, c_systime = 0;
  double p_cputime = 0, p_usertime = 0, p_systime = 0;

#ifdef  HAVE_LIBPTHREAD
  if (CdoDebug::PROCESS)
    MESSAGE("process ",processID," thread %ld", pthread_self());
#endif

  int64_t nvals = processInqNvals(processID);
  nvars = processInqVarNum();
  ntimesteps = processInqTimesteps();

  if (!cdoSilentMode)
    {
      set_text_color(stderr, RESET, GREEN);
      fprintf(stderr, "%s: ", processInqPrompt());
      reset_text_color(stderr);
      if (nvals > 0)
        {
          if (sizeof(int64_t) > sizeof(size_t))
#if defined(_WIN32)
            fprintf(stderr,
                    "Processed %I64d value%s from %d variable%s",
#else
            fprintf(stderr,
                    "Processed %jd value%s from %d variable%s",
#endif
                    (intmax_t) nvals,
                    ADD_PLURAL(nvals),
                    nvars,
                    ADD_PLURAL(nvars));
          else
            fprintf(stderr,
                    "Processed %zu value%s from %d variable%s",
                    (size_t) nvals,
                    ADD_PLURAL(nvals),
                    nvars,
                    ADD_PLURAL(nvars));
        }
      else if (nvars > 0)
        {
          fprintf(stderr, "Processed %d variable%s", nvars, ADD_PLURAL(nvars));
        }

      if (ntimesteps > 0)
        fprintf(stderr, " over %d timestep%s", ntimesteps, ADD_PLURAL(ntimesteps));

      //  fprintf(stderr, ".");
    }
  /*
    fprintf(stderr, "%s: Processed %d variable%s %d timestep%s.",
            processInqPrompt(), nvars, nvars > 1 ? "s" : "",
            ntimesteps, ntimesteps > 1 ? "s" : "");
  */
  processStartTime(&s_utime, &s_stime);
  cdoProcessTime(&e_utime, &e_stime);

  c_usertime = e_utime - s_utime;
  c_systime = e_stime - s_stime;
  c_cputime = c_usertime + c_systime;

#ifdef  HAVE_LIBPTHREAD
  if (getPthreadScope() == PTHREAD_SCOPE_PROCESS)
    {
      c_usertime /= processNums();
      c_systime /= processNums();
      c_cputime /= processNums();
    }
#endif

  processDefCputime(processID, c_cputime);

  processAccuTime(c_usertime, c_systime);

  if (processID == 0)
    {
      size_t memmax = getPeakRSS();
      if (memmax)
        {
          size_t muindex = 0;
          const char *mu[] = { "B", "KB", "MB", "GB", "TB", "PB" };
          const size_t nmu = sizeof(mu)/sizeof(char*);
          while (memmax > 9999 && muindex < nmu-1) { memmax /= 1024; muindex++; }
          snprintf(memstring, sizeof(memstring), " %zu%s", memmax, mu[muindex]);
        }

      processEndTime(&p_usertime, &p_systime);
      p_cputime = p_usertime + p_systime;
    }

#if defined(HAVE_SYS_TIMES_H)
  if (cdoBenchmark)
    fprintf(stderr, " [%.2fs %.2fs %.2fs%s]\n", c_usertime, c_systime, c_cputime, memstring);
  else
    {
      if (!cdoSilentMode)
        fprintf(stderr, " [%.2fs%s]\n", c_cputime, memstring);
    }
  if (cdoBenchmark && processID == 0)
    fprintf(stderr, "total: user %.2fs  sys %.2fs  cpu %.2fs  mem%s\n", p_usertime, p_systime, p_cputime, memstring);
#else
  fprintf(stderr, "\n");
#endif

  processClosePipes();

  processDelete();
}

void
process_t::addChild(process_t *childProcess)
{
  childProcesses.push_back(childProcess);
  nchild = childProcesses.size();
  inputStreams.push_back(create_pstream());
}

void
process_t::addParent(process_t *parentProcess)
{
  parentProcesses.push_back(parentProcess);
}
void
clearProcesses()
{
  Process.clear();
  NumProcess = 0;
  NumProcessActive = 0;
}

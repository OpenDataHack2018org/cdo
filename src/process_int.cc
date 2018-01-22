#include "process.h"
#include "cdoDebugOutput.h"

void
processDefVarNum(int nvars)
{
  process_t &process = processSelf();
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
  process.ntimesteps++;
}

int
processInqTimesteps(void)
{
  return processSelf().ntimesteps;
}

int
operatorArgc(void)
{
  return processSelf().m_oargc;
}

char **
operatorArgv(void)
{
  if(CdoDebug::PROCESS)
  {
      std::string oargv_str = "";
      for( auto entry: processSelf().m_oargv)
      {
          oargv_str += std::string(entry) + " ";
      }
      if(CdoDebug::PROCESS)
      {
        MESSAGE("Getting ",processSelf().m_oargv.size()," operator arguments: ", oargv_str);
      }
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
  process_t &process = processSelf();
  process.inqUserInputForOpArg(enter);
}

int
cdoOperatorAdd(const char *name, int f1, int f2, const char *enter)
{
    process_t &process = processSelf();
    return process.operatorAdd(name, f1, f2, enter);
}

int
cdoOperatorID(void)
{
  process_t &process = processSelf();
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


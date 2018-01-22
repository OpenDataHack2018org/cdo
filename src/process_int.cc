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


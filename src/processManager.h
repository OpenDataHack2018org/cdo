#ifndef PROCESS_MANAGER_H
#define PROCESS_MANAGER_H

// Stdlib includes
#include <map>
#include <vector>

#include <pthread.h>

// cdo includes

// Froward declarations
class ProcessType;

// Error codes
enum class ParseStatus
{
  Ok = 0,
  OpenBracketMissing = -1,
  ClosingBracketMissing = -2,
  UnprocessedInput = -3,
  MissingOutFile = -4,
  MissingObase = -5
};
class ProcessManager
{

private:
  std::map<int, ProcessType> m_processes;
  std::vector<pthread_t> m_threadIDs;

  int m_numProcesses = 0;
  int m_numProcessesActive = 0;

  void handleFirstOperator(int p_argcStart, int argc, const char **argv, ProcessType *p_rootProcess);
  void checkSingleBracketOnly(const char *p_argvEntry, char p_bracketType);
  void handleParseError(ParseStatus p_parseStatus);

  ProcessType *createProcess(const char *p_command);
  ProcessType *addProcess(ProcessType *p_parentProcess, const char *p_command);

public:
  ParseStatus createProcessesFromInput(int argc, const char **argv);
  void createProcesses(int argc, const char **argv);
  void runProcesses();
  void killProcesses();
  void validateProcesses();
  void clearProcesses();
  int processNums();
  int processNumsActive();
  ProcessType &getProcess(int p_processID);
  ProcessType &getProcessByThreadID(pthread_t p_threadID);
};

#endif

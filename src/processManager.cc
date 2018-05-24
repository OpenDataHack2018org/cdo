#include "processManager.h"
#include "process.h"
#include "cdoDebugOutput.h"
#include "util_operatorStrings.h"
#include "cdoOptions.h"
#include "text.h"

#include <stack>
#include <set>

#ifdef HAVE_LIBPTHREAD
pthread_mutex_t processMutex = PTHREAD_MUTEX_INITIALIZER;
#endif

void
ProcessManager::runProcesses()
{
  for (auto &idProcessPair : m_processes)
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
          m_threadIDs.push_back(idProcessPair.second.run());
        }
    }
  m_threadIDs.push_back(pthread_self());
  getProcess(0).m_module.func(&(getProcess(0)));
}

void
ProcessManager::killProcesses()
{
  for (auto threadID : m_threadIDs)
    {
      if (threadID != pthread_self())
        {
          pthread_cancel(threadID);
        }
    }
}
void
ProcessManager::validateProcesses()
{
  for (auto &process : m_processes)
    {
      process.second.validate();
    }
}

void
ProcessManager::clearProcesses()
{
  Cdo_Debug(CdoDebug::PROCESS_MANAGER, "Deleting Processes");
  m_processes.clear();
  m_numProcesses = 0;
  m_numProcessesActive = 0;
}

ProcessType *
ProcessManager::createProcess(const char *command)
{
  Cdo_Debug(CdoDebug::PROCESS, "Creating new process for command: ", command, " with ID: ", m_numProcesses);
  int processID = m_numProcesses++;

  const char *operatorName = get_original(getOperatorName(command));
  auto success = m_processes.insert(std::make_pair(processID, ProcessType(processID, operatorName, command)));
  if (success.second == false)
    {
      ERROR("Process ", processID, " could not be created");
    }

  m_numProcessesActive++;
  Cdo_Debug(CdoDebug::PROCESS, "m_numProcessesActive: ", m_numProcessesActive);

  if (processID >= MAX_PROCESS) ERROR("Limit of ", MAX_PROCESS, " processes reached!");

  return &success.first->second;
}

int
ProcessManager::processNums(void)
{
#ifdef HAVE_LIBPTHREAD
  pthread_mutex_lock(&processMutex);
#endif

  int pnums = m_processes.size();

#ifdef HAVE_LIBPTHREAD
  pthread_mutex_unlock(&processMutex);
#endif

  return pnums;
}

int
ProcessManager::processNumsActive(void)
{
#ifdef HAVE_LIBPTHREAD
  pthread_mutex_lock(&processMutex);
#endif

  int pnums = m_numProcessesActive;

#ifdef HAVE_LIBPTHREAD
  pthread_mutex_unlock(&processMutex);
#endif

  return pnums;
}

/**
 * Handles process Creation for cdo command p_argvEntry
 * Adds \p p_parentProcess as parent to the newly created process.
 * Adds newly created process to \p p_parentProcess.
 * Also checks if \p p_parentProcess accepts another processes streams
 * as input and exits with error message if not.
 */
ProcessType *
ProcessManager::addProcess(ProcessType *p_parentProces, const char *argvEntry)
{
  ProcessType *newProcess = createProcess(argvEntry);
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
ProcessManager::handleFirstOperator(int p_argcStart, int argc, const char **argv, ProcessType *p_rootProcess)
{
  for (int i = p_argcStart; i < argc; i++)
    {
      Cdo_Debug(CdoDebug::PROCESS, "Creating new pstream for output file: ", argv[i]);
      if (strcmp(argv[i], "]") == 0)
        {
          CdoError::Abort(Cdo::progname, "missing output file");
        }
      p_rootProcess->addFileOutStream(argv[i]);
    }
}

void
ProcessManager::checkSingleBracketOnly(const char *p_argvEntry, char p_bracketType)
{
  if (strlen(p_argvEntry) != 1)
    {
      CdoError::Abort(Cdo::progname, "Only single ", p_bracketType, " allowed");
    }
}

void
ProcessManager::createProcesses(int argc, const char **argv)
{
  ParseStatus parseStatus = createProcessesFromInput(argc, argv);
  if (parseStatus != ParseStatus::Ok)
    {
      // Program Exits here
      handleParseError(parseStatus);
    }
  validateProcesses();
}

/* comment FOR DEVELOPERS ONLY (do not include in docu)
 *
 * This is the so to speak parser for cdo console inputs.
 *
 *  This parser runs over every argv that comes after the cdo options.  The
 *  fist thing that is done is processing the first operator, since this
 *  operator is the only one that has output files we can remove the output
 *  file from out argv by limiting the argc. Obase operators are treated as if
 *  they have a single output since the operator itself takes care of writing
 *  and creating the files it needs for its output. We also create the first
 *  process for the operator and push it on out stack.  Our stack will contain
 *  all operators that do not have all in- or output they need.  After the
 *  first operator is handled the parser will go over each other element in
 *  argv.  Here 4 cases can happen. Only one of the 4 will happen in one
 *  iteration.
 *
 *  If an operator is found we create a new process and add this process as
 *  child to the process on top of the stack. Likewise we add the top process
 *  as parent to the new process. Then the new Process is added to the stack.
 *
 *  Does the argv element represent a file (indicated by the missing '-' in
 *  argv[i][0]) we create a file stream and add it to the process at the top of
 *  the stack.
 *
 *  In case of a '[' or ']' we check if there is only one bracket since we
 *  decided to not allow multiple brackets in the same argv entry.Then we add
 *  ('[') or remove (']') the top of the process stack to a set (named
 *  bracketOperators) which will keep track of which operators used a bracket.
 *  This stack allows to 'mark' an operator so that it is only removed in case
 *  of a ']'.  The ']' indicates that the top process should be removed.and
 *  that it SHOULD have the correct number of inputs.
 *
 *  At the end of each iteration we remove all operators that have all their
 *  inputs AND are not contained in out bracket operators set. So a not closed
 *  bracket will cause a wanted miss function of the parser as the operator
 *  will not be removed and more inputs will be added. This will be found later
 *  by our routine (Process::validate) that checks if every process has the
 *  correct number of inputs and outputs and will throw an error.
 */

ParseStatus
ProcessManager::createProcessesFromInput(int argc, const char **argv)
{
  Cdo_Debug(CdoDebug::PROCESS, "== Process Creation Start ==");
  Cdo_Debug(CdoDebug::PROCESS, "operators:  ", CdoDebug::argvToString(argc, argv));

  ProcessType *root_process = createProcess(argv[0]);
  int cntOutFiles = (int) root_process->m_module.streamOutCnt;

  unsigned int maxIdx = argc - cntOutFiles;
  if (cntOutFiles == -1)
    {
      root_process->m_obase = argv[argc - 1];
      cntOutFiles = 1;
      maxIdx = argc - 1;
    }
  else
    {
      if (maxIdx <= 0)
        {
          return ParseStatus::MissingOutFile;
        }
      handleFirstOperator(maxIdx, argc, argv, root_process);
    }

  ProcessType *currentProcess;
  std::stack<ProcessType *> processStack;
  std::set<ProcessType *> bracketOperators;
  const char *argvEntry;
  int unclosedBrackets = 0;
  unsigned int idx = 1;

  processStack.push(root_process);
  while (!processStack.empty() && idx < maxIdx)
    {
      currentProcess = processStack.top();
      Cdo_Debug(CdoDebug::PROCESS, "iteration ", idx, ", current argv: ", argv[idx],
                ",  currentProcess: ", currentProcess->m_operatorCommand);

      argvEntry = argv[idx];
      //------------------------------------------------------
      // case: operator
      if (argvEntry[0] == '-')
        {
          Cdo_Debug(CdoDebug::PROCESS, "Found new Operator: ", argvEntry);
          currentProcess = addProcess(currentProcess, argvEntry);
          processStack.push(currentProcess);
        }
      // - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // case: bracket start
      else if (argvEntry[0] == '[')
        {
          checkSingleBracketOnly(argvEntry, '[');
          bracketOperators.insert(currentProcess);
          unclosedBrackets++;
        }
      // - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // case: bracket end
      else if (argvEntry[0] == ']')
        {
          checkSingleBracketOnly(argvEntry, ']');
          unclosedBrackets--;
          bracketOperators.erase(processStack.top());
          processStack.pop();
        }
      // - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // case: file
      else if (currentProcess->m_module.streamInCnt != 0)
        {
          Cdo_Debug(CdoDebug::PROCESS, "adding in file to ", currentProcess->operatorName);
          currentProcess->addFileInStream(argvEntry);
        }
      // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      // remove finished
      while (!processStack.empty() && processStack.top()->hasAllInputs()
             && bracketOperators.find(processStack.top()) == bracketOperators.end())
        {
          Cdo_Debug(CdoDebug::PROCESS, "Removing ", processStack.top()->operatorName, " from stack");
          processStack.pop();
        }
      //------------------------------------------------------
      idx++;
    }
  //---------------------------------------------------------------
  if (unclosedBrackets > 0)
    {
      return ParseStatus::ClosingBracketMissing;
    }
  if (unclosedBrackets < 0)
    {
      return ParseStatus::OpenBracketMissing;
    }
  if (idx < maxIdx)
    {
      if (argv[idx][0] == ']')
        {
          return ParseStatus::OpenBracketMissing;
        }
      return ParseStatus::UnprocessedInput;
    }

  Cdo_Debug(CdoDebug::PROCESS, "== Process Creation End ==");

  setProcessNum(m_processes.size());
  return ParseStatus::Ok;
}

void
ProcessManager::handleParseError(ParseStatus p_errCode)
{
  switch (p_errCode)
    {
    case ParseStatus::ClosingBracketMissing:
      {
        CdoError::Abort(Cdo::progname, "Missing ']'.");
        break;
      }
    case ParseStatus::OpenBracketMissing:
      {
        CdoError::Abort(Cdo::progname, "Missing '['.");
        break;
      }
    case ParseStatus::UnprocessedInput:
      {
        CdoError::Abort(Cdo::progname, "Unprocessed Input, could not process all Operators/Files");
        break;
      }
    case ParseStatus::MissingOutFile:
      {
        CdoError::Abort(Cdo::progname, "Missing out file for first operator");
        break;
      }
    case ParseStatus::Ok: { return;
      }
    }
}

ProcessType &
ProcessManager::getProcess(int p_processID)
{
  pthread_mutex_lock(&processMutex);
  auto process = m_processes.find(p_processID);
  if (process == m_processes.end())
    {
      ERROR("Process with ID: ", p_processID, " not found");
    }
  pthread_mutex_unlock(&processMutex);
  return process->second;
}

ProcessType &
ProcessManager::getProcessByThreadID(pthread_t p_threadID)
{
  for (auto &id_process_pair : m_processes)
    {
      if (id_process_pair.second.m_isActive)
        {
          if (pthread_equal(id_process_pair.second.threadID, p_threadID))
            {
              return id_process_pair.second;
            }
        }
    }
  CdoError::Abort(Cdo::progname, "could not find process for thread id: ", p_threadID);
}

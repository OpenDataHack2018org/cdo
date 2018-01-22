#include "process.h"


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

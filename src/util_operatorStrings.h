#ifndef UTIL_OPERATORSTRINGS_H
#define UTIL_OPERATORSTRINGS_H

#include <array>
#include <vector>

typedef  std::vector<std::array<std::string,2>> CdoArgv;
typedef  std::array<std::string,2> cdoCommand;

cdoCommand split(std::string p_command);

CdoArgv preProcessArgv(std::vector<std::string> argv);

const char *getOperatorName(const char *operatorCommand);

char *getOperatorArg(const char *p_operatorCommand);

#endif;

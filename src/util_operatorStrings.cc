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
#include <array>
#include <vector>
#include "dmemory.h"
#include "cdoDebugOutput.h"

typedef  std::array<std::vector<std::string>,2> CdoArgv;
typedef  std::array<std::string,2> cdoCommand;
enum class Command{Name = 0, Arg = 1};

cdoCommand split(std::string p_command)
{
    cdoCommand splitted{""};
    int pos =  p_command.find(',');

    splitted[Command::name] += p_command.substr(0, pos);
    if(pos != std::string::npos)
        splitted[CommandArg] += p_command.substr(pos + 1);

    return splitted;
}

CdoArgv preProcessArgv(std::vector<std::string> argv)
{
    CdoArgv preProcessedArgv;
    for(std::string arg : argv)
    {
        std::array<std::string,2> refinedCommand;
        refinedCommand = split(arg);
        preProcessedArgv[Command::Name].push_back(refinedCommand[Command::Name]);
        preProcessedArgv[Command::Arg].push_back(refinedCommand[Command::Arg]);
        MESSAGE("processed: ", refinedCommand[Command::Name], " ", refinedCommand[Command::Arg]);
    }
    return preProcessedArgv;
}

const char *getOperatorName(const char *operatorCommand)
{
  char *operatorName = NULL;

  if ( operatorCommand )
    {
      if ( operatorCommand[0] == '-' ) operatorCommand++;
      char *commapos = (char *)strchr(operatorCommand, ',');
      size_t len = (commapos != NULL) ? (size_t)(commapos - operatorCommand) : strlen(operatorCommand);

      operatorName = (char*) Malloc(len+1);

      memcpy(operatorName, operatorCommand, len);
      operatorName[len] = '\0';
    }

  /*  return operatorName; */
  if(is_alias(operatorName))
    {
      operatorName = get_original(operatorName);
    }

  return operatorName;
}

char *getOperatorArg(const char *p_operatorCommand)
{
  char *operatorCommand = NULL;

  if ( p_operatorCommand )
    {
      char *commapos = (char *)strchr(p_operatorCommand, ',');
      if ( commapos )
        {
          size_t len = strlen(commapos+1);
          if ( len )
            {
              operatorCommand = (char*) Malloc(len+1);
              strcpy(operatorCommand, commapos+1);
            }
        }
    }

  return operatorCommand;
}


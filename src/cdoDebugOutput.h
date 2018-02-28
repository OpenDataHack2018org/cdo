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
#ifndef DEBUG_SWITCHES_H
#define DEBUG_SWITCHES_H
#include <iostream>
#include <sstream>

#include <fstream>
#include <string.h>

namespace CdoLog
{
void StdOut(std::stringstream &message);

template <typename... T>
static void
expand(std::stringstream &p_message, T &&... args)
{
  // for showing that the dummy array is never used
  using expander = int[];
  // creating dummy array with expanding the parameter pack in its
  // initializer list while writing each element of the pack into message
  expander{ 0, (void(p_message << std::forward<T>(args)), 0)... };
  p_message << std::endl;
}

template <typename... T>
void
StdOut(T &&... args)
{
  std::stringstream message;
  expand(message, args...);
  std::cout << message.str();
}

template <typename... T>
void
StdErr(T &&... args)
{
  std::stringstream message;
  expand(message, args...);
  std::cout << message.str();
}
}  // namespace CdoLog

namespace CdoDebug
{
// Debug Switches
extern int cdoDebug;
extern int cdoDebugExt;  //  Debug level for the KNMI extensions
// Subsystem Debug Switches
extern int PSTREAM;
extern int PROCESS;
extern int PIPE;
extern int ARGUMENT;
extern int PTHREAD;

// File switches and streams
extern std::string outfile;
extern bool print_to_seperate_file;
extern std::fstream outfile_stream;

std::string get_padding(const char *p_func);

void CdoStartMessage();
void CdoEndMessage();
void SetDebug(int p_debug_level);
std::string argvToString(int argc, const char **argv);

void printMessage(std::stringstream &p_message, bool both = false);
void printError(std::stringstream &p_message, bool both = false);

template <typename... T>
void
Message_(const char *p_func, T &&... args)
{
  std::stringstream message;
  message << p_func << ": " << get_padding(p_func);
  CdoLog::expand(message, args...);
  printMessage(message);
}

template <typename... T>
void
PrintDebug(const char *p_func, int p_debugScope, T &&... args)
{
  if (p_debugScope > 0)
    {
      std::stringstream message;
      message << p_func << ": " << get_padding(p_func);
      CdoLog::expand(message, args...);
      printMessage(message);
    }
}

template <typename... T>
void
Warning_(T &&... args)
{
  std::stringstream message;
  message << "Warning: ";
  CdoLog::expand(message, args...);
  printMessage(message, true);
}
}  // namespace CdoDebug

namespace CdoError
{
static int _ExitOnError = 1;

template <typename... T>
void
Abort(const char *progname, T &&... args)
{
  std::stringstream message;
  message << "\n" << progname << " (Abort): ";
  CdoLog::expand(message, args...);
  CdoDebug::printError(message, true);
  if (CdoError::_ExitOnError)
    {
      if (CdoDebug::print_to_seperate_file) CdoDebug::outfile_stream.close();
      exit(EXIT_FAILURE);
    }
}

template <typename... T>
void
Error_(const char *p_file, const int p_line, const char *caller, T &&... args)
{
  std::stringstream message;
  message << "Error in: " << p_file << ":" << p_line << " " << caller << " ";
  CdoLog::expand(message, args...);
  CdoDebug::printError(message, true);
  if (CdoError::_ExitOnError)
    {
      if (CdoDebug::print_to_seperate_file) CdoDebug::outfile_stream.close();
      exit(EXIT_FAILURE);
    }
}

template <typename... T>
void
SysError_(const char *p_file, const int p_line, const char *p_func, T &&... args)
{
  int saved_errno = errno;
  std::stringstream message;
  message << "SysError in:" << p_file << std::endl;
  message << "    "
          << "in function: p_func ,line: " << p_line << std::endl;
  CdoLog::StdOut(message, args...);
  if (saved_errno)
    {
      errno = saved_errno;
      perror("Sytem error message");
    }
  exit(EXIT_FAILURE);
}
}  // namespace CdoError
#if defined WITH_CALLER_NAME
#define SYS_ERROR(...) CdoError::SysError_(__FILE__, __LINE__, __func__, __VA_ARGS__)
#define ERROR_C(...) CdoError::Error_(__FILE__, __LINE__, caller, __VA_ARGS__)
#define ERROR(...) CdoError::Error_(__FILE__, __LINE__, __func__, __VA_ARGS__)
#define WARNING(...) CdoError::Warning_(__func__, __VA_ARGS__)
#define MESSAGE_C(...) CdoDebug::Message_(caller, __VA_ARGS__)
#define MESSAGE(...) CdoDebug::Message_(__func__, __VA_ARGS__)
#else
#define SYS_ERROR(...) CdoError::SysError_(__FILE__, __LINE__, "", __VA_ARGS__)
#define ERROR_C(...) CdoError::Error_(__FILE__, __LINE__, "", __VA_ARGS__)
#define ERROR(...) CdoError::Error_(__FILE__, __LINE__, "", __VA_ARGS__)
#define WARNING(...) CdoError::Warning_(__VA_ARGS__)
#define MESSAGE_C(...) CdoDebug::Message_(__VA_ARGS__)
#define MESSAGE(...) CdoDebug::Message_(__func__, __VA_ARGS__)
#endif
#define Cdo_Debug(...) CdoDebug::PrintDebug(__func__, __VA_ARGS__)
#endif

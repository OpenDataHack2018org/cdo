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

#ifndef PROCESS_INT_H
#define PROCESS_INT_H

#include "process.h"

/**
 * Sets the underlying Process threadID and calls omp_set_num_threads
 */
void cdoInitialize(void *process);
void cdoFinish();
void printEndTimes(char *p_memstring);
/**
 * Add operator to process.
 *
 * @param name name of the operator
 * @param f1 value for use in operator(e.g function ID from field.cc)
 * @param f2 value for use in operator
 * @param enter string for inquiring user input
 * @return operatorID
 */
int cdoOperatorAdd(const char *name, int f1, int f2, const char *enter);
int cdoOperatorID(void);
int cdoOperatorF1(int operID);
int cdoOperatorF2(int operID);

/**
 * Returns registered operator name
 */
const char *cdoOperatorName(int operID);
/**
 * Returns operator name
 */
const char *cdoOperatorEnter(int operID);

/**
 * Returns the input stream name as std::string for inStreamID.
 * Aborts Cdo when the ID is larger than the current in stream count
 */
std::string cdoGetInStreamName(int inStreamID);

/**
 * Returns the output stream name as std::string for \p outStreamID.
 * Aborts Cdo when the ID is larger than the current out stream count
 */
std::string cdoGetOutStreamName(int p_outStream);
/**
 * Returns a stream name for given index.
 * if the given index is smaller than the current count of input files
 * the returned name will be a input stream.
 * Otherwise a name for an output stream will be returned.
 * Aborts if the streamID is larger than inStreamCnt + outStreamCnt
 */
std::string cdoGetStreamName(int p_streamIndex);

/**
 *  Returns the basis name for output files.
 *  If no obase was found Cdo will exit
 */
char *cdoGetObase();

/**
 * Returns the number type for the operator used in the current process
 */
int cdoStreamNumber();
/**
 * Returns the current number of in and out streams that were added to this
 * cdoProcess
 */

int cdoStreamCnt(void);
/**
 * Deprecated. Returns parameter cnt
 */
int cdoStreamName(int cnt);

/**
 * Returns the stream ID for \p inStreamIDX. If the input is from another
 * process, a new thread will be started and the found process will run in that
 * thread. Is the ID bigger thatn the count of input streams cdo will exit.
 */
int cdoStreamOpenRead(int inStreamIDX);
int cdoStreamOpenWrite(int outStreamIDX, int filetype);
int cdoStreamOpenWrite(std::string p_filename, int filetype);

/**
 * returns pstreamID for ID \p outStreamIDX.
 *
 */
int cdoStreamOpenAppend(int outStreamIDX);

/**
 * Checks if the output file name behind outStreamIDX exists.
 *  @param[in] outStreamIDX streamID that represents the file
 *  @param[out] return true if file exists else false
 *  true if file exists else false
 */
bool cdoOutFileExists(int outStreamIDX);
/**
 * Checks if the output file behind outStreamIDX exists.
 * [out] true if file exists else false
 */
bool cdoInFileExists(int inStreamIDX);

/**
 * Checks if this process has \p numargs arguments
 */
void operatorCheckArgc(int numargs);

/**
 * Asks user for input.
 * @param enter string that will be displayed when asking for input from the
 * user
 */
void operatorInputArg(const char *enter);
/**
 * Returns pointer to first element of operator arguments
 */
char **operatorArgv(void);

/**
 * Returns the number of operator arguments
 */
int operatorArgc(void);

void processDefVarNum(int nvars);

/**
 * Returns number of variables
 */
int processInqVarNum();
/**
 * Returns c string representation of this process
 */
const char *processInqPrompt(void);

/**
 * Returns number of created processes
 */
int processNums(void);
/**
 * Returns number of running processes.
 */
int processNumsActive();

/**
 * Returns the surrounding process
 */
ProcessType &processSelf(void);
/**
 * Returns process with id \p p_processID
 */
ProcessType *getProcess(int p_processID);

/**
 * Deletes all processes
 */
void clearProcesses();

void processAccuTime(double utime, double stime);
void processStartTime(double *utime, double *stime);

/**
 * Creates processes from argv.
 * @param argc standard argc  minus the options
 * @param argv standard argv after option processing starting with the first
 * operator
 */
void createProcesses(int argc, const char **argv);

/**
 * Creates single process from cdo operator command
 */
ProcessType *processCreate(const char *command);

/**
 * Returns vlistID for pstream with ID = \p pstreamID
 * cdo will exit if no stream for \p pstreamID is found
 */
int cdoStreamInqVlist(int pstreamID);

/**
 * Checks whether pstream with \p pstreamID is a file or a pipe.
 * @return true if pipe otherwise false
 */
bool cdoStreamIsPipe(int pstreamID);

void runProcesses();

void killProcesses();

int cdoStreamInqTimestep(int pstreamID, int tsID);

void cdoValidateProcesses();

#endif

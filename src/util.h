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

#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdbool.h>
#include <string>

/* dummy use of unused parameters to silence compiler warnings */
#define UNUSED(x) (void) x

#undef TRUE
#define TRUE 1
#undef FALSE
#define FALSE 0

#undef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#undef MAX
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#undef SQR
#define SQR(a) ((a) * (a))

#define UNCHANGED_RECORD                                                                                \
  (processSelf().m_ID == 0 && processSelf().inputStreams[0]->ispipe == false && cdoRegulargrid == FALSE \
   && cdoDefaultFileType == -1 && cdoDefaultDataType == -1 && cdoDefaultByteorder == -1)

extern const char *CDO_version;
extern const char *CDO_username;
extern char *cdoGridSearchDir;
extern int CDO_Reduce_Dim;
extern int CDO_Memtype;
extern int CDO_Parallel_Read;
extern int CDO_Append_History;
extern int CDO_Reset_History;

extern int CDO_optind;
extern const char *CDO_optarg;
extern int CDO_opterr;

extern int CDO_flt_digits;
extern int CDO_dbl_digits;

extern bool REMAP_genweights;

extern const char *cdoExpName;

extern int stdin_is_tty;
extern int stdout_is_tty;
extern int stderr_is_tty;

extern int cdoDefaultFileType;
extern int cdoDefaultDataType;
extern int cdoDefaultByteorder;
extern int cdoDefaultTableID;
extern int cdoDefaultInstID;
extern int cdoDefaultTimeType;

extern int cdoCheckDatarange;

extern int cdoOverwriteMode;
extern int cdoRegulargrid;
extern int cdoTimer;
extern int cdoVerbose;
extern int cdoParIO;

extern int cdoChunkType;

extern int CDO_Color;
extern int CDO_Use_FFTW;
extern int CDO_Version_Info;
extern int CDO_CMOR_Mode;
extern int cdoDiag;

extern int cdoNumVarnames;
extern char **cdoVarnames;
extern char CDO_File_Suffix[32];  // refactor: added keyword extern

const char *getProgname(char *string);
const char *getOperatorName(const char *operatorCommand);
char *getOperatorArg(const char *operatorCommand);
const char *cdoComment(void);

char *getFileArg(char *argument);

enum
{
  START_DEC,
  START_JAN
};
int get_season_start(void);
void get_season_name(const char *seas_name[]);
int month_to_season(int month);

void init_is_tty(void);

char *double_to_attstr(int digits, char *str, size_t len, double value);

void progressInit(void);
void progressStatus(double offset, double refval, double curval);

/* convert a CDI datatype to string */
int datatype2str(int datatype, char *datatypestr);
int str2datatype(const char *datatypestr);

int cdoFiletype(void);

void cdoSetNAN(double missval, size_t gridsize, double *array);

int cdoDefineGrid(const char *gridfile);
int cdoDefineZaxis(const char *zaxisfile);

int vlistInqNWPV(int vlistID, int varID);
int vlistIsSzipped(int vlistID);
size_t vlist_check_gridsize(int vlistID);
int vlist_get_psvarid(int vlistID, int zaxisID);
double *vlist_read_vct(int vlistID, int *rzaxisIDh, int *rnvct, int *rnhlev, int *rnhlevf, int *rnhlevh);
void vlist_change_hybrid_zaxis(int vlistID1, int vlistID2, int zaxisID1, int zaxisID2);

void cdoGenFileSuffix(char *filesuffix, size_t maxlen, int filetype, int vlistID, const char *refname);

void writeNCgrid(const char *gridfile, int gridID, int *imask);
void defineZaxis(const char *zaxisarg);

int grid_from_name(const char *gridname);
int zaxisFromName(const char *zaxisname);

/* refactor: moved here from cdo.h */
int cdo_omp_get_thread_num(void);
void cdo_omp_set_num_threads(int nthreads);

/* refactor: moved here from cdo.cc */
void exp_run(int argc, char *argv[], const char *cdoExpName);  // job.cc
void printFeatures(void);                                      // features.cc
void printLibraries(void);                                     // features.cc

int wildcardmatch(const char *w, const char *s);

void cdo_check_round(void);

#endif /* UTIL_H */

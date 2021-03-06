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

#ifndef CDO_INT_H
#define CDO_INT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <cmath>
#include <float.h>
#include <ctype.h>
#include <stdarg.h>

#include "pmlist.h"
#include "listbuf.h"
#include "compare.h"
#include "timebase.h"
#include "field.h"
#include "functs.h"
#include "dmemory.h"
#include "const.h"
#include "util.h"
#include "cdo_vlist.h"
#include "process_int.h"
#include "exception.h"

#ifdef _OPENMP
#define OPENMP3 200805
#define OPENMP4 201307
#define OPENMP45 201511

#if _OPENMP >= OPENMP3
#define HAVE_OPENMP3 1
#endif

#if _OPENMP >= OPENMP4
#define HAVE_OPENMP4 1
#endif

#if _OPENMP >= OPENMP45
#define HAVE_OPENMP45 1
#endif
#endif


#define cmpstr(s1, s2) (strncmp(s1, s2, strlen(s2)))
#define cmpstrlen(s1, s2, len) (strncmp(s1, s2, len = strlen(s2)))

/* sxxxYYYYMMDDhhmm0 */
#define DATE_LEN 31 /* YYYYMMDDhhmmss allocate DTLEN+1 !!!! */
#define SET_DATE(dtstr, date, time) (sprintf(dtstr, "%*lld%*d", DATE_LEN - 6, date, 6, time))
#define DATE_IS_NEQ(dtstr1, dtstr2, len) (memcmp(dtstr1, dtstr2, len) != 0)

enum T_WEIGHT_MODE
{
  WEIGHT_OFF,
  WEIGHT_ON
};
enum T_EIGEN_MODE
{
  JACOBI,
  DANIELSON_LANCZOS
};

#ifndef M_LN10
#define M_LN10 2.30258509299404568401799145468436421 /* log_e 10 */
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288 /* pi */
#endif

template <typename T>
inline void
vectorFree(std::vector<T> &v)
{
  if (v.size())
    {
      v.clear();
      v.shrink_to_fit();
    }
}

#define VECTOR_2D(T, P2D, N, M)       \
  std::vector<std::vector<T>> P2D(N); \
  if ((N))                            \
    for (size_t i = 0; i < (size_t)(N); ++i) P2D[i].resize(M);

#define NEW_2D(T, P2D, N, M)                                              \
  T **P2D = (N) ? new T *[(N)] : nullptr;                                 \
  if ((N))                                                                \
    {                                                                     \
      P2D[0] = (M) ? new T[(N) * (M)] : nullptr;                          \
      for (size_t i = 1; i < (size_t)(N); ++i) P2D[i] = P2D[0] + i * (M); \
    }
#define DELETE_2D(P2D)             \
  if (P2D)                         \
    {                              \
      if (P2D[0]) delete[] P2D[0]; \
      delete[] P2D;                \
      P2D = nullptr;               \
    }

#define IX2D(y, x, nx) ((y) * (nx) + (x))

#define MEMTYPE_DOUBLE 1
#define MEMTYPE_FLOAT 2

void print_pthread_info(void);

void cdoProcessTime(double *utime, double *stime);

int readline(FILE *fp, char *line, int len);

int zaxis2ltype(int zaxisID);

double radius_str_to_deg(const char *string);

void datetime2str(int64_t date, int time, char *datetimestr, int maxlen);
void date2str(int64_t date, char *datestr, int maxlen);
void time2str(int time, char *timestr, int maxlen);

const char *tunit2str(int tunits);
const char *calendar2str(int calendar);

void cdo_set_grids(const char *gridarg);
void defineInstitution(const char *instarg);
int defineTable(const char *tablearg);

void nospec(int vlistID);
void gridWrite(FILE *fp, int gridID);

void openLock(void);
void openUnlock(void);

int cdf_openread(const char *filename);

void printFiletype(int streamID, int vlistID);

off_t fileSize(const char *restrict filename);

const char *parameter2word(const char *string);
double parameter2double(const char *string);
bool parameter2bool(const char *string);
int parameter2int(const char *string);
long parameter2long(const char *string);
size_t parameter2sizet(const char *string);
int parameter2intlist(const char *string);

int referenceToGrid(int gridID1);

void cdo_read_field(const char *name, char *pline, int size, double *field, int *lineno, FILE *fp, const char *dname);

double cdoZaxisInqLevel(int zaxisID, int levelID);
int cdoZaxisInqLevels(int zaxisID, double *levels);

list_t *namelistbuf_to_pmlist(listbuf_t *listbuf);
list_t *namelist_to_pmlist(FILE *fp, const char *name);
list_t *cmortable_to_pmlist(FILE *fp, const char *name);

int literals_find_datatype(int n, char **literals);
int literal_get_datatype(const char *literal);
int literal_to_int(const char *literal);
double literal_to_double(const char *literal);

char *cdoVlistInqVarName(int vlistID, int varID, char *name);

#ifdef __cplusplus
extern "C" {
#endif

void cdiDefTableID(int tableID);

void gridGenXvals(int xsize, double xfirst, double xlast, double xinc, double *xvals);
void gridGenYvals(int gridtype, int ysize, double yfirst, double ylast, double yinc, double *yvals);

void gaussaw(double *restrict pa, double *restrict pw, size_t nlat);

int qu2reg3_double(double *pfield, int *kpoint, int klat, int klon, double msval, int *kret, int omisng, int operio,
                   int oveggy);

void cdoCompareGrids(int gridID1, int gridID2);

#if defined(__cplusplus)
}
#endif

#endif /* CDO_INT_H */

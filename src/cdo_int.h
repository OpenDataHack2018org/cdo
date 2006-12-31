/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2006 Uwe Schulzweida, schulzweida@dkrz.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#ifndef _CDO_INT_H
#define _CDO_INT_H

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "field.h"
#include "functs.h"
#include "dmemory.h"
#include "process.h"

#ifndef strdupx
#ifndef strdup
char *strdup(const char *s);
#endif
#define strdupx  strdup
/*
#define strdupx(s)			          \
({					      	  \
   const char *__old = (s);			  \
   size_t __len = strlen(__old) + 1;		  \
   char *__new = (char *) malloc(__len);	  \
   (char *) memcpy(__new, __old, __len);	  \
})
*/
#endif

#ifndef  DBL_IS_EQUAL
/*
#define  DBL_IS_EQUAL(x,y) (fabs(x - y) <= 2.0*(y*DBL_EPSILON + DBL_MIN))
*/
#define  DBL_IS_EQUAL(x,y) (!(fabs(x - y) > 0))
#endif

#ifndef  M_PI
#define  M_PI		3.14159265358979323846	/* pi */
#endif

#ifndef  M_LN10
#define  M_LN10		2.30258509299404568402	/* log_e 10 */
#endif


#ifndef  RAD2DEG
#define  RAD2DEG  (180./M_PI)   /* conversion for rad to deg */
#endif

#ifndef  DEG2RAD
#define  DEG2RAD  (M_PI/180.)   /* conversion for deg to rad */
#endif

#define  IX2D(y,x,nx)  ((y)*(nx)+(x))


#define CDO_EXP_LOCAL   1
#define CDO_EXP_REMOTE  2

void strtolower(char *str);

void print_pthread_info(void);

void cdoProcessTime(double *utime, double *stime);

void    setCommandLine(int argc, char **argv);
char   *commandLine(void);
int     readline(FILE *fp, char *line, int len);

int ztype2ltype(int zaxistype);
int ltype2ztype(int ltype);

int nlat2ntr(int nlat);
int nlat2ntr_linear(int nlat);
int ntr2nlat(int ntr);
int ntr2nlat_linear(int ntr);
int compNlon(int nlat);

typedef struct {
  int      date;
  int      time;
  int      julval;
}
DATETIME;

void    datetime_avg(int dpy, int ndates, DATETIME *datetime);

void    decode_date(int date, int *year, int *month, int *day);
void    decode_time(int time, int *hour, int *minute);
double  encode_julval(int dpy, int date, int time);
void    decode_julval(int dpy, double value, int *date, int *time);
int     days_per_month(int calendar, int year, int month);
int     days_per_year(int calendar, int year);
int     calendar_dpy(int calendar);

void    defineGrid(const char *gridarg);
void    defineInstitution(char *instarg);
int     defineTable(char *tablearg);

void    cdolog(const char *prompt, double cputime);
void    cdologs(int noper);
void    cdologo(int noper);
void    nospec(int vlistID);
void    gridWrite(FILE *fp, int gridID);

void printFiletype(int streamID, int vlistID);

void job_submit(const char *expname, const char *jobfilename, const char *jobname);

#endif  /* _CDO_INT_H */

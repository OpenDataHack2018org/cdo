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
#ifndef DATETIME_H
#define DATETIME_H

#include <stdio.h>

enum struct TimeStat
{
  UNDEF,
  FIRST,
  LAST,
  MEAN,
  MIDHIGH,
};

struct juldate_t
{
  int64_t julday;
  int secofday;
};

struct CdoDateTime
{
  int64_t date;
  int time;
};

struct DateTimeInfo
{
  CdoDateTime c;     // corrected verification time
  CdoDateTime v;     // verification time
  CdoDateTime b[2];  // time bounds
};

struct DateTimeList
{
  size_t nalloc;
  size_t size;
  int has_bounds;
  int calendar;
  TimeStat stat;
  int timestat_date;
  DateTimeInfo timestat;
  DateTimeInfo *dtinfo;
};

juldate_t juldate_encode(int calendar, int64_t date, int time);
void juldate_decode(int calendar, juldate_t juldate, int64_t *date, int *time);
juldate_t juldate_sub(juldate_t juldate2, juldate_t juldate1);
juldate_t juldate_add_seconds(int64_t seconds, juldate_t juldate);
double juldate_to_seconds(juldate_t juldate);

void datetime_avg(int dpy, int ndates, CdoDateTime *datetime);

DateTimeList *dtlist_new(void);
void dtlist_delete(DateTimeList *dtlist);
void dtlist_shift(DateTimeList *dtlist);
void dtlist_set_stat(DateTimeList *dtlist, TimeStat stat);
void dtlist_set_calendar(DateTimeList *dtlist, int calendar);
int64_t dtlist_get_vdate(DateTimeList *dtlist, int tsID);
int dtlist_get_vtime(DateTimeList *dtlist, int tsID);
void dtlist_taxisInqTimestep(DateTimeList *dtlist, int taxisID, int tsID);
void dtlist_taxisDefTimestep(DateTimeList *dtlist, int taxisID, int tsID);
void dtlist_stat_taxisDefTimestep(DateTimeList *dtlist, int taxisID, int nsteps);

#endif /* DATETIME_H */

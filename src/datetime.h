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

typedef struct
{
  int julday;
  int secofday;
} juldate_t;

typedef struct
{
  int date;
  int time;
} cdo_datetime_t;

typedef struct
{
  cdo_datetime_t c;     // corrected verification time
  cdo_datetime_t v;     // verification time
  cdo_datetime_t b[2];  // time bounds
} dtinfo_type;

typedef struct
{
  size_t nalloc;
  size_t size;
  int has_bounds;
  int calendar;
  TimeStat stat;
  int timestat_date;
  dtinfo_type timestat;
  dtinfo_type *dtinfo;
} dtlist_type;

juldate_t juldate_encode(int calendar, int date, int time);
void juldate_decode(int calendar, juldate_t juldate, int *date, int *time);
juldate_t juldate_sub(juldate_t juldate2, juldate_t juldate1);
juldate_t juldate_add_seconds(int seconds, juldate_t juldate);
double juldate_to_seconds(juldate_t juldate);

void datetime_avg(int dpy, int ndates, cdo_datetime_t *datetime);

dtlist_type *dtlist_new(void);
void dtlist_delete(dtlist_type *dtlist);
void dtlist_shift(dtlist_type *dtlist);
void dtlist_set_stat(dtlist_type *dtlist, TimeStat stat);
void dtlist_set_calendar(dtlist_type *dtlist, int calendar);
int dtlist_get_vdate(dtlist_type *dtlist, int tsID);
int dtlist_get_vtime(dtlist_type *dtlist, int tsID);
void dtlist_taxisInqTimestep(dtlist_type *dtlist, int taxisID, int tsID);
void dtlist_taxisDefTimestep(dtlist_type *dtlist, int taxisID, int tsID);
void dtlist_stat_taxisDefTimestep(dtlist_type *dtlist, int taxisID, int nsteps);

#endif /* DATETIME_H */

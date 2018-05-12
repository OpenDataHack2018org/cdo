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

class DateTimeList
{
 public:
  DateTimeList();
  ~DateTimeList();
  void setStat(TimeStat stat);
  void setCalendar(int calendar);
  int64_t getVdate(int tsID);
  int getVtime(int tsID);
  void taxisInqTimestep(int taxisID, int tsID);
  void taxisDefTimestep(int taxisID, int tsID);
  void statTaxisDefTimestep(int taxisID, int nsteps);
  void shift();

 private:
  size_t nalloc;
  size_t size;
  int has_bounds;
  int calendar;
  TimeStat stat;
  int timestat_date;
  DateTimeInfo timestat;
  std::vector<DateTimeInfo> dtinfo;

  void mean(int nsteps);
  void midhigh(int nsteps);
};

juldate_t juldate_encode(int calendar, int64_t date, int time);
void juldate_decode(int calendar, juldate_t juldate, int64_t *date, int *time);
juldate_t juldate_sub(juldate_t juldate2, juldate_t juldate1);
juldate_t juldate_add_seconds(int64_t seconds, juldate_t juldate);
double juldate_to_seconds(juldate_t juldate);

void datetime_avg(int dpy, int ndates, CdoDateTime *datetime);

#endif /* DATETIME_H */

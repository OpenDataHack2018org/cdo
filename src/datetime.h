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
  void init();

  DateTimeList()
    {
      nalloc = 0;
      size = 0;
      calendar = -1;
      has_bounds = -1;
      stat = TimeStat::LAST;

      init();
    }

  ~DateTimeList()
    {
    }

  void setStat(TimeStat stat)
  {
    this->stat = stat;
  }

  void setCalendar(int calendar)
  {
    this->calendar = calendar;
  }

  int64_t getVdate(int tsID)
  {
    if (tsID < 0 || (size_t) tsID >= this->size) cdoAbort("Internal error; tsID out of bounds!");
  
    return this->dtinfo[tsID].c.date;
  }

  int getVtime(int tsID)
  {
    if (tsID < 0 || (size_t) tsID >= this->size) cdoAbort("Internal error; tsID out of bounds!");
    
    return this->dtinfo[tsID].c.time;
  }

  void shift()
  {
    for (size_t inp = 0; inp < this->size - 1; inp++) this->dtinfo[inp] = this->dtinfo[inp + 1];
  }

  void taxisInqTimestep(int taxisID, int tsID);
  void taxisDefTimestep(int taxisID, int tsID);
  void statTaxisDefTimestep(int taxisID, int nsteps);

 private:
  size_t nalloc;
  size_t size;
  int has_bounds;
  int calendar;
  TimeStat stat;
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

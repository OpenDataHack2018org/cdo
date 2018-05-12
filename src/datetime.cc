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
#include <cdi.h>
#include "cdo_int.h"
#include "datetime.h"
#include "util_string.h"

TimeStat CDO_Timestat_Date = TimeStat::UNDEF;
bool CDO_Timestat_Bounds = false;
static bool dateTimeInit = false;

void
setTimestatDate(const char *optarg)
{
  TimeStat timestatdate = TimeStat::UNDEF;

  // clang-format off
  if      (STR_IS_EQ(optarg, "first"))   timestatdate = TimeStat::FIRST;
  else if (STR_IS_EQ(optarg, "last"))    timestatdate = TimeStat::LAST;
  else if (STR_IS_EQ(optarg, "middle"))  timestatdate = TimeStat::MEAN;
  else if (STR_IS_EQ(optarg, "midhigh")) timestatdate = TimeStat::MIDHIGH;
  // clang-format on

  if (timestatdate == TimeStat::UNDEF) cdoAbort("option --%s: unsupported argument: %s", "timestat_date", optarg);

  CDO_Timestat_Date = timestatdate;
}

static void
getTimestatDate(TimeStat *tstat_date)
{
  char *envstr = getenv("CDO_TIMESTAT_DATE");
  if (envstr == NULL) envstr = getenv("RUNSTAT_DATE");
  if (envstr)
    {
      TimeStat env_date = TimeStat::UNDEF;
      char envstrl[8];

      memcpy(envstrl, envstr, 8);
      envstrl[7] = 0;
      strtolower(envstrl);

      // clang-format off
      if      (memcmp(envstrl, "first", 5) == 0)   env_date = TimeStat::FIRST;
      else if (memcmp(envstrl, "last", 4) == 0)    env_date = TimeStat::LAST;
      else if (memcmp(envstrl, "middle", 6) == 0)  env_date = TimeStat::MEAN;
      else if (memcmp(envstrl, "midhigh", 7) == 0) env_date = TimeStat::MIDHIGH;
      // clang-format on

      if (env_date != TimeStat::UNDEF)
        {
          *tstat_date = env_date;
          if (cdoVerbose) cdoPrint("Set CDO_TIMESTAT_DATE to %s", envstr);
        }
    }
}

void
DateTimeList::init()
{
  if (!dateTimeInit) getTimestatDate(&CDO_Timestat_Date);
  dateTimeInit = true;
}

void
DateTimeList::taxisInqTimestep(int taxisID, int tsID)
{
  size_t NALLOC = 128;

  if ((size_t) tsID >= this->nalloc)
    {
      this->nalloc += NALLOC;
      this->dtinfo.resize(this->nalloc);
    }

  if ((size_t) tsID >= this->size) this->size = (size_t) tsID + 1;

  this->dtinfo[tsID].v.date = taxisInqVdate(taxisID);
  this->dtinfo[tsID].v.time = taxisInqVtime(taxisID);

  this->dtinfo[tsID].c.date = this->dtinfo[tsID].v.date;
  this->dtinfo[tsID].c.time = this->dtinfo[tsID].v.time;

  if (tsID == 0)
    {
      if (this->has_bounds == -1)
        {
          this->has_bounds = 0;
          if (taxisHasBounds(taxisID)) this->has_bounds = 1;
        }

      if (this->calendar == -1)
        {
          this->calendar = taxisInqCalendar(taxisID);
        }
    }

  if (this->has_bounds)
    {
      taxisInqVdateBounds(taxisID, &(this->dtinfo[tsID].b[0].date), &(this->dtinfo[tsID].b[1].date));
      taxisInqVtimeBounds(taxisID, &(this->dtinfo[tsID].b[0].time), &(this->dtinfo[tsID].b[1].time));

      if (CDO_Timestat_Bounds && this->dtinfo[tsID].v.date == this->dtinfo[tsID].b[1].date
          && this->dtinfo[tsID].v.time == this->dtinfo[tsID].b[1].time)
        {
          int calendar = this->calendar;

          int64_t vdate = this->dtinfo[tsID].b[0].date;
          int vtime = this->dtinfo[tsID].b[0].time;
          juldate_t juldate1 = juldate_encode(calendar, vdate, vtime);

          vdate = this->dtinfo[tsID].b[1].date;
          vtime = this->dtinfo[tsID].b[1].time;
          juldate_t juldate2 = juldate_encode(calendar, vdate, vtime);

          // int hour, minute, second;
          // cdiDecodeTime(vtime, &hour, &minute, &second);

          if (vtime == 0 && juldate_to_seconds(juldate1) < juldate_to_seconds(juldate2))
            {
              juldate_t juldate = juldate_add_seconds(-1, juldate2);
              juldate_decode(calendar, juldate, &vdate, &vtime);

              this->dtinfo[tsID].c.date = vdate;
              this->dtinfo[tsID].c.time = vtime;
            }
        }
    }
  else
    {
      this->dtinfo[tsID].b[0].date = 0;
      this->dtinfo[tsID].b[1].date = 0;
      this->dtinfo[tsID].b[0].time = 0;
      this->dtinfo[tsID].b[1].time = 0;
    }
}

void
DateTimeList::taxisDefTimestep(int taxisID, int tsID)
{
  if (tsID < 0 || (size_t) tsID >= this->size) cdoAbort("Internal error; tsID out of bounds!");

  taxisDefVdate(taxisID, this->dtinfo[tsID].v.date);
  taxisDefVtime(taxisID, this->dtinfo[tsID].v.time);
  if (this->has_bounds)
    {
      taxisDefVdateBounds(taxisID, this->dtinfo[tsID].b[0].date, this->dtinfo[tsID].b[1].date);
      taxisDefVtimeBounds(taxisID, this->dtinfo[tsID].b[0].time, this->dtinfo[tsID].b[1].time);
    }
}

void
DateTimeList::mean(int nsteps)
{
  int64_t vdate;
  int vtime;

  if (nsteps % 2 == 0)
    {
      int calendar = this->calendar;

//#define TEST_DTLIST_MEAN 1
#ifdef TEST_DTLIST_MEAN
      vdate = this->dtinfo[0].v.date;
      vtime = this->dtinfo[0].v.time;
      juldate_t juldate0 = juldate_encode(calendar, vdate, vtime);

      juldate_t juldate;
      double seconds = 0;
      for (int i = 1; i < nsteps; ++i)
        {
          vdate = this->dtinfo[i].v.date;
          vtime = this->dtinfo[i].v.time;
          juldate = juldate_encode(calendar, vdate, vtime);

          seconds += juldate_to_seconds(juldate_sub(juldate, juldate0));
        }

      juldate = juldate_add_seconds(lround(seconds / nsteps), juldate0);
      juldate_decode(calendar, juldate, &vdate, &vtime);
#else
      vdate = this->dtinfo[nsteps / 2 - 1].v.date;
      vtime = this->dtinfo[nsteps / 2 - 1].v.time;
      juldate_t juldate1 = juldate_encode(calendar, vdate, vtime);

      vdate = this->dtinfo[nsteps / 2].v.date;
      vtime = this->dtinfo[nsteps / 2].v.time;
      juldate_t juldate2 = juldate_encode(calendar, vdate, vtime);

      double seconds = juldate_to_seconds(juldate_sub(juldate2, juldate1)) / 2;
      juldate_t juldatem = juldate_add_seconds(lround(seconds), juldate1);
      juldate_decode(calendar, juldatem, &vdate, &vtime);
#endif
    }
  else
    {
      vdate = this->dtinfo[nsteps / 2].v.date;
      vtime = this->dtinfo[nsteps / 2].v.time;
    }

  this->timestat.v.date = vdate;
  this->timestat.v.time = vtime;
}

void
DateTimeList::midhigh(int nsteps)
{
  this->timestat.v.date = this->dtinfo[nsteps / 2].v.date;
  this->timestat.v.time = this->dtinfo[nsteps / 2].v.time;
}

void
DateTimeList::statTaxisDefTimestep(int taxisID, int nsteps)
{
  if ((size_t) nsteps > this->size) cdoAbort("Internal error; unexpected nsteps=%d (limit=%ld)!", nsteps, this->size);

  TimeStat stat = this->stat;
  if (CDO_Timestat_Date != TimeStat::UNDEF) stat = CDO_Timestat_Date;

  // clang-format off
  if      (stat == TimeStat::MEAN)    this->mean(nsteps);
  else if (stat == TimeStat::MIDHIGH) this->midhigh(nsteps);
  else if (stat == TimeStat::FIRST)   this->timestat.v = this->dtinfo[0].v;
  else if (stat == TimeStat::LAST)    this->timestat.v = this->dtinfo[nsteps - 1].v;
  else cdoAbort("Internal error; implementation missing for timestat=%d", (int)stat);
  // clang-format on

  if (this->has_bounds)
    {
      this->timestat.b[0] = this->dtinfo[0].b[0];
      this->timestat.b[1] = this->dtinfo[nsteps - 1].b[1];
    }
  else
    {
      this->timestat.b[0] = this->dtinfo[0].v;
      this->timestat.b[1] = this->dtinfo[nsteps - 1].v;
    }

  taxisDefVdate(taxisID, this->timestat.v.date);
  taxisDefVtime(taxisID, this->timestat.v.time);
  // if ( this->has_bounds )
  {
    taxisDefVdateBounds(taxisID, this->timestat.b[0].date, this->timestat.b[1].date);
    taxisDefVtimeBounds(taxisID, this->timestat.b[0].time, this->timestat.b[1].time);
  }
}

void
datetime_avg(int calendar, int ndates, CdoDateTime *datetime)
{
  int64_t vdate;
  int vtime;

  if (ndates % 2 == 0)
    {
      vdate = datetime[ndates / 2 - 1].date;
      vtime = datetime[ndates / 2 - 1].time;
      juldate_t juldate1 = juldate_encode(calendar, vdate, vtime);

      vdate = datetime[ndates / 2].date;
      vtime = datetime[ndates / 2].time;
      juldate_t juldate2 = juldate_encode(calendar, vdate, vtime);

      double seconds = juldate_to_seconds(juldate_sub(juldate2, juldate1)) / 2;
      juldate_t juldatem = juldate_add_seconds(lround(seconds), juldate1);
      juldate_decode(calendar, juldatem, &vdate, &vtime);
    }
  else
    {
      vdate = datetime[ndates / 2].date;
      vtime = datetime[ndates / 2].time;
    }

  datetime[ndates].date = vdate;
  datetime[ndates].time = vtime;
}

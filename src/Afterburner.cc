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
/* ============================================================= */
/*                                                               */
/* postprocessing program for ECHAM data and ECMWF analysis data */
/*                                                               */
/* Luis     Kornblueh   - MPI    Hamburg                         */
/* Uwe      Schulzweida - MPI    Hamburg                         */
/* Arno     Hellbach    - DKRZ   Hamburg                         */
/* Edilbert Kirk        - MI Uni Hamburg                         */
/* Michael  Ponater     - DLR    Oberpfaffenhofen                */
/*                                                               */
/* ============================================================= */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cdi.h>

#ifdef CDO

#include "cdo_int.h"
#include "cdo_task.h"
#include "pstream_int.h"
#define streamOpenWrite cdoStreamOpenWrite
#define streamDefVlist pstreamDefVlist
#define streamDefTimestep pstreamDefTimestep
#endif

#ifdef AFTERBURNER
#include "afterdoc.h"
#endif

#include "afterburner.h"
#include "constants.h"
#include "compare.h"
#include "cdoOptions.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

#if !defined(VERSION)
#define VERSION "0.0.1"
#endif

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC 1000000
#endif

#ifdef AFTERBURNER
static double starttime = 0.0;
#endif

#ifdef AFTERBURNER
void afterInqHistory(int fileID);
void afterDefHistory(int fileID, char *histstring);

long get_nfft(void);
#endif

int scan_par_obsolate(char *namelist, const char *name, int def);
void scan_code(char *namelist, struct Variable *vars, int maxCodes, int *numCodes);
int scan_par(int verbose, char *namelist, const char *name, int def);
int scan_time(int verbose, char *namelist, int *hours, int max_hours);
void scan_darray(char *namelist, const char *name, double *values, int maxValues, int *numValues);

char zaxistypename[CDI_MAX_NAME];

struct RARG
{
  int lana, nrecs;
  struct Variable *vars;
  struct Control *globs;
};

#ifdef AFTERBURNER
int stdin_is_tty = 0;  /* true if stdin  is character device */
int stdout_is_tty = 0; /* true if stdout is character device */
#endif

static bool lstdout = true;

static int Source = 0;

static int ofiletype = -1;

static int DataType = -1;

static char *filename;
static char **ifiles;
static char *ifile = NULL;
static char *ofile = NULL;
#ifdef AFTERBURNER
static char *ofile2 = NULL;
#endif

static int ofileidx = 0;

static int specGridID = -1;
static int gaussGridID = -1;
static int iVertID = -1;
static int oVertID = -1;

static int Lhybrid2pressure = FALSE;

static int TsID;
static bool lparallelread = true;

#define TIMESTEP_INTERVAL -1
#define MONTHLY_INTERVAL 0
#define DAILY_INTERVAL 1
#define UNLIM_INTERVAL 2

#define MaxHours 24
static int nrqh;
static int hours[MaxHours + 1];

static double *LevelFound;

static void
cdiError(int cdiErrno, const char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);

  printf("\n");
  vfprintf(stderr, fmt, args);
  fprintf(stderr, "\n");

  va_end(args);

  fprintf(stderr, "%s\n", cdiStringError(cdiErrno));

  if (_ExitOnError) exit(1);
}

static void
lprintf(FILE *fp)
{
  int num = 67;
  int cval = '-';

  fprintf(fp, " ");
  for (int inum = 0; inum < num; inum++) fprintf(fp, "%c", cval);
  fprintf(fp, "\n");
}

static void
FreeMean(struct Variable *vars)
{
  for (int code = 0; code < MaxCodes; code++)
    if (vars[code].mean)
      {
        Free(vars[code].mean);
        vars[code].mean = NULL;
      }
}

static void
after_PostProcess(struct Control *globs)
{
  if (globs->EndOfInterval)
    {
      if (lstdout)
        {
          if (globs->OutputInterval == DAILY_INTERVAL)
            fprintf(stdout, " Processed Day %2d  Month %2d  Year %04d", globs->OldDate.dy, globs->OldDate.mo, globs->OldDate.yr);
          else if (globs->OutputInterval == MONTHLY_INTERVAL)
            fprintf(stdout, " Processed Month %2d  Year %04d", globs->OldDate.mo, globs->OldDate.yr);
          else if (globs->OutputInterval == UNLIM_INTERVAL)
            fprintf(stdout, " Processed range from %6.4d-%2.2d-%2.2d to %6.4d-%2.2d-%2.2d", globs->StartDate.yr,
                    globs->StartDate.mo, globs->StartDate.dy, globs->OldDate.yr, globs->OldDate.mo, globs->OldDate.dy);

          if (globs->Mean)
            fprintf(stdout, "  (Mean of %3d Terms)\n", globs->MeanCount);
          else
            fprintf(stdout, "   Terms %3d\n", globs->MeanCount);
        }

      globs->EndOfInterval = FALSE;
      globs->MeanCount = 0;
    }
}

/* ================= */
/* switch input file */
/* ================= */
static void
after_SwitchFile(struct Control *globs)
{
  bool echam4 = false;
  int n;
  char y3, y2, y1, y0;
  char m1, m0;
  char d1, d0;

  streamClose(globs->istreamID);

  if (globs->Multi > 0)
    {
      int i = strlen(ifile);
      if (i < 10)
        {
          fprintf(stderr, " Not a valid filename: %s \n", ifile);
          exit(1);
        }

      if (ifile[i - 3] == '.')
        {
          echam4 = true;
          y3 = ifile[i - 9];
          y2 = ifile[i - 8];
          y1 = ifile[i - 7];
          y0 = ifile[i - 6];
          m1 = ifile[i - 5];
          m0 = ifile[i - 4];
          d1 = ifile[i - 2];
          d0 = ifile[i - 1];
        }
      else
        {
          y3 = ifile[i - 6];
          y2 = ifile[i - 5];
          y1 = ifile[i - 4];
          y0 = ifile[i - 3];
          m1 = ifile[i - 2];
          m0 = ifile[i - 1];
          d1 = '0';
          d0 = '1';
        }

      for (n = 0; n < globs->DayIn; n++)
        {
          if (d0 == '9')
            {
              d0 = '0';
              d1++;
            }
          else
            d0++;
          if (d1 == '3' && d0 > '0')
            {
              d1 = '0';
              d0 = '1';
              if (m1 == '0')
                {
                  if (m0 == '9')
                    {
                      m0 = '0';
                      m1 = '1';
                    }
                  else
                    m0++;
                }
              else
                {
                  if (m0 < '2')
                    m0++;
                  else
                    {
                      m1 = '0';
                      m0 = '1';
                      y0++;
                      if (y0 > '9')
                        {
                          y0 = '0';
                          y1++;
                        }
                      if (y1 > '9')
                        {
                          y1 = (char) '0';
                          if (isdigit((int) y2))
                            y2++;
                          else
                            y2 = '1';
                          if (y2 > '9')
                            {
                              y2 = (char) '0';
                              if (isdigit((int) y3))
                                y3++;
                              else
                                y3 = '1';
                            }
                        }
                    }
                }
            }
        }

      if (echam4)
        {
          ifile[i - 9] = y3;
          ifile[i - 8] = y2;
          ifile[i - 7] = y1;
          ifile[i - 6] = y0;
          ifile[i - 5] = m1;
          ifile[i - 4] = m0;
          ifile[i - 2] = d1;
          ifile[i - 1] = d0;
        }
      else
        {
          ifile[i - 6] = y3;
          ifile[i - 5] = y2;
          ifile[i - 4] = y1;
          ifile[i - 3] = y0;
          ifile[i - 2] = m1;
          ifile[i - 1] = m0;
        }

      globs->Multi--;
    }

  if (globs->Nfiles > 0) ifile = ifiles[--globs->Nfiles];

  fprintf(stderr, " Continuation file: %s\n", ifile);

  globs->istreamID = streamOpenRead(ifile);
  if (globs->istreamID < 0) cdiError(globs->istreamID, "Open failed on %s", ifile);

  globs->ivlistID = streamInqVlist(globs->istreamID);
  globs->taxisID = vlistInqTaxis(globs->ivlistID);
}

static int64_t
after_getDate(struct Date datetime)
{
  return cdiEncodeDate(datetime.yr, datetime.mo, datetime.dy);
}

static int
after_getTime(struct Date datetime)
{
  return cdiEncodeTime(datetime.hr, datetime.mn, 0);
}

static void
after_setDateTime(struct Date *datetime, int64_t date, int time)
{
  int sec;
  cdiDecodeDate(date, &datetime->yr, &datetime->mo, &datetime->dy);
  cdiDecodeTime(time, &datetime->hr, &datetime->mn, &sec);
}

static void
after_printProcessStatus(int tsID)
{
  static bool counthead = false;

  if (tsID == -1)
    {
      if (stdout_is_tty)
        {
          fprintf(stdout, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
          fflush(stdout);
        }

      counthead = false;
    }
  else
    {
      if (counthead == false)
        {
          if (stdout_is_tty) fprintf(stdout, " Process timestep :       ");

          counthead = true;
        }

      if (stdout_is_tty)
        {
          fprintf(stdout, "\b\b\b\b\b\b%6d", tsID);
          fflush(stdout);
        }
    }
}

static int
after_setNextDate(struct Control *globs)
{
  int nrecs = 0;

  bool righttime = false;
  while (TRUE)
    {
      nrecs = streamInqTimestep(globs->istreamID, TsID);
      if (nrecs == 0 && (globs->Multi > 0 || globs->Nfiles > 0))
        {
          if (lstdout) after_printProcessStatus(-1);

          after_SwitchFile(globs);

          if (globs->istreamID >= 0)
            {
              TsID = 0;
              nrecs = streamInqTimestep(globs->istreamID, TsID);
            }
        }
      if (nrecs == 0) break;

      int64_t vdate = taxisInqVdate(globs->taxisID);
      int vtime = taxisInqVtime(globs->taxisID);
      after_setDateTime(&globs->NextDate, vdate, vtime);

      for (int i = 0; i < nrqh; i++)
        if (hours[i] < 0 || hours[i] == globs->NextDate.hr)
          {
            righttime = true;
            break;
          }

      if (righttime)
        break;
      else
        TsID += 1;
    }

  return nrecs;
}

static int num_recs = 0;

static void *
after_readTimestep(void *arg)
{
  int varID, gridID, zaxisID, levelID, timeID;
  size_t nmiss;
  RARG *rarg = (RARG *) arg;

  int nrecs = rarg->nrecs;
  int analysisData = rarg->lana;
  struct Variable *vars = rarg->vars;
  struct Control *globs = rarg->globs;

  for (int code = 0; code < MaxCodes; code++) vars[code].nmiss0 = 0;

  for (int recID = 0; recID < nrecs; recID++)
    {
      streamInqRecord(globs->istreamID, &varID, &levelID);

      int code = vlistInqVarCode(globs->ivlistID, varID);
      if (code <= 0 || code >= MaxCodes) continue;

      /* Skip records containing unneeded codes */

      if (!vars[code].needed0) continue;

      vlistInqVar(globs->ivlistID, varID, &gridID, &zaxisID, &timeID);

      int leveltype = zaxisInqType(zaxisID);

      /* Skip records with unselected levels */

      int levelOffset = -1;
      /*
        if ( vars[code].ozaxisID != vars[code].izaxisID && ! Lhybrid2pressure )
      */
      if ((vars[code].ozaxisID != vars[code].izaxisID) && (leveltype == ZAXIS_PRESSURE))
        {
          int level = (int) zaxisInqLevel(zaxisID, levelID);
          for (int i = 0; i < globs->NumLevelRequest; ++i)
            {
              if (IS_EQUAL(globs->LevelRequest[i], level))
                {
                  levelOffset = i;
                  break;
                }
            }

          if (levelOffset < 0) continue;

          zaxisID = vars[code].ozaxisID;
          levelID = levelOffset;
        }

      if (globs->Debug)
        {
          fprintf(stderr, "T%d", globs->Truncation);

          fprintf(stderr, "  Code %3d   Level%6d   %6.4d-%2.2d-%2.2d  %2.2d:%2.2d:00\n", code,
                  (int) zaxisInqLevel(zaxisID, levelID), globs->OldDate.yr, globs->OldDate.mo, globs->OldDate.dy, globs->OldDate.hr,
                  globs->OldDate.mn);
        }

      if (analysisData)
        {
          streamReadRecord(globs->istreamID, globs->Field, &nmiss);
          after_AnalysisAddRecord(globs, vars, code, gridID, zaxisID, levelID, nmiss);
        }
      else
        {
          double *dataptr = after_get_dataptr(vars, code, gridID, zaxisID, levelID);
          streamReadRecord(globs->istreamID, dataptr, &nmiss);
          after_EchamAddRecord(globs, vars, code, gridID, zaxisID, levelID, nmiss);
        }

      if (iVertID != -1 && oVertID != -1 && (vars[code].izaxisID == iVertID)) vars[code].ozaxisID = oVertID;
    }

  TsID++;
  /*
    printf("%3d  date = %d  time = %04d\n", TsID, vdate, vtime);
  */
  num_recs = after_setNextDate(globs);

  return (void *) &num_recs;
}

static void
after_defineNextTimestep(struct Control *globs)
{
  static int otsID = 0;
  int64_t vdate = after_getDate(globs->OldDate);
  int vtime = after_getTime(globs->OldDate);
  taxisDefVdate(globs->taxisID2, vdate);
  taxisDefVtime(globs->taxisID2, vtime);

  if (globs->Mean != 2)
    {
      if (otsID == 0)
        {
          int nvars = vlistNvars(globs->ovlistID);
          if (nvars == 0) Error("No variable selected!");
          vlistDefTaxis(globs->ovlistID, globs->taxisID2);
          streamDefVlist(globs->ostreamID, globs->ovlistID);
        }
      taxisDefNumavg(globs->taxisID2, globs->MeanCount + 1);
      streamDefTimestep(globs->ostreamID, otsID);
    }

#ifdef AFTERBURNER
  if (globs->Mean >= 2)
    {
      if (otsID == 0)
        {
          vlistDefTaxis(globs->ovlistID2, globs->taxisID2);
          streamDefVlist(globs->ostreamID2, globs->ovlistID2);
        }
      taxisDefNumavg(globs->taxisID2, globs->MeanCount + 1);
      streamDefTimestep(globs->ostreamID2, otsID);
    }
#endif

  otsID++;
}

static void
after_setEndOfInterval(struct Control *globs, int nrecs)
{
  if (nrecs == 0)
    {
      globs->EndOfInterval = TRUE;
    }
  else
    {
      if (globs->OutputInterval == DAILY_INTERVAL)
        globs->EndOfInterval = globs->NewDate.dy != globs->OldDate.dy;
      else if (globs->OutputInterval == MONTHLY_INTERVAL)
        globs->EndOfInterval = globs->NewDate.mo != globs->OldDate.mo;
      else if (globs->OutputInterval == UNLIM_INTERVAL)
        globs->EndOfInterval = FALSE;
      else
        Error("output interval %d not implemented!", globs->OutputInterval);
    }
}

static void
after_moveTimestep(struct Variable *vars)
{
  int code;

  for (code = 0; code < MaxCodes; code++) vars[code].nmiss = vars[code].nmiss0;

  for (code = 0; code < MaxCodes; code++)
    if (vars[code].hybrid0)
      {
        vars[code].hybrid = vars[code].hybrid0;
        vars[code].hybrid0 = NULL;
      }

  for (code = 0; code < MaxCodes; code++)
    if (vars[code].spectral0)
      {
        vars[code].spectral = vars[code].spectral0;
        vars[code].spectral0 = NULL;
      }

  for (code = 0; code < MaxCodes; code++)
    if (vars[code].grid0)
      {
        vars[code].grid = vars[code].grid0;
        vars[code].grid0 = NULL;
      }
}

static void
after_check_content(struct Variable *vars, int timestep)
{
  extern int labort_after;
  for (int code = 0; code < 272; code++)
    {
      /*  if ( code == GEOPOTENTIAL ) continue; */
      if (code == SLP) continue;
      if (code == GEOPOTHEIGHT) continue;
      if (code == STREAM) continue;
      if (code == VELOPOT) continue;
      if (code == U_WIND) continue;
      if (code == V_WIND) continue;
      if (code == OMEGA) continue;
      if (code == RHUMIDITY) continue;
      if (code == LOW_CLOUD) continue;
      if (code == MID_CLOUD) continue;
      if (code == HIH_CLOUD) continue;
      if (code == PS) continue;
      if (code == HUMIDITY)
        {
          if (vars[code].needed && !vars[code].selected && vars[code].spectral == NULL && vars[code].hybrid == NULL)
            {
              static bool lwarn = true;
              if (lwarn) Warning("No humidity in data file, set to zero !");
              lwarn = false;
              vars[code].needed = FALSE;
            }
        }
      else
        {
          if (vars[code].needed && !vars[code].comp && vars[code].spectral == NULL && vars[code].hybrid == NULL)
            {
              if (labort_after)
                Error("Code  %3d not found at timestep %d!", code, timestep);
              else
                Warning("Code  %3d not found at timestep %d!", code, timestep);
            }
        }
    }
  /*
  if ( NumLevelRequest > 0 )
    {
      vars[HALF_PRESS].needed = 1;
      vars[FULL_PRESS].needed = 1;
    }

  code = HALF_PRESS;
  if ( vars[code].needed && !vars[code].comp &&
       vars[code].spectral == NULL && vars[code].hybrid == NULL )
    Error( "Hybrid model level not found!");

  code = FULL_PRESS;
  if ( vars[code].needed && !vars[code].comp &&
       vars[code].spectral == NULL && vars[code].hybrid == NULL )
    Error( "Hybrid model level not found!");
  */
}

static void
after_control(struct Control *globs, struct Variable *vars)
{
  int i;
  int nrecs;
  int64_t rdate, vdate;
  int rtime, vtime;
  int code;
  RARG rarg;
  void *statusp = NULL;
  void *read_task = NULL;

  if (lparallelread)
    {
      read_task = cdo_task_new();
      if (read_task == NULL)
        {
          lparallelread = false;
          cdoWarning("CDO tasks not available!");
        }
    }

  for (code = 0; code < MaxCodes; code++) vars[code].needed0 = vars[code].needed;

  TsID = 0;

  bool righttime = false;
  while ((nrecs = streamInqTimestep(globs->istreamID, TsID)) > 0)
    {
      vdate = taxisInqVdate(globs->taxisID);
      vtime = taxisInqVtime(globs->taxisID);
      after_setDateTime(&globs->StartDate, vdate, vtime);
      after_setDateTime(&globs->NewDate, vdate, vtime);

      for (i = 0; i < nrqh; i++)
        if (hours[i] < 0 || hours[i] == globs->NewDate.hr)
          {
            righttime = true;
            break;
          }

      if (righttime)
        break;
      else
        TsID++;
    }

  if (taxisInqType(globs->taxisID) == TAXIS_RELATIVE)
    {
      rdate = taxisInqRdate(globs->taxisID);
      rtime = taxisInqRtime(globs->taxisID);
    }
  else
    {
      rdate = after_getDate(globs->StartDate);
      rtime = after_getTime(globs->StartDate);
    }

  if (ofiletype == CDI_FILETYPE_NC || ofiletype == CDI_FILETYPE_NC2 || ofiletype == CDI_FILETYPE_NC4
      || ofiletype == CDI_FILETYPE_NC4C || ofiletype == CDI_FILETYPE_NC5)
    {
      taxisDefCalendar(globs->taxisID2, CALENDAR_PROLEPTIC);
      taxisDefType(globs->taxisID2, TAXIS_RELATIVE);
      taxisDefTunit(globs->taxisID2, TUNIT_DAY);
      taxisDefRdate(globs->taxisID2, rdate);
      taxisDefRtime(globs->taxisID2, rtime);
    }

  globs->OldDate = globs->NewDate;

  bool tsFirst = true;

  while (nrecs > 0)
    {
      rarg.nrecs = nrecs;
      rarg.lana = globs->AnalysisData;
      rarg.vars = vars;
      rarg.globs = globs;

      if (tsFirst || lparallelread == false)
        {
          if (lparallelread == false)
            {
              statusp = after_readTimestep(&rarg);
            }
          else
            {
              cdo_task_start(read_task, after_readTimestep, &rarg);
            }

          if (tsFirst && globs->Type > 0) after_legini_setup(globs, vars);

          if (lparallelread)
            {
              statusp = cdo_task_wait(read_task);
              if (*(int *) statusp < 0) Error("after_readTimestep error! (status = %d)", *(int *) statusp);
            }
          tsFirst = false;
        }
      else
        {
          statusp = cdo_task_wait(read_task);
          if (*(int *) statusp < 0) Error("after_readTimestep error! (status = %d)", *(int *) statusp);
        }

      nrecs = *(int *) statusp;

      globs->MeanCount0 = globs->MeanCount;
      globs->NewDate = globs->NextDate;

      after_moveTimestep(vars);

      if (nrecs && lparallelread)
        {
          cdo_task_start(read_task, after_readTimestep, &rarg);
        }

      after_setEndOfInterval(globs, nrecs);

      if (lstdout) after_printProcessStatus(TsID);

      if (lstdout && globs->EndOfInterval) after_printProcessStatus(-1);

      if (globs->Mean == 0 || globs->EndOfInterval)
        {
          if (!globs->AnalysisData) after_check_content(vars, globs->TermCount + 1);
          after_defineNextTimestep(globs);
        }

      if (globs->AnalysisData)
        after_processPL(globs, vars);
      else
        after_processML(globs, vars);

      after_PostProcess(globs);

      if (nrecs)
        {
          if (globs->AnalysisData)
            after_AnalysisDependencies(vars, MaxCodes);
          else
            after_EchamDependencies(vars, MaxCodes, globs->Type, Source);
        }

      globs->OldDate = globs->NewDate;
    }

  if (read_task) cdo_task_delete(read_task);
}

static void
after_setLevel(struct Control *globs)
{
  int k, l, found;
  int removeLevel[MaxLevel];
  double level;
  bool checkLevel = true;
  int numplevelDefault; /* default pressure level */
  long plevelDefault[]
      = { 100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000, 20000, 15000, 10000, 7000, 5000, 3000, 2000, 1000 };
  int numhlevelDefault; /* default height level */
  long hlevelDefault[] = { 0, 1000, 2000, 5000, 10000, 15000, 20000, 25000, 30000 };

  numplevelDefault = sizeof(plevelDefault) / sizeof(plevelDefault[0]);
  numhlevelDefault = sizeof(hlevelDefault) / sizeof(hlevelDefault[0]);

  if (iVertID != -1)
    if (zaxisInqType(iVertID) == ZAXIS_HYBRID && globs->Type > 20) Lhybrid2pressure = TRUE;

  if (globs->Verbose) lprintf(stdout);

  if (globs->NumLevelRequest == 0)
    {
      if (iVertID == -1)
        {
          if (globs->Verbose) fprintf(stdout, " No level detected\n");
        }
      else
        {
          if (Lhybrid2pressure)
            {
              if (globs->unitsel == 0)
                {
                  if (globs->Verbose) fprintf(stdout, " Default pressure level selected:\n");
                  globs->NumLevelRequest = numplevelDefault;
                  for (l = 0; l < globs->NumLevelRequest; l++) globs->LevelRequest[l] = plevelDefault[l];
                  oVertID = zaxisCreate(ZAXIS_PRESSURE, globs->NumLevelRequest);
                  zaxisDefLevels(oVertID, globs->LevelRequest);
                }
              else
                {
                  if (globs->Verbose) fprintf(stdout, " Default height level selected:\n");
                  globs->NumLevelRequest = numhlevelDefault;
                  for (l = 0; l < globs->NumLevelRequest; l++) globs->LevelRequest[l] = hlevelDefault[l];
                  oVertID = zaxisCreate(ZAXIS_HEIGHT, globs->NumLevelRequest);
                  zaxisDefLevels(oVertID, globs->LevelRequest);
                }
            }
          else
            {
              if (globs->Verbose)
                {
                  if (zaxisInqType(iVertID) == ZAXIS_HYBRID)
                    fprintf(stdout, " All detected hybrid level selected:\n");
                  else
                    fprintf(stdout, " All detected pressure level selected:\n");
                }
              globs->NumLevelRequest = globs->NumLevelFound;
              for (l = 0; l < globs->NumLevelRequest; l++) globs->LevelRequest[l] = LevelFound[l];
              oVertID = iVertID;
            }
        }
      checkLevel = false;
    }
  else
    {
      if (iVertID == -1)
        {
          if (globs->Verbose) fprintf(stdout, " No level detected\n");
          checkLevel = false;
        }
      else if (globs->NumLevelRequest == 1 && IS_EQUAL(globs->LevelRequest[0], 0))
        {
          if (globs->Verbose) fprintf(stdout, " No level selected\n");
          globs->NumLevelRequest = 0;
          checkLevel = false;
        }
      else if (globs->Verbose)
        {
          if (Lhybrid2pressure)
            {
              if (globs->unitsel == 0)
                fprintf(stdout, " Selected pressure level:\n");
              else
                fprintf(stdout, " Selected height level:\n");
            }
          else
            {
              if (zaxisInqType(iVertID) == ZAXIS_HYBRID)
                fprintf(stdout, " Selected hybrid level:\n");
              else
                {
                  if (globs->unitsel == 0)
                    fprintf(stdout, " Selected pressure level:\n");
                  else
                    fprintf(stdout, " Selected height level:\n");
                }
            }
        }
    }

  if (globs->Verbose && iVertID != -1)
    for (l = 0; l < globs->NumLevelRequest; l++) fprintf(stdout, "  Level %2d = %13.4f\n", l + 1, globs->LevelRequest[l]);

  if (checkLevel)
    {
      for (k = 0; k < globs->NumLevelRequest; k++) removeLevel[k] = FALSE;
      for (k = 0; k < globs->NumLevelRequest; k++)
        {
          level = globs->LevelRequest[k];
          for (l = k + 1; l < globs->NumLevelRequest; l++)
            if (removeLevel[l] == FALSE && IS_EQUAL(level, globs->LevelRequest[l]))
              {
                if (globs->Verbose) fprintf(stdout, "  Level %2d = %13.4f double request\n", l + 1, globs->LevelRequest[l]);
                removeLevel[l] = TRUE;
              }
        }

      l = 0;
      for (k = 0; k < globs->NumLevelRequest; k++)
        if (removeLevel[k] == FALSE) globs->LevelRequest[l++] = globs->LevelRequest[k];

      globs->NumLevelRequest = l;

      if (globs->AnalysisData || globs->Type < 30)
        {
          for (k = 0; k < globs->NumLevelRequest; k++) removeLevel[k] = FALSE;
          for (k = 0; k < globs->NumLevelRequest; k++)
            {
              level = globs->LevelRequest[k];
              found = FALSE;
              for (l = 0; l < globs->NumLevelFound; l++)
                if (IS_EQUAL(level, LevelFound[l])) found = TRUE;

              if (!found)
                {
                  fprintf(stdout, "  Level %2d = %14.4f not in input\n", k + 1, globs->LevelRequest[k]);
                  removeLevel[k] = TRUE;
                }
            }

          l = 0;
          for (k = 0; k < globs->NumLevelRequest; k++)
            if (removeLevel[k] == FALSE) globs->LevelRequest[l++] = globs->LevelRequest[k];

          if (l != globs->NumLevelRequest)
            {
              extern int labort_after;
              if (globs->Verbose) lprintf(stdout);
              if (labort_after)
                Error("Inconsistent or invalid level list!");
              else
                Warning("Inconsistent or invalid level list!");
            }

          globs->NumLevelRequest = l;
        }
    }

  if (globs->Verbose) lprintf(stdout);
}

static void
after_defineLevel(struct Control *globs, struct Variable *vars)
{
  int code, i;

  /* hybrid, pressure, height */

  switch (globs->Type)
    {
    case 0:
    case 10:
    case 11:
    case 20:
      {
        if (iVertID == -1) break;

        if (zaxisInqType(iVertID) == ZAXIS_HYBRID)
          {
            if (oVertID == -1)
              {
                if (globs->NumLevelRequest > globs->NumLevelFound) Error("Too much level requested");

                if (globs->NumLevelFound == globs->NumLevelRequest)
                  {
                    for (i = 0; i < globs->NumLevelRequest; i++)
                      if (IS_NOT_EQUAL(globs->LevelRequest[i], LevelFound[i])) break;

                    if (i == globs->NumLevelRequest) oVertID = iVertID;
                  }

                if (oVertID == -1 && globs->NumLevelRequest > 0)
                  {
                    oVertID = zaxisCreate(ZAXIS_HYBRID, globs->NumLevelRequest);
                    zaxisDefLevels(oVertID, globs->LevelRequest);
                    zaxisDefVct(oVertID, globs->nvct, globs->vct);
                  }
              }

            for (code = 0; code < MaxCodes; code++)
              {
                if (vars[code].selected)
                  {
                    if (vars[code].izaxisID != -1)
                      if (zaxisInqType(vars[code].izaxisID) == ZAXIS_HYBRID
                          && zaxisInqSize(vars[code].izaxisID) >= globs->NumLevelRequest)
                        vars[code].ozaxisID = oVertID;
                  }
              }
          }
        else
          {
            zaxisName(zaxisInqType(iVertID), zaxistypename);
            Error("%s level data unsupported for TYPE %d", zaxistypename, globs->Type);
          }
        break;
      }
    case 30:
    case 40:
    case 41:
    case 50:
    case 60:
    case 61:
    case 70:
      {
        if (iVertID == -1) break;

        if (oVertID == -1)
          {
            if (globs->unitsel == 0)
              oVertID = zaxisCreate(ZAXIS_PRESSURE, globs->NumLevelRequest);
            else
              oVertID = zaxisCreate(ZAXIS_HEIGHT, globs->NumLevelRequest);

            zaxisDefLevels(oVertID, globs->LevelRequest);
          }

        for (code = 0; code < MaxCodes; code++)
          {
            if (vars[code].selected)
              {
                if (vars[code].izaxisID != -1)
                  {
                    int nlev = zaxisInqSize(vars[code].izaxisID);
                    if (zaxisInqType(vars[code].izaxisID) == zaxisInqType(iVertID)
                        && (nlev == globs->NumLevel || nlev == globs->NumLevel + 1) && nlev > 1)
                      vars[code].ozaxisID = oVertID;
                  }
              }
          }

        break;
      }
    default: Error("TYPE %d unsupported", globs->Type);
    }
}

static void
after_defineGrid(struct Control *globs, struct Variable *vars)
{
  int ogridID = -1;
  int code;

  /* spectral, fourier, gauss, zonal mean */

  switch (globs->Type)
    {
    case 0:
    case 50:
      {
        if (specGridID == -1)
          {
            if (globs->DimSP == 0) Error("dim spectral undefined");
            if (globs->Truncation == 0) Error("truncation undefined");

            specGridID = gridCreate(GRID_SPECTRAL, globs->DimSP);
            gridDefTrunc(specGridID, globs->Truncation);
          }

        ogridID = specGridID;
        break;
      }
    case 20:
    case 30:
    case 70:
      {
        if (gaussGridID == -1)
          {
            if (globs->Longitudes == 0) Error("number of longitudes undefined");
            if (globs->Latitudes == 0) Error("number of latitudes undefined");

            gaussGridID = gridCreate(GRID_GAUSSIAN, globs->Longitudes * globs->Latitudes);
            gridDefXsize(gaussGridID, globs->Longitudes);
            gridDefYsize(gaussGridID, globs->Latitudes);
          }

        ogridID = gaussGridID;
        break;
      }
    case 10:
    case 40:
    case 60:
      {
        if (globs->Fouriers == 0) Error("number of fourier coefficients undefined");
        if (globs->Latitudes == 0) Error("number of latitudes undefined");

        ogridID = gridCreate(GRID_FOURIER, globs->Fouriers * globs->Latitudes);
        gridDefXsize(ogridID, globs->Latitudes);
        gridDefYsize(ogridID, globs->Fouriers);
        break;
      }
    case 11:
    case 41:
    case 61:
      {
        if (globs->Latitudes == 0) Error("Number of latitudes undefined");

        ogridID = gridCreate(GRID_GAUSSIAN, globs->Latitudes);
        gridDefXsize(ogridID, 1);
        gridDefYsize(ogridID, globs->Latitudes);
        break;
      }
    default: Error("TYPE %d unsupported", globs->Type);
    }

  if (ogridID != -1)
    for (code = 0; code < MaxCodes; code++)
      {
        if (vars[code].selected)
          {
            vars[code].ogridID = ogridID;
          }
      }

  if (ogridID == -1) Error("out grid undefined");
}

static void
after_setCodes(struct Control *globs, struct Variable *vars, int maxCodes, int numCodes)
{
  if (globs->Verbose) lprintf(stdout);

  if (numCodes == 0)
    {
      if (globs->Verbose) fprintf(stdout, " All detected codes selected:\n");

      for (int code = 0; code < maxCodes; code++)
        if (vars[code].detected) vars[code].selected = 1;
    }
  else if (globs->Verbose)
    fprintf(stdout, " Selected codes:\n");

  if (globs->Verbose)
    {
      fprintf(stdout, "  Table Code Name              Longname\n");
      fprintf(stdout, "  ----- ---- ----              --------\n");
    }

  for (int code = 0; code < maxCodes; code++)
    if (vars[code].selected)
      {
        char name[CDI_MAX_NAME];
        name[0] = 0;
        char longname[CDI_MAX_NAME];
        longname[0] = 0;
        int tableID;
        int table = 0;
        int varID = vars[code].ivarID;

        if (varID == CDI_UNDEFID)
          {
            int modelID = vlistInqVarModel(globs->ivlistID, 0);
            table = 128;
            tableID = tableInq(modelID, table, NULL);

            vars[code].tableID = tableID;
          }
        else
          {
            tableID = vlistInqVarTable(globs->ivlistID, varID);
            table = tableInqNum(tableID);
            vlistInqVarName(globs->ivlistID, varID, name);
            vlistInqVarLongname(globs->ivlistID, varID, longname);
          }

        if (!name[0]) tableInqEntry(tableID, code, -1, name, longname, NULL);

        if (globs->Verbose)
          {
            fprintf(stdout, " %5d", table);
            fprintf(stdout, " %4d", code);
            if (!name[0])
              fprintf(stdout, "  var%d", code);
            else
              {
                fprintf(stdout, "  %-16s", name);
                if (longname[0]) fprintf(stdout, "  %s", longname);
              }
            fprintf(stdout, "\n");
          }
      }
}

static void
after_checkNamelist(struct Control *globs)
{
  if (globs->Mean && globs->Type < 20)
    {
      Error("Mean is only available for TYPE >= 20!");
    }

  if (globs->Extrapolate == FALSE && globs->Type >= 30)
    {
      if (globs->Type > 30) Error("EXTRAPOLATE = 0 is only available for TYPE = 30!");
      if (globs->Mean) Error("EXTRAPOLATE = 0 is only available with MEAN = 0!");
    }
}

#ifdef AFTERBURNER
static void
after_usage(void)
{
  fprintf(stderr, "\nafter [options] <InputFiles> <OutputFile> <VarianceFile>\n");
#if defined(_OPENMP)
  fprintf(stderr, "     option -P <nthreads> : Set number of OpenMP threads\n");
#endif
  fprintf(stderr, "     option -a            : Forces analysis data process\n");
  fprintf(stderr, "     option -c            : Print available codes and names\n");
  fprintf(stderr, "     option -d            : Debug mode\n");
  fprintf(stderr, "     option -v <vctfile>  : Read vct from vctfile\n");
  /*  fprintf(stderr, "     option -h : help (this output)\n"); */
  /*  fprintf(stderr, "     option -p : parallel read on\n"); */
  fprintf(stderr, "  <InputFiles> : ECHAM or ECMWF Ana or ReAna files\n");
  fprintf(stderr, "  <OutputFile> : GRIB, NetCDF or SERVICE format file\n");
  fprintf(stderr, "<VarianceFile> : GRIB, NetCDF or SERVICE format file\n");
  fprintf(stderr, "  namelist is read from <stdin>\n");
  fprintf(stderr, "  output is written to <stdout>\n\n");

  fprintf(stderr, "  default Namelist: \n");
  fprintf(stderr, "  &SELECT\n");
  fprintf(stderr, "    TYPE = 0, CODE = -1, LEVEL = -1, MULTI = 0, DAYIN = 30,\n");
  fprintf(stderr, "    MEAN = 0, TIMESEL = -1, UNITSEL = 0,\n");
  fprintf(stderr, "    FORMAT = 0, PRECISION = 0, SZIP = 0\n");
  fprintf(stderr, "  &END\n");

  exit(1);
}
#endif

static void
after_parini(struct Control *globs, struct Variable *vars)
{
  char namelist[65536];

  if (stdin_is_tty)
    {
#ifdef CDO
      fprintf(stderr, "Default namelist: \n");
      fprintf(stderr, "  TYPE=0, CODE=-1, LEVEL=-1, INTERVAL=0, MEAN=0, EXTRAPOLATE=0\n");
#endif
      fprintf(stdout, "Enter namelist parameter:\n");
    }
  else
    {
      fseek(stdin, 0L, SEEK_END);
      long length = ftell(stdin);
      if (length == 0L)
        {
          fprintf(stderr, "\n stdin not connected\n");
#ifdef AFTERBURNER
          after_usage();
#endif
        }
      fseek(stdin, 0L, SEEK_SET);
    }

  int i = 1;
  namelist[0] = ' ';
  int c = getchar();
  while ((c != EOF) && i < (int) (sizeof(namelist) - 1))
    {
      if ((c >= '0' && c <= '9') || (c == '-' || c == '.'))
        namelist[i++] = c;
      else if (c >= 'a' && c <= 'z')
        namelist[i++] = c;
      else if (c >= 'A' && c <= 'Z')
        namelist[i++] = tolower(c);
      else
        c = ' ';

      if (c == ' ' && namelist[i - 1] != ' ') namelist[i++] = c;
      c = getchar();
    }
  namelist[i] = 0;

  if (globs->Debug)
    {
      lprintf(stderr);
      fprintf(stderr, "  Length of namelist:%4d bytes\n", (int) strlen(namelist));

      for (i = 0; i < (int) strlen(namelist); i += 60) fprintf(stderr, "  namelist[%02d]=%-60.60s\n", i, namelist + i);
      lprintf(stderr);
    }

  if (globs->Verbose)
    {
      lprintf(stdout);
      fprintf(stdout, " Namelist:\n");
    }

  globs->Type = scan_par(globs->Verbose, namelist, "type", 0);
  globs->Multi = scan_par(globs->Verbose, namelist, "multi", 0);
  globs->Mean = scan_par(globs->Verbose, namelist, "mean", 0);
  globs->OutputInterval = scan_par(globs->Verbose, namelist, "interval", MONTHLY_INTERVAL);

#ifdef CDO
  if (globs->Mean >= 2) cdoAbort("Namelist parameter MEAN=%d out of bounds (0:1)", globs->Mean);
#endif

  int fileFormat = scan_par(globs->Verbose, namelist, "format", -1);
  int gribFormat = scan_par_obsolate(namelist, "grib", 0);
  int cdfFormat = scan_par_obsolate(namelist, "netcdf", 0);

  if (gribFormat && cdfFormat) Error("GRIB or NetCDF?");

  switch (fileFormat)
    {
#ifdef CDO
    case -1: ofiletype = -1; break;
#else
    case -1: ofiletype = CDI_FILETYPE_SRV; break;
#endif
    case 0: ofiletype = CDI_FILETYPE_SRV; break;
    case 1: ofiletype = CDI_FILETYPE_GRB; break;
    case 2: ofiletype = CDI_FILETYPE_NC; break;
    case 3: ofiletype = CDI_FILETYPE_EXT; break;
    case 4: ofiletype = CDI_FILETYPE_NC2; break;
    case 5: ofiletype = CDI_FILETYPE_NC5; break;
    case 6: ofiletype = CDI_FILETYPE_NC4; break;
    default: Error("unknown file format %d", fileFormat);
    }

  if (gribFormat) ofiletype = CDI_FILETYPE_GRB;
  if (cdfFormat) ofiletype = CDI_FILETYPE_NC;

  int precision = scan_par(globs->Verbose, namelist, "precision", 0);
  if (precision) switch (precision)
      {
      case 8: DataType = CDI_DATATYPE_PACK8; break;
      case 16: DataType = CDI_DATATYPE_PACK16; break;
      case 24: DataType = CDI_DATATYPE_PACK24; break;
      case 32: DataType = CDI_DATATYPE_FLT32; break;
      case 64: DataType = CDI_DATATYPE_FLT64; break;
      default: Error("unsupported data precision %d", precision);
      }

  globs->unitsel = scan_par(globs->Verbose, namelist, "unitsel", 0);
  globs->DayIn = scan_par(globs->Verbose, namelist, "dayinc", 30);
  globs->Extrapolate = scan_par(globs->Verbose, namelist, "extrapolate", 1);
  globs->Szip = scan_par(globs->Verbose, namelist, "szip", 0);
  int mars = scan_par_obsolate(namelist, "mars", 0);

  if (globs->Multi) --globs->Multi;

  if (mars)
    {
      extern int Mars;
      Mars = 1;
      PlanetRD = C_MARS_RD;
      PlanetGrav = C_MARS_GRAV;
      PlanetRadius = C_MARS_RADIUS;
    }

  nrqh = scan_time(globs->Verbose, namelist, hours, MaxHours);
  scan_code(namelist, vars, MaxCodes, &globs->NumCodesRequest);

  scan_darray(namelist, "level", globs->LevelRequest, MaxLevel, &globs->NumLevelRequest);
  if (globs->NumLevelRequest == 1)
    if (IS_EQUAL(globs->LevelRequest[0], -1)) globs->NumLevelRequest = 0;

  if (globs->Verbose) lprintf(stdout);

  after_checkNamelist(globs);
}

static void
after_dimcalc(struct Control *globs)
{
  if (globs->AnalysisData) globs->NumLevel = globs->NumLevelRequest;

  if (globs->Latitudes == 0)
    {
      globs->Latitudes = 2 * ((globs->Truncation * 3 + 3) / 4);
      if (globs->Truncation == 30) globs->Latitudes = 48;
    }

  if (globs->Longitudes == 0)
    {
      globs->Longitudes = globs->Latitudes * 2;
      if (globs->Truncation == 62) globs->Longitudes = 192;
    }

  globs->Waves = globs->Truncation + 1;
  globs->Fouriers = globs->Waves * 2;
  globs->DimSP = (globs->Truncation + 1) * (globs->Truncation + 2);
  globs->DimFC = globs->Latitudes * globs->Fouriers;
  globs->DimGP = globs->Latitudes * globs->Longitudes;
  globs->Dim3GP = globs->NumLevel * globs->DimGP;
  globs->Dim3FC = globs->NumLevel * globs->DimFC;
  globs->Dim3SP = globs->NumLevel * globs->DimSP;
  globs->HalfLevels = globs->NumLevel + 1;
  globs->DimSP_half = globs->DimSP / 2;

  if (globs->AnalysisData) fprintf(stdout, " Found Ana or Re-Ana Data\n");

  if (globs->Verbose)
    {
      fprintf(stdout, " Dimensions:\n");
      fprintf(stdout, "  Truncation        = %4d\n", globs->Truncation);
      fprintf(stdout, "  Levels            = %4d\n", globs->NumLevel);
      fprintf(stdout, "  Latitudes         = %4d\n", globs->Latitudes);
      fprintf(stdout, "  Longitudes        = %4d\n", globs->Longitudes);
      lprintf(stdout);
    }
}

/* ----------------------------------------------------------- */
/* Extract basic dimension information                         */
/* ----------------------------------------------------------- */
static void
after_precntl(struct Control *globs, struct Variable *vars)
{
  int l;
  int code = 0;
  int gridID, zaxisID, varID, timeID;
  int i, index, leveltype, gridtype;
  int datasize, numlevel;
  int vertfound = 0;
  int nhzaxis = 0;
  int FieldDim = 0;

  int nvars = vlistNvars(globs->ivlistID);
  int ngrids = vlistNgrids(globs->ivlistID);
  int nverts = vlistNzaxis(globs->ivlistID);
  int ntsteps = vlistNtsteps(globs->ivlistID);

  if (globs->Debug)
    {
      Message("nvars      = %d", nvars);
      Message("ngrids     = %d", ngrids);
      Message("nverts     = %d", nverts);
      Message("ntsteps    = %d", ntsteps);
    }

  for (index = 0; index < ngrids; index++)
    {
      gridID = vlistGrid(globs->ivlistID, index);
      gridtype = gridInqType(gridID);
      datasize = gridInqSize(gridID);

      if (datasize > FieldDim) FieldDim = datasize;

      if (gridtype == GRID_SPECTRAL && globs->Truncation == 0)
        {
          specGridID = gridID;
          globs->Truncation = gridInqTrunc(gridID);
        }
      else if (gridtype == GRID_GAUSSIAN && globs->Latitudes == 0)
        {
          gaussGridID = gridID;
          globs->Longitudes = gridInqXsize(gridID);
          globs->Latitudes = gridInqYsize(gridID);
        }
    }

  if (globs->Truncation == 0 && globs->Latitudes == 0) Error("Unsupported file structure (no spectral or Gaussian data found)!");

  if (globs->Truncation == 0)
    {
      if (globs->Latitudes)
        {
          switch (globs->Latitudes)
            {
            case 512: globs->Truncation = 511; break;
            case 320: globs->Truncation = 213; break;
            case 192: globs->Truncation = 127; break;
            case 160: globs->Truncation = 106; break;
            case 128: globs->Truncation = 85; break;
            case 96: globs->Truncation = 63; break;
            case 94: globs->Truncation = 62; break;
            case 64: globs->Truncation = 42; break;
            case 48: globs->Truncation = 31; break;
            case 32: globs->Truncation = 21; break;
            default: fprintf(stderr, "%d Gaussian latitudes not supported.\n", globs->Latitudes);
            }
        }
    }

  for (index = 0; index < nverts; index++)
    {
      zaxisID = vlistZaxis(globs->ivlistID, index);
      leveltype = zaxisInqType(zaxisID);
      numlevel = zaxisInqSize(zaxisID);
      /*
        printf("leveltype : %d %d\n", leveltype, zaxisInqSize(zaxisID));
      */
      if (numlevel > 1)
        {
          if (leveltype == ZAXIS_HYBRID || leveltype == ZAXIS_PRESSURE)
            {
              if (leveltype == ZAXIS_HYBRID && globs->nvct == 0)
                {
                  nhzaxis++;
                  int nvct = zaxisInqVctSize(zaxisID);
                  if (numlevel != (nvct / 2 - 1))
                    {
                      if (nvct == 0)
                        {
                          if (numlevel != 191) Warning("VCT missing for hybrid level data with %d levels!", numlevel);
                        }
                      else
                        {
                          Warning("Skip %d hybrid level data with %d levels!", (nvct / 2 - 1), numlevel);
                        }
                      continue;
                    }
                }
              else if (leveltype == ZAXIS_HYBRID && globs->nvct == zaxisInqVctSize(zaxisID))
                continue;

              if (iVertID != -1) Warning("More than %d different vertical grid structure found!", vertfound);

              vertfound++;

              if (iVertID != -1) continue;

              iVertID = zaxisID;
              globs->NumLevelFound = numlevel;
              LevelFound = (double *) Malloc(globs->NumLevelFound * sizeof(double));
              for (l = 0; l < globs->NumLevelFound; l++) LevelFound[l] = (int) zaxisInqLevel(zaxisID, l);

              if (leveltype == ZAXIS_HYBRID)
                {
                  if (globs->nvct == 0)
                    {
                      if (zaxisInqVctSize(zaxisID))
                        {
                          globs->nvct = zaxisInqVctSize(zaxisID);

                          if (globs->vct == NULL)
                            {
                              globs->vct = (double *) Malloc(globs->nvct * sizeof(double));
                              arrayCopy(globs->nvct, zaxisInqVctPtr(zaxisID), globs->vct);
                            }
                        }
                      else
                        {
                          Error("VCT not defined in inputfile!");
                        }
                    }

                  if (numlevel != (globs->nvct / 2 - 1))
                    Error("Number of hybrid levels %d does not match VCT levels %d", numlevel, globs->nvct / 2 - 1);

                  if (globs->Debug)
                    for (i = 0; i < globs->nvct / 2; i++)
                      fprintf(stderr, " vct: %4d %10.4f %10.4f\n", i, globs->vct[i], globs->vct[i + globs->nvct / 2]);
                }

              if (leveltype == ZAXIS_PRESSURE) globs->AnalysisData = TRUE;
            }
        }
    }

  if (nhzaxis > 0 && globs->nvct == 0) Error("VCT missing!");

  globs->NumLevel = globs->NumLevelFound;

  if (specGridID != -1) globs->Spectral = TRUE;
  if (gaussGridID != -1) globs->Gaussian = TRUE;

  if (globs->Debug) fprintf(stderr, "   T = %3d   L = %2d\n", globs->Truncation, globs->NumLevelFound);

  if (globs->Debug) fprintf(stderr, " CODE CHECK\n");

  if (globs->Verbose)
    {
      int instID = vlistInqVarInstitut(globs->ivlistID, 0);
      int modelID = vlistInqVarModel(globs->ivlistID, 0);

      lprintf(stdout);
      fprintf(stdout, " Institute : ");
      if (instID == CDI_UNDEFID)
        fprintf(stdout, "unknown\n");
      else
        {
          if (institutInqLongnamePtr(instID))
            fprintf(stdout, "%s\n", institutInqLongnamePtr(instID));
          else
            fprintf(stdout, "name unknown\n");
        }

      fprintf(stdout, " Source    : ");
      if (modelID == CDI_UNDEFID)
        fprintf(stdout, "unknown\n");
      else
        {
          if (modelInqNamePtr(modelID))
            {
              if (strncmp(modelInqNamePtr(modelID), "ECHAM5", 6) == 0) Source = S_ECHAM5;
              fprintf(stdout, "%s\n", modelInqNamePtr(modelID));
            }
          else
            fprintf(stdout, "name unknown\n");
        }
    }

  for (varID = 0; varID < nvars; varID++)
    {
      vlistInqVar(globs->ivlistID, varID, &gridID, &zaxisID, &timeID);
      code = vlistInqVarCode(globs->ivlistID, varID);
      if (code <= 0 || code >= MaxCodes)
        {
          Warning("Code number %d out of range, variable ignored!", code);
          continue;
        }
      gridtype = gridInqType(gridID);
      numlevel = zaxisInqSize(zaxisID);
      leveltype = zaxisInqType(zaxisID);

      vars[code].ivarID = varID;
      vars[code].igridID = gridID;
      vars[code].ogridID = gridID;
      vars[code].izaxisID = zaxisID;
      vars[code].ozaxisID = zaxisID;

      vars[code].detected = TRUE;

      if (globs->Debug)
        fprintf(stderr, "Code %3d  Levels = %3d  LevelType = %3d  GridType = %3d\n", code, numlevel, leveltype, gridtype);
    }

  if (globs->Debug) Message("FieldDim = %d", FieldDim);

  globs->Field = (double *) Malloc(FieldDim * sizeof(double));

  if (globs->Debug)
    for (code = 0; code < MaxCodes; code++)
      {
        if (vars[code].detected) fprintf(stderr, " Detected Code %3d with %3d level\n", code, zaxisInqSize(vars[code].izaxisID));
      }
}

/*
 * -----------------------------------------------------------
 * Define output variables
 * -----------------------------------------------------------
 */
static void
after_postcntl(struct Control *globs, struct Variable *vars)
{
  int code = 0;
  int gridID, zaxisID;
  int ovarID, ogridID, ozaxisID;
  int ovarID2;
  int ivarID, instID, modelID, tableID;
  char name[CDI_MAX_NAME], longname[CDI_MAX_NAME], units[CDI_MAX_NAME];
  char histstring[99];
  int datatype;

  snprintf(histstring, sizeof(histstring), "afterburner version %s  type = %d", VERSION, globs->Type);

#ifdef AFTERBURNER
  afterInqHistory(globs->istreamID);
  if (globs->Mean != 2) afterDefHistory(globs->ostreamID, histstring);
  if (globs->Mean >= 2) afterDefHistory(globs->ostreamID2, histstring);
#endif

  if (globs->Debug) lprintf(stdout);
  if (globs->Debug)
    for (code = 0; code < MaxCodes; code++)
      if (vars[code].detected)
        {
          gridID = vars[code].igridID;
          zaxisID = vars[code].izaxisID;
          zaxisName(zaxisInqType(zaxisID), zaxistypename);
          fprintf(stderr, " Detected Code %3d  grid %-8s size %5zu  level %2d %-8s\n", code, gridNamePtr(gridInqType(gridID)),
                  gridInqSize(gridID), zaxisInqSize(zaxisID), zaxistypename);
        }

  if (globs->Debug) lprintf(stdout);
  if (globs->Debug)
    for (code = 0; code < MaxCodes; code++)
      if (vars[code].needed)
        {
          fprintf(stderr, "   Needed Code %3d\n", code);
        }

  for (code = 0; code < MaxCodes; code++)
    if (vars[code].selected)
      {
        name[0] = 0;
        longname[0] = 0;
        units[0] = 0;
        ivarID = vars[code].ivarID;
        ogridID = vars[code].ogridID;
        ozaxisID = vars[code].ozaxisID;

        if (ogridID == -1)
          {
            /*
            Warning( "undefined grid for code %d", code);
            */
            continue;
          }
        if (ozaxisID == -1)
          {
            /*
            Warning( "undefined level for code %d", code);
            */
            continue;
          }

        instID = vlistInqVarInstitut(globs->ivlistID, ivarID);
        modelID = vlistInqVarModel(globs->ivlistID, ivarID);
        tableID = vlistInqVarTable(globs->ivlistID, ivarID);

        vars[code].missval = vlistInqVarMissval(globs->ivlistID, ivarID);
        vars[code].samp = NULL;

        if (DataType != -1)
          datatype = DataType;
        else
          datatype = vlistInqVarDatatype(globs->ivlistID, ivarID);

        if (vars[code].comp)
          {
            tableID = vars[code].tableID;
          }
        else
          {
            vlistInqVarName(globs->ivlistID, ivarID, name);
            vlistInqVarLongname(globs->ivlistID, ivarID, longname);
            vlistInqVarUnits(globs->ivlistID, ivarID, units);
          }

        if (globs->Mean != 2)
          {
            vlistDefTaxis(globs->ovlistID, globs->taxisID2);
            ovarID = vlistDefVar(globs->ovlistID, ogridID, ozaxisID, TIME_VARYING);
            if (globs->Mean) vlistDefVarTsteptype(globs->ovlistID, ovarID, TSTEP_AVG);
            vlistDefVarCode(globs->ovlistID, ovarID, code);
            vars[code].ovarID = ovarID;
            vlistDefVarInstitut(globs->ovlistID, ovarID, instID);
            vlistDefVarModel(globs->ovlistID, ovarID, modelID);
            vlistDefVarTable(globs->ovlistID, ovarID, tableID);
            if (name[0]) vlistDefVarName(globs->ovlistID, ovarID, name);
            if (longname[0]) vlistDefVarLongname(globs->ovlistID, ovarID, longname);
            if (units[0]) vlistDefVarUnits(globs->ovlistID, ovarID, units);
            vlistDefVarDatatype(globs->ovlistID, ovarID, datatype);
            vlistDefVarMissval(globs->ovlistID, ovarID, vars[code].missval);
          }

        if (globs->Mean >= 2)
          {
            vlistDefTaxis(globs->ovlistID2, globs->taxisID2);
            ovarID2 = vlistDefVar(globs->ovlistID2, ogridID, ozaxisID, TIME_VARYING);
            if (globs->Mean) vlistDefVarTsteptype(globs->ovlistID2, ovarID2, TSTEP_AVG);
            vlistDefVarCode(globs->ovlistID2, ovarID2, code);
            vars[code].ovarID2 = ovarID2;
            vlistDefVarInstitut(globs->ovlistID2, ovarID2, instID);
            vlistDefVarModel(globs->ovlistID2, ovarID2, modelID);
            vlistDefVarTable(globs->ovlistID2, ovarID2, tableID);
            if (name[0]) vlistDefVarName(globs->ovlistID2, ovarID2, name);
            if (longname[0]) vlistDefVarLongname(globs->ovlistID2, ovarID2, longname);
            if (units[0]) vlistDefVarUnits(globs->ovlistID2, ovarID2, units);
            vlistDefVarDatatype(globs->ovlistID2, ovarID2, datatype);
            vlistDefVarMissval(globs->ovlistID2, ovarID2, vars[code].missval);
          }
      }

  if (globs->Debug) lprintf(stdout);
  if (globs->Debug)
    for (code = 0; code < MaxCodes; code++)
      if (vars[code].selected)
        {
          gridID = vars[code].ogridID;
          zaxisID = vars[code].ozaxisID;
          zaxisName(zaxisInqType(zaxisID), zaxistypename);
          fprintf(stderr, " Selected Code %3d  grid %-8s size %5zu  level %2d %-8s\n", code, gridNamePtr(gridInqType(gridID)),
                  gridInqSize(gridID), zaxisInqSize(zaxisID), zaxistypename);
        }
}

static void
after_readVct(struct Control *globs, const char *vctfile)
{
  char line[1024];
  int i, n, nlines = 0;
  double va, vb;

  FILE *fp = fopen(vctfile, "r");
  if (fp == NULL) SysError("Open failed on %s", vctfile);

  while (fgets(line, 1023, fp))
    {
      if (line[0] == '#' || line[0] == '\0') continue;
      nlines++;
    }

  globs->nvct = nlines * 2;
  globs->vct = (double *) Malloc(globs->nvct * sizeof(double));

  rewind(fp);

  i = 0;
  while (fgets(line, 1023, fp))
    {
      if (line[0] == '#' || line[0] == '\0') continue;
      sscanf(line, "%d %lg %lg", &n, &va, &vb);
      globs->vct[i] = va;
      globs->vct[i + globs->nvct / 2] = vb;
      i++;
    }
  fprintf(stdout, "  Read VCT with %d hybrid levels from file %s\n", globs->nvct / 2 - 1, vctfile);

  fclose(fp);
}

#ifdef AFTERBURNER
static void
after_version(void)
{
#if defined(COMPILER)
  fprintf(stderr, "Compiler: %s\n", COMPILER);
#endif
#if defined(COMP_VERSION)
  fprintf(stderr, " version: %s\n", COMP_VERSION);
#endif
#if defined(HAVE_LIBSZ) || defined(_OPENMP)
  fprintf(stderr, "    with:");
#if defined(HAVE_LIBSZ)
  fprintf(stderr, " libsz");
#endif
#if defined(_OPENMP)
  fprintf(stderr, " OpenMP");
#endif
  fprintf(stderr, "\n");
#endif
#ifdef SYSTEM_TYPE
  fprintf(stderr, "System: %s\n", SYSTEM_TYPE);
#endif
  cdiPrintVersion();
  fprintf(stderr, "\n");
}
#endif

static void
after_control_init(struct Control *globs)
{
  memset(globs, 0, sizeof(struct Control));

  globs->AnalysisData = 0; /* 0 = ECHAM Data, 1 = ECMWF Spectral Analyses */
  globs->DayIn = 0;        /* day increment of infiles if Multi = TRUE    */
  globs->Debug = FALSE;
  globs->Extrapolate = TRUE;
  globs->Szip = FALSE;

  globs->istreamID = CDI_UNDEFID;
  globs->ostreamID = CDI_UNDEFID;
  globs->ostreamID2 = CDI_UNDEFID;
  globs->ivlistID = CDI_UNDEFID;
  globs->ovlistID = CDI_UNDEFID;
  globs->ovlistID2 = CDI_UNDEFID;
  globs->taxisID = -1;
  globs->taxisID2 = -1;
}

static void
after_variable_init(struct Variable *vars)
{
  memset(vars, 0, sizeof(struct Variable));

  vars->ivarID = -1;
  vars->ovarID = -1;
  vars->ovarID2 = -1;
  vars->izaxisID = -1;
  vars->ozaxisID = -1;
  vars->igridID = -1;
  vars->ogridID = -1;
  vars->tableID = -1;
}

#ifdef AFTERBURNER
static void
after_printCodes(void)
{
  int tableID = tableInq(-1, 128, "echam4");
  int codes[] = { 34, 35, 36, 131, 132, 135, 148, 149, 151, 156, 157, 259, 260, 261, 262, 263, 264, 268, 269, 270, 271, 275 };

  int ncodes = sizeof(codes) / sizeof(codes[0]);

  lprintf(stdout);

  fprintf(stdout, "  Code Name              Longname\n");
  fprintf(stdout, "  ---- ----              --------\n");

  for (int i = 0; i < ncodes; i++)
    {
      int code = codes[i];
      char name[CDI_MAX_NAME];
      name[0] = 0;
      char longname[CDI_MAX_NAME];
      longname[0] = 0;
      tableInqEntry(tableID, code, -1, name, longname, NULL);

      fprintf(stdout, " %4d", code);
      if (name[0])
        {
          fprintf(stdout, "  %-16s", name);
          if (longname[0]) fprintf(stdout, "  %s", longname);
        }
      else
        fprintf(stdout, "  var%d", code);

      fprintf(stdout, "\n");
    }

  lprintf(stdout);
}
#endif

/* =============================================== */
/* procstat   - appends info about memory usage    */
/*              and time consumption               */
/* =============================================== */
#ifdef AFTERBURNER
static void
after_procstat(char *procpath, int truncation)
{
  time_t tp;
  char mtype[12];
  char stat_file[128];

  double CPUTime = ((double) clock() - starttime) / CLOCKS_PER_SEC;

  (void) time(&tp);
  long yy = gmtime(&tp)->tm_year + 1900;
  long mm = gmtime(&tp)->tm_mon + 1;
  long dd = gmtime(&tp)->tm_mday;
  long hh = gmtime(&tp)->tm_hour;
  long mi = gmtime(&tp)->tm_min;
  char *name = getpwuid(getuid())->pw_name;

  char *proc = strrchr(procpath, '/');
  if (proc == 0)
    proc = procpath;
  else
    proc++;

  strcpy(stat_file, "/pf/m/m214003/local/log/after.log");

  double MaxMBytes = (double) memTotal() / 1048576.;

  FILE *sf = fopen(stat_file, "a");
  if (sf)
    {
      char unknown[] = "";
      char *hostname;

      if ((hostname = getenv("HOST")) == NULL) hostname = unknown;

      setvbuf(sf, (char *) NULL, _IONBF, 0);
      fprintf(sf, "%.7s %4.4ld.%2.2ld.%2.2ld %2.2ld:%2.2ld %s %-9.9s %7.1f %7.1f T%3.3d %s\n", name, yy, mm, dd, hh, mi, VERSION,
              proc, MaxMBytes, CPUTime, truncation, hostname);

      fclose(sf);
    }

#if defined(CRAY)
#if defined(_CRAYMPP)
  strcpy(mtype, " CRAYMPP --");
#elif (_MAXVL == 64)
  strcpy(mtype, " CRAYVL64 -");
#elif (_MAXVL == 128)
  strcpy(mtype, " CRAYVL128 ");
#else
  strcpy(mtype, " CRAY -----");
#endif
#elif defined(SX)
  strcpy(mtype, " NECSX ----");
#elif defined(__uxp__)
  strcpy(mtype, " FUJI -----");
#elif defined(sun)
  strcpy(mtype, " SUN ------");
#elif defined(i386)
  strcpy(mtype, " i386 -----");
#elif defined(sgi)
  strcpy(mtype, " sgi ------");
#else
  strcpy(mtype, "-----------");
#endif

  fprintf(stdout, "   NORMAL EXIT\n");
  fprintf(stdout, " ------   End    after  -%-11.11s- %7.1f sec", mtype, CPUTime);
  if (MaxMBytes > 0)
    fprintf(stdout, " --- %7.1f MB ---\n", MaxMBytes);
  else
    fprintf(stdout, " ----------------\n");
}
#endif

static void
after_processing(struct Control *globs, struct Variable *vars)
{
  int i;

  //#if defined(PSTREAM_H)
  //  globs->istreamID = streamOpenRead(cdoStreamName(0));
  //#else
  globs->istreamID = streamOpenRead(ifile);
  if (globs->istreamID < 0) cdiError(globs->istreamID, "Open failed on %s", ifile);
  //#endif
  if (ofiletype == -1) ofiletype = streamInqFiletype(globs->istreamID);

  globs->ivlistID = streamInqVlist(globs->istreamID);
  globs->taxisID = vlistInqTaxis(globs->ivlistID);
  globs->taxisID2 = taxisDuplicate(globs->taxisID);

  if (globs->Mean != 2)
    {
#ifdef CDO
      globs->ostreamID = cdoStreamOpenWrite(cdoStreamName(ofileidx), ofiletype);
#else
      globs->ostreamID = streamOpenWrite(ofile, ofiletype);
      if (globs->ostreamID < 0) cdiError(globs->ostreamID, "Open failed on %s", ofile);
#endif

      if (globs->Szip) streamDefCompType(globs->ostreamID, CDI_COMPRESS_SZIP);

      globs->ovlistID = vlistCreate();
    }
#ifdef AFTERBURNER
  if (globs->Mean >= 2)
    {
      globs->ostreamID2 = streamOpenWrite(ofile2, ofiletype);
      if (globs->ostreamID2 < 0) cdiError(globs->ostreamID2, "Open failed on %s", ofile2);

      if (globs->Szip) streamDefCompType(globs->ostreamID, CDI_COMPRESS_SZIP);

      globs->ovlistID2 = vlistCreate();
    }
#endif

  /* ---------------- */
  /*  pre-processing  */
  /* ---------------- */
  after_precntl(globs, vars);

  /* ----------------- */
  /*  initializations  */
  /* ----------------- */

  after_setCodes(globs, vars, MaxCodes, globs->NumCodesRequest);

  if (globs->unitsel == 2)
    for (i = 0; i < globs->NumLevelRequest; i++) globs->LevelRequest[i] = globs->LevelRequest[i] * 1000;

  if (!globs->AnalysisData)
    for (i = 0; i < globs->NumLevelRequest; i++)
      {
        if ((globs->LevelRequest[i] >= 65535) && globs->unitsel && ofiletype == CDI_FILETYPE_GRB)
          {
            fprintf(stderr, "\n Level %9.2f out of range (max=65535)!\n", globs->LevelRequest[i]);
            exit(1);
          }

        if (!globs->unitsel && globs->Type >= 20 && globs->NumLevelRequest > 1 && IS_EQUAL(globs->LevelRequest[i], 0))
          {
            fprintf(stderr, "\n Level %9.2f illegal for Type %d\n", globs->LevelRequest[i], globs->Type);
            exit(1);
          }
      }

  after_setLevel(globs);

  after_dimcalc(globs);

  globs->rcoslat = (double *) Malloc(globs->Latitudes * sizeof(double));
  globs->coslat = (double *) Malloc(globs->Latitudes * sizeof(double));
  globs->DerivationFactor = (double *) Malloc(globs->Latitudes * sizeof(double));

  if (globs->Type < 50 && globs->AnalysisData)
    {
      fprintf(stderr, " ::::::::::::::::::::::::::::::::::::::::::::::\n");
      fprintf(stderr, " -> Type < 50 is not appropriate for Analysis.\n");
      fprintf(stderr, " -> Please check wether you can use Type >= 50.\n");
      fprintf(stderr, " -> Premature Exit. Sorry.\n");
      exit(1);
    }

  if (globs->Type == 10 || globs->Type == 40 || globs->Type == 60)
    {
      if (ofiletype == CDI_FILETYPE_GRB)
        Error("Can't write fourier coefficients to GRIB!");
      else if (ofiletype == CDI_FILETYPE_NC || ofiletype == CDI_FILETYPE_NC2 || ofiletype == CDI_FILETYPE_NC4
               || ofiletype == CDI_FILETYPE_NC4C || ofiletype == CDI_FILETYPE_NC5)
        Error("Can't write fourier coefficients to NetCDF!");
    }

  filename = strrchr(ifile, '/');
  if (filename == 0)
    filename = ifile;
  else
    filename++;

  if (globs->Type >= 30 && globs->Type < 50 && (vars[DIVERGENCE].selected || vars[VELOPOT].selected || vars[VORTICITY].selected
                                                || vars[STREAM].selected || globs->AnalysisData))
    {
      /*
      int newtype = 0;
      */
      if (globs->Type == 30) globs->Type = 70;
      if (globs->Type == 40) globs->Type = 60;
      if (globs->Type == 41) globs->Type = 61;

      if (globs->AnalysisData)
        fprintf(stderr, "\n TYPE changed to %d (for analysis data)\n", globs->Type);
      else
        fprintf(stderr, "\n TYPE changed to %d (with code %d, %d, %d or %d)\n", globs->Type, DIVERGENCE, VELOPOT, VORTICITY,
                STREAM);
      /*
      if ( globs->Type == 30 ) newtype = 70;
      if ( globs->Type == 40 ) newtype = 60;
      if ( globs->Type == 41 ) newtype = 61;

      if ( globs->AnalysisData )
        fprintf(stderr,"\n Attention: TYPE isn't changed to %d anymore (for
      analysis data)!!!\n", globs->Type); else fprintf(stderr,"\n Attention:
      TYPE isn't changed to %d anymore (with code %d, %d, %d or %d)!!!\n",
                newtype, DIVERGENCE, VELOPOT, VORTICITY, STREAM);
      */
    }

  if (globs->AnalysisData)
    after_AnalysisDependencies(vars, MaxCodes);
  else
    {
      after_EchamDependencies(vars, MaxCodes, globs->Type, Source);
      vars[GEOPOTENTIAL].needed |= globs->Type >= 30 || vars[SLP].comp || vars[GEOPOTHEIGHT].comp;
    }

  /*  if ( vars[U_WIND].needed || vars[V_WIND].needed ) */
  if (vars[U_WIND].comp || vars[V_WIND].comp)
    {
      globs->dv2uv_f1 = (double *) Malloc(globs->DimSP_half * sizeof(double));
      globs->dv2uv_f2 = (double *) Malloc(globs->DimSP_half * sizeof(double));
      geninx(globs->Truncation, globs->dv2uv_f1, globs->dv2uv_f2);
    }

  /* --------- */
  /*  Control  */
  /* --------- */

  after_defineLevel(globs, vars);

  after_defineGrid(globs, vars);

  after_postcntl(globs, vars); /* define output variables */

  after_control(globs, vars);

#ifdef CDO
  if (globs->ostreamID != CDI_UNDEFID) pstreamClose(globs->ostreamID);
#else
  if (globs->ostreamID2 != CDI_UNDEFID) streamClose(globs->ostreamID2);
  if (globs->ostreamID != CDI_UNDEFID) streamClose(globs->ostreamID);
#endif
#ifdef CDO
  processDefVarNum(vlistNvars(globs->ivlistID));
  processSelf().addNvals(streamNvals(globs->istreamID));
#endif
  streamClose(globs->istreamID);

  if (globs->rcoslat) Free(globs->rcoslat);
  if (globs->coslat) Free(globs->coslat);
  if (globs->DerivationFactor) Free(globs->DerivationFactor);

  if (globs->Field) Free(globs->Field);

  if (globs->poli) Free(globs->poli);
  if (globs->pold) Free(globs->pold);
  if (globs->pdev) Free(globs->pdev);
  if (globs->pol2) Free(globs->pol2);
  if (globs->pol3) Free(globs->pol3);
}

extern char *optarg;
extern int optind, opterr, optopt;

#ifdef AFTERBURNER
static int
afterburner(int argc, char *argv[])
{
  int i, code;
  char *proc = argv[0];
  char Line[132];
  int c;
  int fargc0, fargcn;
  FILE *fp;
  int numThreads = 0;
  char *Vctfile = NULL;
  extern int dmemory_ExitOnError;

  dmemory_ExitOnError = 1;

  starttime = (double) clock();

#ifdef AFTERBURNER
  { /* check character device on stdin and stdout */
    struct stat statbuf;
    fstat(0, &statbuf);
    if (S_ISCHR(statbuf.st_mode)) stdin_is_tty = 1;
    fstat(1, &statbuf);
    if (S_ISCHR(statbuf.st_mode)) stdout_is_tty = 1;
  }
#endif

  /* ------------------- */
  /*  print information  */
  /* ------------------- */

  lprintf(stdout);
  fprintf(stdout, "  afterburner version %s\n", VERSION);
  fprintf(stdout, "  ECHAM & analyses postprocessor\n");

  if (sizeof(double) != 8 || sizeof(int) < 4)
    {
      fprintf(stderr, "byte size of type double %d\n", (int) sizeof(double));
      fprintf(stderr, "byte size of type int %d\n", (int) sizeof(int));
      fprintf(stderr, "byte size of type size_t %d\n", (int) sizeof(size_t));
      return 1;
    }

  fp = fopen("/pf/m/m214003/doc/afterburner.doc", "r");
  if (fp)
    {
      do
        {
          fgets(Line, 130, fp);
          fprintf(stdout, "%s", &Line[1]);
        }
      while (!feof(fp) && Line[0] == '#');
      fclose(fp);
    }

  struct Control *globs = (struct Control *) Malloc(sizeof(struct Control));
  after_control_init(globs);

  globs->Verbose = 1;

  /* --------------------- */
  /*  options & filenames  */
  /* --------------------- */
  extern int labort_after;

  while ((c = getopt(argc, argv, "P:b:v:acdgpVw")) != EOF) switch (c)
      {
      case 'a': globs->AnalysisData = 1; break;
      case 'b': Message("option -b not longer needed!"); break;
      case 'c': after_printCodes(); break;
      case 'd': globs->Debug = 1; break;
      case 'p': lparallelread = true; break;
      case 'P': numThreads = atoi(optarg); break;
      case 'V': after_version(); break;
      case 'v': Vctfile = optarg; break;
      case 'w': labort_after = FALSE; break;
      default: Message("option -%c unsupported!", optopt); after_usage();
      }

#if defined(_OPENMP)
  lprintf(stdout);
  if (numThreads <= 0) numThreads = 1;
  omp_set_num_threads(numThreads);
  if (omp_get_max_threads() > omp_get_num_procs())
    fprintf(stdout, " Number of threads is greater than number of Cores=%d!\n", omp_get_num_procs());
  fprintf(stdout, " OpenMP:  num_procs=%d  max_threads=%d\n", omp_get_num_procs(), omp_get_max_threads());
#else
  if (numThreads > 0)
    {
      fprintf(stderr, "Option -P failed, OpenMP support not compiled in!\n");
      return -1;
    }
#endif

  if (lparallelread) fprintf(stdout, " Parallel read enabled\n");

  fargc0 = optind;
  fargcn = argc;

  if (optind < argc) ifile = argv[optind++];
  if (!ifile)
    {
      Message("*** Missing input file ***");
      after_usage();
    }

  struct Variable vars[MaxCodes + 5];
  for (code = 0; code < MaxCodes + 5; code++) after_variable_init(&vars[code]);

  after_parini(globs, vars); /* read namelist parameter */

  fprintf(stdout, "   Input File: %-25s\n", ifile);
  if (globs->Mean >= 2)
    {
      if (fargcn - fargc0 >= 3) ofile2 = argv[--fargcn];
      if (!ofile2)
        {
          Message("*** Missing variance file ***");
          after_usage();
        }
    }

  if (globs->Mean != 2)
    {
      if (optind < argc) ofile = argv[optind++];
      if (fargcn - fargc0 >= 2) ofile = argv[--fargcn];
      if (!ofile)
        {
          Message("*** Missing output file ***");
          after_usage();
        }
      fprintf(stdout, "  Output File: %-25s\n", ofile);
    }

  globs->Nfiles = fargcn - fargc0 - 1;
  if (globs->Nfiles > 0)
    {
      if (globs->Multi > 0) Error("Namelist parameter MULTI works only with one inputfile");

      ifiles = (char **) Malloc(globs->Nfiles * sizeof(char *));
      for (i = 0; i < globs->Nfiles; i++) ifiles[i] = argv[--fargcn];
    }

  if (ofile2) fprintf(stdout, "Variance File: %-25s\n", ofile2);

  if (globs->Debug)
    {
      extern int afterDebug;
      afterDebug = globs->Debug;
      fprintf(stderr, "* Debug on!                              *\n");
      fprintf(stderr, "  Maximum ffts to run in parallel:  %ld\n", get_nfft());
    }

  /* read optional VCT */
  if (Vctfile) after_readVct(globs, Vctfile);

  /* --------------------- */
  /*  open in/output file  */
  /* --------------------- */

  cdiDefGlobal("REGULARGRID", 1);

  after_processing(globs, vars);

  after_procstat(proc, globs->Truncation);

  FreeMean(vars);

  Free(globs);

  return 0;
}
#endif

#ifdef CDO
void *
Afterburner(void *process)
{
  cdoInitialize(process);

  CDO_task = true;

  lstdout = !Options::silentMode;

  struct Control *globs = (struct Control *) Malloc(sizeof(struct Control));
  after_control_init(globs);

  globs->Verbose = cdoVerbose;

  if (operatorArgc() == 1)
    {
      const char *vctfile = operatorArgv()[0];
      after_readVct(globs, vctfile);
    }

  struct Variable vars[MaxCodes + 5];
  for (int code = 0; code < MaxCodes + 5; code++) after_variable_init(&vars[code]);

  after_parini(globs, vars); /* read namelist parameter */

  if (cdoDefaultFileType != CDI_UNDEFID) ofiletype = cdoDefaultFileType;

  int streamCnt = cdoStreamCnt();
  int nfiles = streamCnt - 1;

  ofileidx = nfiles;

  ifile = strdup(cdoGetStreamName(0).c_str());
  ofile = (char *) cdoGetStreamName(nfiles).c_str();

  globs->Nfiles = nfiles - 1;
  if (globs->Nfiles > 0)
    {
      if (globs->Multi > 0) Error("Namelist parameter MULTI works only with one inputfile");

      ifiles = (char **) Malloc(globs->Nfiles * sizeof(char *));
      for (int i = 0; i < globs->Nfiles; ++i) ifiles[i] = (char *) cdoGetStreamName(--nfiles).c_str();
      for (int i = 0; i < globs->Nfiles; ++i) printf("files %d %s\n", i + 1, ifiles[i]);
    }

  after_processing(globs, vars);

  FreeMean(vars);

  Free(globs);
  Free(ifile);

  cdoFinish();

  return 0;
}
#else
int
main(int argc, char *argv[])
{
  return afterburner(argc, argv);
}
#endif

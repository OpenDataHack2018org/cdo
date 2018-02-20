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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <thread>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>
#include <sys/stat.h> /* stat */

#include <map>

#include <cdi.h>
#include "pipe.h"
#include "error.h"
#include "cdoDebugOutput.h"
#include "pstream.h"
#include "timer.h"
#include "exception.h"
#include "cdoOptions.h"
#include "util.h"
#include "commandline.h"
#include "assert.h"
#include "dmemory.h"
#include "compare.h"
#include "cdo_vlist.h"
#include "functs.h"
#include "array.h"
#include "cdoOptions.h"
#include "cdo_history.h"

#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
#include "pthread_debug.h"

// TODO: make threadsafe
static int pthreadScope = 0;
static int processNum = 0;
void
setProcessNum(int p_num)
{
  processNum = p_num;
}
static pthread_mutex_t streamOpenReadMutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t streamOpenWriteMutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t streamMutex = PTHREAD_MUTEX_INITIALIZER;

static std::mutex _pstream_map_mutex;
#define PSTREAM_LOCK() _pstream_map_mutex.lock();
#define PSTREAM_UNLOCK() _pstream_map_mutex.unlock();

#else

#define PSTREAM_LOCK()
#define PSTREAM_UNLOCK()
//-----------
#endif

static std::map<int, PstreamType> _pstream_map;
static int next_pstream_id = 1;
static int createdPstreams = 0;

static void set_comp(int a, int b); /*TEMP*/

PstreamType *
create_pstream()
{
  PSTREAM_LOCK();
  auto new_entry = _pstream_map.insert(std::make_pair(next_pstream_id, PstreamType(next_pstream_id)));
  if (new_entry.second == false)
    {
      ERROR("A Pstream could not be created, ID: ", next_pstream_id);
    }
  if (CdoDebug::PSTREAM)
    {
      MESSAGE("Created new pstream with ID:", next_pstream_id);
    }
  next_pstream_id++;
  createdPstreams++;
  new_entry.first->second.ispipe = true;
  PSTREAM_UNLOCK();

  return &new_entry.first->second;
}
// temporary function: will be replaced by according
// PstreamType::PstreamType(..)
PstreamType *
create_pstream(std::vector<std::string> p_filenameList)
{
  PstreamType *new_entry = create_pstream();
  new_entry->m_mfnames = p_filenameList;
  new_entry->m_name = p_filenameList[0];
  new_entry->ispipe = false;
  return new_entry;
}
// temporary function: will be replaced by according
// PstreamType::PstreamType(..)
PstreamType *
create_pstream(std::string p_filename)
{
  return create_pstream(std::vector<std::string>{ p_filename });
}

// temporary function: will be replaced by according
// PstreamType::PstreamType(..)
PstreamType *
create_pstream(int processID, int pstreamIDX)
{
  PstreamType *new_pstream = create_pstream();
  new_pstream->pipe = std::make_shared<pipe_t>();
  new_pstream->pipe->pipeSetName(processID, pstreamIDX);
  new_pstream->ispipe = true;

  return new_pstream;
}

PstreamType *
pstreamToPointer(int idx)
{
  PSTREAM_LOCK();
  auto pstream_iterator = _pstream_map.find(idx);
  PSTREAM_UNLOCK();
  if (pstream_iterator == _pstream_map.end())
    {
      Error("pstream index %d undefined!", idx);
    }

  return &pstream_iterator->second;
}

void
PstreamType::init()
{
  isopen = true;
  ispipe = false;
  m_fileID = -1;
  m_vlistID = -1;
  tsID = -1;
  m_filetype = -1;
  m_name = "";
  tsID0 = 0;
  mfiles = 0;
  nfiles = 0;
  m_varID = -1;
  m_varlist = NULL;
#ifdef HAVE_LIBPTHREAD
  pipe = NULL;
//  pstreamptr->rthreadID  = 0;
//  pstreamptr->wthreadID  = 0;
#endif
}
PstreamType::PstreamType(int p_id) : self(p_id) { init(); }

bool
PstreamType::isPipe()
{
  return ispipe;
}

int
PstreamType::pstreamOpenReadPipe()
{
#ifdef HAVE_LIBPTHREAD
  rthreadID = pthread_self();

  if (CdoDebug::PSTREAM) MESSAGE("pipe ", pipe->name);

  return self;
#else
  cdoAbort("Cannot use pipes, pthread support not compiled in!");
#endif
}

void
PstreamType::pstreamOpenReadFile(const char *p_args)
{
  std::string filename;

  if (mfiles)
    {
      filename = m_mfnames[0];
      nfiles = 1;
    }
  else
    {
      filename = std::string(p_args);
    }

  if (CdoDebug::PSTREAM) MESSAGE("Opening (r) file: ", filename.c_str());

#ifdef HAVE_LIBPTHREAD
  if (Threading::cdoLockIO)
    pthread_mutex_lock(&streamMutex);
  else
    pthread_mutex_lock(&streamOpenReadMutex);
#endif
  int fileID = streamOpenRead(filename.c_str());
  if (fileID < 0)
    {
      isopen = false;
      cdiOpenError(fileID, "Open failed on >%s<", filename.c_str());
    }

  if (cdoDefaultFileType == CDI_UNDEFID) cdoDefaultFileType = streamInqFiletype(fileID);

  cdoInqHistory(fileID);
#ifdef HAVE_LIBPTHREAD
  if (Threading::cdoLockIO)
    pthread_mutex_unlock(&streamMutex);
  else
    pthread_mutex_unlock(&streamOpenReadMutex);
#endif

  mode = 'r';
  m_name = filename;
  m_fileID = fileID;
}

int
PstreamType::pstreamOpenWritePipe(const char *pipename, int filetype)
{
#ifdef HAVE_LIBPTHREAD
  if (CdoDebug::PSTREAM)
    {
      MESSAGE("pipe ", pipename);
    }
  wthreadID = pthread_self();
  m_filetype = filetype;
  return self;
#else
  return -1;
#endif
}

int
PstreamType::pstreamOpenWriteFile(int filetype)
{
  if (CdoDebug::PSTREAM)
    {
      MESSAGE("Opening (w) file ", m_name);
    }

  if (filetype == CDI_UNDEFID) filetype = CDI_FILETYPE_GRB;

  if (processNum == 1 && Threading::ompNumThreads == 1) timer_start(timer_write);

#ifdef HAVE_LIBPTHREAD
  if (Threading::cdoLockIO)
    pthread_mutex_lock(&streamMutex);
  else
    pthread_mutex_lock(&streamOpenWriteMutex);
#endif

  int fileID = streamOpenWrite(m_name.c_str(), filetype);

#ifdef HAVE_LIBPTHREAD
  if (Threading::cdoLockIO)
    pthread_mutex_unlock(&streamMutex);
  else
    pthread_mutex_unlock(&streamOpenWriteMutex);
#endif

  if (processNum == 1 && Threading::ompNumThreads == 1) timer_stop(timer_write);
  if (fileID < 0) cdiOpenError(fileID, "Open failed on >%s<", m_name.c_str());

  cdoDefHistory(fileID, commandLine());

  if (cdoDefaultByteorder != CDI_UNDEFID) streamDefByteorder(fileID, cdoDefaultByteorder);

  set_comp(fileID, filetype);

  mode = 'w';
  m_fileID = fileID;
  m_filetype = filetype;

  return self;
}

void
PstreamType::openAppend(const char *p_filename)
{
  if (processNum == 1 && Threading::ompNumThreads == 1)
    {
      timer_start(timer_write);
    }
#ifdef HAVE_LIBPTHREAD
  if (Threading::cdoLockIO)
    {
      pthread_mutex_lock(&streamMutex);
    }
  else
    {
      pthread_mutex_lock(&streamOpenReadMutex);
    }
#endif

  int fileID = streamOpenAppend(p_filename);

#ifdef HAVE_LIBPTHREAD
  if (Threading::cdoLockIO)
    {
      pthread_mutex_unlock(&streamMutex);
    }
  else
    {
      pthread_mutex_unlock(&streamOpenReadMutex);
    }
#endif
  if (processNum == 1 && Threading::ompNumThreads == 1)
    {
      timer_stop(timer_write);
    }
  if (fileID < 0)
    {
      cdiOpenError(fileID, "Open failed on >%s<", p_filename);
    }

  int filetype = streamInqFiletype(fileID);
  set_comp(fileID, filetype);

  m_name = p_filename;
  mode = 'a';
  m_fileID = fileID;
}
void
PstreamType::closePipe()
{
  pipe->close();
  pthread_cond_signal(pipe->recInq);

  pthread_mutex_lock(pipe->m_mutex);
  isopen = false;
  pthread_mutex_unlock(pipe->m_mutex);
  pthread_cond_signal(pipe->isclosed);

  pthread_join(wthreadID, NULL);
}

size_t
PstreamType::getNvals()
{
  if (ispipe)
    {
      return pipe->nvals;
    }
  else
    {
      return streamNvals(m_fileID);
    }
}

void
PstreamType::waitForPipe()
{
  pipe->close();
  std::unique_lock<std::mutex> locked_mutex(pipe->m_mutex);
  while (isopen)
    {
      if (CdoDebug::PSTREAM)
        {
          MESSAGE("wait of read close");
        }
      pthread_cond_wait(pipe->isclosed, locked_mutex);
    }
  locked_mutex.unlock();
}

void
PstreamType::close()
{
  if (ispipe)
    {
#ifdef HAVE_LIBPTHREAD
      pthread_t threadID = pthread_self();

      if (CdoDebug::PSTREAM)
        {
          MESSAGE("thID: ", threadID, " rthID: ", rthreadID, " wthID: ", wthreadID);
        }
      if (pthread_equal(threadID, rthreadID))
        {
          closePipe();
        }
      else if (pthread_equal(threadID, wthreadID))
        {
          waitForPipe();
        }
      else
        {
          Error("Internal problem! Close pipe ", m_name.c_str());
        }
#else
      cdoAbort("Cannot use pipes, pthread support not compiled in!");
#endif
    }
  else
    {
      if (CdoDebug::PSTREAM)
        {
          MESSAGE(m_name.c_str(), " fileID ", m_fileID);
        }
#ifdef HAVE_LIBPTHREAD
      if (Threading::cdoLockIO) pthread_mutex_lock(&streamMutex);
#endif
      streamClose(m_fileID);
#ifdef HAVE_LIBPTHREAD
      if (Threading::cdoLockIO) pthread_mutex_unlock(&streamMutex);
#endif

      if (cdoExpMode == CDO_EXP_REMOTE)
        {
          if (mode == 'w')
            {
              extern const char *cdojobfiles;
              FILE *fp = fopen(cdojobfiles, "a");
              fprintf(fp, "%s\n", m_name.c_str());
              fclose(fp);
            }
        }

      if (m_varlist)
        {
          Free(m_varlist);
          m_varlist = NULL;
        }
    }
}

int
PstreamType::inqVlist()
{
  int vlistID = -1;

#ifdef HAVE_LIBPTHREAD
  // read from pipe
  if (ispipe)
    {
      vlistID = pipe->pipeInqVlist(m_vlistID);
      if (vlistID == -1) cdoAbort("Couldn't read data from input stream %s!", m_name.c_str());
    }
  // read from file through cdi streamInqVlist
  else
#endif
    {
      if (processNum == 1 && Threading::ompNumThreads == 1) timer_start(timer_read);
#ifdef HAVE_LIBPTHREAD
      if (Threading::cdoLockIO) pthread_mutex_lock(&streamMutex);
#endif
      vlistID = streamInqVlist(m_fileID);
      if (vlistID == -1)
        {
          cdoAbort("Couldn't read data from input fileID %d!", m_fileID);
        }

#ifdef HAVE_LIBPTHREAD
      if (Threading::cdoLockIO) pthread_mutex_unlock(&streamMutex);
#endif
      if (processNum == 1 && Threading::ompNumThreads == 1) timer_stop(timer_read);

      int nsubtypes = vlistNsubtypes(vlistID);
      if (nsubtypes > 1)
        cdoWarning("Subtypes are unsupported, the processing results are "
                   "possibly wrong!");

      if (cdoDefaultTimeType != CDI_UNDEFID) taxisDefType(vlistInqTaxis(vlistID), cdoDefaultTimeType);

      m_vlistID = vlistID;
    }

  return vlistID;
}

void
PstreamType::defVarList(int p_vlistID)
{
  int filetype = m_filetype;

  if (m_vlistID != -1) cdoAbort("Internal problem, vlist already defined!");

  if (m_varlist != NULL) cdoAbort("Internal problem, varlist already allocated!");

  int nvars = vlistNvars(p_vlistID);
  assert(nvars > 0);

  varlist_t *varlist = (varlist_t *) Malloc(nvars * sizeof(varlist_t));

  for (int varID = 0; varID < nvars; ++varID)
    {
      varlist[varID].gridsize = gridInqSize(vlistInqVarGrid(p_vlistID, varID));
      varlist[varID].datatype = vlistInqVarDatatype(p_vlistID, varID);
      varlist[varID].missval = vlistInqVarMissval(p_vlistID, varID);
      varlist[varID].addoffset = vlistInqVarAddoffset(p_vlistID, varID);
      varlist[varID].scalefactor = vlistInqVarScalefactor(p_vlistID, varID);

      varlist[varID].check_datarange = false;

      int laddoffset = IS_NOT_EQUAL(varlist[varID].addoffset, 0);
      int lscalefactor = IS_NOT_EQUAL(varlist[varID].scalefactor, 1);

      int datatype = varlist[varID].datatype;

      if (filetype == CDI_FILETYPE_NC || filetype == CDI_FILETYPE_NC2 || filetype == CDI_FILETYPE_NC4
          || filetype == CDI_FILETYPE_NC4C || filetype == CDI_FILETYPE_NC5)
        {
          if (datatype == CDI_DATATYPE_UINT8
              && (filetype == CDI_FILETYPE_NC || filetype == CDI_FILETYPE_NC2 || filetype == CDI_FILETYPE_NC5))
            {
              datatype = CDI_DATATYPE_INT16;
              varlist[varID].datatype = datatype;
            }

          if (datatype == CDI_DATATYPE_UINT16
              && (filetype == CDI_FILETYPE_NC || filetype == CDI_FILETYPE_NC2 || filetype == CDI_FILETYPE_NC5))
            {
              datatype = CDI_DATATYPE_INT32;
              varlist[varID].datatype = datatype;
            }

          if (laddoffset || lscalefactor)
            {
              if (datatype == CDI_DATATYPE_INT8 || datatype == CDI_DATATYPE_UINT8 || datatype == CDI_DATATYPE_INT16
                  || datatype == CDI_DATATYPE_UINT16)
                varlist[varID].check_datarange = true;
            }
          else if (cdoCheckDatarange)
            {
              varlist[varID].check_datarange = true;
            }
        }
    }

  m_varlist = varlist;
  m_vlistID = p_vlistID; /* used for -r/-a */
}

void
PstreamType::defVlist(int p_vlistID)
{
#ifdef HAVE_LIBPTHREAD
  if (ispipe)
    {
      if (CdoDebug::PSTREAM) MESSAGE(m_name.c_str(), " pstreamID ", self);
      int vlistIDcp = vlistDuplicate(p_vlistID);
      pipe->pipeDefVlist(m_vlistID, vlistIDcp);
    }
  else
#endif
    {
      if (cdoDefaultDataType != CDI_UNDEFID)
        {
          int varID, nvars = vlistNvars(p_vlistID);

          for (varID = 0; varID < nvars; ++varID)
            vlistDefVarDatatype(p_vlistID, varID, cdoDefaultDataType);

          if (cdoDefaultDataType == CDI_DATATYPE_FLT64 || cdoDefaultDataType == CDI_DATATYPE_FLT32)
            {
              for (varID = 0; varID < nvars; varID++)
                {
                  vlistDefVarAddoffset(p_vlistID, varID, 0.0);
                  vlistDefVarScalefactor(p_vlistID, varID, 1.0);
                }
            }
        }

      if (cdoChunkType != CDI_UNDEFID)
        {
          int varID, nvars = vlistNvars(p_vlistID);

          for (varID = 0; varID < nvars; ++varID)
            vlistDefVarChunkType(p_vlistID, varID, cdoChunkType);
        }

      if (CDO_CMOR_Mode)
        {
          cdo_def_tracking_id(p_vlistID, "tracking_id");
          cdo_def_creation_date(p_vlistID);
        }

      if (CDO_Version_Info) cdiDefAttTxt(p_vlistID, CDI_GLOBAL, "CDO", (int) strlen(cdoComment()), cdoComment());

#ifdef _OPENMP
      if (Threading::ompNumThreads > 1)
        cdiDefAttInt(p_vlistID, CDI_GLOBAL, "cdo_openmp_thread_number", CDI_DATATYPE_INT32, 1,
                     &Threading::ompNumThreads);
#endif
      defVarList(p_vlistID);

      if (processNum == 1 && Threading::ompNumThreads == 1) timer_start(timer_write);
#ifdef HAVE_LIBPTHREAD
      if (Threading::cdoLockIO) pthread_mutex_lock(&streamMutex);
#endif
      streamDefVlist(m_fileID, p_vlistID);
#ifdef HAVE_LIBPTHREAD
      if (Threading::cdoLockIO) pthread_mutex_unlock(&streamMutex);
#endif
      if (processNum == 1 && Threading::ompNumThreads == 1) timer_stop(timer_write);
    }
}

void
PstreamType::inqRecord(int *varID, int *levelID)
{
#ifdef HAVE_LIBPTHREAD
  if (ispipe)
    {
      if (CdoDebug::PSTREAM)
        {
          MESSAGE(pipe->name.c_str(), " pstreamID ", self);
        }
      pipe->pipeInqRecord(varID, levelID);
    }

  else
#endif
    {
      if (processNum == 1 && Threading::ompNumThreads == 1) timer_start(timer_read);
#ifdef HAVE_LIBPTHREAD
      if (Threading::cdoLockIO) pthread_mutex_lock(&streamMutex);
#endif
      streamInqRecord(m_fileID, varID, levelID);
#ifdef HAVE_LIBPTHREAD
      if (Threading::cdoLockIO) pthread_mutex_unlock(&streamMutex);
#endif
      if (processNum == 1 && Threading::ompNumThreads == 1) timer_stop(timer_read);
    }
}

void
PstreamType::defRecord(int varID, int levelID)
{
  m_varID = varID;

#ifdef HAVE_LIBPTHREAD
  if (ispipe)
    {

      if (CdoDebug::PSTREAM)
        {
          MESSAGE(m_name.c_str(), " pstreamid ", self);
        }
      pipe->pipeDefRecord(varID, levelID);
    }
  else
#endif
    {
      if (processNum == 1 && Threading::ompNumThreads == 1) timer_start(timer_write);
#ifdef HAVE_LIBPTHREAD
      if (Threading::cdoLockIO) pthread_mutex_lock(&streamMutex);
#endif
      streamDefRecord(m_fileID, varID, levelID);
#ifdef HAVE_LIBPTHREAD
      if (Threading::cdoLockIO) pthread_mutex_unlock(&streamMutex);
#endif
      if (processNum == 1 && Threading::ompNumThreads == 1) timer_stop(timer_write);
    }
}

void
PstreamType::readRecordF(float *data, size_t *nmiss)
{
#ifdef HAVE_LIBPTHREAD

  if (ispipe)
    {
      cdoAbort("pipeReadRecord not implemented for memtype float!");
      // pipeReadRecord(pstreamptr, data, nmiss);
    }
  else
#endif
    {
      if (processNum == 1 && Threading::ompNumThreads == 1) timer_start(timer_read);
#ifdef HAVE_LIBPTHREAD
      if (Threading::cdoLockIO) pthread_mutex_lock(&streamMutex);
#endif
      streamReadRecordF(m_fileID, data, nmiss);
#ifdef HAVE_LIBPTHREAD
      if (Threading::cdoLockIO) pthread_mutex_unlock(&streamMutex);
#endif
      if (processNum == 1 && Threading::ompNumThreads == 1) timer_stop(timer_read);
    }
}

void
PstreamType::readRecord(double *data, size_t *nmiss)
{
#ifdef HAVE_LIBPTHREAD
  if (ispipe)
    {
      if (CdoDebug::PSTREAM)
        {
          MESSAGE(pipe->name.c_str(), " pstreamID ", self);
        }
      pipe->pipeReadRecord(m_vlistID, data, nmiss);
    }
  else
#endif
    {
      if (processNum == 1 && Threading::ompNumThreads == 1) timer_start(timer_read);
#ifdef HAVE_LIBPTHREAD
      if (Threading::cdoLockIO) pthread_mutex_lock(&streamMutex);
#endif
      streamReadRecord(m_fileID, data, nmiss);
#ifdef HAVE_LIBPTHREAD
      if (Threading::cdoLockIO) pthread_mutex_unlock(&streamMutex);
#endif
      if (processNum == 1 && Threading::ompNumThreads == 1) timer_stop(timer_read);
    }
}

void
PstreamType::checkDatarange(int varID, double *array, size_t nmiss)
{
  long gridsize = m_varlist[varID].gridsize;
  int datatype = m_varlist[varID].datatype;
  double missval = m_varlist[varID].missval;
  double addoffset = m_varlist[varID].addoffset;
  double scalefactor = m_varlist[varID].scalefactor;

  size_t ivals = 0;
  double arrmin, arrmax;

  if (nmiss > 0)
    {
      ivals = arrayMinMaxMV(gridsize, array, missval, &arrmin, &arrmax);
    }
  else
    {
      arrayMinMax(gridsize, array, &arrmin, &arrmax);
      ivals = gridsize;
    }

  if (ivals > 0)
    {
      double smin = (arrmin - addoffset) / scalefactor;
      double smax = (arrmax - addoffset) / scalefactor;

      if (datatype == CDI_DATATYPE_INT8 || datatype == CDI_DATATYPE_UINT8 || datatype == CDI_DATATYPE_INT16
          || datatype == CDI_DATATYPE_UINT16)
        {
          smin = (int) lround(smin);
          smax = (int) lround(smax);
        }

      double vmin = 0, vmax = 0;

      // clang-format off
      if      ( datatype == CDI_DATATYPE_INT8   ) { vmin =        -128.; vmax =        127.; }
      else if ( datatype == CDI_DATATYPE_UINT8  ) { vmin =           0.; vmax =        255.; }
      else if ( datatype == CDI_DATATYPE_INT16  ) { vmin =      -32768.; vmax =      32767.; }
      else if ( datatype == CDI_DATATYPE_UINT16 ) { vmin =           0.; vmax =      65535.; }
      else if ( datatype == CDI_DATATYPE_INT32  ) { vmin = -2147483648.; vmax = 2147483647.; }
      else if ( datatype == CDI_DATATYPE_UINT32 ) { vmin =           0.; vmax = 4294967295.; }
      else if ( datatype == CDI_DATATYPE_FLT32  ) { vmin = -3.40282e+38; vmax = 3.40282e+38; }
      else                                        { vmin =     -1.e+300; vmax =     1.e+300; }
      // clang-format on

      if (smin < vmin || smax > vmax)
        cdoWarning("Some data values (min=%g max=%g) are outside the\n"
                   "    valid range (%g - %g) of the used output precision!\n"
                   "    Use the CDO option%s -b 64 to increase the output precision.",
                   smin, smax, vmin, vmax, (datatype == CDI_DATATYPE_FLT32) ? "" : " -b 32 or");
    }
}

void
PstreamType::writeRecord(double *data, size_t nmiss)
{

#ifdef HAVE_LIBPTHREAD
  if (ispipe)
    {
      if (CdoDebug::PSTREAM)
        {
          MESSAGE(pipe->name.c_str(), " pstreamID ", self);
        }
      pipe->pipeWriteRecord(data, nmiss);
    }
  else
#endif
    {
      int varID = m_varID;
      if (processNum == 1 && Threading::ompNumThreads == 1) timer_start(timer_write);

      if (m_varlist)
        if (m_varlist[varID].check_datarange) checkDatarange(varID, data, nmiss);

#ifdef HAVE_LIBPTHREAD
      if (Threading::cdoLockIO) pthread_mutex_lock(&streamMutex);
#endif
      streamWriteRecord(m_fileID, data, nmiss);
#ifdef HAVE_LIBPTHREAD
      if (Threading::cdoLockIO) pthread_mutex_unlock(&streamMutex);
#endif

      if (processNum == 1 && Threading::ompNumThreads == 1) timer_stop(timer_write);
    }
}

void
PstreamType::writeRecordF(float *data, size_t nmiss)
{
#ifdef HAVE_LIBPTHREAD
  if (ispipe)
    {
      if (CdoDebug::PSTREAM)
        {
          MESSAGE(pipe->name.c_str(), " pstreamID ", self);
        }
      cdoAbort("pipeWriteRecord not implemented for memtype float!");
    }
  else
#endif
    {
      if (processNum == 1 && Threading::ompNumThreads == 1) timer_start(timer_write);
#ifdef HAVE_LIBPTHREAD
      if (Threading::cdoLockIO) pthread_mutex_lock(&streamMutex);
#endif
      streamWriteRecordF(m_fileID, data, nmiss);
#ifdef HAVE_LIBPTHREAD
      if (Threading::cdoLockIO) pthread_mutex_unlock(&streamMutex);
#endif
      if (processNum == 1 && Threading::ompNumThreads == 1) timer_stop(timer_write);
    }
}

int
PstreamType::inqTimestep(int p_tsID)
{
  int nrecs = 0;

  if (processNum == 1 && Threading::ompNumThreads == 1) timer_start(timer_read);
#ifdef HAVE_LIBPTHREAD
  if (Threading::cdoLockIO) pthread_mutex_lock(&streamMutex);
#endif
  nrecs = streamInqTimestep(m_fileID, p_tsID);
#ifdef HAVE_LIBPTHREAD
  if (Threading::cdoLockIO) pthread_mutex_unlock(&streamMutex);
#endif
  if (processNum == 1 && Threading::ompNumThreads == 1) timer_stop(timer_read);

  if (nrecs == 0 && mfiles && (nfiles < mfiles))
    {
      int nfile = nfiles;
      std::string filename;
      int fileID;
      int vlistIDold, vlistIDnew;

      tsID0 += p_tsID;

      vlistIDold = vlistDuplicate(streamInqVlist(m_fileID));
      streamClose(m_fileID);

      filename = m_mfnames[nfile];
      nfiles++;

#ifdef HAVE_LIBPTHREAD
      if (Threading::cdoLockIO)
        pthread_mutex_lock(&streamMutex);
      else
        pthread_mutex_lock(&streamOpenReadMutex);
#endif
      if (cdoVerbose) cdoPrint("Continuation file: %s", filename.c_str());

      if (processNum == 1 && Threading::ompNumThreads == 1) timer_start(timer_read);
      fileID = streamOpenRead(filename.c_str());
      vlistIDnew = streamInqVlist(fileID);
      if (processNum == 1 && Threading::ompNumThreads == 1) timer_stop(timer_read);

      vlistCompare(vlistIDold, vlistIDnew, CMP_HRD);
      vlistDestroy(vlistIDold);
#ifdef HAVE_LIBPTHREAD
      if (Threading::cdoLockIO)
        pthread_mutex_unlock(&streamMutex);
      else
        pthread_mutex_unlock(&streamOpenReadMutex);
#endif
      if (fileID < 0) cdiOpenError(fileID, "Open failed on >%s<", filename.c_str());

      m_name = filename;
      m_fileID = fileID;

      if (processNum == 1 && Threading::ompNumThreads == 1) timer_start(timer_read);
#ifdef HAVE_LIBPTHREAD
      if (Threading::cdoLockIO) pthread_mutex_lock(&streamMutex);
#endif
      nrecs = streamInqTimestep(m_fileID, 0);
#ifdef HAVE_LIBPTHREAD
      if (Threading::cdoLockIO) pthread_mutex_unlock(&streamMutex);
#endif
      if (processNum == 1 && Threading::ompNumThreads == 1) timer_stop(timer_read);
    }

  if (p_tsID == 0 && cdoDefaultTimeType != CDI_UNDEFID) taxisDefType(vlistInqTaxis(m_vlistID), cdoDefaultTimeType);
  return nrecs;
}

void
PstreamType::defTimestep(int p_tsID)
{
#ifdef HAVE_LIBPTHREAD
  if (ispipe)
    {
      if (CdoDebug::PSTREAM)
        {
          MESSAGE(pipe->name.c_str(), " pstreamID ", self);
        }
      pipe->pipeDefTimestep(m_vlistID, p_tsID);
    }
  else
#endif
    {
      if (p_tsID == 0 && cdoDefaultTimeType != CDI_UNDEFID)
        {
          int taxisID, vlistID;
          vlistID = m_vlistID;
          taxisID = vlistInqTaxis(vlistID);
          taxisDefType(taxisID, cdoDefaultTimeType);
        }

      if (processNum == 1 && Threading::ompNumThreads == 1) timer_start(timer_write);
/* don't use sync -> very slow on GPFS */
//  if ( p_tsID > 0 ) streamSync(fileID);
#ifdef HAVE_LIBPTHREAD
      if (Threading::cdoLockIO) pthread_mutex_lock(&streamMutex);
#endif
      streamDefTimestep(m_fileID, p_tsID);
#ifdef HAVE_LIBPTHREAD
      if (Threading::cdoLockIO) pthread_mutex_unlock(&streamMutex);
#endif
      if (processNum == 1 && Threading::ompNumThreads == 1) timer_stop(timer_write);
    }
}

void
PstreamType::copyRecord(PstreamType *src)
{
  if (ispipe || src->ispipe) cdoAbort("This operator can't be combined with other operators!");

#ifdef HAVE_LIBPTHREAD
  if (Threading::cdoLockIO) pthread_mutex_lock(&streamMutex);
#endif
  streamCopyRecord(m_fileID, src->m_fileID);
#ifdef HAVE_LIBPTHREAD
  if (Threading::cdoLockIO) pthread_mutex_unlock(&streamMutex);
#endif
}

int
PstreamType::inqFileType()
{
  int filetype;
#ifdef HAVE_LIBPTHREAD
  if (ispipe)
    {
      filetype = m_filetype;
    }
  else
#endif
    {
      filetype = streamInqFiletype(m_fileID);
    }
  return filetype;
}

int
PstreamType::inqByteorder()
{
  int byteorder;

#ifdef HAVE_LIBPTHREAD
  if (ispipe)
    byteorder = m_filetype;
  else
#endif
    byteorder = streamInqByteorder(m_fileID);

  return byteorder;
}
//- - - - - - - - - - - - - - - - - - - - - - - - has to be moved
void
cdoVlistCopyFlag(int vlistID2, int vlistID1)
{
#ifdef HAVE_LIBPTHREAD
  pthread_mutex_lock(&streamMutex);
#endif

  vlistCopyFlag(vlistID2, vlistID1);

#ifdef HAVE_LIBPTHREAD
  pthread_mutex_unlock(&streamMutex);
#endif
}

void
openLock(void)
{
#ifdef HAVE_LIBPTHREAD
  if (Threading::cdoLockIO)
    pthread_mutex_lock(&streamMutex);
  else
    pthread_mutex_lock(&streamOpenReadMutex);
#endif
}

void
openUnlock(void)
{
#ifdef HAVE_LIBPTHREAD
  if (Threading::cdoLockIO)
    pthread_mutex_unlock(&streamMutex);
  else
    pthread_mutex_unlock(&streamOpenReadMutex);
#endif
}

void
pstreamCloseAll()
{
  for (auto pstream_iter : _pstream_map)
    {
      if (pstream_iter.second.m_fileID != CDI_UNDEFID)
        {
          if (CdoDebug::PSTREAM)
            MESSAGE("Close file ", pstream_iter.second.m_name, " id ", pstream_iter.second.m_fileID);
          streamClose(pstream_iter.second.m_fileID);
        }
    }
  _pstream_map.clear();
}

static void
set_comp(int fileID, int filetype)
{
  if (Options::cdoCompress)
    {
      if (filetype == CDI_FILETYPE_GRB)
        {
          Options::cdoCompType = CDI_COMPRESS_SZIP;
          Options::cdoCompLevel = 0;
        }
      else if (filetype == CDI_FILETYPE_NC4 || filetype == CDI_FILETYPE_NC4C)
        {
          Options::cdoCompType = CDI_COMPRESS_ZIP;
          Options::cdoCompLevel = 1;
        }
    }

  if (Options::cdoCompType != CDI_COMPRESS_NONE)
    {
      streamDefCompType(fileID, Options::cdoCompType);
      streamDefCompLevel(fileID, Options::cdoCompLevel);

      if (Options::cdoCompType == CDI_COMPRESS_SZIP
          && (filetype != CDI_FILETYPE_GRB && filetype != CDI_FILETYPE_GRB2 && filetype != CDI_FILETYPE_NC4
              && filetype != CDI_FILETYPE_NC4C))
        cdoWarning("SZIP compression not available for non GRIB/NetCDF4 data!");

      if (Options::cdoCompType == CDI_COMPRESS_JPEG && filetype != CDI_FILETYPE_GRB2)
        cdoWarning("JPEG compression not available for non GRIB2 data!");

      if (Options::cdoCompType == CDI_COMPRESS_ZIP && (filetype != CDI_FILETYPE_NC4 && filetype != CDI_FILETYPE_NC4C))
        cdoWarning("Deflate compression not available for non NetCDF4 data!");
    }
}

static int
pstreamFindID(const char *p_name)
{
  std::string cur_name;
  for (auto map_pair : _pstream_map)
    {
      cur_name = map_pair.second.m_name;
      if (!(cur_name.empty()))
        {
          if (cur_name.compare(p_name) == 0)
            {
              return map_pair.first;
            }
        }
    }
  return -1;
}

void
pstreamDebug(int debug)
{
  CdoDebug::PSTREAM = debug;
}

static void
pstream_delete_entry(PstreamType *pstreamptr)
{
  int idx = pstreamptr->self;

  PSTREAM_LOCK();

  _pstream_map.erase(idx);

  PSTREAM_UNLOCK();

  if (CdoDebug::PSTREAM) MESSAGE("Removed idx ", idx, " from pstream list");
}

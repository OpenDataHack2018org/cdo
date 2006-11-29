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

#if  defined  (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>

FILE *popen(const char *command, const char *type);
int pclose(FILE *stream);

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "modules.h"
#include "pstream_int.h"
#include "cdo_int.h"
#include "util.h"
#include "process.h"
#include "pipe.h"
#include "error.h"
#include "dmemory.h"


extern int timer_read, timer_write;

static int PSTREAM_Debug = 0;

#define  MAX_PSTREAMS  1024

static int _pstream_max = MAX_PSTREAMS;

static void pstream_initialize(void);

static int _pstream_init = FALSE;

#if  defined  (HAVE_LIBPTHREAD)
#include <pthread.h>
#include "pthread_debug.h"

static pthread_mutex_t streamOpenReadMutex  = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t streamOpenWriteMutex = PTHREAD_MUTEX_INITIALIZER;

static pthread_once_t _pstream_init_thread = PTHREAD_ONCE_INIT;
static pthread_mutex_t _pstream_mutex;

#  define PSTREAM_LOCK           pthread_mutex_lock(&_pstream_mutex);
#  define PSTREAM_UNLOCK         pthread_mutex_unlock(&_pstream_mutex);
#  define PSTREAM_INIT                               \
   if ( _pstream_init == FALSE ) pthread_once(&_pstream_init_thread, pstream_initialize);

#else

#  define PSTREAM_LOCK
#  define PSTREAM_UNLOCK
#  define PSTREAM_INIT                               \
   if ( _pstream_init == FALSE ) pstream_initialize();

#endif


typedef struct _pstreamPtrToIdx {
  int idx;
  PSTREAM *ptr;
  struct _pstreamPtrToIdx *next;
} pstreamPtrToIdx;


static pstreamPtrToIdx *_pstreamList;
static pstreamPtrToIdx *_pstreamAvail = 0;


static void pstream_list_new(void)
{
  static char func[] = "pstream_list_new";

  _pstreamList = (pstreamPtrToIdx *) malloc(_pstream_max*sizeof(pstreamPtrToIdx));
}


static void pstream_list_delete(void)
{
  static char func[] = "pstream_list_delete";

  if ( _pstreamList ) free(_pstreamList);
}


static void pstream_init_pointer(void)
{
  int  i;
  
  for ( i = 0; i < _pstream_max; i++ )
    {
      _pstreamList[i].next = _pstreamList + i + 1;
      _pstreamList[i].idx  = i;
      _pstreamList[i].ptr  = 0;
    }

  _pstreamList[_pstream_max-1].next = 0;

  _pstreamAvail = _pstreamList;
}


static PSTREAM *pstream_to_pointer(int idx)
{
  static char func[] = "pstream_to_pointer";
  PSTREAM *pstreamptr = NULL;

  PSTREAM_INIT

  if ( idx >= 0 && idx < _pstream_max )
    {
      PSTREAM_LOCK

      pstreamptr = _pstreamList[idx].ptr;

      PSTREAM_UNLOCK
    }
  else
    Error(func, "pstream index %d undefined!", idx);

  return (pstreamptr);
}


/* Create an index from a pointer */
static int pstream_from_pointer(PSTREAM *ptr)
{
  static char func[] = "pstream_from_pointer";
  int      idx = -1;
  pstreamPtrToIdx *newptr;

  if ( ptr )
    {
      PSTREAM_LOCK

      if ( _pstreamAvail )
	{
	  newptr        = _pstreamAvail;
	  _pstreamAvail = _pstreamAvail->next;
	  newptr->next  = 0;
	  idx	        = newptr->idx;
	  newptr->ptr   = ptr;
      
	  if ( PSTREAM_Debug )
	    Message(func, "Pointer %p has idx %d from pstream list", ptr, idx);
	}
      else
	Warning(func, "Too many open pstreams (limit is %d)!", _pstream_max);

      PSTREAM_UNLOCK
    }
  else
    Error(func, "Internal problem (pointer %p undefined)", ptr);

  return (idx);
}


static void pstream_init_entry(PSTREAM *pstreamptr)
{
  pstreamptr->self       = pstream_from_pointer(pstreamptr);

  pstreamptr->isopen     = TRUE;
  pstreamptr->ispipe     = FALSE;
  pstreamptr->fileID     = -1;
  pstreamptr->vlistID    = -1;
  pstreamptr->tsID       = -1;
  pstreamptr->filetype   = -1;
  pstreamptr->name       = NULL;
#if  defined  (HAVE_LIBPTHREAD)
  pstreamptr->pipe       = NULL;
  pstreamptr->rthreadID  = 0;
  pstreamptr->wthreadID  = 0;
#endif
  pstreamptr->tsID0      = 0;
  pstreamptr->mfiles     = 0;
  pstreamptr->nfiles     = 0;
  pstreamptr->mfnames    = NULL;
}


static PSTREAM *pstream_new_entry(void)
{
  static char func[] = "pstream_new_entry";
  PSTREAM *pstreamptr;

  pstreamptr = (PSTREAM *) malloc(sizeof(PSTREAM));

  if ( pstreamptr ) pstream_init_entry(pstreamptr);

  return (pstreamptr);
}


static void pstream_delete_entry(PSTREAM *pstreamptr)
{
  static char func[] = "pstream_delete_entry";
  int idx;

  idx = pstreamptr->self;

  PSTREAM_LOCK

  free(pstreamptr);

  _pstreamList[idx].next = _pstreamAvail;
  _pstreamList[idx].ptr  = 0;
  _pstreamAvail   	 = &_pstreamList[idx];

  PSTREAM_UNLOCK

  if ( PSTREAM_Debug )
    Message(func, "Removed idx %d from pstream list", idx);
}


static void pstream_initialize(void)
{
  static char func[] = "pstream_initialize";
  char *env;

#if  defined  (HAVE_LIBPTHREAD)
  /* initialize global API mutex lock */
  pthread_mutex_init(&_pstream_mutex, NULL);
#endif

  env = getenv("PSTREAM_DEBUG");
  if ( env ) PSTREAM_Debug = atoi(env);

  env = getenv("PSTREAM_MAX");
  if ( env ) _pstream_max = atoi(env);

  if ( PSTREAM_Debug )
    Message(func, "PSTREAM_MAX = %d", _pstream_max);

  pstream_list_new();
  atexit(pstream_list_delete);

  pstream_init_pointer();

  _pstream_init = TRUE;
}


static int pstreamFindID(const char *name)
{
  int pstreamID;
  PSTREAM *pstreamptr;

  for ( pstreamID = 0; pstreamID < _pstream_max; pstreamID++ )
    {
      pstreamptr = pstream_to_pointer(pstreamID);

      if ( pstreamptr )
	if ( pstreamptr->name )
	  if ( strcmp(pstreamptr->name, name) == 0 ) break;
    }

  if ( pstreamID == _pstream_max ) pstreamID = -1;

  return (pstreamID);
}


int pstreamIsPipe(int pstreamID)
{
  PSTREAM *pstreamptr;

  pstreamptr = pstream_to_pointer(pstreamID);

  return (pstreamptr->ispipe);
}


int pstreamOpenRead(const char *argument)
{
  static char func[] = "pstreamOpenRead";
  char *operatorArg = NULL;
  char *operatorName = NULL;
  int ispipe = FALSE;
  int fileID;
  int pstreamID;
  PSTREAM *pstreamptr;

  PSTREAM_INIT

  pstreamptr = pstream_new_entry();
  if ( ! pstreamptr ) Error(func, "No memory");

  pstreamID = pstreamptr->self;

  ispipe = argument[0] == '-';

  if ( ispipe )
    {
#if  defined  (HAVE_LIBPTHREAD)
      char *newarg;
      char *pipename = (char *) malloc(10);
      int rval;
      pthread_t thrID;
      pthread_attr_t attr;
      size_t len;
      size_t stacksize;

      operatorArg = getOperator(argument);
      operatorName = getOperatorName(operatorArg);
      free(operatorArg);

      len = strlen(argument);
      newarg = (char *) malloc(len+11);
      strcpy(newarg, argument);
      sprintf(pipename, "(pipe%d.%d)", processSelf() + 1, processInqChildNum() + 1);
      newarg[len] = ' ';
      strcpy(&newarg[len+1], pipename);

      pstreamptr->ispipe = TRUE;
      pstreamptr->name   = pipename;
      pstreamptr->rthreadID = pthread_self();
      pstreamptr->pipe   = pipeNew();
 
      if ( ! cdoSilentMode )
	fprintf(stderr, "%s: Started child process \"%s\".\n", processInqPrompt(), newarg+1);

      /*
      if ( PTHREAD_Debug )
	{
	}
      */
      pthread_attr_init(&attr);
      /*      pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM); */
      pthread_attr_getstacksize(&attr, &stacksize);
      if ( stacksize < 2097152 )
	{
	  stacksize = 2097152;
	  pthread_attr_setstacksize(&attr, stacksize);
	}
      rval = pthread_create(&thrID, &attr, operatorModule(operatorName), newarg);
      if ( rval != 0 )
	{
	  errno = rval;
	  SysError(func, "pthread_create failed for '%s'\n", newarg+1);
	}

      /* free(operatorName); */
      processAddStream(pstreamID);
      /*      pipeInqInfo(pstreamID); */
      if ( PSTREAM_Debug ) Message(func, "pipe %s", pipename);
#else
      cdoAbort("Cannot use pipes, pthread library not available!");
#endif
    }
  else
    {
      extern int cdoDefaultFileType/*, cdoDefaultInstID*/;
      size_t len, i;
      int nfiles = 1, j;
      char *filename = NULL;
      const char *pch;

      len = strlen(argument);

      for ( i = 0; i < len; i++ )
	if ( argument[i] == ':' ) break;

      if ( i < len )
	{
	  pch = &argument[i+1];
	  len -= (i+1);
	  if ( len && ( strncmp(argument, "filelist:", i) == 0 || 
			strncmp(argument, "flist:", i) == 0 ) )
	    {
	      for ( i = 0; i < len; i++ ) if ( pch[i] == ',' ) nfiles++;

	      if ( nfiles == 1 )
		{
		  char line[4096];
		  FILE *fp, *fp2;
		  fp = fopen(pch, "r");
		  if ( fp == NULL ) cdoAbort("Open failed on %s", pch);

		  if ( cdoVerbose )
		    cdoPrint("Reading file names from %s", pch);

		  /* find number of files */
		  nfiles = 0;
		  while ( readline(fp, line, 4096) )
		    {
		      if ( line[0] == '#' || line[0] == '\0' ||
			   line[0] == ' ' ) continue;

		      fp2 = fopen(line, "r" );
		      if ( fp2 == NULL ) cdoAbort("Open failed on %s", line);
		      fclose(fp2);
		      nfiles++;
		      if ( cdoVerbose )
			cdoPrint("File number %d is %s", nfiles, line);
		    }

		  if ( nfiles == 0 ) cdoAbort("No imput file found in %s", pch);

		  pstreamptr->mfiles = nfiles;
		  pstreamptr->mfnames = (char **) malloc(nfiles*sizeof(char *));
		  
		  rewind(fp);

		  nfiles = 0;
		  while ( readline(fp, line, 4096) )
		    {
		      if ( line[0] == '#' || line[0] == '\0' ||
			   line[0] == ' ' ) continue;

		      pstreamptr->mfnames[nfiles] = strdupx(line);
		      nfiles++;
		    }

		  fclose(fp);
		}
	      else
		{
		  char line[65536];

		  pstreamptr->mfiles = nfiles;
		  pstreamptr->mfnames = (char **) malloc(nfiles*sizeof(char *));
		  
		  strcpy(line, pch);
		  for ( i = 0; i < len; i++ ) if ( line[i] == ',' ) line[i] = 0;

		  i = 0;
		  for ( j = 0; j < nfiles; j++ )
		    {
		      pstreamptr->mfnames[j] = strdupx(&line[i]);
		      i += strlen(&line[i]) + 1;
		    }
		}
	    }
	  else if ( len && strncmp(argument, "ls:", i) == 0 )
	    {
	      char line[4096];
	      char command[4096];
	      char *fnames[16384];
	      FILE *pfp;

	      strcpy(command, "ls ");
	      strcat(command, pch);

	      pfp = popen(command, "r");
	      if ( pfp == 0 )
		SysError(func, "popen %s failed", command);

	      nfiles = 0;
	      while ( readline(pfp, line, 4096) )
		{
		  if ( nfiles >= 16384 ) cdoAbort("Too many input files (limit: 16384)");
		  fnames[nfiles++] = strdupx(line);
		}

	      pclose(pfp);

	      pstreamptr->mfiles = nfiles;
	      pstreamptr->mfnames = (char **) malloc(nfiles*sizeof(char *));

	      for ( j = 0; j < nfiles; j++ )
		pstreamptr->mfnames[j] = fnames[j];
	    }
	}

      if ( pstreamptr->mfiles )
	{
	  len = strlen(pstreamptr->mfnames[0]);
	  filename = (char *) malloc(len+1);
	  strcpy(filename, pstreamptr->mfnames[0]);
	  pstreamptr->nfiles = 1;
	}
      else
	{
	  len = strlen(argument);

	  if ( cdoExpMode == CDO_EXP_REMOTE )
	    {
	      char datapath[] = "/scratch/localA/m214003/data/";
	      len += strlen(datapath);

	      filename = (char *) malloc(len+1);

	      strcpy(filename, datapath);
	      strcat(filename, argument);
	    }
	  else
	    {
	      filename = (char *) malloc(len+1);

	      strcpy(filename, argument);
	    }
	}

      if ( PSTREAM_Debug ) Message(func, "file %s", filename);

#if  defined  (HAVE_LIBPTHREAD)
      pthread_mutex_lock(&streamOpenReadMutex);
#endif
      fileID = streamOpenRead(filename);
      if ( fileID < 0 ) cdiError(fileID, "Open failed on >%s<", filename);

      if ( cdoDefaultFileType == CDI_UNDEFID )
	cdoDefaultFileType = streamInqFiletype(fileID);
      /*
      if ( cdoDefaultInstID == CDI_UNDEFID )
	cdoDefaultInstID = streamInqInstID(fileID);
      */
      cdoInqHistory(fileID);
#if  defined  (HAVE_LIBPTHREAD)
      pthread_mutex_unlock(&streamOpenReadMutex);
#endif

      pstreamptr->mode   = 'r';
      pstreamptr->name   = filename;
      pstreamptr->fileID = fileID;
    }

  return (pstreamID);
}


int pstreamOpenWrite(const char *argument, int filetype)
{
  static char func[] = "pstreamOpenWrite";
  int fileID;
  int pstreamID;
  int ispipe;
  PSTREAM *pstreamptr;

  PSTREAM_INIT

  ispipe = strncmp(argument, "(pipe", 5) == 0;

  if ( ispipe )
    {
#if  defined  (HAVE_LIBPTHREAD)
      if ( PSTREAM_Debug ) Message(func, "pipe %s", argument);
      pstreamID = pstreamFindID(argument);
      if ( pstreamID == -1 )
	Error(func, "%s not open", argument);

      pstreamptr = pstream_to_pointer(pstreamID);

      pstreamptr->wthreadID = pthread_self();
      pstreamptr->filetype = filetype;
      processAddStream(pstreamID);
#endif
    }
  else
    {
      /* extern int cdoDefaultInstID; */
      char *filename = (char *) malloc(strlen(argument)+1);

      pstreamptr = pstream_new_entry();
      if ( ! pstreamptr ) Error(func, "No memory");

      pstreamID = pstreamptr->self;
  
      if ( PSTREAM_Debug ) Message(func, "file %s", argument);

      if ( filetype == CDI_UNDEFID ) filetype = FILETYPE_GRB;

#if  defined  (HAVE_LIBPTHREAD)
      pthread_mutex_lock(&streamOpenWriteMutex);
#endif
      if ( cdoTimer ) timer_start(timer_write);
      fileID = streamOpenWrite(argument, filetype);
      if ( cdoTimer ) timer_stop(timer_write);
#if  defined  (HAVE_LIBPTHREAD)
      pthread_mutex_unlock(&streamOpenWriteMutex);
#endif
      if ( fileID < 0 ) cdiError(fileID, "Open failed on %s", argument);

      cdoDefHistory(fileID, commandLine());

      if ( cdoDefaultByteorder != CDI_UNDEFID )
	streamDefByteorder(fileID, cdoDefaultByteorder);
      /*
      if ( cdoDefaultInstID != CDI_UNDEFID )
	streamDefInstID(fileID, cdoDefaultInstID);
      */
      strcpy(filename, argument);

      pstreamptr->mode   = 'w';
      pstreamptr->name   = filename;
      pstreamptr->fileID = fileID;
    }

  return (pstreamID);
}


int pstreamOpenAppend(const char *argument)
{
  static char func[] = "pstreamOpenAppend";
  int fileID;
  int pstreamID = -1;
  int ispipe;
  PSTREAM *pstreamptr;

  ispipe = strncmp(argument, "(pipe", 5) == 0;

  if ( ispipe )
    {
      if ( PSTREAM_Debug ) Message(func, "pipe %s", argument);
      cdoAbort("this operator doesn't work with pipes!");
    }
  else
    {
      char *filename = (char *) malloc(strlen(argument)+1);

      pstreamptr = pstream_new_entry();
      if ( ! pstreamptr ) Error(func, "No memory");

      pstreamID = pstreamptr->self;
  
      if ( PSTREAM_Debug ) Message(func, "file %s", argument);

      if ( cdoTimer ) timer_start(timer_write);
      fileID = streamOpenAppend(argument);
      if ( cdoTimer ) timer_stop(timer_write);
      if ( fileID < 0 ) cdiError(fileID, "Open failed on %s", argument);
      /*
      cdoInqHistory(fileID);
      cdoDefHistory(fileID, commandLine());
      */
      strcpy(filename, argument);
      pstreamptr->name   = filename;
      pstreamptr->fileID = fileID;
    }

  return (pstreamID);
}


void pstreamClose(int pstreamID)
{
  static char func[] = "pstreamClose";
  PSTREAM *pstreamptr;

  pstreamptr = pstream_to_pointer(pstreamID);

  if ( pstreamptr == NULL )
    Error(func, "Internal problem stream %d not open!", pstreamID);

  if ( pstreamptr->ispipe )
    {
#if  defined  (HAVE_LIBPTHREAD)
      PIPE *pipe;
      int lread = FALSE, lwrite = FALSE;
      pthread_t threadID = pthread_self();

      if      ( pthread_equal(threadID, pstreamptr->rthreadID) ) lread  = TRUE;
      else if ( pthread_equal(threadID, pstreamptr->wthreadID) ) lwrite = TRUE;
      else Error(func, "Internal problem! Close pipe %s", pstreamptr->name);

      if ( lread )
	{
	  pipe = pstreamptr->pipe;
	  pthread_mutex_lock(pipe->mutex);
	  pipe->EOP = TRUE;
	  if ( PSTREAM_Debug )
	    Message(func, "%s read closed", pstreamptr->name);
	  pthread_mutex_unlock(pipe->mutex);     
	  pthread_cond_signal(pipe->tsDef);
	  pthread_cond_signal(pipe->tsInq);
	 
	  pthread_cond_signal(pipe->recInq);
	 
	  pthread_mutex_lock(pipe->mutex);
	  pstreamptr->isopen = FALSE;
	  pthread_mutex_unlock(pipe->mutex);     
	  pthread_cond_signal(pipe->isclosed);

	  pthread_join(pstreamptr->wthreadID, NULL);

	  pipeDelete(pipe);
	  pstream_delete_entry(pstreamptr);
	}
      else
	{
	  pipe = pstreamptr->pipe;
	  pthread_mutex_lock(pipe->mutex);
	  pipe->EOP = TRUE;
	  if ( PSTREAM_Debug )
	    Message(func, "%s write closed", pstreamptr->name);
	  pthread_mutex_unlock(pipe->mutex);     
	  pthread_cond_signal(pipe->tsDef);
	  pthread_cond_signal(pipe->tsInq);

	  pthread_mutex_lock(pipe->mutex);
	  while ( pstreamptr->isopen )
	    {
	      if ( PSTREAM_Debug )
		Message(func, "wait of read close");
	      pthread_cond_wait(pipe->isclosed, pipe->mutex);
	    }
	  pthread_mutex_unlock(pipe->mutex);
	}
#else
      cdoAbort("Cannot use pipes, pthread library not available!");
#endif
    }
  else
    {
      if ( PSTREAM_Debug )
	Message(func, "%s fileID %d\n", pstreamptr->name, pstreamptr->fileID);

      if ( pstreamptr->mode == 'r' )
	processAddNvals(streamNvals(pstreamptr->fileID));

      streamClose(pstreamptr->fileID);

      if ( cdoExpMode == CDO_EXP_REMOTE )
	{
	  if ( pstreamptr->mode == 'w' )
	    {
	      extern const char *cdojobfiles;
	      FILE *fp = fopen(cdojobfiles, "a");
	      fprintf(fp, "%s\n", pstreamptr->name);
	      fclose(fp);
	    }
	}

      pstream_delete_entry(pstreamptr);
    }
}


int pstreamInqVlist(int pstreamID)
{
  int vlistID;
  PSTREAM *pstreamptr;

  pstreamptr = pstream_to_pointer(pstreamID);

#if  defined  (HAVE_LIBPTHREAD)
  if ( pstreamptr->ispipe )
    vlistID = pipeInqVlist(pstreamptr);
  else
#endif
    {
      extern int cdoDefaultTimeType;

      if ( cdoTimer ) timer_start(timer_read);
      vlistID = streamInqVlist(pstreamptr->fileID);
      if ( cdoTimer ) timer_stop(timer_read);

      if ( cdoDefaultTimeType != CDI_UNDEFID )
	taxisDefType(vlistInqTaxis(vlistID), cdoDefaultTimeType);
    }

  pstreamptr->vlistID = vlistID;

  processDefVarNum(vlistNvars(vlistID), pstreamID);

  return (vlistID);
}


const char *cdoComment(void)
{
  static char comment[256];
  static int init = 0;
  int size = 0;
  extern char CDO_Version[];

  if ( ! init )
    {
      init = 1;

      size = strlen(CDO_Version);

      strncat(comment, CDO_Version, size);
    }

  return (comment);
}


void pstreamDefVlist(int pstreamID, int vlistID)
{
  PSTREAM *pstreamptr;

  pstreamptr = pstream_to_pointer(pstreamID);

#if  defined  (HAVE_LIBPTHREAD)
  if ( pstreamptr->ispipe )
    {
      /*    pipeDefVlist(pstreamptr, vlistID);*/
      pipeDefVlist(pstreamptr, vlistDuplicate(vlistID));
    }
  else
#endif
    {
      extern int cdoDefaultDataType;
      if ( cdoDefaultDataType != CDI_UNDEFID )
	{
	  int varID, nvars = vlistNvars(vlistID);

	  for ( varID = 0; varID < nvars; varID++ )
	    vlistDefVarDatatype(vlistID, varID, cdoDefaultDataType);

	  if ( cdoDefaultDataType == DATATYPE_FLT64 )
	    {
	      for ( varID = 0; varID < nvars; varID++ )
		{
		  vlistDefVarAddoffset(vlistID, varID, 0.0);
		  vlistDefVarScalefactor(vlistID, varID, 1.0);
		}
	    }
	}
      vlistDefAttribute(vlistID, "CDO", cdoComment());
      if ( cdoTimer ) timer_start(timer_write);
      streamDefVlist(pstreamptr->fileID, vlistID);
      if ( cdoTimer ) timer_stop(timer_write);

      pstreamptr->vlistID = vlistID; /* used for -r/-a */
    }
}


int pstreamInqRecord(int pstreamID, int *varID, int *levelID)
{
  PSTREAM *pstreamptr;

  pstreamptr = pstream_to_pointer(pstreamID);

#if  defined  (HAVE_LIBPTHREAD)
  if ( pstreamptr->ispipe )
    pipeInqRecord(pstreamptr, varID, levelID);
  else
#endif
    {
      if ( cdoTimer ) timer_start(timer_read);
      streamInqRecord(pstreamptr->fileID, varID, levelID);
      if ( cdoTimer ) timer_stop(timer_read);
    }

  return (0);
}


void pstreamDefRecord(int pstreamID, int varID, int levelID)
{
  PSTREAM *pstreamptr;

  pstreamptr = pstream_to_pointer(pstreamID);

#if  defined  (HAVE_LIBPTHREAD)
  if ( pstreamptr->ispipe )
    pipeDefRecord(pstreamptr, varID, levelID);
  else
#endif
    {
      if ( cdoTimer ) timer_start(timer_write);
      streamDefRecord(pstreamptr->fileID, varID, levelID);
      if ( cdoTimer ) timer_stop(timer_write);
    }
}


void pstreamReadRecord(int pstreamID, double *data, int *nmiss)
{
  PSTREAM *pstreamptr;

  pstreamptr = pstream_to_pointer(pstreamID);

#if  defined  (HAVE_LIBPTHREAD)
  if ( pstreamptr->ispipe )
    pipeReadRecord(pstreamptr, data, nmiss);
  else
#endif
    {
      if ( cdoTimer ) timer_start(timer_read);
      streamReadRecord(pstreamptr->fileID, data, nmiss);
      if ( cdoTimer ) timer_stop(timer_read);
    }
}


void pstreamWriteRecord(int pstreamID, double *data, int nmiss)
{
  PSTREAM *pstreamptr;

  pstreamptr = pstream_to_pointer(pstreamID);

#if  defined  (HAVE_LIBPTHREAD)
  if ( pstreamptr->ispipe )
    pipeWriteRecord(pstreamptr, data, nmiss);
  else
#endif
    {
      if ( cdoTimer ) timer_start(timer_write);
      streamWriteRecord(pstreamptr->fileID, data, nmiss);
      if ( cdoTimer ) timer_stop(timer_write);
    }
}


int pstreamInqTimestep(int pstreamID, int tsID)
{
  static char func[] = "pstreamInqTimestep";
  int nrecs = 0;
  PSTREAM *pstreamptr;

  pstreamptr = pstream_to_pointer(pstreamID);

#if  defined  (HAVE_LIBPTHREAD)
  if ( pstreamptr->ispipe )
    nrecs = pipeInqTimestep(pstreamptr, tsID);
  else
#endif
    {
      extern int cdoDefaultTimeType;

      if ( pstreamptr->mfiles ) tsID -= pstreamptr->tsID0;

      if ( cdoTimer ) timer_start(timer_read);
      nrecs = streamInqTimestep(pstreamptr->fileID, tsID);
      if ( cdoTimer ) timer_stop(timer_read);

      if ( nrecs == 0 && pstreamptr->mfiles &&
	   (pstreamptr->nfiles < pstreamptr->mfiles) )
	{
	  size_t len;
	  int nfile = pstreamptr->nfiles;
	  char *filename = NULL;
	  int fileID;
	  int vlistIDold, vlistIDnew;

	  pstreamptr->tsID0 += tsID;

	  vlistIDold = vlistDuplicate(streamInqVlist(pstreamptr->fileID));
	  streamClose(pstreamptr->fileID);

	  len = strlen(pstreamptr->mfnames[nfile]);
	  filename = (char *) malloc(len+1);
	  strcpy(filename, pstreamptr->mfnames[nfile]);
	  pstreamptr->nfiles++;

#if  defined  (HAVE_LIBPTHREAD)
	  pthread_mutex_lock(&streamOpenReadMutex);
#endif
	  if ( cdoVerbose ) cdoPrint("Continuation file: %s", filename);

	  if ( cdoTimer ) timer_start(timer_read);
	  fileID = streamOpenRead(filename);
	  vlistIDnew = streamInqVlist(fileID);
	  if ( cdoTimer ) timer_stop(timer_read);

	  vlistCompare(vlistIDold, vlistIDnew, func_hrd);
	  vlistDestroy(vlistIDold);
#if  defined  (HAVE_LIBPTHREAD)
	  pthread_mutex_unlock(&streamOpenReadMutex);
#endif
	  if ( fileID < 0 ) cdiError(fileID, "Open failed on >%s<", filename);

	  free(pstreamptr->name);

	  pstreamptr->name   = filename;
	  pstreamptr->fileID = fileID;

	  if ( cdoTimer ) timer_start(timer_read);
	  nrecs = streamInqTimestep(pstreamptr->fileID, 0);
	  if ( cdoTimer ) timer_stop(timer_read);
	}

      if ( cdoDefaultTimeType != CDI_UNDEFID )
	taxisDefType(vlistInqTaxis(pstreamptr->vlistID), cdoDefaultTimeType);
    }

  if ( nrecs && tsID != pstreamptr->tsID )
    {
      processDefTimesteps(pstreamID);
      pstreamptr->tsID = tsID;
    }

  return (nrecs);
}


void pstreamDefTimestep(int pstreamID, int tsID)
{
  PSTREAM *pstreamptr;

  pstreamptr = pstream_to_pointer(pstreamID);

#if  defined  (HAVE_LIBPTHREAD)
  if ( pstreamptr->ispipe )
    pipeDefTimestep(pstreamptr, tsID);
  else
#endif
    {
      extern int cdoDefaultTimeType;
      if ( cdoDefaultTimeType != CDI_UNDEFID )
	{
	  int taxisID, vlistID;
	  vlistID = pstreamptr->vlistID;
	  taxisID = vlistInqTaxis(vlistID);
      	  taxisDefType(taxisID, cdoDefaultTimeType);
	}

      if ( cdoTimer ) timer_start(timer_write);
      streamDefTimestep(pstreamptr->fileID, tsID);
      if ( cdoTimer ) timer_stop(timer_write);
    }
}


void pstreamCopyRecord(int pstreamIDdest, int pstreamIDsrc)
{
  static char func[] = "pstreamCopyRecord";
  PSTREAM *pstreamptr_dest, *pstreamptr_src;

  if ( PSTREAM_Debug )
    Message(func, "pstreamIDdest = %d  pstreamIDsrc = %d", pstreamIDdest, pstreamIDsrc);

  pstreamptr_dest = pstream_to_pointer(pstreamIDdest);
  pstreamptr_src  = pstream_to_pointer(pstreamIDsrc);

#if  defined  (HAVE_LIBPTHREAD)
  if ( pstreamptr_dest->ispipe || pstreamptr_src->ispipe )
    pipeCopyRecord(pstreamptr_dest, pstreamptr_src);
  else
#endif
    {
      streamCopyRecord(pstreamptr_dest->fileID, pstreamptr_src->fileID);
    }
}


void pstreamDebug(int debug)
{
  PSTREAM_Debug = debug;
}


void cdoInitialize(void *argument)
{
  static char func[] = "cdoInitialize";
  int processID;

  processID = processCreate();

#if  defined  (HAVE_LIBPTHREAD)
  if ( PSTREAM_Debug )
     Message(func, "process %d  thread %ld", processSelf(), pthread_self());
#endif

  processDefArgument((const char*) argument);
}


void cdoFinish(void)
{
  static char func[] = "cdoFinish";
  int processID = processSelf();
  int sindex, pstreamID;
  int nstream;
  int nvars, ntimesteps;
  char memstring[32] = {""};
  double s_utime, s_stime;
  double e_utime, e_stime;
  double c_cputime = 0, c_usertime = 0, c_systime = 0;
  double p_cputime = 0, p_usertime = 0, p_systime = 0;
  PSTREAM *pstreamptr;

#if  defined  (HAVE_LIBPTHREAD)
  if ( PSTREAM_Debug )
    Message(func, "process %d  thread %ld", processID, pthread_self());
#endif

  nvars = processInqVarNum();
  ntimesteps = processInqTimesteps();

  if ( ! cdoSilentMode )
    fprintf(stderr, "%s: Processed %d variable%s %d timestep%s.",
	    processInqPrompt(), nvars, nvars > 1 ? "s" : "",
	    ntimesteps, ntimesteps > 1 ? "s" : "");

  processStartTime(&s_utime, &s_stime);
  cdoProcessTime(&e_utime, &e_stime);

  c_usertime = e_utime - s_utime;
  c_systime  = e_stime - s_stime;
  c_cputime  = c_usertime + c_systime;

  processAccuTime(c_usertime, c_systime);

  if ( processID == 0 )
    {
      int mu[] = {'B', 'K', 'M', 'G', 'T'};
      int muindex = 0;
      long memmax;
      extern int cdoLogOff;

      memmax = memTotal();
      while ( memmax > 9999 )
	{
	  memmax /= 1024;
	  muindex++;
	}

      if ( memmax )
	sprintf(memstring, " %ld%c ", memmax, mu[muindex]);

      processEndTime(&p_usertime, &p_systime);
      p_cputime  = p_usertime + p_systime;

      if ( cdoLogOff == 0 )
	{
	  cdologs(processNums(), p_cputime, processInqNvals()); 
	  cdolog(processInqPrompt(), p_cputime); 
	}
    }

  if ( cdoBenchmark )
    fprintf(stderr, " ( %.2fs %.2fs %.2fs %s)\n", c_usertime, c_systime, c_cputime, memstring);
  else
    {
      if ( ! cdoSilentMode )
	fprintf(stderr, " ( %.2fs %s)\n", c_cputime, memstring);
    }

  if ( cdoBenchmark && processID == 0 )
    fprintf(stderr, "total: user %.2fs  sys %.2fs  cpu %.2fs  mem%s\n",
	    p_usertime, p_systime, p_cputime, memstring);

  nstream = processInqStreamNum();
  for ( sindex = 0; sindex < nstream; sindex++ )
    {
      pstreamID = processInqStreamID(sindex);
      pstreamptr = pstream_to_pointer(pstreamID);
      if ( PSTREAM_Debug )
	Message(func, "process %d  stream %d  close streamID %d",
		processID, sindex, pstreamID);
      if ( pstreamptr ) pstreamClose(pstreamID);
    }

  processDelete();
}


int pstreamInqFiletype(int pstreamID)
{
  int filetype;
  PSTREAM *pstreamptr;

  pstreamptr = pstream_to_pointer(pstreamID);

#if  defined  (HAVE_LIBPTHREAD)
  if ( pstreamptr->ispipe )
    filetype = pstreamptr->filetype;
  else
#endif
    filetype = streamInqFiletype(pstreamptr->fileID);

  return (filetype);
}

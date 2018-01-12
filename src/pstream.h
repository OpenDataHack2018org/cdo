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

#ifndef PSTREAM_H
#define PSTREAM_H

#ifdef  HAVE_CONFIG_H
#include "config.h" /* _FILE_OFFSET_BITS influence off_t */
#endif

#include "pstream_write.h"
#include "varlist.h"
#include "argument.h"
#include "pipe.h"

#include <sys/types.h> /* off_t */
#include <vector>

class pstream_t
{
public:
  pstream_t(int id);
  int inqVlist();
  int inqFileType();
  void defTimestep(int p_tsID);
  bool isPipe();
  void pstreamOpenReadPipe(const char* pipename);
  void pstreamOpenReadFile(const char *argument);
  int pstreamOpenWriteFile(const char* p_filename, int filetype);
  int pstreamOpenWriteFile(int filetype);
  int pstreamOpenWritePipe(const char* filename, int filetype);
  void openAppend(const char * p_filename);
  void init();
  void defVlist(int p_vlistID);
  void close();
  int self; //aka the id of the pstream
  std::pair<int, int> m_id;
  int mode;
  int m_fileID;
  int m_vlistID;
  int tsID;
  int m_filetype;
  int tsID0;
  int mfiles;
  int nfiles;
  int varID; /* next varID defined with streamDefVar */
  bool ispipe;
  bool isopen;
  std::string m_name;
  std::vector<std::string> m_mfnames;
  varlist_t *m_varlist;
#ifdef  HAVE_LIBPTHREAD
  void *argument;
  std::shared_ptr<pipe_t> pipe;
  pthread_t rthreadID; /* read  thread ID */
  pthread_t wthreadID; /* write thread ID */
private:
   void createFilelist(const char *p_args);
   pstream_t();
   void defVarList(int vlistID);
#endif
};

void pstreamClose(int pstreamID);

int pstreamInqFiletype(int pstreamID);
int pstreamInqByteorder(int pstreamID);

int pstreamInqVlist(int pstreamID);

int pstreamInqTimestep(int pstreamID, int tsID);

int pstreamInqRecord(int pstreamID, int *varID, int *levelID);

void pstreamReadRecord(int pstreamID, double *data, size_t *nmiss);
void pstreamReadRecordF(int pstreamID, float *data, size_t *nmiss);
void pstreamCopyRecord(int pstreamIDdest, int pstreamIDsrc);

void pstreamInqGRIBinfo(int pstreamID, int *intnum, float *fltnum, off_t *bignum);

int pstreamFileID(int pstreamID);

void cdoVlistCopyFlag(int vlistID2, int vlistID1);

const int &getPthreadScope();
pstream_t *create_pstream();
pstream_t *create_pstream(std::vector<std::string> p_filenameList);
pstream_t *create_pstream(std::string p_filename);
pstream_t *create_pstream(int processID, int pstreamIDX);

int get_glob_argc();
void pstreamCloseAll();

#endif /* PSTREAM_H */

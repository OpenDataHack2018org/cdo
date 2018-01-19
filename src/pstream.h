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
#include "pipe.h"

#include <sys/types.h> /* off_t */
#include <vector>

class pstream_t
{
public:
  pstream_t(int id);
  int inqVlist();
  int inqFileType();
  int inqTimestep(int tsID);
  int inqRecord(int *varID, int *levelID);
  bool isPipe();
  int inqByteorder();
  size_t getNvals();

  int pstreamOpenReadPipe();
  int pstreamOpenWritePipe(const char* filename, int filetype);
  int pstreamOpenWriteFile(int filetype);
  void pstreamOpenReadFile(const char* filename);
  void openAppend(const char * p_filename);


  void readRecord(double *data, size_t *nmiss);
  void readRecordF(float *data, size_t *nmiss);
  void copyRecord(pstream_t * dest);

  void defVarList(int vlistID);
  void defTimestep(int p_tsID);
  void defVlist(int p_vlistID);
  void init();
  void close();
  void waitForPipe();
  void closePipe();

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
  std::shared_ptr<pipe_t> pipe;
  pthread_t rthreadID; /* read  thread ID */
  pthread_t wthreadID; /* write thread ID */
private:
   pstream_t();
#endif
};

pstream_t *pstream_to_pointer(int pstreamID);

int pstreamInqRecord(int pstreamID, int *varID, int *levelID);

void cdoVlistCopyFlag(int vlistID2, int vlistID1);

pstream_t *create_pstream();
pstream_t *create_pstream(std::vector<std::string> p_filenameList);
pstream_t *create_pstream(std::string p_filename);
pstream_t *create_pstream(int processID, int pstreamIDX);

int get_glob_argc();
void pstreamCloseAll();
void setProcessNum(int p_num);

#endif /* PSTREAM_H */

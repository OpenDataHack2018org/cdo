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

#ifndef PSTREAM_INT_H
#define PSTREAM_INT_H

#include "pstream.h"
/*clang-format off */

/**
 * Closes pstream and adds its variable count to its process.
 */
void pstreamClose(int pstreamID);

/*
 * Get the variable list
 *
 *    @param  streamID  Stream ID, from a previous call to
 * cdoStreamOpenRead(int) or cdoStreamOpenWrite(int,int).
 *
 * Description
 * he function pstreamInqVlist(int) returns the variable list of a stream.
 *
 * @return
 *  returns an identifier to the variable list.
 *
 */
int pstreamInqVlist(int pstreamID);
/**
 * Define the variable list.
 *
 *     @param  pstreamID Stream ID, from a previous call to
 * cdoStreamOpenRead(int) or cdoStreamOpenWrite(int,int).
 *     @param  vlistID  Variable list ID, from a previous call to
 * pstreamInqVlist(int).
 *
 * The function pstreamDefVlist(int, int) defines the variable list of a stream.
 *
 * To safeguard against errors by modifying the wrong vlist object,
 * this function makes the passed vlist object immutable.
 * All further vlist changes have to use the vlist object returned by
 * pstreamInqVlist().
 *
 */
void pstreamDefVlist(int pstreamID, int vlistID);

/**
 * Inquires record and sets \p varID.
 * @param pstreamID
 * @param varID
 * @param levelID
 * If \p pstreamID represents a pipe \p varID will be -1 if an error occured.
 * If \p pstreamID represents a file either varID will be set or CDO exits with
 * an error message.
 **/
void pstreamInqRecord(int pstreamID, int *varID, int *levelID);
/**
 *Define the next record.
 *The function pstreamDefRecord defines the meta-data of the next record.
 *    @param  pstreamID  Stream ID, from a previous call to cdoStreamOpenRead()
 *    or cdoStreamOpenWrite().
 *    @param  varID     Variable identifier.
 *    @param  levelID   Level identifier.
 */
void pstreamDefRecord(int pstreamID, int varID, int levelID);

int pstreamInqTimestep(PstreamType *p_pstreamptr, int tsID);
void pstreamDefTimestep(int pstreamID, int tsID);

void pstreamInqGRIBinfo(int pstreamID, int *intnum, float *fltnum,
                        off_t *bignum);
int pstreamInqFiletype(int pstreamID);
int pstreamInqByteorder(int pstreamID);
int pstreamFileID(int pstreamID);

void pstreamReadRecord(int pstreamID, double *data, size_t *nmiss);
void pstreamReadRecordF(int pstreamID, float *data, size_t *nmiss);
void pstreamCopyRecord(int pstreamIDdest, int pstreamIDsrc);

void pstreamWriteRecord(int pstreamID, double *data, size_t nmiss);
void pstreamWriteRecordF(int pstreamID, float *data, size_t nmiss);

/*clang-format on */

#endif

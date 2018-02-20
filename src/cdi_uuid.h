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
#ifndef CDI_UUID_H
#define CDI_UUID_H

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include "cdi.h"

#ifdef __cplusplus
extern "C"
{
#endif

  static inline int
  cdiUUIDIsNull(const unsigned char uuid[])
  {
    int isNull = 1;
    for (size_t i = 0; i < CDI_UUID_SIZE; ++i)
      isNull &= (uuid[i] == 0);
    return isNull;
  }

  void cdiCreateUUID(unsigned char uuid[CDI_UUID_SIZE]);

  void cdiUUID2Str(const unsigned char uuid[], char uuidstr[]);
  int cdiStr2UUID(const char *uuidstr, unsigned char uuid[]);

#if defined(__cplusplus)
}
#endif

#endif

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */

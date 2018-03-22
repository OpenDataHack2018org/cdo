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
#ifndef _DMEMORY_H
#define _DMEMORY_H

#include <stdio.h>

// if DEBUG_MEMORY is defined setenv MEMORY_DEBUG to debug memory
#define DEBUG_MEMORY

#ifndef WITH_FUNCTION_NAME
#define WITH_FUNCTION_NAME
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern size_t memTotal(void);
extern void memDebug(int debug);
extern void memExitOnError(void);

#if defined DEBUG_MEMORY

extern void *memRealloc(void *ptr, size_t size, const char *file, const char *functionname, int line);
extern void *memCalloc(size_t nmemb, size_t size, const char *file, const char *functionname, int line);
extern void *memMalloc(size_t size, const char *file, const char *functionname, int line);
extern void memFree(void *ptr, const char *file, const char *functionname, int line);

#if defined(__cplusplus)
}
#endif

#if defined WITH_FUNCTION_NAME
#define Realloc(p, s) memRealloc((p), (s), __FILE__, __func__, __LINE__)
#define Calloc(n, s) memCalloc((n), (s), __FILE__, __func__, __LINE__)
#define Malloc(s) memMalloc((s), __FILE__, __func__, __LINE__)
#define Free(p) memFree((p), __FILE__, __func__, __LINE__)
#else
#define Realloc(p, s) memRealloc((p), (s), __FILE__, (void *) NULL, __LINE__)
#define Calloc(n, s) memCalloc((n), (s), __FILE__, (void *) NULL, __LINE__)
#define Malloc(s) memMalloc((s), __FILE__, (void *) NULL, __LINE__)
#define Free(p) memFree((p), __FILE__, (void *) NULL, __LINE__)
#endif

#endif /* DEBUG_MEMORY */

#endif /* _DMEMORY_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */

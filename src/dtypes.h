#ifndef _DTYPES_H
#define _DTYPES_H

#include <stdio.h>
#include <limits.h>

/* INT32 */

#if ! defined (INT_MAX)
#  error INT_MAX undefined
#endif

#undef  INT32
#if  INT_MAX == 2147483647L
#  define  INT32  int
#elif LONG_MAX == 2147483647L
#  define  INT32  long
#endif

/* INT64 */

#if ! defined (LONG_MAX)
#  error LONG_MAX undefined
#endif

#undef  INT64
#if  LONG_MAX > 2147483647L
#  define  INT64  long
#else
#  define  INT64  long long
#endif

/* REAL4 */

#undef   REAL4
#define  REAL4  float

/* Real8 */

#undef   REAL8
#define  REAL8  double

#endif  /* _DTYPES_H */

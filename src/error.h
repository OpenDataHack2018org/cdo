#ifndef _ERROR_H
#define _ERROR_H

#ifndef  WITH_CALLER_NAME
#define  WITH_CALLER_NAME
#endif

#define  _FATAL     1     /* Error flag: exit on error  */
#define  _VERBOSE   2     /* Error flag: report errors  */
#define  _DEBUG     4     /* Error flag: debug          */

extern int _ExitOnError;  /* If set to 1, exit on error (default 1)       */
extern int _Verbose;      /* If set to 1, errors are reported (default 1) */
extern int _Debug;        /* If set to 1, debuggig (default 0)            */

void SysError(const char *caller, const char *fmt, ...);
void    Error(const char *caller, const char *fmt, ...);
void  Warning(const char *caller, const char *fmt, ...);
void  Message(const char *caller, const char *fmt, ...);

#if  defined  WITH_CALLER_NAME
#  define  SysError(...)  SysError(__func__, __VA_ARGS__)
#  define     Error(...)     Error(__func__, __VA_ARGS__)
#  define    Errorc(...)     Error(caller, __VA_ARGS__)
#  define   Warning(...)   Warning(__func__, __VA_ARGS__)
#  define   Message(...)   Message(__func__, __VA_ARGS__)
#  define  Messagec(...)   Message(caller, __VA_ARGS__)
#else
#  define  SysError(...)  SysError((void *), __VA_ARGS__)
#  define     Error(...)     Error((void *), __VA_ARGS__)
#  define    Errorc(...)     Error((void *), __VA_ARGS__)
#  define   Warning(...)   Warning((void *), __VA_ARGS__)
#  define   Message(...)   Message((void *), __VA_ARGS__)
#  define  Messagec(...)   Message((void *), __VA_ARGS__)
#endif

#endif  /* _ERROR_H */

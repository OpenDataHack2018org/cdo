#ifndef EXCEPTION_H
#define EXCEPTION_H

void    cdiOpenError(int cdiErrno, const char *fmt, const char *path);
void    cdoAbort(const char *fmt, ...);
void    cdoWarning(const char *fmt, ...);
void    cdoPrint(const char *fmt, ...);
void    cdoPrintBlue(const char *fmt, ...);
void    cdoPrintRed(const char *fmt, ...);

#endif

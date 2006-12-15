#if  defined  (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <pwd.h>
#include <unistd.h>  /* write, close */
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>

#include "cdo.h"
#include "dtypes.h"

#if ! defined (VERSION)
#  define  VERSION  "0.0.1"
#endif

#define  MAX_LEN  65536

void cdolog(const char *prompt, double cputime)
{
#if defined (LOGPATH)
#define  XSTRING(x)	#x
#define  STRING(x)	XSTRING(x)
  char logfilename[] = STRING(LOGPATH) "/cdo.log";
  int  logfileno;
  char *username;
  char timestr[30];
  time_t date_and_time_in_sec;
  struct tm *date_and_time;
  int streamID;
  const char *streamName;
  int len, slen, olen, pos;
  int loper;
  char logstring[MAX_LEN];
  int status;
  struct flock mylock;

  memset(logstring, 0, MAX_LEN);

  date_and_time_in_sec = time(NULL);

  if ( date_and_time_in_sec != -1 )
    {
      date_and_time = localtime(&date_and_time_in_sec);
      (void) strftime(timestr, sizeof(timestr), "%d/%m/%Y %H:%M", date_and_time);
    }

  username = getenv("LOGNAME");
  if ( username == NULL )
    {
      username = getenv("USER");
      if ( username == NULL ) username = "unknown";
    }

  slen = sprintf(logstring, "%s %-8s %7.2f %s %-3s",
		 timestr,  username, cputime, VERSION, prompt);

  for ( streamID = 0; streamID < cdoStreamCnt(); streamID++ )
    {
      streamName = cdoStreamName(streamID);
      pos = 0;
      while ( pos < (int) strlen(streamName) )
	{
	  len = 0;
	  loper = 0;
	  if ( streamName[pos++] == '-' ) loper = 1;
	  while ( streamName[pos+len] != ' ' &&  streamName[pos+len] != '\0' ) len++;
	  if ( len && loper )
	    {
	      for ( olen = 1; olen < len; olen++ )
		if ( streamName[pos+olen] == ',' ) break;

	      sprintf(logstring+slen, " %.*s", olen, &streamName[pos]);
	      slen += olen + 1;
	    }
	  pos += len + 1;
	}
    }

  sprintf(logstring+slen, "\n");
  slen++;

  errno = 0;
  logfileno = open(logfilename, O_WRONLY | O_APPEND);
  if ( errno )
    {
      errno = 0;
      return;
    }

  mylock.l_type   = F_WRLCK;
  mylock.l_whence = SEEK_SET;
  mylock.l_start  = 0;
  mylock.l_len    = 0;

  status = fcntl(logfileno, F_SETLKW, &mylock);
  if ( status == 0 )
    {
      write(logfileno, logstring, slen);

      mylock.l_type   = F_UNLCK;
      status = fcntl(logfileno, F_SETLKW, &mylock);
    }

  close(logfileno);

  errno = 0;

  return;

#endif
}


#include <math.h>
/*
 * convert an IMB float to single precision number v1.0
 *
 *                      Wesley Ebisuzaki
 */

static float ibm2flt(unsigned char *ibm) {

	int positive, power;
	unsigned int abspower;
	long int mant;
	double value, exp;

	positive = (ibm[0] & 0x80) == 0;
	mant = (ibm[1] << 16) + (ibm[2] << 8) + ibm[3];
	power = (int) (ibm[0] & 0x7f) - 64;
	abspower = power > 0 ? power : -power;


	/* calc exp */
	exp = 16.0;
	value = 1.0;
	while (abspower) {
		if (abspower & 1) {
			value *= exp;
		}
		exp = exp * exp;
		abspower >>= 1;
	}

	if (power < 0) value = 1.0 / value;
	value = value * mant / 16777216.0;
	if (positive == 0) value = -value;
	return (float)value;
}

/*
 * convert a float to an IBM single precision number v1.0
 *
 *                      Wesley Ebisuzaki
 *
 * doesn't handle subnormal numbers
 */

static int flt2ibm(float x, unsigned char *ibm) {

	int sign, exp, i;
	double mant;

	if ( !(fabs((double)x) > 0) ) {
		ibm[0] = ibm[1] = ibm[2] = ibm[3] = 0;
		return 0;
	}

	/* sign bit */
	if (x < 0.0) {
		sign = 128;
		x = -x;
	}
	else sign = 0;

	mant = frexp((double) x, &exp);

	/* round up by adding 2**-24 */
	/* mant = mant + 1.0/16777216.0; */

	if (mant >= 1.0) {
		mant = 0.5;
		exp++;
	}
	while (exp & 3) {
		mant *= 0.5;
		exp++;
	}
	
	exp = exp/4 + 64;

	if (exp < 0) {
		fprintf(stderr,"underflow in flt2ibm\n");
		ibm[0] = ibm[1] = ibm[2] = ibm[3] = 0;
		return 0;
	}
	if (exp > 127) {
		fprintf(stderr,"overflow in flt2ibm\n");
		ibm[0] = sign | 127;
		ibm[1] = ibm[2] = ibm[3] = 255;
		return -1;
	}

	/* normal number */

	ibm[0] = sign | exp;

	mant = mant * 256.0;
	i = (int) floor(mant);
	mant = mant - i;
	ibm[1] = i;

	mant = mant * 256.0;
	i = (int) floor(mant);
	mant = mant - i;
	ibm[2] = i;

	ibm[3] = (int) floor(mant*256.0);

	return 0;
}


#define  GET_UINT4(xb)        ((int) (((int)xb[0]<<24) + \
                                      ((int)xb[1]<<16) + \
                                      ((int)xb[2]<<8)  + \
                                      ((int)xb[3])))
#define  GET_UINT8(xb)        ((INT64) (((INT64)xb[0]<<56) + \
                                        ((INT64)xb[1]<<48) + \
                                        ((INT64)xb[2]<<40) + \
                                        ((INT64)xb[3]<<32) + \
                                        ((INT64)xb[4]<<24) + \
                                        ((INT64)xb[5]<<16) + \
                                        ((INT64)xb[6]<<8)  + \
					((INT64)xb[7])))

#define  PUT_UINT4(xb, iv)    ((*(xb)   = (iv) >> 24), \
                               (*(xb+1) = (iv) >> 16), \
                               (*(xb+2) = (iv) >>  8), \
                               (*(xb+3) = (iv)))
#define  PUT_UINT8(xb, iv)    ((*(xb)   = (iv) >> 56), \
                               (*(xb+1) = (iv) >> 48), \
                               (*(xb+2) = (iv) >> 40), \
                               (*(xb+3) = (iv) >> 32), \
                               (*(xb+4) = (iv) >> 24), \
                               (*(xb+5) = (iv) >> 16), \
                               (*(xb+6) = (iv) >>  8), \
                               (*(xb+7) = (iv)))

void cdologs(int noper, double cputime, off_t nvals)
{
#if defined (LOGPATH)
#define  LOGSIZE  32
#define  XSTRING(x)	#x
#define  STRING(x)	XSTRING(x)
  char logfilename[] = STRING(LOGPATH) "/cdo.logs";
  int  logfileno;
  char timestr[30];
  time_t date_and_time_in_sec;
  struct tm *date_and_time;
  int status;
  int date = 0, ncdo = 0, nhours = 0;
  int date0 = 0, ncdo0, noper0, nhours0;
  double cputime0;
  INT64 nvals0;
  unsigned char logbuf[LOGSIZE];
  unsigned char *logdate   =  logbuf;
  unsigned char *logncdo   = &logbuf[4];
  unsigned char *lognoper  = &logbuf[8];
  unsigned char *logctime  = &logbuf[12];
  unsigned char *lognvals  = &logbuf[16];
  unsigned char *lognhours = &logbuf[24];
  const size_t logsize = LOGSIZE;
  size_t bufsize;
  struct flock mylock;
  struct stat filestat;

  memset(logbuf, 0, logsize);

  date_and_time_in_sec = time(NULL);

  if ( date_and_time_in_sec != -1 )
    {
      date_and_time = localtime(&date_and_time_in_sec);
      (void) strftime(timestr, sizeof(timestr), "%Y%m%d", date_and_time);
      date = atoi(timestr);
    }

  errno = 0;
  logfileno = open(logfilename, O_RDWR);
  if ( errno )
    {
      errno = 0;
      return;
    }

  mylock.l_type   = F_WRLCK;
  mylock.l_whence = SEEK_SET;
  mylock.l_start  = 0;
  mylock.l_len    = 0;

  status = fcntl(logfileno, F_SETLKW, &mylock);
  errno = 0;
  if ( status != 0 ) goto endlabel;

  status = fstat(logfileno, &filestat);
  errno = 0;
  if ( status != 0 ) goto endlabel;

  bufsize = (size_t) filestat.st_size;

  if ( bufsize > 0 )
    {
      status = (int) lseek(logfileno, (off_t) (bufsize-logsize), SEEK_SET);
      status = (int) read(logfileno, logbuf, logsize);

      date0    = GET_UINT4(logdate);
      ncdo0    = GET_UINT4(logncdo);
      noper0   = GET_UINT4(lognoper);
      cputime0 = (double) ibm2flt(logctime);
      nvals0   = GET_UINT8(lognvals);
      nhours0  = GET_UINT4(lognhours);

      if ( date == date0 )
	{
	  ncdo = ncdo0;
	  nhours = nhours0;
	  noper += noper0;
	  cputime += cputime0;
	  nvals += (off_t) nvals0;
	  status = (int) lseek(logfileno, (off_t) (bufsize-logsize), SEEK_SET);
	}
    }

  if ( date >= date0 )
    {
      ncdo++;

      while ( cputime > 3600 ) { cputime -= 3600; nhours++; }

      PUT_UINT4(logdate, date);
      PUT_UINT4(logncdo, ncdo);
      PUT_UINT4(lognoper, noper);
      flt2ibm((float)cputime, logctime);
      PUT_UINT8(lognvals, nvals);
      PUT_UINT4(lognhours, nhours);
      /*
      mylock.l_type   = F_WRLCK;
      status = fcntl(logfileno, F_SETLKW, &mylock);
      */
      status = (int) write(logfileno, logbuf, logsize);
    }

 endlabel:

  mylock.l_type   = F_UNLCK;
  status = fcntl(logfileno, F_SETLK, &mylock);

  close(logfileno);

  errno = 0;

  return;

#endif
}


void dumplogs(const char *logfilename)
{
#define  LOGSIZE  32
  static const char func[] = "dumplogs";
  int  logfileno;
  int status;
  int date0 = 0, ncdo0, noper0, nhours0;
  int nlogs;
  int i;
  double cputime0;
  INT64 nvals0;
  unsigned char logbuf[LOGSIZE];
  unsigned char *logdate   =  logbuf;
  unsigned char *logncdo   = &logbuf[4];
  unsigned char *lognoper  = &logbuf[8];
  unsigned char *logctime  = &logbuf[12];
  unsigned char *lognvals  = &logbuf[16];
  unsigned char *lognhours = &logbuf[24];
  unsigned char *buffer = NULL;
  const size_t logsize = LOGSIZE;
  size_t bufsize;
  struct stat filestat;

  errno = 0;
  logfileno = open(logfilename, O_RDONLY);
  if ( errno )
    {
      cdoAbort("Open failed on %s", logfilename);
      errno = 0;
      return;
    }

  status = fstat(logfileno, &filestat);
  errno = 0;
  if ( status != 0 ) return;

  bufsize = (size_t) filestat.st_size;

  if ( bufsize > 0 )
    {
      buffer = (unsigned char *) malloc(bufsize);

      status = (int) read(logfileno, buffer, bufsize);

      nlogs = bufsize / logsize;
      for ( i = 0; i < nlogs; i++ )
	{
	  memcpy(logbuf, &buffer[i*logsize], logsize);
	  date0    = GET_UINT4(logdate);
	  ncdo0    = GET_UINT4(logncdo);
	  noper0   = GET_UINT4(lognoper);
	  cputime0 = (double) ibm2flt(logctime);
	  nvals0   = GET_UINT8(lognvals);
	  nhours0  = GET_UINT4(lognhours);

	  if ( sizeof(INT64) > sizeof(long) )
	    fprintf(stdout, "%8d %10d %10d %19lld %10d %8.2f\n", date0, ncdo0, noper0, (long long)nvals0, nhours0, cputime0);
	  else
	    fprintf(stdout, "%8d %10d %10d %19ld %10d %8.2f\n", date0, ncdo0, noper0, (long)nvals0, nhours0, cputime0);
	}

      free(buffer);
    }

  close(logfileno);

  errno = 0;

  return;
}

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

#if ! defined (VERSION)
#  define  VERSION  "0.0.1"
#endif

#define  MAX_LEN  65536

void userlog(const char *prompt, double cputime)
{
#if defined (LOGPATH)
#define  XSTRING(x)	#x
#define  STRING(x)	XSTRING(x)
  char logfilename[] = STRING(LOGPATH) "/cdo.log";
  int   logfileno;
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
      (void) strftime(timestr, 30, "%d/%m/%Y %H:%M", date_and_time);
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

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include <search.h>
#include <ctype.h>

#if defined(HAVE_LIBCMOR)
#include "cmor.h"
#endif

static char *trim(char *s)
{
  int n;
  if (s == NULL)
    return s;
  while ( *s != '\0' && (isspace(*s) || *s == '"') )
    s++;
  n = strlen(s);
  while ( n > 0 && (isspace(s[n - 1]) || s[n - 1] == '"') )
    n--;
  s[n] = '\0';
  return s;
}

static int hinsert(char *kvstr)
{
  char *key, *value;
  ENTRY e;

  key = trim(strtok(kvstr, "="));
  value = trim(strtok(NULL, "="));
  if ( key == NULL || value == NULL )
    return 1;
  e.key = strdup(key);
  e.data = (void *)strdup(value);
  hsearch(e, ENTER);
  return 0;
}

static int parse_kvfile(const char *filename)
{
  FILE *fp;
  char line[1024], *comment;

  hcreate(100);
  fp = fopen(filename, "r");
  if ( fp == NULL )
    return 1;
  while ( fgets(line, sizeof(line), fp) != NULL )
    {
      comment = strchr(line, '#');
      if ( comment )
        *comment = '\0';
      hinsert(line);
    }
  fclose(fp);
  return 0;
}

void *CMOR(void *argument)
{
  cdoInitialize(argument);

#if defined(HAVE_LIBCMOR)
  int nparams = operatorArgc();
  char **params = operatorArgv();
  char tab[CMOR_MAX_STRING];
  char *option, *value;
  char *var, *chunk, *expinfo, *modinfo, *userinfo;
  if ( cdoVerbose )
    for ( int i = 0; i < nparams; ++i )
      printf("param %d: %s\n", i, params[i]);

  if ( nparams < 1 ) cdoAbort("Too few arguments!");
  strncpy(tab, params[0], CMOR_MAX_STRING);
  for ( int i = 1; i < nparams; ++i )
    {
      option = strtok(params[i], "=");
      value = strtok(NULL, "=");
      if ( value == NULL )
        cdoAbort("Missing value!");
      if ( strncasecmp(option, "var", 3) == 0 )
        var = value;
      else if ( strncasecmp(option, "chunk", 5) == 0 )
        chunk = value;
      else if ( strncasecmp(option, "expinfo", 7) == 0 )
        expinfo = value;
      else if ( strncasecmp(option, "modinfo", 7) == 0 )
        modinfo = value;
      else if ( strncasecmp(option, "userinfo", 8) == 0 )
        userinfo = value;
      else
        cdoAbort("Unknown option!");
    }

  parse_kvfile("cmor.rc");

  /*
  if ( cdoVerbose )
    cdoPrint("Using CMOR version %d.%d.%d.", CMOR_VERSION_MAJOR, CMOR_VERSION_MINOR, CMOR_VERSION_PATCH);
  */
  int streamID = streamOpenRead(cdoStreamName(0));

  int vlistID = streamInqVlist(streamID);
  int taxisID = vlistInqTaxis(vlistID);

  int gridsize = vlistGridsizeMax(vlistID);
  double *array = (double*) Malloc(gridsize*sizeof(double));

  int recID, varID, levelID, nrecs, nmiss;
  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID, tsID)) )
    {
      for ( recID = 0; recID < nrecs; recID++ )
        {
          streamInqRecord(streamID, &varID, &levelID);
          streamReadRecord(streamID, array, &nmiss);
        }

      tsID++;
    }

  streamClose(streamID);

  if ( array ) Free(array);
#else
  cdoWarning("CMOR support not compiled in!");
#endif

  hdestroy();
  cdoFinish();

  return 0;
}
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */

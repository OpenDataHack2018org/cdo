#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include <search.h>
#include <ctype.h>

#if defined(HAVE_LIBCMOR)
#include "cmor.h"

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

static ENTRY *hinsert(char *kvstr)
{
  char *key, *value;
  ENTRY e, *ep;

  key = trim(strtok(kvstr, "="));
  value = trim(strtok(NULL, "="));
  if ( key == NULL || value == NULL )
    return NULL;
  e.key = strdup(key);
  e.data = (void *)strdup(value);
  ep = hsearch(e, ENTER);
  return ep;
}

static int parse_kvfile(const char *filename)
{
  FILE *fp;
  char line[1024], *comment;

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

static char *get_val(char *key)
{
  ENTRY e, *ep;

  e.key = key;
  ep = hsearch(e, FIND);
  if ( ep )
    return (char *)ep->data;
  else
    return NULL;
}

static void parse_cmdline_kv(int nparams, char **params)
{
  int i, j, k, size;
  char *p;

  /* Assume key = value pairs. That is, if params[i] contains no '='
   * then treat it as if it belongs to the value of params[i-1],
   * separated by a ','.*/
  i = 0;
  while ( i < nparams )
    {
      j = 1;
      size = strlen(params[i]) + 1;
      while ( i + j < nparams && strchr(params[i + j], '=') == NULL )
        {
          size += strlen(params[i + j]) + 1;
          j++;
        }
      p = (char *) Malloc(size);
      strcpy(p, params[i]);
      for (k = 1; k < j; k++)
        {
          strcat(p, ",");
          strcat(p, params[i + k]);
        }
      hinsert(p);
      free(p);
      i += j;
    }
}
#endif

void *CMOR(void *argument)
{
  cdoInitialize(argument);

#if defined(HAVE_LIBCMOR)
  int nparams = operatorArgc();
  char **params = operatorArgv();
  char *table, *vars;
  char *var_list[CMOR_MAX_VARIABLES];
  int use_n_vars = 0;

  if ( nparams < 1 ) cdoAbort("Too few arguments!");
  hcreate(100);
  parse_kvfile("cmor.rc");
  table = params[0];
  parse_cmdline_kv(nparams - 1, &params[1]);

  vars = get_val("var");
  if ( vars )
    {
      var_list[0] = strtok(strdup(vars), ",");
      do
        use_n_vars++;
      while ( (var_list[use_n_vars] = strtok(NULL, ",")) );
    }

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
  hdestroy();

  if ( array ) Free(array);
#else
  cdoWarning("CMOR support not compiled in!");
#endif

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

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include <search.h>
#include <ctype.h>

#if defined(HAVE_LIBCMOR)
#include "cmor.h"

typedef struct _cmor_var
{
  int cdi_varID;
  int cmor_varID;
  char datatype;
} cmor_var;

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
  ep = hsearch(e, FIND);
  if ( ep )
    ep->data = e.data;
  else
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

static char *get_val(char *key, char *def)
{
  ENTRY e, *ep;

  e.key = key;
  ep = hsearch(e, FIND);
  if ( ep )
    return (char *)ep->data;
  else
    return def;
}

static void append_attr(int vlistID, int varID,
                        char *key, char *value, size_t n)
{
  memset(value, 0, n);
  int status = vlistInqAttTxt(vlistID, varID, key, n - 1, value);
  char *hval = get_val(key, NULL);

  if ( hval != NULL )
    {
      if ( status == CDI_NOERR )
        {
          n -= strlen(value) + 2;
          value += strlen(value);
          if ( n > 0 )
            {
              *value++ = ' ';
              strncpy(value, hval, n);
            }
          }
      else
        {
          strncpy(value, hval, n - 1);
        }
    }
}

static char *substitute(char *word)
{
  ENTRY e, *ep;
  char *key;

  key = (char *) Malloc (strlen(word) + 12);
  sprintf(key, "substitute_%s", word);
  e.key = key;
  ep = hsearch(e, FIND);
  Free(key);
  if ( ep )
    return (char *)ep->data;
  else
    return word;
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
  char *table;
  char *chunk;
  char *logfile;
  int netcdf_file_action, exit_control;
  int set_verbosity;
  int create_subdirectories;

  if ( nparams < 1 ) cdoAbort("Too few arguments!");
  hcreate(100);
  parse_kvfile("cmor.rc");
  table = params[0];
  nparams--;
  params++;
  parse_cmdline_kv(nparams, params);

  chunk = get_val("chunk", "replace");
  if ( strcasecmp(chunk, "replace") == 0 )
    netcdf_file_action = CMOR_REPLACE;
  else if ( strcasecmp(chunk, "append") == 0 )
    netcdf_file_action = CMOR_APPEND;

  set_verbosity = CMOR_NORMAL;
  if ( strcasecmp(get_val("set_verbosity", ""), "CMOR_QUIET") == 0 )
    set_verbosity = CMOR_QUIET;

  exit_control = CMOR_NORMAL;
  if ( strcasecmp(get_val("exit_control", ""), "CMOR_EXIT_ON_MAJOR") == 0 )
    exit_control = CMOR_EXIT_ON_MAJOR;
  if ( strcasecmp(get_val("exit_control", ""), "CMOR_EXIT_ON_WARNING") == 0 )
    exit_control = CMOR_EXIT_ON_WARNING;

  logfile = get_val("logfile", NULL);

  create_subdirectories = atoi(get_val("create_subdirectories", "0"));
  cmor_setup(get_val("inpath", "/usr/share/cmor/"),
             &netcdf_file_action,
             &set_verbosity,
             &exit_control,
             logfile,
             &create_subdirectories);

  int *month_lengths;
  int table_id;
  char *calendar;
  int streamID = streamOpenRead(cdoStreamName(0));
  int vlistID = streamInqVlist(streamID);
  int taxisID = vlistInqTaxis(vlistID);

  switch(taxisInqCalendar(taxisID))
    {
    case CALENDAR_STANDARD:
      calendar = "gregorian";
      break;
    case CALENDAR_PROLEPTIC:
      calendar = "proleptic_gregorian";
      break;
    case CALENDAR_360DAYS:
      calendar = "360_day";
      break;
    case CALENDAR_365DAYS:
      calendar = "noleap";
      break;
    case CALENDAR_366DAYS:
      calendar = "all_leap";
      break;
    default:
      cdoAbort("Unsupported calendar type.");
    }

  if ( get_val("month_lengths", NULL) )
    {
      char *month_lengths_str = strdup(get_val("month_lengths", ""));
      char *month_str = strtok(month_lengths_str, ",");
      int month = 0;
      month_lengths = Malloc (12 * sizeof(int));
      while ( month < 12 && month_str != NULL )
        {
          month_lengths[month++] = atoi(month_str);
          month_str = strtok(NULL, ",");
        }
      if ( month != 12 )
        cdoAbort("Invalid format for month_lengths");
    }
  else
    {
      month_lengths = NULL;
    }

  double branch_time = atof(get_val("branch_time", "0.0"));
  char history[CMOR_MAX_STRING];
  append_attr(vlistID, CDI_GLOBAL, "history", history, CMOR_MAX_STRING);
  cmor_dataset(get_val("outpath", "./"),
               get_val("expinfo", ""),
               get_val("institution", ""),
               get_val("modinfo", ""),
               calendar,
               atoi(get_val("realization", "1")),
               get_val("contact", ""),
               history,
               get_val("comment", ""),
               get_val("references", ""),
               atoi(get_val("leap_year", "0")),
               atoi(get_val("leap_month", "0")),
               month_lengths,
               get_val("model_id", ""),
               get_val("forcing", ""),
               atoi(get_val("initialization_method", "1")),
               atoi(get_val("physics_version", "1")),
               get_val("institute_id", ""),
               get_val("parent_experiment_id", ""),
               &branch_time,
               get_val("parent_experiment_rip", ""));

  cmor_load_table(table, &table_id);
  cmor_set_table(table_id);

  int nvars = vlistNvars(vlistID);
  cmor_var *cmor_var_list = (cmor_var *) Malloc(nvars * sizeof(cmor_var));;
  cmor_var *cmor_var;
  int copy_nvars = 0;
  int recID, varID, levelID, nrecs, nmiss;
  int gridID;
  char name[CDI_MAX_NAME], units[CDI_MAX_NAME];
  int length;
  double *coord_vals, *cell_bounds;
  int ndims;
  double missing_value;
  double tolerance = 1e-4;
  int maxlevels = 0;
  int axis_ids[CMOR_MAX_AXES];
  char *select_vars = get_val("var", NULL);

  for ( varID = 0; varID < nvars; varID++ )
    {
      vlistInqVarName(vlistID, varID, name);
      if ( select_vars == NULL || strstr(select_vars, name) != NULL)
        {
          cmor_var = &cmor_var_list[copy_nvars++];
          cmor_var->cdi_varID = varID;
          gridID = vlistInqVarGrid(vlistID, varID);
          ndims = 0;

          /* Time-Axis */
          int rdate = taxisInqRdate(taxisID);
          int rtime = taxisInqRtime(taxisID);
          int timeunit = taxisInqTunit(taxisID);
          int year, month, day, hour, minute, second;
          cdiDecodeDate(rdate, &year, &month, &day);
          cdiDecodeTime(rtime, &hour, &minute, &second);

          if ( timeunit == TUNIT_QUARTER ||
               timeunit == TUNIT_30MINUTES )
            timeunit = TUNIT_MINUTE;
          if ( timeunit == TUNIT_3HOURS ||
               timeunit == TUNIT_6HOURS ||
               timeunit == TUNIT_12HOURS )
            timeunit = TUNIT_HOUR;

          sprintf(units, "%s since %d-%d-%d %02d:%02d:%02d",
                  tunitNamePtr(timeunit), year, month, day, hour,
                  minute, second);

          cmor_axis(&axis_ids[ndims++],
                    "time",
                    units,
                    0,
                    NULL,
                    0,
                    NULL,
                    0,
                    NULL);

          /* Z-Axis */
          int zaxisID = vlistInqVarZaxis(vlistID, varID);
          zaxisInqName(zaxisID, name);
          zaxisInqUnits(zaxisID, units);
          length = zaxisInqSize(zaxisID);
          if ( length > maxlevels )
            maxlevels = length;
          coord_vals = Malloc(length * sizeof(double));
          zaxisInqLevels(zaxisID, coord_vals);
          cmor_axis(&axis_ids[ndims++],
                    substitute(name),
                    units,
                    length,
                    (void *)coord_vals,
                    'd',
                    NULL,
                    0,
                    NULL);

          /* Y-Axis */
          gridInqYname(gridID, name);
          gridInqYunits(gridID, units);
          length = gridInqYsize(gridID);
          coord_vals = Malloc(length * sizeof(double));
          gridInqYvals(gridID, coord_vals);
          cell_bounds = Malloc(2 * length * sizeof(double));
          gridInqYbounds(gridID, cell_bounds);
          cmor_axis(&axis_ids[ndims++],
                    substitute(name),
                    units,
                    length,
                    (void *)coord_vals,
                    'd',
                    (void *)cell_bounds,
                    2,
                    NULL);

          /* X-Axis */
          gridInqXname(gridID, name);
          gridInqXunits(gridID, units);
          length = gridInqXsize(gridID);
          coord_vals = Malloc(length * sizeof(double));
          gridInqXvals(gridID, coord_vals);
          cell_bounds = Malloc(2 * length * sizeof(double));
          gridInqXbounds(gridID, cell_bounds);
          cmor_axis(&axis_ids[ndims++],
                    substitute(name),
                    units,
                    length,
                    (void *)coord_vals,
                    'd',
                    (void *)cell_bounds,
                    2,
                    NULL);

          /* Variable */
          vlistInqVarUnits(vlistID, varID, units);
          missing_value = vlistInqVarMissval(vlistID, varID);
          if ( vlistInqVarDatatype(vlistID, varID) == DATATYPE_FLT32 )
            cmor_var->datatype = 'f';
          else
            cmor_var->datatype = 'd';
          vlistInqVarName(vlistID, varID, name);
          cmor_variable(&cmor_var->cmor_varID,
                        substitute(name),
                        units,
                        ndims,
                        axis_ids,
                        cmor_var->datatype,
                        &missing_value,
                        &tolerance,
                        NULL, // positive,
                        NULL,
                        NULL,
                        NULL);
        }
      else
        {
          printf("Not found var %s\n", name);
        }
    }

  int gridsize = vlistGridsizeMax(vlistID);
  double *array = (double*) Malloc(gridsize * maxlevels * sizeof(double));
  double time_vals;
  double time_bnds[2];
  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID, tsID)) )
    {
      for ( int i = 0; i < copy_nvars; i++ )
        {
          cmor_var = &cmor_var_list[i];
          time_vals = 0.5 + tsID;
          time_bnds[0] = tsID;
          time_bnds[1] = tsID + 1.0;
          for ( recID = 0; recID < nrecs; recID++ )
            {
              streamInqRecord(streamID, &varID, &levelID);
              if ( cmor_var->datatype == 'f' )
                streamReadRecordF(streamID,
                                  (float *)array + gridsize * levelID, &nmiss);
              else
                streamReadRecord(streamID, array + gridsize * levelID, &nmiss);
            }
          //streamReadVarF(streamID, cmor_var->cdi_varID, array, &nmiss);
          cmor_write(cmor_var->cmor_varID,
                     (void *)array,
                     cmor_var->datatype,
                     NULL, // file_suffix
                     1,    // ntimes_passed
                     &time_vals, // time_vals
                     time_bnds, // time_bnds
                     NULL);// store_with
        }
      tsID++;
    }
  Free (array);
  streamClose(streamID);
  cmor_close();
  hdestroy();

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

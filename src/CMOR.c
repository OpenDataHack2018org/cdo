#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include <search.h>
#include <ctype.h>

#if defined(HAVE_LIBCMOR)
#include "cmor.h"

typedef struct _cc_var_t
{
  int cdi_varID;
  int cmor_varID;
  char datatype;
  void *data;
} cc_var_t;

static cc_var_t *find_var(const int cdi_varID, cc_var_t vars[], int nvars)
{
  for ( int i = 0; i < nvars; i++ )
    if ( cdi_varID == vars[i].cdi_varID )
      return &vars[i];
  return NULL;
}

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

static void setup(int streamID, char *table)
{
  char *chunk;
  char *logfile;
  int netcdf_file_action, exit_control;
  int set_verbosity;
  int create_subdirectories;
  int *month_lengths;
  int table_id;
  int taxisID = vlistInqTaxis(streamInqVlist(streamID));
  char *calendar;
  double branch_time = atof(get_val("branch_time", "0.0"));

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

  switch ( taxisInqCalendar(taxisID) )
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

  cmor_dataset(get_val("outpath", "./"),
               get_val("expinfo", ""),
               get_val("institution", ""),
               get_val("modinfo", ""),
               calendar,
               atoi(get_val("realization", "1")),
               get_val("contact", ""),
               get_val("history", ""),
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
}

static void define_variables(int streamID, cc_var_t vars[], int *nvars)
{
  int vlistID = streamInqVlist(streamID);
  int taxisID = vlistInqTaxis(vlistID);
  size_t gridsize = vlistGridsizeMax(vlistID);
  cc_var_t *var;
  int varID, gridID;
  char name[CDI_MAX_NAME], units[CDI_MAX_NAME];
  int length;
  double *coord_vals, *cell_bounds;
  int ndims, levels;
  char missing_value[sizeof(double)];
  double tolerance = 1e-4;
  int axis_ids[CMOR_MAX_AXES];
  char *select_vars = get_val("var", NULL);
  int year, month, day, hour, minute, second;
  int timeunit = taxisInqTunit(taxisID);
  char taxis_units[CMOR_MAX_STRING];

  cdiDecodeDate(taxisInqRdate(taxisID), &year, &month, &day);
  cdiDecodeTime(taxisInqRtime(taxisID), &hour, &minute, &second);
  if ( timeunit == TUNIT_QUARTER || timeunit == TUNIT_30MINUTES )
    timeunit = TUNIT_MINUTE;
  if ( timeunit == TUNIT_3HOURS ||
       timeunit == TUNIT_6HOURS ||
       timeunit == TUNIT_12HOURS )
    timeunit = TUNIT_HOUR;

  sprintf(taxis_units, "%s since %d-%d-%d %02d:%02d:%02d",
          tunitNamePtr(timeunit), year, month, day, hour,
          minute, second);

  *nvars = 0;
  for ( varID = 0; varID < vlistNvars(vlistID); varID++ )
    {
      vlistInqVarName(vlistID, varID, name);
      if ( select_vars == NULL || strstr(select_vars, name) != NULL)
        {
          var = &vars[(*nvars)++];
          var->cdi_varID = varID;
          gridID = vlistInqVarGrid(vlistID, varID);
          ndims = 0;

          /* Time-Axis */
          cmor_axis(&axis_ids[ndims++],
                    substitute("time"),
                    taxis_units,
                    0,
                    NULL,
                    0,
                    NULL,
                    0,
                    NULL);

          /* Z-Axis */
          int zaxisID = vlistInqVarZaxis(vlistID, varID);
          levels = zaxisInqSize(zaxisID);
          if ( zaxisInqType(zaxisID) != ZAXIS_SURFACE )
            {
              coord_vals = Malloc(levels * sizeof(double));
              zaxisInqLevels(zaxisID, coord_vals);
              zaxisInqName(zaxisID, name);
              zaxisInqUnits(zaxisID, units);
              cmor_axis(&axis_ids[ndims++],
                        substitute(name),
                        units,
                        levels,
                        (void *)coord_vals,
                        'd',
                        NULL,
                        0,
                        NULL);
            }

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
          vlistInqVarName(vlistID, varID, name);
          if ( vlistInqVarDatatype(vlistID, varID) == DATATYPE_FLT32 )
            {
              var->datatype = 'f';
              *(float *) missing_value = vlistInqVarMissval(vlistID, varID);
              var->data = Malloc(gridsize * levels * sizeof(float));
            }
          else
            {
              var->datatype = 'd';
              *(double *) missing_value = vlistInqVarMissval(vlistID, varID);
              var->data = Malloc(gridsize * levels * sizeof(double));
            }
          cmor_variable(&var->cmor_varID,
                        substitute(name),
                        units,
                        ndims,
                        axis_ids,
                        var->datatype,
                        (void *) missing_value,
                        &tolerance,
                        NULL,
                        NULL,
                        NULL,
                        NULL);
        }
    }
}

static void write_variables(int streamID, cc_var_t vars[], int nvars)
{
  cc_var_t *var;
  int vlistID = streamInqVlist(streamID);
  int taxisID = vlistInqTaxis(vlistID);
  size_t gridsize = vlistGridsizeMax(vlistID);
  double time_val;
  double time_bnds[2];
  double *time_bndsp;
  int has_bnds = taxisHasBounds(taxisID);
  int tsID;
  int vdate0b, vdate1b;
  int vtime0b, vtime1b;
  juldate_t juldate, r_juldate;
  int calendar = taxisInqCalendar(taxisID);
  int tunitsec;
  int nrecs;
  int varID, levelID;
  int nmiss;

  switch ( taxisInqTunit(taxisID) )
    {
    case TUNIT_MINUTE: tunitsec = 60; break;
    case TUNIT_HOUR: tunitsec = 3600; break;
    case TUNIT_DAY: tunitsec = 86400; break;
    default: tunitsec = 3600;
    }

  r_juldate = juldate_encode(calendar,
                             taxisInqRdate(taxisID),
                             taxisInqRtime(taxisID));
  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID, tsID++)) )
    {
      juldate = juldate_encode(calendar,
                               taxisInqVdate(taxisID),
                               taxisInqVtime(taxisID));
      time_val = juldate_to_seconds(juldate_sub(juldate, r_juldate))
        / tunitsec;

      if ( has_bnds )
        {
          taxisInqVdateBounds(taxisID, &vdate0b, &vdate1b);
          taxisInqVtimeBounds(taxisID, &vtime0b, &vtime1b);

          juldate = juldate_encode(calendar, vdate0b, vtime0b);
          time_bnds[0] = juldate_to_seconds(juldate_sub(juldate, r_juldate))
            / tunitsec;

          juldate = juldate_encode(calendar, vdate1b, vtime1b);
          time_bnds[1] = juldate_to_seconds(juldate_sub(juldate, r_juldate))
            / tunitsec;
          time_bndsp = time_bnds;
        }
      else
        {
          time_bndsp = NULL;
        }

      while ( nrecs-- )
        {
          streamInqRecord(streamID, &varID, &levelID);
          var = find_var(varID, vars, nvars);
          if ( var->datatype == 'f' )
            streamReadRecordF(streamID,
                              (float *)var->data + gridsize * levelID,
                              &nmiss);
          else
            streamReadRecord(streamID,
                             (double *)var->data + gridsize * levelID,
                             &nmiss);
        }

      for ( int i = 0; i < nvars; i++ )
        cmor_write(vars[i].cmor_varID,
                   vars[i].data,
                   vars[i].datatype,
                   NULL,
                   1,
                   &time_val,
                   time_bndsp,
                   NULL);
    }
}
#endif

void *CMOR(void *argument)
{
  cdoInitialize(argument);

#if defined(HAVE_LIBCMOR)
  int nparams = operatorArgc();
  char **params = operatorArgv();
  int nvars, nvars_max;
  int streamID;
  cc_var_t *vars;

  if ( nparams < 1 )
    cdoAbort("Too few arguments!");

  hcreate(100);
  parse_kvfile("cmor.rc");
  parse_cmdline_kv(nparams - 1, &params[1]);

  streamID = streamOpenRead(cdoStreamName(0));
  nvars_max = vlistNvars(streamInqVlist(streamID));
  vars = (cc_var_t *) Malloc(nvars_max * sizeof(cc_var_t));;

  setup(streamID, params[0]);
  define_variables(streamID, vars, &nvars);
  write_variables(streamID, vars, nvars);

  streamClose(streamID);
  cmor_close();

  hdestroy();
  for ( int i = 0; i < nvars; i++ )
    Free (vars[i].data);
  Free (vars);
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

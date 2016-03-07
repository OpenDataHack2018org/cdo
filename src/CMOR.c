#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

#if defined(HAVE_LIBCMOR)
#include <ctype.h>
#include "uthash.h"
#include "cmor.h"

struct cc_var
{
  int cdi_varID;
  int cmor_varID;
  char datatype;
  void *data;
};

struct kv
{
  char *key;
  char *value;
  UT_hash_handle hh;
};

static struct cc_var *find_var(int cdi_varID, struct cc_var vars[], int nvars)
{
  for ( int i = 0; i < nvars; i++ )
    if ( cdi_varID == vars[i].cdi_varID )
      return &vars[i];
  return NULL;
}

static char *trim(char *s)
{
  if (s == NULL) return s;
  while ( *s != '\0' && (isspace(*s) || *s == '"') )
    s++;
  int n = strlen(s);
  while ( n > 0 && (isspace(s[n - 1]) || s[n - 1] == '"') )
    n--;
  s[n] = '\0';
  return s;
}

static void hinsert(struct kv **ht, const char *key, const char *value)
{
  /* Insert new keys. Do not overwrite values of existing keys. */
  struct kv *e, *s;
  HASH_FIND_STR(*ht, key, s);
  if ( s == NULL)
    {
      e = Malloc(sizeof(struct kv));
      e->key = Malloc(strlen(key) + 1);
      e->value = Malloc(strlen(value) + 1);
      strcpy(e->key, key);
      strcpy(e->value, value);
      HASH_ADD_KEYPTR(hh, *ht, e->key, strlen(e->key), e);
    }
}

static void hreplace(struct kv **ht, const char *key, const char *value)
{
  /* Overwrites values of existing keys. */
  struct kv *s;
  HASH_FIND_STR(*ht, key, s);
  if ( s )
    {
      HASH_DEL(*ht, s);
      Free(s->key);
      Free(s->value);
      Free(s);
    }
  hinsert(ht, key, value);
}

static void parse_kv(struct kv **ht, char *kvstr)
{
  char *key = trim(strtok(kvstr, "="));
  char *value = trim(strtok(NULL, "="));
  if ( key && value )
    hinsert(ht, key, value);
}

static int parse_kv_file(struct kv **ht, const char *filename, int verbose)
{
  FILE *fp = fopen(filename, "r");
  if ( fp == NULL )
    {
      if ( verbose )
        cdoWarning("cannot open '%s'", filename);
      return 1;
    }

  char line[CMOR_MAX_STRING];
  while ( fgets(line, sizeof(line), fp) != NULL )
    {
      char *comment = strchr(line, '#');
      if ( comment ) *comment = '\0';
      parse_kv(ht, line);
    }
  fclose(fp);
  return 0;
}

static void parse_kv_cmdline(struct kv **ht, int nparams, char **params)
{
  /* Assume key = value pairs. That is, if params[i] contains no '='
   * then treat it as if it belongs to the value of params[i-1],
   * separated by a ','.*/
  int i = 0;
  while ( i < nparams )
    {
      int j = 1;
      int size = strlen(params[i]) + 1;
      while ( i + j < nparams && strchr(params[i + j], '=') == NULL )
        {
          size += strlen(params[i + j]) + 1;
          j++;
        }
      char *p = (char *) Malloc(size);
      strcpy(p, params[i]);
      for (int k = 1; k < j; k++)
        {
          strcat(p, ",");
          strcat(p, params[i + k]);
        }
      parse_kv(ht, p);
      Free(p);
      i += j;
    }
}

static char *get_val(struct kv **ht, char *key, char *def)
{
  struct kv *e;
  HASH_FIND_STR(*ht, key, e);
  return e ? e->value : def;
}

static char *substitute(struct kv **ht, char *word)
{
  struct kv *e;
  char *key = (char *) Malloc(strlen(word) + 12);
  sprintf(key, "substitute_%s", word);
  HASH_FIND_STR(*ht, key, e);
  Free(key);
  return e ? e->value : word;
}

static void dump_global_attributes(struct kv **ht, int streamID)
{
  int natts;
  int vlistID = streamInqVlist(streamID);
  vlistInqNatts(vlistID, CDI_GLOBAL, &natts);
  for ( int i = 0; i < natts; i++ )
    {
      char name[CDI_MAX_NAME];
      char *value = NULL;
      char buffer[8];
      int type, len;
      vlistInqAtt(vlistID, CDI_GLOBAL, i, name, &type, &len);
      switch ( type )
        {
        case DATATYPE_TXT:
          value = Malloc(len + 1);
          vlistInqAttTxt(vlistID, CDI_GLOBAL, name, len, value);
          value[len] = '\0';
          break;
        case DATATYPE_INT32:
          value = Malloc(CDI_MAX_NAME);
          vlistInqAttInt(vlistID, CDI_GLOBAL, name, len, (int *)buffer);
          snprintf(value, CDI_MAX_NAME, "%i", *(int *)buffer);
          break;
        case DATATYPE_FLT64:
          value = Malloc(CDI_MAX_NAME);
          vlistInqAttFlt(vlistID, CDI_GLOBAL, name, len, (double *)buffer);
          snprintf(value, CDI_MAX_NAME, "%e", *(double *)buffer);
          break;
        default:
          cdoWarning("Unsupported type %i name %s\n", type, name);
        }
      hinsert(ht, name, value);
      if ( value ) Free(value);
    }
}

static void dump_special_attributes(struct kv **ht, int streamID)
{
  int vlistID = streamInqVlist(streamID);
  /* Any new history will be appended to the existing history. */
  int fileID = pstreamFileID(streamID);
  size_t old_historysize = (size_t) streamInqHistorySize(fileID);
  char *new_history = get_val(ht, "history", "");
  size_t historysize;

  if ( old_historysize )
    {
      historysize = old_historysize;
      if ( new_history )
        historysize += strlen(new_history) + 1;
    }
  else
    {
      historysize = strlen(new_history);
    }

  if ( historysize )
    {
      char *history = Malloc(historysize + 1);
      memset(history, 0, historysize + 1);
      if ( old_historysize )
        {
          streamInqHistoryString(fileID, history);
          if ( new_history )
            {
              strcat(history, " ");
              strcat(history, new_history);
            }
        }
      else
        {
          strcpy(history, new_history);
        }
      hreplace(ht, "history", history);
      Free(history);
    }

  const char *value = institutInqLongnamePtr(vlistInqVarInstitut(vlistID, 0));
  if ( value ) hinsert(ht, "institution", value);

  value = modelInqNamePtr(vlistInqVarModel(vlistID, 0));
  if ( value ) hinsert(ht, "source", value);
}

static void read_config_files(struct kv **ht)
{
  /* Files from info key in command line. */
  char *info = get_val(ht, "info", "");
  char *infoc = Malloc(strlen(info) + 1);
  strcpy(infoc, info);
  char *filename = strtok(infoc, ",");
  while ( filename != NULL )
    {
      parse_kv_file(ht, trim(filename), 1);
      filename = strtok(NULL, ",");
    }
  Free(infoc);

  /* Config file in user's $HOME directory. */
  char *home = getenv("HOME");
  const char *dotconfig = ".cdocmorinfo";
  filename = Malloc(strlen(home) + strlen(dotconfig) + 2);
  sprintf(filename, "%s/%s", home, dotconfig);
  parse_kv_file(ht, filename, 0);
  Free(filename);

  /* System wide configuration. */
  parse_kv_file(ht, "/etc/cdocmor.info", 0);
}

static int in_list(char **list, const char *needle)
{
  while ( *list )
    if ( strcmp(*list++, needle) == 0 )
      return 1;
  return 0;
}

static void setup(struct kv **ht, int streamID, char *table)
{
  int netcdf_file_action;
  char *chunk = get_val(ht, "chunk", "replace");
  if ( strcasecmp(chunk, "replace") == 0 )
    netcdf_file_action = CMOR_REPLACE;
  else if ( strcasecmp(chunk, "append") == 0 )
    netcdf_file_action = CMOR_APPEND;

  int set_verbosity = CMOR_NORMAL;
  if ( strcasecmp(get_val(ht, "set_verbosity", ""), "CMOR_QUIET") == 0 )
    set_verbosity = CMOR_QUIET;

  int exit_control = CMOR_NORMAL;
  if ( strcasecmp(get_val(ht, "exit_control", ""), "CMOR_EXIT_ON_MAJOR") == 0 )
    exit_control = CMOR_EXIT_ON_MAJOR;
  if ( strcasecmp(get_val(ht, "exit_control", ""), "CMOR_EXIT_ON_WARNING")
       == 0 )
    exit_control = CMOR_EXIT_ON_WARNING;

  char *logfile = get_val(ht, "logfile", NULL);
  int create_subdirectories = atoi(get_val(ht, "create_subdirectories", "0"));
  cmor_setup(get_val(ht, "inpath", "/usr/share/cmor/"),
             &netcdf_file_action,
             &set_verbosity,
             &exit_control,
             logfile,
             &create_subdirectories);

  char *calendar;
  int taxisID = vlistInqTaxis(streamInqVlist(streamID));
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

  int *month_lengths;
  char *ml = get_val(ht, "month_lengths", NULL);
  if ( ml )
    {
      char *mlc = Malloc(strlen(ml) + 1);
      char *month_str = strtok(mlc, ",");
      int month = 0;
      month_lengths = Malloc(12 * sizeof(int));
      while ( month < 12 && month_str != NULL )
        {
          month_lengths[month++] = atoi(month_str);
          month_str = strtok(NULL, ",");
        }
      Free(mlc);
      if ( month != 12 )
        cdoAbort("Invalid format for month_lengths");
    }
  else
    {
      month_lengths = NULL;
    }

  double branch_time = atof(get_val(ht, "branch_time", "0.0"));
  cmor_dataset(get_val(ht, "outpath", "./"),
               get_val(ht, "experiment_id", ""),
               get_val(ht, "institution", ""),
               get_val(ht, "source", ""),
               calendar,
               atoi(get_val(ht, "realization", "1")),
               get_val(ht, "contact", ""),
               get_val(ht, "history", ""),
               get_val(ht, "comment", ""),
               get_val(ht, "references", ""),
               atoi(get_val(ht, "leap_year", "0")),
               atoi(get_val(ht, "leap_month", "0")),
               month_lengths,
               get_val(ht, "model_id", ""),
               get_val(ht, "forcing", ""),
               atoi(get_val(ht, "initialization_method", "1")),
               atoi(get_val(ht, "physics_version", "1")),
               get_val(ht, "institute_id", ""),
               get_val(ht, "parent_experiment_id", ""),
               &branch_time,
               get_val(ht, "parent_experiment_rip", ""));

  int table_id;
  cmor_load_table(table, &table_id);
  cmor_set_table(table_id);
}

static void define_variables(struct kv **ht, int streamID,
                             struct cc_var vars[], int *nvars)
{
  int vlistID = streamInqVlist(streamID);
  int taxisID = vlistInqTaxis(vlistID);
  char taxis_units[CMOR_MAX_STRING];
  int timeunit = taxisInqTunit(taxisID);
  int year, month, day, hour, minute, second;
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

  char **name_list;
  char *select_vars = get_val(ht, "var", NULL);
  if ( select_vars )
    {
      name_list = Malloc((strlen(select_vars) + 1) * sizeof(char *));
      char *var_name = strtok(select_vars, ",");
      int i = 0;
      while ( var_name != NULL )
        {
          name_list[i++] = trim(var_name);
          var_name = strtok(NULL, ",");
        }
      name_list[i] = NULL;
    }
  else
    {
      name_list = NULL;
    }

  *nvars = 0;
  for ( int varID = 0; varID < vlistNvars(vlistID); varID++ )
    {
      char name[CDI_MAX_NAME];
      int axis_ids[CMOR_MAX_AXES];
      vlistInqVarName(vlistID, varID, name);
      if ( name_list == NULL || in_list(name_list, name) )
        {
          int ndims = 0;
          /* Time-Axis */
          cmor_axis(&axis_ids[ndims++],
                    substitute(ht, "time"),
                    taxis_units,
                    0,
                    NULL,
                    0,
                    NULL,
                    0,
                    NULL);

          /* Z-Axis */
          int zaxisID = vlistInqVarZaxis(vlistID, varID);
          int levels = zaxisInqSize(zaxisID);
          char units[CDI_MAX_NAME];
          double *coord_vals;
          if ( zaxisInqType(zaxisID) != ZAXIS_SURFACE )
            {
              coord_vals = Malloc(levels * sizeof(double));
              zaxisInqLevels(zaxisID, coord_vals);
              zaxisInqName(zaxisID, name);
              zaxisInqUnits(zaxisID, units);
              cmor_axis(&axis_ids[ndims++],
                        substitute(ht, name),
                        units,
                        levels,
                        (void *)coord_vals,
                        'd',
                        NULL,
                        0,
                        NULL);
            }

          /* Y-Axis */
          int gridID = vlistInqVarGrid(vlistID, varID);
          gridInqYname(gridID, name);
          gridInqYunits(gridID, units);
          int length = gridInqYsize(gridID);
          coord_vals = Malloc(length * sizeof(double));
          gridInqYvals(gridID, coord_vals);
          double *cell_bounds = Malloc(2 * length * sizeof(double));
          int nbounds = gridInqYbounds(gridID, cell_bounds);
          if ( nbounds != 2 * length )
            {
              Free(cell_bounds);
              cell_bounds = NULL;
            }
          cmor_axis(&axis_ids[ndims++],
                    substitute(ht, name),
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
          nbounds = gridInqXbounds(gridID, cell_bounds);
          if ( nbounds != 2 * length )
            {
              Free(cell_bounds);
              cell_bounds = NULL;
            }
          cmor_axis(&axis_ids[ndims++],
                    substitute(ht, name),
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
          char missing_value[sizeof(double)];
          double tolerance = 1e-4;
          size_t gridsize = vlistGridsizeMax(vlistID);
          struct cc_var *var = &vars[(*nvars)++];
          var->cdi_varID = varID;
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
                        substitute(ht, name),
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
  if ( name_list ) Free(name_list);
}

static void write_variables(int streamID, struct cc_var vars[], int nvars)
{
  int vlistID = streamInqVlist(streamID);
  int taxisID = vlistInqTaxis(vlistID);
  int tunitsec;
  switch ( taxisInqTunit(taxisID) )
    {
    case TUNIT_MINUTE: tunitsec = 60; break;
    case TUNIT_HOUR: tunitsec = 3600; break;
    case TUNIT_DAY: tunitsec = 86400; break;
    default: tunitsec = 3600;
    }

  int calendar = taxisInqCalendar(taxisID);
  juldate_t r_juldate = juldate_encode(calendar,
                                       taxisInqRdate(taxisID),
                                       taxisInqRtime(taxisID));
  size_t gridsize = vlistGridsizeMax(vlistID);
  double *buffer = (double *) Malloc(gridsize * sizeof(double));
  int tsID = 0;
  int nrecs;
  while ( (nrecs = streamInqTimestep(streamID, tsID++)) )
    {
      double time_bnds[2];
      double *time_bndsp;
      double time_val;
      juldate_t juldate = juldate_encode(calendar,
                                         taxisInqVdate(taxisID),
                                         taxisInqVtime(taxisID));
      time_val = juldate_to_seconds(juldate_sub(juldate, r_juldate))
        / tunitsec;

      if ( taxisHasBounds(taxisID) )
        {
          int vdate0b, vdate1b, vtime0b, vtime1b;
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
          int varID, levelID;
          streamInqRecord(streamID, &varID, &levelID);
          struct cc_var *var = find_var(varID, vars, nvars);
          if ( var )
            {
              int nmiss;
              if ( var->datatype == 'f' )
                {
                  streamReadRecord(streamID, buffer, &nmiss);
                  for ( size_t i = 0; i < gridsize; i++ )
                    ((float *)var->data)[gridsize * levelID + i] =
                      (float)buffer[i];
                }
              else
                {
                  streamReadRecord(streamID,
                                   (double *)var->data + gridsize * levelID,
                                   &nmiss);
                }
            }
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
  Free(buffer);
}
#endif

void *CMOR(void *argument)
{
  cdoInitialize(argument);

#if defined(HAVE_LIBCMOR)
  int nparams = operatorArgc();
  char **params = operatorArgv();
  struct kv *ht = NULL;
  if ( nparams < 1 ) cdoAbort("Too few arguments!");

  /* Command line config has highest priority. */
  parse_kv_cmdline(&ht, nparams - 1, &params[1]);

  /* Config files are read with descending priority. */
  read_config_files(&ht);

  int streamID = streamOpenRead(cdoStreamName(0));
  /* Existing attributes have lowest priority. */
  dump_global_attributes(&ht, streamID);
  dump_special_attributes(&ht, streamID);

  int nvars;
  int nvars_max = vlistNvars(streamInqVlist(streamID));
  struct cc_var *vars = (struct cc_var *) Malloc(nvars_max
                                                 * sizeof(struct cc_var));
  setup(&ht, streamID, params[0]);
  define_variables(&ht, streamID, vars, &nvars);
  write_variables(streamID, vars, nvars);

  streamClose(streamID);
  cmor_close();

  for ( int i = 0; i < nvars; i++ )
    Free(vars[i].data);
  Free(vars);

  struct kv *s, *tmp;
  HASH_ITER(hh, ht, s, tmp)
    {
      Free(s->key);
      Free(s->value);
      Free(s);
    }
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

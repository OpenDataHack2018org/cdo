#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

#if defined(HAVE_LIBCMOR)
#include <ctype.h>
#include <unistd.h>
#include "uthash.h"
#include "cmor.h"

#define CMOR_UNDEFID (CMOR_MAX_AXES + 1)

struct mapping
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

static char *key_rename(struct kv **ht, char *word)
{
  struct kv *e;
  size_t len = strlen(word) + 12;
  char *key = (char *) Malloc(len);
  snprintf(key, len, "__rename_%s", word);
  HASH_FIND_STR(*ht, key, e);
  Free(key);
  return e ? e->value : word;
}

static void dump_global_attributes(struct kv **ht, int streamID)
{
  int natts;
  int vlistID = streamInqVlist(streamID);
  cdiInqNatts(vlistID, CDI_GLOBAL, &natts);
  for ( int i = 0; i < natts; i++ )
    {
      char name[CDI_MAX_NAME];
      char *value = NULL;
      char buffer[8];
      int type, len;
      cdiInqAtt(vlistID, CDI_GLOBAL, i, name, &type, &len);
      switch ( type )
        {
        case CDI_DATATYPE_TXT:
          value = Malloc(len + 1);
          cdiInqAttTxt(vlistID, CDI_GLOBAL, name, len, value);
          value[len] = '\0';
          break;
        case CDI_DATATYPE_INT32:
          value = Malloc(CDI_MAX_NAME);
          cdiInqAttInt(vlistID, CDI_GLOBAL, name, len, (int *)buffer);
          snprintf(value, CDI_MAX_NAME, "%i", *(int *)buffer);
          break;
        case CDI_DATATYPE_FLT64:
          value = Malloc(CDI_MAX_NAME);
          cdiInqAttFlt(vlistID, CDI_GLOBAL, name, len, (double *)buffer);
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
  int fileID = pstreamFileID(streamID);
  size_t old_historysize;
  char *new_history = get_val(ht, "history", "");
  size_t historysize;
  int natts;
  cdiInqNatts(vlistID, CDI_GLOBAL, &natts);
  if ( natts > 0 )
    old_historysize = (size_t) streamInqHistorySize(fileID);
  else
    old_historysize = 0;

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
  char *info = get_val(ht, "__info", "");
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

static void register_cmor_calendar(struct kv **ht, int calendar)
{
  switch ( calendar )
    {
    case CALENDAR_STANDARD:
      hinsert(ht, "calendar", "gregorian");
      break;
    case CALENDAR_PROLEPTIC:
      hinsert(ht, "calendar", "proleptic_gregorian");
      break;
    case CALENDAR_360DAYS:
      hinsert(ht, "calendar", "360_day");
      break;
    case CALENDAR_365DAYS:
      hinsert(ht, "calendar", "noleap");
      break;
    case CALENDAR_366DAYS:
      hinsert(ht, "calendar", "all_leap");
      break;
    default:
      cdoWarning("Unsupported calendar type %i.\n", calendar);
    }
}

static char *dump_ht_to_json_file(struct kv **ht, char *filename)
{
  int fd = mkstemp(filename);
  dprintf(fd, "{\n");
  for ( struct kv *entry = *ht; entry != NULL; entry = entry->hh.next )
    if ( strncmp(entry->key, "__", 2) != 0 )
      dprintf(fd, "\t\"%s\":\t\t\"%s\",\n", entry->key, entry->value);
  dprintf(fd, "}\n");
  close(fd);
  return filename;
}

static int get_netcdf_file_action(struct kv **ht)
{
  char *chunk = get_val(ht, "__chunk", "replace");
  if ( strcasecmp(chunk, "append") == 0 )
    return CMOR_APPEND;
  else
    return CMOR_REPLACE;
}

static int get_cmor_verbosity(struct kv **ht)
{
  if ( strcasecmp(get_val(ht, "__set_verbosity", ""), "CMOR_QUIET") == 0 )
    return CMOR_QUIET;
  else
    return CMOR_NORMAL;
}

static int get_cmor_exit_control(struct kv **ht)
{
  if ( strcasecmp(get_val(ht, "__exit_control", ""),
                  "CMOR_EXIT_ON_MAJOR") == 0 )
    return CMOR_EXIT_ON_MAJOR;
  else if ( strcasecmp(get_val(ht, "__exit_control", ""),
                       "CMOR_EXIT_ON_WARNING")  == 0 )
    return CMOR_EXIT_ON_WARNING;
  else
    return CMOR_NORMAL;
}

static void feed_json_to_cmor_dataset(struct kv **ht)
{
  char filename[20];
  snprintf(filename, sizeof(filename), "dataset.json_XXXXXX");
  cmor_dataset_json(dump_ht_to_json_file(ht, filename));
  unlink(filename);
}

static void setup_dataset(struct kv **ht, int streamID)
{
  int netcdf_file_action = get_netcdf_file_action(ht);
  int set_verbosity = get_cmor_verbosity(ht);
  int exit_control = get_cmor_exit_control(ht);
  char *logfile = get_val(ht, "__logfile", NULL);
  int create_subdirectories = atoi(get_val(ht, "__create_subdirectories", "0"));

  cmor_setup(get_val(ht, "__inpath", "/usr/share/cmor/"),
             &netcdf_file_action, &set_verbosity, &exit_control,
             logfile, &create_subdirectories);

  int taxisID = vlistInqTaxis(streamInqVlist(streamID));
  register_cmor_calendar(ht, taxisInqCalendar(taxisID));
  feed_json_to_cmor_dataset(ht);
}

static int *new_axis_id(int *axis_ids)
{
  int i;
  for ( i = 0; axis_ids[i] != CMOR_UNDEFID; i++ );
  axis_ids[i + 1] = CMOR_UNDEFID;
  return &axis_ids[i];
}

static int count_axis_ids(int *axis_ids)
{
  int i;
  for ( i = 0; axis_ids[i] != CMOR_UNDEFID; i++ );
  return i;
}

static void register_x_axis(int gridID, char* name, int *axis_ids)
{
  double *coord_vals;
  char units[CDI_MAX_NAME];
  gridInqXunits(gridID, units);
  int length = gridInqXsize(gridID);
  coord_vals = Malloc(length * sizeof(double));
  gridInqXvals(gridID, coord_vals);
  double *cell_bounds = Malloc(2 * length * sizeof(double));
  int nbounds = gridInqXbounds(gridID, cell_bounds);
  if ( nbounds != 2 * length )
    {
      Free(cell_bounds);
      cell_bounds = NULL;
    }
  cmor_axis(new_axis_id(axis_ids), name, units, length, (void *)coord_vals,
            'd', (void *)cell_bounds, 2, NULL);
  Free(coord_vals);
  if ( cell_bounds ) Free(cell_bounds);
}

static void register_y_axis(int gridID, char* name, int *axis_ids)
{
  double *coord_vals;
  char units[CDI_MAX_NAME];
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
  cmor_axis(new_axis_id(axis_ids), name, units, length, (void *)coord_vals,
            'd', (void *)cell_bounds, 2, NULL);
  Free(coord_vals);
  if ( cell_bounds ) Free(cell_bounds);
}

static void register_z_axis(int zaxisID, char *name, int *axis_ids)
{
  int levels = zaxisInqSize(zaxisID);
  char units[CDI_MAX_NAME];
  double *coord_vals;
  if ( zaxisInqType(zaxisID) != ZAXIS_SURFACE )
    {
      coord_vals = Malloc(levels * sizeof(double));
      cdoZaxisInqLevels(zaxisID, coord_vals);
      zaxisInqUnits(zaxisID, units);
      cmor_axis(new_axis_id(axis_ids), name, units, levels, (void *)coord_vals,
                'd', NULL, 0, NULL);
      Free(coord_vals);
    }
}

static void register_xy_only(int gridID, int *axis_ids)
{
  char name[CDI_MAX_NAME];
  name[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_YLONGNAME, CDI_MAX_NAME, name);
  register_y_axis(gridID, name, axis_ids);
  name[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_XLONGNAME, CDI_MAX_NAME, name);
  register_x_axis(gridID, name, axis_ids);
}

static void register_cmor_grid(int gridID, int *axis_ids, int *cmor_grid_id)
{
  int gridsize = gridInqSize(gridID);
  double *latitude = Malloc(gridsize * sizeof(double));
  gridInqYvals(gridID, latitude);
  double *latitude_vertices = Malloc(4 * gridsize * sizeof(double));
  gridInqYbounds(gridID, latitude_vertices);

  double *longitude = Malloc(gridsize * sizeof(double));
  gridInqXvals(gridID, longitude);
  double *longitude_vertices = Malloc(4 * gridsize * sizeof(double));
  gridInqXbounds(gridID, longitude_vertices);
  cmor_grid(cmor_grid_id, count_axis_ids(axis_ids), axis_ids, 'd',
            (void *)latitude, (void *)longitude, 4,
            (void *)latitude_vertices, (void *)longitude_vertices
            );
  Free(latitude);
  Free(latitude_vertices);
  Free(longitude);
  Free(longitude_vertices);
}

static void register_cmor_grid_mapping(int projID, int cmor_grid_id)
{
  char grid_mapping[CDI_MAX_NAME] = "";
  cdiGridInqKeyStr(projID, CDI_KEY_MAPPING, CDI_MAX_NAME, grid_mapping);

  if ( grid_mapping[0] )
    {
      int atttype, attlen;
      char attname[CDI_MAX_NAME];
      int natts, nparameters;
      cdiInqNatts(projID, CDI_GLOBAL, &natts);
      char parameter_names[natts][CDI_MAX_NAME];
      double parameter_values[natts];
      char parameter_units[natts][1];

      nparameters = 0;
      for ( int iatt = 0; iatt < natts; ++iatt )
        {
          cdiInqAtt(projID, CDI_GLOBAL, iatt, attname, &atttype, &attlen);
          if ( atttype == CDI_DATATYPE_FLT32 || atttype == CDI_DATATYPE_FLT64 )
            {
              double attflt[attlen];
              cdiInqAttFlt(projID, CDI_GLOBAL, attname, attlen, attflt);
              strcpy(parameter_names[nparameters], attname);
              parameter_values[nparameters] = attflt[0];
              parameter_units[nparameters][0] = 0;
              nparameters++;
            }
        }
      cmor_set_grid_mapping(cmor_grid_id, grid_mapping, nparameters,
                            (char **)parameter_names, CDI_MAX_NAME,
                            parameter_values,
                            (char **)parameter_units, 1);
    }
}

static void register_projected_grid(int gridID, int *axis_ids)
{
  char name[CDI_MAX_NAME];
  int projID = gridInqProj(gridID);
  int proj_axis_ids[3];
  proj_axis_ids[0] = CMOR_UNDEFID;
  gridInqYstdname(projID, name);
  register_y_axis(projID, name, proj_axis_ids);
  gridInqXstdname(projID, name);
  register_x_axis(projID, name, proj_axis_ids);
  int *cmor_grid_id = new_axis_id(axis_ids);
  register_cmor_grid(gridID, proj_axis_ids, cmor_grid_id);
  register_cmor_grid_mapping(projID, *cmor_grid_id);
}

static void register_xy_and_grid(int gridID, int table_id, int grid_table_id,
                                 int *axis_ids)
{
  if ( gridInqProj(gridID) > 0 )
    {
      if ( grid_table_id == 0 )
        cdoAbort("__grid_table not specified!");
      cmor_set_table(grid_table_id);
      register_projected_grid(gridID, axis_ids);
      cmor_set_table(table_id);
    }
  else
    {
      register_xy_only(gridID, axis_ids);
    }
}

static void get_taxis_units(char *units, int taxisID)
{
  int timeunit = taxisInqTunit(taxisID);
  int year, month, day, hour, minute, second;
  cdiDecodeDate(taxisInqRdate(taxisID), &year, &month, &day);
  cdiDecodeTime(taxisInqRtime(taxisID), &hour, &minute, &second);
  if ( timeunit == TUNIT_QUARTER || timeunit == TUNIT_30MINUTES )
    timeunit = TUNIT_MINUTE;
  if ( timeunit == TUNIT_3HOURS || timeunit == TUNIT_6HOURS ||
       timeunit == TUNIT_12HOURS )
    timeunit = TUNIT_HOUR;

  snprintf(units, CMOR_MAX_STRING, "%s since %d-%d-%d %02d:%02d:%02d", tunitNamePtr(timeunit),
          year, month, day, hour, minute, second);
}

static int get_table_id(char *table)
{
  int table_id;
  cmor_load_table(table, &table_id);
  return table_id;
}

static int get_grid_table_id(struct kv **ht)
{
  if ( get_val(ht, "__grid_table", NULL) )
    {
      return get_table_id(get_val(ht, "__grid_table", NULL));
    }
  else
    {
      cdoAbort("Grid table required but not defined. Set __grid_table!");
      return 0;
    }
}

static char **get_requested_variables(struct kv **ht)
{
  char **name_list = NULL;
  char *select_vars = get_val(ht, "__var", NULL);

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
  return name_list;
}

static void register_variable(int vlistID, int varID, int *axis_ids,
                              struct mapping *var)
{
  char name[CDI_MAX_NAME];
  vlistInqVarName(vlistID, varID, name);
  char units[CDI_MAX_NAME];
  vlistInqVarUnits(vlistID, varID, units);
  char missing_value[sizeof(double)];
  double tolerance = 1e-4;
  size_t gridsize = vlistGridsizeMax(vlistID);
  int levels = zaxisInqSize(vlistInqVarZaxis(vlistID, varID));
  var->cdi_varID = varID;
  if ( vlistInqVarDatatype(vlistID, varID) == CDI_DATATYPE_FLT32 )
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
  cmor_variable(&var->cmor_varID, name, units,
                count_axis_ids(axis_ids), axis_ids,
                var->datatype, (void *) missing_value, &tolerance,
                NULL, NULL, NULL, NULL);
}

static struct mapping *new_var_mapping(struct mapping vars[])
{
  int i;
  for ( i = 0; vars[i].cdi_varID != CDI_UNDEFID; i++ );
  vars[i + 1].cdi_varID = CDI_UNDEFID;
  return &vars[i];
}

static int in_list(char **list, const char *needle)
{
  for ( ; *list; list++ )
    if ( strcmp(*list, needle) == 0 )
      return 1;
  return 0;
}

static void register_axes_and_variables(struct kv **ht, int table_id,
                                        int streamID, struct mapping vars[])
{
  int vlistID = streamInqVlist(streamID);
  int taxisID = vlistInqTaxis(vlistID);
  char taxis_units[CMOR_MAX_STRING];
  get_taxis_units(taxis_units, taxisID);
  char **requested_variables = get_requested_variables(ht);
  int grid_table_id = get_grid_table_id(ht);
  cmor_set_table(table_id);

  for ( int varID = 0; varID < vlistNvars(vlistID); varID++ )
    {
      char var_name[CDI_MAX_NAME];
      int axis_ids[CMOR_MAX_AXES];
      axis_ids[0] = CMOR_UNDEFID;
      vlistInqVarName(vlistID, varID, var_name);
      if ( requested_variables == NULL ||
           in_list(requested_variables, var_name) )
        {
          cmor_axis(new_axis_id(axis_ids), "time", taxis_units, 0,
                    NULL, 0, NULL, 0, NULL);
          int zaxisID = vlistInqVarZaxis(vlistID, varID);
          char z_name[CDI_MAX_NAME];
          zaxisInqName(zaxisID, z_name);
          register_z_axis(zaxisID, key_rename(ht, z_name), axis_ids);
          int gridID = vlistInqVarGrid(vlistID, varID);
          register_xy_and_grid(gridID, table_id, grid_table_id, axis_ids);
          register_variable(vlistID, varID, axis_ids, new_var_mapping(vars));
        }
    }
  if ( requested_variables ) Free(requested_variables);
}

static int time_unit_in_seconds(int taxisID)
{
  switch ( taxisInqTunit(taxisID) )
    {
    case TUNIT_MINUTE: return 60;
    case TUNIT_HOUR: return 3600;
    case TUNIT_DAY: return 86400;
    default: return 3600;
    }
}

static double get_cmor_time_val(int taxisID, juldate_t ref_date)
{
  int tunitsec = time_unit_in_seconds(taxisID);
  int calendar = taxisInqCalendar(taxisID);
  juldate_t juldate = juldate_encode(calendar, taxisInqVdate(taxisID),
                                     taxisInqVtime(taxisID));
  return juldate_to_seconds(juldate_sub(juldate, ref_date)) / tunitsec;
}

static double *get_cmor_time_bounds(int taxisID, juldate_t ref_date,
                                    double *bounds)
{
  if ( taxisHasBounds(taxisID) )
    {
      int vdate0b, vdate1b, vtime0b, vtime1b;
      taxisInqVdateBounds(taxisID, &vdate0b, &vdate1b);
      taxisInqVtimeBounds(taxisID, &vtime0b, &vtime1b);

      int tunitsec = time_unit_in_seconds(taxisID);
      int calendar = taxisInqCalendar(taxisID);
      juldate_t date = juldate_encode(calendar, vdate0b, vtime0b);
      bounds[0] = juldate_to_seconds(juldate_sub(date, ref_date)) / tunitsec;

      date = juldate_encode(calendar, vdate1b, vtime1b);
      bounds[1] = juldate_to_seconds(juldate_sub(date, ref_date)) / tunitsec;
      return bounds;
    }
  else
    {
      return NULL;
    }
}

static struct mapping *map_var(int cdi_varID, struct mapping vars[])
{
  for ( int i = 0; vars[i].cdi_varID != CDI_UNDEFID; i++ )
    if ( cdi_varID == vars[i].cdi_varID )
      return &vars[i];
  return NULL;
}

static void read_record(int streamID, double *buffer, size_t gridsize,
                        struct mapping vars[])
{
  int varID, levelID;
  streamInqRecord(streamID, &varID, &levelID);
  struct mapping *var = map_var(varID, vars);
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
          streamReadRecord(streamID, (double *)var->data + gridsize * levelID,
                           &nmiss);
        }
    }
}

static void write_variables(int streamID, struct mapping vars[])
{
  int vlistID = streamInqVlist(streamID);
  int taxisID = vlistInqTaxis(vlistID);
  int calendar = taxisInqCalendar(taxisID);
  juldate_t ref_date = juldate_encode(calendar, taxisInqRdate(taxisID),
                                      taxisInqRtime(taxisID));
  size_t gridsize = vlistGridsizeMax(vlistID);
  double *buffer = (double *) Malloc(vlistGridsizeMax(vlistID) *
                                     sizeof(double));
  int tsID = 0;
  int nrecs;
  while ( (nrecs = streamInqTimestep(streamID, tsID++)) )
    {
      while ( nrecs-- )
        read_record(streamID, buffer, gridsize, vars);

      double time_val = get_cmor_time_val(taxisID, ref_date);
      double time_bnds[2];
      for ( int i = 0; vars[i].cdi_varID != CDI_UNDEFID; i++ )
        cmor_write(vars[i].cmor_varID, vars[i].data, vars[i].datatype,
                   1, &time_val,
                   get_cmor_time_bounds(taxisID, ref_date, time_bnds), NULL);
    }
  Free(buffer);
}

static void destruct_hash_table(struct kv **ht)
{
  struct kv *s, *tmp;
  HASH_ITER(hh, *ht, s, tmp)
    {
      Free(s->key);
      Free(s->value);
      Free(s);
    }
}

static void destruct_var_mapping(struct mapping vars[])
{
  for ( int i = 0; vars[i].cdi_varID != CDI_UNDEFID; i++ )
    Free(vars[i].data);
  Free(vars);
}

static struct mapping *construct_var_mapping(int streamID)
{
  int nvars_max = vlistNvars(streamInqVlist(streamID));
  struct mapping *vars =
    (struct mapping *) Malloc((nvars_max + 1) * sizeof(struct mapping));
  vars[0].cdi_varID = CDI_UNDEFID;
  return vars;
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

  setup_dataset(&ht, streamID);
  int table_id = get_table_id(params[0]);
  cmor_set_table(table_id);

  struct mapping *vars = construct_var_mapping(streamID);
  register_axes_and_variables(&ht, table_id, streamID, vars);
  write_variables(streamID, vars);
  destruct_var_mapping(vars);
  destruct_hash_table(&ht);

  streamClose(streamID);
  cmor_close();

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

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

#if defined(HAVE_LIBCMOR)
#include <ctype.h>
#include <unistd.h>
#include "uthash.h"
#include "cmor.h"
#include "netcdf.h"

#define CMOR_UNDEFID (CMOR_MAX_AXES + 1)

struct kv
{
  char *key;
  char *value;
  UT_hash_handle hh;
};

struct mapping
{
  int help_var;
  int cdi_varID;
  int cmor_varID;
  char datatype;
  void *data;
};

static struct mapping *map_var(int cdi_varID, struct mapping vars[])
{
  for ( int i = 0; vars[i].cdi_varID != CDI_UNDEFID; i++ )
    if ( cdi_varID == vars[i].cdi_varID )
      return &vars[i];
  return NULL;
}

static struct mapping *new_var_mapping(struct mapping vars[])
{
  int i;
  for ( i = 0; vars[i].cdi_varID != CDI_UNDEFID; i++ );
  vars[i + 1].cdi_varID = CDI_UNDEFID;
  return &vars[i];
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
  char *keylower = trim(strtok(kvstr, "="));

  if ( key )
    {
      int i = 0;
      while( key[i] ) 
       {
         keylower[i]=tolower(key[i]);
         i++;
       }
      if ( value )
        hinsert(ht, keylower, value);
    }
}

static int file_exist(const char *tfilename, int force)
{
  FILE *tfp = fopen(tfilename, "r");
  if ( tfp == NULL && force )
      cdoAbort("cannot open '%s'", tfilename);
  if ( tfp == NULL && !force )
    {
      cdoWarning("cannot open '%s'", tfilename);
      return 0;
    }
  fclose(tfp);
  return 1;
}
  
static int parse_kv_file(struct kv **ht, const char *filename, int verbose)
{
  if ( file_exist(filename, 0) )
    {
      FILE *fp = fopen(filename, "r");

      char line[CMOR_MAX_STRING];
      while ( fgets(line, sizeof(line), fp) != NULL )
        {
          char *comment = strchr(line, '#');
          if ( comment ) *comment = '\0';
          parse_kv(ht, line);
        }
      fclose(fp);
    }
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

static void check_compare_set(char *finalset, char *attribute, char *attname)
{
  if ( strcmp(attribute, "") != 0 )
    if ( strcmp(attribute, finalset) != 0 )
      {
        cdoWarning("%s of variable in input file: '%s' does not agree with configuration attribute %s: '%s'.\nCmor libary is called with attribute unit '%s'.\n", attname, finalset, attname, attribute, attribute);
        strcpy(finalset, attribute);
      }
  else
    if ( !finalset )
      cdoAbort("Required attribute '%s' is not found in the configuration.", attname);
}

static int check_attr(struct kv **ht)
{
  struct kv *s, *project, *pos;
  const char *reqAtt[] = {"institute_id", "institution", "contact", "model_id", "source",
            "experiment_id", "req_time_units", "project_id", NULL};
  const char *reqAttCMIP5[] = {"product", "member", NULL};
  const char *reqAttCORDEX[] = {"product", "member", "cordex_domain", "driving_model_id", NULL};
/* In all Projects needed Attributes are tested first */
  int i = 0;
  while ( reqAtt[i] != NULL )
    {
      HASH_FIND_STR(*ht, reqAtt[i], s);
      if ( !s || strcmp(s->value, "notSet") == 0 ) cdoAbort("Attribute '%s' is required. Either it is missing or notSet", reqAtt[i]);
      printf("Attribute %s is %s \n", reqAtt[i], s->value);
      i++;
    }
/* Set default attributes */
  HASH_FIND_STR(*ht, "references", s);
  if ( !s || strcmp(s->value, "notSet") == 0 )
    {
      char *references = (char *) Malloc(strlen(get_val(ht, "model_id", "")) + 28);
      strcpy(references, "No references available for ");
      strcat(references, get_val(ht, "model_id", ""));
      cdoWarning("Attribute references is set to '%s' ", references);
      hreplace(ht, "references", references);  
    }
/* Special check for CMIP5 or CORDEX projects */
  i=0;
  HASH_FIND_STR(*ht, "project_id", project);
  if ( strcmp(project->value, "CMIP6") == 0 )
    {
      cdoAbort("Not yet possible to create data for project CMIP6 since cmor version 2.9 is used in this operator.\n");
    }
  if ( strcmp(project->value, "CMIP5") == 0 )
    {
      printf("Since the project id is %s further attributes are tested. \n", project->value);
      while ( reqAttCMIP5[i] != NULL )
        {
          HASH_FIND_STR(*ht, reqAttCMIP5[i], s);
          if ( !s || strcmp(s->value, "notSet") == 0 )cdoAbort("Attribute '%s' is required. Either it is missing or notSet", reqAttCMIP5[i]);
          printf("Attribute %s is %s \n", reqAttCMIP5[i], s->value);
          i++;
        }
    }
  else if (strcmp(project->value, "CORDEX") == 0 )
    {
      printf("Since the project id is %s further attributes are tested", project->value);
      i=0;
      while ( reqAttCORDEX[i] != NULL )
        {
          HASH_FIND_STR(*ht, reqAttCORDEX[i], s);
          if ( !s || strcmp(s->value, "notSet") == 0 ) cdoAbort("Attribute '%s' is required. Either it is missing or notSet", reqAttCORDEX[i]);
          printf("Attribute %s is %s \n", reqAttCORDEX[i], s->value);
          i++;
        }
    }
  HASH_FIND_STR(*ht, "positive", pos);
  if ( pos )
    if ( ( strcmp(pos->value, "notSet") == 0 ) || ( strcmp(pos->value, "") != 0 && strcmp(pos->value, "u") != 0 && strcmp(pos->value, "d") != 0 ) )
      {
        cdoWarning("Invalid value '%s' is set for attribute 'positive'. The default (blank) is used.\n", pos-> value);
        hreplace(ht, "positive", "");
      }
  return 0;
} 

static int check_mem(struct kv **ht)
{
  char *member = get_val(ht, "member", "");
  char *project_id = get_val(ht, "project_id", "");
  char *crealiz, *cinitial, *cphysics;
  char workchar[CMOR_MAX_STRING]; 
  int realization, initialization_method, physics_version;
  int ipos=0, ppos=0;

/* Test for the right member, else abort or warn */ 
  if ( strlen(member) >= 6 && member[0] == 'r' )
    {
      crealiz = cinitial = cphysics = (char *) Malloc(strlen(member));
      strcpy(crealiz, &member[1]);
      if ( strtok_r(crealiz, "i", &cinitial) )
        {
          strtok_r(cinitial, "p", &cphysics); 
          realization = strtol(crealiz, NULL, 10);
          initialization_method = strtol(cinitial, NULL, 10);
          physics_version = strtol(cphysics, NULL, 10);
        }
      else cphysics=NULL;
    }
  else crealiz=cinitial=cphysics=NULL;
  if ( realization && initialization_method && physics_version)
    {
      strcpy(workchar, "realization=");
      strcat(workchar, crealiz);
      printf("If no other command line values were set, \n%s \n", workchar);
      parse_kv(ht, workchar); 
      strcpy(workchar, "initialization_method=");
      strcat(workchar, cinitial);
      printf("%s \n", workchar);
      parse_kv(ht, workchar);
      strcpy(workchar, "physics_version=");
      strcat(workchar, cphysics);
      printf("%s \n", workchar);
      parse_kv(ht, workchar);
      return 1;
    }
/* Now abort or warn */ 
  if (strcmp(project_id, "CMIP5") == 0 || strcmp(project_id, "CORDEX") == 0)
    cdoAbort("The member has no RIP format (at least 6 characters and in RIP order)! Found for \n member: %s. This is interpreted as \n Realization: %s \n Initialization: %s \n Physics: %s \n   But three Integers are needed", member, crealiz, cinitial, cphysics);
  else if ( strcmp(member, "notSet") == 0 )
    {
      cdoWarning("The member has no RIP format! We set \n Attribute realization=-1 \n Attribute initialization_method=-1 \n Attribute physics_version=-1 \n");
      hreplace(ht, "realization", "-1"); 
      hreplace(ht, "initialization_method", "-1"); 
      hreplace(ht, "physics_version", "-1"); 
    }
  return 0;
} 


/*
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
*/

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
}

static void read_config_files(struct kv **ht)
{
  /* Files from info key in command line. */
  char *info = get_val(ht, "info", "");
  char *workfile = Malloc(strlen(info) + 1);
  strcpy(workfile, info);
  char *restfile;
  while ( strtok_r(workfile, ",", &restfile) != NULL )
    {
      parse_kv_file(ht, trim(workfile), 1);
      strcpy(workfile, restfile);
    }

  /* Config file in user's $HOME directory. */
  char *home = getenv("HOME");
  const char *dotconfig = ".cdocmorinfo";
  workfile = Malloc(strlen(home) + strlen(dotconfig) + 2);
  sprintf(workfile, "%s/%s", home, dotconfig);
  parse_kv_file(ht, workfile, 0);
  Free(workfile);

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

static int get_netcdf_file_action(struct kv **ht)
{
  char *chunk = get_val(ht, "oflag", "");
  if ( strcmp(chunk, "replace") == 0 )
    return CMOR_REPLACE;
  else if ( strcmp(chunk, "append") == 0 )
    return CMOR_APPEND;
  else if ( strcmp(chunk, "preserve") == 0 )
    return CMOR_APPEND;
  else
    {
      cdoWarning("No valid CMOR output mode! \nAttribute oflag is '%s', but valid are 'append', 'replace' or 'preserve'.\nCMOR output mode is set to: replace.\n", chunk);
      return CMOR_REPLACE;
    }
}

static int get_cmor_verbosity(struct kv **ht)
{
  if ( strcmp(get_val(ht, "set_verbosity", ""), "CMOR_QUIET") == 0 )
    return CMOR_QUIET;
  else
    return CMOR_NORMAL;
}

static int get_cmor_exit_control(struct kv **ht)
{
  if ( strcasecmp(get_val(ht, "exit_control", ""),
                  "CMOR_EXIT_ON_MAJOR") == 0 )
    return CMOR_EXIT_ON_MAJOR;
  else if ( strcasecmp(get_val(ht, "exit_control", ""),
                       "CMOR_EXIT_ON_WARNING")  == 0 )
    return CMOR_EXIT_ON_WARNING;
  else
    return CMOR_NORMAL;
}

static char *get_calendar_ptr(int calendar)
{
  char *calendar_ptr = Malloc(CMOR_MAX_STRING * sizeof(char));
  switch ( calendar )
    {
    case CALENDAR_STANDARD:
      strcpy(calendar_ptr, "gregorian");
    case CALENDAR_PROLEPTIC:
      strcpy(calendar_ptr, "proleptic_gregorian");
    case CALENDAR_360DAYS:
      strcpy(calendar_ptr, "360_day");
    case CALENDAR_365DAYS:
      strcpy(calendar_ptr, "noleap");
    case CALENDAR_366DAYS:
      strcpy(calendar_ptr, "all_leap");
    default:
      strcpy(calendar_ptr, "");
    }
  return calendar_ptr;
}

static int get_calendar_int(char *calendar)
{
  if ( strcmp(calendar, "gregorian") == 0 )
    return CALENDAR_STANDARD;
  else if ( strcmp(calendar, "proleptic_gregorian") == 0 )
    return CALENDAR_PROLEPTIC;
  else if ( strcmp(calendar, "360_day") == 0 )
    return CALENDAR_360DAYS;
  else if ( strcmp(calendar, "noleap") == 0 )
    return  CALENDAR_365DAYS;
  else if ( strcmp(calendar, "all_leap") == 0 )
    return  CALENDAR_366DAYS;
  else
    {
      cdoWarning("Calendar type %s is not supported by CMOR.\n");
      return 0;
    }
}

static void setup_dataset(struct kv **ht, int streamID)
{
  printf("*******Start to process cmor_setup and cmor_dataset.*******\n");
  int netcdf_file_action = get_netcdf_file_action(ht);
  int set_verbosity = get_cmor_verbosity(ht);
  int exit_control = get_cmor_exit_control(ht);

  char *logfile = get_val(ht, "logfile", NULL);
  int create_subdirectories = atoi(get_val(ht, "create_subdirectories", "0"));
  cmor_setup(get_val(ht, "inpath", "/usr/share/cmor/"),
             &netcdf_file_action,
             &set_verbosity,
             &exit_control,
             logfile,
             &create_subdirectories);

  int taxisID = vlistInqTaxis(streamInqVlist(streamID));
  char *attcalendar = get_val(ht, "calendar", "");
  char *calendar = get_calendar_ptr(taxisInqCalendar(taxisID));
  printf("Checking attribute 'calendar' from configuration.\n");
  if ( get_calendar_int(attcalendar) )
    check_compare_set(calendar, attcalendar, "calendar");
  else 
    {
      printf("Try to use Ifile calendar.\n");
      if ( !get_calendar_int(calendar) )
        cdoAbort("No valid configuration and no valid Ifile calendar found.");
      else
        hreplace(ht, "calendar", calendar);
    }
#if defined(CMOR_VERSION_MAJOR)
  int cmor_version_exists = 1;
  if ( CMOR_VERSION_MAJOR == 2 && CMOR_VERSION_MINOR == 9 )
    {
      double branch_time = atof(get_val(ht, "branch_time", "0.0"));
      cmor_dataset(get_val(ht, "outpath", "./"),
               get_val(ht, "experiment_id", ""),
               get_val(ht, "institution", ""),
               get_val(ht, "source", ""),
               calendar,
               atoi(get_val(ht, "realization", "")),
               get_val(ht, "contact", ""),
               get_val(ht, "history", ""),
               get_val(ht, "comment", ""),
               get_val(ht, "references", ""),
               atoi(get_val(ht, "leap_year", "")),
               atoi(get_val(ht, "leap_month", "")),
               NULL,
               get_val(ht, "model_id", ""),
               get_val(ht, "forcing", ""),
               atoi(get_val(ht, "initialization_method", "")),
               atoi(get_val(ht, "physics_version", "")),
               get_val(ht, "institute_id", ""),
               get_val(ht, "parent_experiment_id", ""),
               &branch_time,
               get_val(ht, "parent_experiment_rip", ""));
    }
  else
    cdoAbort("Cmor version %d.%d not yet enabled!\n", (int) CMOR_VERSION_MAJOR, (int) CMOR_VERSION_MINOR);
#endif
  if ( !cmor_version_exists )
    cdoAbort("It is not clear which CMOR version is installed since\nMakros CMOR_VERSION_MAJOR and CMOR_VERSION_MINOR are not available.\n");
  printf("*******Setup finished successfully.*******\n");
}


/*************************/
/* From Invert.c: */
/*************************/

static
void invertLatDataCmor(double *array1, double *array2, int gridID1)
{
  int nlat, nlon;
  int ilat;
  double **field1, **field2;

  nlon = gridInqXsize(gridID1);
  nlat = gridInqYsize(gridID1);

  if ( nlat > 0 )
    {
      field1 = (double **) Malloc(nlat*sizeof(double *));
      field2 = (double **) Malloc(nlat*sizeof(double *));
  
      for ( ilat = 0; ilat < nlat; ilat++ )
	{
	  field1[ilat] = array1 + ilat*nlon;
	  field2[ilat] = array2 + ilat*nlon;
	}

      for ( ilat = 0; ilat < nlat; ilat++ )
	memcpy(field2[nlat-ilat-1], field1[ilat], nlon*sizeof(double));
      
      if ( field1 ) Free(field1);
      if ( field2 ) Free(field2);
    }
  else
    {
      array2[0] = array1[0];
    }
}

/*************************/
/* Until here */
/*************************/


static char *get_time_units(int taxisID)
{
  char *units = Malloc ( CMOR_MAX_STRING * sizeof(char) );
  int timeunit = taxisInqTunit(taxisID);
  int year, month, day, hour, minute, second;
  cdiDecodeDate(taxisInqRdate(taxisID), &year, &month, &day);
  cdiDecodeTime(taxisInqRtime(taxisID), &hour, &minute, &second);
  if ( timeunit == TUNIT_QUARTER || timeunit == TUNIT_30MINUTES )
    timeunit = TUNIT_MINUTE;
  if ( timeunit == TUNIT_3HOURS || timeunit == TUNIT_6HOURS ||
       timeunit == TUNIT_12HOURS )
    timeunit = TUNIT_HOUR;

  sprintf(units, "%s since %d-%d-%d %02d:%02d:%02d", tunitNamePtr(timeunit),
          year, month, day, hour, minute, second);
  return units;
}

static int get_time_step_int(char *time_step)
{
  if ( strcmp(time_step, "hours") == 0  )
    return TUNIT_HOUR;
  else if ( strcmp(time_step, "days") == 0 ) 
    return TUNIT_DAY;
  else if ( strcmp(time_step, "months") ==  0 )
    return TUNIT_MONTH;
  else if ( strcmp(time_step, "years") == 0  )
    return TUNIT_YEAR;
  else
    {
      cdoWarning("Timeunit %s not yet implemented in cmor.\n", time_step);
      return 0;
    }
}

static int check_time_units(char *time_units)
{
  int attyear, attmonth, attday, atthour, attminute, attsecond;
  char time_step[CMOR_MAX_STRING];
  if ( sscanf(time_units, "%s since %d-%d-%d%*1s%02d:%02d:%02d%*1s",
                  time_step, &attyear, &attmonth, &attday, &atthour,
                  &attminute, &attsecond) != 7)
    {
      cdoWarning("Could not read all 7 required time unit values.");
      return 0;
    }
  if ( !get_time_step_int(time_step) )
    return 0;
  return 1;
}

static void get_taxis(char *req_time_units, char *attcalendar, int *sdate, int *stime, int *timeunit, int *calendar)
{
  int attyear, attmonth, attday, atthour, attminute, attsecond;
  char atttimeunit[CMOR_MAX_STRING];

  sscanf(req_time_units, "%s since %d-%d-%d%*1s%02d:%02d:%02d%*1s",
                  atttimeunit, &attyear, &attmonth, &attday, &atthour,
                  &attminute, &attsecond);
  *sdate = cdiEncodeDate(attyear, attmonth, attday);
  *stime = cdiEncodeTime(atthour, attminute, attsecond);
  *timeunit = get_time_step_int(atttimeunit);
  *calendar = get_calendar_int(attcalendar);
}

static char **get_requested_variables(struct kv **ht)
{
  char **name_list = NULL;
  char *select_vars = get_val(ht, "vars", NULL);

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

static void get_time_method(struct kv **ht, int vlistID, int varID, char *cmor_time_name)
{
  char *time_method = Malloc(8192 * sizeof(char));
  cdiInqAttTxt(vlistID, varID, "cell_methods", 8192, time_method);
  char *att_time_method = get_val(ht, "cell_methods", "");
  check_compare_set(time_method, att_time_method, "cell_methods");
  if ( time_method[0] == 'm' || strcmp(time_method, "time: mean") == 0 ) strcpy(cmor_time_name, "time \0");
  else if ( time_method[0] == 'p' ) strcpy(cmor_time_name, "time1\0");
  else if ( time_method[0] == 'c' ) strcpy(cmor_time_name, "time2\0");
  else if ( time_method[0] == 'n' ) strcpy(cmor_time_name, "none\0");
  else
    {
      cdoWarning("Found configuration time cell method '%s' is not valid. Check CF-conventions for allowed time cell methods.\nTime cell method is set to 'mean'. \n", time_method);
      strcpy(cmor_time_name, "time \0");
    }
  hreplace(ht, "time_axis", cmor_time_name); 
}

static void gen_bounds(int n, double *vals, double *bounds)
{
  for ( int i = 0; i < n-1; ++i )
    {
      bounds[2*i+1]   = 0.5*(vals[i] + vals[i+1]);
      bounds[2*(i+1)] = 0.5*(vals[i] + vals[i+1]);
    }

  bounds[0]     = 2*vals[0] - bounds[1];
  bounds[2*n-1] = 2*vals[n-1] - bounds[2*(n-1)];
}

static void get_zcell_bounds(int zaxisID, double *zcell_bounds, double *levels, int zsize)
{
  double *lbounds;
  lbounds = Malloc(zsize * sizeof(double));
  zaxisInqLbounds(zaxisID, lbounds);
  double *ubounds;
  ubounds = Malloc(zsize * sizeof(double));
  zaxisInqUbounds(zaxisID, ubounds);

  if ( !lbounds || !ubounds || ubounds[0] == ubounds[1] || lbounds[0] == lbounds[1] )
    gen_bounds(zsize, levels, zcell_bounds);
  else
    {
      if ( lbounds )
        zcell_bounds[0] = lbounds[0];
      else
        zcell_bounds[0] = 0; 
      for ( int i = 0; i < zsize-1; ++i )
        {
          zcell_bounds[2*i+1]   = ubounds[i];
          zcell_bounds[2*(i+1)] = lbounds[i+1];
        }
      if ( ubounds )
        zcell_bounds[2*zsize-1] = ubounds[zsize-1];
      else
        zcell_bounds[2*zsize-1] = levels[zsize-1] + ( levels[zsize-1] - zcell_bounds[2*zsize-2] );
    }
  Free(lbounds);
  Free(ubounds);
}

static void get_zhybrid(int zaxisID, char *varname, double *p0, double *alev_val, double *alev_bnds, double *b_val, double *b_bnds, double *ap_val, double *ap_bnds)
{
  int zsize = zaxisInqSize(zaxisID);

  int vctsize = zaxisInqVctSize(zaxisID);
  double *vct = Malloc(vctsize * sizeof(double) );
  zaxisInqVct(zaxisID, vct);
  if ( strcmp(varname, "ps") != 0 )
    {
      for ( int i = 0; i<(zsize+1); i++)
        {
          ap_bnds[i] = vct[i];
          b_bnds[i] = vct[zsize+1+i];
        }
      for ( int i = 0; i<zsize; i++)
        {
          ap_val[i] = (ap_bnds[i]+ ap_bnds[i+1]) / 2.0;
          b_val[i] = (b_bnds[i]+ b_bnds[i+1]) / 2.0;
          alev_val[i] = ap_val[i]/p0[0] + b_val[i];
          alev_bnds[i] = ap_bnds[i]/p0[0] + b_bnds[i];
        }
      alev_bnds[zsize] = ap_bnds[zsize]/p0[0] + b_bnds[zsize];
    }
  else
    {
      cdoAbort("Surface Pressure in hybrid z coordinate not yet enabled!");
    }
  Free(vct);
}

static void register_z_axis(struct kv **ht, int zaxisID, char *varname, int *axis_ids, int *zfactor_id)
{

  int zsize = zaxisInqSize(zaxisID);
  char units[CDI_MAX_NAME];
  double *levels;
  if ( zsize > 1)
    {
      levels = Malloc(zsize * sizeof(double));
      zaxisInqLevels(zaxisID, levels);
      double *zcell_bounds;
      zcell_bounds = Malloc( 2*zsize * sizeof(double) );
      get_zcell_bounds(zaxisID, zcell_bounds, levels, zsize);
      if ( zaxisInqType(zaxisID) == ZAXIS_PRESSURE )
        {
          cmor_axis(new_axis_id(axis_ids),
                        "plevs",
                        "Pa",
                        zsize,
                        (void *)levels,
                        'd', NULL, 0, NULL);
        }
      else if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID )
        {
          double *alev_val = Malloc(zsize * sizeof(double));
          double *alev_bnds = Malloc((zsize + 1) * sizeof(double));
          double *ap_val = Malloc(zsize * sizeof(double));
          double *ap_bnds = Malloc((zsize + 1) * sizeof(double));
          double *b_val = Malloc(zsize * sizeof(double));
          double *b_bnds = Malloc((zsize + 1) * sizeof(double));
          double *p0 = Malloc(sizeof(double));
          p0[0] = 101325.0;
          get_zhybrid(zaxisID, varname, p0, alev_val, alev_bnds, b_val, b_bnds, ap_val, ap_bnds);
/*cmor_zfactor (int *zfactor_id,int zaxis_id, char *zfactor_name, char *units, int ndims, int axis_ids[], char type, void *zfactor_values, void *zfactor_bounds)*/
          cmor_axis(new_axis_id(axis_ids),
                        "alternate_hybrid_sigma",
                        "",
                        zsize,
                        (void *)alev_val,
                        'd', alev_bnds,  1, NULL);
          int lev_id = axis_ids[count_axis_ids(axis_ids)-1];
          int lev_id_array[2];
          lev_id_array[0] = lev_id;
          cmor_zfactor(zfactor_id, lev_id, "p0", "Pa", 0, 0, 'd', (void *)p0, NULL);
          cmor_zfactor(zfactor_id, lev_id, "b", "", 1, &lev_id_array[0], 'd', (void *)b_val, (void *)b_bnds);
          cmor_zfactor(zfactor_id, lev_id, "ap", "Pa", 1, &lev_id_array[0], 'd', (void *)ap_val, (void *)ap_bnds);
          cmor_zfactor(zfactor_id, lev_id, "ps", "Pa", count_axis_ids(axis_ids)-1, axis_ids, 'd', NULL, NULL);  
          Free(alev_val);  
          Free(alev_bnds);  
          Free(ap_val);  
          Free(ap_bnds);  
          Free(b_val);  
          Free(b_bnds);  
        }
      else if ( zaxisInqType(zaxisID) == ZAXIS_DEPTH_BELOW_SEA )
        {
          zcell_bounds[0] = (double) 0;
          cmor_axis(new_axis_id(axis_ids),
                        "depth_coord",
                        "m",
                        zsize,
                        (void *)levels,
                        'd', zcell_bounds,  2, NULL);
        }
      else if ( zaxisInqType(zaxisID) == ZAXIS_DEPTH_BELOW_LAND )
        {
          zcell_bounds[0] = (double) 0;
          cmor_axis(new_axis_id(axis_ids),
                        "sdepth",
                        "m",
                        zsize,
                        (void *)levels,
                        'd', zcell_bounds, 2, NULL);
        }
      else if ( zaxisInqType(zaxisID) == ZAXIS_GENERIC || zaxisInqType(zaxisID) == ZAXIS_HEIGHT)
        cdoAbort("Z-axis type %d not yet enabled. \n", zaxisInqType(zaxisID));
      else
        cdoAbort("Invalid Z-axis type %d . \n", zaxisInqType(zaxisID));
      if (zcell_bounds) Free(zcell_bounds);
    }
  if ( zsize == 1 && get_val(ht, "szc", NULL) != NULL )
    {
      char *szc_name = get_val(ht, "szc", NULL);
      char *szc_value;
      strtok_r(szc_name, "_", &szc_value);
      levels = Malloc(sizeof(double));
      levels[0] = (double) atof(szc_value);
      printf("Attribute szc is found.\nScalar z coordinate name is: '%s'\nScalar z coordinate value is: '%f'\n", szc_name, levels[0]);
      cmor_axis(new_axis_id(axis_ids),
                      szc_name,
                      "m",
                      zsize,
                      (void *) levels,
                      'd', NULL, 0, NULL);
   }
}

static void register_character_dimension(int *axis_ids, char *filename)
{
  printf("The grid type is generic and a dimension 'basin' is found.\nTherefore, it is tried to read the character dimension.\n");
  int nc_file_id, nfiledims, nvars, ngatts, unlimdimid;
  nc_type xtypep;
  int varndims, varnattsp;
  int *vardimids;

  char *varname = malloc(36 * sizeof(char));
  char *dimname = malloc(36 * sizeof(char));

  size_t dimlength, dimstrlength;

  nc_open(filename, NC_NOWRITE, &nc_file_id);
  nc_inq(nc_file_id, &nfiledims, &nvars, &ngatts, &unlimdimid);
  vardimids = malloc(nfiledims * sizeof(int));
  void *final_chardim;
  for ( int i = 0; i < nvars; i++ )
    {
      nc_inq_var(nc_file_id, i, varname, &xtypep, &varndims, vardimids, &varnattsp);
      if ( strcmp(varname, "region") == 0 )
        {
          nc_inq_dim(nc_file_id, vardimids[1], dimname, &dimstrlength);
          nc_inq_dim(nc_file_id, vardimids[0], dimname, &dimlength);

          final_chardim = (void *)malloc(dimstrlength * dimlength *sizeof(char));
          nc_get_var(nc_file_id, i, final_chardim);
        }
    }
  nc_close(nc_file_id);
  cmor_axis(new_axis_id(axis_ids), dimname, "", dimlength, final_chardim, 'c',  NULL, dimstrlength, NULL); 
}

static void invert_ygriddes(struct kv **ht, int vlistID, int *gridID, int ylength, double *ycoord_vals, double *ycell_bounds, int *ynbounds)
{
  if ( ( ycoord_vals[0] - ycoord_vals[1] ) > 0 )
    {
      cdoWarning("Latitudes go north to south => reverting latitudes!\n");
      hinsert(ht, "invert_lat", "yes");
          
      int gridID2 = gridDuplicate(*gridID);
      double *yv2;
      double *yb2;

      yv2 = Malloc(ylength * sizeof(double));
      for ( int ilat = 0; ilat < ylength; ilat++ )
        yv2[ylength-ilat-1] = ycoord_vals[ilat];
      gridDefYvals(gridID2, yv2);

      yb2 = Malloc(2* ylength * sizeof(double));
      for ( int ilat = 0; ilat < ylength; ilat++ )
        {
          yb2[ylength*2-ilat*2-1] = ycell_bounds[ilat*2];
          yb2[ylength*2-ilat*2-2] = ycell_bounds[ilat*2+1];
        }
      gridDefYbounds(gridID2, yb2);
      vlistChangeGrid(vlistID, *gridID, gridID2);
      gridInqYvals(*gridID, ycoord_vals);
      *ynbounds = gridInqYbounds(*gridID, ycell_bounds);
    }
}

static void change_grid(char *grid_file, int *gridID, int vlistID)
{
  printf("You configured a grid_info file. For a successfull read, the file probably needs to have at least one variable with ID 0.\n");
  argument_t *fileargument = file_argument_new(grid_file);
  int streamID2 = streamOpenRead(fileargument); 
  int vlistID2 = streamInqVlist(streamID2);
  int gridID2 = vlistInqVarGrid(vlistID2, 0); 
  vlistChangeGrid(vlistID, *gridID, gridID2);
}

static void inquire_vals_and_bounds(int gridID, int *xnbounds, int *ynbounds, double *xcoord_vals, double *ycoord_vals, double *xcell_bounds, double *ycell_bounds)
{
  gridInqYvals(gridID, ycoord_vals);
  gridInqXvals(gridID, xcoord_vals);
  *xnbounds = gridInqXbounds(gridID, xcell_bounds);
  *ynbounds = gridInqYbounds(gridID, ycell_bounds);
}

static void get_cmor_table(struct kv **ht)
{
  int gridtable_id;
  char gridtable[CMOR_MAX_STRING];
  char *mip_table_dir = get_val(ht, "mip_table_dir", NULL);
  char *project_id = get_val(ht, "project_id", NULL);
  if ( mip_table_dir && project_id )
    {
      strcpy(gridtable, mip_table_dir);
      strcat(gridtable, "/");
      strcat(gridtable, project_id);
      strcat(gridtable, "_grids");
      if ( file_exist(gridtable, 1) )  
        {
          cmor_load_table(gridtable, &gridtable_id);
          cmor_set_table(gridtable_id);
        }
    }
  else
    {
      cdoAbort("A project grid table is required for this type of grid but not found in the mip table directory. Check attributes \n mip_table_dir \n and \n project_id !");
    } 
}

static void check_and_gen_bounds(int gridID, int nbounds, int length, double *coord_vals, double *cell_bounds, int x)
{
  if ( nbounds != 2 * length )
    {
      gen_bounds(length, coord_vals, cell_bounds);
      if ( x )
        {
          gridDefNvertex(gridID, 2);
          gridDefXbounds(gridID, cell_bounds);
        }
      else
        gridDefYbounds(gridID, cell_bounds);
    }
}

static void select_and_register_character_dimension(char *grid_file, int *axis_ids)
{
  char *ifile = cdoStreamName(0)->args;
  if ( strcmp(grid_file, "") == 0 )  
    register_character_dimension(axis_ids, ifile);
  else
    register_character_dimension(axis_ids, grid_file);
}

static void register_grid(struct kv **ht, int vlistID, int varID, int *axis_ids, int *grid_ids)
{
  int gridID = vlistInqVarGrid(vlistID, varID);

  char *grid_file = get_val(ht, "ginfo", "");
  if ( strcmp(grid_file, "") != 0 )
    change_grid(grid_file, &gridID, vlistID);

  int type = gridInqType(gridID);
  int ylength = gridInqYsize(gridID);
  int xlength = gridInqXsize(gridID);
  int totalsize = gridInqSize(gridID);

  double *xcoord_vals;
  double *ycoord_vals;
  double *xcell_bounds;
  double *ycell_bounds;
  int xnbounds;
  int ynbounds;

  /* Cmor call per Gridtype */
  if ( type == GRID_GAUSSIAN || type == GRID_LONLAT )
    {
      grid_ids[0] = 0;
      xcoord_vals = Malloc(xlength * sizeof(double));
      ycoord_vals = Malloc(ylength * sizeof(double));
      xcell_bounds = Malloc(2 * xlength * sizeof(double));
      ycell_bounds = Malloc(2 * ylength * sizeof(double));
      inquire_vals_and_bounds(gridID, &xnbounds, &ynbounds, xcoord_vals, ycoord_vals, xcell_bounds, ycell_bounds);

      check_and_gen_bounds(gridID, xnbounds, xlength, xcoord_vals, xcell_bounds, 1);
      check_and_gen_bounds(gridID, ynbounds, ylength, ycoord_vals, ycell_bounds, 0);

      invert_ygriddes(ht, vlistID, &gridID, ylength, ycoord_vals, ycell_bounds, &ynbounds);

      cmor_axis(new_axis_id(axis_ids),    "latitude",    "degrees_north",    ylength,    (void *)ycoord_vals,    'd',    (void *)ycell_bounds,    2,    NULL);
      cmor_axis(new_axis_id(axis_ids),    "longitude",    "degrees_east",    xlength,    (void *)xcoord_vals,    'd', 
   (void *)xcell_bounds,    2,    NULL);

      Free(xcell_bounds);
      Free(ycell_bounds);
    }
  else if ( type == GRID_CURVILINEAR || type == GRID_UNSTRUCTURED)
    {
      xcoord_vals = Malloc(totalsize * sizeof(double));
      ycoord_vals = Malloc(totalsize * sizeof(double));
/* maximal 4 gridbounds per gridcell permitted */
      xcell_bounds = Malloc(4 * totalsize * sizeof(double));
      ycell_bounds = Malloc(4 * totalsize * sizeof(double));
      inquire_vals_and_bounds(gridID, &xnbounds, &ynbounds, xcoord_vals, ycoord_vals, xcell_bounds, ycell_bounds);
      get_cmor_table(ht);
      int grid_axis[2];
      if ( type == GRID_CURVILINEAR )
        {
          double *xncoord_vals;
          double *yncoord_vals;
          xncoord_vals = Malloc(xlength * sizeof(double));
          yncoord_vals = Malloc(ylength * sizeof(double)); 
          for (int j=0; j<ylength; j++) 
            yncoord_vals[j]= (double) j;
          for (int j=0; j<xlength; j++)
            xncoord_vals[j]= (double) j;

          cmor_axis(&grid_axis[0], "j_index",    "1",    ylength,    (void *)yncoord_vals,
    'd', 0, 0, NULL);
          cmor_axis(&grid_axis[1], "i_index",    "1",    xlength,    (void *)xncoord_vals,    'd', 0, 0, NULL);
          cmor_grid(&grid_ids[0],    2,    grid_axis,    'd',    (void *)ycoord_vals,    (void *)xcoord_vals,    4,     (void *)ycell_bounds,    (void *)xcell_bounds);
        }
      /*else
        { 
          cmor_axis(&grid_axis[0],    "grid_longitude",   "degrees",    xlength,    (void *)xcoord_vals,    'd', 0, 0, NULL);
          cmor_axis(&grid_axis[1],    "grid_latitude",    "degrees",    ylength,    (void *)ycoord_vals,    'd', 0, 0, NULL);
          cmor_grid(&grid_ids[0],    2,    grid_axis,    'd',    (void *)ycoord_vals,    (void *)xcoord_vals,    2,     (void *)ycell_bounds,    (void *)xcell_bounds); 
        }*/
      Free(xcell_bounds);
      Free(ycell_bounds);
    }
  else if ( type == GRID_GENERIC && gridInqSize(gridID) > 1 )
    {
      grid_ids[0] = 0;
      xcoord_vals = Malloc(xlength * sizeof(double));
      gridInqXvals(gridID, xcoord_vals);
      ycoord_vals = Malloc(ylength * sizeof(double));
      gridInqYvals(gridID, ycoord_vals);
      char yname[CDI_MAX_NAME];
      cdiGridInqKeyStr(gridID, CDI_KEY_YDIMNAME, CDI_MAX_NAME, yname);
      if ( strcmp(yname, "basin") != 0 )
        {
          invert_ygriddes(ht, vlistID, &gridID, ylength, ycoord_vals, ycell_bounds, &ynbounds);
          ycell_bounds = Malloc(2 * ylength * sizeof(double));
          ynbounds = gridInqYbounds(gridID, ycell_bounds);
          check_and_gen_bounds(gridID, ynbounds, ylength, ycoord_vals, ycell_bounds, 0);
          cmor_axis(new_axis_id(axis_ids),    "latitude",    "degrees_north",    ylength,    (void *)ycoord_vals,    'd',    (void *)ycell_bounds,    2,    NULL);
          Free(ycell_bounds);
        }
      else
        select_and_register_character_dimension(grid_file, axis_ids);

      char xname[CDI_MAX_NAME];
      cdiGridInqKeyStr(gridID, CDI_KEY_XDIMNAME, CDI_MAX_NAME, xname);
      if ( strcmp(xname, "basin") != 0 )
        {
          xcell_bounds = Malloc(2 * xlength * sizeof(double));
          xnbounds = gridInqXbounds(gridID, xcell_bounds);
          check_and_gen_bounds(gridID, xnbounds, xlength, xcoord_vals, xcell_bounds, 1);
          cmor_axis(new_axis_id(axis_ids),    "longitude",    "degrees_east",    xlength,    (void *)xcoord_vals,    'd',    (void *)xcell_bounds,    2,    NULL);
          Free(xcell_bounds);
        }
      else
        select_and_register_character_dimension(grid_file, axis_ids);
    }
  else
    {
      grid_ids[0] = 0;
      cdoWarning("Either the grid type is unknown or a registration is not necessary. However, it is not registered by cdo.\n");
    }
  Free(xcoord_vals);
  Free(ycoord_vals);
}

static void register_variable(struct kv **ht, int vlistID, int varID, int *axis_ids,
                              struct mapping *var, int *grid_ids)
{
  char name[CDI_MAX_NAME];
  vlistInqVarName(vlistID, varID, name);
  char units[CDI_MAX_NAME];
  vlistInqVarUnits(vlistID, varID, units);
  char *attunits = get_val(ht, "units", "" );
  char *attname = get_val(ht, "out_name", "" );
  check_compare_set(name, attname, "out_name");
  check_compare_set(units, attunits, "units");
  char missing_value[sizeof(double)];
  double tolerance = 1e-4;
  size_t gridsize = vlistGridsizeMax(vlistID);
  int zsize = zaxisInqSize(vlistInqVarZaxis(vlistID, varID));
  var->cdi_varID = varID;
  var->help_var = 0;
  if ( vlistInqVarDatatype(vlistID, varID) == DATATYPE_FLT32 )
    {
      var->datatype = 'f';
      *(float *) missing_value = vlistInqVarMissval(vlistID, varID);
      var->data = Malloc(gridsize * zsize * sizeof(float));
    }
  else
    {
      var->datatype = 'd';
      *(double *) missing_value = vlistInqVarMissval(vlistID, varID);
      var->data = Malloc(gridsize * zsize * sizeof(double));
    }
  if ( grid_ids[0] != 0 )
    {
      int *tmp_id = new_axis_id(axis_ids);
      *tmp_id = grid_ids[0];
      cmor_variable(&var->cmor_varID,
            name,units,(count_axis_ids(axis_ids)), axis_ids, var->datatype,
            (void *) missing_value, &tolerance, get_val(ht, "positive", NULL),
                        NULL, NULL, NULL);
    }
  else
    {
      cmor_variable(&var->cmor_varID,
           name, units, count_axis_ids(axis_ids),  axis_ids,   var->datatype,
          (void *) missing_value, &tolerance, get_val(ht, "positive", NULL),
                        NULL, NULL, NULL);
    }
}

static void register_all_dimensions(struct kv **ht, int streamID,
                             struct mapping vars[], int table_id, int *zfactor_id)
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

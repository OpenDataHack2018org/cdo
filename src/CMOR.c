#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

#if defined(HAVE_LIBCMOR)
#include <ctype.h>
#include <unistd.h>
#include "uthash.h"
#include "util.h"
#include "cmor.h"
#include "netcdf.h"
#include "pmlist.h"

#define CMOR_UNDEFID (CMOR_MAX_AXES + 1)

/* */
/* Read Mapping Table */
/* */
int stringToParam(const char *paramstr);

list_t *pmlist_search_kvlist_ventry(list_t *pml, const char *key, const char *value, int nentry, const char **entry);
list_t *pmlist_get_kvlist_ventry(list_t *pml, int nentry, const char **entry);

static
char *readLineFromBuffer(char *buffer, size_t *buffersize, char *line, size_t len)
{
  int ichar;
  size_t ipos = 0;
  while ( *buffersize )
    {
      ichar = *buffer;
      (*buffersize)--;
      buffer++;
      if ( ichar == '\r' )
        {
          if ( *buffersize )
            {
              ichar = *buffer;
              if ( ichar == '\n' )
                {
                  (*buffersize)--;
                  buffer++;
                }
            }
          break;
        }
      if ( ichar == '\n' ) break;
      line[ipos++] = ichar;
      if ( ipos >= len )
        {
          fprintf(stderr, "readLineFromBuffer: end of line not found (maxlen = %ld)!\n", len);
          break;
        }
    }
  line[ipos] = 0;
  if ( *buffersize == 0 && ipos == 0 ) buffer = NULL;
  return buffer;
}

static
char *skipSeparator(char *pline)
{
  while ( isspace((int) *pline) ) pline++;
  if ( *pline == '=' || *pline == ':' ) pline++;
  while ( isspace((int) *pline) ) pline++;

  return pline;
}

static
char *getElementName(char *pline, char *name)
{
  while ( isspace((int) *pline) ) pline++;
  size_t len = strlen(pline);
  size_t pos = 0;
  while ( pos < len && !isspace((int) *(pline+pos)) && *(pline+pos) != '=' && *(pline+pos) != ':' )
    {
      name[pos] = tolower(*(pline+pos));
      pos++;
    }
  name[pos] = 0;

  pline += pos;
  return pline;
}

static void copy_value(char *value, char **values, int *nvalues)
{
  values[*nvalues] = malloc((strlen(value) + 1) * sizeof(values[*nvalues]));
  strncpy(values[*nvalues], value, strlen(value)+1);
  (*nvalues)++;
}

static
char *getElementValues(char *pline, char **values, int *nvalues)
{
  while ( isspace((int) *pline) ) pline++;
  size_t len = strlen(pline);
  *nvalues = 0;
  int i = 0;
  while ( i < len && len )
    {
      if ( *(pline+i) == ',')
        {
          copy_value(pline, values, nvalues);
          *(values[*nvalues-1]+i) = 0;

          i++;
          pline+=i;
          len-=i;
          i=0;
        }
      else if ( *(pline+i) == '"' )
        {
          i++;
          while ( *(pline+i) != '"' )
            {
              i++;
              if ( *(pline+i) == 0 )
                cdoAbort("Found a start quote sign for a value but no end quote sign.\n");
            }
          i++;
        }
      else if ( isspace((int) *(pline+i)) )
        {
          copy_value(pline, values, nvalues);
          if ( *(values[*nvalues-1]+i-1) == '"' || *(values[*nvalues-1]+i-1) == '\'' )
            {
              values[*nvalues-1]++;
              *(values[*nvalues-1]+i-2) = 0;
            }
          else
            *(values[*nvalues-1]+i) = 0;
          i++; 
          pline+=i;
          break;          
        }
      else if ( *(pline+i) == '=' || *(pline+i) == ':' )
        cdoAbort("Found unexpected separator sign in value: '%c'.", *(pline+i) );
      else
        i++;
    }
  if ( i == len && len )
    {
      copy_value(pline, values, nvalues);
      if ( *(values[*nvalues-1]+i-1) == '"' )
        {
          values[*nvalues-1]++;
          *(values[*nvalues-1]+i-2) = 0;
        }
      else
        *(values[*nvalues-1]+i) = 0;
      *pline = 0;
    }
  return pline;
}

static void parse_line_to_list(list_t *list, char *pline, char *kvlname, int checkpml, int lowprior)
{
  char name[256];
  int i = 0, nvalues;
  list_t *kvl = NULL;
  if ( checkpml )
    {
      kvl = list_new(sizeof(keyValues_t *), free_keyval, kvlname);
      list_append(list, &kvl);
    }
  while ( *pline != 0 )
    {
      char **values = malloc( 5 * sizeof(char *) );
      pline = getElementName(pline, name);
      if ( *pline == 0 )
        {
          if ( cdoVerbose )
            printf("Could not find values for: '%s'\n", name);
          break;
        }
      pline = skipSeparator(pline);
      if ( *pline == 0 )
        {
          if ( cdoVerbose )
            printf("Could not find values for: '%s'\n", name);
          break;
        }
      pline = getElementValues(pline, values, &nvalues);
      if ( checkpml )
        kvlist_append(kvl, name, (const char **)values, nvalues);
      else
        {
          if ( lowprior )
            {
              keyValues_t *already = kvlist_search(list, name);
              if ( already )
                continue;
            }
          kvlist_append(list, name, (const char **)values, nvalues);
        }
      while ( isspace((int) *pline) ) pline++;
      if ( *pline == '/' )
        *pline = 0;
      free(values);     
    }
}

void parse_buffer_to_list(list_t *list, size_t buffersize, char *buffer, int checkpml, int lowprior)
{
  char line[4096];
  char name[256];
  char *pline;
  char *listkeys[] = {"axis_entry:", "variable_entry:", "&parameter", NULL};
  int linenumber = 0;
  int listtype = 0;

  while ( (buffer = readLineFromBuffer(buffer, &buffersize, line, sizeof(line))) )
    {
      linenumber++;
      pline = line;
      while ( isspace((int) *pline) ) pline++;
      if ( *pline == '#' || *pline == '!' || *pline == '\0' ) continue;
      //  len = (int) strlen(pline);
      if ( listtype == 0 && *pline == '&' ) listtype = 1;
/* MAXNVALUES*/
      int i = 0;
      while ( listkeys[i] )
        {
          if ( strlen(pline) > strlen(listkeys[i]) )
            if ( strncmp(pline, listkeys[i], strlen(listkeys[i])) == 0 )
              {
	        pline += strlen(listkeys[i]);
 	        listtype = 2;
                break;
	      }
          i++;
        }
      if ( listtype )
        parse_line_to_list(list, pline, listkeys[i], checkpml, lowprior);
      else
        parse_line_to_list(list, pline, "keyvals", checkpml, lowprior);
    }
}


static void kv_insert_a_val(list_t *kvl, const char *key, char *value, int replace)
{
  if ( key )
    {
      keyValues_t *kv = kvlist_search(kvl, key);
      if ( !kv )
        {
          const char *apvalue[] = {value};
          kvlist_append(kvl, key, apvalue, 1);
        }
      else if ( replace )
        {
          free(kv->values[0]);
          kv->values[0] = strdup((const char *) value);
        }
    }
}

static char *kv_get_a_val(list_t *kvl, const char *key, char *replacer)
{
  keyValues_t *kv = kvlist_search(kvl, key);
  if ( kv )
    {
      char *outval = strdup(kv->values[0]);
      return outval;
    }
  return replacer;
}

list_t *cdo_parse_cmor_file(const char *filename)
{
  assert(filename != NULL);

  size_t filesize = fileSize(filename);
  if ( filesize == 0 )
    {
      fprintf(stderr, "Empty table file: %s\n", filename);
      return NULL;
    }

  FILE *fp = fopen(filename, "r");
  if ( fp == NULL )
    {
      fprintf(stderr, "Open failed on %s: %s\n", filename, strerror(errno));
      return NULL;
    }

  char *buffer = (char*) Malloc(filesize);
  size_t nitems = fread(buffer, 1, filesize, fp);

  fclose(fp);

  if ( nitems != filesize )
    {
      fprintf(stderr, "Read failed on %s!\n", filename);
      return NULL;
    }
 
  list_t *pml = list_new(sizeof(list_t *), free_kvlist, filename);

/*  if ( buffer[0] == '{' )
    parse_json_buffer_to_pml(pml, filesize, buffer);
  else */
  parse_buffer_to_list(pml, filesize, buffer, 1, 0);

  Free(buffer);
  
  return pml;
}

static void map_it(list_t *kvl, int vlistID, int varID)
{
  for ( listNode_t *kvnode = kvl->head; kvnode; kvnode = kvnode->next )
    {
      keyValues_t *kv = *(keyValues_t **)kvnode->data;
      const char *key = kv->key;
      const char *value = (kv->nvalues == 1) ? kv->values[0] : NULL;
/*      printf("'%s' = '%s'\n", key, value); */
      if ( !value ) continue;
/* Not necessary because out_name is what we use for renaming :
      else if ( STR_IS_EQ(key, "name")          ) vlistDefVarName(vlistID, varID, parameter2word(value));
*/
      else if ( STR_IS_EQ(key, "out_name")      )
        {
          char name[CDI_MAX_NAME];
          vlistInqVarName(vlistID, varID, name);
          if ( name[0] != 0 )
            cdiDefAttTxt(vlistID, varID, "original_name", (int) strlen(name), parameter2word(name));
          vlistDefVarName(vlistID, varID, parameter2word(value));
         }
      else if ( STR_IS_EQ(key, "units")         ) vlistDefVarUnits(vlistID, varID, value);
      else if ( STR_IS_EQ(key, "cell_methods")  ) cdiDefAttTxt(vlistID, varID, "cell_methods", (int) strlen(value), value);
      else if ( STR_IS_EQ(key, "standard_name") ) vlistDefVarStdname(vlistID, varID, value);
      else if ( STR_IS_EQ(key, "factor")        ) {}
      else if ( STR_IS_EQ(key, "delete")        ) {}
/* Not in mapping table right now: */
      else if ( STR_IS_EQ(key, "long_name")     ) vlistDefVarLongname(vlistID, varID, value);
      else if ( STR_IS_EQ(key, "param")         ) vlistDefVarParam(vlistID, varID, stringToParam(parameter2word(value)));
      else if ( STR_IS_EQ(key, "out_param")     ) vlistDefVarParam(vlistID, varID, stringToParam(parameter2word(value)));
      else if ( STR_IS_EQ(key, "code")          ) {}
              // else if ( STR_IS_EQ(key, "code")          ) vlistDefVarParam(vlistID, varID, cdiEncodeParam(parameter2int(value), ptab, 255));
              // else if ( STR_IS_EQ(key, "out_code")      ) vlistDefVarParam(vlistID, varID, cdiEncodeParam(parameter2int(value), ptab, 255));
      else if ( STR_IS_EQ(key, "comment")       ) cdiDefAttTxt(vlistID, varID, "comment", (int) strlen(value), value);
      else if ( STR_IS_EQ(key, "p")       )
        {
          if ( !isspace(value[0]) )
            cdiDefAttTxt(vlistID, varID, "positive", (int) strlen(value), value);
        }
      else if ( STR_IS_EQ(key, "cell_measures") ) cdiDefAttTxt(vlistID, varID, "cell_measures", (int) strlen(value), value);
      else if ( STR_IS_EQ(key, "convert")       ) {} 
      else if ( STR_IS_EQ(key, "missval")       )   {}  
      else if ( STR_IS_EQ(key, "valid_min")     ){}
      else if ( STR_IS_EQ(key, "valid_max")     ){}
      else if ( STR_IS_EQ(key, "ok_min_mean_abs") ) {}
      else if ( STR_IS_EQ(key, "ok_max_mean_abs") ){}
      else if ( STR_IS_EQ(key, "datatype") || STR_IS_EQ(key, "type") ) {}
      else
        {
          if ( cdoVerbose ) cdoPrint("Key >%s< not supported!", key);
        }
    }
}

static void get_stringcode(int vlistID, int varID, char *varcodestring)
{
  int varcode;
  varcode = vlistInqVarCode(vlistID, varID);
  sprintf(varcodestring, "%03d", varcode);
}

static void get_ifilevalue(char *ifilevalue, const char *key, int vlistID, int varID)
{
  if ( strcmp(key, "name") == 0 )
    vlistInqVarName(vlistID, varID, ifilevalue);
  else if ( strcmp(key, "code") == 0 )
    {
      char varcodestring[CDI_MAX_NAME];
      get_stringcode(vlistID, varID, varcodestring);
      strcpy(ifilevalue, varcodestring);
    }
}

static int map_via_key(list_t *pml, int vlistID, int varID, int nventry, const char **ventry, const char *key)
{
  char ifilevalue[CDI_MAX_NAME];
  get_ifilevalue(ifilevalue, key, vlistID, varID);

  list_t *kvl = pmlist_search_kvlist_ventry(pml, key, ifilevalue, nventry, ventry);
  if ( ifilevalue[0] )
    {
      if ( kvl )
        {
          map_it(kvl, vlistID, varID);
          return 1;
        }
      cdoWarning("Variable with name '%s' could not be mapped via '%s' because no corresponding key '%s' was found in mapping table file.\n", ifilevalue, key, key);
      return 0;
    }
  else
    {
      cdoWarning("Variable with varID '%d' could not be mapped via '%s' because it does not possess a '%s' in Ifile.", varID, key, key);
      return 0;
    }
}

static int map_via_vars_and_key(list_t *kvl_oname, int vlistID, int nvars, char *key)
{
  keyValues_t *kv = kvlist_search(kvl_oname, key);
  if ( kv )
    {
      for ( int varID = 0; varID < nvars; varID++ )
        {
          char ifilevalue[CDI_MAX_NAME];
          get_ifilevalue(ifilevalue, key, vlistID, varID);

          if ( ifilevalue[0] && strcmp(ifilevalue, kv->values[0]) == 0 )
            {
              map_it(kvl_oname, vlistID, varID);
              return 1;
            }
        }
      cdoWarning("The requested variable (via attribute vars) could not be mapped via key '%s' because no Ifile variable '%s' equals '%s'.", key, key, kv->values[0]);
    }
  else
    cdoWarning("The requested variable (via attribute vars) could not be mapped via key '%s' because it possesses no corresponding key '%s' in mapping file.", key, key);
  return 0;
}

static void map_via_vars(list_t *pml, const char **vars, int vlistID, int nvars, int nventry, const char **ventry)
{
  int j = 0;
  while ( vars[j] )
    {
      if ( cdoVerbose )
        printf("*******Try to map requested variable: '%s'********\n", vars[j]);
      list_t *kvl_oname = pmlist_search_kvlist_ventry(pml, "out_name", vars[j], nventry, ventry);
      if ( kvl_oname )
        {
          if ( map_via_vars_and_key(kvl_oname, vlistID, nvars, "name") )
            j++;
          else if ( map_via_vars_and_key(kvl_oname, vlistID, nvars, "code") )
            j++;
          else
            {
              cdoWarning("Could not map requested variable '%s'\n", vars[j]);
              j++;
            }
        }
      else
        {
          cdoWarning("Requested variable '%s' (via attribute vars) is not found in row 'out_name' of mapping table.'\n", vars[j]);
          j++;
        }
    }
}

static
void apply_mapping_table(const char *filename, int nvars, int vlistID, const char **request)
{
  if ( cdoVerbose )
    printf("*******Try to apply mapping table: '%s'*******\n", filename);
  const char *ventry[] = {"&parameter"};
  int nventry = (int) sizeof(ventry)/sizeof(ventry[0]);

  list_t *pml = cdo_parse_cmor_file(filename);
  if ( pml == NULL ) return;
  if ( request )
    map_via_vars(pml, request, vlistID, nvars, nventry, ventry);
  else
    {
      for ( int varID = 0; varID < nvars; varID++ )
        {
          if ( map_via_key(pml, vlistID, varID, nventry, ventry, "name") )
            continue;
          if ( map_via_key(pml, vlistID, varID, nventry, ventry, "code") )
            continue;
          cdoWarning("Could not map variable with id '%d'.", varID);
        }
    }
  list_destroy(pml);
}
/* */
/*... until here */
/* */

struct mapping
{
  int help_var;
  int cdi_varID;
  int cmor_varID;
  int zfactor_id;
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

static int file_exist(const char *tfilename, int force)
{
  assert(tfilename != NULL);
  size_t filesize = fileSize(tfilename);
  if ( filesize == 0 && force)
    {
      fprintf(stderr, "Empty file: %s\n", tfilename);
      return 0;
    }
  else if ( filesize == 0 && !force )
    {
      cdoWarning("cannot open '%s'", tfilename);
      return 0;
    }
  if ( strstr(tfilename, ".nc") || strstr(tfilename, ".grb") )
    return 1;
  FILE *fp = fopen(tfilename, "r");
  if ( fp == NULL && force )
    {
      fprintf(stderr, "Open failed on %s: %s\n", tfilename, strerror(errno));
      return 0;
    }
  else if ( fp == NULL && !force )
    {
      cdoWarning("cannot open '%s'", tfilename);
      return 0;
    }

  fclose(fp);
  return 1;
}
  
static int parse_kv_file(list_t *kvl, const char *filename)
{
  if ( !file_exist(filename, 1) )
    cdoAbort("Configuration failed.\n");

  FILE *fp = fopen(filename, "r");
  size_t filesize = fileSize(filename);
  char *buffer = (char*) Malloc(filesize);
  size_t nitems = fread(buffer, 1, filesize, fp);
  fclose(fp);

  parse_buffer_to_list(kvl, filesize, buffer, 0, 1);
  Free(buffer);
  return 0;
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

static int check_attr(list_t *kvl)
{
  const char *reqAtt[] = {"institute_id", "institution", "contact", "model_id", "source",
            "experiment_id", "req_time_units", "project_id", NULL};
  const char *reqAttCMIP5[] = {"product", "member", NULL};
  const char *reqAttCORDEX[] = {"product", "member", "cordex_domain", "driving_model_id", NULL};
/* In all Projects needed Attributes are tested first */
  int i = 0;

  while ( reqAtt[i] != NULL )
    {
      keyValues_t *kv_reqatt = kvlist_search(kvl, reqAtt[i]);
      
      if ( !kv_reqatt || strcmp(kv_reqatt->values[0], "notSet") == 0 )
        cdoAbort("Attribute '%s' is required. Either it is missing or notSet.", reqAtt[i]);
      i++;
    }
/* Set default attributes */
  keyValues_t *kv = kvlist_search(kvl, "references");

  if ( !kv || strcmp(kv->key, "notSet") == 0 )
    {
      keyValues_t *kv_model_id = kvlist_search(kvl, "model_id");
      char *references = (char *) Malloc(strlen(kv_model_id->values[0]) + 28);
      strcpy(references, "No references available for ");
      strcat(references, kv_model_id->values[0]);
      cdoWarning("Attribute 'references' is set to '%s' ", references);
      kv_insert_a_val(kvl, "references", references, 1);
      Free(references);
    }

/* Special check for CMIP5 or CORDEX projects */
  i=0;
  char *kv_project_id = kv_get_a_val(kvl, "project_id", "");
  if ( strcmp(kv_project_id, "CMIP6") == 0 )
    cdoAbort("Not yet possible to create data for project CMIP6 since cmor version 2.9 is used in this operator.\n");
  if ( strcmp(kv_project_id, "CMIP5") == 0 )
    {
      if ( cdoVerbose )
        printf("Since the project id is %s further attributes are tested. \n", kv_project_id);
      while ( reqAttCMIP5[i] != NULL )
        {
          keyValues_t *kv_reqattCMIP5 = kvlist_search(kvl, reqAttCMIP5[i]);
          if ( !kv_reqattCMIP5 || strcmp(kv_reqattCMIP5->values[0], "notSet") == 0 )
            cdoAbort("Attribute '%s' is required. Either it is missing or notSet", reqAttCMIP5[i]);
          if ( cdoVerbose )
            printf("Attribute '%s' is %s \n", reqAttCMIP5[i], kv_reqattCMIP5->values[0]);
          i++;
        }
    }
  else if (strcmp(kv_project_id, "CORDEX") == 0 )
    {
      if ( cdoVerbose )
        printf("Since the project id is %s further attributes are tested", kv_project_id);
      i=0;
      while ( reqAttCORDEX[i] != NULL )
        {
          keyValues_t *kv_reqattCORDEX = kvlist_search(kvl, reqAttCORDEX[i]);
          if ( !kv_reqattCORDEX || strcmp(kv_reqattCORDEX->values[0], "notSet") == 0 )
            cdoAbort("Attribute '%s' is required. Either it is missing or notSet", reqAttCORDEX[i]);
          if ( cdoVerbose )
            printf("Attribute '%s' is %s \n", reqAttCORDEX[i], kv_reqattCORDEX->values[0]);
          i++;
        }
    }

/* Check for special attributes */
  keyValues_t *kv_pos = kvlist_search(kvl, "p");
  if ( kv_pos )
    if ( ( strcmp(kv_pos->values[0], "notSet") == 0 ) || ( strcmp(kv_pos->values[0], "") != 0 && strcmp(kv_pos->values[0], "u") != 0 && strcmp(kv_pos->values[0], "d") != 0 ) )
      {
        cdoWarning("Invalid value '%s' is set for attribute 'positive'. The default (blank) is used.\n", kv_pos->values[0]);
        free(kv_pos->values[0]);
        kv_pos->values[0] = strdup(" ");
      }
  return 1;
} 

static int check_mem(list_t *kvl)
{
  char *kv_member = kv_get_a_val(kvl, "member", "");
  char *kv_project_id = kv_get_a_val(kvl, "project_id", "");
  char *ripchar[] = {"realization", "initialization", "physics_version"};
  char *crealiz, *cinitial, *cphysics;
  char workchar[CMOR_MAX_STRING]; 
  int realization, initialization_method, physics_version;
  int ipos=0, ppos=0;

/* Test for the right member, else abort or warn */ 
  if ( strlen(kv_member) >= 6 && kv_member[0] == 'r' )
    {
      crealiz = cinitial = cphysics = (char *) Malloc(strlen(kv_member));
      strcpy(crealiz, &kv_member[1]);
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
      char *ripvaluechar[] = {crealiz, cinitial, cphysics};
      for ( int i = 0; i < 3; i++ )
        kv_insert_a_val(kvl, (const char *)ripchar[i], ripvaluechar[i], 1);
      return 1;
    }
  else if ( strcmp(kv_member, "notSet") == 0 )
    {
      cdoWarning("The member has no RIP format! We set \n Attribute realization=-1 \n Attribute initialization_method=-1 \n Attribute physics_version=-1 \n");
      for ( int i = 0; i < 3; i++ )   
        kv_insert_a_val(kvl, (const char *)ripchar[i], "-1", 1);
    }
/* Now abort or warn */ 
  if (strcmp(kv_project_id, "CMIP5") == 0 || strcmp(kv_project_id, "CORDEX") == 0)
    cdoAbort("The member has no RIP format (at least 6 characters and in RIP order)! Found for \n member: %s. This is interpreted as \n Realization: %s \n Initialization: %s \n Physics: %s \n   But three Integers are needed", kv_member, crealiz, cinitial, cphysics);

  return 0;
} 


/*
static void dump_global_attributes(list_t *pml, int streamID)
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

static void dump_special_attributes(list_t *kvl, int streamID)
{
  int vlistID = streamInqVlist(streamID);
  int fileID = pstreamFileID(streamID);
  size_t old_historysize;
  char *new_history = kv_get_a_val(kvl, "history", "");
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
      kv_insert_a_val(kvl, "history", history, 1);
      Free(history);
    }
}

static void read_config_files(list_t *kvl)
{
  /* Files from info key in command line. */
  keyValues_t *info = kvlist_search(kvl, "info");
  int i = 0;
  if ( info )
    while ( i < info->nvalues )
      {
        parse_kv_file(kvl, info->values[i]);
        i++;
      }

  /* Config file in user's $HOME directory. */
  char *home = getenv("HOME");
  const char *dotconfig = ".cdocmorinfo";
  char *workfile = Malloc(strlen(home) + strlen(dotconfig) + 2);
  sprintf(workfile, "%s/%s", home, dotconfig);
  parse_kv_file(kvl, workfile);
  Free(workfile);
  
  if ( i == 0 )
    {
      keyValues_t *info2 = kvlist_search(kvl, "info");
      if ( info2 )
        while ( i < info2->nvalues )
          {
            parse_kv_file(kvl, info2->values[i]);
            i++;
          }
    }
}

static int in_list(char **list, const char *needle)
{
  while ( *list )
    if ( strcmp(*list++, needle) == 0 )
      return 1;
  return 0;
}

static int get_netcdf_file_action(list_t *kvl)
{
  char *chunk = kv_get_a_val(kvl, "oflag", "");
  if ( strcmp(chunk, "r") == 0 )
    return CMOR_REPLACE;
  else if ( strcmp(chunk, "a") == 0 )
    return CMOR_APPEND;
  else if ( strcmp(chunk, "p") == 0 )
    return CMOR_APPEND;
  else
    {
      cdoWarning("No valid CMOR output mode! \nAttribute oflag is '%s', but valid are 'a' for append, 'r' for replace or 'p' for preserve.\nCMOR output mode is set to: replace.\n", chunk);
      return CMOR_REPLACE;
    }
}

static int get_cmor_verbosity(list_t *kvl)
{
  char *verbos = kv_get_a_val(kvl, "set_verbosity", "");
  if ( strcmp(verbos, "CMOR_QUIET") == 0 )
    return CMOR_QUIET;
  else
    return CMOR_NORMAL;
}

static int get_cmor_exit_control(list_t *kvl)
{
  char *exit = kv_get_a_val(kvl, "exit_control", "");
  if ( strcasecmp(exit,"CMOR_EXIT_ON_MAJOR") == 0 )
    return CMOR_EXIT_ON_MAJOR;
  else if ( strcasecmp(exit,"CMOR_EXIT_ON_WARNING")  == 0 )
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

static char *get_txtatt(int vlistID, int varID, char *key)
{
  int natts;
  cdiInqNatts(vlistID, varID, &natts);
  for ( int i = 0; i < natts; i++ )
    {
      char name[CDI_MAX_NAME];
      char *txtatt = NULL;
      char buffer[8];
      int type, len;
      cdiInqAtt(vlistID, varID, i, name, &type, &len);
      if ( strcmp(name, key) == 0 )
        {
          txtatt = Malloc((len + 1) * sizeof(char));
          cdiInqAttTxt(vlistID, varID, name, len, txtatt);
          txtatt[len] = '\0';
          return txtatt;
        }
    }
  return "";
}

static void setup_dataset(list_t *kvl, int streamID)
{
  if ( cdoVerbose )
    printf("*******Start to process cmor_setup and cmor_dataset.*******\n");
  int netcdf_file_action = get_netcdf_file_action(kvl);
  int set_verbosity = get_cmor_verbosity(kvl);
  int exit_control = get_cmor_exit_control(kvl);
  int creat_subs = atol(kv_get_a_val(kvl, "drs", "0"));

  int vlistID = streamInqVlist(streamID);

  cmor_setup(kv_get_a_val(kvl, "inpath", "/usr/share/cmor/"),
             &netcdf_file_action,
             &set_verbosity,
             &exit_control,
             kv_get_a_val(kvl, "logfile", NULL),
             &creat_subs);

  int taxisID = vlistInqTaxis(streamInqVlist(streamID));
  char *attcomment = kv_get_a_val(kvl, "comment", "");
  char *comment = get_txtatt(vlistID, CDI_GLOBAL, "comment");
  check_compare_set(comment, attcomment, "comment");
  
  char *attcalendar = kv_get_a_val(kvl, "calendar", "");
  char *calendar = get_calendar_ptr(taxisInqCalendar(taxisID));
  if ( cdoVerbose )
    printf("Checking attribute 'calendar' from configuration.\n");
  if ( get_calendar_int(attcalendar) )
    check_compare_set(calendar, attcalendar, "calendar");
  else 
    {
      if ( cdoVerbose )
        printf("Try to use Ifile calendar.\n");
      if ( !get_calendar_int(calendar) )
        cdoAbort("No valid configuration and no valid Ifile calendar found.");
      else
        kv_insert_a_val(kvl, "calendar", calendar, 1);
    }
#if defined(CMOR_VERSION_MAJOR)
  int cmor_version_exists = 1;
  if ( CMOR_VERSION_MAJOR == 2 && CMOR_VERSION_MINOR == 9 )
    {
      double branch_time = atof(kv_get_a_val(kvl, "branch_time", "0.0"));
      cmor_dataset(kv_get_a_val(kvl, "odir", "./"),
               kv_get_a_val(kvl, "experiment_id", ""),
               kv_get_a_val(kvl, "institution", ""),
               kv_get_a_val(kvl, "source", ""),
               calendar,
               atoi(kv_get_a_val(kvl, "realization", "")),
               kv_get_a_val(kvl, "contact", ""),
               kv_get_a_val(kvl, "history", ""),
               kv_get_a_val(kvl, "vcomm", ""),
               kv_get_a_val(kvl, "references", ""),
               atoi(kv_get_a_val(kvl, "leap_year", "")),
               atoi(kv_get_a_val(kvl, "leap_month", "")),
               NULL,
               kv_get_a_val(kvl, "model_id", ""),
               kv_get_a_val(kvl, "forcing", ""),
               atoi(kv_get_a_val(kvl, "initialization_method", "")),
               atoi(kv_get_a_val(kvl, "physics_version", "")),
               kv_get_a_val(kvl, "institute_id", ""),
               kv_get_a_val(kvl, "parent_experiment_id", ""),
               &branch_time,
               kv_get_a_val(kvl, "parent_experiment_rip", ""));
    }
  else
    cdoAbort("Cmor version %d.%d not yet enabled!\n", (int) CMOR_VERSION_MAJOR, (int) CMOR_VERSION_MINOR);
#endif
  if ( !cmor_version_exists )
    cdoAbort("It is not clear which CMOR version is installed since\nMakros CMOR_VERSION_MAJOR and CMOR_VERSION_MINOR are not available.\n");
  if ( cdoVerbose )
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
      
      Free(field1);
      Free(field2);
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

static char **get_requested_variables(list_t *kvl)
{
  keyValues_t *kv = kvlist_search(kvl, "vars");
  char **req_vars = NULL;
  if ( kv )
    {
      req_vars = Malloc((kv->nvalues + 1) * sizeof(char *));
      for ( int i = 0; i < kv->nvalues; i++ )
        req_vars[i] = strdup(kv->values[i]);
    }
  if ( req_vars )
    req_vars[kv->nvalues] = NULL;
  return req_vars;
}

static void get_time_method(list_t *kvl, int vlistID, int varID, char *cmor_time_name)
{
  char *mipfreq = kv_get_a_val(kvl, "miptab_freq", "");
  char *project = kv_get_a_val(kvl, "project_id", "");
  if (strcmp(project, "CMIP5") == 0 )
    {
      if ( strcmp(mipfreq, "Oclim") == 0 )
        strcpy(cmor_time_name, "time2");
      if ( strcmp(mipfreq, "Oyr") == 0 )
        strcpy(cmor_time_name, "time");
      if ( strcmp(mipfreq, "day") == 0 )
        strcpy(cmor_time_name, "time");
      if ( strcmp(mipfreq, "6hrPlev") == 0 )
        strcpy(cmor_time_name, "time1");
      if ( strcmp(mipfreq, "6hrLev") == 0 )
      strcpy(cmor_time_name, "time1");
    }
  if ( cmor_time_name[0] != 't' )
    {
      char *time_method = Malloc(8192 * sizeof(char));
      cdiInqAttTxt(vlistID, varID, "cell_methods", 8192, time_method);
      char *att_time_method = kv_get_a_val(kvl, "cell_methods", "");
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
      Free(time_method);
    }
  kv_insert_a_val(kvl, "time_axis", cmor_time_name, 1); 
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
  if ( !lbounds || !ubounds || pow((ubounds[1] - ubounds[0]),2) < 0.001 || pow((lbounds[1] - lbounds[0]), 2) < 0.001 )
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

static void register_z_axis(list_t *kvl, int zaxisID, char *varname, int *axis_ids, int *zfactor_id)
{
  *zfactor_id = 0;
  int zsize = zaxisInqSize(zaxisID);
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
          char *mipfreq = kv_get_a_val(kvl, "miptab_freq", "");
          char *project = kv_get_a_val(kvl, "project_id", "");
          if ( strcmp(project, "CMIP5") != 0 )
            cmor_axis(new_axis_id(axis_ids),
                        "plevs",
                        "Pa",
                        zsize,
                        (void *)levels,
                        'd', NULL, 0, NULL);
          else
            {            
              if ( strcmp(mipfreq, "cfMon") == 0 || strcmp(mipfreq, "cfDay") == 0 )
                {
                  cmor_axis(new_axis_id(axis_ids),
                        "plev7",
                        "Pa",
                        zsize,
                        (void *)levels,
                        'd', NULL, 0, NULL);
                }
              else if ( strcmp(mipfreq, "day") == 0 )
                {
                  cmor_axis(new_axis_id(axis_ids),
                        "plev8",
                        "Pa",
                        zsize,
                        (void *)levels,
                        'd', NULL, 0, NULL);
                }
              else if ( strcmp(mipfreq, "6hrPlev") == 0 )
                {
                  cmor_axis(new_axis_id(axis_ids),
                        "plev3",
                        "Pa",
                        zsize,
                        (void *)levels,
                        'd', NULL, 0, NULL);
                }
              else
                cmor_axis(new_axis_id(axis_ids),
                        "plevs",
                        "Pa",
                        zsize,
                        (void *)levels,
                        'd', NULL, 0, NULL);
            }
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
                        "cm",
                        zsize,
                        (void *)levels,
                        'd', zcell_bounds, 2, NULL);
        }
      else if ( zaxisInqType(zaxisID) == ZAXIS_GENERIC || zaxisInqType(zaxisID) == ZAXIS_HEIGHT)
        {
          char *zaxisname = Malloc(CDI_MAX_NAME * sizeof(char));
          zaxisInqName(zaxisID, zaxisname);
          if ( strcmp(zaxisname, "rho") == 0 )
            {
              char *zaxisunits = Malloc(CDI_MAX_NAME * sizeof(char));
              zaxisInqUnits(zaxisID, zaxisunits);
              if ( strcmp(zaxisunits, "kg m-3") != 0 )
                {
                  cdoAbort("For zaxis with name 'rho' the units must be kg m-3 but are: '%s'", zaxisunits);
                }
              else
                {
                  levels = Malloc(zsize * sizeof(double));
                  zaxisInqLevels(zaxisID, levels);
                  double *zcell_bounds;
                  zcell_bounds = Malloc( 2*zsize * sizeof(double) );
                  get_zcell_bounds(zaxisID, zcell_bounds, levels, zsize);
                  cmor_axis(new_axis_id(axis_ids),
                      "rho",
                      "kg m-3",
                      zsize,
                      (void *) levels,
                      'd', zcell_bounds, 2, NULL);
                }
            }
          else
            cdoAbort("Z-axis type %d with name '%s' not yet enabled.", zaxisInqType(zaxisID), zaxisname);
        }
      else
        cdoAbort("Invalid Z-axis type %d . \n", zaxisInqType(zaxisID));
      Free(zcell_bounds);
      Free(levels);
    }
  char *szc_name = kv_get_a_val(kvl, "szc", "");
  if ( zsize == 1 &&  strcmp(szc_name, "") != 0 )
    {
      char *szc_value;
      strtok_r(szc_name, "_", &szc_value);
      levels = Malloc(sizeof(double));
      levels[0] = (double) atof(szc_value);
      if ( cdoVerbose )
        printf("Attribute szc is found.\nScalar z coordinate name is: '%s'\nScalar z coordinate value is: '%f'\n", szc_name, levels[0]);
      cmor_axis(new_axis_id(axis_ids),
                      szc_name,
                      "m",
                      zsize,
                      (void *) levels,
                      'd', NULL, 0, NULL);
      Free(levels);
    }
}

/*
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
  Free(varname);
  Free(dimname);
  Free(vardimids);
}
*/

static void invert_ygriddes(list_t *kvl, int vlistID, int *gridID, int ylength, double *ycoord_vals, double *ycell_bounds, int *ynbounds)
{
  if ( ( ycoord_vals[0] - ycoord_vals[1] ) > 0 )
    {
      cdoWarning("Latitudes go north to south => reverting latitudes!\n");
      kv_insert_a_val(kvl, "invert_lat", "yes", 0);
          
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
      Free(yv2);
      Free(yb2);
    }
}

static void change_grid(char *grid_file, int *gridID, int vlistID)
{
  if ( cdoVerbose )
    printf("You configured a grid_info file. For a successfull read, the file probably needs to have at least one variable with ID 0.\n");
  argument_t *fileargument = file_argument_new(grid_file);
  int streamID2 = streamOpenRead(fileargument); 
  int vlistID2 = streamInqVlist(streamID2);
  int gridID2 = vlistInqVarGrid(vlistID2, 0); 
  vlistChangeGrid(vlistID, *gridID, gridID2);
}

static void move_lons(double *xcoord_vals, double *xcell_bounds, int xsize, int xboundsize, int xnbounds)
{  
  int testbool = 0;
  for ( int i = 0; i < xsize; i++)
    if ( xcoord_vals[i] < 0.0 )
      {
        testbool = 1;
        break;
      }
  if ( testbool > 0 )
    for ( int i = 0; i < xsize; i++ )
      xcoord_vals[i] += 180.0;
  if ( xnbounds > 1 && testbool > 0 )
    for ( int j = 0; j < xboundsize; j++ )
      xcell_bounds[j] += 180.0;
}

static void inquire_vals_and_bounds(int gridID, int *xnbounds, int *ynbounds, double *xcoord_vals, double *ycoord_vals, double *xcell_bounds, double *ycell_bounds)
{
  gridInqYvals(gridID, ycoord_vals);
  gridInqXvals(gridID, xcoord_vals);
  *xnbounds = gridInqXbounds(gridID, xcell_bounds);
  *ynbounds = gridInqYbounds(gridID, ycell_bounds);
}

static void get_cmor_table(list_t *kvl)
{
  int gridtable_id;
  char gridtable[CMOR_MAX_STRING];
  char *mip_table_dir = kv_get_a_val(kvl, "mip_table_dir", "");
  char *project_id = kv_get_a_val(kvl, "project_id", "");
  if ( mip_table_dir && project_id )
    {
      sprintf(gridtable, "%s/%s_grids\0", mip_table_dir, project_id);
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

static void check_and_gen_bounds_curv(int gridID, int totalsize, int xnbounds, int xlength, double *xcoord_vals, double *xcell_bounds, int ynbounds, int ylength, double *ycoord_vals, double *ycell_bounds)
{
  if ( xnbounds != 4 * totalsize || ynbounds != 4 * totalsize || (xcell_bounds[1] == 0.00 && xcell_bounds[2] == 0.00) || (ycell_bounds[1] == 0.00 && ycell_bounds[2] == 0.00) )
    {
      for ( int j = 1; j < ylength-1; j++ )
        for ( int i = 1; i < xlength-1; i++ )
          {
            ycell_bounds[4*(j*xlength+i)]   = ( ( ycoord_vals[j*xlength+i] + ycoord_vals[xlength*(j+1)+i] ) * 0.5 + ( ycoord_vals[xlength*j+i-1] + ycoord_vals[xlength*(j+1)+i-1] ) * 0.5 ) * 0.5;
            ycell_bounds[4*(j*xlength+i)+1] = ( ( ycoord_vals[j*xlength+i] + ycoord_vals[xlength*(j-1)+i] ) * 0.5 + ( ycoord_vals[xlength*j+i-1] + ycoord_vals[xlength*(j-1)+i-1] ) * 0.5 ) * 0.5;
            ycell_bounds[4*(j*xlength+i)+2] = ( ( ycoord_vals[j*xlength+i] + ycoord_vals[xlength*(j-1)+i] ) * 0.5 + ( ycoord_vals[xlength*j+i+1] + ycoord_vals[xlength*(j-1)+i+1] ) *0.5 ) * 0.5;
            ycell_bounds[4*(j*xlength+i)+3] = ( ( ycoord_vals[j*xlength+i] + ycoord_vals[xlength*(j+1)+i] ) * 0.5 + ( ycoord_vals[xlength*j+i+1] + ycoord_vals[xlength*(j+1)+i+1] ) * 0.5 ) * 0.5;
            xcell_bounds[4*(j*xlength+i)]   = ( ( xcoord_vals[j*xlength+i] + xcoord_vals[xlength*j+i-1] ) * 0.5 + ( xcoord_vals[xlength*(j+1)+i] + xcoord_vals[xlength*(j+1)+i-1] ) * 0.5 ) * 0.5;
            xcell_bounds[4*(j*xlength+i)+1] = ( ( xcoord_vals[j*xlength+i] + xcoord_vals[xlength*j+i-1] ) * 0.5 + ( xcoord_vals[xlength*(j-1)+i] + xcoord_vals[xlength*(j-1)+i-1] ) * 0.5 ) * 0.5;
            xcell_bounds[4*(j*xlength+i)+2] = ( ( xcoord_vals[j*xlength+i] + xcoord_vals[xlength*j+i+1] ) * 0.5 + ( xcoord_vals[xlength*(j-1)+i] + xcoord_vals[xlength*(j-1)+i+1] ) *0.5 ) * 0.5;
            xcell_bounds[4*(j*xlength+i)+3] = ( ( xcoord_vals[j*xlength+i] + xcoord_vals[xlength*j+i+1] ) * 0.5 + ( xcoord_vals[xlength*(j+1)+i] + xcoord_vals[xlength*(j+1)+i+1] ) * 0.5 ) * 0.5;
          }
/* */
/* Bottom left */
/* */
      ycell_bounds[0] = ( ycoord_vals[0] + ycoord_vals[1] ) * 0.5;
      ycell_bounds[1] =   ycoord_vals[0] + ( ycoord_vals[0] - ycoord_vals[1] ) * 0.5;
      ycell_bounds[2] =   ycoord_vals[0] + ( ycoord_vals[0] - ycoord_vals[1] ) * 0.5;
      ycell_bounds[3] = ( ycoord_vals[0] + ycoord_vals[1] ) * 0.5;
      xcell_bounds[0] =   xcoord_vals[0] + ( xcoord_vals[0] - xcoord_vals[xlength] ) * 0.5;
      xcell_bounds[1] =   xcoord_vals[0] + ( xcoord_vals[0] - xcoord_vals[xlength] ) * 0.5;
      xcell_bounds[2] = ( xcoord_vals[0] + xcoord_vals[xlength] ) * 0.5;
      xcell_bounds[3] = ( xcoord_vals[0] + xcoord_vals[xlength] ) * 0.5;

/* */
/* Bottom right */
/* */
      ycell_bounds[4*(xlength-1)]   = ( ycoord_vals[xlength-1] +   ycoord_vals[2*xlength-1] ) * 0.5;
      ycell_bounds[4*(xlength-1)+1] =   ycoord_vals[xlength-1] + ( ycoord_vals[xlength-1] - ycoord_vals[2*xlength-1] ) * 0.5;
      ycell_bounds[4*(xlength-1)+2] =   ycoord_vals[xlength-1] + ( ycoord_vals[xlength-1] - ycoord_vals[2*xlength-1] ) * 0.5;
      ycell_bounds[4*(xlength-1)+3] = ( ycoord_vals[xlength-1] +   ycoord_vals[2*xlength-1] ) * 0.5;
      xcell_bounds[4*(xlength-1)]   = ( xcoord_vals[xlength-1] +   xcoord_vals[xlength-2] ) * 0.5;
      xcell_bounds[4*(xlength-1)+1] = ( xcoord_vals[xlength-1] +   xcoord_vals[xlength-2] ) * 0.5;
      xcell_bounds[4*(xlength-1)+2] =   xcoord_vals[xlength-1] + ( xcoord_vals[xlength-1] - xcoord_vals[xlength-2] ) * 0.5;
      xcell_bounds[4*(xlength-1)+3] =   xcoord_vals[xlength-1] + ( xcoord_vals[xlength-1] - xcoord_vals[xlength-2] ) * 0.5;

/* */
/* Top left */
/* */
      ycell_bounds[4*(totalsize-xlength)]   =   ycoord_vals[totalsize-xlength] + ( ycoord_vals[totalsize-xlength] - ycoord_vals[totalsize-2*xlength] ) * 0.5;
      ycell_bounds[4*(totalsize-xlength)+1] = ( ycoord_vals[totalsize-xlength] + ycoord_vals[totalsize-2*xlength] ) * 0.5;
      ycell_bounds[4*(totalsize-xlength)+2] = ( ycoord_vals[totalsize-xlength] + ycoord_vals[totalsize-2*xlength] ) * 0.5;
      ycell_bounds[4*(totalsize-xlength)+3] =   ycoord_vals[totalsize-xlength] + ( ycoord_vals[totalsize-xlength] - ycoord_vals[totalsize-2*xlength] ) * 0.5;
      xcell_bounds[4*(totalsize-xlength)]   =   xcoord_vals[totalsize-xlength] + ( xcoord_vals[totalsize-xlength] - xcoord_vals[totalsize-xlength+1] ) * 0.5;
      xcell_bounds[4*(totalsize-xlength)+1] =   xcoord_vals[totalsize-xlength] + ( xcoord_vals[totalsize-xlength] - xcoord_vals[totalsize-xlength+1] ) * 0.5;
      xcell_bounds[4*(totalsize-xlength)+2] = ( xcoord_vals[totalsize-xlength] + xcoord_vals[totalsize-xlength+1] ) * 0.5;
      xcell_bounds[4*(totalsize-xlength)+3] = ( xcoord_vals[totalsize-xlength] + xcoord_vals[totalsize-xlength+1] ) * 0.5;

/* */
/* Top right */
/* */
      ycell_bounds[4*totalsize-4] =    ycoord_vals[totalsize-1] + ( ycoord_vals[totalsize-1] - ycoord_vals[totalsize-1-xlength] ) * 0.5;
      ycell_bounds[4*totalsize-3] = (  ycoord_vals[totalsize-1] + ycoord_vals[totalsize-1-xlength] ) * 0.5;
      ycell_bounds[4*totalsize-2] = (  ycoord_vals[totalsize-1] + ycoord_vals[totalsize-1-xlength] ) * 0.5;
      ycell_bounds[4*totalsize-1] =    ycoord_vals[totalsize-1] + ( ycoord_vals[totalsize-1] - ycoord_vals[totalsize-1-xlength] ) * 0.5;
      xcell_bounds[4*totalsize-4] = (  xcoord_vals[totalsize-1] + xcoord_vals[totalsize-2] ) * 0.5;
      xcell_bounds[4*totalsize-3] = (  xcoord_vals[totalsize-1] + xcoord_vals[totalsize-2] ) * 0.5;
      xcell_bounds[4*totalsize-2] =    xcoord_vals[totalsize-1] + ( xcoord_vals[totalsize-1] - xcoord_vals[totalsize-2] ) * 0.5;
      xcell_bounds[4*totalsize-1] =    xcoord_vals[totalsize-1] + ( xcoord_vals[totalsize-1] - xcoord_vals[totalsize-2] ) * 0.5;

      for ( int i = 1; i < xlength-1; i++)
        {
/* */
/*first values: */
/* */
          ycell_bounds[4*i]   = ( ( ycoord_vals[i] + ycoord_vals[i+xlength] ) * 0.5 + ( ycoord_vals[i-1] + ycoord_vals[i+xlength-1] ) * 0.5 ) * 0.5;
          ycell_bounds[4*i+1] =     ycoord_vals[i] + ( ycoord_vals[i] - ycoord_vals[i+xlength] ) * 0.5;
          ycell_bounds[4*i+2] =     ycoord_vals[i] + ( ycoord_vals[i] - ycoord_vals[i+xlength] ) * 0.5;
          ycell_bounds[4*i+3] = ( ( ycoord_vals[i] + ycoord_vals[i+xlength] ) * 0.5 + ( ycoord_vals[i+1] + ycoord_vals[i+xlength+1] ) * 0.5 ) * 0.5;
          xcell_bounds[4*i]   = ( ( xcoord_vals[i] + xcoord_vals[i-1] ) * 0.5 + ( xcoord_vals[i+xlength] + xcoord_vals[i+xlength-1] ) * 0.5 ) * 0.5;
          xcell_bounds[4*i+1] =   ( xcoord_vals[i] + xcoord_vals[i-1] ) * 0.5;
          xcell_bounds[4*i+2] =   ( xcoord_vals[i] + xcoord_vals[i+1] ) * 0.5;
          xcell_bounds[4*i+3] = ( ( xcoord_vals[i] + xcoord_vals[i+1] ) * 0.5 + ( xcoord_vals[i+xlength] + xcoord_vals[i+xlength+1] ) * 0.5 ) * 0.5;
/* */
/*last values: */
/* */
          ycell_bounds[4*(totalsize-xlength+i)]   =     ycoord_vals[totalsize-xlength+i] + ( ycoord_vals[totalsize-xlength+i] - ycoord_vals[totalsize-2*xlength+i] ) * 0.5;
          ycell_bounds[4*(totalsize-xlength+i)+1] = ( ( ycoord_vals[totalsize-xlength+i] + ycoord_vals[totalsize-2*xlength+i] ) * 0.5 + ( ycoord_vals[totalsize-xlength+i-1] + ycoord_vals[totalsize-2*xlength+i-1] ) * 0.5 ) * 0.5;
          ycell_bounds[4*(totalsize-xlength+i)+2] = ( ( ycoord_vals[totalsize-xlength+i] + ycoord_vals[totalsize-2*xlength+i] ) * 0.5 + ( ycoord_vals[totalsize-xlength+i+1] + ycoord_vals[totalsize-2*xlength+i+1] ) * 0.5 ) * 0.5;
          ycell_bounds[4*(totalsize-xlength+i)+3] =     ycoord_vals[totalsize-xlength+i] + ( ycoord_vals[totalsize-xlength+i] - ycoord_vals[totalsize-2*xlength+i] ) * 0.5;

          xcell_bounds[4*(totalsize-xlength+i)]   = ( xcoord_vals[totalsize-xlength+i] + xcoord_vals[totalsize-xlength+i-1] ) * 0.5;
          xcell_bounds[4*(totalsize-xlength+i)+1] = ( ( xcoord_vals[totalsize-xlength+i] + xcoord_vals[totalsize-xlength+i-1] ) * 0.5 + ( xcoord_vals[totalsize-2*xlength+i] + xcoord_vals[totalsize-2*xlength+i-1] ) * 0.5 ) * 0.5;
          xcell_bounds[4*(totalsize-xlength+i)+2] = ( ( xcoord_vals[totalsize-xlength+i] + xcoord_vals[totalsize-xlength+i+1] ) * 0.5 + ( xcoord_vals[totalsize-2*xlength+i] + xcoord_vals[totalsize-2*xlength+i+1] ) * 0.5 ) * 0.5;
          xcell_bounds[4*(totalsize-xlength+i)+3] = ( xcoord_vals[totalsize-xlength+i] + xcoord_vals[totalsize-xlength+i+1] ) * 0.5;
        }

     for ( int j = 1; j < ylength-1; j++)
        {
/* */
/*first values: */
/* */
          ycell_bounds[4*j*xlength]   = (   ycoord_vals[j*xlength] + ycoord_vals[(j+1)*xlength] ) * 0.5;
          ycell_bounds[4*j*xlength+1] = (   ycoord_vals[j*xlength] + ycoord_vals[(j-1)*xlength] ) * 0.5;
          ycell_bounds[4*j*xlength+2] = ( ( ycoord_vals[j*xlength] + ycoord_vals[(j-1)*xlength] ) * 0.5 + ( ycoord_vals[j*xlength+1] + ycoord_vals[(j-1)*xlength+1] ) * 0.5 ) * 0.5;
          ycell_bounds[4*j*xlength+3] = ( ( ycoord_vals[j*xlength] + ycoord_vals[(j+1)*xlength] ) * 0.5 + ( ycoord_vals[j*xlength+1] + ycoord_vals[(j+1)*xlength+1] ) * 0.5 ) * 0.5;
          xcell_bounds[4*j*xlength]   =     xcoord_vals[j*xlength] + ( xcoord_vals[j*xlength] - xcoord_vals[j*xlength+1] ) * 0.5;
          xcell_bounds[4*j*xlength+1] =     xcoord_vals[j*xlength] + ( xcoord_vals[j*xlength] - xcoord_vals[j*xlength+1] ) * 0.5; 
          xcell_bounds[4*j*xlength+2] = ( ( xcoord_vals[j*xlength] + xcoord_vals[j*xlength+1] ) * 0.5 + ( xcoord_vals[(j-1)*xlength] + xcoord_vals[(j-1)*xlength+1] ) * 0.5 ) * 0.5;
          xcell_bounds[4*j*xlength+3] = ( ( xcoord_vals[j*xlength] + xcoord_vals[j*xlength+1] ) * 0.5 + ( xcoord_vals[(j+1)*xlength] + xcoord_vals[(j+1)*xlength+1] ) * 0.5 ) * 0.5;
/* */
/*last values: */
/* */
          ycell_bounds[4*(j+1)*xlength-4] =  ( ( ycoord_vals[(j+1)*xlength-1] + ycoord_vals[(j+2)*xlength-1] ) * 0.5 + ( ycoord_vals[(j+1)*xlength-2] + ycoord_vals[(j+2)*xlength-2] ) * 0.5 ) * 0.5;
          ycell_bounds[4*(j+1)*xlength-3] =  ( ( ycoord_vals[(j+1)*xlength-1] + ycoord_vals[j*xlength-1] )     * 0.5 + ( ycoord_vals[(j+1)*xlength-2] + ycoord_vals[j*xlength-2] )     * 0.5 ) * 0.5;
          ycell_bounds[4*(j+1)*xlength-2] =  (   ycoord_vals[(j+1)*xlength-1] + ycoord_vals[j*xlength-1] ) * 0.5;
          ycell_bounds[4*(j+1)*xlength-1] =  (   ycoord_vals[(j+1)*xlength-1] + ycoord_vals[(j+2)*xlength-1] ) * 0.5;

          xcell_bounds[4*(j+1)*xlength-4] =  ( ( xcoord_vals[(j+1)*xlength-1] + xcoord_vals[(j+1)*xlength-2] ) * 0.5 + ( xcoord_vals[(j+2)*xlength-1] + xcoord_vals[(j+2)*xlength-2] ) * 0.5 ) * 0.5;
          xcell_bounds[4*(j+1)*xlength-3] =  ( ( xcoord_vals[(j+1)*xlength-1] + xcoord_vals[(j+1)*xlength-2] ) * 0.5 + ( xcoord_vals[j*xlength-1] + xcoord_vals[j*xlength-2] ) * 0.5 ) * 0.5;
          xcell_bounds[4*(j+1)*xlength-2] =      xcoord_vals[(j+1)*xlength-1] + ( xcoord_vals[(j+1)*xlength-1] - xcoord_vals[(j+1)*xlength-2] ) * 0.5;
          xcell_bounds[4*(j+1)*xlength-1] =      xcoord_vals[(j+1)*xlength-1] + ( xcoord_vals[(j+1)*xlength-1] - xcoord_vals[(j+1)*xlength-2] ) * 0.5;
        }
      gridDefNvertex(gridID, 4);
      gridDefXbounds(gridID, xcell_bounds);
      gridDefYbounds(gridID, ycell_bounds);
    } 
}

/*
static void select_and_register_character_dimension(char *grid_file, int *axis_ids)
{
  char *ifile = cdoStreamName(0)->args;
  if ( ifile[0] == '-' )
    cdoAbort("Cdo cmor cannot register a character dimension when several cdo operators are piped.");
  if ( strcmp(grid_file, "") == 0 )  
    register_character_dimension(axis_ids, ifile);
  else
    register_character_dimension(axis_ids, grid_file);
}
*/

static void register_grid(list_t *kvl, int vlistID, int varID, int *axis_ids, int *grid_ids)
{
  int gridID = vlistInqVarGrid(vlistID, varID);

  char *grid_file = kv_get_a_val(kvl, "ginfo", "");
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
  double *x2cell_bounds;
  double *y2cell_bounds;
  int xnbounds;
  int ynbounds;

  if ( totalsize > 1 )
  {
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
      
      cmor_axis(new_axis_id(axis_ids),    "latitude",    "degrees_north",    ylength,    (void *)ycoord_vals,    'd',    (void *)ycell_bounds,    2,    NULL);
      cmor_axis(new_axis_id(axis_ids),    "longitude",    "degrees_east",    xlength,    (void *)xcoord_vals,    'd', 
   (void *)xcell_bounds,    2,    NULL);

      Free(xcell_bounds);
      Free(ycell_bounds);
      Free(xcoord_vals);
      Free(ycoord_vals);
    }
  else if ( type == GRID_CURVILINEAR || type == GRID_UNSTRUCTURED)
    {
      xcoord_vals = Malloc(totalsize * sizeof(double));
      ycoord_vals = Malloc(totalsize * sizeof(double));
/* maximal 4 gridbounds per gridcell permitted */
      xcell_bounds = Malloc(4 * totalsize * sizeof(double));
      ycell_bounds = Malloc(4 * totalsize * sizeof(double));
      inquire_vals_and_bounds(gridID, &xnbounds, &ynbounds, xcoord_vals, ycoord_vals, xcell_bounds, ycell_bounds);
      move_lons(xcoord_vals, xcell_bounds, totalsize, 4 * totalsize, xnbounds);   
      x2cell_bounds = Malloc(4 * totalsize * sizeof(double));
      y2cell_bounds = Malloc(4 * totalsize * sizeof(double));
      get_cmor_table(kvl);
      int grid_axis[2];
      check_and_gen_bounds_curv(gridID, totalsize, xnbounds, xlength, xcoord_vals, x2cell_bounds, ynbounds, ylength, ycoord_vals, y2cell_bounds);
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
          cmor_grid(&grid_ids[0],    2,    grid_axis,    'd',    (void *)ycoord_vals,    (void *)xcoord_vals,    4,     (void *)y2cell_bounds,    (void *)x2cell_bounds);
          Free(xncoord_vals);
          Free(yncoord_vals);
        }
      /*else
        { 
          cmor_axis(&grid_axis[0],    "grid_longitude",   "degrees",    xlength,    (void *)xcoord_vals,    'd', 0, 0, NULL);
          cmor_axis(&grid_axis[1],    "grid_latitude",    "degrees",    ylength,    (void *)ycoord_vals,    'd', 0, 0, NULL);
          cmor_grid(&grid_ids[0],    2,    grid_axis,    'd',    (void *)ycoord_vals,    (void *)xcoord_vals,    2,     (void *)ycell_bounds,    (void *)xcell_bounds); 
        }*/
      Free(xcoord_vals);
      Free(ycoord_vals);
      Free(xcell_bounds);
      Free(ycell_bounds);
      Free(x2cell_bounds);
      Free(y2cell_bounds);
    }
/*  else if ( type == GRID_GENERIC )
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
          invert_ygriddes(kvl, vlistID, &gridID, ylength, ycoord_vals, ycell_bounds, &ynbounds);
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
      Free(xcoord_vals);
      Free(ycoord_vals);
    }
*/
  else if ( type == GRID_CHARXY )
    {
      grid_ids[0] = 0;
      char *xname = Malloc(CDI_MAX_NAME * sizeof(char));
      char *yname = Malloc(CDI_MAX_NAME * sizeof(char));
      gridInqXname(gridID, xname);
      gridInqYname(gridID, yname);
      char *xdimname = Malloc(CDI_MAX_NAME * sizeof(char));
      char *ydimname = Malloc(CDI_MAX_NAME * sizeof(char));
      cdiGridInqKeyStr(gridID, 902, CDI_MAX_NAME, xdimname);
      cdiGridInqKeyStr(gridID, 912, CDI_MAX_NAME, ydimname);
      int dimstrlen;   
      if ( dimstrlen = gridInqXIsc(gridID) )
        {
          char **xchars = (char **)Malloc( xlength * dimstrlen * sizeof(char));
          gridInqXCvals(gridID, xchars);
          void *xcharcmor = (void *) Malloc ( xlength * dimstrlen * sizeof(char));
          sprintf((char *)xcharcmor, "%s", xchars[0]);
          for ( int i = 1; i < xlength; i++ )
            {
              sprintf((char *)xcharcmor, "%s%s", (char *)xcharcmor, xchars[i]);
              Free(xchars[i]);
            }
          cmor_axis(new_axis_id(axis_ids), xdimname, "", xlength, (void *)xcharcmor, 'c',  NULL, dimstrlen, NULL); 
          Free(xchars);
        }
      else
        {
          xcoord_vals = Malloc(xlength * sizeof(double));
          gridInqXvals(gridID, xcoord_vals);
          xcell_bounds = Malloc(2 * xlength * sizeof(double));
          xnbounds = gridInqXbounds(gridID, xcell_bounds);
          check_and_gen_bounds(gridID, xnbounds, xlength, xcoord_vals, xcell_bounds, 1);
          cmor_axis(new_axis_id(axis_ids),    "longitude",    "degrees_east",    xlength,    (void *)xcoord_vals,    'd',    (void *)xcell_bounds,    2,    NULL);
          Free(xcell_bounds);
        }
      if ( dimstrlen = gridInqYIsc(gridID) )
        {
          char **ychars = (char **) Malloc( ylength * dimstrlen * sizeof(char));
          gridInqYCvals(gridID, ychars);
          void *ycharcmor = (void *) Malloc ( ylength * dimstrlen * sizeof(char));
          sprintf((char *)ycharcmor, "%s", ychars[0]);
          for ( int i = 1; i < ylength; i++ )
            {
              sprintf((char *)ycharcmor, "%s%s", (char *)ycharcmor, ychars[i]);
              Free(ychars[i]);
            }
          cmor_axis(new_axis_id(axis_ids), ydimname, "", ylength, (void *)ycharcmor, 'c',  NULL, dimstrlen, NULL); 
          Free(ychars);
        }
      else
        {
          ycoord_vals = Malloc(ylength * sizeof(double));
          gridInqYvals(gridID, ycoord_vals);
          invert_ygriddes(kvl, vlistID, &gridID, ylength, ycoord_vals, ycell_bounds, &ynbounds);
          ycell_bounds = Malloc(2 * ylength * sizeof(double));
          ynbounds = gridInqYbounds(gridID, ycell_bounds);
          check_and_gen_bounds(gridID, ynbounds, ylength, ycoord_vals, ycell_bounds, 0);
          cmor_axis(new_axis_id(axis_ids),    "latitude",    "degrees_north",    ylength,    (void *)ycoord_vals,    'd',    (void *)ycell_bounds,    2,    NULL);
          Free(ycell_bounds);  
        }
    }
  else
    {
      grid_ids[0] = 0;
      cdoWarning("Either the grid type is unknown or a registration is not necessary. However, it is not registered by cdo.\n");
    }
  }
  else
    grid_ids[0] = 0;
}

static void register_variable(list_t *kvl, int vlistID, int varID, int *axis_ids,
                              struct mapping *var, int *grid_ids)
{
  char *positive = get_txtatt(vlistID, varID, "positive");
  char name[CDI_MAX_NAME];
  vlistInqVarName(vlistID, varID, name);
  char units[CDI_MAX_NAME];
  vlistInqVarUnits(vlistID, varID, units);
  char *attunits = kv_get_a_val(kvl, "units", "");
  char *attname = kv_get_a_val(kvl, "out_name", "");
  char *attp = kv_get_a_val(kvl, "p", "");
  if ( attp || positive )
    check_compare_set(positive, attp, "positive");
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
            (void *) missing_value, &tolerance, positive,
                        NULL, NULL, NULL);
    }
  else
    {
      cmor_variable(&var->cmor_varID,
           name, units, count_axis_ids(axis_ids),  axis_ids,   var->datatype,
          (void *) missing_value, &tolerance, positive,
                        NULL, NULL, NULL);
    }
}

static void register_all_dimensions(list_t *kvl, int streamID,
                             struct mapping vars[], int table_id)
{
  if ( cdoVerbose )
    printf("\n*******Start to register all dimensions via cmor_axis.******\n");
  int vlistID = streamInqVlist(streamID);
  int taxisID = vlistInqTaxis(vlistID);

  char *time_units = get_time_units(taxisID);
  char *req_time_units = kv_get_a_val(kvl, "req_time_units", "");
  if ( cdoVerbose )
    printf("Checking attribute 'req_time_units' from configuration.\n");
  if ( check_time_units(req_time_units) )
    check_compare_set(time_units, req_time_units, "time_units");
  else 
    cdoAbort("Required Attribute 'req_time_units' from configuration is invalid!");

  char **requested_variables = get_requested_variables(kvl);
  if ( requested_variables == NULL && vlistNvars(vlistID) > 1 )
    cdoWarning("You have not requested any specific variable but there are several in input! Notice that all command line configuration attributes including out_name and units will be used for every variable!\n");
  int foundName = 0;
  int ps_required = 0;
  int ps_in_file = 0;
  for ( int varID = 0; varID < vlistNvars(vlistID); varID++ )
    {
      char name[CDI_MAX_NAME];
      int axis_ids[CMOR_MAX_AXES];
      axis_ids[0] = CMOR_UNDEFID;
      int zaxisID = vlistInqVarZaxis(vlistID, varID);
      vlistInqVarName(vlistID, varID, name);
      if ( requested_variables == NULL || in_list(requested_variables, name) )
        {
          if ( cdoVerbose )
            printf("\n *******Start to define variable with ID: '%d' and name: '%s'*******\n", varID, name);
          if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID )
            {
              if ( cdoVerbose )
                printf("Since the zaxis of variable '%s' is of type HYBRID, surface pressure is required from Ifile.\n It is required to have code nr. 134!\n", name);
              ps_required++;
            }
          foundName++;
          /* Time-Axis */
          char cmor_time_name[CMOR_MAX_STRING];
          get_time_method(kvl, vlistID, varID, cmor_time_name);
          if ( strcmp(cmor_time_name, "none") != 0 )
            cmor_axis(new_axis_id(axis_ids),
                    cmor_time_name,
                    time_units,
                    0,NULL, 0, NULL, 0, NULL);

          /* Grid: */
          int grid_ids[CMOR_MAX_GRIDS];
          register_grid(kvl, vlistID, varID, axis_ids, grid_ids);
          cmor_set_table(table_id);
          /* Z-Axis */
          struct mapping *var = new_var_mapping(vars);
          register_z_axis(kvl, zaxisID, name, axis_ids, &var->zfactor_id);
          /* Variable */
          register_variable(kvl, vlistID, varID, axis_ids, var, grid_ids);          
        }
    }
  if ( ps_required )
    {
      if ( cdoVerbose )
        printf("\n *******Start to find surface pressure.\n");
      for ( int varID = 0; varID < vlistNvars(vlistID); varID++ )
        if ( vlistInqVarCode(vlistID, varID) == 134 )
          {
            ps_in_file++;
            if ( requested_variables == NULL || in_list(requested_variables, "ps") )
              break;
            else
              {
                struct mapping *var = new_var_mapping(vars);
                size_t gridsize = vlistGridsizeMax(vlistID);
                var->cdi_varID = varID;
                var->help_var = 1;
                if ( vlistInqVarDatatype(vlistID, varID) == DATATYPE_FLT32 )
                  {
                    var->datatype = 'f';
                    var->data = Malloc(gridsize * sizeof(float));
                  }
                else
                  {
                    var->datatype = 'd';
                    var->data = Malloc(gridsize * sizeof(double));
                  }
                break;
              }
          }
    }
  if ( ps_required && !ps_in_file )
    cdoAbort("No surface pressure found in Ifile but required for a hybrid sigma pressure z axis!");
  if ( !foundName && requested_variables == NULL )
    cdoAbort("No variables from your table %s found in Ifile.\n");
  if ( !foundName && requested_variables )
    cdoAbort("The given variables to process by attribute vars: '%s' are not found in Ifile.\n", kv_get_a_val(kvl, "vars", ""));
  if ( cdoVerbose )
    printf("*******Called register_all_dimensions for %d variables successfully.*******\n", foundName);
  if ( requested_variables ) Free(requested_variables);
}

static char *get_frequency(list_t *kvl, int streamID, int vlistID, int taxisID)
{
  char *frequency;
  frequency = Malloc(CMOR_MAX_STRING * sizeof(char));
  int ntsteps = vlistNtsteps(vlistID);
  int reccounter = 0;
  int recdummy = 0;

  char *mipfreq = kv_get_a_val(kvl, "miptab_freq", "");
  if ( strstr(mipfreq, "yr") )
    strcpy(frequency, "yr");
  else if ( strstr(mipfreq, "mon") )
    strcpy(frequency, "mon");
  else if ( strstr(mipfreq, "day") )
    strcpy(frequency, "day");
  else if ( strstr(mipfreq, "6hr") )
    strcpy(frequency, "6hr");
  else if ( strstr(mipfreq, "3hr") )
    strcpy(frequency, "3hr");
  else
    {
      if ( cdoStreamName(0)->args[0] == '-' )
        {
          cdoWarning("Cdo cmor cannot check frequency of Ifile recs since you piped several cdo operators.\nIt is tried to use a configuration attribute frequency.");
          strcpy(frequency, kv_get_a_val(kvl, "frequency", ""));
          if ( strcmp(frequency, "") == 0 )
            cdoAbort("No attribute frequency is found.");
          else
            return frequency;
        } 
      
      int streamID2 = streamOpenRead(cdoStreamName(0));
          int vlistID2 = streamInqVlist(streamID2);
      int taxisID2 = vlistInqTaxis(vlistID2);
      if ( ntsteps < 0 )
        {
          while ( recdummy = streamInqTimestep(streamID2, reccounter++) );
          ntsteps = reccounter;
        }    
      ntsteps-=1;
      int fyear, lyear, fmonth, lmonth, dummyone, dummytwo;

      if ( ntsteps > 2 )
        {
          int recfirst = streamInqTimestep(streamID2, 0);
          cdiDecodeDate(taxisInqVdate(taxisID2), &fyear, &fmonth, &dummytwo);
          int reclast = streamInqTimestep(streamID2, ntsteps);    
          cdiDecodeDate(taxisInqVdate(taxisID2), &lyear, &lmonth, &dummytwo);

          double covered_years = lyear-fyear + 1.0;
          if ( DBL_IS_EQUAL(ntsteps / covered_years, 1) )
            strcpy(frequency, "yr");
          else if ( DBL_IS_EQUAL(ntsteps / covered_years, 12) )
            strcpy(frequency, "mon");
          else if ( DBL_IS_EQUAL(ntsteps / covered_years, 365) ||
                    DBL_IS_EQUAL(ntsteps / covered_years, 365.25) ||
                    DBL_IS_EQUAL(ntsteps / covered_years, 366) )
            strcpy(frequency, "day");
          else if ( DBL_IS_EQUAL(ntsteps / covered_years, 365*4) ||
                    DBL_IS_EQUAL(ntsteps / covered_years, 365.25*4) ||
                    DBL_IS_EQUAL(ntsteps / covered_years, 366*4) )
            strcpy(frequency, "6hr");
          else if ( DBL_IS_EQUAL(ntsteps / covered_years, 365*8) ||
                    DBL_IS_EQUAL(ntsteps / covered_years, 365.25*8) ||
                    DBL_IS_EQUAL(ntsteps / covered_years, 366*8) )
            strcpy(frequency, "3hr");
          else 
            {
              int step_per_year = 0;
              reccounter = 0;
              if ( cdoVerbose )
                printf("Frequency is calculated by counting all timesteps in year %d\nin order to calculate time bounds in case they are not given.\n", fyear, fmonth);
              while ( recdummy = streamInqTimestep(streamID2, reccounter++) )
                {
                  int reqyear;
                  cdiDecodeDate(taxisInqVdate(taxisID2), &reqyear, &lmonth, &dummytwo);
                  if ( reqyear == ( fyear + 1 ) )
                    break;
                  step_per_year++;
                } 
              int covered_months = lmonth-fmonth+1;
              if ( step_per_year > 366*8 )
                cdoAbort("Frequency is sub-3hourly! Not yet enabled.");
              else
                {
                  if ( (double)step_per_year / (double)covered_months > 31*8 )
                    cdoAbort("Frequency is sub-3hourly! Not yet enabled.");
                  else if ( (double)step_per_year / (double)covered_months > 31*4 )
                    strcpy(frequency, "3hr");
                  else if ( (double)step_per_year / (double)covered_months > 31 )
                    strcpy(frequency, "6hr");
                  else if ( (double)step_per_year / (double)covered_months > 1 )
                    strcpy(frequency, "day");
                  else
                    strcpy(frequency, "mon");
                }
              if ( cdoVerbose )
                printf("Found %d time steps in year %d.\nTherefore, the frequency is %s.\n", step_per_year, fyear, frequency);
            }
        }
      else
        {
          if ( !taxisHasBounds(taxisID2) && ntsteps > 0 )
            cdoAbort("No time bounds are found in Ifile and for %d found timesteps no frequency can be computed - at least 3 timesteps are required.\nDefine time bounds before cdo cmor.", ntsteps);
          else
            cdoWarning("For %d found timesteps no frequency can be computed - at least 3 timesteps are required.\nTime bounds of the rec are used.\n", ntsteps);
        }
    }
  return frequency;
}

static int get_tunitsec(int tunit)
{
  switch ( tunit )
    {
    case TUNIT_MINUTE: return 60; 
    case TUNIT_HOUR: return 3600; 
    case TUNIT_DAY: return 86400; 
    default: return 3600;
    }
}

static double get_cmor_time_val(int taxisID, juldate_t ref_date, int tunitsec, int calendar)
{
  juldate_t juldate = juldate_encode(calendar, taxisInqVdate(taxisID),
                                     taxisInqVtime(taxisID));
  return juldate_to_seconds(juldate_sub(juldate, ref_date)) / tunitsec;
}

static double *get_time_bounds(int taxisID, char *frequency, juldate_t ref_date, double time_val, int calendar, int tunitsec, double *time_bnds)
{
  int vdate0b, vdate1b, vtime0b, vtime1b;
  int year, month, day;
  cdiDecodeDate(taxisInqVdate(taxisID), &year, &month, &day);
  if ( !taxisHasBounds(taxisID) )
    {
      if ( strcmp(frequency, "yr") == 0 )
        {
          vdate0b = cdiEncodeDate(year,   1, 1);
          vdate1b = cdiEncodeDate(year+1, 1, 1);
        }     
      if ( strcmp(frequency, "mon") == 0 )
        {
          vdate0b = cdiEncodeDate(year, month, 1);
          month++;
          if ( month > 12 ) { month = 1; year++; }
          vdate1b = cdiEncodeDate(year, month, 1);
        }  
      if ( strcmp(frequency, "day") == 0 )
        {
          time_bnds[0] = floor(time_val);
          time_bnds[1] = ceil(time_val);
          return time_bnds;
        }  
      if ( strcmp(frequency, "6hr") == 0 || strcmp(frequency, "3hr") == 0 )
        {
          time_bnds[0] = time_val - 0.125;
          time_bnds[1] = time_val + 0.125;
          return time_bnds;
        }  
      if ( strcmp(frequency, "3hr") == 0 )
        {
          time_bnds[0] = time_val - 0.0625;
          time_bnds[1] = time_val + 0.0625;
          return time_bnds;
        }  
      vtime0b = 0;
      vtime1b = 0;
    }
  else
    {
      taxisInqVdateBounds(taxisID, &vdate0b, &vdate1b);
      taxisInqVtimeBounds(taxisID, &vtime0b, &vtime1b);
    }
  juldate_t juldate = juldate_encode(calendar, vdate0b, vtime0b);
  time_bnds[0] = juldate_to_seconds(juldate_sub(juldate, ref_date))
                / tunitsec;

  juldate = juldate_encode(calendar, vdate1b, vtime1b);
  time_bnds[1] = juldate_to_seconds(juldate_sub(juldate, ref_date))
              / tunitsec;
  return time_bnds;
}


static void read_record(int streamID, double *buffer, size_t gridsize,
                        struct mapping vars[], int invert_lat, int vlistID)
{
  int varID, levelID;
  streamInqRecord(streamID, &varID, &levelID);
  int gridID = vlistInqVarGrid(vlistID, varID);
  int type = gridInqType(gridID);
  struct mapping *var = map_var(varID, vars);
  if ( var )
    {
      int zaxisID = vlistInqVarZaxis(vlistID, varID);
      int latdim = gridInqYsize(gridID);
      int levdim = zaxisInqSize(zaxisID);
      int chardim = gridsize/latdim;
      int nmiss;
      streamReadRecord(streamID, buffer, &nmiss);
      for ( size_t i = 0; i < gridsize; i++ )
        {
/* Wrong:  (lat x basin, lev ) gridsize * levelID + i */
/* Wrong:  (basin x lat, lev) gridsize * levelID + i * chardim - ( int ) floor(i / latdim) * gridsize + ( int ) floor(i/latdim)
/* Wrong:  (basin x lev, lat ) gridsize/latdim * levdim * ( i - ( int ) floor(i/latdim) * latdim ) + ( int ) floor(i/latdim) + gridsize/latdim * levelID; */
/* Wrong:  (lat x lev, basin ) latdim * levdim * ( int ) floor(i/latdim) + ( i - ( int ) floor(i/latdim) * latdim ) + levelID * latdim*/
/* (lev x lat, basin ) */
          int newIndex;
          if ( levdim > 1 && type == GRID_CURVILINEAR )
            newIndex = i + gridsize*levelID;
          else if ( levdim > 1 )
            newIndex = i * levdim + levelID;
          else
            newIndex = i;
          if ( var->datatype == 'f' )
            {
              ((float *)var->data)[newIndex] = (float)buffer[i];
            }
          else
            {
              ((double *)var->data)[newIndex] = (double)buffer[i];
            }
        }
    }
}

static void check_for_sfc_pressure(int *ps_index, struct mapping vars[], int vlistID, int timestep)
{
  int ps_required = 0;
  for ( int j = 0; vars[j].cdi_varID != CDI_UNDEFID; j++ )
    {
      if ( vlistInqVarCode(vlistID, vars[j].cdi_varID) == 134 )
        *ps_index = j;
      else if ( zaxisInqType(vlistInqVarZaxis(vlistID, vars[j].cdi_varID)) == ZAXIS_HYBRID )
        ps_required ++;
    }
  if ( ps_index < 0 && ps_required )
    cdoAbort("No surface pressure found for time step %d but required in Hybrid-sigma-pressure-coordinates! \n", timestep);
}


static int check_append_and_size(list_t *kvl, int vlistID, char *test, int ifreq)
{
  size_t filesize = fileSize(test);
  char old_start_date[CMOR_MAX_STRING];
  char old_end_date[CMOR_MAX_STRING];
  int i = 0, j = 0;
/* Get dates from chunk string */
  while ( *(test+i) != 0 )
    {
      if ( *(test+i) == '_' )
        {
          test+=(i+1);
          i = 0;
        }
      if ( *(test+i) == '-' )
        j = i;
      i++;
    }
  if ( !i || !j || *(test+j+1) == 0 || *(test+2*j) == 0 )
    {
      cdoWarning("Error while checking chunk size for append mode.\nNew data will be appended.");
      return 0;
    }

  strncpy(old_start_date, test, j);
  old_start_date[j] = 0;
  test += (j + 1);
  strncpy(old_end_date, test, j);
  old_end_date[j] = 0;

/* Check frequency of chunk with frequency of file */

  if ( (j == 8 && ifreq !=3) || (ifreq == 3 && j != 8)
    || (j == 6 && ifreq !=2) || (ifreq == 2 && j != 6)
    || (j == 4 && ifreq !=1) || (ifreq == 1 && j != 4) )
    cdoAbort("Frequency of chunk file does not agree with frequency of the working file.");

/* Encode in julseconds depending on frequency */

  int old_start_year, old_start_month = 1, old_start_day = 1;
  int old_end_year, old_end_month = 1, old_end_day = 1;
  int new_end_year, new_end_month = 1, new_end_day = 1;

  switch ( j )
    {
    case ( 8 ):
      sscanf(old_start_date, "%04d%02d%02d", &old_start_year, &old_start_month, &old_start_day);
      sscanf(old_end_date, "%04d%02d%02d", &old_end_year, &old_end_month, &old_end_day);
      break;
    case ( 6 ):
      sscanf(old_start_date, "%04d%02d", &old_start_year, &old_start_month);
      sscanf(old_end_date, "%04d%02d", &old_end_year, &old_end_month);
      break;
    case ( 4 ):
      old_start_year = atol(old_start_date);
      old_end_year = atol(old_end_date);
      break;
    default:
      cdoAbort("Selected chunk to append data has subdaily frequency which is yet not enabled by cdo cmor.\nA new file will be written.");
    }

  char *attcalendar = kv_get_a_val(kvl, "calendar", "");
  int calendar = get_calendar_int(attcalendar);
  int cdi_startdate = cdiEncodeDate(old_start_year, old_start_month, old_start_day);
  int cdi_enddate = cdiEncodeDate(old_end_year, old_end_month, old_end_day);
  int cdi_time = cdiEncodeTime(0, 0, 0);
  juldate_t julostart = juldate_encode(calendar, cdi_startdate, cdi_time);
  juldate_t juloend = juldate_encode(calendar, cdi_enddate, cdi_time);

/* Read in first vdate in case not piped */
  if ( cdoStreamName(0)->args[0] == '-' )
    {
      cdoWarning("Cdo cmor cannot enable append mode since you piped several cdo operators.\nA new file will be written.");
      return 0;
    }
      
  int streamID2 = streamOpenRead(cdoStreamName(0));
  int vlistID2 = streamInqVlist(streamID2);
  int taxisID2 = vlistInqTaxis(vlistID2);
  juldate_t firstdate = juldate_encode(calendar, taxisInqVdate(taxisID2),
                                     taxisInqVtime(taxisID2));

/* Check temporal distance between last chunk date and first file date */
  double append_distance = juldate_to_seconds(juldate_sub(firstdate, juloend)) / 3600.0;

  if ( ( j == 8 && append_distance > 48.0 )
     ||( j == 6 && append_distance/24.0 > 45.0 )
     ||( j == 4 && append_distance/24.0/30.5 > 19.0 ) )
    cdoAbort("A temporal gap is diagnosed between end date of chunk file and first date of working file of: '%f' hours. Maximal valid gaps are:\n47 hours for daily frequency\n45 days for monthly frequency\n19 month for yearly frequency", append_distance);

/* Check file size */
  double old_interval_sec = juldate_to_seconds(juldate_sub(juloend, julostart));
  double size_per_sec = (double) filesize / old_interval_sec;

  int maxsizegb = atol(kv_get_a_val(kvl, "maxsizegb", "2"));
  int maxsizeb = maxsizegb * 1024 * 1024 * 1024;

  int ntsteps = vlistNtsteps(vlistID);
  if ( ntsteps < 0 )
    {
      cdoWarning("Could not check expected file size for exceeding maximal file size: %s gb.\nA new file will be written.", maxsizegb);
      return 0;
    }
  
  double estimated_size;
  switch ( j )
    {
    case ( 8 ):
      estimated_size = ntsteps * 60 * 60 * 24 * size_per_sec + (double) filesize ;
      break;
    case ( 6 ):
      estimated_size = ntsteps * 60 * 60 * 24 * 30.5 * size_per_sec + (double) filesize;
      break;
    case ( 4 ):
      estimated_size = ntsteps * 60 * 60 * 24 * 365.25 * size_per_sec + (double) filesize;
      break;
    default:
      cdoAbort("Selected chunk to append data has subdaily frequency which is yet not enabled by cdo cmor.\nA new file will be written.");
    }

  if ( (unsigned int)estimated_size > (unsigned int) maxsizeb )
    {
      cdoWarning("Estimated file size of appended file is : '%f'gb and exceeds maximal allowed file size: '%d'gb.\nA new file will be written.", estimated_size/1024.0/1024.0/1024.0, maxsizegb);
      return 0;
    }
  return 1;
}

static char *use_chunk_file_des(list_t *kvl, int vlistID, int var_id, char *att_chunk_des_file, int ifreq)
{
  char name[CDI_MAX_NAME];
  vlistInqVarName(vlistID, var_id, name);
  char chunk_des_file[CMOR_MAX_STRING];
  sprintf(chunk_des_file, "APPEND_FILE_%s_%s.txt", name, att_chunk_des_file);
  if ( cdoVerbose )
    printf("It is tried to open a chunk description file named: '%s' where the chunk file name should be noted\n", chunk_des_file); 
  if ( file_exist(chunk_des_file, 0) )
    {
      FILE *fp = fopen(chunk_des_file, "r");
      size_t filesize = fileSize(chunk_des_file);
      char *buffer = (char*) Malloc(filesize);
      char *appendFile = Malloc(4096 * sizeof(char));
      size_t nitems = fread(buffer, 1, filesize, fp);
      buffer = readLineFromBuffer(buffer, &filesize, appendFile, 4096);
      fclose(fp);
      if ( file_exist(appendFile, 0) && check_append_and_size(kvl, vlistID, appendFile, ifreq) )
        return appendFile;
      else
        cdoWarning("Append chunk :'%s' configured via chunk description file could not be opened.\nA new file will be written", appendFile);
    }
  else
    cdoWarning("Chunk description file '%s' could not be opened.\nA new file will be written.", chunk_des_file);
  return " ";  
}

static char **get_append_files(list_t *kvl, struct mapping vars[], int vlistID, int ifreq)
{
  char **append_files;

  char *timeaxis = kv_get_a_val(kvl, "time_axis", "");
  char *att_append_file = kv_get_a_val(kvl, "append_file", " ");
  int create_subs = atol(kv_get_a_val(kvl, "drs", "0"));
/* For chunk description file : */
  char att_chunk_des_file[CMOR_MAX_STRING];
  char *description_atts[] = {"miptab_freq", "model_id", "experiment_id", "member", NULL};
  strcpy(att_chunk_des_file, kv_get_a_val(kvl, description_atts[0], ""));
  int i = 1;
  while ( description_atts[i] != NULL )
    {
      strcat(att_chunk_des_file, "_");
      strcat(att_chunk_des_file, kv_get_a_val(kvl, description_atts[i], ""));
      i++;
    }
  for ( i = 0; vars[i].cdi_varID != CDI_UNDEFID; i++ );
  append_files = Malloc(i * sizeof(char *));

  for ( int j = 0; vars[j].cdi_varID != CDI_UNDEFID; j++ )
    {
      if ( strcmp (kv_get_a_val(kvl, "oflag", ""), "a") == 0 )
        { 
          if ( strcmp(timeaxis, "none") != 0 )
            {
              if ( ( i > 1 && strcmp(att_append_file, " ") == 0 ) || i == 1 )
                {
                  if ( create_subs && strcmp(att_append_file, " ") == 0 )
                    {  
                      if ( cdoVerbose )
                        printf("It is tried to open a chunk description file for varID: '%d'.\n", vars[j].cdi_varID);
                      char *chunk_file=use_chunk_file_des(kvl, vlistID, vars[j].cdi_varID, att_chunk_des_file, ifreq);  
                      append_files[j] = strdup(chunk_file);
                      continue;      
                    }
                  else if ( strcmp(att_append_file, " ") != 0 )
                    {
                      if ( file_exist(att_append_file, 0) && check_append_and_size(kvl, vlistID, att_append_file, ifreq) )
                        {
                          append_files[j] = strdup(att_append_file);
                          continue;
                        }
                      else
                        {
                          cdoWarning("Previous chunk configuration via attribute 'append_file': '%s' could not be opened.", att_append_file);
                        }      
                    }
                  else
                    {
                      cdoWarning("Previous chunk configuration via chunk file description not possible if DRS is not created.\nNew files will be written");
                    }
                }
              else
                {
                  cdoWarning("Previous chunk configuration via attribute 'append_file' not possible if more than one variables are requested.\nNew files will be written");
                }
            }
          else
            {
              cdoWarning("CMOR mode APPEND not possible for time independent variables.\nA new file will be written");
            }
        }
      append_files[j] = strdup(" ");
    }
  return append_files;
}

static void write_variables(list_t *kvl, int streamID, struct mapping vars[])
{
  if ( cdoVerbose )
    printf("\n*******Start to write variables via cmor_write.******\n");
  int vlistID = streamInqVlist(streamID);
  int taxisID = vlistInqTaxis(vlistID);
  int tsID = 0;
  int nrecs;
  int invert_lat = 0;
  if ( strcmp(kv_get_a_val(kvl, "invert_lat", ""), "") != 0 )
    invert_lat = 1;
  size_t gridsize = vlistGridsizeMax(vlistID);
  double *buffer = (double *) Malloc(gridsize * sizeof(double));

  int sdate, stime, time_unit, calendar;
  get_taxis(kv_get_a_val(kvl, "req_time_units", ""), kv_get_a_val(kvl, "calendar", ""), &sdate, &stime, &time_unit, &calendar);
  int tunitsec = get_tunitsec(time_unit);
  juldate_t ref_date = juldate_encode(calendar, sdate, stime);
  char *frequency = NULL;
  char *timeaxis = kv_get_a_val(kvl, "time_axis", "");
  if ( strcmp(timeaxis, "none") != 0 )
    frequency = get_frequency(kvl, streamID, vlistID, taxisID);

  int ifreq = 0;
  if ( frequency )
    {
      if ( strcmp(frequency,"yr") == 0 )
        ifreq = 1;
      if ( strcmp(frequency,"mon") == 0 )
        ifreq = 2;
      if ( strcmp(frequency,"day") == 0 )
        ifreq = 3;
    }
  char **append_files = get_append_files(kvl, vars, vlistID, ifreq);

  while ( (nrecs = streamInqTimestep(streamID, tsID++)) )
    {
      double time_bnds[2];
      double *time_bndsp;
      double time_val;
      if ( strcmp(timeaxis, "none") != 0 )
        {
          time_val = get_cmor_time_val(taxisID, ref_date, tunitsec, calendar);
          time_bndsp = ( strcmp(timeaxis, "time1") != 0 ) ? get_time_bounds(taxisID, frequency, ref_date, time_val, calendar, tunitsec, time_bnds) : 0;
        }
      while ( nrecs-- )
        read_record(streamID, buffer, gridsize, vars, invert_lat, vlistID);

      int ps_index = -1;
      check_for_sfc_pressure(&ps_index, vars, vlistID, tsID);
      for ( int i = 0; vars[i].cdi_varID != CDI_UNDEFID; i++ )
        {
          char name[CDI_MAX_NAME];
          vlistInqVarName(vlistID, vars[i].cdi_varID, name);
          if ( !vars[i].help_var )
            {
              if ( strcmp(timeaxis, "none") != 0 )
                {
                  cmor_write(vars[i].cmor_varID,
                   vars[i].data,
                   vars[i].datatype,
                   append_files[i],
                   1,
                   &time_val,
                   time_bndsp,
                   NULL); 
                  if ( vars[i].zfactor_id > 0 )
                    cmor_write(vars[i].zfactor_id,
                       vars[ps_index].data,
                       vars[ps_index].datatype,
                       append_files[i],
                       1,
                       &time_val,
                       time_bndsp,
                       &vars[i].cmor_varID); 
                }
              else
                cmor_write(vars[i].cmor_varID,
                   vars[i].data,
                   vars[i].datatype,
                   append_files[i], 0, 0, 0, NULL);
            }
        }
    }
  char file_name[CMOR_MAX_STRING];
  for ( int i = 0; vars[i].cdi_varID != CDI_UNDEFID; i++ )
    {
      if ( !vars[i].help_var )
        {
          cmor_close_variable(vars[i].cmor_varID, file_name, NULL);
          printf("*******      File stored in:  '%s' with cmor!*******\n", file_name);
        }
    }
  Free(buffer);
}

static struct mapping *construct_var_mapping(int streamID)
{
  int nvars_max = vlistNvars(streamInqVlist(streamID));
  struct mapping *vars =
    (struct mapping *) Malloc((nvars_max + 1) * sizeof(struct mapping));
  vars[0].cdi_varID = CDI_UNDEFID;
  return vars;
}

static void destruct_var_mapping(struct mapping vars[])
{
  for ( int i = 0; vars[i].cdi_varID != CDI_UNDEFID; i++ )
    Free(vars[i].data);
  Free(vars);
}

static void read_maptab(list_t *kvl, int streamID)
{
  char *maptab = kv_get_a_val(kvl, "mapping_table", "");
  char *maptabdir = kv_get_a_val(kvl, "mapping_table_dir", "");
  keyValues_t *kvvars = kvlist_search(kvl, "vars");

  if ( strcmp(maptab, "") != 0 && strcmp(maptabdir, "") != 0 && maptab[0] != '/' )
    {
      char *maptabbuild = malloc((strlen(maptab)+strlen(maptabdir)+2) * sizeof(char));
      sprintf(maptabbuild, "%s/%s\0", maptabdir, maptab);
      maptab = malloc(strlen(maptabbuild) * sizeof(char));
      strcpy(maptab, maptabbuild);
      free(maptabbuild);
    }
  if ( strcmp(maptab, "") != 0 )
    {
      int vlistID = streamInqVlist(streamID);
      int nvars = vlistNvars(vlistID);
      if ( kvvars )
        {
          const char **vars = Malloc((kvvars->nvalues + 1) * sizeof(char *));
          int i;
          for ( i = 0; i < kvvars->nvalues; i++ )
            vars[i] = strdup(kvvars->values[i]);
          vars[kvvars->nvalues] = NULL;
          apply_mapping_table(maptab, nvars, vlistID, vars);
        }
      else
        apply_mapping_table(maptab, nvars, vlistID, NULL);
      if ( maptab ) Free(maptab);
    }
}

static void parse_cmdline(list_t *pml, char **params, int nparams, char *ventry)
{
  list_t *kvl = NULL;
  kvl = list_new(sizeof(keyValues_t *), free_keyval, ventry);
  list_append(pml, &kvl);

  char *key = NULL;
  char **values = NULL;
  int i = 1, j = 0;
  int MAX_VALUES = 50; 
  while ( params[i] )
    {
      if ( strchr(params[i], '=')  )
        {
          if ( key && values[0] )
            kvlist_append(kvl, (const char *)key, (const char **) values, j);
          else if ( key )
            cdoAbort("Found no value for key '%s'.", key);
          values = malloc(MAX_VALUES * sizeof(char *));
          key = strtok(params[i], "=");
          values[0] = strtok(NULL, "");
          j = 1;
        }
      else
        {
          if ( !key )
            cdoAbort("Found no key for value '%s'.", params[i]);
          else
            {
              values[j] = params[i];
              j++;
            }
        }
      i++;
    }
  if ( key && values )
    {
      kvlist_append(kvl, (const char *)key, (const char **) values, j);
      Free(values);
    }
  else if ( values )
    cdoAbort("Found no key for value '%s'.", params[i-1]);
}

static char *get_mip_table(char *params, list_t *kvl)
{
  if ( !params )
    cdoAbort("A mip table name or path is required as first argument. No first argument found.");
  if ( file_exist(params, 0) )
    return params;
  else
    {
      cdoWarning("Your first argument could not be opened as a file. It is tried to build a path with additional configuration attributes 'mip_table_dir' and 'project_id'");
      char *miptabdir = kv_get_a_val(kvl, "mip_table_dir", "");
      char *projid = kv_get_a_val(kvl, "project_id", "");
      if ( strcmp(miptabdir, "") != 0 && strcmp(projid, "") != 0 )
        {
          char *miptab = malloc((strlen(miptabdir)+strlen(projid)+strlen(params)+3) * sizeof(char));
          sprintf(miptab, "%s/%s_%s\0", miptabdir, projid, params);
          if ( file_exist(miptab, 0) )
            return miptab;
          else
            cdoAbort("Could not open mip table '%s'.", miptab);
        }
      else
        cdoAbort("Could not build a mip table path.");
    }        
}

static void save_miptab_freq(list_t *kvl, char *mip_table)
{
  char *freq = strdup(mip_table);
  int fpos = 0, j = 0;
  while ( *(mip_table + j) )
    {
      j++;
      if ( *(mip_table + j) == '_' && *(mip_table + j + 1) )
        fpos = j + 1;
    }
  freq += fpos;
  if ( freq != NULL && fpos )
    kv_insert_a_val(kvl, "miptab_freq", freq, 0);
}
#endif

void *CMOR(void *argument)
{
  cdoInitialize(argument);

#if defined(HAVE_LIBCMOR)
  int nparams = operatorArgc();
  char **params = operatorArgv();

  /* Definition of pml: */
  list_t *pml = list_new(sizeof(list_t *), free_kvlist, "pml");
  /*Name and size of kvl: */
  const char *ventry[] = {"keyvals"};
  int nentry = (int) sizeof(ventry)/sizeof(ventry[0]);

  if ( nparams < 1 ) cdoAbort("Too few arguments!");

  /* Define kvl and read cmdline */
  parse_cmdline(pml, params, nparams, (char *)ventry[0]);

  /* Get kvl and use it from now on instead of pml */
  list_t *kvl = pmlist_get_kvlist_ventry(pml, nentry, ventry);

  /* Config files are read with descending priority. */
  read_config_files(kvl);

  /* check MIP table */
  char *mip_table = get_mip_table(params[0], kvl);
  save_miptab_freq(kvl, mip_table);

  int streamID = streamOpenRead(cdoStreamName(0));
  /* Existing attributes have lowest priority. */
  dump_special_attributes(kvl, streamID);

  /* Check for attributes and member name */
  printf("*******Start to check attributes.*******\n");
  check_attr(kvl);
  check_mem(kvl);
  printf("*******Succesfully checked global attributes.*******\n");

 /* dump_global_attributes(pml, streamID); */

 /* read mapping table */
  read_maptab(kvl, streamID);

  struct mapping *vars = construct_var_mapping(streamID);

  setup_dataset(kvl, streamID);

  int table_id;
  cmor_load_table(mip_table, &table_id);
  cmor_set_table(table_id);

  int zfactor_id = 0;
  register_all_dimensions(kvl, streamID, vars, table_id, &zfactor_id);
  write_variables(kvl, streamID, vars, &zfactor_id);

  destruct_var_mapping(vars);

  list_destroy(pml);

  streamClose(streamID);
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

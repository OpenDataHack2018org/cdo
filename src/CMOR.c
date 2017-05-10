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

list_t *maptab_search_miptab(list_t *pmlist, const char *cmorname, const char *miptab, char *key)
{
  if ( pmlist && cmorname && miptab )
    {
      listNode_t *node = pmlist->head;
      list_t *listlatest = NULL;
      while ( node )
        {
          if ( node->data )
            {
              list_t *kvlist = *(list_t **)node->data;
              keyValues_t *kvcn = kvlist_search(kvlist, key);
              if ( kvcn && kvcn->nvalues > 0 && *(kvcn->values[0]) == *cmorname && strcmp(kvcn->values[0], cmorname) == 0 )
                {
                  keyValues_t *kvmt = kvlist_search(kvlist, "mip_table");
                  if ( ( kvmt && kvmt->nvalues > 0 && *(kvmt->values[0]) == *miptab && strcmp(kvmt->values[0], miptab) == 0 ) || !kvmt )
                    return kvlist;
                  else
                    listlatest = kvlist;
                }
            }
          node = node->next;
        }
      if ( listlatest )
        {
          printf("No attribute 'mip_table' found in mapping table line for cmorname '%s'.\n The latest line of the mapping table is used.", cmorname);
          return listlatest;
        }
    }

  return NULL;
}

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
  while ( pos < len && *(pline+pos) != '=' )
    {
      if ( isspace((int) *(pline+pos)) )
        cdoAbort("Cannot interpret key of keyvalue because a blank is in: ...'%s'.... Use quotes.\n", pline);
      if ( *(pline+pos) == ',' )
        cdoAbort("Unexpected separator sign ',' in key: ...'%s'.... Use quotes.\n", pline);
      if ( *(pline+pos) == ':' )
        cdoWarning("Separator sign ':' is not supported. Use '=' instead.\n...'%s'...\n", pline);
      name[pos] = tolower(*(pline+pos));
      pos++;
    }
  name[pos] = 0;

  pline += pos;
  return pline;
}

static void copy_value(char *value, char **values, int *nvalues)
{
  if ( *nvalues > 100 )
    cdoAbort("More than 100 values for a key are not supported.");
  values[*nvalues] = strdup(value);
  (*nvalues)++;
  values[*nvalues] = NULL;
}

static void free_array(char **tofree)
{
  int i = 0;
  while ( tofree[i] )
    {
      free(tofree[i]);
      i++;
    }
  Free(tofree);
}

static void quote_replace(char **values, int nvalues, int i)
{
  char *source = values[nvalues];
  char *useful;
  *source++;
  source[i-2] = 0;
  useful = strdup(source);
  free(values[nvalues]);
  values[nvalues] = strdup(useful);
  free(useful);
}

static
char *getElementValues(char *pline, char **values, int *nvalues)
{
  while ( isspace((int) *pline) ) pline++;
  size_t len = strlen(pline);
  while ( isspace((int) *(pline+(int)len)) ) len--;
  *(pline+len) = 0;
  if ( (int)len == 0 )
    cdoAbort("No values found.\n");
  *nvalues = 0;
  int i = 0;
  while ( i < len && len )
    {
      if ( *(pline+i) == ',')
        {
          if ( i == 0 )
            cdoAbort("A value begins with ',': '%s'.\nCheck syntax.", pline);
          copy_value(pline, values, nvalues);
          if ( *(values[*nvalues-1]+i-1) == '"' || *(values[*nvalues-1]+i-1) == '\'' )
            quote_replace(values, *nvalues-1,i);
          else
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
        break;          
      else if ( *(pline+i) == '=' || *(pline+i) == ':' )
        cdoAbort("Found unexpected separator sign in value: '%c'.", *(pline+i) );
      else
        i++;
    }
  copy_value(pline, values, nvalues);
  if ( *(values[*nvalues-1]+i-1) == '"' )
    quote_replace(values, *nvalues-1, i);
  else
    *(values[*nvalues-1]+i) = 0;
  pline+=i;
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
      char **values = (char **) Malloc( 100 * sizeof(char *) );
      
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
      free_array(values);
    }
}

static void remove_space_and_comms(char **pline, char *line)
{
  while ( isspace((int) *(*pline)) ) (*pline)++;
  char *tester = *pline;
  int i = 0;
  while ( *tester != 0 )
    {
      if ( *tester == '#' || *tester == '!' )
        {
          line[i] = '\0';
          break;
        }
      i++; tester++;
    }
}

static int add_lines_tester(char *line)
{
  char *tester = line;
  while ( *tester != 0 ) tester++;
  tester--;
  while ( isspace((int) *tester) ) tester--;
  if ( *tester == ',' )
    return 1;
  return 0;
}

static void add_lines(char *line, char **buffer, size_t *buffersize)
{
  int len = strlen(line);
  char nextline[4096];
  if ( (*buffer = readLineFromBuffer(*buffer, buffersize, nextline, sizeof(nextline))) )
    {
      char *nexttester = nextline;
      remove_space_and_comms(&nexttester, nextline);
      if ( *nexttester != '\0' && *nexttester != '&' )
        {
          if ( strlen(nexttester) + len > 4096 )
            cdoAbort("Line too long!");
          strcat(line, nextline);
        }
      else
        cdoAbort("Found ',' at end of line without information in next line.");
    }
  else
    cdoAbort("Found ',' at end of file.");
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
      remove_space_and_comms(&pline, line);
      if ( *pline == '\0' ) continue;
      while ( add_lines_tester(line) )
        add_lines(line, &buffer, &buffersize);
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
          Free(kv->values[0]);
          kv->values[0] = strdup((const char *) value);
        }
    }
}

static char *kv_get_a_val(list_t *kvl, const char *key, char *replacer)
{
  keyValues_t *kv = kvlist_search(kvl, key);
  if ( kv )
    return kv->values[0];
  return replacer;
}

static char **kv_get_vals(list_t *kvl, const char *key, int *numvals)
{
  keyValues_t *kv = kvlist_search(kvl, key);
  if ( kv )
    { *numvals = kv->nvalues; return kv->values; }
  return NULL;
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

  if ( buffer ) Free(buffer);
  
  return pml;
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

static int getVarIDToMap(int vlistID, int nvars, char *key, char *value)
{
  for ( int varID = 0; varID < nvars; varID++ )
    {
      char ifilevalue[CDI_MAX_NAME];
      get_ifilevalue(ifilevalue, key, vlistID, varID);
      if ( strcmp(ifilevalue, value) == 0 )
        return varID;
    }
  return CDI_UNDEFID;
}

static void map_it(list_t *kvl, int vlistID, int varID)
{
  for ( listNode_t *kvnode = kvl->head; kvnode; kvnode = kvnode->next )
    {
      keyValues_t *kv = *(keyValues_t **)kvnode->data;
      const char *key = kv->key;
      const char *value = kv->values[0];
/*      printf("'%s' = '%s'\n", key, value); */
      if ( !value ) continue;
/* Not necessary because cmor_name is what we use for renaming :
      else if ( STR_IS_EQ(key, "name")          ) vlistDefVarName(vlistID, varID, parameter2word(value));
*/
      else if ( STR_IS_EQ(key, "cmor_name") )
        {
          char name[CDI_MAX_NAME];
          vlistInqVarName(vlistID, varID, name);
          if ( name[0] != 0 )
            cdiDefAttTxt(vlistID, varID, "original_name", (int) strlen(name), parameter2word(name));
          vlistDefVarName(vlistID, varID, parameter2word(value));
         }
      else if ( STR_IS_EQ(key, "units")        ) vlistDefVarUnits(vlistID, varID, value);
      else if ( STR_IS_EQ(key, "cell_methods")  ) cdiDefAttTxt(vlistID, varID, "cell_methods", (int) strlen(value), value);
      else if ( STR_IS_EQ(key, "character_axis")  ) cdiDefAttTxt(vlistID, varID, "character_axis", (int) strlen(value), value);
/*      else if ( STR_IS_EQ(key, "factor")        ) {}
      else if ( STR_IS_EQ(key, "delete")        ) {}
      else if ( STR_IS_EQ(key, "long_name")     ) vlistDefVarLongname(vlistID, varID, value);
      else if ( STR_IS_EQ(key, "param")         ) vlistDefVarParam(vlistID, varID, stringToParam(parameter2word(value)));
      else if ( STR_IS_EQ(key, "out_param")     ) vlistDefVarParam(vlistID, varID, stringToParam(parameter2word(value)));
      else if ( STR_IS_EQ(key, "code")          ) {}
              // else if ( STR_IS_EQ(key, "code")          ) vlistDefVarParam(vlistID, varID, cdiEncodeParam(parameter2int(value), ptab, 255));
              // else if ( STR_IS_EQ(key, "out_code")      ) vlistDefVarParam(vlistID, varID, cdiEncodeParam(parameter2int(value), ptab, 255));
*/
      else if ( STR_IS_EQ(key, "comment")       ) cdiDefAttTxt(vlistID, varID, "comment", (int) strlen(value), value);
      else if ( STR_IS_EQ(key, "p")       )
        {
          if ( !isspace(value[0]) )
            cdiDefAttTxt(vlistID, varID, "positive", (int) strlen(value), value);
        }
/*      else if ( STR_IS_EQ(key, "cell_measures") ) cdiDefAttTxt(vlistID, varID, "cell_measures", (int) strlen(value), value);
      else if ( STR_IS_EQ(key, "convert")       ) {} 
      else if ( STR_IS_EQ(key, "missval")       )   {}  
      else if ( STR_IS_EQ(key, "valid_min")     ){}
      else if ( STR_IS_EQ(key, "valid_max")     ){}
      else if ( STR_IS_EQ(key, "ok_min_mean_abs") ) {}
      else if ( STR_IS_EQ(key, "ok_max_mean_abs") ){}
      else if ( STR_IS_EQ(key, "datatype") || STR_IS_EQ(key, "type") ) {} */
      else
        {
          if ( cdoVerbose ) printf("For Mapping key '%s' is ignored.\n", key);
        }
    }
}


static int change_name_via_name(int vlistID, char *map_name, char *cmor_name)
{
  char name[CDI_MAX_NAME];
  for ( int varID = 0; varID < vlistNvars(vlistID); varID++ )
    {
      vlistInqVarName(vlistID, varID, name);
      if ( strcmp(name, map_name) == 0 )
        {
          vlistDefVarName(vlistID, varID, parameter2word((const char *) cmor_name));
          return 1;
        }
    }
  return 0;
}

static int change_name_via_code(int vlistID, char *map_code, char *cmor_name)
{
  int code;
  char codestring[4];
  for ( int varID = 0; varID < vlistNvars(vlistID); varID++ )
    {
      code = vlistInqVarCode(vlistID, varID);
      sprintf(codestring, "%03d", code);
      if ( strcmp(codestring, map_code) == 0 )
        {
          vlistDefVarName(vlistID, varID, parameter2word((const char *) cmor_name));
          return 1;
        }
    }
  return 0;
}


static int maptab_via_key(list_t *pml, int vlistID, int varID, int nventry, const char **ventry, const char *key, char *miptabfreq)
{
  char ifilevalue[CDI_MAX_NAME];
  get_ifilevalue(ifilevalue, key, vlistID, varID);

  if ( ifilevalue[0] )
    {
      list_t *kvl = maptab_search_miptab(pml, ifilevalue, miptabfreq, (char *)key);
      if ( kvl )
        {
          printf("Started mapping of variable via '%s'.\n", key);
          map_it(kvl, vlistID, varID);
          return 1;
        }
      cdoWarning("Variable named '%s' with varID '%d' could not be mapped via '%s' because no corresponding key '%s' was found in mapping table file.", ifilevalue, varID, key, key);
      return 0;
    }
  else
    {
      cdoWarning("Variable with varID '%d' could not be mapped via '%s' because it does not possess a '%s' in Ifile.", varID, key, key);
      return 0;
    }
}

static int maptab_via_cn_and_key(list_t *kvl_oname, int vlistID, int nvars, char *key)
{
  keyValues_t *kv = kvlist_search(kvl_oname, key);
  if ( kv )
    {
      int varID = getVarIDToMap(vlistID, nvars, key, kv->values[0]);
      if ( varID != CDI_UNDEFID )
        {
          printf("Started mapping of variable via '%s'.\n", key);
          map_it(kvl_oname, vlistID, varID);
          return 1;
        }
      cdoWarning("Could not map via key '%s' because no Ifile variable '%s' equals '%s'.", key, key, kv->values[0]);
    }
  else
    cdoWarning("Could not map via key '%s' because it possesses no corresponding key '%s' in mapping file.", key, key);
  return 0;
}

static void maptab_via_cmd(list_t *pml, char *origValue, int vlistID, int nvars, char *key, char *cmorName, char *miptabfreq)
{
  int varIDToMap = getVarIDToMap(vlistID, nvars, key, origValue);
  if ( varIDToMap == CDI_UNDEFID )
    cdoAbort("Could not find variable with '%s': '%s' in Ifile.", key, origValue);
  list_t *kvl_maptab = maptab_search_miptab(pml, cmorName, miptabfreq, "cmor_name");
  if ( !kvl_maptab )
    {
      cdoWarning("Could not find cmor_name '%s' in mapping table.\n No mapping table is applied.", cmorName);
      vlistDefVarName(vlistID, varIDToMap, parameter2word((const char *) cmorName));
    }
  else
    {
      printf("Started mapping of variable via '%s'.\n", key);
      map_it(kvl_maptab, vlistID, varIDToMap);
    }
}

static void maptab_via_cn(list_t *pml, char **request, int vlistID, int nvars, int numvals, char *miptabfreq)
{
  for ( int j = 0; j<numvals; j++)
    {
      list_t *kvl_oname = maptab_search_miptab(pml, request[j], miptabfreq, "cmor_name");
      if ( kvl_oname )
        {
          if ( maptab_via_cn_and_key(kvl_oname, vlistID, nvars, "name") )
            {
              printf("*******Succesfully mapped variable via name to cmor_name: '%s'.********\n", request[j]); 
              continue;
            }
          else if ( maptab_via_cn_and_key(kvl_oname, vlistID, nvars, "code") )
            {
              printf("*******Succesfully mapped variable via code to cmor_name '%s'.********\n", request[j]); 
              continue;
            }
          else
            {
              cdoWarning("No identification 'name' or 'key' in mapping table line of cmor_name '%s'. Variable not mapped.\n", request[j]);
              continue;
            }
        }
      else
        {
          cdoWarning("Requested cmor_name: '%s' is not found in mapping table.'\n", request[j]);
          continue;
        }
    }
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
  int charvars;
  char datatype;
  void *data;
};

static struct mapping *construct_var_mapping(int streamID)
{
  int nvars_max = vlistNvars(streamInqVlist(streamID));
  struct mapping *vars =
    (struct mapping *) Malloc((nvars_max + 1) * sizeof(struct mapping));
  vars[0].cdi_varID = CDI_UNDEFID;
  vars[0].data = NULL;
  vars[0].charvars = 0;
  return vars;
}

static void destruct_var_mapping(struct mapping vars[])
{
  for ( int i = 0; vars[i].cdi_varID != CDI_UNDEFID; i++ )
    Free(vars[i].data);
  Free(vars);
}

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
  vars[i + 1].data = NULL;
  vars[i + 1].charvars = 0;
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

static void addcharvar(keyValues_t *charvars, int vlistID, char *key, struct mapping vars[])
{
  if ( cdoVerbose )
    printf("*******Start to merge variables to one character coordinate.*******\n");
  int varIDs[charvars->nvalues];
  int nvars = vlistNvars(vlistID);
  for ( int i = 0; i < charvars->nvalues; i++)
    {
      varIDs[i] = getVarIDToMap(vlistID, nvars, key, charvars->values[i]);
      if ( varIDs[i] == CDI_UNDEFID )
        cdoAbort("Could not find '%s' to build a variable with character coordinate.", charvars->values[i]);
    }

  int gridID = vlistInqVarGrid(vlistID, varIDs[0]);
  int subgridID;
  int zaxisID = vlistInqVarZaxis(vlistID, varIDs[0]);
  int subzaxisID;
  int ntsteps = vlistNtsteps(vlistID);

  if ( cdoStreamName(0)->args[0] == '-' )
    cdoAbort("No variables can be merged to one character axis since you piped several cdo operators.");

  int streamID2 = streamOpenRead(cdoStreamName(0));
  if ( ntsteps == -1 )
    {
      ntsteps = 0;
      int dummy;
      while ( ( dummy = streamInqTimestep(streamID2, ntsteps++) ) );
    }

  int axissize[3];
  double *xvals, *yvals, *zvals, *subsvals;

  subsvals = Malloc(charvars->nvalues * sizeof(double));
  for ( int i = 0; i < charvars->nvalues; i++ )
    subsvals[i] = i+1;

  axissize[0] = gridInqXsize(gridID);
  axissize[1] = gridInqYsize(gridID);
  axissize[2] = zaxisInqSize(zaxisID);


  if ( axissize[0] != 1 && axissize[1] != 1 && axissize[2] != 1 )
    cdoAbort("No axis found to merge. One axis may not be allocated with more than one value.");

  int oldgridsize = axissize[0] * axissize[1];
  double *buffer_old = (double *)Malloc(oldgridsize * sizeof(double));

  for ( int i = 1; i < charvars->nvalues; i++)
    {
      gridID = vlistInqVarGrid(vlistID, varIDs[i]);  
      zaxisID = vlistInqVarZaxis(vlistID, varIDs[i]);
      if ( axissize[0] != gridInqXsize(gridID) )
        cdoAbort("Size of x-axis: '%d' of variable '%s'\n differ from x-axis size of variable '%s': '%d'.", gridInqXsize(gridID), charvars->values[i], charvars->values[0], axissize[0]);
      if ( axissize[1] != gridInqYsize(gridID) )
        cdoAbort("Size of y-axis: '%d' of variable '%s'\n differ from y-axis size of variable '%s': '%d'.", gridInqYsize(gridID), charvars->values[i], charvars->values[0], axissize[1]);
      if ( axissize[2] != zaxisInqSize(zaxisID) )
        cdoAbort("Size of z-axis: '%d' of variable '%s'\n differ from z-axis size of variable '%s': '%d'.", zaxisInqSize(zaxisID), charvars->values[i], charvars->values[0], axissize[2]);
    }

  if ( axissize[0] == 1 )
    {
      xvals = subsvals;
      yvals = Malloc(axissize[1] * sizeof(double));
      zvals = Malloc(axissize[2] * sizeof(double));
      gridInqYvals(gridID, yvals);
      zaxisInqLevels(zaxisID, zvals); 
      axissize[0] = charvars->nvalues;
    }
  else if ( axissize[1] == 1 )
    {
      xvals = Malloc(axissize[0] * sizeof(double));
      yvals = subsvals;
      zvals = Malloc(axissize[2] * sizeof(double));
      gridInqXvals(gridID, xvals);
      zaxisInqLevels(zaxisID, zvals); 
      axissize[1] = charvars->nvalues;
    }
  else if ( axissize[2] == 1 )
    {
      xvals = Malloc(axissize[0] * sizeof(double));
      yvals = Malloc(axissize[1] * sizeof(double));
      zvals = subsvals;
      gridInqXvals(gridID, xvals);
      gridInqYvals(gridID, yvals);
      axissize[2] = charvars->nvalues;
    }

  subgridID = gridCreate(GRID_GENERIC, axissize[0]*axissize[1]);
  subzaxisID = zaxisCreate(zaxisInqType(zaxisID), axissize[2]);

  gridDefXsize(subgridID, axissize[0]); 
  gridDefYsize(subgridID, axissize[1]); 
  gridDefXvals(subgridID, xvals); 
  gridDefYvals(subgridID, yvals); 
  zaxisDefLevels(subzaxisID, zvals); 

  struct mapping *var = new_var_mapping(vars);
  var->cdi_varID = vlistDefVar(vlistID, subgridID, subzaxisID,  TSTEP_INSTANT); 
  vlistDefVarName(vlistID, getVarIDToMap(vlistID, nvars+1, key, charvars->values[0]), "ChangedForMap");
  vlistDefVarName(vlistID, var->cdi_varID, charvars->values[0]);
  vlistDefVarDatatype(vlistID, var->cdi_varID,  DATATYPE_FLT64);
  vlistDefVarMissval(vlistID, var->cdi_varID, vlistInqVarMissval(vlistID, varIDs[0]));
  var->datatype = 'd';
  var->data = Malloc(ntsteps*axissize[0]*axissize[1]*axissize[2]*sizeof(double));

  int testzaxisID = vlistInqVarZaxis(vlistID, var->cdi_varID);

  int tsID = 0, nrecs = 0;

  while ( (nrecs = streamInqTimestep(streamID2, tsID)) )
    {
      while ( nrecs-- )
        {
          int varIDrw, levelIDrw, nmiss;
          streamInqRecord(streamID2, &varIDrw, &levelIDrw);
          for ( int i = 0; i < charvars->nvalues; i++ )
            if ( varIDrw == varIDs[i] )
              {
                streamReadRecord(streamID2, buffer_old, &nmiss);
                int newIndex;
                for ( int j = 0; j < oldgridsize; j++ )
                  {
/* (lev x lat, basin ) 
            newIndex = j * levdim + levelID; */
                    if ( oldgridsize == axissize[0]*axissize[1] )
                      newIndex = tsID*axissize[2]*axissize[0]*axissize[1]+i*axissize[0]*axissize[1] + j;
                    else if ( axissize[0] == charvars->nvalues )
                      newIndex = tsID*axissize[2]*axissize[0]*axissize[1]+i*axissize[1]*axissize[2] + j*axissize[2] + levelIDrw;
                    else
                      newIndex = tsID*axissize[2]*axissize[0]*axissize[1]+levelIDrw*axissize[0]*axissize[1] + i*axissize[0]*axissize[1]/oldgridsize + j;
                    ((double *)var->data)[newIndex] = (double) buffer_old[j];
                  }
              }
        }
      tsID++;
    }
  var->charvars = 1;

  streamClose(streamID2);
  Free(buffer_old);

  if ( cdoVerbose )
    printf("*******Successfully merged variables into one character axis. The final variable is called '%s' and has the ID: '%d'*******\n", charvars->values[0], var->cdi_varID);
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

  if ( buffer ) Free(buffer);
  return 0;
}

static void check_compare_set(char **finalset, char *attribute, char *attname, const char *returner)
{
  if ( !(*finalset) )
    {
      if ( !attribute )
        {
          if ( returner )
            *finalset = strdup(returner);
          else
            cdoAbort("Required value for attribute '%s' is neither found in input file nor in the configuration.", attname);
        }
      else 
        *finalset = strdup(attribute);
    }
  else if ( attribute )
    {
      if ( strcmp(attribute, *finalset) != 0 )
        {
          cdoWarning("%s of variable in input file: '%s' does not agree with configuration attribute %s: '%s'.\nCmor libary is called with attribute unit '%s'.\n", attname, *finalset, attname, attribute, attribute);
          strcpy(*finalset, attribute);
        }
    }
}

static int check_attr(list_t *kvl, char *project_id)
{
  const char *longAtt[] = {"required_time_units", "grid_info", "mapping_table", NULL};
  const char *shortAtt[] = {"rtu", "gi", "mt", NULL};

  int i = 0;
  while ( longAtt[i] != NULL )
    {
      keyValues_t *kv_latt = kvlist_search(kvl, longAtt[i]);      
      keyValues_t *kv_satt = kvlist_search(kvl, shortAtt[i]);      
      if ( kv_latt && !kv_satt )
        kv_insert_a_val(kvl, shortAtt[i], kv_latt->values[0], 1);
      else if ( !kv_latt && kv_satt )
        kv_insert_a_val(kvl, longAtt[i], kv_satt->values[0], 1);      
      i++;
    }

/* Project id moved to main void fct */
  const char *reqAtt[] = {"institute_id", "institution", "contact", "model_id", "source",
            "experiment_id", "required_time_units", NULL};
  const char *reqAttCMIP5[] = {"product", "member", NULL};
  const char *reqAttCORDEX[] = {"product", "member", "cordex_domain", "driving_model_id", NULL};
/* In all Projects needed Attributes are tested first */

  while ( reqAtt[i] != NULL )
    {
      keyValues_t *kv_reqatt = kvlist_search(kvl, reqAtt[i]);
      
      if ( !kv_reqatt || strcmp(kv_reqatt->values[0], "notSet") == 0 )
        cdoAbort("Attribute '%s' is required. Either it is missing or notSet.", reqAtt[i]);
      if ( cdoVerbose )
        printf("Attribute '%s' is '%s' \n", reqAtt[i], kv_reqatt->values[0]);
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
  if ( strcmp(project_id, "CMIP6") == 0 )
    cdoAbort("Not yet possible to create data for project CMIP6 since cmor version 2.9 is used in this operator.\n");
  if ( strcmp(project_id, "CMIP5") == 0 )
    {
      if ( cdoVerbose )
        printf("Since the project id is %s further attributes are tested. \n", project_id);
      while ( reqAttCMIP5[i] != NULL )
        {
          keyValues_t *kv_reqattCMIP5 = kvlist_search(kvl, reqAttCMIP5[i]);
          if ( !kv_reqattCMIP5 || strcmp(kv_reqattCMIP5->values[0], "notSet") == 0 )
            cdoAbort("Attribute '%s' is required. Either it is missing or notSet", reqAttCMIP5[i]);
          if ( cdoVerbose )
            printf("Attribute '%s' is '%s' \n", reqAttCMIP5[i], kv_reqattCMIP5->values[0]);
          i++;
        }
    }
  else if (strcmp(project_id, "CORDEX") == 0 )
    {
      if ( cdoVerbose )
        printf("\nSince the project id is %s further attributes are tested\n", project_id);
      i=0;
      while ( reqAttCORDEX[i] != NULL )
        {
          keyValues_t *kv_reqattCORDEX = kvlist_search(kvl, reqAttCORDEX[i]);
          if ( !kv_reqattCORDEX || strcmp(kv_reqattCORDEX->values[0], "notSet") == 0 )
            cdoAbort("Attribute '%s' is required. Either it is missing or notSet", reqAttCORDEX[i]);
          if ( cdoVerbose )
            printf("Attribute '%s' is '%s' \n", reqAttCORDEX[i], kv_reqattCORDEX->values[0]);
          i++;
        }
    }
  return 1;
} 

static int check_mem(list_t *kvl, char *project_id)
{
  char *kv_member = kv_get_a_val(kvl, "member", "");
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
  if (strcmp(project_id, "CMIP5") == 0 || strcmp(project_id, "CORDEX") == 0)
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
  char *new_history = kv_get_a_val(kvl, "history", NULL);
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
  else if ( new_history )
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
  keyValues_t *info = kvlist_search(kvl, "i");
  int i = 0;
  if ( info )
    while ( i < info->nvalues )
      {
        if ( cdoVerbose )
          printf("Try to parse file: '%s' configured with key 'info'.\n", info->values[i]);
        parse_kv_file(kvl, info->values[i]);
        i++;
      }

  /* Config file in user's $cwd directory. */
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  const char *dotconfig = ".cdocmorinfo";
  char *workfile = Malloc(strlen(cwd) + strlen(dotconfig) + 2);
  sprintf(workfile, "%s/%s", cwd, dotconfig);
  if ( cdoVerbose )
    printf("Try to parse default file: '%s'\n", workfile);
  parse_kv_file(kvl, workfile);
  Free(workfile);
  
  if ( i == 0 )
    {
      keyValues_t *info2 = kvlist_search(kvl, "i");
      if ( info2 )
        while ( i < info2->nvalues )
          {
            if ( cdoVerbose )
              printf("Try to parse file: '%s' configured with key 'info' in file '.cdocmorinfo'.\n", info2->values[i]);
            parse_kv_file(kvl, info2->values[i]);
            i++;
          }
    }
}

static int in_list(char **list, const char *needle, int num)
{
  for ( int i = 0; i < num; i++ )
    if ( strcmp(list[i], needle) == 0 )
      return 1;
  return 0;
}

static int get_netcdf_file_action(list_t *kvl)
{
  char *chunk = kv_get_a_val(kvl, "om", "r");
  if ( chunk[0] == 'r' )
    return CMOR_REPLACE;
  else if ( chunk[0] == 'a')
    return CMOR_APPEND;
  else if ( chunk[0] == 'p')
    return CMOR_PRESERVE;
  else
    {
      cdoWarning("No valid CMOR output mode! \nAttribute output_mode is '%s', but valid are 'a' for append ,'r' for replace or 'p' for preserve.\nCMOR output mode is set to: replace.", chunk);
      return CMOR_REPLACE;
    }
}

static int get_cmor_verbosity(list_t *kvl)
{
  char *verbos = kv_get_a_val(kvl, "set_verbosity", NULL);
  if ( !verbos )
    return CMOR_NORMAL;
  if ( strcmp(verbos, "CMOR_QUIET") == 0 )
    return CMOR_QUIET;
  else
    return CMOR_NORMAL;
}

static int get_cmor_exit_control(list_t *kvl)
{
  char *exit = kv_get_a_val(kvl, "exit_control", NULL);
  if ( !exit )
    return CMOR_NORMAL;
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
      strcpy(calendar_ptr, "gregorian"); break;
    case CALENDAR_PROLEPTIC:
      strcpy(calendar_ptr, "proleptic_gregorian"); break;
    case CALENDAR_360DAYS:
      strcpy(calendar_ptr, "360_day"); break;
    case CALENDAR_365DAYS:
      strcpy(calendar_ptr, "noleap"); break;
    case CALENDAR_366DAYS:
      strcpy(calendar_ptr, "all_leap"); break;
    default:
      Free(calendar_ptr); return NULL;
    }
  return calendar_ptr;
}

static int get_calendar_int(char *calendar)
{
  if ( !calendar )
    return 0;
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
      cdoWarning("Calendar type %s is not supported by CMOR.\n", calendar);
      return 0;
    }
}

static char *get_txtatt(int vlistID, int varID, char *key)
{
  int natts;
  cdiInqNatts(vlistID, varID, &natts);
  char *txtatt = NULL;
  for ( int i = 0; i < natts; i++ )
    {
      char name[CDI_MAX_NAME];
      char buffer[8];
      int type, len;
      cdiInqAtt(vlistID, varID, i, name, &type, &len);
      if ( strcmp(name, key) == 0 )
        {
          txtatt = Malloc(CMOR_MAX_STRING * sizeof(char));
          cdiInqAttTxt(vlistID, varID, name, len, txtatt);
          txtatt[len] = '\0';
          return txtatt;
        }
    }
  return txtatt;
}

static void setup_dataset(list_t *kvl, int streamID, int *calendar)
{
  if ( cdoVerbose )
    printf("*******Start to process cmor_setup and cmor_dataset.*******\n");
  int netcdf_file_action = get_netcdf_file_action(kvl);
  int set_verbosity = get_cmor_verbosity(kvl);
  int exit_control = get_cmor_exit_control(kvl);
  int creat_subs = 1;
  char *drs = kv_get_a_val(kvl, "d", "y");
  if ( drs[0] == 'n' )
    creat_subs = 0;
  else if ( drs[0] != 'y' )
    {
      cdoWarning("Unknown value for keyword 'drs' is found: '%s'.\nAllowed are: 'n' or 'y'. DRS is set to 'y'.", drs);
      kv_insert_a_val(kvl, "d", "y", 1);
    }

  int vlistID = streamInqVlist(streamID);

  cmor_setup(kv_get_a_val(kvl, "inpath", "/usr/share/cmor/"),
             &netcdf_file_action,
             &set_verbosity,
             &exit_control,
             kv_get_a_val(kvl, "logfile", NULL),
             &creat_subs);

  int taxisID = vlistInqTaxis(streamInqVlist(streamID));

/*
  char *attcomment = kv_get_a_val(kvl, "comment", NULL);
  char *comment = get_txtatt(vlistID, CDI_GLOBAL, "comment");
*/
  
  char *attcalendar = kv_get_a_val(kvl, "calendar", NULL);
  char *calendarptr = get_calendar_ptr(taxisInqCalendar(taxisID));
  if ( cdoVerbose )
    printf("Checking attribute 'calendar' from configuration.\n");
  if ( *calendar = get_calendar_int(attcalendar) )
    check_compare_set(&calendarptr, attcalendar, "calendar", NULL);
  else 
    {
      if ( cdoVerbose )
        printf("Try to use Ifile calendar.\n");
      if ( !get_calendar_int(calendarptr) )
        cdoAbort("No valid configuration and no valid Ifile calendar found.");
      else
        *calendar = get_calendar_int(calendarptr);
    }

#if defined(CMOR_VERSION_MAJOR)
  int cmor_version_exists = 1;
  if ( CMOR_VERSION_MAJOR == 2 && CMOR_VERSION_MINOR == 9 )
    {
      double branch_time = atof(kv_get_a_val(kvl, "branch_time", "0.0"));
      cmor_dataset(kv_get_a_val(kvl, "dr", "./"),
               kv_get_a_val(kvl, "experiment_id", ""),
               kv_get_a_val(kvl, "institution", ""),
               kv_get_a_val(kvl, "source", ""),
               calendarptr,
               atoi(kv_get_a_val(kvl, "realization", "")),
               kv_get_a_val(kvl, "contact", ""),
               kv_get_a_val(kvl, "history", ""),
/* comment:*/
               kv_get_a_val(kvl, "comment", ""),
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
  Free(calendarptr);
  if ( cdoVerbose )
    printf("*******Setup finished successfully.*******\n");
}


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
/* Required attribute in check_att */
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

static void get_taxis(char *required_time_units, int *sdate, int *stime, int *timeunit)
{
  int attyear, attmonth, attday, atthour, attminute, attsecond;
  char atttimeunit[CMOR_MAX_STRING];

  sscanf(required_time_units, "%s since %d-%d-%d%*1s%02d:%02d:%02d%*1s",
                  atttimeunit, &attyear, &attmonth, &attday, &atthour,
                  &attminute, &attsecond);
  *sdate = cdiEncodeDate(attyear, attmonth, attday);
  *stime = cdiEncodeTime(atthour, attminute, attsecond);
  *timeunit = get_time_step_int(atttimeunit);
}

static void get_time_method(list_t *kvl, int vlistID, int varID, char *cmor_time_name, char *project_id, int miptab_freq, int *time_axis)
{
  if (strcmp(project_id, "CMIP5") == 0 && miptab_freq )
    switch ( miptab_freq )
      {
      case 1: strcpy(cmor_time_name, "time2"); break;
      case 2: strcpy(cmor_time_name, "time"); break;
      case 4: strcpy(cmor_time_name, "time"); break;
      case 5: strcpy(cmor_time_name, "time1"); break;
      case 6: strcpy(cmor_time_name, "time1"); break;
      }
  if ( cmor_time_name[0] != 't' )
    {
      char *time_method = get_txtatt(vlistID, varID, "cell_methods");
      char *att_time_method = kv_get_a_val(kvl, "cm", NULL);
      check_compare_set(&time_method, att_time_method, "cell_methods", " ");
      if ( time_method[0] == 'm' )      { strcpy(cmor_time_name, "time \0"); *time_axis=0; }
      else if ( time_method[0] == 'p' ) { strcpy(cmor_time_name, "time1\0"); *time_axis=1; }
      else if ( time_method[0] == 'c' ) { strcpy(cmor_time_name, "time2\0"); *time_axis=2; }
      else if ( time_method[0] == 'n' ) { strcpy(cmor_time_name, "none\0"); *time_axis=3; }
      else
        {
          cdoWarning("Found configuration time cell method '%s' is not valid. Check CF-conventions for allowed time cell methods.\nTime cell method is set to 'mean'. \n", time_method);
          strcpy(cmor_time_name, "time \0");
        }
      Free(time_method);
    }
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

static void get_zhybrid(int zaxisID, double *p0, double *alev_val, double *alev_bnds, double *b_val, double *b_bnds, double *ap_val, double *ap_bnds)
{
  int zsize = zaxisInqSize(zaxisID);
  int vctsize = zaxisInqVctSize(zaxisID);
  double *vct = Malloc(vctsize * sizeof(double) );
  zaxisInqVct(zaxisID, vct);
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
  Free(vct);
}

static int get_strmaxlen(char **array, int len)
{
  int result = 0, i;
  for (i = 0; i < len; i++)
    if ( result < strlen(array[i]) )
      result = strlen(array[i]);
  return result;     
}

static void register_z_axis(list_t *kvl, int vlistID, int varID, int zaxisID, char *varname, int *axis_ids, int *zfactor_id, char *project_id, int miptab_freq)
{
  *zfactor_id = 0;
  int zsize = zaxisInqSize(zaxisID);
  double *levels;

  char *chardimatt = kv_get_a_val(kvl, "ca", NULL);
  char *chardim = get_txtatt(vlistID, varID, "character_axis");
  check_compare_set(&chardim, chardimatt, "character_axis", "notSet");
  if ( strcmp(chardim, "vegtype") == 0 || strcmp(chardim, "oline") == 0  )
    {
      if ( zsize )
        cdoWarning("You configured a character coordinate '%s' but a zaxis is found with '%d' numerical values. The zaxis attributes are ignored and the '%d' levels are interpreted as the character coordinates in the order they are given for '%s'.", chardim, zsize, zsize, varname);
      int numchar = 0;
      char *charvalstring = Malloc(CMOR_MAX_STRING * sizeof(char));
      sprintf(charvalstring, "char_axis_%s", chardim);
      char **charvals = kv_get_vals(kvl, charvalstring, &numchar);
      Free(charvalstring);
      if ( charvals )
        {
          int maxlen = get_strmaxlen(charvals, numchar);
          void *charcmor = (void *) Malloc ( numchar * maxlen * sizeof(char));
          sprintf((char *)charcmor, "%s", charvals[0]);
          char blanks[maxlen];
          for ( int i = 0; i < maxlen; i++)
            blanks[i] = ' ';
          sprintf((char *)charcmor, "%s%.*s", (char *)charcmor, maxlen-strlen(charvals[0]), blanks);         
          for ( int i = 1; i < numchar; i++ )
            {
              sprintf((char *)charcmor, "%s%s", (char *)charcmor, charvals[i]);
              sprintf((char *)charcmor, "%s%.*s", (char *)charcmor, maxlen-strlen(charvals[i]), blanks);         
            }
          if ( numchar == zsize )
            cmor_axis(new_axis_id(axis_ids), chardim, "", numchar, (void *)charcmor, 'c',  NULL, maxlen, NULL); 
          else
            cdoAbort("The number of registered character coordinates '%d' differ from the number of axis levels '%d'.", numchar, zsize);
          Free(charcmor);
        }
      else
        cdoAbort("You configured a character coordinate '%s' but no values are found! Configure values via attribute 'char_dim_vals'!", chardim);
      Free(chardim);
    }
  else
  {
  if ( zsize > 1)
    {
      levels = Malloc(zsize * sizeof(double));
      zaxisInqLevels(zaxisID, levels);
      double *zcell_bounds;
      zcell_bounds = Malloc( 2*zsize * sizeof(double) );
      get_zcell_bounds(zaxisID, zcell_bounds, levels, zsize);
      if ( zaxisInqType(zaxisID) == ZAXIS_PRESSURE )
        {
          if ( strcmp(project_id, "CMIP5") != 0 )
            cmor_axis(new_axis_id(axis_ids),
                        "plevs",
                        "Pa",
                        zsize,
                        (void *)levels,
                        'd', NULL, 0, NULL);
          else
            {  
              switch ( miptab_freq )
                {
                case 3: cmor_axis(new_axis_id(axis_ids),
                        "plev7",
                        "Pa",
                        zsize,
                        (void *)levels,
                        'd', NULL, 0, NULL); break;
                case 4: cmor_axis(new_axis_id(axis_ids),
                        "plev8",
                        "Pa",
                        zsize,
                        (void *)levels,
                        'd', NULL, 0, NULL); break;
                case 5: cmor_axis(new_axis_id(axis_ids),
                        "plev3",
                        "Pa",
                        zsize,
                        (void *)levels,
                        'd', NULL, 0, NULL); break;
                default: cmor_axis(new_axis_id(axis_ids),
                        "plevs",
                        "Pa",
                        zsize,
                        (void *)levels,
                        'd', NULL, 0, NULL); break;
                }                
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

          char *mtproof = kv_get_a_val(kvl, "mtproof", NULL);
          if ( mtproof )
            {
              if ( cdoVerbose )
                printf("Try to apply mapping table: '%s' for ps.\n", mtproof);
              list_t *pml = cdo_parse_cmor_file(mtproof);
              if ( pml == NULL )
                cdoWarning("Mapping table: '%s' could not be parsed. Infile variable name needs to be 'ps'.", mtproof);
              else
                {
                  char *tempo[] = {"ps"};
                  maptab_via_cn(pml, tempo, vlistID, vlistNvars(vlistID), 1, kv_get_a_val(kvl, "miptab_freq", NULL)); 
                  if ( cdoVerbose )
                    printf("*******Succesfully applied mapping table: '%s' for ps.*******\n", mtproof);
                  list_destroy(pml);
                }
            }
          else
            {
              if ( cdoVerbose )
                printf("Ps needs to be one infile variable name.");
            }
          int psID = getVarIDToMap(vlistID, vlistNvars(vlistID), "name", "ps");
          if ( psID == CDI_UNDEFID )
            cdoAbort("Could not find a surface pressure variable in infile. Cannot register a hybrid zaxis without surface pressure.");

          get_zhybrid(zaxisID, p0, alev_val, alev_bnds, b_val, b_bnds, ap_val, ap_bnds);
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
              Free(zaxisunits);
            }
          else
            cdoAbort("Z-axis type %d with name '%s' not yet enabled.", zaxisInqType(zaxisID), zaxisname);
          Free(zaxisname);
        }
      else
        cdoAbort("Invalid Z-axis type %d . \n", zaxisInqType(zaxisID));
      Free(zcell_bounds);
      Free(levels);
    }
  char *szc_name = kv_get_a_val(kvl, "szc", NULL);
  if ( zsize == 1 &&  szc_name )
    {
      char *szc_value = NULL;
      strtok_r(szc_name, "_", &szc_value);
      if ( !szc_value || !szc_value[0] )
        cdoAbort("Could not find an underscore '_' in szc value '%s' to seperate axis name from axis value", szc_name);
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
}

/*
static void register_character_dimension(int *axis_ids, char *filename)
{
  printf("The grid type is generic and a dimension 'basin' is found.\nTherefore, it is tried to read the character dimension.\n");
  int nc_file_id, nfiledims, nvars, ngatts, unlimdimid;
  nc_type xtypep;
  int varndims, varnattsp;
  int *vardimids;

  char *varname = Malloc(36 * sizeof(char));
  char *dimname = Malloc(36 * sizeof(char));

  size_t dimlength, dimstrlength;

  nc_open(filename, NC_NOWRITE, &nc_file_id);
  nc_inq(nc_file_id, &nfiledims, &nvars, &ngatts, &unlimdimid);
  vardimids = Malloc(nfiledims * sizeof(int));
  void *final_chardim;
  for ( int i = 0; i < nvars; i++ )
    {
      nc_inq_var(nc_file_id, i, varname, &xtypep, &varndims, vardimids, &varnattsp);
      if ( strcmp(varname, "region") == 0 )
        {
          nc_inq_dim(nc_file_id, vardimids[1], dimname, &dimstrlength);
          nc_inq_dim(nc_file_id, vardimids[0], dimname, &dimlength);

          final_chardim = (void *)Malloc(dimstrlength * dimlength *sizeof(char));
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

static void change_grid(char *grid_file, int gridID, int vlistID)
{
  if ( cdoVerbose )
    printf("You configured a grid_info file: '%s'. It is tested for a valid use as substitution.\n");
  argument_t *fileargument = file_argument_new(grid_file);
  int streamID2 = streamOpenRead(fileargument); 
  int vlistID2 = streamInqVlist(streamID2);
  int gridID2 = vlistInqVarGrid(vlistID2, 0); 

  if ( !gridID2 )
    cdoAbort("Could not use grid from file '%s' configured via attribute 'ginfo'\n because of internal problems.", grid_file);

  int a,b;
  a = gridInqSize(gridID);
  b = gridInqSize(gridID2);
  if ( a != b )
    cdoAbort("Could not use grid from file '%s' configured via attribute 'ginfo'\n because total size of $IFILE: '%d' is not identical to total size of ginfo file: '%d'.", grid_file, a, b);

  a = gridInqYsize(gridID);
  b = gridInqYsize(gridID2);
  if ( a != b )
    cdoAbort("Could not use grid from file '%s' configured via attribute 'ginfo'\n because ysize of $IFILE: '%d' is not identical to ysize of ginfo file: '%d'.", grid_file, a, b);

  a = gridInqXsize(gridID);
  b = gridInqXsize(gridID2);
  if ( a != b )
    cdoAbort("Could not use grid from file '%s' configured via attribute 'ginfo'\n because xsize of $IFILE: '%d' is not identical to xsize of ginfo file: '%d'.", grid_file, a, b);

  vlistChangeGrid(vlistID, gridID, gridID2);
  printf("Succesfully substituted grid.\n");

  streamClose(streamID2);
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
      if ( xcoord_vals[i] < 0 )
        xcoord_vals[i] += 360.0;
  if ( xnbounds > 1 && testbool > 0 )
    for ( int j = 0; j < xboundsize; j++ )
      if ( xcell_bounds[j] < 0 )
        xcell_bounds[j] += 360.0;
}

static void inquire_vals_and_bounds(int gridID, int *xnbounds, int *ynbounds, double *xcoord_vals, double *ycoord_vals, double *xcell_bounds, double *ycell_bounds)
{
  gridInqYvals(gridID, ycoord_vals);
  gridInqXvals(gridID, xcoord_vals);
  *xnbounds = gridInqXbounds(gridID, xcell_bounds);
  *ynbounds = gridInqYbounds(gridID, ycell_bounds);
}

static void get_cmor_table(list_t *kvl, char *project_id)
{
  int gridtable_id;
  char gridtable[CMOR_MAX_STRING];
  char *mip_table_dir = kv_get_a_val(kvl, "mip_table_dir", NULL);
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

static double lonbnds_mids_trans_check(double value1, double value2)
{
  if ( abs(value1 - value2) < 180.0 )
    return (value1 + value2) * 0.5;
  else 
    {
      if ( value1 + value2 < 360.0 )
        return (value1 + value2 + 360.0) * 0.5;
      else
        return (value1 + value2 + 360.0) * 0.5 - 360.0;
    }
}

static double lonbnds_bnds_trans_check(double value1, double value2)
{
  if ( abs(value1 - value2) < 180 )
    {
      if ( 2*value1 < value2 )
        return (2*value1 - value2 + 360.0);
      else if ( 2*value1 > value2 + 360.0 )
        return (2*value1 - value2 - 360.0);
      else
        return (2*value1 - value2);
    }
  else if ( value1 - value2 > 180  )
    return (2*value1 - value2 - 360.0);
  else
    return (2*value1 - value2 + 360.0);
}

static void check_and_gen_bounds_curv(int gridID, int totalsize, int xnbounds, int xlength, double *xcoord_vals, double *xcell_bounds, int ynbounds, int ylength, double *ycoord_vals, double *ycell_bounds)
{ 
  if ( xnbounds != 4 * totalsize || ynbounds != 4 * totalsize || (xcell_bounds[1] == 0.00 && xcell_bounds[2] == 0.00) || (ycell_bounds[1] == 0.00 && ycell_bounds[2] == 0.00) )
    {
      double halflons[xlength+1][ylength];
      double halflats[xlength][ylength+1];
      double halflonsOnhalflats[xlength+1][ylength+1];
      double halflatsOnhalflons[xlength+1][ylength+1];

/**/
/*************Half-lons with 360-0 transmission check**************/
/**/
      for ( int j = 0; j < ylength; j++ )
        {
          for ( int i = 1; i < xlength; i++ )
            halflons[i][j] = lonbnds_mids_trans_check(xcoord_vals[i-1+j*xlength], xcoord_vals[i+j*xlength]);
/*left and right boundary: */
          halflons[0][j]       = lonbnds_bnds_trans_check(xcoord_vals[j*xlength], halflons[1][j]);
          halflons[xlength][j] = lonbnds_bnds_trans_check(xcoord_vals[j*xlength-1], halflons[xlength-1][j]);
        }
/**/
/*************Half-lats **************/
/**/
      for ( int i = 0; i < xlength; i++ )
        {
          for ( int j = 1; j < ylength; j++ )
            halflats[i][j] = (ycoord_vals[i+(j-1)*xlength] + ycoord_vals[i+j*xlength]) * 0.5;
/*upper and lower boundary: */
          halflats[i][0]       = 2*ycoord_vals[i] - halflats[i][1];
          halflats[i][ylength] = 2*ycoord_vals[i+(ylength-1)*xlength] - halflats[i][ylength-1];
        }
/**/
/****************Half-lons-on-half-lats with 0-360 transmission check**********/
/****************Half-lats-on-half-lons                              **********/
/**/

      for ( int i = 1; i < xlength; i++ )
        {
          for ( int j = 1; j < ylength; j++ )
            {
              halflonsOnhalflats[i][j] = lonbnds_mids_trans_check(halflons[i][j-1], halflons[i][j]);
              halflatsOnhalflons[i][j] = ( halflats[i-1][j] + halflats[i][j] ) * 0.5;
            }
/*upper and lower boundary: */
          halflonsOnhalflats[i][0]       = lonbnds_bnds_trans_check(halflons[i][0], halflonsOnhalflats[i][1]);
          halflonsOnhalflats[i][ylength] = lonbnds_bnds_trans_check(halflons[i][ylength-1], halflonsOnhalflats[i][ylength-1]);
          halflatsOnhalflons[i][0]       = ( halflats[i-1][0] + halflats[i][0] ) * 0.5;
          halflatsOnhalflons[i][ylength] = ( halflats[i-1][ylength] + halflats[i][ylength] ) * 0.5;
        }      

/*left and right boundary: */
      for ( int j = 1; j < ylength; j++ )
        {
          halflonsOnhalflats[0][j]       = lonbnds_mids_trans_check(halflons[0][j-1], halflons[0][j]);
          halflonsOnhalflats[xlength][j] = lonbnds_mids_trans_check(halflons[xlength][j-1], halflons[xlength][j]);

          halflatsOnhalflons[0][j]       = 2*halflats[0][j] - halflatsOnhalflons[1][j];
          halflatsOnhalflons[xlength][j] = 2*halflats[xlength-1][j] - halflatsOnhalflons[xlength-1][j];
        }
      halflatsOnhalflons[0][0]             = 2*halflats[0][0] - halflatsOnhalflons[1][0];
      halflatsOnhalflons[0][ylength]       = 2*halflats[0][ylength] - halflatsOnhalflons[1][ylength];
      halflatsOnhalflons[xlength][0]       = 2*halflats[xlength-1][0] - halflatsOnhalflons[xlength-1][0];
      halflatsOnhalflons[xlength][ylength] = 2*halflats[xlength-1][ylength] - halflatsOnhalflons[xlength-1][ylength];

      halflonsOnhalflats[0][0]             = lonbnds_bnds_trans_check(halflons[0][0], halflonsOnhalflats[0][1]);
      halflonsOnhalflats[0][ylength]       = lonbnds_bnds_trans_check(halflons[0][ylength-1], halflonsOnhalflats[0][ylength-1]);
      halflonsOnhalflats[xlength][0]       = lonbnds_bnds_trans_check(halflons[xlength][0], halflonsOnhalflats[xlength][1]);
      halflonsOnhalflats[xlength][ylength] = lonbnds_bnds_trans_check(halflons[xlength][ylength-1], halflonsOnhalflats[xlength-1][ylength]);

      for ( int i = 0; i < xlength; i++ )
        for ( int j = 0; j < ylength; j++ )
          {
            xcell_bounds[4*(j*xlength+i)]   = halflonsOnhalflats[i][j+1];
            xcell_bounds[4*(j*xlength+i)+1] = halflonsOnhalflats[i][j];
            xcell_bounds[4*(j*xlength+i)+2] = halflonsOnhalflats[i+1][j];
            xcell_bounds[4*(j*xlength+i)+3] = halflonsOnhalflats[i+1][j+1];
            ycell_bounds[4*(j*xlength+i)]   = halflatsOnhalflons[i][j+1];
            ycell_bounds[4*(j*xlength+i)+1] = halflatsOnhalflons[i][j];
            ycell_bounds[4*(j*xlength+i)+2] = halflatsOnhalflons[i+1][j];
            ycell_bounds[4*(j*xlength+i)+3] = halflatsOnhalflons[i+1][j+1];
          }
      gridDefNvertex(gridID, 4);
      gridDefXbounds(gridID, xcell_bounds);
      gridDefYbounds(gridID, ycell_bounds);
    }
}
/*

  if ( xnbounds != 4 * totalsize || ynbounds != 4 * totalsize || (xcell_bounds[1] == 0.00 && xcell_bounds[2] == 0.00) || (ycell_bounds[1] == 0.00 && ycell_bounds[2] == 0.00) )
    {
      for ( int j = 1; j < ylength-1; j++ )
        for ( int i = 1; i < xlength-1; i++ )
          {
            double *star[9] =
 { &xcoord_vals[(j-1)*xlength+(i-1)], &xcoord_vals[j*xlength+(i-1)], &xcoord_vals[(j+1)*xlength+(i-1)],
   &xcoord_vals[(j-1)*xlength+i],     &xcoord_vals[j*xlength+i],     &xcoord_vals[(j+1)*xlength+i],
   &xcoord_vals[(j-1)*xlength+(i+1)], &xcoord_vals[j*xlength+(i+1)], &xcoord_vals[(j+1)*xlength+i+1] };
            double max = 0, min = 0;
            for ( int k = 0; k < 9; k++ )
              {
                max = (max < *star[k]) ? *star[k] : max;
                min = (min > *star[k]) ? *star[k] : min;
              }
            if ( ( max - min ) > 270 )
              {
                if ( *star[4] < 90 )
                  for ( int l = 0; l < 9; l++)
                    *star[l] = (*star[l] > 270 ) ? *star[l] - 360.0 : *star[l];
                else if ( *star[4] > 270 )
                  for ( int l = 0; l < 9; l++)
                    *star[l] = (*star[l] < 90 ) ? *star[l] + 360.0 : *star[l];
              }

            ycell_bounds[4*(j*xlength+i)]   = ( ( ycoord_vals[j*xlength+i] + ycoord_vals[xlength*(j+1)+i] ) * 0.5 + ( ycoord_vals[xlength*j+i-1] + ycoord_vals[xlength*(j+1)+i-1] ) * 0.5 ) * 0.5;
            ycell_bounds[4*(j*xlength+i)+1] = ( ( ycoord_vals[j*xlength+i] + ycoord_vals[xlength*(j-1)+i] ) * 0.5 + ( ycoord_vals[xlength*j+i-1] + ycoord_vals[xlength*(j-1)+i-1] ) * 0.5 ) * 0.5;
            ycell_bounds[4*(j*xlength+i)+2] = ( ( ycoord_vals[j*xlength+i] + ycoord_vals[xlength*(j-1)+i] ) * 0.5 + ( ycoord_vals[xlength*j+i+1] + ycoord_vals[xlength*(j-1)+i+1] ) *0.5 ) * 0.5;
            ycell_bounds[4*(j*xlength+i)+3] = ( ( ycoord_vals[j*xlength+i] + ycoord_vals[xlength*(j+1)+i] ) * 0.5 + ( ycoord_vals[xlength*j+i+1] + ycoord_vals[xlength*(j+1)+i+1] ) * 0.5 ) * 0.5;
            xcell_bounds[4*(j*xlength+i)]   = ( ( xcoord_vals[j*xlength+i] + xcoord_vals[xlength*j+i-1] ) * 0.5 + ( xcoord_vals[xlength*(j+1)+i] + xcoord_vals[xlength*(j+1)+i-1] ) * 0.5 ) * 0.5;
            xcell_bounds[4*(j*xlength+i)+1] = ( ( xcoord_vals[j*xlength+i] + xcoord_vals[xlength*j+i-1] ) * 0.5 + ( xcoord_vals[xlength*(j-1)+i] + xcoord_vals[xlength*(j-1)+i-1] ) * 0.5 ) * 0.5;
            xcell_bounds[4*(j*xlength+i)+2] = ( ( xcoord_vals[j*xlength+i] + xcoord_vals[xlength*j+i+1] ) * 0.5 + ( xcoord_vals[xlength*(j-1)+i] + xcoord_vals[xlength*(j-1)+i+1] ) *0.5 ) * 0.5;
            xcell_bounds[4*(j*xlength+i)+3] = ( ( xcoord_vals[j*xlength+i] + xcoord_vals[xlength*j+i+1] ) * 0.5 + ( xcoord_vals[xlength*(j+1)+i] + xcoord_vals[xlength*(j+1)+i+1] ) * 0.5 ) * 0.5;
            for ( int m = 0; m < 4; m++ )
              {
                if ( xcell_bounds[4*(j*xlength+i)+m] > 360 ) 
                  xcell_bounds[4*(j*xlength+i)+m] = xcell_bounds[4*(j*xlength+i)+m] - 360.0;
                else if ( xcell_bounds[4*(j*xlength+i)+m] < 0 ) 
                  xcell_bounds[4*(j*xlength+i)+m] = xcell_bounds[4*(j*xlength+i)+m] + 360.0;
              }
            for ( int l = 0; l < 9; l++)
              {
                *star[l] = (*star[l] > 360 ) ? *star[l] - 360.0 : *star[l];
                *star[l] = (*star[l] < 0 )   ? *star[l] + 360.0 : *star[l];
              }
          }

      ycell_bounds[0] = ( ycoord_vals[0] + ycoord_vals[1] ) * 0.5;
      ycell_bounds[1] =   ycoord_vals[0] + ( ycoord_vals[0] - ycoord_vals[1] ) * 0.5;
      ycell_bounds[2] =   ycoord_vals[0] + ( ycoord_vals[0] - ycoord_vals[1] ) * 0.5;
      ycell_bounds[3] = ( ycoord_vals[0] + ycoord_vals[1] ) * 0.5;
      double xcyclicKorr;
      if ( xcyclicKorr =  xcoord_vals[0] - xcoord_vals[xlength] > 270 )
        if ( xcoord_vals[0] < 90 )
           xcyclicKorr -= 
      xcell_bounds[0] =   xcoord_vals[0] + ( xcoord_vals[0] - xcoord_vals[xlength] ) * 0.5;
      xcell_bounds[1] =   xcoord_vals[0] + ( xcoord_vals[0] - xcoord_vals[xlength] ) * 0.5;
      xcell_bounds[2] = ( xcoord_vals[0] + xcoord_vals[xlength] ) * 0.5;
      xcell_bounds[3] = ( xcoord_vals[0] + xcoord_vals[xlength] ) * 0.5;


      ycell_bounds[4*(xlength-1)]   = ( ycoord_vals[xlength-1] +   ycoord_vals[2*xlength-1] ) * 0.5;
      ycell_bounds[4*(xlength-1)+1] =   ycoord_vals[xlength-1] + ( ycoord_vals[xlength-1] - ycoord_vals[2*xlength-1] ) * 0.5;
      ycell_bounds[4*(xlength-1)+2] =   ycoord_vals[xlength-1] + ( ycoord_vals[xlength-1] - ycoord_vals[2*xlength-1] ) * 0.5;
      ycell_bounds[4*(xlength-1)+3] = ( ycoord_vals[xlength-1] +   ycoord_vals[2*xlength-1] ) * 0.5;
      xcell_bounds[4*(xlength-1)]   = ( xcoord_vals[xlength-1] +   xcoord_vals[xlength-2] ) * 0.5;
      xcell_bounds[4*(xlength-1)+1] = ( xcoord_vals[xlength-1] +   xcoord_vals[xlength-2] ) * 0.5;
      xcell_bounds[4*(xlength-1)+2] =   xcoord_vals[xlength-1] + ( xcoord_vals[xlength-1] - xcoord_vals[xlength-2] ) * 0.5;
      xcell_bounds[4*(xlength-1)+3] =   xcoord_vals[xlength-1] + ( xcoord_vals[xlength-1] - xcoord_vals[xlength-2] ) * 0.5;


      ycell_bounds[4*(totalsize-xlength)]   =   ycoord_vals[totalsize-xlength] + ( ycoord_vals[totalsize-xlength] - ycoord_vals[totalsize-2*xlength] ) * 0.5;
      ycell_bounds[4*(totalsize-xlength)+1] = ( ycoord_vals[totalsize-xlength] + ycoord_vals[totalsize-2*xlength] ) * 0.5;
      ycell_bounds[4*(totalsize-xlength)+2] = ( ycoord_vals[totalsize-xlength] + ycoord_vals[totalsize-2*xlength] ) * 0.5;
      ycell_bounds[4*(totalsize-xlength)+3] =   ycoord_vals[totalsize-xlength] + ( ycoord_vals[totalsize-xlength] - ycoord_vals[totalsize-2*xlength] ) * 0.5;
      xcell_bounds[4*(totalsize-xlength)]   =   xcoord_vals[totalsize-xlength] + ( xcoord_vals[totalsize-xlength] - xcoord_vals[totalsize-xlength+1] ) * 0.5;
      xcell_bounds[4*(totalsize-xlength)+1] =   xcoord_vals[totalsize-xlength] + ( xcoord_vals[totalsize-xlength] - xcoord_vals[totalsize-xlength+1] ) * 0.5;
      xcell_bounds[4*(totalsize-xlength)+2] = ( xcoord_vals[totalsize-xlength] + xcoord_vals[totalsize-xlength+1] ) * 0.5;
      xcell_bounds[4*(totalsize-xlength)+3] = ( xcoord_vals[totalsize-xlength] + xcoord_vals[totalsize-xlength+1] ) * 0.5;


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

          ycell_bounds[4*i]   = ( ( ycoord_vals[i] + ycoord_vals[i+xlength] ) * 0.5 + ( ycoord_vals[i-1] + ycoord_vals[i+xlength-1] ) * 0.5 ) * 0.5;
          ycell_bounds[4*i+1] =     ycoord_vals[i] + ( ycoord_vals[i] - ycoord_vals[i+xlength] ) * 0.5;
          ycell_bounds[4*i+2] =     ycoord_vals[i] + ( ycoord_vals[i] - ycoord_vals[i+xlength] ) * 0.5;
          ycell_bounds[4*i+3] = ( ( ycoord_vals[i] + ycoord_vals[i+xlength] ) * 0.5 + ( ycoord_vals[i+1] + ycoord_vals[i+xlength+1] ) * 0.5 ) * 0.5;
          xcell_bounds[4*i]   = ( ( xcoord_vals[i] + xcoord_vals[i-1] ) * 0.5 + ( xcoord_vals[i+xlength] + xcoord_vals[i+xlength-1] ) * 0.5 ) * 0.5;
          xcell_bounds[4*i+1] =   ( xcoord_vals[i] + xcoord_vals[i-1] ) * 0.5;
          xcell_bounds[4*i+2] =   ( xcoord_vals[i] + xcoord_vals[i+1] ) * 0.5;
          xcell_bounds[4*i+3] = ( ( xcoord_vals[i] + xcoord_vals[i+1] ) * 0.5 + ( xcoord_vals[i+xlength] + xcoord_vals[i+xlength+1] ) * 0.5 ) * 0.5;


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

          ycell_bounds[4*j*xlength]   = (   ycoord_vals[j*xlength] + ycoord_vals[(j+1)*xlength] ) * 0.5;
          ycell_bounds[4*j*xlength+1] = (   ycoord_vals[j*xlength] + ycoord_vals[(j-1)*xlength] ) * 0.5;
          ycell_bounds[4*j*xlength+2] = ( ( ycoord_vals[j*xlength] + ycoord_vals[(j-1)*xlength] ) * 0.5 + ( ycoord_vals[j*xlength+1] + ycoord_vals[(j-1)*xlength+1] ) * 0.5 ) * 0.5;
          ycell_bounds[4*j*xlength+3] = ( ( ycoord_vals[j*xlength] + ycoord_vals[(j+1)*xlength] ) * 0.5 + ( ycoord_vals[j*xlength+1] + ycoord_vals[(j+1)*xlength+1] ) * 0.5 ) * 0.5;
          xcell_bounds[4*j*xlength]   =     xcoord_vals[j*xlength] + ( xcoord_vals[j*xlength] - xcoord_vals[j*xlength+1] ) * 0.5;
          xcell_bounds[4*j*xlength+1] =     xcoord_vals[j*xlength] + ( xcoord_vals[j*xlength] - xcoord_vals[j*xlength+1] ) * 0.5; 
          xcell_bounds[4*j*xlength+2] = ( ( xcoord_vals[j*xlength] + xcoord_vals[j*xlength+1] ) * 0.5 + ( xcoord_vals[(j-1)*xlength] + xcoord_vals[(j-1)*xlength+1] ) * 0.5 ) * 0.5;
          xcell_bounds[4*j*xlength+3] = ( ( xcoord_vals[j*xlength] + xcoord_vals[j*xlength+1] ) * 0.5 + ( xcoord_vals[(j+1)*xlength] + xcoord_vals[(j+1)*xlength+1] ) * 0.5 ) * 0.5;


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
} */

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
static void register_lon_axis(int gridID, int xlength, int *axis_ids)
{
  double *xcoord_vals = Malloc(xlength * sizeof(double));
  if ( gridInqXvals(gridID, xcoord_vals) == 0 )
    Free(xcoord_vals);
  else
    {
      double *xcell_bounds = Malloc(2 * xlength * sizeof(double));
      int xnbounds = gridInqXbounds(gridID, xcell_bounds);
      check_and_gen_bounds(gridID, xnbounds, xlength, xcoord_vals, xcell_bounds, 1);
      cmor_axis(new_axis_id(axis_ids),    "longitude",    "degrees_east",    xlength,    (void *)xcoord_vals,    'd',    (void *)xcell_bounds,    2,    NULL);
      if ( xcell_bounds ) Free(xcell_bounds);
      if ( xcoord_vals ) Free(xcoord_vals);
    }
}

static void register_lat_axis(int gridID, int ylength, int *axis_ids)
{
  double *ycoord_vals = Malloc(ylength * sizeof(double));
  if ( gridInqYvals(gridID, ycoord_vals) == 0 )
    Free(ycoord_vals);
  else
    {
      double *ycell_bounds = Malloc(2 * ylength * sizeof(double));
      int ynbounds = gridInqYbounds(gridID, ycell_bounds);
      check_and_gen_bounds(gridID, ynbounds, ylength, ycoord_vals, ycell_bounds, 0);
      cmor_axis(new_axis_id(axis_ids),    "latitude",    "degrees_north",    ylength,    (void *)ycoord_vals,    'd',    (void *)ycell_bounds,    2,    NULL);
      if ( ycell_bounds ) Free(ycell_bounds);
      if ( ycoord_vals ) Free(ycoord_vals);  
    }
}

static void register_char_axis(int numchar, char **charvals, int *axis_ids, char *chardim)
{
  int maxlen = get_strmaxlen(charvals, numchar);
  void *charcmor = (void *) Malloc ( numchar * maxlen * sizeof(char));
  sprintf((char *)charcmor, "%.*s", strlen(charvals[0]), charvals[0]);
  char blanks[maxlen];
  for ( int i = 0; i < maxlen; i++)
    blanks[i] = ' ';
  sprintf((char *)charcmor, "%s%.*s", (char *)charcmor, maxlen-strlen(charvals[0]), blanks);
  for ( int i = 1; i < numchar; i++ )
    {
      sprintf((char *)charcmor, "%s%s", (char *)charcmor, charvals[i]);
      sprintf((char *)charcmor, "%s%.*s", (char *)charcmor, maxlen-strlen(charvals[i]), blanks);   
    }
  cmor_axis(new_axis_id(axis_ids), chardim, "", numchar, (void *)charcmor, 'c',  NULL, maxlen, NULL); 
  Free(charcmor);
}

static void register_grid(list_t *kvl, int vlistID, int varID, int *axis_ids, int *grid_ids, char *project_id)
{
  int gridID = vlistInqVarGrid(vlistID, varID);

  char *grid_file = kv_get_a_val(kvl, "gi", NULL);

  char *chardimatt = kv_get_a_val(kvl, "ca", NULL);
  char *chardim = get_txtatt(vlistID, varID, "character_axis");
  check_compare_set(&chardim, chardimatt, "character_axis", "notSet");

  if ( grid_file )
    {
      change_grid(grid_file, gridID, vlistID);
      gridID = vlistInqVarGrid(vlistID, varID);
    }

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
      get_cmor_table(kvl, project_id);
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
  else if ( type == GRID_GENERIC && ( strcmp(chardim, "oline") == 0 || strcmp(chardim, "basin") == 0 ))
    {
      if ( cdoVerbose )
        printf("*******Start to define a character axis '%s' instead of a grid axis'.******\n", chardim);
      grid_ids[0] = 0;
      int numchar = 0;
      char *charvalstring = Malloc(CMOR_MAX_STRING * sizeof(char));
      sprintf(charvalstring, "char_axis_%s", chardim);
      char **charvals = kv_get_vals(kvl, charvalstring, &numchar);
      Free(charvalstring);
      if ( ( xlength > 0 && xlength != numchar ) && ( ylength > 0 && ylength != numchar ) )
        cdoAbort("You configured a character coordinate '%s' with '%d' string values but you also registered a grid with '%d' numerical values on X axis and '%d' numerical values on Y axis. Both is not supported!", chardim, numchar, xlength, ylength);
      if ( !charvals )
        cdoAbort("You configured a character coordinate '%s' but no values are found! Configure values via attribute 'char_dim_vals'!", chardim);
      if ( charvals && ( xlength == numchar || xlength == 0 ) )
        {
          register_char_axis(numchar, charvals, axis_ids, chardim);
          if ( ylength > 0 )
            register_lat_axis(gridID, ylength, axis_ids);
        }
      else
        {
          register_lon_axis(gridID, xlength, axis_ids);
          register_char_axis(numchar, charvals, axis_ids, chardim);
        }
      if ( cdoVerbose )
        printf("*******Succesfully defined a character axis '%s' instead of a grid axis.******\n", chardim);
      Free(chardim);
    }
/*
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
      if ( strcmp(xdimname, "line") == 0 )
        strcpy(xdimname, "oline");
      int dimstrlen;   
      if ( dimstrlen = gridInqXIsc(gridID) )
        {
          char **xchars = (char **)Malloc( (xlength+1) * sizeof(char *));
          for ( int i = 0; i < xlength; i++ )
            xchars[i] = (char *)Malloc( (dimstrlen+1) * sizeof(char));
          gridInqXCvals(gridID, xchars);
          for ( int j = 0; j < xlength; j++ )
            xchars[j][dimstrlen] = 0;
          xchars[xlength] = NULL;
          register_char_axis(xlength, xchars, axis_ids, xdimname);
          free_array(xchars);
        }
      else if ( xlength)
        register_lon_axis(gridID, xlength, axis_ids);

      if ( dimstrlen = gridInqYIsc(gridID) )
        {
          char **ychars = (char **) Malloc( (ylength + 1) * sizeof(char));
          for ( int i = 0; i < ylength; i++ )
            ychars[i] = (char *)Malloc( (dimstrlen +1) * sizeof(char));
          gridInqYCvals(gridID, ychars);
          for ( int j = 0; j < ylength; j++ )
            ychars[j][dimstrlen] = 0;
          ychars[ylength] = NULL;
          register_char_axis(ylength, ychars, axis_ids, ydimname);
          free_array(ychars);
        }
      else if ( ylength )
        register_lat_axis(gridID, ylength, axis_ids);
      Free(xname); Free(yname); Free(xdimname); Free(ydimname);
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
                              struct mapping *var, int *grid_ids, char *name)
{
  if ( cdoVerbose )
    printf("*******Start to retrieve 'positive' and 'units'.******\n");
  char *positive = get_txtatt(vlistID, varID, "positive");
  char *origname = get_txtatt(vlistID, varID, "original_name");
  char *history = get_txtatt(vlistID, varID, "history");
  char *units = Malloc(CDI_MAX_NAME * sizeof(char));
  vlistInqVarUnits(vlistID, varID, units);
  char *attunits = kv_get_a_val(kvl, "u", NULL);
  char *attp = kv_get_a_val(kvl, "p", NULL);
  char *attorigname = kv_get_a_val(kvl, "original_name", NULL);
  check_compare_set(&positive, attp, "positive", "");
  if ( strcmp(positive, " ") == 0 )
    strcpy(positive, "");
  check_compare_set(&units, attunits, "units", NULL);
  check_compare_set(&origname, attorigname, "original_name", "");
  if ( strcmp(origname, "") == 0 || strstr(origname, "var") )
    origname = NULL;
  if ( cdoVerbose )
    printf("*******Succesfully retrieved 'positive': '%s' and 'units' : '%s'.******\n", positive, units);
  char missing_value[sizeof(double)];
  double tolerance = 1e-4;
  size_t gridsize = vlistGridsizeMax(vlistID);
  int zsize = zaxisInqSize(vlistInqVarZaxis(vlistID, varID));
  var->cdi_varID = varID;
  var->help_var = 0;
  if ( !var->data )
    {
      var->charvars = 0;
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
    }
  else
    *(double *) missing_value = vlistInqVarMissval(vlistID, varID);
   
  if ( cdoVerbose )
    printf("*******Start to call cmor_variable.******\n");
  if ( grid_ids[0] != 0 )
    {
      int *tmp_id = new_axis_id(axis_ids);
      *tmp_id = grid_ids[0];
      cmor_variable(&var->cmor_varID,
            name,units,(count_axis_ids(axis_ids)), axis_ids, var->datatype,
            (void *) missing_value, &tolerance, positive,
                        origname,
                        history,
                        kv_get_a_val(kvl, "vc", NULL));
    }
  else
    {
      cmor_variable(&var->cmor_varID,
           name, units, count_axis_ids(axis_ids),  axis_ids,   var->datatype,
          (void *) missing_value, &tolerance, positive,
                        origname,
                        history,
                        kv_get_a_val(kvl, "vc", NULL));
    }
  if ( cdoVerbose )
    printf("*******Succesfully called cmor_variable.******\n");
  if (positive) Free(positive); 
  if (units) Free(units);
}

static void register_all_dimensions(list_t *kvl, int streamID,
                             struct mapping vars[], int table_id, char *project_id, int miptab_freq, int *time_axis)
{
  int vlistID = streamInqVlist(streamID);
  int taxisID = vlistInqTaxis(vlistID);

  if ( cdoVerbose )
    printf("*******Start to check attribute 'required_time_units'.******\n");
  char *time_units = get_time_units(taxisID);
  char *required_time_units = kv_get_a_val(kvl, "required_time_units", NULL);
  if ( check_time_units(required_time_units) )
    check_compare_set(&time_units, required_time_units, "time_units", NULL);
  else 
    cdoAbort("Required Attribute 'required_time_units' from configuration is invalid!");
  if ( cdoVerbose )
    printf("*******Succesfully checked attribute 'required_time_units'.*******\n");

  if ( cdoVerbose )
    printf("*******Start to retrieve requested variables.******\n");

  int numvals = 0;
  char **cmor_names = kv_get_vals(kvl, "cn", &numvals);

/* Cmdlinemapping: */
  char *mapname, *mapcode;
  if ( kv_get_a_val(kvl, "mt", NULL) && numvals )
    {
      if ( mapname = kv_get_a_val(kvl, "n", NULL) )
        change_name_via_name(vlistID, mapname, cmor_names[0]);
      else if ( mapcode = kv_get_a_val(kvl, "c", NULL) )
        change_name_via_code(vlistID, mapcode, cmor_names[0]);
    }

  if ( cmor_names == NULL && vlistNvars(vlistID) > 1 )
    cdoWarning("You have not requested any specific variable but there are several in input! Notice that all command line configuration attributes including cmor_name and units will be used for every variable!\n");
  if ( cdoVerbose )
    printf("*******Succesfully retrieved requested variables*******\n");
  int foundName = 0;
  int ps_required = 0;
  int ps_in_file = 0;
  for ( int varID = 0; varID < vlistNvars(vlistID); varID++ )
    {
      char name[CDI_MAX_NAME];
      vlistInqVarName(vlistID, varID, name);
      if ( !cmor_names || in_list(cmor_names, name, numvals) )
        {
          struct mapping *var = map_var(varID, vars);
          if ( !var )
            var = new_var_mapping(vars);
          int axis_ids[CMOR_MAX_AXES];
          axis_ids[0] = CMOR_UNDEFID;
          int zaxisID = vlistInqVarZaxis(vlistID, varID);
          if ( cdoVerbose )
            printf("*******Start to define variable with ID: '%d' and name: '%s'*******\n", varID, name);
          if ( zaxisInqType(zaxisID) == ZAXIS_HYBRID )
            {
              if ( cdoVerbose )
                printf("Since the zaxis of variable '%s' is of type HYBRID, surface pressure is required. An infile variable must have the name ps.\n", name);
              ps_required++;
            }
          foundName++;
          /* Time-Axis */
          if ( cdoVerbose )
            printf("*******Start to register a time axis*******\n");
          char cmor_time_name[CMOR_MAX_STRING];
          get_time_method(kvl, vlistID, varID, cmor_time_name, project_id, miptab_freq, time_axis);
          if ( strcmp(cmor_time_name, "none") != 0 )
            cmor_axis(new_axis_id(axis_ids),
                    cmor_time_name,
                    time_units,
                    0,NULL, 0, NULL, 0, NULL);
          if ( cdoVerbose )
            printf("*******Succesfully handled time axis registration*******\n");
          /* Grid: */
          if ( cdoVerbose )
            printf("*******Start to register a grid*******\n");
          int grid_ids[CMOR_MAX_GRIDS];
          register_grid(kvl, vlistID, varID, axis_ids, grid_ids, project_id);
          cmor_set_table(table_id);
          if ( cdoVerbose )
            printf("*******Succesfully handled grid registration*******\n");
          /* Z-Axis */
          if ( cdoVerbose )
            printf("*******Start to register a zaxis*******\n");
          register_z_axis(kvl, vlistID, varID, zaxisID, name, axis_ids, &var->zfactor_id, project_id, miptab_freq);
          if ( cdoVerbose )
            printf("*******Succesfully handled zaxis registration*******\n");
          /* Variable */
          register_variable(kvl, vlistID, varID, axis_ids, var, grid_ids, name);     
          if ( cdoVerbose )
            printf("*******Succesfully defined variable with ID: '%d' and name: '%s'*******\n", varID, name);
        }
    }
  if ( ps_required )
    {
      if ( cdoVerbose )
        printf("\n *******Start to find surface pressure.*******\\n");
      for ( int varID = 0; varID < vlistNvars(vlistID); varID++ )
        if ( vlistInqVarCode(vlistID, varID) == 134 )
          {
            ps_in_file++;
            if ( cmor_names == NULL || in_list(cmor_names, "ps", numvals) )
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
      if ( cdoVerbose )
        printf("\n *******Succesfully registered surface pressure.*******\\n");
    }
  if ( ps_required && !ps_in_file )
    cdoAbort("No surface pressure found in Ifile but required for a hybrid sigma pressure z axis!");
  if ( !foundName && cmor_names == NULL )
    cdoAbort("No variables from your table %s found in Ifile.\n");
  if ( !foundName && cmor_names )
    cdoAbort("None of the given variables to process by attribute 'cmor_name' is found in Ifile.\n");
  if ( time_units) Free(time_units);
  if ( cdoVerbose )
    printf("*******Succesfully registered all dimensions for %d variables successfully.*******\n", foundName);
}

static char *get_frequency(list_t *kvl, int streamID, int vlistID, int taxisID, int miptab_freq)
{
  char *frequency = Malloc(CMOR_MAX_STRING * sizeof(char));
  int ntsteps = vlistNtsteps(vlistID);
  int reccounter = 0;
  int recdummy = 0;

  switch ( miptab_freq )
    {
    case 11: strcpy(frequency, "yr"); break;
    case 2: strcpy(frequency, "yr"); break;
    case 12: strcpy(frequency, "mon"); break;
    case 3: strcpy(frequency, "mon"); break;
    case 13: strcpy(frequency, "day"); break;
    case 4: strcpy(frequency, "day"); break;
    case 14: strcpy(frequency, "6hr"); break;
    case 5: strcpy(frequency, "6hr"); break;
    case 6: strcpy(frequency, "6hr"); break;
    case 15: strcpy(frequency, "3hr"); break;
    default:
    {
      if ( cdoStreamName(0)->args[0] == '-' )
        {
            cdoAbort("No frequency could be determined from MIP-table and cdo cmor cannot check frequency of Ifile recs since you piped several cdo operators.");
/*          char *dummy;
          cdoWarning("Cdo cmor cannot check frequency of Ifile recs since you piped several cdo operators.\nIt is tried to use a configuration attribute frequency.");
          if ( !(dummy = kv_get_a_val(kvl, "frequency", NULL)) )
            cdoAbort("No attribute frequency is found.");
          else
            {
              strcpy(frequency, dummy);
              return frequency;
            }
*/
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
      streamClose(streamID2);
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


static void read_record(int streamID, struct mapping vars[], int vlistID)
{
  int varID, levelID;
  streamInqRecord(streamID, &varID, &levelID);

  int gridID = vlistInqVarGrid(vlistID, varID);
  int type = gridInqType(gridID);
  int gridsize = gridInqSize(gridID);
  double *buffer = (double *) Malloc(gridsize * sizeof(double));

  struct mapping *var = map_var(varID, vars);
  if ( var && var->charvars != 1 )
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
  Free(buffer);
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


static int check_append_and_size(list_t *kvl, int vlistID, char *testIn, int ifreq, int calendar)
{
  char *test = testIn;
  size_t filesize = fileSize((const char *)testIn);
  char old_start_date[CMOR_MAX_STRING];
  char old_end_date[CMOR_MAX_STRING];
  int i = 0, j = 0;
/* Get dates from chunk string */
  if ( cdoVerbose) printf("*******Start to retrieve dates from chunk string.******\n");
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

  if ( cdoVerbose) printf("*******Succesfully retrieved start date: '%s' and end date: '%s' chunk string.******\n", old_start_date, old_end_date);
/* Check frequency of chunk with frequency of file */

  if ( (j == 8 && ifreq !=3) || (ifreq == 3 && j != 8)
    || (j == 6 && ifreq !=2) || (ifreq == 2 && j != 6)
    || (j == 4 && ifreq !=1) || (ifreq == 1 && j != 4) )
    cdoAbort("Frequency of chunk file does not agree with frequency of the working file.");

/* Encode in julseconds depending on frequency */
  if ( cdoVerbose) printf("*******Start to encode dates with frequencies to julseconds.******\n");

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

  int cdi_startdate = cdiEncodeDate(old_start_year, old_start_month, old_start_day);
  int cdi_enddate = cdiEncodeDate(old_end_year, old_end_month, old_end_day);
  int cdi_time = cdiEncodeTime(0, 0, 0);
  juldate_t julostart = juldate_encode(calendar, cdi_startdate, cdi_time);
  juldate_t juloend = juldate_encode(calendar, cdi_enddate, cdi_time);

  if ( cdoVerbose) printf("*******Succesfully calculated juldates.******\n", old_start_date);
/* Read in first vdate in case not piped */
  if ( cdoVerbose) printf("*******Start to calculate temporal gap between chunk and working file.******\n");
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

  if ( ( j == 8 && ( append_distance > 48.0 || append_distance < 0 ) )
     ||( j == 6 && ( append_distance/24.0 > 62.0 || append_distance < 0 ) )
     ||( j == 4 && ( append_distance/24.0/30.5 > 24.0 || append_distance < 0 ) ) )
    {
      cdoWarning("A temporal gap is diagnosed between end date of chunk file and first date of working file of: '%f' hours. Maximal valid gaps are:\n48 hours for daily frequency\n62 days for monthly frequency\n24 month for yearly frequency.\nSwitched to replace mode.", append_distance);
      streamClose(streamID2);
      return 0;
    }

  if ( cdoVerbose) printf("*******Succesfully checked temporal gap.******\n");
/* Check file size */
  if ( cdoVerbose) printf("*******Start to check file size of chunk + working file.******\n");
  double old_interval_sec = juldate_to_seconds(juldate_sub(juloend, julostart));
  double size_per_sec = (double) filesize / old_interval_sec;

  int maxsizegb = atol(kv_get_a_val(kvl, "ms", "2"));
  int maxsizeb = maxsizegb * 1024 * 1024 * 1024;

  int ntsteps = vlistNtsteps(vlistID2);
  if ( ntsteps < 0 )
    {
      ntsteps = 0;
      while ( streamInqTimestep(streamID2, ntsteps++)) ;
      if ( ntsteps == 0 )
        {
          cdoWarning("A mistake occured during timesteps determination.\nSwitched to replace mode.");
          streamClose(streamID2);
          return 0;
        }
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
      {
        cdoWarning("Selected chunk to append data has subdaily frequency which is yet not enabled by cdo cmor.\nA new file will be written.");
        streamClose(streamID2);
        return 0;
      }
    }

  if ( (unsigned int)estimated_size > (unsigned int) maxsizeb )
    {
      cdoWarning("Estimated file size of appended file is : '%f'gb and exceeds maximal allowed file size: '%d'gb.\nA new file will be written.", estimated_size/1024.0/1024.0/1024.0, maxsizegb);
      streamClose(streamID2);
      return 0;
    }
  streamClose(streamID2);
  if ( cdoVerbose) printf("*******Succesfully checked file size of chunk + working file.******\n");
  return 1;
}

static char *use_chunk_des_files(list_t *kvl, int vlistID, int var_id, char *chunk_des_file, int ifreq, int calendar)
{
  char *chunk_file = Malloc(4096 * sizeof(char));
  if ( file_exist(chunk_des_file, 0) )
    {
      FILE *fp = fopen(chunk_des_file, "r");
      size_t filesize = fileSize(chunk_des_file);
      char *buffer = (char*) Malloc(filesize);
      size_t nitems = fread(buffer, 1, filesize, fp);
      buffer = readLineFromBuffer(buffer, &filesize, chunk_file, 4096);
      fclose(fp);
      if ( file_exist(chunk_file, 0) && check_append_and_size(kvl, vlistID, chunk_file, ifreq, calendar) )
        return chunk_file;
      else
        cdoWarning("Chunk '%s' configured via chunk description file could either not be opened or is not suitable to be appended.\nSwitched to replace mode.", chunk_file);
    }
  else
    cdoWarning("Chunk description file '%s' could not be opened.\nSwitched to replace mode.", chunk_des_file);
  strcpy(chunk_file, " \0");
  return chunk_file;
}

static char **empty_array(struct mapping vars[], char ***chunk_files)
{
  for ( int i = 0; vars[i].cdi_varID != CDI_UNDEFID; i++ )
    (*chunk_files)[i] = NULL;
  return *chunk_files;
}

static char **get_chunk_des_files(list_t *kvl, struct mapping vars[], char *miptab_freqptr, int nreq, int vlistID, char *charname)
{
  char **chunk_des_files = Malloc((nreq+1) * sizeof(char *));
  chunk_des_files[nreq] = NULL;

  char *trunk = Malloc(CMOR_MAX_STRING * sizeof(char));
  char *description_atts[] = {"model_id", "experiment_id", "member", NULL};
  strcpy(trunk, miptab_freqptr);
  for ( int i = 0; description_atts[i]; i++ )
    {
      strcat(trunk, "_");
      strcat(trunk, kv_get_a_val(kvl, description_atts[i], ""));
    }

  for ( int j = 0; vars[j].cdi_varID != CDI_UNDEFID; j++)
    {
      char *name = Malloc(CDI_MAX_NAME * sizeof(char));
      if ( charname )
        strcpy(name, charname);
      else
        vlistInqVarName(vlistID, vars[j].cdi_varID, name);
      chunk_des_files[j] = Malloc(CMOR_MAX_STRING * sizeof(char));
      sprintf(chunk_des_files[j], "CHUNK_FILE_%s_%s.txt\0", name, trunk);
      Free(name);
    }
  return chunk_des_files;  
}

static char **get_chunk_files(list_t *kvl, struct mapping vars[], int vlistID, int ifreq, int time_axis, int calendar, char *miptab_freqptr)
{
  int i = 0;
  for ( i = 0; vars[i].cdi_varID != CDI_UNDEFID; i++ );
  char **chunk_files = Malloc((i+1) * sizeof(char *));
  chunk_files[i] = NULL;
  
  char *dummy = kv_get_a_val(kvl, "om", NULL);
  if ( !dummy || strcmp(dummy, "a") != 0 )
    return empty_array(vars, &chunk_files);
  else if ( time_axis == 3 )
    {
      printf("CMOR APPEND mode not possible for time independent variables.\nSwitched to replace mode");
      return empty_array(vars, &chunk_files);
    }

  if ( cdoVerbose )
    printf("\n*******Start to retrieve chunk files to append .******\n");

  int num_aaf = 0;
  char **chunk_att_files = kv_get_vals(kvl, "lc", &num_aaf);
  char **chunk_des_files = NULL;
  if ( num_aaf != i && num_aaf > 0 )
    {
      printf("Number of chunk files '%d' disagree with number of requested variables '%d'.\n Switched to replace mode.\n", num_aaf, i); 
      return empty_array(vars, &chunk_files);
    }  
  else if ( num_aaf == 0 )
    {
      char *nd = kv_get_a_val(kvl, "d", "y");
/* For chunk description file : */
      if ( nd[0] == 'y' )
        chunk_des_files = get_chunk_des_files(kvl, vars, miptab_freqptr, i, vlistID, NULL);
      else if ( cdoVerbose )
        {
          printf("Automatic chunk configuration via file not possible if DRS is not created.\nSwichted to replace mode.");
          return empty_array(vars, &chunk_files);
        }
    }

  for ( int j = 0; vars[j].cdi_varID != CDI_UNDEFID; j++ )
    {
      if ( num_aaf != 0 )
        {
          if ( file_exist(chunk_att_files[j], 0) && check_append_and_size(kvl, vlistID, chunk_att_files[j], ifreq, calendar) )
            chunk_files[j] = strdup(chunk_att_files[j]);
          else
            {
              cdoWarning("Chunk '%s' could not be used.\nSwitched to replace mode for this variable.\n", chunk_att_files[j]);
              chunk_files[j] = strdup(" ");
            }   
        }
      else 
        {
          if ( cdoVerbose )
            printf("It is tried to open a chunk description file for varID: '%d': '%s'.\n", vars[j].cdi_varID, chunk_des_files[j]);
          chunk_files[j] = use_chunk_des_files(kvl, vlistID, vars[j].cdi_varID, chunk_des_files[j], ifreq, calendar);  
        }
      if ( cdoVerbose && strcmp(chunk_files[j], " ") != 0 )
        printf("\n*******Chunk file to append on var with CDI ID %d is: '%s' .******\n", vars[j].cdi_varID, chunk_files[j]);
    }
  if ( chunk_des_files ) free_array(chunk_des_files);
  if ( cdoVerbose )
    printf("\n*******Successfully processed chunk file retrieval.******\n");
  return chunk_files;
}

static void write_variables(list_t *kvl, int *streamID, struct mapping vars[], int miptab_freq, int time_axis, int calendar, char *miptab_freqptr)
{
  int vlistID = streamInqVlist(*streamID);
  int taxisID = vlistInqTaxis(vlistID);
  int tsID = 0;
  int nrecs;
  size_t gridsize = vlistGridsizeMax(vlistID);

  if ( cdoVerbose )
    printf("\n*******Start to retrieve relative start time value from 'required_time_units' and file and start to retrieve frequency.******\n");
  int sdate, stime, time_unit;
  get_taxis(kv_get_a_val(kvl, "rtu", NULL), &sdate, &stime, &time_unit);
  int tunitsec = get_tunitsec(time_unit);
  juldate_t ref_date = juldate_encode(calendar, sdate, stime);
  char *frequency = NULL;
  if ( time_axis != 3 )
    frequency = get_frequency(kvl, *streamID, vlistID, taxisID, miptab_freq);
  if ( cdoVerbose )
    printf("\n*******Succesfully retrieved time value from 'required_time_units' and frequency.******\n");

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

  char **chunk_files = get_chunk_files(kvl, vars, vlistID, ifreq, time_axis, calendar, miptab_freqptr);

  if ( cdoVerbose )
    printf("\n*******Start to write variables via cmor_write.******\n");
  int i = 0;

  int zaxisID, zsize, pscheck = 1;
  char *charname = NULL;
  for ( i = 0; vars[i].cdi_varID != CDI_UNDEFID; i++ )
    if ( vars[i].charvars )
      {
        zaxisID = vlistInqVarZaxis(vlistID, vars[i].cdi_varID);
        zsize = zaxisInqSize(zaxisID);
        charname = Malloc(CDI_MAX_NAME * sizeof(char));
        vlistInqVarName(vlistID, vars[i].cdi_varID, charname);
        
        streamClose(*streamID);
        *streamID = streamOpenRead(cdoStreamName(0));
        pscheck = 0;
        break;
      }
  if ( !pscheck )
    cdoWarning("Since you defined a variable with character coordinate axis you cannot write another variable with zaxis of type ZAXIS_HYBRID.");

  while ( (nrecs = streamInqTimestep(*streamID, tsID++)) )
    { 
      double time_bnds[2];
      double *time_bndsp;
      double time_val;
      if ( time_axis != 3 )
        {
          time_val = get_cmor_time_val(taxisID, ref_date, tunitsec, calendar);
          time_bndsp = ( time_axis != 1 ) ? get_time_bounds(taxisID, frequency, ref_date, time_val, calendar, tunitsec, time_bnds) : 0;
        }
      while ( nrecs-- )
        read_record(*streamID, vars, vlistID);

      int ps_index = -1;
      if ( pscheck )
        check_for_sfc_pressure(&ps_index, vars, vlistID, tsID);
      for ( i = 0; vars[i].cdi_varID != CDI_UNDEFID; i++ )
        {
/*          char name[CDI_MAX_NAME];
          vlistInqVarName(vlistID, vars[i].cdi_varID, name); */
          if ( !vars[i].help_var )
            {
              if ( time_axis != 3 )
                {
                  if ( vars[i].charvars )
                    {
                      void *dataslice = Malloc(gridsize * zsize * sizeof(double));
                      for ( int j = 0; j < gridsize * zsize; j++ )
                        ((double *)dataslice)[j] = ((double *)vars[i].data)[(tsID-1)*gridsize*zsize+j];
                      cmor_write(vars[i].cmor_varID,
                       dataslice,
                       vars[i].datatype,
                       chunk_files[i],
                       1,
                       &time_val,
                       time_bndsp,
                       NULL);
                      Free(dataslice);
                    } 
                  else
                    cmor_write(vars[i].cmor_varID,
                     vars[i].data,
                     vars[i].datatype,
                     chunk_files[i],
                     1,
                     &time_val,
                     time_bndsp,
                     NULL); 
                  if ( vars[i].zfactor_id > 0 )
                    cmor_write(vars[i].zfactor_id,
                       vars[ps_index].data,
                       vars[ps_index].datatype,
                       chunk_files[i],
                       1,
                       &time_val,
                       time_bndsp,
                       &vars[i].cmor_varID);
                }
              else
                cmor_write(vars[i].cmor_varID,
                   vars[i].data,
                   vars[i].datatype,
                   chunk_files[i], 0, 0, 0, NULL);
            }
        }
    }
  if ( cdoVerbose )
    printf("\n*******Succesfully written variables via cmor_write.******\n");
  if ( cdoVerbose )
    printf("\n*******Start to close files and free allocated memory.******\n");
  char **chunkdf = NULL;
  if ( strcmp(kv_get_a_val(kvl, "om", ""), "a") == 0 && strcmp(kv_get_a_val(kvl, "d", "y"), "y") == 0 )
    chunkdf = get_chunk_des_files(kvl, vars, miptab_freqptr, i, vlistID, charname);

  char file_name[CMOR_MAX_STRING];
  for ( i = 0; vars[i].cdi_varID != CDI_UNDEFID; i++ )
    {
      if ( !vars[i].help_var )
        {
          cmor_close_variable(vars[i].cmor_varID, file_name, NULL);
          printf("*******File stored in:  '%s' with cmor!*******\n", file_name);
          if ( chunkdf )
            {
              if ( cdoVerbose )
                printf("*******Start to write a chunk description file.******\n");
              FILE *fp = fopen(chunkdf[i], "w+"); 
              if ( fp )
                fprintf(fp, "%s", file_name);
              else
                {
                  if ( cdoVerbose )
                    printf("Could not open a chunk description file '%s'.\n", chunkdf[i]);
                  continue;
                }
              fclose(fp);  
              if ( cdoVerbose )
                printf("*******Succesfully written a chunk description file '%s'******\n" , chunkdf[i]);            
            }         
        }
    }


  if (frequency) Free(frequency); if (chunk_files) free_array(chunk_files); if (chunkdf) free_array(chunkdf); if (charname) Free(charname);
  if ( cdoVerbose )
    printf("\n*******Succesfully closed files and freed allocated memory.******\n");
}

static list_t *check_for_charvars(list_t *maptab, char *key)
{
  listNode_t *node = maptab->head;
  while ( node )
    {
      if ( node->data )
        {
          list_t *kvlist = *(list_t **)node->data;
          keyValues_t *kvn = NULL;
          if ( key )
            kvn = kvlist_search(kvlist, key);
          else
            {
              kvn = kvlist_search(kvlist, "name");
              if ( !kvn )
                kvn = kvlist_search(kvlist, "code");
            }
          if ( kvn && kvn->nvalues > 1 )
            return kvlist;
          if ( kvn && strstr(kvn->values[0], ",") && kvn->nvalues == 1 )
            {
              char *workchar = strdup(kvn->values[0]);
              Free(kvn->values[0]); Free(kvn->values);
              char *thepoint = workchar;
              int i = 0, j = 0;
              while ( *thepoint != '\0' )
                {
                  thepoint++;
                  if ( *thepoint == ',' )
                    j++;
                }
              j++;
              kvn->nvalues = j;
              kvn->values = (char **) malloc(kvn->nvalues*sizeof(char*)); 

              j = 0; thepoint = workchar;             
              while ( *thepoint != '\0' )
                {
                  if ( *thepoint == ',')
                    {
                      kvn->values[j] = Malloc( (i+1) * sizeof(char) );
                      strncpy(kvn->values[j], workchar, i);
                      kvn->values[j][i] = '\0';
                      j++; thepoint++; workchar+=i+1; i = 0;
                    }
                  else
                    {
                      thepoint++; i++;
                    }
                }
              if ( i > 0 )
                {
                  kvn->values[j] = Malloc( (i+1) * sizeof(char) );
                  strncpy(kvn->values[j], workchar, i);
                  kvn->values[j][i] = '\0';
                  workchar+=i; i = 0; j++; 
                }
              else
                {
                  cdoWarning("Names in String for key '%s' could not be interpreted correctly due to a comma at end of line.");
                  return NULL;
                }
              return kvlist;
            }
        }
      node = node->next;
    } 
  return NULL;
}

static void read_maptab(list_t *kvl, int streamID, char *miptabfreq, struct mapping vars[])
{
  char *maptab = kv_get_a_val(kvl, "mt", NULL);
  char *maptabdir = kv_get_a_val(kvl, "mapping_table_dir", NULL);
  char *maptabbuild = NULL;
  keyValues_t *kvn = kvlist_search(kvl, "n");
  keyValues_t *kvc = kvlist_search(kvl, "c");
  keyValues_t *kvcn = kvlist_search(kvl, "cn");

  if ( maptab && maptabdir ) if ( maptab[0] != '/' )
    {
      maptabbuild = Malloc((strlen(maptab)+strlen(maptabdir)+2) * sizeof(char));
      sprintf(maptabbuild, "%s/%s\0", maptabdir, maptab);
    }
  if ( maptab )
    {
      if ( maptabbuild ) maptab = maptabbuild;
      int vlistID = streamInqVlist(streamID);

      if ( cdoVerbose )
        printf("*******Try to apply mapping table: '%s'*******\n", maptab);
      list_t *pml = cdo_parse_cmor_file(maptab);
      if ( pml == NULL )
        {
          cdoWarning("Mapping table: '%s' could not be parsed. Operator continues.", maptab);
          return;
        }
      const char *ventry[] = {"&parameter"};
      int nventry = (int) sizeof(ventry)/sizeof(ventry[0]);

      list_t *charvarlist = NULL; 
      if ( kvn )
        {
          if ( charvarlist = check_for_charvars(pml, "name") )
            {
              keyValues_t *charkvn = kvlist_search(charvarlist, "name");
              keyValues_t *charkvcn = kvlist_search(charvarlist, "cmor_name");
              if ( !charkvn || !charkvcn );
              else if ( charkvcn == kvcn )
                addcharvar(charkvn, vlistID, "name", vars);
            }
          if ( kvn->nvalues > 1 )
            cdoWarning("Only the first value of variable selection key 'name' is processed.");
          maptab_via_cmd(pml, kvn->values[0], vlistID, vlistNvars(vlistID),  "name", kvcn->values[0], miptabfreq);
          if ( cdoVerbose )
            printf("*******Successfully read mapping '%s' table.*******\n", maptab);
        }
      else if ( kvc )
        {
          if ( charvarlist = check_for_charvars(pml, "code") )
            {
              keyValues_t *charkvc = kvlist_search(charvarlist, "code");
              keyValues_t *charkvcn = kvlist_search(charvarlist, "cmor_name");
              if ( !charkvc || !charkvcn );
              else if ( charkvcn == kvcn )
                addcharvar(charkvc, vlistID, "code", vars);
            }
          if ( kvc->nvalues > 1 )
            cdoWarning("Only the first value of variable selection key 'code' is processed.");
          maptab_via_cmd(pml, kvc->values[0], vlistID, vlistNvars(vlistID), "code", kvcn->values[0], miptabfreq);
          if ( cdoVerbose )
            printf("*******Successfully read mapping '%s' table.*******\n", maptab);
        }
      else if ( kvcn )
        { 
          if ( charvarlist = check_for_charvars(pml, NULL) )
            {
              keyValues_t *charkvn = kvlist_search(charvarlist, "name");
              keyValues_t *charkvcn = kvlist_search(charvarlist, "cmor_name");
              if ( !charkvn || !charkvcn );
              else if ( strcmp(charkvcn->values[0], kvcn->values[0]) == 0 )
                addcharvar(charkvn, vlistID, "name", vars);
            }
          maptab_via_cn(pml, kvcn->values, vlistID, vlistNvars(vlistID), kvcn->nvalues, miptabfreq); 
          if ( cdoVerbose )
            printf("*******Successfully read mapping '%s' table.*******\n", maptab);
        }
      else
        {
          if ( charvarlist = check_for_charvars(pml, NULL) )
            {
              keyValues_t *charkvn = kvlist_search(charvarlist, "name");
              keyValues_t *charkvcn = kvlist_search(charvarlist, "cmor_name");
              if ( !charkvn || !charkvcn );
              else if ( charkvcn == kvcn )
                addcharvar(charkvn, vlistID, "name", vars);
            }
          for ( int varID = 0; varID < vlistNvars(vlistID); varID++ )
            {
              if ( maptab_via_key(pml, vlistID, varID, nventry, ventry, "name", miptabfreq) )
                {
                  printf("*******Successfully mapped varID '%d' via name.*******\n", varID);
                  continue;
                }
              if ( maptab_via_key(pml, vlistID, varID, nventry, ventry, "code", miptabfreq) )
                {
                  printf("*******Successfully mapped varID '%d' via code.*******\n", varID);
                  continue;
                }
              cdoWarning("Could not map variable with id '%d'.", varID);
            }
        }
      kv_insert_a_val(kvl, "mtproof", maptab, 1);
      list_destroy(pml);
      if ( maptabbuild ) Free(maptabbuild);
    }
  else if ( cdoVerbose )
    printf("*******No mapping table found.*******\n");
}

static char *check_short_key(char *key)
{
  char *short_keys[]={"cn", "n", "c", "u", "cm", "vc", "p", "szc", "i", "ca", "gi", "rtu", "mt", "om", "ms", "dr", "d", "lc", NULL};
  char *long_keys[]={"cmor_name", "name", "code", "units", "cell_methods", "variable_comment", "positive", "scalar_z_coordinate", "info", "character_axis", "grid_info", "required_time_units", "mapping_table", "output_mode", "max_size", "drs_root", "drs", "last_chunk", NULL};

  for ( int i = 0; short_keys[i]; i++ )
    if ( strcmp(key, short_keys[i]) == 0 || strcmp(key, long_keys[i]) == 0 )
      return short_keys[i];
/*  if ( strcmp(key, "cmor_name") == 0 ) short_key = strdup("cn");
  else if ( strcmp(key, "name") == 0 ) short_key = strdup("n");
  else if ( strcmp(key, "code") == 0 ) short_key = strdup("c");
  else if ( strcmp(key, "units") == 0 ) short_key = strdup("u");
  else if ( strcmp(key, "cell_methods") == 0 ) short_key = strdup("cm");
  else if ( strcmp(key, "comment") == 0 ) short_key = strdup("k");
  else if ( strcmp(key, "positive") == 0 ) short_key = strdup("p");
  else if ( strcmp(key, "scalar_z_coordinate") == 0 ) short_key = strdup("szc");
  else if ( strcmp(key, "info") == 0 ) short_key = strdup("i");
  else if ( strcmp(key, "character_axis") == 0 ) short_key = strdup("ca");
  else if ( strcmp(key, "grid_info") == 0 ) short_key = strdup("gi");
  else if ( strcmp(key, "required_time_units") == 0 ) short_key = strdup("rtu");
  else if ( strcmp(key, "mapping_table") == 0 ) short_key = strdup("mt");
  else if ( strcmp(key, "calendar") == 0 ) short_key = strdup("l");
  else if ( strcmp(key, "output_mode") == 0 ) short_key = strdup("om");
  else if ( strcmp(key, "max_size") == 0 ) short_key = strdup("ms");
  else if ( strcmp(key, "drs_root") == 0 ) short_key = strdup("dr");
  else if ( strcmp(key, "no_drs") == 0 ) short_key = strdup("d");
  else if ( strcmp(key, "last_chunk") == 0 ) short_key = strdup("lc"); */
  cdoWarning("Unknown commandline keyword: '%s'\n", key);
  return NULL;
}

static void parse_cmdline(list_t *pml, char **params, int nparams, char *ventry)
{
  list_t *kvl = NULL;
  kvl = list_new(sizeof(keyValues_t *), free_keyval, ventry);
  list_append(pml, &kvl);

  char *key = NULL, *eqpos = NULL;
  char **values = NULL;
  int i = 1, j = 0;
  while ( params[i] )
    {
      if ( eqpos = strchr(params[i], '=')  )
        {
          if ( key && values[0] )
            {
              char *short_key = check_short_key(key);
              if ( short_key )
                {
                  if ( strcmp(short_key, key) != 0 )
                    {
                      Free(key);
                      key = strdup(short_key);
                    }
                  kvlist_append(kvl, (const char *)key, (const char **) values, j);
                }
              Free(key);
              free_array(values);
            }
          else if ( key )
            cdoAbort("Found no value for key '%s'.", key);
          if ( strlen(eqpos) == 1 )
            cdoAbort("Could not find values for commandline parameter: '%s'\n", params[i]);
          key = strdup(strtok(params[i], "="));
          values = Malloc(100 * sizeof(char *));
          j = 0;   
          copy_value(strtok(NULL, ""), values, &j);  
        }
      else
        {
          if ( !key )
            cdoAbort("Found no key for value '%s'.", params[i]);
          else
            copy_value(params[i], values, &j);
        }
      i++;
    }
  if ( key && values )
    {
      char *short_key = check_short_key(key);
      if ( short_key )
        {
          if ( strcmp(short_key, key) != 0 )
            {
              Free(key);
              key = strdup(short_key);
            }
          kvlist_append(kvl, (const char *)key, (const char **) values, j);
        }
      Free(key);
      free_array(values);
    }
  else if ( values )
    cdoAbort("Found no key for value '%s'.", params[i-1]);
}

static char *get_mip_table(char *params, list_t *kvl, char *project_id)
{
  if ( !params )
    cdoAbort("A mip table name or path is required as first argument. No first argument found.");
  if ( file_exist(params, 0) )
    return params;
  else
    {
      cdoWarning("Your first argument is not an existing file. It is tried to build a path with additional configuration attributes 'mip_table_dir' and 'project_id'");
      char *miptabdir = kv_get_a_val(kvl, "mip_table_dir", NULL);
      if ( miptabdir && project_id )
        {
          char *miptab = Malloc((strlen(miptabdir)+strlen(project_id)+strlen(params)+3) * sizeof(char));
          sprintf(miptab, "%s/%s_%s\0", miptabdir, project_id, params);
          if ( file_exist(miptab, 0) )
            return miptab;
          else
            cdoAbort("Could not open mip table '%s'.", miptab);
        }
      else
        cdoAbort("Could not build a mip table path.");
    }        
}

static char *freq_from_path(char *mip_table)
{
  char *freq = mip_table;
  int fpos = 0, k = 0, j = 0;
  while ( *(mip_table + j) )
    {
      j++;
      if ( *(mip_table + j) == '/' )
        k = j + 1;
      if ( *(mip_table + j) == '_' && *(mip_table + j + 1) )
        fpos = j + 1;
    }
  freq += k;
  if ( fpos > k )
    freq += fpos-k;
  return freq;
}

static void save_miptab_freq(list_t *kvl, char *mip_table, int *miptab_freq)
{
  char *freq = freq_from_path(mip_table);
  if ( freq != NULL )
    {
      if ( strstr(freq, "yr") || strstr(freq, "Yr") )
        *miptab_freq = 11;
      else if ( strstr(freq, "mon") || strstr(freq, "Mon") )
        *miptab_freq = 12;
      else if ( strstr(freq, "day") || strstr(freq, "Day") )
        *miptab_freq = 13;
      else if ( strstr(freq, "6hr") )
        *miptab_freq = 14;
      else if ( strstr(freq, "3hr") )
        *miptab_freq = 15;

      if ( strcmp(freq, "Oclim") == 0 )
        *miptab_freq = 1;
      else if ( strcmp(freq, "Oyr") == 0 )
        *miptab_freq = 2;
      else if ( strcmp(freq, "cfMon") == 0 )
        *miptab_freq = 3;
      else if ( strcmp(freq, "day") == 0 )
        *miptab_freq = 4;
      else if ( strcmp(freq, "6hrPlev") == 0 )
        *miptab_freq = 5;
      else if ( strcmp(freq, "6hrLev") == 0 )
        *miptab_freq = 6;
    }
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

  if ( nparams < 1 ) cdoAbort("Too few arguments!");

  /* Define kvl and read cmdline */
  parse_cmdline(pml, params, nparams, "cmdline");
  
  const char *pmlistHelper[] = {"cmdline"};
  /* Get kvl and use it from now on instead of pml */
  list_t *kvl = pmlist_get_kvlist_ventry(pml, 1, pmlistHelper);
  char *name = kv_get_a_val(kvl, "n", NULL);
  char *code = kv_get_a_val(kvl, "c", NULL);
  char *cn = kv_get_a_val(kvl, "cn", NULL);
  if ( ( name && code ) )
    cdoAbort("Mapping via command line failed. Only one variable selector of 'name' and 'code' is allowed.");
  if ( ( name && !cn ) || ( code && !cn ) )
    cdoAbort("Mapping via command line failed. A corresponding 'cmor_name' is needed.");

  /* Config files are read with descending priority. */
  if ( cdoVerbose )
    printf("*******Start to read configuration files.*******\n");
  read_config_files(kvl);
  if ( cdoVerbose )
    printf("*******Successfully read configuration files.*******\n");

  /* check MIP table, MIP table frequency and project_id*/
  if ( cdoVerbose )
    printf("*******Start to check MIP table, MIP table frequency and project_id.*******\n");
  int miptab_freq = 0, time_axis = 0, calendar = 0;
  char *project_id, *dummy;
  if ( !(dummy = kv_get_a_val(kvl, "project_id", NULL) ) )
    cdoAbort("Value for attribute 'project_id' is required.");
  else
    project_id = strdup(dummy);
  char *mip_table = get_mip_table(params[0], kvl, project_id);
  save_miptab_freq(kvl, mip_table, &miptab_freq);
  char *miptab_freqptr = strdup(freq_from_path(mip_table));
  kv_insert_a_val(kvl, "miptab_freq", miptab_freqptr, 1);

  if ( cdoVerbose )
    printf("*******Successfully checked MIP table, MIP table frequency and project_id.*******\n");


  int streamID = streamOpenRead(cdoStreamName(0));
  /* Existing attributes have lowest priority. */
  dump_special_attributes(kvl, streamID);

  /* Check for attributes and member name */
  if ( cdoVerbose )
    printf("*******Start to check attributes.*******\n");
  check_attr(kvl, project_id);
  check_mem(kvl, project_id);
  if ( cdoVerbose )
    printf("*******Successfully checked global attributes.*******\n");

 /* dump_global_attributes(pml, streamID); */

  struct mapping *vars = construct_var_mapping(streamID);

 /* read mapping table */
  if ( cdoVerbose )
    printf("*******Start to read mapping table.*******\n");
  read_maptab(kvl, streamID, miptab_freqptr, vars);

  if ( cdoVerbose )
    printf("*******Start to use cmor_setup.*******\n");
  setup_dataset(kvl, streamID, &calendar);
  if ( cdoVerbose )
    printf("*******Succesfully used cmor_setup.*******\n");

  int table_id;
  cmor_load_table(mip_table, &table_id);
  cmor_set_table(table_id);

  register_all_dimensions(kvl, streamID, vars, table_id, project_id, miptab_freq, &time_axis);
  write_variables(kvl, &streamID, vars, miptab_freq, time_axis, calendar, miptab_freqptr);

  destruct_var_mapping(vars);
  Free(mip_table);
  Free(project_id); 
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

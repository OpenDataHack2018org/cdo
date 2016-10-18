/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2016 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

/*
   This module contains the following operators:

*/

#include <errno.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pmlist.h"

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
  while ( pos < len && !isspace((int) *(pline+pos)) && *(pline+pos) != '=' && *(pline+pos) != ':' ) pos++;

  strncpy(name, pline, pos);
  name[pos] = 0;

  pline += pos;
  return pline;
}

static
char *getElementValue(char *pline)
{
  while ( isspace((int) *pline) ) pline++;
  size_t len = strlen(pline);
  while ( isspace((int) *(pline+len-1)) && len ) { *(pline+len-1) = 0; len--;}

  return pline;
}

void pml_parse_buffer(list_t *pml, size_t buffersize, char *buffer)
{
  char line[4096];
  char name[256];
  char *pline;
  char listkey1[] = "axis_entry:";
  char listkey2[] = "variable_entry:";
  int linenumber = 0;
  int listtype = 0;
  list_t *kvl = NULL;

  while ( (buffer = readLineFromBuffer(buffer, &buffersize, line, sizeof(line))) )
    {
      linenumber++;
      pline = line;
      while ( isspace((int) *pline) ) pline++;
      if ( *pline == '#' || *pline == '!' || *pline == '\0' ) continue;
      //  len = (int) strlen(pline);
      if ( listtype == 0 && *pline == '&' )
	{
	  listtype = 1;
	}
      
      if ( strncmp(pline, listkey1, strlen(listkey1)) == 0 )
	{
	  pline += strlen(listkey1);

	  listtype = 2;

          kvl = list_new(sizeof(keyValues_t *), free_keyval, "axis");
          list_append(pml, &kvl);

	  pline = skipSeparator(pline);
	  pline = getElementValue(pline);

	  if ( *pline ) kvlist_append(kvl, "name", (const char **)&pline, 1);
	}
      else if ( strncmp(pline, listkey2, strlen(listkey2)) == 0 )
	{
	  pline += strlen(listkey2);

	  listtype = 2;

          kvl = list_new(sizeof(keyValues_t *), free_keyval, "variable");
          list_append(pml, &kvl);

	  pline = skipSeparator(pline);
	  pline = getElementValue(pline);

	  if ( *pline ) kvlist_append(kvl, "name", (const char **)&pline, 1);
	}
      else
	{
	  pline = getElementName(pline, name);
	  pline = skipSeparator(pline);
	  pline = getElementValue(pline);

	  if ( kvl == NULL )
            {
              kvl = list_new(sizeof(keyValues_t *), free_keyval, "global");
              list_append(pml, &kvl);
            }

	  if ( *pline ) kvlist_append(kvl, name, (const char **)&pline, 1);

	    {
	      //fprintf(stderr, "%d skip line %3d: %s\n", newlist, linenumber, pline);
	    }
	}

      //   printf("%s\n", pline);
    }
}

list_t *pml_parse_cmor_file(const char *filename)
{
  assert(filename != NULL);

  size_t filesize = fileSize(filename);

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

  /*
  if ( buffer[0] == '{' )
    kvlParseBufferJson(kvl);
  else
  */
  pml_parse_buffer(pml, filesize, buffer);
  
  return pml;
}

static
int read_cmor_table(const char *filename)
{
  list_t *pml = pml_parse_cmor_file(filename);
  if ( pml == NULL ) return -1;
  
  printf("# Number of lists: %d\n", list_size(pml));
  int i = 0;
  for ( listNode_t *pmnode = pml->head; pmnode; pmnode = pmnode->next )
    {
      list_t *kvl = *(list_t **)pmnode->data;
      printf("# list ID: %d;   Number of elements: %d\n", i, list_size(kvl));
      printf("&%s\n", list_name(kvl));
      for ( listNode_t *kvnode = kvl->head; kvnode; kvnode = kvnode->next )
        {
          keyValues_t *kv = *(keyValues_t **)kvnode->data;
          if ( kv ) printf("  %s = %s\n", kv->key, kv->values[0]);
        }
      printf("/\n");
      ++i;
    }

  list_destroy(pml);

  return 0;
}

static
int conv_cmor_table(const char *filename)
{
  list_t *pml = pml_parse_cmor_file(filename);
  if ( pml == NULL ) return -1;

  bool hasmissval = false;
  double missval;

  for ( listNode_t *pmnode = pml->head; pmnode; pmnode = pmnode->next )
    {
      list_t *kvl = *(list_t **)pmnode->data;
      const char *listname = list_name(kvl);

      if ( strncmp("global", listname, strlen(listname)) == 0 )
	{
          for ( listNode_t *kvnode = kvl->head; kvnode; kvnode = kvnode->next )
	    {
              keyValues_t *kv = *(keyValues_t **)kvnode->data;
              const char *ename  = kv->key;
	      const char *evalue = kv->values[0];
              size_t len = strlen(ename);

	      if ( strncmp("missing_value", ename, len) == 0 )
		{
		  missval = atof(evalue);
		  hasmissval = true;
		}
	    }
	}
      else if ( strncmp("variable", listname, strlen(listname)) == 0 )
	{
	  printf("&%s\n", "parameter");
          for ( listNode_t *kvnode = kvl->head; kvnode; kvnode = kvnode->next )
	    {
              keyValues_t *kv = *(keyValues_t **)kvnode->data;
              const char *ename  = kv->key;
	      const char *evalue = kv->values[0];
	      int len = strlen(ename);
	      int vlen = strlen(evalue);

	      if ( vlen > 1 && evalue[0] == '"' && evalue[vlen-1] == '"' ) 
		{
		  vlen -= 2;
		  evalue++;
		}

	      char *ovalue = strdup(evalue);
	      for ( int i = 1; i < vlen; ++i )
		{
		  if ( ovalue[i-1] == '"' && ovalue[i] == '"' )
		    {
		      ovalue [i-1] = '\'';
		      for ( int j = i+1; j < vlen; ++j ) ovalue[j-1] = ovalue[j];
		      vlen -= 1;
		    }
		}

	      if ( strncmp("name", ename, len)            == 0 ||
		   strncmp("standard_name", ename, len)   == 0 ||
		   strncmp("out_name", ename, len)        == 0 ||
		   strncmp("type", ename, len)            == 0 ||
		   strncmp("valid_min", ename, len)       == 0 ||
		   strncmp("valid_max", ename, len)       == 0 ||
		   strncmp("ok_min_mean_abs", ename, len) == 0 ||
		   strncmp("ok_max_mean_abs", ename, len) == 0 )
		printf("  %-15s = %s\n", ename, ovalue);
	      else if ( strncmp("long_name", ename, len)  == 0 ||
		   strncmp("units", ename, len)           == 0 ||
		   strncmp("cell_methods", ename, len)    == 0 ||
		   strncmp("cell_measures", ename, len)   == 0 ||
		   strncmp("comment", ename, len)         == 0 )
		printf("  %-15s = \"%.*s\"\n", ename, vlen, ovalue);

	      Free(ovalue);
	    }
	  if ( hasmissval ) printf("  %-15s = %g\n", "missing_value", missval);
	  printf("/\n");
	}
    }

  list_destroy(pml);

  return 0;
}


void *Kvl(void *argument)
{
  cdoInitialize(argument);

  int READ_CMOR_TABLE = cdoOperatorAdd("read_cmor_table",   0,   0, NULL);
  int CONV_CMOR_TABLE = cdoOperatorAdd("conv_cmor_table",   0,   0, NULL);

  int operatorID = cdoOperatorID();

  if ( operatorArgc() != 1 ) cdoAbort("Too few arguments!");
  const char *filename = operatorArgv()[0];

  if ( cdoVerbose ) cdoPrint("Parse file: %s", filename);

  // if      ( operatorID == READ_CMOR_TABLE ) read_cmor_table(filename);
  if      ( operatorID == READ_CMOR_TABLE ) read_cmor_table(filename);
  else if ( operatorID == CONV_CMOR_TABLE ) conv_cmor_table(filename);

  cdoFinish();

  return 0;
}

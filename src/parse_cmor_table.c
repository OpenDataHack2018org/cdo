#include <errno.h>
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

list_t *cdo_parse_cmor_file(const char *filename)
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

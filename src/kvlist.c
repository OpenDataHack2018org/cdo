/* key value list */

/*
  list type 1:
  ============

  &parameter
    name:              ps
    standard_name:     surface_air_pressure
    units:             Pa
    cell_methods:      "time: mean"
    cell_measures:     "area: areacella"
    long_name:         "Surface Air Pressure"
    comment:           "not, in general, the same as mean sea-level pressure"
    valid_min:         4.791e+04
    valid_max:         1.119e+05

  list type 2: one entry on each line
  ============

  variable_entry:    ps
    standard_name:     surface_air_pressure
    units:             Pa
    cell_methods:      time: mean
    cell_measures:     area: areacella
    long_name:         Surface Air Pressure
    comment:           not, in general, the same as mean sea-level pressure
    valid_min:         4.791e+04
    valid_max:         1.119e+05
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <assert.h>


#define MAX_KVLISTS    4096
#define MAX_KVELEMENTS 1024

typedef struct {
  char name[128];
  char *value;
} kvelement_t;

typedef struct {
  char name[128];
  int num_elements;
  kvelement_t elements[MAX_KVELEMENTS];
} kvlist_t;

typedef struct {
  char *filename;
  char *buffer;
  char *bufferp;
  size_t buffersize;
  int num_lists;
  kvlist_t lists[MAX_KVLISTS];
} kvl_t;


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
      if ( ichar == '\r' ) break;
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

  return (buffer);
}

static
void pfree(void *ptr)
{
  if ( ptr ) free(ptr);
}

static
char *skipSeparator(char *pline)
{
  while ( isspace((int) *pline) ) pline++;
  if ( *pline == '=' || *pline == ':' ) pline++;
  while ( isspace((int) *pline) ) pline++;

  return (pline);
}

static
void kvlParseBuffer(kvl_t *kvl)
{
  char line[4096];
  char *listname;
  char *pline;
  char *buffer = kvl->buffer;
  size_t buffersize = kvl->buffersize;
  int len;
  int liststart = 0;
  int linenumber = 0;
  int newlist = 0;
  char listkey2[] = "variable_entry:";
  int listtype = 0;

  while ( (buffer = readLineFromBuffer(buffer, &buffersize, line, sizeof(line))) )
    {
      linenumber++;
      listname = NULL;
      pline = line;
      while ( isspace((int) *pline) ) pline++;
      if ( *pline == '#' || *pline == '!' || *pline == '\0' ) continue;
      //  len = (int) strlen(pline);
      if ( *pline == '&' )
	{
	  // if ( liststart == 1 ) endlist;
	  //	  listname = *pline++;
	}
      
      if ( strncmp(pline, listkey2, strlen(listkey2)) == 0 )
	{
	  pline += strlen(listkey2);

	  if ( newlist == 1 )
	    {
	      //  newlist = 0;
	    }
	  newlist = 1;
	  listtype = 2;

	  kvl->num_lists++;
	  {
	    int nlist = kvl->num_lists-1;
	    strcpy(kvl->lists[nlist].name, "parameter");
	    pline = skipSeparator(pline);
	    while ( isspace((int) *pline) ) pline++;
	    len = strlen(pline);
	    while ( isspace((int) *(pline+len-1)) ) { *(pline+len-1) = 0; len--;}
	    printf(">>>%s< %d\n", pline, len);
	  }
	}
      else
	{
	  if ( newlist == 0 )
	    {
	      //fprintf(stderr, "%d skip line %3d: %s\n", newlist, linenumber, pline);
	    }
	}

      liststart = 1;
      //   printf("%s\n", pline);
    }
}


kvl_t *kvlParseFile(const char *filename)
{
  kvl_t *kvl = NULL;
  FILE *fp;
  char *buffer;
  size_t filesize;
  size_t nitems;

  assert(filename != NULL);

  fp = fopen(filename, "r");
  if ( fp == NULL )
    {
      fprintf(stderr, "Open failed on %s: %s\n", filename, strerror(errno));
      return (kvl);
    }

  /* file size */
  fseek(fp, 0L, SEEK_END);
  filesize = (size_t) ftell(fp);
  fseek(fp, 0L, SEEK_SET);

  buffer = (char *) malloc(filesize);
  nitems = fread(buffer, 1, filesize, fp);

  fclose(fp);

  if ( nitems != filesize )
    {
      fprintf(stderr, "Read failed on %s!\n", filename);
      return (kvl);
    }
 
  kvl = (kvl_t *) calloc(1, sizeof(kvl_t));
  kvl->buffer = buffer;
  kvl->buffersize = filesize;
  kvl->filename = strdup(filename);

  kvlParseBuffer(kvl);
  
  return (kvl);
}


void kvlDelete(kvl_t *kvl)
{
  assert(kvl != NULL);

  pfree(kvl->filename);
  pfree(kvl->buffer);

  free(kvl);
}


int kvlGetNumLists(kvl_t *kvl)
{
  assert(kvl != NULL);

  return(kvl->num_lists);
}


const char *kvlGetListName(kvl_t *kvl, int listID)
{
  char *listname = NULL;
  
  assert(listID < kvl->num_lists);
  
  listname = kvl->lists[listID].name;

  return (listname);
}

int kvlGetListNumElements(kvl_t *kvl, int listID)
{
  int nelements = 0;

  assert(listID < kvl->num_lists);

  nelements = kvl->lists[listID].num_elements;

  return (nelements);
}


const char *kvlGetListElementName(kvl_t *kvl, int listID, int elemID)
{
  char *ename = NULL;

  assert(listID < kvl->num_lists);
  assert(elemID < kvl->lists[listID].num_elements);

  ename = kvl->lists[listID].elements[elemID].name;

  return (ename);
}


const char *kvlGetListElementValue(kvl_t *kvl, int listID, int elemID)
{
  char *evalue = NULL;

  assert(listID < kvl->num_lists);
  assert(elemID < kvl->lists[listID].num_elements);

  evalue = kvl->lists[listID].elements[elemID].value;

  return (evalue);
}


int main(int argc, char *argv[])
{
  char *filename;
  kvl_t *kvlist;
  int nlists, listID;
  int nelements, elemID;
  const char *listname;
  const char *ename;
  const char *evalue;

  if ( argc != 2 ) 
    {
      fprintf(stderr, "usage: kvlist filename\n");
      return (1);
    }

  filename = argv[1];

  printf("Parse kvlist file: %s\n", filename);

  kvlist = kvlParseFile(filename);
  nlists = kvlGetNumLists(kvlist);
  printf("# Number of lists: %d\n", nlists);
  for ( listID = 0; listID < nlists; ++listID )
    {
      listname = kvlGetListName(kvlist, listID);
      nelements = kvlGetListNumElements(kvlist, listID);
      printf("# list ID: %d;   Number of elements: %d\n", listID, nelements);
      printf("&%s\n", listname);
      for ( elemID = 0; elemID < nelements; ++elemID )
	{
	  ename  = kvlGetListElementName(kvlist, listID, elemID);
	  evalue = kvlGetListElementValue(kvlist, listID, elemID);
	  printf("  %s = %s\n", ename, evalue);
	}
      printf("/\n");
    }

  kvlDelete(kvlist);

  return (0);
}

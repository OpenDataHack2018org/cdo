#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "list.h"


bool iterate_int(void *data) 
{
  printf("  Found value: %d\n", *(int *)data);
  return true;
}

bool iterate_string(void *data)
{
  printf("  Found string value: %s\n", *(char **)data);
  return true;
}

bool iterate_list(void *data)
{
  printf("Found %s list with %d strings: \n", list_name(*(list_t **)data), list_size(*(list_t **)data));
  list_for_each(*(list_t **)data, iterate_string);
  return true;
}

void free_string(void *data)
{
  free(*(char **)data);
}

void free_list(void *data)
{
  int n = list_size(*(list_t **)data);
  list_destroy(*(list_t **)data);
  printf("Successfully freed %d strings...\n", n);
}

void list_with_ints()
{
  int numbers = 10;
  printf("Generating list with the first %d positive numbers...\n", numbers);
 
  list_t *list = list_new(sizeof(int), NULL, NULL);
 
  for ( int i = 1; i <= numbers; i++ )
    list_append(list, &i);
 
  list_for_each(list, iterate_int);
 
  list_destroy(list);
  printf("Successfully freed %d numbers...\n", numbers);
}

void list_with_strings()
{
  int numNames = 5;
  const char *names[] = { "Kalle", "Luis", "Ralf", "Rene", "Uwe" };
  printf("Generating list with %d different names...\n", numNames);
 
  list_t *list = list_new(sizeof(char *), free_string, "strings");
 
  for ( int i = 0; i < numNames; i++ )
    {
      char *name = strdup(names[i]);
      list_append(list, &name);
    }
 
  list_for_each(list, iterate_string);
 
  list_destroy(list);
  printf("Successfully freed %d strings...\n", numNames);
}

void list_with_listsofstrings()
{
  int numLists = 0;
  printf("Generating list with lists of names...\n");

  list_t *list = list_new(sizeof(list_t *), free_list, "parameter");

  {
    list_t *strlist = list_new(sizeof(char *), free_string, "p1");
    int numNames = 5;
    const char *names[] = { "Kalle", "Luis", "Ralf", "Rene", "Uwe" };
    for ( int i = 0; i < numNames; i++ )
      {
        char *name = strdup(names[i]);
        list_append(strlist, &name);
      }
    list_append(list, &strlist);
    numLists++;
  }
  {
    list_t *strlist = list_new(sizeof(char *), free_string, "p2");
    int numNames = 3;
    const char *names[] = { "Monika", "Ulrike", "Katja" };
    for ( int i = 0; i < numNames; i++ )
      {
        char *name = strdup(names[i]);
        list_append(strlist, &name);
      }
    list_append(list, &strlist);
    numLists++;
  }
 
  list_for_each(list, iterate_list);
 
  list_destroy(list);
  printf("Successfully freed %d lists...\n", numLists);
}

int main(int argc, char *argv[])
{
  printf("Loading list demo...\n");
  list_with_ints();
  list_with_strings();
  list_with_listsofstrings();
}

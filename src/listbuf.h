#ifndef _LISTBUF_H
#define _LISTBUF_H

#include <stdio.h>


typedef struct
{
  size_t size;
  char *buffer;
  char *name;
} listbuf_t;


listbuf_t *listbuf_new();
void listbuf_destroy(listbuf_t *listbuf);
int listbuf_read(listbuf_t *listbuf, FILE *fp, const char *name);

#endif

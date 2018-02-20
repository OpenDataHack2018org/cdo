/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2018 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/
#ifndef _LIST_H
#define _LIST_H

#include <stdio.h>
#include <stdbool.h>

// a common function used to free malloc'd objects
typedef void (*freeFunction_t)(void *);

typedef bool (*listIterator_t)(void *);

typedef struct _listNode_t
{
  void *data;
  struct _listNode_t *next;
} listNode_t;

typedef struct
{
  size_t logicalLength;
  size_t elementSize;
  char *name;
  listNode_t *head;
  listNode_t *tail;
  freeFunction_t freeFunc;
} list_t;

list_t *list_new(size_t elementSize, freeFunction_t freeFunc, const char *name);
void list_destroy(list_t *list);

void list_prepend(list_t *list, void *element);
void list_append(list_t *list, void *element);
int list_size(list_t *list);
const char *list_name(list_t *list);

void list_for_each(list_t *list, listIterator_t iterator);
void list_head(list_t *list, void *element, bool removeFromList);
void list_tail(list_t *list, void *element);
void *list_entry(list_t *list, int index);

#endif

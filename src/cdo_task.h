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
#ifndef CDO_TASK_H
#define CDO_TASK_H

extern bool CDO_task;

void cdo_task_start(void *task, void *(*task_routine)(void *), void *task_arg);
void *cdo_task_wait(void *task);
void *cdo_task_new();
void cdo_task_delete(void *task);

#endif /* CDO_TASK_H */

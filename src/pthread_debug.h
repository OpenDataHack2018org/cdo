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
#ifndef _PTHREAD_DEBUG_H
#define _PTHREAD_DEBUG_H

#include <mutex>
#include <condition_variable>

void Pthread_debug(int debug);

int Pthread_create(const char *caller, pthread_t *th, pthread_attr_t *attr, void *(*start_routine)(void *), void *arg);

int Pthread_join(const char *caller, pthread_t th, void **thread_return);

void Pthread_mutex_lock(const char *caller, pthread_mutex_t *mutex);
void Pthread_mutex_lock(const char *caller, std::mutex &p_mutex);

void Pthread_mutex_unlock(const char *caller, pthread_mutex_t *mutex);
void Pthread_mutex_unlock(const char *caller, std::mutex &p_mutex);

void Pthread_cond_signal(const char *caller, pthread_cond_t *cond);
void Pthread_cond_signal(const char *caller, std::condition_variable &p_cond_var);
void Pthread_cond_wait(const char *caller, pthread_cond_t *cond, pthread_mutex_t *mutex);
void Pthread_cond_wait(const char *caller, std::condition_variable &p_cond_var, std::unique_lock<std::mutex> &p_mutex);

void print_pthread_attr(const char *caller, pthread_attr_t *attr);
void print_pthread_mutexattr(const char *caller, pthread_mutexattr_t *m_attr);
void print_pthread_condattr(const char *caller, pthread_condattr_t *c_attr);

#define pthread_create(a, b, c, d) Pthread_create(__func__, a, b, c, d)
#define pthread_join(a, b) Pthread_join(__func__, a, b)

#define pthread_mutex_lock(a) Pthread_mutex_lock(__func__, a)
#define pthread_mutex_unlock(a) Pthread_mutex_unlock(__func__, a)

#define pthread_cond_signal(a) Pthread_cond_signal(__func__, a)
#define pthread_cond_wait(a, b) Pthread_cond_wait(__func__, a, b)

#endif /* _PTHREAD_DEBUG_H */

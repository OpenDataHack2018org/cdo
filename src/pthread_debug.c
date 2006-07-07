#if  defined  (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <stdio.h>



#if  defined  (HAVE_LIBPTHREAD)
#include <pthread.h>
#include <errno.h>
#include "error.h"


#define POUT2(caller, x, a, b)     pout2(caller, #x, x, #a, a, #b, b)
#define POUT3(caller, x, a, b, c)  pout3(caller, #x, x, #a, a, #b, b, #c, c)


static void pout2(const char *caller,
		  const char *sval, int ival, 
		  const char *sval1, int oval1,
		  const char *sval2, int oval2)
{
  if ( ival == oval1 )
    fprintf(stderr, "%-18s :  %-14s = %s\n", caller, sval, sval1);
  else if ( ival == oval2 )
    fprintf(stderr, "%-18s :  %-14s = %s\n", caller, sval, sval2);
  else
    fprintf(stderr, "%-18s :  %-14s = %d\n", caller, sval, ival);
}


static void pout3(const char *caller,
		  const char *sval, int ival, 
		  const char *sval1, int oval1,
		  const char *sval2, int oval2,
		  const char *sval3, int oval3)
{
  if ( ival == oval1 )
    fprintf(stderr, "%-18s :  %-14s = %s\n", caller, sval, sval1);
  else if ( ival == oval2 )
    fprintf(stderr, "%-18s :  %-14s = %s\n", caller, sval, sval2);
  else if ( ival == oval3 )
    fprintf(stderr, "%-18s :  %-14s = %s\n", caller, sval, sval3);
  else
    fprintf(stderr, "%-18s :  %-14s = %d\n", caller, sval, ival);
}


void print_pthread_attr(const char *caller, pthread_attr_t *attr)
{
  struct sched_param param;
  int detachstate, policy, inherit, scope, priority;
  size_t stacksize;

  pthread_attr_getdetachstate(attr, &detachstate);
  POUT2(caller, detachstate, PTHREAD_CREATE_JOINABLE, PTHREAD_CREATE_DETACHED);

#if defined (SCHED_FIFO)
  pthread_attr_getschedpolicy(attr, &policy);
  POUT3(caller, policy, SCHED_FIFO, SCHED_RR, SCHED_OTHER);
  pthread_attr_getschedparam(attr, &param);
  priority = param.sched_priority;
  fprintf(stderr, "%-18s :  %-14s = %d\n", caller, "priority", priority);
#endif

#if defined (PTHREAD_INHERIT_SCHED)
  pthread_attr_getinheritsched(attr, &inherit);
  POUT2(caller, inherit, PTHREAD_INHERIT_SCHED, PTHREAD_EXPLICIT_SCHED);
#endif

  pthread_attr_getscope(attr, &scope);
  POUT2(caller, scope, PTHREAD_SCOPE_SYSTEM, PTHREAD_SCOPE_PROCESS);

  pthread_attr_getstacksize(attr, &stacksize);
  fprintf(stderr, "%-18s :  %-14s = %ld\n", caller, "stacksize", (long) stacksize);
}


void print_pthread_mutexattr(const char *caller,  pthread_mutexattr_t *m_attr)
{
  int protocol, kind, pshared;
  /*
    Does not work on cygwin! PTHREAD_PRIO_INHERIT is defined without contents.
#if defined (PTHREAD_PRIO_INHERIT)
  pthread_mutexattr_getprotocol(m_attr, &protocol);
  POUT3(caller, protocol, PTHREAD_PRIO_INHERIT, PTHREAD_PRIO_PROTECT, PTHREAD_PRIO_NONE);
#endif
  */
#if defined (PTHREAD_MUTEX_FAST_NP)
  pthread_mutexattr_getkind_np(m_attr, &kind);
  POUT3(caller, kind, PTHREAD_MUTEX_FAST_NP, PTHREAD_MUTEX_RECURSIVE_NP, PTHREAD_MUTEX_ERRORCHECK_NP);
#endif

  pthread_mutexattr_getpshared(m_attr, &pshared);
  POUT2(caller, pshared, PTHREAD_PROCESS_SHARED, PTHREAD_PROCESS_PRIVATE);
}


void print_pthread_condattr(const char *caller, pthread_condattr_t *c_attr)
{
  int pshared;

  pthread_condattr_getpshared(c_attr, &pshared);
  POUT2(caller, pshared, PTHREAD_PROCESS_SHARED, PTHREAD_PROCESS_PRIVATE);
}


int PTHREAD_Debug = 0;


void Pthread_debug(int debug)
{
  PTHREAD_Debug = debug;
}


int Pthread_create(const char *caller, pthread_t *th,
		   pthread_attr_t *attr, void * (*start_routine)(void *), void *arg)
{
  static char func[] = "Pthread_create";
  int status;

  if ( PTHREAD_Debug ) Message(caller, "+%s", func);

  if ( PTHREAD_Debug )
    {
      Message(caller, "+%s attributes:", func);
      if ( attr )
	print_pthread_attr(func, attr);
      else
	Message(func, "  default attributes");
    }

  status = pthread_create(th, attr, start_routine, arg);

  if ( PTHREAD_Debug ) Message(caller, "-%s (thID = %ld, status = %d)",
			       func, (long) *th, status);

  return (status);
}


int Pthread_join(const char *caller, pthread_t th, void **thread_return)
{
  static char func[] = "Pthread_join";
  int status;

  if ( PTHREAD_Debug ) Message(caller, "+%s (thID = %ld)", func, (void *) th);

  status = pthread_join(th, thread_return);

  if ( PTHREAD_Debug ) Message(caller, "-%s (thID = %ld, status = %d)", func, (void *) th,
			       status);

  return (status);
}


void Pthread_mutex_lock(const char *caller, pthread_mutex_t *mutex)
{
  static char func[] = "Pthread_mutex_lock";
  int status;

  if ( PTHREAD_Debug ) Message(caller, "+%s (mutex = %p)", func, (void *) mutex);

  status = pthread_mutex_lock(mutex);
  if ( status != 0 )
    {
      switch (status)
	{
	case EINVAL:
	  Error(func, "the mutex has not been properly initialized");
	  break;
	case EDEADLK:
	  Error(func, "the mutex is already locked by the calling thread");
	  break;
	default:
	  Error(func, "status %d unknown", status, (void *) mutex);
	}
    }

  if ( PTHREAD_Debug ) Message(caller, "-%s (mutex = %p)", func, (void *) mutex);
}


void Pthread_mutex_unlock(const char *caller, pthread_mutex_t *mutex)
{
  static char func[] = "Pthread_mutex_unlock";
  int status;

  if ( PTHREAD_Debug ) Message(caller, "+%s (mutex = %p)", func, (void *) mutex);

  status = pthread_mutex_unlock(mutex);
  if ( status != 0 )
    {
      switch (status)
	{
	case EINVAL:
	  Error(func, "the mutex has not been properly initialized");
	  break;
	case EPERM:
	  Error(func, "the calling thread does not own the mutex");
	  break;
	default:
	  Error(func, "status %d unknown", status);
	}
    }

  if ( PTHREAD_Debug ) Message(caller, "-%s (mutex = %p)", func, (void *) mutex);
}


void Pthread_cond_signal(const char *caller, pthread_cond_t *cond)
{
  static char func[] = "Pthread_cond_signal";

  if ( PTHREAD_Debug ) Message(caller, "+%s (cond = %p)", func, (void *) cond);

  pthread_cond_signal(cond);

  if ( PTHREAD_Debug ) Message(caller, "-%s (cond = %p)", func, (void *) cond);
}


void Pthread_cond_wait(const char *caller, pthread_cond_t *cond, pthread_mutex_t *mutex)
{
  static char func[] = "Pthread_cond_wait";

  if ( PTHREAD_Debug ) Message(caller, "+%s (cond = %p, mutex =  %p)",
			       func, (void *) cond, (void *) mutex);

  pthread_cond_wait(cond, mutex);

  if ( PTHREAD_Debug ) Message(caller, "-%s (cond = %p, mutex = %p)",
			       func, (void *) cond, (void *) mutex);
}

#endif

#ifndef _CDO_TASK_H
#define _CDO_TASK_H

void cdo_task_start(void *task, void *(*task_routine)(void *), void *task_arg);
void *cdo_task_wait(void *task);
void *cdo_task_new();
void cdo_task_delete(void *task);

#endif  /* _CDO_TASK_H */

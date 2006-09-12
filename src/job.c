#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>

#include "cdo.h"

#if  defined  (HAVE_LIBDRMAA)
#  include "drmaa.h"
#endif


#if  defined  (HAVE_LIBDRMAA)
static drmaa_job_template_t *create_job_template(const char *expname, const char *jobfilename, const char *jobname)
{
  static char func[] = "create_job_template";
  drmaa_job_template_t *job = NULL;

  char error[DRMAA_ERROR_STRING_BUFFER];
  char name[DRMAA_ATTR_BUFFER], value[DRMAA_ATTR_BUFFER];
  char attr[1024];

  long size;
  char *dir, *ptr;

  drmaa_attr_names_t *job_attributes;

  int drmaa_errno;

  char host[1024];

  char *output_path;

  int len, len1, len2;

  /* determine hostname */

  gethostname(host, 255);  

  /* determine current path */

  size = pathconf(".", _PC_PATH_MAX);
  if ((dir = (char *)malloc((size_t)size)) != NULL) {
    ptr = getcwd(dir, (size_t)size);
  }

  /* generate DRMAA conform output path */
  
  len1 = strlen(host);
  len2 = strlen(dir);
  len = len1+len2+1;

  output_path = (char *) malloc(len*sizeof(char));
  strcpy(output_path, host);
  strcat(output_path, ":");
  strcat(output_path, dir);

  /* need to allow chdir on execution host, not thread save! */

  setenv("SGE_DRMAA_ALLOW_CWD", "yes", 1);

  /* allocate job template */

  if (drmaa_allocate_job_template(&job, NULL, 0) != DRMAA_ERRNO_SUCCESS)
    return NULL;

  /* the job's name */
  drmaa_set_attribute(job, DRMAA_JOB_NAME, jobname, NULL, 0);

  /* the job to be run */
  drmaa_set_attribute(job, DRMAA_REMOTE_COMMAND, jobfilename, NULL, 0);

  /* submit state */
  drmaa_set_attribute(job, DRMAA_JS_STATE, "drmaa_active", NULL, 0);

  /* working directory on execution host */
  drmaa_set_attribute(job, DRMAA_WD, dir, NULL, 0);

  /* path for output */
  drmaa_set_attribute(job, DRMAA_OUTPUT_PATH, output_path, NULL, 0);

  /* join output/error file */
  drmaa_set_attribute(job, DRMAA_JOIN_FILES, "y", NULL, 0);

  /* transfer files */
  drmaa_set_attribute(job, DRMAA_TRANSFER_FILES, "ieo", NULL, 0);
  
  /* some native SGE commands necessary */
  sprintf(attr, "-cwd -b n -q %s.q", expname);
  drmaa_set_attribute(job, DRMAA_NATIVE_SPECIFICATION, attr, NULL, 0);  

  /* print out job attributes */
  drmaa_get_attribute_names (&job_attributes, error, DRMAA_ERROR_STRING_BUFFER);

  if ( cdoVerbose )
    while ((drmaa_errno = drmaa_get_next_attr_name(job_attributes, name, DRMAA_ATTR_BUFFER)) == DRMAA_ERRNO_SUCCESS) {
      drmaa_get_attribute (job, name, value, DRMAA_ATTR_BUFFER, error,  DRMAA_ERROR_STRING_BUFFER);

      fprintf (stderr, "name: %-25s \t %s\n", name, value);
    }

  free(dir);

  return job;
}
#endif


#if  defined  (HAVE_LIBDRMAA)
static int drmaa_submit(const char *expname, const char *jobfilename, const char *jobname)
{
  char status[DRMAA_ERROR_STRING_BUFFER];

  char jobid[DRMAA_JOBNAME_BUFFER], jobout[DRMAA_JOBNAME_BUFFER];

  int drmaa_errno, stat;

  drmaa_job_template_t *job;

  int aborted, exited, signaled, exit_status;

  drmaa_attr_values_t *rusage = NULL;

  char usage[DRMAA_ERROR_STRING_BUFFER];

  if (drmaa_init(NULL, status, sizeof(status)-1) != DRMAA_ERRNO_SUCCESS) {
    fprintf(stderr, "drmaa_init() failed: %s\n", status);
    return 1;
  }

  /* submit some sequential jobs */

  if (!(job = create_job_template(expname, jobfilename, jobname))) {
    fprintf(stderr, "create_job_template() failed\n");
    return 1;
  }

  while ((drmaa_errno = drmaa_run_job(jobid, sizeof(jobid)-1, job, status,
	  sizeof(status)-1)) == DRMAA_ERRNO_DRM_COMMUNICATION_FAILURE) {
    fprintf(stderr, "drmaa_run_job() failed - retry: %s\n", status);
    sleep(1);
  }
  if (drmaa_errno != DRMAA_ERRNO_SUCCESS) {
    fprintf(stderr, "drmaa_run_job() failed: %s\n", status);
    return 1;
  }
  fprintf(stderr, "job %s submitted.\n", jobid);

  drmaa_delete_job_template(job, NULL, 0);

  /* wait for job */

  drmaa_errno = drmaa_wait(jobid, jobout, sizeof(jobout)-1, 
			   &stat, DRMAA_TIMEOUT_WAIT_FOREVER, &rusage, status, sizeof(status)-1);

  if (drmaa_errno != DRMAA_ERRNO_SUCCESS) {
    fprintf(stderr, "drmaa_wait(%s) failed: %s\n", jobout, status);
    return 1;
  }
  
  /*
   * report how job finished 
   */
  drmaa_wifaborted(&aborted, stat, NULL, 0);
  if (aborted) {
    printf("job %s never ran\n", jobid);
  } else {
    drmaa_wifexited(&exited, stat, NULL, 0);
    if (exited) {
      drmaa_wexitstatus(&exit_status, stat, NULL, 0);
      printf("job %s finished regularly with exit status %d\n", 
	     jobid, exit_status);
    } else {
      drmaa_wifsignaled(&signaled, stat, NULL, 0);
      if (signaled) {
	char termsig[DRMAA_SIGNAL_BUFFER+1];
	drmaa_wtermsig(termsig, DRMAA_SIGNAL_BUFFER, stat, NULL, 0);
	printf("job %s finished due to signal %s\n", jobid, termsig);

      } else
	printf("job %s finished with unclear conditions\n", jobid);
    }
  }

  if ( cdoVerbose )
    {
      fprintf(stderr, "Job usage:\n");
                
      while (drmaa_get_next_attr_value (rusage, usage, DRMAA_ERROR_STRING_BUFFER) == DRMAA_ERRNO_SUCCESS) {
	fprintf(stderr, "  %s\n", usage);
      }
    }
                
  drmaa_release_attr_values (rusage);


  if (drmaa_exit(status, sizeof(status)-1) != DRMAA_ERRNO_SUCCESS) {
    fprintf(stderr, "drmaa_exit() failed: %s\n", status);
    return 1;
  }

  {
    char commandline[1024];

    sprintf(commandline, "cat %s.o%s | grep -v tty | grep -v cannot | grep -v resize | grep -v shell\n", jobname, jobid);
    system(commandline);

    sprintf(commandline, "rm -f %s.o%s\n", jobname, jobid);
    system(commandline);
  }

  return 0;
}
#endif


void job_submit(const char *expname, const char *jobfilename, const char *jobname)
{
#if  defined  (HAVE_LIBDRMAA)
  int status;

  status = drmaa_submit(expname, jobfilename, jobname);
#else
  fprintf(stderr, "DRMAA library not available!\n");
#endif
}


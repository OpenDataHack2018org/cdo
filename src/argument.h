#ifndef ARGUMENT_H
#define ARGUMENT_H

#include <string>

struct argument_t
{
  //argument_t();
  int process_id;
  int    argc;
  char **argv;
  char  *args;
  std::string operatorName;
  char * operatorArguments;
};

argument_t makeArgument(int argc, char *argv[]);
argument_t *file_argument_new(const char *filename);
void        file_argument_free(argument_t *argument);
argument_t *argument_new(size_t argc, size_t len);
void        argument_free(argument_t *argument);
void        argument_fill(argument_t *argument, int argc, char *argv[]);
void print_argument(argument_t *argument);
#endif
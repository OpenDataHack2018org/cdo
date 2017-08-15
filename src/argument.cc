
#include "argument.h"
#include "dmemory.h"

#include <string>
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

argument_t *file_argument_new(const char *filename)
{
  argument_t *argument = new argument_t();

  argument->argc = 1;
  argument->argv.resize(argument->argc);
  argument->argv[0] = (char *) filename;
  argument->args = (char *) filename;

  return argument;
}

void file_argument_free(argument_t *argument)
{
  if ( argument )
    {
      if ( argument->argc )
        {
          assert(argument->argc == 1);
        }
      delete(argument);
    }
}

argument_t *argument_new(size_t argc, size_t len)
{
  argument_t *argument = new argument_t();

  if ( argc > 0 )
    {
      argument->argc = argc;
      argument->argv.resize(argc);
    }

  if ( len > 0 )
    argument->args = (char*) Calloc(len, sizeof(char));

  return argument;
}

void argument_free(argument_t *argument)
{
  if ( argument )
    {
      if ( argument->argc )
        {
          int argc =  argument->argc;
          for ( int i = 0; i < argc; ++i )
            {
              if ( argument->argv[i] )
                {
                    argument->argv[i] = NULL;
                }
            }

          argument->argc = 0;
        }

      if ( argument->args )
        {
          Free(argument->args);
          argument->args = NULL;
        }

      Free(argument);
    }
}

void argument_fill(argument_t *argument, int argc, char *argv[])
{
  assert(argument->argc == argc);

  for ( int iarg = 0; iarg < argc; ++iarg )
    argument->argv[iarg] = strdup(argv[iarg]);
}

void print_argument(argument_t * p_argument)
{
    std::cout << "argv with " << p_argument->argc << " arguments:" << std::endl;
    for(int i = 0; i < p_argument->argc; i++)
    {
        std::cout << p_argument->argv[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "OperatorName: "<< p_argument->operatorName << std::endl;

    std::cout << "operatorArguments: " << p_argument->operatorArguments << std::endl;
}

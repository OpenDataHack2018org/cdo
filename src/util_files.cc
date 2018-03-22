
#include <sys/stat.h>
#include "util_files.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "cdoOptions.h"
#include "readline.h"
#include "text.h"
#include <ctype.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <sys/stat.h>
bool
fileExists(const char *restrict filename)
{
  bool status = false;
  struct stat buf;

  if (stat(filename, &buf) == 0)
    {
      if (S_ISREG(buf.st_mode) && buf.st_size > 0) status = true;
    }

  return status;
}

bool
userFileOverwrite(const char *restrict filename)
{
  bool status = false;

  if (!Options::silentMode && stdin_is_tty && stderr_is_tty)
    {
      fprintf(stderr, "File %s already exists, overwrite? (yes/no): ", filename);
      char line[1024];
      readline(stdin, line, 1024);
      char *pline = line;
      while (isspace((int) *pline)) pline++;
      int len = (int) strlen(pline);
      if (len == 3)
        {
          if (pline[0] == 'y' && pline[1] == 'e' && pline[2] == 's')
            status = true;
          else if (pline[0] == 'Y' && pline[1] == 'E' && pline[2] == 'S')
            status = true;
        }
      else if (len == 1)
        {
          if (pline[0] == 'y' || pline[0] == 'Y') status = true;
        }
    }

  return status;
}

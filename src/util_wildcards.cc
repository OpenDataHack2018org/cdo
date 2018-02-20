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

#include <wordexp.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <string>
#include <glob.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef strdupx
#ifndef strdup
char *strdup(const char *s);
#endif
#define strdupx strdup
/*
#define strdupx(s)			          \
({					      	  \
   const char *__old = (s);			  \
   size_t __len = strlen(__old) + 1;		  \
   char *__new = Malloc(__len);	  \
   (char *) memcpy(__new, __old, __len);	  \
})
*/
#endif

#if defined(HAVE_GLOB_H)
static int
get_glob_flags(void)
{
  int glob_flags = 0;

#if defined(GLOB_NOCHECK)
  glob_flags |= GLOB_NOCHECK;
#endif
#if defined(GLOB_TILDE)
  glob_flags |= GLOB_TILDE;
#endif

  return glob_flags;
}
#endif

int
find_wildcard(const char *string, size_t len)
{
  int status = 0;

  if (len > 0)
    {
      if (string[0] == '~') status = 1;

      if (status == 0)
        {
          for (size_t i = 0; i < len; ++i)
            if (string[i] == '?' || string[i] == '*' || string[i] == '[')
              {
                status = 1;
                break;
              }
        }
    }

  return status;
}

// used in griddes.cc
char *
expand_filename(const char *string)
{
  char *filename = NULL;

  if (find_wildcard(string, strlen(string)))
    {
#if defined(HAVE_GLOB_H)
      int glob_flags = get_glob_flags();
      glob_t glob_results;

      glob(string, glob_flags, 0, &glob_results);

      if (glob_results.gl_pathc == 1)
        filename = strdupx(glob_results.gl_pathv[0]);

      globfree(&glob_results);
#endif
    }

  return filename;
}
#if defined(HAVE_WORDEXP_H)
/* Expands all input file wildcards and removes the
 * wildcard while inserting all expanded files into argv
 */
std::vector<std::string>
expandWildCards(std::vector<std::string> argv)
{

  int flags = WRDE_UNDEF;
  wordexp_t glob_results;

  for (size_t idx = 1; idx < argv.size(); idx++)
    {
      // if argv[idx] contains wildcard (* or [?]+)
      // multiple ** are ignored
      if (argv[idx][0] != '-'
          && argv[idx].find_first_of("*?") != std::string::npos)
        {
          wordexp(argv[idx].c_str(), &glob_results, flags);
          // range based insert (glob_results.we_wordv is inserted before
          // wildcard
          argv.insert(argv.begin() + idx + 1, glob_results.we_wordv,
                      glob_results.we_wordv + glob_results.we_wordc);
          // delete wildcard
          argv.erase(argv.begin() + idx);
          wordfree(&glob_results);
        }
    }

  return argv;
}
#else
std::vector<std::string>
expandWildCards(std::vector<std::string> argv)
{
  return argv;
}
#endif

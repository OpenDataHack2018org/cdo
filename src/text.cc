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
#include "text.h"

int stdin_is_tty = 0;
int stdout_is_tty = 0;
int stderr_is_tty = 0;

int CDO_Color = 0;

void
set_text_color(FILE *fp, int attr, int fg)
{
  int bg = -1;

  if (fp == stdout && !COLOR_STDOUT) return;
  if (fp == stderr && !COLOR_STDERR) return;

  fprintf(fp, "%c[%d", 0x1B, attr);
  if (fg >= 0)
    {
      fprintf(fp, ";%d", fg + 30);
      if (bg >= 0) fprintf(fp, ";%d", bg + 40);
    }
  fprintf(fp, "m");
}

void
reset_text_color(FILE *fp)
{
  int attr = RESET;

  if (fp == stdout && !COLOR_STDOUT) return;
  if (fp == stderr && !COLOR_STDERR) return;

  fprintf(fp, "%c[%dm", 0x1B, attr);
}

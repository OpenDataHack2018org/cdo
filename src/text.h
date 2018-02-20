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
#ifndef TEXT_H
#define TEXT_H

#include <stdio.h>

enum text_mode
{
  RESET = 0,
  BRIGHT = 1,
  DIM = 2,
  UNDERLINE = 4,
  BLINK = 5,
  REVERSE = 7,
  HIDDEN = 8
};
enum text_color
{
  BLACK = 0,
  RED = 1,
  GREEN = 2,
  YELLOW = 3,
  BLUE = 4,
  MAGENTA = 5,
  CYAN = 6,
  WHITE = 7
};

extern int stdin_is_tty;
extern int stdout_is_tty;
extern int stderr_is_tty;

extern int CDO_Color;

#define COLOR_STDOUT (stdout_is_tty && CDO_Color)
#define COLOR_STDERR (stderr_is_tty && CDO_Color)

void set_text_color(FILE *fp, int attr, int fg);
void reset_text_color(FILE *fp);

#endif /* TEXT_H */

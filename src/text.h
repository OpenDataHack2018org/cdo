#ifndef  TEXT_H
#define  TEXT_H

#include <stdio.h>

enum text_mode {RESET=0, BRIGHT=1, DIM=2, UNDERLINE=4, BLINK=5, REVERSE=7, HIDDEN=8};
enum text_color {BLACK=0, RED=1, GREEN=2, YELLOW=3, BLUE=4, MAGENTA=5, CYAN=6, WHITE=7};

#define COLOR_STDOUT (stdout_is_tty && CDO_Color)
#define COLOR_STDERR (stderr_is_tty && CDO_Color)

void set_text_color(FILE *fp, int attr, int fg);
void reset_text_color(FILE *fp);

#endif  /* TEXT_H */

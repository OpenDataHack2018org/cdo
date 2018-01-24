#ifndef UTIL_STRING_H
#define UTIL_STRING_H

#include "util_string.h"

#define  ADD_PLURAL(n)  ((n)!=1 ? "s" : "")

std::string string2lower(std::string str);
void strtolower(char *str);
void strtoupper(char *str);

#endif

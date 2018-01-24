#include <string>
#include <stdlib.h>

std::string string2lower(std::string str)
{
  std::string lower_case_string = str;
  for(char c : str) c = tolower(c); 
  return lower_case_string;
}

void strtolower(char *str)
{
  if ( str )
    for ( size_t i = 0; str[i]; ++i )
      str[i] = (char)tolower((int)str[i]);
}

void strtoupper(char *str)
{
  if ( str )
    for ( size_t i = 0; str[i]; ++i )
      str[i] = (char)toupper((int)str[i]);
}

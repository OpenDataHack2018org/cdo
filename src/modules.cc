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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "error.h"
#include "modules.h"
#include <cdi.h>
#include "util_string.h"

#include <dirent.h>
#include <dlfcn.h>
#include <regex>
#include <set>
#include <string>
// for std::sort()
#include <algorithm>

/**
 * @param a pointer to a string/substring
 * @param b pointer to a string/substring
 * @param alen length of string a
 * @param blen lenght of string b
 * @retval true if a is similar to b
 * @retval false if a is not similar to b
 *
 * Recursive function for finding substrings of a operator name that match other
 * operators.
 */

static bool
similar(const char *a, const char *b, unsigned long alen, unsigned long blen)
{
  if (alen > 2 && blen > 2 && strstr(b, a)) return true;

  while (*a && *b && *a == *b)
    {
      a++;
      b++;
    }
  if (!*a && !*b) return true;
  /*
    printf("%d %d %s %s\n", alen, blen, a, b);
  */
  if (alen >= 2 && blen >= 1 && *a && similar(a + 1, b, alen - 2, blen - 1)) return true;

  if (alen >= 1 && blen >= 2 && *b && similar(a, b + 1, alen - 1, blen - 2)) return true;

  return false;
}

/**
 * @param original string tested for similarity to \p other
 * @param other string that \p original will be compared to
 * @retval true if original and other are similar
 * @retval false if not
 *
 * Wrapper function for #similar() to parse c++ strings to c strings
 */
static bool
similar(std::string original, std::string other)
{
  return (similar(original.c_str(), other.c_str(), original.size(), other.size()));
}

/**
 * @param operatorName operator name
 * @retval true if #modules_map contains \p operatorName
 * @retval false if not
 */
static bool
operator_name_exists(std::string operatorName)
{
  if (modules_map.find(operatorName) != modules_map.end())
    {
      return true;
    }
  if (aliases.find(operatorName) != aliases.end())
    {
      return true;
    }
  return false;
}

/**
 * @param moduleName module name
 * @retval true if #modules contains \a moduleName
 * @retval false if not
 */
static bool
module_map_contains(std::string moduleName)
{
  if (modules.find(moduleName) != modules.end())
    {
      return true;
    }
  else
    {
      Error("Module %s not found", moduleName.c_str());
    }
  return false;
}

/***
 * function for finding similar operator names for the given string
 * @param operatorName operator name to find similar operators for
 * @returns A string with all found names. The string is seqmented into lines
 * with a max lenght of 75 characters
 */
static std::string
find_similar(std::string operatorName)
{
  std::string found_similar_operators = "";
  unsigned long lines = 1;
  unsigned long line_length = 105;
  if (operatorName != "")
    {
      // Searching for simlar operator names in operator to module map
      for (auto str : modules_map)
        {
          if (similar(string2lower(operatorName), str.first))
            {
              if (found_similar_operators.size() + str.first.size() > lines * line_length)
                {
                  found_similar_operators += "\n";
                  lines++;
                }
              found_similar_operators += str.first;
              found_similar_operators += " ";
            }
        }
      // Searching for similar operator names in aliases to original map
      for (auto str : aliases)
        {
          if (similar(string2lower(operatorName), str.first))
            {
              if (found_similar_operators.size() + str.first.size() > lines * line_length)
                {
                  found_similar_operators += "\n";
                  lines++;
                }
              found_similar_operators += str.first;
              found_similar_operators += " ";
            }
        }
    }
  return found_similar_operators;
}

/**
 * @param operatorName operator name.
 * @retval true if \p operatorName exists.
 * @retval false if \p operatorName is not in #modules_map
 *
 * Checks if given \p operatorName is in #modules_map. Else returns false.

 * Checks if \p operatorName is not a file.

 * If no matching operator is found checks for similar operators using
 find_similar().
 *
 *  \note If \p operatorName is a file name the program will exit.
 */
static bool
check_operator(std::string operatorName)
{
  if (operator_name_exists(operatorName))
    {
      return true;
    }
  else if (operatorName == "")
    Error("Operator name missing!");

  else
    {
      // Checking if the operatorname is an existing file name
      FILE *fp = fopen(operatorName.c_str(), "r");
      if (fp)
        {
          fclose(fp);
          fprintf(stderr, "Use commandline option -h for help.");
          Error("Operator missing, %s is a file on disk!", operatorName.c_str());
        }
      // Operator is no filename
      // Checking for similar operators
      fprintf(stderr, "Operator >%s< not found!\n", operatorName.c_str());
      fprintf(stderr, "Similar operators are:\n");
      std::string found_similar_operators = find_similar(operatorName);

      if (found_similar_operators.size() > 0)
        {
          std::cerr << found_similar_operators << std::endl;
        }
      else
        {
          fprintf(stderr, "(not found)\n");
        }
      exit(EXIT_FAILURE);
    }
  return false;
}

/***
 * Adds a module and its operators to cdo.
 * Adds the module to modules
 * Adds the operators of modules to modules_map
 * @param new_module newly constructed module
 * @note: if an error happens while adding the new module cdo will exit.
 */
void
add_module(std::string module_name, module_t new_module)
{
  if (modules.find(module_name) == modules.end())
    {
      modules[module_name] = new_module;
      for (std::string operatorName : new_module.operators)
        {
          // if the operator name is not already in the map or in the aliases
          if (!operator_name_exists(operatorName))
            {
              modules_map[operatorName] = module_name;
            }
          else
            {
              Error("Tried to add operator but the operator name already exists");
            }
        }
    }
  else
    {
      Error("Module %s name already exists", module_name.c_str());
    }
}

/**
 * adds an key value pair to #modules_map with alias as key and originals name
 * as value
 * @param alias new alias to be added
 * @param original original operator name
 */
int
add_alias(std::string alias, std::string original)
{
  auto iter_original = modules_map.find(original);
  auto iter_alias = aliases.find(alias);

  if (iter_alias != aliases.end())
    {
      Warning("alias %s could not be added: it already exists", alias.c_str());
      return -1;
    }

  if (iter_original == modules_map.end())
    {
      Error("alias %s could not be added: operator %s does not exist", alias.c_str(), original.c_str());
      return -2;
    }
  if (modules_map.find(alias) != modules_map.end())
    {
      Error("alias %s could not be added: alias name already exists as an operator", alias.c_str());
    }
  aliases[alias] = original;

  return 0;
}
/* clang-format on */

/***
 * returns the module of given operator
 *
 * parameters:
 *      std::string operatorName -> name of the operator
 * return value:
 *      std::string -> name of the operators module
 */
static std::string
get_module_name_to(std::string operatorName)
{
  // if not true the programm will exit in function check_operator
  if (check_operator(operatorName))
    {
      if (modules_map.count(operatorName) > 0)
        {
          return modules_map[operatorName];
        }
      else if (aliases.count(operatorName) > 0)
        {
          return modules_map[aliases[operatorName]];
        }
    }
  else
    {
      // TODO: see if depricated since no operator can be added without a
      // module
      Error("Module for >%s< not installed", operatorName.c_str());
    }
  // Only to quell the warning. Function should never reach this.
  return "";
}

/**
 * @fn void *(*operatorModule(const char *operatorName))(void *)

 * returns the function of the module of \a operatorName
 * @param operatorName name of the operator
 * @returns the function of the #module_t of \a operatorName
 */
void *(*operatorModule(const char *operatorName))(void *)
{
  std::string operator_name = std::string(operatorName);
  return modules[get_module_name_to(operator_name)].func;
}

/***
 * returns help for given operator name
 * if there is no help a empty vector will be returned
 * @param operatorName operator name
 * @return vector of strings containing the help
 */
std::vector<std::string>
operatorHelp(std::string operatorName)
{
  std::string operator_name = std::string(operatorName);
  return modules[get_module_name_to(operator_name)].help;
}

/***
 * Returns the number of input streams a operator takes.
 * returns -1 for a unlimited number of input streams.
 * @param operatorName operator name
 * @retval short
 */
int
operatorStreamInCnt(const char *operatorName)
{
  std::string operator_name = std::string(operatorName);
  return modules[get_module_name_to(operator_name)].streamInCnt;
}

/***
 * Returns the number of output streams a given operator has.
 * returns -1 if obase is used
 * @param operatorName operator name
 * @return 1 for CDI_REAL, 2 for CDI_COMP (complex), 3 for CDI_REAL and CDI_COMP
 * @retval short
 */
int
operatorStreamOutCnt(const char *operatorName)
{
  std::string operator_name = std::string(operatorName);
  return modules[get_module_name_to(operator_name)].streamOutCnt;
}

/***
 * Returns the number type this operator works with
 * @param operatorName operator name
 * @reval short
 */
int
operatorStreamNumber(const char *operatorName)
{
  std::string operator_name = std::string(operatorName);
  return modules[get_module_name_to(operator_name)].number;
}

/***
 * Creates a sorted vector with all operator names and alisases excluding all
 * modules that are marked as internal
 * @return a sorted std::vector containing all operator names and aliases
 * excluding all operators which modules are marked as internal
 */
static std::vector<std::string>
get_sorted_operator_name_list()
{

  std::vector<std::string> names;
  for (std::pair<std::string, std::string> operator_module_names_pair : modules_map)
    {
      if (modules[operator_module_names_pair.second].mode == 1)
        {
          names.push_back(operator_module_names_pair.first);
        }
    }
  // adding operators names from alias_map
  for (std::pair<std::string, std::string> alias : aliases)
    {
      names.push_back(alias.first);
    }
  std::sort(names.begin(), names.end());
  return names;
}

std::vector<std::string>
get_no_output_operator_list()
{
  std::vector<std::string> names;
  for (std::pair<std::string, std::string> operator_module_names_pair : modules_map)
    {
      if (modules[operator_module_names_pair.second].mode == 1 && modules[operator_module_names_pair.second].streamOutCnt == 0)
        {
          names.push_back(operator_module_names_pair.first);
        }
    }
  // adding operators names from alias_map
  std::string original;
  for (std::pair<std::string, std::string> alias : aliases)
    {
      original = alias.second;
      if (modules[modules_map[original]].mode == 1 && modules[modules_map[original]].streamOutCnt == 0)
        {
          names.push_back(alias.first);
        }
    }
  std::sort(names.begin(), names.end());
  return names;
}

void
operatorPrintAll(void)
{
  int number_of_chars = 0;
  std::string tab = "   ";
  int tab_width = tab.size();
  // using a set because it sorts the operators alphabetically on its own
  std::vector<std::string> sorted_operator_names = get_sorted_operator_name_list();

  std::cout << tab;
  for (auto operatorName : sorted_operator_names)
    {
      if (number_of_chars > 85)
        {
          number_of_chars = tab_width;
          std::cerr << std::endl << tab;
        }

      std::cerr << " " << operatorName;
      number_of_chars += 1 + operatorName.size();
    }

  std::cerr << std::endl;
}

#ifdef CUSTOM_MODULES
/***
  loads all custom modules in a specified folder
  @param folder_path custom modules folder
*/
#ifdef CUSTOM_MODULES
void
load_custom_modules(std::string folder_path)
{
  DIR *dir = opendir(folder_path.c_str());
  std::string file_path;
  std::regex library_regex("(.*\\.so)");
  if (dir != NULL)
    {
      struct dirent *ent = readdir(dir);
      while (ent != NULL)
        {
          if (std::regex_match(ent->d_name, library_regex))
            {
              file_path = folder_path + "/" + ent->d_name;
              load_custom_module(file_path);
            }
          ent = readdir(dir);
        }
    }
  else
    {
      std::cerr << "Could not find " << folder_path << "for loading custom modules" << std::endl;
    }
}

/***
 * Loads a custom module from given path.
 * Modules must contain a (TODO: rename function) init_custom_module function
 * Program exits if a module could not be loaded.
 * @param module file path
 */
void
load_custom_module(std::string file_path)
{
  void (*init_custom_module)();
  void *lib_handle = dlopen(file_path.c_str(), RTLD_LAZY);
  custom_modules_lib_handles.push_back(lib_handle);
  if (!lib_handle)
    {
      std::cerr << "Cannot open library: " << dlerror() << std::endl;
      return;
    }

  dlerror();
  init_custom_module = (void (*)()) dlsym(lib_handle, "init_custom_module");
  const char *dlsym_error = dlerror();

  if (dlsym_error)
    {
      std::cerr << "Cannot load symbol 'init_custom_module': " << dlsym_error << std::endl;
      dlclose(lib_handle);
      return;
    }

  init_custom_module();
  std::cout << "loaded custom module from '" << file_path << "'" << std::endl;
}
#endif
/***
  closes the handles for the loaded custum modules
*/
void
close_library_handles()
{
  for (unsigned long i = 0; i < custom_modules_lib_handles.size(); i++)
    {
      dlclose(custom_modules_lib_handles[i]);
    }
}
#endif

// helper function for setting the spacing in operatorPrintList
std::string
get_spacing_for(int p_space, std::string str)
{
  std::string spacing = "";
  for (int i = str.size(); i <= p_space; i++)
    {
      spacing += " ";
    }
  return spacing;
}
std::string
get_operator_description(std::string p_current_op_name, std::vector<std::string> help)
{
  std::string description = "";
  unsigned long cur_help_idx;
  bool help_contains_name;
  std::string line;
  unsigned long operator_section = 0;

  cur_help_idx = 0;
  // search for operator section
  while (operator_section == 0 && cur_help_idx < help.size() - 1)
    {
      line = help[++cur_help_idx];
      if (line.find("OPERATORS") != std::string::npos)
        {
          operator_section = cur_help_idx;
        }
    }
  // if no operator section is found
  if (operator_section == 0)
    {
      cur_help_idx = 0;
      line = help[0];
      std::string name_section = help[0];
      help_contains_name = false;
      // search for the operator name in the description
      while (!line.empty())
        {
          line = help[++cur_help_idx];
          if (line.find(p_current_op_name) != std::string::npos)
            {
              help_contains_name = true;
            }
          name_section += line;
        }
      // if the name was found save description for later use
      if (help_contains_name)
        {
          description = name_section.substr(name_section.find_first_of('-') + 2, name_section.size());
        }
    }
  else
    {

      line = help.at(++operator_section);
      // search the operator section for current operator line
      while (line.find(p_current_op_name + " ") == std::string::npos && !line.empty() && operator_section < help.size() - 1)
        {
          line = help.at(++operator_section);
        }
      // if operator line found save description for later use
      if (!line.empty() && line.find(p_current_op_name + " ") != std::string::npos)
        {
          auto op_name_start = line.find_first_not_of(" \t");

          description = line.substr(line.find_first_not_of(" \t", op_name_start + p_current_op_name.size()), line.size());
        }
    }
  return description;
}
/***
 * Prints all operator names and their short descriptions
 * Aliases are listed and point to their original operator name.
 * If the module is not documented the description is empty
 * If a module has only one operator the short module description is listed
 * If the operator is not documented the description is empty
 */
void
operatorPrintList(bool print_no_output)
{
  std::vector<std::string> output_list;
  if (print_no_output)
    {
      output_list = get_no_output_operator_list();
    }
  else
    {
      output_list = get_sorted_operator_name_list();
    }

  unsigned long list_length = output_list.size();
  module_t *current_module;

  // help variables

  for (unsigned long out_list_idx = 0; out_list_idx < list_length; out_list_idx++)
    {
      std::string current_op_name = output_list[out_list_idx];
      current_module = &modules[get_module_name_to(current_op_name)];
      if (aliases.find(current_op_name) != aliases.end())
        {

          output_list[out_list_idx] += std::string(get_spacing_for(16, current_op_name) + "--> " + aliases[current_op_name]);
        }
      else if (!current_module->help.empty())
        {
          // add spaceing and saving output line to the output list
          std::string description = get_operator_description(current_op_name, current_module->help);
          output_list[out_list_idx] += get_spacing_for(16, current_op_name) + description;
        }
      std::string in_out_info
          = "(" + std::to_string(current_module->streamInCnt) + "|" + std::to_string(current_module->streamOutCnt) + ")";
      output_list[out_list_idx] += get_spacing_for(90, output_list[out_list_idx]) + in_out_info;
    }
  // print generated output list
  for (std::string str : output_list)
    {
      std::cout << str << std::endl;
    }
}

bool
is_alias(const char *operatorName)
{
  return (aliases.find(std::string(operatorName)) != aliases.end());
}

const char *
get_original(const char *operatorName)
{
  char *original = NULL;
  if (is_alias(operatorName))
    {
      std::string opName = aliases[std::string(operatorName)];
      original = (char *) malloc(opName.size());
      strcpy(original, opName.c_str());
    }
  else
    {
      return operatorName;
    }
  return original;
}

module_t &
getModule(const std::string &operator_name)
{
  return modules[get_module_name_to(operator_name)];
}

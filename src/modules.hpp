
  /***

   Copyright (C) 2003-2017 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
   See COPYING file for copying and redistribution conditions.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; version 2 of the License.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
 */

#ifndef _MODULES_H
#define _MODULES_H
#define MAX_MOD_OPERATORS 128 // maximum number of operators for a module

#include <iostream>
#include <map>
#include <stdbool.h>
#include <vector>

/***
  type definition for module functions loaded from a custom module
  */
typedef void (*dyn_oper_t)(void *arg);

typedef struct {
    void *(*func)(void *);               // Module
    std::vector<std::string> help;                   // Help
    std::vector<const char *> operators; // Operator names
    short mode;                          // Module mode: 0:intern 1:extern
    short number;                        // Allowed number type
    short streamInCnt;                   // Number of input streams
    short streamOutCnt;                  // Number of output streams
} modules_t;

extern int NumModules;
extern int NumAliases;

/***
  vector for library handles for loaded custom modules
  */
static std::vector<void *> custom_modules_lib_handles;

/**
 * Maps operator names to their module names
 */
static std::map<std::string, std::string> modules_map;

/**
* Contains added modules as values and their names as key
*/
static std::map<std::string, modules_t> modules;

/***
 * Key: operator alias / Value: operator original name
 */
static std::map<std::string, std::string> aliases;

void *(*operatorModule(std::string operatorName))(void *);
void *(*operatorModule(const char *operatorName))(void *);
void init_modules();

void init_aliases();
int add_alias(std::string alias, std::string original);
void add_module(std::string module_name, modules_t new_module);
std::string get_module_name_to(std::string operatorName);

std::vector<std::string> operatorHelp(std::string operatorName);
int operatorStreamInCnt(const char *operatorName);
int operatorStreamOutCnt(const char *operatorName);
int operatorStreamNumber(const char *operatorName);

void operatorPrintAll(void);
void operatorPrintList(bool print_no_output);
#ifdef CUSTOM_MODULES
void load_custom_module(std::string path);
void load_custom_modules(std::string folder_path);
void close_library_handles();
#endif

#endif /* _MODULES_H */

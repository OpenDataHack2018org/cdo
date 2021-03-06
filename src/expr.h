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

#include <stdio.h>
#include <stdbool.h>

#ifndef fileno
int fileno(FILE *stream);
#endif

extern int CDO_parser_errorno;

enum
{
  CTIMESTEP,
  CDATE,
  CTIME,
  CDELTAT,
  CDAY,
  CMONTH,
  CYEAR,
  CSECOND,
  CMINUTE,
  CHOUR,
  CLEN
};

typedef enum { typeCon, typeVar, typeFun, typeFun1c, typeOpr, typeCom } nodeEnum;

// commands
struct comNodeType
{
  char *cname;  // command name
  char *vname;  // variable name
};

// constants
struct conNodeType
{
  double value;  // value of constant
};

// variables
struct varNodeType
{
  char *nm;  // variable name
};

// 1c functions
struct fun1cNodeType
{
  char *name;              // function name
  double value;            // value of constant
  struct nodeTypeTag *op;  // operand
};

// functions
struct funNodeType
{
  char *name;              // function name
  struct nodeTypeTag *op;  // operand
};

// operators
struct oprNodeType
{
  int oper;                   // operator
  int nops;                   // number of operands
  struct nodeTypeTag *op[1];  // operands (expandable)
};

enum {PARAM_VAR, PARAM_CONST};

// parameter
struct paramType
{
  int type;
  bool select;
  bool remove;
  bool lmiss;
  int coord;
  int gridID;
  int zaxisID;
  int datatype;
  int steptype;
  size_t ngp;
  size_t nlat;
  size_t nlev;
  size_t nmiss;
  char *name;
  char *longname;
  char *units;
  double missval;
  double *data;
  double *weight;
};

typedef struct nodeTypeTag
{
  bool ltmpobj;
  paramType param;

  nodeEnum type;  // type of node

  // union must be last entry in nodeType
  // because operNodeType may dynamically increase
  union
  {
    comNodeType com;      // commands
    conNodeType con;      // constants
    varNodeType var;      // variables
    funNodeType fun;      // functions
    fun1cNodeType fun1c;  // functions
    oprNodeType opr;      // operators
  } u;
} nodeType;

struct coordType
{
  bool needed;
  int coord;
  int cdiID;
  size_t size;
  char *units;
  char *longname;
  double *data;
};

struct parseParamType
{
  bool init;
  bool debug;
  bool *needed;
  int maxparams;
  int nparams;
  int nvars1;
  int ncoords;
  int maxcoords;
  int tsID;
  int pointID;
  int zonalID;
  int surfaceID;
  coordType *coords;
  paramType *params;
};

typedef union
{
  double cvalue;   // constant value
  char *varnm;     // variable name
  char *fname;     // function name
  nodeType *nPtr;  // node pointer
} yysType;

#define YYSTYPE yysType
#define YY_EXTRA_TYPE parseParamType *

#define YY_DECL int yylex(YYSTYPE *yylval_param, parseParamType *parse_arg, void *yyscanner)
YY_DECL;

int yyparse(parseParamType *parse_arg, void *);
void yyerror(void *parse_arg, void *scanner, const char *errstr);

int yylex_init(void **);
int yylex_destroy(void *);
void yyset_extra(YY_EXTRA_TYPE, void *);

nodeType *expr_run(nodeType *p, parseParamType *parse_arg);
int params_get_coordID(parseParamType *parse_arg, int coord, int cdiID);

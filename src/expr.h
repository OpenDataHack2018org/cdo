#include <stdio.h>
#include <stdbool.h>

#ifdef __cplusplus
#ifndef register
#define register
#endif
#ifndef fileno
int fileno(FILE *stream);
#endif
#endif


typedef enum { typeCon, typeVar, typeFun, typeOpr } nodeEnum;

// constants
typedef struct {
  double value;               // value of constant
} conNodeType;

// variables
typedef struct {
  char *nm;                   // variable name
} varNodeType;

// functions
typedef struct {
  char *name;                 // function name
  struct nodeTypeTag *op;     // operand
} funNodeType;

// operators
typedef struct {
  int oper;                   // operator             
  int nops;                   // number of operands   
  struct nodeTypeTag *op[1];  // operands (expandable)
} oprNodeType;

// parameter
typedef struct {
  bool    select;
  int     coord;
  int     gridID;
  int     zaxisID;
  int     steptype;
  size_t  ngp;
  size_t  nlev;
  size_t  nmiss;
  char   *name;
  char   *longname;
  char   *units;
  double  missval;
  double *data;
} paramType;


typedef struct nodeTypeTag {
  bool      ltmpobj;
  paramType param;

  nodeEnum  type;             // type of node

  // union must be last entry in nodeType
  // because operNodeType may dynamically increase
  union {
    conNodeType con;          // constants
    varNodeType var;          // variables
    funNodeType fun;          // functions
    oprNodeType opr;          // operators
  } u;
} nodeType;


typedef struct {
  bool       init;
  bool       debug;
  bool      *needed;
  int        maxparams;
  int        nparams;
  int        nvars1;
  int        pointID;
  int        surfaceID;
  paramType  param2;
  paramType *params;
} parse_param_t;


typedef union{
  double    cvalue;           // constant value
  char     *varnm;            // variable name 
  char     *fname;            // function name 
  nodeType *nPtr;             // node pointer  
} stype_t;


#define YYSTYPE        stype_t
#define YY_EXTRA_TYPE  parse_param_t *

#define YY_DECL int yylex(YYSTYPE *yylval_param, parse_param_t *parse_arg, void *yyscanner)
YY_DECL;

int  yyparse(parse_param_t *parse_arg, void*);
void yyerror(void *parse_arg, void *scanner, const char *errstr);

int  yylex_init(void **);
int  yylex_destroy(void *);
void yyset_extra(YY_EXTRA_TYPE, void *);

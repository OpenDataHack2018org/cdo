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
  int gridID;
  int zaxisID;
  int steptype;
  int ngp;
  int nlev;
  int nmiss;
  char *name;
  char *longname;
  char *units;
  double missval;
  double *data;
} paramType;

typedef struct nodeTypeTag {
  bool ltmpvar;
  paramType param;

  nodeEnum type;              // type of node

  // union must be last entry in nodeType
  // because operNodeType may dynamically increase
  union {
    conNodeType con;          // constants
    varNodeType var;          // variables
    funNodeType fun;          // functions
    oprNodeType opr;          // operators
  } u;
} nodeType;

#define MAX_VARS 1024

typedef struct {
  bool   init;
  bool   debug;
  bool  *needed;
  int    maxparams;
  int    nparams;
  int    vlistID2;
  int    nvars1, nvars2;
  int    nmiss[MAX_VARS];
  int    gridID2;
  int    zaxisID2;
  int    tsteptype2;
  double missval2;
  double **vardata2;
  paramType *params;
} parse_param_t;


typedef union{
  double cvalue;              // constant value
  char *varnm;                // variable name 
  char *fname;                // function name 
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

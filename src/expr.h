#include <stdio.h>
int fileno(FILE *stream);

#ifndef strdupx
#define strdupx(s)			          \
({					      	  \
   const char *__old = (s);			  \
   size_t __len = strlen(__old) + 1;		  \
   char *__new = (char *) malloc(__len);	  \
   (char *) memcpy(__new, __old, __len);	  \
})
#endif


typedef enum { typeCon, typeVar, typeFun, typeOpr } nodeEnum;

/* constants */
typedef struct {
  double value;               /* value of constant */
} conNodeType;

/* variables */
typedef struct {
  char *nm;                   /* variable name */
} varNodeType;

/* functions */
typedef struct {
  char *name;                 /* function name */
  struct nodeTypeTag *op;     /* operand */
} funNodeType;

/* operators */
typedef struct {
  int oper;                   /* operator */
  int nops;                   /* number of operands */
  struct nodeTypeTag *op[1];  /* operands (expandable) */
} oprNodeType;

typedef struct nodeTypeTag {
  int tmpvar;
  int gridID, zaxisID;
  int nmiss;
  double missval;
  double *data;
  nodeEnum type;              /* type of node */

  /* union must be last entry in nodeType */
  /* because operNodeType may dynamically increase */
  union {
    conNodeType con;          /* constants   */
    varNodeType var;          /* variables   */
    funNodeType fun;          /* functions   */
    oprNodeType opr;          /* operators   */
  } u;
} nodeType;


typedef struct{ /* prs_sct */
  int    vlistID1, vlistID2;
  int    nvars1, nvars2;
  int    nmiss[1024];
  int    varID[1024];
  int    var_needed[1024];
  char   *var[1024];
  int    init;
  int    debug;
  double **vardata1, **vardata2;
} prs_sct;

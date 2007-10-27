%{
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "expr.h"

/* prototypes */
nodeType *opr(int oper, int nops, ...);
nodeType *var(char *nm);
nodeType *con(double value);
nodeType *fun(char *fname, nodeType *p);

void freeNode(nodeType *p);
int ex(nodeType *p, prs_sct *prs_arg);

void yyerror(char *s);

/* Get YYSTYPE prior to prototyping scanner */
#include "expr_yacc.h" /* expr_yacc.h (y.tab.h) is produced from expr_yacc.y by parser generator */

#define YY_DECL int yylex(YYSTYPE *yylval, prs_sct *prs_arg)

YY_DECL;

/* Bison manual p. 60 describes how to call yyparse() with arguments */
#define YYPARSE_PARAM prs_arg
#define YYLEX_PARAM   &yylval, (prs_sct *) prs_arg

  /* #define YYPURE 1 *//* ??? */

%}

%union {
    double cvalue;              /* constant value */
    char *varnm;                /* variable name */
    char *fname;                /* function name */
    nodeType *nPtr;             /* node pointer */
};

%token <cvalue> CONSTANT
%token <varnm>  VARIABLE
%token <fname>  FUNCTION

%left GE LE EQ NE '>' '<' '='
%left '+' '-'
%left '*' '/'
%right  '^'
%nonassoc UMINUS

%type <nPtr> stmt expr stmt_list

%%

program:
        function                  { return(0); }
        ;

function:
          function stmt           { ex($2, (prs_sct *) prs_arg); freeNode($2); }
        | /* NULL */
        ;

stmt:
          ';'                     { $$ = opr(';', 2, NULL, NULL); }
        | expr ';'                { $$ = $1; }
        | VARIABLE '=' expr ';'   { $$ = opr('=', 2, var($1), $3); }
        | '{' stmt_list '}'       { $$ = $2; }
        ;

stmt_list:
          stmt                    { $$ = $1; }
        | stmt_list stmt          { $$ = opr(';', 2, $1, $2); }
        ;

expr:
          CONSTANT                { $$ = con($1); }
        | VARIABLE                { $$ = var($1); }
        | '-' expr %prec UMINUS   { $$ = opr(UMINUS, 1, $2); }
        | expr '+' expr           { $$ = opr('+', 2, $1, $3); }
        | expr '-' expr           { $$ = opr('-', 2, $1, $3); }
        | expr '*' expr           { $$ = opr('*', 2, $1, $3); }
        | expr '/' expr           { $$ = opr('/', 2, $1, $3); }
        | expr '<' expr           { $$ = opr('<', 2, $1, $3); }
        | expr '>' expr           { $$ = opr('>', 2, $1, $3); }
        | expr '^' expr           { $$ = opr('^', 2, $1, $3); }
        | expr GE expr            { $$ = opr(GE, 2, $1, $3); }
        | expr LE expr            { $$ = opr(LE, 2, $1, $3); }
        | expr NE expr            { $$ = opr(NE, 2, $1, $3); }
        | expr EQ expr            { $$ = opr(EQ, 2, $1, $3); }
        | '(' expr ')'            { $$ = $2; }
        | FUNCTION '(' expr ')'   { $$ = fun($1, $3); }
        ;

%%

#define SIZEOF_NODETYPE ((char *)&p->u.con - (char *)p)

nodeType *con(double value)
{
  nodeType *p = NULL;
  size_t nodeSize;

  /* allocate node */
  nodeSize = SIZEOF_NODETYPE + sizeof(conNodeType);
  if ((p = (nodeType *) malloc(nodeSize)) == NULL)
    yyerror("out of memory");

  /* copy information */
  p->type = typeCon;
  p->u.con.value = value;
    
  return p;
}

nodeType *var(char *nm)
{
  nodeType *p = NULL;
  size_t nodeSize;

  /* allocate node */
  nodeSize = SIZEOF_NODETYPE + sizeof(varNodeType);
  if ((p = (nodeType *) malloc(nodeSize)) == NULL)
    yyerror("out of memory");

  /* copy information */
  p->type = typeVar;
  p->u.var.nm = strdupx(nm);

  return p;
}

nodeType *fun(char *fname, nodeType *op)
{
  nodeType *p = NULL;
  size_t nodeSize;

  /* allocate node */
  nodeSize = SIZEOF_NODETYPE + sizeof(funNodeType);
  if ((p = (nodeType *) malloc(nodeSize)) == NULL)
    yyerror("out of memory");

  /* copy information */
  p->type = typeFun;
  p->u.fun.name = strdupx(fname);
  p->u.fun.op   = op;

  return p;
}

nodeType *opr(int oper, int nops, ...)
{
  va_list ap;
  nodeType *p = NULL;
  size_t nodeSize;
  int i;

  /* allocate node */
  nodeSize = SIZEOF_NODETYPE + sizeof(oprNodeType) + (nops - 1)*sizeof(nodeType*);
  if ((p = (nodeType *) malloc(nodeSize)) == NULL)
    yyerror("out of memory");

  /* copy information */
  p->type = typeOpr;
  p->u.opr.oper = oper;
  p->u.opr.nops = nops;
  va_start(ap, nops);
  for (i = 0; i < nops; i++)
    p->u.opr.op[i] = va_arg(ap, nodeType*);
  va_end(ap);

  return p;
}

void freeNode(nodeType *p)
{
  int i;

  if ( ! p ) return;

  if (p->type == typeOpr)
    {
      for (i = 0; i < p->u.opr.nops; i++)
	freeNode(p->u.opr.op[i]);
    }
  
  free (p);
}

void yyerror(char *s)
{
  fprintf(stdout, "%s\n", s);
}

/*
int main(void)
{
  int i;
  static char fexpr[] = "nvar = q*(geosp+234.56); xx = geosp+999-log(aps);";

  prs_sct prs_arg;

  printf("%s\n", fexpr);

  yy_scan_string(fexpr);

  prs_arg.nvar = 0;
  prs_arg.init = 1;
  prs_arg.debug = 1;

  yyparse((void *)&prs_arg);

  for ( i = 0; i < prs_arg.nvar; i++ )
    printf("vars %d %s\n", i, prs_arg.var[i]);

  yy_scan_string(fexpr);

  prs_arg.init = 0;

  yyparse((void *)&prs_arg);

  for ( i = 0; i < prs_arg.nvar; i++ )
    printf("vars %d %s\n", i, prs_arg.var[i]);

  return 0;
}
*/

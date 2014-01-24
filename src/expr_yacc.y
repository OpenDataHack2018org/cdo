%{
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "expr.h"
#include "expr_yacc.h" /* expr_yacc.h (y.tab.h) is produced from expr_yacc.y by parser generator */

/* Bison manual p. 60 describes how to call yyparse() with arguments */
/* #define YYPARSE_PARAM parse_arg */
/* #define YYLEX_PARAM   ((parse_parm_t *) parse_arg, void *yyscanner) */

  /* #define YYPURE 1 *//* ??? */

/* prototypes */
nodeType *expr_opr(int oper, int nops, ...);
nodeType *expr_var(char *nm);
nodeType *expr_con(double value);
nodeType *expr_fun(char *fname, nodeType *p);

void freeNode(nodeType *p);
int expr_run(nodeType *p, parse_parm_t *parse_arg);

%}

%pure_parser
%parse-param {parse_parm_t *parse_arg}
%parse-param {void *scanner}
%lex-param {parse_parm_t *parse_arg}
%lex-param {yyscan_t *scanner}


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
          function stmt           { expr_run($2, (parse_parm_t *) parse_arg); freeNode($2); }
        | /* NULL */
        ;

stmt:
          ';'                     { $$ = expr_opr(';', 2, NULL, NULL); }
        | expr ';'                { $$ = $1; }
        | VARIABLE '=' expr ';'   { $$ = expr_opr('=', 2, expr_var($1), $3); }
        | VARIABLE ';'            { $$ = expr_opr('=', 2, expr_var($1), expr_var($1)); } /* conflicts: 1 shift/reduce */
        | '{' stmt_list '}'       { $$ = $2; }
        ;

stmt_list:
          stmt                    { $$ = $1; }
        | stmt_list stmt          { $$ = expr_opr(';', 2, $1, $2); }
        ;

expr:
          CONSTANT                { $$ = expr_con($1); }
        | VARIABLE                { $$ = expr_var($1); }
        | '-' expr %prec UMINUS   { $$ = expr_opr(UMINUS, 1, $2); }
        | expr '+' expr           { $$ = expr_opr('+', 2, $1, $3); }
        | expr '-' expr           { $$ = expr_opr('-', 2, $1, $3); }
        | expr '*' expr           { $$ = expr_opr('*', 2, $1, $3); }
        | expr '/' expr           { $$ = expr_opr('/', 2, $1, $3); }
        | expr '<' expr           { $$ = expr_opr('<', 2, $1, $3); }
        | expr '>' expr           { $$ = expr_opr('>', 2, $1, $3); }
        | expr '^' expr           { $$ = expr_opr('^', 2, $1, $3); }
        | expr GE expr            { $$ = expr_opr(GE, 2, $1, $3); }
        | expr LE expr            { $$ = expr_opr(LE, 2, $1, $3); }
        | expr NE expr            { $$ = expr_opr(NE, 2, $1, $3); }
        | expr EQ expr            { $$ = expr_opr(EQ, 2, $1, $3); }
        | '(' expr ')'            { $$ = $2; }
        | FUNCTION '(' expr ')'   { $$ = expr_fun($1, $3); }
        ;

%%

#define SIZEOF_NODETYPE ((char *)&p->u.con - (char *)p)

nodeType *expr_con(double value)
{
  nodeType *p = NULL;
  size_t nodeSize;

  /* allocate node */
  nodeSize = SIZEOF_NODETYPE + sizeof(conNodeType);
  if ((p = (nodeType *) malloc(nodeSize)) == NULL)
    yyerror(NULL, NULL, "Out of memory");

  /* copy information */
  p->type = typeCon;
  p->u.con.value = value;
    
  return p;
}

nodeType *expr_var(char *nm)
{
  nodeType *p = NULL;
  size_t nodeSize;

  /* allocate node */
  nodeSize = SIZEOF_NODETYPE + sizeof(varNodeType);
  if ((p = (nodeType *) malloc(nodeSize)) == NULL)
    yyerror(NULL, NULL, "Out of memory");

  /* copy information */
  p->type = typeVar;
  p->u.var.nm = strdupx(nm);

  return p;
}

nodeType *expr_fun(char *fname, nodeType *op)
{
  nodeType *p = NULL;
  size_t nodeSize;

  /* allocate node */
  nodeSize = SIZEOF_NODETYPE + sizeof(funNodeType);
  if ((p = (nodeType *) malloc(nodeSize)) == NULL)
    yyerror(NULL, NULL, "Out of memory");

  /* copy information */
  p->type = typeFun;
  p->u.fun.name = strdupx(fname);
  p->u.fun.op   = op;

  return p;
}

nodeType *expr_opr(int oper, int nops, ...)
{
  va_list ap;
  nodeType *p = NULL;
  size_t nodeSize;
  int i;

  /* allocate node */
  nodeSize = SIZEOF_NODETYPE + sizeof(oprNodeType) + (nops - 1)*sizeof(nodeType*);
  if ((p = (nodeType *) malloc(nodeSize)) == NULL)
    yyerror(NULL, NULL, "Out of memory");

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

void yyerror(void *parse_arg, void *scanner, char *s)
{
  fprintf(stdout, "%s\n", s);
}
/*
int main(void)
{
  int i;
  static char fexpr[] = "nvar = q*(geosp+234.56); xx = geosp+999-log(aps);";
  void *scanner;
  int yy_scan_string(const char *str, void *scanner);

  parse_parm_t parse_arg;

  printf("%s\n", fexpr);

  yylex_init(&scanner);
  yyset_extra(&parse_arg, scanner);

  yy_scan_string(fexpr, scanner);

  parse_arg.nvars1 = 0;
  parse_arg.init  = 1;
  parse_arg.debug = 1;

  yyparse((void *)&parse_arg, scanner);

  for ( i = 0; i < parse_arg.nvars1; i++ )
    printf("vars %d %s\n", i, parse_arg.var[i]);

  yy_scan_string(fexpr, scanner);

  parse_arg.init = 0;

  yyparse((void *)&parse_arg, scanner);

  for ( i = 0; i < parse_arg.nvars1; i++ )
    printf("vars %d %s\n", i, parse_arg.var[i]);

  return 0;
}
*/

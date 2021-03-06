/* flex -oexpr_lex.c expr_lex.l */
%{
#include <string.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_E
#define M_E		2.7182818284590452354	/* e */
#endif

#ifndef M_PI
#define M_PI		3.14159265358979323846	/* pi */
#endif

#include "expr.h"
#include "expr_yacc.h"

%}

%option reentrant bison-bridge
%option noyywrap
%option nounput

/* Definitions:  LMB92 p. 153
   Definitions are named regular expressions, e.g., DGT [0-9]
   Definitions enclosed in braces in rules section, e.g. {DGT}, are interpreted literally
   DGT [0-9] Digit
   LPH [A-Za-z_] Alphabetic character
   LPHDGT [A-Za-z0-9_] Alphanumeric character
   XPN [eE][+-]?[0-9]+ Real number Exponent */
DGT [0-9]
LPH [A-Za-z_]
LPHDGT [A-Za-z0-9_]
XPN [eE][+-]?[0-9]+

%%

#.*          ; /* ignore comments */

M_PI         {
                yylval->cvalue = M_PI;
                return CONSTANT;
             }

M_E          {
                yylval->cvalue = M_E;
                return CONSTANT;
             }

{DGT}*\.?{DGT}*({XPN})?|{DGT}*({XPN}) {
               /* constant */
               /* NB: Tempted to prepend lexer expressions for floats and doubles with [+-]? 
                so that unary plus/minus is handled in lexer rather than parser.
                However, this has unintended side effects so let parser handle it for now */
                yylval->cvalue = strtod(yytext, (char **)NULL);
                return CONSTANT;
             } /* end constant */

"remove("    { return REMOVE; }

"print("    { return PRINT; }

{LPH}{LPHDGT}*/[ ]*\( {
                yylval->fname = strdup(yytext);
                return FUNCTION;
             } /* end functions */

{LPH}{LPHDGT}*\.{LPH}* {
                /* Recognize coordinates, e.g., var_nm.x (lon of var_nm) */
                yylval->varnm = strdup(yytext);
                return VARIABLE;
             }

{LPH}{LPHDGT}* {
                yylval->varnm = strdup(yytext);
                return VARIABLE;
             }

{LPHDGT}*{LPH}{LPHDGT}* {
                yylval->varnm = strdup(yytext);
                return VARIABLE;
             }


[-()=+*/;,{}^.] {
                return *yytext;
             }

"<=>"           return LEG;
">="            return GE;
"<="            return LE;
">"             return GT;
"<"             return LT;
"=="            return EQ;
"!="            return NE;
"&&"            return AND;
"||"            return OR;
"!"             return NOT;
"?"             return QUESTION;
":"             return COLON;

[ \t\n\r]+      ;       /* ignore whitespace */

.               yyerror(NULL, NULL, "Unknown character");
%%

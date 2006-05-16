%{
#include <string.h>
#include <stdlib.h>

#ifndef M_E
#define M_E		2.7182818284590452354	/* e */
#endif

#ifndef M_PI
#define M_PI		3.14159265358979323846	/* pi */
#endif

#include "expr.h"
#include "expr_yacc.h"

void yyerror(char *);

#define YY_DECL int yylex(YYSTYPE *yylval, prs_sct *prs_arg)

YY_DECL;

%}

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

{DGT}*\.?{DGT}*({XPN})?[LlDd]?|{DGT}*({XPN})[LlDd]? {
               /* constant */
               /* NB: Tempted to prepend lexer expressions for floats and doubles with [+-]? 
                so that unary plus/minus is handled in lexer rather than parser.
                However, this has unintended side effects so let parser handle it for now */
                yylval->cvalue = strtod(yytext, (char **)NULL);
                return CONSTANT;
             } /* end constant */


{LPH}{LPHDGT}*/[ ]*\( {
                yylval->fname = (char *) strdupx(yytext);
                return FUNCTION;
             } /* end functions */


{LPH}{LPHDGT}* {
                yylval->varnm = (char *) strdupx(yytext);
                return VARIABLE;
             }


[-()<>=+*/;{}^.] {
                return *yytext;
             }

">="            return GE;
"<="            return LE;
"=="            return EQ;
"!="            return NE;

[ \t\n]+        ;       /* ignore whitespace */

.               yyerror("Unknown character");
%%
int yywrap(void) {
    return 1;
}

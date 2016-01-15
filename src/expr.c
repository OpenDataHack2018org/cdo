#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <ctype.h>
#include <assert.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "field.h"
#include "expr.h"
#include "expr_yacc.h"

#define    COMPLT(x,y)  ((x) < (y) ? 1 : 0)
#define    COMPGT(x,y)  ((x) > (y) ? 1 : 0)
#define    COMPLE(x,y)  ((x) <= (y) ? 1 : 0)
#define    COMPGE(x,y)  ((x) >= (y) ? 1 : 0)
#define    COMPNE(x,y)  (IS_NOT_EQUAL(x,y) ? 1 : 0)
#define    COMPEQ(x,y)  (IS_EQUAL(x,y) ? 1 : 0)
#define   COMPLEG(x,y)  ((x) < (y) ? -1 : ((x) > (y) ? 1 : 0))
#define   COMPAND(x,y)  (IS_NOT_EQUAL(x,0) && IS_NOT_EQUAL(y,0) ? 1 : 0)
#define    COMPOR(x,y)  (IS_NOT_EQUAL(x,0) || IS_NOT_EQUAL(y,0) ? 1 : 0)
#define  MVCOMPLT(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPLT(x,y))
#define  MVCOMPGT(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPGT(x,y))
#define  MVCOMPLE(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPLE(x,y))
#define  MVCOMPGE(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPGE(x,y))
#define  MVCOMPNE(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPNE(x,y))
#define  MVCOMPEQ(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPEQ(x,y))
#define MVCOMPLEG(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPLEG(x,y))
#define MVCOMPAND(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPAND(x,y))
#define  MVCOMPOR(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPOR(x,y))

static double f_int(double x)  { return ((int)(x)); }
static double f_nint(double x) { return (round(x)); }
static double f_sqr(double x)  { return (x*x);      }

typedef struct {
  int type;
  const char *name;  // function name
  void (*func)();    // pointer to function
}
func_t;

double expr_sum(int n, double *restrict array)
{
  double sum = 0;
  for ( int i = 0; i < n; ++i ) sum += array[i];
  return sum;
}

static func_t fun_sym_tbl[] =
{
  // scalar functions
  {0, "abs",   (void (*)()) fabs},
  {0, "floor", (void (*)()) floor},
  {0, "ceil",  (void (*)()) ceil},
  {0, "int",   (void (*)()) f_int},
  {0, "nint",  (void (*)()) f_nint},
  {0, "sqr",   (void (*)()) f_sqr},
  {0, "sqrt",  (void (*)()) sqrt},
  {0, "exp",   (void (*)()) exp},
  {0, "erf",   (void (*)()) erf},
  {0, "log",   (void (*)()) log},
  {0, "log10", (void (*)()) log10},
  {0, "sin",   (void (*)()) sin},
  {0, "cos",   (void (*)()) cos},
  {0, "tan",   (void (*)()) tan},
  {0, "sinh",  (void (*)()) sinh},
  {0, "cosh",  (void (*)()) cosh},
  {0, "tanh",  (void (*)()) tanh},
  {0, "asin",  (void (*)()) asin},
  {0, "acos",  (void (*)()) acos},
  {0, "atan",  (void (*)()) atan},
  {0, "asinh", (void (*)()) asinh},
  {0, "acosh", (void (*)()) acosh},
  {0, "atanh", (void (*)()) atanh},
  {0, "gamma", (void (*)()) tgamma},

  // cdo functions
  {1, "sum",  (void (*)()) expr_sum},
  // {1, "fldmean",  cdo_fldmean},
  /*
  {1, "min",   min},
  {1, "max",   max},
  {1, "sum",   sum},
  {1, "avg",   avg},
  {1, "mean",  mean},
  {1, "std",   std},
  {1, "var",   var},
  */
};

static int NumFunc = sizeof(fun_sym_tbl) / sizeof(fun_sym_tbl[0]);

static
void param_meta_copy(paramType *out, paramType *in)
{
  out->gridID   = in->gridID;
  out->zaxisID  = in->zaxisID;
  out->steptype = in->steptype;
  out->ngp      = in->ngp;
  out->nlev     = in->nlev;
  out->missval  = in->missval;
  out->nmiss    = 0;
  out->name     = NULL;
}

static
nodeType *expr_con_con(int oper, nodeType *p1, nodeType *p2)
{
  nodeType *p = (nodeType*) Malloc(sizeof(nodeType));

  p->type = typeCon;

  double cval1 = p1->u.con.value;
  double cval2 = p2->u.con.value;

  switch ( oper )
    {
    case '+':  cval1 = cval1 + cval2; break;
    case '-':  cval1 = cval1 - cval2; break;
    case '*':  cval1 = cval1 * cval2; break;
    case '/':  cval1 = cval1 / cval2; break;
    case '^':  cval1 = pow(cval1, cval2); break;
    default:   cdoAbort("%s: operator %c unsupported!", __func__, oper); break;
    }

  p->u.con.value = cval1;

  return p;
}

static
void oper_expr_con_var(int oper, int nmiss, long n, double missval1, double missval2,
                       double *restrict odat, double cval, const double *restrict idat)
{
  long i;

  switch ( oper )
    {
    case '+':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = ADDMN(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] = cval + idat[i];
      break;
    case '-':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = SUBMN(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] = cval - idat[i];
      break;
    case '*':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MULMN(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] = cval * idat[i];
      break;
    case '/':
      for ( i=0; i<n; ++i ) odat[i] = DIVMN(cval, idat[i]);
      break;
    case '^':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = POWMN(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] = pow(cval, idat[i]);
      break;
    case '<':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPLT(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPLT(cval, idat[i]);
      break;
    case '>':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPGT(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPGT(cval, idat[i]);
      break;
    case LE:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPLE(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPLE(cval, idat[i]);
      break;
    case GE:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPGE(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPGE(cval, idat[i]);
      break;
    case NE:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPNE(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPNE(cval, idat[i]);
      break;
    case EQ:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPEQ(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPEQ(cval, idat[i]);
      break;
    case LEG:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPLEG(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPLEG(cval, idat[i]);
      break;
    case AND:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPAND(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPAND(cval, idat[i]);
    case OR:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPOR(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPOR(cval, idat[i]);
      break;
    default:
      cdoAbort("%s: operator %c unsupported!", __func__, oper);
      break;
    }
}

static
void oper_expr_var_con(int oper, int nmiss, long n, double missval1, double missval2,
                       double *restrict odat, const double *restrict idat, double cval)
{
  long i;

  switch ( oper )
    {
    case '+':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = ADDMN(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] = idat[i] + cval;
      break;
    case '-':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = SUBMN(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] = idat[i] - cval;
      break;
    case '*':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MULMN(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] = idat[i] * cval;
      break;
    case '/':
      if ( nmiss || IS_EQUAL(cval, 0) ) for ( i=0; i<n; ++i ) odat[i] = DIVMN(idat[i], cval);
      else                              for ( i=0; i<n; ++i ) odat[i] = idat[i] / cval;
      break;
    case '^':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = POWMN(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] = pow(idat[i], cval);
      break;
    case '<':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPLT(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPLT(idat[i], cval);
      break;
    case '>':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPGT(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPGT(idat[i], cval);
      break;
    case LE:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPLE(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPLE(idat[i], cval);
      break;
    case GE:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPGE(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPGE(idat[i], cval);
      break;
    case NE:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPNE(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPNE(idat[i], cval);
      break;
    case EQ:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPEQ(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPEQ(idat[i], cval);
      break;
    case LEG:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPLEG(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPLEG(idat[i], cval);
      break;
    case AND:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPAND(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPAND(idat[i], cval);
      break;
    case OR:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPOR(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPOR(idat[i], cval);
      break;
    default:
      cdoAbort("%s: operator %c unsupported!", __func__, oper);
      break;
    }
}

static
void oper_expr_var_var(int oper, int nmiss, long ngp, double missval1, double missval2,
                       double *restrict odat, const double *restrict idat1, const double *restrict idat2)
{
  long i;

  switch ( oper )
    {
    case '+':
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = ADDMN(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] = idat1[i] + idat2[i];
      break;
    case '-':
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = SUBMN(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] = idat1[i] - idat2[i];
      break;
    case '*':
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MULMN(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] = idat1[i] * idat2[i];
      break;
    case '/':
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = DIVMN(idat1[i], idat2[i]);
      else
        {
          for ( i = 0; i < ngp; ++i )
            {
              if ( IS_EQUAL(idat2[i], 0.) ) odat[i] = missval1;
              else                          odat[i] = idat1[i] / idat2[i];
            }
        }
      break;
    case '^':
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = POWMN(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] = pow(idat1[i], idat2[i]);
      break;
    case '<':
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPLT(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPLT(idat1[i], idat2[i]);
      break;
    case '>':
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPGT(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPGT(idat1[i], idat2[i]);
      break;
    case LE:
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPLE(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPLE(idat1[i], idat2[i]);
      break;
    case GE:
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPGE(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPGE(idat1[i], idat2[i]);
      break;
    case NE:
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPNE(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPNE(idat1[i], idat2[i]);
      break;
    case EQ:
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPEQ(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPEQ(idat1[i], idat2[i]);
      break;
    case LEG:
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPLEG(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPLEG(idat1[i], idat2[i]);
      break;
    case AND:
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPAND(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPAND(idat1[i], idat2[i]);
      break;
    case OR:
      if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPOR(idat1[i], idat2[i]);
      else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPOR(idat1[i], idat2[i]);
      break;
    default:
      cdoAbort("%s: operator %d (%c) unsupported!", __func__, (int)oper, oper);
      break;
    }
}

static
nodeType *expr_con_var(int oper, nodeType *p1, nodeType *p2)
{
  int gridID   = p2->param.gridID;
  int zaxisID  = p2->param.zaxisID;
  int nmiss    = p2->param.nmiss;
  double missval1 = p2->param.missval;
  double missval2 = p2->param.missval;

  int ngp  = gridInqSize(gridID);
  int nlev = zaxisInqSize(zaxisID);
  long n   = ngp*nlev;

  nodeType *p = (nodeType*) Malloc(sizeof(nodeType));

  p->type     = typeVar;
  p->ltmpvar  = true;
  p->u.var.nm = strdup("tmp");
  param_meta_copy(&p->param, &p2->param);
  p->param.name = p->u.var.nm;

  p->param.data = (double*) Malloc(n*sizeof(double));
  double *restrict odat = p->param.data;
  const double *restrict idat = p2->param.data;
  double cval = p1->u.con.value;

  oper_expr_con_var(oper, nmiss, n, missval1, missval2, odat, cval, idat);

  nmiss = 0;
  for ( long i = 0; i < n; i++ )
    if ( DBL_IS_EQUAL(odat[i], missval1) ) nmiss++;

  p->param.nmiss = nmiss;

  if ( p2->ltmpvar ) Free(p2->param.data);

  return p;
}

static
nodeType *expr_var_con(int oper, nodeType *p1, nodeType *p2)
{
  int gridID   = p1->param.gridID;
  int zaxisID  = p1->param.zaxisID;
  int nmiss    = p1->param.nmiss;
  double missval1 = p1->param.missval;
  double missval2 = p1->param.missval;

  int ngp  = gridInqSize(gridID);
  int nlev = zaxisInqSize(zaxisID);
  long n   = ngp*nlev;

  nodeType *p = (nodeType*) Malloc(sizeof(nodeType));

  p->type     = typeVar;
  p->ltmpvar  = true;
  p->u.var.nm = strdup("tmp");
  param_meta_copy(&p->param, &p1->param);
  p->param.name = p->u.var.nm;

  p->param.data = (double*) Malloc(n*sizeof(double));
  double *restrict odat = p->param.data;
  const double *restrict idat = p1->param.data;
  double cval = p2->u.con.value;

  oper_expr_var_con(oper, nmiss, n, missval1, missval2, odat, idat, cval);

  nmiss = 0;
  for ( long i = 0; i < n; i++ )
    if ( DBL_IS_EQUAL(odat[i], missval1) ) nmiss++;

  p->param.nmiss = nmiss;

  if ( p1->ltmpvar ) Free(p1->param.data);

  return p;
}

static
nodeType *expr_var_var(int oper, nodeType *p1, nodeType *p2)
{
  long i;
  long nlev, k;
  long loff, loff1, loff2;
  int nmiss;

  int nmiss1 = p1->param.nmiss;
  int nmiss2 = p2->param.nmiss;
  double missval1 = p1->param.missval;
  double missval2 = p2->param.missval;

  long ngp1 = gridInqSize(p1->param.gridID);
  long ngp2 = gridInqSize(p2->param.gridID);

  if ( ngp1 != ngp2 && ngp2 != 1 )
    cdoAbort("%s: Number of grid points differ (ngp1 = %ld, ngp2 = %ld)", __func__, ngp1, ngp2);

  long ngp = ngp1;

  long nlev1 = zaxisInqSize(p1->param.zaxisID);
  long nlev2 = zaxisInqSize(p2->param.zaxisID);

  nodeType *p = (nodeType*) Malloc(sizeof(nodeType));

  p->type     = typeVar;
  p->ltmpvar  = true;
  p->u.var.nm = strdup("tmp");

  if ( nlev1 > nlev2 )
    {
      if ( nlev2 != 1 ) cdoAbort("nlev2 = %d must be 1!", nlev2);
      nlev = nlev1;
      param_meta_copy(&p->param, &p1->param);
    }
  else if ( nlev2 > nlev1 )
    {
      if ( nlev1 != 1 ) cdoAbort("nlev1 = %d must be 1!", nlev1);
      nlev = nlev2;
      param_meta_copy(&p->param, &p->param);
    }
  else
    {
      nlev = nlev1;
      param_meta_copy(&p->param, &p1->param);
    }

  p->param.name = p1->param.name;

  p->param.data = (double*) Malloc(ngp*nlev*sizeof(double));

  for ( k = 0; k < nlev; k++ )
    {
      loff = k*ngp;

      if ( nlev1 == 1 ) loff1 = 0;
      else              loff1 = k*ngp1;

      if ( nlev2 == 1 ) loff2 = 0;
      else              loff2 = k*ngp2;

      const double *restrict idat1 = p1->param.data+loff1;
      const double *restrict idat2 = p2->param.data+loff2;
      double *restrict odat = p->param.data+loff;
      int nmiss = nmiss1 > 0 || nmiss2 > 0;

      if ( ngp2 == 1 )
        oper_expr_var_con(oper, nmiss, ngp, missval1, missval2, odat, idat1, idat2[0]);
      else
        oper_expr_var_var(oper, nmiss, ngp, missval1, missval2, odat, idat1, idat2);
    }

  nmiss = 0;
  for ( i = 0; i < ngp*nlev; i++ )
    if ( DBL_IS_EQUAL(p->param.data[i], missval1) ) nmiss++;

  p->param.nmiss = nmiss;

  if ( p1->ltmpvar ) Free(p1->param.data);
  if ( p2->ltmpvar ) Free(p2->param.data);

  return p;
}

static
void ex_copy(nodeType *p2, nodeType *p1)
{
  if ( cdoVerbose ) printf("\tcopy %s\n", p1->u.var.nm);

  int ngp = p1->param.ngp;
  assert(ngp > 0);

  if ( (ngp) != p2->param.ngp )
    cdoAbort("%s: Number of grid points differ (%s[%d] = %s[%d])",
             __func__, p2->param.name, p2->param.ngp, p1->param.name, ngp);

  int nlev = p1->param.nlev;
  assert(nlev > 0);

  if ( nlev != p2->param.nlev )
    cdoAbort("%s: Number of levels differ (nlev1 = %d, nlev2 = %d)",
             __func__, nlev, p2->param.nlev);

  for ( size_t i = 0; i < (size_t)ngp*nlev; ++i ) p2->param.data[i] = p1->param.data[i];

  p2->param.missval = p1->param.missval;
  p2->param.nmiss   = p1->param.nmiss;
}

static
nodeType *expr(int oper, nodeType *p1, nodeType *p2)
{
  nodeType *p = NULL;

  if ( p1->type == typeVar && p2->type == typeVar )
    {
      p = expr_var_var(oper, p1, p2);
      if ( cdoVerbose )
	printf("\t%s %c %s\n", p1->u.var.nm, oper, p2->u.var.nm);
    }
  else if ( p1->type == typeCon && p2->type == typeCon )
    {
      p = expr_con_con(oper, p1, p2);
      if ( cdoVerbose )
	printf("\t%g %c %g\n", p1->u.con.value, oper, p2->u.con.value);
    }
  else if ( p1->type == typeVar && p2->type == typeCon )
    {
      p = expr_var_con(oper, p1, p2);
      if ( cdoVerbose )
	printf("\t%s %c %g\n", p1->u.var.nm, oper, p2->u.con.value);
    }
  else if ( p1->type == typeCon && p2->type == typeVar )
    {
      p = expr_con_var(oper, p1, p2);
      if ( cdoVerbose )
	printf("\t%g %c %s\n", p1->u.con.value, oper, p2->u.var.nm);
    }
  else
    cdoAbort("Internal problem!");

  return p;
}

static
nodeType *ex_fun_con(char *fun, nodeType *p1)
{
  int funcID = -1;

  nodeType *p = (nodeType*) Malloc(sizeof(nodeType));

  p->type = typeCon;

  for ( int i = 0; i < NumFunc; i++)
    if ( fun_sym_tbl[i].type == 0 )
      if ( strcmp(fun, fun_sym_tbl[i].name) == 0 )
	{ 
	  funcID = i;
	  break;
	}

  if ( funcID == -1 )
    cdoAbort("Function >%s< not available!", fun);

  double (*exprfunc)(double) = (double (*)(double)) fun_sym_tbl[funcID].func;
  p->u.con.value = exprfunc(p1->u.con.value);

  return p;
}

static
nodeType *ex_fun_var(char *fun, nodeType *p1)
{
  int funcID = -1;
  for ( int i = 0; i < NumFunc; i++ )
    if ( strcmp(fun, fun_sym_tbl[i].name) == 0 )
      { 
	funcID = i;
	break;
      }

  if ( funcID == -1 )
    cdoAbort("Function >%s< not available!", fun);

  int functype = fun_sym_tbl[funcID].type;

  int gridID  = p1->param.gridID;
  int zaxisID = p1->param.zaxisID;
  int nmiss   = p1->param.nmiss;
  double missval = p1->param.missval;

  long ngp  = gridInqSize(gridID);
  long nlev = zaxisInqSize(zaxisID);

  nodeType *p = (nodeType*) Malloc(sizeof(nodeType));

  p->type     = typeVar;
  p->ltmpvar  = true;
  p->u.var.nm = strdup("tmp");
  p->param.zaxisID  = zaxisID;
  p->param.missval  = missval;

  if ( functype == 0 )
    {
      p->param.gridID = gridID;
    }
  else
    {
      ngp = 1;
      int sgridID = gridCreate(GRID_GENERIC, ngp);
      p->param.gridID = sgridID;
    }

  p->param.data = (double*) Malloc(ngp*nlev*sizeof(double));
  double *restrict pdata = p->param.data;
  double *restrict p1data = p1->param.data;
  
  if ( nmiss > 0 )
    {
      double (*exprfunc)(double) = (double (*)(double)) fun_sym_tbl[funcID].func;
      for ( long i = 0; i < ngp*nlev; i++ )
	{
	  errno = -1;
	  pdata[i] = DBL_IS_EQUAL(p1data[i], missval) ? missval : exprfunc(p1data[i]);
	  if ( errno == EDOM || errno == ERANGE ) pdata[i] = missval;
	  else if ( isnan(pdata[i]) ) pdata[i] = missval;
	}
    }
  else
    {
      if ( fun_sym_tbl[funcID].type == 0 )
        {
          double (*exprfunc)(double) = (double (*)(double)) fun_sym_tbl[funcID].func;
          for ( long i = 0; i < ngp*nlev; i++ )
            {
              errno = -1;
              pdata[i] = exprfunc(p1data[i]);
              if ( errno == EDOM || errno == ERANGE ) pdata[i] = missval;
              else if ( isnan(pdata[i]) ) pdata[i] = missval;
            }
        }
      else
        {
          double (*exprfunc)(int,double*) = (double (*)(int,double*)) fun_sym_tbl[funcID].func;
          for ( int k = 0; k < nlev; k++ )
            pdata[k] = exprfunc(ngp, p1data+k*ngp);
        }
    }

  nmiss = 0;
  for ( long i = 0; i < ngp*nlev; i++ )
    if ( DBL_IS_EQUAL(pdata[i], missval) ) nmiss++;

  p->param.nmiss = nmiss;

  if ( p1->ltmpvar ) Free(p1data);

  return p;
}

static
nodeType *ex_fun(char *fun, nodeType *p1)
{
  nodeType *p = NULL;

  if ( p1->type == typeVar )
    {
      p = ex_fun_var(fun, p1);
      if ( cdoVerbose ) printf("\t%s (%s)\n", fun, p1->u.var.nm);
    }
  else if ( p1->type == typeCon )
    {
      p = ex_fun_con(fun, p1);
      if ( cdoVerbose ) printf("\t%s (%g)\n", fun, p1->u.con.value);
    }
  else
    cdoAbort("Internal problem!");

  return p;
}

static
nodeType *ex_uminus_var(nodeType *p1)
{
  int gridID  = p1->param.gridID;
  int zaxisID = p1->param.zaxisID;
  int nmiss   = p1->param.nmiss;
  double missval = p1->param.missval;

  long ngp  = gridInqSize(gridID);
  long nlev = zaxisInqSize(zaxisID);

  nodeType *p = (nodeType*) Malloc(sizeof(nodeType));

  p->type     = typeVar;
  p->ltmpvar  = true;
  p->u.var.nm = strdup("tmp");
  param_meta_copy(&p->param, &p1->param);
  p->param.name = p->u.var.nm;

  p->param.data = (double*) Malloc(ngp*nlev*sizeof(double));
  double *restrict pdata = p->param.data;
  const double *restrict p1data = p1->param.data;

  if ( nmiss > 0 )
    {
      for ( long i = 0; i < ngp*nlev; i++ )
	pdata[i] = DBL_IS_EQUAL(p1data[i], missval) ? missval : -(p1data[i]);
    }
  else
    {
      for ( long i = 0; i < ngp*nlev; i++ )
	pdata[i] = -(p1data[i]);
    }

  p->param.nmiss = nmiss;
  
  return p;
}

static
nodeType *ex_uminus_con(nodeType *p1)
{
  nodeType *p = (nodeType*) Malloc(sizeof(nodeType));

  p->type = typeCon;

  p->u.con.value = -(p1->u.con.value);

  return p;
}

static
nodeType *ex_uminus(nodeType *p1)
{
  nodeType *p = NULL;

  if ( p1->type == typeVar )
    {
      p = ex_uminus_var(p1);
      if ( cdoVerbose ) printf("\t- (%s)\n", p1->u.var.nm);
    }
  else if ( p1->type == typeCon )
    {
      p = ex_uminus_con(p1);
      if ( cdoVerbose ) printf("\t- (%g)\n", p1->u.con.value);
    }
  else
    cdoAbort("Internal problem!");

  return p;
}

static
nodeType *ex_ifelse(nodeType *p1, nodeType *p2, nodeType *p3)
{
  if ( cdoVerbose ) printf("\t %s ? %s : %s\n", p1->u.var.nm, p2->u.var.nm, p3->u.var.nm);

  if ( p1->type == typeCon ) cdoAbort("expr?expr:expr: First expression is a constant but must be a variable!");

  int nmiss1 = p1->param.nmiss;
  long ngp1 = gridInqSize(p1->param.gridID);
  long nlev1 = zaxisInqSize(p1->param.zaxisID);
  double missval1 = p1->param.missval;
  double *pdata1 = p1->param.data;

  long ngp = ngp1;
  long nlev = nlev1;
  nodeType *px = p1;

  double missval2 = missval1;
  double *pdata2;
  long ngp2 = 1;
  long nlev2 = 1;
  
  if ( p2->type == typeCon )
    {
      pdata2 = &p2->u.con.value;
    }
  else
    {
      ngp2 = gridInqSize(p2->param.gridID);
      nlev2 = zaxisInqSize(p2->param.zaxisID);
      missval2 = p2->param.missval;
      pdata2 = p2->param.data;
      if ( ngp2 > 1 && ngp2 != ngp1 )
	cdoAbort("expr?expr:expr: Number of grid points differ (ngp1 = %ld, ngp2 = %ld)", ngp1, ngp2);
      if ( nlev2 > 1 && nlev2 != nlev )
	{
	  if ( nlev == 1 )
	    {
	      nlev = nlev2;
	      px = p2;
	    }
	  else
	    cdoAbort("expr?expr:expr: Number of levels differ (nlev = %ld, nlev2 = %ld)", nlev, nlev2);
	}
    }

  double missval3 = missval1;
  double *pdata3;
  long ngp3 = 1;
  long nlev3 = 1;
  
  if ( p3->type == typeCon )
    {
      pdata3 = &p3->u.con.value;
    }
  else
    {
      ngp3 = gridInqSize(p3->param.gridID);
      nlev3 = zaxisInqSize(p3->param.zaxisID);
      missval3 = p3->param.missval;
      pdata3 = p3->param.data;
      if ( ngp3 > 1 && ngp3 != ngp1 )
	cdoAbort("expr?expr:expr: Number of grid points differ (ngp1 = %ld, ngp3 = %ld)", ngp1, ngp3);
      if ( nlev3 > 1 && nlev3 != nlev )
	{
	  if ( nlev == 1 )
	    {
	      nlev = nlev3;
	      px = p3;
	    }
	  else
	    cdoAbort("expr?expr:expr: Number of levels differ (nlev = %ld, nlev3 = %ld)", nlev, nlev3);
	}
    }

  nodeType *p = (nodeType*) Malloc(sizeof(nodeType));

  p->type     = typeVar;
  p->ltmpvar  = true;
  p->u.var.nm = strdup("tmp");

  param_meta_copy(&p->param, &px->param);
  p->param.name = p->u.var.nm;

  p->param.data = (double*) Malloc(ngp*nlev*sizeof(double));

  long loff, loff1, loff2, loff3;

  for ( long k = 0; k < nlev; ++k )
    {
      loff = k*ngp;

      if ( nlev1 == 1 ) loff1 = 0;
      else              loff1 = k*ngp;

      if ( nlev2 == 1 ) loff2 = 0;
      else              loff2 = k*ngp;

      if ( nlev3 == 1 ) loff3 = 0;
      else              loff3 = k*ngp;

      const double *restrict idat1 = pdata1+loff1;
      const double *restrict idat2 = pdata2+loff2;
      const double *restrict idat3 = pdata3+loff3;
      double *restrict odat = p->param.data+loff;

      double ival2 = idat2[0];
      double ival3 = idat3[0];
      for ( long i = 0; i < ngp; ++i ) 
	{
	  if ( ngp2 > 1 ) ival2 = idat2[i];
	  if ( ngp3 > 1 ) ival3 = idat3[i];

	  if ( nmiss1 && DBL_IS_EQUAL(idat1[i], missval1) )
	    odat[i] = missval1;
	  else if ( IS_NOT_EQUAL(idat1[i], 0) )
	    odat[i] = DBL_IS_EQUAL(ival2, missval2) ? missval1 : ival2;
	  else
	    odat[i] = DBL_IS_EQUAL(ival3, missval3) ? missval1 : ival3;
	}
    }

  return p;
}
/*
static
int exNode(nodeType *p, parse_param_t *parse_arg)
{
  if ( ! p ) return 0;

  // node is leaf
  if ( p->type == typeCon || p->type == typeVar || p->u.opr.nops == 0 )
    {
      return 0;
    }

  // node has children
  for ( int k = 0; k < p->u.opr.nops; k++ )
    {
      exNode(p->u.opr.op[k], parse_arg);
    }

  return 0;
}
*/

nodeType *expr_run(nodeType *p, parse_param_t *parse_arg)
{
  int gridID1 = -1, zaxisID1 = -1, tsteptype1 = -1;
  double missval = 0;
  int varID;
  nodeType *rnode = NULL;

  if ( ! p ) return rnode;

  /*  if ( ! parse_arg->init ) { exNode(p, parse_arg); return 0; } */

  switch ( p->type )
    {
    case typeCon:       
      if ( parse_arg->init )
	{
	  if ( parse_arg->debug )
	    printf("\tpush const \t%g\n", p->u.con.value);
	}
      else
	{
	  rnode = p;
	}

      break;
    case typeVar:
      /*    if ( parse_arg->init ) */
	{
	  if ( parse_arg->debug )
	    printf("\tpush var \t%s\n", p->u.var.nm);

          int nvars = parse_arg->nparams;
	  for ( varID = nvars-1; varID >= 0; --varID )
	    {
	      if ( strcmp(parse_arg->params[varID].name, p->u.var.nm) == 0 ) break;
	    }
          if ( varID == -1 )
	    {
	      cdoAbort("Variable >%s< not found!", p->u.var.nm);
	    }
	  else
	    {
	      int nlev1, nlev2 = 0;

	      if ( varID < parse_arg->nvars1 && parse_arg->needed[varID] == false )
		{
		  parse_arg->needed[varID] = true;
		}

              gridID1    = parse_arg->params[varID].gridID;
              zaxisID1   = parse_arg->params[varID].zaxisID;
              tsteptype1 = parse_arg->params[varID].steptype;
              missval    = parse_arg->params[varID].missval;
              nlev1      = parse_arg->params[varID].nlev;
              
	      parse_arg->missval2 = missval;

	      if ( parse_arg->gridID2 == -1 )
		parse_arg->gridID2 = gridID1;

              // printf(">> typeVar: zaxisID2 %d %s\n", parse_arg->zaxisID2, parse_arg->params[varID].name);
	      if ( parse_arg->zaxisID2 != -1 ) nlev2 = zaxisInqSize(parse_arg->zaxisID2);

	      if ( parse_arg->zaxisID2 == -1 || (nlev1 > 1 && nlev2 == 1) )
		parse_arg->zaxisID2 = zaxisID1;

	      if ( parse_arg->tsteptype2 == -1 || parse_arg->tsteptype2 == TSTEP_CONSTANT )
		parse_arg->tsteptype2 = tsteptype1;
	    }
	}
	/* else */
	{ 
          param_meta_copy(&p->param, &parse_arg->params[varID]);
          p->param.name = parse_arg->params[varID].name;

 	  if ( parse_arg->debug )
	    printf("var: u.var.nm=%s name=%s gridID=%d zaxisID=%d ngp=%d nlev=%d  varID=%d\n",
                   p->u.var.nm, p->param.name, p->param.gridID, p->param.zaxisID, p->param.ngp, p->param.nlev, varID);
         
	  p->ltmpvar = false;
	  if ( ! parse_arg->init )
	    {
              p->param.data  = parse_arg->params[varID].data;
              p->param.nmiss = parse_arg->params[varID].nmiss;
            }
	  rnode = p;
	}

      break;
    case typeFun:
      if ( parse_arg->init )
	{
	  expr_run(p->u.fun.op, parse_arg);

	  if ( parse_arg->debug ) printf("\tcall \t%s\n", p->u.fun.name); 
	}
      else
	{
	  rnode = ex_fun(p->u.fun.name, expr_run(p->u.fun.op, parse_arg));
	}
      break;
    case typeOpr:
      switch( p->u.opr.oper )
	{
        case '=':
	  parse_arg->gridID2    = -1;
	  parse_arg->zaxisID2   = -1;
          parse_arg->tsteptype2 = -1;

          // printf(">> = init1: zaxisID2 %d %s\n", parse_arg->zaxisID2, parse_arg->params[varID].name);
          rnode = expr_run(p->u.opr.op[1], parse_arg);
          // printf(">> = init2: zaxisID2 %d %s\n", parse_arg->zaxisID2, parse_arg->params[varID].name);

	  if ( parse_arg->init )
	    {
	      if ( parse_arg->debug )
		printf("\tpop  var \t%s\n", p->u.opr.op[0]->u.var.nm);
	      /*
	      if ( p->u.opr.op[1]->type != typeVar )
		cdoAbort("Operand not variable!");
	      */
	      if ( parse_arg->gridID2 == -1 || parse_arg->zaxisID2 == -1 || parse_arg->tsteptype2 == -1 )
		cdoAbort("Operand not variable!");

	      const char *varname2 = p->u.opr.op[0]->u.var.nm;
              int nvars = parse_arg->nparams;
              for ( varID = nvars-1; varID >= 0; --varID )
                {
                  if ( strcmp(parse_arg->params[varID].name, varname2) == 0 ) break;
                }
              if ( varID >= 0 )
                {
                  if ( varID < parse_arg->nvars1 ) parse_arg->params[varID].select = true;
                  else cdoWarning("Variable %s already defined!", varname2);
                }
              else
                {
                  varID = parse_arg->nparams;
                  if ( varID >= parse_arg->maxparams )
                    cdoAbort("Too many parameter (limit=%d)", parse_arg->maxparams);

                  parse_arg->params[varID].name     = strdup(varname2);
                  parse_arg->params[varID].missval  = parse_arg->missval2;
                  parse_arg->params[varID].gridID   = parse_arg->gridID2;
                  parse_arg->params[varID].zaxisID  = parse_arg->zaxisID2;
                  parse_arg->params[varID].steptype = parse_arg->tsteptype2;
                  parse_arg->params[varID].ngp      = gridInqSize(parse_arg->gridID2);
                  parse_arg->params[varID].nlev     = zaxisInqSize(parse_arg->zaxisID2);
                  parse_arg->nparams++;
                }
	    }
	  else
	    {
	      if ( parse_arg->debug )
		printf("\tpop  var\t%s\t%s\n", p->u.opr.op[0]->u.var.nm, rnode->u.var.nm);

              int nvars = parse_arg->nparams;
              for ( varID = nvars-1; varID >= 0; --varID )
                {
                  if ( strcmp(parse_arg->params[varID].name, p->u.opr.op[0]->u.var.nm) == 0 ) break;
                }
              
	      if ( varID < 0 )
		{
		  cdoAbort("Variable >%s< not found!", p->u.opr.op[0]->u.var.nm);
		}
	      else
		{
                  parse_arg->gridID2  = parse_arg->params[varID].gridID;
                  parse_arg->zaxisID2 = parse_arg->params[varID].zaxisID;
                  // printf(">> = run: zaxisID2 %d %s\n", parse_arg->zaxisID2, parse_arg->params[varID].name);
                  parse_arg->tsteptype2 = parse_arg->params[varID].steptype;

                  param_meta_copy(&p->param, &parse_arg->params[varID]);
                  p->param.name = parse_arg->params[varID].name;
                  p->param.data = parse_arg->params[varID].data;
		  p->ltmpvar    = false;

                  //printf(">>copy %s\n", p->param.name);
                  //printf(">>copy %s\n", rnode->param.name);
		  ex_copy(p, rnode);

		  if ( rnode->ltmpvar )
                    {
                      Free(rnode->param.data);
                    }
                }
	    }

	  break;
        case UMINUS:    
	  if ( parse_arg->init )
	    {
	      expr_run(p->u.opr.op[0], parse_arg);

	      if ( parse_arg->debug ) printf("\tneg\n");
	    }
	  else
	    {
	      rnode = ex_uminus(expr_run(p->u.opr.op[0], parse_arg));
	    }

	  break;
        case '?':    
	  if ( parse_arg->init )
	    {
	      expr_run(p->u.opr.op[0], parse_arg);
	      expr_run(p->u.opr.op[1], parse_arg);
	      expr_run(p->u.opr.op[2], parse_arg);

	      if ( parse_arg->debug ) printf("\t?:\n");
	    }
	  else
	    {
	      rnode = ex_ifelse(expr_run(p->u.opr.op[0], parse_arg),
			        expr_run(p->u.opr.op[1], parse_arg),
			        expr_run(p->u.opr.op[2], parse_arg));
	    }

	  break;
        default:
	  if ( parse_arg->init )
	    {
	      expr_run(p->u.opr.op[0], parse_arg);
	      expr_run(p->u.opr.op[1], parse_arg);
	      if ( parse_arg->debug )
		switch( p->u.opr.oper )
		  {
		  case '+':  printf("\tadd\n"); break;
		  case '-':  printf("\tsub\n"); break;
		  case '*':  printf("\tmul\n"); break;
		  case '/':  printf("\tdiv\n"); break;
		  case '<':  printf("\tcompLT\n"); break;
		  case '>':  printf("\tcompGT\n"); break;
		  case LE:   printf("\tcompLE\n"); break;
		  case GE:   printf("\tcompGE\n"); break;
		  case NE:   printf("\tcompNE\n"); break;
		  case EQ:   printf("\tcompEQ\n"); break;
		  case LEG:  printf("\tcompLEG\n"); break;
		  case AND:  printf("\tcompAND\n"); break;
		  case OR:   printf("\tcompOR\n"); break;
		  }
	    }
	  else
	    {
	      rnode = expr(p->u.opr.oper, expr_run(p->u.opr.op[0], parse_arg),
			                  expr_run(p->u.opr.op[1], parse_arg));
	    }
          break;
        }
      break;
    }

  return rnode;
}

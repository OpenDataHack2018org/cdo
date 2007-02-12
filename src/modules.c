/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2007 Uwe Schulzweida, schulzweida@dkrz.de
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
#include <string.h>
#include <stdlib.h>
#include "cdo.h"
#include "operator_help.h"
#include "modules.h"
#include "error.h"


#define  MAX_MOD_OPERATORS  99        /* maximum number of operators for a module */

typedef struct {
  void *(*func)(void *);              /* Module                   */
  char **help;                        /* Help                     */
  char *operators[MAX_MOD_OPERATORS]; /* Operator names           */
  int  streamInCnt;                   /* number of input streams  */
  int  streamOutCnt;                  /* number of output streams */
}
MODULES;


void *Arith(void *argument);
void *Arithc(void *argument);
void *Arithdays(void *argument);
void *Arithlat(void *argument);
void *Cat(void *argument);
void *Change(void *argument);
void *Comp(void *argument);
void *Compc(void *argument);
void *Cond(void *argument);
void *Cond2(void *argument);
void *Condc(void *argument);
void *Copy(void *argument);
void *Detrend(void *argument);
void *Diff(void *argument);
void *Echam5ini(void *argument);
void *Enlarge(void *argument);
void *Enlargegrid(void *argument);
void *Ensstat(void *argument);
void *Expr(void *argument);
void *Filedes(void *argument);
void *Fillmiss(void *argument);
void *Fldrms(void *argument);
void *Fldstat(void *argument);
void *Gradsdes(void *argument);
void *Histogram(void *argument);
void *Info(void *argument);
void *Input(void *argument);
void *Intgrid(void *argument);
void *Intgridtraj(void *argument);
void *Inttime(void *argument);
void *Intntime(void *argument);
void *Intyear(void *argument);
void *Invert(void *argument);
void *Log(void *argument);
void *Maskbox(void *argument);
void *Mastrfu(void *argument);
void *Math(void *argument);
void *Merge(void *argument);
void *Mergegrid(void *argument);
void *Mergetime(void *argument);
void *Merstat(void *argument);
void *Mrotuv(void *argument);
void *Ninfo(void *argument);
void *Nmltest(void *argument);
void *Output(void *argument);
void *Outputgmt(void *argument);
void *Pinfo(void *argument);
void *Remap(void *argument);
void *Replace(void *argument);
void *Rotuv(void *argument);
/* RQ */
void *Runpctl(void *argument);
/* QR */
void *Runstat(void *argument);
/* RQ */
void *Seaspctl(void *argument);
/* QR */
void *Seasstat(void *argument);
void *Selbox(void *argument);
void *Select(void *argument);
void *Seloperator(void *argument);
void *Selrec(void *argument);
/* RQ */
void *Selpctl(void *argument);
/* QR */
void *Selstat(void *argument);
void *Seltime(void *argument);
void *Set(void *argument);
void *Setbox(void *argument);
void *Setgatt(void *argument);
void *Setgrid(void *argument);
void *Sethalo(void *argument);
void *Setmiss(void *argument);
void *Setrcaname(void *argument);
void *Settime(void *argument);
void *Setzaxis(void *argument);
void *Showinfo(void *argument);
void *Sinfo(void *argument);
void *Sort(void *argument);
void *Specinfo(void *argument);
void *Spectral(void *argument);
void *Split(void *argument);
void *Splitrec(void *argument);
void *Splittime(void *argument);
void *Splityear(void *argument);
void *Subtrend(void *argument);
void *Template1(void *argument);
void *Template2(void *argument);
void *Test(void *argument);
void *Test2(void *argument);
void *Timsort(void *argument);
/* RQ */
void *Timpctl(void *argument);
/* QR */
void *Timstat(void *argument);
void *Trend(void *argument);
void *Trms(void *argument);
void *Vardup(void *argument);
void *Vargen(void *argument);
void *Varrms(void *argument);
void *Vertint(void *argument);
void *Vertstat(void *argument);
void *Wind(void *argument);
void *Writegrid(void *argument);
void *Writerandom(void *argument);
/* RQ */
void *Ydaypctl(void *argument);
/* QR */
void *Ydaystat(void *argument);
/* RQ */
void *Ydrunpctl(void *argument);
void *Ydrunstat(void *argument);
/* QR */
void *Ymonarith(void *argument);
/* RQ */
void *Ymonpctl(void *argument);
/* QR */
void *Ymonstat(void *argument);
/* RQ */
void *Yseaspctl(void *argument);
/* QR */
void *Yseasstat(void *argument);
void *Zonstat(void *argument);
/* RQ */
void *EcaCfd(void *argument);
void *EcaCsu(void *argument);
void *EcaCwdi(void *argument);
void *EcaCwfi(void *argument);
void *EcaEtr(void *argument);
void *EcaFd(void *argument);
void *EcaGsl(void *argument);
void *EcaHd(void *argument);
void *EcaHwdi(void *argument);
void *EcaHwfi(void *argument);
void *EcaId(void *argument);
void *EcaSu(void *argument);
void *EcaTr(void *argument);
void *EcaTg10p(void *argument);
void *EcaTg90p(void *argument);
void *EcaTn10p(void *argument);
void *EcaTn90p(void *argument);
void *EcaTx10p(void *argument);
void *EcaTx90p(void *argument);

void *EcaCdd(void *argument);
void *EcaCwd(void *argument);
void *EcaRr1(void *argument);
void *EcaR10mm(void *argument);
void *EcaR20mm(void *argument);
void *EcaR75p(void *argument);
void *EcaR75ptot(void *argument);
void *EcaR90p(void *argument);
void *EcaR90ptot(void *argument);
void *EcaR95p(void *argument);
void *EcaR95ptot(void *argument);
void *EcaR99p(void *argument);
void *EcaR99ptot(void *argument);
void *EcaRx1day(void *argument);
void *EcaRx5day(void *argument);
void *EcaSdii(void *argument);

void *EcaFdns(void *argument);
void *EcaStrwin(void *argument);
void *EcaStrbre(void *argument);
void *EcaStrgal(void *argument);
void *EcaHurr(void *argument);

void *Hi(void *argument);
void *Wct(void *argument);
/* QR */


#define  ArithOperators         {"add",  "sub",  "mul",  "div", "min", "max", "atan2"}
#define  ArithcOperators        {"addc", "subc", "mulc", "divc"}
#define  ArithdaysOperators     {"muldpm", "divdpm", "muldpy", "divdpy"}
#define  ArithlatOperators      {"mulcoslat", "divcoslat"}
#define  CatOperators           {"cat"}
#define  ChangeOperators        {"chcode", "chvar", "chlevel", "chlevelc", "chlevelv", "chltype"}
#define  CompOperators          {"eq",  "ne",  "le",  "lt",  "ge",  "gt"}
#define  CompcOperators         {"eqc", "nec", "lec", "ltc", "gec", "gtc"}
#define  CondOperators          {"ifthen",  "ifnotthen"}
#define  Cond2Operators         {"ifthenelse"}
#define  CondcOperators         {"ifthenc", "ifnotthenc"}
#define  CopyOperators          {"copy", "selall"}
#define  DetrendOperators       {"detrend"}
#define  DiffOperators          {"diff", "diffv"}
#define  Echam5iniOperators     {"read_e5ini", "write_e5ini"}
#define  EnlargeOperators       {"enlarge"}
#define  EnlargegridOperators   {"enlargegrid"}
#define  EnsstatOperators       {"ensmin", "ensmax", "enssum", "ensmean", "ensavg", "ensvar", "ensstd", "enspctl"}
#define  ExprOperators          {"expr", "exprf"}
#define  FiledesOperators       {"filedes", "griddes", "griddes2", "zaxisdes", "vct", "vardes", "taxisdes", "vlist", "partab"}
#define  FillmissOperators      {"fillmiss"}
#define  FldrmsOperators        {"fldrms"}
#define  FldstatOperators       {"fldmin", "fldmax", "fldsum", "fldmean", "fldavg", "fldvar", "fldstd", "fldpctl"}
#define  GradsdesOperators      {"gradsdes1", "gradsdes2", "dumpmap"}
#define  HistogramOperators     {"histcount", "histsum", "histmean"}
#define  InfoOperators          {"info", "infov", "map"}
#define  InputOperators         {"input", "inputsrv", "inputext"}
#define  IntgridOperators       {"intgridbil", "intpoint", "interpolate", "intarea"}
#define  IntgridtrajOperators   {"intgridtraj"}
#define  InttimeOperators       {"inttime"}
#define  IntntimeOperators      {"intntime"}
#define  IntyearOperators       {"intyear"}
#define  InvertOperators        {"invertlat", "invertlon", "invertlatdes", "invertlondes", \
                                 "invertlatdata", "invertlondata"}
#define  LogOperators           {"dumplogs", "daylogs", "monlogs", "dumplogo"}
#define  MaskboxOperators       {"masklonlatbox", "maskindexbox"}
#define  MastrfuOperators       {"mastrfu"}
#define  MathOperators          {"abs", "int", "nint", "sqr", "sqrt", "exp", "ln", "log10", "sin", "cos", "tan", "asin", "acos", "atan"}
#define  MergeOperators         {"merge"}
#define  MergegridOperators     {"mergegrid"}
#define  MergetimeOperators     {"mergetime"}
#define  MerstatOperators       {"mermin", "mermax", "mersum", "mermean", "meravg", "mervar", "merstd", "merpctl"}
#define  MrotuvOperators        {"mrotuvb"}
#define  NinfoOperators         {"nyear", "nmon", "ndate", "ntime", "ncode", "nvar", "nlevel"}
#define  NmltestOperators       {"nmltest"}
#define  OutputOperators        {"output", "outputint", "outputsrv", "outputext", "outputf", "outputts", "outputfld", "outputarr"}
#define  OutputgmtOperators     {"outputcenter", "outputcentercpt", "outputbounds", "outputboundscpt", "outputvector"}
#define  PinfoOperators         {"pinfo", "pinfov"}
#define  RemapOperators         {"remap"}
#define    RemapgridOperators   {"remapcon", "remapbil", "remapbic", "remapdis", "remapdis1", "remapcon1"}
#define    GenweightsOperators  {"gencon", "genbil", "genbic", "gendis"}
#define  ReplaceOperators       {"replace"}
#define  RotuvOperators         {"rotuvb"}
/* RQ */
#define  RunpctlOperators       {"runpctl"}
/* QR */
#define  RunstatOperators       {"runmin",  "runmax",  "runsum",  "runmean",  "runavg",  "runvar",  "runstd"}
/* RQ */
#define  SeaspctlOperators      {"seaspctl"}
/* QR */
#define  SeasstatOperators      {"seasmin",  "seasmax",  "seassum",  "seasmean",  "seasavg",  "seasvar",  "seasstd"}
#define  SelboxOperators        {"sellonlatbox", "selindexbox"}
#define  SelectOperators        {"selcode", "selvar", "selstdname", "sellevel", "selgrid", "selgridname", \
                                 "selzaxis", "selzaxisname", "seltabnum", "delcode", "delvar", "selltype"}
#define  SeloperatorOperators   {"seloperator"}
/* RQ */
#define  SelpctlOperators       {"selpctl"}
/* QR */
#define  SelstatOperators       {"selmin",  "selmax",  "selsum",  "selmean",  "selavg",  "selvar",  "selstd"}
#define  SelrecOperators        {"selrec"}
#define  SeltimeOperators       {"seltimestep", "selyear", "selseas", "selmon", "selday", "selhour", "seldate", "seltime", "selsmon"}
#define  SetOperators           {"setpartab", "setpartabv", "setcode", "setvar", "setlevel", "setltype"}
#define  SetboxOperators        {"setclonlatbox", "setcindexbox"}
#define  SetgattOperators       {"setgatt", "setgatts"}
#define  SetgridOperators       {"setgrid", "setgridtype", "setgridarea"}
#define  SethaloOperators       {"sethalo"}
#define  SetmissOperators       {"setmissval", "setctomiss", "setmisstoc", "setrtomiss"}
#define  SetrcanameOperators    {"setrcaname"}
#define  SettimeOperators       {"setyear", "setmon", "setday", "setdate", "settime", "settunits", \
                                 "settaxis", "setreftime", "setcalendar", "shifttime"}
#define  SetzaxisOperators      {"setzaxis"}
#define  ShowinfoOperators      {"showyear", "showmon", "showdate", "showtime", "showcode", "showvar", \
                                 "showstdname", "showlevel", "showltype", "showformat"}
#define  SinfoOperators         {"sinfo", "sinfov", "sinfop"}
#define  SortOperators          {"sortcode", "sortvar", "sortlevel"}
#define  SpecinfoOperators      {"specinfo"}
#define  SpectralOperators      {"gp2sp", "gp2spl", "sp2gp", "sp2gpl", "sp2sp", "spcut"}
#define  SplitOperators         {"splitcode", "splitvar", "splitlevel", "splitgrid", "splitzaxis"}
#define  SplitrecOperators      {"splitrec"}
#define  SplittimeOperators     {"splithour", "splitday", "splitmon", "splitseas"}
#define  SplityearOperators     {"splityear"}
#define  SubtrendOperators      {"subtrend"}
#define  Template1Operators     {"template1"}
#define  Template2Operators     {"template2"}
#define  TestOperators          {"test"}
#define  Test2Operators         {"test2"}
#define  TimsortOperators       {"timsort"}
/* RQ */
#define  TimpctlOperators       {"timpctl"}
#define    YearpctlOperators    {"yearpctl"}
#define    MonpctlOperators     {"monpctl"}
#define    DaypctlOperators     {"daypctl"}
#define    HourpctlOperators    {"hourpctl"}
/* QR */
#define  TimstatOperators       {"timmin",  "timmax",  "timsum",  "timmean",  "timavg",  "timvar",  "timstd"}
#define    YearstatOperators    {"yearmin", "yearmax", "yearsum", "yearmean", "yearavg", "yearvar", "yearstd"}
#define    MonstatOperators     {"monmin",  "monmax",  "monsum",  "monmean",  "monavg",  "monvar",  "monstd"}
#define    DaystatOperators     {"daymin",  "daymax",  "daysum",  "daymean",  "dayavg",  "dayvar",  "daystd"}
#define    HourstatOperators    {"hourmin", "hourmax", "hoursum", "hourmean", "houravg", "hourvar", "hourstd"}
#define  TrendOperators         {"trend"}
#define  TrmsOperators          {"trms"}
#define  VardupOperators        {"vardup", "varmul"}
#define  VargenOperators        {"random", "const", "topo"}
#define  VarrmsOperators        {"varrms"}
#define  VertintOperators       {"ml2pl", "ml2hl"}
#define  VertstatOperators      {"vertmin", "vertmax", "vertsum", "vertmean", "vertavg", "vertvar", "vertstd"}
#define  WindOperators          {"uv2dv", "uv2dvl", "dv2uv", "dv2uvl", "dv2ps"}
#define  WritegridOperators     {"writegrid", "gridarea"}
#define  WriterandomOperators   {"writerandom"}
/* RQ */
#define  YdaypctlOperators      {"ydaypctl"}
/* QR */
#define  YdaystatOperators      {"ydaymin", "ydaymax", "ydaysum", "ydaymean", "ydayavg", "ydayvar", "ydaystd"}
/* RQ */
#define  YdrunpctlOperators     {"ydrunpctl"}
#define  YdrunstatOperators     {"ydrunmin", "ydrunmax", "ydrunsum", "ydrunmean", "ydrunavg", "ydrunvar", "ydrunstd"}
/* QR */
#define  YmonarithOperators     {"ymonadd", "ymonsub", "ymonmul", "ymondiv"}
/* RQ */
#define  YmonpctlOperators      {"ymonpctl"}
/* QR */
#define  YmonstatOperators      {"ymonmin", "ymonmax", "ymonsum", "ymonmean", "ymonavg", "ymonvar", "ymonstd"}
/* RQ */
#define  YseaspctlOperators     {"yseaspctl"}
/* QR */
#define  YseasstatOperators     {"yseasmin", "yseasmax", "yseassum", "yseasmean", "yseasavg", "yseasvar", "yseasstd"}
#define  ZonstatOperators       {"zonmin", "zonmax", "zonsum", "zonmean", "zonavg", "zonstd", "zonpctl"}

/* RQ */
#define  EcaCfdOperators        {"eca_cfd"}
#define  EcaCsuOperators        {"eca_csu"}
#define  EcaCwfiOperators       {"eca_cwfi"}
#define  EcaHwdiOperators       {"eca_hwdi"}
#define  EcaEtrOperators        {"eca_etr"}
#define  EcaFdOperators         {"eca_fd"}
#define  EcaGslOperators        {"eca_gsl"}
#define  EcaHdOperators         {"eca_hd"}
#define  EcaCwdiOperators       {"eca_cwdi"}
#define  EcaHwfiOperators       {"eca_hwfi"}
#define  EcaIdOperators         {"eca_id"}
#define  EcaSuOperators         {"eca_su"}
#define  EcaTrOperators         {"eca_tr"}
#define  EcaTg10pOperators      {"eca_tg10p"}
#define  EcaTg90pOperators      {"eca_tg90p"}
#define  EcaTn10pOperators      {"eca_tn10p"}
#define  EcaTn90pOperators      {"eca_tn90p"}
#define  EcaTx10pOperators      {"eca_tx10p"}
#define  EcaTx90pOperators      {"eca_tx90p"}

#define  EcaCddOperators        {"eca_cdd"}
#define  EcaCwdOperators        {"eca_cwd"}
#define  EcaRr1Operators        {"eca_rr1"}
#define  EcaR10mmOperators      {"eca_r10mm"}
#define  EcaR20mmOperators      {"eca_r20mm"}
#define  EcaR75pOperators       {"eca_r75p"}
#define  EcaR75ptotOperators    {"eca_r75ptot"}
#define  EcaR90pOperators       {"eca_r90p"}
#define  EcaR90ptotOperators    {"eca_r90ptot"}
#define  EcaR95pOperators       {"eca_r95p"}
#define  EcaR95ptotOperators    {"eca_r95ptot"}
#define  EcaR99pOperators       {"eca_r99p"}
#define  EcaR99ptotOperators    {"eca_r99ptot"}
#define  EcaRx1dayOperators     {"eca_rx1day"}
#define  EcaRx5dayOperators     {"eca_rx5day"}
#define  EcaSdiiOperators       {"eca_sdii"}

#define  EcaFdnsOperators       {"eca_fdns"}

#define  EcaStrwinOperators     {"eca_strwin"}
#define  EcaStrbreOperators     {"eca_strbre"}
#define  EcaStrgalOperators     {"eca_strgal"}
#define  EcaHurrOperators       {"eca_hurr"}

#define  HiOperators            {"hi"}
#define  WctOperators           {"wct"}
/* QR */


static MODULES Modules[] =
{
  /*
    function        help function      operator names          num streams
                                                               in  out
  */
  { Arith,          ArithHelp,         ArithOperators,          2,  1 },
  { Arithc,         ArithcHelp,        ArithcOperators,         1,  1 },
  { Arithdays,      ArithdaysHelp,     ArithdaysOperators,      1,  1 },
  { Arithlat,       NULL,              ArithlatOperators,       1,  1 },
  { Cat,            CopyHelp,          CatOperators,           -1,  1 },
  { Change,         ChangeHelp,        ChangeOperators,         1,  1 },
  { Comp,           CompHelp,          CompOperators,           2,  1 },
  { Compc,          CompcHelp,         CompcOperators,          1,  1 },
  { Cond,           CondHelp,          CondOperators,           2,  1 },
  { Cond2,          Cond2Help,         Cond2Operators,          3,  1 },
  { Condc,          CondcHelp,         CondcOperators,          1,  1 },
  { Copy,           CopyHelp,          CopyOperators,          -1,  1 },
  { Detrend,        DetrendHelp,       DetrendOperators,        1,  1 },
  { Diff,           DiffHelp,          DiffOperators,           2,  0 },
  { Echam5ini,      NULL,              Echam5iniOperators,      1,  1 },
  { Enlarge,        EnlargeHelp,       EnlargeOperators,        1,  1 },
  { Enlargegrid,    NULL,              EnlargegridOperators,    1,  1 },
  { Ensstat,        EnsstatHelp,       EnsstatOperators,       -1,  1 },
  { Expr,           ExprHelp,          ExprOperators,           1,  1 },
  { Filedes,        FiledesHelp,       FiledesOperators,        1,  0 },
  { Fillmiss,       NULL,              FillmissOperators,       1,  1 },
  { Fldrms,         NULL,              FldrmsOperators,         2,  1 },
  { Fldstat,        FldstatHelp,       FldstatOperators,        1,  1 },
  { Gradsdes,       GradsdesHelp,      GradsdesOperators,       1,  0 },
  { Histogram,      NULL,              HistogramOperators,      1,  1 },
  { Info,           InfoHelp,          InfoOperators,          -1,  0 },
  { Input,          InputHelp,         InputOperators,          0,  1 },
  { Intgrid,        IntgridHelp,       IntgridOperators,        1,  1 },
  { Intgridtraj,    NULL,              IntgridtrajOperators,    1,  1 },
  { Inttime,        InttimeHelp,       InttimeOperators,        1,  1 },
  { Intntime,       InttimeHelp,       IntntimeOperators,       1,  1 },
  { Intyear,        IntyearHelp,       IntyearOperators,        2,  1 },
  { Invert,         InvertHelp,        InvertOperators,         1,  1 },
  { Log,            NULL,              LogOperators,            1,  0 },
  { Maskbox,        MaskboxHelp,       MaskboxOperators,        1,  1 },
  { Mastrfu,        MastrfuHelp,       MastrfuOperators,        1,  1 },
  { Math,           MathHelp,          MathOperators,           1,  1 },
  { Merge,          MergeHelp,         MergeOperators,         -1,  1 },
  { Mergegrid,      NULL,              MergegridOperators,      2,  1 },
  { Mergetime,      MergeHelp,         MergetimeOperators,     -1,  1 },
  { Merstat,        MerstatHelp,       MerstatOperators,        1,  1 },
  { Mrotuv,         NULL,              MrotuvOperators,         2,  1 },
  { Ninfo,          NinfoHelp,         NinfoOperators,          1,  0 },
  { Nmltest,        NULL,              NmltestOperators,        0,  0 },
  { Output,         OutputHelp,        OutputOperators,        -1,  0 },
  { Outputgmt,      NULL,              OutputgmtOperators,      1,  0 },
  { Pinfo,          NULL,              PinfoOperators,          1,  1 },
  { Remap,          RemapHelp,         RemapOperators,          1,  1 },
  { Remap,          RemapgridHelp,     RemapgridOperators,      1,  1 },
  { Remap,          GenweightsHelp,    GenweightsOperators,     1,  1 },
  { Replace,        ReplaceHelp,       ReplaceOperators,        2,  1 },
  { Rotuv,          RotuvHelp,         RotuvOperators,          1,  1 },
  /* RQ */
  { Runpctl,        RunpctlHelp,       RunpctlOperators,        1,  1 },
  /* QR */
  { Runstat,        RunstatHelp,       RunstatOperators,        1,  1 },
  /* RQ */
  { Seaspctl,       SeaspctlHelp,      SeaspctlOperators,       3,  1 },
  /* QR */
  { Seasstat,       SeasstatHelp,      SeasstatOperators,       1,  1 },
  { Selbox,         SelboxHelp,        SelboxOperators,         1,  1 },
  { Select,         SelectHelp,        SelectOperators,         1,  1 },
  { Seloperator,    NULL,              SeloperatorOperators,    1,  1 },
  { Selrec,         SelectHelp,        SelrecOperators,         1,  1 },
  /* RQ */
  { Selpctl,        SelpctlHelp,       SelpctlOperators,        3,  1 },
  /* QR */
  { Selstat,        SelstatHelp,       SelstatOperators,        1,  1 },
  { Seltime,        SeltimeHelp,       SeltimeOperators,        1,  1 },
  { Set,            SetHelp,           SetOperators,            1,  1 },
  { Setbox,         SetboxHelp,        SetboxOperators,         1,  1 },
  { Setgatt,        SetgattHelp,       SetgattOperators,        1,  1 },
  { Setgrid,        SetgridHelp,       SetgridOperators,        1,  1 },
  { Sethalo,        NULL,              SethaloOperators,        1,  1 },
  { Setrcaname,     NULL,              SetrcanameOperators,     1,  1 },
  { Setmiss,        SetmissHelp,       SetmissOperators,        1,  1 },
  { Settime,        SettimeHelp,       SettimeOperators,        1,  1 },
  { Setzaxis,       SetzaxisHelp,      SetzaxisOperators,       1,  1 },
  { Showinfo,       ShowinfoHelp,      ShowinfoOperators,       1,  0 },
  { Sinfo,          SinfoHelp,         SinfoOperators,         -1,  0 },
  { Sort,           NULL,              SortOperators,           1,  1 },
  { Specinfo,       NULL,              SpecinfoOperators,       0,  0 },
  { Spectral,       SpectralHelp,      SpectralOperators,       1,  1 },
  { Split,          SplitHelp,         SplitOperators,          1,  1 },
  { Splitrec,       SplitHelp,         SplitrecOperators,       1,  1 },
  { Splittime,      SplittimeHelp,     SplittimeOperators,      1,  1 },
  { Splityear,      SplittimeHelp,     SplityearOperators,      1,  1 },
  { Subtrend,       SubtrendHelp,      SubtrendOperators,       3,  1 },
  { Template1,      NULL,              Template1Operators,      1,  1 },
  { Template2,      NULL,              Template2Operators,      1,  1 },
  { Test,           NULL,              TestOperators,           1,  1 },
  { Test2,          NULL,              Test2Operators,          2,  1 },
  /* RQ */
  { Timpctl,        TimpctlHelp,       TimpctlOperators,        3,  1 },
  { Timpctl,        YearpctlHelp,      YearpctlOperators,       3,  1 },
  { Timpctl,        MonpctlHelp,       MonpctlOperators,        3,  1 },
  { Timpctl,        DaypctlHelp,       DaypctlOperators,        3,  1 },
  { Timpctl,        HourpctlHelp,      HourpctlOperators,       3,  1 },
  /* QR */
  { Timsort,        TimsortHelp,       TimsortOperators,        1,  1 },
  { Timstat,        TimstatHelp,       TimstatOperators,        1,  1 },
  { Timstat,        YearstatHelp,      YearstatOperators,       1,  1 },
  { Timstat,        MonstatHelp,       MonstatOperators,        1,  1 },
  { Timstat,        DaystatHelp,       DaystatOperators,        1,  1 },
  { Timstat,        HourstatHelp,      HourstatOperators,       1,  1 },
  { Trend,          TrendHelp,         TrendOperators,          1,  2 },
  { Trms,           NULL,              TrmsOperators,           2,  1 },
  { Vardup,         VardupHelp,        VardupOperators,         1,  1 },
  { Vargen,         VargenHelp,        VargenOperators,         0,  1 },
  { Varrms,         NULL,              VarrmsOperators,         2,  1 },
  { Vertint,        IntvertHelp,       VertintOperators,        1,  1 },
  { Vertstat,       VertstatHelp,      VertstatOperators,       1,  1 },
  { Wind,           WindHelp,          WindOperators,           1,  1 },
  { Writegrid,      NULL,              WritegridOperators,      1,  1 },  /* no cdi output */
  { Writerandom,    NULL,              WriterandomOperators,    1,  1 },
  /* RQ */
  { Ydaypctl,       YdaypctlHelp,      YdaypctlOperators,       3,  1 },
  /* QR */
  { Ydaystat,       YdaystatHelp,      YdaystatOperators,       1,  1 },
  /* RQ */
  { Ydrunpctl,      YdrunpctlHelp,     YdrunpctlOperators,      3,  1 },
  { Ydrunstat,      YdrunstatHelp,     YdrunstatOperators,      1,  1 },
  /* QR */
  { Ymonarith,      YmonarithHelp,     YmonarithOperators,      2,  1 },
  /* RQ */
  { Ymonpctl,       YmonpctlHelp,      YmonpctlOperators,       3,  1 },
  /* QR */
  { Ymonstat,       YmonstatHelp,      YmonstatOperators,       1,  1 },
  /* RQ */
  { Yseaspctl,      YseaspctlHelp,     YseaspctlOperators,      3,  1 },
  /* QR */
  { Yseasstat,      YseasstatHelp,     YseasstatOperators,      1,  1 },
  { Zonstat,        ZonstatHelp,       ZonstatOperators,        1,  1 },
  /* RQ */
  { EcaCfd,         EcaCfdHelp,        EcaCfdOperators,         1,  1 },
  { EcaCsu,         EcaCsuHelp,        EcaCsuOperators,         1,  1 },
  { EcaCwdi,        EcaCwdiHelp,       EcaCwdiOperators,        2,  1 },
  { EcaCwfi,        EcaCwfiHelp,       EcaCwfiOperators,        2,  1 },
  { EcaEtr,         EcaEtrHelp,        EcaEtrOperators,         2,  1 },
  { EcaFd,          EcaFdHelp,         EcaFdOperators,          1,  1 },
  { EcaGsl,         EcaGslHelp,        EcaGslOperators,         1,  1 },
  { EcaHd,          EcaHdHelp,         EcaHdOperators,          1,  1 },
  { EcaHwdi,        EcaHwdiHelp,       EcaHwdiOperators,        2,  1 },
  { EcaHwfi,        EcaHwfiHelp,       EcaHwfiOperators,        2,  1 },
  { EcaId,          EcaIdHelp,         EcaIdOperators,          1,  1 },
  { EcaSu,          EcaSuHelp,         EcaSuOperators,          1,  1 },
  { EcaTr,          EcaTrHelp,         EcaTrOperators,          1,  1 },
  { EcaTg10p,       EcaTg10pHelp,      EcaTg10pOperators,       2,  1 },
  { EcaTg90p,       EcaTg90pHelp,      EcaTg90pOperators,       2,  1 },
  { EcaTn10p,       EcaTn10pHelp,      EcaTn10pOperators,       2,  1 },
  { EcaTn90p,       EcaTn90pHelp,      EcaTn90pOperators,       2,  1 },
  { EcaTx10p,       EcaTx10pHelp,      EcaTx10pOperators,       2,  1 },
  { EcaTx90p,       EcaTx90pHelp,      EcaTx90pOperators,       2,  1 },
  { EcaCdd,         EcaCddHelp,        EcaCddOperators,         1,  1 },
  { EcaCwd,         EcaCwdHelp,        EcaCwdOperators,         1,  1 },
  { EcaRr1,         EcaRr1Help,        EcaRr1Operators,         1,  1 },
  { EcaR10mm,       EcaR10mmHelp,      EcaR10mmOperators,       1,  1 },
  { EcaR20mm,       EcaR20mmHelp,      EcaR20mmOperators,       1,  1 },
  { EcaR75p,        EcaR75pHelp,       EcaR75pOperators,        2,  1 },
  { EcaR75ptot,     EcaR75ptotHelp,    EcaR75ptotOperators,     2,  1 },
  { EcaR90p,        EcaR90pHelp,       EcaR90pOperators,        2,  1 },
  { EcaR90ptot,     EcaR90ptotHelp,    EcaR90ptotOperators,     2,  1 },
  { EcaR95p,        EcaR95pHelp,       EcaR95pOperators,        2,  1 },
  { EcaR95ptot,     EcaR95ptotHelp,    EcaR95ptotOperators,     2,  1 },
  { EcaR99p,        EcaR99pHelp,       EcaR99pOperators,        2,  1 },
  { EcaR99ptot,     EcaR99ptotHelp,    EcaR99ptotOperators,     2,  1 },
  { EcaRx1day,      EcaRx1dayHelp,     EcaRx1dayOperators,      1,  1 },
  { EcaRx5day,      EcaRx5dayHelp,     EcaRx5dayOperators,      1,  1 },
  { EcaSdii,        EcaSdiiHelp,       EcaSdiiOperators,        1,  1 },
  { EcaFdns,        EcaFdnsHelp,       EcaFdnsOperators,        2,  1 },
  { EcaStrwin,      EcaStrwinHelp,     EcaStrwinOperators,      1,  1 },
  { EcaStrbre,      EcaStrbreHelp,     EcaStrbreOperators,      1,  1 },
  { EcaStrgal,      EcaStrgalHelp,     EcaStrgalOperators,      1,  1 },
  { EcaHurr,        EcaHurrHelp,       EcaHurrOperators,        1,  1 },

  { Hi,             HiHelp,            HiOperators,             3,  1 },
  { Wct,            WctHelp,           WctOperators,            2,  1 },
  /* QR */
};

static int NumModules = sizeof(Modules) / sizeof(Modules[0]);

static char *opalias[][2] =
{
  {"anomaly",             "ymonsub"    },
  {"ggstat",              "info"       },
  {"ggstats",             "sinfo"      },
  {"globavg",             "fldavg"     },
  {"gradsdes",            "gradsdes2"  },
  {"gp2sp2",              "gp2spl"     },
  {"sp2gp2",              "sp2gpl"     },
  {"infos",               "sinfo"      },
  {"intgrid",             "intgridbil" },
  {"log",                 "ln"         },
  {"lmean",               "ymonmean"   },
  {"lmmean",              "ymonmean"   },
  {"lmavg",               "ymonavg"    },
  {"lmstd",               "ymonstd"    },
  {"lsmean",              "yseasmean"  },
  {"remaplin",            "remapbil"   },
  {"remapcub",            "remapbic"   },
  {"remapbilinear",       "remapbil"   },
  {"remapbicubic",        "remapbic"   },
  {"remapconservative",   "remapcon"   },
  {"remapdistance",       "remapdis"   },
  {"sort",                "timsort"    },
  {"vinfos",              "sinfov"     },
  /* RQ */
  {"eca_r1mm",            "eca_rr1"    },
  /* QR */  
};

static int nopalias = sizeof(opalias) / (2*sizeof(opalias[0][0]));

static int similar(char *a, char *b, int level_a, int level_b)
{
  while ( *a && *b && *a == *b )
    { 
      a++ ;
      b++ ;
    }
  if ( !*a && !*b )
    return TRUE ;
  /*
  printf("%d %d %s %s\n", level_a, level_b, a, b);
  */

  if ( level_a >= 2 && level_b >= 1 && *a && 
       similar(a+1, b, level_a-2, level_b-1) )
    return TRUE ;

  if ( level_a >= 1 && level_b >= 2 && *b && 
       similar(a, b+1, level_a-1, level_b-2 ) )
    return TRUE ;

  return FALSE ; 
}

char *operatorAlias(char *operatorName)
{
  char *operatorNameNew;
  int i;

  operatorNameNew = operatorName;

  for ( i = 0; i < nopalias; i++ )
    {
      /*   printf("%d %d %s %s\n", nopalias, i, opalias[i][0], opalias[i][1]); */
      if ( strcmp(operatorName, opalias[i][0]) == 0 ) break;
    }

  if ( i < nopalias )
    {
      /* fprintf(stdout, "%s is an alias for %s\n", operatorName, opalias[i][1]); */
      operatorNameNew = opalias[i][1];
    }

  return (operatorNameNew);
}

static int operatorInqModID(char *operatorName)
{
  static char func[] = "operatorInqModID";
  int i, j, modID = -1;

  if ( operatorName )
    {
      for ( i = 0; i < NumModules; i++ )
	{
	  j = 0;
	  for ( j = 0; j < MAX_MOD_OPERATORS; j++ )
	    {
	      if ( Modules[i].operators[j] == NULL ) break;

	      if ( operatorName[0] == Modules[i].operators[j][0] )
		{
		  if ( strcmp(operatorName, Modules[i].operators[j]) == 0 )
		    {
		      modID = i;
		      break;
		    }
		}
	    }
	  if ( modID != -1 ) break;
	}
    }
  
  if ( modID == -1 )
    {
      FILE *fp;
      int nbyte;
      int error = TRUE;
      fp = fopen(operatorName, "r");
      if ( fp )
	{
	  fclose(fp);
	  fprintf(stderr, "Use commandline option -h for help.");
	  Error(func, "operator missing! %s is a file on disk!", operatorName);
	}
      fprintf(stderr, "Operator >%s< not found!\n", operatorName);
      fprintf(stderr, "Similar operators are:\n");
      nbyte = fprintf(stderr, "   ");
      if ( operatorName )
	for ( i = 0; i < NumModules; i++ )
	  {
	    if ( Modules[i].help == NULL ) continue;
	    j = 0;
	    while ( Modules[i].operators[j] )
	      {
		if( similar(operatorName, Modules[i].operators[j],
			    strlen(operatorName), strlen(Modules[i].operators[j])) )
		  {
		    if ( nbyte > 75 )
		      {
			fprintf(stdout, "\n");
			nbyte = fprintf(stderr, "   ");
		      }
		    nbyte += fprintf(stderr, " %s", Modules[i].operators[j]);
		    error = FALSE ;
		  }
		j++;
	      }
	  }
      if ( error )
	fprintf(stderr, "(not found)\n") ;
      else
	fprintf(stderr, "\n");

      exit(EXIT_FAILURE);
    }

  if ( modID != -1 )
    if ( ! Modules[modID].func )
      Error(func, "Module for operator >%s< not installed!", operatorName);

  return (modID);
}

void *(*operatorModule(char *operatorName))(void *)
{
  int modID;
  modID = operatorInqModID(operatorName);
  return (Modules[modID].func);
}

char **operatorHelp(char *operatorName)
{
  int modID;
  modID = operatorInqModID(operatorName);
  return (Modules[modID].help);
}

int operatorStreamInCnt(char *operatorName)
{
  int modID;
  modID = operatorInqModID(operatorAlias(operatorName));
  return (Modules[modID].streamInCnt);
}

int operatorStreamOutCnt(char *operatorName)
{
  int modID;
  modID = operatorInqModID(operatorAlias(operatorName));
  return (Modules[modID].streamOutCnt);
}

int cmpname(const void *s1, const void *s2)
{
  char **c1 = (char **) s1;
  char **c2 = (char **) s2;

  return (strcmp((const char *)*c1, (const char *)*c2));
}

void operatorPrintAll(void)
{
  int i, j, nbyte, nop = 0;
  char *opernames[4096];

  for ( i = 0; i < NumModules; i++ )
    {
      if ( Modules[i].help == NULL ) continue;
      j = 0;
      while ( Modules[i].operators[j] )
	{
	  opernames[nop++] = Modules[i].operators[j++];
	}
    }

  qsort(opernames, nop, sizeof(char *), cmpname);

  nbyte = fprintf(stderr, "   ");
  for ( i = 0; i < nop; i++ )
    {
      if ( nbyte > 65 )
	{
	  fprintf(stdout, "\n");
	  nbyte = fprintf(stderr, "   ");
	}
      nbyte += fprintf(stderr, " %s", opernames[i]);
    }
  fprintf(stderr, "\n");
}

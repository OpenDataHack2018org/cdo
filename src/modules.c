/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2009 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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
#include "cdi.h"
#include "cdo.h"
#include "operator_help.h"
#include "modules.h"
#include "error.h"


#define  MAX_MOD_OPERATORS  128         /* maximum number of operators for a module */

typedef struct {
  void  *(*func)(void *);               /* Module                   */
  char **help;                          /* Help                     */
  char  *operators[MAX_MOD_OPERATORS];  /* Operator names           */
  short  number;                        /* Allowed number type      */
  short  streamInCnt;                   /* Number of input streams  */
  short  streamOutCnt;                  /* Number of output streams */
}
modules_t;


void *Arith(void *argument);
void *Arithc(void *argument);
void *Arithdays(void *argument);
void *Arithlat(void *argument);
void *Cat(void *argument);
void *CDItest(void *argument);
void *Change(void *argument);
void *Change_e5slm(void *argument);
void *Comp(void *argument);
void *Compc(void *argument);
void *Cond(void *argument);
void *Cond2(void *argument);
void *Condc(void *argument);
void *Copy(void *argument);
void *Del29feb(void *argument);
void *Detrend(void *argument);
void *Diff(void *argument);
void *Echam5ini(void *argument);
void *Enlarge(void *argument);
void *Enlargegrid(void *argument);
void *Ensstat(void *argument);
void *Eofcoeff(void *argument);
void *Timeof(void *argument);
void *Expr(void *argument);
void *Filedes(void *argument);
void *Fillmiss(void *argument);
void *Filter(void *argument);
void *Fldrms(void *argument);
void *Fldstat(void *argument);
void *Gengrid(void *argument);
void *Gradsdes(void *argument);
void *Gridboxstat(void *argument);
void *Gridcell(void *argument);
void *Harmonic(void *argument);
void *Histogram(void *argument);
void *Importamsr(void *argument);
void *Importbinary(void *argument);
void *Importcmsaf(void *argument);
void *Importobs(void *argument);
void *Info(void *argument);
void *Input(void *argument);
void *Intgrid(void *argument);
void *Intgridtraj(void *argument);
void *Intlevel(void *argument);
void *Inttime(void *argument);
void *Intntime(void *argument);
void *Intyear(void *argument);
void *Invert(void *argument);
void *Invertlev(void *argument);
void *Log(void *argument);
void *Maskbox(void *argument);
void *Mastrfu(void *argument);
void *Math(void *argument);
void *Merge(void *argument);
void *Mergegrid(void *argument);
void *Mergetime(void *argument);
void *Merstat(void *argument);
void *Monarith(void *argument);
void *Mrotuv(void *argument);
void *Mrotuvb(void *argument);
void *Ninfo(void *argument);
void *Nmltest(void *argument);
void *Output(void *argument);
void *Outputgmt(void *argument);
void *Pressure(void *argument);
void *Pinfo(void *argument);
void *Regres(void *argument);
void *Remap(void *argument);
void *Remapeta(void *argument);
void *Replace(void *argument);
void *Rotuv(void *argument);
void *Runpctl(void *argument);
void *Runstat(void *argument);
void *Seascount(void *argument);
void *Seaspctl(void *argument);
void *Seasstat(void *argument);
void *Selbox(void *argument);
void *Select(void *argument);
void *Selvar(void *argument);
void *Seloperator(void *argument);
void *Selrec(void *argument);
void *Seltime(void *argument);
void *Set(void *argument);
void *Setbox(void *argument);
void *Setgatt(void *argument);
void *Setgrid(void *argument);
void *Sethalo(void *argument);
void *Setmiss(void *argument);
void *Setrange(void *argument);
void *Setrcaname(void *argument);
void *Settime(void *argument);
void *Setzaxis(void *argument);
void *Showinfo(void *argument);
void *Sinfo(void *argument);
void *Smooth9(void *argument);
void *Sort(void *argument);
void *Sorttimestamp(void *argument);
void *Specinfo(void *argument);
void *Spectral(void *argument);
void *Spectrum(void *argument);
void *Split(void *argument);
void *Splitrec(void *argument);
void *Splitsel(void *argument);
void *Splittime(void *argument);
void *Splityear(void *argument);
void *Subtrend(void *argument);
void *Template1(void *argument);
void *Template2(void *argument);
void *Test(void *argument);
void *Test2(void *argument);
void *Testdata(void *argument);
void *Tests(void *argument);
void *Timsort(void *argument);
void *Timcount(void *argument);
void *Timpctl(void *argument);
void *Timselpctl(void *argument);
void *Timselstat(void *argument);
void *Timstat(void *argument);
void *Timstat2(void *argument);
void *Timstat3(void *argument);
void *Tinfo(void *argument);
void *Tocomplex(void *argument);
void *Transpose(void *argument);
void *Trend(void *argument);
void *Trms(void *argument);
void *Vardup(void *argument);
void *Vargen(void *argument);
void *Varrms(void *argument);
void *Vertint(void *argument);
void *Vertstat(void *argument);
void *Vertwind(void *argument);
void *Wind(void *argument);
void *Writegrid(void *argument);
void *Writerandom(void *argument);
void *Yhourstat(void *argument);
void *Ydaypctl(void *argument);
void *Ydaystat(void *argument);
void *Ydrunpctl(void *argument);
void *Ydrunstat(void *argument);
void *Ymonarith(void *argument);
void *Ymonpctl(void *argument);
void *Ymonstat(void *argument);
void *Yseaspctl(void *argument);
void *Yseasstat(void *argument);
void *Zonstat(void *argument);
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

void *Fdns(void *argument);
void *Strwin(void *argument);
void *Strbre(void *argument);
void *Strgal(void *argument);
void *Hurr(void *argument);

void *Hi(void *argument);
void *Wct(void *argument);


#define  ArithOperators         {"add",  "sub",  "mul",  "div", "min", "max", "atan2"}
#define  ArithcOperators        {"addc", "subc", "mulc", "divc"}
#define  ArithdaysOperators     {"muldpm", "divdpm", "muldpy", "divdpy"}
#define  ArithlatOperators      {"mulcoslat", "divcoslat"}
#define  CatOperators           {"cat"}
#define  CDItestOperators       {"ncopy"}
#define  ChangeOperators        {"chcode", "chname", "chlevel", "chlevelc", "chlevelv", "chltype"}
#define  Change_e5slmOperators  {"change_e5slm", "change_e5lsm", "change_e5mask"}
#define  CompOperators          {"eq",  "ne",  "le",  "lt",  "ge",  "gt"}
#define  CompcOperators         {"eqc", "nec", "lec", "ltc", "gec", "gtc"}
#define  CondOperators          {"ifthen",  "ifnotthen"}
#define  Cond2Operators         {"ifthenelse"}
#define  CondcOperators         {"ifthenc", "ifnotthenc"}
#define  CopyOperators          {"copy", "selall"}
#define  Del29febOperators      {"del29feb"}
#define  DetrendOperators       {"detrend"}
#define  DiffOperators          {"diff", "diffp", "diffv", "sdiff"}
#define  Echam5iniOperators     {"import_e5ml", "import_e5res", \
                                 "export_e5ml", "export_e5res"}
#define  EnlargeOperators       {"enlarge"}
#define  EnlargegridOperators   {"enlargegrid"}
#define  EnsstatOperators       {"ensmin", "ensmax", "enssum", "ensmean", "ensavg", "ensvar", "ensstd", "enspctl"}
#define  EofcoeffOperators      {"eofcoeff"}
#define  TimeofOperators        {"eof", "eofspatial", "eoftime"}
#define  ExprOperators          {"expr", "exprf", "aexpr", "aexprf"}
#define  FiledesOperators       {"filedes", "griddes", "griddes2", "zaxisdes", "vct", "vct2", "pardes", \
                                 "vlist", "partab", "partab2"}
#define  FillmissOperators      {"fillmiss"}
#define  FilterOperators        {"bandpass", "highpass", "lowpass"}
#define  FldrmsOperators        {"fldrms"}
#define  FldstatOperators       {"fldmin", "fldmax", "fldsum", "fldmean", "fldavg", "fldvar", "fldstd", "fldpctl"}
#define  GengridOperators       {"gengrid"}
#define  GradsdesOperators      {"gradsdes1", "gradsdes2", "dumpmap"}
#define  GridboxstatOperators   {"gridboxmin", "gridboxmax", "gridboxsum", "gridboxmean", "gridboxavg", \
                                 "gridboxvar", "gridboxstd"}
#define  GridcellOperators      {"gridarea", "gridweights"}
#define  HarmonicOperators      {"harmonic"}
#define  HistogramOperators     {"histcount", "histsum", "histmean", "histfreq"}
#define  ImportamsrOperators    {"import_amsr"}
#define  ImportbinaryOperators  {"import_binary", "import_grads"}
#define  ImportcmsafOperators   {"import_cmsaf"}
#define  ImportobsOperators     {"import_obs"}
#define  InfoOperators          {"info", "infov", "infop", "map"}
#define  InputOperators         {"input", "inputsrv", "inputext"}
#define  IntgridOperators       {"intgridbil", "intpoint", "interpolate", "intarea", "boxavg", "thinout"}
#define  IntgridtrajOperators   {"intgridtraj"}
#define  IntlevelOperators      {"intlevel", "intlevelx"}
#define  InttimeOperators       {"inttime"}
#define  IntntimeOperators      {"intntime"}
#define  IntyearOperators       {"intyear"}
#define  InvertOperators        {"invertlat", "invertlon", "invertlatdes", "invertlondes", \
                                 "invertlatdata", "invertlondata"}
#define  InvertlevOperators     {"invertlev"}
#define  LogOperators           {"dumplogs", "daylogs", "monlogs", "dumplogo", \
                                 "snamelogo", "scalllogo", "smemlogo", "stimelogo", "sperclogo"}
#define  MaskboxOperators       {"masklonlatbox", "maskindexbox"}
#define  MaskregionOperators    {"maskregion"}
#define  MastrfuOperators       {"mastrfu"}
#define  MathOperators          {"abs", "int", "nint", "sqr", "sqrt", "exp", "ln", "log10", "sin", \
                                 "cos", "tan", "asin", "acos", "atan", "pow", "reci"}
#define  MergeOperators         {"merge"}
#define  MergegridOperators     {"mergegrid"}
#define  MergetimeOperators     {"mergetime"}
#define  MerstatOperators       {"mermin", "mermax", "mersum", "mermean", "meravg", "mervar", "merstd", "merpctl"}
#define  MonarithOperators      {"monadd", "monsub", "monmul", "mondiv"}
#define  MrotuvOperators        {"mrotuv"}
#define  MrotuvbOperators       {"mrotuvb"}
#define  NinfoOperators         {"nyear", "nmon", "ndate", "ntime", "ncode", "npar", "nlevel"}
#define  NmltestOperators       {"nmltest"}
#define  OutputOperators        {"output", "outputint", "outputsrv", "outputext", "outputf", "outputts", \
                                 "outputfld", "outputarr", "outputxyz", "outputkey"}
#define  OutputgmtOperators     {"gridverify", "outputcenter", "outputcenter2", "outputcentercpt", "outputbounds", \
                                 "outputboundscpt", "outputvector", "outputtri"}
#define  PressureOperators      {"fpressure", "hpressure"}
#define  PinfoOperators         {"pinfo", "pinfov"}
#define  RegresOperators        {"regres"}
#define  RemapOperators         {"remap"}
#define    RemapgridOperators   {"remapcon", "remapbil", "remapbic", "remapdis", "remapnn", "remaplaf", "remapcon2", "remapsum"}
#define    GenweightsOperators  {"gencon", "genbil", "genbic", "gendis", "gennn", "genlaf", "gencon2"}
#define  RemapetaOperators      {"remapeta"}
#define  ReplaceOperators       {"replace"}
#define  RotuvOperators         {"rotuvb"}
#define  RunpctlOperators       {"runpctl"}
#define  RunstatOperators       {"runmin",  "runmax",  "runsum",  "runmean",  "runavg",  "runvar",  "runstd"}
#define  SeascountOperators     {"seascount"}
#define  SeaspctlOperators      {"seaspctl"}
#define  SeasstatOperators      {"seasmin",  "seasmax",  "seassum",  "seasmean",  "seasavg",  "seasvar",  "seasstd"}
#define  SelboxOperators        {"sellonlatbox", "selindexbox"}
#define  SelectOperators        {"select"}
#define  SelvarOperators        {"selparam", "selcode", "selname", "selstdname", "sellevel", "sellevidx", "selgrid", \
                                 "selzaxis", "seltabnum", "delparam", "delcode", "delname", "selltype"}
#define  SeloperatorOperators   {"seloperator"}
#define  SelrecOperators        {"selrec"}
#define  SeltimeOperators       {"seltimestep", "selyear", "selseas", "selmon", "selday", "selhour", "seldate", \
                                 "seltime", "selsmon"}
#define  SetOperators           {"setpartab", "setpartabv", "setcode", "setname", "setlevel", "setltype", "settabnum"}
#define  SetboxOperators        {"setclonlatbox", "setcindexbox"}
#define  SetgattOperators       {"setgatt", "setgatts"}
#define  SetgridOperators       {"setgrid", "setgridtype", "setgridarea"}
#define  SethaloOperators       {"sethalo", "tpnhalo"}
#define  SetmissOperators       {"setmissval", "setctomiss", "setmisstoc", "setrtomiss", "setvrange"}
#define  SetrangeOperators      {"setrtoc", "setrtoc2"}
#define  SetrcanameOperators    {"setrcaname"}
#define  SettimeOperators       {"setyear", "setmon", "setday", "setdate", "settime", "settunits", \
                                 "settaxis", "setreftime", "setcalendar", "shifttime"}
#define  SetzaxisOperators      {"setzaxis"}
#define  ShowinfoOperators      {"showyear", "showmon", "showdate", "showtime", "showtimestamp", "showcode", \
                                 "showparam", "showname", "showstdname", "showlevel", "showltype", "showformat"}
#define  SinfoOperators         {"sinfo", "sinfov", "sinfop"}
#define  Smooth9Operators       {"smooth9"}
#define  SortOperators          {"sortcode", "sortname", "sortlevel"}
#define  SorttimestampOperators {"sorttimestamp", "sorttaxis"}
#define  SpecinfoOperators      {"specinfo"}
#define  SpectralOperators      {"gp2sp", "gp2spl", "sp2gp", "sp2gpl", "sp2sp", "spcut"}
#define  SpectrumOperators      {"spectrum"}
#define  SplitOperators         {"splitcode", "splitname", "splitlevel", "splitgrid", "splitzaxis", "splittabnum"}
#define  SplitrecOperators      {"splitrec"}
#define  SplitselOperators      {"splitsel"}
#define  SplittimeOperators     {"splithour", "splitday", "splitmon", "splitseas"}
#define  SplityearOperators     {"splityear"}
#define  SubtrendOperators      {"subtrend"}
#define  Template1Operators     {"template1"}
#define  Template2Operators     {"template2"}
#define  TestOperators          {"test"}
#define  Test2Operators         {"test2"}
#define  TestdataOperators      {"testdata"}
#define  TestsOperators         {"normal", "studentt", "chisquare", "beta", "fisher"}
#define  TimsortOperators       {"timsort"}
#define  TimcountOperators      {"timcount"}
#define    YearcountOperators   {"yearcount"}
#define    MoncountOperators    {"moncount"}
#define    DaycountOperators    {"daycount"}
#define    HourcountOperators   {"hourcount"}
#define  TimpctlOperators       {"timpctl"}
#define    YearpctlOperators    {"yearpctl"}
#define    MonpctlOperators     {"monpctl"}
#define    DaypctlOperators     {"daypctl"}
#define    HourpctlOperators    {"hourpctl"}
#define  TimselpctlOperators    {"timselpctl"}
#define  TimselstatOperators    {"timselmin", "timselmax", "timselsum", "timselmean", "timselavg", "timselvar", "timselstd"}
#define  TimstatOperators       {"timmin",  "timmax",  "timsum",  "timmean",  "timavg",  "timvar",  "timstd"}
#define    YearstatOperators    {"yearmin", "yearmax", "yearsum", "yearmean", "yearavg", "yearvar", "yearstd"}
#define    MonstatOperators     {"monmin",  "monmax",  "monsum",  "monmean",  "monavg",  "monvar",  "monstd"}
#define    DaystatOperators     {"daymin",  "daymax",  "daysum",  "daymean",  "dayavg",  "dayvar",  "daystd"}
#define    HourstatOperators    {"hourmin", "hourmax", "hoursum", "hourmean", "houravg", "hourvar", "hourstd"}
#define  Timstat2Operators      {"timcor"}
#define  Timstat3Operators      {"meandiff2test", "varquot2test"}
#define  TinfoOperators         {"tinfo"}
#define  TocomplexOperators     {"retocomplex", "imtocomplex"}
#define  TransposeOperators     {"transxy"}
#define  TrendOperators         {"trend"}
#define  TrmsOperators          {"trms"}
#define  VardupOperators        {"pardup", "parmul"}
#define  VargenOperators        {"random", "const", "topo", "for"}
#define  VarrmsOperators        {"varrms"}
#define  VertintOperators       {"ml2pl", "ml2hl", "ml2plx", "ml2hlx", \
                                 "ml2pl_lp", "ml2hl_lp", "ml2plx_lp", "ml2hlx_lp"}
#define  VertstatOperators      {"vertmin", "vertmax", "vertsum", "vertmean", "vertavg", "vertvar", "vertstd"}
#define  VertwindOperators      {"vertwind"}
#define  WindOperators          {"uv2dv", "uv2dvl", "dv2uv", "dv2uvl", "dv2ps"}
#define  WritegridOperators     {"writegrid"}
#define  WriterandomOperators   {"writerandom"}
#define  YhourstatOperators     {"yhourmin", "yhourmax", "yhoursum", "yhourmean", "yhouravg", "yhourvar", "yhourstd"}
#define  YdaypctlOperators      {"ydaypctl"}
#define  YdaystatOperators      {"ydaymin", "ydaymax", "ydaysum", "ydaymean", "ydayavg", "ydayvar", "ydaystd"}
#define  YdrunpctlOperators     {"ydrunpctl"}
#define  YdrunstatOperators     {"ydrunmin", "ydrunmax", "ydrunsum", "ydrunmean", "ydrunavg", "ydrunvar", "ydrunstd"}
#define  YmonarithOperators     {"ymonadd", "ymonsub", "ymonmul", "ymondiv"}
#define  YmonpctlOperators      {"ymonpctl"}
#define  YmonstatOperators      {"ymonmin", "ymonmax", "ymonsum", "ymonmean", "ymonavg", "ymonvar", "ymonstd"}
#define  YseaspctlOperators     {"yseaspctl"}
#define  YseasstatOperators     {"yseasmin", "yseasmax", "yseassum", "yseasmean", "yseasavg", "yseasvar", "yseasstd"}
#define  ZonstatOperators       {"zonmin", "zonmax", "zonsum", "zonmean", "zonavg", "zonvar", "zonstd", "zonpctl"}

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

#define  FdnsOperators          {"fdns"}

#define  StrwinOperators        {"strwin"}
#define  StrbreOperators        {"strbre"}
#define  StrgalOperators        {"strgal"}
#define  HurrOperators          {"hurr"}

#define  HiOperators            {"hi"}
#define  WctOperators           {"wct"}


static modules_t Modules[] =
{
  /*
    function        help function      operator names          number     num streams
                                                               type       in  out
  */
  { Arith,          ArithHelp,         ArithOperators,         CDI_REAL,  2,  1 },
  { Arithc,         ArithcHelp,        ArithcOperators,        CDI_REAL,  1,  1 },
  { Arithdays,      ArithdaysHelp,     ArithdaysOperators,     CDI_REAL,  1,  1 },
  { Arithlat,       NULL,              ArithlatOperators,      CDI_REAL,  1,  1 },
  { Cat,            CopyHelp,          CatOperators,           CDI_REAL, -1,  1 },
  { CDItest,        NULL,              CDItestOperators,       CDI_REAL,  1,  1 },
  { Change,         ChangeHelp,        ChangeOperators,        CDI_REAL,  1,  1 },
  { Change_e5slm,   NULL,              Change_e5slmOperators,  CDI_REAL,  1,  1 },
  { Comp,           CompHelp,          CompOperators,          CDI_REAL,  2,  1 },
  { Compc,          CompcHelp,         CompcOperators,         CDI_REAL,  1,  1 },
  { Cond,           CondHelp,          CondOperators,          CDI_REAL,  2,  1 },
  { Cond2,          Cond2Help,         Cond2Operators,         CDI_REAL,  3,  1 },
  { Condc,          CondcHelp,         CondcOperators,         CDI_REAL,  1,  1 },
  { Copy,           CopyHelp,          CopyOperators,          CDI_REAL, -1,  1 },
  { Del29feb,       NULL,              Del29febOperators,      CDI_REAL,  1,  1 },
  { Detrend,        DetrendHelp,       DetrendOperators,       CDI_REAL,  1,  1 },
  { Diff,           DiffHelp,          DiffOperators,          CDI_REAL,  2,  0 },
  { Echam5ini,      NULL,              Echam5iniOperators,     CDI_REAL,  1,  1 },
  { Enlarge,        EnlargeHelp,       EnlargeOperators,       CDI_REAL,  1,  1 },
  { Enlargegrid,    NULL,              EnlargegridOperators,   CDI_REAL,  1,  1 },
  { Ensstat,        EnsstatHelp,       EnsstatOperators,       CDI_REAL, -1,  1 },
  { Eofcoeff,       NULL,              EofcoeffOperators,      CDI_REAL,  2,  1 },
  { Timeof,         NULL,              TimeofOperators,        CDI_REAL,  1,  2 },
  { Expr,           ExprHelp,          ExprOperators,          CDI_REAL,  1,  1 },
  { Filedes,        FiledesHelp,       FiledesOperators,       CDI_REAL,  1,  0 },
  { Fillmiss,       NULL,              FillmissOperators,      CDI_REAL,  1,  1 },
  { Filter,         NULL,              FilterOperators,        CDI_REAL,  1,  1 },
  { Fldrms,         NULL,              FldrmsOperators,        CDI_REAL,  2,  1 },
  { Fldstat,        FldstatHelp,       FldstatOperators,       CDI_REAL,  1,  1 },
  { Gengrid,        NULL,              GengridOperators,       CDI_REAL,  2,  1 },
  { Gradsdes,       GradsdesHelp,      GradsdesOperators,      CDI_REAL,  1,  0 },
  { Gridboxstat,    NULL,              GridboxstatOperators,   CDI_REAL,  1,  1 },
  { Gridcell,       GridcellHelp,      GridcellOperators,      CDI_REAL,  1,  1 },
  { Harmonic,       NULL,              HarmonicOperators,      CDI_REAL,  1,  1 },
  { Histogram,      HistogramHelp,     HistogramOperators,     CDI_REAL,  1,  1 },
  { Importamsr,     ImportamsrHelp,    ImportamsrOperators,    CDI_REAL,  1,  1 },
  { Importbinary,   ImportbinaryHelp,  ImportbinaryOperators,  CDI_REAL,  1,  1 },
  { Importcmsaf,    ImportcmsafHelp,   ImportcmsafOperators,   CDI_REAL,  1,  1 },
  { Importobs,      NULL,              ImportobsOperators,     CDI_REAL,  1,  1 },
  { Info,           InfoHelp,          InfoOperators,          CDI_BOTH, -1,  0 },
  { Input,          InputHelp,         InputOperators,         CDI_REAL,  0,  1 },
  { Intgrid,        IntgridHelp,       IntgridOperators,       CDI_REAL,  1,  1 },
  { Intgridtraj,    NULL,              IntgridtrajOperators,   CDI_REAL,  1,  1 },
  { Intlevel,       IntlevelHelp,      IntlevelOperators,      CDI_REAL,  1,  1 },
  { Inttime,        InttimeHelp,       InttimeOperators,       CDI_REAL,  1,  1 },
  { Intntime,       InttimeHelp,       IntntimeOperators,      CDI_REAL,  1,  1 },
  { Intyear,        IntyearHelp,       IntyearOperators,       CDI_REAL,  2,  1 },
  { Invert,         InvertHelp,        InvertOperators,        CDI_REAL,  1,  1 },
  { Invertlev,      InvertlevHelp,     InvertlevOperators,     CDI_REAL,  1,  1 },
  { Log,            NULL,              LogOperators,           CDI_REAL,  1,  0 },
  { Maskbox,        MaskboxHelp,       MaskboxOperators,       CDI_REAL,  1,  1 },
  { Maskbox,        MaskregionHelp,    MaskregionOperators,    CDI_REAL,  1,  1 },
  { Mastrfu,        MastrfuHelp,       MastrfuOperators,       CDI_REAL,  1,  1 },
  { Math,           MathHelp,          MathOperators,          CDI_REAL,  1,  1 },
  { Merge,          MergeHelp,         MergeOperators,         CDI_REAL, -1,  1 },
  { Mergegrid,      NULL,              MergegridOperators,     CDI_REAL,  2,  1 },
  { Mergetime,      MergeHelp,         MergetimeOperators,     CDI_REAL, -1,  1 },
  { Merstat,        MerstatHelp,       MerstatOperators,       CDI_REAL,  1,  1 },
  { Monarith,       MonarithHelp,      MonarithOperators,      CDI_REAL,  2,  1 },
  { Mrotuv,         NULL,              MrotuvOperators,        CDI_REAL,  1,  2 },
  { Mrotuvb,        NULL,              MrotuvbOperators,       CDI_REAL,  2,  1 },
  { Ninfo,          NinfoHelp,         NinfoOperators,         CDI_REAL,  1,  0 },
  { Nmltest,        NULL,              NmltestOperators,       CDI_REAL,  0,  0 },
  { Output,         OutputHelp,        OutputOperators,        CDI_REAL, -1,  0 },
  { Outputgmt,      NULL,              OutputgmtOperators,     CDI_REAL,  1,  0 },
  { Pressure,       NULL,              PressureOperators,      CDI_REAL,  1,  1 },
  { Pinfo,          NULL,              PinfoOperators,         CDI_REAL,  1,  1 },
  { Regres,         RegresHelp,        RegresOperators,        CDI_REAL,  1,  1 },
  { Remap,          RemapHelp,         RemapOperators,         CDI_REAL,  1,  1 },
  { Remap,          RemapgridHelp,     RemapgridOperators,     CDI_REAL,  1,  1 },
  { Remap,          GenweightsHelp,    GenweightsOperators,    CDI_REAL,  1,  1 },
  { Remapeta,       RemapetaHelp,      RemapetaOperators,      CDI_REAL,  1,  1 },
  { Replace,        ReplaceHelp,       ReplaceOperators,       CDI_REAL,  2,  1 },
  { Rotuv,          RotuvHelp,         RotuvOperators,         CDI_REAL,  1,  1 },
  { Runpctl,        RunpctlHelp,       RunpctlOperators,       CDI_REAL,  1,  1 },
  { Runstat,        RunstatHelp,       RunstatOperators,       CDI_REAL,  1,  1 },
  { Seascount,      NULL,              SeascountOperators,     CDI_REAL,  1,  1 },
  { Seaspctl,       SeaspctlHelp,      SeaspctlOperators,      CDI_REAL,  3,  1 },
  { Seasstat,       SeasstatHelp,      SeasstatOperators,      CDI_REAL,  1,  1 },
  { Selbox,         SelboxHelp,        SelboxOperators,        CDI_REAL,  1,  1 },
  { Select,         NULL,              SelectOperators,        CDI_REAL, -1,  1 },
  { Selvar,         SelvarHelp,        SelvarOperators,        CDI_REAL,  1,  1 },
  { Selrec,         SelvarHelp,        SelrecOperators,        CDI_REAL,  1,  1 },
  { Seloperator,    NULL,              SeloperatorOperators,   CDI_REAL,  1,  1 },
  { Seltime,        SeltimeHelp,       SeltimeOperators,       CDI_REAL,  1,  1 },
  { Set,            SetHelp,           SetOperators,           CDI_REAL,  1,  1 },
  { Setbox,         SetboxHelp,        SetboxOperators,        CDI_REAL,  1,  1 },
  { Setgatt,        SetgattHelp,       SetgattOperators,       CDI_REAL,  1,  1 },
  { Setgrid,        SetgridHelp,       SetgridOperators,       CDI_REAL,  1,  1 },
  { Sethalo,        SethaloHelp,       SethaloOperators,       CDI_REAL,  1,  1 },
  { Setrange,       SetrangeHelp,      SetrangeOperators,      CDI_REAL,  1,  1 },
  { Setrcaname,     NULL,              SetrcanameOperators,    CDI_REAL,  1,  1 },
  { Setmiss,        SetmissHelp,       SetmissOperators,       CDI_REAL,  1,  1 },
  { Settime,        SettimeHelp,       SettimeOperators,       CDI_REAL,  1,  1 },
  { Setzaxis,       SetzaxisHelp,      SetzaxisOperators,      CDI_REAL,  1,  1 },
  { Showinfo,       ShowinfoHelp,      ShowinfoOperators,      CDI_REAL,  1,  0 },
  { Sinfo,          SinfoHelp,         SinfoOperators,         CDI_BOTH, -1,  0 },
  { Smooth9,        Smooth9Help,       Smooth9Operators,       CDI_REAL,  1,  1 },
  { Sort,           NULL,              SortOperators,          CDI_REAL,  1,  1 },
  { Sorttimestamp,  NULL,              SorttimestampOperators, CDI_REAL, -1,  1 },
  { Specinfo,       NULL,              SpecinfoOperators,      CDI_REAL,  0,  0 },
  { Spectral,       SpectralHelp,      SpectralOperators,      CDI_REAL,  1,  1 },
  { Spectrum,       NULL,              SpectrumOperators,      CDI_REAL,  1,  1 },
  { Split,          SplitHelp,         SplitOperators,         CDI_REAL,  1,  1 },
  { Splitrec,       SplitHelp,         SplitrecOperators,      CDI_REAL,  1,  1 },
  { Splitsel,       SplitselHelp,      SplitselOperators,      CDI_REAL,  1,  1 },
  { Splittime,      SplittimeHelp,     SplittimeOperators,     CDI_REAL,  1,  1 },
  { Splityear,      SplittimeHelp,     SplityearOperators,     CDI_REAL,  1,  1 },
  { Subtrend,       SubtrendHelp,      SubtrendOperators,      CDI_REAL,  3,  1 },
  { Template1,      NULL,              Template1Operators,     CDI_REAL,  1,  1 },
  { Template2,      NULL,              Template2Operators,     CDI_REAL,  1,  1 },
  { Test,           NULL,              TestOperators,          CDI_REAL,  1,  1 },
  { Test2,          NULL,              Test2Operators,         CDI_REAL,  2,  1 },
  { Testdata,       NULL,              TestdataOperators,      CDI_REAL,  1,  1 },
  { Tests,          NULL,              TestsOperators,         CDI_REAL,  1,  1 },
  { Timcount,       NULL,              TimcountOperators,      CDI_REAL,  1,  1 },
  { Timcount,       NULL,              YearcountOperators,     CDI_REAL,  1,  1 },
  { Timcount,       NULL,              MoncountOperators,      CDI_REAL,  1,  1 },
  { Timcount,       NULL,              DaycountOperators,      CDI_REAL,  1,  1 },
  { Timcount,       NULL,              HourcountOperators,     CDI_REAL,  1,  1 },
  { Timpctl,        TimpctlHelp,       TimpctlOperators,       CDI_REAL,  3,  1 },
  { Timpctl,        YearpctlHelp,      YearpctlOperators,      CDI_REAL,  3,  1 },
  { Timpctl,        MonpctlHelp,       MonpctlOperators,       CDI_REAL,  3,  1 },
  { Timpctl,        DaypctlHelp,       DaypctlOperators,       CDI_REAL,  3,  1 },
  { Timpctl,        HourpctlHelp,      HourpctlOperators,      CDI_REAL,  3,  1 },
  { Timselpctl,     TimselpctlHelp,    TimselpctlOperators,    CDI_REAL,  3,  1 },
  { Timsort,        TimsortHelp,       TimsortOperators,       CDI_REAL,  1,  1 },
  { Timselstat,     TimselstatHelp,    TimselstatOperators,    CDI_REAL,  1,  1 },
  { Timstat,        TimstatHelp,       TimstatOperators,       CDI_REAL,  1,  1 },
  { Timstat,        YearstatHelp,      YearstatOperators,      CDI_REAL,  1,  1 },
  { Timstat,        MonstatHelp,       MonstatOperators,       CDI_REAL,  1,  1 },
  { Timstat,        DaystatHelp,       DaystatOperators,       CDI_REAL,  1,  1 },
  { Timstat,        HourstatHelp,      HourstatOperators,      CDI_REAL,  1,  1 },
  { Timstat2,       NULL,              Timstat2Operators,      CDI_REAL,  2,  1 },
  { Timstat3,       NULL,              Timstat3Operators,      CDI_REAL,  2,  1 },
  { Tinfo,          NULL,              TinfoOperators,         CDI_REAL,  1,  0 },
  { Tocomplex,      NULL,              TocomplexOperators,     CDI_REAL,  1,  1 },
  { Transpose,      NULL,              TransposeOperators,     CDI_REAL,  1,  1 },
  { Trend,          TrendHelp,         TrendOperators,         CDI_REAL,  1,  2 },
  { Trms,           NULL,              TrmsOperators,          CDI_REAL,  2,  1 },
  { Vardup,         NULL,              VardupOperators,        CDI_REAL,  1,  1 },
  { Vargen,         VargenHelp,        VargenOperators,        CDI_REAL,  0,  1 },
  { Varrms,         NULL,              VarrmsOperators,        CDI_REAL,  2,  1 },
  { Vertint,        IntvertHelp,       VertintOperators,       CDI_REAL,  1,  1 },
  { Vertstat,       VertstatHelp,      VertstatOperators,      CDI_REAL,  1,  1 },
  { Vertwind,       NULL,              VertwindOperators,      CDI_REAL,  1,  1 },
  { Wind,           WindHelp,          WindOperators,          CDI_REAL,  1,  1 },
  { Writegrid,      NULL,              WritegridOperators,     CDI_REAL,  1,  1 },  /* no cdi output */
  { Writerandom,    NULL,              WriterandomOperators,   CDI_REAL,  1,  1 },
  { Yhourstat,      YhourstatHelp,     YhourstatOperators,     CDI_REAL,  1,  1 },
  { Ydaypctl,       YdaypctlHelp,      YdaypctlOperators,      CDI_REAL,  3,  1 },
  { Ydaystat,       YdaystatHelp,      YdaystatOperators,      CDI_REAL,  1,  1 },
  { Ydrunpctl,      YdrunpctlHelp,     YdrunpctlOperators,     CDI_REAL,  3,  1 },
  { Ydrunstat,      YdrunstatHelp,     YdrunstatOperators,     CDI_REAL,  1,  1 },
  { Ymonarith,      YmonarithHelp,     YmonarithOperators,     CDI_REAL,  2,  1 },
  { Ymonpctl,       YmonpctlHelp,      YmonpctlOperators,      CDI_REAL,  3,  1 },
  { Ymonstat,       YmonstatHelp,      YmonstatOperators,      CDI_REAL,  1,  1 },
  { Yseaspctl,      YseaspctlHelp,     YseaspctlOperators,     CDI_REAL,  3,  1 },
  { Yseasstat,      YseasstatHelp,     YseasstatOperators,     CDI_REAL,  1,  1 },
  { Zonstat,        ZonstatHelp,       ZonstatOperators,       CDI_REAL,  1,  1 },
  { EcaCfd,         EcaCfdHelp,        EcaCfdOperators,        CDI_REAL,  1,  1 },
  { EcaCsu,         EcaCsuHelp,        EcaCsuOperators,        CDI_REAL,  1,  1 },
  { EcaCwdi,        EcaCwdiHelp,       EcaCwdiOperators,       CDI_REAL,  2,  1 },
  { EcaCwfi,        EcaCwfiHelp,       EcaCwfiOperators,       CDI_REAL,  2,  1 },
  { EcaEtr,         EcaEtrHelp,        EcaEtrOperators,        CDI_REAL,  2,  1 },
  { EcaFd,          EcaFdHelp,         EcaFdOperators,         CDI_REAL,  1,  1 },
  { EcaGsl,         EcaGslHelp,        EcaGslOperators,        CDI_REAL,  2,  1 },
  { EcaHd,          EcaHdHelp,         EcaHdOperators,         CDI_REAL,  1,  1 },
  { EcaHwdi,        EcaHwdiHelp,       EcaHwdiOperators,       CDI_REAL,  2,  1 },
  { EcaHwfi,        EcaHwfiHelp,       EcaHwfiOperators,       CDI_REAL,  2,  1 },
  { EcaId,          EcaIdHelp,         EcaIdOperators,         CDI_REAL,  1,  1 },
  { EcaSu,          EcaSuHelp,         EcaSuOperators,         CDI_REAL,  1,  1 },
  { EcaTr,          EcaTrHelp,         EcaTrOperators,         CDI_REAL,  1,  1 },
  { EcaTg10p,       EcaTg10pHelp,      EcaTg10pOperators,      CDI_REAL,  2,  1 },
  { EcaTg90p,       EcaTg90pHelp,      EcaTg90pOperators,      CDI_REAL,  2,  1 },
  { EcaTn10p,       EcaTn10pHelp,      EcaTn10pOperators,      CDI_REAL,  2,  1 },
  { EcaTn90p,       EcaTn90pHelp,      EcaTn90pOperators,      CDI_REAL,  2,  1 },
  { EcaTx10p,       EcaTx10pHelp,      EcaTx10pOperators,      CDI_REAL,  2,  1 },
  { EcaTx90p,       EcaTx90pHelp,      EcaTx90pOperators,      CDI_REAL,  2,  1 },
  { EcaCdd,         EcaCddHelp,        EcaCddOperators,        CDI_REAL,  1,  1 },
  { EcaCwd,         EcaCwdHelp,        EcaCwdOperators,        CDI_REAL,  1,  1 },
  { EcaRr1,         EcaRr1Help,        EcaRr1Operators,        CDI_REAL,  1,  1 },
  { EcaR10mm,       EcaR10mmHelp,      EcaR10mmOperators,      CDI_REAL,  1,  1 },
  { EcaR20mm,       EcaR20mmHelp,      EcaR20mmOperators,      CDI_REAL,  1,  1 },
  { EcaR75p,        EcaR75pHelp,       EcaR75pOperators,       CDI_REAL,  2,  1 },
  { EcaR75ptot,     EcaR75ptotHelp,    EcaR75ptotOperators,    CDI_REAL,  2,  1 },
  { EcaR90p,        EcaR90pHelp,       EcaR90pOperators,       CDI_REAL,  2,  1 },
  { EcaR90ptot,     EcaR90ptotHelp,    EcaR90ptotOperators,    CDI_REAL,  2,  1 },
  { EcaR95p,        EcaR95pHelp,       EcaR95pOperators,       CDI_REAL,  2,  1 },
  { EcaR95ptot,     EcaR95ptotHelp,    EcaR95ptotOperators,    CDI_REAL,  2,  1 },
  { EcaR99p,        EcaR99pHelp,       EcaR99pOperators,       CDI_REAL,  2,  1 },
  { EcaR99ptot,     EcaR99ptotHelp,    EcaR99ptotOperators,    CDI_REAL,  2,  1 },
  { EcaRx1day,      EcaRx1dayHelp,     EcaRx1dayOperators,     CDI_REAL,  1,  1 },
  { EcaRx5day,      EcaRx5dayHelp,     EcaRx5dayOperators,     CDI_REAL,  1,  1 },
  { EcaSdii,        EcaSdiiHelp,       EcaSdiiOperators,       CDI_REAL,  1,  1 },
  { Fdns,           FdnsHelp,          FdnsOperators,          CDI_REAL,  2,  1 },
  { Strwin,         StrwinHelp,        StrwinOperators,        CDI_REAL,  1,  1 },
  { Strbre,         StrbreHelp,        StrbreOperators,        CDI_REAL,  1,  1 },
  { Strgal,         StrgalHelp,        StrgalOperators,        CDI_REAL,  1,  1 },
  { Hurr,           HurrHelp,          HurrOperators,          CDI_REAL,  1,  1 },
  /*  { Hi,             NULL,              HiOperators,        CDI_REAL,  3,  1 }, */
  { Wct,            WctHelp,           WctOperators,           CDI_REAL,  2,  1 },
};							       
							       
static int NumModules = sizeof(Modules) / sizeof(Modules[0]);

static char *opalias[][2] =
{
  {"anomaly",             "ymonsub"    },
  {"ggstat",              "info"       },
  {"ggstats",             "sinfo"      },
  {"globavg",             "fldavg"     },
  {"gradsdes",            "gradsdes2"  },
  {"infos",               "sinfo"      },
  {"intgrid",             "intgridbil" },
  {"log",                 "ln"         },
  {"lmean",               "ymonmean"   },
  {"lmmean",              "ymonmean"   },
  {"lmavg",               "ymonavg"    },
  {"lmstd",               "ymonstd"    },
  {"lsmean",              "yseasmean"  },
  {"chvar",               "chname"     },
  {"ncode",               "npar"       },
  {"nvar",                "npar"       },
  {"vardes",              "pardes"     },
  {"delvar",              "delname"    },
  {"vardup",              "pardup"     },
  {"varmul",              "parmul"     },
  {"read_e5ml",           "import_e5ml"},
  {"remapcon1",           "remaplaf"   },
  {"remapdis1",           "remapnn"    },
  {"showvar",             "showname"   },
  {"selgridname",         "selgrid"    },
  {"selzaxisname",        "selzaxis"   },
  {"selvar",              "selname"    },
  {"setvar",              "setname"    },
  {"sortvar",             "sortname"   },
  {"splitvar",            "splitname"  },
  {"sort",                "timsort"    },
  {"write_e5ml",          "export_e5ml"},
  {"eca_r1mm",            "eca_rr1"    },
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

int operatorStreamNumber(char *operatorName)
{
  int modID;
  modID = operatorInqModID(operatorAlias(operatorName));
  return (Modules[modID].number);
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

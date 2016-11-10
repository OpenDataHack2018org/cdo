/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2016 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <cdi.h>
#include "cdo_int.h"
#include "operator_help.h"
#include "modules.h"
#include "error.h"


#define  MAX_MOD_OPERATORS  128               // maximum number of operators for a module

typedef struct {
  void  *(*func)(void *);                     // Module
  const char **help;                          // Help
  const char  *operators[MAX_MOD_OPERATORS];  // Operator names
  short  mode;                                // Module mode: 0:intern 1:extern
  short  number;                              // Allowed number type
  short  streamInCnt;                         // Number of input streams
  short  streamOutCnt;                        // Number of output streams
}
modules_t;


void *Adisit(void *argument);
void *Afterburner(void *argument);
void *Arith(void *argument);
void *Arithc(void *argument);
void *Arithdays(void *argument);
void *Arithlat(void *argument);
void *Cat(void *argument);
void *CDItest(void *argument);
void *CDIread(void *argument);
void *CDIwrite(void *argument);
void *Change(void *argument);
void *Change_e5slm(void *argument);
void *Cloudlayer(void *argument);
void *CMOR(void *argument);
void *CMOR_lite(void *argument);
void *CMOR_table(void *argument);
void *Collgrid(void *argument);
void *Command(void *argument);
void *Comp(void *argument);
void *Compc(void *argument);
void *Complextorect(void *argument);
void *Cond(void *argument);
void *Cond2(void *argument);
void *Condc(void *argument);
void *Consecstat(void *argument);
void *Copy(void *argument);
void *Deltat(void *argument);
void *Deltime(void *argument);
void *Derivepar(void *argument);
void *Detrend(void *argument);
void *Diff(void *argument);
void *Distgrid(void *argument);
void *Duplicate(void *argument);
void *Echam5ini(void *argument);
void *Enlarge(void *argument);
void *Enlargegrid(void *argument);
void *Ensstat(void *argument);
void *Ensstat3(void *argument);
void *Ensval(void *argument);
void *Eofcoeff(void *argument);
void *Eofcoeff3d(void *argument);
void *EOFs(void *argument);
void *EOF3d(void *argument);
void *Expr(void *argument);
void *FC(void *argument);
void *Filedes(void *argument);
void *Fillmiss(void *argument);
void *Filter(void *argument);
void *Fldrms(void *argument);
void *Fldstat(void *argument);
void *Fldstat2(void *argument);
void *Fourier(void *argument);
void *Gengrid(void *argument);
void *Gradsdes(void *argument);
void *Gridboxstat(void *argument);
void *Gridcell(void *argument);
void *Gridsearch(void *argument);
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
void *Intlevel3d(void *argument);
void *Inttime(void *argument);
void *Intntime(void *argument);
void *Intyear(void *argument);
void *Invert(void *argument);
void *Invertlev(void *argument);
void *Isosurface(void *argument);
void *Log(void *argument);
void *MapReduce(void *argument);
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
void *Nmldump(void *argument);
void *Output(void *argument);
void *Outputgmt(void *argument);
void *Pack(void *argument);
void *Pardup(void *argument);
void *Pinfo(void *argument);
void *Pressure(void *argument);
void *Regres(void *argument);
void *Remap(void *argument);
void *Remapeta(void *argument);
void *Replace(void *argument);
void *Replacevalues(void *argument);
void *Rotuv(void *argument);
void *Rhopot(void *argument);
void *Runpctl(void *argument);
void *Runstat(void *argument);
void *Seascount(void *argument);
void *Seaspctl(void *argument);
void *Seasstat(void *argument);
void *Selbox(void *argument);
void *Selindex(void *argument);
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
void *Setpartab(void *argument);
void *Setrcaname(void *argument);
void *Settime(void *argument);
void *Setzaxis(void *argument);
void *Shiftxy(void *argument);
void *Showinfo(void *argument);
void *Sinfo(void *argument);
void *Smooth(void *argument);
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
void *Tee(void *argument);
void *Template1(void *argument);
void *Template2(void *argument);
void *Test(void *argument);
void *Test2(void *argument);
void *Testdata(void *argument);
void *Tests(void *argument);
void *Timsort(void *argument);
void *Timcount(void *argument);
void *Timcumsum(void *argument);
void *Timpctl(void *argument);
void *Timselpctl(void *argument);
void *Timselstat(void *argument);
void *XTimstat(void *argument);
void *Timstat(void *argument);
void *Timstat2(void *argument);
void *Timstat3(void *argument);
void *Tinfo(void *argument);
void *Tocomplex(void *argument);
void *Transpose(void *argument);
void *Trend(void *argument);
void *Trms(void *argument);
void *Tstepcount(void *argument);
void *Vargen(void *argument);
void *Varrms(void *argument);
void *Vertintml(void *argument);
void *Vertintap(void *argument);
void *Vertstat(void *argument);
void *Vertcum(void *argument);
void *Vertwind(void *argument);
void *Verifygrid(void *argument);
void *Wind(void *argument);
void *Writegrid(void *argument);
void *Writerandom(void *argument);
void *YAR(void *argument);
void *Yearmonstat(void *argument);
void *Ydayarith(void *argument);
void *Ydaypctl(void *argument);
void *Ydaystat(void *argument);
void *Ydrunpctl(void *argument);
void *Ydrunstat(void *argument);
void *Yhourarith(void *argument);
void *Yhourstat(void *argument);
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
void *EcaPd(void *argument);
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

//void *Hi(void *argument);
void *Wct(void *argument);

void *Magplot(void *argument);
void *Magvector(void *argument);
void *Maggraph(void *argument);


#define  AdisitOperators        {"adisit", "adipot"}
#define  AfterburnerOperators   {"after"}
#define  ArithOperators         {"add",  "sub",  "mul",  "div", "min", "max", "atan2"}
#define  ArithcOperators        {"addc", "subc", "mulc", "divc", "mod"}
#define  ArithdaysOperators     {"muldpm", "divdpm", "muldpy", "divdpy", "muldoy"}
#define  ArithlatOperators      {"mulcoslat", "divcoslat"}
#define  CatOperators           {"cat"}
#define  CDItestOperators       {"ncopy"}
#define  CDIreadOperators       {"cdiread"}
#define  CDIwriteOperators      {"cdiwrite"}
#define  ChangeOperators        {"chcode", "chtabnum", "chparam", "chname", "chunit", "chlevel", "chlevelc", "chlevelv", "chltype"}
#define  Change_e5slmOperators  {"change_e5slm", "change_e5lsm", "change_e5mask"}
#define  CloudlayerOperators    {"cloudlayer"}
#define  CMOROperators          {"cmor"}
#define  CMORliteOperators      {"cmorlite"}
#define  CMORtableOperators     {"dump_cmor_table", "conv_cmor_table"}
#define  CollgridOperators      {"collgrid"}
#define  CommandOperators       {"command", "com", "cmd"}
#define  CompOperators          {"eq",  "ne",  "le",  "lt",  "ge",  "gt"}
#define  CompcOperators         {"eqc", "nec", "lec", "ltc", "gec", "gtc"}
#define  ComplextorectOperators {"complextorect"}
#define  CondOperators          {"ifthen",  "ifnotthen"}
#define  Cond2Operators         {"ifthenelse"}
#define  CondcOperators         {"ifthenc", "ifnotthenc"}
#define  ConsecstatOperators    {"consects", "consecsum"}
#define  CopyOperators          {"copy", "selall", "szip"}
#define  DeltatOperators        {"deltat"}
#define  DeltimeOperators       {"delday", "del29feb"}
#define  DeriveparOperators     {"gheight", "sealevelpressure"}
#define  DetrendOperators       {"detrend"}
#define  DiffOperators          {"diff", "diffp", "diffn", "diffc"}
#define  DistgridOperators      {"distgrid"}
#define  DuplicateOperators     {"duplicate"}
#define  Echam5iniOperators     {"import_e5ml", "import_e5res", \
                                 "export_e5ml", "export_e5res"}
#define  EnlargeOperators       {"enlarge"}
#define  EnlargegridOperators   {"enlargegrid"}
#define  EnsstatOperators       {"ensmin", "ensmax", "enssum", "ensmean", "ensavg", "ensvar", "ensvar1", "ensstd", "ensstd1", "enspctl"}
#define  Ensstat3Operators      {"ensrkhistspace", "ensrkhisttime", "ensroc"}
#define  EnsvalOperators        {"enscrps", "ensbrs"}
#define  EofcoeffOperators      {"eofcoeff"}
#define  Eofcoeff3dOperators    {"eofcoeff3d"}
#define  EOFsOperators          {"eof", "eofspatial", "eoftime"}
#define  EOF3dOperators         {"eof3d","eof3dspatial","eof3dtime"}
#define  ExprOperators          {"expr", "exprf", "aexpr", "aexprf"}
#define  FCOperators            {"fc2sp", "sp2fc", "fc2gp", "gp2fc"}
#define  FiledesOperators       {"filedes", "griddes", "griddes2", "zaxisdes", "vct", "vct2", "codetab", \
                                 "vlist", "partab", "partab2", "spartab"}
#define  FillmissOperators      {"fillmiss", "fillmiss2"}
#define  FilterOperators        {"bandpass", "highpass", "lowpass"}
#define  FldrmsOperators        {"fldrms"}
#define  FldstatOperators       {"fldmin", "fldmax", "fldsum", "fldmean", "fldavg", "fldstd", "fldstd1", "fldvar", "fldvar1", "fldpctl"}
#define  FldcorOperators        {"fldcor"}
#define  FldcovarOperators      {"fldcovar"}
#define  FourierOperators       {"fourier"}
#define  GengridOperators       {"gengrid"}
#define  GradsdesOperators      {"gradsdes", "dumpmap"}
#define  GridboxstatOperators   {"gridboxmin", "gridboxmax", "gridboxsum", "gridboxmean", "gridboxavg", "gridboxstd", "gridboxstd1", "gridboxvar", "gridboxvar1"}
#define  GridcellOperators      {"gridarea", "gridweights", "gridmask", "griddx", "griddy"}
#define  GridsearchOperators    {"testpointsearch", "testcellsearch"}
#define  HarmonicOperators      {"harmonic"}
#define  HistogramOperators     {"histcount", "histsum", "histmean", "histfreq"}
#define  ImportamsrOperators    {"import_amsr"}
#define  ImportbinaryOperators  {"import_binary"}
#define  ImportcmsafOperators   {"import_cmsaf"}
#define  ImportobsOperators     {"import_obs"}
#define  InfoOperators          {"info", "infop", "infon", "infoc", "map"}
#define  InputOperators         {"input", "inputsrv", "inputext"}
#define  IntgridOperators       {"intgridbil", "intpoint", "interpolate", "boxavg", "thinout"}
#define  IntgridtrajOperators   {"intgridtraj"}
#define  IntlevelOperators      {"intlevel", "intlevelx"}
#define  Intlevel3dOperators    {"intlevel3d", "intlevelx3d"}
#define  InttimeOperators       {"inttime"}
#define  IntntimeOperators      {"intntime"}
#define  IntyearOperators       {"intyear"}
#define  InvertOperators        {"invertlat", "invertlon", "invertlatdes", "invertlondes", "invertlatdata", "invertlondata"}
#define  InvertlevOperators     {"invertlev"}
#define  IsosurfaceOperators    {"isosurface"}
#define  LogOperators           {"dumplogs", "daylogs", "monlogs", "dumplogo", "snamelogo", "scalllogo", "smemlogo", "stimelogo", "sperclogo"}
#define  MapReduceOperators     {"reducegrid"}
#define  MaskboxOperators       {"masklonlatbox", "maskindexbox"}
#define  MaskregionOperators    {"maskregion"}
#define  MastrfuOperators       {"mastrfu"}
#define  MathOperators          {"abs", "int", "nint", "sqr", "sqrt", "exp", "ln", "log10", "sin", \
                                 "cos", "tan", "asin", "acos", "atan", "pow", "reci"}
#define  MergeOperators         {"merge"}
#define  MergegridOperators     {"mergegrid"}
#define  MergetimeOperators     {"mergetime"}
#define  MerstatOperators       {"mermin", "mermax", "mersum", "mermean", "meravg", "merstd", "merstd1", "mervar", "mervar1", "merpctl"}
#define  MonarithOperators      {"monadd", "monsub", "monmul", "mondiv"}
#define  MrotuvOperators        {"mrotuv"}
#define  MrotuvbOperators       {"mrotuvb"}
#define  NinfoOperators         {"nyear", "nmon", "ndate", "ntime", "ncode", "npar", "nlevel", "ngridpoints", "ngrids"}
#define  NmldumpOperators       {"nmldump"}
#define  OutputOperators        {"output", "outputint", "outputsrv", "outputext", "outputf", "outputts", \
                                 "outputfld", "outputarr", "outputxyz"}
#define  OutputtabOperators     {"outputtab"}
#define  OutputgmtOperators     {"gmtxyz", "gmtcells", "outputcenter2", "outputcentercpt", \
                                 "outputboundscpt", "outputvector", "outputtri", "outputvrml"}
#define  PackOperators          {"pack"}
#define  PardupOperators        {"pardup", "parmul"}
#define  PinfoOperators         {"pinfo", "pinfov"}
#define  PressureOperators      {"pressure_fl", "pressure_hl", "deltap"}
#define  RegresOperators        {"regres"}
#define  RemapOperators         {"remap"}
#define  RemapbilOperators      {"remapbil", "genbil"}
#define  RemapbicOperators      {"remapbic", "genbic"}
#define  RemapnnOperators       {"remapnn", "gennn"}
#define  RemapdisOperators      {"remapdis", "gendis"}
#define  RemapyconOperators     {"remapycon", "genycon"}
#define  RemapconOperators      {"remapcon", "gencon"}
#define  Remapcon2Operators     {"remapcon2", "gencon2"}
#define  RemaplafOperators      {"remaplaf", "genlaf"}
#define    RemapgridOperators   {"remapsum"}
#define  RemapetaOperators      {"remapeta", "remapeta_s", "remapeta_z"}
#define  ReplaceOperators       {"replace"}
#define  ReplacevaluesOperators {"setvals", "setrtoc", "setrtoc2"}
#define  RhopotOperators        {"rhopot"}
#define  RotuvOperators         {"rotuvb"}
#define  RunpctlOperators       {"runpctl"}
#define  RunstatOperators       {"runmin", "runmax", "runsum", "runmean", "runavg", "runstd", "runstd1", "runvar", "runvar1"}
#define  SeascountOperators     {"seascount"}
#define  SeaspctlOperators      {"seaspctl"}
#define  SeasstatOperators      {"seasmin", "seasmax", "seassum", "seasmean", "seasavg", "seasstd", "seasstd1", "seasvar", "seasvar1"}
#define  SelboxOperators        {"sellonlatbox", "selindexbox"}
#define  SelindexOperators      {"selindex"}
#define  SelectOperators        {"select", "delete"}
#define  SelvarOperators        {"selparam", "selcode", "selname", "selstdname", "sellevel", "sellevidx", "selgrid", \
                                 "selzaxis", "selzaxisname", "seltabnum", "delparam", "delcode", "delname", "selltype"}
#define  SeloperatorOperators   {"seloperator"}
#define  SelrecOperators        {"selrec"}
#define  SeltimeOperators       {"seltimestep", "selyear", "selseason", "selmonth", "selday", "selhour", "seldate", \
                                 "seltime", "selsmon"}
#define  SetOperators           {"setcode", "setparam", "setname", "setunit", "setlevel", "setltype", "settabnum"}
#define  SetboxOperators        {"setclonlatbox", "setcindexbox"}
#define  SetgattOperators       {"setgatt", "setgatts"}
#define  SetgridOperators       {"setgrid", "setgridtype", "setgridarea", "setgridmask", "unsetgridmask", "setgridnumber", "setgriduri"}
#define  SethaloOperators       {"sethalo", "tpnhalo"}
#define  SetmissOperators       {"setmissval", "setctomiss", "setmisstoc", "setrtomiss", "setvrange"}
#define  SetmisstonnOperators   {"setmisstonn", "setmisstodis"}
#define  SetcodetabOperators    {"setcodetab"}
#define  SetpartabOperators     {"setpartabc", "setpartabp", "setpartabn"}
#define  SetrcanameOperators    {"setrcaname"}
#define  SettimeOperators       {"setyear", "setmon", "setday", "setdate", "settime", "settunits", \
                                 "settaxis", "settbounds", "setreftime", "setcalendar", "shifttime"}
#define  SetzaxisOperators      {"setzaxis", "genlevelbounds"}
#define  ShiftxyOperators       {"shiftx", "shifty"}
#define  ShowinfoOperators      {"showyear", "showmon", "showdate", "showtime", "showtimestamp", "showcode", "showunit", \
                                 "showparam", "showname", "showstdname", "showlevel", "showltype", "showformat"}
#define  SinfoOperators         {"sinfo", "sinfop", "sinfon", "sinfoc", "seinfo", "seinfop", "seinfon", "seinfoc"}
#define  SmoothOperators        {"smooth", "smooth9"}
#define  SortOperators          {"sortcode", "sortparam", "sortname", "sortlevel"}
#define  SorttimestampOperators {"sorttimestamp", "sorttaxis"}
#define  SpecinfoOperators      {"specinfo"}
#define  SpectralOperators      {"gp2sp", "gp2spl", "sp2gp", "sp2gpl", "sp2sp", "spcut"}
#define  SpectrumOperators      {"spectrum"}
#define  SplitOperators         {"splitcode", "splitparam", "splitname", "splitlevel", "splitgrid", "splitzaxis", "splittabnum"}
#define  SplitrecOperators      {"splitrec"}
#define  SplitselOperators      {"splitsel"}
#define  SplittimeOperators     {"splithour", "splitday", "splitmon", "splitseas"}
#define  SplityearOperators     {"splityear", "splityearmon"}
#define  SubtrendOperators      {"subtrend"}
#define  TeeOperators           {"tee"}
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
#define  TimcumsumOperators     {"timcumsum"}
#define  TimpctlOperators       {"timpctl"}
#define    YearpctlOperators    {"yearpctl"}
#define    MonpctlOperators     {"monpctl"}
#define    DaypctlOperators     {"daypctl"}
#define    HourpctlOperators    {"hourpctl"}
#define  TimselpctlOperators    {"timselpctl"}
#define  TimselstatOperators    {"timselmin", "timselmax", "timselsum", "timselmean", "timselavg", "timselvar", "timselvar1", "timselstd", "timselstd1"}
#define  XTimstatOperators      {"xtimmin",  "xtimmax",  "xtimsum",  "xtimmean",  "xtimavg",  "xtimvar",  "xtimvar1",  "xtimstd",  "xtimstd1", \
                                 "xyearmin", "xyearmax", "xyearsum", "xyearmean", "xyearavg", "xyearvar", "xyearvar1", "xyearstd", "xyearstd1", \
                                 "xmonmin",  "xmonmax",  "xmonsum",  "xmonmean",  "xmonavg",  "xmonvar",  "xmonvar1",  "xmonstd",  "xmonstd1"}
#define  TimstatOperators       {"timmin",  "timmax",  "timsum",  "timmean",  "timavg",  "timvar",  "timvar1",  "timstd",  "timstd1"}
#define    YearstatOperators    {"yearmin", "yearmax", "yearsum", "yearmean", "yearavg", "yearvar", "yearvar1", "yearstd", "yearstd1"}
#define    MonstatOperators     {"monmin",  "monmax",  "monsum",  "monmean",  "monavg",  "monvar",  "monvar1",  "monstd",  "monstd1"}
#define    DaystatOperators     {"daymin",  "daymax",  "daysum",  "daymean",  "dayavg",  "dayvar",  "dayvar1",  "daystd",  "daystd1"}
#define    HourstatOperators    {"hourmin", "hourmax", "hoursum", "hourmean", "houravg", "hourvar", "hourvar1", "hourstd", "hourstd1"}
#define  TimcorOperators        {"timcor"}
#define  TimcovarOperators      {"timcovar"}
#define  Timstat3Operators      {"meandiff2test", "varquot2test"}
#define  TinfoOperators         {"tinfo"}
#define  TocomplexOperators     {"retocomplex", "imtocomplex"}
#define  TransposeOperators     {"transxy"}
#define  TrendOperators         {"trend"}
#define  TrmsOperators          {"trms"}
#define  TstepcountOperators    {"tstepcount"}
#define  VargenOperators        {"random", "const", "sincos", "coshill", "for", "topo", "temp", "mask", "stdatm"}
#define  VarrmsOperators        {"varrms"}
#define  VertintmlOperators     {"ml2pl", "ml2hl", "ml2plx", "ml2hlx", "ml2pl_lp", "ml2hl_lp", "ml2plx_lp", "ml2hlx_lp"}
#define  VertintapOperators     {"ap2pl", "ap2plx", "ap2pl_lp", "ap2plx_lp", "ap2hl", "ap2hlx"}
#define  VertstatOperators      {"vertmin", "vertmax", "vertsum", "vertint", "vertmean", "vertavg", "vertstd", "vertstd1", "vertvar", "vertvar1"}
#define  VertcumOperators       {"vertcum", "vertcumhl"}
#define  VertwindOperators      {"vertwind"}
#define  VerifygridOperators    {"verifygrid"}
#define  WindOperators          {"uv2dv", "uv2dvl", "dv2uv", "dv2uvl", "dv2ps"}
#define  WritegridOperators     {"writegrid"}
#define  WriterandomOperators   {"writerandom"}
#define  YAROperators           {"yarbil", "yarnn", "yarcon"}
#define  YearmonstatOperators   {"yearmonmean", "yearmonavg"}
#define  YdayarithOperators     {"ydayadd", "ydaysub", "ydaymul", "ydaydiv"}
#define  YdaypctlOperators      {"ydaypctl"}
#define  YdaystatOperators      {"ydaymin", "ydaymax", "ydaysum", "ydaymean", "ydayavg", "ydaystd", "ydaystd1", "ydayvar", "ydayvar1"}
#define  YdrunpctlOperators     {"ydrunpctl"}
#define  YdrunstatOperators     {"ydrunmin", "ydrunmax", "ydrunsum", "ydrunmean", "ydrunavg", "ydrunstd", "ydrunstd1", "ydrunvar", "ydrunvar1"}
#define  YhourarithOperators    {"yhouradd", "yhoursub", "yhourmul", "yhourdiv"}
#define  YhourstatOperators     {"yhourmin", "yhourmax", "yhoursum", "yhourmean", "yhouravg", "yhourstd", "yhourstd1", "yhourvar", "yhourvar1"}
#define  YmonarithOperators     {"ymonadd", "ymonsub", "ymonmul", "ymondiv"}
#define  YseasarithOperators    {"yseasadd", "yseassub", "yseasmul", "yseasdiv"}
#define  YmonpctlOperators      {"ymonpctl"}
#define  YmonstatOperators      {"ymonmin", "ymonmax", "ymonsum", "ymonmean", "ymonavg", "ymonstd", "ymonstd1", "ymonvar", "ymonvar1"}
#define  YseaspctlOperators     {"yseaspctl"}
#define  YseasstatOperators     {"yseasmin", "yseasmax", "yseassum", "yseasmean", "yseasavg", "yseasstd", "yseasstd1", "yseasvar", "yseasvar1"}
#define  ZonstatOperators       {"zonmin", "zonmax", "zonrange", "zonsum", "zonmean", "zonavg", "zonstd", "zonstd1", "zonvar", "zonvar1", "zonpctl"}

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
/*
#define  EcaR10mmOperators      {"eca_r10mm"}
#define  EcaR20mmOperators      {"eca_r20mm"}
*/
#define  EcaPdOperators         {"eca_pd", "eca_r10mm", "eca_r20mm"}
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

#define  MagplotOperators       {"contour", "shaded", "grfill"}
#define  MagvectorOperators     {"vector"}
#define  MaggraphOperators      {"graph"}

static modules_t Modules[] =
{
  // stream in  -1 means: unlimited number of input streams
  // stream out -1 means: usage of obase
  /*
    function        help function      operator names          mode number     num streams
                                                                    type       in  out
  */
  { Adisit,         AdisitHelp,        AdisitOperators,        1,   CDI_REAL,  1,  1 },
  { Afterburner,    AfterburnerHelp,   AfterburnerOperators,   1,   CDI_REAL, -1,  1 },
  { Arith,          ArithHelp,         ArithOperators,         1,   CDI_REAL,  2,  1 },
  { Arithc,         ArithcHelp,        ArithcOperators,        1,   CDI_REAL,  1,  1 },
  { Arithdays,      ArithdaysHelp,     ArithdaysOperators,     1,   CDI_REAL,  1,  1 },
  { Arithlat,       NULL,              ArithlatOperators,      1,   CDI_REAL,  1,  1 },
  { Cat,            CopyHelp,          CatOperators,           1,   CDI_REAL, -1,  1 },
  { CDItest,        NULL,              CDItestOperators,       1,   CDI_REAL,  1,  1 },
  { CDIread,        NULL,              CDIreadOperators,       1,   CDI_REAL,  1,  0 },
  { CDIwrite,       NULL,              CDIwriteOperators,      1,   CDI_REAL,  0,  1 },
  { Change,         ChangeHelp,        ChangeOperators,        1,   CDI_REAL,  1,  1 },
  { Change_e5slm,   NULL,              Change_e5slmOperators,  0,   CDI_REAL,  1,  1 },
  { Cloudlayer,     NULL,              CloudlayerOperators,    1,   CDI_REAL,  1,  1 },
  { CMOR,           CMORHelp,          CMOROperators,          1,   CDI_REAL,  1,  0 },
  { CMOR_lite,      NULL,              CMORliteOperators,      1,   CDI_REAL,  1,  1 },
  { CMOR_table,     NULL,              CMORtableOperators,     1,   CDI_REAL,  0,  0 },
  { Collgrid,       CollgridHelp,      CollgridOperators,      1,   CDI_REAL, -1,  1 },
  { Command,        NULL,              CommandOperators,       0,   CDI_REAL,  1,  0 },
  { Comp,           CompHelp,          CompOperators,          1,   CDI_REAL,  2,  1 },
  { Compc,          CompcHelp,         CompcOperators,         1,   CDI_REAL,  1,  1 },
  { Complextorect,  NULL,              ComplextorectOperators, 1,   CDI_COMP,  1,  2 },
  { Cond,           CondHelp,          CondOperators,          1,   CDI_REAL,  2,  1 },
  { Cond2,          Cond2Help,         Cond2Operators,         1,   CDI_REAL,  3,  1 },
  { Condc,          CondcHelp,         CondcOperators,         1,   CDI_REAL,  1,  1 },
  { Consecstat,     ConsecstatHelp,    ConsecstatOperators,    1,   CDI_REAL,  1,  1 },
  { Copy,           CopyHelp,          CopyOperators,          1,   CDI_REAL, -1,  1 },
  { Deltat,         NULL,              DeltatOperators,        1,   CDI_REAL,  1,  1 },
  { Deltime,        NULL,              DeltimeOperators,       1,   CDI_REAL,  1,  1 },
  { Derivepar,      DeriveparHelp,     DeriveparOperators,     1,   CDI_REAL,  1,  1 },
  { Detrend,        DetrendHelp,       DetrendOperators,       1,   CDI_REAL,  1,  1 },
  { Diff,           DiffHelp,          DiffOperators,          1,   CDI_REAL,  2,  0 },
  { Distgrid,       DistgridHelp,      DistgridOperators,      1,   CDI_REAL,  1,  1 },
  { Duplicate,      DuplicateHelp,     DuplicateOperators,     1,   CDI_REAL,  1,  1 },
  { Echam5ini,      NULL,              Echam5iniOperators,     1,   CDI_REAL,  1,  1 },
  { Enlarge,        EnlargeHelp,       EnlargeOperators,       1,   CDI_REAL,  1,  1 },
  { Enlargegrid,    NULL,              EnlargegridOperators,   0,   CDI_REAL,  1,  1 },
  { Ensstat,        EnsstatHelp,       EnsstatOperators,       1,   CDI_REAL, -1,  1 },
  { Ensstat3,       Ensstat2Help,      Ensstat3Operators,      1,   CDI_REAL, -1,  1 },
  { Ensval,         EnsvalHelp,        EnsvalOperators,        1,   CDI_REAL, -1,  1 },
  { Eofcoeff,       EofcoeffHelp,      EofcoeffOperators,      1,   CDI_REAL,  2, -1 },
  { Eofcoeff3d,     EofcoeffHelp,      Eofcoeff3dOperators,    1,   CDI_REAL,  2, -1 },
  { EOFs,           EOFsHelp,          EOFsOperators,          1,   CDI_REAL,  1,  2 },
  { EOF3d,          EOFsHelp,          EOF3dOperators,         1,   CDI_REAL,  1,  2 },
  { Expr,           ExprHelp,          ExprOperators,          1,   CDI_REAL,  1,  1 },
  { FC,             NULL,              FCOperators,            1,   CDI_REAL,  1,  1 },
  { Filedes,        FiledesHelp,       FiledesOperators,       1,   CDI_BOTH,  1,  0 },
  { Fillmiss,       NULL,              FillmissOperators,      1,   CDI_REAL,  1,  1 },
  { Filter,         FilterHelp,        FilterOperators,        1,   CDI_REAL,  1,  1 },
  { Fldrms,         NULL,              FldrmsOperators,        1,   CDI_REAL,  2,  1 },
  { Fldstat,        FldstatHelp,       FldstatOperators,       1,   CDI_REAL,  1,  1 },
  { Fldstat2,       FldcorHelp,        FldcorOperators,        1,   CDI_REAL,  2,  1 },
  { Fldstat2,       FldcovarHelp,      FldcovarOperators,      1,   CDI_REAL,  2,  1 },
  { Fourier,        NULL,              FourierOperators,       1,   CDI_COMP,  1,  1 },
  { Gengrid,        NULL,              GengridOperators,       1,   CDI_REAL,  2,  1 },
  { Gradsdes,       GradsdesHelp,      GradsdesOperators,      1,   CDI_REAL,  1,  0 },
  { Gridboxstat,    GridboxstatHelp,   GridboxstatOperators,   1,   CDI_REAL,  1,  1 },
  { Gridcell,       GridcellHelp,      GridcellOperators,      1,   CDI_REAL,  1,  1 },
  { Gridsearch,     NULL,              GridsearchOperators,    0,   CDI_REAL,  0,  0 },
  { Harmonic,       NULL,              HarmonicOperators,      1,   CDI_REAL,  1,  1 },
  { Histogram,      HistogramHelp,     HistogramOperators,     1,   CDI_REAL,  1,  1 },
  { Importamsr,     ImportamsrHelp,    ImportamsrOperators,    1,   CDI_REAL,  1,  1 },
  { Importbinary,   ImportbinaryHelp,  ImportbinaryOperators,  1,   CDI_REAL,  1,  1 },
  { Importcmsaf,    ImportcmsafHelp,   ImportcmsafOperators,   1,   CDI_REAL,  1,  1 },
  { Importobs,      NULL,              ImportobsOperators,     1,   CDI_REAL,  1,  1 },
  { Info,           InfoHelp,          InfoOperators,          1,   CDI_BOTH, -1,  0 },
  { Input,          InputHelp,         InputOperators,         1,   CDI_REAL,  0,  1 },
  { Intgrid,        NULL,              IntgridOperators,       1,   CDI_REAL,  1,  1 },
  { Intgridtraj,    NULL,              IntgridtrajOperators,   1,   CDI_REAL,  1,  1 },
  { Intlevel,       IntlevelHelp,      IntlevelOperators,      1,   CDI_REAL,  1,  1 },
  { Intlevel3d,     Intlevel3dHelp,    Intlevel3dOperators,    1,   CDI_REAL,  2,  1 },
  { Inttime,        InttimeHelp,       InttimeOperators,       1,   CDI_REAL,  1,  1 },
  { Intntime,       InttimeHelp,       IntntimeOperators,      1,   CDI_REAL,  1,  1 },
  { Intyear,        IntyearHelp,       IntyearOperators,       1,   CDI_REAL,  2, -1 },
  { Invert,         InvertHelp,        InvertOperators,        1,   CDI_REAL,  1,  1 },
  { Invertlev,      InvertlevHelp,     InvertlevOperators,     1,   CDI_REAL,  1,  1 },
  { Isosurface,     NULL,              IsosurfaceOperators,    1,   CDI_REAL,  1,  1 },
  { Log,            NULL,              LogOperators,           0,   CDI_REAL,  1,  0 },
  { MapReduce,      MapReduceHelp,     MapReduceOperators,     1,   CDI_REAL,  1,  1 },
  { Maskbox,        MaskboxHelp,       MaskboxOperators,       1,   CDI_REAL,  1,  1 },
  { Maskbox,        MaskregionHelp,    MaskregionOperators,    1,   CDI_REAL,  1,  1 },
  { Mastrfu,        MastrfuHelp,       MastrfuOperators,       1,   CDI_REAL,  1,  1 },
  { Math,           MathHelp,          MathOperators,          1,   CDI_BOTH,  1,  1 },
  { Merge,          MergeHelp,         MergeOperators,         1,   CDI_REAL, -1,  1 },
  { Mergetime,      MergeHelp,         MergetimeOperators,     1,   CDI_REAL, -1,  1 },
  { Mergegrid,      MergegridHelp,     MergegridOperators,     1,   CDI_REAL,  2,  1 },
  { Merstat,        MerstatHelp,       MerstatOperators,       1,   CDI_REAL,  1,  1 },
  { Monarith,       MonarithHelp,      MonarithOperators,      1,   CDI_REAL,  2,  1 },
  { Mrotuv,         NULL,              MrotuvOperators,        1,   CDI_REAL,  1,  2 },
  { Mrotuvb,        NULL,              MrotuvbOperators,       1,   CDI_REAL,  2,  1 },
  { Ninfo,          NinfoHelp,         NinfoOperators,         1,   CDI_BOTH,  1,  0 },
  { Nmldump,        NULL,              NmldumpOperators,       0,   CDI_REAL,  0,  0 },
  { Output,         OutputHelp,        OutputOperators,        1,   CDI_REAL, -1,  0 },
  { Output,         OutputtabHelp,     OutputtabOperators,     1,   CDI_REAL, -1,  0 },
  { Outputgmt,      OutputgmtHelp,     OutputgmtOperators,     1,   CDI_REAL,  1,  0 },
  { Pack,           NULL,              PackOperators,          1,   CDI_REAL,  1,  1 },
  { Pardup,         NULL,              PardupOperators,        1,   CDI_REAL,  1,  1 },
  { Pinfo,          NULL,              PinfoOperators,         1,   CDI_REAL,  1,  1 },
  { Pressure,       NULL,              PressureOperators,      1,   CDI_REAL,  1,  1 },
  { Regres,         RegresHelp,        RegresOperators,        1,   CDI_REAL,  1,  1 },
  { Remap,          RemapHelp,         RemapOperators,         1,   CDI_REAL,  1,  1 },
  { Remap,          RemapbilHelp,      RemapbilOperators,      1,   CDI_REAL,  1,  1 },
  { Remap,          RemapbicHelp,      RemapbicOperators,      1,   CDI_REAL,  1,  1 },
  { Remap,          RemapnnHelp,       RemapnnOperators,       1,   CDI_REAL,  1,  1 },
  { Remap,          RemapdisHelp,      RemapdisOperators,      1,   CDI_REAL,  1,  1 },
  { Remap,          RemapyconHelp,     RemapyconOperators,     1,   CDI_REAL,  1,  1 },
  { Remap,          RemapconHelp,      RemapconOperators,      1,   CDI_REAL,  1,  1 },
  { Remap,          Remapcon2Help,     Remapcon2Operators,     1,   CDI_REAL,  1,  1 },
  { Remap,          RemaplafHelp,      RemaplafOperators,      1,   CDI_REAL,  1,  1 },
  { Remap,          NULL,              RemapgridOperators,     1,   CDI_REAL,  1,  1 },
  { Remapeta,       RemapetaHelp,      RemapetaOperators,      1,   CDI_REAL,  1,  1 },
  { Replace,        ReplaceHelp,       ReplaceOperators,       1,   CDI_REAL,  2,  1 },
  { Replacevalues,  ReplacevaluesHelp, ReplacevaluesOperators, 1,   CDI_REAL,  1,  1 },
  { Rhopot,         RhopotHelp,        RhopotOperators,        1,   CDI_REAL,  1,  1 },
  { Rotuv,          RotuvbHelp,        RotuvOperators,         1,   CDI_REAL,  1,  1 },
  { Runpctl,        RunpctlHelp,       RunpctlOperators,       1,   CDI_REAL,  1,  1 },
  { Runstat,        RunstatHelp,       RunstatOperators,       1,   CDI_REAL,  1,  1 },
  { Seascount,      NULL,              SeascountOperators,     1,   CDI_BOTH,  1,  1 },
  { Seaspctl,       SeaspctlHelp,      SeaspctlOperators,      1,   CDI_REAL,  3,  1 },
  { Seasstat,       SeasstatHelp,      SeasstatOperators,      1,   CDI_REAL,  1,  1 },
  { Selbox,         SelboxHelp,        SelboxOperators,        1,   CDI_BOTH,  1,  1 },
  { Selindex,       NULL,              SelindexOperators,      1,   CDI_BOTH,  1,  1 },
  { Select,         SelectHelp,        SelectOperators,        1,   CDI_BOTH, -1,  1 },
  { Selvar,         SelvarHelp,        SelvarOperators,        1,   CDI_BOTH,  1,  1 },
  { Selrec,         SelvarHelp,        SelrecOperators,        1,   CDI_BOTH,  1,  1 },
  { Seloperator,    NULL,              SeloperatorOperators,   1,   CDI_REAL,  1,  1 },
  { Seltime,        SeltimeHelp,       SeltimeOperators,       1,   CDI_BOTH,  1,  1 },
  { Set,            SetHelp,           SetOperators,           1,   CDI_BOTH,  1,  1 },
  { Setbox,         SetboxHelp,        SetboxOperators,        1,   CDI_REAL,  1,  1 },
  { Setgatt,        SetgattHelp,       SetgattOperators,       1,   CDI_BOTH,  1,  1 },
  { Setgrid,        SetgridHelp,       SetgridOperators,       1,   CDI_BOTH,  1,  1 },
  { Sethalo,        SethaloHelp,       SethaloOperators,       1,   CDI_REAL,  1,  1 },
  { Setmiss,        SetmissHelp,       SetmissOperators,       1,   CDI_REAL,  1,  1 },
  { Fillmiss,       SetmissHelp,       SetmisstonnOperators,   1,   CDI_REAL,  1,  1 },
  { Setpartab,      SetHelp,           SetcodetabOperators,    1,   CDI_REAL,  1,  1 },
  { Setpartab,      SetpartabHelp,     SetpartabOperators,     1,   CDI_REAL,  1,  1 },
  { Setrcaname,     NULL,              SetrcanameOperators,    1,   CDI_REAL,  1,  1 },
  { Settime,        SettimeHelp,       SettimeOperators,       1,   CDI_BOTH,  1,  1 },
  { Setzaxis,       SetzaxisHelp,      SetzaxisOperators,      1,   CDI_BOTH,  1,  1 },
  { Shiftxy,        NULL,              ShiftxyOperators,       1,   CDI_REAL,  1,  1 },
  { Showinfo,       ShowinfoHelp,      ShowinfoOperators,      1,   CDI_BOTH,  1,  0 },
  { Sinfo,          SinfoHelp,         SinfoOperators,         1,   CDI_BOTH, -1,  0 },
  { Smooth,         SmoothHelp,        SmoothOperators,        1,   CDI_REAL,  1,  1 },
  { Sort,           NULL,              SortOperators,          1,   CDI_REAL,  1,  1 },
  { Sorttimestamp,  NULL,              SorttimestampOperators, 1,   CDI_REAL, -1,  1 },
  { Specinfo,       NULL,              SpecinfoOperators,      1,   CDI_REAL,  0,  0 },
  { Spectral,       SpectralHelp,      SpectralOperators,      1,   CDI_REAL,  1,  1 },
  { Spectrum,       NULL,              SpectrumOperators,      1,   CDI_REAL,  1,  1 },
  { Split,          SplitHelp,         SplitOperators,         1,   CDI_BOTH,  1, -1 },
  { Splitrec,       SplitHelp,         SplitrecOperators,      1,   CDI_BOTH,  1, -1 },
  { Splitsel,       SplitselHelp,      SplitselOperators,      1,   CDI_BOTH,  1, -1 },
  { Splittime,      SplittimeHelp,     SplittimeOperators,     1,   CDI_BOTH,  1, -1 },
  { Splityear,      SplittimeHelp,     SplityearOperators,     1,   CDI_BOTH,  1, -1 },
  { Subtrend,       SubtrendHelp,      SubtrendOperators,      1,   CDI_REAL,  3,  1 },
  { Tee,            NULL,              TeeOperators,           1,   CDI_REAL,  2,  1 },
  { Template1,      NULL,              Template1Operators,     0,   CDI_REAL,  1,  1 },
  { Template2,      NULL,              Template2Operators,     0,   CDI_REAL,  1,  1 },
  { Test,           NULL,              TestOperators,          0,   CDI_REAL,  1,  1 },
  { Test2,          NULL,              Test2Operators,         0,   CDI_REAL,  2,  1 },
  { Testdata,       NULL,              TestdataOperators,      0,   CDI_REAL,  1,  1 },
  { Tests,          NULL,              TestsOperators,         0,   CDI_REAL,  1,  1 },
  { Timcount,       NULL,              TimcountOperators,      1,   CDI_BOTH,  1,  1 },
  { Timcount,       NULL,              YearcountOperators,     1,   CDI_BOTH,  1,  1 },
  { Timcount,       NULL,              MoncountOperators,      1,   CDI_BOTH,  1,  1 },
  { Timcount,       NULL,              DaycountOperators,      1,   CDI_BOTH,  1,  1 },
  { Timcount,       NULL,              HourcountOperators,     1,   CDI_BOTH,  1,  1 },
  { Timcumsum,      NULL,              TimcumsumOperators,     1,   CDI_BOTH,  1,  1 },
  { Timpctl,        TimpctlHelp,       TimpctlOperators,       1,   CDI_REAL,  3,  1 },
  { Timpctl,        YearpctlHelp,      YearpctlOperators,      1,   CDI_REAL,  3,  1 },
  { Timpctl,        MonpctlHelp,       MonpctlOperators,       1,   CDI_REAL,  3,  1 },
  { Timpctl,        DaypctlHelp,       DaypctlOperators,       1,   CDI_REAL,  3,  1 },
  { Timpctl,        HourpctlHelp,      HourpctlOperators,      1,   CDI_REAL,  3,  1 },
  { Timselpctl,     TimselpctlHelp,    TimselpctlOperators,    1,   CDI_REAL,  3,  1 },
  { Timsort,        TimsortHelp,       TimsortOperators,       1,   CDI_REAL,  1,  1 },
  { Timselstat,     TimselstatHelp,    TimselstatOperators,    1,   CDI_REAL,  1,  1 },
  { XTimstat,       NULL,              XTimstatOperators,      0,   CDI_BOTH,  1,  1 },
  { Timstat,        TimstatHelp,       TimstatOperators,       1,   CDI_BOTH,  1,  1 },
  { Timstat,        YearstatHelp,      YearstatOperators,      1,   CDI_BOTH,  1,  1 },
  { Timstat,        MonstatHelp,       MonstatOperators,       1,   CDI_BOTH,  1,  1 },
  { Timstat,        DaystatHelp,       DaystatOperators,       1,   CDI_BOTH,  1,  1 },
  { Timstat,        HourstatHelp,      HourstatOperators,      1,   CDI_BOTH,  1,  1 },
  { Timstat2,       TimcorHelp,        TimcorOperators,        1,   CDI_REAL,  2,  1 },
  { Timstat2,       TimcovarHelp,      TimcovarOperators,      1,   CDI_REAL,  2,  1 },
  { Timstat3,       NULL,              Timstat3Operators,      1,   CDI_REAL,  2,  1 },
  { Tinfo,          NULL,              TinfoOperators,         1,   CDI_BOTH,  1,  0 },
  { Tocomplex,      NULL,              TocomplexOperators,     1,   CDI_REAL,  1,  1 },
  { Transpose,      NULL,              TransposeOperators,     1,   CDI_REAL,  1,  1 },
  { Trend,          TrendHelp,         TrendOperators,         1,   CDI_REAL,  1,  2 },
  { Trms,           NULL,              TrmsOperators,          0,   CDI_REAL,  2,  1 },
  { Tstepcount,     NULL,              TstepcountOperators,    1,   CDI_REAL,  1,  1 },
  { Vargen,         VargenHelp,        VargenOperators,        1,   CDI_REAL,  0,  1 },
  { Varrms,         NULL,              VarrmsOperators,        0,   CDI_REAL,  2,  1 },
  { Vertintml,      VertintmlHelp,     VertintmlOperators,     1,   CDI_REAL,  1,  1 },
  { Vertintap,      VertintapHelp,     VertintapOperators,     1,   CDI_REAL,  1,  1 },
  { Vertstat,       VertstatHelp,      VertstatOperators,      1,   CDI_REAL,  1,  1 },
  { Vertcum,        NULL,              VertcumOperators,       1,   CDI_REAL,  1,  1 },
  { Vertwind,       NULL,              VertwindOperators,      1,   CDI_REAL,  1,  1 },
  { Verifygrid,     NULL,              VerifygridOperators,    1,   CDI_REAL,  1,  0 },
  { Wind,           WindHelp,          WindOperators,          1,   CDI_REAL,  1,  1 },
  { Writegrid,      NULL,              WritegridOperators,     1,   CDI_REAL,  1,  1 },  /* no cdi output */
  { Writerandom,    NULL,              WriterandomOperators,   1,   CDI_REAL,  1,  1 },
  { YAR,            NULL,              YAROperators,           0,   CDI_REAL,  1,  1 },
  { Yearmonstat,    YearmonstatHelp,   YearmonstatOperators,   1,   CDI_REAL,  1,  1 },
  { Ydayarith,      YdayarithHelp,     YdayarithOperators,     1,   CDI_REAL,  2,  1 },
  { Ydaypctl,       YdaypctlHelp,      YdaypctlOperators,      1,   CDI_REAL,  3,  1 },
  { Ydaystat,       YdaystatHelp,      YdaystatOperators,      1,   CDI_REAL,  1,  1 },
  { Ydrunpctl,      YdrunpctlHelp,     YdrunpctlOperators,     1,   CDI_REAL,  3,  1 },
  { Ydrunstat,      YdrunstatHelp,     YdrunstatOperators,     1,   CDI_REAL,  1,  1 },
  { Yhourarith,     YhourarithHelp,    YhourarithOperators,    1,   CDI_REAL,  2,  1 },
  { Yhourstat,      YhourstatHelp,     YhourstatOperators,     1,   CDI_REAL,  1,  1 },
  { Ymonarith,      YmonarithHelp,     YmonarithOperators,     1,   CDI_REAL,  2,  1 },
  { Ymonarith,      YseasarithHelp,    YseasarithOperators,    1,   CDI_REAL,  2,  1 },
  { Ymonpctl,       YmonpctlHelp,      YmonpctlOperators,      1,   CDI_REAL,  3,  1 },
  { Ymonstat,       YmonstatHelp,      YmonstatOperators,      1,   CDI_REAL,  1,  1 },
  { Yseaspctl,      YseaspctlHelp,     YseaspctlOperators,     1,   CDI_REAL,  3,  1 },
  { Yseasstat,      YseasstatHelp,     YseasstatOperators,     1,   CDI_REAL,  1,  1 },
  { Zonstat,        ZonstatHelp,       ZonstatOperators,       1,   CDI_REAL,  1,  1 },
  { EcaCfd,         EcaCfdHelp,        EcaCfdOperators,        1,   CDI_REAL,  1,  1 },
  { EcaCsu,         EcaCsuHelp,        EcaCsuOperators,        1,   CDI_REAL,  1,  1 },
  { EcaCwdi,        EcaCwdiHelp,       EcaCwdiOperators,       1,   CDI_REAL,  2,  1 },
  { EcaCwfi,        EcaCwfiHelp,       EcaCwfiOperators,       1,   CDI_REAL,  2,  1 },
  { EcaEtr,         EcaEtrHelp,        EcaEtrOperators,        1,   CDI_REAL,  2,  1 },
  { EcaFd,          EcaFdHelp,         EcaFdOperators,         1,   CDI_REAL,  1,  1 },
  { EcaGsl,         EcaGslHelp,        EcaGslOperators,        1,   CDI_REAL,  2,  1 },
  { EcaHd,          EcaHdHelp,         EcaHdOperators,         1,   CDI_REAL,  1,  1 },
  { EcaHwdi,        EcaHwdiHelp,       EcaHwdiOperators,       1,   CDI_REAL,  2,  1 },
  { EcaHwfi,        EcaHwfiHelp,       EcaHwfiOperators,       1,   CDI_REAL,  2,  1 },
  { EcaId,          EcaIdHelp,         EcaIdOperators,         1,   CDI_REAL,  1,  1 },
  { EcaSu,          EcaSuHelp,         EcaSuOperators,         1,   CDI_REAL,  1,  1 },
  { EcaTr,          EcaTrHelp,         EcaTrOperators,         1,   CDI_REAL,  1,  1 },
  { EcaTg10p,       EcaTg10pHelp,      EcaTg10pOperators,      1,   CDI_REAL,  2,  1 },
  { EcaTg90p,       EcaTg90pHelp,      EcaTg90pOperators,      1,   CDI_REAL,  2,  1 },
  { EcaTn10p,       EcaTn10pHelp,      EcaTn10pOperators,      1,   CDI_REAL,  2,  1 },
  { EcaTn90p,       EcaTn90pHelp,      EcaTn90pOperators,      1,   CDI_REAL,  2,  1 },
  { EcaTx10p,       EcaTx10pHelp,      EcaTx10pOperators,      1,   CDI_REAL,  2,  1 },
  { EcaTx90p,       EcaTx90pHelp,      EcaTx90pOperators,      1,   CDI_REAL,  2,  1 },
  { EcaCdd,         EcaCddHelp,        EcaCddOperators,        1,   CDI_REAL,  1,  1 },
  { EcaCwd,         EcaCwdHelp,        EcaCwdOperators,        1,   CDI_REAL,  1,  1 },
  { EcaRr1,         EcaRr1Help,        EcaRr1Operators,        1,   CDI_REAL,  1,  1 },
  { EcaPd,          EcaPdHelp,         EcaPdOperators,         1,   CDI_REAL,  1,  1 },
  { EcaR75p,        EcaR75pHelp,       EcaR75pOperators,       1,   CDI_REAL,  2,  1 },
  { EcaR75ptot,     EcaR75ptotHelp,    EcaR75ptotOperators,    1,   CDI_REAL,  2,  1 },
  { EcaR90p,        EcaR90pHelp,       EcaR90pOperators,       1,   CDI_REAL,  2,  1 },
  { EcaR90ptot,     EcaR90ptotHelp,    EcaR90ptotOperators,    1,   CDI_REAL,  2,  1 },
  { EcaR95p,        EcaR95pHelp,       EcaR95pOperators,       1,   CDI_REAL,  2,  1 },
  { EcaR95ptot,     EcaR95ptotHelp,    EcaR95ptotOperators,    1,   CDI_REAL,  2,  1 },
  { EcaR99p,        EcaR99pHelp,       EcaR99pOperators,       1,   CDI_REAL,  2,  1 },
  { EcaR99ptot,     EcaR99ptotHelp,    EcaR99ptotOperators,    1,   CDI_REAL,  2,  1 },
  { EcaRx1day,      EcaRx1dayHelp,     EcaRx1dayOperators,     1,   CDI_REAL,  1,  1 },
  { EcaRx5day,      EcaRx5dayHelp,     EcaRx5dayOperators,     1,   CDI_REAL,  1,  1 },
  { EcaSdii,        EcaSdiiHelp,       EcaSdiiOperators,       1,   CDI_REAL,  1,  1 },
  { Fdns,           FdnsHelp,          FdnsOperators,          1,   CDI_REAL,  2,  1 },
  { Strwin,         StrwinHelp,        StrwinOperators,        1,   CDI_REAL,  1,  1 },
  { Strbre,         StrbreHelp,        StrbreOperators,        1,   CDI_REAL,  1,  1 },
  { Strgal,         StrgalHelp,        StrgalOperators,        1,   CDI_REAL,  1,  1 },
  { Hurr,           HurrHelp,          HurrOperators,          1,   CDI_REAL,  1,  1 },
  /*  { Hi,             NULL,              HiOperators,        1,   CDI_REAL,  3,  1 }, */
  { Wct,            WctHelp,           WctOperators,           1,   CDI_REAL,  2,  1 },
  { Magplot,        MagplotHelp,       MagplotOperators,       1,   CDI_REAL,  1,  1 },
  { Magvector,      MagvectorHelp,     MagvectorOperators,     1,   CDI_REAL,  1,  1 },
  { Maggraph,       MaggraphHelp,      MaggraphOperators,      1,   CDI_REAL, -1,  1 },
};							       
							       
static int NumModules = sizeof(Modules) / sizeof(Modules[0]);

static const char *opalias[][2] =
{
  {"afterburner",         "after"},
  {"anomaly",             "ymonsub"},
  {"deltap_fl",           "deltap"},
  {"diffv",               "diffn"},
  {"covar0",              "timcovar"},
  {"covar0r",             "fldcovar"},
  {"gather",              "collgrid"},
  {"geopotheight",        "gheight"},
  {"ggstat",              "info"},
  {"ggstats",             "sinfo"},
  {"globavg",             "fldavg"},
  {"import_grads",        "import_binary"},
  {"infos",               "sinfo"},
  {"infov",               "infon"},
  {"intgrid",             "intgridbil"},
  {"log",                 "ln"},
  {"lmean",               "ymonmean"},
  {"lmmean",              "ymonmean"},
  {"lmavg",               "ymonavg"},
  {"lmstd",               "ymonstd"},
  {"lsmean",              "yseasmean"},
  {"chvar",               "chname"},
  {"ncode",               "npar"},
  {"nvar",                "npar"},
  {"outputkey",           "outputtab"},
  {"vardes",              "codetab"},
  {"pardes",              "codetab"},
  {"selvar",              "selname"},
  {"delvar",              "delname"},
  {"remapcon1",           "remaplaf"},
  {"remapdis1",           "remapnn"},
  {"scatter",             "distgrid"},
  {"showvar",             "showname"},
  {"selgridname",         "selgrid"},
  {"setvar",              "setname"},
  {"setpartabv",          "setpartabn"},
  {"setpartab",           "setcodetab"},
  {"sinfov",              "sinfon"},
  {"sortvar",             "sortname"},
  {"splitvar",            "splitname"},
  {"sort",                "timsort"},
  {"eca_r1mm",            "eca_rr1"},
  {"fpressure",           "pressure_fl"},
  {"hpressure",           "pressure_hl"},
  {"ensrkhist_space",     "ensrkhistspace"},
  {"ensrkhist_time",      "ensrkhisttime"},
  {"gridverify",          "verifygrid"},
  {"outputcenter",        "gmtxyz"},
  {"outputbounds",        "gmtcells"},
  {"selseas",             "selseason"},
  {"selmon",              "selmonth"},
};

static int nopalias = sizeof(opalias) / (2*sizeof(opalias[0][0]));


static
bool similar(const char *a, const char *b, int alen, int blen)
{
  if ( alen > 2 && blen > 2 && strstr(b, a) )
    return true;

  while ( *a && *b && *a == *b )
    { 
      a++;
      b++;
    }
  if ( !*a && !*b )
    return true;
  /*
    printf("%d %d %s %s\n", alen, blen, a, b);
  */
  if ( alen >= 2 && blen >= 1 && *a && similar(a+1, b, alen-2, blen-1) )
    return true;

  if ( alen >= 1 && blen >= 2 && *b && similar(a, b+1, alen-1, blen-2) )
    return true;

  return false; 
}


const char *operatorAlias(const char *operatorName)
{
  int i;

  const char *operatorNameNew = operatorName;

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

  return operatorNameNew;
}

static
int operatorInqModID(const char *operatorName)
{
  int i, j, modID = -1;

  if ( operatorName )
    {
      for ( i = 0; i < NumModules; i++ )
	{
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

  if ( modID == -1 && *operatorName == 0 )
    Error("Operator name missing!");

  if ( modID == -1 )
    {
      bool lfound = false;
      FILE *fp = fopen(operatorName, "r");
      if ( fp )
	{
	  fclose(fp);
	  fprintf(stderr, "Use commandline option -h for help.");
	  Error("Operator missing, %s is a file on disk!", operatorName);
	}

      fprintf(stderr, "Operator >%s< not found!\n", operatorName);
      fprintf(stderr, "Similar operators are:\n");
      int nbyte = fprintf(stderr, "   ");
      if ( operatorName )
        {
          int len = strlen(operatorName);
          char *opname = strdup(operatorName);
          strtolower(opname);
          for ( i = 0; i < NumModules; i++ )
            {
              if ( Modules[i].help == NULL ) continue;
              j = 0;
              while ( Modules[i].operators[j] )
                {
                  if ( similar(opname, Modules[i].operators[j],
                               len, strlen(Modules[i].operators[j])) )
                    {
                      if ( nbyte > 75 )
                        {
                          fprintf(stdout, "\n");
                          nbyte = fprintf(stderr, "   ");
                        }
                      nbyte += fprintf(stderr, " %s", Modules[i].operators[j]);
                      lfound = true;
                    }
                  j++;
                }
            }
          free(opname);
        }

      if ( !lfound )
	fprintf(stderr, "(not found)\n") ;
      else
	fprintf(stderr, "\n");

      exit(EXIT_FAILURE);
    }

  if ( modID != -1 )
    if ( ! Modules[modID].func )
      Error("Module for operator >%s< not installed!", operatorName);

  return modID;
}

void *(*operatorModule(const char *operatorName))(void *)
{
  int modID = operatorInqModID(operatorName);
  return Modules[modID].func;
}

const char **operatorHelp(const char *operatorName)
{
  int modID = operatorInqModID(operatorName);
  return Modules[modID].help;
}

int operatorStreamInCnt(const char *operatorName)
{
  int modID = operatorInqModID(operatorAlias(operatorName));
  return Modules[modID].streamInCnt;
}

int operatorStreamOutCnt(const char *operatorName)
{
  int modID = operatorInqModID(operatorAlias(operatorName));
  return Modules[modID].streamOutCnt;
}

int operatorStreamNumber(const char *operatorName)
{
  int modID = operatorInqModID(operatorAlias(operatorName));
  return Modules[modID].number;
}

int cmpname(const void *s1, const void *s2)
{
  char **c1 = (char **) s1;
  char **c2 = (char **) s2;

  return strcmp((const char *)*c1, (const char *)*c2);
}

void operatorPrintAll(void)
{
  int i, j, nbyte, nop = 0;
  const char *opernames[4096];
  FILE *pout = stderr;

  for ( i = 0; i < NumModules; i++ )
    {
      if ( Modules[i].mode )
        {
          j = 0;
          while ( Modules[i].operators[j] )
            opernames[nop++] = Modules[i].operators[j++];
	}
    }

  // Add operator aliases
  for ( i = 0; i < nopalias; i++ )
    {
      opernames[nop++] = opalias[i][0];
    }

  qsort(opernames, nop, sizeof(char *), cmpname);

  nbyte = fprintf(pout, "   ");
  for ( i = 0; i < nop; i++ )
    {
      if ( nbyte > 85 )
	{
	  fprintf(pout, "\n");
	  nbyte = fprintf(pout, "   ");
	}
      nbyte += fprintf(pout, " %s", opernames[i]);
    }
  fprintf(pout, "\n");
}


void operatorPrintList(bool print_no_output)
{
  int i, j, nbyte, nop = 0;
  const char *opernames[4096];
  FILE *pout = stdout;

  for ( i = 0; i < NumModules; i++ )
    {
      if ( Modules[i].mode )
        {
          j = 0;
          while ( Modules[i].operators[j] )
            opernames[nop++] = Modules[i].operators[j++];
        }
    }

  // Add operator aliases
  for ( i = 0; i < nopalias; i++ ) opernames[nop++] = opalias[i][0];

  qsort(opernames, nop, sizeof(char *), cmpname);

  for ( i = 0; i < nop; i++ )
    {
      int iname = -1;
      int ialias = -1;
      size_t operlen = strlen(opernames[i]);
      const char *pdes = NULL;
      const char **phelp = NULL;
      const char **help = NULL;

      for ( j = 0; j < nopalias; j++ )
        {
          if ( strcmp(opernames[i], opalias[j][0]) == 0 )
            {
              ialias = j;
              break;
            }
        }

      if ( ialias == -1 ) help = operatorHelp((char *)opernames[i]);

      phelp = help;
      if ( phelp )
        {
          int loper = FALSE;
          int index = 0;
          while ( *phelp )
            {
              if ( loper )
                {
                  if ( *phelp[0] == '\0' ) break;
                  const char *ph = *phelp;
                  while ( *ph != '\0' )
                    {
                      if ( *ph != ' ' ) break;
                      ph++;
                    }
                  if ( *ph != '\0' && strncmp(ph, opernames[i], operlen) == 0 && ph[operlen] == ' ' )
                    {
                      ph += operlen;
                      while ( *ph != '\0' )
                        {
                          if ( *ph != ' ' ) break;
                          ph++;
                        }
                      if ( *ph != '\0' ) pdes = ph;
                    }
                }
              
              if ( strcmp(*phelp, "NAME") == 0 ) iname = index+1;
              if ( strcmp(*phelp, "OPERATORS") == 0 ) loper = TRUE;
              index++;
              phelp++;
            }
        }

      phelp = help;
      if ( phelp && pdes == NULL && iname >= 0 )
        {
          const char *ph = phelp[iname];
          while ( *ph != '\0' )
            {
              if ( *ph != ' ' ) break;
              ph++;
            }
          if ( *ph != '\0' && strncmp(ph, opernames[i], operlen) == 0 && ph[operlen] == ' ' )
            {
              ph += operlen;
              while ( *ph != '\0' )
                {
                  if ( *ph != ' ' && *ph != '-' ) break;
                  ph++;
                }
              if ( *ph != '\0' ) pdes = ph;
            }      
        }

      if ( !print_no_output || operatorStreamOutCnt(opernames[i]) == 0 )
        {
          nbyte = fprintf(pout, "%s ", opernames[i]);
          for ( int i = nbyte; i <= 16; ++i ) fprintf(pout, " ");
          if ( pdes ) fprintf(pout, "%s", pdes);
          else if ( ialias >= 0 )  fprintf(pout, "--> %s", opalias[ialias][1]);
          fprintf(pout, "\n");
        }
    }
}

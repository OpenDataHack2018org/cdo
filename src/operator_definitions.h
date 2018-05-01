/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2018 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/
#ifndef OPERATOR_DEFINITIONS_H
#define OPERATOR_DEFINITIONS_H

/* clang-format off */
#define  AdisitOperators        {"adisit", "adipot"}
#define  AfterburnerOperators   {"after"}
#define  ArithOperators         {"add",  "sub",  "mul",  "div", "min", "max", "atan2", "setmiss"}
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
#define  ComplextorectOperators {"complextorect", "complextopol"}
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
#define  Echam5iniOperators     {"import_e5ml", "import_e5res", "export_e5ml", "export_e5res"}
#define  EnlargeOperators       {"enlarge"}
#define  EnlargegridOperators   {"enlargegrid"}
#define  EnsstatOperators       {"ensrange", "ensmin", "ensmax", "enssum", "ensmean", "ensavg", "ensvar", "ensvar1", "ensstd", "ensstd1", "enspctl", "ensskew", "enskurt"}
#define  Ensstat3Operators      {"ensrkhistspace", "ensrkhisttime", "ensroc"}
#define  EnsvalOperators        {"enscrps", "ensbrs"}
#define  EofcoeffOperators      {"eofcoeff"}
#define  Eofcoeff3dOperators    {"eofcoeff3d"}
#define  EOFsOperators          {"eof", "eofspatial", "eoftime"}
#define  EOF3dOperators         {"eof3d","eof3dspatial","eof3dtime"}
#define  EstFreqOperators       {"estfreq"}
#define  ExprOperators          {"expr", "exprf", "aexpr", "aexprf"}
#define  FCOperators            {"fc2sp", "sp2fc", "fc2gp", "gp2fc"}
#define  FiledesOperators       {"filedes", "griddes", "griddes2", "zaxisdes", "vct", "vct2", "codetab", \
                                 "vlist", "partab", "partab2", "spartab"}
#define  FillmissOperators      {"fillmiss", "fillmiss2"}
#define  FilterOperators        {"bandpass", "highpass", "lowpass"}
#define  FldrmsOperators        {"fldrms"}
#define  FldstatOperators       {"fldrange", "fldmin", "fldmax", "fldsum", "fldmean", "fldavg", "fldstd", "fldstd1", "fldvar", "fldvar1", "fldpctl"}
#define  FldcorOperators        {"fldcor"}
#define  FldcovarOperators      {"fldcovar"}
#define  FourierOperators       {"fourier"}
#define  GengridOperators       {"gengrid"}
#define  GradsdesOperators      {"gradsdes", "dumpmap"}
#define  GridboxstatOperators   {"gridboxrange", "gridboxmin", "gridboxmax", "gridboxsum", "gridboxmean", "gridboxavg", "gridboxstd", "gridboxstd1", "gridboxvar", "gridboxvar1"}
#define  GridcellOperators      {"gridarea", "gridweights", "gridmask", "griddx", "griddy", "gridcellidx"}
#define  GridsearchOperators    {"testpointsearch", "testcellsearch"}
#define  HarmonicOperators      {"harmonic"}
#define  HistogramOperators     {"histcount", "histsum", "histmean", "histfreq"}
#define  ImportamsrOperators    {"import_amsr"}
#define  ImportbinaryOperators  {"import_binary"}
#define  ImportcmsafOperators   {"import_cmsaf"}
#define  ImportobsOperators     {"import_obs"}
#define  InfoOperators          {"info", "infop", "infon", "infoc", "xinfon", "map"}
#define  InputOperators         {"input", "inputsrv", "inputext"}
#define  IntgridOperators       {"intgridbil", "intgriddis", "intgridnn", "interpolate", "boxavg", "thinout"}
#define  IntgridtrajOperators   {"intgridtraj"}
#define  IntlevelOperators      {"intlevel", "intlevelx"}
#define  Intlevel3dOperators    {"intlevel3d", "intlevelx3d"}
#define  InttimeOperators       {"inttime"}
#define  IntntimeOperators      {"intntime"}
#define  IntyearOperators       {"intyear"}
#define  InvertOperators        {"invertlat", "invertlon", "invertlatdes", "invertlondes", "invertlatdata", "invertlondata"}
#define  InvertlevOperators     {"invertlev"}
#define  IsosurfaceOperators    {"isosurface"}
#define  MapReduceOperators     {"reducegrid"}
#define  MaskboxOperators       {"masklonlatbox", "maskindexbox"}
#define  MaskregionOperators    {"maskregion"}
#define  MastrfuOperators       {"mastrfu"}
#define  MathOperators          {"abs", "int", "nint", "sqr", "sqrt", "exp", "ln", "log10", "sin", "cos", "tan", "asin", "acos", "atan", "pow", "reci", "not", \
                                 "conj", "re", "im", "arg"}
#define  MergeOperators         {"merge"}
#define  MergegridOperators     {"mergegrid"}
#define  MergetimeOperators     {"mergetime"}
#define  MerstatOperators       {"merrange", "mermin", "mermax", "mersum", "mermean", "meravg", "merstd", "merstd1", "mervar", "mervar1", "merpctl"}
#define  MonarithOperators      {"monadd", "monsub", "monmul", "mondiv"}
#define  MrotuvOperators        {"mrotuv"}
#define  MrotuvbOperators       {"mrotuvb"}
#define  NCL_windOperators      {"uv2dv_cfd", "uv2vr_cfd"}
#define  NinfoOperators         {"nyear", "nmon", "ndate", "ntime", "ncode", "npar", "nlevel", "ngridpoints", "ngrids"}
#define  NmldumpOperators       {"nmldump", "kvldump"}
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
#define  RunstatOperators       {"runrange", "runmin", "runmax", "runsum", "runmean", "runavg", "runstd", "runstd1", "runvar", "runvar1"}
#define  SamplegridiconOperators {"samplegridicon"}
#define  SeascountOperators     {"seascount"}
#define  SeaspctlOperators      {"seaspctl"}
#define  SeasstatOperators      {"seasrange", "seasmin", "seasmax", "seassum", "seasmean", "seasavg", "seasstd", "seasstd1", "seasvar", "seasvar1"}
#define  SelboxOperators        {"sellonlatbox", "selindexbox"}
#define  SelgridcellOperators   {"selgridcell", "delgridcell"}
#define  SelectOperators        {"select", "delete"}
#define  SelvarOperators        {"selparam", "selcode", "selname", "selstdname", "sellevel", "sellevidx", "selgrid", \
                                 "selzaxis", "selzaxisname", "seltabnum", "delparam", "delcode", "delname", "selltype"}
#define  SeloperatorOperators   {"seloperator"}
#define  SelrecOperators        {"selrec"}
#define  SeltimeOperators       {"seltimestep", "selyear", "selseason", "selmonth", "selday", "selhour", "seldate", \
                                 "seltime", "selsmon"}
#define  SelyearidxOperators    {"selyearidx", "seltimeidx"}
#define  SetOperators           {"setcode", "setparam", "setname", "setunit", "setlevel", "setltype", "settabnum"}
#define  SetattributeOperators  {"setattribute"}
#define  SetboxOperators        {"setclonlatbox", "setcindexbox"}
#define  SetgattOperators       {"setgatt", "setgatts"}
#define  SetgridOperators       {"setgrid", "setgridtype", "setgridarea", "setgridmask", "unsetgridmask", "setgridnumber", "setgriduri", "usegridnumber"}
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
                                 "showparam", "showname", "showstdname", "showlevel", "showltype", "showformat", "showgrid", "showatts", "showattsglob"}
#define  ShowattributeOperators {"showattribute", "showattsvar"}
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
#define  TimselstatOperators    {"timselrange", "timselmin", "timselmax", "timselsum", "timselmean", "timselavg", "timselvar", "timselvar1", "timselstd", "timselstd1"}
#define  XTimstatOperators      {"xtimmin",  "xtimmax",  "xtimsum",  "xtimmean",  "xtimavg",  "xtimvar",  "xtimvar1",  "xtimstd",  "xtimstd1", \
                                 "xyearmin", "xyearmax", "xyearsum", "xyearmean", "xyearavg", "xyearvar", "xyearvar1", "xyearstd", "xyearstd1", \
                                 "xmonmin",  "xmonmax",  "xmonsum",  "xmonmean",  "xmonavg",  "xmonvar",  "xmonvar1",  "xmonstd",  "xmonstd1"}
#define  TimstatOperators       {"timrange",  "timmin",  "timmax",  "timsum",  "timmean",  "timavg",  "timvar",  "timvar1",  "timstd",  "timstd1"}
#define    YearstatOperators    {"yearrange", "yearmin", "yearmax", "yearsum", "yearmean", "yearavg", "yearvar", "yearvar1", "yearstd", "yearstd1", "yearminidx", "yearmaxidx"}
#define    MonstatOperators     {"monrange",  "monmin",  "monmax",  "monsum",  "monmean",  "monavg",  "monvar",  "monvar1",  "monstd",  "monstd1"}
#define    DaystatOperators     {"dayrange",  "daymin",  "daymax",  "daysum",  "daymean",  "dayavg",  "dayvar",  "dayvar1",  "daystd",  "daystd1"}
#define    HourstatOperators    {"hourrange", "hourmin", "hourmax", "hoursum", "hourmean", "houravg", "hourvar", "hourvar1", "hourstd", "hourstd1"}
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
#define  VertstatOperators      {"vertrange", "vertmin", "vertmax", "vertsum", "vertint", "vertmean", "vertavg", "vertstd", "vertstd1", "vertvar", "vertvar1"}
#define  VertcumOperators       {"vertcum", "vertcumhl"}
#define  VertwindOperators      {"vertwind"}
#define  VerifygridOperators    {"verifygrid"}
#define  WindOperators          {"uv2dv", "uv2dvl", "dv2uv", "dv2uvl", "dv2ps"}
#define  WritegridOperators     {"writegrid"}
#define  WriterandomOperators   {"writerandom"}
#define  YearmonstatOperators   {"yearmonmean", "yearmonavg"}
#define  YdayarithOperators     {"ydayadd", "ydaysub", "ydaymul", "ydaydiv"}
#define  YdaypctlOperators      {"ydaypctl"}
#define  YdaystatOperators      {"ydayrange", "ydaymin", "ydaymax", "ydaysum", "ydaymean", "ydayavg", "ydaystd", "ydaystd1", "ydayvar", "ydayvar1"}
#define  YdrunpctlOperators     {"ydrunpctl"}
#define  YdrunstatOperators     {"ydrunmin", "ydrunmax", "ydrunsum", "ydrunmean", "ydrunavg", "ydrunstd", "ydrunstd1", "ydrunvar", "ydrunvar1"}
#define  YhourarithOperators    {"yhouradd", "yhoursub", "yhourmul", "yhourdiv"}
#define  YhourstatOperators     {"yhourrange", "yhourmin", "yhourmax", "yhoursum", "yhourmean", "yhouravg", "yhourstd", "yhourstd1", "yhourvar", "yhourvar1"}
#define  YmonarithOperators     {"ymonadd", "ymonsub", "ymonmul", "ymondiv"}
#define  YseasarithOperators    {"yseasadd", "yseassub", "yseasmul", "yseasdiv"}
#define  YmonpctlOperators      {"ymonpctl"}
#define  YmonstatOperators      {"ymonrange", "ymonmin", "ymonmax", "ymonsum", "ymonmean", "ymonavg", "ymonstd", "ymonstd1", "ymonvar", "ymonvar1"}
#define  YseaspctlOperators     {"yseaspctl"}
#define  YseasstatOperators     {"yseasrange", "yseasmin", "yseasmax", "yseassum", "yseasmean", "yseasavg", "yseasstd", "yseasstd1", "yseasvar", "yseasvar1"}
#define  ZonstatOperators       {"zonrange", "zonmin", "zonmax", "zonsum", "zonmean", "zonavg", "zonstd", "zonstd1", "zonvar", "zonvar1", "zonpctl"}

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

// HIRLAM_EXTENSIONS
#define  SelmultiOperators      {"selmulti", "delmulti", "changemulti"}
#define  WindTransOperators     {"uvDestag", "rotuvN", "rotuvNorth", "projuvLatLon"}
#define  SamplegridOperators    {"samplegrid", "subgrid"}

/* clang-format on */
#endif

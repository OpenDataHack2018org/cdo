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
#ifndef MODULE_DEFINITIONS_H
#define MODULE_DEFINITIONS_H
/* \cond */
#include "operator_definitions.h"
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
void *EstFreq(void *argument);
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
void *NCL_wind(void *argument);
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
void *Samplegridicon(void *argument);
void *Seascount(void *argument);
void *Seaspctl(void *argument);
void *Seasstat(void *argument);
void *Selbox(void *argument);
void *Selgridcell(void *argument);
void *Select(void *argument);
void *Selvar(void *argument);
void *Seloperator(void *argument);
void *Selrec(void *argument);
void *Seltime(void *argument);
void *Selyearidx(void *argument);
void *Set(void *argument);
void *Setattribute(void *argument);
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
void *Showattribute(void *argument);
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

// void *Hi(void *argument);
void *Wct(void *argument);

void *Magplot(void *argument);
void *Magvector(void *argument);
void *Maggraph(void *argument);

// HIRLAM_EXTENSIONS
void *Selmulti(void *argument);    // "selmulti", "delmulti"
void *WindTrans(void *argument);   // "uvDestag", "rotuvN", "rotuvNorth", "projuvLatLon"
void *Samplegrid(void *argument);  // "samplegrid", "subgrid"

/* \endcond */
#endif

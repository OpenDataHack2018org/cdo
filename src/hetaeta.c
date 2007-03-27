
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

const double ap0 = 100000.0;
const double apr = 101325.0;
const double aipr = 1.0/101325.0;

const double epsilon = 0.622;
const double rair = 287.04;
const double cpair = 1004.6;

const double eta_pbl = 0.8;

const double g = 9.81;

FILE *old, *new;

void interpolate_linear(int n, double *xa, double *ya, double x, double *y)
{
  int klo, khi, k;
  double h, d;

  klo = 0;
  khi = n-1;

  while (khi-klo > 1) {
    k = (khi+klo)/2;
      if (xa[k] > x) {
	khi = k;
      } else {
	klo = k;
      }
  }
  h = xa[khi]-xa[klo];
  d = ya[khi]-ya[klo];
  *y = ya[klo]+(x-xa[klo])/h*d;

  return;
}

double esat(double temperature)
{
  double zes;
  double es;
  double tc;

  tc = temperature-273.16;
  if (tc < 0.0) {
    zes = 21.8745584*tc / (temperature-7.66);
  } else {
    zes =  17.2693882*tc / (temperature-35.86);
  }
  es = 610.78*exp(zes);

  return es;
}


void hetaeta(int in_nlev, double *in_ah, double *in_bh,
             double in_fis, double in_ps, 
             double *in_t, double *in_q, double *in_u,   double *in_v,  double *in_cl,double *in_ci, double *in_cc,
             int nlev, double *ah, double *bh, 
             double fis, double *ps, 
             double *t, double *q, double *u, double *v, 
	     double *cl, double *ci, double *cc,
	     double *tscor, double *pscor, double *secor)
{
  double epsm1i, zdff, zdffl, ztv, zb, zbb, zc, zps;
  double zsump, zsumpp, zsumt, zsumtp;
  double dfi, fiadj, dteta;
  double pbl_lim, pbl_lim_need;
  int jblt, jjblt;
  int k;
  int jlev, jlevr, jnop;

  double *in_etah, *in_ph, *in_lnph, *in_fi;
  double *in_af, *in_bf, *in_etaf, *in_pf, *in_lnpf;
  double *in_tv, *in_theta, *in_rh;
  double *etah, *ph, *lnph, *fi;
  double *af, *bf, *etaf, *pf, *lnpf;

  double *rhpbl, *upbl, *vpbl;
  double *thetapbl, *clpbl, *cipbl, *ccpbl;

  double *rh;

  double *w1, *w2;

  int *jl1, *jl2;
  
  int in_nlevp1;
  int nlevp1;

  old = fopen("old.dat","w");
  new = fopen("new.dat","w");

  in_nlevp1 = in_nlev+1;
  nlevp1    = nlev+1;

  in_etah  = (double *) malloc(in_nlevp1*sizeof(double));			   
  in_ph    = (double *) malloc(in_nlevp1*sizeof(double));			   
  in_lnph  = (double *) malloc(in_nlevp1*sizeof(double));			   
  in_fi    = (double *) malloc(in_nlevp1*sizeof(double));			   
										   
  in_af    = (double *) malloc(in_nlev*sizeof(double));				   
  in_bf    = (double *) malloc(in_nlev*sizeof(double));				   
										   
  in_etaf  = (double *) malloc(in_nlev*sizeof(double));				   
  in_pf    = (double *) malloc(in_nlev*sizeof(double));				   
  in_lnpf  = (double *) malloc(in_nlev*sizeof(double));				   
										   
  in_tv    = (double *) malloc(in_nlev*sizeof(double));				   
  in_theta = (double *) malloc(in_nlev*sizeof(double));				   
  in_rh    = (double *) malloc(in_nlev*sizeof(double));				   
  										   
  etah     = (double *) malloc(nlevp1*sizeof(double));				   
  ph       = (double *) malloc(nlevp1*sizeof(double));				   
  lnph     = (double *) malloc(nlevp1*sizeof(double));				   
  fi       = (double *) malloc(nlevp1*sizeof(double));				   
										   
  af       = (double *) malloc(nlev*sizeof(double));				   
  bf       = (double *) malloc(nlev*sizeof(double));				   
										   
  etaf     = (double *) malloc(nlev*sizeof(double));				   
  pf       = (double *) malloc(nlev*sizeof(double));				   
  lnpf     = (double *) malloc(nlev*sizeof(double));				   
										   
										   
  rhpbl    = (double *) malloc(nlev*sizeof(double));				   
  upbl     = (double *) malloc(nlev*sizeof(double));				   
  vpbl     = (double *) malloc(nlev*sizeof(double));				   
										   
  thetapbl = (double *) malloc(nlev*sizeof(double));				   
  clpbl    = (double *) malloc(nlev*sizeof(double));				   
  cipbl    = (double *) malloc(nlev*sizeof(double));				   
  ccpbl    = (double *) malloc(nlev*sizeof(double));				   
										   
  rh       = (double *) malloc(nlev*sizeof(double));				   
										   
  w1       = (double *) malloc(nlev*sizeof(double));				   
  w2       = (double *) malloc(nlev*sizeof(double));				   
										   
  jl1      = (int *)    malloc(nlev*sizeof(int));				   
  jl2      = (int *)    malloc(nlev*sizeof(int));                                  
  
  for (k = 0; k < in_nlev; k++) {
    in_af[k] = 0.5*(in_ah[k]+in_ah[k+1]);
    in_bf[k] = 0.5*(in_bh[k]+in_bh[k+1]);
  }

  in_etah[in_nlev] = in_ah[in_nlev]*aipr+in_bh[in_nlev];
  for (k = 0; k < in_nlev; k++) { 
    in_etah[k] = in_ah[k]*aipr+in_bh[k];
    in_etaf[k] = in_af[k]*aipr+in_bf[k];
  }

  for (k = 0; k < nlev; k++) {
    af[k] = 0.5*(ah[k]+ah[k+1]);
    bf[k] = 0.5*(bh[k]+bh[k+1]);
  }

  etah[nlev] = ah[nlev]*aipr+bh[nlev];
  jblt = nlev;
  for (k = nlev-1; k >= 0; k--) { 
    etah[k] = ah[k]*aipr+bh[k];
    etaf[k] = af[k]*aipr+bf[k];
    if (etah[k] > eta_pbl) jblt = k;
  }

  for (k = 0; k < nlev; k++) { 
    if (etaf[k] <= in_etaf[0]) {
      jl1[k] = 0;
      jl2[k] = 1;
      w2[k]  = 0.0;
    } else if ( etaf[k] >= in_etaf[in_nlev-1]) {
      jl1[k] = in_nlev-2;
      jl2[k] = in_nlev-1;
      w2[k]  = 1.0;
    } else {
      for (jlev = in_nlev-2; jlev >= 1; jlev--) {
	jl1[k] = jlev;
	if (etaf[k] > in_etaf[jlev]) break;
      }
      jl2[k] = jl1[k]+1;
      w2[k] = log(etaf[k]/in_etaf[jl1[k]])
	/log(in_etaf[jl2[k]]/in_etaf[jl1[k]]);
    }
    w1[k] = 1.0-w2[k];
  }

  in_ph[0]       =  0.0;
  in_lnph[0]     = -1.0; 
  for (k = 1; k < in_nlevp1; k++) {
    in_ph[k]   = in_ah[k]+in_bh[k]*in_ps;
    in_lnph[k] = log(in_ph[k]);
  } 

  for (k = 0; k < in_nlev; k++) {
    in_pf[k]     = in_af[k]+in_bf[k]*in_ps;
    in_lnpf[k]   = log(in_pf[k]);
  }

  epsm1i = 1.0/epsilon-1.0;
  for (k = 0; k < in_nlev; k++) {
    in_tv[k]    = (1.0+epsm1i*in_q[k])*in_t[k];
    in_rh[k]    = in_q[k]*in_pf[k]/(epsilon*esat(in_t[k]));
    in_theta[k] = in_t[k]*pow(apr/in_pf[k],rair/cpair);
  };

  in_fi[0] = 0.0;
  in_fi[in_nlev] = in_fis;
  for (k = in_nlev-1; k > 0; k--) {
    in_fi[k] = in_fi[k+1]+rair*in_tv[k]*(in_lnph[k+1]-in_lnph[k]);
  }

  for (k = in_nlev-1; k >= 0; k--) { 
    fprintf(old, "%3d %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f\n", 
	   k, in_fi[k]/g, in_pf[k], in_t[k], in_q[k], in_u[k], in_v[k],
              in_cl[k], in_ci[k], in_cc[k]);
  }


  for (k = in_nlev-2; k > 0; k--) {
    jlev = k;
    if (fis < in_fi[k]) break;
  }
  zdff = in_fi[jlev+1]-fis;

  for (k = jlev-1; k > 0; k--) {
    jlevr = k;
    zdffl = in_fi[k]-in_fi[jlev+1];
    if (zdffl >= zdff) break;
  }

  jnop = jlev+1-jlevr+1;

  zsumt  = 0.0;
  zsump  = 0.0;
  zsumpp = 0.0;
  zsumtp = 0.0;

  for (k = jlevr; k <= jlev+1; k++) {
    zsumt  = zsumt  + in_tv[k];
    zsump  = zsump  + in_lnpf[k];
    zsumpp = zsumpp + in_lnpf[k]*in_lnpf[k];
    zsumtp = zsumtp + in_tv[k]*in_lnpf[k];
  }

  zb = jnop*zsumpp - zsump*zsump;
  zc = (zsumt*zsumpp-zsump*zsumtp)/zb;
  zb = (jnop*zsumtp-zsump*zsumt)/zb;

  zps = in_lnph[jlev];

  if (fabs(zb) < 1.0e-20) {
    *ps = exp(zps+(in_fi[jlev]-fis)/(zc*rair));
  } else {
    zbb = zc*zc + zb*(zps*(zb*zps+2.0*zc)+2.0*(in_fi[jlev]-fis)/rair);
    *ps = exp((sqrt(zbb)-zc)/zb);
  }

  ph[0]       =  0.0;
  lnph[0]     = -1.0; 
  for (k = 1; k < nlevp1; k++) {
    ph[k]   = ah[k]+bh[k]* *ps;
    lnph[k] = log(ph[k]);
  } 

  for (k = 0; k < nlev; k++) {
    pf[k]     = af[k]+bf[k]* *ps;
    lnpf[k]   = log(pf[k]);
  }

  for (k = 1; k < in_nlevp1; k++) {
    jlev = k;
    if (in_ph[k] > 40000.0) break;
  }

  fiadj = in_fi[jlev]+(in_fi[jlev-1]-in_fi[jlev])*
    log(in_ph[jlev]/40000.0)/log(in_ph[jlev]/in_ph[jlev-1]);

  pbl_lim = in_ps*eta_pbl;
  jjblt = nlev-1;
  for (k = nlev-1; k > 0; k--) {
    pbl_lim_need = *ps *etah[k];
    if (pbl_lim > pbl_lim_need) break;
    jjblt = jjblt-1;
  }

  if (jblt < jjblt) jjblt = jblt;

  for (k = jjblt; k < nlev; k++) {

    thetapbl[k] = w1[k]*in_theta[jl1[k]]+w2[k]*in_theta[jl2[k]];

    rhpbl[k] = w1[k]*in_rh[jl1[k]]+w2[k]*in_rh[jl2[k]];

    upbl[k] = w1[k]*in_u[jl1[k]]+w2[k]*in_u[jl2[k]];
    vpbl[k] = w1[k]*in_v[jl1[k]]+w2[k]*in_v[jl2[k]];

    clpbl[k] = w1[k]*in_cl[jl1[k]]+ w2[k]*in_cl[jl2[k]];
    cipbl[k] = w1[k]*in_ci[jl1[k]]+ w2[k]*in_ci[jl2[k]];
    ccpbl[k] = w1[k]*in_cc[jl1[k]]+ w2[k]*in_cc[jl2[k]];
  }

  for (k = 0; k <= jjblt; k++) {
    interpolate_linear (in_nlev, in_pf, in_t,  pf[k], &t[k]); 
    interpolate_linear (in_nlev, in_pf, in_rh, pf[k], &rh[k]); 
    interpolate_linear (in_nlev, in_pf, in_u,  pf[k], &u[k]); 
    interpolate_linear (in_nlev, in_pf, in_v,  pf[k], &v[k]); 
    interpolate_linear (in_nlev, in_pf, in_cl, pf[k], &cl[k]); 
    interpolate_linear (in_nlev, in_pf, in_ci, pf[k], &ci[k]); 
    interpolate_linear (in_nlev, in_pf, in_cc, pf[k], &cc[k]); 
  }

  dteta = t[jjblt]*pow(apr/pf[jjblt],rair/cpair)-thetapbl[jjblt];

  rh[jjblt] = 0.5*(rh[jjblt]+rhpbl[jjblt]);
  u [jjblt] = 0.5*(u [jjblt]+upbl [jjblt]);
  v [jjblt] = 0.5*(v [jjblt]+vpbl [jjblt]);
  cl[jjblt] = 0.5*(cl[jjblt]+clpbl[jjblt]);
  ci[jjblt] = 0.5*(ci[jjblt]+cipbl[jjblt]);
  cc[jjblt] = 0.5*(cc[jjblt]+ccpbl[jjblt]);

  for (k = jjblt+1; k < nlev; k++) {
    t [k] = (thetapbl[k]+dteta)*pow(pf[k]/apr,rair/cpair);
    rh[k] = rhpbl[k];
    u [k] = upbl [k];
    v [k] = vpbl [k];
    cl[k] = clpbl[k];
    ci[k] = cipbl[k];
    cc[k] = ccpbl[k];
  }

  for (k = 0; k < nlev; k++) {
    q[k] = rh[k]*epsilon*esat(t[k])/pf[k];
  }

  fi[nlev] = fis;
  fi[0]    = -1.0;
    
  for (k = nlev-1; k > 0; k--) {
    fi[k] = fi[k+1]+rair*t[k]*(lnph[k+1]-lnph[k])*(1.0+epsm1i*q[k]);
  }

  for (k = nlev-1; k > 0; k--) { 
    jlev = k; 
    if (ph[k] < 40000.0) break;
  }

  dfi = fiadj-(fi[jlev+1]+(fi[jlev]-fi[jlev+1])* 
	       log(ph[jlev+1]/40000.0)/log(ph[jlev+1]/ph[jlev]));
  ztv     = (1.0+epsm1i*q[nlev-1])*t[nlev-1];
  *ps = *ps *exp(dfi/(rair*ztv));

  for (k = 0; k < nlev; k++) {
    pf[k] = af[k] + bf[k]* *ps;
  }

  for (k = 0; k < nlev; k++) {
    q[k] = rh[k]*epsilon*esat(t[k])/pf[k];
  }

  for (k = nlev-1; k >= 0; k--) { 
    fprintf(new, "%3d %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f\n", 
	   k, fi[k]/g, pf[k], t[k], q[k], u[k], v[k], cl[k], ci[k], cc[k]);
  }

  *tscor = dteta*pow(*ps/apr,rair/cpair);
  *pscor = pow(*ps/in_ps,rair/cpair);
  *secor = in_tv[in_nlev-1]*
    (cpair+rair*(1.0-in_ph[in_nlev-1]
		 /(in_ps-in_ph[in_nlev-1])
		 *log(in_ps/in_ph[in_nlev-1])));

  fclose(old);
  fclose(new);

  free(jl2);      
  free(jl1);      

  free(w2);       
  free(w1);       

  free(rh);      

  free(ccpbl);    
  free(cipbl);    
  free(clpbl);     
  free(thetapbl); 

  free(vpbl);     
  free(upbl);     
  free(rhpbl);    

  free(lnpf);     
  free(pf);       
  free(etaf);     

  free(bf);       
  free(af);       

  free(fi);       
  free(lnph);     
  free(ph);       
  free(etah);    

  free(in_rh);    
  free(in_theta); 
  free(in_tv);    

  free(in_lnpf);  
  free(in_pf);    
  free(in_etaf);  

  free(in_bf);    
  free(in_af);    

  free(in_fi);    
  free(in_lnph);  
  free(in_ph);    
  free(in_etah);      

  return;
}

/*
int main (int argc, char *argv[])
{
  double a1[41] = {
       0.00000000000000000,       2000.00000000000000000,       4000.00000000000000000,
    6000.00000000000000000,       8000.00000000000000000,       9976.13671875000000000,
   11902.14453125000000000,      13722.03125000000000000,      15379.80468750000000000,
   16819.47265625000000000,      18045.18359375000000000,      19027.69531250000000000,
   19755.10937500000000000,      20222.20312500000000000,      20429.86328125000000000,
   20384.48046875000000000,      20097.40234375000000000,      19584.32812500000000000,
   18864.75000000000000000,      17961.35937500000000000,      16899.46875000000000000,
   15706.44921875000000000,      14411.12500000000000000,      13043.21875000000000000,
   11632.75781250000000000,      10209.50000000000000000,       8802.35546875000000000,
    7438.80468750000000000,       6144.31640625000000000,       4941.77734375000000000,
    3850.91333007812500000,       2887.69653320312500000,       2063.77978515625000000,
    1385.91259765625000000,        855.36181640625000000,        467.33349609375000000,
     210.39390563964843750,         65.88919067382812500,          7.36769962310791016,
       0.00000000000000000,          0.00000000000000000 };

  double b1[41] = {
    0.00000000000000000,    0.00000000000000000,    0.00000000000000000,
    0.00000000000000000,    0.00000000000000000,    0.00039085815660655,
    0.00182679994031787,    0.00513499975204468,    0.01114289835095406,
    0.02067789807915688,    0.03412120044231415,    0.05169039964675903,
    0.07353377342224121,    0.09967470169067383,    0.13002246618270874,
    0.16438430547714233,    0.20247590541839600,    0.24393308162689209,
    0.28832298517227173,    0.33515489101409912,    0.38389205932617188,
    0.43396288156509399,    0.48477149009704590,    0.53570991754531860,
    0.58616840839385986,    0.63554751873016357,    0.68326860666275024,
    0.72878581285476685,    0.77159661054611206,    0.81125342845916748,
    0.84737491607666016,    0.87965691089630127,    0.90788388252258301,
    0.93194031715393066,    0.95182150602340698,    0.96764522790908813,
    0.97966271638870239,    0.98827010393142700,    0.99401938915252686,
    0.99763011932373047,    1.00000000000000000 };

  double a2[20] = {
       0.00000000000000000,       2000.00000000000000000,       4000.00000000000000000,   
    6046.10937500000000000,       8267.92968750000000000,      10609.51171875000000000,   
   12851.10156250000000000,      14698.50000000000000000,      15861.12890625000000000,   
   16116.23828125000000000,      15356.92187500000000000,      13621.46093750000000000,  
   11101.55859375000000000,       8127.14453125000000000,       5125.14062500000000000,   
    2549.96899414062500000,        783.19506835937500000,          0.00000000000000000,   
       0.00000000000000000,          0.00000000000000000 };  

  double b2[20] = {
    0.00000000000000000,    0.00000000000000000,    0.00000000000000000,
    0.00033899326808751,    0.00335718691349030,    0.01307003945112228,
    0.03407714888453484,    0.07064980268478394,    0.12591671943664551,
    0.20119541883468628,    0.29551959037780762,    0.40540921688079834,
    0.52493220567703247,    0.64610791206359863,    0.75969839096069336,
    0.85643762350082397,    0.92874687910079956,    0.97298520803451538,
    0.99228149652481079,    1.00000000000000000 };

  double iu[19] = {
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0 };

  double iv[19] = {
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0 };

  double it[19] = {
    224.2257,  212.7505,  209.553,  208.7785,  211.7619,  220.2336,  221.2698,
    220.3876,  227.1461,  237.6735, 248.1776,  258.1013,  264.4792,  269.1322,
    271.9017,  275.6761,  279.819,  282.2512,  284.141 };

  double iq[19] = {
    2.512447e-06,  2.176736e-06,  2.170464e-06,  2.01653e-06,  1.805185e-06,
    1.726813e-06,  3.75322e-06,   8.901303e-06,  3.285719e-05, 0.0001270178,
    0.0003347051,  0.0007223329,  0.001228461,   0.001733165,  0.002967748,
    0.004558741,   0.004706143,   0.004668835,   0.004677606 };

  double icl[19] = {
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  -4.987459e-40,  -4.791847e-39,  -3.970467e-23,
    -1.902515e-23,  -1.694066e-21,  -3.705769e-22,  -1.799945e-21,  -4.632211e-22,
    2.072752e-05,  0.000149563,  -1.482308e-20,  -2.541099e-21,  5.033612e-05 };

  double ici[19] = {
    -4.408104e-37, 0.0,  0.0,  -2.003328e-25,  -9.305782e-24,  -2.15067e-23,  -9.926167e-23,
    -1.958764e-21,  -8.735027e-22,  -2.779327e-22,  -2.117582e-21,  -1.323489e-21,  -8.470329e-22,
    -4.102816e-22,  -1.429368e-21,  -2.646978e-21,  -5.029258e-22,  -8.205632e-22,  -1.588187e-21 };

  double icc[19] = {
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.4445496,  0.0,  0.0,  0.001098633} ;

  double ifis = 9.121094, ips = 102511.8;

  double fis = 9.121094, ps;
 
  double *t, *q, *u, *v, *cl, *ci, *cc;

  double tscor, pscor, secor; 

  t  = (double *) malloc(40*sizeof(double));
  q  = (double *) malloc(40*sizeof(double));
  u  = (double *) malloc(40*sizeof(double));
  v  = (double *) malloc(40*sizeof(double));
  cl = (double *) malloc(40*sizeof(double));
  ci = (double *) malloc(40*sizeof(double));
  cc = (double *) malloc(40*sizeof(double));

  hetaeta(19, a2, b2,
          ifis, ips,
          it, iq, iu, iv, icl, ici, icc,
	  40, a1, b1,
          fis, &ps,
          t, q, u, v, cl, ci, cc,
	  &tscor, &pscor, &secor);

  return 0;
}
*/

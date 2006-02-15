
typedef struct {
  int nlon;
  int nlat;
  int trunc;
  int poldim;
  int ifax[10];
  double *trig;
  double *poli;
  double *pold;
}
SPTRANS;


SPTRANS *sptrans_new(int nlon, int nlat, int trunc);
void sptrans_delete(SPTRANS *sptrans);

void grid2spec(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut);
void spec2grid(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut);

void spec2spec(int gridIDin, double *arrayIn, int gridIDout, double *arrayOut);
void speccut(int gridIDin, double *arrayIn, double *arrayOut, int waves[]);

void sp2fc(double *sa, double *fa, double *poli, int nlev, int nlat, int nfc, int nt);
void fc2sp(double *fa, double *sa, double *poli, int nlev, int nlat, int nfc, int nt);

void fc2gp(double *trig, int *ifax, double *fc, double *gp, int nlat, int nlon, int nlev, int nfc);
void gp2fc(double *trig, int *ifax, double *gp, double *fc, int nlat, int nlon, int nlev, int nfc);

void sp2sp(double *arrayIn, int truncIn, double *arrayOut, int truncOut);
void spcut(double *arrayIn, double *arrayOut, int trunc, int waves[]);

void phcs(double *pnm, double *hnm, int waves, double pmu,
	  double *ztemp1, double *ztemp2);
void jspleg1(double *pleg, double plat, int ktrunc, double *work);

void fft_set(double *trigs, int *ifax, int n);

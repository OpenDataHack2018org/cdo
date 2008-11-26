#ifndef _GRID_H
#define _GRID_H

int gridToZonal(int gridID);
int gridToMeridional(int gridID);
int gridToCell(int gridID);
int gridToCurvilinear(int gridID);

/* GME grid */
struct cart {
  double x[3];
};

struct geo {
  double lon;
  double lat;
};

double areas(struct cart *dv1, struct cart *dv2, struct cart *dv3);
struct cart gc2cc(struct geo *position);
void factorni(int kni, int *kni2, int *kni3);
void gme_grid_restore(double *p, int ni, int nd);
void gme_grid(int gridsize, double *rlon, double *rlat,
	      double *blon, double *blat, int *imask,
              int ni, int nd, int ni2, int ni3);

/* Rotated grid */
double rls_to_rl(double phis, double rlas, double polphi, double pollam);
double phs_to_ph(double phis, double rlas, double polphi);
double rl_to_rls(double phi, double rla, double polphi, double pollam);
double ph_to_phs(double phi, double rla, double polphi, double pollam);
void usvs_to_uv(double us, double vs, double phi, double rla,
		double polphi, double pollam, double *u, double *v);

/* Lambert Conformal grid */
int W3FB12(double xi, double xj, double alat1, double elon1, double dx,
	   double elonv, double alatan, double *alat, double *elon);

#endif  /* _GRID_H */

#ifndef _VINTERP_H
#define _VINTERP_H

#if defined(__cplusplus)
extern "C" {
#endif

void h2p(double *phlev, double *hlev, int nphlev);

void presh(double * restrict fullp, double * halfp, const double *vct, const double *ps, long nhlev, long ngp);

void genind(int *nx, const double *plev, const double *fullp, long ngp, long nplev, long nhlev);
void genindmiss(int *nx, double *plev, int ngp, int nplev, double *ps_prog, int *pnmiss);

void extra_P(double *slp, double *halfp, double *fullp, double *geop, double *temp, long ngp);

void interp_T(const double *geop, const double *gt, double *pt, const double *fullp, const double *halfp,
              const int *nx, const double *plev, long nplev, long ngp, long nhlev, double missval);
void interp_Z(const double *geop, const double *gz, double *pz, const double *fullp, const double *halfp,
	      const int *nx, const double *gt, const double *plev, long nplev, long ngp, long nhlev, double missval);
void interp_X(const double *gt, double *pt, const double *hyb_press, const int *nx, const double *plev, long nplev,
	      long ngp, long nhlev, double missval);

#if defined(__cplusplus)
}
#endif

#endif  /* _VINTERP_H */

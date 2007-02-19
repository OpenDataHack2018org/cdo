#ifndef _STATISTIC_H

#define _STATISTIC_H

void eigen_solution_of_symmetric_matrix (double **a, double *eig_val,
					 int n, char *prompt);
int solution_of_linear_equation (double **a, double *b, int n);
int inverse_of_matrix (double **a, double **b, int n);
void fft (double *real, double *imag, int n, int sign);
void ft (double *real, double *imag, int n, int sign, char *prompt);
double lngamma (double x);
double beta (double a, double b, char *prompt);
double incomplete_gamma (double a, double x, char *prompt);
double incomplete_beta (double a, double b, double x, char *prompt);
double normal_density (double x);
double normal (double x, char *prompt);
double normal_inv (double p, char *prompt);
double student_t_density (double n, double x, char *prompt);
double student_t (double n, double x, char *prompt);
double student_t_inv (double n, double p, char *prompt);
double chi_square_density (double n, double x, char *prompt);
double chi_square (double n, double x, char *prompt);
double chi_square_inv (double n, double p, char *prompt);
void chi_square_constants (double n, double p, double *c1, double *c2, char *prompt);
double beta_distr_density (double a, double b, double x, char *prompt);
double beta_distr (double a, double b, double x, char *prompt);
double beta_distr_inv (double a, double b, double p, char *prompt);
void beta_distr_constants (double a, double b, double p, double *c1, double *c2, char *prompt);
double fisher (double m, double n, double x, char *prompt);

#endif

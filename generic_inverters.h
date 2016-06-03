// Wed Jan  6 09:49:39 EST 2016
// Evan S Weinberg
// Include file for generic inverters.

#ifndef ESW_INVERTER
#define ESW_INVERTER

// Uncomment to print warnings. 
//#define VERBOSE_WARN

#include <complex>

using std::complex; 

// Struct that contains information that
// all matrix functions return.
// Return for various matrix functions.
struct inversion_info
{
  double resSq; // squared residual.
  int iter; // number of iterations.
  bool success; // did we reach residual?
  std::string name; // name of algorithm.
};

// Performs gaussian elimination.
int gaussian_elimination(double* x, double* b, double** matrix, int size);
int gaussian_elimination(complex<double>* x, complex<double>* b, complex<double>** matrix, int size);

// Solves lhs = A^(-1) rhs with conjugate gradient
// Requires matrix to be symmetric (Hermitian) positive definite 
inversion_info minv_vector_cg(double  *phi, double  *phi0, int size, int max_iter, double res, void (*matrix_vector)(double*,double*,void*), void* extra_info);

inversion_info minv_vector_cg(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double res, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info);

// Solves lhs = A^(-1) rhs using generalized conjugate residual
// Makes no assumptions about the matrix. 
inversion_info minv_vector_gcr(double  *phi, double  *phi0, int size, int max_iter, double res, void (*matrix_vector)(double*,double*,void*), void* extra_info);

inversion_info minv_vector_gcr(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double res, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info);

// Solves lhs = A^(-1) with GCR(m), where m is restart_freq. 
inversion_info minv_vector_gcr_restart(double  *phi, double  *phi0, int size, int max_iter, double res, int restart_freq, void (*matrix_vector)(double*,double*,void*), void* extra_info);

inversion_info minv_vector_gcr_restart(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double res, int restart_freq, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info);

// Solves lhs = A^(-1) rhs with bicgstab
// Makes no assumptions about the matrix.
inversion_info minv_vector_bicgstab(double  *phi, double  *phi0, int size, int max_iter, double res, void (*matrix_vector)(double*,double*,void*), void* extra_info);

inversion_info minv_vector_bicgstab(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double res, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info);

// Solves lhs = A^(-1) rhs with GMRES, no restarts.
// Makes no assumption about matrix.
inversion_info minv_vector_gmres_norestart(double  *phi, double  *phi0, int size, int max_iter, double res, void (*matrix_vector)(double*,double*,void*), void* extra_info);

inversion_info minv_vector_gmres_norestart(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double res, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info);

// Solves lhs = A^(-1) rhs with GMRES with restarts.
// Makes no assumption about matrix.
inversion_info minv_vector_gmres_restart(double  *phi, double  *phi0, int size, int max_iter, double res, int restart_freq, void (*matrix_vector)(double*,double*,void*), void* extra_info);

inversion_info minv_vector_gmres_restart(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double res, int restart_freq, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info);

// Solves lhs = A^(-1) rhs with SOR.
// Makes no assumption about matrix, but omega has to be properly set for convergence. 
inversion_info minv_vector_sor(double  *phi, double  *phi0, int size, int max_iter, double eps, double omega, void (*matrix_vector)(double*,double*,void*), void* extra_info);

inversion_info minv_vector_sor(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, double omega, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info);

// Solves lhs = A^(-1) rhs with MinRes.
// Assumes the symmetric part of the matrix is positive definite.
inversion_info minv_vector_minres(double  *phi, double  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(double*,double*,void*), void* extra_info);

inversion_info minv_vector_minres(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info);


// Could also want (as a smoother): 
// Gauss-Seidel: https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method
// Jacobi iterations: https://en.wikipedia.org/wiki/Jacobi_method


#endif // define ESW_INVERTER

// Wed Aug 24 15:46:50 EDT 2016
// Evan S Weinberg
// Include file for Preconditioned BiCGStab

#ifndef ESW_INVERTER_BICGSTAB_PRECOND
#define ESW_INVERTER_BICGSTAB_PRECOND

#include <string>
#include <complex>
using std::complex;

#include "inverter_struct.h"
#include "verbosity.h"

// Solves lhs = A^(-1) rhs with bicgstab
// Makes no assumptions about the matrix.
inversion_info minv_vector_bicgstab_precond(double  *phi, double  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(double*,double*,void*), void* extra_info, void (*precond_matrix_vector)(double*,double*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verb);
//inversion_info minv_vector_bicgstab_restart(double  *phi, double  *phi0, int size, int max_iter, double res, int restart_freq, void (*matrix_vector)(double*,double*,void*), void* extra_info, inversion_verbose_struct* verbosity = 0);

//inversion_info minv_vector_bicgstab(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double res, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verbosity = 0);
//inversion_info minv_vector_bicgstab_restart(complex<double>  *phi, complex<double> *phi0, int size, int max_iter, double res, int restart_freq, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verbosity = 0);

#endif
// Wed Aug 24 15:46:50 EDT 2016
// Evan S Weinberg
// Include file for Flexibly Preconditioned CG

#ifndef ESW_INVERTER_FPCG
#define ESW_INVERTER_FPCG

#include <string>
#include <complex>
using std::complex;

#include "inverter_struct.h"
#include "verbosity.h"

// Flexibly Preconditioned Conjugate Gradient:
inversion_info minv_vector_cg_flex_precond(double  *phi, double  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(double*,double*,void*), void* extra_info, void (*precond_matrix_vector)(double*,double*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verbosity = 0);

inversion_info minv_vector_cg_flex_precond(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, void (*precond_matrix_vector)(complex<double>*,complex<double>*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verbosity = 0);

// Restarted Flexibly Preconditioned Conjugate Gradient:
inversion_info minv_vector_cg_flex_precond_restart(double  *phi, double  *phi0, int size, int max_iter, double eps, int restart_freq, void (*matrix_vector)(double*,double*,void*), void* extra_info, void (*precond_matrix_vector)(double*,double*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verbosity = 0);

inversion_info minv_vector_cg_flex_precond_restart(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, int restart_freq, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, void (*precond_matrix_vector)(complex<double>*,complex<double>*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verbosity = 0);
#endif
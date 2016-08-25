// Wed Aug 24 15:46:50 EDT 2016
// Evan S Weinberg
// Include file for MinRes

#ifndef ESW_INVERTER_MINRES
#define ESW_INVERTER_MINRES

#include <string>
#include <complex>
using std::complex;

#include "inverter_struct.h"
#include "verbosity.h"

// Solves lhs = A^(-1) rhs with MinRes
// Assumes the symmetric part of the matrix is positive definite.
inversion_info minv_vector_minres(double  *phi, double  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(double*,double*,void*), void* extra_info, inversion_verbose_struct* verbosity = 0);
inversion_info minv_vector_minres(double  *phi, double  *phi0, int size, int max_iter, double eps, double omega, void (*matrix_vector)(double*,double*,void*), void* extra_info, inversion_verbose_struct* verbosity = 0);

inversion_info minv_vector_minres(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verbosity = 0);
inversion_info minv_vector_minres(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, double omega, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verbosity = 0);



#endif
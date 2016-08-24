// Wed Aug 24 15:46:50 EDT 2016
// Evan S Weinberg
// Include file for SOR

#ifndef ESW_INVERTER_SOR
#define ESW_INVERTER_SOR

#include <string>
#include <complex>
using std::complex;

#include "inverter_struct.h"
#include "verbosity.h"

// Solves lhs = A^(-1) rhs with SOR.
// Makes no assumption about matrix, but omega has to be properly set for convergence. 
inversion_info minv_vector_sor(double  *phi, double  *phi0, int size, int max_iter, double eps, double omega, void (*matrix_vector)(double*,double*,void*), void* extra_info, inversion_verbose_struct* verbosity = 0);

inversion_info minv_vector_sor(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, double omega, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verbosity = 0);

#endif
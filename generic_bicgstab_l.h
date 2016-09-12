// Wed Aug 24 15:46:50 EDT 2016
// Evan S Weinberg
// Include file for BiCGStab-l

#ifndef ESW_INVERTER_BICGSTAB_L
#define ESW_INVERTER_BICGSTAB_L

#include <string>
#include <complex>
using std::complex;

#include "inverter_struct.h"
#include "verbosity.h"

// Solves lhs = A^(-1) rhs with bicgstab
// Makes no assumptions about the matrix.
inversion_info minv_vector_bicgstab_l(double  *phi, double  *phi0, int size, int max_iter, double res, int l, void (*matrix_vector)(double*,double*,void*), void* extra_info, inversion_verbose_struct* verbosity = 0);
inversion_info minv_vector_bicgstab_l_restart(double  *phi, double  *phi0, int size, int max_iter, double res, int restart_freq, int l, void (*matrix_vector)(double*,double*,void*), void* extra_info, inversion_verbose_struct* verbosity = 0);

inversion_info minv_vector_bicgstab_l(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double res, int l, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verbosity = 0);
inversion_info minv_vector_bicgstab_l_restart(complex<double>  *phi, complex<double> *phi0, int size, int max_iter, double res, int restart_freq, int l, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verbosity = 0);

#endif
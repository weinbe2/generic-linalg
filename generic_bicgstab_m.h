// Wed Aug 24 15:46:50 EDT 2016
// Evan S Weinberg
// Include file for BICGSTAB-M

#ifndef ESW_INVERTER_BICGSTAB_M
#define ESW_INVERTER_BICGSTAB_M

#include <string>
#include <complex>
using std::complex;

#include "inverter_struct.h"
#include "verbosity.h"

// Multishift BiCGstab defined in http://arxiv.org/pdf/hep-lat/9612014.pdf
inversion_info minv_vector_bicgstab_m(double **phi, double *phi0, int n_shift, int size, int resid_freq_check, int max_iter, double eps, double* shifts, void (*matrix_vector)(double*,double*,void*), void* extra_info, bool worst_first = false, inversion_verbose_struct* verbosity = 0);
inversion_info minv_vector_bicgstab_m(complex<double> **phi, complex<double>  *phi0, int n_shift, int size, int resid_freq_check, int max_iter, double eps, double* shifts, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, bool worst_first = false, inversion_verbose_struct* verbosity = 0);

// multishift BICGSTAB starts with 0 initial guess, so a restarted version wouldn't work. 

#endif
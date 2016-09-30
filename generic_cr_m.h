// Wed Aug 24 15:46:50 EDT 2016
// Evan S Weinberg
// Include file for CR-M

// WARNING!!!
// MULTISHIFT CR IS A PARTICULARLY UNSTABLE ALGORITHM.
// TECHNICALLY, ONE SHOULD MONITOR CONVERGENCE BY LOOKING FOR \zeta_n < 1.
// On a more mathematical level, CR appears to generate a Krylov space in
// Ab, instead of in b (which is what CG does). This follows from noting
// the solution is A-orthogonal to the Krylov space K(A, b), as opposed to
// just orthogonal.
// Under a shift, a Krylov space in Ab -> a shift in (A+\sigma)b.
// Thus, the shifted operator doesn't generate the same Krylov space
// as the original operator!

#ifndef ESW_INVERTER_CR_M
#define ESW_INVERTER_CR_M

#include <string>
#include <complex>
using std::complex;

#include "inverter_struct.h"
#include "verbosity.h"

// Multishift CR defined in http://arxiv.org/pdf/hep-lat/9612014.pdf
inversion_info minv_vector_cr_m(double **phi, double *phi0, int n_shift, int size, int resid_freq_check, int max_iter, double eps, double* shifts, void (*matrix_vector)(double*,double*,void*), void* extra_info, inversion_verbose_struct* verbosity = 0);
inversion_info minv_vector_cr_m(complex<double> **phi, complex<double>  *phi0, int n_shift, int size, int resid_freq_check, int max_iter, double eps, double* shifts, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verbosity = 0);

// multishift CR starts with 0 initial guess, so a restarted version wouldn't work. 

#endif
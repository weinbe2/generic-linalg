// Wed Aug 24 15:46:50 EDT 2016
// Evan S Weinberg
// Include file for CG

#ifndef ESW_INVERTER_CG
#define ESW_INVERTER_CG

#include <string>
#include <complex>
using std::complex;

#include "inverter_struct.h"
#include "verbosity.h"

// CG with verbosity setting.
/*template <typename T, typename Treal> // = typename RealType<T>::Type >
inversion_info minv_vector_cg(T  *phi, T *phi0, int size, int max_iter, Treal eps, void (*matrix_vector)(T*,T*,void*), void* extra_info, inversion_verbose_struct* verbosity = 0);*/
inversion_info minv_vector_cg(double  *phi, double  *phi0, int size, int max_iter, double res, void (*matrix_vector)(double*,double*,void*), void* extra_info, inversion_verbose_struct* verbosity = 0);
inversion_info minv_vector_cg(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double res, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verbosity = 0);

#endif
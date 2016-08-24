// Wed Jan  6 09:49:39 EST 2016
// Evan S Weinberg
// Include file for generic inverters.

#ifndef ESW_INVERTER
#define ESW_INVERTER

// Uncomment to print warnings. 
//#define VERBOSE_WARN

#include <string>
#include <complex>

using std::complex; 

// Struct that contains information that
// all matrix functions return.
// Return for various matrix functions.
#include "inverter_struct.h"


#include "generic_traits.h"
#include "verbosity.h"

// Performs gaussian elimination.
#include "generic_gelim.h"

// Conjugate Gradient!

// Solves lhs = A^(-1) rhs with conjugate gradient
// Requires matrix to be symmetric (Hermitian) positive definite 

// CG with verbosity setting.
#include "generic_cg.h"

// Solves lhs = A^(-1) rhs with conjugate residual
// Requires matrix to be symmetric (Hermitian) 
#include "generic_cr.h"

// Solves lhs = A^(-1) rhs using generalized conjugate residual
// Makes no assumptions about the matrix. 
#include "generic_gcr.h"

// Solves lhs = A^(-1) rhs with bicgstab
// Makes no assumptions about the matrix.
#include "generic_bicgstab.h"

// Solves lhs = A^(-1) rhs with GMRES, no restarts.
// Makes no assumption about matrix.
#include "generic_gmres.h"

// Solves lhs = A^(-1) rhs with SOR.
// Makes no assumption about matrix, but omega has to be properly set for convergence. 
#include "generic_sor.h"

// Solves lhs = A^(-1) rhs with MinRes
// Assumes the symmetric part of the matrix is positive definite.
#include "generic_minres.h"

// MinRes with relaxation parameter:
inversion_info minv_vector_minres(double  *phi, double  *phi0, int size, int max_iter, double eps, double omega, void (*matrix_vector)(double*,double*,void*), void* extra_info, inversion_verbose_struct* verbosity = 0);

inversion_info minv_vector_minres(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, double omega, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verbosity = 0);


// Could also want (as a smoother): 
// Gauss-Seidel: https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method
// Jacobi iterations: https://en.wikipedia.org/wiki/Jacobi_method


#endif // define ESW_INVERTER

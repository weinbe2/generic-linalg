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

// Solves lhs = A^(-1) rhs with bicgstab-l
// Makes no assumptions about the matrix.
#include "generic_bicgstab_l.h"

// Solves lhs = A^(-1) rhs with GMRES, no restarts.
// Makes no assumption about matrix.
#include "generic_gmres.h"

// Solves lhs = A^(-1) rhs with SOR.
// Makes no assumption about matrix, but omega has to be properly set for convergence. 
#include "generic_sor.h"

// Solves lhs = A^(-1) rhs with MinRes
// Assumes the symmetric part of the matrix is positive definite.
#include "generic_minres.h"

// Could also want (as a smoother): 
// Gauss-Seidel: https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method
// Jacobi iterations: https://en.wikipedia.org/wiki/Jacobi_method

// Enumerates all current solvers without preconditioners.
// Add new inverters here.
enum minv_inverter
{
    MINV_CG = 0,
    MINV_CR = 1,
    MINV_GCR = 2,
    MINV_BICGSTAB = 3,
    MINV_BICGSTAB_L = 4,
    MINV_GMRES = 5,
    MINV_SOR = 6,
    MINV_MINRES = 7,
    MINV_INVALID = -1,
};

// A structure which holds all properties shared by different inverters, and
// custom properties that are only relevant for some inverters (i.e., are 
// ignored by others.
struct minv_inverter_params
{
    // General properties.
    double tol;
    int max_iters;
    bool restart;
    int restart_freq;
    
    // Solver-specific properties.
    
    // SOR:
    double sor_omega;
    
    // MinRes:
    double minres_omega; 
    
    // BiCGstab-L:
    int bicgstabl_l;
};

// A general function that takes an 'minv_inverter_params' and an 'minv_inverter' enum, then calls
// the right inverter. Not as efficient as it could be (due to overhead of conditional statements),
// but it's a bit more compact than other methods.
inversion_info minv_unpreconditioned(double* lhs, double* rhs, int size, minv_inverter type, minv_inverter_params& params, void (*matrix_vector)(double*,double*,void*), void* extra_info, inversion_verbose_struct* verbosity = 0);

inversion_info minv_unpreconditioned(complex<double>* lhs, complex<double>* rhs, int size, minv_inverter type, minv_inverter_params& params, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verbosity = 0);


#endif // define ESW_INVERTER

// Wed Jan  6 09:49:39 EST 2016
// Evan S Weinberg
// Include file for generic inverters.

#ifndef ESW_INVERTER_PRECOND
#define ESW_INVERTER_PRECOND

// Uncomment to print warnings. 
//#define VERBOSE_WARN

#include "verbosity.h"

// Preconditioned Conjugate Gradient:
#include "generic_cg_precond.h"

// Flexibly Preconditioned Conjugate Gradient:
#include "generic_cg_flex_precond.h"

// Variably Preconditioned GCR:
#include "generic_gcr_var_precond.h"

// (Flexibily) Preconditioned BiCGStab:
#include "generic_bicgstab_precond.h"

// Pre-implemented preconditioners: 
// Identity, MinRes, GCR.
#include "generic_precond.h"


// Enumerates all current solvers without preconditioners.
// Add new inverters here.
enum minv_inverter_precond
{
    MINV_PRE_CG = 0,
    MINV_PRE_FPCG = 1,
    MINV_PRE_VPGCR = 2,
    MINV_PRE_BICGSTAB = 3,
    MINV_PRE_INVALID = -1,
};

// A structure which holds all properties shared by different inverters, and
// custom properties that are only relevant for some inverters (i.e., are 
// ignored by others.
struct minv_inverter_precond_params
{
    // General properties.
    double tol;
    int max_iters;
    bool restart;
    int restart_freq;
    
    // Solver-specific properties.
    
    // None so far!
};

// A general function that takes an 'minv_inverter_precond_params' and an 'minv_inverter_precond' enum,
// then calls the right inverter. 
inversion_info minv_preconditioned(double* lhs, double* rhs, int size, minv_inverter_precond type, minv_inverter_precond_params& params, void (*matrix_vector)(double*,double*,void*), void* extra_info, void (*precond_matrix_vector)(double*,double*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verbosity = 0);

inversion_info minv_preconditioned(complex<double>* lhs, complex<double>* rhs, int size, minv_inverter_precond type, minv_inverter_precond_params& params, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, void (*precond_matrix_vector)(complex<double>*,complex<double>*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verbosity = 0);

#endif // define ESW_INVERTER_PRECOND

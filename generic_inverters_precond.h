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

#endif // define ESW_INVERTER_PRECOND

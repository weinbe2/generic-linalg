// Wed Aug 24 15:46:50 EDT 2016
// Evan S Weinberg
// Include file for various inverter preconditioners

#ifndef ESW_INVERTER_PREC
#define ESW_INVERTER_PREC

#include <string>
#include <complex>
using std::complex;

#include "inverter_struct.h"
#include "verbosity.h"

// Preconditioning function.
void identity_preconditioner(double* lhs, double* rhs, int size, void* extra_data, inversion_verbose_struct* verb = 0); 

void identity_preconditioner(complex<double>* lhs, complex<double>* rhs, int size, void* extra_data, inversion_verbose_struct* verb = 0);

/* BEGIN MinRes PRECONDITIONER */

// Preconditioning function.
void minres_preconditioner(double* lhs, double* rhs, int size, void* extra_data, inversion_verbose_struct* verb = 0);

void minres_preconditioner(complex<double>* lhs, complex<double>* rhs, int size, void* extra_data, inversion_verbose_struct* verb = 0); 

// Preconditioning struct. 
struct minres_precond_struct_real
{
    int n_step; 
    double rel_res; 
    void (*matrix_vector)(double*, double*, void*);
    void* matrix_extra_data; 
};

// Preconditioning struct. 
struct minres_precond_struct_complex
{
    int n_step; 
    double rel_res; 
    void (*matrix_vector)(complex<double>*, complex<double>*, void*);
    void* matrix_extra_data; 
};

/* END MinRes PRECONDITIONER */

/* BEGIN GCR PRECONDITIONER */

// Preconditioning function.
void gcr_preconditioner(double* lhs, double* rhs, int size, void* extra_data, inversion_verbose_struct* verb = 0); 

void gcr_preconditioner(complex<double>* lhs, complex<double>* rhs, int size, void* extra_data, inversion_verbose_struct* verb = 0); 

// Preconditioning struct. 
struct gcr_precond_struct_real
{
    int n_step; 
    double rel_res; 
    void (*matrix_vector)(double*, double*, void*);
    void* matrix_extra_data; 
};

// Preconditioning struct. 
struct gcr_precond_struct_complex
{
    int n_step; 
    double rel_res; 
    void (*matrix_vector)(complex<double>*, complex<double>*, void*);
    void* matrix_extra_data; 
};

/* END GCR PRECONDITIONER */


#endif
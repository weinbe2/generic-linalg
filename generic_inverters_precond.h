// Wed Jan  6 09:49:39 EST 2016
// Evan S Weinberg
// Include file for generic inverters.

#ifndef ESW_INVERTER_PRECOND
#define ESW_INVERTER_PRECOND

// Uncomment to print warnings. 
//#define VERBOSE_WARN

#include "generic_inverters.h"

// Preconditioned Conjugate Gradient:
inversion_info minv_vector_cg_precond(double  *phi, double  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(double*,double*,void*), void* extra_info, void (*precond_matrix_vector)(double*,double*,int,void*), void* precond_info);

inversion_info minv_vector_cg_precond(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, void (*precond_matrix_vector)(complex<double>*,complex<double>*,int,void*), void* precond_info);

// Variably Preconditioned GCR:
inversion_info minv_vector_gcr_var_precond(double  *phi, double  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(double*,double*,void*), void* extra_info, void (*precond_matrix_vector)(double*,double*,int,void*), void* precond_info);

inversion_info minv_vector_gcr_var_precond(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, void (*precond_matrix_vector)(complex<double>*,complex<double>*,int,void*), void* precond_info);

// Pre-implemented preconditioners: 

// Preconditioning function.
void identity_preconditioner(double* lhs, double* rhs, int size, void* extra_data); 

void identity_preconditioner(complex<double>* lhs, complex<double>* rhs, int size, void* extra_data);

// Preconditioning function.
void minres_preconditioner(double* lhs, double* rhs, int size, void* extra_data); 

void minres_preconditioner(complex<double>* lhs, complex<double>* rhs, int size, void* extra_data); 

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

#endif // define ESW_INVERTER_PRECOND

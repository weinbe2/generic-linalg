// Wed Jan  6 09:49:39 EST 2016
// Evan S Weinberg
// Include file for generic inverters.

#ifndef ESW_INVERTER_PRECOND
#define ESW_INVERTER_PRECOND

// Uncomment to print warnings. 
//#define VERBOSE_WARN

#include "generic_inverters.h"
#include "verbosity.h"

// Preconditioned Conjugate Gradient:
inversion_info minv_vector_cg_precond(double  *phi, double  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(double*,double*,void*), void* extra_info, void (*precond_matrix_vector)(double*,double*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verbosity = 0);

inversion_info minv_vector_cg_precond(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, void (*precond_matrix_vector)(complex<double>*,complex<double>*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verbosity = 0);

// Flexibly Preconditioned Conjugate Gradient:
inversion_info minv_vector_cg_flex_precond(double  *phi, double  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(double*,double*,void*), void* extra_info, void (*precond_matrix_vector)(double*,double*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verbosity = 0);

inversion_info minv_vector_cg_flex_precond(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, void (*precond_matrix_vector)(complex<double>*,complex<double>*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verbosity = 0);

// Restarted Flexibly Preconditioned Conjugate Gradient:
inversion_info minv_vector_cg_flex_precond_restart(double  *phi, double  *phi0, int size, int max_iter, double eps, int restart_freq, void (*matrix_vector)(double*,double*,void*), void* extra_info, void (*precond_matrix_vector)(double*,double*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verbosity = 0);

inversion_info minv_vector_cg_flex_precond_restart(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, int restart_freq, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, void (*precond_matrix_vector)(complex<double>*,complex<double>*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verbosity = 0);

// Variably Preconditioned GCR:
inversion_info minv_vector_gcr_var_precond(double  *phi, double  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(double*,double*,void*), void* extra_info, void (*precond_matrix_vector)(double*,double*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verbosity = 0);

inversion_info minv_vector_gcr_var_precond(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, void (*precond_matrix_vector)(complex<double>*,complex<double>*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verbosity = 0);

// Restarted Variably Preconditioned GCR:
inversion_info minv_vector_gcr_var_precond_restart(double  *phi, double  *phi0, int size, int max_iter, double eps, int restart_freq, void (*matrix_vector)(double*,double*,void*), void* extra_info, void (*precond_matrix_vector)(double*,double*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verbosity = 0);

inversion_info minv_vector_gcr_var_precond_restart(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, int restart_freq, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, void (*precond_matrix_vector)(complex<double>*,complex<double>*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verbosity = 0);

// Pre-implemented preconditioners: 

// Preconditioning function.
void identity_preconditioner(double* lhs, double* rhs, int size, void* extra_data, inversion_verbose_struct* verb = 0); 

void identity_preconditioner(complex<double>* lhs, complex<double>* rhs, int size, void* extra_data, inversion_verbose_struct* verb = 0);

/* BEGIN MINRES PRECONDITIONER */

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

/* END MINRES PRECONDITIONER */

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

#endif // define ESW_INVERTER_PRECOND

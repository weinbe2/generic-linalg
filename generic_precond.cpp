// Fri Apr 29 11:54:22 EDT 2016
// Evan S Weinberg
// C++ file for CG inverter.

// To do:
// 1. Template to support float, double. 

#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <complex>

#include "generic_inverters_precond.h"
#include "generic_vector.h"

using namespace std;

// Preconditioning function.
void identity_preconditioner(double* lhs, double* rhs, int size, void* extra_data)
{
    int i = 0; 
    
    for (i = 0; i < size; i++)
    {
        lhs[i] = rhs[i];
    }
}

void identity_preconditioner(complex<double>* lhs, complex<double>* rhs, int size, void* extra_data)
{
    int i = 0; 
    
    for (i = 0; i < size; i++)
    {
        lhs[i] = rhs[i];
    }
}

// Minres preconditioning function. 
void minres_preconditioner(double* lhs, double* rhs, int size, void* extra_data)
{
    minres_precond_struct_real* mps = (minres_precond_struct_real*)extra_data; 
    
    // Run mps->nstep iterations of minres. 
    minv_vector_minres(lhs, rhs, size, mps->n_step, mps->rel_res, mps->matrix_vector, mps->matrix_extra_data);
}

void minres_preconditioner(complex<double>* lhs, complex<double>* rhs, int size, void* extra_data)
{
    minres_precond_struct_complex* mps = (minres_precond_struct_complex*)extra_data; 
    
    // Run mps->nstep iterations of minres. 
    minv_vector_minres(lhs, rhs, size, mps->n_step, mps->rel_res, mps->matrix_vector, mps->matrix_extra_data);
}

// GCR preconditioning function. 
void gcr_preconditioner(double* lhs, double* rhs, int size, void* extra_data)
{
    gcr_precond_struct_real* gps = (gcr_precond_struct_real*)extra_data; 
    
    // Run mps->nstep iterations of minres. 
    minv_vector_gcr(lhs, rhs, size, gps->n_step, gps->rel_res, gps->matrix_vector, gps->matrix_extra_data);
}

void gcr_preconditioner(complex<double>* lhs, complex<double>* rhs, int size, void* extra_data)
{
    gcr_precond_struct_complex* gps = (gcr_precond_struct_complex*)extra_data; 
    
    // Run mps->nstep iterations of minres. 
    minv_vector_gcr(lhs, rhs, size, gps->n_step, gps->rel_res, gps->matrix_vector, gps->matrix_extra_data);
}


// Wed Oct 26 11:47:24 EDT 2016
// Evan S Weinberg
// C++ file for calling a generic preconditioned inverter.

// To do:
// 1. Template to support float, double. 

#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <complex>

#include "generic_inverters_precond.h"

using namespace std;

inversion_info minv_preconditioned(double* lhs, double* rhs, int size, minv_inverter_precond type, minv_inverter_precond_params& params, void (*matrix_vector)(double*,double*,void*), void* extra_info, void (*precond_matrix_vector)(double*,double*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verb)
{
  if (type == MINV_PRE_CG)
  {
    // There's no restarted preconditioned CG. There probably should be.
    return minv_vector_cg_precond(lhs, rhs, size, params.max_iters, params.tol, matrix_vector, extra_info, precond_matrix_vector, precond_info, verb);
  }
  else if (type == MINV_PRE_FPCG)
  {
    if (params.restart)
    {
      return minv_vector_cg_flex_precond_restart(lhs, rhs, size, params.max_iters, params.tol, params.restart_freq, matrix_vector, extra_info, precond_matrix_vector, precond_info, verb);
    }
    else
    {
      return minv_vector_cg_flex_precond(lhs, rhs, size, params.max_iters, params.tol, matrix_vector, extra_info, precond_matrix_vector, precond_info, verb);
    }
  }
  else if (type == MINV_PRE_VPGCR)
  {
    if (params.restart) 
    {
      return minv_vector_gcr_var_precond_restart(lhs, rhs, size, params.max_iters, params.tol, params.restart_freq, matrix_vector, extra_info, precond_matrix_vector, precond_info, verb);
    }
    else
    {
      return minv_vector_gcr_var_precond(lhs, rhs, size, params.max_iters, params.tol, matrix_vector, extra_info, precond_matrix_vector, precond_info, verb);
    }
  }
  else if (type == MINV_PRE_BICGSTAB)
  {
    if (params.restart)
    {
      return minv_vector_bicgstab_precond_restart(lhs, rhs, size, params.max_iters, params.tol, params.restart_freq, matrix_vector, extra_info, precond_matrix_vector, precond_info, verb);
    }
    else
    {
      return minv_vector_bicgstab_precond(lhs, rhs, size, params.max_iters, params.tol, matrix_vector, extra_info, precond_matrix_vector, precond_info, verb);
    }
  }
  else if (type == MINV_PRE_INVALID)
  {
		return inversion_info();
  }
  
	return inversion_info();
}

inversion_info minv_preconditioned(complex<double>* lhs, complex<double>* rhs, int size, minv_inverter_precond type, minv_inverter_precond_params& params, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, void (*precond_matrix_vector)(complex<double>*,complex<double>*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verb)
{
  if (type == MINV_PRE_CG)
  {
    // There's no restarted preconditioned CG. There probably should be.
    return minv_vector_cg_precond(lhs, rhs, size, params.max_iters, params.tol, matrix_vector, extra_info, precond_matrix_vector, precond_info, verb);
  }
  else if (type == MINV_PRE_FPCG)
  {
    if (params.restart)
    {
      return minv_vector_cg_flex_precond_restart(lhs, rhs, size, params.max_iters, params.tol, params.restart_freq, matrix_vector, extra_info, precond_matrix_vector, precond_info, verb);
    }
    else
    {
      return minv_vector_cg_flex_precond(lhs, rhs, size, params.max_iters, params.tol, matrix_vector, extra_info, precond_matrix_vector, precond_info, verb);
    }
  }
  else if (type == MINV_PRE_VPGCR)
  {
    if (params.restart) 
    {
      return minv_vector_gcr_var_precond_restart(lhs, rhs, size, params.max_iters, params.tol, params.restart_freq, matrix_vector, extra_info, precond_matrix_vector, precond_info, verb);
    }
    else
    {
      return minv_vector_gcr_var_precond(lhs, rhs, size, params.max_iters, params.tol, matrix_vector, extra_info, precond_matrix_vector, precond_info, verb);
    }
  }
  else if (type == MINV_PRE_BICGSTAB)
  {
    if (params.restart)
    {
      return minv_vector_bicgstab_precond_restart(lhs, rhs, size, params.max_iters, params.tol, params.restart_freq, matrix_vector, extra_info, precond_matrix_vector, precond_info, verb);
    }
    else
    {
      return minv_vector_bicgstab_precond(lhs, rhs, size, params.max_iters, params.tol, matrix_vector, extra_info, precond_matrix_vector, precond_info, verb);
    }
  }
  else if (type == MINV_PRE_INVALID)
  {
		return inversion_info();
  }
  
	return inversion_info();
}



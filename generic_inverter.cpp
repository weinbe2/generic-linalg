// Wed Oct 26 11:47:24 EDT 2016
// Evan S Weinberg
// C++ file for calling a generic inverter.

// To do:
// 1. Template to support float, double. 

#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <complex>

#include "generic_inverters.h"

using namespace std;

inversion_info minv_unpreconditioned(double* lhs, double* rhs, int size, minv_inverter type, minv_inverter_params& params, void (*matrix_vector)(double*,double*,void*), void* extra_info, inversion_verbose_struct* verb)
{
  if (type == MINV_CG)
  {
   	if (params.restart) // why would you do this I don't know it's CG come on
    {
      return minv_vector_cg_restart(lhs, rhs, size, params.max_iters, params.tol, params.restart_freq, matrix_vector, extra_info, verb);
    }
    else
    {
      return minv_vector_cg(lhs, rhs, size, params.max_iters, params.tol, matrix_vector, extra_info, verb);
    }
  }
  else if (type == MINV_CR)
  {
    if (params.restart)
    {
      return minv_vector_cr_restart(lhs, rhs, size, params.max_iters, params.tol, params.restart_freq, matrix_vector, extra_info, verb);
    }
    else
    {
      return minv_vector_cr(lhs, rhs, size, params.max_iters, params.tol, matrix_vector, extra_info, verb);
    }
  }
  else if (type == MINV_GCR)
  {
    if (params.restart) // why would you do this I don't know it's CG come on
    {
      return minv_vector_gcr_restart(lhs, rhs, size, params.max_iters, params.tol, params.restart_freq, matrix_vector, extra_info, verb);
    }
    else
    {
      return minv_vector_gcr(lhs, rhs, size, params.max_iters, params.tol, matrix_vector, extra_info, verb);
    }
  }
  else if (type == MINV_BICGSTAB)
  {
    if (params.restart)
    {
      return minv_vector_bicgstab_restart(lhs, rhs, size, params.max_iters, params.tol, params.restart_freq, matrix_vector, extra_info, verb);
    }
    else
    {
      return minv_vector_bicgstab(lhs, rhs, size, params.max_iters, params.tol, matrix_vector, extra_info, verb);
    }
  }
  else if (type == MINV_BICGSTAB_L)
  {
    if (params.restart)
    {
      return minv_vector_bicgstab_l_restart(lhs, rhs, size, params.max_iters, params.tol, params.restart_freq, params.bicgstabl_l, matrix_vector, extra_info, verb);
    }
    else
    {
      return minv_vector_bicgstab_l(lhs, rhs, size, params.max_iters, params.tol, params.bicgstabl_l, matrix_vector, extra_info, verb);
    }
  }
  else if (type == MINV_GMRES)
  {
    if (params.restart)
    {
      return minv_vector_gmres_restart(lhs, rhs, size, params.max_iters, params.tol, params.restart_freq, matrix_vector, extra_info, verb);
    }
    else
    {
      return minv_vector_gmres(lhs, rhs, size, params.max_iters, params.tol, matrix_vector, extra_info, verb);
    }
  }
  else if (type == MINV_SOR)
  {
    // Restarting doesn't make sense for SOR. 
			return minv_vector_sor(lhs, rhs, size, params.max_iters, params.tol, params.sor_omega, matrix_vector, extra_info, verb);
  }
  else if (type == MINV_MINRES)
  {
    // Restarting doesn't make sense for MinRes. 
    return minv_vector_minres(lhs, rhs, size, params.max_iters, params.tol, params.minres_omega, matrix_vector, extra_info, verb);
  }
  else if (type == MINV_INVALID)
  {
		return inversion_info();
  }
  
	return inversion_info();
}

inversion_info minv_unpreconditioned(complex<double>* lhs, complex<double>* rhs, int size, minv_inverter type, minv_inverter_params& params, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verb)
{
  if (type == MINV_CG)
  {
   	if (params.restart) // why would you do this I don't know it's CG come on
    {
      return minv_vector_cg_restart(lhs, rhs, size, params.max_iters, params.tol, params.restart_freq, matrix_vector, extra_info, verb);
    }
    else
    {
      return minv_vector_cg(lhs, rhs, size, params.max_iters, params.tol, matrix_vector, extra_info, verb);
    }
  }
  else if (type == MINV_CR)
  {
    if (params.restart)
    {
      return minv_vector_cr_restart(lhs, rhs, size, params.max_iters, params.tol, params.restart_freq, matrix_vector, extra_info, verb);
    }
    else
    {
      return minv_vector_cr(lhs, rhs, size, params.max_iters, params.tol, matrix_vector, extra_info, verb);
    }
  }
  else if (type == MINV_GCR)
  {
    if (params.restart) // why would you do this I don't know it's CG come on
    {
      return minv_vector_gcr_restart(lhs, rhs, size, params.max_iters, params.tol, params.restart_freq, matrix_vector, extra_info, verb);
    }
    else
    {
      return minv_vector_gcr(lhs, rhs, size, params.max_iters, params.tol, matrix_vector, extra_info, verb);
    }
  }
  else if (type == MINV_BICGSTAB)
  {
    if (params.restart)
    {
      return minv_vector_bicgstab_restart(lhs, rhs, size, params.max_iters, params.tol, params.restart_freq, matrix_vector, extra_info, verb);
    }
    else
    {
      return minv_vector_bicgstab(lhs, rhs, size, params.max_iters, params.tol, matrix_vector, extra_info, verb);
    }
  }
  else if (type == MINV_BICGSTAB_L)
  {
    if (params.restart)
    {
      return minv_vector_bicgstab_l_restart(lhs, rhs, size, params.max_iters, params.tol, params.restart_freq, params.bicgstabl_l, matrix_vector, extra_info, verb);
    }
    else
    {
      return minv_vector_bicgstab_l(lhs, rhs, size, params.max_iters, params.tol, params.bicgstabl_l, matrix_vector, extra_info, verb);
    }
  }
  else if (type == MINV_GMRES)
  {
    if (params.restart)
    {
      return minv_vector_gmres_restart(lhs, rhs, size, params.max_iters, params.tol, params.restart_freq, matrix_vector, extra_info, verb);
    }
    else
    {
      return minv_vector_gmres(lhs, rhs, size, params.max_iters, params.tol, matrix_vector, extra_info, verb);
    }
  }
  else if (type == MINV_SOR)
  {
    // Restarting doesn't make sense for SOR. 
			return minv_vector_sor(lhs, rhs, size, params.max_iters, params.tol, params.sor_omega, matrix_vector, extra_info, verb);
  }
  else if (type == MINV_MINRES)
  {
    // Restarting doesn't make sense for MinRes. 
    return minv_vector_minres(lhs, rhs, size, params.max_iters, params.tol, params.minres_omega, matrix_vector, extra_info, verb);
  }
  else if (type == MINV_INVALID)
  {
		return inversion_info();
  }
  
	return inversion_info();
}



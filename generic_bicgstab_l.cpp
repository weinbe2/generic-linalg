// Fri Apr 29 11:54:22 EDT 2016
// Evan S Weinberg
// C++ file for BiCGStab inverter.

// To do:
// 1. Template to support float, double. 

#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <complex>

#include "generic_traits.h"
#include "generic_vector.h"

#include "generic_bicgstab.h"

using namespace std;

// Solves lhs = A^(-1) rhs using bicgstab-l.
// Defined in the paper "BICGSTAB(L) for linear equations involving unsymmetric matrices with complex spectrum"
// G. Sleijpen, D. Fokkema, 1993. 
// Based on Kate Clark's implementation in CPS, src file src/util/dirac_op/d_op_wilson_types/bicgstab.C
inversion_info minv_vector_bicgstab_l(double  *phi, double  *phi0, int size, int max_iter, double eps, int l, void (*matrix_vector)(double*,double*,void*), void* extra_info, inversion_verbose_struct* verb)
{

  // Initialize vectors.
  double *r0, **r, **u;
  double rho0, rho1, alpha, omega, beta;
  double *sigma, *gamma, *gamma_prime, *gamma_prime_prime, **tau; 
  double rsq, bsqrt, truersq; 
  int k,i,j,h;
  inversion_info invif;
  
  // Prepare verbosity.
  stringstream ss;
  ss << "BiCGStab-" << l;

  // Allocate memory.
  r0 = new double[size];
  r = new double*[l+1];
  u = new double*[l+1];
  for (i = 0; i < l+1; i++)
  {
    r[i] = new double[size];
    u[i] = new double[size];
  }
  
  sigma = new double[l+1];
  gamma = new double[l+1];
  gamma_prime = new double[l+1];
  gamma_prime_prime = new double[l+1];
  tau = new double*[l+1];
  for (i = 0; i < l+1; i++)
  {
    tau[i] = new double[l+1];
  }
  
  // Zero vectors. 
  zero<double>(r0, size);
  for (i = 0; i < l+1; i++)
  {
    zero<double>(r[i], size);
    zero<double>(u[i], size);
  }
  
  zero<double>(sigma, l+1);
  zero<double>(gamma, l+1);
  zero<double>(gamma_prime, l+1);
  zero<double>(gamma_prime_prime, l+1);
  for (i = 0; i < l+1; i++)
  {
    zero<double>(tau[i], l+1);
  }

  // Initialize values.
  rsq = 0.0; bsqrt = 0.0; truersq = 0.0;
  rho0 = 1;
  alpha = 0;
  omega = 1; 
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq<double>(phi0, size));

  // 1. r[0]0 = b - Ax. r0 = r[0]. sigma[0] = ||r[0]||^2
  // Take advantage of initial guess in phi. Use u[0] as a tmp.
  zero<double>(u[0], size); 
  (*matrix_vector)(u[0], phi, extra_info); invif.ops_count++;
  for(i = 0; i<size; i++) {
    r[0][i] = phi0[i] - u[0][i]; // 1. r0 = b-Ax0
  }
  sigma[0] = norm2sq<double>(r[0], size);
  copy<double>(r0, r[0], size);
  zero<double>(u[0], size);
  
  
  // 2. iterate till convergence
  // This was written with a bit of a headache, so my conventions in my comments are
  // all over the place... sorry about that...
  for(k = 0; k < max_iter; k+=l) {
    
    // rho0 = -omega*rho0;
    rho0 *= -omega; 
    
    // BiCG part.
    for (j = 0; j < l; j++)
    {
      // rho1 = <r0, r_j>, beta = alpha*rho1/rho0, rho0 = rho1
      rho1 = dot<double>(r0, r[j], size);
      beta = alpha*rho1/rho0;
      rho0 = rho1;
      // for i = 0 .. j, u[i] = r[i] - beta*u[i];
      for (i = 0; i <= j; i++)
      {
        for (h = 0; h < size; h++)
        {
          u[i][h] = r[i][h] - beta*u[i][h];
        }
      }
      // u[j+1] = A u[j];
      zero<double>(u[j+1], size); 
      (*matrix_vector)(u[j+1], u[j], extra_info); invif.ops_count++;
      
      // alpha = rho0/<r0, u[j+1]>
      alpha = rho0/dot<double>(r0, u[j+1], size);
      
      // for i = 0 .. j, r[i] = r[i] - alpha u[i+1]
      for (i = 0; i <= j; i++)
      {
        for (h = 0; h < size; h++)
        {
          r[i][h] = r[i][h] - alpha*u[i+1][h];
        }
      }
      
      // r[j+1] = A r[j], x = x + alpha*u[0]
      (*matrix_vector)(r[j+1], r[j], extra_info); invif.ops_count++;
      for (h = 0; h < size; h++)
      {
        phi[h] = phi[h] + alpha*u[0][h];
      }
    } // End BiCG part.
    
    // MR part. Really just modified Gram-Schmidt.
    // I could probably write this in terms of the normalize, orthonormalize functions I have, but eh.
    // The algorithm definition uses the byproducts of the Gram-Schmidt to update x, etc,
    // Don't want to disentangle it at the moment,
    for (j = 1; j <= l; j++)
    {
      for (i = 1; i < j; i++)
      {
        // tau_ij = <r_i,r_j>/sigma_i
        tau[i][j] = dot<double>(r[i], r[j], size)/sigma[i];
        
        // r_j = r_j - tau_ij r[i];
        for (h = 0; h < size; h++)
        {
          r[j][h] = r[j][h] - tau[i][j]*r[i][h];
        }
      }
        
      // sigma_j = r_j^2, gamma'_j = <r_0, r_j>/sigma_j
      sigma[j] = norm2sq<double>(r[j], size);
      gamma_prime[j] = dot<double>(r[j],r[0], size)/sigma[j];
    }
        
    // gamma[l] = gamma'_l, omega = gamma[l]
    gamma[l] = gamma_prime[l];
    omega = gamma[l];
    
    // gamma = T^(-1) gamma_prime. Check paper for defn of T.
    for (j = l-1; j > 0; j--)
    {
      // Internal def: gamma[j] = gamma'_j - \sum_{i = j+1 to l} tau_ji gamma_i
      gamma[j] = gamma_prime[j];
      for (i = j+1; i <= l; i++)
      {
        gamma[j] = gamma[j] - tau[j][i]*gamma[i];
      }
    }
    
    // gamma'' = T S gamma. Check paper for defn of S.
    for (j = 1; j < l; j++)
    {
      gamma_prime_prime[j] = gamma[j+1];
      for (i = j+1; i < l; i++)
      {
        gamma_prime_prime[j] = gamma_prime_prime[j] + tau[j][i]*gamma[i+1];
      }
    }
    
    // Update x, r, u.
    // x = x + gamma_1 r_0, r_0 = r_0 - gamma'_l r_l, u_0 = u_0 - gamma_l u_l
    for (h = 0; h < size; h++)
    {
      phi[h] = phi[h] + gamma[1]*r[0][h];
      u[0][h] = u[0][h] - gamma[l]*u[l][h];
      r[0][h] = r[0][h] - gamma_prime[l]*r[l][h];
    }
           
    // for j = 1 .. l-1: u[0] -= gamma_j u[j], phi += gamma''_j r[j], r[0] -= gamma'_j r[j].
    for (j = 1; j < l; j++)
    {
      for (h = 0; h < size; h++)
      {
        u[0][h] = u[0][h] - gamma[j]*u[j][h];
        phi[h] = phi[h] + gamma_prime_prime[j]*r[j][h];
        r[0][h] = r[0][h] - gamma_prime[j]*r[j][h];
      }
    }
    
    // sigma[0] = r_0^2. This is rsq in my other codes. 
    sigma[0] = norm2sq<double>(r[0], size);
    print_verbosity_resid(verb, ss.str(), k+l, invif.ops_count, sqrt(sigma[0])/bsqrt);
    
    // Check for convergence.
    if (sqrt(sigma[0]) < eps*bsqrt)
    {
      rsq = sigma[0];
      break;
    }
  }
  
  // Weird counting conventions give me headaches.
  if(k >= max_iter-1) {
    //printf("CG: Failed to converge iter = %d, rsq = %e\n", k,rsq);
    invif.success = false;
  }
  else
  {
     //printf("CG: Converged in %d iterations.\n", k);
     invif.success = true;
  }
	k++; 
  
  // Check the true residual. Use u[0] as a tmp.
  (*matrix_vector)(u[0],phi,extra_info); invif.ops_count++; 
  truersq = diffnorm2sq<double>(u[0], phi0, size);
  
  // Free all the things!
  delete[] r0;
  for (i = 0; i < l+1; i++)
  {
    delete[] r[i];
    delete[] u[i];
  }
  delete[] r;
  delete[] u;
  
  delete[] sigma; 
  delete[] gamma;
  delete[] gamma_prime;
  delete[] gamma_prime_prime;
  for (i = 0; i < l+1; i++)
  {
    delete[] tau[i];
  }
  delete[] tau;
  
  print_verbosity_summary(verb, ss.str(), invif.success, k, invif.ops_count, sqrt(truersq)/bsqrt);

  //  printf("# CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = ss.str();
  return invif; // Convergence 

} 

// Performs BiCGStab-l(restart_freq) with restarts when restart_freq is hit.
inversion_info minv_vector_bicgstab_l_restart(double  *phi, double  *phi0, int size, int max_iter, double res, int restart_freq, int l, void (*matrix_vector)(double*,double*,void*), void* extra_info, inversion_verbose_struct* verb)
{
  int iter; // counts total number of iterations.
  int ops_count;
  inversion_info invif;
  double bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  inversion_verbose_struct verb_rest;
  shuffle_verbosity_restart(&verb_rest, verb);
  
  stringstream ss;
  ss << "BiCGStab-" << l << "(" << restart_freq << ")";
  
  iter = 0; ops_count = 0; 
  do
  {
    invif = minv_vector_bicgstab_l(phi, phi0, size, restart_freq, res, l, matrix_vector, extra_info, &verb_rest);
    iter += invif.iter; ops_count += invif.ops_count;
    
    print_verbosity_restart(verb, ss.str(), iter, ops_count, sqrt(invif.resSq)/bsqrt);
  }
  while (iter < max_iter && invif.success == false && sqrt(invif.resSq)/bsqrt > res);
  
  invif.iter = iter; invif.ops_count = ops_count; 
  
  print_verbosity_summary(verb, ss.str(), invif.success, iter, invif.ops_count, sqrt(invif.resSq)/bsqrt);
  
  invif.name = ss.str();
  // invif.resSq is good.
  if (sqrt(invif.resSq)/bsqrt > res)
  {
    invif.success = false;
  }
  else
  {
    invif.success = true;
  }
  
  return invif;
}

// Solves lhs = A^(-1) rhs using bicgstab-l.
// Defined in the paper "BICGSTAB(L) for linear equations involving unsymmetric matrices with complex spectrum"
// G. Sleijpen, D. Fokkema, 1993. 
// Based on Kate Clark's implementation in CPS, src file src/util/dirac_op/d_op_wilson_types/bicgstab.C
inversion_info minv_vector_bicgstab_l(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, int l, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verb)
{

  // Initialize vectors.
  complex<double> *r0, **r, **u;
  complex<double> rho0, rho1, alpha, omega, beta;
  complex<double> *gamma, *gamma_prime, *gamma_prime_prime, **tau; 
  double *sigma, rsq, bsqrt, truersq; 
  int k,i,j,h;
  inversion_info invif;
  
  // Prepare verbosity.
  stringstream ss;
  ss << "BiCGStab-" << l;

  // Allocate memory.
  r0 = new complex<double>[size];
  r = new complex<double>*[l+1];
  u = new complex<double>*[l+1];
  for (i = 0; i < l+1; i++)
  {
    r[i] = new complex<double>[size];
    u[i] = new complex<double>[size];
  }
  
  sigma = new double[l+1];
  gamma = new complex<double>[l+1];
  gamma_prime = new complex<double>[l+1];
  gamma_prime_prime = new complex<double>[l+1];
  tau = new complex<double>*[l+1];
  for (i = 0; i < l+1; i++)
  {
    tau[i] = new complex<double>[l+1];
  }
  
  // Zero vectors. 
  zero<double>(r0, size);
  for (i = 0; i < l+1; i++)
  {
    zero<double>(r[i], size);
    zero<double>(u[i], size);
  }
  
  zero<double>(sigma, l+1);
  zero<double>(gamma, l+1);
  zero<double>(gamma_prime, l+1);
  zero<double>(gamma_prime_prime, l+1);
  for (i = 0; i < l+1; i++)
  {
    zero<double>(tau[i], l+1);
  }

  // Initialize values.
  rsq = 0.0; bsqrt = 0.0; truersq = 0.0;
  rho0 = 1;
  alpha = 0;
  omega = 1; 
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq<double>(phi0, size));

  // 1. r[0]0 = b - Ax. r0 = r[0]. sigma[0] = ||r[0]||^2
  // Take advantage of initial guess in phi. Use u[0] as a tmp.
  zero<double>(u[0], size); 
  (*matrix_vector)(u[0], phi, extra_info); invif.ops_count++;
  for(i = 0; i<size; i++) {
    r[0][i] = phi0[i] - u[0][i]; // 1. r0 = b-Ax0
  }
  sigma[0] = norm2sq<double>(r[0], size);
  copy<double>(r0, r[0], size);
  zero<double>(u[0], size);
  
  
  // 2. iterate till convergence
  // This was written with a bit of a headache, so my conventions in my comments are
  // all over the place... sorry about that...
  for(k = 0; k < max_iter; k+=l) {
    
    // rho0 = -omega*rho0;
    rho0 *= -omega; 
    
    // BiCG part.
    for (j = 0; j < l; j++)
    {
      // rho1 = <r0, r_j>, beta = alpha*rho1/rho0, rho0 = rho1
      rho1 = dot<double>(r0, r[j], size);
      beta = alpha*rho1/rho0;
      rho0 = rho1;
      // for i = 0 .. j, u[i] = r[i] - beta*u[i];
      for (i = 0; i <= j; i++)
      {
        for (h = 0; h < size; h++)
        {
          u[i][h] = r[i][h] - beta*u[i][h];
        }
      }
      // u[j+1] = A u[j];
      zero<double>(u[j+1], size); 
      (*matrix_vector)(u[j+1], u[j], extra_info); invif.ops_count++;
      
      // alpha = rho0/<r0, u[j+1]>
      alpha = rho0/dot<double>(r0, u[j+1], size);
      
      // for i = 0 .. j, r[i] = r[i] - alpha u[i+1]
      for (i = 0; i <= j; i++)
      {
        for (h = 0; h < size; h++)
        {
          r[i][h] = r[i][h] - alpha*u[i+1][h];
        }
      }
      
      // r[j+1] = A r[j], x = x + alpha*u[0]
      (*matrix_vector)(r[j+1], r[j], extra_info); invif.ops_count++;
      for (h = 0; h < size; h++)
      {
        phi[h] = phi[h] + alpha*u[0][h];
      }
    } // End BiCG part.
    
    // MR part. Really just modified Gram-Schmidt.
    // I could probably write this in terms of the normalize, orthonormalize functions I have, but eh.
    // The algorithm definition uses the byproducts of the Gram-Schmidt to update x, etc,
    // Don't want to disentangle it at the moment,
    for (j = 1; j <= l; j++)
    {
      for (i = 1; i < j; i++)
      {
        // tau_ij = <r_i,r_j>/sigma_i
        tau[i][j] = dot<double>(r[i], r[j], size)/sigma[i];
        
        // r_j = r_j - tau_ij r[i];
        for (h = 0; h < size; h++)
        {
          r[j][h] = r[j][h] - tau[i][j]*r[i][h];
        }
      }
        
      // sigma_j = r_j^2, gamma'_j = <r_0, r_j>/sigma_j
      sigma[j] = norm2sq<double>(r[j], size);
      gamma_prime[j] = dot<double>(r[j],r[0], size)/sigma[j];
    }
        
    // gamma[l] = gamma'_l, omega = gamma[l]
    gamma[l] = gamma_prime[l];
    omega = gamma[l];
    
    // gamma = T^(-1) gamma_prime. Check paper for defn of T.
    for (j = l-1; j > 0; j--)
    {
      // Internal def: gamma[j] = gamma'_j - \sum_{i = j+1 to l} tau_ji gamma_i
      gamma[j] = gamma_prime[j];
      for (i = j+1; i <= l; i++)
      {
        gamma[j] = gamma[j] - tau[j][i]*gamma[i];
      }
    }
    
    // gamma'' = T S gamma. Check paper for defn of S.
    for (j = 1; j < l; j++)
    {
      gamma_prime_prime[j] = gamma[j+1];
      for (i = j+1; i < l; i++)
      {
        gamma_prime_prime[j] = gamma_prime_prime[j] + tau[j][i]*gamma[i+1];
      }
    }
    
    // Update x, r, u.
    // x = x + gamma_1 r_0, r_0 = r_0 - gamma'_l r_l, u_0 = u_0 - gamma_l u_l
    for (h = 0; h < size; h++)
    {
      phi[h] = phi[h] + gamma[1]*r[0][h];
      u[0][h] = u[0][h] - gamma[l]*u[l][h];
      r[0][h] = r[0][h] - gamma_prime[l]*r[l][h];
    }
           
    // for j = 1 .. l-1: u[0] -= gamma_j u[j], phi += gamma''_j r[j], r[0] -= gamma'_j r[j].
    for (j = 1; j < l; j++)
    {
      for (h = 0; h < size; h++)
      {
        u[0][h] = u[0][h] - gamma[j]*u[j][h];
        phi[h] = phi[h] + gamma_prime_prime[j]*r[j][h];
        r[0][h] = r[0][h] - gamma_prime[j]*r[j][h];
      }
    }
    
    // sigma[0] = r_0^2. This is rsq in my other codes. 
    sigma[0] = norm2sq<double>(r[0], size);
    print_verbosity_resid(verb, ss.str(), k+l, invif.ops_count, sqrt(sigma[0])/bsqrt);
    
    // Check for convergence.
    if (sqrt(sigma[0]) < eps*bsqrt)
    {
      rsq = sigma[0];
      break;
    }
  }
  
  // Weird counting conventions give me headaches.
  if(k >= max_iter-1) {
    //printf("CG: Failed to converge iter = %d, rsq = %e\n", k,rsq);
    invif.success = false;
  }
  else
  {
     //printf("CG: Converged in %d iterations.\n", k);
     invif.success = true;
  }
	k++; 
  
  // Check the true residual. Use u[0] as a tmp.
  (*matrix_vector)(u[0],phi,extra_info); invif.ops_count++; 
  truersq = diffnorm2sq<double>(u[0], phi0, size);
  
  // Free all the things!
  delete[] r0;
  for (i = 0; i < l+1; i++)
  {
    delete[] r[i];
    delete[] u[i];
  }
  delete[] r;
  delete[] u;
  
  delete[] sigma; 
  delete[] gamma;
  delete[] gamma_prime;
  delete[] gamma_prime_prime;
  for (i = 0; i < l+1; i++)
  {
    delete[] tau[i];
  }
  delete[] tau;
  
  print_verbosity_summary(verb, ss.str(), invif.success, k, invif.ops_count, sqrt(truersq)/bsqrt);

  //  printf("# CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = ss.str();
  return invif; // Convergence 

} 

// Performs BiCGStab-l(restart_freq) with restarts when restart_freq is hit.
inversion_info minv_vector_bicgstab_l_restart(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double res, int restart_freq, int l, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verb)
{
  int iter; // counts total number of iterations.
  int ops_count;
  inversion_info invif;
  double bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  inversion_verbose_struct verb_rest;
  shuffle_verbosity_restart(&verb_rest, verb);
  
  stringstream ss;
  ss << "BiCGStab-" << l << "(" << restart_freq << ")";
  
  iter = 0; ops_count = 0; 
  do
  {
    invif = minv_vector_bicgstab_l(phi, phi0, size, restart_freq, res, l, matrix_vector, extra_info, &verb_rest);
    iter += invif.iter; ops_count += invif.ops_count;
    
    print_verbosity_restart(verb, ss.str(), iter, ops_count, sqrt(invif.resSq)/bsqrt);
  }
  while (iter < max_iter && invif.success == false && sqrt(invif.resSq)/bsqrt > res);
  
  invif.iter = iter; invif.ops_count = ops_count; 
  
  print_verbosity_summary(verb, ss.str(), invif.success, iter, invif.ops_count, sqrt(invif.resSq)/bsqrt);
  
  invif.name = ss.str();
  // invif.resSq is good.
  if (sqrt(invif.resSq)/bsqrt > res)
  {
    invif.success = false;
  }
  else
  {
    invif.success = true;
  }
  
  return invif;
}


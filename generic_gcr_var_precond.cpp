// Fri Apr 29 11:54:22 EDT 2016
// Evan S Weinberg
// C++ file for variably preconditioned GCR inverter.
// Based on VPGCR defined in:
//   "A VARIABLE PRECONDITIONING USING THE SOR METHOD FOR GCR-LIKE METHODS"
// And on the equivalent rewrite in:
//   "A variable preconditioned GCR(m) method using the GSOR method for singular and rectangular linear systems"
//   However, this should work for a generic variable preconditioner. 

// To do:
// 1. Template various functions to support float, double,
//    as well as complex< > variants.

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <complex>

#include "generic_vector.h"

#include "generic_gcr_var_precond.h"

using namespace std; 

    

// Solves lhs = A^(-1) rhs using GCR.
// Taken from section 6.9 of Saad, 2nd Edition.
inversion_info minv_vector_gcr_var_precond(double  *phi, double  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(double*,double*,void*), void* extra_info, void (*precond_matrix_vector)(double*,double*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verb)
{
  // Initialize vectors.
  double *x, *r, *Ar, *p, *Ap, *z, *Az; 
  vector<double*> p_store, Ap_store; // GCR requires explicit reorthogonalization against old search vectors.
  double alpha, beta_ij, rsq, bsqrt, truersq;
  int k,i,ii;
  inversion_info invif;
  
  // For preconditioning verbosity.
  inversion_verbose_struct verb_prec;
  shuffle_verbosity_precond(&verb_prec, verb);

  // Allocate memory.
  x = new double[size];
  r = new double[size];
  Ar = new double[size];
  p = new double[size];
  Ap = new double[size];
  z = new double[size];
  Az = new double[size];
  
  // Zero vectors. 
  zero<double>(p, size);  zero<double>(r, size);
  zero<double>(Ap, size); zero<double>(Ar, size);
  zero<double>(z, size);  zero<double>(Az, size); 

  // Initialize values.
  rsq = 0.0; bsqrt = 0.0; truersq = 0.0;
  
  // Copy initial guess into solution.
  copy<double>(x, phi, size);
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  // 1. r_0 = b - Ax_0. x is phi, the initial guess.
  (*matrix_vector)(p, x, extra_info); invif.ops_count++; // Put Ax_0 into p, temp.
  for (i = 0; i < size; i++)
  {
    r[i] = phi0[i] - p[i]; // r_0 = b - Ax_0
  }
  
  // 2. z_0 = M^(-1) r_0.
  (*precond_matrix_vector)(z, r, size, precond_info, &verb_prec); 
  
  // 3. p_0 = z_0, Compute A p_0 (called q_0 in paper)
  copy<double>(p, z, size); 
  (*matrix_vector)(Ap, p, extra_info); invif.ops_count++;

  // iterate until convergence
  for(k = 0; k< max_iter; k++) {
    
    // If we've hit here, push the latest p, Ap to storage.
    p_store.push_back(p);
    Ap_store.push_back(Ap);
    
    // 4. alpha = <Ap_k, r>/<Ap_k, Ap_k>
    alpha = dot<double>(Ap, r, size)/norm2sq<double>(Ap, size);
    
    // 5. x = x + alpha p_k
    // 6. r = r - alpha Ap_k
    for (i = 0; i < size; i++)
    {
      x[i] = x[i] + alpha*p[i];
      r[i] = r[i] - alpha*Ap[i];
    }
    
    // Compute norm.
    rsq = norm2sq<double>(r, size);
    
    print_verbosity_resid(verb, "VPGCR", k+1, invif.ops_count, sqrt(rsq)/bsqrt); 
    
    //printf("Rel residual: %.8e\n", sqrt(rsq)/bsqrt); fflush(stdout);
    
    // Check convergence. 
    if (sqrt(rsq) < eps*bsqrt || k==max_iter-1) {
      //        printf("Final rsq = %g\n", rsqNew);
      break;
    }
    
    // 7. z = M^(-1) r;
    zero<double>(z, size); 
    (*precond_matrix_vector)(z, r, size, precond_info, &verb_prec); 
    
    // 8. Compute Az.
    zero<double>(Az, size);
    (*matrix_vector)(Az, z, extra_info); invif.ops_count++;
    
    // 8. b_ij = -<Az_{j+1}, Ap_i>/<Ap_i, Ap_i> for i = 0, ..., j
    // 9. p_{j+1} = z_{j+1} + sum_i=0^j b_ij p_i
    // 10. Ap_{j+1} = Az_{j+1} + sum_i=0^j b_ij Ap_i
    p = new double[size];  copy<double>(p, z, size);
    Ap = new double[size]; copy<double>(Ap, Az, size);
    for (ii = 0; ii <= k; ii++)
    {
      beta_ij = -dot<double>(Ap_store[ii], Az, size)/norm2sq<double>(Ap_store[ii], size);
      for (i = 0; i < size; i++)
      {
        p[i] += beta_ij*p_store[ii][i];
        Ap[i] += beta_ij*Ap_store[ii][i]; 
      }
    }
  } 
    
  if(k == max_iter-1) {
    //printf("CG: Failed to converge iter = %d, rsq = %e\n", k,rsq);
    invif.success = false;
    //return 0;// Failed convergence 
  }
  else
  {
     invif.success = true;
     //printf("CG: Converged in %d iterations.\n", k);
  }
	k++;
  
  // Check true residual.
  zero<double>(p,size);
  (*matrix_vector)(p,x,extra_info); invif.ops_count++;
  for(i=0; i < size; i++) truersq += (p[i] - phi0[i])*(p[i] - phi0[i]);
  
  // Copy solution into phi.
  copy<double>(phi, x, size);
  
  // Free all the things!
  delete[] x;
  delete[] r;
  delete[] Ar;
  ii = p_store.size();
  for (i = 0; i < ii; i++)
  {
    delete[] p_store[i];
    delete[] Ap_store[i];
  }
  delete[] z; 
  delete[] Az; 

  
  print_verbosity_summary(verb, "VPGCR", invif.success, k, invif.ops_count, sqrt(truersq)/bsqrt);
  
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "Variably Preconditioned GCR";
  return invif; // Convergence 
} 

// Performs VPGCR(restart_freq) with restarts when restart_freq is hit.
// This may be sloppy, but it works.
inversion_info minv_vector_gcr_var_precond_restart(double  *phi, double  *phi0, int size, int max_iter, double res, int restart_freq, void (*matrix_vector)(double*,double*,void*), void* extra_info, void (*precond_matrix_vector)(double*,double*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verb)
{
  int iter; // counts total number of iterations.
  int ops_count; 
  inversion_info invif;
  double bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  stringstream ss;
  ss << "Variably Preconditioned Restarted GCR(" << restart_freq << ")";

  inversion_verbose_struct verb_rest;
  shuffle_verbosity_restart(&verb_rest, verb);
  
  iter = 0; ops_count = 0; 
  do
  {
    invif = minv_vector_gcr_var_precond(phi, phi0, size, restart_freq, res, matrix_vector, extra_info, precond_matrix_vector, precond_info, &verb_rest);
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

inversion_info minv_vector_gcr_var_precond(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, void (*precond_matrix_vector)(complex<double>*,complex<double>*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verb)
{

  // Initialize vectors.
  complex<double> *x, *r, *Ar, *p, *Ap,  *z, *Az; 
  vector<complex<double>*> p_store, Ap_store; // GCR requires explicit reorthogonalization against old search vectors.
  double rsq, bsqrt, truersq;
  complex<double> alpha, beta_ij;
  int k,i,ii;
  inversion_info invif;
  
  // For preconditioning verbosity.
  inversion_verbose_struct verb_prec;
  shuffle_verbosity_precond(&verb_prec, verb);

  // Allocate memory.
  x = new complex<double>[size];
  r = new complex<double>[size];
  Ar = new complex<double>[size];
  p = new complex<double>[size];
  Ap = new complex<double>[size];
  z = new complex<double>[size];
  Az = new complex<double>[size];
  
  // Zero vectors. 
  zero<double>(p, size);  zero<double>(r, size);
  zero<double>(Ap, size); zero<double>(Ar, size);
  zero<double>(z, size);  zero<double>(Az, size); 

  // Initialize values.
  rsq = 0.0; bsqrt = 0.0; truersq = 0.0;
  
  // Copy initial guess into solution.
  copy<double>(x, phi, size);
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  // 1. r_0 = b - Ax_0. x is phi, the initial guess.
  (*matrix_vector)(p, x, extra_info); invif.ops_count++; // Put Ax_0 into p, temp.
  for (i = 0; i < size; i++)
  {
    r[i] = phi0[i] - p[i]; // r_0 = b - Ax_0
  }
  
  // 2. z_0 = M^(-1) r_0.
  (*precond_matrix_vector)(z, r, size, precond_info, &verb_prec); 
  
  // 3. p_0 = z_0, Compute A p_0 (called q_0 in paper)
  copy<double>(p, z, size); 
  (*matrix_vector)(Ap, p, extra_info); invif.ops_count++;

  // iterate until convergence
  for(k = 0; k< max_iter; k++) {
    
    // If we've hit here, push the latest p, Ap to storage.
    p_store.push_back(p);
    Ap_store.push_back(Ap);
    
    // 4. alpha = <Ap_k, r>/<Ap_k, Ap_k>
    alpha = dot<double>(Ap, r, size)/norm2sq<double>(Ap, size);
    
    // 5. x = x + alpha p_k
    // 6. r = r - alpha Ap_k
    for (i = 0; i < size; i++)
    {
      x[i] = x[i] + alpha*p[i];
      r[i] = r[i] - alpha*Ap[i];
    }
    
    // Compute norm.
    rsq = norm2sq<double>(r, size);
    
    print_verbosity_resid(verb, "VPGCR", k+1, invif.ops_count, sqrt(rsq)/bsqrt); 
    
    // Check convergence. 
    if (sqrt(rsq) < eps*bsqrt || k==max_iter-1) {
      //        printf("Final rsq = %g\n", rsqNew);
      break;
    }
    
    // 7. z = M^(-1) r;
    zero<double>(z, size); 
    (*precond_matrix_vector)(z, r, size, precond_info, &verb_prec); 
    
    // 8. Compute Az.
    zero<double>(Az, size);
    (*matrix_vector)(Az, z, extra_info); invif.ops_count++;
    
    // 8. b_ij = -<Az_{j+1}, Ap_i>/<Ap_i, Ap_i> for i = 0, ..., j
    // 9. p_{j+1} = z_{j+1} + sum_i=0^j b_ij p_i
    // 10. Ap_{j+1} = Az_{j+1} + sum_i=0^j b_ij Ap_i
    p = new complex<double>[size];  copy<double>(p, z, size);
    Ap = new complex<double>[size]; copy<double>(Ap, Az, size);
    for (ii = 0; ii <= k; ii++)
    {
      beta_ij = -dot<double>(Ap_store[ii], Az, size)/norm2sq<double>(Ap_store[ii], size);
      for (i = 0; i < size; i++)
      {
        p[i] += beta_ij*p_store[ii][i];
        Ap[i] += beta_ij*Ap_store[ii][i]; 
      }
    }
  } 
    
  if(k == max_iter-1) {
    //printf("CG: Failed to converge iter = %d, rsq = %e\n", k,rsq);
    invif.success = false;
    //return 0;// Failed convergence 
  }
  else
  {
     invif.success = true;
     //printf("CG: Converged in %d iterations.\n", k);
  }
	k++;
  
  // Check true residual.
  zero<double>(p,size);
  (*matrix_vector)(p,x,extra_info); invif.ops_count++;
  for(i=0; i < size; i++) truersq += real(conj(p[i] - phi0[i])*(p[i] - phi0[i]));
  
  // Copy solution into phi.
  copy<double>(phi, x, size);
  
  // Free all the things!
  delete[] x;
  delete[] r;
  delete[] Ar;
  k = p_store.size();
  for (i = 0; i < k; i++)
  {
    delete[] p_store[i];
    delete[] Ap_store[i];
  }
  delete[] z; 
  delete[] Az; 

  
  print_verbosity_summary(verb, "VPGCR", invif.success, k, invif.ops_count, sqrt(truersq)/bsqrt);
  
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "Variably Preconditioned GCR";
  return invif; // Convergence 
} 

// Performs VPGCR(restart_freq) with restarts when restart_freq is hit.
// This may be sloppy, but it works.
inversion_info minv_vector_gcr_var_precond_restart(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double res, int restart_freq, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, void (*precond_matrix_vector)(complex<double>*,complex<double>*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verb)
{
  int iter; // counts total number of iterations.
  int ops_count; 
  inversion_info invif;
  double bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  stringstream ss;
  ss << "Variably Preconditioned Restarted GCR(" << restart_freq << ")";

  inversion_verbose_struct verb_rest;
  shuffle_verbosity_restart(&verb_rest, verb);
  
  iter = 0; ops_count = 0; 
  do
  {
    invif = minv_vector_gcr_var_precond(phi, phi0, size, restart_freq, res, matrix_vector, extra_info, precond_matrix_vector, precond_info, &verb_rest);
    iter += invif.iter; ops_count += invif.ops_count; 
    
    print_verbosity_restart(verb, ss.str(), iter, ops_count, sqrt(invif.resSq)/bsqrt);
  }
  while (iter < max_iter && invif.success == false && sqrt(invif.resSq)/bsqrt > res);
  
  print_verbosity_summary(verb, ss.str(), invif.success, iter, invif.ops_count, sqrt(invif.resSq)/bsqrt);
  
  invif.iter = iter;
  invif.ops_count = ops_count; 
  
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



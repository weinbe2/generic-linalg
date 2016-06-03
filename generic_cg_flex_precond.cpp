// Fri Apr 29 11:54:22 EDT 2016
// Evan S Weinberg
// C++ file for preconditioned CG inverter.

// To do:
// 1. Template to support float, double. 

#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <complex>
#include <vector>

#include "generic_inverters_precond.h"
#include "generic_vector.h"

using namespace std;

// Based on Minv_phi_3d in fem.c
// Solves lhs = A^(-1) rhs using flexibly preconditioned Conjugate gradient. 
// Defined in FLEXIBLE CONJUGATE GRADIENTS by YVAN NOTAY. 
// Currently implemented without truncation. 
// Uses the preconditioner given in precond_matrix_vector. 
inversion_info minv_vector_cg_flex_precond(double  *phi, double  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(double*,double*,void*), void* extra_info, void (*precond_matrix_vector)(double*,double*,int,void*), void* precond_info)
{
//  see http://en.wikipedia.org/wiki/Conjugate_gradient_method

  // Initialize vectors.
  double *r, *p, *Ap, *z, *Az;
  vector<double*> p_store, Ap_store; // Flexible CG requires explicit reorthogonalization against old search vectors. 
  double alpha, beta_ij, rsq, bsqrt, truersq;
  int k,i,ii;
  inversion_info invif;

  // Allocate memory.
  r = new double[size];
  p = new double[size];
  Ap = new double[size];
  z = new double[size];
  Az = new double[size]; 

  // Initialize values.
  rsq = bsqrt = 0.0;

  // Zero vectors;
  zero<double>(r, size); zero<double>(z, size);
  zero<double>(p, size); zero<double>(Ap, size);
  zero<double>(Az, size); 
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  // 1. Compute r = b - Ax
  (*matrix_vector)(p, phi, extra_info);
  for (i = 0; i < size; i++)
  {
    r[i] = phi0[i] - p[i];
  }
  
  // 2. z = M^(-1) r
  (*precond_matrix_vector)(z, r, size, precond_info); 
  
  // 3. p_0 = z_0.
  copy<double>(p, z, size);
  
  // Compute Ap.
  zero<double>(Ap, size);
  (*matrix_vector)(Ap, p, extra_info);

  // iterate till convergence
  for(k = 0; k< max_iter; k++) {
    
    // If we've hit here, push the latest p, Ap to storage.
    p_store.push_back(p); 
    Ap_store.push_back(Ap);
    
    // 4. alpha = z dot r / p A p 
    alpha = dot<double>(p, r, size)/dot<double>(p, Ap, size);

    // 5. phi = phi + alpha p
    // 6. r = r - alpha A p
    for (i = 0; i < size; i++)
    {
      phi[i] = phi[i] + alpha*p[i];
      r[i] = r[i] - alpha*Ap[i];
    }
    
    // Exit if new residual is small enough
    rsq = norm2sq<double>(r, size);

    if (sqrt(rsq) < eps*bsqrt) {
      //        printf("Final rsq = %g\n", rsqNew);
      break;
    }
    
    // 7. z = M^(-1) r
    zero<double>(z, size);
    (*precond_matrix_vector)(z, r, size, precond_info); 
    
    // 8. Compute Az.
    zero<double>(Az, size);
    (*matrix_vector)(Az, z, extra_info); 
    
    // 9. b_ij = -<p_i, z>/<p_i, A p_i>;
    // 10. p_{j+1} = z_{j+1} + sum_i=0^j b_ij p_i
    // 11. Ap_{j+1} = Az_{j+1} + sum_i=0^j b_ij Ap_i
    p = new double[size];  copy<double>(p, z, size);
    Ap = new double[size]; copy<double>(Ap, Az, size);
    for (ii = 0; ii <= k; ii++)
    {
      beta_ij = -dot<double>(Ap_store[ii], z, size)/dot<double>(p_store[ii], Ap_store[ii], size); 
      for (i = 0; i < size; i++)
      {
        p[i] += beta_ij*p_store[ii][i];
        Ap[i] += beta_ij*Ap_store[ii][i];
      }
    }
    
  } 
    
  if(k == max_iter) {
    //printf("CG: Failed to converge iter = %d, rsq = %e\n", k,rsq);
    invif.success = false;
    //return 0;// Failed convergence 
  }
  else
  {
     invif.success = true;
     //printf("CG: Converged in %d iterations.\n", k);
  }
  
  
  zero<double>(Ap, size); 
  (*matrix_vector)(Ap,phi,extra_info);
  for(i=0; i < size; i++) truersq += (Ap[i] - phi0[i])*(Ap[i] - phi0[i]);
  
  // Free all the things!
  delete[] r;
  delete[] p;
  delete[] Ap;
  delete[] z; 
  delete[] Az; 
  ii = p_store.size();
  for (i = 0; i < ii-1; i++) // not sure why it's ii-1...
  {
    delete[] p_store[i];
    delete[] Ap_store[i];
  }

  
  //  printf("# CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "Flexibly Preconditioned CG";
  return invif; // Convergence 

} 


inversion_info minv_vector_cg_flex_precond(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, void (*precond_matrix_vector)(complex<double>*,complex<double>*,int,void*), void* precond_info)
{
  // Initialize vectors.
  complex<double> *r, *p, *Ap, *z, *Az;
  vector<complex<double>*> p_store, Ap_store; // Flexible CG requires explicit reorthogonalization against old search vectors. 
  complex<double> alpha, beta_ij;
  double rsq, bsqrt, truersq;
  int k,i,ii;
  inversion_info invif;

  // Allocate memory.
  r = new complex<double>[size];
  p = new complex<double>[size];
  Ap = new complex<double>[size];
  z = new complex<double>[size];
  Az = new complex<double>[size]; 

  // Initialize values.
  rsq = bsqrt = 0.0;

  // Zero vectors;
  zero<double>(r, size); zero<double>(z, size);
  zero<double>(p, size); zero<double>(Ap, size);
  zero<double>(Az, size); 
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  // 1. Compute r = b - Ax
  (*matrix_vector)(p, phi, extra_info);
  for (i = 0; i < size; i++)
  {
    r[i] = phi0[i] - p[i];
  }
  
  // 2. z = M^(-1) r
  (*precond_matrix_vector)(z, r, size, precond_info); 
  
  // 3. p_0 = z_0.
  copy<double>(p, z, size);
  
  // Compute Ap.
  zero<double>(Ap, size);
  (*matrix_vector)(Ap, p, extra_info);

  // iterate till convergence
  for(k = 0; k< max_iter; k++) {
    
    // If we've hit here, push the latest p, Ap to storage.
    p_store.push_back(p); 
    Ap_store.push_back(Ap);
    
    // 4. alpha = z dot r / p A p 
    alpha = dot<double>(p, r, size)/dot<double>(p, Ap, size);

    // 5. phi = phi + alpha p
    // 6. r = r - alpha A p
    for (i = 0; i < size; i++)
    {
      phi[i] = phi[i] + alpha*p[i];
      r[i] = r[i] - alpha*Ap[i];
    }
    
    // Exit if new residual is small enough
    rsq = norm2sq<double>(r, size);

    if (sqrt(rsq) < eps*bsqrt) {
      //        printf("Final rsq = %g\n", rsqNew);
      break;
    }
    
    // 7. z = M^(-1) r
    zero<double>(z, size);
    (*precond_matrix_vector)(z, r, size, precond_info); 
    
    // 8. Compute Az.
    zero<double>(Az, size);
    (*matrix_vector)(Az, z, extra_info); 
    
    // 9. b_ij = -<p_i, z>/<p_i, A p_i>;
    // 10. p_{j+1} = z_{j+1} + sum_i=0^j b_ij p_i
    // 11. Ap_{j+1} = Az_{j+1} + sum_i=0^j b_ij Ap_i
    p = new complex<double>[size];  copy<double>(p, z, size);
    Ap = new complex<double>[size]; copy<double>(Ap, Az, size);
    for (ii = 0; ii <= k; ii++)
    {
      beta_ij = -dot<double>(Ap_store[ii], z, size)/dot<double>(p_store[ii], Ap_store[ii], size); 
      for (i = 0; i < size; i++)
      {
        p[i] += beta_ij*p_store[ii][i];
        Ap[i] += beta_ij*Ap_store[ii][i];
      }
    }
    
  } 
    
  if(k == max_iter) {
    //printf("CG: Failed to converge iter = %d, rsq = %e\n", k,rsq);
    invif.success = false;
    //return 0;// Failed convergence 
  }
  else
  {
     invif.success = true;
     //printf("CG: Converged in %d iterations.\n", k);
  }
  
  
  zero<double>(Ap, size); 
  (*matrix_vector)(Ap,phi,extra_info);
  for(i=0; i < size; i++) truersq += real(conj(Ap[i] - phi0[i])*(Ap[i] - phi0[i]));
  
  // Free all the things!
  delete[] r;
  delete[] p;
  delete[] Ap;
  delete[] z; 
  delete[] Az; 
  ii = p_store.size();
  for (i = 0; i < ii-1; i++) // not sure why it's ii-1...
  {
    delete[] p_store[i];
    delete[] Ap_store[i];
  }

  
  //  printf("# CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "Flexibly Preconditioned CG";
  return invif; // Convergence 


}

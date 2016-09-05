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

#include "generic_vector.h"

#include "generic_cg_precond.h"

using namespace std;

// Based on Minv_phi_3d in fem.c
// Solves lhs = A^(-1) rhs using Conjugate gradient. 
// Uses the preconditioner given in precond_matrix_vector. 
inversion_info minv_vector_cg_precond(double  *phi, double  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(double*,double*,void*), void* extra_info, void (*precond_matrix_vector)(double*,double*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verb)
{
/// CG solutions to Mphi = b 
//  see http://en.wikipedia.org/wiki/Conjugate_gradient_method

  // Initialize vectors.
  double *r, *p, *Ap, *z;
  double alpha, beta, zdotr, zdotr_new, rsq, bsqrt, truersq;
  int k,i;
  inversion_info invif;
  
  // For preconditioning verbosity.
  inversion_verbose_struct verb_prec;
  shuffle_verbosity_precond(&verb_prec, verb);
    

  // Allocate memory.
  r = new double[size];
  p = new double[size];
  Ap = new double[size];
  z = new double[size];

  // Initialize values.
  rsq = zdotr = zdotr_new = bsqrt = truersq = 0.0;

  // Zero vectors;
  zero<double>(r, size); zero<double>(z, size);
  zero<double>(p, size); zero<double>(Ap, size);
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  // 1. Compute r = b - Ax
  (*matrix_vector)(p, phi, extra_info); invif.ops_count++;
  for (i = 0; i < size; i++)
  {
    r[i] = phi0[i] - p[i];
  }
  
  // 2. z = M^(-1) r
  (*precond_matrix_vector)(z, r, size, precond_info, &verb_prec); 
  
  // 3. p_0 = z_0.
  copy<double>(p, z, size);
  
  // Compute Ap.
  zero<double>(Ap, size);
  (*matrix_vector)(Ap, p, extra_info); invif.ops_count++;
  
  // Compute zdotr.
  zdotr = dot<double>(z, r, size);

  // iterate till convergence
  for(k = 0; k< max_iter; k++) {
    
    // 4. alpha = z dot r / p A p 
    alpha = zdotr/dot<double>(p, Ap, size);

    // 5. phi = phi + alpha p
    // 6. r = r - alpha A p
    for (i = 0; i < size; i++)
    {
      phi[i] = phi[i] + alpha*p[i];
      r[i] = r[i] - alpha*Ap[i];
    }
    
    // Exit if new residual is small enough
    rsq = norm2sq<double>(r, size);
    
    print_verbosity_resid(verb, "PCG", k+1, invif.ops_count, sqrt(rsq)/bsqrt); 

    if (sqrt(rsq) < eps*bsqrt) {
      //        printf("Final rsq = %g\n", rsqNew);
      break;
    }
    
    // 7. z = M^(-1) r
    zero<double>(z, size);
    (*precond_matrix_vector)(z, r, size, precond_info, &verb_prec); 
    
    // 8. beta = r dot z (new) / r dot z
    zdotr_new = dot<double>(r, z, size);
    beta = zdotr_new / zdotr;
    
    zdotr = zdotr_new; 
    
    for (i = 0; i < size; i++) {
      p[i] = z[i] + beta * p[i];
    }
    
    // Compute the new Ap.
    (*matrix_vector)(Ap, p, extra_info); invif.ops_count++;
  } 
    
  if(k == max_iter) {
    //printf("CG: Failed to converge iter = %d, rsq = %e\n", k,rsq);
    invif.success = false;
    //return 0;// Failed convergence 
  }
  else
  {
     invif.success = true;
    k++;
     //printf("CG: Converged in %d iterations.\n", k);
  }
  
  zero<double>(Ap, size); 
  (*matrix_vector)(Ap,phi,extra_info); invif.ops_count++;
  for(i=0; i < size; i++) truersq += (Ap[i] - phi0[i])*(Ap[i] - phi0[i]);
  
  // Free all the things!
  delete[] r;
  delete[] p;
  delete[] Ap;
  delete[] z; 

  
  print_verbosity_summary(verb, "PCG", invif.success, k, invif.ops_count, sqrt(truersq)/bsqrt);
  
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "Preconditioned CG";
  return invif; // Convergence 

} 


inversion_info minv_vector_cg_precond(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, void (*precond_matrix_vector)(complex<double>*,complex<double>*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verb)
{
/// CG solutions to Mphi = b 
//  see http://en.wikipedia.org/wiki/Conjugate_gradient_method

  // Initialize vectors.
  complex<double> *r, *p, *Ap, *z;
  complex<double> alpha, beta, zdotr, zdotr_new;
  double rsq, bsqrt, truersq;
  int k,i;
  inversion_info invif;

  // For preconditioning verbosity.
  inversion_verbose_struct verb_prec;
  shuffle_verbosity_precond(&verb_prec, verb);
  
  // Allocate memory.
  r = new complex<double>[size];
  p = new complex<double>[size];
  Ap = new complex<double>[size];
  z = new complex<double>[size];

  // Initialize values.
  zdotr = zdotr_new = alpha = beta = truersq = 0.0;
  rsq = bsqrt = 0.0;

  // Zero vectors;
  zero<double>(r, size); zero<double>(z, size);
  zero<double>(p, size); zero<double>(Ap, size);
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  // 1. Compute r = b - Ax
  (*matrix_vector)(p, phi, extra_info); invif.ops_count++;
  for (i = 0; i < size; i++)
  {
    r[i] = phi0[i] - p[i];
  }
  
  // 2. z = M^(-1) r
  (*precond_matrix_vector)(z, r, size, precond_info, &verb_prec); 
  
  // 3. p_0 = z_0.
  copy<double>(p, z, size);
  
  // Compute Ap.
  zero<double>(Ap, size);
  (*matrix_vector)(Ap, p, extra_info); invif.ops_count++;
  
  // Compute zdotr.
  zdotr = dot<double>(z, r, size);
  
  printf("Starting loop!\n"); fflush(stdout); 

  // iterate till convergence
  for(k = 0; k< max_iter; k++) {
    
    // 4. alpha = z dot r / p A p 
    alpha = zdotr/dot<double>(p, Ap, size);

    // 5. phi = phi + alpha p
    // 6. r = r - alpha A p
    for (i = 0; i < size; i++)
    {
      phi[i] = phi[i] + alpha*p[i];
      r[i] = r[i] - alpha*Ap[i];
    }
    
    // Exit if new residual is small enough
    rsq = norm2sq<double>(r, size);
    
    print_verbosity_resid(verb, "PCG", k+1, invif.ops_count, sqrt(rsq)/bsqrt); 

    if (sqrt(rsq) < eps*bsqrt) {
      //        printf("Final rsq = %g\n", rsqNew);
      break;
    }
    
    // 7. z = M^(-1) r
    zero<double>(z, size);
    (*precond_matrix_vector)(z, r, size, precond_info, &verb_prec); 
    
    // 8. beta = r dot z (new) / r dot z
    zdotr_new = dot<double>(r, z, size);
    beta = zdotr_new / zdotr;
    
    zdotr = zdotr_new; 
    
    for (i = 0; i < size; i++) {
      p[i] = z[i] + beta * p[i];
    }
    
    // Compute the new Ap.
    (*matrix_vector)(Ap, p, extra_info); invif.ops_count++;
  } 
    
  if(k == max_iter) {
    //printf("CG: Failed to converge iter = %d, rsq = %e\n", k,rsq);
    invif.success = false;
    //return 0;// Failed convergence 
  }
  else
  {
     invif.success = true;
     k++; // Fix a counting issue...
  }
  
  zero<double>(Ap, size); 
  (*matrix_vector)(Ap,phi,extra_info); invif.ops_count++;
  for(i=0; i < size; i++) truersq += real(conj(Ap[i] - phi0[i])*(Ap[i] - phi0[i]));
  
  // Free all the things!
  delete[] r;
  delete[] p;
  delete[] Ap;
  delete[] z; 

  
  print_verbosity_summary(verb, "PCG", invif.success, k, invif.ops_count, sqrt(truersq)/bsqrt);
  
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "Preconditioned CG";
  return invif; // Convergence 
}

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

#include "generic_inverters.h"
#include "generic_vector.h"

using namespace std;

// Solves lhs = A^(-1) rhs
inversion_info minv_vector_cg(double  *phi, double  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(double*,double*,void*), void* extra_info, inversion_verbose_struct* verb)
{
// CG solutions to Mphi = b 
//  see http://en.wikipedia.org/wiki/Conjugate_gradient_method
  
  // Initialize vectors.
  double *r, *p, *Ap;
  double alpha, beta, rsq, rsqNew, bsqrt, truersq;
  int k,i;
  inversion_info invif;

  // Allocate memory.
  r = new double[size];
  p = new double[size];
  Ap = new double[size];

  // Initialize values.
  rsq = 0.0; rsqNew = 0.0; bsqrt = 0.0; truersq = 0.0; k=0;

  // Zero vectors;
  zero<double>(r, size); 
  zero<double>(p, size); zero<double>(Ap, size);
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  // 1. Compute r = b - Ax
  (*matrix_vector)(p, phi, extra_info);
  for (i = 0; i < size; i++)
  {
    r[i] = phi0[i] - p[i];
  }
  
  // 2. p_0 = r_0.
  copy<double>(p, r, size);
  
  // Compute Ap.
  zero<double>(Ap, size);
  (*matrix_vector)(Ap, p, extra_info);
  
  // Compute rsq.
  rsq = norm2sq<double>(r, size);

  // iterate till convergence
  for(k = 0; k< max_iter; k++) {
    
    alpha = rsq/dot<double>(p, Ap, size);

    for (i = 0; i < size; i++)
    {
      phi[i] = phi[i] + alpha*p[i];
      r[i] = r[i] - alpha*Ap[i];
    }
    
    // Exit if new residual is small enough
    rsqNew = norm2sq<double>(r, size);
    
    print_verbosity_resid(verb, "CG", k+1, sqrt(rsqNew)/bsqrt); 

    if (sqrt(rsqNew) < eps*bsqrt) {
      //        printf("Final rsq = %g\n", rsqNew);
      break;
    }
  
    // Update vec using new residual
    beta = rsqNew / rsq;
    rsq = rsqNew; 
    
    for (i = 0; i < size; i++) {
      p[i] = r[i] + beta * p[i];
    }
    
    // Compute the new Ap.
    (*matrix_vector)(Ap, p, extra_info);
  } 
    
  if(k == max_iter) {
    invif.success = false;
    //return 0;// Failed convergence 
  }
  else
  {
     invif.success = true;
    k++; 
  }
  
  (*matrix_vector)(Ap,phi,extra_info);
  for(i=0; i < size; i++) truersq += (Ap[i] - phi0[i])*(Ap[i] - phi0[i]);
  
  // Free all the things!
  delete[] r;
  delete[] p;
  delete[] Ap;

  print_verbosity_summary(verb, "CG", invif.success, k, sqrt(truersq)/bsqrt);
  
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "CG";
  return invif; // Convergence 
} 


inversion_info minv_vector_cg(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verb)
{
// CG solutions to Mphi = b 
//  see http://en.wikipedia.org/wiki/Conjugate_gradient_method

  
  // Initialize vectors.
  complex<double> *r, *p, *Ap;
  complex<double> alpha, beta, denom;
  double rsq, rsqNew, bsqrt, truersq;
  int k,i;
  inversion_info invif;

  // Allocate memory.
  r = new complex<double>[size];
  p = new complex<double>[size];
  Ap = new complex<double>[size];

  // Initialize values.
  rsq = 0.0; rsqNew = 0.0; bsqrt = 0.0; truersq = 0.0; k=0;

  // Zero vectors;
  zero<double>(r, size); 
  zero<double>(p, size); zero<double>(Ap, size);
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  // 1. Compute r = b - Ax
  (*matrix_vector)(p, phi, extra_info);
  for (i = 0; i < size; i++)
  {
    r[i] = phi0[i] - p[i];
  }
  
  // 2. p_0 = r_0.
  copy<double>(p, r, size);
  
  // Compute Ap.
  zero<double>(Ap, size);
  (*matrix_vector)(Ap, p, extra_info);
  
  // Compute rsq.
  rsq = norm2sq<double>(r, size);
  
  print_verbosity_resid(verb, "CG", k+1, sqrt(rsqNew)/bsqrt); 

  // iterate till convergence
  for(k = 0; k< max_iter; k++) {
    
    alpha = rsq/dot<double>(p, Ap, size);

    for (i = 0; i < size; i++)
    {
      phi[i] = phi[i] + alpha*p[i];
      r[i] = r[i] - alpha*Ap[i];
    }
    
    // Exit if new residual is small enough
    rsqNew = norm2sq<double>(r, size);

    if (sqrt(rsqNew) < eps*bsqrt) {
      //        printf("Final rsq = %g\n", rsqNew);
      break;
    }
  
    // Update vec using new residual
    beta = rsqNew / rsq;
    rsq = rsqNew; 
    
    for (i = 0; i < size; i++) {
      p[i] = r[i] + beta * p[i];
    }
    
    // Compute the new Ap.
    (*matrix_vector)(Ap, p, extra_info);
  } 
    
  if(k == max_iter) {
    invif.success = false;
  }
  else
  {
     invif.success = true;
    k++;
  }
  
  (*matrix_vector)(Ap,phi,extra_info);
  for(i=0; i < size; i++) truersq += real(conj(Ap[i] - phi0[i])*(Ap[i] - phi0[i]));
  
  // Free all the things!
  delete[] r;
  delete[] p;
  delete[] Ap;

  
  print_verbosity_summary(verb, "CG", invif.success, k, sqrt(truersq)/bsqrt);
  
  
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "CG";
  return invif; // Convergence 
} 

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

// Solves lhs = A^(-1) rhs using bicgstab
inversion_info minv_vector_bicgstab(double  *phi, double  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(double*,double*,void*), void* extra_info, inversion_verbose_struct* verb)
{
// BICGSTAB solutions to Mphi = b 
//  see https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method

  // Initialize vectors.
  double *r, *r0, *v, *p, *s, *t; 
  double rho, rhoNew, alpha, beta, omega;
  double ssq, bsqrt, truersq; 
  int k,i;
  inversion_info invif;

  // Allocate memory.
  r = new double[size];
  r0 = new double[size];
  v = new double[size];
  p = new double[size];
  s = new double[size];
  t = new double[size];
  
  // Zero vectors. 
  zero<double>(r, size);
  zero<double>(r0, size);
  zero<double>(v, size);
  zero<double>(p, size);
  zero<double>(s, size);
  zero<double>(t, size);

  // Initialize values.
  ssq = 0.0; bsqrt = 0.0; truersq = 0.0;
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq<double>(phi0, size));

  // 1. r = b - Ax. 
  // Take advantage of initial guess in phi.
  (*matrix_vector)(v, phi, extra_info);
  for(i = 0; i<size; i++) {
    r[i] = phi0[i] - v[i]; // 1. r0 = b-Ax0
  }
  
  // 2. rhat0 = r.
  copy<double>(r0, r, size);
  
  // 3. rho = alpha = omega = 1.0.
  rho = alpha = omega = 1.0;
  
  // 4. v = p = 0 (already done).
  
  // iterate till convergence
  for(k = 0; k< max_iter; k++) {
    
    // 5.1. rhoNew = <rhat0, ri-1>
    rhoNew = dot<double>(r0, r, size);

    
    // 5.2. beta = (rhoNew/rho)(alpha/omega_i-1)
    beta = (rhoNew/rho)*(alpha/omega);
    rho = rhoNew;

    // 5.3. p = r + beta(p - omega v)
    for (i = 0; i < size; i++) {
      p[i] = r[i] + beta*p[i] - beta*omega*v[i];
    }
    
    // 5.4. v = Ap
    (*matrix_vector)(v, p, extra_info);
    
    // 5.5. alpha = rho/<rhat0, v>
    alpha = rho/dot<double>(r0, v, size);
    
    // 5.6. s = r - alpha v
    for (i = 0; i < size; i++) {
      s[i] = r[i] - alpha*v[i];
    }
    
    // 5.7. If ||s|| is sufficiently small, x = x+alpha p, quit.
    ssq = norm2sq<double>(s, size);
    
    print_verbosity_resid(verb, "BiCGStab", k+1, sqrt(ssq)/bsqrt);
    
    if (sqrt(ssq) < eps*bsqrt)
    {
      // printf("Final rsq = %g\n", ssq);
      for (i = 0; i < size; i++)
      {
        phi[i] = phi[i] + alpha*p[i];
      }
      break;
    }
    
    // 5.8. t = As
    (*matrix_vector)(t, s, extra_info);
    
    // 4.9. omega = <t, s>/<t, t>;
    omega = dot<double>(t, s, size)/norm2sq<double>(t, size);
    
    // 4.10. x = x + alpha p + omega s
    for (i = 0; i < size; i++) {
      phi[i] = phi[i] + alpha*p[i] + omega*s[i];
    }
    
    // 4.11. If x_i is accurate enough, then quit.
    // We'll ignore this for now.
    
    // 4.12. r = s - omega t;
    for (i = 0; i < size; i++) {
      r[i] = s[i] - omega*t[i];
    }
    
  }
  
  if(k == max_iter) {
    //printf("CG: Failed to converge iter = %d, rsq = %e\n", k,rsq);
    invif.success = false;
  }
  else
  {
     //printf("CG: Converged in %d iterations.\n", k);
    k++; // Fix a counting issue. 
     invif.success = true;
  }
  
  // Check the true residual. 
  truersq = 0.0;
  (*matrix_vector)(v,phi,extra_info);
  for(i=0; i < size; i++) truersq += real(conj(v[i] - phi0[i])*(v[i] - phi0[i]));
  
  // Free all the things!
  delete[] r;
  delete[] r0;
  delete[] v;
  delete[] p;
  delete[] s;
  delete[] t;
  
  print_verbosity_summary(verb, "BiCGStab", invif.success, k, sqrt(truersq)/bsqrt);

  //  printf("# CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "BiCGStab";
  return invif; // Convergence 

} 

inversion_info minv_vector_bicgstab(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verb)
{
// BICGSTAB solutions to Mphi = b 
//  see https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method

  // Initialize vectors.
  complex<double> *r, *r0, *v, *p, *s, *t; 
  complex<double> rho, rhoNew, alpha, beta, omega;
  double ssq, bsqrt, truersq; 
  int k,i;
  inversion_info invif;

  // Allocate memory.
  r = new complex<double>[size];
  r0 = new complex<double>[size];
  v = new complex<double>[size];
  p = new complex<double>[size];
  s = new complex<double>[size];
  t = new complex<double>[size];
  
  // Zero vectors. 
  zero<double>(r, size);
  zero<double>(r0, size);
  zero<double>(v, size);
  zero<double>(p, size);
  zero<double>(s, size);
  zero<double>(t, size);

  // Initialize values.
  ssq = 0.0; bsqrt = 0.0; truersq = 0.0;
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq<double>(phi0, size));

  // 1. r = b - Ax. 
  // Take advantage of initial guess in phi.
  (*matrix_vector)(v, phi, extra_info);
  for(i = 0; i<size; i++) {
    r[i] = phi0[i] - v[i]; // 1. r0 = b-Ax0
  }
  
  // 2. rhat0 = r.
  copy<double>(r0, r, size);
  
  // 3. rho = alpha = omega = 1.0.
  rho = alpha = omega = 1.0;
  
  // 4. v = p = 0 (already done).
  
  // iterate till convergence
  for(k = 0; k< max_iter; k++) {
    
    // 5.1. rhoNew = <rhat0, ri-1>
    rhoNew = dot<double>(r0, r, size);

    
    // 5.2. beta = (rhoNew/rho)(alpha/omega_i-1)
    beta = (rhoNew/rho)*(alpha/omega);
    rho = rhoNew;

    // 5.3. p = r + beta(p - omega v)
    for (i = 0; i < size; i++) {
      p[i] = r[i] + beta*p[i] - beta*omega*v[i];
    }
    
    // 5.4. v = Ap
    (*matrix_vector)(v, p, extra_info);
    
    // 5.5. alpha = rho/<rhat0, v>
    alpha = rho/dot<double>(r0, v, size);
    
    // 5.6. s = r - alpha v
    for (i = 0; i < size; i++) {
      s[i] = r[i] - alpha*v[i];
    }
    
    // 5.7. If ||s|| is sufficiently small, x = x+alpha p, quit.
    ssq = norm2sq<double>(s, size);
    
    print_verbosity_resid(verb, "BiCGStab", k+1, sqrt(ssq)/bsqrt);
    
    if (sqrt(ssq) < eps*bsqrt)
    {
      // printf("Final rsq = %g\n", ssq);
      for (i = 0; i < size; i++)
      {
        phi[i] = phi[i] + alpha*p[i];
      }
      break;
    }
    
    // 5.8. t = As
    (*matrix_vector)(t, s, extra_info);
    
    // 4.9. omega = <t, s>/<t, t>;
    omega = dot<double>(t, s, size)/norm2sq<double>(t, size);
    
    // 4.10. x = x + alpha p + omega s
    for (i = 0; i < size; i++) {
      phi[i] = phi[i] + alpha*p[i] + omega*s[i];
    }
    
    // 4.11. If x_i is accurate enough, then quit.
    // We'll ignore this for now.
    
    // 4.12. r = s - omega t;
    for (i = 0; i < size; i++) {
      r[i] = s[i] - omega*t[i];
    }
    
  }
  
  if(k == max_iter) {
    //printf("CG: Failed to converge iter = %d, rsq = %e\n", k,rsq);
    invif.success = false;
  }
  else
  {
     //printf("CG: Converged in %d iterations.\n", k);
    k++; // Fix a counting issue. 
     invif.success = true;
  }
  
  // Check the true residual. 
  truersq = 0.0;
  (*matrix_vector)(v,phi,extra_info);
  for(i=0; i < size; i++) truersq += real(conj(v[i] - phi0[i])*(v[i] - phi0[i]));
  
  // Free all the things!
  delete[] r;
  delete[] r0;
  delete[] v;
  delete[] p;
  delete[] s;
  delete[] t;
  
  print_verbosity_summary(verb, "BiCGStab", invif.success, k, sqrt(truersq)/bsqrt);

  //  printf("# CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "BiCGStab";
  return invif; // Convergence 
} 

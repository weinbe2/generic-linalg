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
inversion_info minv_vector_bicgstab_precond(double  *phi, double  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(double*,double*,void*), void* extra_info, void (*precond_matrix_vector)(double*,double*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verb)
{
// BICGSTAB solutions to Mphi = b 
  //  see www.mcs.anl.gov/papers/P3039-0912.pdf
  // "Analysis and Practical Use of Flexible BiCGStab

  // Initialize vectors.
  double *r, *r0, *p, *s, *ptilde, *Aptilde, *stilde, *Astilde;
  double rho, rhoNew, alpha, beta, omega;
  double rsq, bsqrt, truersq; 
  int k,i;
  inversion_info invif;
  
  // For preconditioning verbosity.
  inversion_verbose_struct verb_prec;
  shuffle_verbosity_precond(&verb_prec, verb);

  // Allocate memory.
  r = new double[size];
  r0 = new double[size];
  p = new double[size];
  s = new double[size];
  ptilde = new double[size];
  Aptilde = new double[size];
  stilde = new double[size];
  Astilde = new double[size];
  
  // Zero vectors. 
  zero<double>(r, size);
  zero<double>(r0, size);
  zero<double>(p, size);
  zero<double>(s, size);
  zero<double>(ptilde, size);
  zero<double>(Aptilde, size);
  zero<double>(stilde, size);
  zero<double>(Astilde, size);

  // Initialize values.
  rsq = 0.0; bsqrt = 0.0; truersq = 0.0;
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq<double>(phi0, size));

  // 1. r = b - Ax. , r0 = arbitrary (use r).
  // Take advantage of initial guess in phi.
  (*matrix_vector)(Aptilde, phi, extra_info);
  for(i = 0; i<size; i++) {
    r[i] = phi0[i] - Aptilde[i]; // 1. r0 = b-Ax0
  }
  copy<double>(r0, r, size);
  
  // 2. p = r
  copy<double>(p, r, size); 
  
  // 2a. Initialize rho = <r, r0>.
  rho = dot<double>(r0, r, size);
  
  
  // 3. iterate till convergence
  for(k = 0; k< max_iter; k++) {
    
    // 4. ptilde = M^(-1) p.
    zero<double>(ptilde, size);
    (*precond_matrix_vector)(ptilde, p, size, precond_info, &verb_prec);
    
    // 4a. Construct A ptilde.
    zero<double>(Aptilde, size); 
    (*matrix_vector)(Aptilde, ptilde, extra_info);
    
    // 5. alpha = <r0, r>/<r0, Aptilde>
    alpha = rho/dot<double>(r0, Aptilde, size);
    
    // 6. s = r - alpha Aptilde
    for (i = 0; i < size; i++)
    {
      s[i] = r[i] - alpha*Aptilde[i];
    }
    
    // 7. stilde = M^(-1) s
    zero<double>(stilde, size);
    (*precond_matrix_vector)(stilde, s, size, precond_info, &verb_prec);
    
    // 8. Compute Astilde, w = <stilde, Astilde>/(Astilde, Astilde)
    zero<double>(Astilde, size);
    (*matrix_vector)(Astilde, stilde, extra_info);
    //zero<double>(Aptilde, size);
    //(*precond_matrix_vector)(Aptilde, Astilde, size, precond_info, &verb_prec);
    //omega = dot<double>(Astilde, Aptilde, size)/dot<double>(Aptilde, Aptilde, size);
    omega = dot<double>(Astilde, stilde, size)/dot<double>(Astilde, Astilde, size);
    
    // 9. Update phi = phi + alpha*ptilde + omega*stilde
    for (i = 0; i < size; i++)
    {
      phi[i] = phi[i] + alpha*ptilde[i] + omega*stilde[i];
    }
    
    // 10. Update r = s - omega*Astilde
    for (i = 0; i < size; i++)
    {
      r[i] = s[i] - omega*Astilde[i];
    }
    
    // 10a. If ||r|| is sufficiently small, quit.
    rsq = norm2sq<double>(r, size);
    print_verbosity_resid(verb, "Preconditioned BiCGStab", k+1, sqrt(rsq)/bsqrt);
    
    if (sqrt(rsq) < eps*bsqrt)
    {
      break;
    }
    
    // 11. rhoNew = <r0, r>.
    rhoNew = dot<double>(r0, r, size);
    beta = rhoNew/rho*(alpha/omega);
    rho = rhoNew;
    
    // 12. Update p = r + beta*p - omega*beta*Aptilde
    for (i = 0; i < size; i++)
    {
      p[i] = r[i] + beta*(p[i] - omega*Aptilde[i]);
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
  (*matrix_vector)(Aptilde,phi,extra_info);
  truersq = diffnorm2sq<double>(Aptilde, phi0, size);
  
  // Free all the things!
  delete[] r;
  delete[] r0;
  delete[] p;
  delete[] s;
  delete[] ptilde;
  delete[] Aptilde;
  delete[] stilde;
  delete[] Astilde; 
  
  print_verbosity_summary(verb, "Preconditioned BiCGStab", invif.success, k, sqrt(truersq)/bsqrt);

  //  printf("# CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "BiCGStab";
  return invif; // Convergence 

} 

// Performs PBiCGStab(restart_freq) with restarts when restart_freq is hit.
// This may be sloppy, but it works.
inversion_info minv_vector_bicgstab_precond_restart(double  *phi, double  *phi0, int size, int max_iter, double res, int restart_freq, void (*matrix_vector)(double*,double*,void*), void* extra_info, void (*precond_matrix_vector)(double*,double*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verb)
{
  int iter; // counts total number of iterations.
  inversion_info invif;
  double bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  stringstream ss;
  ss << "Preconditioned Restarted BiCGStab(" << restart_freq << ")";

  inversion_verbose_struct verb_rest;
  shuffle_verbosity_restart(&verb_rest, verb);
  
  iter = 0;  
  do
  {
    invif = minv_vector_bicgstab_precond(phi, phi0, size, restart_freq, res, matrix_vector, extra_info, precond_matrix_vector, precond_info, &verb_rest);
    iter += invif.iter;
    
    print_verbosity_restart(verb, ss.str(), iter, sqrt(invif.resSq)/bsqrt);
  }
  while (iter < max_iter && invif.success == false && sqrt(invif.resSq)/bsqrt > res);
  
  invif.iter = iter;
  
  print_verbosity_summary(verb, ss.str(), invif.success, iter, sqrt(invif.resSq)/bsqrt);
  
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


// Solves lhs = A^(-1) rhs using bicgstab
inversion_info minv_vector_bicgstab_precond(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, void (*precond_matrix_vector)(complex<double>*,complex<double>*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verb)
{
// BICGSTAB solutions to Mphi = b 
  //  see www.mcs.anl.gov/papers/P3039-0912.pdf
  // "Analysis and Practical Use of Flexible BiCGStab

  // Initialize vectors.
  complex<double> *r, *r0, *p, *s, *ptilde, *Aptilde, *stilde, *Astilde;
  complex<double> rho, rhoNew, alpha, beta, omega;
  double rsq, bsqrt, truersq; 
  int k,i;
  inversion_info invif;
  
  // For preconditioning verbosity.
  inversion_verbose_struct verb_prec;
  shuffle_verbosity_precond(&verb_prec, verb);

  // Allocate memory.
  r = new complex<double>[size];
  r0 = new complex<double>[size];
  p = new complex<double>[size];
  s = new complex<double>[size];
  ptilde = new complex<double>[size];
  Aptilde = new complex<double>[size];
  stilde = new complex<double>[size];
  Astilde = new complex<double>[size];
  
  // Zero vectors. 
  zero<complex<double>>(r, size);
  zero<complex<double>>(r0, size);
  zero<complex<double>>(p, size);
  zero<complex<double>>(s, size);
  zero<complex<double>>(ptilde, size);
  zero<complex<double>>(Aptilde, size);
  zero<complex<double>>(stilde, size);
  zero<complex<double>>(Astilde, size);

  // Initialize values.
  rsq = 0.0; bsqrt = 0.0; truersq = 0.0;
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq<double>(phi0, size));

  // 1. r = b - Ax. , r0 = arbitrary (use r).
  // Take advantage of initial guess in phi.
  (*matrix_vector)(Aptilde, phi, extra_info);
  for(i = 0; i<size; i++) {
    r[i] = phi0[i] - Aptilde[i]; // 1. r0 = b-Ax0
  }
  copy<double>(r0, r, size);
  
  // 2. p = r
  copy<double>(p, r, size); 
  
  // 2a. Initialize rho = <r, r0>.
  rho = dot<double>(r0, r, size);
  
  
  // 3. iterate till convergence
  for(k = 0; k< max_iter; k++) {
    
    // 4. ptilde = M^(-1) p.
    zero<double>(ptilde, size);
    (*precond_matrix_vector)(ptilde, p, size, precond_info, &verb_prec);
    
    // 4a. Construct A ptilde.
    zero<double>(Aptilde, size); 
    (*matrix_vector)(Aptilde, ptilde, extra_info);
    
    // 5. alpha = <r0, r>/<r0, Aptilde>
    alpha = rho/dot<double>(r0, Aptilde, size);
    
    // 6. s = r - alpha Aptilde
    for (i = 0; i < size; i++)
    {
      s[i] = r[i] - alpha*Aptilde[i];
    }
    
    // 7. stilde = M^(-1) s
    zero<double>(stilde, size);
    (*precond_matrix_vector)(stilde, s, size, precond_info, &verb_prec);
    
    // 8. Compute Astilde, w = <stilde, Astilde>/(Astilde, Astilde)
    zero<double>(Astilde, size);
    (*matrix_vector)(Astilde, stilde, extra_info);
    //zero<double>(Aptilde, size);
    //(*precond_matrix_vector)(Aptilde, Astilde, size, precond_info, &verb_prec);
    //omega = dot<double>(Astilde, Aptilde, size)/dot<double>(Aptilde, Aptilde, size);
    omega = dot<double>(Astilde, stilde, size)/dot<double>(Astilde, Astilde, size);
    
    // 9. Update phi = phi + alpha*ptilde + omega*stilde
    for (i = 0; i < size; i++)
    {
      phi[i] = phi[i] + alpha*ptilde[i] + omega*stilde[i];
    }
    
    // 10. Update r = s - omega*Astilde
    for (i = 0; i < size; i++)
    {
      r[i] = s[i] - omega*Astilde[i];
    }
    
    // 10a. If ||r|| is sufficiently small, quit.
    rsq = norm2sq<double>(r, size);
    print_verbosity_resid(verb, "Preconditioned BiCGStab", k+1, sqrt(rsq)/bsqrt);
    
    if (sqrt(rsq) < eps*bsqrt)
    {
      break;
    }
    
    // 11. rhoNew = <r0, r>.
    rhoNew = dot<double>(r0, r, size);
    beta = rhoNew/rho*(alpha/omega);
    rho = rhoNew;
    
    // 12. Update p = r + beta*p - omega*beta*Aptilde
    for (i = 0; i < size; i++)
    {
      p[i] = r[i] + beta*(p[i] - omega*Aptilde[i]);
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
  (*matrix_vector)(Aptilde,phi,extra_info);
  truersq = diffnorm2sq<double>(Aptilde, phi0, size);
  
  // Free all the things!
  delete[] r;
  delete[] r0;
  delete[] p;
  delete[] s;
  delete[] ptilde;
  delete[] Aptilde;
  delete[] stilde;
  delete[] Astilde; 
  
  print_verbosity_summary(verb, "Preconditioned BiCGStab", invif.success, k, sqrt(truersq)/bsqrt);

  //  printf("# CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "BiCGStab";
  return invif; // Convergence 

} 

// Performs PBiCGStab(restart_freq) with restarts when restart_freq is hit.
// This may be sloppy, but it works.
inversion_info minv_vector_bicgstab_precond_restart(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double res, int restart_freq, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, void (*precond_matrix_vector)(complex<double>*,complex<double>*,int,void*,inversion_verbose_struct*), void* precond_info, inversion_verbose_struct* verb)
{
  int iter; // counts total number of iterations.
  inversion_info invif;
  double bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  stringstream ss;
  ss << "Preconditioned Restarted BiCGStab(" << restart_freq << ")";

  inversion_verbose_struct verb_rest;
  shuffle_verbosity_restart(&verb_rest, verb);
  
  iter = 0;  
  do
  {
    invif = minv_vector_bicgstab_precond(phi, phi0, size, restart_freq, res, matrix_vector, extra_info, precond_matrix_vector, precond_info, &verb_rest);
    iter += invif.iter;
    
    print_verbosity_restart(verb, ss.str(), iter, sqrt(invif.resSq)/bsqrt);
  }
  while (iter < max_iter && invif.success == false && sqrt(invif.resSq)/bsqrt > res);
  
  invif.iter = iter;
  
  print_verbosity_summary(verb, ss.str(), invif.success, iter, sqrt(invif.resSq)/bsqrt);
  
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

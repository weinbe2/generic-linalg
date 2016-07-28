// Fri Apr 29 11:54:22 EDT 2016
// Evan S Weinberg
// C++ file for CR inverter.

// To do:
// 1. Template various functions to support float, double,
//    as well as complex< > variants.

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <complex>

#include "generic_inverters.h"
#include "generic_vector.h"

using namespace std; 

    

// Solves lhs = A^(-1) rhs using CR.
// Taken from section 6.8 of Saad, 2nd Edition.
inversion_info minv_vector_cr(double  *phi, double  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(double*,double*,void*), void* extra_info, inversion_verbose_struct* verb)
{

  // Initialize vectors.
  double *x, *r, *Ar, *p, *Ap;
  double alpha, beta, rsq, bsqrt, truersq;
  int k,i;
  inversion_info invif;

  // Allocate memory.
  x = new double[size];
  r = new double[size];
  Ar = new double[size];
  p = new double[size];
  Ap = new double[size];
  
  // Zero vectors. 
  zero<double>(p, size);  zero<double>(r, size);
  zero<double>(Ap, size); zero<double>(Ar, size);

  // Initialize values.
  rsq = 0.0; bsqrt = 0.0; truersq = 0.0;
  
  // Copy initial guess into solution.
  copy<double>(x, phi, size);
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  // 1. r_0 = b - Ax_0. x is phi, the initial guess.
  (*matrix_vector)(p, x, extra_info); // Put Ax_0 into p, temp.
  for (i = 0; i < size; i++)
  {
    r[i] = phi0[i] - p[i]; // r_0 = b - Ax_0
  }
  
  // 2. p_0 = r_0.
  copy<double>(p, r, size);
  
  // 3. Compute A p_0 = A r_0, presave beta.
  (*matrix_vector)(Ap, p, extra_info); 
  
  copy<double>(Ar, Ap, size);
  beta = dot<double>(Ar, r, size);
  

  // iterate until convergence
  for(k = 0; k< max_iter; k++) {
    
    // 4. alpha = <Ap_k, r>/<Ap_k, Ap_k>
    alpha = dot<double>(Ar, r, size)/norm2sq<double>(Ap, size);

    // 5. x = x + alpha p_k
    // 6. r = r - alpha Ap_k
    for (i = 0; i < size; i++)
    {
      x[i] = x[i] + alpha*p[i];
      r[i] = r[i] - alpha*Ap[i];
    }
    
    // Compute norm.
    rsq = norm2sq<double>(r, size);
    
    print_verbosity_resid(verb, "CR", k+1, sqrt(rsq)/bsqrt); 
    
    // Check convergence. 
    if (sqrt(rsq) < eps*bsqrt) {
      //        printf("Final rsq = %g\n", rsqNew);
      break;
    }
    
    // 7. Compute Ar.
    zero<double>(Ar, size);
    (*matrix_vector)(Ar, r, extra_info);
    
    // 8. b_j = <Ar_{j+1}, Ar_{j+1}>/<Ar_j, r_j> (Update beta)
    beta = dot<double>(Ar, r, size)/beta; // Might be unstable, there's a way to correct this.
    
    // 9. p_{j+1} = r_{j+1} + b_j p_j
    // 10. Ap_{j+1} = Ar_{j+1} + b_j Ap_j
    for (i = 0; i < size; i++)
    {
      p[i] = r[i] + beta*p[i];
      Ap[i] = Ar[i] + beta*Ap[i]; 
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
    k++; // Fix a counting issue...
     //printf("CG: Converged in %d iterations.\n", k);
  }
  
  // Check true residual.
  zero<double>(p,size);
  (*matrix_vector)(p,x,extra_info);
  for(i=0; i < size; i++) truersq += (p[i] - phi0[i])*(p[i] - phi0[i]);
  
  // Copy solution into phi.
  copy<double>(phi, x, size);
  
  // Free all the things!
  delete[] x;
  delete[] r;
  delete[] Ar;
  delete[] p;
  delete[] Ap;

  print_verbosity_summary(verb, "CR", invif.success, k, sqrt(truersq)/bsqrt);
  
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "CR";
  return invif; // Convergence 
} 

// Performs CR(restart_freq) with restarts when restart_freq is hit.
// This may be sloppy, but it works.
inversion_info minv_vector_cr_restart(double  *phi, double  *phi0, int size, int max_iter, double res, int restart_freq, void (*matrix_vector)(double*,double*,void*), void* extra_info, inversion_verbose_struct* verb)
{
  int iter; // counts total number of iterations.
  inversion_info invif;
  double bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  inversion_verbose_struct verb_rest;
  shuffle_verbosity_restart(&verb_rest, verb);
  
  stringstream ss;
  ss << "CR(" << restart_freq << ")";
  
  iter = 0;  
  do
  {
    invif = minv_vector_cr(phi, phi0, size, restart_freq, res, matrix_vector, extra_info, &verb_rest);
    iter += invif.iter;
    
    print_verbosity_restart(verb, ss.str(), iter, sqrt(invif.resSq)/bsqrt);
  }
  while (iter < max_iter && invif.success == false && sqrt(invif.resSq)/bsqrt > res);
  
  invif.iter = iter;
  
  print_verbosity_summary(verb, ss.str(), invif.success, iter, sqrt(invif.resSq)/bsqrt);
  
  invif.name = ss.str();
  // invif.resSq is good.
  if (sqrt(invif.resSq) > res)
  {
    invif.success = false;
  }
  else
  {
    invif.success = true;
  }
  
  return invif;
}


inversion_info minv_vector_cr(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verb)
{

  // Initialize vectors.
  complex<double> *x, *r, *Ar, *p, *Ap;
  double rsq, bsqrt, truersq;
  complex<double> alpha, beta;
  int k,i;
  inversion_info invif;

  // Allocate memory.
  x = new complex<double>[size];
  r = new complex<double>[size];
  Ar = new complex<double>[size];
  p = new complex<double>[size];
  Ap = new complex<double>[size];
  
  // Zero vectors. 
  zero<double>(p, size);  zero<double>(r, size);
  zero<double>(Ap, size); zero<double>(Ar, size);

  // Initialize values.
  rsq = 0.0; bsqrt = 0.0; truersq = 0.0;
  
  // Copy initial guess into solution.
  copy<double>(x, phi, size);
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  // 1. r_0 = b - Ax_0. x is phi, the initial guess.
  (*matrix_vector)(p, x, extra_info); // Put Ax_0 into p, temp.
  for (i = 0; i < size; i++)
  {
    r[i] = phi0[i] - p[i]; // r_0 = b - Ax_0
  }
  
  // 2. p_0 = r_0.
  copy<double>(p, r, size);
  
  // 3. Compute A p_0 = A r_0, presave beta.
  (*matrix_vector)(Ap, p, extra_info); 
  
  copy<double>(Ar, Ap, size);
  beta = dot<double>(Ar, r, size);
  

  // iterate until convergence
  for(k = 0; k< max_iter; k++) {
    
    // 4. alpha = <Ap_k, r>/<Ap_k, Ap_k>
    alpha = dot<double>(Ar, r, size)/norm2sq<double>(Ap, size);

    // 5. x = x + alpha p_k
    // 6. r = r - alpha Ap_k
    for (i = 0; i < size; i++)
    {
      x[i] = x[i] + alpha*p[i];
      r[i] = r[i] - alpha*Ap[i];
    }
    
    // Compute norm.
    rsq = norm2sq<double>(r, size);
    
    print_verbosity_resid(verb, "CR", k+1, sqrt(rsq)/bsqrt); 
    
    // Check convergence. 
    if (sqrt(rsq) < eps*bsqrt) {
      //        printf("Final rsq = %g\n", rsqNew);
      break;
    }
    
    // 7. Compute Ar.
    zero<double>(Ar, size);
    (*matrix_vector)(Ar, r, extra_info);
    
    // 8. b_j = <Ar_{j+1}, Ar_{j+1}>/<Ar_j, r_j> (Update beta)
    beta = dot<double>(Ar, r, size)/beta; // Might be unstable, there's a way to correct this.
    
    // 9. p_{j+1} = r_{j+1} + b_j p_j
    // 10. Ap_{j+1} = Ar_{j+1} + b_j Ap_j
    for (i = 0; i < size; i++)
    {
      p[i] = r[i] + beta*p[i];
      Ap[i] = Ar[i] + beta*Ap[i]; 
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
    k++; // Fix a counting issue...
     //printf("CG: Converged in %d iterations.\n", k);
  }
  
  // Check true residual.
  zero<double>(p,size);
  (*matrix_vector)(p,x,extra_info);
  for(i=0; i < size; i++) truersq += real(conj(p[i] - phi0[i])*(p[i] - phi0[i]));
  
  // Copy solution into phi.
  copy<double>(phi, x, size);
  
  // Free all the things!
  delete[] x;
  delete[] r;
  delete[] Ar;
  delete[] p;
  delete[] Ap;

  print_verbosity_summary(verb, "CR", invif.success, k, sqrt(truersq)/bsqrt);
  
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "CR";
  return invif; // Convergence 
  
}

// Performs CR(restart_freq) with restarts when restart_freq is hit.
// This may be sloppy, but it works.
inversion_info minv_vector_cr_restart(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double res, int restart_freq, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verb)
{
  int iter; // counts total number of iterations.
  inversion_info invif;
  double bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  inversion_verbose_struct verb_rest;
  shuffle_verbosity_restart(&verb_rest, verb);
  
  stringstream ss;
  ss << "CR(" << restart_freq << ")";
  
  iter = 0;  
  do
  {
    invif = minv_vector_cr(phi, phi0, size, restart_freq, res, matrix_vector, extra_info, &verb_rest);
    iter += invif.iter;
    
    print_verbosity_restart(verb, ss.str(), iter, sqrt(invif.resSq)/bsqrt);
  }
  while (iter < max_iter && invif.success == false && sqrt(invif.resSq)/bsqrt > res);
  
  invif.iter = iter;
  
  print_verbosity_summary(verb, ss.str(), invif.success, iter, sqrt(invif.resSq)/bsqrt);
  
  invif.name = ss.str();
  // invif.resSq is good.
  if (sqrt(invif.resSq) > res)
  {
    invif.success = false;
  }
  else
  {
    invif.success = true;
  }
  
  return invif;
}

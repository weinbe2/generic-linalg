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

#include "generic_traits.h"
#include "generic_vector.h"

#include "generic_cg.h"

using namespace std;

/*
// Solves phi = A^(-1) phi0.
template <typename T, typename Treal = typename RealType<T>::Type >
inversion_info minv_vector_cg(T  *phi, T *phi0, int size, int max_iter, Treal eps, void (*matrix_vector)(T*,T*,void*), void* extra_info, inversion_verbose_struct* verb)
{
  // CG solutions to Mphi = b 
//  see http://en.wikipedia.org/wiki/Conjugate_gradient_method
  //typedef typename Tr::Type Treal; // Type with same precision as type passed in, i.e., if T is complex<double>, Treal is a double, if T is double, Treal is also double.
  
  // Initialize vectors.
  T *r, *p, *Ap; // complex or real
  T alpha, beta, denom; // complex or real
  Treal rsq, rsqNew, bsqrt, truersq; 
  int k,i;
  inversion_info invif;

  // Allocate memory.
  r = new T[size];
  p = new T[size];
  Ap = new T[size];

  // Initialize values.
  rsq = 0.0; rsqNew = 0.0; bsqrt = 0.0; truersq = 0.0; k=0;

  // Zero vectors;
  zero<Treal>(r, size); 
  zero<Treal>(p, size); zero<Treal>(Ap, size);
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq<Treal>(phi0, size));
  
  // 1. Compute r = b - Ax
  (*matrix_vector)(p, phi, extra_info);
  for (i = 0; i < size; i++)
  {
    r[i] = phi0[i] - p[i];
  }
  
  // 2. p_0 = r_0.
  copy<Treal>(p, r, size);
  
  // Compute Ap.
  zero<Treal>(Ap, size);
  (*matrix_vector)(Ap, p, extra_info);
  
  // Compute rsq.
  rsq = norm2sq<Treal>(r, size);
  
  print_verbosity_resid(verb, "CG", k+1, sqrt(rsqNew)/bsqrt); 

  // iterate till convergence
  for(k = 0; k< max_iter; k++) {
    
    alpha = rsq/dot<Treal>(p, Ap, size);

    for (i = 0; i < size; i++)
    {
      phi[i] = phi[i] + alpha*p[i];
      r[i] = r[i] - alpha*Ap[i];
    }
    
    // Exit if new residual is small enough
    rsqNew = norm2sq<Treal>(r, size);
      
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
  }
  else
  {
     invif.success = true;
    k++;
  }
  
  (*matrix_vector)(Ap,phi,extra_info);
  truersq = diffnorm2sq<Treal>(Ap, phi0, size);
  
  // Free all the things!
  delete[] r;
  delete[] p;
  delete[] Ap;

  
  print_verbosity_summary(verb, "CG", invif.success, k, sqrt(truersq)/bsqrt);
  
  
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "CG";
  return invif; // Convergence 
}*/

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
  (*matrix_vector)(p, phi, extra_info); invif.ops_count++;
  for (i = 0; i < size; i++)
  {
    r[i] = phi0[i] - p[i];
  }
  
  // 2. p_0 = r_0.
  copy<double>(p, r, size);
  
  // Compute Ap.
  zero<double>(Ap, size);
  (*matrix_vector)(Ap, p, extra_info); invif.ops_count++;
  
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
    
    print_verbosity_resid(verb, "CG", k+1, invif.ops_count, sqrt(rsqNew)/bsqrt); 

    if (sqrt(rsqNew) < eps*bsqrt || k == max_iter-1) {
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
    (*matrix_vector)(Ap, p, extra_info); invif.ops_count++;
  } 
    
  if(k == max_iter-1) {
    invif.success = false;
  }
  else
  {
     invif.success = true;
  }
  k++;
  
  (*matrix_vector)(Ap,phi,extra_info); invif.ops_count++;
  truersq = diffnorm2sq<double>(Ap, phi0, size);
  
  // Free all the things!
  delete[] r;
  delete[] p;
  delete[] Ap;

  print_verbosity_summary(verb, "CG", invif.success, k, invif.ops_count, sqrt(truersq)/bsqrt);
  
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "CG";
  return invif; // Convergence 
} 

// Performs CG(restart_freq) with restarts when restart_freq is hit.
// This may be sloppy, but it works.
inversion_info minv_vector_cg_restart(double  *phi, double  *phi0, int size, int max_iter, double res, int restart_freq, void (*matrix_vector)(double*,double*,void*), void* extra_info, inversion_verbose_struct* verb)
{
  int iter; // counts total number of iterations.
  int ops_count; 
  inversion_info invif;
  double bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  inversion_verbose_struct verb_rest;
  shuffle_verbosity_restart(&verb_rest, verb);
  
  stringstream ss;
  ss << "CG(" << restart_freq << ")";
  
  iter = 0; ops_count = 0; 
  do
  {
    invif = minv_vector_cg(phi, phi0, size, restart_freq, res, matrix_vector, extra_info, &verb_rest);
    iter += invif.iter;
    ops_count += invif.ops_count; 
    
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
  (*matrix_vector)(p, phi, extra_info); invif.ops_count++;
  for (i = 0; i < size; i++)
  {
    r[i] = phi0[i] - p[i];
  }
  
  // 2. p_0 = r_0.
  copy<double>(p, r, size);
  
  // Compute Ap.
  zero<double>(Ap, size);
  (*matrix_vector)(Ap, p, extra_info); invif.ops_count++;
  
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
      
    print_verbosity_resid(verb, "CG", k+1, invif.ops_count, sqrt(rsqNew)/bsqrt);

    if (sqrt(rsqNew) < eps*bsqrt || k == max_iter - 1) {
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
    (*matrix_vector)(Ap, p, extra_info); invif.ops_count++;
  } 
    
  if(k == max_iter-1) {
    invif.success = false;
  }
  else
  {
     invif.success = true;
  }
	
  k++; 
  
  (*matrix_vector)(Ap,phi,extra_info); invif.ops_count++;
  truersq = diffnorm2sq<double>(Ap, phi0, size);
  
  // Free all the things!
  delete[] r;
  delete[] p;
  delete[] Ap;

  
  print_verbosity_summary(verb, "CG", invif.success, k, invif.ops_count, sqrt(truersq)/bsqrt);
  
  
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "CG";
  return invif; // Convergence 
} 


// Performs CG(restart_freq) with restarts when restart_freq is hit.
// This may be sloppy, but it works.
inversion_info minv_vector_cg_restart(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double res, int restart_freq, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verb)
{
  int iter; // counts total number of iterations.
  int ops_count; 
  inversion_info invif;
  double bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  inversion_verbose_struct verb_rest;
  shuffle_verbosity_restart(&verb_rest, verb);
  
  stringstream ss;
  ss << "CG(" << restart_freq << ")";
  
  iter = 0; ops_count = 0; 
  do
  {
    invif = minv_vector_cg(phi, phi0, size, restart_freq, res, matrix_vector, extra_info, &verb_rest);
    iter += invif.iter;
    ops_count += invif.ops_count; 
    
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


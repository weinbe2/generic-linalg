// Mon May  2 14:34:46 EDT 2016
// Evan S Weinberg
// C++ file for successive overrelaxation. 

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

// Solves lhs = A^(-1) rhs using SOR. omega is the OR parameter. 
// If lambda(p) are the eigenvalues, this only converges if |1 - omega*lambda(p)| < 1 is true
// for all eigenvalues. 
inversion_info minv_vector_sor(double  *phi, double  *phi0, int size, int max_iter, double eps, double omega, void (*matrix_vector)(double*,double*,void*), void* extra_info, inversion_verbose_struct* verb)
{
  // Initialize vectors.
  double *Ax, *x, *xnew, *check; 
  double bsqrt, rsq, truersq;
  int i,k;
  inversion_info invif;
  
  stringstream ss;
  ss << "SOR_" << omega;

  // Allocate memory.
  Ax = new double[size];
  x = new double[size];
  xnew = new double[size];
  check = new double[size];
  
  // Copy initial guess phi into x.
  for (i = 0; i < size; i++) {
    x[i] = phi[i];
  }
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  // Initialize values.
  rsq = 0.0; truersq = 0.0;
  
  // Iterate until convergence: x_{n+1} = x_n + omega(b - Ax_n)
  for (k = 0; k < max_iter; k++) {
    // Apply A to x_n.
    (*matrix_vector)(Ax, x, extra_info);
    
    // Update x_new = x + omega(b - Ax)
    for (i = 0; i < size; i++)
    {
      xnew[i] = x[i] + omega*(phi0[i] - Ax[i]);
      check[i] = Ax[i] - phi0[i]; 
    }
    
    // Compute norm.
    rsq = norm2sq<double>(check, size);
    
    print_verbosity_resid(verb, ss.str(), k+1, sqrt(rsq)/bsqrt); 
    
    // Check convergence. 
    if (sqrt(rsq) < eps*bsqrt) {
      //        printf("Final rsq = %g\n", rsqNew);
      break;
    }
    
    // Copy xnew back into x.
    for (i = 0; i < size; i++)
    {
      x[i] = xnew[i]; 
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
    k++;
     //printf("CG: Converged in %d iterations.\n", k);
  }
  
  // Check true residual. 
  (*matrix_vector)(Ax,x,extra_info);
  for(i=0; i < size; i++) truersq += (Ax[i] - phi0[i])*(Ax[i] - phi0[i]);
  
  // Copy solution into phi.
  for (i = 0; i < size; i++) { phi[i] = x[i]; }
  
  // Free all the things!
  delete[] Ax;
  delete[] xnew;
  delete[] x;
  delete[] check;
  
  print_verbosity_summary(verb, ss.str(), invif.success, k, sqrt(invif.resSq)/bsqrt);
  
  invif.resSq = truersq;
  invif.iter = k;
  stringstream ss2;
  ss2 << "SOR omega=" << omega;
  invif.name = ss2.str();
  return invif; // Convergence 
} 

// Solves lhs = A^(-1) rhs using SOR. omega is the OR parameter. 
// If lambda(p) are the eigenvalues, this only converges if |1 - omega*lambda(p)| < 1 is true
// for all eigenvalues. 
inversion_info minv_vector_sor(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, double omega, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verb)
{
  // Initialize vectors.
  complex<double> *Ax, *x, *xnew, *check; 
  double bsqrt, rsq, truersq;
  int i,k;
  inversion_info invif;

  stringstream ss;
  ss << "SOR_" << omega;
  
  // Allocate memory.
  Ax = new complex<double>[size];
  x = new complex<double>[size];
  xnew = new complex<double>[size];
  check = new complex<double>[size];
  
  // Copy initial guess phi into x.
  for (i = 0; i < size; i++) {
    x[i] = phi[i];
  }
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  // Initialize values.
  rsq = 0.0; truersq = 0.0;
  
  // Iterate until convergence: x_{n+1} = x_n + omega(b - Ax_n)
  for (k = 0; k < max_iter; k++) {
    // Apply A to x_n.
    (*matrix_vector)(Ax, x, extra_info);
    
    // Update x_new = x + omega(b - Ax)
    for (i = 0; i < size; i++)
    {
      xnew[i] = x[i] + omega*(phi0[i] - Ax[i]);
      check[i] = Ax[i] - phi0[i]; 
    }
    
    // Compute norm.
    rsq = norm2sq<double>(check, size);
    
    print_verbosity_resid(verb, ss.str(), k+1, sqrt(rsq)/bsqrt); 
    
    // Check convergence. 
    if (sqrt(rsq) < eps*bsqrt) {
      //        printf("Final rsq = %g\n", rsqNew);
      break;
    }
    
    // Copy xnew back into x.
    for (i = 0; i < size; i++)
    {
      x[i] = xnew[i]; 
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
  
  // Check true residual. 
  (*matrix_vector)(Ax,x,extra_info);
  for(i=0; i < size; i++) truersq += real(conj(Ax[i] - phi0[i])*(Ax[i] - phi0[i]));
  
  // Copy solution into phi.
  for (i = 0; i < size; i++) { phi[i] = x[i]; }
  
  // Free all the things!
  delete[] Ax;
  delete[] xnew;
  delete[] x;
  delete[] check;
  
  print_verbosity_summary(verb, ss.str(), invif.success, k, sqrt(invif.resSq)/bsqrt);
  
  invif.resSq = truersq;
  invif.iter = k;
  stringstream ss2;
  ss2 << "SOR omega=" << omega;
  invif.name = ss2.str();
  return invif; // Convergence 
} 

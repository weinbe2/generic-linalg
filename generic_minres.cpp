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

#include "generic_vector.h"

#include "generic_minres.h"

using namespace std;

// Assumes the symmetric part of the matrix is positive definite.
// Taken from section 5.3.2 of Saad, 2nd Edition.
inversion_info minv_vector_minres(double  *phi, double  *phi0, int size, int max_iter, double eps, double omega, void (*matrix_vector)(double*,double*,void*), void* extra_info, inversion_verbose_struct* verbosity)
{
  // Initialize vectors.
  double *p, *r, *x;
  double alpha, rsq, truersq, bsqrt;
  int i,k;
  inversion_info invif;
  
  stringstream ss;
  ss << "MR_" << omega;
  

  // Allocate memory.
  p = new double[size];
  r = new double[size];
  x = new double[size];
  
  // Zero vectors. 
  zero<double>(p, size); zero<double>(r, size);
  
  // Copy initial guess into solution.
  copy<double>(x, phi, size);
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  // 1. r = b - Ax. x is phi, the initial guess. 
  (*matrix_vector)(p,phi,extra_info); // p = Ax, temporarily.
  for (i = 0; i < size; i++) {
    r[i] = phi0[i] - p[i]; // r = b - Ax
  }
  
  // Initialize values.
  alpha = 0.0; rsq = 0.0; truersq = 0.0;
  
  // Iterate until convergence or max_iter is reached.
  for (k = 0; k < max_iter; k++) {
    
    // 2. p = Ar.
    zero<double>(p, size); 
    (*matrix_vector)(p, r, extra_info); // p = Ar
    
    // 3. alpha = <p,r>/<p,p>.
    alpha = dot<double>(p, r, size)/norm2sq<double>(p, size); 
    
    // 3a. Over/underrelaxation alpha *= omega.
    alpha *= omega; 
    
    // 4. x = x + alpha r
    // 5. r = r - alpha p
    for (i = 0; i < size; i++)
    {
      x[i] = x[i] + alpha*r[i];
      r[i] = r[i] - alpha*p[i];
    }
    
    // Compute norm.
    rsq = norm2sq<double>(r, size);
    
    print_verbosity_resid(verbosity, ss.str(), k+1, sqrt(rsq)/bsqrt); 
    
    // Check convergence. 
    if (sqrt(rsq) < eps*bsqrt) {
      //        printf("Final rsq = %g\n", rsqNew);
      break;
    }
    
  }
   
  if(k == max_iter) {
    invif.success = false;
  }
  else
  {
     invif.success = true;
     k++;
  }
  
  // Check true residual. 
  zero<double>(p, size);
  (*matrix_vector)(p,x,extra_info);
  for(i=0; i < size; i++) truersq += (p[i] - phi0[i])*(p[i] - phi0[i]);
  
  // Copy solution into phi.
  for (i = 0; i < size; i++) { phi[i] = x[i]; }
  
  // Free all the things!
  delete[] p;
  delete[] r;
  delete[] x;
  
  print_verbosity_summary(verbosity, ss.str(), invif.success, k, sqrt(truersq)/bsqrt);
  
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "Minimum Residual (MinRes)";
  return invif; // Convergence 
} 

// Version without overrelaxation parameter. 
inversion_info minv_vector_minres(double  *phi, double  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(double*,double*,void*), void* extra_info, inversion_verbose_struct* verbosity)
{
  return minv_vector_minres(phi, phi0, size, max_iter, eps, 1.0, matrix_vector, extra_info, verbosity);
}


// Assumes the Hermitian part of the matrix is positive definite.
// Taken from section 5.3.2 of Saad, 2nd Edition.
inversion_info minv_vector_minres(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, double omega, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verbosity)
{
  // Initialize vectors.
  complex<double> *p, *r, *x;
  double rsq, truersq, bsqrt;
  complex<double> alpha;
  int i,k;
  inversion_info invif;
  
  stringstream ss;
  ss << "MR_" << omega;

  // Allocate memory.
  p = new complex<double>[size];
  r = new complex<double>[size];
  x = new complex<double>[size];
  
  // Zero vectors. 
  zero<double>(p, size); zero<double>(r, size);
  
  // Copy initial guess into solution.
  copy<double>(x, phi, size);
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  // 1. r = b - Ax. x is phi, the initial guess. 
  (*matrix_vector)(p,phi,extra_info); // p = Ax, temporarily.
  for (i = 0; i < size; i++) {
    r[i] = phi0[i] - p[i]; // r = b - Ax
  }
  
  // Initialize values.
  alpha = 0.0; rsq = 0.0; truersq = 0.0;
  
  // Iterate until convergence or max_iter is reached.
  for (k = 0; k < max_iter; k++) {
    
    // 2. p = Ar.
    zero<double>(p, size); 
    (*matrix_vector)(p, r, extra_info); // p = Ar
    
    // 3. alpha = <p,r>/<p,p>.
    alpha = dot<double>(p, r, size)/norm2sq<double>(p, size); 
    
    // 3a. Over/underrelaxation alpha *= omega.
    alpha *= omega; 
    
    // 4. x = x + alpha r
    // 5. r = r - alpha p
    for (i = 0; i < size; i++)
    {
      x[i] = x[i] + alpha*r[i];
      r[i] = r[i] - alpha*p[i];
    }
    
    // Compute norm.
    rsq = norm2sq<double>(r, size);
    
    print_verbosity_resid(verbosity, ss.str(), k+1, sqrt(rsq)/bsqrt); 
    
    // Check convergence. 
    if (sqrt(rsq) < eps*bsqrt) {
      //        printf("Final rsq = %g\n", rsqNew);
      break;
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
  zero<double>(p, size); 
  (*matrix_vector)(p,x,extra_info);
  for(i=0; i < size; i++) truersq += real(conj(p[i] - phi0[i])*(p[i] - phi0[i]));
  
  // Copy solution into phi.
  for (i = 0; i < size; i++) { phi[i] = x[i]; }
  
  // Free all the things!
  delete[] p;
  delete[] r;
  delete[] x;
  
  print_verbosity_summary(verbosity, ss.str(), invif.success, k, sqrt(truersq)/bsqrt);
  
  //  printf("# CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "Minimum Residual (MinRes)";
  return invif; // Convergence 
} 

// Version without overrelaxation parameter. 
inversion_info minv_vector_minres(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verbosity)
{
  return minv_vector_minres(phi, phi0, size, max_iter, eps, 1.0, matrix_vector, extra_info, verbosity);
}



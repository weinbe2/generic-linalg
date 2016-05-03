// Tue May  3 16:52:46 EDT 2016
// Evan S Weinberg
// C++ file for power iterations. 

// To do:
// 1. Add complex. 
// 2. Make poweriter return largest eigenvector, as well. 
// 2. Template to support float, double. 

#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <complex>

#include "generic_eigenvalues.h"
#include "generic_vector.h"

using namespace std;

// Performs power iterations to find the largest eigenvalue. 
// Put the result into eig, use phi0 as an initial guess. 
eigenvalue_info eig_vector_poweriter(double* eig, double* phi0, int size, int max_iter, double relres, void (*matrix_vector)(double*,double*,void*), void* extra_info)
{
  // Initialize variables.
  double *x, *q;
  double beta, beta_new; 
  int i,k;
  eigenvalue_info eigif; 
  
  // Allocate memory.
  x = new double[size];
  q = new double[size];
  
  // Copy phi0 into x.
  copy(x, phi0, size); 
  
  // Find the norm of x.
  beta_new = sqrt(norm2sq<double>(x, size)); 
  
  // Normalize x.
  for (i = 0; i < size; i++)
  {
    q[i] = x[i]/beta_new;
  }
  
  // iterate until convergence.
  for (k = 0; k < max_iter; k++)
  {
    // Save the old beta.
    beta = beta_new;
    
    // x = A*q
    (*matrix_vector)(x, q, extra_info);
    
    // beta_n is the norm of the new vector.
    beta_new = sqrt(norm2sq<double>(x, size));
    
    // if beta hasn't changed by much, quit.
    if (fabs(beta - beta_new) < relres) { break; }
    
    // Normalize x.
    for (i = 0; i < size; i++)
    {
      q[i] = x[i]/beta_new;
    }

  }
  
  if (k == max_iter)
  {
    eigif.success = false;
  }
  else
  {
    eigif.success = true;
  }
  
  // Clean up!
  delete[] x;
  delete[] q;
  
  // Return
  *eig = beta_new; 
  eigif.relative_diff = fabs(beta - beta_new);
  eigif.iter = k;
  eigif.name = "Power Iteration";
  return eigif; // Convergence 
}

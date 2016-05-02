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
inversion_info minv_vector_sor(double  *phi, double  *phi0, int size, int max_iter, double eps, double omega, void (*matrix_vector)(double*,double*,void*), void* extra_info)
{
  // Initialize vectors.
  double *Ax, *x, *xnew, *check; 
  double bsqrt, rsq, truersq;
  int i,k;
  inversion_info invif;

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
  for(i=0; i < size; i++) truersq += (Ax[i] - phi0[i])*(Ax[i] - phi0[i]);
  
  // Copy solution into phi.
  for (i = 0; i < size; i++) { phi[i] = x[i]; }
  
  // Free all the things!
  delete[] Ax;
  delete[] xnew;
  delete[] x;
  delete[] check;
  
  //  printf("# CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  
  invif.resSq = truersq;
  invif.iter = k;
  stringstream ss;
  ss << "SOR omega=" << omega;
  invif.name = ss.str();
  return invif; // Convergence 
} 

/*
inversion_info minv_vector_cg(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info)
{
// CG solutions to Mphi = b 
//  see http://en.wikipedia.org/wiki/Conjugate_gradient_method

  // Initialize vectors.
  complex<double> *res, *resNew, *pvec, *Apvec, *pvec_tmp;
  complex<double> alpha, beta, denom;
  double rsq, rsqNew, bsqrt, truersq;
  int k,i;
  inversion_info invif;

  // Allocate memory.
  res = new complex<double>[size];
  resNew = new complex<double>[size];
  pvec = new complex<double>[size];
  Apvec = new complex<double>[size];
  pvec_tmp = new complex<double>[size];
  //res = (double*)malloc(size*sizeof(double));
  //resNew = (double*)malloc(size*sizeof(double));
  //pvec = (double*)malloc(size*sizeof(double));
  //Apvec = (double*)malloc(size*sizeof(double));
  //pvec_tmp = (double*)malloc(size*sizeof(double));

  // Initialize values.
  rsq = 0.0; rsqNew = 0.0; bsqrt = 0.0; truersq = 0.0;

  // Take advantage of initial guess in phi.
  (*matrix_vector)(Apvec, phi, extra_info);
  for(i = 0; i<size; i++) {
    res[i] = phi0[i] - Apvec[i];  
    resNew[i] = 0.0;
    pvec[i] = res[i];
    Apvec[i] = 0.0; // We don't need this component anymore.
    //bsqrt += phi0[i]*phi0[i];
    pvec_tmp[i] = 0.0;
  }
  bsqrt = sqrt(norm2sq<double>(phi0, size));

  // iterate till convergence
  for(k = 0; k< max_iter; k++) {
    rsq = norm2sq<double>(res, size);
    //rsq = 0;
    //for (i = 0; i < size; i++) rsq += res[i]*res[i];

    (*matrix_vector)(Apvec, pvec, extra_info);
    denom = dot<double>(pvec, Apvec, size);
    //denom = 0;
    //for(i=0; i< size; i++) denom += pvec[i]*Apvec[i];
    alpha = rsq/denom;

    for(i=0; i < size; i++)  phi[i] +=  alpha * pvec[i];
    for(i=0; i < size; i++) resNew[i] = res[i]- alpha*Apvec[i];
    
    // Exit if new residual is small enough
    rsqNew = norm2sq<double>(resNew, size);
    //rsqNew = 0;
    //for (i = 0; i < size; i++) rsqNew += resNew[i]*resNew[i];
    
    if (sqrt(rsqNew) < eps*bsqrt) {
      //        printf("Final rsq = %g\n", rsqNew);
      break;
    }
  
    // Update vec using new residual
    beta = rsqNew / rsq;
    for (i = 0; i < size; i++) {
      pvec[i] = resNew[i] + beta * pvec[i];
      res[i] = resNew[i];
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
  
  (*matrix_vector)(Apvec,phi,extra_info);
  for(i=0; i < size; i++) truersq += pow(abs(Apvec[i] - phi0[i]),2);
  
  // Free all the things!
  delete[] res;
  delete[] resNew;
  delete[] pvec;
  delete[] Apvec;
  delete[] pvec_tmp;
  //free(res);
  //free(resNew);
  //free(pvec);
  //free(Apvec);
  //free(pvec_tmp);

  
  //  printf("# CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "CG";
  return invif; // Convergence 
} */

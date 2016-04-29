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

// Based on Minv_phi_3d in fem.c
// Solves lhs = A^(-1) rhs
inversion_info minv_vector_cg(double  *phi, double  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(double*,double*,void*), void* extra_info)
{
// CG solutions to Mphi = b 
//  see http://en.wikipedia.org/wiki/Conjugate_gradient_method

  // Initialize vectors.
  double *res, *resNew, *pvec, *Apvec, *pvec_tmp;
  double alpha, beta, denom, rsq, rsqNew, bsqrt, truersq;
  int k,i;
  inversion_info invif;

  // Allocate memory.
  res = new double[size];
  resNew = new double[size];
  pvec = new double[size];
  Apvec = new double[size];
  pvec_tmp = new double[size];
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
  for(i=0; i < size; i++) truersq += (Apvec[i] - phi0[i])*(Apvec[i] - phi0[i]);
  
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
} 


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
} 

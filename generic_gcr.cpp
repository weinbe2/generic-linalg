// Fri Apr 29 11:54:22 EDT 2016
// Evan S Weinberg
// C++ file for CR/GCR inverter.

// To do:
// 1. Template various functions to support float, double,
//    as well as complex< > variants.
// 2. Correct CR/GCR. 

#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <complex>

#include "generic_inverters.h"
#include "generic_vector.h"

using namespace std; 

    

// Solves lhs = A^(-1) rhs
// Currently does not work!
/*inversion_info minv_vector_gcr(double  *phi, double  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(double*,double*,void*), void* extra_info)
{
// GCR solutions to Mphi = b 
//  see https://en.wikipedia.org/wiki/Conjugate_residual_method

  // Initialize vectors.
  double *res, *resNew, *Ares, *pvec, *Apvec, *pvec_tmp;
  double alpha, beta, denom, rAr, rArNew, rsq, bsqrt, truersq, Apsq; 
  int k,i;
  inversion_info invif;

  // Allocate memory.
  res = (double*)malloc(size*sizeof(double));
  resNew = (double*)malloc(size*sizeof(double));
  Ares = (double*)malloc(size*sizeof(double));
  pvec = (double*)malloc(size*sizeof(double));
  Apvec = (double*)malloc(size*sizeof(double));
  pvec_tmp = (double*)malloc(size*sizeof(double));

  // Initialize values.
  rAr = 0.0; rArNew = 0.0; rsq = 0.0; bsqrt = 0.0; truersq = 0.0, Apsq = 0.0;

  // Take advantage of initial guess in phi.
  (*matrix_vector)(Apvec, phi, extra_info);
  for(i = 0; i<size; i++) {
    res[i] = phi0[i] - Apvec[i];  
    resNew[i] = 0.0;
	Ares[i] = 0.0;
    pvec[i] = res[i];
    Apvec[i] = 0.0; // We don't need this component anymore.
    bsqrt += phi0[i]*phi0[i];
    pvec_tmp[i] = 0.0;
  }
  bsqrt = sqrt(bsqrt);
  
  // Compute A r_0, A p_0 ( = A r_0)
  (*matrix_vector)(Ares, res, extra_info);
  for (i = 0; i < size; i++) {
	  Apvec[i] = Ares[i];
	  Apsq = Apsq + (Apvec[i]*Apvec[i]);
	  rAr += Ares[i]*res[i];
  }
  
  //cout << "Iteration 0 rAr " << rAr << endl;
  //cout << "Iteration 0 Apsq " << Apsq << endl;
  
  if (abs(rAr) < 1e-20)
  {
    printf("GCR: Failed to converge iter = 0 because rAr is zero.\n");
    invif.resSq = rAr;
    invif.iter = 0;
    invif.name = "GCR";
    invif.success = false;
    return invif;
  }

  // iterate till convergence.
  // This can be optimised to only require one matrix
  // multiply per iteration.
  for(k = 0; k< max_iter; k++) {
	  
	// alpha = r_k A r_k / (A p_k)^2 = rAr/Apsq
	alpha = rAr/Apsq; 
	
	// x_{k+1} = x_k + alpha_k p_k
	// r_{k+1} = r_k - alpha_k A p_k
	for(i=0;i < size; i++) {
		phi[i] += alpha * pvec[i];
		res[i] -= alpha * Apvec[i];
    }
	
	// Exit if new residual is small enough
    rsq = 0.0;
    for (i = 0; i < size; i++) rsq += res[i]*res[i];
    
	if (sqrt(rsq) < eps*bsqrt) {
      //        printf("Final rsq = %g\n", rsqNew);
      break;
    }

	// Else, update rAr for beta = rAr_{new}/rAr_{old}
	
    (*matrix_vector)(Ares, res, extra_info);
    rArNew = 0.0;
    for(i=0; i< size; i++) rArNew += res[i]*Ares[i];
    beta = rArNew/rAr;
	 rAr = rArNew;
	 
	 if (abs(rAr) < 1e-20)
    {
      printf("GCR: Failed to converge iter = %d because rAr is zero.\n", k+1);
      invif.resSq = rsq;
      invif.iter = k;
      invif.name = "GCR";
      invif.success = false;
      return invif;
    }
	 
	 //cout << "Iteration " << (k+1) << " rAr " << rAr << endl;
	
	// Update p, Ap, Apsq.
	// p_{k+1} = r_{k+1} + beta_k p_k
	// A p_{k+1} = A r_{k+1} + beta_k A p_k
	Apsq = 0.0;
	for (i = 0; i < size; i++) {
	  pvec[i] = res[i] + beta*pvec[i];
	  Apvec[i] = Ares[i] + beta*Apvec[i];
	  Apsq += Apvec[i]*Apvec[i];
	}
	
   //cout << "Iteration " << (k+1) << " Apsq " << Apsq << endl;
	
  } 
    
  if(k == max_iter) {
    invif.success = false;
    //printf("GCR: Failed to converge iter = %d, rsq = %e\n", k,rsq);
    //return 0;// Failed convergence 
  }
  else
  {
    invif.success = true;
     //printf("GCR: Converged in %d iterations.\n", k);
  }
  
  truersq = 0.0;
  (*matrix_vector)(Apvec,phi,extra_info);
  for(i=0; i < size; i++) truersq += (Apvec[i] - phi0[i])*(Apvec[i] - phi0[i]);
  
  // Free all the things!
  free(res);
  free(resNew);
  free(Ares);
  free(pvec);
  free(Apvec);
  free(pvec_tmp);


  //  printf("# CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "GCR";
  return invif; // Convergence 
} */




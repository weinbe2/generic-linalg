// Fri Apr 29 11:54:22 EDT 2016
// Evan S Weinberg
// C++ file for GCR inverter.

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

    

// Solves lhs = A^(-1) rhs using GCR.
// Taken from section 6.9 of Saad, 2nd Edition.
inversion_info minv_vector_gcr(double  *phi, double  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(double*,double*,void*), void* extra_info)
{

  // Initialize vectors.
  double *x, *r, *Ar, *p, *Ap;
  vector<double*> p_store, Ap_store; // GCR requires explicit reorthogonalization against old search vectors.
  double alpha, beta_ij, rsq, bsqrt, truersq;
  int k,i,ii;
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
  
  // 3. Compute A p_0.
  (*matrix_vector)(Ap, p, extra_info); 

  // iterate until convergence
  for(k = 0; k< max_iter; k++) {
    
    // If we've hit here, push the latest p, Ap to storage.
    p_store.push_back(p);
    Ap_store.push_back(Ap);
    
    // 4. alpha = <Ap_k, r>/<Ap_k, Ap_k>
    alpha = dot<double>(Ap, r, size)/norm2sq<double>(Ap, size);
    
    // 5. x = x + alpha p_k
    // 6. r = r - alpha Ap_k
    for (i = 0; i < size; i++)
    {
      x[i] = x[i] + alpha*p[i];
      r[i] = r[i] - alpha*Ap[i];
    }
    
    // Compute norm.
    rsq = norm2sq<double>(r, size);
    
    // Check convergence. 
    if (sqrt(rsq) < eps*bsqrt) {
      //        printf("Final rsq = %g\n", rsqNew);
      break;
    }
    
    // 7. Compute Ar.
    zero<double>(Ar, size);
    (*matrix_vector)(Ar, r, extra_info);
    
    // 8. b_ij = -<Ar_{j+1}, Ap_i>/<Ap_i, Ap_i> for i = 0, ..., j
    // 9. p_{j+1} = r_{j+1} + sum_i=0^j b_ij p_i
    // 10. Ap_{j+1} = Ar_{j+1} + sum_i=0^j b_ij Ap_i
    p = new double[size];  copy<double>(p, r, size);
    Ap = new double[size]; copy<double>(Ap, Ar, size);
    for (ii = 0; ii <= k; ii++)
    {
      beta_ij = -dot<double>(Ap_store[ii], Ar, size)/norm2sq<double>(Ap_store[ii], size);
      for (i = 0; i < size; i++)
      {
        p[i] += beta_ij*p_store[ii][i];
        Ap[i] += beta_ij*Ap_store[ii][i]; 
      }
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
  k = p_store.size();
  for (i = 0; i < k; i++)
  {
    delete[] p_store[i];
    delete[] Ap_store[i];
  }

  
  //  printf("# CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "GCR";
  return invif; // Convergence 
} 

// Performs GCR(restart_freq) with restarts when restart_freq is hit.
// This may be sloppy, but it works.
inversion_info minv_vector_gcr_restart(double  *phi, double  *phi0, int size, int max_iter, double res, int restart_freq, void (*matrix_vector)(double*,double*,void*), void* extra_info)
{
  int iter; // counts total number of iterations.
  inversion_info invif;

  iter = 0;  
  do
  {
    invif = minv_vector_gcr(phi, phi0, size, restart_freq, res, matrix_vector, extra_info);
    iter += invif.iter;
  }
  while (iter < max_iter && invif.success == false && sqrt(invif.resSq) > res);
  
  invif.iter = iter;
  stringstream ss;
  ss << "GCR(" << restart_freq << ")";
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


inversion_info minv_vector_gcr(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info)
{

  // Initialize vectors.
  complex<double> *x, *r, *Ar, *p, *Ap;
  vector<complex<double>*> p_store, Ap_store; // GCR requires explicit reorthogonalization against old search vectors.
  double rsq, bsqrt, truersq;
  complex<double> alpha, beta_ij;
  int k,i,ii;
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
  
  // 3. Compute A p_0.
  (*matrix_vector)(Ap, p, extra_info); 

  // iterate until convergence
  for(k = 0; k< max_iter; k++) {
    
    // If we've hit here, push the latest p, Ap to storage.
    p_store.push_back(p);
    Ap_store.push_back(Ap);
    
    // 4. alpha = <r, Ap_k>/<Ap_k, Ap_k>
    alpha = dot<double>(Ap, r, size)/norm2sq<double>(Ap, size);
    
    // 5. x = x + alpha p_k
    // 6. r = r - alpha Ap_k
    for (i = 0; i < size; i++)
    {
      x[i] = x[i] + alpha*p[i];
      r[i] = r[i] - alpha*Ap[i];
    }
    
    // Compute norm.
    rsq = norm2sq<double>(r, size);
    
    // Check convergence. 
    if (sqrt(rsq) < eps*bsqrt) {
      //        printf("Final rsq = %g\n", rsqNew);
      break;
    }
    
    // 7. Compute Ar.
    zero<double>(Ar, size);
    (*matrix_vector)(Ar, r, extra_info);
    
    // 8. b_ij = -<Ap_i, Ar_{j+1}>/<Ap_i, Ap_i> for i = 0, ..., j
    // 9. p_{j+1} = r_{j+1} + sum_i=0^j b_ij p_i
    // 10. Ap_{j+1} = Ar_{j+1} + sum_i=0^j b_ij Ap_i
    p = new complex<double>[size];  copy<double>(p, r, size);
    Ap = new complex<double>[size]; copy<double>(Ap, Ar, size);
    for (ii = 0; ii <= k; ii++)
    {
      beta_ij = -dot<double>(Ap_store[ii], Ar, size)/norm2sq<double>(Ap_store[ii], size);
      for (i = 0; i < size; i++)
      {
        p[i] += beta_ij*p_store[ii][i];
        Ap[i] += beta_ij*Ap_store[ii][i]; 
      }
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
  k = p_store.size();
  for (i = 0; i < k; i++)
  {
    delete[] p_store[i];
    delete[] Ap_store[i];
  }

  
  //  printf("# CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "GCR";
  return invif; // Convergence 
} 

// Performs GCR with restarts when restart_freq is hit.
// This may be sloppy, but it works.
inversion_info minv_vector_gcr_restart(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double res, int restart_freq, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info)
{
  int iter; // counts total number of iterations.
  inversion_info invif;

  iter = 0;  
  do
  {
    invif = minv_vector_gcr(phi, phi0, size, restart_freq, res, matrix_vector, extra_info);
    iter += invif.iter;
  }
  while (iter < max_iter && invif.success == false && sqrt(invif.resSq) > res);
  
  invif.iter = iter;
  stringstream ss;
  ss << "GCR(" << restart_freq << ")";
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

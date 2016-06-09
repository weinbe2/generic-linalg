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

#include "generic_inverters.h"
#include "generic_vector.h"

using namespace std;

// Solves lhs = A^(-1) rhs using bicgstab
inversion_info minv_vector_bicgstab(double  *phi, double  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(double*,double*,void*), void* extra_info)
{
// BICGSTAB solutions to Mphi = b 
//  see https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method

  // Initialize vectors.
  double *r, *r0, *v, *p, *s, *t; 
  double rho, rhoNew, alpha, beta, omega, tmp; 
  double ssq, bsqrt, truersq, ts, tt; 
  int k,i;
  inversion_info invif;

  // Allocate memory.
  r = new double[size];
  r0 = new double[size];
  v = new double[size];
  p = new double[size];
  s = new double[size];
  t = new double[size];
  //r = (double*)malloc(size*sizeof(double));
  //r0 = (double*)malloc(size*sizeof(double));
  //v = (double*)malloc(size*sizeof(double));
  //p = (double*)malloc(size*sizeof(double));
  //s = (double*)malloc(size*sizeof(double));
  //t = (double*)malloc(size*sizeof(double));

  // Initialize values.
  ssq = 0.0; bsqrt = 0.0; truersq = 0.0;

  // Take advantage of initial guess in phi.
  (*matrix_vector)(v, phi, extra_info);
  for(i = 0; i<size; i++) {
    r[i] = phi0[i] - v[i]; // 1. r0 = b-Ax0
    r0[i] = r[i]; // 2. Assign rhat0 = r0.
    rho = alpha = omega = 1.0; // 3. Assign initial values.
    v[i] = p[i] = 0.0; // 4. v0 = p0 = 0.
    //bsqrt += phi0[i]*phi0[i]; // Used to check if residual is small.
  }
  bsqrt = sqrt(norm2sq<double>(phi0, size)); // Used to check if residual is small.
  //bsqrt = sqrt(bsqrt);

  // iterate till convergence
  for(k = 0; k< max_iter; k++) {
    // 5.1. rhoNew = <rhat0, ri-1>
    rhoNew = dot<double>(r0, r, size);
    //rhoNew = 0.0;
    //for (i = 0; i < size; i++) { rhoNew += r0[i]*r[i]; }
    
    // 5.2. beta = (rhoNew/rho)(alpha/omega_i-1)
    beta = (rhoNew/rho)*(alpha/omega);
    rho = rhoNew;

    // 5.3. p = r + beta(p - omega v)
    for (i = 0; i < size; i++) {
      p[i] = r[i] + beta*p[i] - beta*omega*v[i];
    }
    
    // 5.4. v = Ap
    (*matrix_vector)(v, p, extra_info);
    
    // 5.5. alpha = rho/<rhat0, v>
    tmp = dot<double>(r0, v, size);
    //tmp = 0.0;
    //for (i = 0; i < size; i++) {
    //  tmp += r0[i]*v[i];
    //}
    alpha = rho/tmp;
    
    // 5.6. s = r - alpha v
    for (i = 0; i < size; i++) {
      s[i] = r[i] - alpha*v[i];
    }
    
    // 5.7. If ||s|| is sufficiently small, x = x+alpha p, quit.
    ssq = norm2sq<double>(s, size);
    //ssq = 0.0;
    //for (i = 0; i < size; i++) {
    //  ssq += s[i]*s[i];
    //}
    if (sqrt(ssq) < eps*bsqrt)
    {
      // printf("Final rsq = %g\n", ssq);
      phi[i] = phi[i] + alpha*p[i];
      break;
    }
    
    // 5.8. t = As
    (*matrix_vector)(t, s, extra_info);
    
    // 4.9. omega = <t, s>/<t, t>;
    ts = dot<double>(t, s, size);
    tt = norm2sq<double>(t, size);
    //ts = tt = 0.0;
    //for (i = 0; i < size; i++) {
    //  ts += t[i]*s[i];
    //  tt += t[i]*t[i];
    //}
    omega = ts/tt;
    
    // 4.10. x = x + alpha p + omega s
    for (i = 0; i < size; i++) {
      phi[i] = phi[i] + alpha*p[i] + omega*s[i];
    }
    
    // 4.11. If x_i is accurate enough, then quit.
    // We'll ignore this for now.
    
    // 4.12. r = s - omega t;
    for (i = 0; i < size; i++) {
      r[i] = s[i] - omega*t[i];
    }
    
  } 
    
  if(k == max_iter) {
    //printf("CG: Failed to converge iter = %d, rsq = %e\n", k,rsq);
    invif.success = false;
  }
  else
  {
     //printf("CG: Converged in %d iterations.\n", k);
     invif.success = true;
  }
  
  (*matrix_vector)(v,phi,extra_info);
  for(i=0; i < size; i++) truersq += (v[i] - phi0[i])*(v[i] - phi0[i]);
  
  // Free all the things!
  delete[] r;
  delete[] r0;
  delete[] v;
  delete[] p;
  delete[] s;
  delete[] t;
  //free(r);
  //free(r0);
  //free(v);
  //free(p);
  //free(s);
  //free(t);


  //  printf("# CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "BiCGStab";
  return invif; // Convergence 
} 

inversion_info minv_vector_bicgstab(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info)
{
// BICGSTAB solutions to Mphi = b 
//  see https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method

  // Initialize vectors.
  complex<double> *r, *r0, *v, *p, *s, *t; 
  complex<double> rho, rhoNew, alpha, beta, omega, tmp, ts; 
  double ssq, bsqrt, truersq, tt; 
  int k,i;
  inversion_info invif;

  // Allocate memory.
  r = new complex<double>[size];
  r0 = new complex<double>[size];
  v = new complex<double>[size];
  p = new complex<double>[size];
  s = new complex<double>[size];
  t = new complex<double>[size];

  // Initialize values.
  ssq = 0.0; bsqrt = 0.0; truersq = 0.0;

  // Take advantage of initial guess in phi.
  (*matrix_vector)(v, phi, extra_info);
  for(i = 0; i<size; i++) {
    r[i] = phi0[i] - v[i]; // 1. r0 = b-Ax0
    r0[i] = r[i]; // 2. Assign rhat0 = r0.
    rho = alpha = omega = 1.0; // 3. Assign initial values.
    v[i] = p[i] = 0.0; // 4. v0 = p0 = 0.
    //bsqrt += phi0[i]*phi0[i]; // Used to check if residual is small.
  }
  bsqrt = sqrt(norm2sq<double>(phi0, size)); // Used to check if residual is small.
  //bsqrt = sqrt(bsqrt);
  // iterate till convergence
  for(k = 0; k< max_iter; k++) {
    // 5.1. rhoNew = <rhat0, ri-1>
    rhoNew = dot<double>(r0, r, size);
    //rhoNew = 0.0;
    //for (i = 0; i < size; i++) { rhoNew += r0[i]*r[i]; }
    
    // 5.2. beta = (rhoNew/rho)(alpha/omega_i-1)
    beta = (rhoNew/rho)*(alpha/omega);
    rho = rhoNew;

    // 5.3. p = r + beta(p - omega v)
    for (i = 0; i < size; i++) {
      p[i] = r[i] + beta*p[i] - beta*omega*v[i];
    }
    
    // 5.4. v = Ap
    (*matrix_vector)(v, p, extra_info);
    
    // 5.5. alpha = rho/<rhat0, v>
    tmp = dot<double>(r0, v, size);
    //tmp = 0.0;
    //for (i = 0; i < size; i++) {
    //  tmp += r0[i]*v[i];
    //}
    alpha = rho/tmp;
    
    // 5.6. s = r - alpha v
    for (i = 0; i < size; i++) {
      s[i] = r[i] - alpha*v[i];
    }
    
    // 5.7. If ||s|| is sufficiently small, x = x+alpha p, quit.
    ssq = norm2sq<double>(s, size);
    //ssq = 0.0;
    //for (i = 0; i < size; i++) {
    //  ssq += s[i]*s[i];
    //}
    if (sqrt(ssq) < eps*bsqrt)
    {
      // printf("Final rsq = %g\n", ssq);
      phi[i] = phi[i] + alpha*p[i];
      break;
    }
    
    // 5.8. t = As
    (*matrix_vector)(t, s, extra_info);
    
    // 4.9. omega = <t, s>/<t, t>;
    ts = dot<double>(t, s, size);
    tt = norm2sq<double>(t, size);
    //ts = tt = 0.0;
    //for (i = 0; i < size; i++) {
    //  ts += t[i]*s[i];
    //  tt += t[i]*t[i];
    //}
    omega = ts/tt;
    
    // 4.10. x = x + alpha p + omega s
    for (i = 0; i < size; i++) {
      phi[i] = phi[i] + alpha*p[i] + omega*s[i];
    }
    
    // 4.11. If x_i is accurate enough, then quit.
    // We'll ignore this for now.
    
    // 4.12. r = s - omega t;
    for (i = 0; i < size; i++) {
      r[i] = s[i] - omega*t[i];
    }
    
  }
  
  if(k == max_iter) {
    //printf("CG: Failed to converge iter = %d, rsq = %e\n", k,rsq);
    invif.success = false;
  }
  else
  {
     //printf("CG: Converged in %d iterations.\n", k);
     invif.success = true;
  }
  
  (*matrix_vector)(v,phi,extra_info);
  for(i=0; i < size; i++) truersq += real(conj(v[i] - phi0[i])*(v[i] - phi0[i]));
  
  // Free all the things!
  delete[] r;
  delete[] r0;
  delete[] v;
  delete[] p;
  delete[] s;
  delete[] t;

  //  printf("# CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "BiCGStab";
  return invif; // Convergence 
} 

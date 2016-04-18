// Wed Jan  6 09:49:39 EST 2016
// Evan S Weinberg
// C++ file for generic inverters.

// To do:
// 1. Template various functions to support float, double,
//    as well as complex< > variants.
// 2. Improve GMRES to get residual from reduced matrix.
//    This saves one matrix op per iteration. (We currently
//    lazily use the matrix op again to get the residual.)

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


// Solves lhs = A^(-1) rhs using bicgstab
inversion_info minv_vector_bicgstab(double  *phi, double  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(double*,double*,void*), void* extra_info)
{
// BICGSTAB solutions to Mphi = b 
//  see https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method

  // Initialize vectors.
  double *r, *r0, *v, *p, *s, *t; 
  double rho, rhoNew, alpha, beta, omega, tmp; 
  double rsq, ssq, bsqrt, truersq, ts, tt; 
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
  rsq = 0.0; ssq = 0.0; bsqrt = 0.0; truersq = 0.0;

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
  double rsq, ssq, bsqrt, truersq, tt; 
  int k,i;
  inversion_info invif;

  // Allocate memory.
  r = new complex<double>[size];
  r0 = new complex<double>[size];
  v = new complex<double>[size];
  p = new complex<double>[size];
  s = new complex<double>[size];
  t = new complex<double>[size];
  //r = (double*)malloc(size*sizeof(double));
  //r0 = (double*)malloc(size*sizeof(double));
  //v = (double*)malloc(size*sizeof(double));
  //p = (double*)malloc(size*sizeof(double));
  //s = (double*)malloc(size*sizeof(double));
  //t = (double*)malloc(size*sizeof(double));

  // Initialize values.
  rsq = 0.0; ssq = 0.0; bsqrt = 0.0; truersq = 0.0;

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
  for(i=0; i < size; i++) truersq += pow(abs(v[i] - phi0[i]),2);
  
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

// Solves lhs = A^(-1) rhs using gmres

// At one point during development, we solved the
// subspace normal equations within gmres with CG. 
// These functions and structs support this.
// At this point, we just use Gaussian elimination,
// but these are here for legacy.
/*
struct gmres_struct
{
  double** h;
  int size;
};
void gmres_hTh(double* phi, double* phi0, void* extra_data);
*/

inversion_info minv_vector_gmres_norestart(double  *phi, double  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(double*,double*,void*), void* extra_info)
{
// GMRES solutions to Mphi = b 
//  see https://en.wikipedia.org/wiki/Generalized_minimal_residual_method
// A clearer article: 
// http://www.math.iit.edu/~fass/477577_Chapter_14.pdf
  
  // Solve a counting issue...
  max_iter++;
  
  // Initialize vectors.
  double **q, **h; // Q matrix, Hessenberg matrix.
  double **hTh; // stores normal eqn matrix for hTh. 
  double *res, *tmp, *tmp2;
  double *y, *bhTy; // Holds the reduced solution. 
  double bres, localres; // norm of b, curr residual.
  double beta; // beta is the norm of res = b-Ax_0.
  //gmres_struct gmstr; // Passed to solve lls problem, no longer needed: we use gaussian elim.
  int i, j, k, iter;
  inversion_info invif;
  
  // Prepare q, h. These get larger on each iteration.
  q = new double*[max_iter];
  h = new double*[max_iter-1];
  hTh = new double*[max_iter-1];
  //q = (double**)malloc(max_iter*sizeof(double*));
  //h = (double**)malloc((max_iter-1)*sizeof(double*));
  //hTh = (double**)malloc((max_iter-1)*sizeof(double*));
  //gmstr.h = h; // Depreciated, no longer use CG for normal subspace.
  
  // Prepare space for the residual, and a temporary vector.
  res = new double[size];
  tmp = new double[size];
  tmp2 = new double[size];
  y = new double[max_iter];
  bhTy = new double [max_iter];
  //res = (double*)malloc(size*sizeof(double));
  //tmp = (double*)malloc(size*sizeof(double));
  //tmp2 = (double*)malloc(size*sizeof(double));
  //y = (double*)malloc(max_iter*sizeof(double));
  //bhTy = (double*)malloc(max_iter*sizeof(double));
  
  // Compute the residual. Ax' = res gives the solution phi+x'.
  (*matrix_vector)(tmp, phi, extra_info);
  for (i=0;i<size;i++) // res = b - Ax_0
  {
    res[i] = phi0[i]-tmp[i]; 
  }
  
  // Initialize the first q vector. q is res/||res||
  q[0] = new double[size];
  //q[0] = (double*)malloc(size*sizeof(double));
  
  beta = 0.0; bres = 0.0;
  for (i=0;i<size;i++)
  {
    beta = beta + res[i]*res[i];
    bres = bres + phi0[i]*phi0[i];
  }
  beta = sqrt(beta);
  bres = sqrt(bres);
  
  for (i=0;i<size;i++)
  {
    q[0][i] = res[i]/beta;
    //printf("%.8f ", q[0][i]);
  }
  //printf("\n");
  
  
  // Begin Arnoldi iterations. Taken from
  // https://en.wikipedia.org/wiki/Arnoldi_iteration
  for (iter = 1; iter < max_iter; iter++)
  {
    // Allocate a new q, h.
    q[iter] = new double[size];
    //q[iter] = (double*)malloc(size*sizeof(double));
    // Column 'i' of the Hessenberg matrix has
    // 'i+1' non-zero components.
    // To be fancier, we could do malloc((iter+1)...),
    // but then we have to do book-keeping everywhere.
    // That's why I kept having problems.
    h[iter-1] = new double[size];
    hTh[iter-1] = new double[size];
    //h[iter-1] = (double*)malloc(size*sizeof(double));
    //hTh[iter-1] = (double*)malloc(size*sizeof(double));
    for (j=0;j<size;j++)
    {
      q[iter][j] = h[iter-1][j] = hTh[iter-1][j] = 0.0;
    }
    
    // Compute q_{iter} = Aq_{iter-1}
    (*matrix_vector)(q[iter],q[iter-1],extra_info);
    
    // Perform an Arnoldi iteration. In some regards, this is just
    // a Gram-Schmidt process.
    for (j=0;j<iter;j++)
    {
      // Compute q[j] dot q[iter].
      h[iter-1][j] = 0.0;
      for (i=0;i<size;i++)
      {
        h[iter-1][j] = h[iter-1][j] + q[j][i]*q[iter][i]; //q* q
      }
      
      // Subtract off part of q[iter] along q[j].
      for (i=0;i<size;i++)
      {
        q[iter][i] = q[iter][i] - h[iter-1][j]*q[j][i];
      }
    } // Go to next existing q vector.
    
    // Normalize the new q[iter].
    h[iter-1][iter] = 0.0;
    for(i=0;i<size;i++)
    {
      h[iter-1][iter] = h[iter-1][iter] + q[iter][i]*q[iter][i];
    }
    h[iter-1][iter] = sqrt(h[iter-1][iter]);
    
    for(i=0;i<size;i++)
    {
      q[iter][i] = q[iter][i]/h[iter-1][iter];
      //printf("%.8f ", q[iter][i]);
    }
    //printf("\n");
    
    // This is formally the end of the Arnoldi iteration.
    // Right now, H is an (iter+1)x(iter) matrix stored
    // as its transpose.
    // h is iter x iter+1
    // For GMRES, we minimize || H y - beta e_1 ||_2 for y.
    // This becomes the normal equation: H^T H y = beta H^T e_1.
    // We solve this with gaussian elimination, but QR is apparently
    // a better way to do this.
    
    // Compute bhTy = beta H^T e_1, or the first column of H^T.
    //printf("bhTy: ");
    for(i=0;i<iter;i++)
    {
      bhTy[i] = beta*h[i][0];
      //printf("%.8f ", bhTy[i]);
    }
    //printf("\n");
    
    // Need to compute hTh, which is the iter x iter
    // normal matrix.
    // Recall h is iter x (iter+1)
    for (i=0;i<iter;i++)
    {
      for (j=0;j<iter;j++)
      {
        hTh[i][j] = 0.0;
        for (k=0;k<iter+1;k++)
        {
          hTh[i][j] = hTh[i][j] + h[i][k]*h[j][k];
        }
        //printf("%.8f ", hTh[i][j]);
      }
      //printf("\n");
    }
    gaussian_elimination(y, bhTy, hTh, iter);
    
    //printf("Exit gaussian elimination.\n"); fflush(stdout);
    
    // Alternatively, we can solve it with CG.
    // This is bad because we want an "exact" solution, which we
    // can only get in exact arithmetic. In practice, I didn't
    // get any benefit from using CG since the matrix is
    // so dense.
    
    //gmstr.size = iter;
    //minv_vector_cg(y, bhTy, iter, iter, eps*0.001, &gmres_hTh, (void*)(&gmstr));
    
    //double* meh = (double*)malloc(iter*sizeof(double));
    //gmres_hTh(meh, y, (void*)(&gmstr));
    //printf("meh: ");
    //for(i=0;i<iter;i++)
    //{
    //  printf("%.8f ", meh[i]);
    //}
    //printf("\n");
    
    // Use the q's to prolong y to the solution, put in tmp.
    for (i=0;i<size;i++)
    {
      tmp[i] = 0.0;
    }
    
    for (j=0;j<iter;j++)
    {
      for (i=0;i<size;i++)
      {
        tmp[i] = tmp[i]+q[j][i]*y[j];
      }
    }
    
    for (i=0;i<size;i++)
    {
      tmp2[i] = phi[i]+tmp[i];
    }
    
    // Check the residual. In principle, we can do this by looking at
    // the subspace norm, but I'm not sure how to combine that with
    // the fact we have an initial guess. Have to think about that. 
    (*matrix_vector)(res,tmp2,extra_info);
    
    localres = 0.0;
    for (i=0;i<size;i++)
    {
      localres = localres+(phi0[i]-res[i])*(phi0[i]-res[i]);
    }
    localres=sqrt(localres);
    
    if (localres < eps*bres)
    {
      break;
    }
    
    //printf("GMRES: Iter %d Resid %.8e\n", iter, localres/bres);
    
    
  }
  
  for (i=0;i<size;i++)
  {
    phi[i] = tmp2[i];
  }
  
  if(iter == max_iter) {
    //printf("GMRES: Failed to converge iter = %d, rsq = %e\n", iter,localres);
    invif.success = false;
    //return 0;// Failed convergence 
  }
  else
  {
     invif.success = true;
     //printf("GMRES: Converged in %d iterations.\n", iter);
  }
  
  // Check q's. Passes.
  /*
  for (i=0;i<iter;i++)
  {
    for (j=0;j<iter;j++)
    {
      bres = 0.0;
      for (k=0;k<size;k++)
      {
        bres = bres + q[i][k]*q[j][k];
      }
      printf("%.8f ", bres);
    }
    printf("\n");
  }
  */
  
  // Check A Q_n = Q_{n+1} \tilde H_n.
  // This requires a size x (iter-1) space.
  // We can check this column by column (of size 'size').
  // Works!
  /*for (i=0;i<iter-1;i++)
  {
    // Column of A Q_n
    (*matrix_vector)(res, q[i], extra_info);
    
    for (j=0;j<size;j++)
    {
      for (k=0;k<iter;k++)
      {
        res[j] = res[j] - q[k][j]*h[i][k];
      }
    }
    
    bres = 0.0;
    // Find magnitude of res.
    for (j=0;j<size;j++)
    {
      bres = bres + res[j]*res[j];
    }
    printf("Mag column %d: %.8f\n", i, sqrt(bres));
  }*/
    
    
  
  
  // Clean up.
  for (i=0;i<iter-1;i++)
  {
    delete[] h[i];
    delete[] q[i];
    delete[] hTh[i];
    //free(h[i]);
    //free(q[i]);
    //free(hTh[i]);
  }
  delete[] q[iter-1];
  delete[] h;
  delete[] q;
  delete[] res;
  delete[] tmp;
  delete[] tmp2;
  delete[] y;
  delete[] bhTy;
  //free(q[iter-1]);
  //free(h);
  //free(q);
  //free(res);
  //free(tmp);
  //free(tmp2);
  //free(y);
  //free(bhTy);
  
  invif.resSq = localres*localres;
  invif.iter = iter;
  if (invif.success == false) { invif.iter--; } // For loop has an extra incr. at the end.
  invif.name = "GMRES";
  return invif; // Convergence 

} 


// Performs GMRES with restarts when restart_freq is hit.
// This may be sloppy, but it works.
inversion_info minv_vector_gmres_restart(double  *phi, double  *phi0, int size, int max_iter, double res, int restart_freq, void (*matrix_vector)(double*,double*,void*), void* extra_info)
{
  int iter; // counts total number of iterations.
  inversion_info invif;

  iter = 0;  
  do
  {
    invif = minv_vector_gmres_norestart(phi, phi0, size, restart_freq, res, matrix_vector, extra_info);
    iter += invif.iter;
  }
  while (iter < max_iter && invif.success == false && sqrt(invif.resSq) > res);
  
  invif.iter = iter;
  stringstream ss;
  ss << "GMRES(" << restart_freq << ")";
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

inversion_info minv_vector_gmres_norestart(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info)
{
// GMRES solutions to Mphi = b 
//  see https://en.wikipedia.org/wiki/Generalized_minimal_residual_method
// A clearer article: 
// http://www.math.iit.edu/~fass/477577_Chapter_14.pdf
  
  // Solve a counting issue...
  max_iter++;
  
  // Initialize vectors.
  complex<double> **q, **h; // Q matrix, Hessenberg matrix.
  complex<double> **hTh; // stores normal eqn matrix for hTh. 
  complex<double> *res, *tmp, *tmp2;
  complex<double> *y, *bhTy; // Holds the reduced solution. 
  double bres, localres; // norm of b, curr residual.
  double beta; // beta is the norm of res = b-Ax_0.
  //gmres_struct gmstr; // Passed to solve lls problem, no longer needed: we use gaussian elim.
  int i, j, k, iter;
  inversion_info invif;
  
  // Prepare q, h. These get larger on each iteration.
  q = new complex<double>*[max_iter];
  h = new complex<double>*[max_iter-1];
  hTh = new complex<double>*[max_iter-1];
  //q = (double**)malloc(max_iter*sizeof(double*));
  //h = (double**)malloc((max_iter-1)*sizeof(double*));
  //hTh = (double**)malloc((max_iter-1)*sizeof(double*));
  //gmstr.h = h; // Depreciated, no longer use CG for normal subspace.
  
  // Prepare space for the residual, and a temporary vector.
  res = new complex<double>[size];
  tmp = new complex<double>[size];
  tmp2 = new complex<double>[size];
  y = new complex<double>[max_iter];
  bhTy = new complex<double>[max_iter];
  //res = (double*)malloc(size*sizeof(double));
  //tmp = (double*)malloc(size*sizeof(double));
  //tmp2 = (double*)malloc(size*sizeof(double));
  //y = (double*)malloc(max_iter*sizeof(double));
  //bhTy = (double*)malloc(max_iter*sizeof(double));
  
  // Compute the residual. Ax' = res gives the solution phi+x'.
  (*matrix_vector)(tmp, phi, extra_info);
  for (i=0;i<size;i++) // res = b - Ax_0
  {
    res[i] = phi0[i]-tmp[i]; 
  }
  
  // Initialize the first q vector. q is res/||res||
  q[0] = new complex<double>[size];
  //q[0] = (double*)malloc(size*sizeof(double));
  
  beta = sqrt(norm2sq<double>(res, size));
  bres = sqrt(norm2sq<double>(phi0, size));
  
  for (i=0;i<size;i++)
  {
    q[0][i] = res[i]/beta;
    //printf("%.8f ", q[0][i]);
  }
  //printf("\n");
  
  
  // Begin Arnoldi iterations. Taken from
  // https://en.wikipedia.org/wiki/Arnoldi_iteration
  for (iter = 1; iter < max_iter; iter++)
  {
    // Allocate a new q, h.
    q[iter] = new complex<double>[size];
    //q[iter] = (double*)malloc(size*sizeof(double));
    // Column 'i' of the Hessenberg matrix has
    // 'i+1' non-zero components.
    // To be fancier, we could do malloc((iter+1)...),
    // but then we have to do book-keeping everywhere.
    // That's why I kept having problems.
    h[iter-1] = new complex<double>[size];
    hTh[iter-1] = new complex<double>[size];
    //h[iter-1] = (double*)malloc(size*sizeof(double));
    //hTh[iter-1] = (double*)malloc(size*sizeof(double));
    for (j=0;j<size;j++)
    {
      q[iter][j] = h[iter-1][j] = hTh[iter-1][j] = 0.0;
    }
    
    // Compute q_{iter} = Aq_{iter-1}
    (*matrix_vector)(q[iter],q[iter-1],extra_info);
    
    // Perform an Arnoldi iteration. In some regards, this is just
    // a Gram-Schmidt process.
    for (j=0;j<iter;j++)
    {
      // Compute q[j] dot q[iter].
      h[iter-1][j] = dot<double>(q[j],q[iter], size);
      
      // Subtract off part of q[iter] along q[j].
      for (i=0;i<size;i++)
      {
        q[iter][i] = q[iter][i] - h[iter-1][j]*q[j][i];
      }
    } // Go to next existing q vector.
    
    // Normalize the new q[iter].
    h[iter-1][iter] = sqrt(norm2sq<double>(q[iter], size));
    for(i=0;i<size;i++)
    {
      q[iter][i] = q[iter][i]/h[iter-1][iter];
      //printf("%.8f ", q[iter][i]);
    }
    //printf("\n");
    
    // This is formally the end of the Arnoldi iteration.
    // Right now, H is an (iter+1)x(iter) matrix stored
    // as its transpose.
    // h is iter x iter+1
    // For GMRES, we minimize || H y - beta e_1 ||_2 for y.
    // This becomes the normal equation: H^T H y = beta H^T e_1.
    // We solve this with gaussian elimination, but QR is apparently
    // a better way to do this.
    
    // Compute bhTy = beta H^T e_1, or the first column of H^T.
    //printf("bhTy: ");
    for(i=0;i<iter;i++)
    {
      bhTy[i] = beta*h[i][0];
      //printf("%.8f ", bhTy[i]);
    }
    //printf("\n");
    
    // Need to compute hTh, which is the iter x iter
    // normal matrix.
    // Recall h is iter x (iter+1)
    for (i=0;i<iter;i++)
    {
      for (j=0;j<iter;j++)
      {
        hTh[i][j] = dot<double>(h[i],h[j], iter+1);
        //printf("%.8f ", hTh[i][j]);
      }
      //printf("\n");
    }
    gaussian_elimination(y, bhTy, hTh, iter);
    
    //printf("Exit gaussian elimination.\n"); fflush(stdout);
    
    // Alternatively, we can solve it with CG.
    // This is bad because we want an "exact" solution, which we
    // can only get in exact arithmetic. In practice, I didn't
    // get any benefit from using CG since the matrix is
    // so dense.
    
    //gmstr.size = iter;
    //minv_vector_cg(y, bhTy, iter, iter, eps*0.001, &gmres_hTh, (void*)(&gmstr));
    
    //double* meh = (double*)malloc(iter*sizeof(double));
    //gmres_hTh(meh, y, (void*)(&gmstr));
    //printf("meh: ");
    //for(i=0;i<iter;i++)
    //{
    //  printf("%.8f ", meh[i]);
    //}
    //printf("\n");
    
    // Use the q's to prolong y to the solution, put in tmp.
    for (i=0;i<size;i++)
    {
      tmp[i] = 0.0;
    }
    
    for (j=0;j<iter;j++)
    {
      for (i=0;i<size;i++)
      {
        tmp[i] = tmp[i]+q[j][i]*y[j];
      }
    }
    
    for (i=0;i<size;i++)
    {
      tmp2[i] = phi[i]+tmp[i];
    }
    
    // Check the residual. In principle, we can do this by looking at
    // the subspace norm, but I'm not sure how to combine that with
    // the fact we have an initial guess. Have to think about that. 
    (*matrix_vector)(res,tmp2,extra_info);
    
    localres = 0.0;
    for (i=0;i<size;i++)
    {
      localres = localres+pow(abs(phi0[i]-res[i]),2);
    }
    localres=sqrt(localres);
    
    if (localres < eps*bres)
    {
      break;
    }
    
    //printf("GMRES: Iter %d Resid %.8e\n", iter, localres/bres);
    
    
  }
  
  for (i=0;i<size;i++)
  {
    phi[i] = tmp2[i];
  }
  
  if(iter == max_iter) {
    //printf("GMRES: Failed to converge iter = %d, rsq = %e\n", iter,localres);
    invif.success = false;
    //return 0;// Failed convergence 
  }
  else
  {
     invif.success = true;
     //printf("GMRES: Converged in %d iterations.\n", iter);
  }
  
  // Check q's. Passes.
  /*
  complex<double> tmp3;
  for (i=0;i<iter;i++)
  {
    for (j=0;j<iter;j++)
    {
      tmp3 = 0.0;
      for (k=0;k<size;k++)
      {
        tmp3 = tmp3 + conj(q[i][k])*q[j][k];
      }
      cout << tmp3 << " ";
    }
    cout << "\n";
  }*/
  
  
  // Check A Q_n = Q_{n+1} \tilde H_n.
  // This requires a size x (iter-1) space.
  // We can check this column by column (of size 'size').
  // Works!
  /*
  for (i=0;i<iter-1;i++)
  {
    // Column of A Q_n
    (*matrix_vector)(res, q[i], extra_info);
    
    for (j=0;j<size;j++)
    {
      for (k=0;k<iter;k++)
      {
        res[j] = res[j] - q[k][j]*h[i][k];
      }
    }
    
    complex<double> tmp3 = 0.0;
    // Find magnitude of res.
    for (j=0;j<size;j++)
    {
      tmp3 = tmp3 + conj(res[j])*res[j];
    }
    cout << "Mag column " << i << ": " << sqrt(tmp3) << "\n";
  }
  */
    
  
  
  // Clean up.
  for (i=0;i<iter-1;i++)
  {
    delete[] h[i];
    delete[] q[i];
    delete[] hTh[i];
    //free(h[i]);
    //free(q[i]);
    //free(hTh[i]);
  }
  delete[] q[iter-1];
  delete[] h;
  delete[] q;
  delete[] res;
  delete[] tmp;
  delete[] tmp2;
  delete[] y;
  delete[] bhTy;
  //free(q[iter-1]);
  //free(h);
  //free(q);
  //free(res);
  //free(tmp);
  //free(tmp2);
  //free(y);
  //free(bhTy);
  
  invif.resSq = localres*localres;
  invif.iter = iter;
  if (invif.success == false) { invif.iter--; } // For loop has an extra incr. at the end.
  invif.name = "GMRES";
  return invif; // Convergence 

} 


// Performs GMRES with restarts when restart_freq is hit.
// This may be sloppy, but it works.
inversion_info minv_vector_gmres_restart(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double res, int restart_freq, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info)
{
  int iter; // counts total number of iterations.
  inversion_info invif;

  iter = 0;  
  do
  {
    invif = minv_vector_gmres_norestart(phi, phi0, size, restart_freq, res, matrix_vector, extra_info);
    iter += invif.iter;
  }
  while (iter < max_iter && invif.success == false && sqrt(invif.resSq) > res);
  
  invif.iter = iter;
  stringstream ss;
  ss << "GMRES(" << restart_freq << ")";
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


// This is no longer necessary, we now use
// explicit gaussian elimination instead of CG.
// Solves H^T H phi = phi0. 
// extra_data is really a gmres_struct
// containing h = H^T which is n by (n+1),
// and the size.
/*void gmres_hTh(double* phi, double* phi0, void* extra_data)
{
  gmres_struct gmstr = *(gmres_struct*)(extra_data);
  double* yinter; // Holds the intermediate of H^T H.
  int i,j;
  
  yinter = (double*)malloc((gmstr.size+1)*sizeof(double));
  for (i=0;i<gmstr.size+1;i++)
  {
    yinter[i] = 0.0;
  }
  for (i=0;i<gmstr.size;i++)
  {
    phi[i] = 0.0;
  }
  
  // Do the H phi0 multiply.
  for (i=0;i<gmstr.size+1;i++)
  {
    for (j=0;j<gmstr.size;j++)
    {
      yinter[i] = yinter[i] + gmstr.h[j][i]*phi0[j];
    }
  }
  
  // Do the H^T yinter multiply.
  for (i=0;i<gmstr.size;i++)
  {
    for (j=0;j<gmstr.size+1;j++)
    {
      phi[i] = phi[i] + gmstr.h[i][j]*yinter[j];
    }
  }
  
  free(yinter);
}*/


// Perform gaussian elimination on a matrix to solve Ax=b
void gaussian_elimination(double* x, double* b, double** matrix, int size)
{
  // Initialize variables.
  int i,j,k;
  double** grown_matrix;
  double pivot_space; // For pivot.
  double max_val = 0.0;
  int max_index = -1;
  
  // Declare size of matrix.
  grown_matrix = new double*[size];
  for (i=0;i<size;i++)
  {
    grown_matrix[i] = new double[size+1];
  }
  
  // Copy things over.
  for (i=0;i<size;i++)
  {
    for (j=0;j<size;j++)
    {
      grown_matrix[i][j] = matrix[i][j];
    }
    grown_matrix[i][size] = b[i];
  }
  
  // Cool! Begin gaussian elimination.
  // Iterate over all rows.
  for (i=0;i<size;i++)
  {
    max_val = 0.0; max_index = -1;
    // First, pivot over self and all lower rows.
    for (j=i;j<size;j++)
    {
      if (abs(grown_matrix[j][i]) > max_val)
      {
        max_index = j;
        max_val = abs(grown_matrix[j][i]);
      }
    }
    
    //printf("Max: %d %.8f\n", max_index, max_val); fflush(stdout);
    if (max_index == -1) // everything's 0!
    {
      return;
    }
    
    // Put maximal row in pivot location.
    // Switch i, max_index.
    if (max_index != i)
    {
      for(j=i;j<size+1;j++)
      {
        pivot_space = grown_matrix[i][j];
        grown_matrix[i][j] = grown_matrix[max_index][j];
        grown_matrix[max_index][j] = pivot_space;
      }
      //printf("Pivoted.\n"); fflush(stdout);
    }
    
    
    // Good! We've safely pivoted. Now, normalize the top row.
    for (j=i+1;j<size+1;j++)
    {
      grown_matrix[i][j] = grown_matrix[i][j]/grown_matrix[i][i];
    }
    grown_matrix[i][i] = 1.0;
    
    // Eliminate the top row from all other rows.
    // This part can get parallelized!
    for (j=0;j<size;j++)
    {
      if (j == i) continue;
      for (k=i+1;k<size+1;k++)
      {
        grown_matrix[j][k] = grown_matrix[j][k] - grown_matrix[j][i]*grown_matrix[i][k];
      }
      grown_matrix[j][i] = 0.0;
    }
    
  } // end loop over rows = i.
  //printf("Exited reduction.\n"); fflush(stdout);
  
  // Copy result in!
  for (i=0;i<size;i++)
  {
    x[i] = grown_matrix[i][size];
  }
  
  //printf("Copied solution.\n"); fflush(stdout);
  
  // Free
  for (i=0;i<size;i++)
  {
    delete[] grown_matrix[i];
  }
  delete[] grown_matrix;
  
}

void gaussian_elimination(complex<double>* x, complex<double>* b, complex<double>** matrix, int size)
{
  // Initialize variables.
  int i,j,k;
  complex<double>** grown_matrix;
  complex<double> pivot_space; // For pivot.
  double max_val = 0.0;
  int max_index = -1;
  
  // Declare size of matrix.
  grown_matrix = new complex<double>*[size];
  for (i=0;i<size;i++)
  {
    grown_matrix[i] = new complex<double>[size+1];
  }
  
  // Copy things over.
  for (i=0;i<size;i++)
  {
    for (j=0;j<size;j++)
    {
      grown_matrix[i][j] = matrix[i][j];
    }
    grown_matrix[i][size] = b[i];
  }
  
  // Cool! Begin gaussian elimination.
  // Iterate over all rows.
  for (i=0;i<size;i++)
  {
    max_val = 0.0; max_index = -1;
    // First, pivot over self and all lower rows.
    for (j=i;j<size;j++)
    {
      if (abs(grown_matrix[j][i]) > max_val)
      {
        max_index = j;
        max_val = abs(grown_matrix[j][i]);
      }
    }
    
    //printf("Max: %d %.8f\n", max_index, max_val); fflush(stdout);
    if (max_index == -1) // everything's 0!
    {
      return;
    }
    
    // Put maximal row in pivot location.
    // Switch i, max_index.
    if (max_index != i)
    {
      for(j=i;j<size+1;j++)
      {
        pivot_space = grown_matrix[i][j];
        grown_matrix[i][j] = grown_matrix[max_index][j];
        grown_matrix[max_index][j] = pivot_space;
      }
      //printf("Pivoted.\n"); fflush(stdout);
    }
    
    
    // Good! We've safely pivoted. Now, normalize the top row.
    for (j=i+1;j<size+1;j++)
    {
      grown_matrix[i][j] = grown_matrix[i][j]/grown_matrix[i][i];
    }
    grown_matrix[i][i] = 1.0;
    
    // Eliminate the top row from all other rows.
    // This part can get parallelized!
    for (j=0;j<size;j++)
    {
      if (j == i) continue;
      for (k=i+1;k<size+1;k++)
      {
        grown_matrix[j][k] = grown_matrix[j][k] - grown_matrix[j][i]*grown_matrix[i][k];
      }
      grown_matrix[j][i] = 0.0;
    }
    
  } // end loop over rows = i.
  //printf("Exited reduction.\n"); fflush(stdout);
  
  // Copy result in!
  for (i=0;i<size;i++)
  {
    x[i] = grown_matrix[i][size];
  }
  
  //printf("Copied solution.\n"); fflush(stdout);
  
  // Free
  for (i=0;i<size;i++)
  {
    free(grown_matrix[i]);
  }
  free(grown_matrix);
  
}


// This code tests the gaussian elimination routine.
/*
void test_gaussian()
{
   int m,n;
   int lsize = 20;
   double** test;
   double* b;
   double* x;
   double* b2;
   
   test = (double**)malloc(lsize*sizeof(double*));
   for (m=0;m<lsize;m++)
   {
     test[m] = (double*)malloc(lsize*sizeof(double));
   }
   b = (double*)malloc(lsize*sizeof(double));
   b2 = (double*)malloc(lsize*sizeof(double));
   x = (double*)malloc(lsize*sizeof(double));
   
   for (m=0;m<lsize;m++)
   {
     test[m][m] = 3;
     test[m][(m+1)%lsize] = -1;
     test[m][(m+3)%lsize] = -6;
     test[m][(m-1+lsize)%lsize] = -1;
     test[m][(m-3+lsize)%lsize] = -6;
   }
   
   b[0] = 1;
   gaussian_elimination(x,b,test,lsize);
   
   for (m=0;m<lsize;m++)
   {
     printf("%f ", x[m]);
   }
   printf("\n");
   
   for (m=0;m<lsize;m++)
   {
     b2[m] = 0.0;
     for (n=0;n<lsize;n++)
     {
       b2[m] = b2[m] + test[m][n]*x[n];
     }
   }
   
   for (m=0;m<lsize;m++)
   {
     printf("%f ", b[m]-b2[m]);
   }
   printf("\n");
   

   return;
}
*/


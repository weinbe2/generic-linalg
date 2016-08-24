
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>

using std::complex;

#include "arpack_interface.h"

// Various maxima. This is used to allocate more than
// enough memory.

// Maximum dimension of the matrix.
#define MAXN 1024

// Maximum number of eigenvalues requested.
#define MAXNEV 256

// Maximum number of basis vectors for Arnoldi iterations.
// Remember: MAXNEV+2 < MAXNCV < MAXN
#define MAXNCV 512

// Unique to this problem: define size of 1D laplace chain.
#define NDIM 200
#define MASS 0.1
#define PI 3.14159265358979323846
typedef struct {
   double ndim;
   double mass;
} laplace_1d_t;
   

// This function does the multiply!
void mat_vec(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Solves lhs = (A-sigma*I)^(-1) rhs with bicgstab
double minv_vector_bicgstab_shift(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double res, void (*matrix_vector)(complex<double>*,complex<double>*,void*), complex<double> sigma, void* extra_info);

int main(int argc, char** argv)
{
   // Declare variables.
   int i, j; // Iterators.
   int n, nev, ncv; // matrix dimension, eigenvals, internal vals.
   int maxitr; // maximum number of iterations.
   char* which; // Sets which part of spectrum to get.
   double tol; // Set tolerance. 0 = machine.
   complex<double> sigma; // Set to zero for now.
   laplace_1d_t sys_info; // Set the mass term and dim of matrix.
   arpack_solve_t info_solve; // Hold info on the solve.
   complex<double> *evec, *eval; // Will hold eigenvalues, eigenvectors.
   complex<double> *resid; // Hold residuals of eigenvalues.
   complex<double> *ax, *r; // Temporary space to hold vectors to compute
                   // residuals.
   int maxitr_cg; // Max number of iterations for cg.
   double tol_cg; // Tolerance for cg.
                   
   arpack_dcn_t* ar_strc; // Holds arpack info.
   
   // Set dimension of matrix, number of eigenvalues, etc.
   n = NDIM;
   nev = 20; // nev < n
   ncv = 40; // nev+2 < ncv < n.
   sys_info.mass = MASS;
   sys_info.ndim = NDIM;
   
   // Configure calculating eigenvectors.
   maxitr = 4000; // 4000 iterations!
   which = "SM"; // This could be SM (smallest magnitude).
                 // This could also be LM (largest magnitude)
                 // SR (smallest real), SI (smallest imaginary),
                 // and similar for largest.
   tol = 1e-7; // 0 = machine precision.
   sigma = MASS*0.99; // Zero for now.

   // Allocate space.
   evec = new complex<double>[n*nev]; // eigenvals
   eval = new complex<double>[nev]; // eigenvecs
   resid = new complex<double>[nev]; // residuals
   ax = new complex<double>[n]; // temp vector.
   r = new complex<double>[n]; // temp vector.
   
   printf("Allocated.\n"); fflush(stdout);
   
   // Initialize arpack. Largely just allocates internal memory.
   ar_strc = arpack_dcn_init(MAXN, MAXNEV, MAXNCV);
   
   printf("Init arpack.\n"); fflush(stdout);
   
   // Get some eigenvalues and vectors!
   info_solve = arpack_dcn_getev(ar_strc, eval, evec, n, nev, ncv,
                              maxitr, which, tol, sigma,
                              &mat_vec, (void*)(&sys_info));
   
   printf("Got evals/evecs.\n"); fflush(stdout);         
   
   if (info_solve.is_error == 0 && info_solve.nconv > 0)
   {
      // Print info about the eigensolve.
      printf("Number of converged eigenvalues: %d\n", info_solve.nconv);
      printf("Number of iteration steps:       %d\n", info_solve.niter);
      printf("Number of matrix multiplies:     %d\n", info_solve.nops);
      
      // Get the residual norm:
      // || A*x - lambda*x ||
      // for the nconv accurately computed eigenvals/vecs.
      // Put the residual here:
      
      // Loop over all good eigenvalues.
      for (i=0;i<info_solve.nconv;i++)
      {
         // Get pointer to start of i'th eigenvector.
         complex<double> *evec_local = (evec+i*n);
         
         // If you want to print it out, uncomment this.
         //for (int j=0;j<NDIM;j++)
         //{
         //   printf("%d\t%.8e\n", i, evec_local[j]);
         //}
         
         // Get A*x, where x is an eigenvector.
         mat_vec((complex<double>*)ax, evec_local, (void*)(&sys_info));
         
         resid[i] = 0;
         

         // Get the norm squared.
         for (j=0;j<NDIM;j++)
         {
			 resid[i] += real(conj(ax[j]-eval[i]*evec_local[j])*(ax[j]-eval[i]*evec_local[j]));
         }
         resid[i] = sqrt(resid[i])/(abs(eval[i])*abs(eval[i]));
      }
      
      printf("Got residuals.\n"); fflush(stdout);
   
      // Print eigenvalues, residual.
      // Also, print what mode we get for 1D Laplace
      for (i=0;i<info_solve.nconv;i++)
      {
         // We compute the mode. This can overflow on some
         // ops, though, so we're careful of that.
         double tmp,tmp2; 
         tmp = (real(eval[i])-MASS)/4; 
         tmp2 = (tmp > 0 ? sqrt(tmp) : 0); // make sure tmp > 0
         tmp = (tmp2 < 1 ? (tmp2 > -1 ? asin(tmp2) : asin(-1.0)) : asin(1.0))*NDIM/PI; // make sure -1 < tmp2 < 1
         // tmp = NDIM/PI*asin(sqrt((d[i]-MASS)/4));
         
         printf("%d\t%.8e\t%.8e\t%.8e\t%.8e\n", i+1, real(eval[i]), imag(eval[i]), real(resid[i]), tmp);
      }
   
   }
   
   // Now redo it with shift invert!
   which = "LM"; // This could be SM (smallest magnitude).
                 // This could also be LM (largest magnitude)
                 // SR (smallest real), SI (smallest imaginary),
                 // and similar for largest.
   maxitr_cg = 4000;
   tol_cg = 1e-9;
   
   // Get some eigenvalues and vectors!
   info_solve = arpack_dcn_getev_sinv(ar_strc, eval, evec, n, nev, ncv,
                              maxitr, which, tol, sigma,
                              &mat_vec, &minv_vector_bicgstab_shift, maxitr_cg, tol_cg, (void*)(&sys_info));
   
   printf("Got evals/evecs.\n"); fflush(stdout);         
   
   if (info_solve.is_error == 0 && info_solve.nconv > 0)
   {
      // Print info about the eigensolve.
      printf("Number of converged eigenvalues: %d\n", info_solve.nconv);
      printf("Number of iteration steps:       %d\n", info_solve.niter);
      printf("Number of matrix multiplies:     %d\n", info_solve.nops);
      
      // Get the residual norm:
      // || A*x - lambda*x ||
      // for the nconv accurately computed eigenvals/vecs.
      // Put the residual here:
      
      // Loop over all good eigenvalues.
      for (i=0;i<info_solve.nconv;i++)
      {
         // Get pointer to start of i'th eigenvector.
         complex<double> *evec_local = (evec+i*n);
         
         // If you want to print it out, uncomment this.
         //for (int j=0;j<NDIM;j++)
         //{
         //   printf("%d\t%.8e\n", i, evec_local[j]);
         //}
         
         // Get A*x, where x is an eigenvector.
         mat_vec((complex<double>*)ax, evec_local, (void*)(&sys_info));
         
         resid[i] = 0;
         

         // Get the norm squared.
         for (j=0;j<NDIM;j++)
         {
			 resid[i] += real(conj(ax[j]-eval[i]*evec_local[j])*(ax[j]-eval[i]*evec_local[j]));
         }
         resid[i] = sqrt(resid[i])/(abs(eval[i])*abs(eval[i]));
      }
      
      printf("Got residuals.\n"); fflush(stdout);
   
      // Print eigenvalues, residual.
      // Also, print what mode we get for 1D Laplace
      for (i=0;i<info_solve.nconv;i++)
      {
         // We compute the mode. This can overflow on some
         // ops, though, so we're careful of that.
         double tmp,tmp2; 
         tmp = (real(eval[i])-MASS)/4; 
         tmp2 = (tmp > 0 ? sqrt(tmp) : 0); // make sure tmp > 0
         tmp = (tmp2 < 1 ? (tmp2 > -1 ? asin(tmp2) : asin(-1.0)) : asin(1.0))*NDIM/PI; // make sure -1 < tmp2 < 1
         // tmp = NDIM/PI*asin(sqrt((d[i]-MASS)/4));
         
         printf("%d\t%.8e\t%.8e\t%.8e\t%.8e\n", i+1, real(eval[i]), imag(eval[i]), real(resid[i]), tmp);
      }
   
   }
   
   // Clean up!
   delete[] evec;
	delete[] eval;
	delete[] resid;
	delete[] ax;
	delete[] r;
   arpack_dcn_free(&ar_strc);
   
}

// matrix multiplication function!
// Return lhs = A*rhs;

void mat_vec(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
   // Apply the 1D Laplace eqn with a mass!
   
   // Declare variables.
   int ndim, i;
   double mass;
   laplace_1d_t* system_info;
   
   // Get input data.
   system_info = (laplace_1d_t*)extra_data;
   mass = system_info->mass;
   ndim = system_info->ndim;
         
   // Non boundary case!
   for (i=1;i<(ndim-1);i++)
   {
      lhs[i] = (2+mass)*rhs[i]-rhs[i-1]-rhs[i+1];
   }
   // Boundary.
   lhs[0] = (2+mass)*rhs[0]-rhs[1]-rhs[ndim-1];
   lhs[ndim-1] = (2+mass)*rhs[ndim-1]-rhs[0]-rhs[ndim-2];
}

// Solves lhs = (A-sigma*I)^(-1) rhs using bicgstab
double minv_vector_bicgstab_shift(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(complex<double>*,complex<double>*,void*), complex<double> sigma, void* extra_info)
{
// BICGSTAB solutions to Mphi = b 
//  see https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method

  // Initialize vectors.
  complex<double> *r, *r0, *v, *p, *s, *t; 
  complex<double> rho, rhoNew, alpha, beta, omega, tmp, ts; 
  double rsq, ssq, bsqrt, truersq, tt; 
  int k,i;

  // Allocate memory.
  r = new complex<double>[size*sizeof(complex<double>)];
  r0 = new complex<double>[size*sizeof(complex<double>)];
  v = new complex<double>[size*sizeof(complex<double>)];
  p = new complex<double>[size*sizeof(complex<double>)];
  s = new complex<double>[size*sizeof(complex<double>)];
  t = new complex<double>[size*sizeof(complex<double>)];

  // Initialize values.
  rsq = 0.0; ssq = 0.0; bsqrt = 0.0; truersq = 0.0;

  // Take advantage of initial guess in phi.
  (*matrix_vector)(v, phi, extra_info);
  for(i = 0; i<size; i++) {
    r[i] = phi0[i] - v[i] + sigma*phi[i]; // 1. r0 = b-Ax0
    r0[i] = r[i]; // 2. Assign rhat0 = r0.
    rho = alpha = omega = 1.0; // 3. Assign initial values.
    v[i] = p[i] = 0.0; // 4. v0 = p0 = 0.
    bsqrt += real(conj(phi0[i])*phi0[i]); // Used to check if residual is small.
  }
  bsqrt = sqrt(bsqrt);

  // iterate till convergence
  for(k = 0; k< max_iter; k++) {
    // 5.1. rhoNew = <rhat0, ri-1>
    rhoNew = 0.0;
    for (i = 0; i < size; i++) { rhoNew += conj(r0[i])*r[i]; }
    
    // 5.2. beta = (rhoNew/rho)(alpha/omega_i-1)
    beta = (rhoNew/rho)*(alpha/omega);
    rho = rhoNew;

    // 5.3. p = r + beta(p - omega v)
    for (i = 0; i < size; i++) {
      p[i] = r[i] + beta*p[i] - beta*omega*v[i];
    }
    
    // 5.4. v = (A+sigma*I)p
    (*matrix_vector)(v, p, extra_info);
    for (i = 0; i < size; i++) {
      v[i] -= sigma*p[i];
    }
    
    // 5.5. alpha = rho/<rhat0, v>
    tmp = 0.0;
    for (i = 0; i < size; i++) {
      tmp += conj(r0[i])*v[i];
    }
    alpha = rho/tmp;
    
    // 5.6. s = r - alpha v
    for (i = 0; i < size; i++) {
      s[i] = r[i] - alpha*v[i];
    }
    
    // 5.7. If ||s|| is sufficiently small, x = x+alpha p, quit.
    ssq = 0.0;
    for (i = 0; i < size; i++) {
      ssq += real(conj(s[i])*s[i]);
    }
    if (sqrt(ssq) < eps*bsqrt)
    {
      // printf("Final rsq = %g\n", ssq);
		for (i = 0; i < size; i++)
		{
      		phi[i] = phi[i] + alpha*p[i];
		}
      break;
    }
    
    // 5.8. t = (A+sigma*I)s
    (*matrix_vector)(t, s, extra_info);
    for (i = 0; i < size; i++) {
      t[i] -= sigma*s[i];
    }
    
    // 4.9. omega = <t, s>/<t, t>;
    ts = tt = 0.0;
    for (i = 0; i < size; i++) {
      ts += conj(t[i])*s[i];
      tt += real(conj(t[i])*t[i]);
    }
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
    printf("CG: Failed to converge iter = %d, rsq = %e\n", k,rsq);
    return 0;// Failed convergence 
  }
  
  (*matrix_vector)(v,phi,extra_info);
  for(i=0; i < size; i++) truersq += real(conj(v[i] - sigma*phi[i] - phi0[i])*(v[i] - sigma*phi[i] - phi0[i]));
  
  // Free all the things!
  free(r);
  free(r0);
  free(v);
  free(p);
  free(s);
  free(t);


  //  printf("# CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  return truersq; // Convergence 
} 



